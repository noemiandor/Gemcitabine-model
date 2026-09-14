#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "Matrix", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05h_existing_neutrophil_states", "mouse_with_existing_neutrophil_states.rds")
out_dir <- file.path(results_root, "05i_myeloid_extraction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05h input: ", input_rds, call. = FALSE)

full_obj <- readRDS(input_rds)
if (!inherits(full_obj, "Seurat")) stop("05h input is not a Seurat object.", call. = FALSE)
required_meta <- c("sample", "barcode_raw", "mouse_final_cell_type", "mouse_fraction", "scDblFinder.class", "seurat_clusters")
missing_meta <- setdiff(required_meta, colnames(full_obj@meta.data))
if (length(missing_meta) > 0L) stop("05h object lacks metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)

source_labels <- as.character(unlist(cfg$myeloid$source_final_labels, use.names = FALSE))
selected <- as.character(full_obj$mouse_final_cell_type) %in% source_labels
selected_cells <- colnames(full_obj)[selected]
if (length(selected_cells) == 0L) stop("No mouse myeloid cells matched configured source labels.", call. = FALSE)
if (any(as.character(full_obj$scDblFinder.class[selected]) != "singlet")) stop("Selected myeloid compartment contains non-singlets.", call. = FALSE)
if (any(full_obj$mouse_fraction[selected] < as.numeric(cfg$species$mouse_fraction_min), na.rm = TRUE)) stop("Selected myeloid compartment violates mouse threshold.", call. = FALSE)
if (length(unique(as.character(full_obj$sample[selected]))) != length(cfg$sample_sources)) stop("Myeloid extraction does not retain all configured tissue samples.", call. = FALSE)

counts <- Seurat::GetAssayData(full_obj, assay = "RNA", slot = "counts")[, selected_cells, drop = FALSE]
if (sum(startsWith(rownames(counts), as.character(cfg$species$human_prefix))) != 0L) stop("Human features detected in myeloid counts.", call. = FALSE)
metadata <- full_obj@meta.data[selected_cells, , drop = FALSE]
metadata$mouse_myeloid_source_label <- as.character(metadata$mouse_final_cell_type)
metadata$mouse_parent_cluster_05d <- as.character(metadata$seurat_clusters)

# Build a clean RNA-only object. This prevents the full-object integrated graph,
# PCA, and UMAP from leaking into the independent myeloid re-analysis.
myeloid <- Seurat::CreateSeuratObject(
  counts = counts,
  project = "SUM159_mouse_myeloid",
  assay = "RNA",
  min.cells = 0,
  min.features = 0,
  meta.data = metadata
)
if (!identical(colnames(myeloid), selected_cells)) stop("Cell order changed while creating the clean myeloid object.", call. = FALSE)

selection <- data.frame(
  cell = colnames(full_obj), sample = as.character(full_obj$sample), barcode_raw = as.character(full_obj$barcode_raw),
  source_final_cell_type = as.character(full_obj$mouse_final_cell_type),
  selected_myeloid = selected, selection_reason = ifelse(selected, "configured_mouse_myeloid_lineage", "not_configured_myeloid_lineage"),
  stringsAsFactors = FALSE
)
utils::write.csv(selection, file.path(out_dir, "mouse_myeloid_selection_decisions.csv"), row.names = FALSE)

sample_lineage_counts <- as.data.frame(table(
  sample = as.character(myeloid$sample), source_lineage = as.character(myeloid$mouse_myeloid_source_label)
), stringsAsFactors = FALSE)
colnames(sample_lineage_counts)[3] <- "n_cells"
sample_lineage_counts <- sample_lineage_counts[sample_lineage_counts$n_cells > 0L, , drop = FALSE]
utils::write.csv(sample_lineage_counts, file.path(out_dir, "source_myeloid_counts_by_sample.csv"), row.names = FALSE)

source_colors <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
source_short <- sub(paste0("^", as.character(cfg$immune_composition$final_label_prefix), "\\s*"), "", as.character(myeloid$mouse_myeloid_source_label))
myeloid$mouse_myeloid_source_lineage <- factor(source_short, levels = intersect(c("Monocyte", "Macrophage", "Neutrophil"), unique(source_short)))
parent_umap <- Seurat::Embeddings(full_obj, "umap")[selected_cells, , drop = FALSE]
myeloid[["parent_umap_05d"]] <- Seurat::CreateDimReducObject(embeddings = parent_umap, key = "PARENTUMAP_", assay = "RNA")
p <- Seurat::DimPlot(
  myeloid, reduction = "parent_umap_05d", group.by = "mouse_myeloid_source_lineage",
  cols = unname(source_colors[levels(myeloid$mouse_myeloid_source_lineage)]), raster = FALSE
) + ggplot2::labs(title = "Extracted mouse myeloid compartment", subtitle = "Coordinates from the full mouse-cell 05d UMAP")
ggplot2::ggsave(file.path(out_dir, "parent_umap_extracted_mouse_myeloid.pdf"), p, width = 10, height = 8)
ggplot2::ggsave(file.path(out_dir, "parent_umap_extracted_mouse_myeloid.png"), p, width = 10, height = 8, dpi = 300)

myeloid@misc$mouse05_myeloid_extraction_contract <- list(
  source_object = normalizePath(input_rds, mustWork = TRUE), source_labels = source_labels,
  doublet_rule = "05b scDblFinder.class equals singlet", mouse_fraction_min = as.numeric(cfg$species$mouse_fraction_min),
  retained_assays = "RNA only", inherited_reductions = "parent_umap_05d for visualization only",
  cells = ncol(myeloid), samples = sort(unique(as.character(myeloid$sample)))
)
saveRDS(myeloid, file.path(out_dir, "mouse_myeloid_rna_only_preintegration.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05i.txt"))
writeLines(c(
  "05i mouse myeloid extraction completed.", paste0("Cells: ", ncol(myeloid)),
  paste0("Samples: ", length(unique(myeloid$sample))), paste0("Source labels: ", paste(source_labels, collapse = "; ")),
  "Only 05b QC-passing scDblFinder singlets were retained."
), file.path(out_dir, "completion.txt"))
message("[05i] Completed. Output: ", out_dir)
