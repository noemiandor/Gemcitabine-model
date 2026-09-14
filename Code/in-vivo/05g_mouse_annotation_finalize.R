#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
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
required_packages <- c("yaml", "Seurat", "Matrix", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05f_immune_annotation", "mouse_broad_and_immune_annotated.rds")
feature_map_file <- file.path(results_root, "05a_species_extraction", "GRCm39_feature_mapping.csv")
out_dir <- file.path(results_root, "05g_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05f input: ", input_rds, call. = FALSE)

obj <- readRDS(input_rds)
feature_map <- utils::read.csv(feature_map_file, stringsAsFactors = FALSE, check.names = FALSE)
if (sum(startsWith(rownames(obj[["RNA"]]), cfg$species$human_prefix)) != 0L) stop("Human-prefixed features detected in final object.", call. = FALSE)
if (any(feature_map$source_genome != "GRCm39") || any(startsWith(feature_map$source_feature_id, cfg$species$human_prefix))) {
  stop("Final feature provenance is not exclusively GRCm39.", call. = FALSE)
}
if (any(obj$mouse_fraction < as.numeric(cfg$species$mouse_fraction_min), na.rm = TRUE)) {
  stop("A final cell violates the >=90% mouse-transcript rule.", call. = FALSE)
}

required_annotation_columns <- c(
  "mouse_broad_cell_type", "mouse_broad_confidence", "mouse_broad_signature",
  "mouse_broad_overlap", "mouse_broad_f1", "mouse_broad_overlap_genes",
  "mouse_immune_candidate", "mouse_immune_cell_type", "mouse_immune_confidence",
  "mouse_immune_signature", "mouse_immune_overlap", "mouse_immune_f1", "mouse_immune_overlap_genes"
)
missing_annotation_columns <- setdiff(required_annotation_columns, colnames(obj@meta.data))
if (length(missing_annotation_columns) > 0L) {
  stop("05f object lacks cluster-DEG annotation columns: ", paste(missing_annotation_columns, collapse = ", "), call. = FALSE)
}

immune_usable <- obj$mouse_immune_candidate &
  !is.na(obj$mouse_immune_cell_type) &
  obj$mouse_immune_cell_type != "Ambiguous"

final_type <- as.character(obj$mouse_broad_cell_type)
final_type[immune_usable] <- paste0("Immune: ", obj$mouse_immune_cell_type[immune_usable])
final_type[obj$mouse_immune_candidate & !immune_usable & obj$mouse_broad_cell_type == "Immune"] <- "Immune: Ambiguous"
final_type[is.na(final_type) | !nzchar(final_type)] <- "Ambiguous"

final_confidence <- as.character(obj$mouse_broad_confidence)
final_confidence[immune_usable] <- as.character(obj$mouse_immune_confidence[immune_usable])
final_confidence[obj$mouse_immune_candidate & !immune_usable & obj$mouse_broad_cell_type == "Immune"] <- "ambiguous"

evidence <- paste0(
  "broad_signature=", ifelse(is.na(obj$mouse_broad_signature), "NA", obj$mouse_broad_signature),
  ";broad_overlap=", ifelse(is.na(obj$mouse_broad_overlap), "NA", obj$mouse_broad_overlap),
  ";broad_f1=", ifelse(is.na(obj$mouse_broad_f1), "NA", signif(obj$mouse_broad_f1, 5)),
  ";broad_overlap_genes=", ifelse(is.na(obj$mouse_broad_overlap_genes), "NA", obj$mouse_broad_overlap_genes),
  ";immune_signature=", ifelse(is.na(obj$mouse_immune_signature), "NA", obj$mouse_immune_signature),
  ";immune_overlap=", ifelse(is.na(obj$mouse_immune_overlap), "NA", obj$mouse_immune_overlap),
  ";immune_f1=", ifelse(is.na(obj$mouse_immune_f1), "NA", signif(obj$mouse_immune_f1, 5)),
  ";immune_overlap_genes=", ifelse(is.na(obj$mouse_immune_overlap_genes), "NA", obj$mouse_immune_overlap_genes)
)

obj$mouse_final_cell_type <- final_type
obj$mouse_final_confidence <- final_confidence
obj$mouse_annotation_evidence <- evidence

# The final label is also strictly cluster-level. Mixed labels within a cluster
# indicate an upstream contract violation and are not allowed to pass silently.
final_labels_per_cluster <- tapply(final_type, as.character(obj$seurat_clusters), function(x) length(unique(x)))
if (any(final_labels_per_cluster != 1L)) stop("Final annotations are not uniform within every cluster.", call. = FALSE)

# The config is the single source of truth for final UMAP colors. Immune labels
# use the exact same subtype colors as all downstream 05z figures.
immune_colors <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
nonimmune_colors <- unlist(cfg$annotation$final_nonimmune_colors, use.names = TRUE)
if (is.null(names(immune_colors)) || any(!nzchar(names(immune_colors)))) {
  stop("immune_composition.subtype_colors must be a named mapping.", call. = FALSE)
}
if (is.null(names(nonimmune_colors)) || any(!nzchar(names(nonimmune_colors)))) {
  stop("annotation.final_nonimmune_colors must be a named mapping.", call. = FALSE)
}
immune_prefix <- as.character(cfg$immune_composition$final_label_prefix)
immune_final_colors <- stats::setNames(
  unname(immune_colors),
  paste0(immune_prefix, " ", names(immune_colors))
)
configured_final_colors <- c(nonimmune_colors, immune_final_colors)
final_type_levels <- sort(unique(final_type))
missing_final_colors <- setdiff(final_type_levels, names(configured_final_colors))
if (length(missing_final_colors) > 0L) {
  stop("Final UMAP palette is missing labels: ", paste(missing_final_colors, collapse = ", "), call. = FALSE)
}
final_palette <- configured_final_colors[final_type_levels]
obj$mouse_final_cell_type <- factor(final_type, levels = final_type_levels)
final_palette_df <- data.frame(
  final_cell_type = final_type_levels,
  color_hex = unname(final_palette),
  stringsAsFactors = FALSE
)
utils::write.csv(final_palette_df, file.path(out_dir, "mouse_final_cell_type_palette.csv"), row.names = FALSE)

cell_df <- data.frame(
  cell = colnames(obj),
  sample = obj$sample,
  barcode_raw = obj$barcode_raw,
  cluster = as.character(obj$seurat_clusters),
  human_umi = obj$human_umi,
  mouse_umi = obj$mouse_umi,
  mouse_fraction = obj$mouse_fraction,
  mouse_percent = obj$mouse_percent,
  nCount_RNA = obj$nCount_RNA,
  nFeature_RNA = obj$nFeature_RNA,
  percent.mt = obj$percent.mt,
  broad_cell_type = obj$mouse_broad_cell_type,
  broad_confidence = obj$mouse_broad_confidence,
  broad_signature = obj$mouse_broad_signature,
  broad_overlap = obj$mouse_broad_overlap,
  broad_f1 = obj$mouse_broad_f1,
  broad_overlap_genes = obj$mouse_broad_overlap_genes,
  immune_candidate = obj$mouse_immune_candidate,
  immune_cell_type = obj$mouse_immune_cell_type,
  immune_subtype = obj$mouse_immune_subtype,
  immune_confidence = obj$mouse_immune_confidence,
  immune_signature = obj$mouse_immune_signature,
  immune_overlap = obj$mouse_immune_overlap,
  immune_f1 = obj$mouse_immune_f1,
  immune_overlap_genes = obj$mouse_immune_overlap_genes,
  final_cell_type = obj$mouse_final_cell_type,
  final_confidence = obj$mouse_final_confidence,
  annotation_evidence = obj$mouse_annotation_evidence,
  stringsAsFactors = FALSE
)
utils::write.csv(cell_df, file.path(out_dir, "mouse_cell_annotations_final.csv"), row.names = FALSE, na = "NA")

cluster_summary <- as.data.frame(table(cluster = cell_df$cluster, final_cell_type = cell_df$final_cell_type), stringsAsFactors = FALSE)
cluster_summary <- cluster_summary[cluster_summary$Freq > 0, , drop = FALSE]
colnames(cluster_summary)[3] <- "n_cells"
cluster_totals <- tapply(cluster_summary$n_cells, cluster_summary$cluster, sum)
cluster_summary$fraction_within_cluster <- cluster_summary$n_cells / unname(cluster_totals[cluster_summary$cluster])
utils::write.csv(cluster_summary, file.path(out_dir, "mouse_cluster_annotation_summary.csv"), row.names = FALSE)

sample_summary <- as.data.frame(table(sample = cell_df$sample, final_cell_type = cell_df$final_cell_type), stringsAsFactors = FALSE)
sample_summary <- sample_summary[sample_summary$Freq > 0, , drop = FALSE]
colnames(sample_summary)[3] <- "n_cells"
utils::write.csv(sample_summary, file.path(out_dir, "mouse_sample_annotation_summary.csv"), row.names = FALSE)

Seurat::DefaultAssay(obj) <- "RNA"
expr <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
canonical_markers <- intersect(
  c("Ptprc", "Cd3d", "Cd3e", "Trac", "Nkg7", "Klrd1", "Cd79a", "Ms4a1", "Jchain", "Lyz2", "Adgre1", "Csf1r", "Ly6c2", "S100a8", "S100a9", "Ly6g", "Itgax", "Flt3", "Kit", "Cpa3", "Pecam1", "Vwf", "Col1a1", "Dcn", "Rgs5", "Pdgfrb", "Epcam", "Krt8", "Adipoq", "Plin1"),
  rownames(expr)
)
marker_audit <- do.call(rbind, lapply(sort(unique(final_type)), function(ct) {
  cells <- which(final_type == ct)
  if (length(cells) == 0L || length(canonical_markers) == 0L) return(NULL)
  data.frame(
    final_cell_type = ct,
    gene = canonical_markers,
    mean_log_normalized_expression = as.numeric(Matrix::rowMeans(expr[canonical_markers, cells, drop = FALSE])),
    fraction_detected = as.numeric(Matrix::rowMeans(expr[canonical_markers, cells, drop = FALSE] > 0)),
    n_cells = length(cells),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(marker_audit, file.path(out_dir, "final_annotation_marker_audit.csv"), row.names = FALSE, na = "NA")

p <- Seurat::DimPlot(
  obj,
  reduction = "umap",
  group.by = "mouse_final_cell_type",
  cols = unname(final_palette),
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggplot2::labs(title = "Final mouse cell-type annotation")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_final_cell_type.pdf"), p, width = 12, height = 9)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_final_cell_type.png"), p, width = 12, height = 9, dpi = 300)

obj@misc$mouse05_final_contract <- list(
  mouse_fraction_min = as.numeric(cfg$species$mouse_fraction_min),
  human_features_retained = 0L,
  mouse_feature_provenance = "GRCm39 only",
  annotation_level = "cluster",
  annotation_method = "top positive cluster DEG to curated signature overlap, ranked by F1 then overlap then Jaccard",
  cell_level_classifier_used = FALSE,
  final_cells = ncol(obj),
  final_cell_types = final_type_levels,
  final_cell_type_palette = as.list(final_palette)
)
final_rds <- file.path(out_dir, "mouse_integrated_clustered_annotated_final.rds")
final_rds_tmp <- file.path(out_dir, paste0(".mouse_integrated_clustered_annotated_final.rds.tmp.", Sys.getpid()))
on.exit(unlink(final_rds_tmp, force = TRUE), add = TRUE)
saveRDS(obj, final_rds_tmp, compress = FALSE)
if (!file.rename(final_rds_tmp, final_rds)) {
  stop("Failed to atomically replace the final annotated RDS.", call. = FALSE)
}

output_files <- c(
  final_rds,
  file.path(out_dir, "mouse_cell_annotations_final.csv"),
  file.path(out_dir, "mouse_cluster_annotation_summary.csv"),
  file.path(out_dir, "mouse_sample_annotation_summary.csv"),
  file.path(out_dir, "final_annotation_marker_audit.csv"),
  file.path(out_dir, "mouse_final_cell_type_palette.csv"),
  file.path(out_dir, "umap_mouse_final_cell_type.pdf"),
  file.path(out_dir, "umap_mouse_final_cell_type.png")
)
fi <- file.info(output_files)
manifest <- data.frame(
  path = normalizePath(output_files, mustWork = TRUE),
  size_bytes = fi$size,
  sha256 = vapply(output_files, function(x) digest::digest(x, algo = "sha256", file = TRUE, serialize = FALSE), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "final_output_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05g.txt"))
writeLines(
  c(
    "05g final cluster-level annotation completed.",
    paste0("Final cells: ", ncol(obj)),
    paste0("Final cell types: ", length(unique(final_type))),
    "Final UMAP colors are fixed in 05_mouse_cell_analysis_config.yaml.",
    "Immune colors are shared with all downstream 05z figures.",
    "No cell-level classifier was used.",
    paste0("Minimum observed mouse fraction: ", min(obj$mouse_fraction, na.rm = TRUE)),
    "Human features retained: 0"
  ),
  file.path(out_dir, "completion.txt")
)
message("[05g] Completed. Output: ", out_dir)
