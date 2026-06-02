#!/usr/bin/env Rscript

resolve_in_vivo_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  candidate_files <- character(0)
  if (length(file_match) > 0) {
    candidate_files <- c(candidate_files, sub("--file=", "", file_match[1]))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  candidate_files <- c(candidate_files, frame_files[!is.na(frame_files)])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_path)) candidate_files <- c(candidate_files, active_path)
  }

  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cwd_parts <- strsplit(cwd, .Platform$file.sep, fixed = TRUE)[[1]]
  parent_dirs <- vapply(
    seq_along(cwd_parts),
    function(i) {
      paste(c(cwd_parts[seq_len(length(cwd_parts) - i + 1)]), collapse = .Platform$file.sep)
    },
    character(1)
  )
  parent_dirs <- parent_dirs[nzchar(parent_dirs)]
  parent_dirs <- if (grepl("^/", cwd)) paste0("/", sub("^/+", "", parent_dirs)) else parent_dirs
  candidate_dirs <- unique(c(
    candidate_dirs,
    cwd,
    file.path(cwd, "Code", "in-vivo"),
    parent_dirs,
    file.path(parent_dirs, "Code", "in-vivo")
  ))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

required_packages <- c(
  "Seurat",
  "dplyr",
  "ggplot2",
  "readr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  KMP_DUPLICATE_LIB_OK = "TRUE",
  KMP_INIT_AT_FORK = "FALSE"
)

set.seed(1234)

input_rds <- file.path(
  results_root,
  "02d_manual_cluster_merge",
  "objects",
  "integrated_sct_cca_seurat_cluster_refine_manual_merge_test.rds"
)
output_root <- file.path(results_root, "03_final_cluster")

filter_cluster_col <- "seurat_cluster_refine"
retained_cluster_col <- "seurat_cluster_refine_retained"
original_cluster_col <- "seurat_clusters_orig"
new_cluster_col <- "seurat_clusters"
sample_col <- "sample"
remove_cells_from_clusters <- c("9", "4", "3", "9c")
remove_metadata_cols <- c("cluster_refine_focus_status", "manual_merge_test")

integrated_assay <- "integrated"
rna_assay <- "RNA"
nfeatures_integration <- 3000L
npcs <- 30L
cluster_resolution <- 0.6
future_globals_maxsize_gb <- 60

deg_min_pct <- 0.10
deg_logfc_threshold <- 0
deg_padj_cutoff <- 0.05
top_markers_n <- 30L
deg_parallel_workers <- as.integer(Sys.getenv("DEG_PARALLEL_WORKERS", unset = "8"))

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_deg <- .ensure_dir(file.path(output_root, "01_reclustered_cluster_vs_rest_DEG"))
out_plots <- .ensure_dir(file.path(output_root, "02_plots"))
out_objects <- .ensure_dir(file.path(output_root, "03_objects"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}

configure_future_for_seurat(max_size_gb = future_globals_maxsize_gb)

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

assert_required_column(obj@meta.data, filter_cluster_col)
assert_required_column(obj@meta.data, new_cluster_col)
assert_required_column(obj@meta.data, sample_col)

if (!(rna_assay %in% names(obj@assays))) {
  stop("RNA assay is missing: ", rna_assay, call. = FALSE)
}

filter_clusters <- as.character(obj@meta.data[[filter_cluster_col]])
if (any(is.na(filter_clusters) | filter_clusters == "")) {
  stop("Missing values detected in ", filter_cluster_col, ".", call. = FALSE)
}

filter_cluster_levels <- if (is.factor(obj@meta.data[[filter_cluster_col]])) {
  as.character(levels(obj@meta.data[[filter_cluster_col]]))
} else {
  sort_maybe_numeric(filter_clusters)
}

missing_remove_clusters <- setdiff(remove_cells_from_clusters, unique(filter_clusters))
if (length(missing_remove_clusters) > 0) {
  stop(
    "Configured clusters to remove were not found in ",
    filter_cluster_col,
    ": ",
    paste(missing_remove_clusters, collapse = ", "),
    call. = FALSE
  )
}

keep_flag <- !(filter_clusters %in% remove_cells_from_clusters)
keep_cells <- rownames(obj@meta.data)[keep_flag]
removed_cells <- rownames(obj@meta.data)[!keep_flag]
retained_filter_levels <- setdiff(filter_cluster_levels, remove_cells_from_clusters)

filter_summary_df <- data.frame(
  metric = c("input_cells", "removed_cells", "retained_cells", "removed_cluster_n", "retained_cluster_n"),
  value = c(
    ncol(obj),
    length(removed_cells),
    length(keep_cells),
    length(remove_cells_from_clusters),
    length(retained_filter_levels)
  ),
  stringsAsFactors = FALSE
)
write_table_csv(filter_summary_df, file.path(out_summary, "filter_summary.csv"))

cluster_counts_before <- data.frame(
  cluster = names(table(filter_clusters)),
  n_cells = as.integer(table(filter_clusters)),
  status = ifelse(names(table(filter_clusters)) %in% remove_cells_from_clusters, "removed", "retained"),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, filter_cluster_levels))
write_table_csv(cluster_counts_before, file.path(out_summary, "seurat_cluster_refine_counts_before_filter.csv"))

removed_cell_df <- obj@meta.data[!keep_flag, , drop = FALSE] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    !!filter_cluster_col := as.character(.data[[filter_cluster_col]])
  ) %>%
  dplyr::select(dplyr::any_of(c("cell", filter_cluster_col, "seurat_clusters", "sample", "Dose")))
write_table_csv(removed_cell_df, file.path(out_summary, "removed_cells.csv"))

message("Subsetting retained cells.")
obj_filtered <- subset(obj, cells = keep_cells)
obj_filtered@meta.data[[retained_cluster_col]] <- factor(
  as.character(obj_filtered@meta.data[[filter_cluster_col]]),
  levels = retained_filter_levels
)

cluster_counts_after_filter <- data.frame(
  cluster = names(table(obj_filtered@meta.data[[retained_cluster_col]])),
  n_cells = as.integer(table(obj_filtered@meta.data[[retained_cluster_col]])),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, retained_filter_levels))
write_table_csv(cluster_counts_after_filter, file.path(out_summary, "seurat_cluster_refine_counts_after_filter.csv"))

obj_filtered@meta.data[[original_cluster_col]] <- obj_filtered@meta.data[[new_cluster_col]]
obj_filtered@meta.data[[new_cluster_col]] <- NULL

if ("integrated_snn_res.0.6" %in% colnames(obj_filtered@meta.data)) {
  obj_filtered@meta.data[["integrated_snn_res.0.6_orig"]] <- obj_filtered@meta.data[["integrated_snn_res.0.6"]]
}

for (meta_col in remove_metadata_cols) {
  if (meta_col %in% colnames(obj_filtered@meta.data)) {
    obj_filtered@meta.data[[meta_col]] <- NULL
  }
}

message("Re-running SCT-CCA integration on retained cells.")
DefaultAssay(obj_filtered) <- rna_assay
obj_filtered <- maybe_join_layers(obj_filtered, assay = rna_assay)
if (!("percent.mt" %in% colnames(obj_filtered@meta.data))) {
  obj_filtered[["percent.mt"]] <- PercentageFeatureSet(obj_filtered, assay = rna_assay, pattern = "^MT-")
}

sample_values <- as.character(obj_filtered@meta.data[[sample_col]])
if (any(is.na(sample_values) | sample_values == "")) {
  stop("Missing sample values detected in metadata column: ", sample_col, call. = FALSE)
}
sample_levels <- sort(unique(sample_values))
if (length(sample_levels) < 2L) {
  stop("At least 2 non-empty samples are required for SCT-CCA integration.", call. = FALSE)
}

retained_sample_count_df <- data.frame(
  sample = names(table(sample_values)),
  n_cells = as.integer(table(sample_values)),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(sample, sample_levels))
write_table_csv(retained_sample_count_df, file.path(out_summary, "retained_sample_cell_counts.csv"))

message("Building retained-cell Seurat list by sample.")
seurat_list <- setNames(vector("list", length(sample_levels)), sample_levels)
for (sample_id in sample_levels) {
  sample_cells <- rownames(obj_filtered@meta.data)[sample_values == sample_id]
  sample_obj <- subset(obj_filtered, cells = sample_cells)
  DefaultAssay(sample_obj) <- rna_assay
  sample_obj <- maybe_join_layers(sample_obj, assay = rna_assay)
  counts <- get_assay_data_slot(sample_obj, assay = rna_assay, slot_name = "counts")
  if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
    stop("RNA counts slot is empty for sample: ", sample_id, call. = FALSE)
  }

  sample_meta <- sample_obj@meta.data[colnames(counts), , drop = FALSE]
  auto_meta_cols <- intersect(
    c("orig.ident", paste0("nCount_", rna_assay), paste0("nFeature_", rna_assay)),
    colnames(sample_meta)
  )
  sample_meta <- sample_meta[, setdiff(colnames(sample_meta), auto_meta_cols), drop = FALSE]

  fresh_obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  fresh_obj <- AddMetaData(fresh_obj, metadata = sample_meta[colnames(fresh_obj), , drop = FALSE])
  seurat_list[[sample_id]] <- fresh_obj
}

message("Running SCTransform.")
seurat_list <- lapply(
  seurat_list,
  function(x) {
    SCTransform(
      object = x,
      assay = rna_assay,
      new.assay.name = "SCT",
      vars.to.regress = "percent.mt",
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
  }
)

message("Selecting features and preparing SCT integration.")
integration_features <- SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = nfeatures_integration
)
configure_future_for_seurat(max_size_gb = future_globals_maxsize_gb)
invisible(gc())
seurat_list <- PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = integration_features,
  verbose = FALSE
)

message("Finding SCT-CCA anchors and integrating retained cells.")
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "cca",
  verbose = FALSE
)
obj_filtered <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)
DefaultAssay(obj_filtered) <- integrated_assay

message("Running PCA / neighbors / clusters / UMAP.")
obj_filtered <- RunPCA(obj_filtered, npcs = max(50L, npcs), verbose = FALSE)
dims_use <- seq_len(min(npcs, ncol(Embeddings(obj_filtered, reduction = "pca"))))
obj_filtered <- FindNeighbors(obj_filtered, dims = dims_use, verbose = FALSE)
obj_filtered <- FindClusters(obj_filtered, resolution = cluster_resolution, verbose = FALSE)
obj_filtered <- RunUMAP(obj_filtered, dims = dims_use, verbose = FALSE)
message("Reintegration and reclustering finished.")

new_clusters <- as.character(obj_filtered@meta.data[[new_cluster_col]])
new_cluster_levels <- sort_maybe_numeric(new_clusters)
obj_filtered@meta.data[[new_cluster_col]] <- factor(new_clusters, levels = new_cluster_levels)
Idents(obj_filtered) <- obj_filtered@meta.data[[new_cluster_col]]

recluster_parameters_df <- data.frame(
  parameter = c(
    "input_rds",
    "output_root",
    "filter_cluster_col",
    "retained_cluster_col",
    "original_cluster_col",
    "new_cluster_col",
    "remove_cells_from_clusters",
    "remove_metadata_cols",
    "sample_col",
    "integrated_assay",
    "rna_assay",
    "nfeatures_integration",
    "npcs",
    "dims_used",
    "cluster_resolution",
    "deg_min_pct",
    "deg_logfc_threshold",
    "deg_parallel_workers"
  ),
  value = c(
    input_rds,
    output_root,
    filter_cluster_col,
    retained_cluster_col,
    original_cluster_col,
    new_cluster_col,
    paste(remove_cells_from_clusters, collapse = ","),
    paste(remove_metadata_cols, collapse = ","),
    sample_col,
    integrated_assay,
    rna_assay,
    as.character(nfeatures_integration),
    as.character(npcs),
    paste(dims_use, collapse = ","),
    as.character(cluster_resolution),
    as.character(deg_min_pct),
    as.character(deg_logfc_threshold),
    as.character(deg_parallel_workers)
  ),
  stringsAsFactors = FALSE
)
write_table_csv(recluster_parameters_df, file.path(out_summary, "recluster_parameters.csv"))

new_cluster_count_df <- data.frame(
  cluster = names(table(obj_filtered@meta.data[[new_cluster_col]])),
  n_cells = as.integer(table(obj_filtered@meta.data[[new_cluster_col]])),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, new_cluster_levels))
write_table_csv(new_cluster_count_df, file.path(out_summary, "new_seurat_cluster_cell_counts.csv"))

if ("umap" %in% names(obj_filtered@reductions)) {
  p_refine <- DimPlot(
    obj_filtered,
    reduction = "umap",
    group.by = retained_cluster_col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.30,
    raster = FALSE
  ) +
    labs(title = "seurat_cluster_refine on final UMAP")
  save_plot_pdf_png(
    plot_obj = p_refine,
    file_stub = file.path(out_plots, "umap_retained_seurat_cluster_refine"),
    width = 9,
    height = 7
  )

  p_new <- DimPlot(
    obj_filtered,
    reduction = "umap",
    group.by = new_cluster_col,
    label = TRUE,
    repel = TRUE,
    pt.size = 0.30,
    raster = FALSE
  ) +
    labs(title = "Final seurat_clusters")
  save_plot_pdf_png(
    plot_obj = p_new,
    file_stub = file.path(out_plots, "umap_new_seurat_clusters"),
    width = 9,
    height = 7
  )
}

message("Preparing RNA assay for DEG.")
obj_filtered <- prepare_obj_for_markers(obj_filtered, assay = rna_assay)
Idents(obj_filtered) <- obj_filtered@meta.data[[new_cluster_col]]

message("Checking existing one-vs-rest DEG results for new seurat_clusters.")
marker_results <- stats::setNames(vector("list", length(new_cluster_levels)), new_cluster_levels)
clusters_to_run_deg <- character(0)
for (cluster_id in new_cluster_levels) {
  deg_dir_cluster <- .ensure_dir(file.path(out_deg, paste0("cluster_", cluster_id)))
  deg_markers_file <- file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_markers.csv"))
  if (file.exists(deg_markers_file)) {
    message("Reusing existing DEG file: ", deg_markers_file)
    de_cached <- readr::read_csv(deg_markers_file, show_col_types = FALSE)
    marker_results[[cluster_id]] <- list(
      cluster = as.character(cluster_id),
      markers = as.data.frame(de_cached, stringsAsFactors = FALSE),
      error = NULL,
      cached = TRUE
    )
  } else {
    clusters_to_run_deg <- c(clusters_to_run_deg, cluster_id)
  }
}

if (length(clusters_to_run_deg) > 0) {
  message("Running missing one-vs-rest DEG for new seurat_clusters with parallel Seurat::FindMarkers.")
  marker_results_new <- run_seurat_cluster_markers_parallel(
    obj = obj_filtered,
    cluster_levels = clusters_to_run_deg,
    assay = rna_assay,
    min_pct = deg_min_pct,
    logfc_threshold = deg_logfc_threshold,
    workers = deg_parallel_workers
  )
  deg_workers_used <- attr(marker_results_new, "workers")
  for (cluster_id in clusters_to_run_deg) {
    marker_results[[cluster_id]] <- marker_results_new[[cluster_id]]
    marker_results[[cluster_id]]$cached <- FALSE
  }
} else {
  message("All new seurat_clusters DEG files already exist; skipping DEG analysis.")
  deg_workers_used <- 0L
}

deg_summary_list <- list()
top_positive_list <- list()
all_markers_list <- list()

for (cluster_id in new_cluster_levels) {
  deg_dir_cluster <- .ensure_dir(file.path(out_deg, paste0("cluster_", cluster_id)))
  cluster_result <- marker_results[[cluster_id]]
  de_cluster <- cluster_result$markers

  if (!is.null(cluster_result$error)) {
    msg <- paste0("FindMarkers failed for seurat_clusters ", cluster_id, ": ", cluster_result$error)
    warning(msg)
    writeLines(msg, con = file.path(deg_dir_cluster, "FindMarkers_error.txt"))
  }

  if (is.null(de_cluster) || nrow(de_cluster) == 0) {
    deg_summary_list[[cluster_id]] <- data.frame(
      cluster = as.character(cluster_id),
      n_cells = new_cluster_count_df$n_cells[match(cluster_id, new_cluster_count_df$cluster)],
      tested_genes = 0L,
      sig_genes = 0L,
      top_positive_gene = NA_character_,
      stringsAsFactors = FALSE
    )
    top_positive_list[[cluster_id]] <- data.frame()
    all_markers_list[[cluster_id]] <- data.frame()
    next
  }

  deg_markers_file <- file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_vs_rest_markers.csv"))
  if (!isTRUE(cluster_result$cached) || !file.exists(deg_markers_file)) {
    write_table_csv(de_cluster, deg_markers_file)
  }
  lfc_col <- resolve_lfc_col(de_cluster)

  top_positive <- de_cluster %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < deg_padj_cutoff, .data[[lfc_col]] > 0) %>%
    dplyr::arrange(dplyr::desc(.data[[lfc_col]]), p_val_adj, gene) %>%
    dplyr::slice_head(n = top_markers_n) %>%
    dplyr::mutate(rank_within_cluster = dplyr::row_number())
  write_table_csv(top_positive, file.path(deg_dir_cluster, paste0("cluster_", cluster_id, "_top_positive_markers.csv")))

  deg_summary_list[[cluster_id]] <- data.frame(
    cluster = as.character(cluster_id),
    n_cells = new_cluster_count_df$n_cells[match(cluster_id, new_cluster_count_df$cluster)],
    tested_genes = nrow(de_cluster),
    sig_genes = sum(!is.na(de_cluster$p_val_adj) & de_cluster$p_val_adj < deg_padj_cutoff),
    top_positive_gene = ifelse(nrow(top_positive) == 0, NA_character_, top_positive$gene[1]),
    stringsAsFactors = FALSE
  )

  top_positive_list[[cluster_id]] <- top_positive
  all_markers_list[[cluster_id]] <- de_cluster
}

deg_summary_df <- bind_rows(deg_summary_list) %>%
  dplyr::arrange(match(cluster, new_cluster_levels))
write_table_csv(deg_summary_df, file.path(out_summary, "new_seurat_cluster_deg_summary.csv"))

all_markers <- if (length(all_markers_list) > 0) bind_rows(all_markers_list) else data.frame()
write_table_csv(all_markers, file.path(out_summary, "new_seurat_clusters_all_markers.csv"))

top_positive_df <- if (length(top_positive_list) > 0) bind_rows(top_positive_list) else data.frame()
write_table_csv(top_positive_df, file.path(out_summary, "new_seurat_cluster_top_positive_markers_all.csv"))

message("Saving final Seurat object.")
saveRDS(
  obj_filtered,
  file.path(out_objects, "integrated_sct_cca_seurat_final_reclustered.rds")
)

summary_lines <- c(
  "03 final cluster completed.",
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Filter cluster column: ", filter_cluster_col),
  paste0("Removed cells from clusters: ", paste(remove_cells_from_clusters, collapse = ", ")),
  paste0("Removed metadata columns: ", paste(remove_metadata_cols, collapse = ", ")),
  paste0("Cells removed: ", length(removed_cells)),
  paste0("Cells retained: ", length(keep_cells)),
  paste0("Retained ", filter_cluster_col, " labels: ", paste(retained_filter_levels, collapse = ", ")),
  paste0("Integration workflow: SCTransform -> SCT-CCA anchors -> IntegrateData"),
  paste0("Integration features: ", nfeatures_integration),
  paste0("New seurat_clusters count: ", length(new_cluster_levels)),
  paste0("New seurat_clusters labels: ", paste(new_cluster_levels, collapse = ", ")),
  paste0("PCA/neighbor dims used: ", paste(dims_use, collapse = ",")),
  paste0("Cluster resolution: ", cluster_resolution),
  paste0("DEG backend: Seurat::FindMarkers via parallel::mclapply fork"),
  paste0("DEG method: Seurat::FindMarkers default"),
  paste0("DEG parallel workers requested: ", deg_parallel_workers),
  paste0("DEG parallel workers used: ", deg_workers_used)
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("03_final_cluster finished.")
message("Output root: ", output_root)
