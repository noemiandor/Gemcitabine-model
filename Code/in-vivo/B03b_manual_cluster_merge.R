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

required_utils_functions <- c(
  "build_manual_merge_labels",
  "make_manual_merge_mapping_table",
  "make_internal_pair_df_from_merge_map",
  "run_markers_vs_rest_default",
  "run_seurat_cluster_markers_parallel",
  "run_pairwise_markers_default",
  "build_manual_merge_cosine_summary",
  "cosine_summary_to_matrix",
  "plot_cosine_heatmap"
)
missing_utils_functions <- required_utils_functions[
  !vapply(required_utils_functions, exists, logical(1), mode = "function")
]
if (length(missing_utils_functions) > 0) {
  stop(
    "Required function(s) were not loaded from Utils.R: ",
    paste(missing_utils_functions, collapse = ", "),
    ". Run the script from the top, or first run: source('",
    file.path(script_dir, "Utils.R"),
    "').",
    call. = FALSE
  )
}

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
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
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
  "03_final_cluster",
  "03_objects",
  "integrated_sct_cca_seurat_final_reclustered.rds"
)
output_root <- file.path(results_root, "03b_manual_cluster_merge")

base_cluster_col <- "seurat_clusters"
analysis_cluster_col <- "cluster_final"
pre_renumber_cluster_col <- "cluster_final_pre_renumber"
expected_merged_cluster_count <- 9L

manual_merge_map <- list(
  "8" = c("8", "9", "11"),
  "1" = c("1", "10"),
  "3" = c("3", "2", "7"),
  "0" = c("0", "12")
)

marker_min_pct <- 0.10
marker_logfc_threshold <- 0.25
deg_min_pct <- 0.10
deg_logfc_threshold <- 0
deg_padj_cutoff <- 0.05
deg_parallel_workers <- as.integer(Sys.getenv("DEG_PARALLEL_WORKERS", unset = "8"))
top_markers_n <- 20L
top_dotplot_n <- 5L
cosine_top_n_each <- 50L

message("Preparing output directories.")
.ensure_dir(output_root)
out_plots <- .ensure_dir(file.path(output_root, "plots"))
out_markers_original <- .ensure_dir(file.path(output_root, "markers_original"))
out_markers_merged <- .ensure_dir(file.path(output_root, "markers_merged"))
out_markers_internal <- .ensure_dir(file.path(output_root, "markers_internal_pairs"))
out_summary <- .ensure_dir(file.path(output_root, "summaries"))
out_objects <- .ensure_dir(file.path(output_root, "objects"))

if (!file.exists(input_rds)) {
  stop("Input RDS does not exist: ", input_rds, call. = FALSE)
}

configure_future_for_seurat(max_size_gb = 30)

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

assert_required_column(obj@meta.data, base_cluster_col)
if (!("umap" %in% names(obj@reductions))) {
  stop("UMAP reduction not found in Seurat object.", call. = FALSE)
}

obj <- prepare_obj_for_markers(obj, assay = "RNA")

cluster_vec <- as.character(obj@meta.data[[base_cluster_col]])
if (any(is.na(cluster_vec) | cluster_vec == "")) {
  stop("Missing values detected in metadata column: ", base_cluster_col, call. = FALSE)
}

cluster_levels <- if (is.factor(obj@meta.data[[base_cluster_col]])) {
  as.character(levels(obj@meta.data[[base_cluster_col]]))
} else {
  sort_maybe_numeric(cluster_vec)
}
obj@meta.data[[base_cluster_col]] <- factor(cluster_vec, levels = cluster_levels)

configured_sources <- unique(as.character(unlist(manual_merge_map, use.names = FALSE)))
missing_sources <- setdiff(configured_sources, unique(cluster_vec))
if (length(missing_sources) > 0) {
  stop(
    "Manual merge source cluster(s) not found in ",
    base_cluster_col,
    ": ",
    paste(missing_sources, collapse = ", "),
    call. = FALSE
  )
}

obj@meta.data[[pre_renumber_cluster_col]] <- build_manual_merge_labels(
  cluster_vec,
  manual_merge_map,
  cluster_order = cluster_levels
)
pre_renumber_levels <- levels(obj@meta.data[[pre_renumber_cluster_col]])
renumber_source_levels <- sort_maybe_numeric(pre_renumber_levels)
renumber_map <- stats::setNames(
  as.character(seq_along(renumber_source_levels) - 1L),
  renumber_source_levels
)
manual_merge_map_renumbered <- manual_merge_map
names(manual_merge_map_renumbered) <- unname(renumber_map[names(manual_merge_map)])
manual_merge_map_renumbered <- manual_merge_map_renumbered[sort_maybe_numeric(names(manual_merge_map_renumbered))]
obj@meta.data[[analysis_cluster_col]] <- factor(
  unname(renumber_map[as.character(obj@meta.data[[pre_renumber_cluster_col]])]),
  levels = as.character(seq_along(renumber_source_levels) - 1L)
)
merged_levels <- levels(obj@meta.data[[analysis_cluster_col]])
if (length(merged_levels) != expected_merged_cluster_count) {
  stop(
    "Expected ",
    expected_merged_cluster_count,
    " merged clusters, but found ",
    length(merged_levels),
    ": ",
    paste(merged_levels, collapse = ", "),
    call. = FALSE
  )
}
Idents(obj) <- obj@meta.data[[analysis_cluster_col]]

merge_map_df <- make_manual_merge_mapping_table(cluster_levels, manual_merge_map)
merge_map_df$cluster_final_pre_renumber <- merge_map_df$merged_group
merge_map_df$cluster_final <- unname(renumber_map[merge_map_df$merged_group])
merge_map_df <- merge_map_df %>%
  dplyr::select(
    original_cluster,
    cluster_final_pre_renumber,
    cluster_final,
    merge_status
  )
write_table_csv(merge_map_df, file.path(out_summary, "manual_merge_mapping.csv"))

renumber_map_df <- data.frame(
  cluster_final_pre_renumber = names(renumber_map),
  cluster_final = unname(renumber_map),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(as.numeric(cluster_final))
write_table_csv(renumber_map_df, file.path(out_summary, "cluster_final_renumber_mapping.csv"))

manual_merge_membership_df <- data.frame(
  cluster_final_pre_renumber = names(manual_merge_map),
  cluster_final = unname(renumber_map[names(manual_merge_map)]),
  source_clusters = vapply(manual_merge_map, function(x) paste(x, collapse = ","), character(1)),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(as.numeric(cluster_final))
write_table_csv(manual_merge_membership_df, file.path(out_summary, "manual_merge_membership_renumbered.csv"))

cluster_count_df <- data.frame(
  cluster = names(table(obj@meta.data[[base_cluster_col]])),
  n_cells = as.integer(table(obj@meta.data[[base_cluster_col]])),
  stringsAsFactors = FALSE
)
cluster_count_df <- cluster_count_df %>%
  dplyr::arrange(match(cluster, cluster_levels))
write_table_csv(cluster_count_df, file.path(out_summary, "original_cluster_cell_counts.csv"))

merged_count_df <- data.frame(
  cluster = names(table(obj@meta.data[[analysis_cluster_col]])),
  n_cells = as.integer(table(obj@meta.data[[analysis_cluster_col]])),
  stringsAsFactors = FALSE
)
merged_count_df <- merged_count_df %>%
  dplyr::arrange(match(cluster, merged_levels))
write_table_csv(merged_count_df, file.path(out_summary, "cluster_final_cell_counts.csv"))

message("Writing UMAP plots.")
p_umap_original <- DimPlot(
  obj,
  reduction = "umap",
  group.by = base_cluster_col,
  label = FALSE,
  repel = FALSE,
  pt.size = 0.30,
  raster = FALSE
) +
  labs(title = paste("Original clusters:", base_cluster_col))
save_plot_pdf_png(p_umap_original, file.path(out_plots, "umap_original_clusters"), width = 9, height = 7)

p_umap_merged <- DimPlot(
  obj,
  reduction = "umap",
  group.by = analysis_cluster_col,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.30,
  raster = FALSE
) +
  labs(title = "Manual merged clusters")
save_plot_pdf_png(p_umap_merged, file.path(out_plots, "umap_cluster_final"), width = 9, height = 7)

focus_clusters <- configured_sources
focus_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[base_cluster_col]]) %in% focus_clusters]
focus_obj <- subset(obj, cells = focus_cells)

p_focus_original <- DimPlot(
  focus_obj,
  reduction = "umap",
  group.by = base_cluster_col,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.35,
  raster = FALSE
) +
  labs(title = "Focus UMAP: original clusters used for manual merge")
save_plot_pdf_png(p_focus_original, file.path(out_plots, "umap_focus_original_clusters"), width = 8, height = 6)

p_focus_merged <- DimPlot(
  focus_obj,
  reduction = "umap",
  group.by = analysis_cluster_col,
  label = TRUE,
  repel = TRUE,
  pt.size = 0.35,
  raster = FALSE
) +
  labs(title = "Focus UMAP: manual merged clusters")
save_plot_pdf_png(p_focus_merged, file.path(out_plots, "umap_focus_cluster_final"), width = 8, height = 6)

message("Running original-cluster markers.")
original_marker_groups <- sort_maybe_numeric(configured_sources)
original_markers <- run_markers_vs_rest_default(
  obj = obj,
  ident_col = base_cluster_col,
  groups = original_marker_groups,
  out_dir = out_markers_original,
  min_pct = marker_min_pct,
  logfc_threshold = marker_logfc_threshold,
  assay = "RNA",
  top_n = top_markers_n
)
write_table_csv(original_markers$full, file.path(out_markers_original, "markers_original_focus_clusters_all.csv"))
write_table_csv(original_markers$top_positive, file.path(out_markers_original, "markers_original_focus_clusters_top_positive.csv"))

message("Running merged-group markers.")
merged_marker_groups <- merged_levels
message("Checking existing one-vs-rest DEG results for cluster_final.")
marker_results <- stats::setNames(vector("list", length(merged_marker_groups)), merged_marker_groups)
clusters_to_run_deg <- character(0)
for (cluster_id in merged_marker_groups) {
  deg_markers_file <- file.path(out_markers_merged, paste0("markers_", cluster_id, "_vs_rest.csv"))
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
  message("Running missing one-vs-rest DEG for cluster_final with parallel Seurat::FindMarkers.")
  marker_results_new <- run_seurat_cluster_markers_parallel(
    obj = obj,
    cluster_levels = clusters_to_run_deg,
    assay = "RNA",
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
  message("All cluster_final DEG files already exist; skipping DEG analysis.")
  deg_workers_used <- 0L
}

merged_all_markers <- list()
merged_top_markers <- list()
merged_marker_tables <- list()
merged_deg_summary <- list()

for (cluster_id in merged_marker_groups) {
  cluster_result <- marker_results[[cluster_id]]
  de <- cluster_result$markers

  if (!is.null(cluster_result$error)) {
    msg <- paste0("FindMarkers failed for cluster_final ", cluster_id, ": ", cluster_result$error)
    warning(msg)
    writeLines(msg, con = file.path(out_markers_merged, paste0("FindMarkers_error_cluster_", cluster_id, ".txt")))
  }

  if (is.null(de) || nrow(de) == 0) {
    merged_deg_summary[[cluster_id]] <- data.frame(
      cluster = as.character(cluster_id),
      n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
      tested_genes = 0L,
      sig_genes = 0L,
      top_positive_gene = NA_character_,
      stringsAsFactors = FALSE
    )
    merged_top_markers[[cluster_id]] <- data.frame()
    merged_all_markers[[cluster_id]] <- data.frame()
    merged_marker_tables[[cluster_id]] <- list(full = data.frame(), top_positive = data.frame(), lfc_col = NA_character_)
    next
  }

  lfc_col <- resolve_lfc_col(de)
  de$group <- as.character(cluster_id)
  de$comparison <- paste0(cluster_id, "_vs_rest")
  de <- de[, c("group", "comparison", setdiff(colnames(de), c("group", "comparison"))), drop = FALSE]
  deg_markers_file <- file.path(out_markers_merged, paste0("markers_", cluster_id, "_vs_rest.csv"))
  if (!isTRUE(cluster_result$cached) || !file.exists(deg_markers_file)) {
    write_table_csv(de, deg_markers_file)
  }

  top_pos <- de %>%
    dplyr::filter(!is.na(p_val_adj), p_val_adj < deg_padj_cutoff, .data[[lfc_col]] > 0) %>%
    dplyr::arrange(dplyr::desc(.data[[lfc_col]]), p_val_adj, gene) %>%
    dplyr::slice_head(n = top_markers_n) %>%
    dplyr::mutate(rank_within_group = dplyr::row_number())
  write_table_csv(top_pos, file.path(out_markers_merged, paste0("top_positive_markers_", cluster_id, "_vs_rest.csv")))

  merged_deg_summary[[cluster_id]] <- data.frame(
    cluster = as.character(cluster_id),
    n_cells = merged_count_df$n_cells[match(cluster_id, merged_count_df$cluster)],
    tested_genes = nrow(de),
    sig_genes = sum(!is.na(de$p_val_adj) & de$p_val_adj < deg_padj_cutoff),
    top_positive_gene = ifelse(nrow(top_pos) == 0, NA_character_, top_pos$gene[1]),
    stringsAsFactors = FALSE
  )
  merged_all_markers[[cluster_id]] <- de
  merged_top_markers[[cluster_id]] <- top_pos
  merged_marker_tables[[cluster_id]] <- list(full = de, top_positive = top_pos, lfc_col = lfc_col)
}

merged_markers <- list(
  full = dplyr::bind_rows(merged_all_markers),
  top_positive = dplyr::bind_rows(merged_top_markers),
  by_group = merged_marker_tables
)
merged_deg_summary_df <- dplyr::bind_rows(merged_deg_summary) %>%
  dplyr::arrange(match(cluster, merged_marker_groups))
write_table_csv(merged_deg_summary_df, file.path(out_markers_merged, "cluster_final_deg_summary.csv"))
write_table_csv(merged_markers$full, file.path(out_markers_merged, "markers_manual_merged_groups_all.csv"))
write_table_csv(merged_markers$top_positive, file.path(out_markers_merged, "markers_manual_merged_groups_top_positive.csv"))

message("Running internal pairwise marker checks.")
pair_df <- make_internal_pair_df_from_merge_map(manual_merge_map)
write_table_csv(pair_df, file.path(out_summary, "internal_pairwise_marker_pairs.csv"))
internal_pairwise <- run_pairwise_markers_default(
  obj = focus_obj,
  ident_col = base_cluster_col,
  pair_df = pair_df,
  out_dir = out_markers_internal,
  min_pct = marker_min_pct,
  logfc_threshold = marker_logfc_threshold,
  assay = "RNA"
)

message("Computing top50_each signed-LFC cosine similarities.")
cosine_feature_spec <- list(
  type = "top_n",
  padj_max = 0.05,
  abs_log2fc_min = 0.25,
  abs_delta_pct_min = 0.05,
  top_n_up = cosine_top_n_each,
  top_n_down = cosine_top_n_each
)
cosine_summary <- build_manual_merge_cosine_summary(
  original_markers = original_markers,
  merged_markers = merged_markers,
  original_groups = original_marker_groups,
  merged_groups = merged_marker_groups,
  merge_map = manual_merge_map_renumbered,
  spec = cosine_feature_spec
)
write_table_csv(cosine_summary, file.path(out_summary, "manual_merge_vs_original_top50_each_cosine_similarity.csv"))

cosine_mat <- cosine_summary_to_matrix(
  cosine_df = cosine_summary,
  merged_groups = merged_marker_groups,
  original_groups = original_marker_groups
)
write_matrix_csv(
  cosine_mat,
  file.path(out_summary, "manual_merge_vs_original_top50_each_cosine_matrix.csv"),
  row_id = "merged_group"
)
plot_cosine_heatmap(
  mat = cosine_mat,
  title_text = "Merged groups vs original clusters | top50_each signed-LFC cosine",
  file_stub = file.path(out_plots, "heatmap_manual_merge_vs_original_top50_each_cosine")
)

marker_overlap_summary <- build_marker_overlap_summary(
  original_markers = original_markers,
  merged_markers = merged_markers,
  merge_map = manual_merge_map_renumbered,
  n_top = cosine_top_n_each
)
write_table_csv(marker_overlap_summary, file.path(out_summary, "manual_merge_top_positive_marker_overlap.csv"))

message("Writing dotplots.")
original_dotplot_features <- extract_features_for_dotplot(
  top_marker_df = original_markers$top_positive,
  group_col = "group",
  n_per_group = top_dotplot_n
)
plot_dotplot_safe(
  obj = focus_obj,
  features = original_dotplot_features,
  group.by = base_cluster_col,
  title_text = "Top original-cluster markers for manual merge sources",
  file_stub = file.path(out_plots, "dotplot_original_focus_clusters"),
  assay = "RNA"
)

merged_dotplot_features <- extract_features_for_dotplot(
  top_marker_df = merged_markers$top_positive,
  group_col = "group",
  n_per_group = top_dotplot_n
)
focus_merged_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[analysis_cluster_col]]) %in% merged_marker_groups]
focus_merged_obj <- subset(obj, cells = focus_merged_cells)
plot_dotplot_safe(
  obj = focus_merged_obj,
  features = merged_dotplot_features,
  group.by = analysis_cluster_col,
  title_text = "Top merged-cluster markers after manual merge",
  file_stub = file.path(out_plots, "dotplot_manual_merged_groups"),
  assay = "RNA"
)

message("Saving merged Seurat object.")
output_rds <- file.path(out_objects, "integrated_sct_cca_seurat_final_manual_merge.rds")
saveRDS(obj, output_rds)

summary_lines <- c(
  "03b manual cluster merge completed.",
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Base cluster column: ", base_cluster_col),
  paste0("Analysis cluster column: ", analysis_cluster_col),
  paste0("Pre-renumber cluster column: ", pre_renumber_cluster_col),
  paste0("Merged cluster count: ", length(merged_levels)),
  paste0("Merged cluster labels: ", paste(merged_levels, collapse = ", ")),
  paste0("Pre-renumber merged labels: ", paste(pre_renumber_levels, collapse = ", ")),
  paste0("Original-cluster marker check settings: min.pct=", marker_min_pct, ", logfc.threshold=", marker_logfc_threshold, ", Seurat default test."),
  paste0("cluster_final DEG settings: min.pct=", deg_min_pct, ", logfc.threshold=", deg_logfc_threshold, ", Seurat default test."),
  paste0("cluster_final DEG parallel workers requested: ", deg_parallel_workers),
  paste0("cluster_final DEG parallel workers used: ", deg_workers_used),
  paste0(
    "Cosine feature rule: padj < ", cosine_feature_spec$padj_max,
    ", abs(log2FC) >= ", cosine_feature_spec$abs_log2fc_min,
    ", abs(pct.1 - pct.2) >= ", cosine_feature_spec$abs_delta_pct_min,
    ", top ", cosine_feature_spec$top_n_up, " up + top ", cosine_feature_spec$top_n_down, " down genes."
  ),
  "",
  "Merged cluster membership:",
  paste0(
    "  final ",
    manual_merge_membership_df$cluster_final,
    " (pre ",
    manual_merge_membership_df$cluster_final_pre_renumber,
    "): ",
    manual_merge_membership_df$source_clusters
  ),
  "",
  paste0("Saved merged object: ", output_rds)
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("03b manual cluster merge finished.")
message("Outputs written to: ", normalizePath(output_root, mustWork = FALSE))
