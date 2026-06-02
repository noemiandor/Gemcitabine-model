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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

set.seed(12345)

input_rds <- file.path(
  results_root,
  "03_final_cluster",
  "03_objects",
  "integrated_sct_cca_seurat_final_reclustered.rds"
)
cluster_col <- "seurat_clusters"
input_root <- file.path(results_root, "03_final_cluster", "01_reclustered_cluster_vs_rest_DEG")
output_root <- file.path(results_root, "03a_refine_cluster")

feature_specs <- default_deg_similarity_feature_specs()
metric_names <- c("directional_jaccard", "signed_lfc_cosine", "signed_lfc_spearman")
top_n_positive_highlight <- 5L

message("Preparing output directories.")
.ensure_dir(output_root)
matrix_dir <- .ensure_dir(file.path(output_root, "matrices"))
plot_dir <- .ensure_dir(file.path(output_root, "plots"))
summary_dir <- .ensure_dir(file.path(output_root, "summaries"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}
if (!dir.exists(input_root)) {
  stop("Input DEG directory does not exist: ", input_root, call. = FALSE)
}

message("Reading final Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input RDS is not a Seurat object.", call. = FALSE)
}
if (!(cluster_col %in% colnames(obj@meta.data))) {
  stop("Cluster column is missing from Seurat metadata: ", cluster_col, call. = FALSE)
}

cluster_values <- as.character(obj@meta.data[[cluster_col]])
if (any(is.na(cluster_values) | cluster_values == "")) {
  stop("Missing cluster values detected in metadata column: ", cluster_col, call. = FALSE)
}
cluster_ids <- sort_maybe_numeric(cluster_values)

cluster_count_df <- data.frame(
  cluster = names(table(cluster_values)),
  n_cells = as.integer(table(cluster_values)),
  stringsAsFactors = FALSE
) %>%
  arrange(match(cluster, cluster_ids))
readr::write_csv(cluster_count_df, file.path(summary_dir, "seurat_cluster_cell_counts.csv"))
rm(obj)
gc(verbose = FALSE)

threshold_metadata <- bind_rows(lapply(feature_specs, function(spec) {
  data.frame(
    threshold = spec$name,
    type = spec$type,
    description = spec$description,
    stringsAsFactors = FALSE
  )
}))
readr::write_csv(threshold_metadata, file.path(summary_dir, "threshold_metadata.csv"))

cluster_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)
cluster_dirs <- cluster_dirs[grepl("^cluster_[0-9]+$", basename(cluster_dirs))]
if (length(cluster_dirs) == 0) {
  stop("No cluster directories found in: ", input_root, call. = FALSE)
}

cluster_dirs <- file.path(input_root, paste0("cluster_", cluster_ids))
missing_cluster_dirs <- cluster_dirs[!dir.exists(cluster_dirs)]
if (length(missing_cluster_dirs) > 0) {
  stop("Missing DEG cluster directories: ", paste(missing_cluster_dirs, collapse = ", "), call. = FALSE)
}

extra_deg_clusters <- setdiff(sub("^cluster_", "", basename(list.dirs(input_root, recursive = FALSE, full.names = TRUE))), cluster_ids)
extra_deg_clusters <- extra_deg_clusters[extra_deg_clusters != basename(input_root)]
if (length(extra_deg_clusters) > 0) {
  message("Ignoring DEG directories not present in ", cluster_col, ": ", paste(sort_maybe_numeric(extra_deg_clusters), collapse = ", "))
}

message("Reading cluster DEG tables.")
deg_tables <- setNames(lapply(cluster_dirs, read_cluster_deg), cluster_ids)

message("Building feature sets under multiple thresholds.")
feature_sizes <- list()
pairwise_all <- list()
nearest_neighbor_all <- list()
summary_all <- list()

for (spec in feature_specs) {
  threshold_name <- spec$name
  message("Processing threshold: ", threshold_name)

  feature_map <- setNames(vector("list", length(cluster_ids)), cluster_ids)
  for (cluster_id in cluster_ids) {
    feature_map[[cluster_id]] <- make_feature_object(deg_tables[[cluster_id]], spec)
  }

  feature_sizes[[threshold_name]] <- bind_rows(lapply(cluster_ids, function(cluster_id) {
    fo <- feature_map[[cluster_id]]
    data.frame(
      threshold = threshold_name,
      cluster = cluster_id,
      n_features = nrow(fo$table),
      n_up = length(fo$up),
      n_down = length(fo$down),
      stringsAsFactors = FALSE
    )
  }))

  pair_df <- compute_pair_metrics(cluster_ids, feature_map, threshold_name)
  pairwise_all[[threshold_name]] <- pair_df

  for (metric_name in metric_names) {
    mat <- pair_to_matrix(pair_df, cluster_ids, metric_name)
    write_matrix_csv(mat, file.path(matrix_dir, paste0("similarity_matrix__", metric_name, "__", threshold_name, ".csv")))
    plot_similarity_heatmap(
      mat = mat,
      title = paste0(metric_name, " | ", threshold_name),
      file_pdf = file.path(plot_dir, paste0("heatmap__", metric_name, "__", threshold_name, ".pdf")),
      na_to_zero = TRUE
    )

    nearest_neighbor_all[[paste(metric_name, threshold_name, sep = "__")]] <- bind_rows(lapply(cluster_ids, function(cluster_id) {
      values <- mat[cluster_id, , drop = TRUE]
      values <- values[names(values) != cluster_id]
      if (length(values) == 0 || all(is.na(values))) {
        return(data.frame(
          threshold = threshold_name,
          metric = metric_name,
          cluster = cluster_id,
          nearest_cluster = NA_character_,
          similarity = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      best_idx <- which.max(values)
      data.frame(
        threshold = threshold_name,
        metric = metric_name,
        cluster = cluster_id,
        nearest_cluster = names(values)[best_idx],
        similarity = as.numeric(values[best_idx]),
        stringsAsFactors = FALSE
      )
    }))

    off_diag_values <- mat[row(mat) != col(mat)]
    off_diag_values <- as.numeric(off_diag_values[is.finite(off_diag_values)])
    summary_all[[paste(metric_name, threshold_name, sep = "__")]] <- data.frame(
      threshold = threshold_name,
      metric = metric_name,
      mean_offdiag_similarity = if (length(off_diag_values) > 0) mean(off_diag_values) else NA_real_,
      median_offdiag_similarity = if (length(off_diag_values) > 0) median(off_diag_values) else NA_real_,
      sd_offdiag_similarity = if (length(off_diag_values) > 1) stats::sd(off_diag_values) else NA_real_,
      max_offdiag_similarity = if (length(off_diag_values) > 0) max(off_diag_values) else NA_real_,
      min_offdiag_similarity = if (length(off_diag_values) > 0) min(off_diag_values) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}

feature_size_df <- bind_rows(feature_sizes)
readr::write_csv(feature_size_df, file.path(summary_dir, "feature_set_sizes_by_threshold.csv"))

pairwise_long <- bind_rows(pairwise_all) %>%
  tidyr::pivot_longer(
    cols = all_of(metric_names),
    names_to = "metric",
    values_to = "similarity"
  ) %>%
  filter(cluster_1 != cluster_2) %>%
  mutate(
    pair_id = paste(pmin(cluster_1, cluster_2), pmax(cluster_1, cluster_2), sep = "__"),
    threshold = factor(threshold, levels = vapply(feature_specs, `[[`, character(1), "name"))
  ) %>%
  distinct(threshold, metric, pair_id, .keep_all = TRUE)

readr::write_csv(pairwise_long, file.path(summary_dir, "pairwise_similarity_long.csv"))

nearest_neighbor_df <- bind_rows(nearest_neighbor_all)
readr::write_csv(nearest_neighbor_df, file.path(summary_dir, "nearest_neighbor_by_threshold_metric.csv"))

summary_df <- bind_rows(summary_all)
readr::write_csv(summary_df, file.path(summary_dir, "threshold_metric_summary.csv"))

pair_rank_summary <- pairwise_long %>%
  group_by(metric, pair_id) %>%
  summarise(
    n_thresholds = sum(!is.na(similarity)),
    mean_similarity = mean(similarity, na.rm = TRUE),
    median_similarity = median(similarity, na.rm = TRUE),
    max_similarity = max(similarity, na.rm = TRUE),
    min_similarity = min(similarity, na.rm = TRUE),
    range_similarity = max_similarity - min_similarity,
    sd_similarity = if (sum(!is.na(similarity)) > 1) stats::sd(similarity, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  separate(pair_id, into = c("cluster_1", "cluster_2"), sep = "__", remove = FALSE)

readr::write_csv(pair_rank_summary, file.path(summary_dir, "pair_similarity_across_thresholds_summary.csv"))

stable_positive_summary <- pairwise_long %>%
  group_by(metric, pair_id) %>%
  summarise(
    n_thresholds = sum(!is.na(similarity)),
    min_similarity = min(similarity, na.rm = TRUE),
    max_similarity = max(similarity, na.rm = TRUE),
    mean_similarity = mean(similarity, na.rm = TRUE),
    median_similarity = median(similarity, na.rm = TRUE),
    range_similarity = max_similarity - min_similarity,
    sd_similarity = if (sum(!is.na(similarity)) > 1) stats::sd(similarity, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    always_positive = !is.na(min_similarity) & min_similarity > 0
  ) %>%
  group_by(metric) %>%
  mutate(
    mean_positive_median = {
      positive_means <- mean_similarity[always_positive]
      if (length(positive_means) > 0) stats::median(positive_means, na.rm = TRUE) else NA_real_
    },
    stability_score = if_else(always_positive, mean_similarity / pmax(range_similarity, 1e-6), NA_real_),
    eligible_highlight = always_positive & !is.na(mean_positive_median) & mean_similarity >= mean_positive_median
  ) %>%
  arrange(desc(eligible_highlight), desc(stability_score), desc(mean_similarity), range_similarity, pair_id, .by_group = TRUE) %>%
  mutate(
    highlight_rank = if_else(eligible_highlight, cumsum(eligible_highlight), NA_integer_),
    highlight_pair = eligible_highlight & highlight_rank <= top_n_positive_highlight
  ) %>%
  ungroup() %>%
  separate(pair_id, into = c("cluster_1", "cluster_2"), sep = "__", remove = FALSE)

readr::write_csv(stable_positive_summary, file.path(summary_dir, "pair_similarity_stable_positive_summary.csv"))

stable_positive_highlight_df <- stable_positive_summary %>%
  filter(highlight_pair) %>%
  select(
    metric,
    pair_id,
    cluster_1,
    cluster_2,
    mean_similarity,
    min_similarity,
    max_similarity,
    range_similarity,
    sd_similarity,
    stability_score,
    highlight_rank
  ) %>%
  arrange(metric, highlight_rank, pair_id)

readr::write_csv(stable_positive_highlight_df, file.path(summary_dir, "stable_positive_high_similarity_pairs.csv"))

stable_positive_plot_df <- pairwise_long %>%
  left_join(
    stable_positive_summary %>%
      select(metric, pair_id, highlight_pair),
    by = c("metric", "pair_id")
  ) %>%
  mutate(
    highlight_pair = dplyr::coalesce(highlight_pair, FALSE)
  )

trend_plot <- ggplot(pairwise_long, aes(x = threshold, y = similarity, group = pair_id, color = pair_id)) +
  geom_line(alpha = 0.35, linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 1.0) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    x = "Feature-set threshold",
    y = "Similarity",
    title = "Cluster-pair similarity trends across DEG feature definitions"
  )

ggsave(file.path(plot_dir, "pairwise_similarity_trends.pdf"), trend_plot, width = 11, height = 10)

stable_positive_trend_plot <- ggplot() +
  geom_line(
    data = stable_positive_plot_df %>% filter(!highlight_pair),
    aes(x = threshold, y = similarity, group = pair_id),
    color = "grey80",
    alpha = 0.6,
    linewidth = 0.35
  ) +
  geom_point(
    data = stable_positive_plot_df %>% filter(!highlight_pair),
    aes(x = threshold, y = similarity, group = pair_id),
    color = "grey80",
    alpha = 0.5,
    size = 0.8
  ) +
  geom_line(
    data = stable_positive_plot_df %>% filter(highlight_pair),
    aes(x = threshold, y = similarity, group = pair_id, color = pair_id),
    alpha = 0.9,
    linewidth = 0.7
  ) +
  geom_point(
    data = stable_positive_plot_df %>% filter(highlight_pair),
    aes(x = threshold, y = similarity, group = pair_id, color = pair_id),
    alpha = 0.95,
    size = 1.5
  ) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y") +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(
    x = "Feature-set threshold",
    y = "Similarity",
    color = "Stable positive pair",
    title = "Cluster-pair similarity trends with stable positive pairs highlighted",
    subtitle = "Highlighted pairs remain positive across thresholds, exceed the within-metric positive median, and rank highest by mean/range"
  )

ggsave(file.path(plot_dir, "pairwise_similarity_trends__stable_positive_highlighted.pdf"), stable_positive_trend_plot, width = 11, height = 10)

feature_size_plot <- ggplot(feature_size_df, aes(x = threshold, y = n_features, group = cluster, color = cluster)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.2) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Feature-set threshold",
    y = "Number of retained feature genes",
    title = "Feature-set size by cluster across DEG thresholds"
  )

ggsave(file.path(plot_dir, "feature_set_sizes_by_threshold.pdf"), feature_size_plot, width = 10, height = 6)

integrated_text <- c(
  "03a refine cluster similarity by existing DEGs workflow completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Cluster column: ", cluster_col),
  paste0("Input root: ", input_root),
  paste0("Output root: ", output_root),
  paste0("Number of clusters processed: ", length(cluster_ids)),
  paste0("Thresholds: ", paste(vapply(feature_specs, `[[`, character(1), "name"), collapse = ", ")),
  paste0("Metrics: ", paste(metric_names, collapse = ", "))
)
writeLines(integrated_text, con = file.path(summary_dir, "run_summary.txt"))

message("Cluster similarity workflow completed.")
message("Output root: ", output_root)
