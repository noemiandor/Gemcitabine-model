#!/usr/bin/env Rscript

script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_match <- grep(file_arg, cmd_args, value = TRUE)
if (length(file_match) > 0) {
  script_path <- normalizePath(sub(file_arg, "", file_match[1]), mustWork = FALSE)
}
if (is.null(script_path) || !nzchar(script_path)) {
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = FALSE) else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    script_path <- frame_files[length(frame_files)]
  }
}
script_dir <- if (!is.null(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(pheatmap)
  library(tibble)
})

set.seed(12345)

input_root <- file.path(results_root, "022_clusterEnrich", "01_cluster_vs_rest_DEG_ORA")
output_root <- file.path(results_root, "02a_cluster_similarity_by_DEGs")

feature_specs <- list(
  list(
    name = "all_ranked",
    type = "all",
    description = "All genes retained with signed avg_log2FC."
  ),
  list(
    name = "lenient",
    type = "filtered",
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    description = "padj < 0.05, abs(avg_log2FC) >= 0.25, abs(pct.1 - pct.2) >= 0.05."
  ),
  list(
    name = "moderate",
    type = "filtered",
    padj_max = 0.01,
    abs_log2fc_min = 0.50,
    abs_delta_pct_min = 0.10,
    description = "padj < 0.01, abs(avg_log2FC) >= 0.50, abs(pct.1 - pct.2) >= 0.10."
  ),
  list(
    name = "strict",
    type = "filtered",
    padj_max = 0.001,
    abs_log2fc_min = 1.00,
    abs_delta_pct_min = 0.20,
    description = "padj < 0.001, abs(avg_log2FC) >= 1.00, abs(pct.1 - pct.2) >= 0.20."
  ),
  list(
    name = "top50_each",
    type = "top_n",
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    top_n_up = 50L,
    top_n_down = 50L,
    description = "After lenient filtering, keep top 50 up and top 50 down genes by avg_log2FC."
  ),
  list(
    name = "top30_each",
    type = "top_n",
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    top_n_up = 30L,
    top_n_down = 30L,
    description = "After lenient filtering, keep top 30 up and top 30 down genes by avg_log2FC."
  ),
  list(
    name = "top20_each",
    type = "top_n",
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    top_n_up = 20L,
    top_n_down = 20L,
    description = "After lenient filtering, keep top 20 up and top 20 down genes by avg_log2FC."
  ),
  list(
    name = "top10_each",
    type = "top_n",
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    top_n_up = 10L,
    top_n_down = 10L,
    description = "After lenient filtering, keep top 10 up and top 10 down genes by avg_log2FC."
  )
)

metric_names <- c("directional_jaccard", "signed_lfc_cosine", "signed_lfc_spearman")
top_n_positive_highlight <- 5L

normalize_gene_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^GRCh[0-9]+[-_]", "", x, ignore.case = TRUE)
  x <- sub("^GRCm39[-_]", "", x, ignore.case = TRUE)
  x <- sub("^hg38[-_]", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

resolve_gene_label <- function(df) {
  gene_symbol <- if ("gene_symbol" %in% colnames(df)) as.character(df$gene_symbol) else rep(NA_character_, nrow(df))
  gene <- if ("gene" %in% colnames(df)) as.character(df$gene) else rep(NA_character_, nrow(df))
  out <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)
  out[is.na(out) | out == ""] <- gene[is.na(out) | out == ""]
  out
}

collapse_deg_table <- function(df) {
  df$gene_label <- resolve_gene_label(df)
  df$gene_key <- normalize_gene_key(df$gene_label)
  df$delta_pct <- as.numeric(df$pct.1) - as.numeric(df$pct.2)
  df$abs_log2fc <- abs(as.numeric(df$avg_log2FC))
  df$abs_delta_pct <- abs(df$delta_pct)
  df$p_val_adj_num <- as.numeric(df$p_val_adj)

  df <- df %>%
    filter(!is.na(gene_key), !is.na(avg_log2FC), !is.na(p_val_adj_num)) %>%
    arrange(desc(abs_log2fc), p_val_adj_num)

  df %>%
    group_by(gene_key) %>%
    slice(1) %>%
    ungroup()
}

read_cluster_deg <- function(cluster_dir) {
  cluster_name <- basename(cluster_dir)
  cluster_id <- sub("^cluster_", "", cluster_name)
  deg_path <- file.path(cluster_dir, "DEG", paste0(cluster_name, "_vs_rest_markers.csv"))

  if (!file.exists(deg_path)) {
    stop("Missing DEG file: ", deg_path, call. = FALSE)
  }

  df <- readr::read_csv(deg_path, show_col_types = FALSE)
  required_cols <- c("avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required DEG columns in ", deg_path, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- collapse_deg_table(df)
  df$cluster <- cluster_id
  df
}

filter_feature_table <- function(df, spec) {
  if (spec$type == "all") {
    keep <- rep(TRUE, nrow(df))
  } else {
    keep <- df$p_val_adj_num < spec$padj_max &
      df$abs_log2fc >= spec$abs_log2fc_min &
      df$abs_delta_pct >= spec$abs_delta_pct_min
  }

  filtered <- df[keep, , drop = FALSE]

  if (spec$type == "top_n") {
    up_tbl <- filtered %>%
      filter(avg_log2FC > 0) %>%
      arrange(desc(avg_log2FC), p_val_adj_num) %>%
      slice_head(n = spec$top_n_up)
    down_tbl <- filtered %>%
      filter(avg_log2FC < 0) %>%
      arrange(avg_log2FC, p_val_adj_num) %>%
      slice_head(n = spec$top_n_down)
    filtered <- bind_rows(up_tbl, down_tbl) %>%
      arrange(desc(abs_log2fc), p_val_adj_num) %>%
      distinct(gene_key, .keep_all = TRUE)
  }

  filtered
}

make_feature_object <- function(df, spec) {
  feature_tbl <- filter_feature_table(df, spec)
  up_tbl <- feature_tbl %>% filter(avg_log2FC > 0)
  down_tbl <- feature_tbl %>% filter(avg_log2FC < 0)

  signed_vec <- feature_tbl$avg_log2FC
  names(signed_vec) <- feature_tbl$gene_key

  list(
    table = feature_tbl,
    up = unique(up_tbl$gene_key),
    down = unique(down_tbl$gene_key),
    signed_vec = signed_vec
  )
}

safe_directional_jaccard <- function(f1, f2) {
  same_up <- safe_jaccard(f1$up, f2$up)
  same_down <- safe_jaccard(f1$down, f2$down)
  cross_up_down <- safe_jaccard(f1$up, f2$down)
  cross_down_up <- safe_jaccard(f1$down, f2$up)

  if (all(is.na(c(same_up, same_down, cross_up_down, cross_down_up)))) {
    return(NA_real_)
  }

  0.5 * dplyr::coalesce(same_up, 0) +
    0.5 * dplyr::coalesce(same_down, 0) -
    0.5 * dplyr::coalesce(cross_up_down, 0) -
    0.5 * dplyr::coalesce(cross_down_up, 0)
}

make_union_vectors <- function(v1, v2) {
  genes <- union(names(v1), names(v2))
  if (length(genes) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }
  x <- setNames(rep(0, length(genes)), genes)
  y <- setNames(rep(0, length(genes)), genes)
  x[names(v1)] <- v1
  y[names(v2)] <- v2
  list(x = unname(x), y = unname(y))
}

safe_cosine <- function(v1, v2) {
  vv <- make_union_vectors(v1, v2)
  x <- vv$x
  y <- vv$y
  if (length(x) == 0) return(NA_real_)
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(x * y) / denom
}

safe_spearman <- function(v1, v2) {
  vv <- make_union_vectors(v1, v2)
  x <- vv$x
  y <- vv$y
  if (length(x) < 3) return(NA_real_)
  if (length(unique(x)) < 2 || length(unique(y)) < 2) return(NA_real_)
  suppressWarnings(stats::cor(x, y, method = "spearman", use = "pairwise.complete.obs"))
}

compute_pair_metrics <- function(cluster_ids, feature_map, threshold_name) {
  out <- list()
  idx <- 1L

  for (i in seq_len(length(cluster_ids))) {
    for (j in seq(i, length(cluster_ids))) {
      c1 <- cluster_ids[i]
      c2 <- cluster_ids[j]

      if (c1 == c2) {
        directional_jaccard <- 1
        signed_lfc_cosine <- 1
        signed_lfc_spearman <- 1
      } else {
        f1 <- feature_map[[c1]]
        f2 <- feature_map[[c2]]
        directional_jaccard <- safe_directional_jaccard(f1, f2)
        signed_lfc_cosine <- safe_cosine(f1$signed_vec, f2$signed_vec)
        signed_lfc_spearman <- safe_spearman(f1$signed_vec, f2$signed_vec)
      }

      out[[idx]] <- data.frame(
        threshold = threshold_name,
        cluster_1 = c1,
        cluster_2 = c2,
        directional_jaccard = directional_jaccard,
        signed_lfc_cosine = signed_lfc_cosine,
        signed_lfc_spearman = signed_lfc_spearman,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  bind_rows(out)
}

pair_to_matrix <- function(pair_df, cluster_ids, metric_name) {
  mat <- matrix(NA_real_, nrow = length(cluster_ids), ncol = length(cluster_ids), dimnames = list(cluster_ids, cluster_ids))
  for (i in seq_len(nrow(pair_df))) {
    c1 <- as.character(pair_df$cluster_1[i])
    c2 <- as.character(pair_df$cluster_2[i])
    value <- as.numeric(pair_df[[metric_name]][i])
    mat[c1, c2] <- value
    mat[c2, c1] <- value
  }
  diag(mat) <- 1
  mat
}

write_matrix_csv <- function(mat, file_path) {
  df <- as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column("cluster")
  readr::write_csv(df, file_path)
}

plot_similarity_heatmap <- function(mat, title, file_pdf) {
  pdf(file_pdf, width = 8, height = 7)
  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(c("#1f4e79", "white", "#b03a2e"))(101),
    breaks = seq(-1, 1, length.out = 102),
    border_color = NA,
    na_col = "grey90",
    main = title
  )
  dev.off()
}

message("Preparing output directories.")
.ensure_dir(output_root)
matrix_dir <- .ensure_dir(file.path(output_root, "matrices"))
plot_dir <- .ensure_dir(file.path(output_root, "plots"))
summary_dir <- .ensure_dir(file.path(output_root, "summaries"))

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

cluster_ids <- sub("^cluster_", "", basename(cluster_dirs))
cluster_ids <- sort_maybe_numeric(cluster_ids)
cluster_dirs <- file.path(input_root, paste0("cluster_", cluster_ids))

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
      file_pdf = file.path(plot_dir, paste0("heatmap__", metric_name, "__", threshold_name, ".pdf"))
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
  "Cluster similarity by DEGs workflow completed.",
  paste0("Input root: ", input_root),
  paste0("Output root: ", output_root),
  paste0("Number of clusters processed: ", length(cluster_ids)),
  paste0("Thresholds: ", paste(vapply(feature_specs, `[[`, character(1), "name"), collapse = ", ")),
  paste0("Metrics: ", paste(metric_names, collapse = ", "))
)
writeLines(integrated_text, con = file.path(summary_dir, "run_summary.txt"))

message("Cluster similarity workflow completed.")
message("Output root: ", output_root)
