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

# Cell-cycle diagnostic workflow is implemented for an existing Seurat object.

required_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "pheatmap",
  "patchwork"
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
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(pheatmap)
  library(patchwork)
})

set.seed(12345)

input_rds <- file.path(results_root, "01_data", "integrated_sct_cca_seurat.rds")
output_dir <- file.path(results_root, "01a_cell_cycle")
updated_rds_path <- file.path(output_dir, "integrated_sct_cca_seurat_cell_cycle.rds")
methods_long_path <- file.path(output_dir, "candidate_cell_cycle_cluster_methods_long.csv")
integrated_candidate_path <- file.path(output_dir, "candidate_cell_cycle_cluster_integrated_summary.csv")

phase_levels <- c("G1", "S", "G2M")

.ensure_dir(output_dir)

save_plot_both <- function(plot_obj, file_prefix, out_dir, width = 8, height = 6, dpi = 300) {
  pdf_path <- file.path(out_dir, paste0(file_prefix, ".pdf"))
  png_path <- file.path(out_dir, paste0(file_prefix, ".png"))
  ggplot2::ggsave(filename = pdf_path, plot = plot_obj, width = width, height = height, units = "in")
  ggplot2::ggsave(filename = png_path, plot = plot_obj, width = width, height = height, units = "in", dpi = dpi)
}

na_to_group <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "NA"
  x
}

stable_levels <- function(x) {
  unique(as.character(x))
}

normalize_feature_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^GRCh[0-9]+-", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

extract_feature_metadata <- function(assay_obj) {
  assay_slots <- methods::slotNames(assay_obj)

  if ("meta.features" %in% assay_slots) {
    return(as.data.frame(assay_obj@meta.features, stringsAsFactors = FALSE))
  }

  if ("meta.data" %in% assay_slots) {
    return(as.data.frame(assay_obj@meta.data, stringsAsFactors = FALSE))
  }

  data.frame(row.names = rownames(assay_obj))
}

build_feature_lookup <- function(feature_ids, alias_values) {
  alias_keys <- normalize_feature_key(alias_values)
  keep <- !is.na(alias_keys) & alias_keys != "" & !duplicated(alias_keys)
  stats::setNames(feature_ids[keep], alias_keys[keep])
}

fill_feature_matches <- function(current_matches, target_keys, lookup) {
  proposed <- unname(lookup[target_keys])
  fill_idx <- is.na(current_matches) & !is.na(proposed)
  current_matches[fill_idx] <- proposed[fill_idx]
  current_matches
}

resolve_gene_set_features <- function(seu, assay_name, gene_symbols, set_name, min_features = 5) {
  assay_obj <- seu[[assay_name]]
  feature_ids <- rownames(assay_obj)

  if (length(feature_ids) == 0) {
    stop("No features are present in assay '", assay_name, "'.", call. = FALSE)
  }

  gene_keys <- normalize_feature_key(gene_symbols)
  matched_features <- rep(NA_character_, length(gene_symbols))

  lookup_list <- list(
    feature_ids = build_feature_lookup(feature_ids, feature_ids)
  )

  if (any(grepl("\\|", feature_ids))) {
    lookup_list$feature_before_pipe <- build_feature_lookup(feature_ids, sub("\\|.*$", "", feature_ids))
    lookup_list$feature_after_pipe <- build_feature_lookup(feature_ids, sub("^.*\\|", "", feature_ids))
  }

  if (any(grepl("_", feature_ids))) {
    lookup_list$feature_before_underscore <- build_feature_lookup(feature_ids, sub("_.*$", "", feature_ids))
    lookup_list$feature_after_underscore <- build_feature_lookup(feature_ids, sub("^.*_", "", feature_ids))
  }

  feature_meta <- extract_feature_metadata(assay_obj)
  candidate_meta_cols <- intersect(
    c(
      "gene_name",
      "gene",
      "symbol",
      "gene_symbol",
      "gene_symbols",
      "SYMBOL",
      "GENE",
      "Gene",
      "GeneName",
      "GeneSymbol",
      "feature_name",
      "feature",
      "features"
    ),
    colnames(feature_meta)
  )

  if (length(candidate_meta_cols) > 0) {
    for (col_name in candidate_meta_cols) {
      lookup_list[[paste0("meta_", col_name)]] <- build_feature_lookup(feature_ids, feature_meta[[col_name]])
    }
  }

  for (lookup in lookup_list) {
    matched_features <- fill_feature_matches(matched_features, gene_keys, lookup)
  }

  matched_unique <- unique(matched_features[!is.na(matched_features)])
  missing_symbols <- gene_symbols[is.na(matched_features)]

  message(
    set_name,
    ": ",
    length(matched_unique),
    "/",
    length(gene_symbols),
    " genes matched in assay '",
    assay_name,
    "'."
  )

  if (length(matched_unique) < min_features) {
    stop(
      set_name,
      " gene set has only ",
      length(matched_unique),
      " matched features in assay '",
      assay_name,
      "'. Feature names may not be stored as standard gene symbols.",
      call. = FALSE
    )
  }

  list(
    features = matched_unique,
    missing_symbols = missing_symbols
  )
}

write_phase_table <- function(df, id_col, file_path) {
  out_df <- as.data.frame(df, stringsAsFactors = FALSE)
  readr::write_csv(out_df, file_path)
}

make_phase_tables <- function(group_values, phase_values, group_levels, group_name) {
  tmp_df <- data.frame(
    group = factor(as.character(group_values), levels = group_levels),
    Phase = factor(as.character(phase_values), levels = phase_levels),
    stringsAsFactors = FALSE
  )

  count_long <- dplyr::count(tmp_df, group, Phase, .drop = FALSE, name = "n")
  count_wide <- tidyr::pivot_wider(
    count_long,
    names_from = Phase,
    values_from = n,
    values_fill = 0
  )
  count_wide <- as.data.frame(count_wide, stringsAsFactors = FALSE)
  colnames(count_wide)[1] <- group_name

  count_mat <- as.matrix(count_wide[, phase_levels, drop = FALSE])
  row_totals <- rowSums(count_mat)
  denom <- ifelse(row_totals == 0, 1, row_totals)
  prop_mat <- sweep(count_mat, 1, denom, "/")
  prop_mat[row_totals == 0, ] <- 0

  prop_wide <- data.frame(
    count_wide[[group_name]],
    prop_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  colnames(prop_wide)[1] <- group_name

  prop_long <- tidyr::pivot_longer(
    prop_wide,
    cols = all_of(phase_levels),
    names_to = "Phase",
    values_to = "proportion"
  )

  list(
    counts = count_wide,
    props = prop_wide,
    prop_long = prop_long
  )
}

make_score_distribution_plot <- function(group_values, score_values, group_levels, x_label, y_label) {
  plot_df <- data.frame(
    group = factor(as.character(group_values), levels = group_levels),
    value = as.numeric(score_values),
    stringsAsFactors = FALSE
  )

  ggplot(plot_df, aes(x = group, y = value, fill = group)) +
    geom_violin(scale = "width", trim = TRUE) +
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = x_label, y = y_label)
}

message("Loading Seurat object.")
seu <- readRDS(input_rds)

if (!inherits(seu, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

if (!("seurat_clusters" %in% colnames(seu@meta.data))) {
  stop("Metadata column 'seurat_clusters' is missing.", call. = FALSE)
}

if (!("sample" %in% colnames(seu@meta.data))) {
  stop("Metadata column 'sample' is missing.", call. = FALSE)
}

if (!("Dose" %in% colnames(seu@meta.data))) {
  stop("Metadata column 'Dose' is missing.", call. = FALSE)
}

if (!("orig.ident" %in% colnames(seu@meta.data))) {
  stop("Metadata column 'orig.ident' is missing.", call. = FALSE)
}

if (!("RNA" %in% names(seu@assays))) {
  stop("Assay 'RNA' is missing.", call. = FALSE)
}

if (!("pca" %in% names(seu@reductions))) {
  stop("Reduction 'pca' is missing.", call. = FALSE)
}

if (!("umap" %in% names(seu@reductions))) {
  stop("Reduction 'umap' is missing.", call. = FALSE)
}

message("Running cell-cycle scoring on the RNA assay.")
DefaultAssay(seu) <- "RNA"
cc_genes <- Seurat::cc.genes.updated.2019
s_gene_info <- resolve_gene_set_features(seu, "RNA", cc_genes$s.genes, "S.Score")
g2m_gene_info <- resolve_gene_set_features(seu, "RNA", cc_genes$g2m.genes, "G2M.Score")
seu <- Seurat::CellCycleScoring(
  object = seu,
  s.features = s_gene_info$features,
  g2m.features = g2m_gene_info$features,
  set.ident = FALSE
)

meta_df <- seu@meta.data
meta_df$cluster <- as.character(meta_df$seurat_clusters)
meta_df$sample_group <- na_to_group(meta_df$sample)
meta_df$dose_group <- na_to_group(meta_df$Dose)
meta_df$Phase <- factor(as.character(meta_df$Phase), levels = phase_levels)

cluster_levels <- stable_levels(meta_df$cluster)
sample_levels <- stable_levels(meta_df$sample_group)
dose_levels <- stable_levels(meta_df$dose_group)

meta_df$cluster <- factor(meta_df$cluster, levels = cluster_levels)
meta_df$sample_group <- factor(meta_df$sample_group, levels = sample_levels)
meta_df$dose_group <- factor(meta_df$dose_group, levels = dose_levels)

message("Building summary tables.")

cluster_summary <- meta_df %>%
  group_by(cluster) %>%
  summarise(
    n_cells = n(),
    mean_S.Score = mean(S.Score, na.rm = TRUE),
    median_S.Score = median(S.Score, na.rm = TRUE),
    mean_G2M.Score = mean(G2M.Score, na.rm = TRUE),
    median_G2M.Score = median(G2M.Score, na.rm = TRUE),
    frac_G1 = mean(Phase == "G1", na.rm = TRUE),
    frac_S = mean(Phase == "S", na.rm = TRUE),
    frac_G2M = mean(Phase == "G2M", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cluster = as.character(cluster))

readr::write_csv(cluster_summary, file.path(output_dir, "cluster_cell_cycle_summary.csv"))

cluster_phase_tables <- make_phase_tables(meta_df$cluster, meta_df$Phase, cluster_levels, "cluster")
write_phase_table(cluster_phase_tables$counts, "cluster", file.path(output_dir, "cluster_phase_counts.csv"))
write_phase_table(cluster_phase_tables$props, "cluster", file.path(output_dir, "cluster_phase_row_prop.csv"))

sample_phase_tables <- make_phase_tables(meta_df$sample_group, meta_df$Phase, sample_levels, "sample")
write_phase_table(sample_phase_tables$counts, "sample", file.path(output_dir, "sample_phase_counts.csv"))
write_phase_table(sample_phase_tables$props, "sample", file.path(output_dir, "sample_phase_row_prop.csv"))

dose_phase_tables <- make_phase_tables(meta_df$dose_group, meta_df$Phase, dose_levels, "Dose")
write_phase_table(dose_phase_tables$counts, "Dose", file.path(output_dir, "Dose_phase_counts.csv"))
write_phase_table(dose_phase_tables$props, "Dose", file.path(output_dir, "Dose_phase_row_prop.csv"))

cluster_by_sample_summary <- meta_df %>%
  group_by(cluster, sample = sample_group) %>%
  summarise(
    n_cells = n(),
    mean_S.Score = mean(S.Score, na.rm = TRUE),
    mean_G2M.Score = mean(G2M.Score, na.rm = TRUE),
    frac_G1 = mean(Phase == "G1", na.rm = TRUE),
    frac_S = mean(Phase == "S", na.rm = TRUE),
    frac_G2M = mean(Phase == "G2M", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cluster = as.character(cluster), sample = as.character(sample))
readr::write_csv(cluster_by_sample_summary, file.path(output_dir, "cluster_by_sample_cell_cycle_summary.csv"))

cluster_by_dose_summary <- meta_df %>%
  group_by(cluster, Dose = dose_group) %>%
  summarise(
    n_cells = n(),
    mean_S.Score = mean(S.Score, na.rm = TRUE),
    mean_G2M.Score = mean(G2M.Score, na.rm = TRUE),
    frac_G1 = mean(Phase == "G1", na.rm = TRUE),
    frac_S = mean(Phase == "S", na.rm = TRUE),
    frac_G2M = mean(Phase == "G2M", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cluster = as.character(cluster), Dose = as.character(Dose))
readr::write_csv(cluster_by_dose_summary, file.path(output_dir, "cluster_by_Dose_cell_cycle_summary.csv"))

message("Computing PCA correlations with cell-cycle scores.")
pca_embed <- Embeddings(seu, "pca")
if (ncol(pca_embed) < 20) {
  stop("At least 20 PCs are required in the PCA reduction.", call. = FALSE)
}

pca_use <- pca_embed[, 1:20, drop = FALSE]
if (is.null(colnames(pca_use)) || any(colnames(pca_use) == "")) {
  colnames(pca_use) <- paste0("PC", seq_len(ncol(pca_use)))
}

cor_s <- apply(
  pca_use,
  2,
  function(x) stats::cor(x, meta_df$S.Score, method = "pearson", use = "pairwise.complete.obs")
)
cor_g2m <- apply(
  pca_use,
  2,
  function(x) stats::cor(x, meta_df$G2M.Score, method = "pearson", use = "pairwise.complete.obs")
)

pc_cell_cycle_correlations <- data.frame(
  PC = colnames(pca_use),
  check.names = FALSE,
  cor_with_S.Score = as.numeric(cor_s),
  cor_with_G2M.Score = as.numeric(cor_g2m),
  abs_cor_with_S.Score = abs(as.numeric(cor_s)),
  abs_cor_with_G2M.Score = abs(as.numeric(cor_g2m))
)
readr::write_csv(pc_cell_cycle_correlations, file.path(output_dir, "pc_cell_cycle_correlations.csv"))

pc_corr_plot_df <- tidyr::pivot_longer(
  pc_cell_cycle_correlations[, c("PC", "cor_with_S.Score", "cor_with_G2M.Score")],
  cols = c("cor_with_S.Score", "cor_with_G2M.Score"),
  names_to = "score_metric",
  values_to = "correlation"
)
pc_corr_plot_df$score_metric <- factor(
  pc_corr_plot_df$score_metric,
  levels = c("cor_with_S.Score", "cor_with_G2M.Score"),
  labels = c("S.Score", "G2M.Score")
)
pc_corr_plot_df$PC <- factor(pc_corr_plot_df$PC, levels = rev(colnames(pca_use)))

message("Building diagnostic plots.")

plot_umap_cluster <- DimPlot(seu, reduction = "umap", group.by = "seurat_clusters", raster = FALSE)
save_plot_both(plot_umap_cluster, "umap_by_cluster", output_dir)

plot_umap_phase <- DimPlot(seu, reduction = "umap", group.by = "Phase", raster = FALSE)
save_plot_both(plot_umap_phase, "umap_by_phase", output_dir)

plot_umap_s <- FeaturePlot(seu, reduction = "umap", features = "S.Score", raster = FALSE)
save_plot_both(plot_umap_s, "umap_S_Score", output_dir)

plot_umap_g2m <- FeaturePlot(seu, reduction = "umap", features = "G2M.Score", raster = FALSE)
save_plot_both(plot_umap_g2m, "umap_G2M_Score", output_dir)

plot_vln_s_cluster <- VlnPlot(seu, features = "S.Score", group.by = "seurat_clusters", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_both(plot_vln_s_cluster, "vln_S_Score_by_cluster", output_dir)

plot_vln_g2m_cluster <- VlnPlot(seu, features = "G2M.Score", group.by = "seurat_clusters", pt.size = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot_both(plot_vln_g2m_cluster, "vln_G2M_Score_by_cluster", output_dir)

cluster_mean_scores_plot_df <- tidyr::pivot_longer(
  cluster_summary[, c("cluster", "mean_S.Score", "mean_G2M.Score")],
  cols = c("mean_S.Score", "mean_G2M.Score"),
  names_to = "metric",
  values_to = "mean_score"
)
cluster_mean_scores_plot_df$metric <- factor(
  cluster_mean_scores_plot_df$metric,
  levels = c("mean_S.Score", "mean_G2M.Score"),
  labels = c("mean S.Score", "mean G2M.Score")
)
cluster_mean_scores_plot_df$cluster <- factor(cluster_mean_scores_plot_df$cluster, levels = rev(cluster_levels))

cluster_mean_scores_plot <- ggplot(cluster_mean_scores_plot_df, aes(x = metric, y = cluster, fill = mean_score)) +
  geom_tile() +
  theme_bw(base_size = 11) +
  labs(x = "metric", y = "cluster", fill = "mean score")
save_plot_both(cluster_mean_scores_plot, "cluster_mean_scores_heatmap", output_dir)

plot_s_by_sample <- make_score_distribution_plot(meta_df$sample_group, meta_df$S.Score, sample_levels, "sample", "S.Score")
save_plot_both(plot_s_by_sample, "S_Score_by_sample", output_dir)

plot_g2m_by_sample <- make_score_distribution_plot(meta_df$sample_group, meta_df$G2M.Score, sample_levels, "sample", "G2M.Score")
save_plot_both(plot_g2m_by_sample, "G2M_Score_by_sample", output_dir)

plot_s_by_dose <- make_score_distribution_plot(meta_df$dose_group, meta_df$S.Score, dose_levels, "Dose", "S.Score")
save_plot_both(plot_s_by_dose, "S_Score_by_Dose", output_dir)

plot_g2m_by_dose <- make_score_distribution_plot(meta_df$dose_group, meta_df$G2M.Score, dose_levels, "Dose", "G2M.Score")
save_plot_both(plot_g2m_by_dose, "G2M_Score_by_Dose", output_dir)

phase_by_sample_plot <- ggplot(
  sample_phase_tables$prop_long,
  aes(x = factor(sample, levels = sample_levels), y = proportion, fill = Phase)
) +
  geom_col() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "sample", y = "proportion")
save_plot_both(phase_by_sample_plot, "phase_composition_by_sample", output_dir)

phase_by_dose_plot <- ggplot(
  dose_phase_tables$prop_long,
  aes(x = factor(Dose, levels = dose_levels), y = proportion, fill = Phase)
) +
  geom_col() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Dose", y = "proportion")
save_plot_both(phase_by_dose_plot, "phase_composition_by_Dose", output_dir)

pc_corr_heatmap_plot <- ggplot(pc_corr_plot_df, aes(x = score_metric, y = PC, fill = correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_bw(base_size = 11) +
  labs(x = "score", y = "PC", fill = "correlation")
save_plot_both(pc_corr_heatmap_plot, "pc_cell_cycle_correlation_heatmap", output_dir)

message("Calling candidate cell-cycle-driven clusters.")

safe_sd <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  if (is.na(out)) 0 else out
}

score_threshold_s <- mean(cluster_summary$mean_S.Score, na.rm = TRUE) + safe_sd(cluster_summary$mean_S.Score)
score_threshold_g2m <- mean(cluster_summary$mean_G2M.Score, na.rm = TRUE) + safe_sd(cluster_summary$mean_G2M.Score)

method_score_outlier <- bind_rows(
  data.frame(
    cluster = cluster_summary$cluster,
    method_name = "score_outlier",
    metric_name = "mean_S.Score",
    observed_value = cluster_summary$mean_S.Score,
    threshold_or_rule = sprintf("mean_S.Score > %.6f (global_mean + 1 SD)", score_threshold_s),
    how_judged = "cluster mean S.Score was compared against across-cluster threshold",
    result_flag = cluster_summary$mean_S.Score > score_threshold_s,
    optional_note = ifelse(cluster_summary$mean_S.Score > score_threshold_s, "above threshold", "below threshold"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cluster = cluster_summary$cluster,
    method_name = "score_outlier",
    metric_name = "mean_G2M.Score",
    observed_value = cluster_summary$mean_G2M.Score,
    threshold_or_rule = sprintf("mean_G2M.Score > %.6f (global_mean + 1 SD)", score_threshold_g2m),
    how_judged = "cluster mean G2M.Score was compared against across-cluster threshold",
    result_flag = cluster_summary$mean_G2M.Score > score_threshold_g2m,
    optional_note = ifelse(cluster_summary$mean_G2M.Score > score_threshold_g2m, "above threshold", "below threshold"),
    stringsAsFactors = FALSE
  )
)

flag_score_outlier <- method_score_outlier %>%
  group_by(cluster) %>%
  summarise(flag_score_outlier = any(result_flag), .groups = "drop")

phase_threshold_s <- mean(cluster_summary$frac_S, na.rm = TRUE) + safe_sd(cluster_summary$frac_S)
phase_threshold_g2m <- mean(cluster_summary$frac_G2M, na.rm = TRUE) + safe_sd(cluster_summary$frac_G2M)

method_phase_enrichment <- bind_rows(
  data.frame(
    cluster = cluster_summary$cluster,
    method_name = "phase_enrichment",
    metric_name = "frac_S",
    observed_value = cluster_summary$frac_S,
    threshold_or_rule = sprintf("frac_S > %.6f (global_mean + 1 SD)", phase_threshold_s),
    how_judged = "cluster S-phase fraction was compared against across-cluster threshold",
    result_flag = cluster_summary$frac_S > phase_threshold_s,
    optional_note = ifelse(cluster_summary$frac_S > phase_threshold_s, "above threshold", "below threshold"),
    stringsAsFactors = FALSE
  ),
  data.frame(
    cluster = cluster_summary$cluster,
    method_name = "phase_enrichment",
    metric_name = "frac_G2M",
    observed_value = cluster_summary$frac_G2M,
    threshold_or_rule = sprintf("frac_G2M > %.6f (global_mean + 1 SD)", phase_threshold_g2m),
    how_judged = "cluster G2M-phase fraction was compared against across-cluster threshold",
    result_flag = cluster_summary$frac_G2M > phase_threshold_g2m,
    optional_note = ifelse(cluster_summary$frac_G2M > phase_threshold_g2m, "above threshold", "below threshold"),
    stringsAsFactors = FALSE
  )
)

flag_phase_enrichment <- method_phase_enrichment %>%
  group_by(cluster) %>%
  summarise(flag_phase_enrichment = any(result_flag), .groups = "drop")

pc_selection_score <- pmax(
  pc_cell_cycle_correlations[["abs_cor_with_S.Score"]],
  pc_cell_cycle_correlations[["abs_cor_with_G2M.Score"]]
)
selected_pc_table <- pc_cell_cycle_correlations[
  pc_cell_cycle_correlations[["abs_cor_with_S.Score"]] >= 0.30 |
    pc_cell_cycle_correlations[["abs_cor_with_G2M.Score"]] >= 0.30,
  ,
  drop = FALSE
]

pc_how_judged <- "selected because abs(cor with S.Score) >= 0.30 or abs(cor with G2M.Score) >= 0.30"
if (nrow(selected_pc_table) == 0) {
  selected_pc_table <- pc_cell_cycle_correlations[which.max(pc_selection_score), , drop = FALSE]
  pc_how_judged <- "fallback highest-correlation PC used"
}

selected_pcs <- selected_pc_table$PC

pc_cluster_df <- as.data.frame(pca_use[, selected_pcs, drop = FALSE], check.names = FALSE)
pc_cluster_df$cluster <- as.character(meta_df$cluster)
pc_cluster_long <- tidyr::pivot_longer(
  pc_cluster_df,
  cols = all_of(selected_pcs),
  names_to = "metric_name",
  values_to = "pc_value"
)
pc_cluster_mean <- pc_cluster_long %>%
  group_by(cluster, metric_name) %>%
  summarise(observed_value = mean(pc_value, na.rm = TRUE), .groups = "drop")

pc_cluster_stats <- pc_cluster_mean %>%
  group_by(metric_name) %>%
  summarise(
    across_cluster_mean = mean(observed_value, na.rm = TRUE),
    across_cluster_sd = safe_sd(observed_value),
    .groups = "drop"
  )

method_pc_support <- pc_cluster_mean %>%
  left_join(pc_cluster_stats, by = "metric_name") %>%
  mutate(
    method_name = "pc_separation_support",
    threshold_or_rule = "abs(cluster_mean_PC - across_cluster_mean) > 1 SD on selected cell-cycle-associated PC",
    how_judged = pc_how_judged,
    result_flag = abs(observed_value - across_cluster_mean) > across_cluster_sd,
    optional_note = ifelse(result_flag, "cluster mean PC outlier", "within 1 SD")
  ) %>%
  select(
    cluster,
    method_name,
    metric_name,
    observed_value,
    threshold_or_rule,
    how_judged,
    result_flag,
    optional_note
  )

flag_pc_support <- method_pc_support %>%
  group_by(cluster) %>%
  summarise(flag_pc_separation_support = any(result_flag), .groups = "drop")

candidate_methods_long <- bind_rows(
  method_score_outlier,
  method_phase_enrichment,
  method_pc_support
) %>%
  mutate(cluster = as.character(cluster))

readr::write_csv(candidate_methods_long, methods_long_path)

candidate_integrated <- data.frame(cluster = cluster_levels, stringsAsFactors = FALSE) %>%
  left_join(flag_score_outlier, by = "cluster") %>%
  left_join(flag_phase_enrichment, by = "cluster") %>%
  left_join(flag_pc_support, by = "cluster") %>%
  mutate(
    flag_score_outlier = ifelse(is.na(flag_score_outlier), FALSE, flag_score_outlier),
    flag_phase_enrichment = ifelse(is.na(flag_phase_enrichment), FALSE, flag_phase_enrichment),
    flag_pc_separation_support = ifelse(is.na(flag_pc_separation_support), FALSE, flag_pc_separation_support),
    n_methods_flagged =
      as.integer(flag_score_outlier) +
      as.integer(flag_phase_enrichment) +
      as.integer(flag_pc_separation_support),
    final_candidate_flag = n_methods_flagged >= 2,
    optional_note = dplyr::case_when(
      final_candidate_flag ~ "flagged by >=2 methods",
      n_methods_flagged == 1 ~ "only 1 method flagged",
      TRUE ~ "not flagged"
    )
  ) %>%
  select(
    cluster,
    flag_score_outlier,
    flag_phase_enrichment,
    flag_pc_separation_support,
    n_methods_flagged,
    final_candidate_flag,
    optional_note
  )

readr::write_csv(candidate_integrated, integrated_candidate_path)

cluster_annotation_map <- candidate_integrated %>%
  mutate(
    cluster_cell_cycle_annotation = ifelse(
      final_candidate_flag,
      "cell_cycle_candidate",
      "not_cell_cycle_candidate"
    )
  ) %>%
  select(cluster, cluster_cell_cycle_annotation)

seu$cluster_cell_cycle_annotation <- cluster_annotation_map$cluster_cell_cycle_annotation[
  match(as.character(seu$seurat_clusters), cluster_annotation_map$cluster)
]

message("Saving the Seurat object with cluster-level cell-cycle annotation.")
saveRDS(seu, updated_rds_path)

final_candidate_count <- sum(candidate_integrated$final_candidate_flag, na.rm = TRUE)

message("Cell-cycle diagnostic workflow completed.")
message("Number of cells: ", ncol(seu))
message("Number of clusters: ", length(cluster_levels))
message("Number of final candidate clusters: ", final_candidate_count)
message("Updated RDS saved to: ", updated_rds_path)
message("Candidate methods table saved to: ", methods_long_path)
message("Integrated candidate summary saved to: ", integrated_candidate_path)
