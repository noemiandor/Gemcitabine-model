#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

cmd_args <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("pseudotimeAssociations.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))

parse_args <- function(args) {
  out <- list(input = NULL, output_dir = NULL)
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--input=")) {
      out$input <- sub("^--input=", "", arg)
    } else if (arg == "--input") {
      i <- i + 1
      out$input <- args[[i]]
    } else if (startsWith(arg, "--output-dir=")) {
      out$output_dir <- sub("^--output-dir=", "", arg)
    } else if (arg == "--output-dir") {
      i <- i + 1
      out$output_dir <- args[[i]]
    } else {
      stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
    }
    i <- i + 1
  }
  out
}

args <- parse_args(commandArgs(TRUE))
default_output_dir <- file.path(
  repo_root,
  "Results",
  "in-vivo",
  "pseudotime_associations",
  "runs",
  paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_pseudotime_associations")
)
input_csv <- if (!is.null(args$input)) args$input else file.path(
  repo_root,
  "Data",
  "in-vivo",
  "CellCycelCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
)
out_dir <- if (!is.null(args$output_dir)) args$output_dir else default_output_dir
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

df <- read.csv(input_csv, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c(
  "sample_id",
  "initial_ploidy",
  "gemcitabine_dose",
  "gemcitabine_dose_mg_per_kg",
  "pseudotime",
  "cell_ploidy",
  "TGI_percent_auc"
)
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

df <- df[is.finite(df$pseudotime), ]
sample_ids <- sort(unique(df$sample_id))
sample_meta <- do.call(
  rbind,
  lapply(sample_ids, function(sample_id) {
    idx <- df$sample_id == sample_id
    meta <- unique(df[idx, c(
      "sample_id",
      "initial_ploidy",
      "gemcitabine_dose",
      "gemcitabine_dose_mg_per_kg",
      "TGI_percent_auc"
    )])
    if (nrow(meta) != 1) {
      stop("Sample has inconsistent metadata: ", sample_id)
    }
    data.frame(
      sample_id = sample_id,
      initial_ploidy = meta$initial_ploidy,
      dose = meta$gemcitabine_dose,
      dose_mg = meta$gemcitabine_dose_mg_per_kg,
      TGI_percent_auc = meta$TGI_percent_auc,
      mean_pseudotime = mean(df$pseudotime[idx]),
      mean_cell_ploidy = mean(df$cell_ploidy[idx], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
rownames(sample_meta) <- sample_meta$sample_id
sample_meta <- sample_meta[order(sample_meta$dose_mg, sample_meta$initial_ploidy, sample_meta$sample_id), ]

dose_levels <- unique(sample_meta$dose[order(sample_meta$dose_mg)])
dose_cols <- c("0mg/kg" = "#8c8c8c", "30mg/kg" = "#fdae61", "120mg/kg" = "#7b3294")
dose_cols <- dose_cols[names(dose_cols) %in% dose_levels]
plot_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

ecdf_grid <- sort(unique(df$pseudotime))

ecdf_distance_test <- function(comparison_df, group_col, group_a, group_b, stratify_by_ploidy = FALSE) {
  local_samples <- sort(unique(comparison_df$sample_id))
  local_meta <- unique(comparison_df[, c("sample_id", "initial_ploidy", group_col)])
  groups <- setNames(as.character(local_meta[[group_col]]), local_meta$sample_id)[local_samples]
  ploidy <- setNames(as.character(local_meta$initial_ploidy), local_meta$sample_id)[local_samples]

  local_grid <- sort(unique(comparison_df$pseudotime))
  local_ecdf_mat <- t(vapply(local_samples, function(sample_id) {
    stats::ecdf(comparison_df$pseudotime[comparison_df$sample_id == sample_id])(local_grid)
  }, numeric(length(local_grid))))

  stat_fun <- function(groups_now) {
    mean((
      colMeans(local_ecdf_mat[groups_now == group_a, , drop = FALSE]) -
        colMeans(local_ecdf_mat[groups_now == group_b, , drop = FALSE])
    )^2)
  }

  observed <- stat_fun(groups)
  if (!stratify_by_ploidy) {
    choices <- combn(seq_along(local_samples), sum(groups == group_a))
    perm_stats <- apply(choices, 2, function(group_a_idx) {
      perm_groups <- rep(group_b, length(local_samples))
      perm_groups[group_a_idx] <- group_a
      stat_fun(perm_groups)
    })
  } else {
    choices_by_ploidy <- lapply(split(seq_along(local_samples), ploidy), function(idx) {
      combn(idx, sum(groups[idx] == group_a), simplify = FALSE)
    })
    choice_grid <- expand.grid(lapply(choices_by_ploidy, seq_along), KEEP.OUT.ATTRS = FALSE)
    perm_stats <- apply(choice_grid, 1, function(choice_idx) {
      group_a_idx <- unlist(
        mapply(
          function(choices, idx) choices[[idx]],
          choices_by_ploidy,
          as.integer(choice_idx),
          SIMPLIFY = FALSE
        ),
        use.names = FALSE
      )
      perm_groups <- rep(group_b, length(local_samples))
      perm_groups[group_a_idx] <- group_a
      stat_fun(perm_groups)
    })
  }

  list(statistic = observed, p_value = mean(perm_stats >= observed - 1e-15))
}

df$combined_dose_group <- ifelse(df$gemcitabine_dose_mg_per_kg == 0, "0mg/kg", "30+120mg/kg")
combined_test <- ecdf_distance_test(df, "combined_dose_group", "0mg/kg", "30+120mg/kg")
combined_stratified_test <- ecdf_distance_test(df, "combined_dose_group", "0mg/kg", "30+120mg/kg", TRUE)

pairwise_dose_test <- function(dose_a, dose_b) {
  comparison_df <- df[df$gemcitabine_dose_mg_per_kg %in% c(dose_a, dose_b), ]
  comparison_df$comparison_group <- as.character(comparison_df$gemcitabine_dose)
  ecdf_distance_test(comparison_df, "comparison_group", paste0(dose_a, "mg/kg"), paste0(dose_b, "mg/kg"))
}
test_0_vs_30 <- pairwise_dose_test(0, 30)
test_0_vs_120 <- pairwise_dose_test(0, 120)
test_30_vs_120 <- pairwise_dose_test(30, 120)

shift_metrics <- do.call(
  rbind,
  lapply(sample_ids, function(sample_id) {
    sample_ploidy <- sample_meta[sample_id, "initial_ploidy"]
    sample_dose <- sample_meta[sample_id, "dose_mg"]
    controls <- sample_meta$sample_id[
      sample_meta$initial_ploidy == sample_ploidy &
        sample_meta$dose_mg == 0
    ]
    reference_samples <- if (sample_dose == 0) setdiff(controls, sample_id) else controls
    sample_pseudotime <- df$pseudotime[df$sample_id == sample_id]
    reference_pseudotime <- df$pseudotime[df$sample_id %in% reference_samples]
    ecdf_delta <- stats::ecdf(sample_pseudotime)(ecdf_grid) -
      stats::ecdf(reference_pseudotime)(ecdf_grid)
    data.frame(
      sample_id = sample_id,
      initial_ploidy = sample_ploidy,
      dose = sample_meta[sample_id, "dose"],
      dose_mg = sample_dose,
      TGI_percent_auc = sample_meta[sample_id, "TGI_percent_auc"],
      mean_cell_ploidy = sample_meta[sample_id, "mean_cell_ploidy"],
      ecdf_rmse = sqrt(mean(ecdf_delta^2)),
      stringsAsFactors = FALSE
    )
  })
)
shift_metrics <- shift_metrics[order(shift_metrics$dose_mg, shift_metrics$initial_ploidy, shift_metrics$sample_id), ]

treated_shift_metrics <- shift_metrics[shift_metrics$dose_mg > 0, ]
auc_pearson <- suppressWarnings(stats::cor.test(
  treated_shift_metrics$ecdf_rmse,
  treated_shift_metrics$TGI_percent_auc,
  method = "pearson"
))
auc_spearman <- suppressWarnings(stats::cor.test(
  treated_shift_metrics$ecdf_rmse,
  treated_shift_metrics$TGI_percent_auc,
  method = "spearman",
  exact = FALSE
))
mean_ploidy_all <- suppressWarnings(stats::cor.test(
  shift_metrics$ecdf_rmse,
  shift_metrics$mean_cell_ploidy,
  method = "pearson"
))
mean_ploidy_treated <- suppressWarnings(stats::cor.test(
  treated_shift_metrics$ecdf_rmse,
  treated_shift_metrics$mean_cell_ploidy,
  method = "pearson"
))
initial_ploidy_all <- suppressWarnings(stats::wilcox.test(
  ecdf_rmse ~ initial_ploidy,
  data = shift_metrics,
  exact = FALSE
))
initial_ploidy_treated <- suppressWarnings(stats::wilcox.test(
  ecdf_rmse ~ initial_ploidy,
  data = treated_shift_metrics,
  exact = FALSE
))

reported_stats <- data.frame(
  statistic = c(
    "0_vs_30plus120_sample_aware_ECDF_p",
    "0_vs_30plus120_ploidy_stratified_ECDF_p",
    "0_vs_30_sample_aware_ECDF_p",
    "0_vs_120_sample_aware_ECDF_p",
    "30_vs_120_sample_aware_ECDF_p",
    "treated_TGI_AUC_vs_ECDF_RMSE_Pearson_r",
    "treated_TGI_AUC_vs_ECDF_RMSE_Pearson_p",
    "treated_TGI_AUC_vs_ECDF_RMSE_Spearman_rho",
    "treated_TGI_AUC_vs_ECDF_RMSE_Spearman_p",
    "all_samples_mean_cell_ploidy_vs_ECDF_RMSE_Pearson_p",
    "treated_mean_cell_ploidy_vs_ECDF_RMSE_Pearson_p",
    "all_samples_initial_ploidy_vs_ECDF_RMSE_Wilcoxon_p",
    "treated_initial_ploidy_vs_ECDF_RMSE_Wilcoxon_p"
  ),
  value = c(
    combined_test$p_value,
    combined_stratified_test$p_value,
    test_0_vs_30$p_value,
    test_0_vs_120$p_value,
    test_30_vs_120$p_value,
    unname(auc_pearson$estimate),
    auc_pearson$p.value,
    unname(auc_spearman$estimate),
    auc_spearman$p.value,
    mean_ploidy_all$p.value,
    mean_ploidy_treated$p.value,
    initial_ploidy_all$p.value,
    initial_ploidy_treated$p.value
  )
)
write.csv(
  reported_stats,
  file.path(table_dir, "pseudotimeAssociations_reported_stats.csv"),
  row.names = FALSE
)

df$dose <- factor(df$gemcitabine_dose, levels = dose_levels)
df$initial_ploidy <- factor(df$initial_ploidy, levels = unique(sample_meta$initial_ploidy))
sample_meta$dose <- factor(sample_meta$dose, levels = dose_levels)
shift_metrics$dose <- factor(shift_metrics$dose, levels = dose_levels)
treated_shift_metrics$dose <- factor(treated_shift_metrics$dose, levels = dose_levels[dose_levels != "0mg/kg"])

density_plot <- ggplot(df, aes(x = pseudotime, color = dose, fill = dose)) +
  geom_density(alpha = 0.12, linewidth = 0.8) +
  facet_wrap(~ initial_ploidy, ncol = 1) +
  scale_color_manual(values = dose_cols, name = "dose") +
  scale_fill_manual(values = dose_cols, name = "dose") +
  labs(
    title = "Pseudotime distributions by dose and initial ploidy",
    x = "pseudotime",
    y = "density"
  ) +
  plot_theme

ecdf_plot_df <- do.call(
  rbind,
  lapply(sample_ids, function(sample_id) {
    data.frame(
      sample_id = sample_id,
      initial_ploidy = sample_meta[sample_id, "initial_ploidy"],
      dose = sample_meta[sample_id, "dose"],
      pseudotime = ecdf_grid,
      ecdf = stats::ecdf(df$pseudotime[df$sample_id == sample_id])(ecdf_grid),
      stringsAsFactors = FALSE
    )
  })
)
ecdf_plot_df$dose <- factor(ecdf_plot_df$dose, levels = dose_levels)
ecdf_plot <- ggplot(ecdf_plot_df, aes(x = pseudotime, y = ecdf, color = dose, group = sample_id)) +
  geom_line(alpha = 0.55, linewidth = 0.55) +
  facet_wrap(~ initial_ploidy, ncol = 1) +
  scale_color_manual(values = dose_cols, name = "dose") +
  labs(
    title = "Sample ECDFs by dose and initial ploidy",
    x = "pseudotime",
    y = "ECDF"
  ) +
  plot_theme

mean_plot <- ggplot(sample_meta, aes(x = dose, y = mean_pseudotime, color = dose)) +
  geom_boxplot(aes(group = dose), outlier.shape = NA, alpha = 0, color = "grey55", linewidth = 0.45) +
  geom_point(aes(shape = initial_ploidy), size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
  scale_color_manual(values = dose_cols, name = "dose") +
  labs(
    title = "Sample mean pseudotime by dose",
    x = "dose",
    y = "sample mean pseudotime"
  ) +
  plot_theme

treated_cols <- dose_cols[names(dose_cols) != "0mg/kg"]
ecdf_auc_plot <- ggplot(
  treated_shift_metrics,
  aes(x = ecdf_rmse, y = TGI_percent_auc, color = dose)
) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  geom_point(aes(shape = initial_ploidy), size = 2.8) +
  geom_text(aes(label = sample_id), size = 2.5, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
  scale_color_manual(values = treated_cols, name = "dose") +
  labs(
    title = "TGI AUC versus pseudotime distance from ploidy-matched controls",
    x = "ECDF RMSE from ploidy-matched 0 mg/kg reference",
    y = "TGI AUC (%)"
  ) +
  plot_theme

save_plot <- function(filename, plot, width, height) {
  suppressMessages(ggsave(file.path(fig_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300))
  suppressMessages(ggsave(file.path(fig_dir, paste0(filename, ".pdf")), plot, width = width, height = height))
}

save_plot("CellCycleCells_dose_group_comparisons_density_by_dose_and_ploidy", density_plot, 8, 5.5)
save_plot("CellCycleCells_dose_group_comparisons_sample_ecdf_by_dose_and_ploidy", ecdf_plot, 8, 5.5)
save_plot("CellCycleCells_dose_group_comparisons_sample_mean_pseudotime_by_dose", mean_plot, 6.5, 4.8)
save_plot("CellCycleCells_pseudotime_shift_associations_TGI_AUC_ecdf_distance_vs_TGI_AUC", ecdf_auc_plot, 8, 5.5)
