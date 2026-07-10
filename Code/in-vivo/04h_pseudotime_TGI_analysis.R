#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

cmd_args <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]))
} else {
  normalizePath("04h_pseudotime_TGI_analysis.R", mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
data_dir <- file.path(repo_root, "Data", "in-vivo")
results_root <- Sys.getenv(
  "RESULTS_ROOT",
  "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results"
)
output_root <- file.path(results_root, "04h_pseudotime_TGI_analysis")

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  write.csv(x, path, row.names = FALSE)
  invisible(path)
}

write_text <- function(lines, path) {
  ensure_dir(dirname(path))
  writeLines(lines, path)
  invisible(path)
}

save_plot <- function(plot, output_dir, filename, width, height) {
  ensure_dir(output_dir)
  suppressMessages(ggsave(file.path(output_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300))
  suppressMessages(ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, width = width, height = height))
  invisible(filename)
}

format_num <- function(x, digits = 3) {
  ifelse(is.finite(x), trimws(formatC(x, digits = digits, format = "fg")), "NA")
}

format_p <- function(x) {
  ifelse(is.finite(x), ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "fg")), "NA")
}

dose_cols_all <- c("0mg/kg" = "#8c8c8c", "30mg/kg" = "#fdae61", "120mg/kg" = "#7b3294")
plot_theme <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank())

compartments <- list(
  CellCycle = list(
    label = "Cell-cycle-associated tumor cells",
    prefix = "CellCycleCells",
    input = file.path(data_dir, "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    output_dir = file.path(output_root, "CellCycle")
  ),
  NonCellCycle = list(
    label = "Non-cell-cycle-associated tumor cells",
    prefix = "NonCellCycleCells",
    input = file.path(data_dir, "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    output_dir = file.path(output_root, "NonCellCycle")
  )
)

required_cols <- c(
  "sample_id",
  "cluster",
  "initial_ploidy",
  "gemcitabine_dose",
  "gemcitabine_dose_mg_per_kg",
  "pseudotime",
  "cell_ploidy",
  "TGI_percent_auc"
)

endpoint_tgi_cols <- function(nms) {
  cols <- grep("^TGI_percent_Day_[0-9]+$", nms, value = TRUE)
  days <- suppressWarnings(as.numeric(sub("^TGI_percent_Day_", "", cols)))
  cols[order(days)]
}

endpoint_day_num <- function(tgi_col) {
  suppressWarnings(as.numeric(sub("^TGI_percent_Day_", "", tgi_col)))
}

day_cols_by_prefix <- function(nms, prefix) {
  cols <- grep(paste0("^", prefix, "Day_[0-9]+$"), nms, value = TRUE)
  days <- suppressWarnings(as.numeric(sub(paste0("^", prefix, "Day_"), "", cols)))
  cols[order(days)]
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
    return(data.frame(n = sum(ok), estimate = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  test <- if (identical(method, "spearman")) {
    suppressWarnings(stats::cor.test(x[ok], y[ok], method = method, exact = FALSE))
  } else {
    suppressWarnings(stats::cor.test(x[ok], y[ok], method = method))
  }
  data.frame(
    n = sum(ok),
    estimate = unname(test$estimate),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

safe_wilcox <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  group <- as.character(group)
  if (sum(ok) < 3 || length(unique(group[ok])) != 2) {
    return(data.frame(n = sum(ok), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  test <- suppressWarnings(stats::wilcox.test(value[ok] ~ group[ok], exact = FALSE))
  data.frame(
    n = sum(ok),
    statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

safe_paired_test <- function(x, y, method) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) {
    return(data.frame(n = sum(ok), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  test <- if (identical(method, "paired_t")) {
    suppressWarnings(stats::t.test(x[ok], y[ok], paired = TRUE))
  } else {
    suppressWarnings(stats::wilcox.test(x[ok], y[ok], paired = TRUE, exact = FALSE))
  }
  data.frame(
    n = sum(ok),
    statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

read_compartment_data <- function(compartment_name, cfg) {
  if (!file.exists(cfg$input)) {
    stop("Missing input CSV for ", compartment_name, ": ", cfg$input, call. = FALSE)
  }
  df <- read.csv(cfg$input, check.names = FALSE, stringsAsFactors = FALSE)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in ", cfg$input, ": ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  df$pseudotime <- suppressWarnings(as.numeric(df$pseudotime))
  df$cell_ploidy <- suppressWarnings(as.numeric(df$cell_ploidy))
  df$gemcitabine_dose_mg_per_kg <- suppressWarnings(as.numeric(df$gemcitabine_dose_mg_per_kg))
  for (col in grep("^(TGI_percent|tumor_volume|matched_control|auc_n_days)", names(df), value = TRUE)) {
    if (!grepl("_day$|_sample_id|_harvest|_ploidy|_dose", col)) {
      df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
    }
  }
  df <- df[is.finite(df$pseudotime), , drop = FALSE]
  df$compartment <- compartment_name
  df
}

build_sample_meta <- function(df) {
  sample_ids <- sort(unique(df$sample_id))
  sample_level_cols <- setdiff(names(df), c("cell_id", "cluster", "pseudotime", "cell_ploidy", "compartment"))
  rows <- lapply(sample_ids, function(sample_id) {
    sample_df <- df[df$sample_id == sample_id, , drop = FALSE]
    row <- sample_df[1, sample_level_cols, drop = FALSE]
    row$n_cells <- nrow(sample_df)
    row$clusters <- paste(sort(unique(as.character(sample_df$cluster))), collapse = ";")
    row$mean_pseudotime <- mean(sample_df$pseudotime, na.rm = TRUE)
    row$median_pseudotime <- stats::median(sample_df$pseudotime, na.rm = TRUE)
    row$mean_cell_ploidy_exported <- mean(sample_df$cell_ploidy, na.rm = TRUE)
    row$median_cell_ploidy_exported <- stats::median(sample_df$cell_ploidy, na.rm = TRUE)
    row$ploidy_p90 <- as.numeric(stats::quantile(sample_df$cell_ploidy, probs = 0.90, na.rm = TRUE, names = FALSE, type = 7))
    row
  })
  sample_meta <- do.call(rbind, rows)
  rownames(sample_meta) <- sample_meta$sample_id
  sample_meta <- sample_meta[order(sample_meta$gemcitabine_dose_mg_per_kg, sample_meta$initial_ploidy, sample_meta$sample_id), , drop = FALSE]
  sample_meta$dose <- sample_meta$gemcitabine_dose
  sample_meta$dose_mg <- sample_meta$gemcitabine_dose_mg_per_kg
  sample_meta
}

ecdf_distance_test <- function(comparison_df, group_col, group_a, group_b, stratify_by_ploidy = FALSE) {
  comparison_df <- comparison_df[comparison_df[[group_col]] %in% c(group_a, group_b), , drop = FALSE]
  local_samples <- sort(unique(comparison_df$sample_id))
  if (length(local_samples) < 3) {
    return(list(statistic = NA_real_, p_value = NA_real_, n_permutations = NA_integer_, n_samples = length(local_samples)))
  }
  local_meta <- unique(comparison_df[, c("sample_id", "initial_ploidy", group_col), drop = FALSE])
  groups <- setNames(as.character(local_meta[[group_col]]), local_meta$sample_id)[local_samples]
  ploidy <- setNames(as.character(local_meta$initial_ploidy), local_meta$sample_id)[local_samples]
  if (!all(c(group_a, group_b) %in% groups)) {
    return(list(statistic = NA_real_, p_value = NA_real_, n_permutations = NA_integer_, n_samples = length(local_samples)))
  }

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
    choices <- combn(seq_along(local_samples), sum(groups == group_a), simplify = FALSE)
    perm_stats <- vapply(choices, function(group_a_idx) {
      perm_groups <- rep(group_b, length(local_samples))
      perm_groups[group_a_idx] <- group_a
      stat_fun(perm_groups)
    }, numeric(1))
  } else {
    choices_by_ploidy <- lapply(split(seq_along(local_samples), ploidy), function(idx) {
      k <- sum(groups[idx] == group_a)
      combn(idx, k, simplify = FALSE)
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

  list(
    statistic = observed,
    p_value = mean(perm_stats >= observed - 1e-15),
    n_permutations = length(perm_stats),
    n_samples = length(local_samples)
  )
}

run_dose_tests <- function(df) {
  out <- list()
  df$combined_dose_group <- ifelse(df$gemcitabine_dose_mg_per_kg == 0, "0mg/kg", "30+120mg/kg")
  combined <- ecdf_distance_test(df, "combined_dose_group", "0mg/kg", "30+120mg/kg", FALSE)
  combined_strat <- ecdf_distance_test(df, "combined_dose_group", "0mg/kg", "30+120mg/kg", TRUE)
  out[[length(out) + 1]] <- data.frame(
    comparison = "0_vs_30plus120",
    group_a = "0mg/kg",
    group_b = "30+120mg/kg",
    stratified_by_initial_ploidy = FALSE,
    statistic = combined$statistic,
    p_value = combined$p_value,
    n_samples = combined$n_samples,
    n_permutations = combined$n_permutations,
    stringsAsFactors = FALSE
  )
  out[[length(out) + 1]] <- data.frame(
    comparison = "0_vs_30plus120",
    group_a = "0mg/kg",
    group_b = "30+120mg/kg",
    stratified_by_initial_ploidy = TRUE,
    statistic = combined_strat$statistic,
    p_value = combined_strat$p_value,
    n_samples = combined_strat$n_samples,
    n_permutations = combined_strat$n_permutations,
    stringsAsFactors = FALSE
  )

  dose_values <- sort(unique(df$gemcitabine_dose_mg_per_kg))
  if (length(dose_values) >= 2) {
    dose_pairs <- combn(dose_values, 2)
    for (i in seq_len(ncol(dose_pairs))) {
      dose_a <- dose_pairs[1, i]
      dose_b <- dose_pairs[2, i]
      comparison_df <- df[df$gemcitabine_dose_mg_per_kg %in% c(dose_a, dose_b), , drop = FALSE]
      dose_a_label <- unique(comparison_df$gemcitabine_dose[comparison_df$gemcitabine_dose_mg_per_kg == dose_a])[1]
      dose_b_label <- unique(comparison_df$gemcitabine_dose[comparison_df$gemcitabine_dose_mg_per_kg == dose_b])[1]
      comparison_df$comparison_group <- comparison_df$gemcitabine_dose
      test <- ecdf_distance_test(comparison_df, "comparison_group", dose_a_label, dose_b_label, FALSE)
      out[[length(out) + 1]] <- data.frame(
        comparison = paste0(dose_a, "_vs_", dose_b),
        group_a = dose_a_label,
        group_b = dose_b_label,
        stratified_by_initial_ploidy = FALSE,
        statistic = test$statistic,
        p_value = test$p_value,
        n_samples = test$n_samples,
        n_permutations = test$n_permutations,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

calculate_shift_metrics <- function(df, sample_meta) {
  ecdf_grid <- sort(unique(df$pseudotime))
  sample_ids <- rownames(sample_meta)
  rows <- lapply(sample_ids, function(sample_id) {
    sample_ploidy <- sample_meta[sample_id, "initial_ploidy"]
    sample_dose <- sample_meta[sample_id, "dose_mg"]
    controls <- sample_meta$sample_id[
      sample_meta$initial_ploidy == sample_ploidy &
        sample_meta$dose_mg == 0
    ]
    reference_samples <- if (sample_dose == 0) setdiff(controls, sample_id) else controls
    sample_pseudotime <- df$pseudotime[df$sample_id == sample_id]
    reference_pseudotime <- df$pseudotime[df$sample_id %in% reference_samples]
    if (length(reference_pseudotime) == 0 || length(sample_pseudotime) == 0) {
      ecdf_rmse <- NA_real_
      ecdf_ks <- NA_real_
      ecdf_mean_abs <- NA_real_
      signed_mean_shift <- NA_real_
    } else {
      ecdf_delta <- stats::ecdf(sample_pseudotime)(ecdf_grid) -
        stats::ecdf(reference_pseudotime)(ecdf_grid)
      ecdf_rmse <- sqrt(mean(ecdf_delta^2))
      ecdf_ks <- max(abs(ecdf_delta))
      ecdf_mean_abs <- mean(abs(ecdf_delta))
      signed_mean_shift <- mean(sample_pseudotime) - mean(reference_pseudotime)
    }
    data.frame(
      sample_id = sample_id,
      initial_ploidy = sample_ploidy,
      dose = sample_meta[sample_id, "dose"],
      dose_mg = sample_dose,
      TGI_percent_auc = sample_meta[sample_id, "TGI_percent_auc"],
      mean_cell_ploidy = sample_meta[sample_id, "mean_cell_ploidy_exported"],
      ploidy_p90 = sample_meta[sample_id, "ploidy_p90"],
      n_cells = sample_meta[sample_id, "n_cells"],
      reference_n_samples = length(reference_samples),
      ecdf_rmse = ecdf_rmse,
      ecdf_ks = ecdf_ks,
      ecdf_mean_abs = ecdf_mean_abs,
      signed_mean_shift = signed_mean_shift,
      stringsAsFactors = FALSE
    )
  })
  shift_metrics <- do.call(rbind, rows)
  for (col in endpoint_tgi_cols(names(sample_meta))) {
    shift_metrics[[col]] <- sample_meta[shift_metrics$sample_id, col]
  }
  shift_metrics[order(shift_metrics$dose_mg, shift_metrics$initial_ploidy, shift_metrics$sample_id), , drop = FALSE]
}

run_tgi_associations <- function(shift_metrics) {
  treated <- shift_metrics[shift_metrics$dose_mg > 0, , drop = FALSE]
  tgi_cols <- c("TGI_percent_auc", endpoint_tgi_cols(names(shift_metrics)))
  shift_cols <- c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")
  rows <- list()
  for (tgi_col in tgi_cols) {
    for (shift_col in shift_cols) {
      for (method in c("pearson", "spearman")) {
        test <- safe_cor(treated[[shift_col]], treated[[tgi_col]], method = method)
        rows[[length(rows) + 1]] <- data.frame(
          sample_set = "treated",
          shift_metric = shift_col,
          tgi_measure = tgi_col,
          tgi_day = ifelse(grepl("^TGI_percent_Day_", tgi_col), endpoint_day_num(tgi_col), NA_real_),
          method = method,
          n = test$n,
          estimate = test$estimate,
          p_value = test$p_value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

run_ploidy_associations <- function(shift_metrics) {
  rows <- list()
  for (sample_set in c("all", "treated")) {
    local <- if (identical(sample_set, "treated")) shift_metrics[shift_metrics$dose_mg > 0, , drop = FALSE] else shift_metrics
    for (shift_col in c("ecdf_rmse", "signed_mean_shift")) {
      for (ploidy_col in c("mean_cell_ploidy", "ploidy_p90")) {
        for (method in c("pearson", "spearman")) {
          test <- safe_cor(local[[shift_col]], local[[ploidy_col]], method = method)
          rows[[length(rows) + 1]] <- data.frame(
            sample_set = sample_set,
            shift_metric = shift_col,
            ploidy_measure = ploidy_col,
            test = method,
            n = test$n,
            estimate = test$estimate,
            p_value = test$p_value,
            stringsAsFactors = FALSE
          )
        }
      }
      wt <- safe_wilcox(local[[shift_col]], local$initial_ploidy)
      rows[[length(rows) + 1]] <- data.frame(
        sample_set = sample_set,
        shift_metric = shift_col,
        ploidy_measure = "initial_ploidy",
        test = "wilcoxon",
        n = wt$n,
        estimate = wt$statistic,
        p_value = wt$p_value,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

make_growth_curve_long <- function(sample_meta) {
  raw_cols <- day_cols_by_prefix(names(sample_meta), "tumor_volume_")
  delta_cols <- day_cols_by_prefix(names(sample_meta), "tumor_volume_delta_")
  tgi_cols <- endpoint_tgi_cols(names(sample_meta))
  all_days <- sort(unique(c(
    suppressWarnings(as.numeric(sub("^tumor_volume_Day_", "", raw_cols))),
    suppressWarnings(as.numeric(sub("^tumor_volume_delta_Day_", "", delta_cols))),
    suppressWarnings(as.numeric(sub("^TGI_percent_Day_", "", tgi_cols)))
  )))
  rows <- lapply(seq_len(nrow(sample_meta)), function(i) {
    sample_row <- sample_meta[i, , drop = FALSE]
    do.call(rbind, lapply(all_days, function(day) {
      day_label <- paste0("Day_", day)
      raw_col <- paste0("tumor_volume_", day_label)
      delta_col <- paste0("tumor_volume_delta_", day_label)
      tgi_col <- paste0("TGI_percent_", day_label)
      control_delta_col <- paste0("matched_control_mean_delta_", day_label)
      data.frame(
        sample_id = sample_row$sample_id,
        initial_ploidy = sample_row$initial_ploidy,
        dose = sample_row$dose,
        dose_mg = sample_row$dose_mg,
        day = day,
        day_label = day_label,
        tumor_volume = if (raw_col %in% names(sample_row)) sample_row[[raw_col]] else NA_real_,
        tumor_volume_delta = if (delta_col %in% names(sample_row)) sample_row[[delta_col]] else NA_real_,
        matched_control_mean_delta = if (control_delta_col %in% names(sample_row)) sample_row[[control_delta_col]] else NA_real_,
        TGI_percent_endpoint = if (tgi_col %in% names(sample_row)) sample_row[[tgi_col]] else NA_real_,
        TGI_percent_auc = sample_row$TGI_percent_auc,
        tumor_volume_auc_delta = if ("tumor_volume_auc_delta" %in% names(sample_row)) sample_row$tumor_volume_auc_delta else NA_real_,
        matched_control_mean_auc_delta = if ("matched_control_mean_auc_delta" %in% names(sample_row)) sample_row$matched_control_mean_auc_delta else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  out$dose <- factor(out$dose, levels = unique(sample_meta$dose[order(sample_meta$dose_mg)]))
  out
}

make_auc_calculation_table <- function(sample_meta, growth_long) {
  control_curve <- aggregate(
    cbind(tumor_volume_delta, tumor_volume) ~ initial_ploidy + day + day_label,
    growth_long[growth_long$dose_mg == 0, , drop = FALSE],
    function(x) mean(x, na.rm = TRUE)
  )
  names(control_curve)[names(control_curve) == "tumor_volume_delta"] <- "control_mean_delta_curve"
  names(control_curve)[names(control_curve) == "tumor_volume"] <- "control_mean_volume_curve"

  auc_table <- sample_meta[, c(
    "sample_id", "initial_ploidy", "dose", "dose_mg",
    "tumor_volume_baseline_day", "tumor_volume_baseline",
    "tumor_volume_auc", "tumor_volume_auc_delta",
    "auc_start_day", "auc_end_day", "auc_n_days",
    "matched_control_mean_auc_delta", "matched_control_n_auc",
    "TGI_percent_auc"
  ), drop = FALSE]
  auc_table$auc_formula <- "100 * (1 - tumor_volume_auc_delta / matched_control_mean_auc_delta)"
  list(auc_table = auc_table, control_curve = control_curve)
}

run_endpoint_vs_auc_stats <- function(shift_metrics) {
  treated <- shift_metrics[shift_metrics$dose_mg > 0, , drop = FALSE]
  tgi_cols <- endpoint_tgi_cols(names(treated))
  rows <- lapply(tgi_cols, function(col) {
    day <- endpoint_day_num(col)
    pearson_auc <- safe_cor(treated[[col]], treated$TGI_percent_auc, "pearson")
    spearman_auc <- safe_cor(treated[[col]], treated$TGI_percent_auc, "spearman")
    pearson_shift <- safe_cor(treated[[col]], treated$ecdf_rmse, "pearson")
    data.frame(
      tgi_measure = col,
      tgi_day = day,
      endpoint_vs_auc_pearson_r = pearson_auc$estimate,
      endpoint_vs_auc_pearson_p = pearson_auc$p_value,
      endpoint_vs_auc_spearman_rho = spearman_auc$estimate,
      endpoint_vs_auc_spearman_p = spearman_auc$p_value,
      endpoint_vs_ecdf_rmse_pearson_r = pearson_shift$estimate,
      endpoint_vs_ecdf_rmse_pearson_p = pearson_shift$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_original_reported_stats <- function(dose_tests, tgi_stats, ploidy_stats) {
  get_p <- function(comparison, stratified = FALSE) {
    hit <- dose_tests$comparison == comparison & dose_tests$stratified_by_initial_ploidy == stratified
    if (any(hit)) dose_tests$p_value[which(hit)[1]] else NA_real_
  }
  get_tgi <- function(method, value_col) {
    hit <- tgi_stats$sample_set == "treated" &
      tgi_stats$shift_metric == "ecdf_rmse" &
      tgi_stats$tgi_measure == "TGI_percent_auc" &
      tgi_stats$method == method
    if (!any(hit)) return(NA_real_)
    tgi_stats[which(hit)[1], value_col]
  }
  get_ploidy_p <- function(sample_set, ploidy_measure, test) {
    hit <- ploidy_stats$sample_set == sample_set &
      ploidy_stats$shift_metric == "ecdf_rmse" &
      ploidy_stats$ploidy_measure == ploidy_measure &
      ploidy_stats$test == test
    if (!any(hit)) return(NA_real_)
    ploidy_stats$p_value[which(hit)[1]]
  }
  data.frame(
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
      get_p("0_vs_30plus120", FALSE),
      get_p("0_vs_30plus120", TRUE),
      get_p("0_vs_30", FALSE),
      get_p("0_vs_120", FALSE),
      get_p("30_vs_120", FALSE),
      get_tgi("pearson", "estimate"),
      get_tgi("pearson", "p_value"),
      get_tgi("spearman", "estimate"),
      get_tgi("spearman", "p_value"),
      get_ploidy_p("all", "mean_cell_ploidy", "pearson"),
      get_ploidy_p("treated", "mean_cell_ploidy", "pearson"),
      get_ploidy_p("all", "initial_ploidy", "wilcoxon"),
      get_ploidy_p("treated", "initial_ploidy", "wilcoxon")
    ),
    stringsAsFactors = FALSE
  )
}

plot_compartment_results <- function(df, sample_meta, shift_metrics, tgi_stats, cfg) {
  fig_dir <- ensure_dir(file.path(cfg$output_dir, "figures"))
  original_dir <- ensure_dir(file.path(fig_dir, "original_analysis"))
  auc_dir <- ensure_dir(file.path(fig_dir, "auc_growth_curves"))
  association_dir <- ensure_dir(file.path(fig_dir, "tgi_associations"))
  dose_levels <- unique(sample_meta$dose[order(sample_meta$dose_mg)])
  dose_cols <- dose_cols_all[names(dose_cols_all) %in% dose_levels]
  treated_cols <- dose_cols[names(dose_cols) != "0mg/kg"]

  df$dose <- factor(df$gemcitabine_dose, levels = dose_levels)
  df$initial_ploidy <- factor(df$initial_ploidy, levels = unique(sample_meta$initial_ploidy))
  sample_meta$dose <- factor(sample_meta$dose, levels = dose_levels)
  shift_metrics$dose <- factor(shift_metrics$dose, levels = dose_levels)
  treated_shift <- shift_metrics[shift_metrics$dose_mg > 0, , drop = FALSE]
  treated_shift$dose <- factor(treated_shift$dose, levels = dose_levels[dose_levels != "0mg/kg"])

  density_plot <- ggplot(df, aes(x = pseudotime, color = dose, fill = dose)) +
    geom_density(alpha = 0.12, linewidth = 0.8) +
    facet_wrap(~ initial_ploidy, ncol = 1) +
    scale_color_manual(values = dose_cols, name = "dose") +
    scale_fill_manual(values = dose_cols, name = "dose") +
    labs(title = paste0(cfg$label, ": pseudotime distributions"), x = "pseudotime", y = "density") +
    plot_theme
  save_plot(density_plot, fig_dir, paste0(cfg$prefix, "_density_by_dose_and_ploidy"), 8, 5.5)
  save_plot(density_plot, original_dir, paste0(cfg$prefix, "_dose_group_comparisons_density_by_dose_and_ploidy"), 8, 5.5)

  ecdf_grid <- sort(unique(df$pseudotime))
  sample_ids <- rownames(sample_meta)
  ecdf_plot_df <- do.call(rbind, lapply(sample_ids, function(sample_id) {
    data.frame(
      sample_id = sample_id,
      initial_ploidy = sample_meta[sample_id, "initial_ploidy"],
      dose = sample_meta[sample_id, "dose"],
      pseudotime = ecdf_grid,
      ecdf = stats::ecdf(df$pseudotime[df$sample_id == sample_id])(ecdf_grid),
      stringsAsFactors = FALSE
    )
  }))
  ecdf_plot_df$dose <- factor(ecdf_plot_df$dose, levels = dose_levels)
  ecdf_plot <- ggplot(ecdf_plot_df, aes(x = pseudotime, y = ecdf, color = dose, group = sample_id)) +
    geom_line(alpha = 0.55, linewidth = 0.55) +
    facet_wrap(~ initial_ploidy, ncol = 1) +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(title = paste0(cfg$label, ": sample ECDFs"), x = "pseudotime", y = "ECDF") +
    plot_theme
  save_plot(ecdf_plot, fig_dir, paste0(cfg$prefix, "_sample_ecdf_by_dose_and_ploidy"), 8, 5.5)
  save_plot(ecdf_plot, original_dir, paste0(cfg$prefix, "_dose_group_comparisons_sample_ecdf_by_dose_and_ploidy"), 8, 5.5)

  mean_plot <- ggplot(sample_meta, aes(x = dose, y = mean_pseudotime, color = dose)) +
    geom_boxplot(aes(group = dose), outlier.shape = NA, alpha = 0, color = "grey55", linewidth = 0.45) +
    geom_point(aes(shape = initial_ploidy), size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(title = paste0(cfg$label, ": sample mean pseudotime"), x = "dose", y = "sample mean pseudotime") +
    plot_theme
  save_plot(mean_plot, fig_dir, paste0(cfg$prefix, "_sample_mean_pseudotime_by_dose"), 6.5, 4.8)
  save_plot(mean_plot, original_dir, paste0(cfg$prefix, "_dose_group_comparisons_sample_mean_pseudotime_by_dose"), 6.5, 4.8)

  tgi_auc_plot <- ggplot(treated_shift, aes(x = ecdf_rmse, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    geom_point(aes(shape = initial_ploidy), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.5, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = treated_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": TGI AUC versus pseudotime distance"),
      x = "ECDF RMSE from ploidy-matched 0 mg/kg reference",
      y = "TGI AUC (%)"
    ) +
    plot_theme
  save_plot(tgi_auc_plot, fig_dir, paste0(cfg$prefix, "_TGI_AUC_vs_ecdf_rmse"), 8, 5.5)
  save_plot(tgi_auc_plot, original_dir, paste0(cfg$prefix, "_pseudotime_shift_associations_TGI_AUC_ecdf_distance_vs_TGI_AUC"), 8, 5.5)
  save_plot(tgi_auc_plot, association_dir, paste0(cfg$prefix, "_TGI_AUC_vs_ecdf_rmse"), 8, 5.5)

  signed_shift_plot <- ggplot(treated_shift, aes(x = signed_mean_shift, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    geom_point(aes(shape = initial_ploidy), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.5, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = treated_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": TGI AUC versus signed pseudotime shift"),
      x = "sample mean pseudotime minus ploidy-matched 0 mg/kg reference",
      y = "TGI AUC (%)"
    ) +
    plot_theme
  save_plot(signed_shift_plot, fig_dir, paste0(cfg$prefix, "_TGI_AUC_vs_signed_mean_shift"), 8, 5.5)
  save_plot(signed_shift_plot, association_dir, paste0(cfg$prefix, "_TGI_AUC_vs_signed_mean_shift"), 8, 5.5)

  endpoint_stats <- tgi_stats[
    tgi_stats$shift_metric == "ecdf_rmse" &
      tgi_stats$method == "pearson" &
      grepl("^TGI_percent_Day_", tgi_stats$tgi_measure),
    ,
    drop = FALSE
  ]
  if (nrow(endpoint_stats) > 0) {
    endpoint_stats <- endpoint_stats[order(endpoint_stats$tgi_day), , drop = FALSE]
    endpoint_stats$neg_log10_p <- -log10(pmax(endpoint_stats$p_value, .Machine$double.xmin))
    endpoint_plot <- ggplot(endpoint_stats, aes(x = tgi_day, y = estimate)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
      geom_line(color = "black", linewidth = 0.6) +
      geom_point(aes(size = neg_log10_p), color = "#2166ac") +
      scale_size_continuous(name = "-log10(p)", range = c(1.8, 4.5)) +
      labs(
        title = paste0(cfg$label, ": endpoint TGI correlation over tumor-growth days"),
        x = "tumor-growth day",
        y = "Pearson r: ECDF RMSE vs endpoint TGI"
      ) +
      plot_theme
    save_plot(endpoint_plot, fig_dir, paste0(cfg$prefix, "_endpoint_TGI_correlations_over_days"), 7, 4.6)
    save_plot(endpoint_plot, association_dir, paste0(cfg$prefix, "_endpoint_TGI_correlations_over_days"), 7, 4.6)
  }

  ploidy_auc_plot <- ggplot(treated_shift, aes(x = ploidy_p90, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
    geom_point(aes(shape = initial_ploidy), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.5, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = treated_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": TGI AUC versus high-end cell ploidy"),
      x = "90th pct cell ploidy",
      y = "TGI AUC (%)"
    ) +
    plot_theme
  save_plot(ploidy_auc_plot, fig_dir, paste0(cfg$prefix, "_TGI_AUC_vs_ploidy_p90"), 7, 5)
  save_plot(ploidy_auc_plot, association_dir, paste0(cfg$prefix, "_TGI_AUC_vs_ploidy_p90"), 7, 5)

  if ("TGI_percent_Day_24" %in% names(treated_shift)) {
    ploidy_day24_plot <- ggplot(treated_shift, aes(x = ploidy_p90, y = TGI_percent_Day_24, color = dose)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
      geom_point(aes(shape = initial_ploidy), size = 2.8) +
      geom_text(aes(label = sample_id), size = 2.5, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
      scale_color_manual(values = treated_cols, name = "dose") +
      labs(
        title = paste0(cfg$label, ": Day 24 TGI versus high-end cell ploidy"),
        x = "90th pct cell ploidy",
        y = "TGI percent: Day 24"
      ) +
      plot_theme
    save_plot(ploidy_day24_plot, fig_dir, paste0(cfg$prefix, "_TGI_Day24_vs_ploidy_p90"), 7, 5)
    save_plot(ploidy_day24_plot, original_dir, paste0(cfg$prefix, "_ploidy_vs_TGI_Day24_p90_treated_by_dose"), 7, 5)
    save_plot(ploidy_day24_plot, association_dir, paste0(cfg$prefix, "_TGI_Day24_vs_ploidy_p90"), 7, 5)
  }

  growth_long <- make_growth_curve_long(sample_meta)
  growth_long$dose <- factor(growth_long$dose, levels = dose_levels)
  auc_tables <- make_auc_calculation_table(sample_meta, growth_long)
  control_curve <- auc_tables$control_curve

  raw_growth_plot <- ggplot(growth_long, aes(x = day, y = tumor_volume, color = dose, group = sample_id)) +
    geom_line(alpha = 0.62, linewidth = 0.55) +
    geom_point(size = 1.35, alpha = 0.8) +
    facet_wrap(~ initial_ploidy, ncol = 1, scales = "free_y") +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": tumor growth curves"),
      x = "tumor-growth day",
      y = "tumor volume"
    ) +
    plot_theme
  save_plot(raw_growth_plot, auc_dir, paste0(cfg$prefix, "_tumor_growth_raw_curves_by_dose_and_ploidy"), 8, 5.5)

  baseline_adjusted_plot <- ggplot(growth_long, aes(x = day, y = tumor_volume_delta, color = dose, group = sample_id)) +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.35) +
    geom_line(alpha = 0.62, linewidth = 0.55) +
    geom_point(size = 1.35, alpha = 0.8) +
    facet_wrap(~ initial_ploidy, ncol = 1, scales = "free_y") +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": baseline-adjusted tumor growth curves for AUC TGI"),
      x = "tumor-growth day",
      y = "tumor volume minus Day 0"
    ) +
    plot_theme
  save_plot(baseline_adjusted_plot, auc_dir, paste0(cfg$prefix, "_tumor_growth_baseline_adjusted_auc_curves"), 8, 5.5)

  mean_auc_curve <- aggregate(
    tumor_volume_delta ~ initial_ploidy + dose + dose_mg + day,
    growth_long,
    function(x) mean(x, na.rm = TRUE)
  )
  mean_auc_curve$dose <- factor(mean_auc_curve$dose, levels = dose_levels)
  mean_auc_area_plot <- ggplot(mean_auc_curve, aes(x = day, y = tumor_volume_delta, color = dose, fill = dose, group = dose)) +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.35) +
    geom_area(alpha = 0.14, position = "identity", color = NA) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~ initial_ploidy, ncol = 1, scales = "free_y") +
    scale_color_manual(values = dose_cols, name = "dose") +
    scale_fill_manual(values = dose_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": mean baseline-adjusted AUC curves by dose"),
      x = "tumor-growth day",
      y = "mean tumor volume minus Day 0"
    ) +
    plot_theme
  save_plot(mean_auc_area_plot, auc_dir, paste0(cfg$prefix, "_mean_baseline_adjusted_AUC_area_curves_by_dose"), 8, 5.5)

  control_curve$dose <- "matched 0mg/kg mean"
  control_curve$dose_mg <- 0
  control_curve$sample_id <- paste0(control_curve$initial_ploidy, "_matched_control_mean")
  control_curve$curve_type <- "matched 0mg/kg mean"
  treated_curve <- growth_long[growth_long$dose_mg > 0, , drop = FALSE]
  treated_curve$curve_type <- "treated sample"
  treated_curve$control_mean_delta_curve <- NA_real_
  auc_overlay <- rbind(
    data.frame(
      sample_id = treated_curve$sample_id,
      initial_ploidy = treated_curve$initial_ploidy,
      dose = as.character(treated_curve$dose),
      dose_mg = treated_curve$dose_mg,
      day = treated_curve$day,
      day_label = treated_curve$day_label,
      curve_type = treated_curve$curve_type,
      baseline_adjusted_volume = treated_curve$tumor_volume_delta,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sample_id = control_curve$sample_id,
      initial_ploidy = control_curve$initial_ploidy,
      dose = control_curve$dose,
      dose_mg = control_curve$dose_mg,
      day = control_curve$day,
      day_label = control_curve$day_label,
      curve_type = control_curve$curve_type,
      baseline_adjusted_volume = control_curve$control_mean_delta_curve,
      stringsAsFactors = FALSE
    )
  )
  auc_overlay$curve_type <- factor(auc_overlay$curve_type, levels = c("matched 0mg/kg mean", "treated sample"))
  auc_overlay$dose <- factor(auc_overlay$dose, levels = c("matched 0mg/kg mean", setdiff(dose_levels, "0mg/kg")))
  overlay_cols <- c("matched 0mg/kg mean" = "black", treated_cols)
  auc_overlay_plot <- ggplot(auc_overlay, aes(x = day, y = baseline_adjusted_volume, color = dose, group = interaction(sample_id, curve_type), linetype = curve_type)) +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.35) +
    geom_line(alpha = 0.75, linewidth = 0.65) +
    facet_wrap(~ initial_ploidy, ncol = 1, scales = "free_y") +
    scale_color_manual(values = overlay_cols, name = "curve") +
    scale_linetype_manual(values = c("matched 0mg/kg mean" = "dashed", "treated sample" = "solid"), name = "curve type") +
    labs(
      title = paste0(cfg$label, ": treated AUC curves versus ploidy-matched controls"),
      x = "tumor-growth day",
      y = "tumor volume minus Day 0"
    ) +
    plot_theme
  save_plot(auc_overlay_plot, auc_dir, paste0(cfg$prefix, "_treated_sample_auc_vs_matched_control_auc"), 8, 5.5)

  tgi_auc_boxplot <- ggplot(sample_meta, aes(x = dose, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_boxplot(aes(group = dose), outlier.shape = NA, alpha = 0, color = "grey55", linewidth = 0.45) +
    geom_point(aes(shape = initial_ploidy), size = 2.8, position = position_jitter(width = 0.08, height = 0)) +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": AUC-based TGI by dose"),
      x = "dose",
      y = "TGI AUC (%)"
    ) +
    plot_theme
  save_plot(tgi_auc_boxplot, auc_dir, paste0(cfg$prefix, "_TGI_AUC_by_dose_and_ploidy"), 6.5, 4.8)

  endpoint_curve <- growth_long[is.finite(growth_long$TGI_percent_endpoint), , drop = FALSE]
  endpoint_curve <- endpoint_curve[endpoint_curve$dose_mg > 0, , drop = FALSE]
  endpoint_tgi_plot <- ggplot(endpoint_curve, aes(x = day, y = TGI_percent_endpoint, color = dose, group = sample_id)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_line(alpha = 0.7, linewidth = 0.6) +
    geom_point(size = 1.4, alpha = 0.85) +
    facet_wrap(~ initial_ploidy, ncol = 1) +
    scale_color_manual(values = treated_cols, name = "dose") +
    labs(
      title = paste0(cfg$label, ": endpoint TGI over tumor-growth days"),
      x = "tumor-growth day",
      y = "endpoint TGI (%)"
    ) +
    plot_theme
  save_plot(endpoint_tgi_plot, auc_dir, paste0(cfg$prefix, "_endpoint_TGI_curves_over_days"), 8, 5.5)

  endpoint_vs_auc_stats <- run_endpoint_vs_auc_stats(shift_metrics)
  if (nrow(endpoint_vs_auc_stats) > 0) {
    endpoint_vs_auc_plot <- ggplot(endpoint_vs_auc_stats, aes(x = tgi_day, y = endpoint_vs_auc_pearson_r)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
      geom_line(color = "black", linewidth = 0.6) +
      geom_point(aes(size = -log10(pmax(endpoint_vs_auc_pearson_p, .Machine$double.xmin))), color = "#b2182b") +
      scale_size_continuous(name = "-log10(p)", range = c(1.8, 4.5)) +
      labs(
        title = paste0(cfg$label, ": endpoint TGI agreement with AUC TGI"),
        x = "tumor-growth day",
        y = "Pearson r: endpoint TGI vs AUC TGI"
      ) +
      plot_theme
    save_plot(endpoint_vs_auc_plot, auc_dir, paste0(cfg$prefix, "_endpoint_vs_auc_TGI_summary"), 7, 4.6)
  }

  invisible(list(
    growth_long = growth_long,
    auc_calculation = auc_tables$auc_table,
    control_curve = control_curve,
    endpoint_vs_auc_stats = endpoint_vs_auc_stats
  ))
}

write_compartment_summary <- function(cfg, df, dose_tests, tgi_stats, ploidy_stats) {
  get_p <- function(comparison, stratified = FALSE) {
    hit <- dose_tests$comparison == comparison & dose_tests$stratified_by_initial_ploidy == stratified
    if (any(hit)) dose_tests$p_value[which(hit)[1]] else NA_real_
  }
  auc_row <- tgi_stats[
    tgi_stats$sample_set == "treated" &
      tgi_stats$shift_metric == "ecdf_rmse" &
      tgi_stats$tgi_measure == "TGI_percent_auc" &
      tgi_stats$method == "pearson",
    ,
    drop = FALSE
  ]
  auc_r <- if (nrow(auc_row) > 0) auc_row$estimate[1] else NA_real_
  auc_p <- if (nrow(auc_row) > 0) auc_row$p_value[1] else NA_real_
  lines <- c(
    paste0("# 04h Summary: ", cfg$label),
    "",
    paste0("Input cells: ", nrow(df), ". Samples: ", length(unique(df$sample_id)), "."),
    "",
    "## Pseudotime dose comparisons",
    "",
    paste0("- 0 mg/kg vs 30+120 mg/kg ECDF permutation p = ", format_p(get_p("0_vs_30plus120", FALSE)), "."),
    paste0("- Ploidy-stratified 0 mg/kg vs 30+120 mg/kg ECDF permutation p = ", format_p(get_p("0_vs_30plus120", TRUE)), "."),
    paste0("- 0 vs 30 mg/kg p = ", format_p(get_p("0_vs_30", FALSE)), "."),
    paste0("- 0 vs 120 mg/kg p = ", format_p(get_p("0_vs_120", FALSE)), "."),
    paste0("- 30 vs 120 mg/kg p = ", format_p(get_p("30_vs_120", FALSE)), "."),
    "",
    "## TGI association",
    "",
    paste0("- Treated-only ECDF RMSE vs TGI AUC Pearson r = ", format_num(auc_r), ", p = ", format_p(auc_p), "."),
    "",
    "## Output files",
    "",
    "- `data/sample_metrics.csv`",
    "- `data/shift_metrics.csv`",
    "- `stats/dose_group_tests.csv`",
    "- `stats/tgi_association_stats.csv`",
    "- `stats/ploidy_association_stats.csv`",
    "- `stats/endpoint_vs_auc_tgi_stats.csv`",
    "- `data/auc_growth_curve_long.csv`",
    "- `data/auc_calculation_by_sample.csv`",
    "- `figures/original_analysis/*.png` and `*.pdf`",
    "- `figures/auc_growth_curves/*.png` and `*.pdf`",
    "- `figures/tgi_associations/*.png` and `*.pdf`"
  )
  write_text(lines, file.path(cfg$output_dir, "analysis_summary.md"))
}

run_compartment <- function(compartment_name, cfg) {
  ensure_dir(cfg$output_dir)
  df <- read_compartment_data(compartment_name, cfg)
  sample_meta <- build_sample_meta(df)
  dose_tests <- run_dose_tests(df)
  shift_metrics <- calculate_shift_metrics(df, sample_meta)
  tgi_stats <- run_tgi_associations(shift_metrics)
  ploidy_stats <- run_ploidy_associations(shift_metrics)
  original_reported_stats <- make_original_reported_stats(dose_tests, tgi_stats, ploidy_stats)

  write_csv(sample_meta, file.path(cfg$output_dir, "data", "sample_metrics.csv"))
  write_csv(shift_metrics, file.path(cfg$output_dir, "data", "shift_metrics.csv"))
  write_csv(dose_tests, file.path(cfg$output_dir, "stats", "dose_group_tests.csv"))
  write_csv(tgi_stats, file.path(cfg$output_dir, "stats", "tgi_association_stats.csv"))
  write_csv(ploidy_stats, file.path(cfg$output_dir, "stats", "ploidy_association_stats.csv"))
  write_csv(original_reported_stats, file.path(cfg$output_dir, "stats", paste0(cfg$prefix, "_pseudotimeAssociations_reported_stats.csv")))
  plot_outputs <- plot_compartment_results(df, sample_meta, shift_metrics, tgi_stats, cfg)
  write_csv(plot_outputs$growth_long, file.path(cfg$output_dir, "data", "auc_growth_curve_long.csv"))
  write_csv(plot_outputs$auc_calculation, file.path(cfg$output_dir, "data", "auc_calculation_by_sample.csv"))
  write_csv(plot_outputs$control_curve, file.path(cfg$output_dir, "data", "ploidy_matched_control_mean_growth_curve.csv"))
  write_csv(plot_outputs$endpoint_vs_auc_stats, file.path(cfg$output_dir, "stats", "endpoint_vs_auc_tgi_stats.csv"))
  write_compartment_summary(cfg, df, dose_tests, tgi_stats, ploidy_stats)

  list(
    df = df,
    sample_meta = sample_meta,
    dose_tests = dose_tests,
    shift_metrics = shift_metrics,
    tgi_stats = tgi_stats,
    ploidy_stats = ploidy_stats,
    original_reported_stats = original_reported_stats,
    growth_long = plot_outputs$growth_long,
    endpoint_vs_auc_stats = plot_outputs$endpoint_vs_auc_stats
  )
}

paired_comparison <- function(results) {
  cfg <- list(output_dir = file.path(output_root, "Comparison"))
  ensure_dir(cfg$output_dir)
  fig_dir <- ensure_dir(file.path(cfg$output_dir, "figures"))
  data_dir_out <- ensure_dir(file.path(cfg$output_dir, "data"))
  stats_dir_out <- ensure_dir(file.path(cfg$output_dir, "stats"))

  cc <- results$CellCycle$shift_metrics
  nc <- results$NonCellCycle$shift_metrics
  keep_cols <- c(
    "sample_id", "initial_ploidy", "dose", "dose_mg", "TGI_percent_auc",
    "mean_cell_ploidy", "ploidy_p90", "n_cells", "ecdf_rmse", "ecdf_ks",
    "ecdf_mean_abs", "signed_mean_shift"
  )
  comp <- merge(cc[, keep_cols, drop = FALSE], nc[, keep_cols, drop = FALSE], by = "sample_id", suffixes = c("_CellCycle", "_NonCellCycle"))
  comp$dose <- comp$dose_CellCycle
  comp$dose_mg <- comp$dose_mg_CellCycle
  comp$initial_ploidy <- comp$initial_ploidy_CellCycle
  comp$TGI_percent_auc <- comp$TGI_percent_auc_CellCycle
  comp$ecdf_rmse_delta_CellCycle_minus_NonCellCycle <- comp$ecdf_rmse_CellCycle - comp$ecdf_rmse_NonCellCycle
  comp$signed_mean_shift_delta_CellCycle_minus_NonCellCycle <- comp$signed_mean_shift_CellCycle - comp$signed_mean_shift_NonCellCycle
  comp <- comp[order(comp$dose_mg, comp$initial_ploidy, comp$sample_id), , drop = FALSE]
  write_csv(comp, file.path(data_dir_out, "cellcycle_vs_noncellcycle_sample_shift_comparison.csv"))

  rows <- list()
  for (sample_set in c("all", "treated")) {
    local <- if (identical(sample_set, "treated")) comp[comp$dose_mg > 0, , drop = FALSE] else comp
    for (metric in c("ecdf_rmse", "signed_mean_shift")) {
      x <- local[[paste0(metric, "_CellCycle")]]
      y <- local[[paste0(metric, "_NonCellCycle")]]
      for (method in c("paired_t", "paired_wilcox")) {
        test <- safe_paired_test(x, y, method)
        rows[[length(rows) + 1]] <- data.frame(
          sample_set = sample_set,
          comparison = paste0(metric, "_CellCycle_vs_NonCellCycle"),
          test = method,
          n = test$n,
          estimate = test$statistic,
          p_value = test$p_value,
          stringsAsFactors = FALSE
        )
      }
      for (method in c("pearson", "spearman")) {
        test <- safe_cor(x, y, method)
        rows[[length(rows) + 1]] <- data.frame(
          sample_set = sample_set,
          comparison = paste0(metric, "_CellCycle_cor_NonCellCycle"),
          test = method,
          n = test$n,
          estimate = test$estimate,
          p_value = test$p_value,
          stringsAsFactors = FALSE
        )
      }
    }
    for (metric in c("ecdf_rmse_delta_CellCycle_minus_NonCellCycle", "signed_mean_shift_delta_CellCycle_minus_NonCellCycle")) {
      for (method in c("pearson", "spearman")) {
        test <- safe_cor(local[[metric]], local$TGI_percent_auc, method)
        rows[[length(rows) + 1]] <- data.frame(
          sample_set = sample_set,
          comparison = paste0(metric, "_vs_TGI_AUC"),
          test = method,
          n = test$n,
          estimate = test$estimate,
          p_value = test$p_value,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  comparison_stats <- do.call(rbind, rows)
  write_csv(comparison_stats, file.path(stats_dir_out, "comparison_stats.csv"))

  long_shift <- rbind(
    data.frame(
      sample_id = comp$sample_id,
      compartment = "CellCycle",
      dose = comp$dose,
      dose_mg = comp$dose_mg,
      initial_ploidy = comp$initial_ploidy,
      TGI_percent_auc = comp$TGI_percent_auc,
      ecdf_rmse = comp$ecdf_rmse_CellCycle,
      signed_mean_shift = comp$signed_mean_shift_CellCycle,
      stringsAsFactors = FALSE
    ),
    data.frame(
      sample_id = comp$sample_id,
      compartment = "NonCellCycle",
      dose = comp$dose,
      dose_mg = comp$dose_mg,
      initial_ploidy = comp$initial_ploidy,
      TGI_percent_auc = comp$TGI_percent_auc,
      ecdf_rmse = comp$ecdf_rmse_NonCellCycle,
      signed_mean_shift = comp$signed_mean_shift_NonCellCycle,
      stringsAsFactors = FALSE
    )
  )
  long_shift$compartment <- factor(long_shift$compartment, levels = c("CellCycle", "NonCellCycle"))
  long_shift$dose <- factor(long_shift$dose, levels = unique(comp$dose[order(comp$dose_mg)]))
  write_csv(long_shift, file.path(data_dir_out, "cellcycle_noncellcycle_shift_long.csv"))

  dose_cols <- dose_cols_all[names(dose_cols_all) %in% levels(long_shift$dose)]
  paired_plot <- ggplot(long_shift, aes(x = compartment, y = ecdf_rmse, group = sample_id, color = dose)) +
    geom_line(alpha = 0.45, linewidth = 0.45) +
    geom_point(aes(shape = initial_ploidy), size = 2.7) +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(
      title = "Pseudotime shift magnitude by compartment",
      x = NULL,
      y = "ECDF RMSE from ploidy-matched 0 mg/kg reference"
    ) +
    plot_theme
  save_plot(paired_plot, fig_dir, "CellCycle_vs_NonCellCycle_paired_ecdf_rmse", 6.8, 5)

  scatter_plot <- ggplot(comp, aes(x = ecdf_rmse_CellCycle, y = ecdf_rmse_NonCellCycle, color = dose)) +
    geom_abline(slope = 1, intercept = 0, color = "grey70", linewidth = 0.4) +
    geom_point(aes(shape = initial_ploidy), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.5, nudge_y = 0.006, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = dose_cols, name = "dose") +
    labs(
      title = "Cell-cycle versus non-cell-cycle pseudotime shift",
      x = "CellCycle ECDF RMSE",
      y = "NonCellCycle ECDF RMSE"
    ) +
    plot_theme
  save_plot(scatter_plot, fig_dir, "CellCycle_vs_NonCellCycle_ecdf_rmse_scatter", 6.2, 5.4)

  treated_long <- long_shift[long_shift$dose_mg > 0, , drop = FALSE]
  treated_cols <- dose_cols[names(dose_cols) != "0mg/kg"]
  tgi_shift_plot <- ggplot(treated_long, aes(x = ecdf_rmse, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    geom_point(aes(shape = initial_ploidy), size = 2.6) +
    facet_wrap(~ compartment, nrow = 1) +
    scale_color_manual(values = treated_cols, name = "dose") +
    labs(
      title = "TGI AUC versus pseudotime shift by compartment",
      x = "ECDF RMSE from ploidy-matched 0 mg/kg reference",
      y = "TGI AUC (%)"
    ) +
    plot_theme
  save_plot(tgi_shift_plot, fig_dir, "CellCycle_vs_NonCellCycle_TGI_AUC_vs_ecdf_rmse", 8.8, 4.7)

  endpoint_stats <- rbind(
    data.frame(compartment = "CellCycle", results$CellCycle$tgi_stats, stringsAsFactors = FALSE),
    data.frame(compartment = "NonCellCycle", results$NonCellCycle$tgi_stats, stringsAsFactors = FALSE)
  )
  endpoint_stats <- endpoint_stats[
    endpoint_stats$shift_metric == "ecdf_rmse" &
      endpoint_stats$method == "pearson" &
      grepl("^TGI_percent_Day_", endpoint_stats$tgi_measure),
    ,
    drop = FALSE
  ]
  endpoint_stats <- endpoint_stats[order(endpoint_stats$compartment, endpoint_stats$tgi_day), , drop = FALSE]
  write_csv(endpoint_stats, file.path(stats_dir_out, "endpoint_TGI_correlation_comparison.csv"))
  if (nrow(endpoint_stats) > 0) {
    endpoint_plot <- ggplot(endpoint_stats, aes(x = tgi_day, y = estimate, color = compartment)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
      geom_line(linewidth = 0.65) +
      geom_point(size = 2.3) +
      scale_color_manual(values = c("CellCycle" = "#2166ac", "NonCellCycle" = "#1b9e77"), name = "compartment") +
      labs(
        title = "Endpoint TGI correlations across tumor-growth days",
        x = "tumor-growth day",
        y = "Pearson r: ECDF RMSE vs endpoint TGI"
      ) +
      plot_theme
    save_plot(endpoint_plot, fig_dir, "CellCycle_vs_NonCellCycle_endpoint_TGI_correlations_over_days", 7.3, 4.8)
  }

  cc_auc <- results$CellCycle$tgi_stats[
    results$CellCycle$tgi_stats$shift_metric == "ecdf_rmse" &
      results$CellCycle$tgi_stats$tgi_measure == "TGI_percent_auc" &
      results$CellCycle$tgi_stats$method == "pearson",
    ,
    drop = FALSE
  ]
  nc_auc <- results$NonCellCycle$tgi_stats[
    results$NonCellCycle$tgi_stats$shift_metric == "ecdf_rmse" &
      results$NonCellCycle$tgi_stats$tgi_measure == "TGI_percent_auc" &
      results$NonCellCycle$tgi_stats$method == "pearson",
    ,
    drop = FALSE
  ]
  lines <- c(
    "# 04h CellCycle vs NonCellCycle Summary",
    "",
    paste0("Paired samples: ", nrow(comp), ". Treated paired samples: ", sum(comp$dose_mg > 0), "."),
    "",
    "## TGI AUC association",
    "",
    paste0("- CellCycle ECDF RMSE vs TGI AUC: Pearson r = ", format_num(cc_auc$estimate[1]), ", p = ", format_p(cc_auc$p_value[1]), "."),
    paste0("- NonCellCycle ECDF RMSE vs TGI AUC: Pearson r = ", format_num(nc_auc$estimate[1]), ", p = ", format_p(nc_auc$p_value[1]), "."),
    "",
    "## Files",
    "",
    "- `data/cellcycle_vs_noncellcycle_sample_shift_comparison.csv`",
    "- `stats/comparison_stats.csv`",
    "- `stats/endpoint_TGI_correlation_comparison.csv`",
    "- `figures/*.png` and `figures/*.pdf`"
  )
  write_text(lines, file.path(cfg$output_dir, "analysis_summary.md"))

  list(comparison = comp, comparison_stats = comparison_stats, endpoint_stats = endpoint_stats)
}

write_tex_ready_summary <- function(results, comparison, output_root) {
  get_dose_p <- function(result, comparison_name, stratified = FALSE) {
    x <- result$dose_tests
    hit <- x$comparison == comparison_name & x$stratified_by_initial_ploidy == stratified
    if (any(hit)) x$p_value[which(hit)[1]] else NA_real_
  }
  get_tgi_auc <- function(result, shift_metric = "ecdf_rmse", method = "pearson") {
    x <- result$tgi_stats
    hit <- x$sample_set == "treated" &
      x$shift_metric == shift_metric &
      x$tgi_measure == "TGI_percent_auc" &
      x$method == method
    if (!any(hit)) return(data.frame(estimate = NA_real_, p_value = NA_real_))
    x[which(hit)[1], c("estimate", "p_value"), drop = FALSE]
  }
  get_best_endpoint <- function(result) {
    x <- result$tgi_stats
    x <- x[
      x$sample_set == "treated" &
        x$shift_metric == "ecdf_rmse" &
        x$method == "pearson" &
        grepl("^TGI_percent_Day_", x$tgi_measure),
      ,
      drop = FALSE
    ]
    if (nrow(x) == 0) return(data.frame(tgi_day = NA_real_, estimate = NA_real_, p_value = NA_real_))
    x <- x[order(x$p_value), , drop = FALSE]
    x[1, c("tgi_day", "estimate", "p_value"), drop = FALSE]
  }
  get_ploidy_p <- function(result, sample_set, ploidy_measure, test) {
    x <- result$ploidy_stats
    hit <- x$sample_set == sample_set &
      x$shift_metric == "ecdf_rmse" &
      x$ploidy_measure == ploidy_measure &
      x$test == test
    if (any(hit)) x$p_value[which(hit)[1]] else NA_real_
  }
  comparison_p <- function(comparison_name, sample_set = "treated", test = "paired_t") {
    x <- comparison$comparison_stats
    hit <- x$sample_set == sample_set & x$comparison == comparison_name & x$test == test
    if (any(hit)) x$p_value[which(hit)[1]] else NA_real_
  }

  cc_auc <- get_tgi_auc(results$CellCycle)
  nc_auc <- get_tgi_auc(results$NonCellCycle)
  cc_auc_abs <- get_tgi_auc(results$CellCycle, "ecdf_mean_abs")
  cc_best_endpoint <- get_best_endpoint(results$CellCycle)
  nc_best_endpoint <- get_best_endpoint(results$NonCellCycle)

  lines <- c(
    "# Tex-Ready 04h Result Summary",
    "",
    "## Main question",
    "",
    "Does gemcitabine treatment shift the pseudotime distribution of cell-cycle-associated tumor cells relative to untreated controls, and does the magnitude of that shift relate to tumor-growth inhibition or ploidy?",
    "",
    "## Primary conclusions",
    "",
    paste0(
      "1. In cell-cycle-associated tumor cells, pseudotime distributions differ between untreated and treated tumors ",
      "(sample-aware ECDF permutation p = ", format_p(get_dose_p(results$CellCycle, "0_vs_30plus120", FALSE)),
      "; ploidy-stratified p = ", format_p(get_dose_p(results$CellCycle, "0_vs_30plus120", TRUE)), ")."
    ),
    paste0(
      "2. Both 30 mg/kg and 120 mg/kg differ from untreated controls in the cell-cycle compartment ",
      "(0 vs 30 p = ", format_p(get_dose_p(results$CellCycle, "0_vs_30", FALSE)),
      "; 0 vs 120 p = ", format_p(get_dose_p(results$CellCycle, "0_vs_120", FALSE)),
      "), while 30 and 120 mg/kg do not clearly separate from each other ",
      "(p = ", format_p(get_dose_p(results$CellCycle, "30_vs_120", FALSE)), ")."
    ),
    paste0(
      "3. Among treated tumors, larger cell-cycle pseudotime shifts are associated with higher AUC-based TGI ",
      "(ECDF RMSE Pearson r = ", format_num(cc_auc$estimate[1]),
      ", p = ", format_p(cc_auc$p_value[1]),
      "; ECDF mean absolute distance Pearson r = ", format_num(cc_auc_abs$estimate[1]),
      ", p = ", format_p(cc_auc_abs$p_value[1]), ")."
    ),
    paste0(
      "4. Endpoint TGI associations are time-dependent. The strongest cell-cycle endpoint association occurs around Day ",
      format_num(cc_best_endpoint$tgi_day[1], 3),
      " (Pearson r = ", format_num(cc_best_endpoint$estimate[1]),
      ", p = ", format_p(cc_best_endpoint$p_value[1]),
      "), supporting use of full-curve AUC rather than relying only on the final endpoint."
    ),
    paste0(
      "5. Non-cell-cycle tumor cells do not show the same untreated-versus-treated pseudotime shift ",
      "(sample-aware ECDF permutation p = ", format_p(get_dose_p(results$NonCellCycle, "0_vs_30plus120", FALSE)),
      "; ploidy-stratified p = ", format_p(get_dose_p(results$NonCellCycle, "0_vs_30plus120", TRUE)),
      ") and do not show a positive TGI AUC association ",
      "(Pearson r = ", format_num(nc_auc$estimate[1]),
      ", p = ", format_p(nc_auc$p_value[1]), ")."
    ),
    paste0(
      "6. Cell-cycle pseudotime-shift magnitude is not detectably explained by mean cell ploidy or initial ploidy ",
      "(mean ploidy all-sample p = ", format_p(get_ploidy_p(results$CellCycle, "all", "mean_cell_ploidy", "pearson")),
      "; treated-only p = ", format_p(get_ploidy_p(results$CellCycle, "treated", "mean_cell_ploidy", "pearson")),
      "; initial ploidy all-sample p = ", format_p(get_ploidy_p(results$CellCycle, "all", "initial_ploidy", "wilcoxon")),
      ")."
    ),
    "",
    "## AUC/TGI figures added in 04h",
    "",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_tumor_growth_raw_curves_by_dose_and_ploidy.*`: raw tumor-growth curves.",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_tumor_growth_baseline_adjusted_auc_curves.*`: baseline-adjusted curves used for AUC TGI.",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_mean_baseline_adjusted_AUC_area_curves_by_dose.*`: mean baseline-adjusted AUC curves with area fill.",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_treated_sample_auc_vs_matched_control_auc.*`: treated curves overlaid with ploidy-matched 0 mg/kg control mean curves.",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_TGI_AUC_by_dose_and_ploidy.*`: sample-level AUC-based TGI by dose.",
    "- `CellCycle/figures/auc_growth_curves/CellCycleCells_endpoint_TGI_curves_over_days.*`: endpoint TGI trajectories across tumor-growth days.",
    "- `Comparison/figures/CellCycle_vs_NonCellCycle_endpoint_TGI_correlations_over_days.*`: endpoint-day correlation comparison between compartments.",
    "",
    "## Original-analysis compatibility outputs",
    "",
    "- `CellCycle/figures/original_analysis/CellCycleCells_dose_group_comparisons_density_by_dose_and_ploidy.*`",
    "- `CellCycle/figures/original_analysis/CellCycleCells_dose_group_comparisons_sample_ecdf_by_dose_and_ploidy.*`",
    "- `CellCycle/figures/original_analysis/CellCycleCells_dose_group_comparisons_sample_mean_pseudotime_by_dose.*`",
    "- `CellCycle/figures/original_analysis/CellCycleCells_pseudotime_shift_associations_TGI_AUC_ecdf_distance_vs_TGI_AUC.*`",
    "- `CellCycle/stats/CellCycleCells_pseudotimeAssociations_reported_stats.csv`",
    "",
    "## Suggested manuscript-style wording",
    "",
    paste0(
      "Gemcitabine treatment was associated with a significant, sample-aware shift in the pseudotime distribution of cell-cycle-associated tumor cells relative to untreated controls ",
      "(ECDF permutation p = ", format_p(get_dose_p(results$CellCycle, "0_vs_30plus120", FALSE)),
      "), and this remained significant after stratifying permutations by initial ploidy ",
      "(p = ", format_p(get_dose_p(results$CellCycle, "0_vs_30plus120", TRUE)),
      "). Both 30 and 120 mg/kg differed from untreated controls, but the two treated dose groups did not clearly differ from each other. Among treated tumors, the magnitude of the cell-cycle pseudotime shift was positively associated with AUC-based TGI ",
      "(Pearson r = ", format_num(cc_auc$estimate[1]),
      ", p = ", format_p(cc_auc$p_value[1]),
      "). By contrast, non-cell-cycle tumor cells did not show a comparable untreated-versus-treated pseudotime shift or a positive TGI AUC association. Baseline-adjusted tumor-growth curves and endpoint TGI trajectories support the AUC analysis by showing that response information is distributed across the growth curve rather than captured solely by the final endpoint."
    ),
    "",
    "## Comparison notes",
    "",
    paste0(
      "- Treated paired CellCycle vs NonCellCycle signed mean shift differs by paired t-test p = ",
      format_p(comparison_p("signed_mean_shift_CellCycle_vs_NonCellCycle", "treated", "paired_t")),
      "."
    ),
    paste0(
      "- Treated CellCycle-minus-NonCellCycle ECDF RMSE delta versus TGI AUC Pearson p = ",
      format_p(comparison_p("ecdf_rmse_delta_CellCycle_minus_NonCellCycle_vs_TGI_AUC", "treated", "pearson")),
      "."
    ),
    paste0(
      "- Best NonCellCycle endpoint association is Day ",
      format_num(nc_best_endpoint$tgi_day[1], 3),
      " with Pearson r = ", format_num(nc_best_endpoint$estimate[1]),
      ", p = ", format_p(nc_best_endpoint$p_value[1]),
      "; interpret as exploratory."
    )
  )
  write_text(lines, file.path(output_root, "tex_ready_summary.md"))
}

ensure_dir(output_root)
code_dir <- ensure_dir(file.path(output_root, "code"))
if (file.exists(script_path)) {
  invisible(file.copy(script_path, file.path(code_dir, basename(script_path)), overwrite = TRUE))
}

results <- list()
for (compartment_name in names(compartments)) {
  message("Running 04h analysis for ", compartment_name)
  results[[compartment_name]] <- run_compartment(compartment_name, compartments[[compartment_name]])
}

comparison <- paired_comparison(results)
write_tex_ready_summary(results, comparison, output_root)

all_reported_stats <- rbind(
  data.frame(compartment = "CellCycle", result_type = "dose_group_tests", results$CellCycle$dose_tests, stringsAsFactors = FALSE),
  data.frame(compartment = "NonCellCycle", result_type = "dose_group_tests", results$NonCellCycle$dose_tests, stringsAsFactors = FALSE)
)
write_csv(all_reported_stats, file.path(output_root, "reported_dose_group_tests.csv"))

all_tgi_stats <- rbind(
  data.frame(compartment = "CellCycle", results$CellCycle$tgi_stats, stringsAsFactors = FALSE),
  data.frame(compartment = "NonCellCycle", results$NonCellCycle$tgi_stats, stringsAsFactors = FALSE)
)
write_csv(all_tgi_stats, file.path(output_root, "reported_tgi_association_stats.csv"))

root_summary <- c(
  "# 04h Pseudotime and TGI Analysis",
  "",
  "This analysis consolidates the existing simple cell-cycle-associated pseudotime/TGI analysis, preserves original-analysis-compatible outputs, adds AUC/growth-curve visualizations for TGI, and applies the same logic to non-cell-cycle-associated tumor cells.",
  "",
  "The output is organized by compartment:",
  "",
  "- `CellCycle/`: cell-cycle-associated tumor cells from clusters 4c, 6, and 10.",
  "- `NonCellCycle/`: tumor cells from all clusters except 4c, 6, and 10.",
  "- `Comparison/`: paired sample-level comparison between the two compartments.",
  "",
  "The main distribution test is a sample-aware ECDF permutation test. The main TGI association uses the treated-only correlation between pseudotime-shift magnitude and `TGI_percent_auc`.",
  "",
  "AUC/TGI curve outputs are under each compartment's `figures/auc_growth_curves/` folder. Original pull-in analysis compatibility outputs are under each compartment's `figures/original_analysis/` folder.",
  "",
  "Key combined tables:",
  "",
  "- `reported_dose_group_tests.csv`",
  "- `reported_tgi_association_stats.csv`",
  "- `tex_ready_summary.md`"
)
write_text(root_summary, file.path(output_root, "README.md"))

message("04h pseudotime/TGI analysis completed.")
message("Output root: ", output_root)
