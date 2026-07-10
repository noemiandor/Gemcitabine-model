#!/usr/bin/env Rscript

# Unified 04h workflow. This entry point keeps the validated V4/V5 analysis
# bodies intact, gives every grouping method a concrete peer-level result
# directory, and runs the post-analysis interpretation/comparison stage.

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    x <- args[[i]]
    if (grepl("^--[^=]+=", x)) {
      key <- sub("^--([^=]+)=.*$", "\\1", x)
      out[[key]] <- sub("^--[^=]+=", "", x)
      i <- i + 1L
    } else if (grepl("^--", x)) {
      key <- sub("^--", "", x)
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- "TRUE"
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

cmd_args <- parse_args(commandArgs(trailingOnly = TRUE))
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = FALSE)
} else {
  normalizePath("Code/in-vivo/04h_pseudotime_TGI.R", mustWork = FALSE)
}
script_dir <- dirname(script_path)
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

arg_or <- function(name, default) {
  value <- cmd_args[[name]]
  if (!is.null(value) && nzchar(value)) value else default
}

normalize_choice <- function(value) tolower(gsub("-", "_", trimws(value)))

grouping <- normalize_choice(arg_or("grouping", "all"))
if (grouping == "initial") grouping <- "initial_ploidy"
if (!grouping %in% c("initial_ploidy", "etp", "all")) {
  stop("--grouping must be one of: initial_ploidy, etp, all", call. = FALSE)
}

etp_methods <- c(
  fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
  boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
  reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24"
)
etp_aliases <- c(
  all = "all",
  etp_fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
  etp_boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
  etp_reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24",
  etp_methods
)
requested_etp <- normalize_choice(arg_or("etp_threshold", "all"))
selected_etp <- unname(etp_aliases[requested_etp])
if (length(selected_etp) != 1L || is.na(selected_etp)) {
  stop(
    "--etp_threshold must be one of: all, fixed_threshold_2_25, ",
    "boundary_stress_threshold_2_375, reference_balanced_threshold_2_24",
    call. = FALSE
  )
}
if (identical(selected_etp, "all")) selected_etp <- unname(etp_methods)

comparison_mode <- normalize_choice(arg_or("comparison", "auto"))
if (!comparison_mode %in% c("auto", "skip", "only")) {
  stop("--comparison must be one of: auto, skip, only", call. = FALSE)
}

input_root <- normalizePath(
  arg_or("input_root", file.path(repo_root, "Data", "in-vivo")),
  mustWork = FALSE
)
results_root <- normalizePath(
  arg_or(
    "results_root",
    "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI"
  ),
  mustWork = FALSE
)
seed <- as.integer(arg_or("seed", "1"))
n_perm <- as.integer(arg_or("n_perm", "10000"))
n_boot <- as.integer(arg_or("n_boot", "1000"))
n_downsample <- as.integer(arg_or("n_downsample", "200"))
run_optional_methods <- normalize_choice(arg_or("run_optional_methods", "true")) %in%
  c("true", "t", "1", "yes", "y")

dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
workflow_log_dir <- file.path(results_root, "workflow_logs")
dir.create(workflow_log_dir, recursive = TRUE, showWarnings = FALSE)

initial_method <- "initial_ploidy"
selected_methods <- character()
if (grouping %in% c("initial_ploidy", "all")) selected_methods <- c(selected_methods, initial_method)
if (grouping %in% c("etp", "all")) selected_methods <- c(selected_methods, selected_etp)
method_workers <- as.integer(arg_or("method_workers", as.character(length(selected_methods))))
if (is.na(method_workers) || method_workers < 1L) {
  stop("--method_workers must be a positive integer", call. = FALSE)
}
method_workers <- min(method_workers, max(1L, length(selected_methods)))

method_registry <- data.frame(
  method = c(initial_method, unname(etp_methods)),
  grouping = c("initial_ploidy", rep("ETP", length(etp_methods))),
  etp_threshold = c(NA_real_, 2.25, 2.375, 2.24),
  threshold_scheme = c(NA_character_, names(etp_methods)),
  output_dir = file.path(results_root, c(initial_method, unname(etp_methods))),
  stringsAsFactors = FALSE
)
write.csv(method_registry, file.path(results_root, "method_registry.csv"), row.names = FALSE, na = "NA")

rscript <- file.path(R.home("bin"), "Rscript")
v4_script <- file.path(script_dir, "04h_pseudotime_TGI_analysis_v4.R")
v5_script <- file.path(script_dir, "04h_pseudotime_TGI_analysis_v5.R")
comparison_script <- file.path(script_dir, "04h_v4_results_interpretation_and_comparison.R")
required_scripts <- c(v4_script, v5_script, comparison_script)
missing_scripts <- required_scripts[!file.exists(required_scripts)]
if (length(missing_scripts) > 0L) {
  stop("Missing workflow scripts:\n", paste(missing_scripts, collapse = "\n"), call. = FALSE)
}

quote_args <- function(args) vapply(args, shQuote, character(1L), USE.NAMES = FALSE)

run_child <- function(label, script, args) {
  log_path <- file.path(workflow_log_dir, paste0(label, ".log"))
  message("Running ", label)
  message("  output log: ", log_path)
  status <- system2(
    rscript,
    args = c(shQuote(script), quote_args(args)),
    stdout = log_path,
    stderr = log_path,
    wait = TRUE
  )
  data.frame(
    method = label,
    exit_status = as.integer(status),
    log = log_path,
    completed_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
}

psock_lapply <- function(items, fun, workers, exports) {
  if (workers <= 1L || length(items) <= 1L) return(lapply(items, fun))
  cluster <- parallel::makePSOCKcluster(workers)
  on.exit(parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterExport(cluster, exports, envir = .GlobalEnv)
  parallel::parLapply(cluster, items, fun)
}

# -----------------------------------------------------------------------------
# Unified mean-ETP/TGI post-processing
# -----------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The unified mean-ETP analysis requires the ggplot2 package.", call. = FALSE)
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

endpoint_tgi_cols <- function(nms) {
  cols <- grep("^TGI_percent_Day_[0-9]+$", nms, value = TRUE)
  days <- safe_num(sub("^TGI_percent_Day_", "", cols))
  cols[order(days)]
}

safe_cor_result <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(
      n = length(x), estimate = NA_real_, p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  data.frame(
    n = length(x), estimate = unname(ct$estimate), p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}

.perm_cache <- new.env(parent = emptyenv())

perm_index_matrix <- function(n) {
  key <- as.character(n)
  if (exists(key, envir = .perm_cache, inherits = FALSE)) {
    return(get(key, envir = .perm_cache, inherits = FALSE))
  }
  total <- as.integer(round(factorial(n)))
  mat <- matrix(NA_integer_, nrow = total, ncol = n)
  j <- 0L
  recurse <- function(prefix, rest) {
    if (length(rest) == 0L) {
      j <<- j + 1L
      mat[j, ] <<- prefix
    } else {
      for (k in seq_along(rest)) recurse(c(prefix, rest[k]), rest[-k])
    }
  }
  recurse(integer(0), seq_len(n))
  assign(key, mat, envir = .perm_cache)
  mat
}

permutation_cor_p <- function(x, y, method = "pearson", n_perm = 10000L, exact = FALSE) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(
      n = n,
      estimate = NA_real_,
      asymptotic_p = NA_real_,
      permutation_p_two_sided = NA_real_,
      permutation_p_positive = NA_real_,
      n_permutations = NA_integer_,
      permutation_mode = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  observed <- suppressWarnings(stats::cor(x, y, method = method))
  asym <- safe_cor_result(x, y, method)
  if (isTRUE(exact) && n <= 9L && factorial(n) <= 400000) {
    idx <- perm_index_matrix(n)
    perm_est <- vapply(
      seq_len(nrow(idx)),
      function(i) suppressWarnings(stats::cor(x, y[idx[i, ]], method = method)),
      numeric(1L)
    )
    mode <- "exact_TGI_label_enumeration"
  } else {
    perm_est <- replicate(
      n_perm,
      suppressWarnings(stats::cor(x, sample(y), method = method))
    )
    mode <- "monte_carlo_TGI_label_permutation"
  }
  data.frame(
    n = n,
    estimate = observed,
    asymptotic_p = asym$p_value,
    permutation_p_two_sided = mean(abs(perm_est) >= abs(observed) - 1e-15, na.rm = TRUE),
    permutation_p_positive = mean(perm_est >= observed - 1e-15, na.rm = TRUE),
    n_permutations = length(perm_est),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

dose_cols_all <- c(
  "0mg/kg" = "#666666",
  "30mg/kg" = "#d95f02",
  "120mg/kg" = "#1b9e77",
  "treated" = "#377eb8"
)

mean_etp_plot_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
    legend.key.height = grid::unit(0.45, "cm"),
    aspect.ratio = 1
  )

save_mean_etp_plot <- function(plot, output_dir, filename, width = 6.8, height = 6.8) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  suppressMessages(ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 300
  ))
  suppressMessages(ggplot2::ggsave(
    file.path(output_dir, paste0(filename, ".pdf")),
    plot,
    width = width,
    height = height
  ))
  invisible(file.path(output_dir, paste0(filename, ".pdf")))
}

method_group_spec <- function(method) {
  if (identical(method, initial_method)) {
    list(column = "initial_ploidy", label = "initial_ploidy")
  } else {
    list(column = "end_timepoint_ploidy_group", label = "end_timepoint_ploidy_group")
  }
}

read_mean_etp_source <- function(method, compartment) {
  path <- file.path(
    results_root,
    method,
    "A_Pseudotime",
    "A1_TreatmentDistributionShift",
    compartment,
    "tables",
    paste0(compartment, "_shift_metrics_equal_sample_ref.csv")
  )
  if (!file.exists(path)) stop("Missing mean-ETP source table: ", path, call. = FALSE)
  dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "sample_id", "dose", "dose_mg", "sample_mean_endpoint_ploidy", "TGI_percent_auc"
  )
  group_spec <- method_group_spec(method)
  required <- c(required, group_spec$column)
  missing <- setdiff(required, names(dat))
  if (length(missing) > 0L) {
    stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(dat$sample_id)) {
    stop("Duplicate sample_id values in ", path, call. = FALSE)
  }
  dat
}

run_mean_etp_tgi_associations <- function(sample_table, method, compartment) {
  tgi_cols <- unique(c("TGI_percent_auc", endpoint_tgi_cols(names(sample_table))))
  rows <- vector("list", length(tgi_cols))
  for (i in seq_along(tgi_cols)) {
    tgi_col <- tgi_cols[[i]]
    use_exact <- identical(tgi_col, "TGI_percent_auc")
    pearson <- permutation_cor_p(
      sample_table$sample_mean_endpoint_ploidy,
      sample_table[[tgi_col]],
      "pearson",
      n_perm,
      exact = use_exact
    )
    spearman <- permutation_cor_p(
      sample_table$sample_mean_endpoint_ploidy,
      sample_table[[tgi_col]],
      "spearman",
      n_perm,
      exact = use_exact
    )
    rows[[i]] <- data.frame(
      method = method,
      compartment = compartment,
      sample_set = "treated",
      predictor = "sample_mean_endpoint_ploidy",
      predictor_label = "sample_mean_ETP",
      tgi_measure = tgi_col,
      tgi_day = ifelse(
        grepl("^TGI_percent_Day_", tgi_col),
        safe_num(sub("^TGI_percent_Day_", "", tgi_col)),
        NA_real_
      ),
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p_asymptotic = pearson$asymptotic_p,
      pearson_p_permutation_two_sided = pearson$permutation_p_two_sided,
      pearson_p_permutation_positive = pearson$permutation_p_positive,
      pearson_n_permutations = pearson$n_permutations,
      pearson_permutation_mode = pearson$permutation_mode,
      spearman_rho = spearman$estimate,
      spearman_p_asymptotic = spearman$asymptotic_p,
      spearman_p_permutation_two_sided = spearman$permutation_p_two_sided,
      spearman_p_permutation_positive = spearman$permutation_p_positive,
      spearman_n_permutations = spearman$n_permutations,
      spearman_permutation_mode = spearman$permutation_mode,
      multiplicity_family = paste0(compartment, "_mean_ETP_TGI"),
      analysis_role = ifelse(
        compartment == "CellCycle" && tgi_col == "TGI_percent_auc",
        "parallel_mean_ETP_AUC",
        "exploratory_mean_ETP"
      ),
      pre_specified_primary = FALSE,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  out$pearson_q_exploratory <- stats::p.adjust(
    out$pearson_p_permutation_two_sided,
    method = "BH"
  )
  out$spearman_q_exploratory <- stats::p.adjust(
    out$spearman_p_permutation_two_sided,
    method = "BH"
  )
  out
}

plot_mean_etp_tgi_association <- function(sample_table, method, compartment, figure_dir) {
  local <- sample_table[
    sample_table$analysis_included_AUC & !is.na(sample_table$method_group),
    ,
    drop = FALSE
  ]
  dose_levels <- unique(local$dose[order(local$dose_mg)])
  group_spec <- method_group_spec(method)
  compartment_label <- if (identical(compartment, "CellCycle")) {
    "Cell-cycle-associated tumor cells"
  } else {
    "Non-cell-cycle-associated tumor cells"
  }
  p <- ggplot2::ggplot(
    local,
    ggplot2::aes(x = sample_mean_endpoint_ploidy, y = TGI_percent_auc, color = dose)
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_smooth(
      ggplot2::aes(group = 1),
      method = "lm",
      se = TRUE,
      color = "black",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(ggplot2::aes(shape = method_group), size = 2.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = sample_id),
      size = 2.4,
      nudge_y = 2.5,
      check_overlap = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = dose_cols_all[names(dose_cols_all) %in% dose_levels],
      name = "dose"
    ) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.10, 0.10))
    ) +
    ggplot2::labs(
      title = paste0(compartment_label, ": AUC TGI vs sample mean ETP"),
      x = "Sample mean ETP",
      y = "AUC-based TGI (%)",
      shape = group_spec$label
    ) +
    mean_etp_plot_theme
  save_mean_etp_plot(
    p,
    figure_dir,
    paste0(compartment, "_TGI_AUC_vs_mean_ETP")
  )
}

make_mean_etp_sample_table <- function(source, method, compartment) {
  group_spec <- method_group_spec(method)
  treated <- source[source$dose_mg > 0, , drop = FALSE]
  endpoint_cols <- endpoint_tgi_cols(names(treated))
  keep <- unique(c(
    "sample_id", "dose", "dose_mg", group_spec$column,
    "sample_mean_endpoint_ploidy", "TGI_percent_auc", endpoint_cols
  ))
  out <- treated[, keep, drop = FALSE]
  out$method <- method
  out$compartment <- compartment
  out$sample_set <- "treated"
  out$method_group_type <- group_spec$label
  out$method_group <- out[[group_spec$column]]
  out$analysis_included_AUC <- is.finite(out$sample_mean_endpoint_ploidy) &
    is.finite(out$TGI_percent_auc)
  out$exclusion_reason_AUC <- ifelse(
    out$analysis_included_AUC,
    "included",
    ifelse(
      !is.finite(out$sample_mean_endpoint_ploidy),
      "missing_mean_ETP",
      "missing_AUC_TGI"
    )
  )
  out$n_available_TGI_measures <- rowSums(is.finite(as.matrix(out[, c("TGI_percent_auc", endpoint_cols), drop = FALSE])))
  out <- out[order(out$sample_id), , drop = FALSE]
  rownames(out) <- NULL
  leading <- c(
    "method", "compartment", "sample_set", "sample_id", "dose", "dose_mg",
    "method_group_type", "method_group", "sample_mean_endpoint_ploidy",
    "TGI_percent_auc", "analysis_included_AUC", "exclusion_reason_AUC",
    "n_available_TGI_measures"
  )
  out[, c(leading, setdiff(names(out), leading)), drop = FALSE]
}

leave_one_out_mean_etp <- function(sample_table, method) {
  local <- sample_table[sample_table$analysis_included_AUC, , drop = FALSE]
  full <- safe_cor_result(local$sample_mean_endpoint_ploidy, local$TGI_percent_auc, "pearson")
  rows <- lapply(local$sample_id, function(sample_id) {
    dat <- local[local$sample_id != sample_id, , drop = FALSE]
    p <- safe_cor_result(dat$sample_mean_endpoint_ploidy, dat$TGI_percent_auc, "pearson")
    s <- safe_cor_result(dat$sample_mean_endpoint_ploidy, dat$TGI_percent_auc, "spearman")
    data.frame(
      method = method,
      predictor = "sample_mean_endpoint_ploidy",
      omitted_sample_id = sample_id,
      n = p$n,
      pearson_r = p$estimate,
      pearson_p = p$p_value,
      spearman_rho = s$estimate,
      spearman_p = s$p_value,
      sign_matches_full_pearson = sign(p$estimate) == sign(full$estimate),
      delta_pearson_r_from_full = p$estimate - full$estimate,
      influential_flag = is.finite(p$estimate) && is.finite(full$estimate) &&
        (sign(p$estimate) != sign(full$estimate) || abs(p$estimate - full$estimate) > 0.20),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bootstrap_mean_etp <- function(sample_table, method) {
  local <- sample_table[sample_table$analysis_included_AUC, , drop = FALSE]
  rows <- lapply(seq_len(n_boot), function(i) {
    idx <- sample(seq_len(nrow(local)), replace = TRUE)
    dat <- local[idx, , drop = FALSE]
    p <- safe_cor_result(dat$sample_mean_endpoint_ploidy, dat$TGI_percent_auc, "pearson")
    s <- safe_cor_result(dat$sample_mean_endpoint_ploidy, dat$TGI_percent_auc, "spearman")
    data.frame(
      method = method,
      predictor = "sample_mean_endpoint_ploidy",
      iteration = i,
      n = p$n,
      pearson_r = p$estimate,
      pearson_p = p$p_value,
      spearman_rho = s$estimate,
      spearman_p = s$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_mean_etp_robustness <- function(loo, boot, auc_stats, method) {
  finite_quantile <- function(x, p) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return(NA_real_)
    as.numeric(stats::quantile(x, p, names = FALSE, na.rm = TRUE))
  }
  data.frame(
    method = method,
    predictor = "sample_mean_endpoint_ploidy",
    full_n = auc_stats$n[[1L]],
    full_pearson_r = auc_stats$pearson_r[[1L]],
    full_pearson_asymptotic_p = auc_stats$pearson_p_asymptotic[[1L]],
    full_pearson_perm_p = auc_stats$pearson_p_permutation_two_sided[[1L]],
    full_pearson_q = auc_stats$pearson_q_exploratory[[1L]],
    full_spearman_rho = auc_stats$spearman_rho[[1L]],
    full_spearman_asymptotic_p = auc_stats$spearman_p_asymptotic[[1L]],
    full_spearman_perm_p = auc_stats$spearman_p_permutation_two_sided[[1L]],
    full_spearman_q = auc_stats$spearman_q_exploratory[[1L]],
    leave_one_out_all_same_sign = all(loo$sign_matches_full_pearson, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(loo$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_r = max(loo$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(loo$pearson_p, na.rm = TRUE),
    leave_one_out_max_abs_delta_r = max(abs(loo$delta_pearson_r_from_full), na.rm = TRUE),
    n_influential_flags = sum(loo$influential_flag, na.rm = TRUE),
    n_bootstrap_requested = n_boot,
    n_bootstrap_finite_pearson = sum(is.finite(boot$pearson_r)),
    n_bootstrap_finite_spearman = sum(is.finite(boot$spearman_rho)),
    bootstrap_pearson_r_ci025 = finite_quantile(boot$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = finite_quantile(boot$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = finite_quantile(boot$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = finite_quantile(boot$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

plot_mean_etp_robustness <- function(loo, boot, robustness_dir) {
  p1 <- ggplot2::ggplot(
    loo,
    ggplot2::aes(
      x = stats::reorder(omitted_sample_id, pearson_r),
      y = pearson_r,
      color = influential_flag
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    ggplot2::geom_point(size = 2.6) +
    ggplot2::coord_flip() +
    ggplot2::scale_color_manual(
      values = c("FALSE" = "#377eb8", "TRUE" = "#e41a1c"),
      name = "influential"
    ) +
    ggplot2::labs(
      title = "CellCycle mean ETP-TGI association: leave-one-out Pearson r",
      x = "omitted sample",
      y = "Pearson r"
    ) +
    mean_etp_plot_theme
  save_mean_etp_plot(
    p1,
    robustness_dir,
    "CellCycle_mean_ETP_TGI_leave_one_out_forest"
  )

  p2 <- ggplot2::ggplot(boot, ggplot2::aes(x = pearson_r)) +
    ggplot2::geom_histogram(bins = 35, fill = "#377eb8", color = "white", alpha = 0.85) +
    ggplot2::geom_vline(xintercept = 0, color = "grey40", linewidth = 0.35) +
    ggplot2::labs(
      title = "CellCycle mean ETP-TGI association: bootstrap Pearson r",
      x = "bootstrap Pearson r",
      y = "iterations"
    ) +
    mean_etp_plot_theme
  save_mean_etp_plot(
    p2,
    robustness_dir,
    "CellCycle_mean_ETP_TGI_bootstrap_distribution"
  )
}

postprocess_mean_etp_method <- function(method) {
  method_root <- file.path(results_root, method)
  compartment_results <- list()
  for (compartment in c("CellCycle", "NonCellCycle")) {
    source <- read_mean_etp_source(method, compartment)
    sample_table <- make_mean_etp_sample_table(source, method, compartment)
    set.seed(seed + ifelse(compartment == "CellCycle", 1101L, 1201L))
    stats_table <- run_mean_etp_tgi_associations(sample_table, method, compartment)
    compartment_root <- file.path(
      method_root,
      "A_Pseudotime",
      "A1_TreatmentDistributionShift",
      compartment
    )
    write_csv(
      sample_table,
      file.path(
        compartment_root,
        "tables",
        paste0(compartment, "_mean_ETP_TGI_sample_table.csv")
      )
    )
    write_csv(
      stats_table,
      file.path(
        compartment_root,
        "stats",
        paste0(compartment, "_TGI_associations_mean_ETP.csv")
      )
    )
    plot_mean_etp_tgi_association(
      sample_table,
      method,
      compartment,
      file.path(compartment_root, "figures")
    )
    compartment_results[[compartment]] <- list(
      samples = sample_table,
      stats = stats_table
    )
  }

  robustness_dir <- file.path(method_root, "C_SupportingAnalyses", "C1_Robustness")
  set.seed(seed + 2101L)
  loo <- leave_one_out_mean_etp(compartment_results$CellCycle$samples, method)
  set.seed(seed + 2201L)
  boot <- bootstrap_mean_etp(compartment_results$CellCycle$samples, method)
  auc_stats <- compartment_results$CellCycle$stats[
    compartment_results$CellCycle$stats$tgi_measure == "TGI_percent_auc",
    ,
    drop = FALSE
  ]
  robustness <- summarize_mean_etp_robustness(loo, boot, auc_stats, method)
  write_csv(loo, file.path(robustness_dir, "CellCycle_mean_ETP_TGI_leave_one_out.csv"))
  write_csv(boot, file.path(robustness_dir, "CellCycle_mean_ETP_TGI_bootstrap.csv"))
  write_csv(
    robustness,
    file.path(robustness_dir, "CellCycle_mean_ETP_TGI_robustness_summary.csv")
  )
  write_csv(
    data.frame(
      method = method,
      analysis = "mean_ETP_cell_count_downsampling",
      status = "not_run_not_applicable_to_existing_ECDF_downsampling_design",
      reason = paste(
        "The existing downsampling routine targets pseudotime ECDF shift.",
        "A cell-resampled mean-ETP estimand would require a separate raw-cell resampling design."
      ),
      stringsAsFactors = FALSE
    ),
    file.path(robustness_dir, "CellCycle_mean_ETP_TGI_downsampling_status.csv")
  )
  plot_mean_etp_robustness(loo, boot, robustness_dir)
  write_csv(
    rbind(
      compartment_results$CellCycle$stats[
        compartment_results$CellCycle$stats$tgi_measure == "TGI_percent_auc",
        ,
        drop = FALSE
      ],
      compartment_results$NonCellCycle$stats[
        compartment_results$NonCellCycle$stats$tgi_measure == "TGI_percent_auc",
        ,
        drop = FALSE
      ]
    ),
    file.path(method_root, "mean_ETP_TGI_analysis_summary.csv")
  )
  invisible(TRUE)
}

run_mean_etp_postprocess <- function(method) {
  tryCatch({
    message("Running unified mean-ETP post-processing: ", method)
    postprocess_mean_etp_method(method)
    data.frame(
      method = method,
      exit_status = 0L,
      error = NA_character_,
      completed_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      method = method,
      exit_status = 1L,
      error = conditionMessage(e),
      completed_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
  })
}

run_mean_etp_postprocess_child <- function(method) {
  log_path <- file.path(workflow_log_dir, paste0("mean_ETP_postprocess_", method, ".log"))
  status <- system2(
    rscript,
    args = c(
      shQuote(script_path),
      quote_args(c(
        paste0("--postprocess_method=", method),
        paste0("--results_root=", results_root),
        paste0("--seed=", seed),
        paste0("--n_perm=", n_perm),
        paste0("--n_boot=", n_boot)
      ))
    ),
    stdout = log_path,
    stderr = log_path,
    wait = TRUE
  )
  data.frame(
    method = method,
    exit_status = as.integer(status),
    error = ifelse(status == 0L, NA_character_, paste0("See ", log_path)),
    completed_at = as.character(Sys.time()),
    stringsAsFactors = FALSE
  )
}

mean_etp_required_files <- function(method) {
  method_root <- file.path(results_root, method)
  compartment_files <- unlist(lapply(c("CellCycle", "NonCellCycle"), function(compartment) {
    root <- file.path(
      method_root,
      "A_Pseudotime",
      "A1_TreatmentDistributionShift",
      compartment
    )
    c(
      file.path(root, "figures", paste0(compartment, "_TGI_AUC_vs_mean_ETP.pdf")),
      file.path(root, "figures", paste0(compartment, "_TGI_AUC_vs_mean_ETP.png")),
      file.path(root, "tables", paste0(compartment, "_mean_ETP_TGI_sample_table.csv")),
      file.path(root, "stats", paste0(compartment, "_TGI_associations_mean_ETP.csv"))
    )
  }), use.names = FALSE)
  robustness_root <- file.path(method_root, "C_SupportingAnalyses", "C1_Robustness")
  b1_root <- file.path(method_root, "B_TGI_GrowthResponse", "B1_AUC_TGI_Association", "CellCycle")
  c(
    compartment_files,
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_leave_one_out.csv"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_bootstrap.csv"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_robustness_summary.csv"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_downsampling_status.csv"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_leave_one_out_forest.pdf"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_leave_one_out_forest.png"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_bootstrap_distribution.pdf"),
    file.path(robustness_root, "CellCycle_mean_ETP_TGI_bootstrap_distribution.png"),
    file.path(b1_root, "stats", "CellCycle_AUC_TGI_mean_ETP_ploidy_dose_models.csv"),
    file.path(b1_root, "stats", "CellCycle_AUC_TGI_mean_ETP_ploidy_dose_model_summaries.csv"),
    file.path(method_root, "mean_ETP_TGI_analysis_summary.csv")
  )
}

consolidate_mean_etp_comparison <- function(methods) {
  comparison_dir <- file.path(results_root, "method_comparison")
  dir.create(comparison_dir, recursive = TRUE, showWarnings = FALSE)

  association_rows <- list()
  robustness_rows <- list()
  sample_tables <- list()
  for (method in methods) {
    for (compartment in c("CellCycle", "NonCellCycle")) {
      compartment_root <- file.path(
        results_root,
        method,
        "A_Pseudotime",
        "A1_TreatmentDistributionShift",
        compartment
      )
      association_rows[[length(association_rows) + 1L]] <- read.csv(
        file.path(
          compartment_root,
          "stats",
          paste0(compartment, "_TGI_associations_mean_ETP.csv")
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      key <- paste(method, compartment, sep = "::")
      sample_tables[[key]] <- read.csv(
        file.path(
          compartment_root,
          "tables",
          paste0(compartment, "_mean_ETP_TGI_sample_table.csv")
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
    robustness_rows[[length(robustness_rows) + 1L]] <- read.csv(
      file.path(
        results_root,
        method,
        "C_SupportingAnalyses",
        "C1_Robustness",
        "CellCycle_mean_ETP_TGI_robustness_summary.csv"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
  associations <- do.call(rbind, association_rows)
  robustness <- do.call(rbind, robustness_rows)
  write_csv(
    associations,
    file.path(comparison_dir, "mean_ETP_TGI_associations_all_methods.csv")
  )
  write_csv(
    robustness,
    file.path(comparison_dir, "mean_ETP_TGI_robustness_all_methods.csv")
  )

  audit_rows <- list()
  for (compartment in c("CellCycle", "NonCellCycle")) {
    reference <- sample_tables[[paste(initial_method, compartment, sep = "::")]]
    reference <- reference[, c("sample_id", "sample_mean_endpoint_ploidy", "TGI_percent_auc")]
    names(reference)[-1L] <- paste0(names(reference)[-1L], "_reference")
    for (method in methods) {
      candidate <- sample_tables[[paste(method, compartment, sep = "::")]]
      candidate <- candidate[, c("sample_id", "sample_mean_endpoint_ploidy", "TGI_percent_auc")]
      merged <- merge(reference, candidate, by = "sample_id", all = TRUE)
      same_ids <- setequal(reference$sample_id, candidate$sample_id) &&
        nrow(merged) == nrow(reference) && nrow(merged) == nrow(candidate)
      mean_etp_delta <- abs(
        merged$sample_mean_endpoint_ploidy_reference - merged$sample_mean_endpoint_ploidy
      )
      tgi_delta <- abs(merged$TGI_percent_auc_reference - merged$TGI_percent_auc)
      max_mean_etp_delta <- if (all(is.na(mean_etp_delta))) NA_real_ else max(mean_etp_delta, na.rm = TRUE)
      max_tgi_delta <- if (all(is.na(tgi_delta))) NA_real_ else max(tgi_delta, na.rm = TRUE)
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        method = method,
        compartment = compartment,
        reference_method = initial_method,
        n_reference = nrow(reference),
        n_candidate = nrow(candidate),
        sample_ids_identical = same_ids,
        max_abs_mean_ETP_difference = max_mean_etp_delta,
        max_abs_AUC_TGI_difference = max_tgi_delta,
        invariant_values_pass = same_ids &&
          is.finite(max_mean_etp_delta) && max_mean_etp_delta < 1e-12 &&
          is.finite(max_tgi_delta) && max_tgi_delta < 1e-12,
        stringsAsFactors = FALSE
      )
    }
  }
  audit <- do.call(rbind, audit_rows)
  write_csv(
    audit,
    file.path(comparison_dir, "mean_ETP_TGI_sample_invariance_audit.csv")
  )

  metric_specs <- list(
    c("cellcycle_mean_ETP_AUC_n", "CellCycle mean ETP vs AUC-TGI n", "assoc", "n"),
    c("cellcycle_mean_ETP_AUC_pearson_r", "CellCycle mean ETP vs AUC-TGI Pearson r", "assoc", "pearson_r"),
    c("cellcycle_mean_ETP_AUC_pearson_asymptotic_p", "CellCycle mean ETP vs AUC-TGI Pearson asymptotic p", "assoc", "pearson_p_asymptotic"),
    c("cellcycle_mean_ETP_AUC_pearson_perm_two_sided_p", "CellCycle mean ETP vs AUC-TGI Pearson TGI-label permutation p", "assoc", "pearson_p_permutation_two_sided"),
    c("cellcycle_mean_ETP_AUC_pearson_perm_positive_p", "CellCycle mean ETP vs AUC-TGI Pearson positive-direction permutation p", "assoc", "pearson_p_permutation_positive"),
    c("cellcycle_mean_ETP_AUC_pearson_q", "CellCycle mean ETP vs AUC-TGI Pearson BH q", "assoc", "pearson_q_exploratory"),
    c("cellcycle_mean_ETP_AUC_spearman_rho", "CellCycle mean ETP vs AUC-TGI Spearman rho", "assoc", "spearman_rho"),
    c("cellcycle_mean_ETP_AUC_spearman_asymptotic_p", "CellCycle mean ETP vs AUC-TGI Spearman asymptotic p", "assoc", "spearman_p_asymptotic"),
    c("cellcycle_mean_ETP_AUC_spearman_perm_two_sided_p", "CellCycle mean ETP vs AUC-TGI Spearman TGI-label permutation p", "assoc", "spearman_p_permutation_two_sided"),
    c("cellcycle_mean_ETP_AUC_spearman_perm_positive_p", "CellCycle mean ETP vs AUC-TGI Spearman positive-direction permutation p", "assoc", "spearman_p_permutation_positive"),
    c("cellcycle_mean_ETP_AUC_spearman_q", "CellCycle mean ETP vs AUC-TGI Spearman BH q", "assoc", "spearman_q_exploratory"),
    c("cellcycle_mean_ETP_LOO_min_pearson_r", "CellCycle mean ETP leave-one-out minimum Pearson r", "robust", "leave_one_out_min_pearson_r"),
    c("cellcycle_mean_ETP_LOO_max_pearson_r", "CellCycle mean ETP leave-one-out maximum Pearson r", "robust", "leave_one_out_max_pearson_r"),
    c("cellcycle_mean_ETP_LOO_max_pearson_p", "CellCycle mean ETP leave-one-out maximum Pearson p", "robust", "leave_one_out_max_pearson_p"),
    c("cellcycle_mean_ETP_LOO_influential_n", "CellCycle mean ETP leave-one-out influential sample count", "robust", "n_influential_flags"),
    c("cellcycle_mean_ETP_bootstrap_pearson_ci025", "CellCycle mean ETP bootstrap Pearson r CI lower", "robust", "bootstrap_pearson_r_ci025"),
    c("cellcycle_mean_ETP_bootstrap_pearson_ci975", "CellCycle mean ETP bootstrap Pearson r CI upper", "robust", "bootstrap_pearson_r_ci975"),
    c("cellcycle_mean_ETP_bootstrap_spearman_ci025", "CellCycle mean ETP bootstrap Spearman rho CI lower", "robust", "bootstrap_spearman_rho_ci025"),
    c("cellcycle_mean_ETP_bootstrap_spearman_ci975", "CellCycle mean ETP bootstrap Spearman rho CI upper", "robust", "bootstrap_spearman_rho_ci975")
  )
  new_metrics <- do.call(rbind, lapply(metric_specs, function(spec) {
    source_name <- spec[[3L]]
    value_col <- spec[[4L]]
    row <- data.frame(
      metric_id = spec[[1L]],
      metric_label = spec[[2L]],
      analysis_family = ifelse(
        source_name == "assoc",
        "A1 continuous mean-ETP-TGI association",
        "C1 mean-ETP-TGI robustness"
      ),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    for (method in methods) {
      source <- if (source_name == "assoc") {
        associations[
          associations$method == method &
            associations$compartment == "CellCycle" &
            associations$tgi_measure == "TGI_percent_auc",
          ,
          drop = FALSE
        ]
      } else {
        robustness[robustness$method == method, , drop = FALSE]
      }
      row[[method]] <- if (nrow(source) > 0L) source[[value_col]][[1L]] else NA_real_
    }
    row$comparability_note <- if (source_name == "assoc") {
      "Directly comparable: treated sample IDs, continuous mean ETP, and AUC-TGI are invariant across methods."
    } else {
      "Directly comparable sample-level robustness analysis; method grouping is not used."
    }
    row
  }))

  master_path <- file.path(comparison_dir, "comparison_master_metrics.csv")
  if (!file.exists(master_path)) {
    stop("Missing comparison master metrics after comparison stage: ", master_path, call. = FALSE)
  }
  master <- read.csv(master_path, stringsAsFactors = FALSE, check.names = FALSE)
  master <- master[!master$metric_id %in% new_metrics$metric_id, , drop = FALSE]
  new_metrics <- new_metrics[, names(master), drop = FALSE]
  write_csv(rbind(master, new_metrics), master_path)

  auc_reference <- associations[
    associations$method == initial_method &
      associations$compartment == "CellCycle" &
      associations$tgi_measure == "TGI_percent_auc",
    ,
    drop = FALSE
  ]
  summary_lines <- c(
    "# Continuous mean-ETP/TGI comparison",
    "",
    "This unified post-processing analysis uses treated biological samples and sample mean endpoint ploidy as the continuous predictor.",
    "It is parallel to the A1 ECDF-RMSE/TGI association but does not use an untreated ECDF reference.",
    "",
    paste0("- CellCycle n: ", auc_reference$n[[1L]]),
    paste0("- Pearson r: ", format(auc_reference$pearson_r[[1L]], digits = 6)),
    paste0("- Pearson asymptotic p: ", format(auc_reference$pearson_p_asymptotic[[1L]], digits = 6)),
    paste0("- Pearson exact TGI-label permutation p: ", format(auc_reference$pearson_p_permutation_two_sided[[1L]], digits = 6)),
    paste0("- Spearman rho: ", format(auc_reference$spearman_rho[[1L]], digits = 6)),
    paste0("- Spearman exact TGI-label permutation p: ", format(auc_reference$spearman_p_permutation_two_sided[[1L]], digits = 6)),
    paste0("- All sample invariance checks passed: ", all(audit$invariant_values_pass)),
    "",
    "Unadjusted associations and robustness outputs are expected to be identical across methods because sample IDs, mean ETP, and TGI are identical. Method-specific grouping remains visible in plot shapes and in the existing B1 adjusted/interaction models.",
    "",
    "Complete outputs:",
    "",
    "- `mean_ETP_TGI_associations_all_methods.csv`",
    "- `mean_ETP_TGI_robustness_all_methods.csv`",
    "- `mean_ETP_TGI_sample_invariance_audit.csv`",
    "- `comparison_master_metrics.csv`"
  )
  writeLines(
    summary_lines,
    file.path(comparison_dir, "mean_ETP_TGI_comparison_summary.md")
  )
  invisible(TRUE)
}

postprocess_method_arg <- cmd_args[["postprocess_method"]]
if (!is.null(postprocess_method_arg) && nzchar(postprocess_method_arg)) {
  if (!postprocess_method_arg %in% method_registry$method) {
    stop(
      "--postprocess_method must name a registered method: ",
      paste(method_registry$method, collapse = ", "),
      call. = FALSE
    )
  }
  status <- run_mean_etp_postprocess(postprocess_method_arg)
  if (status$exit_status != 0L) {
    stop(status$error, call. = FALSE)
  }
  quit(save = "no", status = 0L)
}

common_args <- c(
  paste0("--input_root=", input_root),
  paste0("--seed=", seed),
  paste0("--n_perm=", n_perm),
  paste0("--n_boot=", n_boot),
  paste0("--n_downsample=", n_downsample),
  paste0("--run_optional_methods=", ifelse(run_optional_methods, "TRUE", "FALSE"))
)

run_method <- function(method) {
  tryCatch({
    if (identical(method, initial_method)) {
      method_args <- c(
        common_args,
        paste0("--method_name=", initial_method),
        paste0("--results_root=", file.path(results_root, initial_method))
      )
      run_child(method, v5_script, method_args)
    } else {
      method_args <- c(
        common_args,
        paste0("--etp_threshold=", method),
        paste0("--results_root=", results_root)
      )
      run_child(method, v4_script, method_args)
    }
  }, error = function(e) {
    data.frame(
      method = method,
      exit_status = 1L,
      log = paste0("wrapper error: ", conditionMessage(e)),
      completed_at = as.character(Sys.time()),
      stringsAsFactors = FALSE
    )
  })
}

status_rows <- list()
if (identical(comparison_mode, "only")) {
  prior_status_path <- file.path(results_root, "workflow_run_status.csv")
  prior_status <- if (file.exists(prior_status_path)) {
    read.csv(prior_status_path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    data.frame()
  }
  if (nrow(prior_status) > 0L) {
    prior_status <- prior_status[
      prior_status$method %in% method_registry$method & prior_status$exit_status == 0L,
      ,
      drop = FALSE
    ]
  }
  missing_status_methods <- setdiff(method_registry$method, prior_status$method)
  reconstructed <- lapply(missing_status_methods, function(method) {
    summary_path <- file.path(results_root, method, "analysis_summary.json")
    if (!file.exists(summary_path)) return(NULL)
    data.frame(
      method = method,
      exit_status = 0L,
      log = file.path(workflow_log_dir, paste0(method, ".log")),
      completed_at = as.character(file.info(summary_path)$mtime),
      stringsAsFactors = FALSE
    )
  })
  reconstructed <- reconstructed[!vapply(reconstructed, is.null, logical(1L))]
  if (length(reconstructed) > 0L) {
    prior_status <- rbind(prior_status, do.call(rbind, reconstructed))
  }
  if (nrow(prior_status) > 0L) {
    prior_status <- prior_status[
      match(method_registry$method, prior_status$method, nomatch = 0L),
      ,
      drop = FALSE
    ]
    status_rows <- lapply(seq_len(nrow(prior_status)), function(i) {
      prior_status[i, , drop = FALSE]
    })
  }
}
if (!identical(comparison_mode, "only")) {
  message(
    "Running ", length(selected_methods), " method task(s) with ",
    method_workers, " parallel worker(s)"
  )
  status_rows <- if (method_workers > 1L && length(selected_methods) > 1L) {
    psock_lapply(
      selected_methods,
      run_method,
      workers = method_workers,
      exports = c(
        "run_method", "run_child", "quote_args", "rscript", "workflow_log_dir",
        "initial_method", "common_args", "results_root", "v5_script", "v4_script"
      )
    )
  } else {
    lapply(selected_methods, run_method)
  }
  status_table <- do.call(rbind, status_rows)
  write.csv(status_table, file.path(results_root, "workflow_run_status.csv"), row.names = FALSE)
  failures <- status_table[status_table$exit_status != 0L, , drop = FALSE]
  if (nrow(failures) > 0L) {
    stop(
      "Method task(s) failed: ", paste(failures$method, collapse = ", "),
      ". See workflow_run_status.csv and method logs.",
      call. = FALSE
    )
  }
}

all_methods <- method_registry$method
base_ready <- file.exists(file.path(results_root, all_methods, "analysis_summary.json"))
names(base_ready) <- all_methods

if (identical(comparison_mode, "only") && !all(base_ready)) {
  missing <- names(base_ready)[!base_ready]
  stop(
    "Comparison-only mode is missing completed methods: ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

methods_to_postprocess <- if (identical(comparison_mode, "skip")) {
  selected_methods
} else {
  unique(c(selected_methods, names(base_ready)[base_ready]))
}
methods_to_postprocess <- methods_to_postprocess[
  file.exists(file.path(results_root, methods_to_postprocess, "analysis_summary.json"))
]

if (length(methods_to_postprocess) > 0L) {
  postprocess_workers <- min(method_workers, length(methods_to_postprocess))
  message(
    "Running unified mean-ETP post-processing for ",
    length(methods_to_postprocess),
    " method(s) with ",
    postprocess_workers,
    " parallel worker(s)"
  )
  postprocess_rows <- if (postprocess_workers > 1L) {
    psock_lapply(
      methods_to_postprocess,
      run_mean_etp_postprocess_child,
      workers = postprocess_workers,
      exports = c(
        "run_mean_etp_postprocess_child", "quote_args", "rscript", "script_path",
        "workflow_log_dir", "results_root", "seed", "n_perm", "n_boot"
      )
    )
  } else {
    lapply(methods_to_postprocess, run_mean_etp_postprocess_child)
  }
  valid_rows <- vapply(postprocess_rows, is.data.frame, logical(1L))
  if (!all(valid_rows)) {
    missing_methods <- methods_to_postprocess[!valid_rows]
    postprocess_rows[!valid_rows] <- lapply(missing_methods, function(method) {
      data.frame(
        method = method,
        exit_status = 1L,
        error = "Post-processing child did not return a status row.",
        completed_at = as.character(Sys.time()),
        stringsAsFactors = FALSE
      )
    })
  }
  postprocess_status <- do.call(rbind, postprocess_rows)
  write.csv(
    postprocess_status,
    file.path(results_root, "mean_ETP_postprocess_status.csv"),
    row.names = FALSE,
    na = "NA"
  )
  postprocess_failures <- postprocess_status[
    postprocess_status$exit_status != 0L,
    ,
    drop = FALSE
  ]
  if (nrow(postprocess_failures) > 0L) {
    stop(
      "Unified mean-ETP post-processing failed for: ",
      paste(postprocess_failures$method, collapse = ", "),
      ". See mean_ETP_postprocess_status.csv.",
      call. = FALSE
    )
  }
}

base_ready <- file.exists(file.path(results_root, all_methods, "analysis_summary.json"))
mean_etp_ready <- vapply(
  all_methods,
  function(method) all(file.exists(mean_etp_required_files(method))),
  logical(1L)
)
names(base_ready) <- names(mean_etp_ready) <- all_methods
all_method_outputs_ready <- all(base_ready & mean_etp_ready)
run_comparison <- FALSE
if (identical(comparison_mode, "only")) {
  if (!all_method_outputs_ready) {
    missing_base <- all_methods[!base_ready]
    missing_mean <- unlist(lapply(all_methods[!mean_etp_ready], function(method) {
      required <- mean_etp_required_files(method)
      required[!file.exists(required)]
    }), use.names = FALSE)
    stop(
      "Comparison-only mode is missing required outputs:\n",
      paste(c(missing_base, missing_mean), collapse = "\n"),
      call. = FALSE
    )
  }
  run_comparison <- TRUE
} else if (identical(comparison_mode, "auto") && all_method_outputs_ready) {
  run_comparison <- TRUE
} else if (identical(comparison_mode, "auto")) {
  message("Comparison skipped because the complete four-method result set is not available.")
}

if (run_comparison) {
  comparison_status <- run_child(
    "method_comparison",
    comparison_script,
    c(paste0("--results_root=", results_root))
  )
  if (comparison_status$exit_status == 0L) {
    consolidation_error <- tryCatch({
      consolidate_mean_etp_comparison(all_methods)
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(consolidation_error)) {
      comparison_status$exit_status <- 1L
      writeLines(
        consolidation_error,
        file.path(workflow_log_dir, "mean_ETP_comparison_consolidation_error.log")
      )
    }
  }
  status_rows[[length(status_rows) + 1L]] <- comparison_status
  if (comparison_status$exit_status != 0L) {
    status_table <- do.call(rbind, status_rows)
    write.csv(status_table, file.path(results_root, "workflow_run_status.csv"), row.names = FALSE)
    stop("Method comparison failed. See ", comparison_status$log, call. = FALSE)
  }
}

status_table <- if (length(status_rows) > 0L) do.call(rbind, status_rows) else data.frame()
if (nrow(status_table) > 0L) {
  write.csv(status_table, file.path(results_root, "workflow_run_status.csv"), row.names = FALSE)
}

ppt_files <- list.files(results_root, pattern = "\\.(ppt|pptx)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
if (length(ppt_files) > 0L) {
  stop("Unexpected PPT/PPTX outputs were produced:\n", paste(ppt_files, collapse = "\n"), call. = FALSE)
}

message("Completed unified 04h pseudotime/TGI workflow")
message("Results root: ", results_root)
message("Methods: ", paste(selected_methods, collapse = ", "))
message("Comparison: ", ifelse(run_comparison, "completed", "not run"))
