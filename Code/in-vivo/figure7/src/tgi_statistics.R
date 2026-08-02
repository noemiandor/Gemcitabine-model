# Sample-aware ECDF and exact-permutation statistics for Figure 7.

figure7_grid <- function(data, n = 501L) {
  limits <- range(data$pseudotime, na.rm = TRUE)
  if (!all(is.finite(limits)) || diff(limits) == 0) figure7_stop("Pseudotime must have a finite nonzero range")
  seq(limits[[1L]], limits[[2L]], length.out = n)
}

figure7_sample_ecdfs <- function(data, sample_ids, grid) {
  values <- t(vapply(sample_ids, function(id) {
    x <- data$pseudotime[data$sample_id == id]
    if (!length(x)) figure7_stop("No CellCycle pseudotime values for ", id)
    stats::ecdf(x)(grid)
  }, numeric(length(grid))))
  rownames(values) <- sample_ids
  values
}

figure7_shift_metrics <- function(data, samples, config, group_column = NULL) {
  if (is.null(group_column)) {
    group_column <- as.character(config$tgi$matched_control_group)
  }
  if (!identical(group_column, "initial_ploidy")) {
    figure7_stop(
      "Panel 7D/J ECDF references must be matched by injected initial_ploidy; ",
      "run-confounded endpoint-CN-score groups are not permitted"
    )
  }
  if (!group_column %in% names(samples)) {
    figure7_stop("Missing ECDF reference group column: ", group_column)
  }
  tgi_measure <- figure7_tgi_measure(config)
  grid <- figure7_grid(data)
  ids <- samples$sample_id
  matrix <- figure7_sample_ecdfs(data, ids, grid)
  rows <- lapply(ids, function(id) {
    meta <- samples[samples$sample_id == id, , drop = FALSE]
    group <- as.character(meta[[group_column]])
    controls <- samples$sample_id[samples$dose_mg == 0 & as.character(samples[[group_column]]) == group]
    if (!length(controls)) figure7_stop("No untreated ECDF reference for group ", group)
    reference <- colMeans(matrix[controls, , drop = FALSE])
    delta <- matrix[id, ] - reference
    row <- data.frame(
      sample_id = id, initial_ploidy = meta$initial_ploidy, dose = meta$dose,
      dose_mg = meta$dose_mg, etp_group = meta$etp_group,
      sample_mean_endpoint_ploidy = meta$sample_mean_endpoint_ploidy,
      ecdf_rmse = sqrt(mean(delta^2)), ecdf_ks = max(abs(delta)),
      reference_n_samples = length(controls),
      reference_sample_ids = paste(sort(controls), collapse = ";"),
      reference_type = "untreated_equal_sample_mean_ecdf",
      reference_group_column = group_column,
      reference_group_value = group,
      stringsAsFactors = FALSE
    )
    row[[tgi_measure]] <- meta[[tgi_measure]]
    figure7_add_tgi_metadata(row, config)
  })
  do.call(rbind, rows)
}

figure7_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y); x <- x[keep]; y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(c(n = length(x), estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  c(n = length(x), estimate = unname(test$estimate), p_value = test$p.value)
}

figure7_permutation_matrix <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (exists(key, cache, inherits = FALSE)) return(get(key, cache, inherits = FALSE))
    out <- matrix(NA_integer_, nrow = factorial(n), ncol = n); row <- 0L
    visit <- function(prefix, remaining) {
      if (!length(remaining)) { row <<- row + 1L; out[row, ] <<- prefix; return() }
      for (i in seq_along(remaining)) visit(c(prefix, remaining[[i]]), remaining[-i])
    }
    visit(integer(), seq_len(n)); assign(key, out, cache); out
  }
})

figure7_stratified_permutation_matrix <- function(strata) {
  groups <- split(seq_along(strata), as.character(strata))
  local <- lapply(groups, function(index) {
    permutations <- figure7_permutation_matrix(length(index))
    matrix(index[permutations], nrow = nrow(permutations))
  })
  out <- list()
  visit <- function(i, current) {
    if (i > length(local)) { out[[length(out) + 1L]] <<- current; return() }
    index <- groups[[i]]
    for (row in seq_len(nrow(local[[i]]))) {
      next_index <- current; next_index[index] <- local[[i]][row, ]; visit(i + 1L, next_index)
    }
  }
  visit(1L, seq_along(strata)); do.call(rbind, out)
}

figure7_exact_cor <- function(x, y, strata = NULL, method = "pearson") {
  observed <- suppressWarnings(stats::cor(x, y, method = method))
  index <- if (is.null(strata)) figure7_permutation_matrix(length(y)) else figure7_stratified_permutation_matrix(strata)
  permuted <- vapply(seq_len(nrow(index)), function(i) {
    suppressWarnings(stats::cor(x, y[index[i, ]], method = method))
  }, numeric(1L))
  asymptotic <- figure7_safe_cor(x, y, method)
  data.frame(
    n = length(x), estimate = observed, asymptotic_p = asymptotic[["p_value"]],
    permutation_p_two_sided = mean(abs(permuted) >= abs(observed) - 1e-15),
    n_permutations = length(permuted),
    permutation_mode = if (is.null(strata)) "exact_TGI_label_enumeration" else "exact_within_stratum_TGI_label_enumeration",
    stringsAsFactors = FALSE
  )
}

figure7_group_assignments <- function(groups, low, high, strata) {
  blocks <- split(seq_along(groups), as.character(strata))
  choices <- lapply(blocks, function(index) utils::combn(index, sum(groups[index] == low), simplify = FALSE))
  out <- list()
  visit <- function(i, selected) {
    if (i > length(choices)) {
      x <- rep(high, length(groups)); x[unlist(selected, use.names = FALSE)] <- low
      out[[length(out) + 1L]] <<- x; return()
    }
    for (choice in choices[[i]]) visit(i + 1L, c(selected, list(choice)))
  }
  visit(1L, list()); out
}

figure7_dose_adjusted_effect <- function(y, groups, dose, high) {
  design <- stats::model.matrix(~factor(dose) + I(groups == high))
  unname(stats::lm.fit(design, y)$coefficients[["I(groups == high)TRUE"]])
}

figure7_panel_c_test <- function(samples, config) {
  tgi_measure <- figure7_tgi_measure(config)
  x <- samples[samples$dose_mg > 0, , drop = FALSE]
  x <- x[order(x$dose_mg, x$initial_ploidy, x$sample_id), , drop = FALSE]
  outcome <- x[[tgi_measure]]
  groups <- x$initial_ploidy; observed <- figure7_dose_adjusted_effect(outcome, groups, x$dose_mg, "4N")
  assignments <- figure7_group_assignments(groups, "2N", "4N", x$dose_mg)
  effects <- vapply(assignments, function(g) figure7_dose_adjusted_effect(outcome, g, x$dose_mg, "4N"), numeric(1L))
  result <- data.frame(
    sample_set = "treated", group_low = "2N", group_high = "4N",
    n = nrow(x), n_group_low = sum(groups == "2N"), n_group_high = sum(groups == "4N"),
    mean_group_low = mean(outcome[groups == "2N"]),
    mean_group_high = mean(outcome[groups == "4N"]),
    dose_adjusted_difference_high_minus_low = observed,
    permutation_p_two_sided = mean(abs(effects) >= abs(observed) - 1e-15),
    n_permutations = length(effects), permutation_mode = "exact_group_label_enumeration_within_dose",
    comparison_design = "independent_tumors_with_group_labels_permuted_within_dose",
    pairing_status = "not_paired_no_one_to_one_mouse_key",
    stringsAsFactors = FALSE
  )
  figure7_add_tgi_metadata(result, config)
}

figure7_panel_d <- function(shifts, config) {
  tgi_measure <- figure7_tgi_measure(config)
  x <- shifts[shifts$dose_mg > 0, , drop = FALSE]
  x$shift_centered <- x$ecdf_rmse - ave(x$ecdf_rmse, x$dose_mg, FUN = mean)
  x$tgi_centered <- x[[tgi_measure]] - ave(x[[tgi_measure]], x$dose_mg, FUN = mean)
  if (!identical(unique(as.character(x$reference_group_column)), "initial_ploidy")) {
    figure7_stop("Panel 7D/J requires initial-ploidy-matched ECDF references")
  }
  test <- figure7_exact_cor(x$shift_centered, x$tgi_centered, x$dose_mg)
  test$ecdf_reference_group_column <- "initial_ploidy"
  test$centering <- "within_dose"
  test$permutation_strata <- "dose_mg"
  list(data = x, test = test)
}

figure7_endpoint_ploidy_association <- function(data, outcome_column) {
  required <- c(
    "sample_id", "initial_ploidy", "dose", "dose_mg",
    "endpoint_ploidy_file", "sample_mean_endpoint_ploidy",
    "n_endpoint_ploidy_cells", "endpoint_ploidy_source_total_cells",
    "endpoint_ploidy_score_universe_total_cells",
    "endpoint_ploidy_source_file_count", "endpoint_ploidy_source_sha256",
    "endpoint_ploidy_score_policy", "endpoint_ploidy_mapping_policy",
    outcome_column
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    figure7_stop(
      "Panel 7E/K endpoint-ploidy analysis is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  data <- data[order(data$initial_ploidy, data$dose_mg, data$sample_id), , drop = FALSE]
  data$etp_group <- NULL
  if (nrow(data) != 8L || any(!is.finite(data[[outcome_column]])) ||
      any(!is.finite(data$sample_mean_endpoint_ploidy)) ||
      any(!is.finite(data$n_endpoint_ploidy_cells)) ||
      any(data$n_endpoint_ploidy_cells <= 0) ||
      any(data$endpoint_ploidy_source_total_cells != 14125L) ||
      any(data$endpoint_ploidy_score_universe_total_cells != 9832L) ||
      any(data$endpoint_ploidy_source_file_count != 16L) ||
      length(unique(data$endpoint_ploidy_source_sha256)) != 1L ||
      any(!grepl("^[0-9a-f]{64}$", data$endpoint_ploidy_source_sha256)) ||
      length(unique(data$endpoint_ploidy_score_policy)) != 1L ||
      length(unique(data$endpoint_ploidy_mapping_policy)) != 1L) {
    figure7_stop(
      "Panel 7E/K requires complete finite mouse-level endpoint ploidy ",
      "and TGI outcomes"
    )
  }
  test <- figure7_exact_cor(
    data$sample_mean_endpoint_ploidy,
    data[[outcome_column]]
  )
  if (!is.finite(test$estimate) || !is.finite(test$asymptotic_p) ||
      !is.finite(test$permutation_p_two_sided) ||
      test$n_permutations != factorial(8L)) {
    figure7_stop("Panel 7E/K raw endpoint-ploidy association is not estimable")
  }
  test$permutation_strata <- "none"
  test$association_type <- "unadjusted_mouse_level_pearson"
  test$score_variable <- "sample_mean_qc_passed_curated_cbs_cell_ploidy"
  test$score_source_sha256 <- unique(data$endpoint_ploidy_source_sha256)
  test$score_source_n_cells <-
    unique(data$endpoint_ploidy_score_universe_total_cells)
  test$score_inventory_n_cells <- unique(data$endpoint_ploidy_source_total_cells)
  test$score_source_n_files <- unique(data$endpoint_ploidy_source_file_count)
  test$treated_score_n_cells <- sum(data$n_endpoint_ploidy_cells)
  test$score_aggregation_policy <- unique(data$endpoint_ploidy_score_policy)
  test$sample_mapping_policy <- unique(data$endpoint_ploidy_mapping_policy)
  test$score_standardization <- "none"
  test$adjustment_terms <- "none"
  test$outcome_variable <- outcome_column
  test$plot_x <- "sample_mean_endpoint_ploidy"
  test$plot_y <- outcome_column
  list(data = data, test = test)
}

figure7_panel_e <- function(samples, config) {
  tgi_measure <- figure7_tgi_measure(config)
  x <- samples[samples$dose_mg > 0, , drop = FALSE]
  if (nrow(x) != 8L || sum(x$n_endpoint_ploidy_cells) != 5335L) {
    figure7_stop(
      "Panel 7E/K requires the exact 5,335 QC-passed curated CBS cells ",
      "from the eight treated tumors"
    )
  }
  figure7_endpoint_ploidy_association(x, tgi_measure)
}

figure7_mean_ecdf <- function(data, ids, grid) {
  colMeans(figure7_sample_ecdfs(data, ids, grid))
}

figure7_density_localization_contract <- function(config) {
  local <- config$density_localization
  required <- c(
    "analysis_id", "cell_universe", "expected_cells", "expected_samples",
    "expected_vehicle_samples", "expected_treated_samples",
    "expected_samples_per_treatment_within_origin", "estimator",
    "bandwidth_method", "bandwidth", "grid_start", "grid_end",
    "grid_points", "mouse_weighting", "treatment_contrast",
    "permutation_strata", "expected_permutations", "pointwise_alpha",
    "simultaneous_alpha", "simultaneous_method",
    "expected_pointwise_start", "expected_pointwise_end",
    "expected_simultaneous_start", "expected_simultaneous_end",
    "expected_simultaneous_critical",
    "expected_raw_excess_peak_pseudotime",
    "expected_raw_excess_peak_density_difference",
    "expected_max_abs_t_pseudotime", "expected_global_max_abs_t",
    "expected_global_max_t_p_two_sided"
  )
  if (is.null(local) || !all(required %in% names(local))) {
    figure7_stop(
      "Source Figure 7B/SI4I density-localization config is incomplete: ",
      paste(setdiff(required, names(local)), collapse = ", ")
    )
  }
  expected_strings <- c(
    analysis_id = "equal_mouse_kde_exact_origin_stratified_max_t_v1",
    cell_universe = "reviewed_qc_retained_cellcycle_2881",
    estimator = "stats_density_gaussian",
    bandwidth_method = "pooled_label_invariant_bw_nrd0",
    mouse_weighting = "equal",
    treatment_contrast = "treated_minus_vehicle",
    permutation_strata = "initial_ploidy",
    simultaneous_method = "studentized_max_abs_t"
  )
  observed_strings <- vapply(
    names(expected_strings),
    function(name) as.character(local[[name]]),
    character(1L)
  )
  if (!identical(unname(observed_strings), unname(expected_strings))) {
    figure7_stop(
      "Source Figure 7B/SI4I density-localization method must remain the reviewed ",
      "equal-mouse, injected-origin-stratified contract"
    )
  }
  numeric_fields <- c(
    "expected_cells", "expected_samples", "expected_vehicle_samples",
    "expected_treated_samples", "expected_samples_per_treatment_within_origin",
    "bandwidth", "grid_start", "grid_end", "grid_points",
    "expected_permutations", "pointwise_alpha", "simultaneous_alpha",
    "expected_pointwise_start", "expected_pointwise_end",
    "expected_simultaneous_start", "expected_simultaneous_end",
    "expected_simultaneous_critical",
    "expected_raw_excess_peak_pseudotime",
    "expected_raw_excess_peak_density_difference",
    "expected_max_abs_t_pseudotime", "expected_global_max_abs_t",
    "expected_global_max_t_p_two_sided"
  )
  values <- vapply(
    numeric_fields,
    function(name) suppressWarnings(as.numeric(local[[name]])),
    numeric(1L)
  )
  if (any(!is.finite(values)) || values[["expected_cells"]] != 2881 ||
      values[["expected_samples"]] != 16 ||
      values[["expected_vehicle_samples"]] != 8 ||
      values[["expected_treated_samples"]] != 8 ||
      values[["expected_samples_per_treatment_within_origin"]] != 4 ||
      values[["grid_start"]] != 0 || values[["grid_end"]] != 1 ||
      values[["grid_points"]] != 501 ||
      values[["expected_permutations"]] != 4900 ||
      values[["pointwise_alpha"]] != 0.05 ||
      values[["simultaneous_alpha"]] != 0.05) {
    figure7_stop("Source Figure 7B/SI4I density-localization numeric contract is invalid")
  }
  list(raw = local, numeric = values)
}

figure7_supported_intervals <- function(grid, supported, support_type, alpha) {
  if (length(grid) != length(supported) || anyNA(supported)) {
    figure7_stop("Density-localization support mask is invalid")
  }
  runs <- rle(as.logical(supported))
  end_index <- cumsum(runs$lengths)
  start_index <- end_index - runs$lengths + 1L
  keep <- which(runs$values)
  if (!length(keep)) {
    return(data.frame(
      support_type = character(), start = numeric(), end = numeric(),
      width = numeric(), alpha = numeric(), stringsAsFactors = FALSE
    ))
  }
  data.frame(
    support_type = rep(support_type, length(keep)),
    start = grid[start_index[keep]],
    end = grid[end_index[keep]],
    width = grid[end_index[keep]] - grid[start_index[keep]],
    alpha = rep(alpha, length(keep)),
    stringsAsFactors = FALSE
  )
}

figure7_density_localization <- function(data, samples, config) {
  contract <- figure7_density_localization_contract(config)
  local <- contract$raw
  values <- contract$numeric
  required_data <- c("sample_id", "pseudotime")
  required_samples <- c("sample_id", "initial_ploidy", "dose_mg")
  if (!all(required_data %in% names(data)) ||
      !all(required_samples %in% names(samples))) {
    figure7_stop("Source Figure 7B/SI4I density localization is missing required inputs")
  }
  samples <- samples[order(samples$sample_id), , drop = FALSE]
  rownames(samples) <- NULL
  if (anyDuplicated(samples$sample_id) ||
      nrow(data) != as.integer(values[["expected_cells"]]) ||
      nrow(samples) != as.integer(values[["expected_samples"]]) ||
      !setequal(unique(as.character(data$sample_id)), samples$sample_id) ||
      any(!is.finite(data$pseudotime)) ||
      any(data$pseudotime < values[["grid_start"]] |
          data$pseudotime > values[["grid_end"]])) {
    figure7_stop(
      "Source Figure 7B/SI4I density localization requires the exact 2,881-cell, ",
      "16-mouse reviewed QC universe"
    )
  }
  treatment <- ifelse(samples$dose_mg == 0, "vehicle", "treated")
  if (sum(treatment == "vehicle") != values[["expected_vehicle_samples"]] ||
      sum(treatment == "treated") != values[["expected_treated_samples"]]) {
    figure7_stop("Source Figure 7B/SI4I density-localization treatment balance is invalid")
  }
  balance <- table(samples$initial_ploidy, treatment)
  if (!identical(sort(unique(as.character(samples$initial_ploidy))), c("2N", "4N")) ||
      !all(balance == values[["expected_samples_per_treatment_within_origin"]])) {
    figure7_stop(
      "Source Figure 7B/SI4I density localization requires four vehicle and four ",
      "treated mice within each injected-origin stratum"
    )
  }
  bandwidth <- stats::bw.nrd0(data$pseudotime)
  if (!isTRUE(all.equal(
    bandwidth,
    values[["bandwidth"]],
    tolerance = 1e-13,
    check.attributes = FALSE
  ))) {
    figure7_stop(
      "Pooled label-invariant density bandwidth disagrees with the reviewed ",
      "source Figure 7B/SI4I contract"
    )
  }
  grid <- seq(
    values[["grid_start"]], values[["grid_end"]],
    length.out = as.integer(values[["grid_points"]])
  )
  density_matrix <- t(vapply(samples$sample_id, function(id) {
    x <- data$pseudotime[as.character(data$sample_id) == id]
    if (!length(x)) figure7_stop("No pseudotime values for density sample ", id)
    stats::density(
      x,
      bw = bandwidth,
      kernel = "gaussian",
      from = values[["grid_start"]],
      to = values[["grid_end"]],
      n = as.integer(values[["grid_points"]]),
      na.rm = TRUE
    )$y
  }, numeric(length(grid))))
  rownames(density_matrix) <- samples$sample_id
  vehicle_mean <- colMeans(density_matrix[treatment == "vehicle", , drop = FALSE])
  treated_mean <- colMeans(density_matrix[treatment == "treated", , drop = FALSE])
  observed <- treated_mean - vehicle_mean

  assignments <- figure7_group_assignments(
    treatment, "vehicle", "treated", samples$initial_ploidy
  )
  if (length(assignments) != as.integer(values[["expected_permutations"]])) {
    figure7_stop("Source Figure 7B/SI4I density-localization permutation count is invalid")
  }
  permuted <- vapply(assignments, function(group) {
    colMeans(density_matrix[group == "treated", , drop = FALSE]) -
      colMeans(density_matrix[group == "vehicle", , drop = FALSE])
  }, numeric(length(grid)))
  permutation_sd <- apply(permuted, 1L, stats::sd)
  if (any(!is.finite(permutation_sd)) || any(permutation_sd <= 0)) {
    figure7_stop("Source Figure 7B/SI4I density-localization permutation variance is invalid")
  }
  observed_studentized <- observed / permutation_sd
  permutation_studentized <- sweep(permuted, 1L, permutation_sd, "/")
  max_abs_permutation <- apply(abs(permutation_studentized), 2L, max)
  simultaneous_critical <- as.numeric(stats::quantile(
    max_abs_permutation,
    probs = 1 - values[["simultaneous_alpha"]],
    type = 1,
    names = FALSE
  ))
  tolerance <- 1e-15
  pointwise_p <- rowMeans(
    abs(permuted) >= abs(observed) - tolerance
  )
  max_t_adjusted_p <- vapply(abs(observed_studentized), function(value) {
    mean(max_abs_permutation >= value - tolerance)
  }, numeric(1L))
  pointwise_supported <- observed > 0 &
    pointwise_p <= values[["pointwise_alpha"]]
  simultaneous_supported <- observed_studentized > simultaneous_critical
  pointwise_intervals <- figure7_supported_intervals(
    grid, pointwise_supported, "positive_pointwise_two_sided",
    values[["pointwise_alpha"]]
  )
  simultaneous_intervals <- figure7_supported_intervals(
    grid, simultaneous_supported, "positive_simultaneous_max_abs_t",
    values[["simultaneous_alpha"]]
  )
  intervals <- rbind(pointwise_intervals, simultaneous_intervals)
  expected_intervals <- data.frame(
    support_type = c(
      "positive_pointwise_two_sided", "positive_simultaneous_max_abs_t"
    ),
    start = c(
      values[["expected_pointwise_start"]],
      values[["expected_simultaneous_start"]]
    ),
    end = c(
      values[["expected_pointwise_end"]],
      values[["expected_simultaneous_end"]]
    ),
    stringsAsFactors = FALSE
  )
  observed_intervals <- intervals[, c("support_type", "start", "end"), drop = FALSE]
  rownames(observed_intervals) <- NULL
  if (!identical(
        observed_intervals$support_type,
        expected_intervals$support_type
      ) ||
      !isTRUE(all.equal(
        observed_intervals$start,
        expected_intervals$start,
        tolerance = 1e-12,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        observed_intervals$end,
        expected_intervals$end,
        tolerance = 1e-12,
        check.attributes = FALSE
      ))) {
    figure7_stop(
      "Source Figure 7B/SI4I density-localization support intervals disagree with ",
      "the reviewed 2,881-cell analysis"
    )
  }
  raw_excess_peak_index <- which.max(observed)
  max_abs_t_index <- which.max(abs(observed_studentized))
  global_max_t <- max(abs(observed_studentized))
  global_p <- mean(max_abs_permutation >= global_max_t - tolerance)
  reviewed_values <- c(
    simultaneous_critical = simultaneous_critical,
    raw_excess_peak_pseudotime = grid[[raw_excess_peak_index]],
    raw_excess_peak_density_difference = observed[[raw_excess_peak_index]],
    max_abs_t_pseudotime = grid[[max_abs_t_index]],
    global_max_abs_t = global_max_t,
    global_max_t_p_two_sided = global_p
  )
  expected_reviewed_values <- c(
    simultaneous_critical = values[["expected_simultaneous_critical"]],
    raw_excess_peak_pseudotime = values[["expected_raw_excess_peak_pseudotime"]],
    raw_excess_peak_density_difference = values[["expected_raw_excess_peak_density_difference"]],
    max_abs_t_pseudotime = values[["expected_max_abs_t_pseudotime"]],
    global_max_abs_t = values[["expected_global_max_abs_t"]],
    global_max_t_p_two_sided = values[["expected_global_max_t_p_two_sided"]]
  )
  if (!isTRUE(all.equal(
    reviewed_values,
    expected_reviewed_values,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    figure7_stop(
      "Source Figure 7B/SI4I density-localization peak or global inference disagrees ",
      "with the reviewed 2,881-cell analysis"
    )
  }
  grid_data <- data.frame(
    pseudotime = grid,
    vehicle_mean_density = vehicle_mean,
    treated_mean_density = treated_mean,
    treated_minus_vehicle_density = observed,
    permutation_sd = permutation_sd,
    observed_studentized = observed_studentized,
    pointwise_p_two_sided = pointwise_p,
    max_t_adjusted_p_two_sided = max_t_adjusted_p,
    simultaneous_critical = rep(simultaneous_critical, length(grid)),
    simultaneous_lower_envelope = -simultaneous_critical * permutation_sd,
    simultaneous_upper_envelope = simultaneous_critical * permutation_sd,
    pointwise_positive_supported = pointwise_supported,
    simultaneous_positive_supported = simultaneous_supported,
    frozen_state_interval = grid >= config$intervals$primary_accumulated_state$start &
      grid <= config$intervals$primary_accumulated_state$end,
    stringsAsFactors = FALSE
  )
  test <- data.frame(
    analysis_id = as.character(local$analysis_id),
    cell_universe = as.character(local$cell_universe),
    n_cells = nrow(data),
    n_samples = nrow(samples),
    n_vehicle_samples = sum(treatment == "vehicle"),
    n_treated_samples = sum(treatment == "treated"),
    sample_weighting = "equal_mouse",
    density_estimator = "stats::density Gaussian kernel",
    bandwidth_method = "pooled label-invariant stats::bw.nrd0",
    bandwidth = bandwidth,
    grid_start = min(grid),
    grid_end = max(grid),
    grid_points = length(grid),
    contrast = "treated_minus_vehicle",
    permutation_strata = "initial_ploidy",
    n_permutations = length(assignments),
    pointwise_test = "two-sided exact permutation at each grid point",
    pointwise_alpha = values[["pointwise_alpha"]],
    pointwise_start = pointwise_intervals$start,
    pointwise_end = pointwise_intervals$end,
    simultaneous_test = "studentized max-absolute-T exact permutation null envelope",
    simultaneous_alpha = values[["simultaneous_alpha"]],
    simultaneous_critical = simultaneous_critical,
    simultaneous_start = simultaneous_intervals$start,
    simultaneous_end = simultaneous_intervals$end,
    raw_excess_peak_pseudotime = grid[[raw_excess_peak_index]],
    raw_excess_peak_density_difference = observed[[raw_excess_peak_index]],
    max_abs_t_pseudotime = grid[[max_abs_t_index]],
    max_abs_t_observed_statistic = observed_studentized[[max_abs_t_index]],
    global_max_abs_t = global_max_t,
    global_max_t_p_two_sided = global_p,
    stringsAsFactors = FALSE
  )
  list(grid = grid_data, intervals = intervals, test = test)
}

figure7_panel_b <- function(data, samples, config) {
  grid <- figure7_grid(data)
  definitions <- list(
    list(id = 1L, label = "1. 0 vs treated", subset = rep(TRUE, nrow(samples))),
    list(id = 8L, label = "8. 4N: 0 vs treated", subset = samples$initial_ploidy == "4N"),
    list(id = 9L, label = "9. 2N: 0 vs treated", subset = samples$initial_ploidy == "2N")
  )
  curve_rows <- list(); test_rows <- list()
  for (definition in definitions) {
    meta <- samples[definition$subset, , drop = FALSE]
    groups <- ifelse(meta$dose_mg == 0, "0mg/kg", "treated")
    ids_a <- meta$sample_id[groups == "0mg/kg"]; ids_b <- meta$sample_id[groups == "treated"]
    delta <- figure7_mean_ecdf(data, ids_a, grid) - figure7_mean_ecdf(data, ids_b, grid)
    assignments <- figure7_group_assignments(groups, "0mg/kg", "treated", rep("all", nrow(meta)))
    ecdf_matrix <- figure7_sample_ecdfs(data, meta$sample_id, grid)
    permuted <- vapply(assignments, function(g) {
      d <- colMeans(ecdf_matrix[g == "0mg/kg", , drop = FALSE]) - colMeans(ecdf_matrix[g == "treated", , drop = FALSE])
      sqrt(mean(d^2))
    }, numeric(1L))
    observed <- sqrt(mean(delta^2))
    test_rows[[length(test_rows) + 1L]] <- data.frame(
      comparison_id = definition$id, panel = definition$label, group_a = "0mg/kg", group_b = "treated",
      n_group_a = length(ids_a), n_group_b = length(ids_b), observed_ecdf_rmse = observed,
      observed_ecdf_ks = max(abs(delta)), p_ecdf_rmse = mean(permuted >= observed - 1e-15),
      n_permutations = length(permuted), permutation_mode = "exact_label_enumeration", stringsAsFactors = FALSE
    )
    for (group in c("0mg/kg", "treated")) {
      ids <- meta$sample_id[groups == group]
      curve_rows[[length(curve_rows) + 1L]] <- data.frame(
        comparison_id = definition$id, panel = definition$label, pseudotime = grid,
        mean_ecdf = figure7_mean_ecdf(data, ids, grid), curve_label = group,
        color_group = group,
        line_group = if (definition$id == 8L) "4N" else if (definition$id == 9L) "2N" else "All",
        n_samples = length(ids), stringsAsFactors = FALSE
      )
    }
  }
  tests <- do.call(rbind, test_rows)
  tests$annotation <- sprintf("%s n=%d; %s n=%d\nRMSE=%.4f, P=%.3g\nKS=%.4f",
                              tests$group_a, tests$n_group_a, tests$group_b, tests$n_group_b,
                              tests$observed_ecdf_rmse, tests$p_ecdf_rmse, tests$observed_ecdf_ks)
  localization <- figure7_density_localization(data, samples, config)
  list(
    data = do.call(rbind, curve_rows),
    tests = tests,
    localization_grid = localization$grid,
    localization_intervals = localization$intervals,
    localization_test = localization$test
  )
}
