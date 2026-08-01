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

figure7_within_group_z <- function(values, groups, value_label) {
  values <- as.numeric(values)
  groups <- as.character(groups)
  if (length(values) != length(groups) || !length(values) ||
      any(!is.finite(values)) || any(is.na(groups)) || any(!nzchar(groups))) {
    figure7_stop("Cannot standardize ", value_label, ": values and groups must be complete")
  }
  standardized <- rep(NA_real_, length(values))
  for (group in unique(groups)) {
    index <- which(groups == group)
    spread <- stats::sd(values[index])
    if (length(index) < 2L || !is.finite(spread) || spread <= 0) {
      figure7_stop(
        "Cannot standardize ", value_label, " within injected origin ",
        group, ": at least two distinct finite values are required"
      )
    }
    standardized[index] <- (values[index] - mean(values[index])) / spread
  }
  standardized
}

figure7_adjusted_cn_association <- function(data, outcome_column) {
  required <- c(
    "sample_id", "initial_ploidy", "dose", "dose_mg",
    "endpoint_ploidy_file", "sample_mean_endpoint_ploidy",
    "n_endpoint_ploidy_cells", "endpoint_ploidy_source_total_cells",
    "endpoint_ploidy_source_file_count", "endpoint_ploidy_source_sha256",
    "endpoint_ploidy_score_policy", "endpoint_ploidy_mapping_policy",
    outcome_column
  )
  missing <- setdiff(required, names(data))
  if (length(missing)) {
    figure7_stop(
      "Panel 7E/K adjusted CN-score analysis is missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
  data <- data[order(data$initial_ploidy, data$dose_mg, data$sample_id), , drop = FALSE]
  data$etp_group <- NULL
  if (nrow(data) < 5L || any(!is.finite(data[[outcome_column]])) ||
      any(!is.finite(data$sample_mean_endpoint_ploidy)) ||
      any(!is.finite(data$n_endpoint_ploidy_cells)) ||
      any(data$n_endpoint_ploidy_cells <= 0) ||
      any(data$endpoint_ploidy_source_total_cells != 14125L) ||
      any(data$endpoint_ploidy_source_file_count != 16L) ||
      length(unique(data$endpoint_ploidy_source_sha256)) != 1L ||
      any(!grepl("^[0-9a-f]{64}$", data$endpoint_ploidy_source_sha256)) ||
      length(unique(data$endpoint_ploidy_score_policy)) != 1L ||
      length(unique(data$endpoint_ploidy_mapping_policy)) != 1L) {
    figure7_stop("Panel 7E/K requires complete finite tumor-level CN scores and TGI outcomes")
  }
  data$terminal_postprocessed_cn_score <- data$sample_mean_endpoint_ploidy
  data$terminal_cn_score_within_origin_z <- figure7_within_group_z(
    data$terminal_postprocessed_cn_score,
    data$initial_ploidy,
    "terminal postprocessed CN score"
  )
  data$permutation_stratum <- paste(data$initial_ploidy, data$dose_mg, sep = "|")

  initial_ploidy_factor <- factor(data$initial_ploidy)
  dose_factor <- factor(data$dose_mg)
  nuisance <- stats::model.matrix(~ initial_ploidy_factor + dose_factor)
  design <- stats::model.matrix(
    ~ terminal_cn_score_within_origin_z + initial_ploidy_factor + dose_factor,
    data = data
  )
  if (qr(nuisance)$rank != ncol(nuisance) || qr(design)$rank != ncol(design)) {
    figure7_stop("Panel 7E/K origin- and dose-adjusted design is rank deficient")
  }

  outcome <- data[[outcome_column]]
  fit <- stats::lm.fit(design, outcome)
  coefficient_name <- "terminal_cn_score_within_origin_z"
  effect <- unname(fit$coefficients[[coefficient_name]])
  data$terminal_cn_score_nuisance_residual <- unname(
    stats::lm.fit(nuisance, data$terminal_cn_score_within_origin_z)$residuals
  )
  data$tgi_origin_dose_residual <- unname(
    stats::lm.fit(nuisance, outcome)$residuals
  )
  partial_correlation <- suppressWarnings(stats::cor(
    data$terminal_cn_score_nuisance_residual,
    data$tgi_origin_dose_residual
  ))
  if (!is.finite(effect) || !is.finite(partial_correlation)) {
    figure7_stop("Panel 7E/K adjusted CN-score association is not estimable")
  }

  permutation_index <- figure7_stratified_permutation_matrix(data$permutation_stratum)
  permuted_effects <- vapply(seq_len(nrow(permutation_index)), function(i) {
    unname(stats::lm.fit(design, outcome[permutation_index[i, ]])$coefficients[[coefficient_name]])
  }, numeric(1L))
  model <- stats::lm(
    outcome ~ terminal_cn_score_within_origin_z + initial_ploidy_factor + dose_factor,
    data = data
  )
  asymptotic_p <- unname(summary(model)$coefficients[coefficient_name, "Pr(>|t|)"])
  test <- data.frame(
    n = nrow(data),
    estimate = partial_correlation,
    partial_correlation = partial_correlation,
    effect_per_within_origin_sd = effect,
    asymptotic_p = asymptotic_p,
    permutation_p_two_sided = mean(
      abs(permuted_effects) >= abs(effect) - 1e-15
    ),
    n_permutations = nrow(permutation_index),
    permutation_mode = "exact_TGI_label_enumeration_within_initial_ploidy_x_dose",
    permutation_strata = "initial_ploidy:dose_mg",
    score_variable =
      "sample_mean_all_canonical_cbs_cell_ploidy",
    score_source_sha256 = unique(data$endpoint_ploidy_source_sha256),
    score_source_n_cells = unique(data$endpoint_ploidy_source_total_cells),
    score_source_n_files = unique(data$endpoint_ploidy_source_file_count),
    treated_score_n_cells = sum(data$n_endpoint_ploidy_cells),
    score_aggregation_policy = unique(data$endpoint_ploidy_score_policy),
    sample_mapping_policy = unique(data$endpoint_ploidy_mapping_policy),
    score_standardization = "z_score_within_initial_ploidy",
    adjustment_terms = "initial_ploidy+dose_mg",
    outcome_variable = outcome_column,
    plot_x = "terminal_cn_score_nuisance_residual",
    plot_y = "tgi_origin_dose_residual",
    stringsAsFactors = FALSE
  )
  list(data = data, test = test)
}

figure7_panel_e <- function(samples, config) {
  tgi_measure <- figure7_tgi_measure(config)
  x <- samples[samples$dose_mg > 0, , drop = FALSE]
  if (nrow(x) != 8L || sum(x$n_endpoint_ploidy_cells) != 7623L) {
    figure7_stop(
      "Panel 7E/K requires all 7,623 canonical CBS cells from the exact ",
      "eight treated tumors"
    )
  }
  figure7_adjusted_cn_association(x, tgi_measure)
}

figure7_mean_ecdf <- function(data, ids, grid) {
  colMeans(figure7_sample_ecdfs(data, ids, grid))
}

figure7_panel_b <- function(data, samples) {
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
  list(data = do.call(rbind, curve_rows), tests = tests)
}
