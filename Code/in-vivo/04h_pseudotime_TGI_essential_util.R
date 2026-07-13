essential_parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      out[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (grepl("^--", token)) {
      key <- sub("^--", "", token)
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

essential_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

essential_bool <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

essential_safe_numeric <- function(x) suppressWarnings(as.numeric(x))

essential_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

essential_write_csv <- function(x, path) {
  essential_ensure_dir(dirname(path))
  write.csv(x, path, row.names = FALSE, na = "NA")
  invisible(path)
}

essential_save_pdf <- function(plot, path, width, height) {
  essential_ensure_dir(dirname(path))
  suppressMessages(ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    device = grDevices::cairo_pdf
  ))
  invisible(path)
}

essential_format_number <- function(x, digits = 3L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "fg"), "NA")
}

essential_format_p <- function(x) {
  ifelse(
    is.finite(x),
    ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "fg")),
    "NA"
  )
}

essential_endpoint_tgi_columns <- function(names_vector) {
  columns <- grep("^TGI_percent_Day_[0-9]+$", names_vector, value = TRUE)
  days <- essential_safe_numeric(sub("^TGI_percent_Day_", "", columns))
  columns[order(days)]
}

essential_endpoint_day <- function(column) {
  essential_safe_numeric(sub("^TGI_percent_Day_", "", column))
}

essential_method_specs <- function() {
  list(
    initial_ploidy = list(
      method = "initial_ploidy",
      type = "initial_ploidy",
      threshold = NA_real_,
      group_levels = c("2N", "4N"),
      group_label = "Initial ploidy",
      group_column = "initial_ploidy",
      group_numeric_column = "initial_ploidy_numeric"
    ),
    ETP_fixed_threshold_2_25 = list(
      method = "ETP_fixed_threshold_2_25",
      type = "ETP",
      threshold = 2.25,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    ),
    ETP_boundary_stress_threshold_2_375 = list(
      method = "ETP_boundary_stress_threshold_2_375",
      type = "ETP",
      threshold = 2.375,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    ),
    ETP_reference_balanced_threshold_2_24 = list(
      method = "ETP_reference_balanced_threshold_2_24",
      type = "ETP",
      threshold = 2.24,
      group_levels = c("ETP-lower", "ETP-higher"),
      group_label = "EndTimePoint ploidy",
      group_column = "end_timepoint_ploidy_group",
      group_numeric_column = "end_timepoint_ploidy_group_numeric"
    )
  )
}

essential_required_input_columns <- function() {
  c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "TGI_percent_auc"
  )
}

essential_read_input <- function(path, compartment) {
  if (!file.exists(path)) stop("Missing input CSV: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(essential_required_input_columns(), names(data))
  if (length(missing) > 0L) {
    stop("Input is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  numeric_columns <- unique(c(
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "TGI_percent_auc",
    essential_endpoint_tgi_columns(names(data))
  ))
  for (column in intersect(numeric_columns, names(data))) {
    data[[column]] <- essential_safe_numeric(data[[column]])
  }
  data$compartment <- compartment
  data
}

essential_unique_sample_value <- function(data, column, sample_id, numeric = FALSE) {
  values <- data[data$sample_id == sample_id, column]
  if (numeric) {
    values <- essential_safe_numeric(values)
    values <- unique(values[is.finite(values)])
  } else {
    values <- unique(as.character(values[!is.na(values)]))
  }
  if (length(values) > 1L) {
    stop("Inconsistent sample-level column ", column, " for sample ", sample_id, call. = FALSE)
  }
  if (length(values) == 0L) {
    if (numeric) NA_real_ else NA_character_
  } else {
    values[[1L]]
  }
}

essential_derive_assignments <- function(cellcycle, noncellcycle) {
  needed <- c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "cell_ploidy"
  )
  union_cells <- rbind(
    cellcycle[, needed, drop = FALSE],
    noncellcycle[, needed, drop = FALSE]
  )
  union_cells$cell_ploidy <- essential_safe_numeric(union_cells$cell_ploidy)
  union_cells$gemcitabine_dose_mg_per_kg <- essential_safe_numeric(union_cells$gemcitabine_dose_mg_per_kg)
  union_cells <- union_cells[is.finite(union_cells$cell_ploidy), , drop = FALSE]
  if (nrow(union_cells) == 0L) stop("No finite cell ploidy values were found", call. = FALSE)

  duplicated_ids <- unique(union_cells$cell_id[duplicated(union_cells$cell_id)])
  for (cell_id in duplicated_ids) {
    local <- union_cells[union_cells$cell_id == cell_id, , drop = FALSE]
    if (length(unique(local$sample_id)) != 1L || length(unique(local$cell_ploidy)) != 1L) {
      stop("Conflicting duplicated cell_id: ", cell_id, call. = FALSE)
    }
  }
  union_cells <- union_cells[!duplicated(union_cells$cell_id), , drop = FALSE]

  sample_ids <- sort(unique(union_cells$sample_id))
  rows <- lapply(sample_ids, function(sample_id) {
    local <- union_cells[union_cells$sample_id == sample_id, , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      initial_ploidy = essential_unique_sample_value(local, "initial_ploidy", sample_id),
      dose = essential_unique_sample_value(local, "gemcitabine_dose", sample_id),
      dose_mg = essential_unique_sample_value(local, "gemcitabine_dose_mg_per_kg", sample_id, numeric = TRUE),
      n_endpoint_ploidy_cells = nrow(local),
      sample_mean_endpoint_ploidy = mean(local$cell_ploidy),
      sample_median_endpoint_ploidy = median(local$cell_ploidy),
      sample_max_endpoint_ploidy = max(local$cell_ploidy),
      stringsAsFactors = FALSE
    )
  })
  assignments <- do.call(rbind, rows)
  unexpected <- setdiff(unique(assignments$initial_ploidy), c("2N", "4N"))
  if (length(unexpected) > 0L) {
    stop("Unexpected initial ploidy labels: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }
  assignments
}

essential_prepare_cellcycle <- function(cellcycle, assignments, spec) {
  data <- cellcycle[is.finite(cellcycle$pseudotime), , drop = FALSE]
  assignment_index <- match(data$sample_id, assignments$sample_id)
  if (anyNA(assignment_index)) {
    stop("Missing assignment for CellCycle samples", call. = FALSE)
  }
  joined <- assignments[assignment_index, , drop = FALSE]
  if (any(as.character(data$initial_ploidy) != as.character(joined$initial_ploidy))) {
    stop("Initial ploidy metadata disagree between inputs", call. = FALSE)
  }
  data$dose <- as.character(data$gemcitabine_dose)
  data$dose_mg <- essential_safe_numeric(data$gemcitabine_dose_mg_per_kg)
  data$sample_mean_endpoint_ploidy <- joined$sample_mean_endpoint_ploidy
  data$sample_median_endpoint_ploidy <- joined$sample_median_endpoint_ploidy
  data$sample_max_endpoint_ploidy <- joined$sample_max_endpoint_ploidy
  data$initial_ploidy <- factor(as.character(data$initial_ploidy), levels = c("2N", "4N"))
  data$initial_ploidy_numeric <- ifelse(as.character(data$initial_ploidy) == "4N", 1, 0)
  if (identical(spec$type, "ETP")) {
    etp_group <- ifelse(
      data$sample_mean_endpoint_ploidy > spec$threshold,
      "ETP-higher",
      "ETP-lower"
    )
    data$end_timepoint_ploidy_group <- factor(etp_group, levels = spec$group_levels)
    data$end_timepoint_ploidy_group_numeric <- ifelse(etp_group == "ETP-higher", 1, 0)
    data$analysis_group <- data$end_timepoint_ploidy_group
    data$analysis_group_numeric <- data$end_timepoint_ploidy_group_numeric
  } else {
    data$analysis_group <- factor(as.character(data$initial_ploidy), levels = spec$group_levels)
    data$analysis_group_numeric <- data$initial_ploidy_numeric
  }
  data
}

essential_build_sample_meta <- function(data) {
  sample_ids <- sort(unique(data$sample_id))
  tgi_columns <- unique(c("TGI_percent_auc", essential_endpoint_tgi_columns(names(data))))
  rows <- lapply(sample_ids, function(sample_id) {
    local <- data[data$sample_id == sample_id, , drop = FALSE]
    base <- data.frame(
      sample_id = sample_id,
      initial_ploidy = as.character(local$initial_ploidy[[1L]]),
      initial_ploidy_numeric = local$initial_ploidy_numeric[[1L]],
      analysis_group = as.character(local$analysis_group[[1L]]),
      analysis_group_numeric = local$analysis_group_numeric[[1L]],
      dose = as.character(local$dose[[1L]]),
      dose_mg = local$dose_mg[[1L]],
      sample_mean_endpoint_ploidy = local$sample_mean_endpoint_ploidy[[1L]],
      n_cells = nrow(local),
      mean_pseudotime = mean(local$pseudotime),
      median_pseudotime = median(local$pseudotime),
      mean_cell_ploidy = mean(local$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = median(local$cell_ploidy, na.rm = TRUE),
      p90_cell_ploidy = as.numeric(stats::quantile(local$cell_ploidy, 0.90, na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
    for (column in tgi_columns) {
      base[[column]] <- essential_unique_sample_value(local, column, sample_id, numeric = TRUE)
    }
    base
  })
  meta <- do.call(rbind, rows)
  meta <- meta[order(meta$dose_mg, meta$analysis_group_numeric, meta$sample_id), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  meta
}

essential_make_grid <- function(data, n = 501L) {
  range_values <- range(data$pseudotime, na.rm = TRUE)
  if (!all(is.finite(range_values)) || diff(range_values) == 0) {
    return(sort(unique(data$pseudotime)))
  }
  seq(range_values[[1L]], range_values[[2L]], length.out = n)
}

essential_sample_ecdf_matrix <- function(data, sample_ids, grid) {
  matrix_values <- t(vapply(sample_ids, function(sample_id) {
    values <- data$pseudotime[data$sample_id == sample_id]
    if (length(values) == 0L) rep(NA_real_, length(grid)) else stats::ecdf(values)(grid)
  }, numeric(length(grid))))
  rownames(matrix_values) <- sample_ids
  matrix_values
}

essential_calculate_shift_metrics <- function(data, sample_meta, pooled_reference = FALSE) {
  grid <- essential_make_grid(data)
  sample_ids <- rownames(sample_meta)
  ecdf_matrix <- essential_sample_ecdf_matrix(data, sample_ids, grid)
  rows <- lapply(sample_ids, function(sample_id) {
    group <- sample_meta[sample_id, "analysis_group"]
    dose_mg <- sample_meta[sample_id, "dose_mg"]
    controls <- sample_meta$sample_id[
      sample_meta$analysis_group == group & sample_meta$dose_mg == 0
    ]
    reference_samples <- if (dose_mg == 0) setdiff(controls, sample_id) else controls
    if (length(reference_samples) == 0L) {
      delta <- rep(NA_real_, length(grid))
      reference_mean <- NA_real_
      reference_cells <- 0L
    } else if (!pooled_reference) {
      reference_ecdf <- colMeans(ecdf_matrix[reference_samples, , drop = FALSE], na.rm = TRUE)
      delta <- ecdf_matrix[sample_id, ] - reference_ecdf
      reference_mean <- mean(sample_meta[reference_samples, "mean_pseudotime"])
      reference_cells <- sum(sample_meta[reference_samples, "n_cells"])
    } else {
      sample_values <- data$pseudotime[data$sample_id == sample_id]
      reference_values <- data$pseudotime[data$sample_id %in% reference_samples]
      delta <- stats::ecdf(sample_values)(grid) - stats::ecdf(reference_values)(grid)
      reference_mean <- mean(reference_values)
      reference_cells <- length(reference_values)
    }
    finite_delta <- delta[is.finite(delta)]
    row <- data.frame(
      sample_id = sample_id,
      reference_type = if (pooled_reference) "sensitivity_pooled_cell_reference" else "primary_equal_sample_reference",
      initial_ploidy = sample_meta[sample_id, "initial_ploidy"],
      initial_ploidy_numeric = sample_meta[sample_id, "initial_ploidy_numeric"],
      analysis_group = group,
      analysis_group_numeric = sample_meta[sample_id, "analysis_group_numeric"],
      dose = sample_meta[sample_id, "dose"],
      dose_mg = dose_mg,
      TGI_percent_auc = sample_meta[sample_id, "TGI_percent_auc"],
      cell_count = sample_meta[sample_id, "n_cells"],
      sample_mean_endpoint_ploidy = sample_meta[sample_id, "sample_mean_endpoint_ploidy"],
      mean_cell_ploidy = sample_meta[sample_id, "mean_cell_ploidy"],
      median_cell_ploidy = sample_meta[sample_id, "median_cell_ploidy"],
      p90_cell_ploidy = sample_meta[sample_id, "p90_cell_ploidy"],
      mean_pseudotime = sample_meta[sample_id, "mean_pseudotime"],
      median_pseudotime = sample_meta[sample_id, "median_pseudotime"],
      reference_n_samples = length(reference_samples),
      reference_cell_count = reference_cells,
      ecdf_rmse = if (length(finite_delta) > 0L) sqrt(mean(finite_delta^2)) else NA_real_,
      ecdf_ks = if (length(finite_delta) > 0L) max(abs(finite_delta)) else NA_real_,
      ecdf_mean_abs = if (length(finite_delta) > 0L) mean(abs(finite_delta)) else NA_real_,
      signed_mean_shift = sample_meta[sample_id, "mean_pseudotime"] - reference_mean,
      stringsAsFactors = FALSE
    )
    for (column in essential_endpoint_tgi_columns(names(sample_meta))) {
      row[[column]] <- sample_meta[sample_id, column]
    }
    row
  })
  output <- do.call(rbind, rows)
  output <- output[order(output$dose_mg, output$analysis_group_numeric, output$sample_id), , drop = FALSE]
  rownames(output) <- NULL
  output
}

essential_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(n = length(x), estimate = NA_real_, p_value = NA_real_))
  }
  test <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  data.frame(n = length(x), estimate = unname(test$estimate), p_value = test$p.value)
}

essential_permutation_indices <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(get(key, envir = cache, inherits = FALSE))
    }
    total <- as.integer(round(factorial(n)))
    matrix_values <- matrix(NA_integer_, nrow = total, ncol = n)
    row_index <- 0L
    recurse <- function(prefix, remaining) {
      if (length(remaining) == 0L) {
        row_index <<- row_index + 1L
        matrix_values[row_index, ] <<- prefix
      } else {
        for (i in seq_along(remaining)) {
          recurse(c(prefix, remaining[[i]]), remaining[-i])
        }
      }
    }
    recurse(integer(0), seq_len(n))
    assign(key, matrix_values, envir = cache)
    matrix_values
  }
})

essential_stratified_permutation_indices <- function(strata) {
  split_indices <- split(seq_along(strata), as.character(strata))
  local_permutations <- lapply(split_indices, function(indices) {
    local <- essential_permutation_indices(length(indices))
    t(apply(local, 1L, function(order_index) indices[order_index]))
  })
  output <- list()
  recurse <- function(i, current) {
    if (i > length(local_permutations)) {
      output[[length(output) + 1L]] <<- current
      return(invisible(NULL))
    }
    for (row in seq_len(nrow(local_permutations[[i]]))) {
      next_index <- current
      indices <- split_indices[[i]]
      next_index[indices] <- local_permutations[[i]][row, ]
      recurse(i + 1L, next_index)
    }
    invisible(NULL)
  }
  recurse(1L, seq_along(strata))
  do.call(rbind, output)
}

essential_permutation_cor <- function(
  x,
  y,
  method = "pearson",
  n_perm = 10000L,
  exact = FALSE,
  strata = NULL
) {
  keep <- is.finite(x) & is.finite(y)
  if (!is.null(strata)) keep <- keep & !is.na(strata)
  x <- x[keep]
  y <- y[keep]
  if (!is.null(strata)) strata <- as.character(strata[keep])
  n <- length(x)
  if (n < 3L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(data.frame(
      n = n, estimate = NA_real_, asymptotic_p = NA_real_,
      permutation_p_two_sided = NA_real_, permutation_p_positive = NA_real_,
      n_permutations = NA_integer_, permutation_mode = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  observed <- suppressWarnings(stats::cor(x, y, method = method))
  asymptotic <- essential_safe_cor(x, y, method)
  exact_count <- if (is.null(strata)) factorial(n) else prod(vapply(table(strata), factorial, numeric(1L)))
  if (isTRUE(exact) && is.finite(exact_count) && exact_count <= 400000) {
    index_matrix <- if (is.null(strata)) {
      essential_permutation_indices(n)
    } else {
      essential_stratified_permutation_indices(strata)
    }
    permuted <- vapply(seq_len(nrow(index_matrix)), function(i) {
      suppressWarnings(stats::cor(x, y[index_matrix[i, ]], method = method))
    }, numeric(1L))
    mode <- if (is.null(strata)) "exact_TGI_label_enumeration" else "exact_within_stratum_TGI_label_enumeration"
    p_two_sided <- mean(abs(permuted) >= abs(observed) - 1e-15, na.rm = TRUE)
    p_positive <- mean(permuted >= observed - 1e-15, na.rm = TRUE)
  } else {
    permuted <- replicate(n_perm, {
      permuted_y <- y
      if (is.null(strata)) {
        permuted_y <- sample(y)
      } else {
        for (level in unique(strata)) {
          indices <- which(strata == level)
          permuted_y[indices] <- sample(y[indices])
        }
      }
      suppressWarnings(stats::cor(x, permuted_y, method = method))
    })
    mode <- if (is.null(strata)) "monte_carlo_TGI_label_permutation" else "monte_carlo_within_stratum_TGI_label_permutation"
    exceed_two <- sum(abs(permuted) >= abs(observed) - 1e-15, na.rm = TRUE)
    exceed_positive <- sum(permuted >= observed - 1e-15, na.rm = TRUE)
    p_two_sided <- (exceed_two + 1) / (sum(is.finite(permuted)) + 1)
    p_positive <- (exceed_positive + 1) / (sum(is.finite(permuted)) + 1)
  }
  data.frame(
    n = n,
    estimate = observed,
    asymptotic_p = asymptotic$p_value,
    permutation_p_two_sided = p_two_sided,
    permutation_p_positive = p_positive,
    n_permutations = length(permuted),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

essential_choice_lists <- function(groups, group_a, strata) {
  split_indices <- split(seq_along(groups), strata)
  lapply(split_indices, function(indices) {
    k <- sum(groups[indices] == group_a)
    if (k == 0L) return(list(integer(0)))
    if (k == length(indices)) return(list(indices))
    utils::combn(indices, k, simplify = FALSE)
  })
}

essential_exact_group_vectors <- function(groups, group_a, group_b, strata) {
  choices <- essential_choice_lists(groups, group_a, strata)
  output <- list()
  recurse <- function(i, selected) {
    if (i > length(choices)) {
      vector <- rep(group_b, length(groups))
      vector[unlist(selected, use.names = FALSE)] <- group_a
      output[[length(output) + 1L]] <<- vector
    } else {
      for (choice in choices[[i]]) recurse(i + 1L, c(selected, list(choice)))
    }
  }
  recurse(1L, list())
  output
}

essential_sample_group_vector <- function(groups, group_a, group_b, strata) {
  output <- rep(group_b, length(groups))
  for (indices in split(seq_along(groups), strata)) {
    k <- sum(groups[indices] == group_a)
    if (k > 0L) output[sample(indices, k)] <- group_a
  }
  output
}

essential_ecdf_distance <- function(ecdf_matrix, groups, group_a, group_b) {
  if (!any(groups == group_a) || !any(groups == group_b)) {
    return(c(ecdf_rmse = NA_real_, ecdf_ks = NA_real_, ecdf_mean_abs = NA_real_))
  }
  delta <- colMeans(ecdf_matrix[groups == group_a, , drop = FALSE], na.rm = TRUE) -
    colMeans(ecdf_matrix[groups == group_b, , drop = FALSE], na.rm = TRUE)
  c(
    ecdf_rmse = sqrt(mean(delta^2, na.rm = TRUE)),
    ecdf_ks = max(abs(delta), na.rm = TRUE),
    ecdf_mean_abs = mean(abs(delta), na.rm = TRUE)
  )
}

essential_ecdf_group_test <- function(data, group_column, group_a, group_b, stratify, n_perm) {
  local <- data[data[[group_column]] %in% c(group_a, group_b), , drop = FALSE]
  sample_ids <- sort(unique(local$sample_id))
  meta_columns <- unique(c("sample_id", "analysis_group", group_column))
  meta <- local[!duplicated(local$sample_id), meta_columns, drop = FALSE]
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
  groups <- as.character(meta[[group_column]])
  not_estimable <- length(sample_ids) < 3L || !all(c(group_a, group_b) %in% groups)
  if (not_estimable) {
    return(data.frame(
      group_a = group_a, group_b = group_b, stratified_by_analysis_group = stratify,
      n_samples = length(sample_ids), n_group_a = sum(groups == group_a), n_group_b = sum(groups == group_b),
      observed_ecdf_rmse = NA_real_, observed_ecdf_ks = NA_real_, observed_ecdf_mean_abs = NA_real_,
      p_ecdf_rmse = NA_real_, p_ecdf_ks = NA_real_, p_ecdf_mean_abs = NA_real_,
      n_permutations = NA_integer_, permutation_mode = "not_estimable", stringsAsFactors = FALSE
    ))
  }
  strata <- if (stratify) as.character(meta$analysis_group) else rep("all", length(groups))
  grid <- essential_make_grid(local)
  ecdf_matrix <- essential_sample_ecdf_matrix(local, sample_ids, grid)
  observed <- essential_ecdf_distance(ecdf_matrix, groups, group_a, group_b)
  choices <- essential_choice_lists(groups, group_a, strata)
  n_exact <- prod(vapply(choices, length, integer(1L)))
  if (is.finite(n_exact) && n_exact <= n_perm) {
    permuted_groups <- essential_exact_group_vectors(groups, group_a, group_b, strata)
    mode <- "exact_label_enumeration"
  } else {
    permuted_groups <- replicate(
      n_perm,
      essential_sample_group_vector(groups, group_a, group_b, strata),
      simplify = FALSE
    )
    mode <- "monte_carlo_label_permutation"
  }
  permuted_stats <- t(vapply(permuted_groups, function(vector) {
    essential_ecdf_distance(ecdf_matrix, vector, group_a, group_b)
  }, numeric(3L)))
  data.frame(
    group_a = group_a,
    group_b = group_b,
    stratified_by_analysis_group = stratify,
    n_samples = length(sample_ids),
    n_group_a = sum(groups == group_a),
    n_group_b = sum(groups == group_b),
    observed_ecdf_rmse = observed[["ecdf_rmse"]],
    observed_ecdf_ks = observed[["ecdf_ks"]],
    observed_ecdf_mean_abs = observed[["ecdf_mean_abs"]],
    p_ecdf_rmse = mean(permuted_stats[, "ecdf_rmse"] >= observed[["ecdf_rmse"]] - 1e-15, na.rm = TRUE),
    p_ecdf_ks = mean(permuted_stats[, "ecdf_ks"] >= observed[["ecdf_ks"]] - 1e-15, na.rm = TRUE),
    p_ecdf_mean_abs = mean(permuted_stats[, "ecdf_mean_abs"] >= observed[["ecdf_mean_abs"]] - 1e-15, na.rm = TRUE),
    n_permutations = nrow(permuted_stats),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

essential_run_dose_tests <- function(data, n_perm) {
  local <- data
  local$combined_treatment_group <- ifelse(local$dose_mg == 0, "0mg/kg", "treated")
  output <- list(
    cbind(
      comparison = "0_vs_30plus120",
      essential_ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", FALSE, n_perm)
    ),
    cbind(
      comparison = "0_vs_30plus120",
      essential_ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", TRUE, n_perm)
    )
  )
  dose_values <- sort(unique(local$dose_mg))
  for (pair in split(utils::combn(dose_values, 2L), col(utils::combn(dose_values, 2L)))) {
    dose_a <- pair[[1L]]
    dose_b <- pair[[2L]]
    subset <- local[local$dose_mg %in% c(dose_a, dose_b), , drop = FALSE]
    group_a <- unique(subset$dose[subset$dose_mg == dose_a])[[1L]]
    group_b <- unique(subset$dose[subset$dose_mg == dose_b])[[1L]]
    output[[length(output) + 1L]] <- cbind(
      comparison = paste0(dose_a, "_vs_", dose_b),
      essential_ecdf_group_test(subset, "dose", group_a, group_b, FALSE, n_perm)
    )
  }
  do.call(rbind, output)
}

essential_run_tgi_associations <- function(equal_shift, pooled_shift, method, n_perm) {
  references <- list(
    primary_equal_sample_reference = equal_shift,
    sensitivity_pooled_cell_reference = pooled_shift
  )
  shift_columns <- c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")
  tgi_columns <- unique(c("TGI_percent_auc", essential_endpoint_tgi_columns(names(equal_shift))))
  rows <- list()
  for (reference_name in names(references)) {
    treated <- references[[reference_name]][references[[reference_name]]$dose_mg > 0, , drop = FALSE]
    for (shift_column in shift_columns) {
      for (tgi_column in tgi_columns) {
        exact <- identical(tgi_column, "TGI_percent_auc")
        pearson <- essential_permutation_cor(treated[[shift_column]], treated[[tgi_column]], "pearson", n_perm, exact)
        spearman <- essential_permutation_cor(treated[[shift_column]], treated[[tgi_column]], "spearman", n_perm, exact)
        rows[[length(rows) + 1L]] <- data.frame(
          method = method,
          compartment = "CellCycle",
          sample_set = "treated",
          reference_type = reference_name,
          shift_metric = shift_column,
          tgi_measure = tgi_column,
          tgi_day = if (grepl("^TGI_percent_Day_", tgi_column)) essential_endpoint_day(tgi_column) else NA_real_,
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
          pre_specified_primary = reference_name == "primary_equal_sample_reference" &&
            shift_column == "ecdf_rmse" && tgi_column == "TGI_percent_auc",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  output <- do.call(rbind, rows)
  output$pearson_q_exploratory <- stats::p.adjust(output$pearson_p_permutation_two_sided, method = "BH")
  output$spearman_q_exploratory <- stats::p.adjust(output$spearman_p_permutation_two_sided, method = "BH")
  output
}

essential_run_mean_etp_associations <- function(equal_shift, method, n_perm) {
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  tgi_columns <- unique(c("TGI_percent_auc", essential_endpoint_tgi_columns(names(treated))))
  rows <- lapply(tgi_columns, function(tgi_column) {
    exact <- identical(tgi_column, "TGI_percent_auc")
    pearson <- essential_permutation_cor(
      treated$sample_mean_endpoint_ploidy,
      treated[[tgi_column]],
      "pearson",
      n_perm,
      exact
    )
    spearman <- essential_permutation_cor(
      treated$sample_mean_endpoint_ploidy,
      treated[[tgi_column]],
      "spearman",
      n_perm,
      exact
    )
    data.frame(
      method = method,
      compartment = "CellCycle",
      sample_set = "treated",
      predictor = "sample_mean_endpoint_ploidy",
      predictor_label = "sample_mean_ETP",
      tgi_measure = tgi_column,
      tgi_day = if (grepl("^TGI_percent_Day_", tgi_column)) essential_endpoint_day(tgi_column) else NA_real_,
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
      multiplicity_family = "CellCycle_mean_ETP_TGI",
      analysis_role = if (exact) "parallel_mean_ETP_AUC" else "exploratory_mean_ETP",
      pre_specified_primary = FALSE,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  output$pearson_q_exploratory <- stats::p.adjust(output$pearson_p_permutation_two_sided, method = "BH")
  output$spearman_q_exploratory <- stats::p.adjust(output$spearman_p_permutation_two_sided, method = "BH")
  output
}

essential_leave_one_out <- function(data, predictor) {
  full_pearson <- essential_safe_cor(data[[predictor]], data$TGI_percent_auc, "pearson")
  rows <- lapply(data$sample_id, function(sample_id) {
    local <- data[data$sample_id != sample_id, , drop = FALSE]
    pearson <- essential_safe_cor(local[[predictor]], local$TGI_percent_auc, "pearson")
    spearman <- essential_safe_cor(local[[predictor]], local$TGI_percent_auc, "spearman")
    data.frame(
      omitted_sample_id = sample_id,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p = pearson$p_value,
      spearman_rho = spearman$estimate,
      spearman_p = spearman$p_value,
      sign_matches_full_pearson = sign(pearson$estimate) == sign(full_pearson$estimate),
      delta_pearson_r_from_full = pearson$estimate - full_pearson$estimate,
      influential_flag = is.finite(pearson$estimate) && is.finite(full_pearson$estimate) &&
        (sign(pearson$estimate) != sign(full_pearson$estimate) || abs(pearson$estimate - full_pearson$estimate) > 0.20),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

essential_bootstrap <- function(data, predictor, n_boot) {
  rows <- lapply(seq_len(n_boot), function(iteration) {
    indices <- sample(seq_len(nrow(data)), replace = TRUE)
    local <- data[indices, , drop = FALSE]
    pearson <- essential_safe_cor(local[[predictor]], local$TGI_percent_auc, "pearson")
    spearman <- essential_safe_cor(local[[predictor]], local$TGI_percent_auc, "spearman")
    data.frame(
      iteration = iteration,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p = pearson$p_value,
      spearman_rho = spearman$estimate,
      spearman_p = spearman$p_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

essential_quantile <- function(x, probability) {
  finite <- x[is.finite(x)]
  if (length(finite) == 0L) return(NA_real_)
  as.numeric(stats::quantile(finite, probability, names = FALSE))
}

essential_primary_robustness_summary <- function(association, leave_one_out, bootstrap) {
  primary <- association[association$pre_specified_primary, , drop = FALSE]
  data.frame(
    full_n = primary$n[[1L]],
    full_pearson_r = primary$pearson_r[[1L]],
    full_pearson_perm_p = primary$pearson_p_permutation_two_sided[[1L]],
    full_spearman_rho = primary$spearman_rho[[1L]],
    full_spearman_perm_p = primary$spearman_p_permutation_two_sided[[1L]],
    leave_one_out_all_positive = all(leave_one_out$pearson_r > 0, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(leave_one_out$pearson_p, na.rm = TRUE),
    n_influential_flags = sum(leave_one_out$influential_flag, na.rm = TRUE),
    bootstrap_pearson_r_ci025 = essential_quantile(bootstrap$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = essential_quantile(bootstrap$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = essential_quantile(bootstrap$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = essential_quantile(bootstrap$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

essential_mean_etp_robustness_summary <- function(method, association, leave_one_out, bootstrap, n_boot) {
  primary <- association[association$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  data.frame(
    method = method,
    predictor = "sample_mean_endpoint_ploidy",
    full_n = primary$n[[1L]],
    full_pearson_r = primary$pearson_r[[1L]],
    full_pearson_asymptotic_p = primary$pearson_p_asymptotic[[1L]],
    full_pearson_perm_p = primary$pearson_p_permutation_two_sided[[1L]],
    full_pearson_q = primary$pearson_q_exploratory[[1L]],
    full_spearman_rho = primary$spearman_rho[[1L]],
    full_spearman_asymptotic_p = primary$spearman_p_asymptotic[[1L]],
    full_spearman_perm_p = primary$spearman_p_permutation_two_sided[[1L]],
    full_spearman_q = primary$spearman_q_exploratory[[1L]],
    leave_one_out_all_same_sign = all(leave_one_out$sign_matches_full_pearson, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_r = max(leave_one_out$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(leave_one_out$pearson_p, na.rm = TRUE),
    leave_one_out_max_abs_delta_r = max(abs(leave_one_out$delta_pearson_r_from_full), na.rm = TRUE),
    n_influential_flags = sum(leave_one_out$influential_flag, na.rm = TRUE),
    n_bootstrap_requested = n_boot,
    n_bootstrap_finite_pearson = sum(is.finite(bootstrap$pearson_r)),
    n_bootstrap_finite_spearman = sum(is.finite(bootstrap$spearman_rho)),
    bootstrap_pearson_r_ci025 = essential_quantile(bootstrap$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = essential_quantile(bootstrap$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = essential_quantile(bootstrap$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = essential_quantile(bootstrap$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

essential_lm_terms <- function(fit, analysis, model_label, n) {
  coefficient_matrix <- summary(fit)$coefficients
  data.frame(
    compartment = "CellCycle",
    analysis = analysis,
    model = model_label,
    term = rownames(coefficient_matrix),
    estimate = coefficient_matrix[, 1L],
    std_error = coefficient_matrix[, 2L],
    statistic = coefficient_matrix[, 3L],
    p_value = coefficient_matrix[, 4L],
    n = n,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

essential_fit_shift_models <- function(treated, spec) {
  local <- treated[is.finite(treated$ecdf_rmse) & is.finite(treated$TGI_percent_auc), , drop = FALSE]
  group_column <- spec$group_column
  local[[group_column]] <- factor(as.character(local$analysis_group), levels = spec$group_levels)
  local$dose_mg_factor <- factor(local$dose_mg, levels = sort(unique(local$dose_mg)))
  formulas <- c(
    shift_only = "TGI_percent_auc ~ scale(ecdf_rmse)",
    shift_ploidy_dose_adjusted = paste0(
      "TGI_percent_auc ~ scale(ecdf_rmse) + ", group_column, " + dose_mg_factor"
    ),
    shift_by_ploidy_dose_adjusted = paste0(
      "TGI_percent_auc ~ scale(ecdf_rmse) * ", group_column, " + dose_mg_factor"
    )
  )
  models <- list()
  summaries <- list()
  for (analysis in names(formulas)) {
    fit <- stats::lm(stats::as.formula(formulas[[analysis]]), data = local)
    models[[analysis]] <- essential_lm_terms(fit, analysis, formulas[[analysis]], nrow(local))
    summary_fit <- summary(fit)
    summaries[[analysis]] <- data.frame(
      compartment = "CellCycle",
      analysis = analysis,
      model = formulas[[analysis]],
      n = stats::nobs(fit),
      residual_df = stats::df.residual(fit),
      r_squared = unname(summary_fit$r.squared),
      adjusted_r_squared = unname(summary_fit$adj.r.squared),
      sigma = unname(summary_fit$sigma),
      aic = stats::AIC(fit),
      stringsAsFactors = FALSE
    )
  }
  model_table <- do.call(rbind, models)
  model_table$pre_specified_primary <- model_table$analysis == "shift_ploidy_dose_adjusted" &
    model_table$term == "scale(ecdf_rmse)"
  list(models = model_table, summaries = do.call(rbind, summaries))
}

essential_confounding <- function(equal_shift, spec, n_perm) {
  data <- equal_shift
  data[[spec$group_numeric_column]] <- data$analysis_group_numeric
  ploidy_columns <- c(
    "mean_cell_ploidy", "median_cell_ploidy", "p90_cell_ploidy", spec$group_numeric_column
  )
  rows <- list()
  for (sample_set in c("all", "treated")) {
    local <- if (sample_set == "treated") data[data$dose_mg > 0, , drop = FALSE] else data
    for (shift_metric in c("ecdf_rmse", "signed_mean_shift")) {
      for (ploidy_column in ploidy_columns) {
        for (method in c("pearson", "spearman")) {
          asymptotic <- essential_safe_cor(local[[shift_metric]], local[[ploidy_column]], method)
          permutation <- essential_permutation_cor(
            local[[shift_metric]], local[[ploidy_column]], method, n_perm, exact = FALSE
          )
          rows[[length(rows) + 1L]] <- data.frame(
            compartment = "CellCycle",
            sample_set = sample_set,
            shift_metric = shift_metric,
            ploidy_measure = ploidy_column,
            method = method,
            n = asymptotic$n,
            estimate = asymptotic$estimate,
            p_value = asymptotic$p_value,
            p_permutation_two_sided = permutation$permutation_p_two_sided,
            n_permutations = permutation$n_permutations,
            permutation_mode = permutation$permutation_mode,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  ploidy_tests <- do.call(rbind, rows)

  treated <- data[data$dose_mg > 0, , drop = FALSE]
  centered <- treated
  centered$shift_centered <- centered$ecdf_rmse - ave(centered$ecdf_rmse, centered$dose_mg, FUN = mean)
  centered$tgi_centered <- centered$TGI_percent_auc - ave(centered$TGI_percent_auc, centered$dose_mg, FUN = mean)
  centered_cor <- essential_safe_cor(centered$shift_centered, centered$tgi_centered, "pearson")
  centered_perm <- essential_permutation_cor(
    centered$shift_centered,
    centered$tgi_centered,
    "pearson",
    n_perm,
    exact = TRUE,
    strata = centered$dose_mg
  )

  adjusted <- treated[is.finite(treated$ecdf_rmse) & is.finite(treated$TGI_percent_auc), , drop = FALSE]
  shift_residual <- stats::resid(stats::lm(
    ecdf_rmse ~ factor(dose_mg) + analysis_group_numeric,
    data = adjusted
  ))
  tgi_residual <- stats::resid(stats::lm(
    TGI_percent_auc ~ factor(dose_mg) + analysis_group_numeric,
    data = adjusted
  ))
  adjusted_cor <- essential_safe_cor(shift_residual, tgi_residual, "pearson")
  adjusted_perm <- essential_permutation_cor(shift_residual, tgi_residual, "pearson", n_perm, exact = TRUE)
  residualized <- rbind(
    data.frame(
      compartment = "CellCycle",
      analysis = paste0("residualized_against_dose_and_", spec$group_column),
      method = "pearson", n = adjusted_cor$n, estimate = adjusted_cor$estimate,
      p_value = adjusted_cor$p_value,
      p_permutation_two_sided = adjusted_perm$permutation_p_two_sided,
      n_permutations = adjusted_perm$n_permutations,
      permutation_mode = adjusted_perm$permutation_mode,
      stringsAsFactors = FALSE
    ),
    data.frame(
      compartment = "CellCycle", analysis = "within_dose_centered", method = "pearson",
      n = centered_cor$n, estimate = centered_cor$estimate, p_value = centered_cor$p_value,
      p_permutation_two_sided = centered_perm$permutation_p_two_sided,
      n_permutations = centered_perm$n_permutations,
      permutation_mode = centered_perm$permutation_mode,
      stringsAsFactors = FALSE
    )
  )
  list(ploidy_tests = ploidy_tests, residualized = residualized, centered = centered)
}

essential_plot_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      legend.key.height = grid::unit(0.45, "cm"),
      aspect.ratio = 1
    )
}

essential_dose_colors <- function() {
  c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77", "treated" = "#377eb8")
}

essential_group_colors <- function(spec) {
  if (identical(spec$type, "ETP")) {
    c("ETP-lower" = "#4C78A8", "ETP-higher" = "#E45756")
  } else {
    c("2N" = "#4C78A8", "4N" = "#E45756")
  }
}

essential_association_annotation <- function(association) {
  paste0(
    "Pearson r = ", essential_format_number(association$estimate[[1L]], 3),
    "\nPermutation P = ", essential_format_p(association$permutation_p_two_sided[[1L]]),
    "\nn = ", association$n[[1L]]
  )
}

essential_attach_association <- function(data, analysis_id, pearson, spearman) {
  data$analysis_id <- analysis_id
  data$pearson_r <- pearson$estimate[[1L]]
  data$pearson_p_asymptotic <- pearson$asymptotic_p[[1L]]
  data$pearson_p_permutation_two_sided <- pearson$permutation_p_two_sided[[1L]]
  data$pearson_p_permutation_positive <- pearson$permutation_p_positive[[1L]]
  data$pearson_n_permutations <- pearson$n_permutations[[1L]]
  data$pearson_permutation_mode <- pearson$permutation_mode[[1L]]
  data$spearman_rho <- spearman$estimate[[1L]]
  data$spearman_p_asymptotic <- spearman$asymptotic_p[[1L]]
  data$spearman_p_permutation_two_sided <- spearman$permutation_p_two_sided[[1L]]
  data$spearman_p_permutation_positive <- spearman$permutation_p_positive[[1L]]
  data$spearman_n_permutations <- spearman$n_permutations[[1L]]
  data$spearman_permutation_mode <- spearman$permutation_mode[[1L]]
  data$plot_annotation <- essential_association_annotation(pearson)
  data
}

essential_mean_ecdf_curves <- function(data, panel, grid) {
  rows <- lapply(unique(as.character(data$curve_label)), function(curve_label) {
    local <- data[as.character(data$curve_label) == curve_label, , drop = FALSE]
    sample_ids <- unique(local$sample_id)
    if (length(sample_ids) == 0L) return(NULL)
    matrix_values <- t(vapply(sample_ids, function(sample_id) {
      stats::ecdf(local$pseudotime[local$sample_id == sample_id])(grid)
    }, numeric(length(grid))))
    data.frame(
      panel = panel,
      pseudotime = grid,
      mean_ecdf = colMeans(matrix_values, na.rm = TRUE),
      curve_label = curve_label,
      color_group = as.character(local$color_group[[1L]]),
      line_group = as.character(local$line_group[[1L]]),
      n_samples = length(sample_ids),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

essential_direct_ecdf <- function(data, dose_tests, spec, n_perm) {
  data$treatment_group <- ifelse(data$dose_mg == 0, "0mg/kg", "treated")
  grid <- essential_make_grid(data)
  curve_rows <- list()
  add_curves <- function(local, panel, curve_label, color_group, line_group) {
    if (nrow(local) == 0L) return(invisible(NULL))
    local$curve_label <- curve_label
    local$color_group <- color_group
    local$line_group <- line_group
    curve_rows[[length(curve_rows) + 1L]] <<- essential_mean_ecdf_curves(local, panel, grid)
    invisible(NULL)
  }
  group_low <- spec$group_levels[[1L]]
  group_high <- spec$group_levels[[2L]]
  stratified_panel <- paste0("2. 0 vs treated, ", if (identical(spec$type, "ETP")) "ETP" else "initial-ploidy", "-stratified")
  panel_six <- paste0("6. Treated: ", group_low, " vs ", group_high)
  panel_seven <- paste0("7. Control: ", group_low, " vs ", group_high)
  panel_eight <- paste0("8. ", group_high, ": 0 vs treated")
  panel_nine <- paste0("9. ", group_low, ": 0 vs treated")

  add_curves(data, "1. 0 vs treated", data$treatment_group, data$treatment_group, "All")
  add_curves(
    data,
    stratified_panel,
    paste(data$treatment_group, data$analysis_group, sep = " | "),
    data$treatment_group,
    as.character(data$analysis_group)
  )
  dose_panels <- list(
    list(panel = "3. 0 vs 30 mg/kg", doses = c(0, 30)),
    list(panel = "4. 0 vs 120 mg/kg", doses = c(0, 120)),
    list(panel = "5. 30 vs 120 mg/kg", doses = c(30, 120))
  )
  for (definition in dose_panels) {
    local <- data[data$dose_mg %in% definition$doses, , drop = FALSE]
    add_curves(local, definition$panel, local$dose, local$dose, "All")
  }
  treated <- data[data$dose_mg > 0, , drop = FALSE]
  control <- data[data$dose_mg == 0, , drop = FALSE]
  add_curves(treated, panel_six, treated$analysis_group, "treated", treated$analysis_group)
  add_curves(control, panel_seven, control$analysis_group, "0mg/kg", control$analysis_group)
  high <- data[data$analysis_group == group_high, , drop = FALSE]
  low <- data[data$analysis_group == group_low, , drop = FALSE]
  add_curves(high, panel_eight, high$treatment_group, high$treatment_group, group_high)
  add_curves(low, panel_nine, low$treatment_group, low$treatment_group, group_low)
  curves <- do.call(rbind, curve_rows)

  get_test <- function(comparison, stratified) {
    hit <- dose_tests[
      dose_tests$comparison == comparison & dose_tests$stratified_by_analysis_group == stratified,
      ,
      drop = FALSE
    ]
    if (nrow(hit) == 0L) data.frame() else hit[1L, , drop = FALSE]
  }
  between_group <- function(local) {
    essential_ecdf_group_test(local, "analysis_group", group_low, group_high, FALSE, n_perm)
  }
  treatment_within_group <- function(group) {
    local <- data[data$analysis_group == group, , drop = FALSE]
    essential_ecdf_group_test(local, "treatment_group", "0mg/kg", "treated", FALSE, n_perm)
  }
  test_rows <- list(
    cbind(panel = "1. 0 vs treated", test_design = "sample_label_permutation", get_test("0_vs_30plus120", FALSE)),
    cbind(panel = stratified_panel, test_design = "treatment_label_permutation_within_group", get_test("0_vs_30plus120", TRUE)),
    cbind(panel = "3. 0 vs 30 mg/kg", test_design = "sample_label_permutation", get_test("0_vs_30", FALSE)),
    cbind(panel = "4. 0 vs 120 mg/kg", test_design = "sample_label_permutation", get_test("0_vs_120", FALSE)),
    cbind(panel = "5. 30 vs 120 mg/kg", test_design = "sample_label_permutation", get_test("30_vs_120", FALSE)),
    cbind(panel = panel_six, test_design = "unpaired_group_label_permutation", between_group(treated)),
    cbind(panel = panel_seven, test_design = "unpaired_group_label_permutation", between_group(control)),
    cbind(panel = panel_eight, test_design = "treatment_label_permutation_within_group", treatment_within_group(group_high)),
    cbind(panel = panel_nine, test_design = "treatment_label_permutation_within_group", treatment_within_group(group_low))
  )
  all_columns <- unique(unlist(lapply(test_rows, names), use.names = FALSE))
  test_rows <- lapply(test_rows, function(row) {
    for (column in setdiff(all_columns, names(row))) row[[column]] <- NA
    row[, all_columns, drop = FALSE]
  })
  tests <- do.call(rbind, test_rows)
  tests$method <- spec$method
  tests$q_ecdf_rmse_exploratory <- stats::p.adjust(tests$p_ecdf_rmse, method = "BH")
  tests$q_ecdf_ks_exploratory <- stats::p.adjust(tests$p_ecdf_ks, method = "BH")
  tests$q_ecdf_mean_abs_exploratory <- stats::p.adjust(tests$p_ecdf_mean_abs, method = "BH")
  tests$threshold <- spec$threshold
  tests$estimable <- is.finite(tests$p_ecdf_rmse)
  tests$annotation <- paste0(
    tests$group_a, " n=", tests$n_group_a, "; ", tests$group_b, " n=", tests$n_group_b,
    "\nRMSE=", essential_format_number(tests$observed_ecdf_rmse, 4), ", P=", essential_format_p(tests$p_ecdf_rmse),
    "\nKS=", essential_format_number(tests$observed_ecdf_ks, 4), ", P=", essential_format_p(tests$p_ecdf_ks)
  )
  panel_levels <- c(
    "1. 0 vs treated", stratified_panel, "3. 0 vs 30 mg/kg", "4. 0 vs 120 mg/kg",
    "5. 30 vs 120 mg/kg", panel_six, panel_seven, panel_eight, panel_nine
  )
  curves$panel <- factor(curves$panel, levels = panel_levels)
  curves$color_group <- factor(curves$color_group, levels = names(essential_dose_colors()))
  curves$line_group <- factor(curves$line_group, levels = c("All", spec$group_levels))
  tests$panel <- factor(tests$panel, levels = panel_levels)
  plot_data <- merge(curves, tests, by = "panel", all.x = TRUE, sort = FALSE)
  plot_data$method <- spec$method
  line_values <- c("All" = "solid")
  line_values[spec$group_levels] <- c("solid", "22")
  plot <- ggplot2::ggplot(
    curves,
    ggplot2::aes(x = pseudotime, y = mean_ecdf, color = color_group, linetype = line_group, group = curve_label)
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_text(
      data = tests,
      ggplot2::aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.35,
      lineheight = 0.92
    ) +
    ggplot2::facet_wrap(~panel, ncol = 3L, drop = FALSE) +
    ggplot2::scale_color_manual(values = essential_dose_colors(), name = "Group") +
    ggplot2::scale_linetype_manual(values = line_values, name = spec$group_label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    ggplot2::labs(
      title = paste0("CellCycle: direct group mean ECDF comparisons (", spec$method, ")"),
      subtitle = "Equal-sample mean ECDFs",
      x = "Pseudotime",
      y = "Mean ECDF"
    ) +
    essential_plot_theme() +
    ggplot2::theme(legend.position = "bottom", aspect.ratio = NULL)
  list(plot = plot, tests = tests, plot_data = plot_data)
}

essential_scatter_plot <- function(
  plot_data,
  x_column,
  y_column,
  color_column,
  shape_column,
  color_values,
  title,
  x_label,
  y_label,
  color_label,
  shape_label,
  annotation,
  label_nudge_y = 2.5
) {
  label_layer <- ggrepel::geom_text_repel(
    data = plot_data,
    ggplot2::aes(
      x = .data[[x_column]],
      y = .data[[y_column]],
      label = sample_id,
      color = .data[[color_column]]
    ),
    inherit.aes = FALSE,
    nudge_y = label_nudge_y,
    size = 2.4,
    seed = 1,
    box.padding = 0.22,
    point.padding = 0.12,
    min.segment.length = 0,
    segment.color = "grey65",
    segment.size = 0.25,
    max.overlaps = Inf,
    show.legend = FALSE
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]], color = .data[[color_column]], shape = .data[[shape_column]])
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    ggplot2::geom_smooth(
      data = plot_data,
      ggplot2::aes(x = .data[[x_column]], y = .data[[y_column]], group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = TRUE,
      color = "black",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(size = 2.9) +
    label_layer +
    ggplot2::annotate(
      "label",
      x = Inf,
      y = Inf,
      label = annotation,
      hjust = 1.05,
      vjust = 1.1,
      size = 3.0,
      linewidth = 0.2
    ) +
    ggplot2::scale_color_manual(values = color_values, name = color_label) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.15, 0.12))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.15))) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(title = title, x = x_label, y = y_label, shape = shape_label) +
    essential_plot_theme() +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 14, 8, 14))
}

essential_figure_filenames <- function() {
  c(
    "CellCycle_direct_group_ecdf_comparisons.pdf",
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf",
    "CellCycle_TGI_AUC_vs_mean_ETP.pdf",
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf",
    "CellCycle_TGI_association_within_dose_centered.pdf",
    "CellCycle_ecdf_rmse_vs_ploidy.pdf"
  )
}

essential_read_figure_table <- function(path, required_columns, label) {
  if (!file.exists(path)) stop("Missing figure-only table: ", path, call. = FALSE)
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required_columns, names(data))
  if (length(missing) > 0L) {
    stop(
      label, " is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(data) == 0L) stop(label, " contains no rows: ", path, call. = FALSE)
  data
}

essential_require_single_row <- function(data, keep, label) {
  keep[is.na(keep)] <- FALSE
  selected <- data[keep, , drop = FALSE]
  if (nrow(selected) != 1L) {
    stop(label, " must select exactly one statistics row; found ", nrow(selected), call. = FALSE)
  }
  selected
}

essential_prefixed_association <- function(row, prefix) {
  estimate_column <- if (identical(prefix, "pearson")) "pearson_r" else "spearman_rho"
  required <- c(
    "n", estimate_column, paste0(prefix, "_p_asymptotic"),
    paste0(prefix, "_p_permutation_two_sided"),
    paste0(prefix, "_p_permutation_positive"),
    paste0(prefix, "_n_permutations"), paste0(prefix, "_permutation_mode")
  )
  missing <- setdiff(required, names(row))
  if (length(missing) > 0L) {
    stop("Statistics row is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    n = essential_safe_numeric(row$n[[1L]]),
    estimate = essential_safe_numeric(row[[estimate_column]][[1L]]),
    asymptotic_p = essential_safe_numeric(row[[paste0(prefix, "_p_asymptotic")]][[1L]]),
    permutation_p_two_sided = essential_safe_numeric(
      row[[paste0(prefix, "_p_permutation_two_sided")]][[1L]]
    ),
    permutation_p_positive = essential_safe_numeric(
      row[[paste0(prefix, "_p_permutation_positive")]][[1L]]
    ),
    n_permutations = essential_safe_numeric(row[[paste0(prefix, "_n_permutations")]][[1L]]),
    permutation_mode = as.character(row[[paste0(prefix, "_permutation_mode")]][[1L]]),
    stringsAsFactors = FALSE
  )
}

essential_generic_association <- function(row) {
  required <- c(
    "n", "estimate", "p_value", "p_permutation_two_sided",
    "n_permutations", "permutation_mode"
  )
  missing <- setdiff(required, names(row))
  if (length(missing) > 0L) {
    stop("Statistics row is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  data.frame(
    n = essential_safe_numeric(row$n[[1L]]),
    estimate = essential_safe_numeric(row$estimate[[1L]]),
    asymptotic_p = essential_safe_numeric(row$p_value[[1L]]),
    permutation_p_two_sided = essential_safe_numeric(row$p_permutation_two_sided[[1L]]),
    permutation_p_positive = NA_real_,
    n_permutations = essential_safe_numeric(row$n_permutations[[1L]]),
    permutation_mode = as.character(row$permutation_mode[[1L]]),
    stringsAsFactors = FALSE
  )
}

essential_plot_data_association <- function(data, prefix) {
  row <- data[1L, , drop = FALSE]
  row$n <- nrow(data)
  essential_prefixed_association(row, prefix)
}

essential_restore_scatter_factors <- function(data, spec) {
  required <- c("sample_id", "dose", "dose_mg", "analysis_group")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop("Scatter plot data are missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  dose_order <- unique(as.character(data$dose[order(essential_safe_numeric(data$dose_mg))]))
  data$dose <- factor(as.character(data$dose), levels = dose_order)
  data$analysis_group <- factor(as.character(data$analysis_group), levels = spec$group_levels)
  data
}

essential_direct_plot_from_tables <- function(curves, tests, spec) {
  curve_columns <- c("panel", "pseudotime", "mean_ecdf", "curve_label", "color_group", "line_group")
  test_columns <- c("panel", "annotation")
  curve_missing <- setdiff(curve_columns, names(curves))
  test_missing <- setdiff(test_columns, names(tests))
  if (length(curve_missing) > 0L) {
    stop("Direct ECDF plot data are missing columns: ", paste(curve_missing, collapse = ", "), call. = FALSE)
  }
  if (length(test_missing) > 0L) {
    stop("Direct ECDF statistics are missing columns: ", paste(test_missing, collapse = ", "), call. = FALSE)
  }
  panel_levels <- unique(as.character(tests$panel))
  curves$panel <- factor(as.character(curves$panel), levels = panel_levels)
  curves$color_group <- factor(as.character(curves$color_group), levels = names(essential_dose_colors()))
  curves$line_group <- factor(as.character(curves$line_group), levels = c("All", spec$group_levels))
  tests$panel <- factor(as.character(tests$panel), levels = panel_levels)
  line_values <- c("All" = "solid")
  line_values[spec$group_levels] <- c("solid", "22")
  ggplot2::ggplot(
    curves,
    ggplot2::aes(
      x = pseudotime, y = mean_ecdf, color = color_group,
      linetype = line_group, group = curve_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_text(
      data = tests,
      ggplot2::aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.35,
      lineheight = 0.92
    ) +
    ggplot2::facet_wrap(~panel, ncol = 3L, drop = FALSE) +
    ggplot2::scale_color_manual(values = essential_dose_colors(), name = "Group") +
    ggplot2::scale_linetype_manual(values = line_values, name = spec$group_label) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    ggplot2::labs(
      title = paste0("CellCycle: direct group mean ECDF comparisons (", spec$method, ")"),
      subtitle = "Equal-sample mean ECDFs",
      x = "Pseudotime",
      y = "Mean ECDF"
    ) +
    essential_plot_theme() +
    ggplot2::theme(legend.position = "bottom", aspect.ratio = NULL)
}

essential_regenerate_figures_from_tables <- function(tables_root, output_root, spec) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)

  plot_data_dir <- file.path(tables_root, "plot_data", spec$method)
  stats_dir <- file.path(tables_root, "stats", spec$method)
  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", spec$method))
  read_plot_data <- function(stem, required) {
    essential_read_figure_table(
      file.path(plot_data_dir, paste0(stem, "_plot_data.csv")),
      required,
      paste0(stem, " plot data")
    )
  }
  read_stats <- function(filename, required) {
    essential_read_figure_table(
      file.path(stats_dir, filename), required, paste0(filename, " statistics")
    )
  }

  association_stats <- read_stats(
    "CellCycle_TGI_associations_ecdf_rmse.csv",
    c(
      "reference_type", "shift_metric", "tgi_measure", "n", "pearson_r",
      "pearson_p_asymptotic", "pearson_p_permutation_two_sided",
      "pearson_p_permutation_positive", "pearson_n_permutations",
      "pearson_permutation_mode", "spearman_rho", "spearman_p_asymptotic",
      "spearman_p_permutation_two_sided", "spearman_p_permutation_positive",
      "spearman_n_permutations", "spearman_permutation_mode"
    )
  )
  primary_stats <- essential_require_single_row(
    association_stats,
    association_stats$reference_type == "primary_equal_sample_reference" &
      association_stats$shift_metric == "ecdf_rmse" &
      association_stats$tgi_measure == "TGI_percent_auc",
    "Primary ECDF-TGI association"
  )
  primary_pearson <- essential_prefixed_association(primary_stats, "pearson")
  primary_spearman <- essential_prefixed_association(primary_stats, "spearman")

  primary_data <- read_plot_data(
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    c("sample_id", "ecdf_rmse", "TGI_percent_auc", "dose", "dose_mg", "analysis_group")
  )
  primary_data <- essential_restore_scatter_factors(primary_data, spec)
  primary_data <- essential_attach_association(
    primary_data,
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    primary_pearson,
    primary_spearman
  )
  primary_plot <- essential_scatter_plot(
    primary_data,
    "ecdf_rmse",
    "TGI_percent_auc",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "Cell-cycle-associated tumor cells: AUC TGI vs sample-equal ECDF shift",
    "ECDF RMSE from group-matched equal-sample 0 mg/kg reference",
    "AUC-based TGI (%)",
    "Dose",
    spec$group_label,
    primary_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    primary_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf"),
    6.8,
    6.8
  )

  direct_data <- read_plot_data(
    "CellCycle_direct_group_ecdf_comparisons",
    c("panel", "pseudotime", "mean_ecdf", "curve_label", "color_group", "line_group")
  )
  direct_stats <- read_stats(
    "CellCycle_direct_group_ecdf_comparisons_9panel_tests.csv",
    c("panel", "annotation", "p_ecdf_rmse", "p_ecdf_ks")
  )
  direct_plot <- essential_direct_plot_from_tables(direct_data, direct_stats, spec)
  essential_save_pdf(
    direct_plot,
    file.path(figures_dir, "CellCycle_direct_group_ecdf_comparisons.pdf"),
    15,
    11.5
  )

  mean_stats_all <- read_stats(
    "CellCycle_TGI_associations_mean_ETP.csv",
    c(
      "tgi_measure", "n", "pearson_r", "pearson_p_asymptotic",
      "pearson_p_permutation_two_sided", "pearson_p_permutation_positive",
      "pearson_n_permutations", "pearson_permutation_mode", "spearman_rho",
      "spearman_p_asymptotic", "spearman_p_permutation_two_sided",
      "spearman_p_permutation_positive", "spearman_n_permutations",
      "spearman_permutation_mode"
    )
  )
  mean_stats <- essential_require_single_row(
    mean_stats_all,
    mean_stats_all$tgi_measure == "TGI_percent_auc",
    "Mean ETP-TGI association"
  )
  mean_data <- read_plot_data(
    "CellCycle_TGI_AUC_vs_mean_ETP",
    c(
      "sample_id", "sample_mean_endpoint_ploidy", "TGI_percent_auc",
      "dose", "dose_mg", "analysis_group"
    )
  )
  mean_data <- essential_restore_scatter_factors(mean_data, spec)
  mean_data <- essential_attach_association(
    mean_data,
    "CellCycle_TGI_AUC_vs_mean_ETP",
    essential_prefixed_association(mean_stats, "pearson"),
    essential_prefixed_association(mean_stats, "spearman")
  )
  mean_plot <- essential_scatter_plot(
    mean_data,
    "sample_mean_endpoint_ploidy",
    "TGI_percent_auc",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "Cell-cycle-associated tumor cells: AUC TGI vs sample mean ETP",
    "Sample mean ETP",
    "AUC-based TGI (%)",
    "Dose",
    spec$group_label,
    mean_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    mean_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_mean_ETP.pdf"),
    6.8,
    6.8
  )

  model_data <- read_plot_data(
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose",
    c("sample_id", "ecdf_rmse", "TGI_percent_auc", "dose", "dose_mg", "analysis_group")
  )
  model_data <- essential_restore_scatter_factors(model_data, spec)
  model_data <- essential_attach_association(
    model_data,
    "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose",
    primary_pearson,
    primary_spearman
  )
  model_plot <- essential_scatter_plot(
    model_data,
    "ecdf_rmse",
    "TGI_percent_auc",
    "analysis_group",
    "dose",
    essential_group_colors(spec),
    "CellCycle: AUC-TGI association with pseudotime shift and ploidy",
    "ECDF RMSE from group-matched 0 mg/kg reference",
    "AUC-based TGI (%)",
    spec$group_label,
    "Dose",
    model_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    model_plot,
    file.path(figures_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf"),
    6.8,
    6.8
  )

  centered_stats_all <- read_stats(
    "residualized_TGI_associations.csv",
    c(
      "analysis", "method", "n", "estimate", "p_value",
      "p_permutation_two_sided", "n_permutations", "permutation_mode"
    )
  )
  centered_stats <- essential_require_single_row(
    centered_stats_all,
    centered_stats_all$analysis == "within_dose_centered" & centered_stats_all$method == "pearson",
    "Within-dose-centered TGI association"
  )
  centered_data <- read_plot_data(
    "CellCycle_TGI_association_within_dose_centered",
    c(
      "sample_id", "shift_centered", "tgi_centered", "dose", "dose_mg",
      "analysis_group", "spearman_rho", "spearman_p_asymptotic",
      "spearman_p_permutation_two_sided", "spearman_p_permutation_positive",
      "spearman_n_permutations", "spearman_permutation_mode"
    )
  )
  centered_data <- essential_restore_scatter_factors(centered_data, spec)
  centered_data <- essential_attach_association(
    centered_data,
    "CellCycle_TGI_association_within_dose_centered",
    essential_generic_association(centered_stats),
    essential_plot_data_association(centered_data, "spearman")
  )
  centered_plot <- essential_scatter_plot(
    centered_data,
    "shift_centered",
    "tgi_centered",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle TGI association after within-dose centering",
    "Dose-centered ECDF RMSE",
    "Dose-centered AUC TGI",
    "Dose",
    spec$group_label,
    centered_data$plot_annotation[[1L]],
    label_nudge_y = 1.0
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35)
  essential_save_pdf(
    centered_plot,
    file.path(figures_dir, "CellCycle_TGI_association_within_dose_centered.pdf"),
    6.6,
    6.6
  )

  ploidy_stats_all <- read_stats(
    "ploidy_confounding_tests.csv",
    c(
      "sample_set", "shift_metric", "ploidy_measure", "method", "n",
      "estimate", "p_value", "p_permutation_two_sided", "n_permutations",
      "permutation_mode"
    )
  )
  ploidy_primary <- ploidy_stats_all[
    ploidy_stats_all$sample_set == "all" &
      ploidy_stats_all$shift_metric == "ecdf_rmse" &
      ploidy_stats_all$ploidy_measure == "mean_cell_ploidy",
    ,
    drop = FALSE
  ]
  ploidy_pearson <- essential_require_single_row(
    ploidy_primary, ploidy_primary$method == "pearson", "Ploidy Pearson association"
  )
  ploidy_spearman <- essential_require_single_row(
    ploidy_primary, ploidy_primary$method == "spearman", "Ploidy Spearman association"
  )
  ploidy_data <- read_plot_data(
    "CellCycle_ecdf_rmse_vs_ploidy",
    c("sample_id", "mean_cell_ploidy", "ecdf_rmse", "dose", "dose_mg", "analysis_group")
  )
  ploidy_data <- ploidy_data[
    is.finite(essential_safe_numeric(ploidy_data$mean_cell_ploidy)) &
      is.finite(essential_safe_numeric(ploidy_data$ecdf_rmse)),
    ,
    drop = FALSE
  ]
  ploidy_data <- essential_restore_scatter_factors(ploidy_data, spec)
  ploidy_data <- essential_attach_association(
    ploidy_data,
    "CellCycle_ecdf_rmse_vs_ploidy",
    essential_generic_association(ploidy_pearson),
    essential_generic_association(ploidy_spearman)
  )
  ploidy_plot <- essential_scatter_plot(
    ploidy_data,
    "mean_cell_ploidy",
    "ecdf_rmse",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle ECDF RMSE vs mean cell ploidy",
    "Mean cell ploidy",
    "ECDF RMSE",
    "Dose",
    spec$group_label,
    ploidy_data$plot_annotation[[1L]],
    label_nudge_y = 0.004
  )
  essential_save_pdf(
    ploidy_plot,
    file.path(figures_dir, "CellCycle_ecdf_rmse_vs_ploidy.pdf"),
    6.4,
    6.4
  )

  invisible(list(method = spec$method, figures = length(essential_figure_filenames())))
}

essential_validate_figures_only_inventory <- function(output_root, methods) {
  expected <- sort(essential_figure_filenames())
  for (method in methods) {
    method_dir <- file.path(output_root, "Figures", method)
    actual <- sort(list.files(method_dir, pattern = "[.]pdf$", full.names = FALSE))
    other <- list.files(method_dir, all.files = TRUE, no.. = TRUE, full.names = FALSE)
    other <- setdiff(other, actual)
    if (!identical(actual, expected)) {
      stop("Unexpected figure-only PDF inventory for method: ", method, call. = FALSE)
    }
    if (length(other) > 0L) {
      stop("Unexpected non-PDF figure-only output for method: ", method, call. = FALSE)
    }
  }
  invisible(list(figures = length(expected) * length(methods)))
}

essential_write_readme <- function(output_root) {
  lines <- c(
    "# Essential pseudotime-TGI analysis",
    "",
    "This directory was generated directly from the two 04h cell-level input CSV files.",
    "",
    "The standalone workflow uses base R plus the `ggplot2` and `ggrepel` packages.",
    "",
    "## Contents",
    "",
    "- `Figures/<method>/`: six CellCycle PDF figures per method.",
    "- `stats/<method>/`: thirteen statistical CSV files per method.",
    "- `plot_data/<method>/`: one plotting-data CSV for each figure.",
    "",
    "## Methods",
    "",
    "- `initial_ploidy`: original 2N/4N sample grouping.",
    "- `ETP_fixed_threshold_2_25`: sample mean ETP threshold 2.25.",
    "- `ETP_boundary_stress_threshold_2_375`: sample mean ETP threshold 2.375.",
    "- `ETP_reference_balanced_threshold_2_24`: sample mean ETP threshold 2.24.",
    "",
    "## Figure data",
    "",
    "Each plotting-data CSV contains the exact rows used by its PDF. Scatter-plot tables also contain the Pearson and Spearman statistics, permutation P values, permutation mode, and the annotation text printed in the figure.",
    "",
    "Figures can be regenerated without the cell-level inputs and without rerunning statistical tests by using `--figures_only=TRUE --tables_root=<existing-output-root>`. In this mode the workflow reads only `plot_data/` and the required files in `stats/`, and writes only `Figures/`.",
    "",
    "The NonCellCycle input is used only when deriving sample mean end-timepoint ploidy. All figures and reported associations use CellCycle cells."
  )
  writeLines(lines, file.path(output_root, "README.md"))
  invisible(file.path(output_root, "README.md"))
}

essential_run_method <- function(cellcycle_path, noncellcycle_path, output_root, spec, seed, n_perm, n_boot) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required", call. = FALSE)
  if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required", call. = FALSE)
  set.seed(seed)
  figures_dir <- essential_ensure_dir(file.path(output_root, "Figures", spec$method))
  stats_dir <- essential_ensure_dir(file.path(output_root, "stats", spec$method))
  plot_data_dir <- essential_ensure_dir(file.path(output_root, "plot_data", spec$method))

  cellcycle <- essential_read_input(cellcycle_path, "CellCycle")
  noncellcycle <- essential_read_input(noncellcycle_path, "NonCellCycle")
  assignments <- essential_derive_assignments(cellcycle, noncellcycle)
  data <- essential_prepare_cellcycle(cellcycle, assignments, spec)
  sample_meta <- essential_build_sample_meta(data)
  equal_shift <- essential_calculate_shift_metrics(data, sample_meta, pooled_reference = FALSE)
  pooled_shift <- essential_calculate_shift_metrics(data, sample_meta, pooled_reference = TRUE)
  dose_tests <- essential_run_dose_tests(data, n_perm)

  association <- essential_run_tgi_associations(equal_shift, pooled_shift, spec$method, n_perm)
  essential_write_csv(association, file.path(stats_dir, "CellCycle_TGI_associations_ecdf_rmse.csv"))
  primary <- association[association$pre_specified_primary, , drop = FALSE]
  primary_pearson <- data.frame(
    n = primary$n, estimate = primary$pearson_r, asymptotic_p = primary$pearson_p_asymptotic,
    permutation_p_two_sided = primary$pearson_p_permutation_two_sided,
    permutation_p_positive = primary$pearson_p_permutation_positive,
    n_permutations = primary$pearson_n_permutations,
    permutation_mode = primary$pearson_permutation_mode
  )
  primary_spearman <- data.frame(
    n = primary$n, estimate = primary$spearman_rho, asymptotic_p = primary$spearman_p_asymptotic,
    permutation_p_two_sided = primary$spearman_p_permutation_two_sided,
    permutation_p_positive = primary$spearman_p_permutation_positive,
    n_permutations = primary$spearman_n_permutations,
    permutation_mode = primary$spearman_permutation_mode
  )
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  treated$analysis_group <- factor(treated$analysis_group, levels = spec$group_levels)
  treated$dose <- factor(treated$dose, levels = unique(treated$dose[order(treated$dose_mg)]))
  primary_plot_data <- essential_attach_association(
    treated,
    "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref",
    primary_pearson,
    primary_spearman
  )
  essential_write_csv(
    primary_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref_plot_data.csv")
  )
  primary_plot <- essential_scatter_plot(
    primary_plot_data,
    "ecdf_rmse",
    "TGI_percent_auc",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "Cell-cycle-associated tumor cells: AUC TGI vs sample-equal ECDF shift",
    "ECDF RMSE from group-matched equal-sample 0 mg/kg reference",
    "AUC-based TGI (%)",
    "Dose",
    spec$group_label,
    primary_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    primary_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_ecdf_rmse_equal_sample_ref.pdf"),
    6.8,
    6.8
  )

  direct <- essential_direct_ecdf(data, dose_tests, spec, n_perm)
  essential_write_csv(
    direct$tests,
    file.path(stats_dir, "CellCycle_direct_group_ecdf_comparisons_9panel_tests.csv")
  )
  essential_write_csv(
    direct$plot_data,
    file.path(plot_data_dir, "CellCycle_direct_group_ecdf_comparisons_plot_data.csv")
  )
  essential_save_pdf(
    direct$plot,
    file.path(figures_dir, "CellCycle_direct_group_ecdf_comparisons.pdf"),
    15,
    11.5
  )

  set.seed(seed + 101L)
  mean_association <- essential_run_mean_etp_associations(equal_shift, spec$method, n_perm)
  essential_write_csv(mean_association, file.path(stats_dir, "CellCycle_TGI_associations_mean_ETP.csv"))
  mean_primary <- mean_association[mean_association$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  mean_pearson <- data.frame(
    n = mean_primary$n, estimate = mean_primary$pearson_r, asymptotic_p = mean_primary$pearson_p_asymptotic,
    permutation_p_two_sided = mean_primary$pearson_p_permutation_two_sided,
    permutation_p_positive = mean_primary$pearson_p_permutation_positive,
    n_permutations = mean_primary$pearson_n_permutations,
    permutation_mode = mean_primary$pearson_permutation_mode
  )
  mean_spearman <- data.frame(
    n = mean_primary$n, estimate = mean_primary$spearman_rho, asymptotic_p = mean_primary$spearman_p_asymptotic,
    permutation_p_two_sided = mean_primary$spearman_p_permutation_two_sided,
    permutation_p_positive = mean_primary$spearman_p_permutation_positive,
    n_permutations = mean_primary$spearman_n_permutations,
    permutation_mode = mean_primary$spearman_permutation_mode
  )
  mean_plot_data <- essential_attach_association(
    treated,
    "CellCycle_TGI_AUC_vs_mean_ETP",
    mean_pearson,
    mean_spearman
  )
  essential_write_csv(
    mean_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_AUC_vs_mean_ETP_plot_data.csv")
  )
  mean_plot <- essential_scatter_plot(
    mean_plot_data,
    "sample_mean_endpoint_ploidy",
    "TGI_percent_auc",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "Cell-cycle-associated tumor cells: AUC TGI vs sample mean ETP",
    "Sample mean ETP",
    "AUC-based TGI (%)",
    "Dose",
    spec$group_label,
    mean_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    mean_plot,
    file.path(figures_dir, "CellCycle_TGI_AUC_vs_mean_ETP.pdf"),
    6.8,
    6.8
  )

  set.seed(seed + 201L)
  primary_loo <- essential_leave_one_out(treated, "ecdf_rmse")
  primary_boot <- essential_bootstrap(treated, "ecdf_rmse", n_boot)
  primary_summary <- essential_primary_robustness_summary(association, primary_loo, primary_boot)
  essential_write_csv(primary_loo, file.path(stats_dir, "CellCycle_primary_TGI_leave_one_out.csv"))
  essential_write_csv(primary_boot, file.path(stats_dir, "CellCycle_primary_TGI_bootstrap.csv"))
  essential_write_csv(primary_summary, file.path(stats_dir, "CellCycle_primary_TGI_robustness_summary.csv"))

  set.seed(seed + 301L)
  mean_loo <- essential_leave_one_out(treated, "sample_mean_endpoint_ploidy")
  mean_boot <- essential_bootstrap(treated, "sample_mean_endpoint_ploidy", n_boot)
  mean_summary <- essential_mean_etp_robustness_summary(
    spec$method, mean_association, mean_loo, mean_boot, n_boot
  )
  essential_write_csv(mean_loo, file.path(stats_dir, "CellCycle_mean_ETP_TGI_leave_one_out.csv"))
  essential_write_csv(mean_boot, file.path(stats_dir, "CellCycle_mean_ETP_TGI_bootstrap.csv"))
  essential_write_csv(mean_summary, file.path(stats_dir, "CellCycle_mean_ETP_TGI_robustness_summary.csv"))

  model_results <- essential_fit_shift_models(treated, spec)
  essential_write_csv(
    model_results$models,
    file.path(stats_dir, "CellCycle_AUC_TGI_shift_ploidy_dose_models.csv")
  )
  essential_write_csv(
    model_results$summaries,
    file.path(stats_dir, "CellCycle_AUC_TGI_shift_ploidy_dose_model_summaries.csv")
  )
  model_plot_data <- primary_plot_data
  model_plot_data$analysis_id <- "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose"
  essential_write_csv(
    model_plot_data,
    file.path(plot_data_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose_plot_data.csv")
  )
  model_plot <- essential_scatter_plot(
    model_plot_data,
    "ecdf_rmse",
    "TGI_percent_auc",
    "analysis_group",
    "dose",
    essential_group_colors(spec),
    "CellCycle: AUC-TGI association with pseudotime shift and ploidy",
    "ECDF RMSE from group-matched 0 mg/kg reference",
    "AUC-based TGI (%)",
    spec$group_label,
    "Dose",
    model_plot_data$plot_annotation[[1L]]
  )
  essential_save_pdf(
    model_plot,
    file.path(figures_dir, "CellCycle_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose.pdf"),
    6.8,
    6.8
  )

  set.seed(seed + 401L)
  confounding <- essential_confounding(equal_shift, spec, n_perm)
  essential_write_csv(confounding$ploidy_tests, file.path(stats_dir, "ploidy_confounding_tests.csv"))
  essential_write_csv(confounding$residualized, file.path(stats_dir, "residualized_TGI_associations.csv"))

  centered_result <- confounding$residualized[confounding$residualized$analysis == "within_dose_centered", , drop = FALSE]
  centered_pearson <- data.frame(
    n = centered_result$n,
    estimate = centered_result$estimate,
    asymptotic_p = centered_result$p_value,
    permutation_p_two_sided = centered_result$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = centered_result$n_permutations,
    permutation_mode = centered_result$permutation_mode
  )
  centered_spearman <- essential_permutation_cor(
    confounding$centered$shift_centered,
    confounding$centered$tgi_centered,
    "spearman",
    n_perm,
    exact = TRUE,
    strata = confounding$centered$dose_mg
  )
  centered_plot_data <- essential_attach_association(
    confounding$centered,
    "CellCycle_TGI_association_within_dose_centered",
    centered_pearson,
    centered_spearman
  )
  centered_plot_data$dose <- factor(
    centered_plot_data$dose,
    levels = unique(centered_plot_data$dose[order(centered_plot_data$dose_mg)])
  )
  centered_plot_data$analysis_group <- factor(
    centered_plot_data$analysis_group,
    levels = spec$group_levels
  )
  essential_write_csv(
    centered_plot_data,
    file.path(plot_data_dir, "CellCycle_TGI_association_within_dose_centered_plot_data.csv")
  )
  centered_plot <- essential_scatter_plot(
    centered_plot_data,
    "shift_centered",
    "tgi_centered",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle TGI association after within-dose centering",
    "Dose-centered ECDF RMSE",
    "Dose-centered AUC TGI",
    "Dose",
    spec$group_label,
    centered_plot_data$plot_annotation[[1L]],
    label_nudge_y = 1.0
  ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35)
  essential_save_pdf(
    centered_plot,
    file.path(figures_dir, "CellCycle_TGI_association_within_dose_centered.pdf"),
    6.6,
    6.6
  )

  ploidy_primary <- confounding$ploidy_tests[
    confounding$ploidy_tests$sample_set == "all" &
      confounding$ploidy_tests$shift_metric == "ecdf_rmse" &
      confounding$ploidy_tests$ploidy_measure == "mean_cell_ploidy",
    ,
    drop = FALSE
  ]
  pearson_row <- ploidy_primary[ploidy_primary$method == "pearson", , drop = FALSE]
  spearman_row <- ploidy_primary[ploidy_primary$method == "spearman", , drop = FALSE]
  ploidy_pearson <- data.frame(
    n = pearson_row$n,
    estimate = pearson_row$estimate,
    asymptotic_p = pearson_row$p_value,
    permutation_p_two_sided = pearson_row$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = pearson_row$n_permutations,
    permutation_mode = pearson_row$permutation_mode
  )
  ploidy_spearman <- data.frame(
    n = spearman_row$n,
    estimate = spearman_row$estimate,
    asymptotic_p = spearman_row$p_value,
    permutation_p_two_sided = spearman_row$p_permutation_two_sided,
    permutation_p_positive = NA_real_,
    n_permutations = spearman_row$n_permutations,
    permutation_mode = spearman_row$permutation_mode
  )
  ploidy_source <- equal_shift[
    is.finite(equal_shift$mean_cell_ploidy) & is.finite(equal_shift$ecdf_rmse),
    ,
    drop = FALSE
  ]
  ploidy_source$dose <- factor(
    ploidy_source$dose,
    levels = unique(ploidy_source$dose[order(ploidy_source$dose_mg)])
  )
  ploidy_source$analysis_group <- factor(ploidy_source$analysis_group, levels = spec$group_levels)
  ploidy_plot_data <- essential_attach_association(
    ploidy_source,
    "CellCycle_ecdf_rmse_vs_ploidy",
    ploidy_pearson,
    ploidy_spearman
  )
  essential_write_csv(
    ploidy_plot_data,
    file.path(plot_data_dir, "CellCycle_ecdf_rmse_vs_ploidy_plot_data.csv")
  )
  ploidy_plot <- essential_scatter_plot(
    ploidy_plot_data,
    "mean_cell_ploidy",
    "ecdf_rmse",
    "dose",
    "analysis_group",
    essential_dose_colors(),
    "CellCycle ECDF RMSE vs mean cell ploidy",
    "Mean cell ploidy",
    "ECDF RMSE",
    "Dose",
    spec$group_label,
    ploidy_plot_data$plot_annotation[[1L]],
    label_nudge_y = 0.004
  )
  essential_save_pdf(
    ploidy_plot,
    file.path(figures_dir, "CellCycle_ecdf_rmse_vs_ploidy.pdf"),
    6.4,
    6.4
  )

  invisible(list(method = spec$method, figures = 6L, stats = 13L, plot_data = 6L))
}

essential_validate_inventory <- function(output_root, methods) {
  os_metadata <- list.files(
    output_root,
    pattern = "^[.]DS_Store$",
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE
  )
  if (length(os_metadata) > 0L) unlink(os_metadata, force = TRUE)
  expected_figures <- 6L * length(methods)
  expected_stats <- 13L * length(methods)
  expected_plot_data <- 6L * length(methods)
  figures <- list.files(file.path(output_root, "Figures"), pattern = "[.]pdf$", recursive = TRUE, full.names = TRUE)
  stats <- list.files(file.path(output_root, "stats"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
  plot_data <- list.files(file.path(output_root, "plot_data"), pattern = "[.]csv$", recursive = TRUE, full.names = TRUE)
  if (length(figures) != expected_figures) {
    stop("Unexpected PDF count: ", length(figures), "; expected ", expected_figures, call. = FALSE)
  }
  if (length(stats) != expected_stats) {
    stop("Unexpected statistics CSV count: ", length(stats), "; expected ", expected_stats, call. = FALSE)
  }
  if (length(plot_data) != expected_plot_data) {
    stop("Unexpected plot-data CSV count: ", length(plot_data), "; expected ", expected_plot_data, call. = FALSE)
  }
  unexpected <- list.files(
    output_root,
    pattern = "[.](png|ppt|pptx)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(unexpected) > 0L) {
    stop("Unexpected output files: ", paste(unexpected, collapse = ", "), call. = FALSE)
  }
  all_files <- list.files(
    output_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE
  )
  allowed <- grepl("[.](pdf|csv)$", all_files, ignore.case = TRUE) |
    basename(all_files) == "README.md"
  if (any(!allowed)) {
    stop("Unexpected essential output: ", paste(all_files[!allowed], collapse = ", "), call. = FALSE)
  }
  invisible(list(figures = length(figures), stats = length(stats), plot_data = length(plot_data)))
}
