#!/usr/bin/env Rscript

required_packages <- c("drc", "ggplot2", "scales")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(ggplot2)
})

command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("Code/in-vitro/drug_response/plot_gemcitabine_ploidy_auc_association.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    input_file = NULL,
    output_dir = NULL,
    ploidy_file = NULL,
    used_named_output_dir = FALSE,
    positional = character()
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (startsWith(arg, "--input=")) {
      out$input_file <- sub("^--input=", "", arg)
    } else if (arg == "--input") {
      i <- i + 1
      out$input_file <- args[[i]]
    } else if (startsWith(arg, "--output-dir=")) {
      out$output_dir <- sub("^--output-dir=", "", arg)
      out$used_named_output_dir <- TRUE
    } else if (arg == "--output-dir") {
      i <- i + 1
      out$output_dir <- args[[i]]
      out$used_named_output_dir <- TRUE
    } else if (startsWith(arg, "--ploidy-file=")) {
      out$ploidy_file <- sub("^--ploidy-file=", "", arg)
    } else if (arg == "--ploidy-file") {
      i <- i + 1
      out$ploidy_file <- args[[i]]
    } else if (startsWith(arg, "--")) {
      stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
    } else {
      out$positional <- c(out$positional, arg)
    }
    i <- i + 1
  }
  if (is.null(out$input_file) && length(out$positional) >= 1) {
    out$input_file <- out$positional[[1]]
  }
  if (is.null(out$output_dir) && length(out$positional) >= 2) {
    out$output_dir <- out$positional[[2]]
  }
  if (is.null(out$ploidy_file) && length(out$positional) >= 3) {
    out$ploidy_file <- out$positional[[3]]
  }
  out
}

parsed_args <- parse_args(args)
input_file <- if (!is.null(parsed_args$input_file)) {
  parsed_args$input_file
} else {
  file.path(repo_root, "Data/in-vitro/drug_response/Gemcitabine.txt")
}
out_dir <- if (!is.null(parsed_args$output_dir)) {
  parsed_args$output_dir
} else {
  file.path(repo_root, "Figs")
}
ploidy_file <- if (!is.null(parsed_args$ploidy_file)) {
  parsed_args$ploidy_file
} else {
  file.path(repo_root, "Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv")
}
fig_dir <- if (isTRUE(parsed_args$used_named_output_dir)) file.path(out_dir, "figures") else out_dir
table_dir <- if (isTRUE(parsed_args$used_named_output_dir)) file.path(out_dir, "tables") else out_dir
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

common_auc_min_uM <- 0.005
common_auc_max_uM <- 0.9
common_auc_label <- sprintf("%.3g-%.3g uM", common_auc_min_uM, common_auc_max_uM)
common_auc_id <- "0p005_0p9_uM"

# Gemcitabine.txt contains only response columns; low/high pair labels and ploidy
# values are carried here from the original Figure 3H analysis metadata.
sample_metadata <- data.frame(
  group = c("HGC", "HGC", "MDA", "MDA", "SNU", "SNU", "SUM", "SUM",
            "SUMOld", "SUMOld", "MDAOld", "MDAOld"),
  ploidy_class = rep(c("Low", "High"), 6),
  sample = c(
    "HGC-27_MRCA_harvest",
    "HGC-27_G1_A10_harvesT3",
    "MDA-MB-231_O1_A23-K_harvest",
    "MDA-MB-231_G2_A3_harvesT8",
    "SNU-668_P1_A19kT_harvest",
    "SNU-668_r2_GFP_VT_A7_harvesT3",
    "SUM159_NLS_2N_O1_A7K_harvest",
    "SUM-159_NLS_2N_O2_A24_seedT2",
    "SUM-159_2N",
    "SUM-159_4N",
    "MDAMB231_2N",
    "MDAMB231_4N"
  ),
  embedded_ploidy = c(
    2.594,
    2.814,
    3.430,
    4.249,
    2.462,
    4.726,
    1.996,
    3.726,
    2.009975,
    3.49593528659252,
    4.15675859608933,
    6.5
  ),
  stringsAsFactors = FALSE
)

sample_metadata_base <- sample_metadata
sample_order <- sample_metadata_base$sample

ploidy_stat_specs <- data.frame(
  ploidy_stat = c("mean", "median", "q90", "q10", "q75", "max", "min"),
  ploidy_label = c("mean", "median", "90th percentile", "10th percentile", "75th percentile", "maximum", "minimum"),
  ploidy_column = c(
    "mean_population_ploidy",
    "median_population_ploidy",
    "q90_population_ploidy",
    "q10_population_ploidy",
    "q75_population_ploidy",
    "max_population_ploidy",
    "min_population_ploidy"
  ),
  output_suffix = c("mean", "median", "q90", "q10", "q75", "max", "min"),
  stringsAsFactors = FALSE
)

load_cloneid_ploidy <- function(path, stat_specs) {
  if (!file.exists(path)) {
    return(NULL)
  }

  cloneid_ploidy <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_ploidy_columns <- unique(c(
    "dose_response_label",
    stat_specs$ploidy_column,
    "analysis_cloneid_origin",
    "origin_resolution_status",
    "origin_resolution_distance_edges",
    "origin_resolution_mrca",
    "ploidy_status"
  ))
  missing_ploidy_columns <- setdiff(required_ploidy_columns, names(cloneid_ploidy))
  if (length(missing_ploidy_columns) > 0) {
    stop(
      "CLONEID ploidy table is missing columns: ",
      paste(missing_ploidy_columns, collapse = ", "),
      ". Regenerate it with query_cloneid_fig3h_ploidy.R.",
      call. = FALSE
    )
  }
  cloneid_ploidy[, required_ploidy_columns, drop = FALSE]
}

build_sample_metadata_for_ploidy_stat <- function(base_metadata, cloneid_ploidy, stat_spec) {
  sample_metadata <- base_metadata
  sample_metadata$ploidy <- sample_metadata$embedded_ploidy
  sample_metadata$ploidy_stat <- stat_spec$ploidy_stat
  sample_metadata$ploidy_stat_label <- stat_spec$ploidy_label
  sample_metadata$ploidy_source <- "embedded_fallback"
  sample_metadata$analysis_cloneid_origin <- NA_character_
  sample_metadata$origin_resolution_status <- NA_character_
  sample_metadata$origin_resolution_distance_edges <- NA_real_
  sample_metadata$origin_resolution_mrca <- NA_character_
  sample_metadata$ploidy_status <- NA_character_

  if (!is.null(cloneid_ploidy)) {
    stat_column <- stat_spec$ploidy_column
    cloneid_stat <- cloneid_ploidy[, c(
      "dose_response_label",
      stat_column,
      "analysis_cloneid_origin",
      "origin_resolution_status",
      "origin_resolution_distance_edges",
      "origin_resolution_mrca",
      "ploidy_status"
    ), drop = FALSE]
    names(cloneid_stat)[names(cloneid_stat) == stat_column] <- "cloneid_ploidy"

    sample_metadata <- merge(
      sample_metadata,
      cloneid_stat,
      by.x = "sample",
      by.y = "dose_response_label",
      all.x = TRUE,
      sort = FALSE,
      suffixes = c("", "_cloneid")
    )
    sample_metadata <- sample_metadata[match(base_metadata$sample, sample_metadata$sample), , drop = FALSE]

    use_cloneid_ploidy <- is.finite(sample_metadata$cloneid_ploidy) &
      !is.na(sample_metadata$ploidy_status_cloneid) &
      sample_metadata$ploidy_status_cloneid == "ok"
    sample_metadata$ploidy[use_cloneid_ploidy] <- sample_metadata$cloneid_ploidy[use_cloneid_ploidy]
    sample_metadata$ploidy_source[use_cloneid_ploidy] <- "cloneid_genome_perspective"
    sample_metadata$analysis_cloneid_origin <- sample_metadata$analysis_cloneid_origin_cloneid
    sample_metadata$origin_resolution_status <- sample_metadata$origin_resolution_status_cloneid
    sample_metadata$origin_resolution_distance_edges <- sample_metadata$origin_resolution_distance_edges_cloneid
    sample_metadata$origin_resolution_mrca <- sample_metadata$origin_resolution_mrca_cloneid
    sample_metadata$ploidy_status <- sample_metadata$ploidy_status_cloneid
    sample_metadata$cloneid_ploidy <- NULL
    sample_metadata$analysis_cloneid_origin_cloneid <- NULL
    sample_metadata$origin_resolution_status_cloneid <- NULL
    sample_metadata$origin_resolution_distance_edges_cloneid <- NULL
    sample_metadata$origin_resolution_mrca_cloneid <- NULL
    sample_metadata$ploidy_status_cloneid <- NULL
  }

  sample_metadata
}

read_gemcitabine_table <- function(path) {
  raw <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  if (ncol(raw) < 2) {
    stop("Expected at least one dose column and one sample column in ", path, call. = FALSE)
  }

  dose_uM <- suppressWarnings(as.numeric(raw[[1]]))
  sample_headers <- names(raw)[-1]

  long <- do.call(rbind, lapply(seq_along(sample_headers), function(index) {
    data.frame(
      dose_uM = dose_uM,
      sample = sample_headers[[index]],
      replicate_column = index,
      response = suppressWarnings(as.numeric(raw[[index + 1]])),
      stringsAsFactors = FALSE
    )
  }))
  long[is.finite(long$dose_uM) & is.finite(long$response), , drop = FALSE]
}

trapezoid_auc <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  unique_x <- !duplicated(x)
  x <- x[unique_x]
  y <- y[unique_x]
  if (length(x) < 2) {
    return(NA_real_)
  }
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

trapezoid_auc_over_dose_interval <- function(dose_uM, y, lower_uM, upper_uM) {
  keep <- is.finite(dose_uM) & dose_uM > 0 & is.finite(y)
  dose_uM <- dose_uM[keep]
  y <- y[keep]
  if (length(dose_uM) < 2 || lower_uM >= upper_uM) {
    return(NA_real_)
  }

  log_dose <- log10(dose_uM)
  ord <- order(log_dose)
  log_dose <- log_dose[ord]
  y <- y[ord]
  unique_log_dose <- !duplicated(log_dose)
  log_dose <- log_dose[unique_log_dose]
  y <- y[unique_log_dose]

  log_lower <- log10(lower_uM)
  log_upper <- log10(upper_uM)
  if (log_lower < min(log_dose) || log_upper > max(log_dose)) {
    return(NA_real_)
  }

  auc_grid <- sort(unique(c(
    log_lower,
    log_dose[log_dose > log_lower & log_dose < log_upper],
    log_upper
  )))
  auc_y <- approx(log_dose, y, xout = auc_grid, ties = mean)$y
  trapezoid_auc(auc_grid, auc_y)
}

safe_log10 <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x) & x > 0
  out[ok] <- log10(x[ok])
  out
}

estimate_fitted_response_crossing <- function(fit, target_response, lower_uM, upper_uM) {
  if (inherits(fit, "try-error") || lower_uM <= 0 || upper_uM <= lower_uM) {
    return(NA_real_)
  }

  objective <- function(log_dose) {
    as.numeric(predict(fit, newdata = data.frame(dose_uM = 10^log_dose))) - target_response
  }
  lower_log <- log10(lower_uM)
  upper_log <- log10(upper_uM)
  lower_value <- try(objective(lower_log), silent = TRUE)
  upper_value <- try(objective(upper_log), silent = TRUE)
  if (inherits(lower_value, "try-error") || inherits(upper_value, "try-error") ||
      !is.finite(lower_value) || !is.finite(upper_value)) {
    return(NA_real_)
  }
  if (lower_value == 0) {
    return(lower_uM)
  }
  if (upper_value == 0) {
    return(upper_uM)
  }
  if (sign(lower_value) == sign(upper_value)) {
    return(NA_real_)
  }

  root <- try(stats::uniroot(objective, lower = lower_log, upper = upper_log)$root, silent = TRUE)
  if (inherits(root, "try-error") || !is.finite(root)) {
    return(NA_real_)
  }
  10^root
}

estimate_observed_response_crossing <- function(dose_uM, y, target_response, lower_uM, upper_uM) {
  keep <- is.finite(dose_uM) & dose_uM > 0 & is.finite(y) & dose_uM >= lower_uM & dose_uM <= upper_uM
  dose_uM <- dose_uM[keep]
  y <- y[keep]
  if (length(dose_uM) < 2) {
    return(NA_real_)
  }

  log_dose <- log10(dose_uM)
  ord <- order(log_dose)
  log_dose <- log_dose[ord]
  y <- y[ord]
  for (i in seq_len(length(log_dose) - 1L)) {
    y0 <- y[[i]] - target_response
    y1 <- y[[i + 1L]] - target_response
    if (y0 == 0) {
      return(10^log_dose[[i]])
    }
    if (y1 == 0) {
      return(10^log_dose[[i + 1L]])
    }
    if (sign(y0) != sign(y1)) {
      fraction <- abs(y0) / (abs(y0) + abs(y1))
      return(10^(log_dose[[i]] + fraction * (log_dose[[i + 1L]] - log_dose[[i]])))
    }
  }
  NA_real_
}

fit_sample_response <- function(sample_name, response_long) {
  sample_data <- response_long[response_long$sample == sample_name, , drop = FALSE]
  if (nrow(sample_data) == 0) {
    stop("No response data found for sample: ", sample_name, call. = FALSE)
  }

  mean_response <- aggregate(response ~ dose_uM, sample_data, mean)
  control_mean <- mean(mean_response$response[mean_response$dose_uM == 0], na.rm = TRUE)
  if (!is.finite(control_mean) || control_mean <= 0) {
    stop("Missing or invalid untreated control mean for sample: ", sample_name, call. = FALSE)
  }

  mean_response$response_norm <- mean_response$response / control_mean
  fit_data <- mean_response[mean_response$dose_uM > 0 & is.finite(mean_response$response_norm), , drop = FALSE]
  if (nrow(fit_data) < 4) {
    stop("Too few positive-dose observations to fit sample: ", sample_name, call. = FALSE)
  }

  dose_min <- min(fit_data$dose_uM)
  dose_max <- max(fit_data$dose_uM)
  log_dose_min <- log10(dose_min)
  log_dose_max <- log10(dose_max)
  if (common_auc_min_uM < dose_min || common_auc_max_uM > dose_max) {
    stop(
      sprintf(
        "Common AUC interval %s is outside the fitted positive-dose range %.4g-%.4g uM for sample %s.",
        common_auc_label,
        dose_min,
        dose_max,
        sample_name
      ),
      call. = FALSE
    )
  }

  fit <- try(
    suppressWarnings(drc::drm(
      response_norm ~ dose_uM,
      data = fit_data,
      fct = drc::LL.4(names = c("Slope", "LowerLimit", "UpperLimit", "IC50")),
      control = drc::drmc(noMessage = TRUE)
    )),
    silent = TRUE
  )

  model_type <- "LL.4"
  coefficients <- rep(NA_real_, 4)
  names(coefficients) <- c("slope", "lower_limit", "upper_limit", "ic50_uM")
  ec50_uM <- NA_real_
  ic50_uM <- NA_real_
  auc_norm_log10dose_full_observed_range <- NA_real_
  auc_norm_log10dose_common_range <- NA_real_

  if (!inherits(fit, "try-error")) {
    coefficients[] <- as.numeric(coef(fit)[seq_len(4)])
    ec50_uM <- coefficients[["ic50_uM"]]
    ic50_uM <- estimate_fitted_response_crossing(
      fit,
      target_response = 0.5,
      lower_uM = dose_min,
      upper_uM = dose_max
    )
    predicted_auc_full <- try(
      integrate(
        function(log_dose) {
          as.numeric(predict(fit, newdata = data.frame(dose_uM = 10^log_dose)))
        },
        lower = log_dose_min,
        upper = log_dose_max,
        subdivisions = 200,
        rel.tol = 1e-6
      )$value,
      silent = TRUE
    )
    if (!inherits(predicted_auc_full, "try-error") && is.finite(predicted_auc_full)) {
      auc_norm_log10dose_full_observed_range <- as.numeric(predicted_auc_full)
    }

    predicted_auc_common <- try(
      integrate(
        function(log_dose) {
          as.numeric(predict(fit, newdata = data.frame(dose_uM = 10^log_dose)))
        },
        lower = log10(common_auc_min_uM),
        upper = log10(common_auc_max_uM),
        subdivisions = 200,
        rel.tol = 1e-6
      )$value,
      silent = TRUE
    )
    if (!inherits(predicted_auc_common, "try-error") && is.finite(predicted_auc_common)) {
      auc_norm_log10dose_common_range <- as.numeric(predicted_auc_common)
    }
  }

  if (!is.finite(auc_norm_log10dose_full_observed_range) || !is.finite(auc_norm_log10dose_common_range)) {
    model_type <- "observed_trapezoid"
    auc_norm_log10dose_full_observed_range <- trapezoid_auc(log10(fit_data$dose_uM), fit_data$response_norm)
    auc_norm_log10dose_common_range <- trapezoid_auc_over_dose_interval(
      fit_data$dose_uM,
      fit_data$response_norm,
      common_auc_min_uM,
      common_auc_max_uM
    )
  }
  if (!is.finite(ic50_uM)) {
    ic50_uM <- estimate_observed_response_crossing(
      fit_data$dose_uM,
      fit_data$response_norm,
      target_response = 0.5,
      lower_uM = dose_min,
      upper_uM = dose_max
    )
  }

  if (!is.finite(auc_norm_log10dose_common_range)) {
    stop("Could not compute common-interval AUC for sample: ", sample_name, call. = FALSE)
  }

  data.frame(
    sample = sample_name,
    model_type = model_type,
    control_mean = control_mean,
    slope = coefficients[["slope"]],
    lower_limit = coefficients[["lower_limit"]],
    upper_limit = coefficients[["upper_limit"]],
    ec50_uM = ec50_uM,
    ic50_uM = ic50_uM,
    log10_ec50_uM = safe_log10(ec50_uM),
    log10_ic50_uM = safe_log10(ic50_uM),
    auc_metric = paste0("normalized_response_auc_log10dose_common_", common_auc_id),
    auc_norm_log10dose = auc_norm_log10dose_common_range,
    auc_norm_log10dose_common_range = auc_norm_log10dose_common_range,
    auc_norm_log10dose_full_observed_range = auc_norm_log10dose_full_observed_range,
    auc_integration_min_uM = common_auc_min_uM,
    auc_integration_max_uM = common_auc_max_uM,
    full_auc_integration_min_uM = dose_min,
    full_auc_integration_max_uM = dose_max,
    positive_dose_min_uM = dose_min,
    positive_dose_max_uM = dose_max,
    n_positive_doses = length(unique(fit_data$dose_uM)),
    n_observations = nrow(sample_data),
    stringsAsFactors = FALSE
  )
}

build_normalized_fit_plot_data <- function(sample_name, response_long) {
  sample_data <- response_long[response_long$sample == sample_name, , drop = FALSE]
  mean_response <- aggregate(response ~ dose_uM, sample_data, mean)
  control_mean <- mean(mean_response$response[mean_response$dose_uM == 0], na.rm = TRUE)
  if (!is.finite(control_mean) || control_mean <= 0) {
    stop("Missing or invalid untreated control mean for sample: ", sample_name, call. = FALSE)
  }

  raw_points <- sample_data
  raw_points$response_norm <- raw_points$response / control_mean

  fit_data <- mean_response[mean_response$dose_uM > 0 & is.finite(mean_response$response), , drop = FALSE]
  fit_data$response_norm <- fit_data$response / control_mean
  fit <- try(
    suppressWarnings(drc::drm(
      response_norm ~ dose_uM,
      data = fit_data,
      fct = drc::LL.4(names = c("Slope", "LowerLimit", "UpperLimit", "IC50")),
      control = drc::drmc(noMessage = TRUE)
    )),
    silent = TRUE
  )

  predicted_curve <- data.frame()
  if (!inherits(fit, "try-error")) {
    dose_grid <- 10^seq(log10(min(fit_data$dose_uM)), log10(max(fit_data$dose_uM)), length.out = 250)
    predicted_curve <- data.frame(
      sample = sample_name,
      dose_uM = dose_grid,
      response_norm = as.numeric(predict(fit, newdata = data.frame(dose_uM = dose_grid))),
      stringsAsFactors = FALSE
    )
  }

  list(raw_points = raw_points, predicted_curve = predicted_curve)
}

format_p <- function(p_value) {
  if (!is.finite(p_value)) {
    return("NA")
  }
  if (p_value < 0.001) {
    return("<0.001")
  }
  sprintf("%.3f", p_value)
}

correlation_table <- function(x, y, count_name) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n >= 3 && length(unique(x[ok])) > 1 && length(unique(y[ok])) > 1) {
    pearson <- cor.test(x[ok], y[ok], method = "pearson")
    r <- unname(pearson$estimate)
    p <- pearson$p.value
  } else {
    r <- NA_real_
    p <- NA_real_
  }
  out <- data.frame(
    drug = "Gemcitabine",
    pearson_r = r,
    p_value = p,
    stringsAsFactors = FALSE
  )
  out[[count_name]] <- n
  out[, c("drug", count_name, "pearson_r", "p_value"), drop = FALSE]
}

format_ploidy_for_legend <- function(x) {
  if (!is.finite(x)) {
    return("NA")
  }
  sprintf("%.1fN", x)
}

build_group_pair_labels <- function(fit_results) {
  groups <- unique(as.character(fit_results$group))
  labels <- vapply(groups, function(group_name) {
    rows <- fit_results[fit_results$group == group_name, , drop = FALSE]
    low <- rows[rows$ploidy_class == "Low", , drop = FALSE]
    high <- rows[rows$ploidy_class == "High", , drop = FALSE]
    cell_line <- unname(group_display_names[group_name])
    if (is.na(cell_line) || !nzchar(cell_line)) {
      cell_line <- group_name
    }
    sprintf(
      "%s: %s vs. %s",
      cell_line,
      format_ploidy_for_legend(low$ploidy[[1]]),
      format_ploidy_for_legend(high$ploidy[[1]])
    )
  }, character(1))
  names(labels) <- groups
  labels
}

add_group_pair_labels <- function(fit_results, group_pair_labels) {
  fit_results$group_pair_label <- unname(group_pair_labels[as.character(fit_results$group)])
  fit_results
}

group_color_scale <- function(groups, group_pair_labels) {
  legend_breaks <- intersect(names(group_colors), unique(as.character(groups)))
  legend_labels <- group_pair_labels[legend_breaks]
  missing_labels <- is.na(legend_labels) | !nzchar(legend_labels)
  legend_labels[missing_labels] <- legend_breaks[missing_labels]
  scale_color_manual(
    values = group_colors,
    breaks = legend_breaks,
    labels = legend_labels,
    name = "Cell line (ploidy)"
  )
}

build_pair_stats <- function(fit_results, stat_spec) {
  pair_stats <- do.call(rbind, lapply(unique(fit_results$group), function(group_name) {
    rows <- fit_results[fit_results$group == group_name, , drop = FALSE]
    low <- rows[rows$ploidy_class == "Low", , drop = FALSE]
    high <- rows[rows$ploidy_class == "High", , drop = FALSE]
    data.frame(
      drug = "Gemcitabine",
      ploidy_stat = stat_spec$ploidy_stat,
      ploidy_stat_label = stat_spec$ploidy_label,
      group = group_name,
      group_pair_label = low$group_pair_label,
      low_sample = low$sample,
      high_sample = high$sample,
      low_ploidy = low$ploidy,
      high_ploidy = high$ploidy,
      low_ploidy_source = low$ploidy_source,
      high_ploidy_source = high$ploidy_source,
      low_analysis_cloneid_origin = low$analysis_cloneid_origin,
      high_analysis_cloneid_origin = high$analysis_cloneid_origin,
      low_origin_resolution_status = low$origin_resolution_status,
      high_origin_resolution_status = high$origin_resolution_status,
      low_origin_resolution_distance_edges = low$origin_resolution_distance_edges,
      high_origin_resolution_distance_edges = high$origin_resolution_distance_edges,
      low_origin_resolution_mrca = low$origin_resolution_mrca,
      high_origin_resolution_mrca = high$origin_resolution_mrca,
      delta_ploidy = high$ploidy - low$ploidy,
      auc_metric = low$auc_metric,
      auc_integration_min_uM = low$auc_integration_min_uM,
      auc_integration_max_uM = low$auc_integration_max_uM,
      low_auc = low$auc_norm_log10dose,
      high_auc = high$auc_norm_log10dose,
      delta_auc_high_minus_low = high$auc_norm_log10dose - low$auc_norm_log10dose,
      low_auc_full_observed_range = low$auc_norm_log10dose_full_observed_range,
      high_auc_full_observed_range = high$auc_norm_log10dose_full_observed_range,
      delta_auc_full_observed_range_high_minus_low =
        high$auc_norm_log10dose_full_observed_range - low$auc_norm_log10dose_full_observed_range,
      low_full_auc_integration_min_uM = low$full_auc_integration_min_uM,
      low_full_auc_integration_max_uM = low$full_auc_integration_max_uM,
      high_full_auc_integration_min_uM = high$full_auc_integration_min_uM,
      high_full_auc_integration_max_uM = high$full_auc_integration_max_uM,
      low_ec50_uM = low$ec50_uM,
      high_ec50_uM = high$ec50_uM,
      low_ic50_uM = low$ic50_uM,
      high_ic50_uM = high$ic50_uM,
      log2_fc_ec50_high_vs_low = log2(high$ec50_uM / low$ec50_uM),
      log2_fc_ic50_high_vs_low = log2(high$ic50_uM / low$ic50_uM),
      stringsAsFactors = FALSE
    )
  }))
  pair_stats
}

build_absolute_auc_stats <- function(fit_results, stat_spec) {
  absolute_auc_stats <- fit_results[, c(
    "sample", "group", "group_pair_label", "ploidy_class", "ploidy", "ploidy_source",
    "analysis_cloneid_origin", "origin_resolution_status",
    "origin_resolution_distance_edges", "origin_resolution_mrca",
    "auc_metric", "auc_norm_log10dose", "auc_norm_log10dose_common_range",
    "auc_norm_log10dose_full_observed_range", "auc_integration_min_uM",
    "auc_integration_max_uM", "full_auc_integration_min_uM",
    "full_auc_integration_max_uM", "ec50_uM", "ic50_uM",
    "log10_ec50_uM", "log10_ic50_uM", "model_type"
  ), drop = FALSE]
  absolute_auc_stats$drug <- "Gemcitabine"
  absolute_auc_stats$ploidy_stat <- stat_spec$ploidy_stat
  absolute_auc_stats$ploidy_stat_label <- stat_spec$ploidy_label
  absolute_auc_stats$ploidy <- as.numeric(absolute_auc_stats$ploidy)
  absolute_auc_stats$auc_norm_log10dose <- as.numeric(absolute_auc_stats$auc_norm_log10dose)
  absolute_auc_stats[, c(
    "drug", "ploidy_stat", "ploidy_stat_label",
    setdiff(names(absolute_auc_stats), c("drug", "ploidy_stat", "ploidy_stat_label"))
  ), drop = FALSE]
}

make_delta_plot <- function(pair_stats, correlation_summary, stat_label, group_pair_labels) {
  annotation <- data.frame(
    drug = "Gemcitabine",
    label = sprintf(
      "R = %.2f, p = %s",
      correlation_summary$pearson_r,
      format_p(correlation_summary$p_value)
    ),
    stringsAsFactors = FALSE
  )

  ggplot(
    pair_stats,
    aes(x = delta_ploidy, y = delta_auc_high_minus_low, color = group)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dotted", color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dotted", color = "black") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, linetype = "longdash", color = "black") +
    geom_point(size = 3.8) +
    geom_text(
      data = annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.35,
      size = 4.2
    ) +
    facet_wrap(~ drug) +
    group_color_scale(pair_stats$group, group_pair_labels) +
    expand_limits(x = 0, y = 0) +
    labs(
      title = sprintf("Relative: Delta AUC (High - Low), %s ploidy", stat_label),
      x = sprintf("Delta ploidy (%s)", stat_label),
      y = sprintf("Delta standardized AUC, %s (High - Low)", common_auc_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, hjust = 0, face = "plain"),
      strip.background = element_rect(fill = "grey85", color = "grey35", linewidth = 0.4),
      strip.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.3),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

make_absolute_plot <- function(absolute_auc_stats, absolute_correlation_summary, stat_label, group_pair_labels) {
  absolute_auc_ok <- is.finite(absolute_auc_stats$ploidy) & is.finite(absolute_auc_stats$auc_norm_log10dose)
  absolute_annotation <- data.frame(
    label = sprintf(
      "R = %.2f, p = %s",
      absolute_correlation_summary$pearson_r,
      format_p(absolute_correlation_summary$p_value)
    ),
    stringsAsFactors = FALSE
  )

  ggplot(
    absolute_auc_stats[absolute_auc_ok, , drop = FALSE],
    aes(x = ploidy, y = auc_norm_log10dose, color = group, shape = ploidy_class)
  ) +
    geom_smooth(
      data = absolute_auc_stats[absolute_auc_ok, , drop = FALSE],
      aes(x = ploidy, y = auc_norm_log10dose, group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = FALSE,
      linewidth = 0.9,
      linetype = "longdash",
      color = "black"
    ) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text(
      data = absolute_annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.35,
      size = 4.2
    ) +
    group_color_scale(absolute_auc_stats$group[absolute_auc_ok], group_pair_labels) +
    scale_shape_manual(values = c(Low = 16, High = 17), name = "Ploidy class") +
    labs(
      title = sprintf("Absolute %s ploidy versus normalized gemcitabine AUC", stat_label),
      x = sprintf("Ploidy (%s)", stat_label),
      y = sprintf("Standardized AUC of normalized response, %s", common_auc_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, hjust = 0, face = "plain"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.3),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
}

concentration_metric_specs <- data.frame(
  response_metric = c("ec50", "ic50"),
  response_metric_label = c("EC50", "IC50"),
  transformed_value_column = c("log10_ec50_uM", "log10_ic50_uM"),
  raw_value_column = c("ec50_uM", "ic50_uM"),
  delta_base = c(
    "gemcitabine_delta_ec50_vs_delta_ploidy",
    "gemcitabine_delta_ic50_vs_delta_ploidy"
  ),
  absolute_base = c(
    "gemcitabine_absolute_ec50_vs_ploidy",
    "gemcitabine_absolute_ic50_vs_ploidy"
  ),
  stringsAsFactors = FALSE
)

build_pair_concentration_metric_stats <- function(fit_results, stat_spec, metric_spec) {
  value_column <- metric_spec$transformed_value_column
  raw_column <- metric_spec$raw_value_column
  pair_stats <- do.call(rbind, lapply(unique(fit_results$group), function(group_name) {
    rows <- fit_results[fit_results$group == group_name, , drop = FALSE]
    low <- rows[rows$ploidy_class == "Low", , drop = FALSE]
    high <- rows[rows$ploidy_class == "High", , drop = FALSE]
    low_value <- low[[value_column]]
    high_value <- high[[value_column]]
    low_raw <- low[[raw_column]]
    high_raw <- high[[raw_column]]
    data.frame(
      drug = "Gemcitabine",
      response_metric = metric_spec$response_metric,
      response_metric_label = metric_spec$response_metric_label,
      response_metric_scale = "log10_uM",
      ploidy_stat = stat_spec$ploidy_stat,
      ploidy_stat_label = stat_spec$ploidy_label,
      group = group_name,
      group_pair_label = low$group_pair_label,
      low_sample = low$sample,
      high_sample = high$sample,
      low_ploidy = low$ploidy,
      high_ploidy = high$ploidy,
      low_ploidy_source = low$ploidy_source,
      high_ploidy_source = high$ploidy_source,
      low_analysis_cloneid_origin = low$analysis_cloneid_origin,
      high_analysis_cloneid_origin = high$analysis_cloneid_origin,
      low_origin_resolution_status = low$origin_resolution_status,
      high_origin_resolution_status = high$origin_resolution_status,
      low_origin_resolution_distance_edges = low$origin_resolution_distance_edges,
      high_origin_resolution_distance_edges = high$origin_resolution_distance_edges,
      low_origin_resolution_mrca = low$origin_resolution_mrca,
      high_origin_resolution_mrca = high$origin_resolution_mrca,
      delta_ploidy = high$ploidy - low$ploidy,
      low_metric_value = low_value,
      high_metric_value = high_value,
      delta_metric_high_minus_low = high_value - low_value,
      low_metric_uM = low_raw,
      high_metric_uM = high_raw,
      log2_fc_metric_high_vs_low = log2(high_raw / low_raw),
      stringsAsFactors = FALSE
    )
  }))
  pair_stats
}

build_absolute_concentration_metric_stats <- function(fit_results, stat_spec, metric_spec) {
  value_column <- metric_spec$transformed_value_column
  raw_column <- metric_spec$raw_value_column
  absolute_stats <- fit_results[, c(
    "sample", "group", "group_pair_label", "ploidy_class", "ploidy", "ploidy_source",
    "analysis_cloneid_origin", "origin_resolution_status",
    "origin_resolution_distance_edges", "origin_resolution_mrca",
    raw_column, value_column, "model_type"
  ), drop = FALSE]
  names(absolute_stats)[names(absolute_stats) == raw_column] <- "metric_uM"
  names(absolute_stats)[names(absolute_stats) == value_column] <- "metric_value"
  absolute_stats$drug <- "Gemcitabine"
  absolute_stats$response_metric <- metric_spec$response_metric
  absolute_stats$response_metric_label <- metric_spec$response_metric_label
  absolute_stats$response_metric_scale <- "log10_uM"
  absolute_stats$ploidy_stat <- stat_spec$ploidy_stat
  absolute_stats$ploidy_stat_label <- stat_spec$ploidy_label
  absolute_stats$ploidy <- as.numeric(absolute_stats$ploidy)
  absolute_stats$metric_uM <- as.numeric(absolute_stats$metric_uM)
  absolute_stats$metric_value <- as.numeric(absolute_stats$metric_value)
  absolute_stats[, c(
    "drug", "response_metric", "response_metric_label", "response_metric_scale",
    "ploidy_stat", "ploidy_stat_label",
    setdiff(
      names(absolute_stats),
      c(
        "drug", "response_metric", "response_metric_label", "response_metric_scale",
        "ploidy_stat", "ploidy_stat_label"
      )
    )
  ), drop = FALSE]
}

make_delta_concentration_metric_plot <- function(pair_stats, correlation_summary, stat_label, metric_label, group_pair_labels) {
  plot_data <- pair_stats[is.finite(pair_stats$delta_ploidy) & is.finite(pair_stats$delta_metric_high_minus_low), , drop = FALSE]
  annotation <- data.frame(
    drug = "Gemcitabine",
    label = sprintf(
      "R = %.2f, p = %s",
      correlation_summary$pearson_r,
      format_p(correlation_summary$p_value)
    ),
    stringsAsFactors = FALSE
  )

  plot_obj <- ggplot(
    plot_data,
    aes(x = delta_ploidy, y = delta_metric_high_minus_low, color = group)
  ) +
    geom_hline(yintercept = 0, linewidth = 0.45, linetype = "dotted", color = "black") +
    geom_vline(xintercept = 0, linewidth = 0.45, linetype = "dotted", color = "black") +
    geom_point(size = 3.8) +
    geom_text(
      data = annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.35,
      size = 4.2
    ) +
    facet_wrap(~ drug) +
    group_color_scale(plot_data$group, group_pair_labels) +
    expand_limits(x = 0, y = 0) +
    labs(
      title = sprintf("Relative: Delta %s (High - Low), %s ploidy", metric_label, stat_label),
      x = sprintf("Delta ploidy (%s)", stat_label),
      y = sprintf("Delta log10 %s (uM), High - Low", metric_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, hjust = 0, face = "plain"),
      strip.background = element_rect(fill = "grey85", color = "grey35", linewidth = 0.4),
      strip.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.3),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  if (nrow(plot_data) >= 3) {
    plot_obj <- plot_obj +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.9, linetype = "longdash", color = "black")
  }
  plot_obj
}

make_absolute_concentration_metric_plot <- function(absolute_stats, correlation_summary, stat_label, metric_label, group_pair_labels) {
  plot_data <- absolute_stats[is.finite(absolute_stats$ploidy) & is.finite(absolute_stats$metric_value), , drop = FALSE]
  annotation <- data.frame(
    label = sprintf(
      "R = %.2f, p = %s",
      correlation_summary$pearson_r,
      format_p(correlation_summary$p_value)
    ),
    stringsAsFactors = FALSE
  )

  plot_obj <- ggplot(
    plot_data,
    aes(x = ploidy, y = metric_value, color = group, shape = ploidy_class)
  ) +
    geom_point(size = 3.5, alpha = 0.9) +
    geom_text(
      data = annotation,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.02,
      vjust = 1.35,
      size = 4.2
    ) +
    group_color_scale(plot_data$group, group_pair_labels) +
    scale_shape_manual(values = c(Low = 16, High = 17), name = "Ploidy class") +
    labs(
      title = sprintf("Absolute %s ploidy versus %s", stat_label, metric_label),
      x = sprintf("Ploidy (%s)", stat_label),
      y = sprintf("log10 %s (uM)", metric_label)
    ) +
    theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, hjust = 0, face = "plain"),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_line(color = "grey94", linewidth = 0.3),
      legend.position = "right",
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11)
    )
  if (nrow(plot_data) >= 3) {
    plot_obj <- plot_obj +
      geom_smooth(
        data = plot_data,
        aes(x = ploidy, y = metric_value, group = 1),
        inherit.aes = FALSE,
        method = "lm",
        se = FALSE,
        linewidth = 0.9,
        linetype = "longdash",
        color = "black"
      )
  }
  plot_obj
}

write_table_to_bases <- function(table, bases, suffix) {
  for (base in unique(bases)) {
    write.table(
      table,
      file = file.path(table_dir, paste0(base, suffix, ".tsv")),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

save_plot_to_bases <- function(plot_obj, bases, width, height) {
  for (base in unique(bases)) {
    ggsave(
      filename = file.path(fig_dir, paste0(base, ".png")),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 300
    )
    ggsave(
      filename = file.path(fig_dir, paste0(base, ".pdf")),
      plot = plot_obj,
      width = width,
      height = height
    )
  }
}

group_colors <- c(
  HGC = "#F8766D",
  MDA = "#B79F00",
  MDAOld = "#00BA38",
  SNU = "#00BFC4",
  SUM = "#619CFF",
  SUMOld = "#F564E3"
)
group_display_names <- c(
  HGC = "HGC-27",
  MDA = "MDA-MB-231",
  MDAOld = "MDA-MB-231",
  SNU = "SNU-668",
  SUM = "SUM-159",
  SUMOld = "SUM-159"
)

response_long <- read_gemcitabine_table(input_file)
missing_samples <- setdiff(sample_order, unique(response_long$sample))
if (length(missing_samples) > 0) {
  stop("Input table is missing expected samples: ", paste(missing_samples, collapse = ", "), call. = FALSE)
}

cloneid_ploidy <- load_cloneid_ploidy(ploidy_file, ploidy_stat_specs)
fit_results_base <- do.call(rbind, lapply(sample_order, fit_sample_response, response_long = response_long))

mean_spec <- ploidy_stat_specs[ploidy_stat_specs$ploidy_stat == "mean", , drop = FALSE]
mean_metadata <- build_sample_metadata_for_ploidy_stat(sample_metadata_base, cloneid_ploidy, mean_spec[1, , drop = FALSE])
fit_results_mean <- merge(mean_metadata, fit_results_base, by = "sample", sort = FALSE)
fit_results_mean <- fit_results_mean[match(sample_order, fit_results_mean$sample), , drop = FALSE]
mean_group_pair_labels <- build_group_pair_labels(fit_results_mean)
fit_results_mean <- add_group_pair_labels(fit_results_mean, mean_group_pair_labels)

fit_plot_data <- lapply(sample_order, build_normalized_fit_plot_data, response_long = response_long)
fit_raw_points <- do.call(rbind, lapply(fit_plot_data, `[[`, "raw_points"))
fit_predicted_curves <- do.call(rbind, lapply(fit_plot_data, `[[`, "predicted_curve"))

plot_metadata <- sample_metadata_base
plot_metadata$facet_label <- paste0(plot_metadata$group, " ", plot_metadata$ploidy_class)
fit_raw_points <- merge(fit_raw_points, plot_metadata, by = "sample", sort = FALSE)
fit_predicted_curves <- merge(fit_predicted_curves, plot_metadata, by = "sample", sort = FALSE)
fit_raw_points$group <- factor(fit_raw_points$group, levels = unique(plot_metadata$group))
fit_predicted_curves$group <- factor(fit_predicted_curves$group, levels = unique(plot_metadata$group))

dose_breaks <- c(0, 0.005, 0.03, 0.1, 0.3, 1, 5)
dose_response_plot <- ggplot() +
  geom_hline(yintercept = 1, linewidth = 0.35, linetype = "dotted", color = "grey35") +
  geom_point(
    data = fit_raw_points,
    aes(x = dose_uM, y = response_norm, color = ploidy_class),
    alpha = 0.7,
    size = 1.7
  ) +
  geom_line(
    data = fit_predicted_curves,
    aes(x = dose_uM, y = response_norm, color = ploidy_class),
    linewidth = 0.85
  ) +
  facet_wrap(~ group, ncol = 3) +
  scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.001),
    breaks = dose_breaks,
    labels = c("0", "0.005", "0.03", "0.1", "0.3", "1", "5")
  ) +
  scale_color_manual(values = c(Low = "#404040", High = "#D55E00"), name = "Ploidy class") +
  labs(
    title = "Gemcitabine paired normalized dose-response fits",
    x = "Gemcitabine (uM)",
    y = "Response normalized to zero-dose mean"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(size = 14, hjust = 0),
    strip.background = element_rect(fill = "grey88", color = "grey45", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 8),
    legend.position = "bottom"
  )

write.table(
  fit_results_mean,
  file = file.path(table_dir, "gemcitabine_dose_response_fit_parameters.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

delta_pair_stats_list <- list()
delta_correlation_list <- list()
absolute_auc_stats_list <- list()
absolute_correlation_list <- list()
concentration_delta_pair_stats_list <- list()
concentration_delta_correlation_list <- list()
concentration_absolute_stats_list <- list()
concentration_absolute_correlation_list <- list()

for (i in seq_len(nrow(ploidy_stat_specs))) {
  stat_spec <- ploidy_stat_specs[i, , drop = FALSE]
  stat_metadata <- build_sample_metadata_for_ploidy_stat(sample_metadata_base, cloneid_ploidy, stat_spec)
  fit_results <- merge(stat_metadata, fit_results_base, by = "sample", sort = FALSE)
  fit_results <- fit_results[match(sample_order, fit_results$sample), , drop = FALSE]
  group_pair_labels <- build_group_pair_labels(fit_results)
  fit_results <- add_group_pair_labels(fit_results, group_pair_labels)

  pair_stats <- build_pair_stats(fit_results, stat_spec)
  correlation_summary <- correlation_table(
    pair_stats$delta_ploidy,
    pair_stats$delta_auc_high_minus_low,
    "n_pairs"
  )
  correlation_summary$ploidy_stat <- stat_spec$ploidy_stat
  correlation_summary$ploidy_stat_label <- stat_spec$ploidy_label
  correlation_summary <- correlation_summary[, c(
    "drug", "ploidy_stat", "ploidy_stat_label",
    "n_pairs", "pearson_r", "p_value"
  ), drop = FALSE]

  absolute_auc_stats <- build_absolute_auc_stats(fit_results, stat_spec)
  absolute_auc_ok <- is.finite(absolute_auc_stats$ploidy) & is.finite(absolute_auc_stats$auc_norm_log10dose)
  absolute_correlation_summary <- correlation_table(
    absolute_auc_stats$ploidy[absolute_auc_ok],
    absolute_auc_stats$auc_norm_log10dose[absolute_auc_ok],
    "n_samples"
  )
  absolute_correlation_summary$ploidy_stat <- stat_spec$ploidy_stat
  absolute_correlation_summary$ploidy_stat_label <- stat_spec$ploidy_label
  absolute_correlation_summary <- absolute_correlation_summary[, c(
    "drug", "ploidy_stat", "ploidy_stat_label",
    "n_samples", "pearson_r", "p_value"
  ), drop = FALSE]

  delta_plot <- make_delta_plot(pair_stats, correlation_summary, stat_spec$ploidy_label, group_pair_labels)
  absolute_auc_plot <- make_absolute_plot(
    absolute_auc_stats,
    absolute_correlation_summary,
    stat_spec$ploidy_label,
    group_pair_labels
  )

  delta_bases <- paste0("gemcitabine_delta_auc_vs_delta_ploidy_", stat_spec$output_suffix)
  absolute_bases <- paste0("gemcitabine_absolute_auc_vs_ploidy_", stat_spec$output_suffix)
  if (identical(stat_spec$output_suffix, "mean")) {
    delta_bases <- c("gemcitabine_delta_auc_vs_delta_ploidy", delta_bases)
    absolute_bases <- c("gemcitabine_absolute_auc_vs_ploidy", absolute_bases)
  }

  write_table_to_bases(pair_stats, delta_bases, "_pair_stats")
  write_table_to_bases(correlation_summary, delta_bases, "_correlation")
  write_table_to_bases(absolute_auc_stats, absolute_bases, "_sample_stats")
  write_table_to_bases(absolute_correlation_summary, absolute_bases, "_correlation")
  save_plot_to_bases(delta_plot, delta_bases, width = 6.8, height = 4.9)
  save_plot_to_bases(absolute_auc_plot, absolute_bases, width = 6.8, height = 4.9)

  delta_pair_stats_list[[stat_spec$ploidy_stat]] <- pair_stats
  delta_correlation_list[[stat_spec$ploidy_stat]] <- correlation_summary
  absolute_auc_stats_list[[stat_spec$ploidy_stat]] <- absolute_auc_stats
  absolute_correlation_list[[stat_spec$ploidy_stat]] <- absolute_correlation_summary

  for (metric_i in seq_len(nrow(concentration_metric_specs))) {
    metric_spec <- concentration_metric_specs[metric_i, , drop = FALSE]
    metric_key <- paste(metric_spec$response_metric, stat_spec$ploidy_stat, sep = "_")
    metric_pair_stats <- build_pair_concentration_metric_stats(fit_results, stat_spec, metric_spec)
    metric_delta_correlation <- correlation_table(
      metric_pair_stats$delta_ploidy,
      metric_pair_stats$delta_metric_high_minus_low,
      "n_pairs"
    )
    metric_delta_correlation$response_metric <- metric_spec$response_metric
    metric_delta_correlation$response_metric_label <- metric_spec$response_metric_label
    metric_delta_correlation$response_metric_scale <- "log10_uM"
    metric_delta_correlation$ploidy_stat <- stat_spec$ploidy_stat
    metric_delta_correlation$ploidy_stat_label <- stat_spec$ploidy_label
    metric_delta_correlation <- metric_delta_correlation[, c(
      "drug", "response_metric", "response_metric_label", "response_metric_scale",
      "ploidy_stat", "ploidy_stat_label", "n_pairs", "pearson_r", "p_value"
    ), drop = FALSE]

    metric_absolute_stats <- build_absolute_concentration_metric_stats(fit_results, stat_spec, metric_spec)
    metric_absolute_ok <- is.finite(metric_absolute_stats$ploidy) & is.finite(metric_absolute_stats$metric_value)
    metric_absolute_correlation <- correlation_table(
      metric_absolute_stats$ploidy[metric_absolute_ok],
      metric_absolute_stats$metric_value[metric_absolute_ok],
      "n_samples"
    )
    metric_absolute_correlation$response_metric <- metric_spec$response_metric
    metric_absolute_correlation$response_metric_label <- metric_spec$response_metric_label
    metric_absolute_correlation$response_metric_scale <- "log10_uM"
    metric_absolute_correlation$ploidy_stat <- stat_spec$ploidy_stat
    metric_absolute_correlation$ploidy_stat_label <- stat_spec$ploidy_label
    metric_absolute_correlation <- metric_absolute_correlation[, c(
      "drug", "response_metric", "response_metric_label", "response_metric_scale",
      "ploidy_stat", "ploidy_stat_label", "n_samples", "pearson_r", "p_value"
    ), drop = FALSE]

    metric_delta_plot <- make_delta_concentration_metric_plot(
      metric_pair_stats,
      metric_delta_correlation,
      stat_spec$ploidy_label,
      metric_spec$response_metric_label,
      group_pair_labels
    )
    metric_absolute_plot <- make_absolute_concentration_metric_plot(
      metric_absolute_stats,
      metric_absolute_correlation,
      stat_spec$ploidy_label,
      metric_spec$response_metric_label,
      group_pair_labels
    )

    metric_delta_bases <- paste0(metric_spec$delta_base, "_", stat_spec$output_suffix)
    metric_absolute_bases <- paste0(metric_spec$absolute_base, "_", stat_spec$output_suffix)
    if (identical(stat_spec$output_suffix, "mean")) {
      metric_delta_bases <- c(metric_spec$delta_base, metric_delta_bases)
      metric_absolute_bases <- c(metric_spec$absolute_base, metric_absolute_bases)
    }

    write_table_to_bases(metric_pair_stats, metric_delta_bases, "_pair_stats")
    write_table_to_bases(metric_delta_correlation, metric_delta_bases, "_correlation")
    write_table_to_bases(metric_absolute_stats, metric_absolute_bases, "_sample_stats")
    write_table_to_bases(metric_absolute_correlation, metric_absolute_bases, "_correlation")
    save_plot_to_bases(metric_delta_plot, metric_delta_bases, width = 6.8, height = 4.9)
    save_plot_to_bases(metric_absolute_plot, metric_absolute_bases, width = 6.8, height = 4.9)

    concentration_delta_pair_stats_list[[metric_key]] <- metric_pair_stats
    concentration_delta_correlation_list[[metric_key]] <- metric_delta_correlation
    concentration_absolute_stats_list[[metric_key]] <- metric_absolute_stats
    concentration_absolute_correlation_list[[metric_key]] <- metric_absolute_correlation
  }
}

combined_delta_pair_stats <- do.call(rbind, delta_pair_stats_list)
combined_delta_correlations <- do.call(rbind, delta_correlation_list)
combined_absolute_auc_stats <- do.call(rbind, absolute_auc_stats_list)
combined_absolute_correlations <- do.call(rbind, absolute_correlation_list)

write.table(
  combined_delta_pair_stats,
  file = file.path(table_dir, "gemcitabine_delta_auc_vs_delta_ploidy_all_ploidy_stats_pair_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_delta_correlations,
  file = file.path(table_dir, "gemcitabine_delta_auc_vs_delta_ploidy_all_ploidy_stats_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_absolute_auc_stats,
  file = file.path(table_dir, "gemcitabine_absolute_auc_vs_ploidy_all_ploidy_stats_sample_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_absolute_correlations,
  file = file.path(table_dir, "gemcitabine_absolute_auc_vs_ploidy_all_ploidy_stats_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_concentration_delta_pair_stats <- do.call(rbind, concentration_delta_pair_stats_list)
combined_concentration_delta_correlations <- do.call(rbind, concentration_delta_correlation_list)
combined_concentration_absolute_stats <- do.call(rbind, concentration_absolute_stats_list)
combined_concentration_absolute_correlations <- do.call(rbind, concentration_absolute_correlation_list)

for (metric_i in seq_len(nrow(concentration_metric_specs))) {
  metric_spec <- concentration_metric_specs[metric_i, , drop = FALSE]
  metric_delta_rows <- combined_concentration_delta_pair_stats$response_metric == metric_spec$response_metric
  metric_delta_correlation_rows <- combined_concentration_delta_correlations$response_metric == metric_spec$response_metric
  metric_absolute_rows <- combined_concentration_absolute_stats$response_metric == metric_spec$response_metric
  metric_absolute_correlation_rows <- combined_concentration_absolute_correlations$response_metric == metric_spec$response_metric

  write.table(
    combined_concentration_delta_pair_stats[metric_delta_rows, , drop = FALSE],
    file = file.path(table_dir, paste0(metric_spec$delta_base, "_all_ploidy_stats_pair_stats.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    combined_concentration_delta_correlations[metric_delta_correlation_rows, , drop = FALSE],
    file = file.path(table_dir, paste0(metric_spec$delta_base, "_all_ploidy_stats_correlation.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    combined_concentration_absolute_stats[metric_absolute_rows, , drop = FALSE],
    file = file.path(table_dir, paste0(metric_spec$absolute_base, "_all_ploidy_stats_sample_stats.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(
    combined_concentration_absolute_correlations[metric_absolute_correlation_rows, , drop = FALSE],
    file = file.path(table_dir, paste0(metric_spec$absolute_base, "_all_ploidy_stats_correlation.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

write.table(
  combined_concentration_delta_pair_stats,
  file = file.path(table_dir, "gemcitabine_delta_ec50_ic50_vs_delta_ploidy_all_ploidy_stats_pair_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_concentration_delta_correlations,
  file = file.path(table_dir, "gemcitabine_delta_ec50_ic50_vs_delta_ploidy_all_ploidy_stats_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_concentration_absolute_stats,
  file = file.path(table_dir, "gemcitabine_absolute_ec50_ic50_vs_ploidy_all_ploidy_stats_sample_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  combined_concentration_absolute_correlations,
  file = file.path(table_dir, "gemcitabine_absolute_ec50_ic50_vs_ploidy_all_ploidy_stats_correlation.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ggsave(
  filename = file.path(fig_dir, "gemcitabine_normalized_dose_response_fits.png"),
  plot = dose_response_plot,
  width = 10.5,
  height = 7.2,
  dpi = 300
)
ggsave(
  filename = file.path(fig_dir, "gemcitabine_normalized_dose_response_fits.pdf"),
  plot = dose_response_plot,
  width = 10.5,
  height = 7.2
)

mean_delta <- combined_delta_correlations[combined_delta_correlations$ploidy_stat == "mean", , drop = FALSE]
mean_absolute <- combined_absolute_correlations[combined_absolute_correlations$ploidy_stat == "mean", , drop = FALSE]
message(sprintf(
  "Wrote Gemcitabine delta and absolute AUC, EC50, and IC50 plots/tables for %d ploidy statistics to %s (mean-ploidy AUC delta R = %.3f, p = %.4g; mean-ploidy AUC absolute R = %.3f, p = %.4g).",
  nrow(ploidy_stat_specs),
  out_dir,
  mean_delta$pearson_r,
  mean_delta$p_value,
  mean_absolute$pearson_r,
  mean_absolute$p_value
))
