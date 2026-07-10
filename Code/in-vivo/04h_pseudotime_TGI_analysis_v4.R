#!/usr/bin/env Rscript

# EndTimePoint-ploidy (ETP) workflow. The public method names describe the
# actual threshold design rather than the historical V4.1/V4.2/V4.3 labels:
# fixed_threshold_2_25, boundary_stress_threshold_2_375, and
# reference_balanced_threshold_2_24. All thresholds are applied to every
# sample regardless of original ploidy.

suppressPackageStartupMessages({
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    x <- args[[i]]
    if (grepl("^--[^=]+=", x)) {
      key <- sub("^--([^=]+)=.*$", "\\1", x)
      val <- sub("^--[^=]+=", "", x)
      out[[key]] <- val
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
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE)
} else {
  normalizePath("Code/in-vivo/04h_pseudotime_TGI_analysis_v4.R", mustWork = FALSE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)

arg_or <- function(name, default) {
  if (!is.null(cmd_args[[name]]) && nzchar(cmd_args[[name]])) cmd_args[[name]] else default
}

input_root <- normalizePath(
  arg_or("input_root", file.path(repo_root, "Data", "in-vivo")),
  mustWork = FALSE
)
results_base <- arg_or(
  "results_root",
  Sys.getenv(
    "RESULTS_ROOT",
    "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI_analysis_v4"
  )
)
results_base <- normalizePath(results_base, mustWork = FALSE)
etp_thresholds <- c(
  ETP_fixed_threshold_2_25 = 2.25,
  ETP_boundary_stress_threshold_2_375 = 2.375,
  ETP_reference_balanced_threshold_2_24 = 2.24
)
etp_threshold_sources <- c(
  ETP_fixed_threshold_2_25 = "manual_fixed_2.25",
  ETP_boundary_stress_threshold_2_375 = "boundary_stress_test_2.375",
  ETP_reference_balanced_threshold_2_24 = "untreated_reference_balanced_2.24"
)
etp_threshold_aliases <- c(
  all = "all",
  fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
  boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
  reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24",
  ETP_fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
  ETP_boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
  ETP_reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24",
  `V4.1` = "ETP_fixed_threshold_2_25",
  `V4.2` = "ETP_boundary_stress_threshold_2_375",
  `V4.3` = "ETP_reference_balanced_threshold_2_24"
)
requested_threshold <- arg_or("etp_threshold", arg_or("version", "all"))
analysis_version <- unname(etp_threshold_aliases[requested_threshold])
if (length(analysis_version) != 1L || is.na(analysis_version)) {
  stop(
    "--etp_threshold must be one of: all, fixed_threshold_2_25, ",
    "boundary_stress_threshold_2_375, reference_balanced_threshold_2_24",
    call. = FALSE
  )
}
seed <- as.integer(arg_or("seed", "1"))
n_perm <- as.integer(arg_or("n_perm", "10000"))
n_boot <- as.integer(arg_or("n_boot", "1000"))
n_downsample <- as.integer(arg_or("n_downsample", "200"))
run_optional_methods <- tolower(arg_or("run_optional_methods", "TRUE")) %in% c("true", "t", "1", "yes", "y")

if (identical(analysis_version, "all")) {
  dir.create(results_base, recursive = TRUE, showWarnings = FALSE)
  rscript <- file.path(R.home("bin"), "Rscript")
  child_status <- vapply(names(etp_thresholds), function(method) {
    message("Running ", method, " under ", results_base)
    args <- c(
      shQuote(script_path),
      paste0("--etp_threshold=", method),
      shQuote(paste0("--input_root=", input_root)),
      shQuote(paste0("--results_root=", results_base)),
      paste0("--seed=", seed),
      paste0("--n_perm=", n_perm),
      paste0("--n_boot=", n_boot),
      paste0("--n_downsample=", n_downsample),
      paste0("--run_optional_methods=", ifelse(run_optional_methods, "TRUE", "FALSE"))
    )
    system2(rscript, args = args)
  }, integer(1))
  names(child_status) <- names(etp_thresholds)
  write.csv(
    data.frame(method = names(child_status), exit_status = unname(child_status), stringsAsFactors = FALSE),
    file.path(results_base, "all_mode_run_status.csv"),
    row.names = FALSE
  )
  if (any(child_status != 0L)) quit(save = "no", status = 1L)
  quit(save = "no", status = 0L)
}

results_root <- file.path(results_base, analysis_version)
etp_group_levels <- c("ETP-lower", "ETP-higher")
endpoint_threshold <- unname(etp_thresholds[[analysis_version]])
endpoint_threshold_source <- unname(etp_threshold_sources[[analysis_version]])

set.seed(seed)

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

save_plot <- function(plot, output_dir, filename, width = 7, height = 5) {
  ensure_dir(output_dir)
  suppressMessages(ggsave(file.path(output_dir, paste0(filename, ".png")), plot, width = width, height = height, dpi = 300))
  suppressMessages(ggsave(file.path(output_dir, paste0(filename, ".pdf")), plot, width = width, height = height))
  invisible(file.path(output_dir, paste0(filename, ".png")))
}

square_panel <- function(plot) {
  plot + theme(aspect.ratio = 1)
}

save_square_plot <- function(plot, output_dir, filename, width = 6.5, height = 6.5) {
  save_plot(square_panel(plot), output_dir, filename, width, height)
}

format_num <- function(x, digits = 3) {
  ifelse(is.finite(x), trimws(formatC(x, digits = digits, format = "fg")), "NA")
}

format_p <- function(x) {
  ifelse(is.finite(x), ifelse(x < 0.001, formatC(x, digits = 2, format = "e"), formatC(x, digits = 3, format = "fg")), "NA")
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

bool_text <- function(x) ifelse(isTRUE(x), "TRUE", "FALSE")

plot_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", color = "grey70"),
    legend.key.height = unit(0.45, "cm")
  )

dose_cols_all <- c("0mg/kg" = "#666666", "30mg/kg" = "#d95f02", "120mg/kg" = "#1b9e77", "treated" = "#377eb8")

analysis_roots <- list(
  A = file.path(results_root, "A_Pseudotime"),
  B = file.path(results_root, "B_TGI_GrowthResponse"),
  C = file.path(results_root, "C_SupportingAnalyses")
)

dirs <- list(
  CellCycle = file.path(analysis_roots$A, "A1_TreatmentDistributionShift", "CellCycle"),
  NonCellCycle = file.path(analysis_roots$A, "A1_TreatmentDistributionShift", "NonCellCycle"),
  Pseudotime = file.path(analysis_roots$A, "A2_EndTimePointPloidyEffectModification"),
  AUCTGIAssociation = file.path(analysis_roots$B, "B1_AUC_TGI_Association"),
  DoseSpecific = file.path(analysis_roots$B, "B2_DoseSpecific_TGI"),
  PloidyResponse = file.path(analysis_roots$B, "B3_EndTimePointPloidy_TGI_Response"),
  GrowthCurveMixedModel = file.path(analysis_roots$B, "B4_LongitudinalGrowth_MixedModel"),
  GrowthModel = file.path(analysis_roots$B, "B5_GrowthModel_Sensitivity"),
  Robustness = file.path(analysis_roots$C, "C1_Robustness"),
  Comparison = file.path(analysis_roots$C, "C2_CellCycle_vs_NonCellCycle_Comparison"),
  Composition = file.path(analysis_roots$C, "C3_Composition"),
  GenePrograms = file.path(analysis_roots$C, "C4_GenePrograms"),
  Pseudobulk = file.path(analysis_roots$C, "C5_Pseudobulk"),
  NeighborhoodDA = file.path(analysis_roots$C, "C6_NeighborhoodDA"),
  EndpointExploration = file.path(analysis_roots$C, "C7_EndTimePointPloidy_Exploration"),
  Manuscript = file.path(results_root, "Manuscript"),
  logs = file.path(results_root, "logs")
)
invisible(lapply(dirs, ensure_dir))

write_result_structure <- function() {
  write_text(c(
    "# 04h Result Folder Structure",
    "",
    "The result folder follows the two-branch analysis logic used in the manuscript and slides.",
    "",
    "## A_Pseudotime",
    "",
    "- `A1_TreatmentDistributionShift/CellCycle`: CellCycle pseudotime distribution shift, sample-level ECDF tests, shift metrics, and primary CellCycle shift-TGI scatter.",
    "- `A1_TreatmentDistributionShift/NonCellCycle`: the same pseudotime distribution-shift outputs for NonCellCycle cells.",
    "- `A2_EndTimePointPloidyEffectModification`: formal ETP-lower-vs-ETP-higher pseudotime treatment-effect modification using Delta ECDF difference-in-differences curves.",
    "",
    "## B_TGI_GrowthResponse",
    "",
    "- `B1_AUC_TGI_Association`: AUC-TGI models asking whether overall tumor-growth inhibition is associated with pseudotime shift and EndTimePoint ploidy.",
    "- `B2_DoseSpecific_TGI`: dose-specific AUC/endpoint TGI comparisons and within-dose shift-TGI association checks.",
    "- `B3_EndTimePointPloidy_TGI_Response`: treated-sample AUC/endpoint TGI by EndTimePoint ploidy plus supporting response-shift summaries.",
    "- `B4_LongitudinalGrowth_MixedModel`: repeated-measures spline/linear mixed-model growth-response analysis.",
    "- `B5_GrowthModel_Sensitivity`: supportive growth-curve summaries and model-based sensitivity outputs.",
    "",
    "## C_SupportingAnalyses",
    "",
    "- `C1_Robustness`: leave-one-out, bootstrap, downsampling, and confounding checks.",
    "- `C2_CellCycle_vs_NonCellCycle_Comparison`: paired CellCycle-vs-NonCellCycle shift/TGI-comparison outputs.",
    "- `C3_Composition`: CellCycle fraction, cluster-composition, cluster-specific pseudotime, and TGI association outputs.",
    "- `C4_GenePrograms`: pseudotime direction and gene-program support analyses.",
    "- `C5_Pseudobulk`: targeted pseudobulk differential-expression support analyses.",
    "- `C6_NeighborhoodDA`: pseudotime-bin fallback neighborhood differential-abundance analyses.",
    "- `C7_EndTimePointPloidy_Exploration`: dose-colored continuous EndTimePoint-ploidy/TGI and ploidy/pseudotime scatter plots.",
    "",
    "## Manuscript and Logs",
    "",
    "- `Manuscript`: final report, TeX conclusion, and slide deck.",
    "- `logs`: run parameters and R session information.",
    "",
    "Legacy flat folders from earlier runs can be preserved under `Legacy_FlatLayout_before_AB_structure` after a manual cleanup step."
  ), file.path(results_root, "RESULT_STRUCTURE.md"))
}
write_result_structure()

run_parameters <- list(
  script_path = script_path,
  repo_root = repo_root,
  input_root = input_root,
  results_root = results_root,
  analysis_version = analysis_version,
  endpoint_ploidy_threshold = endpoint_threshold,
  endpoint_ploidy_rule = "sample_mean_ploidy > threshold is ETP-higher; otherwise ETP-lower",
  threshold_source = endpoint_threshold_source,
  seed = seed,
  n_perm = n_perm,
  n_boot = n_boot,
  n_downsample = n_downsample,
  run_optional_methods = run_optional_methods,
  run_time = as.character(Sys.time())
)

write_yaml_like <- function(x, path) {
  lines <- unlist(Map(function(k, v) paste0(k, ": ", as.character(v)), names(x), x), use.names = FALSE)
  write_text(lines, path)
}

write_json_like <- function(x, path) {
  ensure_dir(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  } else {
    lines <- c("{", paste0('  "', names(x), '": "', gsub('"', '\\"', as.character(x)), '"', collapse = ",\n"), "}")
    writeLines(lines, path)
  }
  invisible(path)
}

if (requireNamespace("yaml", quietly = TRUE)) {
  yaml::write_yaml(run_parameters, file.path(dirs$logs, "run_parameters.yaml"))
} else {
  write_yaml_like(run_parameters, file.path(dirs$logs, "run_parameters.yaml"))
}

capture.output(sessionInfo(), file = file.path(dirs$logs, "sessionInfo.txt"))

compartments <- list(
  CellCycle = list(
    label = "Cell-cycle-associated tumor cells",
    prefix = "CellCycle",
    input = file.path(input_root, "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    dir = dirs$CellCycle
  ),
  NonCellCycle = list(
    label = "Non-cell-cycle-associated tumor cells",
    prefix = "NonCellCycle",
    input = file.path(input_root, "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    dir = dirs$NonCellCycle
  )
)

derive_sample_etp_assignments <- function(compartment_cfg, threshold) {
  raw_list <- lapply(names(compartment_cfg), function(name) {
    path <- compartment_cfg[[name]]$input
    if (!file.exists(path)) stop("Missing ETP assignment input: ", path, call. = FALSE)
    raw <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
    needed <- c(
      "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
      "gemcitabine_dose_mg_per_kg", "cell_ploidy"
    )
    missing <- setdiff(needed, names(raw))
    if (length(missing) > 0) {
      stop("ETP assignment input is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    raw <- raw[, needed, drop = FALSE]
    raw$source_compartment <- name
    raw$cell_ploidy <- safe_num(raw$cell_ploidy)
    raw$gemcitabine_dose_mg_per_kg <- safe_num(raw$gemcitabine_dose_mg_per_kg)
    raw
  })
  union_cells <- do.call(rbind, raw_list)
  union_cells <- union_cells[is.finite(union_cells$cell_ploidy), , drop = FALSE]
  if (nrow(union_cells) == 0) stop("No finite cell_ploidy values are available for ETP grouping", call. = FALSE)

  duplicated_ids <- unique(union_cells$cell_id[duplicated(union_cells$cell_id)])
  if (length(duplicated_ids) > 0) {
    for (cell_id in duplicated_ids) {
      local <- union_cells[union_cells$cell_id == cell_id, , drop = FALSE]
      if (length(unique(local$sample_id)) != 1L || length(unique(local$cell_ploidy)) != 1L) {
        stop("Conflicting duplicated cell_id in ETP inputs: ", cell_id, call. = FALSE)
      }
    }
    union_cells <- union_cells[!duplicated(union_cells$cell_id), , drop = FALSE]
  }

  sample_rows <- lapply(split(union_cells, union_cells$sample_id), function(local) {
    original_ploidy <- unique(as.character(local$initial_ploidy))
    dose <- unique(as.character(local$gemcitabine_dose))
    dose_mg <- unique(local$gemcitabine_dose_mg_per_kg)
    if (length(original_ploidy) != 1L || length(dose) != 1L || length(dose_mg) != 1L) {
      stop("Inconsistent sample metadata for ", local$sample_id[[1]], call. = FALSE)
    }
    data.frame(
      sample_id = local$sample_id[[1]],
      original_initial_ploidy = original_ploidy,
      dose = dose,
      dose_mg = dose_mg,
      n_endpoint_ploidy_cells = nrow(local),
      sample_mean_ploidy = mean(local$cell_ploidy),
      sample_median_ploidy = median(local$cell_ploidy),
      sample_max_ploidy = max(local$cell_ploidy),
      stringsAsFactors = FALSE
    )
  })
  assignments <- do.call(rbind, sample_rows)
  assignments$end_timepoint_ploidy_group <- factor(
    ifelse(assignments$sample_mean_ploidy > threshold, "ETP-higher", "ETP-lower"),
    levels = etp_group_levels
  )
  assignments$threshold <- threshold
  assignments$threshold_source <- endpoint_threshold_source
  assignments$version <- analysis_version
  assignments <- assignments[
    order(assignments$dose_mg, assignments$end_timepoint_ploidy_group, assignments$sample_id),
    ,
    drop = FALSE
  ]
  rownames(assignments) <- NULL

  threshold_definition <- data.frame(
    version = analysis_version,
    threshold = threshold,
    threshold_source = endpoint_threshold_source,
    assignment_population = "all samples regardless of original initial ploidy",
    sample_mean_source = "deduplicated union of CellCycle and NonCellCycle cell_ploidy",
    comparison_rule = "sample_mean_ploidy > threshold is ETP-higher; otherwise ETP-lower",
    n_samples = nrow(assignments),
    n_cells = nrow(union_cells),
    original_2N_mean_of_sample_means = mean(assignments$sample_mean_ploidy[assignments$original_initial_ploidy == "2N"]),
    original_4N_mean_of_sample_means = mean(assignments$sample_mean_ploidy[assignments$original_initial_ploidy == "4N"]),
    stringsAsFactors = FALSE
  )
  write_csv(threshold_definition, file.path(results_root, "threshold_definition.csv"))
  write_csv(assignments, file.path(results_root, "sample_end_timepoint_ploidy_assignments.csv"))

  observed <- aggregate(
    sample_id ~ dose + dose_mg + end_timepoint_ploidy_group,
    assignments,
    function(x) length(unique(x))
  )
  names(observed)[names(observed) == "sample_id"] <- "n_samples"
  full <- expand.grid(
    dose = unique(assignments$dose[order(assignments$dose_mg)]),
    end_timepoint_ploidy_group = etp_group_levels,
    stringsAsFactors = FALSE
  )
  dose_map <- unique(assignments[, c("dose", "dose_mg")])
  full$dose_mg <- dose_map$dose_mg[match(full$dose, dose_map$dose)]
  coverage <- merge(
    full,
    observed,
    by = c("dose", "dose_mg", "end_timepoint_ploidy_group"),
    all.x = TRUE
  )
  coverage$n_samples[is.na(coverage$n_samples)] <- 0L
  coverage$empty_group <- coverage$n_samples == 0L
  coverage$version <- analysis_version
  coverage$threshold <- threshold
  coverage <- coverage[order(coverage$dose_mg, match(coverage$end_timepoint_ploidy_group, etp_group_levels)), ]
  write_csv(coverage, file.path(results_root, "ETP_group_dose_coverage_audit.csv"))

  observed_initial <- aggregate(
    sample_id ~ original_initial_ploidy + dose + dose_mg + end_timepoint_ploidy_group,
    assignments,
    function(x) length(unique(x))
  )
  names(observed_initial)[names(observed_initial) == "sample_id"] <- "n_samples"
  full_initial <- expand.grid(
    original_initial_ploidy = sort(unique(assignments$original_initial_ploidy)),
    dose = unique(assignments$dose[order(assignments$dose_mg)]),
    end_timepoint_ploidy_group = etp_group_levels,
    stringsAsFactors = FALSE
  )
  full_initial$dose_mg <- dose_map$dose_mg[match(full_initial$dose, dose_map$dose)]
  coverage_initial <- merge(
    full_initial,
    observed_initial,
    by = c("original_initial_ploidy", "dose", "dose_mg", "end_timepoint_ploidy_group"),
    all.x = TRUE
  )
  coverage_initial$n_samples[is.na(coverage_initial$n_samples)] <- 0L
  coverage_initial$empty_group <- coverage_initial$n_samples == 0L
  coverage_initial$version <- analysis_version
  coverage_initial$threshold <- threshold
  coverage_initial$assignment_population <- "all samples regardless of original initial ploidy"
  coverage_initial <- coverage_initial[
    order(
      match(coverage_initial$original_initial_ploidy, c("2N", "4N")),
      coverage_initial$dose_mg,
      match(coverage_initial$end_timepoint_ploidy_group, etp_group_levels)
    ),
    ,
    drop = FALSE
  ]
  write_csv(
    coverage_initial,
    file.path(results_root, "ETP_group_initial_ploidy_dose_coverage_audit.csv")
  )
  assignments
}

etp_assignments <- derive_sample_etp_assignments(compartments, endpoint_threshold)

required_cols <- c(
  "sample_id", "cluster", "initial_ploidy", "gemcitabine_dose",
  "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy", "TGI_percent_auc"
)

endpoint_tgi_cols <- function(nms) {
  cols <- grep("^TGI_percent_Day_[0-9]+$", nms, value = TRUE)
  days <- safe_num(sub("^TGI_percent_Day_", "", cols))
  cols[order(days)]
}

day_cols_by_prefix <- function(nms, prefix) {
  cols <- grep(paste0("^", prefix, "Day_[0-9]+$"), nms, value = TRUE)
  days <- safe_num(sub(paste0("^", prefix, "Day_"), "", cols))
  cols[order(days)]
}

endpoint_day_num <- function(col) safe_num(sub("^TGI_percent_Day_", "", col))

day_label_num <- function(x) safe_num(gsub("[^0-9.]+", "", as.character(x)))

ploidy_numeric <- function(x) {
  x <- as.character(x)
  ifelse(x == "ETP-lower", 0, ifelse(x == "ETP-higher", 1, NA_real_))
}

safe_cor_result <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.frame(n = length(x), estimate = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  ct <- suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE))
  data.frame(n = length(x), estimate = unname(ct$estimate), p_value = ct$p.value, stringsAsFactors = FALSE)
}

safe_wilcox <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]
  group <- as.character(group[ok])
  if (length(value) < 3 || length(unique(group)) != 2) {
    return(data.frame(n = length(value), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  wt <- try(suppressWarnings(stats::wilcox.test(value ~ group, exact = FALSE)), silent = TRUE)
  if (inherits(wt, "try-error")) {
    return(data.frame(n = length(value), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(n = length(value), statistic = unname(wt$statistic), p_value = wt$p.value, stringsAsFactors = FALSE)
}

safe_ttest <- function(value, group) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]
  group <- as.character(group[ok])
  group_counts <- table(group)
  if (length(value) < 4 || length(group_counts) != 2 || any(group_counts < 2)) {
    return(data.frame(n = length(value), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  tt <- try(suppressWarnings(stats::t.test(value ~ group)), silent = TRUE)
  if (inherits(tt, "try-error")) {
    return(data.frame(n = length(value), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(n = length(value), statistic = unname(tt$statistic), p_value = tt$p.value, stringsAsFactors = FALSE)
}

safe_paired <- function(x, y, method = "paired_t") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) {
    return(data.frame(n = length(x), statistic = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
  }
  if (identical(method, "paired_t")) {
    tst <- suppressWarnings(stats::t.test(x, y, paired = TRUE))
  } else {
    tst <- suppressWarnings(stats::wilcox.test(x, y, paired = TRUE, exact = FALSE))
  }
  data.frame(n = length(x), statistic = unname(tst$statistic), p_value = tst$p.value, stringsAsFactors = FALSE)
}

read_compartment_data <- function(name, cfg) {
  if (!file.exists(cfg$input)) stop("Missing input CSV for ", name, ": ", cfg$input, call. = FALSE)
  df <- read.csv(cfg$input, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) stop("Missing required columns in ", cfg$input, ": ", paste(missing, collapse = ", "), call. = FALSE)
  numeric_cols <- grep("^(pseudotime|cell_ploidy|average_ploidy|median_ploidy|sample_|n_|target_|gemcitabine_dose_mg_per_kg|TGI_percent|tumor_volume|matched_control|auc_)", names(df), value = TRUE)
  numeric_cols <- setdiff(numeric_cols, c("sample_id", "growth_curve_sample_id_raw", "growth_curve_harvest"))
  for (col in numeric_cols) df[[col]] <- safe_num(df[[col]])
  df <- df[is.finite(df$pseudotime), , drop = FALSE]
  assignment_index <- match(df$sample_id, etp_assignments$sample_id)
  if (anyNA(assignment_index)) {
    stop("Missing ETP assignment for samples: ", paste(unique(df$sample_id[is.na(assignment_index)]), collapse = ", "), call. = FALSE)
  }
  df$original_initial_ploidy <- as.character(df$initial_ploidy)
  expected_initial <- etp_assignments$original_initial_ploidy[assignment_index]
  if (any(df$original_initial_ploidy != expected_initial)) {
    stop("Original initial-ploidy metadata disagree with the ETP assignment table", call. = FALSE)
  }
  df$end_timepoint_ploidy_group <- factor(
    as.character(etp_assignments$end_timepoint_ploidy_group[assignment_index]),
    levels = etp_group_levels
  )
  df$sample_mean_endpoint_ploidy <- etp_assignments$sample_mean_ploidy[assignment_index]
  df$sample_median_endpoint_ploidy <- etp_assignments$sample_median_ploidy[assignment_index]
  df$sample_max_endpoint_ploidy <- etp_assignments$sample_max_ploidy[assignment_index]
  df$compartment <- name
  df$dose <- df$gemcitabine_dose
  df$dose_mg <- df$gemcitabine_dose_mg_per_kg
  df
}

build_sample_meta <- function(df) {
  sample_ids <- sort(unique(df$sample_id))
  sample_level_cols <- setdiff(names(df), c("cell_id", "cluster", "pseudotime", "cell_ploidy", "compartment", "dose", "dose_mg"))
  rows <- lapply(sample_ids, function(sample_id) {
    sub <- df[df$sample_id == sample_id, , drop = FALSE]
    row <- sub[1, sample_level_cols, drop = FALSE]
    row$n_cells <- nrow(sub)
    row$clusters <- paste(sort(unique(as.character(sub$cluster))), collapse = ";")
    row$mean_pseudotime <- mean(sub$pseudotime, na.rm = TRUE)
    row$median_pseudotime <- median(sub$pseudotime, na.rm = TRUE)
    row$mean_cell_ploidy <- mean(sub$cell_ploidy, na.rm = TRUE)
    row$median_cell_ploidy <- median(sub$cell_ploidy, na.rm = TRUE)
    row$max_cell_ploidy <- max(sub$cell_ploidy, na.rm = TRUE)
    row$p90_cell_ploidy <- as.numeric(quantile(sub$cell_ploidy, 0.90, na.rm = TRUE, names = FALSE))
    row
  })
  meta <- do.call(rbind, rows)
  names(meta) <- make.unique(names(meta), sep = "_dup")
  rownames(meta) <- meta$sample_id
  meta[["dose"]] <- meta[["gemcitabine_dose"]]
  meta[["dose_mg"]] <- safe_num(meta[["gemcitabine_dose_mg_per_kg"]])
  meta[["end_timepoint_ploidy_group_numeric"]] <- ploidy_numeric(meta[["end_timepoint_ploidy_group"]])
  meta <- meta[order(meta[["dose_mg"]], meta[["end_timepoint_ploidy_group_numeric"]], meta[["sample_id"]]), , drop = FALSE]
  rownames(meta) <- meta$sample_id
  meta
}

make_grid <- function(df, n = 501L) {
  rng <- range(df$pseudotime, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(sort(unique(df$pseudotime)))
  seq(rng[1], rng[2], length.out = n)
}

make_sample_ecdf_matrix <- function(df, sample_ids, grid) {
  mat <- t(vapply(sample_ids, function(sample_id) {
    x <- df$pseudotime[df$sample_id == sample_id]
    if (length(x) == 0) rep(NA_real_, length(grid)) else stats::ecdf(x)(grid)
  }, numeric(length(grid))))
  rownames(mat) <- sample_ids
  colnames(mat) <- paste0("grid_", seq_along(grid))
  mat
}

mean_row <- function(mat, sample_ids) {
  sample_ids <- intersect(sample_ids, rownames(mat))
  if (length(sample_ids) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[sample_ids, , drop = FALSE], na.rm = TRUE)
}

calculate_shift_metrics <- function(df, sample_meta, reference_type = c("primary_equal_sample_reference", "sensitivity_pooled_cell_reference")) {
  reference_type <- match.arg(reference_type)
  grid <- make_grid(df)
  sample_ids <- rownames(sample_meta)
  ecdf_mat <- make_sample_ecdf_matrix(df, sample_ids, grid)
  sample_mean_pt <- setNames(sample_meta$mean_pseudotime, sample_meta$sample_id)
  rows <- lapply(sample_ids, function(sample_id) {
    sample_ploidy <- sample_meta[sample_id, "end_timepoint_ploidy_group"]
    sample_dose <- sample_meta[sample_id, "dose_mg"]
    controls <- sample_meta$sample_id[sample_meta$end_timepoint_ploidy_group == sample_ploidy & sample_meta$dose_mg == 0]
    reference_samples <- if (sample_dose == 0) setdiff(controls, sample_id) else controls
    if (length(reference_samples) == 0) {
      delta <- rep(NA_real_, length(grid))
      reference_mean <- NA_real_
      reference_cells <- 0L
    } else if (identical(reference_type, "primary_equal_sample_reference")) {
      ref_ecdf <- mean_row(ecdf_mat, reference_samples)
      delta <- ecdf_mat[sample_id, ] - ref_ecdf
      reference_mean <- mean(sample_mean_pt[reference_samples], na.rm = TRUE)
      reference_cells <- sum(sample_meta[reference_samples, "n_cells"], na.rm = TRUE)
    } else {
      sample_pseudotime <- df$pseudotime[df$sample_id == sample_id]
      reference_pseudotime <- df$pseudotime[df$sample_id %in% reference_samples]
      delta <- stats::ecdf(sample_pseudotime)(grid) - stats::ecdf(reference_pseudotime)(grid)
      reference_mean <- mean(reference_pseudotime, na.rm = TRUE)
      reference_cells <- length(reference_pseudotime)
    }
    finite_delta <- delta[is.finite(delta)]
    ecdf_rmse <- if (length(finite_delta) > 0) sqrt(mean(finite_delta^2)) else NA_real_
    ecdf_ks <- if (length(finite_delta) > 0) max(abs(finite_delta)) else NA_real_
    ecdf_mean_abs <- if (length(finite_delta) > 0) mean(abs(finite_delta)) else NA_real_
    data.frame(
      sample_id = sample_id,
      reference_type = reference_type,
      end_timepoint_ploidy_group = sample_ploidy,
      end_timepoint_ploidy_group_numeric = sample_meta[sample_id, "end_timepoint_ploidy_group_numeric"],
      dose = sample_meta[sample_id, "dose"],
      dose_mg = sample_dose,
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
      ecdf_rmse = ecdf_rmse,
      ecdf_ks = ecdf_ks,
      ecdf_mean_abs = ecdf_mean_abs,
      signed_mean_shift = sample_meta[sample_id, "mean_pseudotime"] - reference_mean,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  for (col in c("TGI_percent", endpoint_tgi_cols(names(sample_meta)))) {
    if (col %in% names(sample_meta)) out[[col]] <- sample_meta[out$sample_id, col]
  }
  out <- out[order(out$dose_mg, out$end_timepoint_ploidy_group_numeric, out$sample_id), , drop = FALSE]
  rownames(out) <- NULL
  out
}

choice_lists_for_strata <- function(groups, group_a, strata) {
  split_idx <- split(seq_along(groups), strata)
  lapply(split_idx, function(idx) {
    k <- sum(groups[idx] == group_a)
    if (k == 0) return(list(integer(0)))
    if (k == length(idx)) return(list(idx))
    combn(idx, k, simplify = FALSE)
  })
}

expand_exact_group_vectors <- function(groups, group_a, group_b, strata) {
  choices_by_strata <- choice_lists_for_strata(groups, group_a, strata)
  out <- list()
  recurse <- function(i, selected) {
    if (i > length(choices_by_strata)) {
      g <- rep(group_b, length(groups))
      g[unlist(selected, use.names = FALSE)] <- group_a
      out[[length(out) + 1L]] <<- g
    } else {
      for (choice in choices_by_strata[[i]]) recurse(i + 1L, c(selected, list(choice)))
    }
  }
  recurse(1L, list())
  out
}

sample_group_vector <- function(groups, group_a, group_b, strata) {
  g <- rep(group_b, length(groups))
  split_idx <- split(seq_along(groups), strata)
  for (idx in split_idx) {
    k <- sum(groups[idx] == group_a)
    if (k > 0) g[sample(idx, k)] <- group_a
  }
  g
}

ecdf_distance_from_groups <- function(ecdf_mat, groups, group_a, group_b) {
  if (!any(groups == group_a) || !any(groups == group_b)) {
    return(c(ecdf_rmse = NA_real_, ecdf_ks = NA_real_, ecdf_mean_abs = NA_real_))
  }
  delta <- colMeans(ecdf_mat[groups == group_a, , drop = FALSE], na.rm = TRUE) -
    colMeans(ecdf_mat[groups == group_b, , drop = FALSE], na.rm = TRUE)
  c(ecdf_rmse = sqrt(mean(delta^2, na.rm = TRUE)), ecdf_ks = max(abs(delta), na.rm = TRUE), ecdf_mean_abs = mean(abs(delta), na.rm = TRUE))
}

ecdf_group_test <- function(df, group_col, group_a, group_b, stratify_by_end_timepoint_ploidy_group = FALSE, n_perm = 10000L) {
  local <- df[df[[group_col]] %in% c(group_a, group_b), , drop = FALSE]
  sample_ids <- sort(unique(local$sample_id))
  if (length(sample_ids) < 3) {
    small_meta <- local[!duplicated(local$sample_id), c("sample_id", group_col), drop = FALSE]
    small_groups <- as.character(small_meta[[group_col]])
    return(data.frame(
      group_a = group_a, group_b = group_b, stratified_by_end_timepoint_ploidy_group = stratify_by_end_timepoint_ploidy_group,
      n_samples = length(sample_ids), n_group_a = sum(small_groups == group_a), n_group_b = sum(small_groups == group_b),
      observed_ecdf_rmse = NA_real_, observed_ecdf_ks = NA_real_, observed_ecdf_mean_abs = NA_real_,
      p_ecdf_rmse = NA_real_, p_ecdf_ks = NA_real_, p_ecdf_mean_abs = NA_real_,
      n_permutations = NA_integer_, permutation_mode = "not_estimable_insufficient_samples", stringsAsFactors = FALSE
    ))
  }
  meta <- local[!duplicated(local$sample_id), c("sample_id", "end_timepoint_ploidy_group", group_col), drop = FALSE]
  meta <- meta[match(sample_ids, meta$sample_id), , drop = FALSE]
  groups <- as.character(meta[[group_col]])
  strata <- if (stratify_by_end_timepoint_ploidy_group) as.character(meta$end_timepoint_ploidy_group) else rep("all", length(groups))
  if (!all(c(group_a, group_b) %in% groups)) {
    return(data.frame(
      group_a = group_a, group_b = group_b, stratified_by_end_timepoint_ploidy_group = stratify_by_end_timepoint_ploidy_group,
      n_samples = length(sample_ids), n_group_a = sum(groups == group_a), n_group_b = sum(groups == group_b),
      observed_ecdf_rmse = NA_real_, observed_ecdf_ks = NA_real_, observed_ecdf_mean_abs = NA_real_,
      p_ecdf_rmse = NA_real_, p_ecdf_ks = NA_real_, p_ecdf_mean_abs = NA_real_,
      n_permutations = NA_integer_, permutation_mode = "not_estimable_missing_group", stringsAsFactors = FALSE
    ))
  }
  grid <- make_grid(local)
  ecdf_mat <- make_sample_ecdf_matrix(local, sample_ids, grid)
  observed <- ecdf_distance_from_groups(ecdf_mat, groups, group_a, group_b)
  choices_by_strata <- choice_lists_for_strata(groups, group_a, strata)
  n_exact <- prod(vapply(choices_by_strata, length, integer(1)))
  if (is.finite(n_exact) && n_exact <= n_perm) {
    perm_groups <- expand_exact_group_vectors(groups, group_a, group_b, strata)
    permutation_mode <- "exact_label_enumeration"
  } else {
    perm_groups <- replicate(n_perm, sample_group_vector(groups, group_a, group_b, strata), simplify = FALSE)
    permutation_mode <- "monte_carlo_label_permutation"
  }
  perm_stats <- t(vapply(perm_groups, function(g) ecdf_distance_from_groups(ecdf_mat, g, group_a, group_b), numeric(3)))
  data.frame(
    group_a = group_a,
    group_b = group_b,
    stratified_by_end_timepoint_ploidy_group = stratify_by_end_timepoint_ploidy_group,
    n_samples = length(sample_ids),
    n_group_a = sum(groups == group_a),
    n_group_b = sum(groups == group_b),
    observed_ecdf_rmse = observed[["ecdf_rmse"]],
    observed_ecdf_ks = observed[["ecdf_ks"]],
    observed_ecdf_mean_abs = observed[["ecdf_mean_abs"]],
    p_ecdf_rmse = mean(perm_stats[, "ecdf_rmse"] >= observed[["ecdf_rmse"]] - 1e-15, na.rm = TRUE),
    p_ecdf_ks = mean(perm_stats[, "ecdf_ks"] >= observed[["ecdf_ks"]] - 1e-15, na.rm = TRUE),
    p_ecdf_mean_abs = mean(perm_stats[, "ecdf_mean_abs"] >= observed[["ecdf_mean_abs"]] - 1e-15, na.rm = TRUE),
    n_permutations = nrow(perm_stats),
    permutation_mode = permutation_mode,
    stringsAsFactors = FALSE
  )
}

run_dose_tests <- function(df, n_perm = 10000L) {
  local <- df
  local$combined_treatment_group <- ifelse(local$dose_mg == 0, "0mg/kg", "treated")
  out <- list(
    cbind(comparison = "0_vs_30plus120", ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", FALSE, n_perm)),
    cbind(comparison = "0_vs_30plus120", ecdf_group_test(local, "combined_treatment_group", "0mg/kg", "treated", TRUE, n_perm))
  )
  dose_values <- sort(unique(local$dose_mg))
  if (length(dose_values) >= 2) {
    pairs <- combn(dose_values, 2)
    for (i in seq_len(ncol(pairs))) {
      da <- pairs[1, i]
      db <- pairs[2, i]
      sub <- local[local$dose_mg %in% c(da, db), , drop = FALSE]
      group_a <- unique(sub$dose[sub$dose_mg == da])[1]
      group_b <- unique(sub$dose[sub$dose_mg == db])[1]
      cmp <- paste0(da, "_vs_", db)
      out[[length(out) + 1L]] <- cbind(comparison = cmp, ecdf_group_test(sub, "dose", group_a, group_b, FALSE, n_perm))
    }
  }
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans
}

.perm_cache <- new.env(parent = emptyenv())

perm_index_matrix <- function(n) {
  key <- as.character(n)
  if (exists(key, envir = .perm_cache, inherits = FALSE)) return(get(key, envir = .perm_cache))
  total <- as.integer(round(factorial(n)))
  mat <- matrix(NA_integer_, nrow = total, ncol = n)
  j <- 0L
  recurse <- function(prefix, rest) {
    if (length(rest) == 0) {
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

permutation_cor_p <- function(x, y, method = "pearson", n_perm = 10000L, positive_direction = TRUE, exact = FALSE) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(data.frame(
      n = n, estimate = NA_real_, asymptotic_p = NA_real_,
      permutation_p_two_sided = NA_real_, permutation_p_positive = NA_real_,
      n_permutations = NA_integer_, permutation_mode = NA_character_, stringsAsFactors = FALSE
    ))
  }
  observed <- suppressWarnings(cor(x, y, method = method))
  asym <- safe_cor_result(x, y, method)
  if (isTRUE(exact) && n <= 9 && factorial(n) <= 400000) {
    idx <- perm_index_matrix(n)
    perm_est <- vapply(seq_len(nrow(idx)), function(i) suppressWarnings(cor(x, y[idx[i, ]], method = method)), numeric(1))
    mode <- "exact_TGI_label_enumeration"
  } else {
    perm_est <- replicate(n_perm, suppressWarnings(cor(x, sample(y), method = method)))
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

run_tgi_associations <- function(equal_shift, pooled_shift, compartment_name, n_perm = 10000L, exact_auc = TRUE) {
  refs <- list(primary_equal_sample_reference = equal_shift, sensitivity_pooled_cell_reference = pooled_shift)
  shift_cols <- c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")
  tgi_cols <- unique(c("TGI_percent_auc", endpoint_tgi_cols(names(equal_shift))))
  rows <- list()
  for (ref_name in names(refs)) {
    treated <- refs[[ref_name]][refs[[ref_name]]$dose_mg > 0, , drop = FALSE]
    for (shift_col in shift_cols) {
      for (tgi_col in tgi_cols) {
        use_exact <- isTRUE(exact_auc) && identical(tgi_col, "TGI_percent_auc")
        pearson <- permutation_cor_p(treated[[shift_col]], treated[[tgi_col]], "pearson", n_perm, exact = use_exact)
        spearman <- permutation_cor_p(treated[[shift_col]], treated[[tgi_col]], "spearman", n_perm, exact = use_exact)
        rows[[length(rows) + 1L]] <- data.frame(
          compartment = compartment_name,
          sample_set = "treated",
          reference_type = ref_name,
          shift_metric = shift_col,
          tgi_measure = tgi_col,
          tgi_day = ifelse(grepl("^TGI_percent_Day_", tgi_col), endpoint_day_num(tgi_col), NA_real_),
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
          pre_specified_primary = compartment_name == "CellCycle" &&
            ref_name == "primary_equal_sample_reference" &&
            shift_col == "ecdf_rmse" &&
            tgi_col == "TGI_percent_auc",
          stringsAsFactors = FALSE
        )
      }
    }
  }
  out <- do.call(rbind, rows)
  out$pearson_q_exploratory <- p.adjust(out$pearson_p_permutation_two_sided, method = "BH")
  out$spearman_q_exploratory <- p.adjust(out$spearman_p_permutation_two_sided, method = "BH")
  out
}

plot_sample_ecdf <- function(df, sample_meta, cfg) {
  fig_dir <- ensure_dir(file.path(cfg$dir, "figures"))
  grid <- make_grid(df)
  sample_ids <- rownames(sample_meta)
  mat <- make_sample_ecdf_matrix(df, sample_ids, grid)
  plot_df <- do.call(rbind, lapply(sample_ids, function(sample_id) {
    data.frame(
      sample_id = sample_id,
      pseudotime = grid,
      ecdf = mat[sample_id, ],
      dose = sample_meta[sample_id, "dose"],
      end_timepoint_ploidy_group = sample_meta[sample_id, "end_timepoint_ploidy_group"],
      stringsAsFactors = FALSE
    )
  }))
  dose_levels <- unique(sample_meta$dose[order(sample_meta$dose_mg)])
  plot_df$dose <- factor(plot_df$dose, levels = dose_levels)
  p <- ggplot(plot_df, aes(x = pseudotime, y = ecdf, color = dose, group = sample_id)) +
    geom_line(alpha = 0.55, linewidth = 0.55) +
    facet_wrap(~ end_timepoint_ploidy_group, ncol = 1) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = paste0(cfg$label, ": sample-equal ECDFs by dose and ploidy"), x = "pseudotime", y = "sample ECDF") +
    plot_theme
  save_square_plot(p, fig_dir, paste0(cfg$prefix, "_sample_equal_ecdf_by_dose_and_ploidy"), 7.2, 9.4)
}

plot_tgi_association <- function(shift, cfg) {
  fig_dir <- ensure_dir(file.path(cfg$dir, "figures"))
  treated <- shift[shift$dose_mg > 0, , drop = FALSE]
  dose_levels <- unique(treated$dose[order(treated$dose_mg)])
  p <- ggplot(treated, aes(x = ecdf_rmse, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
    geom_point(aes(shape = end_timepoint_ploidy_group), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.4, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(
      title = paste0(cfg$label, ": AUC TGI vs sample-equal ECDF shift"),
      x = "ECDF RMSE from ploidy-matched equal-sample 0 mg/kg reference",
      y = "AUC-based TGI (%)"
    ) +
    plot_theme
  save_square_plot(p, fig_dir, paste0(cfg$prefix, "_TGI_AUC_vs_ecdf_rmse_equal_sample_ref"), 6.8, 6.8)
}

plot_basic_compartment <- function(df, sample_meta, equal_shift, cfg) {
  fig_dir <- ensure_dir(file.path(cfg$dir, "figures"))
  dose_levels <- unique(sample_meta$dose[order(sample_meta$dose_mg)])
  df$dose_factor <- factor(df$dose, levels = dose_levels)
  df$end_timepoint_ploidy_group_factor <- factor(df$end_timepoint_ploidy_group, levels = etp_group_levels)
  p_density <- ggplot(df, aes(x = pseudotime, color = dose_factor, fill = dose_factor)) +
    geom_density(alpha = 0.12, linewidth = 0.8) +
    facet_wrap(~ end_timepoint_ploidy_group_factor, ncol = 1) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    scale_fill_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = paste0(cfg$label, ": cell-level pseudotime density"), x = "pseudotime", y = "density") +
    plot_theme
  save_plot(p_density, fig_dir, paste0(cfg$prefix, "_density_by_dose_and_ploidy"), 8, 5.5)
  shift_plot <- ggplot(equal_shift, aes(x = dose, y = ecdf_rmse, color = dose)) +
    geom_boxplot(aes(group = dose), outlier.shape = NA, color = "grey55", linewidth = 0.45) +
    geom_point(aes(shape = end_timepoint_ploidy_group), position = position_jitter(width = 0.08, height = 0), size = 2.6) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = paste0(cfg$label, ": sample-equal ECDF shift by dose"), x = "dose", y = "ECDF RMSE") +
    plot_theme
  save_square_plot(shift_plot, fig_dir, paste0(cfg$prefix, "_ecdf_rmse_by_dose_equal_sample_ref"), 6.4, 6.4)
}

run_compartment <- function(name, cfg, n_perm) {
  ensure_dir(file.path(cfg$dir, "tables"))
  ensure_dir(file.path(cfg$dir, "stats"))
  ensure_dir(file.path(cfg$dir, "figures"))
  df <- read_compartment_data(name, cfg)
  sample_meta <- build_sample_meta(df)
  equal_shift <- calculate_shift_metrics(df, sample_meta, "primary_equal_sample_reference")
  pooled_shift <- calculate_shift_metrics(df, sample_meta, "sensitivity_pooled_cell_reference")
  dose_tests <- run_dose_tests(df, n_perm)
  tgi_stats <- run_tgi_associations(equal_shift, pooled_shift, name, n_perm)
  write_csv(sample_meta, file.path(cfg$dir, "tables", paste0(cfg$prefix, "_sample_metadata.csv")))
  write_csv(equal_shift, file.path(cfg$dir, "tables", paste0(cfg$prefix, "_shift_metrics_equal_sample_ref.csv")))
  write_csv(pooled_shift, file.path(cfg$dir, "tables", paste0(cfg$prefix, "_shift_metrics_pooled_ref_sensitivity.csv")))
  write_csv(dose_tests, file.path(cfg$dir, "stats", "dose_group_ecdf_tests_v4.csv"))
  write_csv(tgi_stats, file.path(cfg$dir, "stats", paste0(cfg$prefix, "_TGI_associations_v4.csv")))
  plot_sample_ecdf(df, sample_meta, cfg)
  plot_tgi_association(equal_shift, cfg)
  plot_basic_compartment(df, sample_meta, equal_shift, cfg)
  list(df = df, sample_meta = sample_meta, equal_shift = equal_shift, pooled_shift = pooled_shift, dose_tests = dose_tests, tgi_stats = tgi_stats)
}

plot_endpoint_continuous_exploration <- function(name, result) {
  comp_dir <- ensure_dir(file.path(dirs$EndpointExploration, name))
  table_dir <- ensure_dir(file.path(comp_dir, "tables"))
  figure_dir <- ensure_dir(file.path(comp_dir, "figures"))
  meta <- result$sample_meta
  meta$end_timepoint_ploidy_group <- factor(meta$end_timepoint_ploidy_group, levels = etp_group_levels)
  meta$dose <- factor(meta$dose, levels = unique(meta$dose[order(meta$dose_mg)]))

  ploidy_metrics <- c("mean_cell_ploidy", "median_cell_ploidy", "max_cell_ploidy")
  ploidy_long <- do.call(rbind, lapply(ploidy_metrics, function(metric) {
    data.frame(
      meta,
      endpoint_ploidy_metric = metric,
      endpoint_ploidy_value = meta[[metric]],
      stringsAsFactors = FALSE
    )
  }))

  tgi_cols <- intersect(c("TGI_percent_auc", "TGI_percent"), names(meta))
  if (length(tgi_cols) > 0) {
    tgi_long <- do.call(rbind, lapply(tgi_cols, function(metric) {
      data.frame(
        ploidy_long,
        tgi_measure = metric,
        tgi_value = ploidy_long[[metric]],
        stringsAsFactors = FALSE
      )
    }))
    write_csv(tgi_long, file.path(table_dir, paste0(name, "_endpoint_ploidy_mean_median_max_vs_TGI_plot_data.csv")))
    p_tgi <- ggplot(
      tgi_long,
      aes(x = endpoint_ploidy_value, y = tgi_value, color = dose, shape = end_timepoint_ploidy_group)
    ) +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_smooth(
        data = tgi_long,
        aes(x = endpoint_ploidy_value, y = tgi_value, group = 1),
        inherit.aes = FALSE,
        method = "lm",
        se = FALSE,
        color = "black",
        linewidth = 0.45
      ) +
      geom_point(size = 2.7) +
      geom_text(aes(label = sample_id), size = 2.1, nudge_y = 2, check_overlap = TRUE, show.legend = FALSE) +
      facet_grid(tgi_measure ~ endpoint_ploidy_metric, scales = "free") +
      scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% levels(meta$dose)], drop = FALSE, name = "dose") +
      scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
      scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
      coord_cartesian(clip = "off") +
      labs(
        title = paste0(name, ": sample EndTimePoint ploidy summaries vs TGI"),
        subtitle = paste0(analysis_version, "; fixed ETP threshold = ", format_num(endpoint_threshold, 6)),
        x = "EndTimePoint ploidy",
        y = "TGI (%)",
        shape = "EndTimePoint ploidy group"
      ) +
      plot_theme +
      theme(plot.margin = margin(8, 24, 8, 24))
    save_plot(p_tgi, figure_dir, paste0(name, "_endpoint_ploidy_mean_median_max_vs_TGI"), 12, 7.2)
  }

  pseudotime_metrics <- c("mean_pseudotime", "median_pseudotime")
  pseudotime_long <- do.call(rbind, lapply(pseudotime_metrics, function(metric) {
    data.frame(
      ploidy_long,
      pseudotime_metric = metric,
      pseudotime_value = ploidy_long[[metric]],
      stringsAsFactors = FALSE
    )
  }))
  write_csv(pseudotime_long, file.path(table_dir, paste0(name, "_endpoint_ploidy_mean_median_max_vs_pseudotime_plot_data.csv")))
  p_pt <- ggplot(
    pseudotime_long,
    aes(x = endpoint_ploidy_value, y = pseudotime_value, color = dose, shape = end_timepoint_ploidy_group)
  ) +
    geom_smooth(
      data = pseudotime_long,
      aes(x = endpoint_ploidy_value, y = pseudotime_value, group = 1),
      inherit.aes = FALSE,
      method = "lm",
      se = FALSE,
      color = "black",
      linewidth = 0.45
    ) +
    geom_point(size = 2.7) +
    geom_text(aes(label = sample_id), size = 2.1, nudge_y = 0.015, check_overlap = TRUE, show.legend = FALSE) +
    facet_grid(pseudotime_metric ~ endpoint_ploidy_metric, scales = "free") +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% levels(meta$dose)], drop = FALSE, name = "dose") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0.18))) +
    scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(name, ": sample EndTimePoint ploidy summaries vs pseudotime"),
      subtitle = paste0(analysis_version, "; fixed ETP threshold = ", format_num(endpoint_threshold, 6)),
      x = "EndTimePoint ploidy",
      y = "sample pseudotime",
      shape = "EndTimePoint ploidy group"
    ) +
    plot_theme +
    theme(plot.margin = margin(8, 24, 8, 24))
  save_plot(p_pt, figure_dir, paste0(name, "_endpoint_ploidy_mean_median_max_vs_pseudotime"), 12, 7.2)

  cell_df <- result$df
  cell_df$dose <- factor(cell_df$dose, levels = unique(meta$dose))
  cell_df$end_timepoint_ploidy_group <- factor(cell_df$end_timepoint_ploidy_group, levels = etp_group_levels)
  cell_plot_data <- cell_df[, c(
    "cell_id", "sample_id", "original_initial_ploidy", "dose", "dose_mg",
    "end_timepoint_ploidy_group", "cell_ploidy", "pseudotime"
  )]
  write_csv(cell_plot_data, file.path(table_dir, paste0(name, "_single_cell_ploidy_vs_single_cell_pseudotime_plot_data.csv")))
  p_cell <- ggplot(cell_plot_data, aes(x = cell_ploidy, y = pseudotime, color = dose)) +
    geom_vline(xintercept = endpoint_threshold, linetype = "dashed", color = "black", linewidth = 0.45) +
    geom_point(size = 0.45, alpha = 0.18) +
    facet_wrap(~ end_timepoint_ploidy_group, ncol = 1, drop = FALSE) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% levels(meta$dose)], drop = FALSE, name = "dose") +
    labs(
      title = paste0(name, ": single-cell ploidy vs single-cell pseudotime"),
      subtitle = paste0("Dashed line: fixed sample-level ETP threshold = ", format_num(endpoint_threshold, 6)),
      x = "single-cell ploidy",
      y = "single-cell pseudotime"
    ) +
    plot_theme
  save_plot(p_cell, figure_dir, paste0(name, "_single_cell_ploidy_vs_single_cell_pseudotime_by_dose"), 8.4, 7.2)
}

sample_level_mean_difference_test <- function(
  data,
  value_col,
  group_col,
  group_a,
  group_b,
  n_perm = 10000L
) {
  value <- safe_num(data[[value_col]])
  group <- as.character(data[[group_col]])
  keep <- is.finite(value) & !is.na(group) & group %in% as.character(c(group_a, group_b))
  value <- value[keep]
  group <- group[keep]
  group_a <- as.character(group_a)
  group_b <- as.character(group_b)
  n_a <- sum(group == group_a)
  n_b <- sum(group == group_b)
  mean_a <- if (n_a > 0) mean(value[group == group_a]) else NA_real_
  mean_b <- if (n_b > 0) mean(value[group == group_b]) else NA_real_
  observed <- if (n_a > 0 && n_b > 0) mean_b - mean_a else NA_real_

  base <- data.frame(
    group_a = group_a,
    group_b = group_b,
    n_group_a = n_a,
    n_group_b = n_b,
    mean_group_a = mean_a,
    mean_group_b = mean_b,
    observed_mean_diff_group_b_minus_a = observed,
    permutation_p_two_sided = NA_real_,
    n_permutations = NA_integer_,
    permutation_mode = NA_character_,
    estimable = FALSE,
    non_estimable_reason = NA_character_,
    stringsAsFactors = FALSE
  )
  if (n_a == 0 || n_b == 0) {
    base$permutation_mode <- "not_estimable_missing_group"
    base$non_estimable_reason <- "at least one comparison group has no biological sample"
    return(base)
  }
  if (n_a < 2 || n_b < 2) {
    base$permutation_mode <- "not_estimable_insufficient_biological_samples"
    base$non_estimable_reason <- "at least two biological samples are required in each comparison group"
    return(base)
  }

  n <- length(value)
  n_exact <- choose(n, n_a)
  if (is.finite(n_exact) && n_exact <= n_perm) {
    group_a_indices <- combn(seq_len(n), n_a, simplify = FALSE)
    permuted_differences <- vapply(group_a_indices, function(index) {
      permuted_group <- rep(group_b, n)
      permuted_group[index] <- group_a
      mean(value[permuted_group == group_b]) - mean(value[permuted_group == group_a])
    }, numeric(1))
    mode <- "exact_sample_label_enumeration"
  } else {
    permuted_differences <- replicate(n_perm, {
      index <- sample(seq_len(n), n_a)
      permuted_group <- rep(group_b, n)
      permuted_group[index] <- group_a
      mean(value[permuted_group == group_b]) - mean(value[permuted_group == group_a])
    })
    mode <- "monte_carlo_sample_label_permutation"
  }
  base$permutation_p_two_sided <- mean(
    abs(permuted_differences) >= abs(observed) - 1e-15,
    na.rm = TRUE
  )
  base$n_permutations <- length(permuted_differences)
  base$permutation_mode <- mode
  base$estimable <- TRUE
  base$non_estimable_reason <- ""
  base
}

format_signed_num <- function(x, digits = 3) {
  ifelse(is.finite(x), sprintf(paste0("%+.", digits, "f"), x), "NA")
}

format_plot_test_p <- function(p, estimable) {
  ifelse(estimable & is.finite(p), format_p(p), "NE")
}

plot_sample_level_etp_ploidy_violin <- function(name, result, cfg, n_perm) {
  df <- result$df
  assignment_index <- match(as.character(df$sample_id), as.character(etp_assignments$sample_id))
  if (anyNA(assignment_index)) {
    stop("Plot data contain samples without a global ETP assignment", call. = FALSE)
  }
  global_etp_group <- as.character(etp_assignments$end_timepoint_ploidy_group[assignment_index])
  if (any(as.character(df$end_timepoint_ploidy_group) != global_etp_group)) {
    stop("Plot data ETP labels disagree with the global all-sample assignment table", call. = FALSE)
  }
  # Reapply the single global sample assignment explicitly before plotting.
  # No ETP threshold is recalculated within initial-ploidy or dose strata.
  df$end_timepoint_ploidy_group <- factor(global_etp_group, levels = etp_group_levels)
  dose_lookup <- unique(df[, c("dose", "dose_mg")])
  dose_lookup <- dose_lookup[order(dose_lookup$dose_mg), , drop = FALSE]
  dose_levels <- as.character(dose_lookup$dose)
  initial_levels <- c("2N", "4N")
  etp_offsets <- c("ETP-lower" = -0.20, "ETP-higher" = 0.20)
  df$dose_plot <- factor(as.character(df$dose), levels = dose_levels)
  df$original_initial_ploidy <- factor(
    as.character(df$original_initial_ploidy),
    levels = initial_levels
  )
  df$initial_etp_x <- match(as.character(df$dose), dose_levels) +
    unname(etp_offsets[as.character(df$end_timepoint_ploidy_group)])

  sample_summary <- do.call(rbind, lapply(split(df, as.character(df$sample_id)), function(local) {
    data.frame(
      sample_id = as.character(local$sample_id[[1]]),
      original_initial_ploidy = as.character(local$original_initial_ploidy[[1]]),
      dose = as.character(local$dose[[1]]),
      dose_mg = safe_num(local$dose_mg[[1]]),
      end_timepoint_ploidy_group = as.character(local$end_timepoint_ploidy_group[[1]]),
      sample_mean_cell_ploidy = mean(local$cell_ploidy, na.rm = TRUE),
      sample_median_cell_ploidy = median(local$cell_ploidy, na.rm = TRUE),
      ETP_classification_sample_mean_ploidy = safe_num(local$sample_mean_endpoint_ploidy[[1]]),
      n_cells = nrow(local),
      stringsAsFactors = FALSE
    )
  }))
  rownames(sample_summary) <- NULL
  sample_summary$dose_plot <- factor(sample_summary$dose, levels = dose_levels)
  sample_summary$end_timepoint_ploidy_group <- factor(
    sample_summary$end_timepoint_ploidy_group,
    levels = etp_group_levels
  )
  sample_summary$original_initial_ploidy <- factor(
    sample_summary$original_initial_ploidy,
    levels = initial_levels
  )
  sample_summary$initial_etp_x <- match(sample_summary$dose, dose_levels) +
    unname(etp_offsets[as.character(sample_summary$end_timepoint_ploidy_group)])
  sample_summary$version <- analysis_version
  sample_summary$threshold <- endpoint_threshold
  sample_summary$analysis_unit <- "biological_sample"
  write_csv(
    sample_summary,
    file.path(
      cfg$dir,
      "tables",
      paste0(cfg$prefix, "_ploidy_violin_dose_end_timepoint_ploidy_group_sample_medians.csv")
    )
  )

  dose_pairs <- if (nrow(dose_lookup) >= 2) {
    combn(seq_len(nrow(dose_lookup)), 2, simplify = FALSE)
  } else {
    list()
  }
  dose_test_rows <- list()
  for (etp_group in etp_group_levels) {
    local <- sample_summary[
      as.character(sample_summary$end_timepoint_ploidy_group) == etp_group,
      ,
      drop = FALSE
    ]
    for (pair in dose_pairs) {
      dose_a <- dose_lookup$dose_mg[pair[[1]]]
      dose_b <- dose_lookup$dose_mg[pair[[2]]]
      test <- sample_level_mean_difference_test(
        local,
        "sample_mean_cell_ploidy",
        "dose_mg",
        dose_a,
        dose_b,
        n_perm
      )
      dose_test_rows[[length(dose_test_rows) + 1L]] <- data.frame(
        version = analysis_version,
        compartment = name,
        comparison_family = "pairwise_dose_within_ETP_group",
        end_timepoint_ploidy_group = etp_group,
        dose_a = as.character(dose_lookup$dose[pair[[1]]]),
        dose_b = as.character(dose_lookup$dose[pair[[2]]]),
        test,
        analysis_unit = "biological_sample_mean_cell_ploidy",
        stringsAsFactors = FALSE
      )
    }
  }
  dose_tests <- if (length(dose_test_rows) > 0) do.call(rbind, dose_test_rows) else data.frame()
  if (nrow(dose_tests) > 0) {
    dose_tests$permutation_q_BH_exploratory <- p.adjust(
      dose_tests$permutation_p_two_sided,
      method = "BH"
    )
    dose_tests$interpretation_caveat <- paste0(
      "Dose comparisons are conditional on sample ETP groups defined from endpoint ploidy; ",
      "treat as exploratory."
    )
  }
  write_csv(
    dose_tests,
    file.path(
      cfg$dir,
      "stats",
      paste0(cfg$prefix, "_ploidy_violin_ETP_facets_pairwise_dose_mean_tests.csv")
    )
  )

  first_dose_factor <- factor(dose_levels[[1]], levels = dose_levels)
  y_range <- range(df$cell_ploidy, na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  dose_annotations <- do.call(rbind, lapply(etp_group_levels, function(etp_group) {
    local <- dose_tests[
      as.character(dose_tests$end_timepoint_ploidy_group) == etp_group,
      ,
      drop = FALSE
    ]
    lines <- if (nrow(local) > 0) vapply(seq_len(nrow(local)), function(i) {
      paste0(
        local$dose_a[[i]], " vs ", local$dose_b[[i]],
        ": d=", format_signed_num(local$observed_mean_diff_group_b_minus_a[[i]]),
        "; p/q=", format_plot_test_p(local$permutation_p_two_sided[[i]], local$estimable[[i]]),
        "/", format_plot_test_p(local$permutation_q_BH_exploratory[[i]], local$estimable[[i]]),
        "; n=", local$n_group_a[[i]], "/", local$n_group_b[[i]]
      )
    }, character(1)) else "No dose comparison available"
    data.frame(
      end_timepoint_ploidy_group = factor(etp_group, levels = etp_group_levels),
      dose_plot = first_dose_factor,
      annotation_y = y_range[[2]] + 0.30 * y_span,
      stat_label = paste(lines, collapse = "\n"),
      stringsAsFactors = FALSE
    )
  }))

  etp_cols <- c("ETP-lower" = "#4C78A8", "ETP-higher" = "#E45756")
  p_etp_facets <- ggplot(
    df,
    aes(
      x = dose_plot,
      y = cell_ploidy,
      fill = end_timepoint_ploidy_group,
      group = dose_plot
    )
  ) +
    geom_violin(
      width = 0.78,
      scale = "width",
      trim = FALSE,
      alpha = 0.58,
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_hline(yintercept = endpoint_threshold, linetype = "dashed", color = "grey20", linewidth = 0.5) +
    geom_point(
      data = sample_summary,
      aes(
        x = dose_plot,
        y = sample_mean_cell_ploidy,
        color = end_timepoint_ploidy_group
      ),
      position = position_jitter(width = 0.08, height = 0),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    geom_label(
      data = dose_annotations,
      aes(x = dose_plot, y = annotation_y, label = stat_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.8,
      lineheight = 1.05,
      linewidth = 0.2,
      label.padding = unit(0.14, "lines"),
      fill = "white"
    ) +
    facet_wrap(~ end_timepoint_ploidy_group, nrow = 1, scales = "fixed", drop = FALSE) +
    scale_fill_manual(values = etp_cols, limits = etp_group_levels, drop = FALSE, name = "ETP group") +
    scale_color_manual(values = etp_cols, limits = etp_group_levels, drop = FALSE, guide = "none") +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.05))) +
    labs(
      title = paste0(cfg$label, ": ploidy by dose within sample-level ETP groups"),
      subtitle = paste0(
        analysis_version,
        "; dashed line = fixed sample-mean ETP threshold ",
        format_num(endpoint_threshold, 6),
        "; violins = cells; points = biological-sample means"
      ),
      caption = paste0(
        "Exact sample-label tests compare sample mean cell ploidy between doses. ",
        "\n",
        "d = second dose minus first; n = first/second; p = raw; q = BH-adjusted; ",
        "NE = not estimable (<2 biological samples in either group)."
      ),
      x = "gemcitabine dose",
      y = "single-cell ploidy"
    ) +
    plot_theme +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = 8.5),
      plot.margin = margin(8, 10, 8, 10)
    )
  save_plot(
    p_etp_facets,
    file.path(cfg$dir, "figures"),
    paste0(cfg$prefix, "_ploidy_violin_dose_end_timepoint_ploidy_group"),
    12.0,
    6.6
  )

  initial_etp_test_rows <- list()
  for (initial_group in initial_levels) {
    for (dose_index in seq_len(nrow(dose_lookup))) {
      dose_value <- dose_lookup$dose_mg[[dose_index]]
      local <- sample_summary[
        as.character(sample_summary$original_initial_ploidy) == initial_group &
          sample_summary$dose_mg == dose_value,
        ,
        drop = FALSE
      ]
      test <- sample_level_mean_difference_test(
        local,
        "sample_mean_cell_ploidy",
        "end_timepoint_ploidy_group",
        "ETP-lower",
        "ETP-higher",
        n_perm
      )
      initial_etp_test_rows[[length(initial_etp_test_rows) + 1L]] <- data.frame(
        version = analysis_version,
        compartment = name,
        comparison_family = "ETP_group_within_initial_ploidy_and_dose",
        original_initial_ploidy = initial_group,
        dose = as.character(dose_lookup$dose[[dose_index]]),
        dose_mg = dose_value,
        test,
        analysis_unit = "biological_sample_mean_cell_ploidy",
        stringsAsFactors = FALSE
      )
    }
  }
  initial_etp_tests <- do.call(rbind, initial_etp_test_rows)
  initial_etp_tests$permutation_q_BH_exploratory <- p.adjust(
    initial_etp_tests$permutation_p_two_sided,
    method = "BH"
  )
  initial_etp_tests$interpretation_caveat <- paste0(
    "ETP groups are defined from endpoint ploidy, so ploidy contrasts between ETP groups ",
    "are descriptive and selection-linked rather than independent validation."
  )
  write_csv(
    initial_etp_tests,
    file.path(
      cfg$dir,
      "stats",
      paste0(cfg$prefix, "_ploidy_violin_initial_ploidy_dose_ETP_mean_tests.csv")
    )
  )
  write_csv(
    sample_summary,
    file.path(
      cfg$dir,
      "tables",
      paste0(cfg$prefix, "_ploidy_violin_initial_ploidy_dose_endpoint_group_sample_summary.csv")
    )
  )

  initial_annotations <- data.frame(
    original_initial_ploidy = factor(
      initial_etp_tests$original_initial_ploidy,
      levels = initial_levels
    ),
    dose_x = match(initial_etp_tests$dose, dose_levels),
    annotation_y = y_range[[2]] + 0.23 * y_span,
    stat_label = vapply(seq_len(nrow(initial_etp_tests)), function(i) {
      paste0(
        "H-L=", format_signed_num(initial_etp_tests$observed_mean_diff_group_b_minus_a[[i]]),
        "\np/q=", format_plot_test_p(
          initial_etp_tests$permutation_p_two_sided[[i]],
          initial_etp_tests$estimable[[i]]
        ),
        "/", format_plot_test_p(
          initial_etp_tests$permutation_q_BH_exploratory[[i]],
          initial_etp_tests$estimable[[i]]
        ),
        "; n=", initial_etp_tests$n_group_a[[i]], "/", initial_etp_tests$n_group_b[[i]]
      )
    }, character(1)),
    stringsAsFactors = FALSE
  )

  initial_slot_labels <- expand.grid(
    original_initial_ploidy = initial_levels,
    dose_x = seq_along(dose_levels),
    end_timepoint_ploidy_group = etp_group_levels,
    stringsAsFactors = FALSE
  )
  initial_slot_labels$slot_x <- initial_slot_labels$dose_x +
    unname(etp_offsets[initial_slot_labels$end_timepoint_ploidy_group])
  initial_slot_labels$slot_y <- y_range[[1]] - 0.08 * y_span
  initial_slot_labels$slot_label <- ifelse(
    initial_slot_labels$end_timepoint_ploidy_group == "ETP-lower",
    "L",
    "H"
  )
  initial_slot_labels$original_initial_ploidy <- factor(
    initial_slot_labels$original_initial_ploidy,
    levels = initial_levels
  )
  initial_slot_labels$end_timepoint_ploidy_group <- factor(
    initial_slot_labels$end_timepoint_ploidy_group,
    levels = etp_group_levels
  )

  p_initial_etp <- ggplot(
    df,
    aes(
      x = initial_etp_x,
      y = cell_ploidy,
      fill = end_timepoint_ploidy_group,
      group = interaction(dose_plot, end_timepoint_ploidy_group)
    )
  ) +
    geom_violin(
      width = 0.36,
      scale = "width",
      trim = TRUE,
      alpha = 0.58,
      color = "grey35",
      linewidth = 0.35
    ) +
    geom_hline(yintercept = endpoint_threshold, linetype = "dashed", color = "grey20", linewidth = 0.5) +
    geom_point(
      data = sample_summary,
      aes(
        x = initial_etp_x,
        y = sample_mean_cell_ploidy,
        color = end_timepoint_ploidy_group
      ),
      position = position_jitter(width = 0.035, height = 0),
      inherit.aes = FALSE,
      size = 2.5
    ) +
    geom_text(
      data = initial_slot_labels,
      aes(
        x = slot_x,
        y = slot_y,
        label = slot_label,
        color = end_timepoint_ploidy_group
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 3.0,
      show.legend = FALSE
    ) +
    geom_label(
      data = initial_annotations,
      aes(x = dose_x, y = annotation_y, label = stat_label),
      inherit.aes = FALSE,
      hjust = 0.5,
      vjust = 1,
      size = 2.65,
      lineheight = 1.02,
      linewidth = 0.2,
      label.padding = unit(0.12, "lines"),
      fill = "white"
    ) +
    facet_wrap(~ original_initial_ploidy, nrow = 1, scales = "fixed", drop = FALSE) +
    scale_fill_manual(values = etp_cols, limits = etp_group_levels, drop = FALSE, name = "ETP group") +
    scale_color_manual(values = etp_cols, limits = etp_group_levels, drop = FALSE, guide = "none") +
    scale_x_continuous(
      breaks = seq_along(dose_levels),
      labels = dose_levels,
      limits = c(0.52, length(dose_levels) + 0.48),
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.05))) +
    labs(
      title = paste0(cfg$label, ": ETP-group ploidy differences within initial-ploidy strata"),
      subtitle = paste0(
        analysis_version,
        "; facets = original initial ploidy; dashed line = fixed ETP threshold ",
        format_num(endpoint_threshold, 6),
        "; points = biological-sample means"
      ),
      caption = paste0(
        "Within each initial-ploidy x dose stratum, exact sample-label tests compare ",
        "ETP-higher minus ETP-lower sample mean cell ploidy.",
        "\n",
        "H-L = ETP-higher minus ETP-lower; n = lower/higher; ",
        "p = raw; q = BH-adjusted; NE = not estimable (<2 biological samples/group or a missing group).",
        "\n",
        "Within each dose, fixed L/H slots show ETP-lower/ETP-higher; an empty slot means n=0. ",
        "ETP is defined from endpoint ploidy, so H-L contrasts are descriptive and selection-linked."
      ),
      x = "gemcitabine dose",
      y = "single-cell ploidy"
    ) +
    plot_theme +
    theme(
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, size = 8.5),
      plot.margin = margin(8, 10, 8, 10)
    )
  save_plot(
    p_initial_etp,
    file.path(cfg$dir, "figures"),
    paste0(cfg$prefix, "_ploidy_violin_initial_ploidy_dose_endpoint_group"),
    12.0,
    6.6
  )
}

mean_sample_ecdf_curve_data <- function(df, panel, grid) {
  rows <- lapply(unique(as.character(df$curve_label)), function(curve_label) {
    local <- df[as.character(df$curve_label) == curve_label, , drop = FALSE]
    sample_ids <- unique(as.character(local$sample_id))
    if (length(sample_ids) == 0) return(NULL)
    mat <- t(vapply(sample_ids, function(sample_id) {
      values <- local$pseudotime[as.character(local$sample_id) == sample_id]
      stats::ecdf(values)(grid)
    }, numeric(length(grid))))
    data.frame(
      panel = panel,
      pseudotime = grid,
      mean_ecdf = colMeans(mat, na.rm = TRUE),
      curve_label = curve_label,
      color_group = as.character(local$color_group[[1]]),
      line_group = as.character(local$line_group[[1]]),
      n_samples = length(sample_ids),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0) return(data.frame())
  do.call(rbind, rows)
}

plot_direct_group_ecdf_comparisons <- function(name, result, cfg, n_perm) {
  df <- result$df
  df$treatment_group <- ifelse(df$dose_mg == 0, "0mg/kg", "treated")
  grid <- make_grid(df)
  curve_rows <- list()
  add_curves <- function(local, panel, curve_label, color_group, line_group) {
    if (nrow(local) == 0) return(invisible(NULL))
    local$curve_label <- curve_label
    local$color_group <- color_group
    local$line_group <- line_group
    curve_rows[[length(curve_rows) + 1L]] <<- mean_sample_ecdf_curve_data(local, panel, grid)
    invisible(NULL)
  }

  add_curves(df, "1. 0 vs treated", df$treatment_group, df$treatment_group, "All")
  add_curves(
    df,
    "2. 0 vs treated, ETP-stratified",
    paste(df$treatment_group, df$end_timepoint_ploidy_group, sep = " | "),
    df$treatment_group,
    as.character(df$end_timepoint_ploidy_group)
  )
  dose_specs <- list(
    list(panel = "3. 0 vs 30 mg/kg", doses = c(0, 30)),
    list(panel = "4. 0 vs 120 mg/kg", doses = c(0, 120)),
    list(panel = "5. 30 vs 120 mg/kg", doses = c(30, 120))
  )
  for (spec in dose_specs) {
    local <- df[df$dose_mg %in% spec$doses, , drop = FALSE]
    add_curves(local, spec$panel, as.character(local$dose), as.character(local$dose), "All")
  }
  treated <- df[df$dose_mg > 0, , drop = FALSE]
  add_curves(
    treated,
    "6. Treated: ETP-lower vs ETP-higher",
    as.character(treated$end_timepoint_ploidy_group),
    "treated",
    as.character(treated$end_timepoint_ploidy_group)
  )
  control <- df[df$dose_mg == 0, , drop = FALSE]
  add_curves(
    control,
    "7. Control: ETP-lower vs ETP-higher",
    as.character(control$end_timepoint_ploidy_group),
    "0mg/kg",
    as.character(control$end_timepoint_ploidy_group)
  )
  high <- df[df$end_timepoint_ploidy_group == "ETP-higher", , drop = FALSE]
  add_curves(high, "8. ETP-higher: 0 vs treated", high$treatment_group, high$treatment_group, "ETP-higher")
  low <- df[df$end_timepoint_ploidy_group == "ETP-lower", , drop = FALSE]
  add_curves(low, "9. ETP-lower: 0 vs treated", low$treatment_group, low$treatment_group, "ETP-lower")
  curves <- do.call(rbind, curve_rows)

  overall_tests <- result$dose_tests
  get_overall_test <- function(comparison, stratified_flag = FALSE) {
    hit <- overall_tests[
      overall_tests$comparison == comparison &
        overall_tests$stratified_by_end_timepoint_ploidy_group == stratified_flag,
      ,
      drop = FALSE
    ]
    if (nrow(hit) == 0) return(data.frame())
    hit[1, , drop = FALSE]
  }
  etp_between <- function(local) {
    ecdf_group_test(local, "end_timepoint_ploidy_group", "ETP-lower", "ETP-higher", FALSE, n_perm)
  }
  treatment_within_etp <- function(group) {
    local <- df[df$end_timepoint_ploidy_group == group, , drop = FALSE]
    ecdf_group_test(local, "treatment_group", "0mg/kg", "treated", FALSE, n_perm)
  }
  test_rows <- list(
    cbind(panel = "1. 0 vs treated", test_design = "sample_label_permutation", get_overall_test("0_vs_30plus120", FALSE)),
    cbind(panel = "2. 0 vs treated, ETP-stratified", test_design = "treatment_label_permutation_within_ETP_group", get_overall_test("0_vs_30plus120", TRUE)),
    cbind(panel = "3. 0 vs 30 mg/kg", test_design = "sample_label_permutation", get_overall_test("0_vs_30", FALSE)),
    cbind(panel = "4. 0 vs 120 mg/kg", test_design = "sample_label_permutation", get_overall_test("0_vs_120", FALSE)),
    cbind(panel = "5. 30 vs 120 mg/kg", test_design = "sample_label_permutation", get_overall_test("30_vs_120", FALSE)),
    cbind(panel = "6. Treated: ETP-lower vs ETP-higher", test_design = "unpaired_sample_ETP_label_permutation", etp_between(treated)),
    cbind(panel = "7. Control: ETP-lower vs ETP-higher", test_design = "unpaired_sample_ETP_label_permutation", etp_between(control)),
    cbind(panel = "8. ETP-higher: 0 vs treated", test_design = "sample_treatment_label_permutation_within_ETP_group", treatment_within_etp("ETP-higher")),
    cbind(panel = "9. ETP-lower: 0 vs treated", test_design = "sample_treatment_label_permutation_within_ETP_group", treatment_within_etp("ETP-lower"))
  )
  all_test_columns <- unique(unlist(lapply(test_rows, names), use.names = FALSE))
  test_rows <- lapply(test_rows, function(x) {
    missing <- setdiff(all_test_columns, names(x))
    for (column in missing) x[[column]] <- NA
    x[, all_test_columns, drop = FALSE]
  })
  tests <- do.call(rbind, test_rows)
  tests$q_ecdf_rmse_exploratory <- p.adjust(tests$p_ecdf_rmse, method = "BH")
  tests$q_ecdf_ks_exploratory <- p.adjust(tests$p_ecdf_ks, method = "BH")
  tests$q_ecdf_mean_abs_exploratory <- p.adjust(tests$p_ecdf_mean_abs, method = "BH")
  tests$version <- analysis_version
  tests$threshold <- endpoint_threshold
  tests$estimable <- is.finite(tests$p_ecdf_rmse)
  tests$annotation <- paste0(
    tests$group_a, " n=", tests$n_group_a, "; ", tests$group_b, " n=", tests$n_group_b,
    "\nRMSE=", format_num(tests$observed_ecdf_rmse, 4), ", p=", format_p(tests$p_ecdf_rmse),
    "\nKS=", format_num(tests$observed_ecdf_ks, 4), ", p=", format_p(tests$p_ecdf_ks),
    "\n", ifelse(tests$estimable, tests$permutation_mode, "not estimable")
  )
  write_csv(curves, file.path(cfg$dir, "tables", paste0(cfg$prefix, "_direct_group_ecdf_comparisons_9panel_curves.csv")))
  write_csv(tests, file.path(cfg$dir, "stats", paste0(cfg$prefix, "_direct_group_ecdf_comparisons_9panel_tests.csv")))

  panel_levels <- paste0(seq_len(9), ". ", c(
    "0 vs treated", "0 vs treated, ETP-stratified", "0 vs 30 mg/kg", "0 vs 120 mg/kg", "30 vs 120 mg/kg",
    "Treated: ETP-lower vs ETP-higher", "Control: ETP-lower vs ETP-higher",
    "ETP-higher: 0 vs treated", "ETP-lower: 0 vs treated"
  ))
  curves$panel <- factor(curves$panel, levels = panel_levels)
  tests$panel <- factor(tests$panel, levels = panel_levels)
  direct_colors <- c("0mg/kg" = "#666666", "treated" = "#1b9e77", "30mg/kg" = "#d95f02", "120mg/kg" = "#7570b3")
  p <- ggplot(curves, aes(x = pseudotime, y = mean_ecdf, color = color_group, linetype = line_group, group = curve_label)) +
    geom_line(linewidth = 0.8) +
    geom_text(
      data = tests,
      aes(x = 0.98, y = 0.05, label = annotation),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0,
      size = 2.35,
      lineheight = 0.92
    ) +
    facet_wrap(~ panel, ncol = 3, drop = FALSE) +
    scale_color_manual(values = direct_colors, name = "group") +
    scale_linetype_manual(values = c("All" = "solid", "ETP-lower" = "solid", "ETP-higher" = "22"), name = "EndTimePoint ploidy") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1.05), clip = "off") +
    labs(
      title = paste0(name, ": direct group mean ECDF comparisons (", analysis_version, ")"),
      subtitle = paste0("Equal-sample mean ECDFs; fixed sample-level ETP threshold = ", format_num(endpoint_threshold, 6)),
      x = "pseudotime",
      y = "mean ECDF"
    ) +
    plot_theme +
    theme(legend.position = "bottom")
  save_plot(p, file.path(cfg$dir, "figures"), paste0(cfg$prefix, "_direct_group_ecdf_comparisons"), 15, 11.5)
  list(curves = curves, tests = tests)
}

leave_one_out_primary <- function(equal_shift) {
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  full_p <- safe_cor_result(treated$ecdf_rmse, treated$TGI_percent_auc, "pearson")
  full_s <- safe_cor_result(treated$ecdf_rmse, treated$TGI_percent_auc, "spearman")
  rows <- lapply(treated$sample_id, function(sample_id) {
    local <- treated[treated$sample_id != sample_id, , drop = FALSE]
    p <- safe_cor_result(local$ecdf_rmse, local$TGI_percent_auc, "pearson")
    s <- safe_cor_result(local$ecdf_rmse, local$TGI_percent_auc, "spearman")
    data.frame(
      omitted_sample_id = sample_id,
      n = p$n,
      pearson_r = p$estimate,
      pearson_p = p$p_value,
      spearman_rho = s$estimate,
      spearman_p = s$p_value,
      sign_matches_full_pearson = sign(p$estimate) == sign(full_p$estimate),
      delta_pearson_r_from_full = p$estimate - full_p$estimate,
      influential_flag = is.finite(p$estimate) && is.finite(full_p$estimate) &&
        (sign(p$estimate) != sign(full_p$estimate) || abs(p$estimate - full_p$estimate) > 0.20),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bootstrap_primary <- function(equal_shift, n_boot = 1000L) {
  treated <- equal_shift[equal_shift$dose_mg > 0, , drop = FALSE]
  rows <- lapply(seq_len(n_boot), function(i) {
    idx <- sample(seq_len(nrow(treated)), replace = TRUE)
    local <- treated[idx, , drop = FALSE]
    p <- safe_cor_result(local$ecdf_rmse, local$TGI_percent_auc, "pearson")
    s <- safe_cor_result(local$ecdf_rmse, local$TGI_percent_auc, "spearman")
    data.frame(iteration = i, n = p$n, pearson_r = p$estimate, pearson_p = p$p_value, spearman_rho = s$estimate, spearman_p = s$p_value, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

summarize_bootstrap <- function(boot, loo, full_assoc) {
  q <- function(x, p) as.numeric(quantile(x[is.finite(x)], probs = p, na.rm = TRUE, names = FALSE))
  primary <- full_assoc[full_assoc$pre_specified_primary & full_assoc$reference_type == "primary_equal_sample_reference", , drop = FALSE]
  data.frame(
    full_n = if (nrow(primary) > 0) primary$n[1] else NA_integer_,
    full_pearson_r = if (nrow(primary) > 0) primary$pearson_r[1] else NA_real_,
    full_pearson_perm_p = if (nrow(primary) > 0) primary$pearson_p_permutation_two_sided[1] else NA_real_,
    full_spearman_rho = if (nrow(primary) > 0) primary$spearman_rho[1] else NA_real_,
    full_spearman_perm_p = if (nrow(primary) > 0) primary$spearman_p_permutation_two_sided[1] else NA_real_,
    leave_one_out_all_positive = all(loo$pearson_r > 0, na.rm = TRUE),
    leave_one_out_min_pearson_r = min(loo$pearson_r, na.rm = TRUE),
    leave_one_out_max_pearson_p = max(loo$pearson_p, na.rm = TRUE),
    n_influential_flags = sum(loo$influential_flag, na.rm = TRUE),
    bootstrap_pearson_r_ci025 = q(boot$pearson_r, 0.025),
    bootstrap_pearson_r_ci975 = q(boot$pearson_r, 0.975),
    bootstrap_spearman_rho_ci025 = q(boot$spearman_rho, 0.025),
    bootstrap_spearman_rho_ci975 = q(boot$spearman_rho, 0.975),
    stringsAsFactors = FALSE
  )
}

plot_robustness <- function(loo, boot) {
  p1 <- ggplot(loo, aes(x = reorder(omitted_sample_id, pearson_r), y = pearson_r, color = influential_flag)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(size = 2.6) +
    coord_flip() +
    scale_color_manual(values = c("FALSE" = "#377eb8", "TRUE" = "#e41a1c"), name = "influential") +
    labs(title = "CellCycle primary TGI association: leave-one-out Pearson r", x = "omitted sample", y = "Pearson r") +
    plot_theme
  save_square_plot(p1, dirs$Robustness, "CellCycle_primary_TGI_leave_one_out_forest", 6.8, 6.8)
  p2 <- ggplot(boot, aes(x = pearson_r)) +
    geom_histogram(bins = 35, fill = "#377eb8", color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, color = "grey40", linewidth = 0.35) +
    labs(title = "CellCycle primary TGI association: bootstrap Pearson r", x = "bootstrap Pearson r", y = "iterations") +
    plot_theme
  save_square_plot(p2, dirs$Robustness, "CellCycle_primary_TGI_bootstrap_distribution", 6.6, 6.6)
}

downsample_sensitivity <- function(name, df, sample_meta, n_iter = 200L, n_perm_iter = 1000L) {
  counts <- table(df$sample_id)
  target_n <- as.integer(min(counts))
  retained <- names(counts)[counts >= target_n]
  exclusions <- data.frame(
    analysis = paste0("downsampling_", name),
    sample_id = names(counts),
    reason = ifelse(counts < target_n, paste0("cell_count_below_target_", target_n), "retained"),
    cell_count = as.integer(counts),
    target_n = target_n,
    stringsAsFactors = FALSE
  )
  iter_rows <- lapply(seq_len(n_iter), function(i) {
    sampled <- do.call(rbind, lapply(retained, function(sample_id) {
      sub <- df[df$sample_id == sample_id, , drop = FALSE]
      sub[sample(seq_len(nrow(sub)), target_n, replace = FALSE), , drop = FALSE]
    }))
    smeta <- build_sample_meta(sampled)
    sh <- calculate_shift_metrics(sampled, smeta, "primary_equal_sample_reference")
    dt <- run_dose_tests(sampled, n_perm_iter)
    treated <- sh[sh$dose_mg > 0, , drop = FALSE]
    pearson <- permutation_cor_p(treated$ecdf_rmse, treated$TGI_percent_auc, "pearson", n_perm_iter, exact = FALSE)
    spearman <- permutation_cor_p(treated$ecdf_rmse, treated$TGI_percent_auc, "spearman", n_perm_iter, exact = FALSE)
    combined <- dt[dt$comparison == "0_vs_30plus120" & !dt$stratified_by_end_timepoint_ploidy_group, , drop = FALSE]
    data.frame(
      iteration = i,
      target_cells_per_sample = target_n,
      n_samples_retained = length(retained),
      treatment_p_ecdf_rmse = if (nrow(combined) > 0) combined$p_ecdf_rmse[1] else NA_real_,
      pearson_r = pearson$estimate,
      pearson_p_permutation = pearson$permutation_p_two_sided,
      spearman_rho = spearman$estimate,
      stringsAsFactors = FALSE
    )
  })
  iterations <- do.call(rbind, iter_rows)
  summary <- data.frame(
    compartment = name,
    target_cells_per_sample = target_n,
    n_iterations = n_iter,
    n_samples_retained = length(retained),
    median_treatment_p_ecdf_rmse = median(iterations$treatment_p_ecdf_rmse, na.rm = TRUE),
    prop_treatment_p_lt_0_05 = mean(iterations$treatment_p_ecdf_rmse < 0.05, na.rm = TRUE),
    median_pearson_r = median(iterations$pearson_r, na.rm = TRUE),
    pearson_r_ci025 = as.numeric(quantile(iterations$pearson_r, 0.025, na.rm = TRUE)),
    pearson_r_ci975 = as.numeric(quantile(iterations$pearson_r, 0.975, na.rm = TRUE)),
    median_spearman_rho = median(iterations$spearman_rho, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  p_r <- ggplot(iterations, aes(x = pearson_r)) +
    geom_histogram(bins = 35, fill = "#377eb8", color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0, color = "grey45", linewidth = 0.35) +
    labs(title = paste0(name, ": downsampled TGI Pearson r"), x = "Pearson r", y = "iterations") +
    plot_theme
  save_square_plot(p_r, dirs$Robustness, paste0("downsampling_", name, "_TGI_r_distribution"), 6.4, 6.4)
  p_p <- ggplot(iterations, aes(x = treatment_p_ecdf_rmse)) +
    geom_histogram(bins = 35, fill = "#4daf4a", color = "white", alpha = 0.85) +
    geom_vline(xintercept = 0.05, color = "#e41a1c", linewidth = 0.45) +
    labs(title = paste0(name, ": downsampled treatment ECDF p-values"), x = "ECDF RMSE permutation p", y = "iterations") +
    plot_theme
  save_square_plot(p_p, dirs$Robustness, paste0("downsampling_", name, "_treatment_p_distribution"), 6.4, 6.4)
  list(iterations = iterations, summary = summary, exclusions = exclusions)
}

run_confounding <- function(results) {
  ploidy_rows <- list()
  dose_rows <- list()
  resid_rows <- list()
  for (name in names(results)) {
    sh <- results[[name]]$equal_shift
    for (sample_set in c("all", "treated")) {
      local <- if (sample_set == "treated") sh[sh$dose_mg > 0, , drop = FALSE] else sh
      for (metric in c("ecdf_rmse", "signed_mean_shift")) {
        for (ploidy_col in c("mean_cell_ploidy", "median_cell_ploidy", "p90_cell_ploidy", "end_timepoint_ploidy_group_numeric")) {
          for (method in c("pearson", "spearman")) {
            cr <- safe_cor_result(local[[metric]], local[[ploidy_col]], method)
            ploidy_rows[[length(ploidy_rows) + 1L]] <- data.frame(
              compartment = name, sample_set = sample_set, shift_metric = metric, ploidy_measure = ploidy_col,
              method = method, n = cr$n, estimate = cr$estimate, p_value = cr$p_value, stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    treated <- sh[sh$dose_mg > 0, , drop = FALSE]
    for (value_col in c("ecdf_rmse", "signed_mean_shift", "TGI_percent_auc")) {
      wt <- safe_wilcox(treated[[value_col]], treated$dose)
      tt <- safe_ttest(treated[[value_col]], treated$dose)
      dose_rows[[length(dose_rows) + 1L]] <- data.frame(compartment = name, value = value_col, comparison = "30_vs_120", test = "wilcoxon", n = wt$n, statistic = wt$statistic, p_value = wt$p_value, stringsAsFactors = FALSE)
      dose_rows[[length(dose_rows) + 1L]] <- data.frame(compartment = name, value = value_col, comparison = "30_vs_120", test = "t_test", n = tt$n, statistic = tt$statistic, p_value = tt$p_value, stringsAsFactors = FALSE)
    }
    if (nrow(treated) >= 6) {
      local <- treated[is.finite(treated$ecdf_rmse) & is.finite(treated$TGI_percent_auc), , drop = FALSE]
      if (nrow(local) >= 6 && length(unique(local$dose_mg)) > 1) {
        shift_resid <- resid(lm(ecdf_rmse ~ factor(dose_mg) + end_timepoint_ploidy_group_numeric, data = local))
        tgi_resid <- resid(lm(TGI_percent_auc ~ factor(dose_mg) + end_timepoint_ploidy_group_numeric, data = local))
        cr <- safe_cor_result(shift_resid, tgi_resid, "pearson")
        resid_rows[[length(resid_rows) + 1L]] <- data.frame(compartment = name, analysis = "residualized_against_dose_and_end_timepoint_ploidy_group", method = "pearson", n = cr$n, estimate = cr$estimate, p_value = cr$p_value, stringsAsFactors = FALSE)
      }
      centered <- local
      centered$shift_centered <- centered$ecdf_rmse - ave(centered$ecdf_rmse, centered$dose_mg, FUN = function(x) mean(x, na.rm = TRUE))
      centered$tgi_centered <- centered$TGI_percent_auc - ave(centered$TGI_percent_auc, centered$dose_mg, FUN = function(x) mean(x, na.rm = TRUE))
      cr <- safe_cor_result(centered$shift_centered, centered$tgi_centered, "pearson")
      resid_rows[[length(resid_rows) + 1L]] <- data.frame(compartment = name, analysis = "within_dose_centered", method = "pearson", n = cr$n, estimate = cr$estimate, p_value = cr$p_value, stringsAsFactors = FALSE)
    }
  }
  ploidy <- do.call(rbind, ploidy_rows)
  dose <- do.call(rbind, dose_rows)
  resid <- do.call(rbind, resid_rows)
  write_csv(ploidy, file.path(dirs$Robustness, "ploidy_confounding_tests.csv"))
  write_csv(dose, file.path(dirs$Robustness, "dose_confounding_tests.csv"))
  write_csv(resid, file.path(dirs$Robustness, "residualized_TGI_associations.csv"))
  cc <- results$CellCycle$equal_shift
  p1 <- ggplot(cc, aes(x = mean_cell_ploidy, y = ecdf_rmse, color = dose)) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.5) +
    geom_point(aes(shape = end_timepoint_ploidy_group), size = 2.6) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% unique(cc$dose)], name = "dose") +
    labs(title = "CellCycle ECDF RMSE vs mean cell ploidy", x = "mean cell ploidy", y = "ECDF RMSE") +
    plot_theme
  save_square_plot(p1, dirs$Robustness, "CellCycle_ecdf_rmse_vs_ploidy", 6.4, 6.4)
  treated <- cc[cc$dose_mg > 0, , drop = FALSE]
  treated$shift_centered <- treated$ecdf_rmse - ave(treated$ecdf_rmse, treated$dose_mg, FUN = function(x) mean(x, na.rm = TRUE))
  treated$tgi_centered <- treated$TGI_percent_auc - ave(treated$TGI_percent_auc, treated$dose_mg, FUN = function(x) mean(x, na.rm = TRUE))
  p2 <- ggplot(treated, aes(x = shift_centered, y = tgi_centered, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
    geom_point(aes(shape = end_timepoint_ploidy_group), size = 2.7) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% unique(treated$dose)], name = "dose") +
    labs(title = "CellCycle TGI association after within-dose centering", x = "dose-centered ECDF RMSE", y = "dose-centered AUC TGI") +
    plot_theme
  save_square_plot(p2, dirs$Robustness, "CellCycle_TGI_association_within_dose_centered", 6.6, 6.6)
  list(ploidy = ploidy, dose = dose, resid = resid)
}

permutation_mean_diff_p <- function(value, group, group_a, group_b, n_perm = 10000L, exact = TRUE) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]
  group <- as.character(group[ok])
  group_a <- as.character(group_a)
  group_b <- as.character(group_b)
  keep <- group %in% c(group_a, group_b)
  value <- value[keep]
  group <- group[keep]
  n <- length(value)
  n_a <- sum(group == group_a)
  n_b <- sum(group == group_b)
  if (n < 3 || n_a < 2 || n_b < 2) {
    return(data.frame(
      n = n, n_group_a = n_a, n_group_b = n_b, observed_mean_diff_group_b_minus_a = NA_real_,
      permutation_p_two_sided = NA_real_, permutation_p_group_b_greater = NA_real_,
      n_permutations = NA_integer_, permutation_mode = NA_character_, stringsAsFactors = FALSE
    ))
  }
  observed <- mean(value[group == group_b], na.rm = TRUE) - mean(value[group == group_a], na.rm = TRUE)
  n_exact <- choose(n, n_a)
  if (isTRUE(exact) && is.finite(n_exact) && n_exact <= n_perm) {
    combos <- combn(seq_len(n), n_a, simplify = FALSE)
    perm <- vapply(combos, function(idx) {
      g <- rep(group_b, n)
      g[idx] <- group_a
      mean(value[g == group_b], na.rm = TRUE) - mean(value[g == group_a], na.rm = TRUE)
    }, numeric(1))
    mode <- "exact_dose_label_enumeration"
  } else {
    perm <- replicate(n_perm, {
      idx <- sample(seq_len(n), n_a)
      g <- rep(group_b, n)
      g[idx] <- group_a
      mean(value[g == group_b], na.rm = TRUE) - mean(value[g == group_a], na.rm = TRUE)
    })
    mode <- "monte_carlo_dose_label_permutation"
  }
  data.frame(
    n = n,
    n_group_a = n_a,
    n_group_b = n_b,
    observed_mean_diff_group_b_minus_a = observed,
    permutation_p_two_sided = mean(abs(perm) >= abs(observed) - 1e-15, na.rm = TRUE),
    permutation_p_group_b_greater = mean(perm >= observed - 1e-15, na.rm = TRUE),
    n_permutations = length(perm),
    permutation_mode = mode,
    stringsAsFactors = FALSE
  )
}

dose_value_summary <- function(value, group, group_value) {
  x <- value[as.character(group) == as.character(group_value)]
  data.frame(
    n = sum(is.finite(x)),
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_endpoint_tgi_long <- function(treated, endpoint_cols) {
  if (length(endpoint_cols) == 0) return(data.frame())
  do.call(rbind, lapply(endpoint_cols, function(col) {
    data.frame(
      sample_id = treated$sample_id,
      end_timepoint_ploidy_group = treated$end_timepoint_ploidy_group,
      dose = treated$dose,
      dose_mg = treated$dose_mg,
      day = endpoint_day_num(col),
      tgi_measure = col,
      TGI_percent_endpoint = treated[[col]],
      stringsAsFactors = FALSE
    )
  }))
}

plot_dose_specific_tgi <- function(name, treated, comp_dir, tgi_comparisons) {
  fig_dir <- ensure_dir(file.path(comp_dir, "figures"))
  dose_levels <- unique(treated$dose[order(treated$dose_mg)])
  comp_prefix <- ifelse(name == "CellCycle", "CellCycle", "NonCellCycle")

  p_auc <- ggplot(treated, aes(x = dose, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_boxplot(aes(group = dose), outlier.shape = NA, color = "grey55", linewidth = 0.45) +
    geom_point(aes(shape = end_timepoint_ploidy_group), position = position_jitter(width = 0.08, height = 0), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.35, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = paste0(comp_prefix, ": AUC TGI by gemcitabine dose"), x = "dose", y = "AUC-based TGI (%)") +
    plot_theme
  save_square_plot(p_auc, fig_dir, paste0(comp_prefix, "_TGI_AUC_by_treated_dose"), 6.6, 6.6)

  p_shift <- ggplot(treated, aes(x = ecdf_rmse, y = TGI_percent_auc, color = dose)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
    geom_point(aes(shape = end_timepoint_ploidy_group), size = 2.8) +
    geom_text(aes(label = sample_id), size = 2.35, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
    facet_wrap(~ dose, nrow = 1) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(
      title = paste0(comp_prefix, ": shift-TGI association within each treated dose"),
      x = "ECDF RMSE from ploidy-matched 0 mg/kg reference",
      y = "AUC-based TGI (%)"
    ) +
    plot_theme
  save_square_plot(p_shift, fig_dir, paste0(comp_prefix, "_shift_TGI_AUC_by_treated_dose_faceted"), 9.2, 5.8)

  endpoint_cols <- endpoint_tgi_cols(names(treated))
  endpoint_long <- make_endpoint_tgi_long(treated, endpoint_cols)
  if (nrow(endpoint_long) > 0) {
    p_endpoint <- ggplot(endpoint_long, aes(x = day, y = TGI_percent_endpoint, color = dose, group = sample_id)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_line(alpha = 0.6, linewidth = 0.55) +
      geom_point(aes(shape = end_timepoint_ploidy_group), size = 2.1) +
      scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
      labs(title = paste0(comp_prefix, ": endpoint TGI across all measured days"), x = "day", y = "endpoint TGI (%)") +
      plot_theme
    save_plot(p_endpoint, fig_dir, paste0(comp_prefix, "_endpoint_TGI_by_treated_dose_over_days"), 7.4, 4.8)
    write_csv(endpoint_long, file.path(comp_dir, "tables", paste0(comp_prefix, "_endpoint_TGI_long_by_sample.csv")))
  }
  invisible(TRUE)
}

run_dose_specific_tgi_exploration <- function(results, n_perm = 10000L) {
  ensure_dir(dirs$DoseSpecific)
  all_comparisons <- list()
  all_correlations <- list()
  all_interactions <- list()
  all_slopes <- list()

  for (name in names(results)) {
    sh <- results[[name]]$equal_shift
    treated <- sh[sh$dose_mg > 0, , drop = FALSE]
    comp_dir <- ensure_dir(file.path(dirs$DoseSpecific, name))
    ensure_dir(file.path(comp_dir, "tables"))
    ensure_dir(file.path(comp_dir, "stats"))
    if (nrow(treated) == 0) next

    dose_values <- sort(unique(treated$dose_mg[is.finite(treated$dose_mg)]))
    tgi_cols <- unique(c("TGI_percent_auc", endpoint_tgi_cols(names(treated))))
    shift_cols <- c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")

    comp_comparisons <- list()
    if (length(dose_values) >= 2) {
      dose_pairs <- combn(dose_values, 2)
      for (i in seq_len(ncol(dose_pairs))) {
        dose_a <- dose_pairs[1, i]
        dose_b <- dose_pairs[2, i]
        local_pair <- treated[treated$dose_mg %in% c(dose_a, dose_b), , drop = FALSE]
        dose_a_label <- unique(local_pair$dose[local_pair$dose_mg == dose_a])[1]
        dose_b_label <- unique(local_pair$dose[local_pair$dose_mg == dose_b])[1]
        for (tgi_col in tgi_cols) {
          a_sum <- dose_value_summary(local_pair[[tgi_col]], local_pair$dose_mg, dose_a)
          b_sum <- dose_value_summary(local_pair[[tgi_col]], local_pair$dose_mg, dose_b)
          wt <- safe_wilcox(local_pair[[tgi_col]], local_pair$dose_mg)
          tt <- safe_ttest(local_pair[[tgi_col]], local_pair$dose_mg)
          perm <- permutation_mean_diff_p(local_pair[[tgi_col]], local_pair$dose_mg, dose_a, dose_b, n_perm, exact = TRUE)
          comp_comparisons[[length(comp_comparisons) + 1L]] <- data.frame(
            compartment = name,
            comparison = paste0(dose_a, "_vs_", dose_b),
            dose_a = dose_a_label,
            dose_b = dose_b_label,
            dose_a_mg = dose_a,
            dose_b_mg = dose_b,
            tgi_measure = tgi_col,
            tgi_day = ifelse(grepl("^TGI_percent_Day_", tgi_col), endpoint_day_num(tgi_col), NA_real_),
            n_a = a_sum$n,
            n_b = b_sum$n,
            mean_a = a_sum$mean,
            mean_b = b_sum$mean,
            median_a = a_sum$median,
            median_b = b_sum$median,
            sd_a = a_sum$sd,
            sd_b = b_sum$sd,
            mean_diff_b_minus_a = perm$observed_mean_diff_group_b_minus_a,
            wilcoxon_p = wt$p_value,
            t_test_p = tt$p_value,
            permutation_p_two_sided = perm$permutation_p_two_sided,
            permutation_p_group_b_greater = perm$permutation_p_group_b_greater,
            n_permutations = perm$n_permutations,
            permutation_mode = perm$permutation_mode,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    comp_comparisons <- if (length(comp_comparisons) > 0) do.call(rbind, comp_comparisons) else data.frame()
    if (nrow(comp_comparisons) > 0) {
      comp_comparisons$permutation_q_exploratory <- p.adjust(comp_comparisons$permutation_p_two_sided, method = "BH")
      write_csv(comp_comparisons, file.path(comp_dir, "stats", paste0(name, "_TGI_AUC_endpoint_30_vs_120_comparisons.csv")))
      all_comparisons[[length(all_comparisons) + 1L]] <- comp_comparisons
    }

    comp_correlations <- list()
    for (dose_value in dose_values) {
      local <- treated[treated$dose_mg == dose_value, , drop = FALSE]
      dose_label <- unique(local$dose)[1]
      for (shift_col in shift_cols) {
        for (tgi_col in tgi_cols) {
          pearson <- permutation_cor_p(local[[shift_col]], local[[tgi_col]], "pearson", n_perm, exact = TRUE)
          spearman <- permutation_cor_p(local[[shift_col]], local[[tgi_col]], "spearman", n_perm, exact = TRUE)
          comp_correlations[[length(comp_correlations) + 1L]] <- data.frame(
            compartment = name,
            dose = dose_label,
            dose_mg = dose_value,
            reference_type = "primary_equal_sample_reference",
            shift_metric = shift_col,
            tgi_measure = tgi_col,
            tgi_day = ifelse(grepl("^TGI_percent_Day_", tgi_col), endpoint_day_num(tgi_col), NA_real_),
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
            pre_specified_primary = name == "CellCycle" && shift_col == "ecdf_rmse" && tgi_col == "TGI_percent_auc",
            stringsAsFactors = FALSE
          )
        }
      }
    }
    comp_correlations <- if (length(comp_correlations) > 0) do.call(rbind, comp_correlations) else data.frame()
    if (nrow(comp_correlations) > 0) {
      comp_correlations$pearson_q_exploratory <- p.adjust(comp_correlations$pearson_p_permutation_two_sided, method = "BH")
      comp_correlations$spearman_q_exploratory <- p.adjust(comp_correlations$spearman_p_permutation_two_sided, method = "BH")
      write_csv(comp_correlations, file.path(comp_dir, "stats", paste0(name, "_dose_specific_shift_TGI_correlations.csv")))
      all_correlations[[length(all_correlations) + 1L]] <- comp_correlations
    }

    local_primary <- treated[is.finite(treated$ecdf_rmse) & is.finite(treated$TGI_percent_auc), , drop = FALSE]
    if (nrow(local_primary) >= 6 && length(unique(local_primary$dose_mg)) > 1 && length(unique(local_primary$ecdf_rmse)) > 1) {
      fit <- try(lm(TGI_percent_auc ~ scale(ecdf_rmse) * factor(dose_mg), data = local_primary), silent = TRUE)
      if (!inherits(fit, "try-error")) {
        cf <- coef(summary(fit))
        term_rows <- data.frame(
          compartment = name,
          analysis = "treated_only_AUC_TGI_by_shift_dose_interaction",
          model = "TGI_percent_auc ~ scale(ecdf_rmse) * factor(dose_mg)",
          term = rownames(cf),
          estimate = cf[, "Estimate"],
          std_error = cf[, "Std. Error"],
          statistic = cf[, "t value"],
          p_value = cf[, "Pr(>|t|)"],
          n = nrow(local_primary),
          stringsAsFactors = FALSE
        )
        write_csv(term_rows, file.path(comp_dir, "stats", paste0(name, "_shift_TGI_AUC_dose_interaction_lm.csv")))
        all_interactions[[length(all_interactions) + 1L]] <- term_rows
      }
    }

    comp_slopes <- list()
    for (dose_value in dose_values) {
      local <- local_primary[local_primary$dose_mg == dose_value, , drop = FALSE]
      dose_label <- unique(local$dose)[1]
      if (nrow(local) >= 3 && length(unique(local$ecdf_rmse)) > 1) {
        fit <- try(lm(TGI_percent_auc ~ ecdf_rmse, data = local), silent = TRUE)
        if (!inherits(fit, "try-error")) {
          cf <- coef(summary(fit))
          if ("ecdf_rmse" %in% rownames(cf)) {
            comp_slopes[[length(comp_slopes) + 1L]] <- data.frame(
              compartment = name,
              dose = dose_label,
              dose_mg = dose_value,
              analysis = "within_dose_AUC_TGI_slope",
              model = "TGI_percent_auc ~ ecdf_rmse",
              n = nrow(local),
              slope = cf["ecdf_rmse", "Estimate"],
              std_error = cf["ecdf_rmse", "Std. Error"],
              statistic = cf["ecdf_rmse", "t value"],
              p_value = cf["ecdf_rmse", "Pr(>|t|)"],
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    comp_slopes <- if (length(comp_slopes) > 0) do.call(rbind, comp_slopes) else data.frame()
    if (nrow(comp_slopes) > 0) {
      write_csv(comp_slopes, file.path(comp_dir, "stats", paste0(name, "_within_dose_shift_TGI_AUC_slopes.csv")))
      all_slopes[[length(all_slopes) + 1L]] <- comp_slopes
    }

    plot_dose_specific_tgi(name, treated, comp_dir, comp_comparisons)
    summary_lines <- c(
      "# Dose-specific TGI Exploration",
      "",
      "This exploratory module separates treated tumors by gemcitabine dose. It asks whether AUC/endpoint TGI differs between treated doses and whether the shift-TGI association has the same direction within each dose stratum.",
      "",
      "The inference unit is the biological sample. Per-dose correlations have four samples per dose in the current dataset, so they should be interpreted as descriptive direction checks rather than definitive dose-specific effects.",
      "",
      "Primary outputs:",
      "",
      paste0("- `stats/", name, "_TGI_AUC_endpoint_30_vs_120_comparisons.csv`: AUC and all available endpoint TGI 30 vs 120 mg/kg comparisons."),
      paste0("- `stats/", name, "_dose_specific_shift_TGI_correlations.csv`: within-dose shift-TGI correlations for AUC and all endpoint TGI measures."),
      paste0("- `stats/", name, "_shift_TGI_AUC_dose_interaction_lm.csv`: treated-only AUC TGI interaction sensitivity model."),
      paste0("- `figures/", name, "_TGI_AUC_by_treated_dose.*`: AUC TGI by dose."),
      paste0("- `figures/", name, "_shift_TGI_AUC_by_treated_dose_faceted.*`: shift-TGI relationship within each dose."),
      paste0("- `figures/", name, "_endpoint_TGI_by_treated_dose_over_days.*`: endpoint TGI over all measured days.")
    )
    write_text(summary_lines, file.path(comp_dir, "dose_specific_TGI_summary.md"))
  }

  comparisons <- if (length(all_comparisons) > 0) do.call(rbind, all_comparisons) else data.frame()
  correlations <- if (length(all_correlations) > 0) do.call(rbind, all_correlations) else data.frame()
  interactions <- if (length(all_interactions) > 0) do.call(rbind, all_interactions) else data.frame()
  slopes <- if (length(all_slopes) > 0) do.call(rbind, all_slopes) else data.frame()
  if (nrow(comparisons) > 0) write_csv(comparisons, file.path(dirs$DoseSpecific, "all_compartments_TGI_AUC_endpoint_dose_comparisons.csv"))
  if (nrow(correlations) > 0) write_csv(correlations, file.path(dirs$DoseSpecific, "all_compartments_dose_specific_shift_TGI_correlations.csv"))
  if (nrow(interactions) > 0) write_csv(interactions, file.path(dirs$DoseSpecific, "all_compartments_shift_TGI_AUC_dose_interaction_lm.csv"))
  if (nrow(slopes) > 0) write_csv(slopes, file.path(dirs$DoseSpecific, "all_compartments_within_dose_shift_TGI_AUC_slopes.csv"))
  list(comparisons = comparisons, correlations = correlations, interactions = interactions, slopes = slopes)
}

run_auc_tgi_association_analysis <- function(results) {
  ensure_dir(dirs$AUCTGIAssociation)
  all_models <- list()
  all_summaries <- list()
  all_mean_etp_models <- list()
  all_mean_etp_summaries <- list()
  all_samples <- list()

  for (name in names(results)) {
    comp_dir <- ensure_dir(file.path(dirs$AUCTGIAssociation, name))
    ensure_dir(file.path(comp_dir, "tables"))
    ensure_dir(file.path(comp_dir, "stats"))
    ensure_dir(file.path(comp_dir, "figures"))

    treated <- results[[name]]$equal_shift[results[[name]]$equal_shift$dose_mg > 0, , drop = FALSE]
    keep <- c(
      "sample_id", "end_timepoint_ploidy_group", "end_timepoint_ploidy_group_numeric", "dose", "dose_mg",
      "sample_mean_endpoint_ploidy",
      "ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift",
      "TGI_percent_auc", endpoint_tgi_cols(names(treated))
    )
    keep <- intersect(keep, names(treated))
    treated <- treated[, keep, drop = FALSE]
    treated$compartment <- name
    treated$end_timepoint_ploidy_group <- factor(treated$end_timepoint_ploidy_group, levels = unique(treated$end_timepoint_ploidy_group[order(treated$end_timepoint_ploidy_group_numeric)]))
    treated$dose <- factor(treated$dose, levels = unique(treated$dose[order(treated$dose_mg)]))
    treated$dose_mg_factor <- factor(treated$dose_mg, levels = sort(unique(treated$dose_mg)))
    write_csv(treated, file.path(comp_dir, "tables", paste0(name, "_treated_AUC_TGI_shift_ploidy_sample_table.csv")))
    all_samples[[length(all_samples) + 1L]] <- treated

    local <- treated[is.finite(treated$TGI_percent_auc) & is.finite(treated$ecdf_rmse), , drop = FALSE]
    local_mean_etp <- treated[
      is.finite(treated$TGI_percent_auc) & is.finite(treated$sample_mean_endpoint_ploidy),
      ,
      drop = FALSE
    ]
    model_rows <- list()
    model_summaries <- list()
    fits <- list(
      shift_only = try(lm(TGI_percent_auc ~ scale(ecdf_rmse), data = local), silent = TRUE),
      shift_ploidy_dose_adjusted = try(lm(TGI_percent_auc ~ scale(ecdf_rmse) + end_timepoint_ploidy_group + dose_mg_factor, data = local), silent = TRUE),
      shift_by_ploidy_dose_adjusted = try(lm(TGI_percent_auc ~ scale(ecdf_rmse) * end_timepoint_ploidy_group + dose_mg_factor, data = local), silent = TRUE)
    )
    labels <- c(
      shift_only = "TGI_percent_auc ~ scale(ecdf_rmse)",
      shift_ploidy_dose_adjusted = "TGI_percent_auc ~ scale(ecdf_rmse) + end_timepoint_ploidy_group + factor(dose_mg)",
      shift_by_ploidy_dose_adjusted = "TGI_percent_auc ~ scale(ecdf_rmse) * end_timepoint_ploidy_group + factor(dose_mg)"
    )
    for (analysis in names(fits)) {
      fit <- fits[[analysis]]
      if (!inherits(fit, "try-error") && !is.null(fit)) {
        model_rows[[length(model_rows) + 1L]] <- lm_terms_table(fit, name, analysis, labels[[analysis]], nrow(local))
        sm <- summary(fit)
        model_summaries[[length(model_summaries) + 1L]] <- data.frame(
          compartment = name,
          analysis = analysis,
          model = labels[[analysis]],
          n = nobs(fit),
          residual_df = stats::df.residual(fit),
          r_squared = unname(sm$r.squared),
          adjusted_r_squared = unname(sm$adj.r.squared),
          sigma = unname(sm$sigma),
          aic = AIC(fit),
          stringsAsFactors = FALSE
        )
      }
    }
    models <- if (length(model_rows) > 0) do.call(rbind, model_rows) else data.frame()
    summaries <- if (length(model_summaries) > 0) do.call(rbind, model_summaries) else data.frame()
    if (nrow(models) > 0) {
      models$pre_specified_primary <- models$compartment == "CellCycle" &
        models$analysis == "shift_ploidy_dose_adjusted" &
        models$term == "scale(ecdf_rmse)"
      write_csv(models, file.path(comp_dir, "stats", paste0(name, "_AUC_TGI_shift_ploidy_dose_models.csv")))
      all_models[[length(all_models) + 1L]] <- models
    }
    if (nrow(summaries) > 0) {
      write_csv(summaries, file.path(comp_dir, "stats", paste0(name, "_AUC_TGI_shift_ploidy_dose_model_summaries.csv")))
      all_summaries[[length(all_summaries) + 1L]] <- summaries
    }

    mean_etp_model_rows <- list()
    mean_etp_model_summaries <- list()
    mean_etp_fits <- list(
      mean_etp_only = try(lm(TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy), data = local_mean_etp), silent = TRUE),
      mean_etp_ploidy_dose_adjusted = try(lm(TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) + end_timepoint_ploidy_group + dose_mg_factor, data = local_mean_etp), silent = TRUE),
      mean_etp_by_ploidy_dose_adjusted = try(lm(TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) * end_timepoint_ploidy_group + dose_mg_factor, data = local_mean_etp), silent = TRUE)
    )
    mean_etp_labels <- c(
      mean_etp_only = "TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy)",
      mean_etp_ploidy_dose_adjusted = "TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) + end_timepoint_ploidy_group + factor(dose_mg)",
      mean_etp_by_ploidy_dose_adjusted = "TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) * end_timepoint_ploidy_group + factor(dose_mg)"
    )
    for (analysis in names(mean_etp_fits)) {
      fit <- mean_etp_fits[[analysis]]
      if (!inherits(fit, "try-error") && !is.null(fit)) {
        mean_etp_model_rows[[length(mean_etp_model_rows) + 1L]] <- lm_terms_table(
          fit,
          name,
          analysis,
          mean_etp_labels[[analysis]],
          nrow(local_mean_etp)
        )
        sm <- summary(fit)
        mean_etp_model_summaries[[length(mean_etp_model_summaries) + 1L]] <- data.frame(
          compartment = name,
          analysis = analysis,
          model = mean_etp_labels[[analysis]],
          n = nobs(fit),
          residual_df = stats::df.residual(fit),
          r_squared = unname(sm$r.squared),
          adjusted_r_squared = unname(sm$adj.r.squared),
          sigma = unname(sm$sigma),
          aic = AIC(fit),
          stringsAsFactors = FALSE
        )
      }
    }
    mean_etp_models <- if (length(mean_etp_model_rows) > 0) do.call(rbind, mean_etp_model_rows) else data.frame()
    mean_etp_summaries <- if (length(mean_etp_model_summaries) > 0) do.call(rbind, mean_etp_model_summaries) else data.frame()
    if (nrow(mean_etp_models) > 0) {
      mean_etp_models$pre_specified_primary <- FALSE
      write_csv(mean_etp_models, file.path(comp_dir, "stats", paste0(name, "_AUC_TGI_mean_ETP_ploidy_dose_models.csv")))
      all_mean_etp_models[[length(all_mean_etp_models) + 1L]] <- mean_etp_models
    }
    if (nrow(mean_etp_summaries) > 0) {
      write_csv(mean_etp_summaries, file.path(comp_dir, "stats", paste0(name, "_AUC_TGI_mean_ETP_ploidy_dose_model_summaries.csv")))
      all_mean_etp_summaries[[length(all_mean_etp_summaries) + 1L]] <- mean_etp_summaries
    }

    if (nrow(local) > 0) {
      p <- ggplot(local, aes(x = ecdf_rmse, y = TGI_percent_auc, color = end_timepoint_ploidy_group, shape = dose)) +
        geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
        geom_point(size = 3) +
        geom_text(aes(label = sample_id), size = 2.4, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
        labs(
          title = paste0(name, ": AUC-TGI association with pseudotime shift and ploidy"),
          x = "ECDF RMSE from ploidy-matched 0 mg/kg reference",
          y = "AUC-based TGI (%)",
          color = "EndTimePoint ploidy",
          shape = "dose"
        ) +
        plot_theme
      save_square_plot(p, file.path(comp_dir, "figures"), paste0(name, "_AUC_TGI_vs_ecdf_rmse_by_ploidy_dose"), 6.8, 6.8)
    }

    if (nrow(local_mean_etp) > 0) {
      p_mean_etp <- ggplot(
        local_mean_etp,
        aes(
          x = sample_mean_endpoint_ploidy,
          y = TGI_percent_auc,
          color = end_timepoint_ploidy_group,
          shape = dose
        )
      ) +
        geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linewidth = 0.55) +
        geom_point(size = 3) +
        geom_text(aes(label = sample_id), size = 2.4, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
        labs(
          title = paste0(name, ": AUC-TGI association with sample mean ETP and ploidy"),
          x = "Sample mean ETP",
          y = "AUC-based TGI (%)",
          color = "EndTimePoint ploidy",
          shape = "dose"
        ) +
        plot_theme
      save_square_plot(
        p_mean_etp,
        file.path(comp_dir, "figures"),
        paste0(name, "_AUC_TGI_vs_mean_ETP_by_ploidy_dose"),
        6.8,
        6.8
      )
    }

    summary_lines <- c(
      "# AUC-TGI Association Summary",
      "",
      "This module addresses Question 1: whether overall tumor-growth inhibition is associated with pseudotime shift and EndTimePoint ploidy.",
      "",
      "The inference unit is the biological sample. The primary outcome is `TGI_percent_auc`, and the primary pseudotime-shift metric is sample-equal ploidy-matched `ecdf_rmse`.",
      "",
      "Primary models:",
      "",
      "- `TGI_percent_auc ~ scale(ecdf_rmse) + end_timepoint_ploidy_group + factor(dose_mg)`",
      "- `TGI_percent_auc ~ scale(ecdf_rmse) * end_timepoint_ploidy_group + factor(dose_mg)`",
      "",
      "The first model asks whether pseudotime shift and ploidy explain overall AUC-TGI after accounting for treated dose. The second model asks whether the shift-TGI relationship differs by EndTimePoint ploidy.",
      "",
      "Parallel exploratory mean-ETP models:",
      "",
      "- `TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy)`",
      "- `TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) + end_timepoint_ploidy_group + factor(dose_mg)`",
      "- `TGI_percent_auc ~ scale(sample_mean_endpoint_ploidy) * end_timepoint_ploidy_group + factor(dose_mg)`",
      "",
      paste0("Mean-ETP model terms: `stats/", name, "_AUC_TGI_mean_ETP_ploidy_dose_models.csv`."),
      paste0("Mean-ETP model summaries: `stats/", name, "_AUC_TGI_mean_ETP_ploidy_dose_model_summaries.csv`.")
    )
    write_text(summary_lines, file.path(comp_dir, "AUC_TGI_association_summary.md"))
  }

  models <- if (length(all_models) > 0) do.call(rbind, all_models) else data.frame()
  summaries <- if (length(all_summaries) > 0) do.call(rbind, all_summaries) else data.frame()
  mean_etp_models <- if (length(all_mean_etp_models) > 0) do.call(rbind, all_mean_etp_models) else data.frame()
  mean_etp_summaries <- if (length(all_mean_etp_summaries) > 0) do.call(rbind, all_mean_etp_summaries) else data.frame()
  samples <- if (length(all_samples) > 0) do.call(rbind, all_samples) else data.frame()
  if (nrow(models) > 0) write_csv(models, file.path(dirs$AUCTGIAssociation, "all_compartments_AUC_TGI_shift_ploidy_dose_models.csv"))
  if (nrow(summaries) > 0) write_csv(summaries, file.path(dirs$AUCTGIAssociation, "all_compartments_AUC_TGI_shift_ploidy_dose_model_summaries.csv"))
  if (nrow(mean_etp_models) > 0) write_csv(mean_etp_models, file.path(dirs$AUCTGIAssociation, "all_compartments_AUC_TGI_mean_ETP_ploidy_dose_models.csv"))
  if (nrow(mean_etp_summaries) > 0) write_csv(mean_etp_summaries, file.path(dirs$AUCTGIAssociation, "all_compartments_AUC_TGI_mean_ETP_ploidy_dose_model_summaries.csv"))
  if (nrow(samples) > 0) write_csv(samples, file.path(dirs$AUCTGIAssociation, "all_compartments_treated_AUC_TGI_shift_ploidy_sample_table.csv"))
  list(
    models = models,
    summaries = summaries,
    mean_etp_models = mean_etp_models,
    mean_etp_summaries = mean_etp_summaries,
    samples = samples
  )
}

mean_row_allow_duplicate <- function(mat, sample_ids) {
  sample_ids <- sample_ids[sample_ids %in% rownames(mat)]
  if (length(sample_ids) == 0) return(rep(NA_real_, ncol(mat)))
  colMeans(mat[sample_ids, , drop = FALSE], na.rm = TRUE)
}

pseudotime_delta_from_ids <- function(ecdf_mat, sample_meta, treated_ids, control_ids) {
  treated_ecdf <- mean_row_allow_duplicate(ecdf_mat, treated_ids)
  control_ecdf <- mean_row_allow_duplicate(ecdf_mat, control_ids)
  delta <- treated_ecdf - control_ecdf
  signed_mean_shift <- mean(sample_meta[treated_ids, "mean_pseudotime"], na.rm = TRUE) -
    mean(sample_meta[control_ids, "mean_pseudotime"], na.rm = TRUE)
  list(
    treated_ecdf = treated_ecdf,
    control_ecdf = control_ecdf,
    delta = delta,
    signed_mean_shift = signed_mean_shift,
    n_treated = length(treated_ids),
    n_control = length(control_ids)
  )
}

pseudotime_did_stats_from_ids <- function(ecdf_mat, sample_meta, group_a, group_b, ids) {
  a <- pseudotime_delta_from_ids(ecdf_mat, sample_meta, ids$treated_a, ids$control_a)
  b <- pseudotime_delta_from_ids(ecdf_mat, sample_meta, ids$treated_b, ids$control_b)
  did <- b$delta - a$delta
  c(
    did_ecdf_rmse = sqrt(mean(did^2, na.rm = TRUE)),
    did_ecdf_ks = max(abs(did), na.rm = TRUE),
    did_ecdf_mean_abs = mean(abs(did), na.rm = TRUE),
    did_signed_mean_shift = b$signed_mean_shift - a$signed_mean_shift
  )
}

pseudotime_did_stats_from_labels <- function(ecdf_mat, sample_meta, labels, group_a, group_b) {
  sample_ids <- rownames(sample_meta)
  treatment_group <- ifelse(sample_meta$dose_mg == 0, "control", "treated")
  ids <- list(
    control_a = sample_ids[treatment_group == "control" & labels == group_a],
    treated_a = sample_ids[treatment_group == "treated" & labels == group_a],
    control_b = sample_ids[treatment_group == "control" & labels == group_b],
    treated_b = sample_ids[treatment_group == "treated" & labels == group_b]
  )
  pseudotime_did_stats_from_ids(ecdf_mat, sample_meta, group_a, group_b, ids)
}

ploidy_label_permutations_within_treatment <- function(sample_meta, group_a, group_b, n_perm = 10000L) {
  treatment_group <- ifelse(sample_meta$dose_mg == 0, "control", "treated")
  choices_by_group <- lapply(split(seq_len(nrow(sample_meta)), treatment_group), function(idx) {
    n_a <- sum(as.character(sample_meta$end_timepoint_ploidy_group[idx]) == group_a)
    if (n_a == 0) return(list(integer(0)))
    if (n_a == length(idx)) return(list(idx))
    combn(idx, n_a, simplify = FALSE)
  })
  n_exact <- prod(vapply(choices_by_group, length, integer(1)))
  if (is.finite(n_exact) && n_exact <= n_perm) {
    out <- list()
    recurse <- function(i, selected) {
      if (i > length(choices_by_group)) {
        labels <- rep(group_b, nrow(sample_meta))
        labels[unlist(selected, use.names = FALSE)] <- group_a
        out[[length(out) + 1L]] <<- labels
      } else {
        for (choice in choices_by_group[[i]]) recurse(i + 1L, c(selected, list(choice)))
      }
    }
    recurse(1L, list())
    list(labels = out, mode = "exact_ploidy_label_enumeration_within_treatment_group")
  } else {
    split_idx <- split(seq_len(nrow(sample_meta)), treatment_group)
    observed_n_a <- vapply(split_idx, function(idx) sum(as.character(sample_meta$end_timepoint_ploidy_group[idx]) == group_a), integer(1))
    out <- replicate(n_perm, {
      labels <- rep(group_b, nrow(sample_meta))
      for (nm in names(split_idx)) {
        idx <- split_idx[[nm]]
        if (observed_n_a[[nm]] > 0) labels[sample(idx, observed_n_a[[nm]])] <- group_a
      }
      labels
    }, simplify = FALSE)
    list(labels = out, mode = "monte_carlo_ploidy_label_permutation_within_treatment_group")
  }
}

bootstrap_pseudotime_did <- function(ecdf_mat, sample_meta, group_a, group_b, n_boot = 1000L) {
  sample_ids <- rownames(sample_meta)
  treatment_group <- ifelse(sample_meta$dose_mg == 0, "control", "treated")
  strata <- list(
    control_a = sample_ids[treatment_group == "control" & as.character(sample_meta$end_timepoint_ploidy_group) == group_a],
    treated_a = sample_ids[treatment_group == "treated" & as.character(sample_meta$end_timepoint_ploidy_group) == group_a],
    control_b = sample_ids[treatment_group == "control" & as.character(sample_meta$end_timepoint_ploidy_group) == group_b],
    treated_b = sample_ids[treatment_group == "treated" & as.character(sample_meta$end_timepoint_ploidy_group) == group_b]
  )
  if (any(vapply(strata, length, integer(1)) < 2)) return(data.frame())
  boot <- replicate(n_boot, {
    ids <- lapply(strata, function(x) sample(x, length(x), replace = TRUE))
    pseudotime_did_stats_from_ids(ecdf_mat, sample_meta, group_a, group_b, ids)
  })
  boot <- t(boot)
  q <- function(metric, prob) as.numeric(quantile(boot[, metric], prob, na.rm = TRUE, names = FALSE))
  data.frame(
    did_ecdf_rmse_ci025 = q("did_ecdf_rmse", 0.025),
    did_ecdf_rmse_ci975 = q("did_ecdf_rmse", 0.975),
    did_ecdf_ks_ci025 = q("did_ecdf_ks", 0.025),
    did_ecdf_ks_ci975 = q("did_ecdf_ks", 0.975),
    did_ecdf_mean_abs_ci025 = q("did_ecdf_mean_abs", 0.025),
    did_ecdf_mean_abs_ci975 = q("did_ecdf_mean_abs", 0.975),
    did_signed_mean_shift_ci025 = q("did_signed_mean_shift", 0.025),
    did_signed_mean_shift_ci975 = q("did_signed_mean_shift", 0.975),
    n_boot = n_boot,
    stringsAsFactors = FALSE
  )
}

run_pseudotime_ploidy_effect_modification <- function(results, n_perm = 10000L, n_boot = 1000L) {
  out_dir <- ensure_dir(dirs$Pseudotime)
  ensure_dir(file.path(out_dir, "stats"))
  ensure_dir(file.path(out_dir, "tables"))
  ensure_dir(file.path(out_dir, "figures"))
  all_tests <- list()
  all_curves <- list()

  for (name in names(results)) {
    df <- results[[name]]$df
    sample_meta <- results[[name]]$sample_meta
    ploidies <- etp_group_levels[etp_group_levels %in% as.character(sample_meta$end_timepoint_ploidy_group)]
    if (length(ploidies) < 2) next
    group_a <- ploidies[1]
    group_b <- ploidies[2]
    keep_samples <- rownames(sample_meta)[as.character(sample_meta$end_timepoint_ploidy_group) %in% c(group_a, group_b)]
    local_meta <- sample_meta[keep_samples, , drop = FALSE]
    local_df <- df[df$sample_id %in% keep_samples, , drop = FALSE]
    grid <- make_grid(local_df)
    ecdf_mat <- make_sample_ecdf_matrix(local_df, rownames(local_meta), grid)
    treatment_group <- ifelse(local_meta$dose_mg == 0, "control", "treated")
    original_ids <- list(
      control_a = rownames(local_meta)[treatment_group == "control" & as.character(local_meta$end_timepoint_ploidy_group) == group_a],
      treated_a = rownames(local_meta)[treatment_group == "treated" & as.character(local_meta$end_timepoint_ploidy_group) == group_a],
      control_b = rownames(local_meta)[treatment_group == "control" & as.character(local_meta$end_timepoint_ploidy_group) == group_b],
      treated_b = rownames(local_meta)[treatment_group == "treated" & as.character(local_meta$end_timepoint_ploidy_group) == group_b]
    )
    delta_a <- pseudotime_delta_from_ids(ecdf_mat, local_meta, original_ids$treated_a, original_ids$control_a)
    delta_b <- pseudotime_delta_from_ids(ecdf_mat, local_meta, original_ids$treated_b, original_ids$control_b)
    did_curve <- delta_b$delta - delta_a$delta
    observed <- pseudotime_did_stats_from_ids(ecdf_mat, local_meta, group_a, group_b, original_ids)

    perms <- ploidy_label_permutations_within_treatment(local_meta, group_a, group_b, n_perm)
    perm_stats <- t(vapply(perms$labels, function(labels) {
      pseudotime_did_stats_from_labels(ecdf_mat, local_meta, labels, group_a, group_b)
    }, numeric(4)))
    boot_ci <- bootstrap_pseudotime_did(ecdf_mat, local_meta, group_a, group_b, n_boot)
    if (nrow(boot_ci) == 0) {
      boot_ci <- data.frame(
        did_ecdf_rmse_ci025 = NA_real_, did_ecdf_rmse_ci975 = NA_real_,
        did_ecdf_ks_ci025 = NA_real_, did_ecdf_ks_ci975 = NA_real_,
        did_ecdf_mean_abs_ci025 = NA_real_, did_ecdf_mean_abs_ci975 = NA_real_,
        did_signed_mean_shift_ci025 = NA_real_, did_signed_mean_shift_ci975 = NA_real_,
        n_boot = NA_integer_
      )
    }
    test <- data.frame(
      compartment = name,
      comparison = paste0(group_b, "_minus_", group_a, "_treatment_effect"),
      group_a = group_a,
      group_b = group_b,
      n_control_a = length(original_ids$control_a),
      n_treated_a = length(original_ids$treated_a),
      n_control_b = length(original_ids$control_b),
      n_treated_b = length(original_ids$treated_b),
      did_ecdf_rmse = observed[["did_ecdf_rmse"]],
      did_ecdf_ks = observed[["did_ecdf_ks"]],
      did_ecdf_mean_abs = observed[["did_ecdf_mean_abs"]],
      did_signed_mean_shift = observed[["did_signed_mean_shift"]],
      p_did_ecdf_rmse = mean(perm_stats[, "did_ecdf_rmse"] >= observed[["did_ecdf_rmse"]] - 1e-15, na.rm = TRUE),
      p_did_ecdf_ks = mean(perm_stats[, "did_ecdf_ks"] >= observed[["did_ecdf_ks"]] - 1e-15, na.rm = TRUE),
      p_did_ecdf_mean_abs = mean(perm_stats[, "did_ecdf_mean_abs"] >= observed[["did_ecdf_mean_abs"]] - 1e-15, na.rm = TRUE),
      p_did_signed_mean_shift_two_sided = mean(abs(perm_stats[, "did_signed_mean_shift"]) >= abs(observed[["did_signed_mean_shift"]]) - 1e-15, na.rm = TRUE),
      n_permutations = nrow(perm_stats),
      permutation_mode = perms$mode,
      stringsAsFactors = FALSE
    )
    test <- cbind(test, boot_ci)
    all_tests[[length(all_tests) + 1L]] <- test

    curves <- data.frame(
      compartment = name,
      pseudotime = grid,
      group_a = group_a,
      group_b = group_b,
      delta_ecdf_group_a = delta_a$delta,
      delta_ecdf_group_b = delta_b$delta,
      did_delta_ecdf_group_b_minus_a = did_curve,
      control_ecdf_group_a = delta_a$control_ecdf,
      treated_ecdf_group_a = delta_a$treated_ecdf,
      control_ecdf_group_b = delta_b$control_ecdf,
      treated_ecdf_group_b = delta_b$treated_ecdf,
      stringsAsFactors = FALSE
    )
    all_curves[[length(all_curves) + 1L]] <- curves

    curve_long <- rbind(
      data.frame(compartment = name, pseudotime = grid, curve_type = "treatment_delta_ecdf", group = group_a, value = delta_a$delta, stringsAsFactors = FALSE),
      data.frame(compartment = name, pseudotime = grid, curve_type = "treatment_delta_ecdf", group = group_b, value = delta_b$delta, stringsAsFactors = FALSE),
      data.frame(compartment = name, pseudotime = grid, curve_type = "ploidy_difference_in_treatment_delta", group = paste0(group_b, " - ", group_a), value = did_curve, stringsAsFactors = FALSE)
    )
    p <- ggplot(curve_long, aes(x = pseudotime, y = value, color = group)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_line(linewidth = 0.8) +
      facet_wrap(~ curve_type, ncol = 1, scales = "free_y") +
      labs(
        title = paste0(name, ": pseudotime treatment effect by EndTimePoint ploidy"),
        x = "pseudotime",
        y = "ECDF difference",
        color = "curve"
      ) +
      plot_theme
    save_square_plot(p, file.path(out_dir, "figures"), paste0(name, "_delta_ECDF_ETP_lower_vs_ETP_higher_effect_modification"), 7.2, 9.4)
  }

  tests <- if (length(all_tests) > 0) do.call(rbind, all_tests) else data.frame()
  curves <- if (length(all_curves) > 0) do.call(rbind, all_curves) else data.frame()
  if (nrow(tests) > 0) write_csv(tests, file.path(out_dir, "stats", "pseudotime_ploidy_effect_modification_tests.csv"))
  if (nrow(curves) > 0) write_csv(curves, file.path(out_dir, "tables", "pseudotime_delta_ecdf_curves_by_ploidy.csv"))
  summary_lines <- c(
    "# Pseudotime EndTimePoint-Ploidy Effect-Modification Summary",
    "",
    "This module addresses pseudotime question A2: whether the gemcitabine-induced pseudotime distribution shift differs between ETP-lower and ETP-higher sample groups.",
    "",
    "The estimand is a difference-in-differences ECDF curve:",
    "",
    "- `Delta ECDF_ETP-lower(t) = mean ECDF(treated ETP-lower) - mean ECDF(control ETP-lower)`",
    "- `Delta ECDF_ETP-higher(t) = mean ECDF(treated ETP-higher) - mean ECDF(control ETP-higher)`",
    "- `DID(t) = Delta ECDF_ETP-higher(t) - Delta ECDF_ETP-lower(t)`",
    "",
    "Permutation p values shuffle ETP-lower/ETP-higher labels within treatment groups, preserving the observed number of both ETP groups among control and treated samples. Bootstrap CIs resample samples within each treatment-by-ETP stratum.",
    "",
    "Primary outputs:",
    "",
    "- `stats/pseudotime_ploidy_effect_modification_tests.csv`: DID ECDF RMSE/KS/mean-absolute/signed-mean statistics, permutation p values, and bootstrap CIs.",
    "- `tables/pseudotime_delta_ecdf_curves_by_ploidy.csv`: Delta ECDF curves by ETP group plus the ETP-higher-minus-ETP-lower DID curve.",
    "- `figures/*_delta_ECDF_ETP_lower_vs_ETP_higher_effect_modification.*`: treatment-effect and DID ECDF curves."
  )
  write_text(summary_lines, file.path(out_dir, "pseudotime_ploidy_effect_modification_summary.md"))
  list(tests = tests, curves = curves)
}

bootstrap_mean_diff_ci <- function(value, group, group_a, group_b, n_boot = 1000L, margin = NA_real_) {
  ok <- is.finite(value) & !is.na(group)
  value <- value[ok]
  group <- as.character(group[ok])
  group_a <- as.character(group_a)
  group_b <- as.character(group_b)
  xa <- value[group == group_a]
  xb <- value[group == group_b]
  if (length(xa) < 2 || length(xb) < 2) {
    return(data.frame(
      n_group_a = length(xa), n_group_b = length(xb),
      observed_mean_diff_group_b_minus_a = NA_real_,
      ci025 = NA_real_, ci975 = NA_real_, margin = margin,
      ci_within_margin = NA, n_boot = NA_integer_, stringsAsFactors = FALSE
    ))
  }
  observed <- mean(xb, na.rm = TRUE) - mean(xa, na.rm = TRUE)
  boot <- replicate(n_boot, mean(sample(xb, replace = TRUE), na.rm = TRUE) - mean(sample(xa, replace = TRUE), na.rm = TRUE))
  ci <- as.numeric(quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
  data.frame(
    n_group_a = length(xa),
    n_group_b = length(xb),
    observed_mean_diff_group_b_minus_a = observed,
    ci025 = ci[1],
    ci975 = ci[2],
    margin = margin,
    ci_within_margin = if (is.finite(margin)) ci[1] >= -abs(margin) && ci[2] <= abs(margin) else NA,
    n_boot = n_boot,
    stringsAsFactors = FALSE
  )
}

lm_terms_table <- function(fit, compartment, analysis, model_label, n) {
  if (inherits(fit, "try-error") || is.null(fit)) return(data.frame())
  cf <- coef(summary(fit))
  data.frame(
    compartment = compartment,
    analysis = analysis,
    model = model_label,
    term = rownames(cf),
    estimate = cf[, "Estimate"],
    std_error = cf[, "Std. Error"],
    statistic = cf[, "t value"],
    p_value = cf[, "Pr(>|t|)"],
    n = n,
    stringsAsFactors = FALSE
  )
}

run_ploidy_response_analysis <- function(results, n_perm = 10000L, n_boot = 1000L) {
  ensure_dir(dirs$PloidyResponse)
  all_ecdf <- list()
  all_shift <- list()
  all_models <- list()

  for (name in names(results)) {
    comp_dir <- ensure_dir(file.path(dirs$PloidyResponse, name))
    ensure_dir(file.path(comp_dir, "stats"))
    ensure_dir(file.path(comp_dir, "figures"))
    df <- results[[name]]$df
    sh <- results[[name]]$equal_shift
    df$treatment_group <- ifelse(df$dose_mg == 0, "0mg/kg", "treated")
    ploidies <- etp_group_levels[etp_group_levels %in% as.character(sh$end_timepoint_ploidy_group)]

    ecdf_rows <- list()
    for (ploidy in ploidies) {
      local <- df[as.character(df$end_timepoint_ploidy_group) == ploidy, , drop = FALSE]
      tst <- ecdf_group_test(local, "treatment_group", "0mg/kg", "treated", FALSE, n_perm)
      ecdf_rows[[length(ecdf_rows) + 1L]] <- cbind(compartment = name, end_timepoint_ploidy_group = ploidy, comparison = "0_vs_treated_within_ploidy", tst)
    }
    ecdf_tests <- if (length(ecdf_rows) > 0) do.call(rbind, ecdf_rows) else data.frame()
    if (nrow(ecdf_tests) > 0) write_csv(ecdf_tests, file.path(comp_dir, "stats", paste0(name, "_ploidy_stratified_treatment_ecdf_tests.csv")))
    all_ecdf[[length(all_ecdf) + 1L]] <- ecdf_tests

    treated <- sh[sh$dose_mg > 0, , drop = FALSE]
    shift_rows <- list()
    if (length(ploidies) >= 2 && nrow(treated) > 0) {
      group_a <- ploidies[1]
      group_b <- ploidies[2]
      for (metric in c("ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift")) {
        a_sum <- dose_value_summary(treated[[metric]], treated$end_timepoint_ploidy_group, group_a)
        b_sum <- dose_value_summary(treated[[metric]], treated$end_timepoint_ploidy_group, group_b)
        wt <- safe_wilcox(treated[[metric]], treated$end_timepoint_ploidy_group)
        tt <- safe_ttest(treated[[metric]], treated$end_timepoint_ploidy_group)
        perm <- permutation_mean_diff_p(treated[[metric]], treated$end_timepoint_ploidy_group, group_a, group_b, n_perm, exact = TRUE)
        margin <- if (metric == "ecdf_rmse") 0.01 else NA_real_
        eq <- bootstrap_mean_diff_ci(treated[[metric]], treated$end_timepoint_ploidy_group, group_a, group_b, n_boot, margin)
        shift_rows[[length(shift_rows) + 1L]] <- data.frame(
          compartment = name,
          sample_set = "treated",
          comparison = paste0(group_a, "_vs_", group_b),
          group_a = group_a,
          group_b = group_b,
          response_metric = metric,
          n_a = a_sum$n,
          n_b = b_sum$n,
          mean_a = a_sum$mean,
          mean_b = b_sum$mean,
          median_a = a_sum$median,
          median_b = b_sum$median,
          mean_diff_b_minus_a = perm$observed_mean_diff_group_b_minus_a,
          wilcoxon_p = wt$p_value,
          t_test_p = tt$p_value,
          permutation_p_two_sided = perm$permutation_p_two_sided,
          bootstrap_ci025 = eq$ci025,
          bootstrap_ci975 = eq$ci975,
          equivalence_margin = eq$margin,
          ci_within_equivalence_margin = eq$ci_within_margin,
          stringsAsFactors = FALSE
        )
      }
    }
    shift_tests <- if (length(shift_rows) > 0) do.call(rbind, shift_rows) else data.frame()
    if (nrow(shift_tests) > 0) {
      shift_tests$permutation_q_exploratory <- p.adjust(shift_tests$permutation_p_two_sided, method = "BH")
      write_csv(shift_tests, file.path(comp_dir, "stats", paste0(name, "_treated_shift_by_ploidy_tests.csv")))
      all_shift[[length(all_shift) + 1L]] <- shift_tests
    }

    plot_data <- treated
    if (nrow(plot_data) > 0) {
      p_shift <- ggplot(plot_data, aes(x = end_timepoint_ploidy_group, y = ecdf_rmse, color = end_timepoint_ploidy_group)) +
        geom_boxplot(aes(group = end_timepoint_ploidy_group), outlier.shape = NA, color = "grey55", linewidth = 0.45) +
        geom_point(aes(shape = dose), position = position_jitter(width = 0.08, height = 0), size = 2.8) +
        labs(title = paste0(name, ": treated ECDF RMSE by EndTimePoint ploidy"), x = "EndTimePoint ploidy", y = "ECDF RMSE") +
        plot_theme +
        theme(legend.position = "right")
      save_square_plot(p_shift, file.path(comp_dir, "figures"), paste0(name, "_treated_ecdf_rmse_by_end_timepoint_ploidy_group"), 6.6, 6.6)
    }

    model_rows <- list()
    local_all <- sh[is.finite(sh$ecdf_rmse), , drop = FALSE]
    if (nrow(local_all) >= 6 && length(unique(local_all$end_timepoint_ploidy_group)) > 1 && length(unique(local_all$dose_mg > 0)) > 1) {
      local_all$treatment_group <- ifelse(local_all$dose_mg == 0, "control", "treated")
      fit <- try(lm(ecdf_rmse ~ treatment_group * end_timepoint_ploidy_group, data = local_all), silent = TRUE)
      model_rows[[length(model_rows) + 1L]] <- lm_terms_table(fit, name, "shift_treatment_by_ploidy_interaction", "ecdf_rmse ~ treatment_group * end_timepoint_ploidy_group", nrow(local_all))
    }
    local_treated <- treated[is.finite(treated$ecdf_rmse), , drop = FALSE]
    if (nrow(local_treated) >= 6 && length(unique(local_treated$end_timepoint_ploidy_group)) > 1) {
      fit <- try(lm(ecdf_rmse ~ factor(dose_mg) + end_timepoint_ploidy_group, data = local_treated), silent = TRUE)
      model_rows[[length(model_rows) + 1L]] <- lm_terms_table(fit, name, "treated_shift_by_dose_and_ploidy", "ecdf_rmse ~ factor(dose_mg) + end_timepoint_ploidy_group", nrow(local_treated))
    }
    local_tgi <- treated[is.finite(treated$TGI_percent_auc), , drop = FALSE]
    if (nrow(local_tgi) >= 6 && length(unique(local_tgi$end_timepoint_ploidy_group)) > 1) {
      fit <- try(lm(TGI_percent_auc ~ factor(dose_mg) + end_timepoint_ploidy_group, data = local_tgi), silent = TRUE)
      model_rows[[length(model_rows) + 1L]] <- lm_terms_table(fit, name, "treated_TGI_by_dose_and_ploidy", "TGI_percent_auc ~ factor(dose_mg) + end_timepoint_ploidy_group", nrow(local_tgi))
      fit <- try(lm(TGI_percent_auc ~ scale(ecdf_rmse) * end_timepoint_ploidy_group + factor(dose_mg), data = local_tgi), silent = TRUE)
      model_rows[[length(model_rows) + 1L]] <- lm_terms_table(fit, name, "treated_TGI_shift_by_ploidy_interaction", "TGI_percent_auc ~ scale(ecdf_rmse) * end_timepoint_ploidy_group + factor(dose_mg)", nrow(local_tgi))
    }
    models <- if (length(model_rows) > 0) do.call(rbind, model_rows) else data.frame()
    if (nrow(models) > 0) {
      write_csv(models, file.path(comp_dir, "stats", paste0(name, "_ploidy_response_interaction_models.csv")))
      all_models[[length(all_models) + 1L]] <- models
    }
  }

  tgi_dir <- ensure_dir(file.path(dirs$PloidyResponse, "TGI"))
  ensure_dir(file.path(tgi_dir, "stats"))
  ensure_dir(file.path(tgi_dir, "figures"))
  treated_samples <- results$CellCycle$equal_shift[results$CellCycle$equal_shift$dose_mg > 0, , drop = FALSE]
  ploidies <- etp_group_levels[etp_group_levels %in% as.character(treated_samples$end_timepoint_ploidy_group)]
  tgi_cols <- unique(c("TGI_percent_auc", endpoint_tgi_cols(names(treated_samples))))
  tgi_rows <- list()
  if (length(ploidies) >= 2) {
    group_a <- ploidies[1]
    group_b <- ploidies[2]
    for (tgi_col in tgi_cols) {
      a_sum <- dose_value_summary(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group, group_a)
      b_sum <- dose_value_summary(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group, group_b)
      wt <- safe_wilcox(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group)
      tt <- safe_ttest(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group)
      perm <- permutation_mean_diff_p(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group, group_a, group_b, n_perm, exact = TRUE)
      margin <- if (tgi_col == "TGI_percent_auc") 10 else NA_real_
      eq <- bootstrap_mean_diff_ci(treated_samples[[tgi_col]], treated_samples$end_timepoint_ploidy_group, group_a, group_b, n_boot, margin)
      tgi_rows[[length(tgi_rows) + 1L]] <- data.frame(
        comparison = paste0(group_a, "_vs_", group_b),
        group_a = group_a,
        group_b = group_b,
        tgi_measure = tgi_col,
        tgi_day = ifelse(grepl("^TGI_percent_Day_", tgi_col), endpoint_day_num(tgi_col), NA_real_),
        n_a = a_sum$n,
        n_b = b_sum$n,
        mean_a = a_sum$mean,
        mean_b = b_sum$mean,
        median_a = a_sum$median,
        median_b = b_sum$median,
        mean_diff_b_minus_a = perm$observed_mean_diff_group_b_minus_a,
        wilcoxon_p = wt$p_value,
        t_test_p = tt$p_value,
        permutation_p_two_sided = perm$permutation_p_two_sided,
        bootstrap_ci025 = eq$ci025,
        bootstrap_ci975 = eq$ci975,
        equivalence_margin = eq$margin,
        ci_within_equivalence_margin = eq$ci_within_margin,
        stringsAsFactors = FALSE
      )
    }
  }
  tgi_tests <- if (length(tgi_rows) > 0) do.call(rbind, tgi_rows) else data.frame()
  if (nrow(tgi_tests) > 0) {
    tgi_tests$permutation_q_exploratory <- p.adjust(tgi_tests$permutation_p_two_sided, method = "BH")
    write_csv(tgi_tests, file.path(tgi_dir, "stats", "treated_TGI_auc_endpoint_by_ploidy_tests.csv"))
  }
  endpoint_long <- make_endpoint_tgi_long(treated_samples, endpoint_tgi_cols(names(treated_samples)))
  if (nrow(endpoint_long) > 0) {
    write_csv(endpoint_long, file.path(tgi_dir, "stats", "treated_endpoint_TGI_by_ploidy_long.csv"))
    p_endpoint <- ggplot(endpoint_long, aes(x = day, y = TGI_percent_endpoint, color = end_timepoint_ploidy_group, group = sample_id)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_line(alpha = 0.6, linewidth = 0.55) +
      geom_point(aes(shape = dose), size = 2.1) +
      labs(title = "Endpoint TGI over days by EndTimePoint ploidy", x = "day", y = "endpoint TGI (%)") +
      plot_theme
    save_plot(p_endpoint, file.path(tgi_dir, "figures"), "treated_endpoint_TGI_by_end_timepoint_ploidy_group_over_days", 7.4, 4.8)
  }
  if (nrow(treated_samples) > 0) {
    p_auc <- ggplot(treated_samples, aes(x = end_timepoint_ploidy_group, y = TGI_percent_auc, color = end_timepoint_ploidy_group)) +
      geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
      geom_boxplot(aes(group = end_timepoint_ploidy_group), outlier.shape = NA, color = "grey55", linewidth = 0.45) +
      geom_point(aes(shape = dose), position = position_jitter(width = 0.08, height = 0), size = 2.8) +
      geom_text(aes(label = sample_id), size = 2.35, nudge_y = 2.5, check_overlap = TRUE, show.legend = FALSE) +
      labs(title = "AUC TGI by EndTimePoint ploidy", x = "EndTimePoint ploidy", y = "AUC-based TGI (%)") +
      plot_theme
    save_square_plot(p_auc, file.path(tgi_dir, "figures"), "treated_TGI_AUC_by_end_timepoint_ploidy_group", 6.6, 6.6)
  }

  ecdf <- if (length(all_ecdf) > 0) do.call(rbind, all_ecdf) else data.frame()
  shift <- if (length(all_shift) > 0) do.call(rbind, all_shift) else data.frame()
  models <- if (length(all_models) > 0) do.call(rbind, all_models) else data.frame()
  if (nrow(ecdf) > 0) write_csv(ecdf, file.path(dirs$PloidyResponse, "all_compartments_ploidy_stratified_treatment_ecdf_tests.csv"))
  if (nrow(shift) > 0) write_csv(shift, file.path(dirs$PloidyResponse, "all_compartments_treated_shift_by_ploidy_tests.csv"))
  if (nrow(models) > 0) write_csv(models, file.path(dirs$PloidyResponse, "all_compartments_ploidy_response_interaction_models.csv"))

  summary_lines <- c(
    "# EndTimePoint-Ploidy Response Summary",
    "",
    "This module treats EndTimePoint ploidy as a potential response modifier rather than only as a confounder.",
    "",
    "Difference tests ask whether ETP-lower and ETP-higher responses are detectably different. Equivalence-style rows use pre-specified descriptive margins: ECDF RMSE +/-0.01 and AUC TGI +/-10 percentage points. If the bootstrap CI is not contained inside the margin, equivalence is not established.",
    "",
    "Primary outputs:",
    "",
    "- `all_compartments_ploidy_stratified_treatment_ecdf_tests.csv`: 0 vs treated ECDF tests within each ploidy and compartment.",
    "- `all_compartments_treated_shift_by_ploidy_tests.csv`: treated-sample shift differences between ETP-lower and ETP-higher.",
    "- `TGI/stats/treated_TGI_auc_endpoint_by_ploidy_tests.csv`: treated-sample AUC and endpoint TGI differences between ETP-lower and ETP-higher.",
    "- `all_compartments_ploidy_response_interaction_models.csv`: ploidy interaction and adjusted sensitivity models."
  )
  write_text(summary_lines, file.path(dirs$PloidyResponse, "ploidy_response_summary.md"))

  list(ecdf = ecdf, shift = shift, tgi = tgi_tests, models = models)
}

compartment_comparison <- function(results) {
  cc <- results$CellCycle$equal_shift
  nc <- results$NonCellCycle$equal_shift
  keep <- c("sample_id", "end_timepoint_ploidy_group", "dose", "dose_mg", "TGI_percent_auc", "ecdf_rmse", "ecdf_ks", "ecdf_mean_abs", "signed_mean_shift", "cell_count")
  comp <- merge(cc[, keep], nc[, keep], by = "sample_id", suffixes = c("_CellCycle", "_NonCellCycle"))
  comp$dose <- comp$dose_CellCycle
  comp$dose_mg <- comp$dose_mg_CellCycle
  comp$end_timepoint_ploidy_group <- comp$end_timepoint_ploidy_group_CellCycle
  comp$TGI_percent_auc <- comp$TGI_percent_auc_CellCycle
  comp$ecdf_rmse_delta <- comp$ecdf_rmse_CellCycle - comp$ecdf_rmse_NonCellCycle
  comp$signed_mean_shift_delta <- comp$signed_mean_shift_CellCycle - comp$signed_mean_shift_NonCellCycle
  ensure_dir(file.path(dirs$Comparison, "stats"))
  ensure_dir(file.path(dirs$Comparison, "figures"))
  write_csv(comp, file.path(dirs$Comparison, "stats", "cellcycle_vs_noncellcycle_sample_shift_table.csv"))
  rows <- list()
  for (sample_set in c("all", "treated")) {
    local <- if (sample_set == "treated") comp[comp$dose_mg > 0, , drop = FALSE] else comp
    for (metric in c("ecdf_rmse", "signed_mean_shift")) {
      x <- local[[paste0(metric, "_CellCycle")]]
      y <- local[[paste0(metric, "_NonCellCycle")]]
      pt <- safe_paired(x, y, "paired_t")
      pw <- safe_paired(x, y, "wilcoxon")
      rows[[length(rows) + 1L]] <- data.frame(sample_set = sample_set, metric = metric, test = "paired_t", n = pt$n, statistic = pt$statistic, p_value = pt$p_value, mean_delta_CellCycle_minus_NonCellCycle = mean(x - y, na.rm = TRUE), stringsAsFactors = FALSE)
      rows[[length(rows) + 1L]] <- data.frame(sample_set = sample_set, metric = metric, test = "wilcoxon_signed_rank", n = pw$n, statistic = pw$statistic, p_value = pw$p_value, mean_delta_CellCycle_minus_NonCellCycle = mean(x - y, na.rm = TRUE), stringsAsFactors = FALSE)
    }
  }
  paired <- do.call(rbind, rows)
  write_csv(paired, file.path(dirs$Comparison, "stats", "paired_compartment_shift_tests_v4.csv"))
  treated <- comp[comp$dose_mg > 0, , drop = FALSE]
  obs_cc <- safe_cor_result(treated$ecdf_rmse_CellCycle, treated$TGI_percent_auc, "pearson")$estimate
  obs_nc <- safe_cor_result(treated$ecdf_rmse_NonCellCycle, treated$TGI_percent_auc, "pearson")$estimate
  obs_delta <- obs_cc - obs_nc
  n <- nrow(treated)
  swap_mat <- expand.grid(rep(list(c(FALSE, TRUE)), n))
  delta_perm <- apply(swap_mat, 1, function(sw) {
    x_cc <- ifelse(sw, treated$ecdf_rmse_NonCellCycle, treated$ecdf_rmse_CellCycle)
    x_nc <- ifelse(sw, treated$ecdf_rmse_CellCycle, treated$ecdf_rmse_NonCellCycle)
    safe_cor_result(x_cc, treated$TGI_percent_auc, "pearson")$estimate - safe_cor_result(x_nc, treated$TGI_percent_auc, "pearson")$estimate
  })
  perm_table <- data.frame(
    observed_r_CellCycle = obs_cc,
    observed_r_NonCellCycle = obs_nc,
    observed_delta_r = obs_delta,
    n = n,
    n_permutations = length(delta_perm),
    p_two_sided = mean(abs(delta_perm) >= abs(obs_delta) - 1e-15, na.rm = TRUE),
    p_positive_delta = mean(delta_perm >= obs_delta - 1e-15, na.rm = TRUE),
    permutation_mode = "within_sample_compartment_label_swap",
    stringsAsFactors = FALSE
  )
  write_csv(perm_table, file.path(dirs$Comparison, "stats", "compartment_TGI_correlation_difference_permutation.csv"))
  assoc_compare <- rbind(
    results$CellCycle$tgi_stats[results$CellCycle$tgi_stats$reference_type == "primary_equal_sample_reference" & results$CellCycle$tgi_stats$shift_metric == "ecdf_rmse" & results$CellCycle$tgi_stats$tgi_measure == "TGI_percent_auc", ],
    results$NonCellCycle$tgi_stats[results$NonCellCycle$tgi_stats$reference_type == "primary_equal_sample_reference" & results$NonCellCycle$tgi_stats$shift_metric == "ecdf_rmse" & results$NonCellCycle$tgi_stats$tgi_measure == "TGI_percent_auc", ]
  )
  write_csv(assoc_compare, file.path(dirs$Comparison, "stats", "CellCycle_vs_NonCellCycle_TGI_association_comparison.csv"))
  long <- rbind(
    data.frame(sample_id = comp$sample_id, dose = comp$dose, end_timepoint_ploidy_group = comp$end_timepoint_ploidy_group, compartment = "CellCycle", ecdf_rmse = comp$ecdf_rmse_CellCycle, signed_mean_shift = comp$signed_mean_shift_CellCycle, TGI_percent_auc = comp$TGI_percent_auc, stringsAsFactors = FALSE),
    data.frame(sample_id = comp$sample_id, dose = comp$dose, end_timepoint_ploidy_group = comp$end_timepoint_ploidy_group, compartment = "NonCellCycle", ecdf_rmse = comp$ecdf_rmse_NonCellCycle, signed_mean_shift = comp$signed_mean_shift_NonCellCycle, TGI_percent_auc = comp$TGI_percent_auc, stringsAsFactors = FALSE)
  )
  p1 <- ggplot(long, aes(x = compartment, y = ecdf_rmse, group = sample_id)) +
    geom_line(color = "grey65", linewidth = 0.4) +
    geom_point(aes(color = dose, shape = end_timepoint_ploidy_group), size = 2.5) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% unique(long$dose)], name = "dose") +
    labs(title = "Paired CellCycle vs NonCellCycle ECDF RMSE", x = NULL, y = "ECDF RMSE") +
    plot_theme
  save_square_plot(p1, file.path(dirs$Comparison, "figures"), "CellCycle_vs_NonCellCycle_paired_ecdf_rmse_v4", 6.4, 6.4)
  p2 <- ggplot(long, aes(x = compartment, y = signed_mean_shift, group = sample_id)) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_line(color = "grey65", linewidth = 0.4) +
    geom_point(aes(color = dose, shape = end_timepoint_ploidy_group), size = 2.5) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% unique(long$dose)], name = "dose") +
    labs(title = "Paired CellCycle vs NonCellCycle signed shift", x = NULL, y = "signed mean shift") +
    plot_theme
  save_plot(p2, file.path(dirs$Comparison, "figures"), "CellCycle_vs_NonCellCycle_paired_signed_shift_v4", 6.3, 4.8)
  p3 <- ggplot(data.frame(delta_r = delta_perm), aes(x = delta_r)) +
    geom_histogram(bins = 35, fill = "#999999", color = "white") +
    geom_vline(xintercept = obs_delta, color = "#e41a1c", linewidth = 0.8) +
    labs(title = "Permutation null for CellCycle minus NonCellCycle TGI correlation", x = "delta r", y = "permutations") +
    plot_theme
  save_square_plot(p3, file.path(dirs$Comparison, "figures"), "compartment_TGI_delta_r_permutation", 6.4, 6.4)
  p4 <- ggplot(long[long$dose != "0mg/kg", ], aes(x = ecdf_rmse, y = TGI_percent_auc, color = compartment, shape = dose)) +
    geom_smooth(aes(group = compartment), method = "lm", se = FALSE, linewidth = 0.55) +
    geom_point(size = 2.6) +
    labs(title = "Treated samples: TGI AUC vs ECDF RMSE by compartment", x = "ECDF RMSE", y = "AUC TGI (%)") +
    plot_theme
  save_square_plot(p4, file.path(dirs$Comparison, "figures"), "CellCycle_vs_NonCellCycle_TGI_association", 6.8, 6.8)
  list(table = comp, paired = paired, delta = perm_table, assoc_compare = assoc_compare)
}

permutation_prop_test <- function(prop, group, n_perm = 10000L) {
  ok <- is.finite(prop) & !is.na(group)
  prop <- prop[ok]
  group <- as.character(group[ok])
  if (length(prop) < 3 || length(unique(group)) != 2) return(c(statistic = NA_real_, p_value = NA_real_, n_permutations = NA_real_))
  lev <- unique(group)
  obs <- abs(mean(prop[group == lev[2]], na.rm = TRUE) - mean(prop[group == lev[1]], na.rm = TRUE))
  n1 <- sum(group == lev[1])
  total <- choose(length(group), n1)
  if (total <= n_perm) {
    choices <- combn(seq_along(group), n1, simplify = FALSE)
    perm <- vapply(choices, function(idx) {
      gp <- rep(lev[2], length(group))
      gp[idx] <- lev[1]
      abs(mean(prop[gp == lev[2]], na.rm = TRUE) - mean(prop[gp == lev[1]], na.rm = TRUE))
    }, numeric(1))
  } else {
    perm <- replicate(n_perm, {
      gp <- sample(group)
      abs(mean(prop[gp == lev[2]], na.rm = TRUE) - mean(prop[gp == lev[1]], na.rm = TRUE))
    })
  }
  c(statistic = obs, p_value = mean(perm >= obs - 1e-15, na.rm = TRUE), n_permutations = length(perm))
}

composition_analysis <- function(results, n_perm = 10000L) {
  cc <- results$CellCycle$df
  nc <- results$NonCellCycle$df
  cc_counts <- as.data.frame(table(sample_id = cc$sample_id), stringsAsFactors = FALSE)
  names(cc_counts)[2] <- "CellCycle_cells"
  nc_counts <- as.data.frame(table(sample_id = nc$sample_id), stringsAsFactors = FALSE)
  names(nc_counts)[2] <- "NonCellCycle_cells"
  comp <- merge(cc_counts, nc_counts, by = "sample_id", all = TRUE)
  comp[is.na(comp)] <- 0
  meta <- results$CellCycle$sample_meta[, c("sample_id", "end_timepoint_ploidy_group", "dose", "dose_mg", "TGI_percent_auc"), drop = FALSE]
  comp <- merge(meta, comp, by = "sample_id", all.x = TRUE)
  comp$total_tumor_cells <- comp$CellCycle_cells + comp$NonCellCycle_cells
  comp$CellCycle_fraction <- comp$CellCycle_cells / comp$total_tumor_cells
  comp$NonCellCycle_fraction <- comp$NonCellCycle_cells / comp$total_tumor_cells
  cluster_counts <- aggregate(cell_id ~ sample_id + cluster, rbind(cc, nc), length)
  names(cluster_counts)[3] <- "cells"
  cluster_totals <- aggregate(cells ~ sample_id, cluster_counts, sum)
  names(cluster_totals)[2] <- "total_cells_cluster_table"
  cluster_counts <- merge(cluster_counts, cluster_totals, by = "sample_id")
  cluster_counts$fraction_all_tumor <- cluster_counts$cells / cluster_counts$total_cells_cluster_table
  write_csv(comp, file.path(dirs$Composition, "sample_celltype_cluster_composition.csv"))
  group <- ifelse(comp$dose_mg == 0, "control", "treated")
  tests <- list()
  for (value in c("CellCycle_fraction", "NonCellCycle_fraction")) {
    wt <- safe_wilcox(comp[[value]], group)
    perm <- permutation_prop_test(comp[[value]], group, n_perm)
    glm_p <- NA_real_
    if (all(is.finite(comp$CellCycle_cells)) && length(unique(group)) == 2) {
      tmp <- comp
      tmp$treatment_group <- factor(group)
      fit <- try(suppressWarnings(glm(cbind(CellCycle_cells, NonCellCycle_cells) ~ treatment_group + end_timepoint_ploidy_group, family = quasibinomial(), data = tmp)), silent = TRUE)
      if (!inherits(fit, "try-error") && "treatment_grouptreated" %in% rownames(summary(fit)$coefficients)) {
        glm_p <- summary(fit)$coefficients["treatment_grouptreated", "Pr(>|t|)"]
      }
    }
    tests[[length(tests) + 1L]] <- data.frame(fraction = value, comparison = "control_vs_treated", wilcox_p = wt$p_value, permutation_p = perm[["p_value"]], n_permutations = perm[["n_permutations"]], quasibinomial_p = glm_p, stringsAsFactors = FALSE)
  }
  fraction_tests <- do.call(rbind, tests)
  write_csv(fraction_tests, file.path(dirs$Composition, "cellcycle_fraction_treatment_tests.csv"))
  cl_tests <- lapply(split(cluster_counts, cluster_counts$cluster), function(tab) {
    tab <- merge(meta[, c("sample_id", "dose_mg", "end_timepoint_ploidy_group")], tab, by = "sample_id", all.x = TRUE)
    tab$cells[is.na(tab$cells)] <- 0
    tab$fraction_all_tumor[is.na(tab$fraction_all_tumor)] <- 0
    gp <- ifelse(tab$dose_mg == 0, "control", "treated")
    wt <- safe_wilcox(tab$fraction_all_tumor, gp)
    data.frame(cluster = as.character(tab$cluster[which(!is.na(tab$cluster))[1]]), comparison = "control_vs_treated", n = wt$n, wilcox_p = wt$p_value, treated_minus_control_mean_fraction = mean(tab$fraction_all_tumor[gp == "treated"], na.rm = TRUE) - mean(tab$fraction_all_tumor[gp == "control"], na.rm = TRUE), stringsAsFactors = FALSE)
  })
  cl_tests <- do.call(rbind, cl_tests)
  write_csv(cl_tests, file.path(dirs$Composition, "cluster_fraction_treatment_tests.csv"))
  cluster_shift_rows <- list()
  cluster_tgi_rows <- list()
  for (cl in sort(unique(cc$cluster))) {
    sub <- cc[cc$cluster == cl, , drop = FALSE]
    if (length(unique(sub$sample_id)) < 6) next
    sm <- build_sample_meta(sub)
    sh <- calculate_shift_metrics(sub, sm, "primary_equal_sample_reference")
    dt <- run_dose_tests(sub, min(n_perm, 2000L))
    treated_sh <- sh[sh$dose_mg > 0, , drop = FALSE]
    pearson <- permutation_cor_p(treated_sh$ecdf_rmse, treated_sh$TGI_percent_auc, "pearson", min(n_perm, 2000L), exact = FALSE)
    spearman <- permutation_cor_p(treated_sh$ecdf_rmse, treated_sh$TGI_percent_auc, "spearman", min(n_perm, 2000L), exact = FALSE)
    primary <- data.frame(
      compartment = paste0("CellCycle_cluster_", cl),
      sample_set = "treated",
      reference_type = "primary_equal_sample_reference",
      shift_metric = "ecdf_rmse",
      tgi_measure = "TGI_percent_auc",
      tgi_day = NA_real_,
      n = pearson$n,
      pearson_r = pearson$estimate,
      pearson_p_asymptotic = pearson$asymptotic_p,
      pearson_p_permutation_two_sided = pearson$permutation_p_two_sided,
      pearson_p_permutation_positive = pearson$permutation_p_positive,
      spearman_rho = spearman$estimate,
      spearman_p_asymptotic = spearman$asymptotic_p,
      spearman_p_permutation_two_sided = spearman$permutation_p_two_sided,
      spearman_p_permutation_positive = spearman$permutation_p_positive,
      stringsAsFactors = FALSE
    )
    cluster_shift_rows[[length(cluster_shift_rows) + 1L]] <- data.frame(cluster = cl, dt, stringsAsFactors = FALSE)
    cluster_tgi_rows[[length(cluster_tgi_rows) + 1L]] <- data.frame(cluster = cl, primary, stringsAsFactors = FALSE)
  }
  cluster_shift <- if (length(cluster_shift_rows) > 0) do.call(rbind, cluster_shift_rows) else data.frame()
  cluster_tgi <- if (length(cluster_tgi_rows) > 0) do.call(rbind, cluster_tgi_rows) else data.frame()
  write_csv(cluster_shift, file.path(dirs$Composition, "cellcycle_cluster_specific_pseudotime_tests.csv"))
  write_csv(cluster_tgi, file.path(dirs$Composition, "cellcycle_cluster_specific_TGI_associations.csv"))
  p1 <- ggplot(comp, aes(x = dose, y = CellCycle_fraction, color = dose)) +
    geom_boxplot(aes(group = dose), outlier.shape = NA, color = "grey55", linewidth = 0.45) +
    geom_point(aes(shape = end_timepoint_ploidy_group), position = position_jitter(width = 0.08, height = 0), size = 2.6) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% unique(comp$dose)], name = "dose") +
    labs(title = "CellCycle fraction among tumor cells by dose", x = "dose", y = "CellCycle fraction") +
    plot_theme
  save_plot(p1, dirs$Composition, "CellCycle_fraction_by_dose", 6.4, 4.6)
  p2 <- ggplot(cluster_counts, aes(x = cluster, y = fraction_all_tumor, fill = cluster)) +
    geom_boxplot(outlier.shape = NA, color = "grey45", linewidth = 0.35) +
    facet_wrap(~ sample_id, ncol = 4) +
    labs(title = "Cluster fractions by sample", x = "cluster", y = "fraction of tumor cells") +
    plot_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  save_plot(p2, dirs$Composition, "CellCycle_cluster_fractions_by_dose", 9, 7)
  if (nrow(cluster_shift) > 0) {
    csp <- cluster_shift[cluster_shift$comparison == "0_vs_30plus120" & !cluster_shift$stratified_by_end_timepoint_ploidy_group, , drop = FALSE]
    p3 <- ggplot(csp, aes(x = as.character(cluster), y = observed_ecdf_rmse)) +
      geom_col(fill = "#377eb8", alpha = 0.85) +
      labs(title = "CellCycle cluster-specific treated-vs-control pseudotime shift", x = "cluster", y = "ECDF RMSE") +
      plot_theme
    save_square_plot(p3, dirs$Composition, "CellCycle_cluster_specific_pseudotime_shift", 6.4, 6.4)
  }
  if (nrow(cluster_tgi) > 0) {
    p4 <- ggplot(cluster_tgi, aes(x = as.character(cluster), y = pearson_r)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
      geom_col(fill = "#1b9e77", alpha = 0.85) +
      labs(title = "CellCycle cluster-specific TGI association", x = "cluster", y = "Pearson r") +
      plot_theme
    save_square_plot(p4, dirs$Composition, "CellCycle_cluster_specific_TGI_associations", 6.4, 6.4)
  }
  list(sample_composition = comp, fraction_tests = fraction_tests, cluster_tests = cl_tests, cluster_shift = cluster_shift, cluster_tgi = cluster_tgi)
}

make_growth_curve_long <- function(sample_meta) {
  raw_cols <- day_cols_by_prefix(names(sample_meta), "tumor_volume_")
  delta_cols <- day_cols_by_prefix(names(sample_meta), "tumor_volume_delta_")
  tgi_cols <- endpoint_tgi_cols(names(sample_meta))
  baseline_days <- if ("tumor_volume_baseline_day" %in% names(sample_meta)) day_label_num(sample_meta$tumor_volume_baseline_day) else numeric(0)
  days <- sort(unique(c(
    baseline_days[is.finite(baseline_days)],
    safe_num(sub("^tumor_volume_Day_", "", raw_cols)),
    safe_num(sub("^tumor_volume_delta_Day_", "", delta_cols)),
    safe_num(sub("^TGI_percent_Day_", "", tgi_cols))
  )))
  rows <- lapply(seq_len(nrow(sample_meta)), function(i) {
    sr <- sample_meta[i, , drop = FALSE]
    baseline_day <- if ("tumor_volume_baseline_day" %in% names(sr)) day_label_num(sr$tumor_volume_baseline_day) else 0
    baseline_volume <- if ("tumor_volume_baseline" %in% names(sr)) safe_num(sr$tumor_volume_baseline) else NA_real_
    do.call(rbind, lapply(days, function(day) {
      dl <- paste0("Day_", day)
      is_baseline <- is.finite(baseline_day) && day == baseline_day
      tumor_volume <- if (is_baseline) {
        baseline_volume
      } else if (paste0("tumor_volume_", dl) %in% names(sr)) {
        sr[[paste0("tumor_volume_", dl)]]
      } else {
        NA_real_
      }
      tumor_volume_delta <- if (is_baseline) {
        0
      } else if (paste0("tumor_volume_delta_", dl) %in% names(sr)) {
        sr[[paste0("tumor_volume_delta_", dl)]]
      } else {
        NA_real_
      }
      data.frame(
        sample_id = sr$sample_id,
        end_timepoint_ploidy_group = sr$end_timepoint_ploidy_group,
        dose = sr$dose,
        dose_mg = sr$dose_mg,
        day = day,
        day_label = dl,
        tumor_volume_baseline = baseline_volume,
        tumor_volume = tumor_volume,
        tumor_volume_delta = tumor_volume_delta,
        relative_volume = if (is.finite(tumor_volume) && is.finite(baseline_volume) && baseline_volume > 0) tumor_volume / baseline_volume else NA_real_,
        log_relative_volume = if (is.finite(tumor_volume) && is.finite(baseline_volume) && baseline_volume > 0 && tumor_volume > 0) log(tumor_volume / baseline_volume) else NA_real_,
        TGI_percent_endpoint = if (is_baseline) 0 else if (paste0("TGI_percent_", dl) %in% names(sr)) sr[[paste0("TGI_percent_", dl)]] else NA_real_,
        TGI_percent_auc = sr$TGI_percent_auc,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

extract_model_coefs <- function(fit, model_name) {
  if (inherits(fit, "try-error") || is.null(fit)) return(data.frame())
  if (inherits(fit, "lme")) {
    tt <- summary(fit)$tTable
    data.frame(model = model_name, term = rownames(tt), estimate = tt[, "Value"], std_error = tt[, "Std.Error"], statistic = tt[, "t-value"], p_value = tt[, "p-value"], stringsAsFactors = FALSE)
  } else {
    tt <- coef(summary(fit))
    data.frame(model = model_name, term = rownames(tt), estimate = tt[, "Estimate"], std_error = tt[, "Std. Error"], statistic = tt[, "t value"], p_value = tt[, "Pr(>|t|)"], stringsAsFactors = FALSE)
  }
}

growth_model_analysis <- function(results) {
  meta <- results$CellCycle$sample_meta
  growth <- make_growth_curve_long(meta)
  growth <- merge(growth, results$CellCycle$equal_shift[, c("sample_id", "ecdf_rmse", "signed_mean_shift")], by = "sample_id", all.x = TRUE)
  growth$log_volume <- log1p(growth$tumor_volume)
  growth$dose <- factor(growth$dose, levels = unique(meta$dose[order(meta$dose_mg)]))
  growth$end_timepoint_ploidy_group <- factor(growth$end_timepoint_ploidy_group, levels = etp_group_levels)
  fit1 <- NULL
  fit2 <- NULL
  mode1 <- "lm"
  mode2 <- "lm"
  if (requireNamespace("nlme", quietly = TRUE) && requireNamespace("splines", quietly = TRUE)) {
    fit1 <- try(nlme::lme(log_volume ~ splines::ns(day, df = 3) * dose + end_timepoint_ploidy_group, random = ~ 1 | sample_id, data = growth, na.action = na.omit, control = nlme::lmeControl(msMaxIter = 100)), silent = TRUE)
    fit2 <- try(nlme::lme(log_volume ~ splines::ns(day, df = 2) * scale(ecdf_rmse) + dose + end_timepoint_ploidy_group, random = ~ 1 | sample_id, data = growth[growth$dose_mg > 0, ], na.action = na.omit, control = nlme::lmeControl(msMaxIter = 100)), silent = TRUE)
    if (inherits(fit1, "try-error")) {
      fit1 <- NULL
    } else {
      mode1 <- "nlme_lme_random_intercept"
    }
    if (inherits(fit2, "try-error")) {
      fit2 <- NULL
    } else {
      mode2 <- "nlme_lme_random_intercept"
    }
  }
  if (is.null(fit1) && requireNamespace("splines", quietly = TRUE)) fit1 <- lm(log_volume ~ splines::ns(day, df = 3) * dose + end_timepoint_ploidy_group, data = growth)
  if (is.null(fit2) && requireNamespace("splines", quietly = TRUE)) fit2 <- lm(log_volume ~ splines::ns(day, df = 2) * scale(ecdf_rmse) + dose + end_timepoint_ploidy_group, data = growth[growth$dose_mg > 0, ])
  coefs <- rbind(extract_model_coefs(fit1, paste0("dose_growth_model_", mode1)), extract_model_coefs(fit2, paste0("treated_shift_growth_model_", mode2)))
  write_csv(coefs, file.path(dirs$GrowthModel, "growth_model_coefficients.csv"))
  auc_table <- unique(meta[, c("sample_id", "end_timepoint_ploidy_group", "dose", "dose_mg", "TGI_percent_auc", "tumor_volume_auc", "tumor_volume_auc_delta", "matched_control_mean_auc_delta", "auc_start_day", "auc_end_day", "auc_n_days"), drop = FALSE])
  auc_table <- merge(auc_table, results$CellCycle$equal_shift[, c("sample_id", "ecdf_rmse", "signed_mean_shift")], by = "sample_id", all.x = TRUE)
  write_csv(auc_table, file.path(dirs$GrowthModel, "model_based_TGI_or_AUC.csv"))
  dose_levels <- unique(meta$dose[order(meta$dose_mg)])
  p1 <- ggplot(growth, aes(x = day, y = tumor_volume, color = dose, group = sample_id)) +
    geom_line(alpha = 0.65, linewidth = 0.55) +
    geom_point(size = 1.4, alpha = 0.85) +
    facet_wrap(~ end_timepoint_ploidy_group, ncol = 1, scales = "free_y") +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = "Tumor growth curves by dose and EndTimePoint ploidy", x = "day", y = "tumor volume") +
    plot_theme
  save_plot(p1, dirs$GrowthModel, "growth_curves_by_dose", 8, 5.6)
  treated <- growth[growth$dose_mg > 0, , drop = FALSE]
  cut <- median(unique(auc_table$ecdf_rmse[auc_table$dose_mg > 0]), na.rm = TRUE)
  treated$shift_group <- ifelse(treated$ecdf_rmse >= cut, "high shift", "low shift")
  p2 <- ggplot(treated, aes(x = day, y = tumor_volume, color = shift_group, group = sample_id)) +
    geom_line(alpha = 0.65, linewidth = 0.55) +
    geom_point(size = 1.4, alpha = 0.85) +
    labs(title = "Treated tumor growth curves by CellCycle pseudotime-shift group", x = "day", y = "tumor volume", color = "shift group") +
    plot_theme
  save_plot(p2, dirs$GrowthModel, "growth_curves_by_shift_high_low", 7.2, 5)
  if (!is.null(fit2)) {
    pred_grid <- expand.grid(day = seq(min(treated$day, na.rm = TRUE), max(treated$day, na.rm = TRUE), length.out = 60), ecdf_rmse = quantile(treated$ecdf_rmse, c(0.25, 0.75), na.rm = TRUE), dose = unique(treated$dose)[1], end_timepoint_ploidy_group = unique(treated$end_timepoint_ploidy_group)[1])
    pred_grid$shift_group <- rep(c("low shift", "high shift"), each = 60)
    pred <- try(predict(fit2, newdata = pred_grid, level = 0), silent = TRUE)
    if (!inherits(pred, "try-error")) {
      pred_grid$predicted_volume <- expm1(pred)
      p3 <- ggplot(pred_grid, aes(x = day, y = predicted_volume, color = shift_group)) +
        geom_line(linewidth = 0.9) +
        labs(title = "Model-predicted treated growth by pseudotime shift", x = "day", y = "predicted tumor volume", color = "shift group") +
        plot_theme
      save_plot(p3, dirs$GrowthModel, "model_predicted_growth_by_shift", 6.6, 4.6)
    }
  }
  summary_lines <- c(
    "# Growth Model Summary",
    "",
    paste0("- Growth rows: ", nrow(growth), ". Samples: ", length(unique(growth$sample_id)), "."),
    paste0("- Dose growth model mode: ", mode1, "."),
    paste0("- Treated shift growth model mode: ", mode2, "."),
    "- These models are used as sensitivity support for AUC-based TGI, not as definitive mechanistic growth models."
  )
  write_text(summary_lines, file.path(dirs$GrowthModel, "growth_model_summary.md"))
  list(growth = growth, coefs = coefs, auc_table = auc_table)
}

fixed_effect_table <- function(fit, analysis, model_label, time_basis, n_rows, n_samples) {
  if (inherits(fit, "try-error") || is.null(fit)) return(data.frame())
  tt <- coef(summary(fit))
  p_col <- grep("^Pr\\(", colnames(tt), value = TRUE)[1]
  stat_col <- intersect(c("t value", "z value", "t-value"), colnames(tt))[1]
  se_col <- intersect(c("Std. Error", "Std.Error"), colnames(tt))[1]
  df_col <- intersect(c("df", "DF"), colnames(tt))[1]
  data.frame(
    analysis = analysis,
    model = model_label,
    time_basis = time_basis,
    term = rownames(tt),
    estimate = tt[, "Estimate"],
    std_error = if (!is.na(se_col)) tt[, se_col] else NA_real_,
    df = if (!is.na(df_col)) tt[, df_col] else NA_real_,
    statistic = if (!is.na(stat_col)) tt[, stat_col] else NA_real_,
    p_value = if (!is.na(p_col)) tt[, p_col] else NA_real_,
    n_rows = n_rows,
    n_samples = n_samples,
    singular_fit = if (inherits(fit, "merMod")) lme4::isSingular(fit, tol = 1e-4) else NA,
    stringsAsFactors = FALSE
  )
}

lrt_row <- function(reduced_fit, full_fit, analysis, reduced_model, full_model, time_basis, n_rows, n_samples) {
  if (inherits(reduced_fit, "try-error") || inherits(full_fit, "try-error") || is.null(reduced_fit) || is.null(full_fit)) {
    return(data.frame(
      analysis = analysis,
      time_basis = time_basis,
      reduced_model = reduced_model,
      full_model = full_model,
      n_rows = n_rows,
      n_samples = n_samples,
      reduced_npar = NA_real_, full_npar = NA_real_,
      reduced_aic = NA_real_, full_aic = NA_real_,
      reduced_bic = NA_real_, full_bic = NA_real_,
      reduced_logLik = NA_real_, full_logLik = NA_real_,
      chisq = NA_real_, df = NA_real_, p_value = NA_real_,
      stringsAsFactors = FALSE
    ))
  }
  tab <- as.data.frame(anova(reduced_fit, full_fit))
  p_col <- grep("^Pr\\(", names(tab), value = TRUE)[1]
  chisq_col <- grep("Chisq|L.Ratio", names(tab), value = TRUE)[1]
  df_col <- intersect(c("Df", "Chi Df"), names(tab))[1]
  data.frame(
    analysis = analysis,
    time_basis = time_basis,
    reduced_model = reduced_model,
    full_model = full_model,
    n_rows = n_rows,
    n_samples = n_samples,
    reduced_npar = if ("npar" %in% names(tab)) tab$npar[1] else NA_real_,
    full_npar = if ("npar" %in% names(tab)) tab$npar[2] else NA_real_,
    reduced_aic = if ("AIC" %in% names(tab)) tab$AIC[1] else NA_real_,
    full_aic = if ("AIC" %in% names(tab)) tab$AIC[2] else NA_real_,
    reduced_bic = if ("BIC" %in% names(tab)) tab$BIC[1] else NA_real_,
    full_bic = if ("BIC" %in% names(tab)) tab$BIC[2] else NA_real_,
    reduced_logLik = if ("logLik" %in% names(tab)) tab$logLik[1] else NA_real_,
    full_logLik = if ("logLik" %in% names(tab)) tab$logLik[2] else NA_real_,
    chisq = if (!is.na(chisq_col)) tab[[chisq_col]][2] else NA_real_,
    df = if (!is.na(df_col)) tab[[df_col]][2] else NA_real_,
    p_value = if (!is.na(p_col)) tab[[p_col]][2] else NA_real_,
    stringsAsFactors = FALSE
  )
}

emm_prediction_and_contrasts <- function(fit, growth, observed_days, control_label) {
  if (!requireNamespace("emmeans", quietly = TRUE)) {
    return(list(predictions = data.frame(), contrasts = data.frame()))
  }
  emm <- try(emmeans::emmeans(
    fit,
    ~ dose * end_timepoint_ploidy_group | day,
    at = list(day = observed_days),
    lmer.df = "asymptotic"
  ), silent = TRUE)
  if (inherits(emm, "try-error")) {
    return(list(predictions = data.frame(), contrasts = data.frame()))
  }
  pred <- as.data.frame(summary(emm, infer = c(TRUE, TRUE)))
  lcl_col <- intersect(c("asymp.LCL", "lower.CL"), names(pred))[1]
  ucl_col <- intersect(c("asymp.UCL", "upper.CL"), names(pred))[1]
  pred$predicted_relative_volume <- exp(pred$emmean)
  pred$predicted_relative_volume_lcl <- if (!is.na(lcl_col)) exp(pred[[lcl_col]]) else NA_real_
  pred$predicted_relative_volume_ucl <- if (!is.na(ucl_col)) exp(pred[[ucl_col]]) else NA_real_

  grid <- emm@grid
  linfct <- emm@linfct
  beta <- lme4::fixef(fit)
  vc <- as.matrix(stats::vcov(fit))
  common <- intersect(colnames(linfct), names(beta))
  linfct <- linfct[, common, drop = FALSE]
  beta <- beta[common]
  vc <- vc[common, common, drop = FALSE]
  dose_levels <- levels(growth$dose)
  ploidy_levels <- levels(growth$end_timepoint_ploidy_group)
  dose_map <- unique(growth[, c("dose", "dose_mg"), drop = FALSE])

  row_id <- function(day_value, dose_value, ploidy_value) {
    hit <- which(
      abs(safe_num(as.character(grid$day)) - day_value) < 1e-8 &
        as.character(grid$dose) == as.character(dose_value) &
        as.character(grid$end_timepoint_ploidy_group) == as.character(ploidy_value)
    )
    if (length(hit) == 0) NA_integer_ else hit[1]
  }

  estimate_l <- function(l_vec) {
    estimate <- as.numeric(sum(l_vec * beta))
    se <- as.numeric(sqrt(t(l_vec) %*% vc %*% l_vec))
    z <- estimate / se
    data.frame(
      estimate = estimate,
      std_error = se,
      statistic = z,
      p_value = 2 * stats::pnorm(-abs(z)),
      ci025 = estimate + stats::qnorm(0.025) * se,
      ci975 = estimate + stats::qnorm(0.975) * se,
      stringsAsFactors = FALSE
    )
  }

  rows <- list()
  treated_doses <- setdiff(dose_levels, control_label)
  for (day_value in observed_days) {
    for (ploidy_value in ploidy_levels) {
      control_idx <- row_id(day_value, control_label, ploidy_value)
      if (!is.finite(control_idx)) next
      for (dose_value in treated_doses) {
        treated_idx <- row_id(day_value, dose_value, ploidy_value)
        if (!is.finite(treated_idx)) next
        l_vec <- linfct[treated_idx, ] - linfct[control_idx, ]
        est <- estimate_l(l_vec)
        rows[[length(rows) + 1L]] <- data.frame(
          contrast_type = "treated_minus_control_within_ploidy",
          comparison = paste0(dose_value, "_minus_", control_label),
          dose = dose_value,
          dose_mg = dose_map$dose_mg[match(as.character(dose_value), as.character(dose_map$dose))],
          end_timepoint_ploidy_group = ploidy_value,
          day = day_value,
          estimate_log_relative_volume_diff = est$estimate,
          std_error = est$std_error,
          statistic = est$statistic,
          p_value = est$p_value,
          ci025 = est$ci025,
          ci975 = est$ci975,
          growth_ratio_treated_vs_control = exp(est$estimate),
          percent_inhibition_vs_control = 100 * (1 - exp(est$estimate)),
          stringsAsFactors = FALSE
        )
      }
    }
    if (length(ploidy_levels) >= 2) {
      ploidy_a <- ploidy_levels[1]
      ploidy_b <- ploidy_levels[2]
      control_a <- row_id(day_value, control_label, ploidy_a)
      control_b <- row_id(day_value, control_label, ploidy_b)
      if (is.finite(control_a) && is.finite(control_b)) {
        for (dose_value in treated_doses) {
          treated_a <- row_id(day_value, dose_value, ploidy_a)
          treated_b <- row_id(day_value, dose_value, ploidy_b)
          if (!is.finite(treated_a) || !is.finite(treated_b)) next
          l_vec <- (linfct[treated_b, ] - linfct[control_b, ]) - (linfct[treated_a, ] - linfct[control_a, ])
          est <- estimate_l(l_vec)
          rows[[length(rows) + 1L]] <- data.frame(
            contrast_type = "ploidy_difference_in_treatment_effect",
            comparison = paste0("(", dose_value, "_minus_", control_label, ")_", ploidy_b, "_minus_", ploidy_a),
            dose = dose_value,
            dose_mg = dose_map$dose_mg[match(as.character(dose_value), as.character(dose_map$dose))],
            end_timepoint_ploidy_group = paste0(ploidy_b, "_minus_", ploidy_a),
            day = day_value,
            estimate_log_relative_volume_diff = est$estimate,
            std_error = est$std_error,
            statistic = est$statistic,
            p_value = est$p_value,
            ci025 = est$ci025,
            ci975 = est$ci975,
            growth_ratio_treated_vs_control = exp(est$estimate),
            percent_inhibition_vs_control = 100 * (1 - exp(est$estimate)),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  contrasts <- if (length(rows) > 0) do.call(rbind, rows) else data.frame()
  list(predictions = pred, contrasts = contrasts)
}

run_growth_curve_mixed_model <- function(results, spline_df = 3L) {
  ensure_dir(dirs$GrowthCurveMixedModel)
  ensure_dir(file.path(dirs$GrowthCurveMixedModel, "tables"))
  ensure_dir(file.path(dirs$GrowthCurveMixedModel, "stats"))
  ensure_dir(file.path(dirs$GrowthCurveMixedModel, "figures"))

  meta <- results$CellCycle$sample_meta
  growth <- make_growth_curve_long(meta)
  growth <- merge(
    growth,
    results$CellCycle$equal_shift[, c("sample_id", "ecdf_rmse", "signed_mean_shift"), drop = FALSE],
    by = "sample_id",
    all.x = TRUE
  )
  growth <- growth[is.finite(growth$day) & is.finite(growth$tumor_volume) & growth$tumor_volume > 0 &
    is.finite(growth$log_relative_volume), , drop = FALSE]
  dose_levels <- unique(meta$dose[order(meta$dose_mg)])
  ploidy_levels <- unique(meta$end_timepoint_ploidy_group[order(meta$end_timepoint_ploidy_group_numeric)])
  growth$dose <- factor(growth$dose, levels = dose_levels)
  growth$end_timepoint_ploidy_group <- factor(growth$end_timepoint_ploidy_group, levels = ploidy_levels)
  growth$sample_id <- factor(growth$sample_id)
  write_csv(growth, file.path(dirs$GrowthCurveMixedModel, "tables", "tumor_growth_long_for_mixed_model.csv"))

  sample_design <- unique(growth[, c("sample_id", "dose", "dose_mg", "end_timepoint_ploidy_group")])
  observed_design <- aggregate(
    sample_id ~ dose + dose_mg + end_timepoint_ploidy_group,
    sample_design,
    function(x) length(unique(x))
  )
  names(observed_design)[names(observed_design) == "sample_id"] <- "n_samples"
  full_design <- expand.grid(
    dose = dose_levels,
    end_timepoint_ploidy_group = ploidy_levels,
    stringsAsFactors = FALSE
  )
  dose_map <- unique(sample_design[, c("dose", "dose_mg")])
  full_design$dose_mg <- dose_map$dose_mg[match(as.character(full_design$dose), as.character(dose_map$dose))]
  design_coverage <- merge(
    full_design,
    observed_design,
    by = c("dose", "dose_mg", "end_timepoint_ploidy_group"),
    all.x = TRUE
  )
  design_coverage$n_samples[is.na(design_coverage$n_samples)] <- 0L
  design_coverage$empty_combination <- design_coverage$n_samples == 0L
  design_coverage$version <- analysis_version
  write_csv(design_coverage, file.path(dirs$GrowthCurveMixedModel, "stats", "growth_curve_dose_ETP_design_coverage.csv"))

  if (!requireNamespace("lmerTest", quietly = TRUE) || !requireNamespace("lme4", quietly = TRUE) || !requireNamespace("splines", quietly = TRUE)) {
    write_text(c(
      "# Growth-Curve Mixed Model Summary",
      "",
      "Required packages `lmerTest`, `lme4`, or `splines` were unavailable, so the mixed-model analysis was not run."
    ), file.path(dirs$GrowthCurveMixedModel, "growth_curve_mixed_model_summary.md"))
    return(list(growth = growth, fixed_effects = data.frame(), lrt = data.frame(), predictions = data.frame(), contrasts = data.frame()))
  }

  n_rows <- nrow(growth)
  n_samples <- length(unique(growth$sample_id))
  control_label <- as.character(dose_levels[which(safe_num(gsub("[^0-9.]+", "", dose_levels)) == 0)[1]])
  if (!nzchar(control_label) || is.na(control_label)) control_label <- as.character(dose_levels[1])
  full_formula_label <- paste0("log_relative_volume ~ ns(day, df = ", spline_df, ") * dose * end_timepoint_ploidy_group + (1 | sample_id)")
  reduced_formula_label <- paste0("log_relative_volume ~ ns(day, df = ", spline_df, ") * dose + ns(day, df = ", spline_df, ") * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group + (1 | sample_id)")
  full_formula <- as.formula(paste0("log_relative_volume ~ splines::ns(day, df = ", spline_df, ") * dose * end_timepoint_ploidy_group + (1 | sample_id)"))
  reduced_formula <- as.formula(paste0("log_relative_volume ~ splines::ns(day, df = ", spline_df, ") * dose + splines::ns(day, df = ", spline_df, ") * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group + (1 | sample_id)"))
  linear_full_formula <- log_relative_volume ~ day * dose * end_timepoint_ploidy_group + (1 | sample_id)
  linear_reduced_formula <- log_relative_volume ~ day * dose + day * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group + (1 | sample_id)
  spline_fixed_full <- as.formula(paste0("~ splines::ns(day, df = ", spline_df, ") * dose * end_timepoint_ploidy_group"))
  spline_fixed_reduced <- as.formula(paste0("~ splines::ns(day, df = ", spline_df, ") * dose + splines::ns(day, df = ", spline_df, ") * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group"))
  linear_fixed_full <- ~ day * dose * end_timepoint_ploidy_group
  linear_fixed_reduced <- ~ day * dose + day * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group
  matrix_audit <- function(formula, analysis) {
    mm <- model.matrix(formula, data = growth)
    rank <- qr(mm)$rank
    data.frame(
      analysis = analysis,
      n_rows = nrow(mm),
      n_columns = ncol(mm),
      rank = rank,
      full_rank = rank == ncol(mm),
      empty_dose_ETP_combinations = sum(design_coverage$empty_combination),
      stringsAsFactors = FALSE
    )
  }
  design_audit <- rbind(
    matrix_audit(spline_fixed_full, "spline_full"),
    matrix_audit(spline_fixed_reduced, "spline_reduced"),
    matrix_audit(linear_fixed_full, "linear_full"),
    matrix_audit(linear_fixed_reduced, "linear_reduced")
  )
  design_audit$version <- analysis_version
  write_csv(design_audit, file.path(dirs$GrowthCurveMixedModel, "stats", "growth_curve_mixed_model_estimability_audit.csv"))
  ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000))

  fit_full <- try(lmerTest::lmer(full_formula, data = growth, REML = FALSE, control = ctrl), silent = TRUE)
  fit_reduced <- try(lmerTest::lmer(reduced_formula, data = growth, REML = FALSE, control = ctrl), silent = TRUE)
  fit_linear_full <- try(lmerTest::lmer(linear_full_formula, data = growth, REML = FALSE, control = ctrl), silent = TRUE)
  fit_linear_reduced <- try(lmerTest::lmer(linear_reduced_formula, data = growth, REML = FALSE, control = ctrl), silent = TRUE)

  fixed_effects <- rbind(
    fixed_effect_table(fit_full, "spline_full_time_dose_ploidy", full_formula_label, "natural_spline_df3", n_rows, n_samples),
    fixed_effect_table(fit_reduced, "spline_reduced_without_time_dose_ploidy", reduced_formula_label, "natural_spline_df3", n_rows, n_samples),
    fixed_effect_table(fit_linear_full, "linear_full_time_dose_ploidy", "log_relative_volume ~ day * dose * end_timepoint_ploidy_group + (1 | sample_id)", "linear_day", n_rows, n_samples),
    fixed_effect_table(fit_linear_reduced, "linear_reduced_without_time_dose_ploidy", "log_relative_volume ~ day * dose + day * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group + (1 | sample_id)", "linear_day", n_rows, n_samples)
  )
  if (nrow(fixed_effects) > 0) write_csv(fixed_effects, file.path(dirs$GrowthCurveMixedModel, "stats", "growth_curve_mixed_model_fixed_effects.csv"))

  lrt <- rbind(
    lrt_row(fit_reduced, fit_full, "primary_test_time_dose_ploidy_interaction", reduced_formula_label, full_formula_label, "natural_spline_df3", n_rows, n_samples),
    lrt_row(fit_linear_reduced, fit_linear_full, "sensitivity_linear_time_dose_ploidy_interaction", "log_relative_volume ~ day * dose + day * end_timepoint_ploidy_group + dose * end_timepoint_ploidy_group + (1 | sample_id)", "log_relative_volume ~ day * dose * end_timepoint_ploidy_group + (1 | sample_id)", "linear_day", n_rows, n_samples)
  )
  lrt$design_full_rank <- c(
    design_audit$full_rank[design_audit$analysis == "spline_full"] && design_audit$full_rank[design_audit$analysis == "spline_reduced"],
    design_audit$full_rank[design_audit$analysis == "linear_full"] && design_audit$full_rank[design_audit$analysis == "linear_reduced"]
  )
  lrt$estimable <- lrt$design_full_rank & sum(design_coverage$empty_combination) == 0L
  lrt$estimability_note <- ifelse(
    lrt$estimable,
    "complete dose-by-ETP design and full-rank fixed-effect matrices",
    "not estimable as the complete prespecified interaction because at least one dose-by-ETP combination is empty or the fixed-effect matrix is rank deficient"
  )
  lrt$p_value[!lrt$estimable] <- NA_real_
  write_csv(lrt, file.path(dirs$GrowthCurveMixedModel, "stats", "growth_curve_mixed_model_likelihood_ratio_tests.csv"))

  observed_days <- sort(unique(growth$day))
  pred_con <- if (!inherits(fit_full, "try-error") && !is.null(fit_full)) {
    emm_prediction_and_contrasts(fit_full, growth, observed_days, control_label)
  } else {
    list(predictions = data.frame(), contrasts = data.frame())
  }
  predictions <- pred_con$predictions
  contrasts <- pred_con$contrasts
  if (nrow(predictions) > 0) write_csv(predictions, file.path(dirs$GrowthCurveMixedModel, "tables", "growth_curve_mixed_model_predicted_relative_volume.csv"))
  if (nrow(contrasts) > 0) write_csv(contrasts, file.path(dirs$GrowthCurveMixedModel, "stats", "growth_curve_mixed_model_timepoint_contrasts.csv"))

  p_obs <- ggplot(growth, aes(x = day, y = relative_volume, color = dose, group = sample_id)) +
    geom_hline(yintercept = 1, color = "grey75", linewidth = 0.35) +
    geom_line(alpha = 0.65, linewidth = 0.55) +
    geom_point(size = 1.5, alpha = 0.85) +
    facet_wrap(~ end_timepoint_ploidy_group, ncol = 1) +
    scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
    labs(title = "Observed relative tumor-growth curves", x = "day", y = "relative tumor volume vs Day 0") +
    plot_theme
  save_plot(p_obs, file.path(dirs$GrowthCurveMixedModel, "figures"), "observed_relative_growth_curves_by_dose_ploidy", 8, 5.6)

  if (nrow(predictions) > 0) {
    p_pred <- ggplot(predictions, aes(x = day, y = predicted_relative_volume, color = dose, fill = dose)) +
      geom_hline(yintercept = 1, color = "grey75", linewidth = 0.35) +
      geom_ribbon(aes(ymin = predicted_relative_volume_lcl, ymax = predicted_relative_volume_ucl), alpha = 0.12, color = NA) +
      geom_line(linewidth = 0.9) +
      facet_wrap(~ end_timepoint_ploidy_group, ncol = 1) +
      scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
      scale_fill_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
      labs(title = "Spline mixed-model predicted relative tumor growth", x = "day", y = "model-predicted relative volume") +
      plot_theme
    save_plot(p_pred, file.path(dirs$GrowthCurveMixedModel, "figures"), "spline_mixed_model_predicted_relative_growth_by_dose_ploidy", 8, 5.6)
  }

  if (nrow(contrasts) > 0) {
    did <- contrasts[contrasts$contrast_type == "ploidy_difference_in_treatment_effect", , drop = FALSE]
    if (nrow(did) > 0) {
      p_did <- ggplot(did, aes(x = day, y = estimate_log_relative_volume_diff, color = dose, fill = dose)) +
        geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
        geom_ribbon(aes(ymin = ci025, ymax = ci975), alpha = 0.12, color = NA) +
        geom_line(linewidth = 0.85) +
        geom_point(size = 1.8) +
        scale_color_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
        scale_fill_manual(values = dose_cols_all[names(dose_cols_all) %in% dose_levels], name = "dose") +
        labs(
          title = "Ploidy difference in treatment effect over time",
          x = "day",
          y = "log growth effect difference: ETP-higher minus ETP-lower"
        ) +
        plot_theme
      save_plot(p_did, file.path(dirs$GrowthCurveMixedModel, "figures"), "ploidy_difference_in_treatment_effect_over_time", 7.4, 5)
    }
  }

  primary_lrt_p <- lrt$p_value[lrt$analysis == "primary_test_time_dose_ploidy_interaction"][1]
  final_day <- max(observed_days, na.rm = TRUE)
  final_did <- if (nrow(contrasts) > 0) {
    contrasts[contrasts$contrast_type == "ploidy_difference_in_treatment_effect" & contrasts$day == final_day, , drop = FALSE]
  } else data.frame()
  summary_lines <- c(
    "# Growth-Curve Mixed Model Summary",
    "",
    "This module addresses Question 2: whether longitudinal tumor-growth response to gemcitabine differs by EndTimePoint ploidy.",
    "",
    "Response variable:",
    "",
    "- `log_relative_volume = log(tumor_volume / tumor_volume_Day_0)`",
    "",
    "Primary model:",
    "",
    paste0("- `", full_formula_label, "`"),
    "",
    "Primary likelihood-ratio test:",
    "",
    paste0("- Full spline model vs reduced model without `ns(day):dose:end_timepoint_ploidy_group`; p = ", format_p(primary_lrt_p), "."),
    paste0("- Complete prespecified interaction estimable: ", bool_text(lrt$estimable[lrt$analysis == "primary_test_time_dose_ploidy_interaction"][1]), "."),
    paste0("- Empty dose-by-ETP combinations: ", sum(design_coverage$empty_combination), "."),
    "",
    "Timepoint contrasts:",
    "",
    "- `treated_minus_control_within_ploidy`: model-estimated treated-vs-control log relative-volume difference within ETP-lower or ETP-higher.",
    "- `ploidy_difference_in_treatment_effect`: `(ETP-higher treated - ETP-higher control) - (ETP-lower treated - ETP-lower control)` on the log-relative-volume scale.",
    "",
    "Interpretation note:",
    "",
    "- ETP is a post-treatment sample attribute. These models describe endpoint-defined association/effect modification and should not be interpreted as a baseline predictive biomarker.",
    "",
    "Interpretation rule: a nonsignificant difference test means no clear detected ploidy difference in longitudinal response; it does not prove equality unless an equivalence margin is separately specified."
  )
  if (nrow(final_did) > 0) {
    summary_lines <- c(
      summary_lines,
      "",
      paste0("Final observed day: Day ", final_day, "."),
      paste0(
        "- Day ", final_day, " ploidy difference-in-differences: ",
        paste0(final_did$dose, " estimate = ", format_num(final_did$estimate_log_relative_volume_diff), ", p = ", format_p(final_did$p_value), collapse = "; "),
        "."
      )
    )
  }
  write_text(summary_lines, file.path(dirs$GrowthCurveMixedModel, "growth_curve_mixed_model_summary.md"))

  list(growth = growth, fixed_effects = fixed_effects, lrt = lrt, predictions = predictions, contrasts = contrasts)
}

parse_group_ids <- function(ids) {
  parts <- strsplit(ids, "|", fixed = TRUE)
  sample_ids <- vapply(parts, `[`, character(1), 1)
  assignment_index <- match(sample_ids, etp_assignments$sample_id)
  if (anyNA(assignment_index)) {
    stop("Pseudobulk group IDs contain samples without ETP assignments: ", paste(unique(sample_ids[is.na(assignment_index)]), collapse = ", "), call. = FALSE)
  }
  data.frame(
    group_id = ids,
    sample_id = sample_ids,
    original_initial_ploidy = vapply(parts, `[`, character(1), 2),
    end_timepoint_ploidy_group = as.character(etp_assignments$end_timepoint_ploidy_group[assignment_index]),
    dose = vapply(parts, `[`, character(1), 3),
    group_feature = vapply(parts, function(x) paste(x[-seq_len(3)], collapse = "|"), character(1)),
    stringsAsFactors = FALSE
  )
}

read_pseudobulk <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  meta <- parse_group_ids(x$group_id)
  mat <- as.matrix(x[, setdiff(names(x), "group_id"), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- x$group_id
  list(meta = meta, mat = mat, path = path)
}

logcpm <- function(mat) {
  lib <- rowSums(mat, na.rm = TRUE)
  lib[!is.finite(lib) | lib <= 0] <- 1
  log2(sweep(mat, 1, lib / 1e6, "/") + 1)
}

module_definitions <- function(genes) {
  g <- function(x) intersect(paste0("GRCh38_", x), genes)
  list(
    E2F_S_phase_DNA_replication = g(c("MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "TYMS", "TK1", "DHFR", "CDC6", "ORC1")),
    G2M_checkpoint = g(c("CDK1", "CCNB1", "CCNB2", "CDC20", "CDC25C", "PLK1", "AURKA", "AURKB", "TOP2A", "UBE2C", "BIRC5")),
    DNA_damage_replication_stress = g(c("RPA1", "RPA2", "CLSPN", "ATR", "ATRIP", "CHEK1", "CHEK2", "WEE1", "TIMELESS", "TIPIN", "RAD51", "BRCA1", "BRCA2", "PARP1", "TP53BP1", "MDC1")),
    p53_apoptosis_senescence = g(c("GADD45A", "GADD45B", "GADD45G", "BAX", "BAK1", "BCL2", "BCL2L1", "BBC3", "PMAIP1", "CASP3", "CASP7", "CASP8", "CASP9", "CDKN1A", "CDKN2A", "MDM2", "SERPINE1", "LMNB1", "TP53I3", "BTG2", "DDB2")),
    Gemcitabine_transport_metabolism = g(c("SLC29A1", "SLC29A2", "SLC28A1", "DCK", "CDA", "DCTD", "CMPK1", "NME1", "NME2", "NT5C2", "NT5C3A", "RRM1", "RRM2", "RRM2B"))
  )
}

run_optional_expression_modules <- function(results) {
  external_base <- dirname(results_base)
  paths <- list(
    CellCycle_cluster = file.path(external_base, "04f_cell_cycle_results", "ROOT_6_END_NULL", "data", "pseudobulk_counts_by_cluster.tsv"),
    CellCycle_bin = file.path(external_base, "04f_cell_cycle_results", "ROOT_6_END_NULL", "data", "pseudobulk_counts_by_pseudotime_bin.tsv"),
    NonCellCycle_cluster = file.path(external_base, "04g_nonCellCycle_results", "ROOT_6_END_NULL", "data", "pseudobulk_counts_by_cluster.tsv"),
    NonCellCycle_bin = file.path(external_base, "04g_nonCellCycle_results", "ROOT_6_END_NULL", "data", "pseudobulk_counts_by_pseudotime_bin.tsv")
  )
  notes <- character()
  if (!run_optional_methods) {
    msg <- "# Optional Expression Modules\n\nOptional methods were disabled by `--run_optional_methods FALSE`."
    write_text(msg, file.path(dirs$GenePrograms, "pseudotime_direction_validation_summary.md"))
    write_text(msg, file.path(dirs$Pseudobulk, "pseudobulk_DE_summary.md"))
    return(list(notes = msg))
  }
  pb <- lapply(paths, read_pseudobulk)
  missing <- names(pb)[vapply(pb, is.null, logical(1))]
  if (length(missing) > 0) notes <- c(notes, paste0("Missing pseudobulk inputs: ", paste(missing, collapse = ", "), "."))
  if (is.null(pb$CellCycle_bin)) {
    msg <- c("# Pseudotime Direction Validation", "", "Skipped: CellCycle pseudotime-bin pseudobulk counts were not found.", notes)
    write_text(msg, file.path(dirs$GenePrograms, "pseudotime_direction_validation_summary.md"))
  } else {
    lcp <- logcpm(pb$CellCycle_bin$mat)
    modules <- module_definitions(colnames(lcp))
    scores <- do.call(cbind, lapply(modules, function(genes) {
      if (length(genes) == 0) rep(NA_real_, nrow(lcp)) else rowMeans(lcp[, genes, drop = FALSE], na.rm = TRUE)
    }))
    score_df <- cbind(pb$CellCycle_bin$meta, as.data.frame(scores, check.names = FALSE))
    score_df$bin_numeric <- safe_num(sub("^bin_", "", score_df$group_feature))
    score_df$treatment_group <- ifelse(score_df$dose == "0mg/kg", "control", "treated")
    score_df <- merge(score_df, results$CellCycle$equal_shift[, c("sample_id", "TGI_percent_auc", "ecdf_rmse")], by = "sample_id", all.x = TRUE)
    score_df$aggregation_level <- "sample_pseudotime_bin"
    write_csv(score_df, file.path(dirs$GenePrograms, "cellcycle_module_scores_by_cell.csv"))
    trend_rows <- list()
    for (module in names(modules)) {
      for (group in c("control", "treated")) {
        local <- score_df[score_df$treatment_group == group, , drop = FALSE]
        cr <- safe_cor_result(local$bin_numeric, local[[module]], "spearman")
        trend_rows[[length(trend_rows) + 1L]] <- data.frame(module = module, analysis = paste0("score_vs_pseudotime_", group), n = cr$n, estimate = cr$estimate, p_value = cr$p_value, stringsAsFactors = FALSE)
      }
      fit <- try(lm(score_df[[module]] ~ score_df$bin_numeric * score_df$treatment_group), silent = TRUE)
      p_int <- if (!inherits(fit, "try-error")) {
        cf <- coef(summary(fit))
        idx <- grep(":", rownames(cf))
        if (length(idx) > 0) cf[idx[1], "Pr(>|t|)"] else NA_real_
      } else NA_real_
      trend_rows[[length(trend_rows) + 1L]] <- data.frame(module = module, analysis = "pseudotime_by_treatment_interaction_lm", n = nrow(score_df), estimate = NA_real_, p_value = p_int, stringsAsFactors = FALSE)
    }
    trends <- do.call(rbind, trend_rows)
    write_csv(trends, file.path(dirs$GenePrograms, "module_score_pseudotime_trends.csv"))
    write_csv(trends[grep("score_vs_pseudotime|interaction", trends$analysis), ], file.path(dirs$GenePrograms, "pseudotime_gene_trends.csv"))
    write_csv(trends[grep("treatment_interaction", trends$analysis), ], file.path(dirs$GenePrograms, "pseudotime_treatment_interaction_trends.csv"))
    tgi_rows <- list()
    for (module in names(modules)) {
      local <- score_df[score_df$treatment_group == "treated", , drop = FALSE]
      fit <- try(lm(local[[module]] ~ local$bin_numeric * local$TGI_percent_auc), silent = TRUE)
      p_int <- if (!inherits(fit, "try-error")) {
        cf <- coef(summary(fit))
        idx <- grep(":", rownames(cf))
        if (length(idx) > 0) cf[idx[1], "Pr(>|t|)"] else NA_real_
      } else NA_real_
      tgi_rows[[length(tgi_rows) + 1L]] <- data.frame(module = module, analysis = "pseudotime_by_TGI_interaction_lm", n = nrow(local), p_value = p_int, stringsAsFactors = FALSE)
    }
    tgi_trends <- do.call(rbind, tgi_rows)
    write_csv(tgi_trends, file.path(dirs$GenePrograms, "pseudotime_TGI_interaction_trends.csv"))
    heat <- reshape(score_df[, c("sample_id", "dose", "bin_numeric", names(modules))], varying = names(modules), v.names = "score", timevar = "module", times = names(modules), direction = "long")
    heat_summary <- aggregate(score ~ dose + bin_numeric + module, heat, function(x) mean(x, na.rm = TRUE))
    p_heat <- ggplot(heat_summary, aes(x = bin_numeric, y = module, fill = score)) +
      geom_tile(color = "white", linewidth = 0.2) +
      facet_wrap(~ dose) +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "score") +
      labs(title = "CellCycle module scores along pseudotime", x = "pseudotime bin", y = NULL) +
      plot_theme
    save_plot(p_heat, dirs$GenePrograms, "module_scores_along_pseudotime_heatmap", 8, 4.8)
    p_line <- ggplot(heat[heat$module %in% c("E2F_S_phase_DNA_replication", "G2M_checkpoint", "DNA_damage_replication_stress"), ], aes(x = bin_numeric, y = score, color = dose)) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
      facet_wrap(~ module, scales = "free_y") +
      labs(title = "Cell-cycle and DNA-replication module scores along pseudotime", x = "pseudotime bin", y = "mean logCPM module score") +
      plot_theme
    save_plot(p_line, dirs$GenePrograms, "E2F_G2M_DNAreplication_scores_along_pseudotime", 8, 4.8)
    p_gem <- ggplot(heat[heat$module == "Gemcitabine_transport_metabolism", ], aes(x = bin_numeric, y = score, color = dose)) +
      stat_summary(fun = mean, geom = "line", linewidth = 0.8) +
      labs(title = "Gemcitabine metabolism module along pseudotime", x = "pseudotime bin", y = "mean logCPM module score") +
      plot_theme
    save_plot(p_gem, dirs$GenePrograms, "gemcitabine_metabolism_genes_along_pseudotime", 6.6, 4.4)
    write_text(c(
      "# Pseudotime Direction Validation",
      "",
      "Raw cell-level expression was not provided inside `--input_root`; module trends were therefore computed from sample x pseudotime-bin pseudobulk counts when available.",
      paste0("- Available modules: ", paste(names(modules), collapse = ", "), "."),
      paste0("- Pseudobulk source: ", pb$CellCycle_bin$path),
      "",
      "Interpretation should be limited to aggregated pseudobulk/module trends, not single-cell gene-level dynamics."
    ), file.path(dirs$GenePrograms, "pseudotime_direction_validation_summary.md"))
  }
  run_pseudobulk_de(pb, results, notes)
}

limma_de <- function(count_mat_by_sample, meta, group_col, coef_name = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) return(data.frame())
  lcp <- t(logcpm(count_mat_by_sample))
  meta <- meta[match(colnames(lcp), meta$sample_id), , drop = FALSE]
  design <- model.matrix(as.formula(paste0("~ ", group_col)), data = meta)
  fit <- limma::eBayes(limma::lmFit(lcp, design))
  coef_idx <- if (!is.null(coef_name) && coef_name %in% colnames(design)) coef_name else colnames(design)[2]
  tt <- limma::topTable(fit, coef = coef_idx, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt[, c("gene", setdiff(names(tt), "gene"))]
}

aggregate_pseudobulk_by_sample <- function(pb) {
  mat <- pb$mat
  meta <- pb$meta
  sample_mat <- rowsum(mat, group = meta$sample_id, reorder = FALSE)
  sample_meta <- meta[!duplicated(meta$sample_id), c("sample_id", "end_timepoint_ploidy_group", "dose"), drop = FALSE]
  sample_meta$dose_mg <- safe_num(sub("mg/kg", "", sample_meta$dose))
  sample_meta$treatment_group <- ifelse(sample_meta$dose_mg == 0, "control", "treated")
  list(mat = sample_mat, meta = sample_meta)
}

run_pseudobulk_de <- function(pb, results, notes = character()) {
  if (is.null(pb$CellCycle_cluster) || is.null(pb$NonCellCycle_cluster)) {
    write_text(c("# Pseudobulk DE Summary", "", "Skipped or partial: required cluster-level pseudobulk counts were not all found.", notes), file.path(dirs$Pseudobulk, "pseudobulk_DE_summary.md"))
    return(invisible(NULL))
  }
  cc <- aggregate_pseudobulk_by_sample(pb$CellCycle_cluster)
  nc <- aggregate_pseudobulk_by_sample(pb$NonCellCycle_cluster)
  cc$meta <- merge(cc$meta, results$CellCycle$equal_shift[, c("sample_id", "TGI_percent_auc")], by = "sample_id", all.x = TRUE)
  nc$meta <- merge(nc$meta, results$NonCellCycle$equal_shift[, c("sample_id", "TGI_percent_auc")], by = "sample_id", all.x = TRUE)
  write_csv(rbind(data.frame(compartment = "CellCycle", cc$meta), data.frame(compartment = "NonCellCycle", nc$meta)), file.path(dirs$Pseudobulk, "pseudobulk_sample_metadata.csv"))
  cc_de <- limma_de(cc$mat, cc$meta, "treatment_group")
  nc_de <- limma_de(nc$mat, nc$meta, "treatment_group")
  treated_meta <- cc$meta[cc$meta$treatment_group == "treated", , drop = FALSE]
  treated_meta$TGI_group <- ifelse(treated_meta$TGI_percent_auc >= median(treated_meta$TGI_percent_auc, na.rm = TRUE), "highTGI", "lowTGI")
  cc_high <- limma_de(cc$mat[treated_meta$sample_id, , drop = FALSE], treated_meta, "TGI_group")
  write_csv(cc_de, file.path(dirs$Pseudobulk, "CellCycle_treated_vs_control_DE.csv"))
  write_csv(cc_high, file.path(dirs$Pseudobulk, "CellCycle_highTGI_vs_lowTGI_DE.csv"))
  write_csv(nc_de, file.path(dirs$Pseudobulk, "NonCellCycle_treated_vs_control_DE.csv"))
  modules <- module_definitions(colnames(cc$mat))
  enrichment <- do.call(rbind, lapply(names(modules), function(module) {
    genes <- intersect(modules[[module]], cc_de$gene)
    if (length(genes) < 2) return(data.frame(module = module, n_genes = length(genes), mean_t = NA_real_, p_value = NA_real_))
    in_set <- cc_de$t[match(genes, cc_de$gene)]
    out_set <- cc_de$t[!cc_de$gene %in% genes]
    tt <- suppressWarnings(t.test(in_set, out_set))
    data.frame(module = module, n_genes = length(genes), mean_t = mean(in_set, na.rm = TRUE), p_value = tt$p.value, stringsAsFactors = FALSE)
  }))
  enrichment$q_value <- p.adjust(enrichment$p_value, method = "BH")
  write_csv(enrichment, file.path(dirs$Pseudobulk, "pathway_enrichment_results.csv"))
  volcano <- function(de, filename, title) {
    if (nrow(de) == 0) return(invisible(NULL))
    de$neg_log10_p <- -log10(pmax(de$P.Value, .Machine$double.xmin))
    p <- ggplot(de, aes(x = logFC, y = neg_log10_p)) +
      geom_point(alpha = 0.75, color = "#377eb8") +
      geom_hline(yintercept = -log10(0.05), color = "grey60", linewidth = 0.35) +
      labs(title = title, x = "logFC", y = "-log10 p") +
      plot_theme
    save_plot(p, dirs$Pseudobulk, filename, 6.4, 4.6)
  }
  volcano(cc_de, "CellCycle_treated_vs_control_volcano", "CellCycle pseudobulk: treated vs control")
  volcano(cc_high, "CellCycle_highTGI_vs_lowTGI_volcano", "CellCycle pseudobulk: high vs low TGI")
  p_en <- ggplot(enrichment, aes(x = mean_t, y = reorder(module, mean_t), size = n_genes, color = -log10(pmax(q_value, .Machine$double.xmin)))) +
    geom_point() +
    scale_color_gradient(low = "#2166ac", high = "#b2182b", name = "-log10 q") +
    labs(title = "Module-level enrichment summary from CellCycle pseudobulk DE", x = "mean moderated t", y = NULL) +
    plot_theme
  save_plot(p_en, dirs$Pseudobulk, "pathway_enrichment_dotplot", 7, 4.6)
  gem_genes <- intersect(module_definitions(colnames(cc$mat))$Gemcitabine_transport_metabolism, colnames(cc$mat))
  if (length(gem_genes) > 0) {
    lcp <- logcpm(cc$mat)[, gem_genes, drop = FALSE]
    heat <- cbind(cc$meta, as.data.frame(lcp[cc$meta$sample_id, , drop = FALSE], check.names = FALSE))
    heat_l <- reshape(heat[, c("sample_id", "dose", gem_genes)], varying = gem_genes, v.names = "logcpm", timevar = "gene", times = gem_genes, direction = "long")
    p_h <- ggplot(heat_l, aes(x = sample_id, y = gene, fill = logcpm)) +
      geom_tile() +
      facet_grid(. ~ dose, scales = "free_x", space = "free_x") +
      scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", name = "logCPM") +
      labs(title = "Gemcitabine-related genes in CellCycle pseudobulk", x = "sample", y = NULL) +
      plot_theme + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    save_plot(p_h, dirs$Pseudobulk, "gemcitabine_gene_heatmap", 8, 4.8)
  }
  write_text(c(
    "# Pseudobulk DE Summary",
    "",
    "Pseudobulk DE was run on sample-level sums of the 04f/04g targeted pseudobulk count tables.",
    "The available pseudobulk tables contain a targeted gene panel, not a genome-wide matrix, so pathway results are module summaries over available genes.",
    paste0("- CellCycle genes tested: ", nrow(cc_de), "."),
    paste0("- NonCellCycle genes tested: ", nrow(nc_de), ".")
  ), file.path(dirs$Pseudobulk, "pseudobulk_DE_summary.md"))
  invisible(list(cc_de = cc_de, nc_de = nc_de, enrichment = enrichment))
}

neighborhood_fallback <- function(results) {
  combined <- rbind(
    results$CellCycle$df[, c("cell_id", "sample_id", "cluster", "end_timepoint_ploidy_group", "dose", "dose_mg", "pseudotime", "compartment")],
    results$NonCellCycle$df[, c("cell_id", "sample_id", "cluster", "end_timepoint_ploidy_group", "dose", "dose_mg", "pseudotime", "compartment")]
  )
  breaks <- unique(quantile(combined$pseudotime, probs = seq(0, 1, length.out = 11), na.rm = TRUE))
  combined$neighborhood <- cut(combined$pseudotime, breaks = breaks, include.lowest = TRUE, labels = paste0("pt_bin_", seq_len(length(breaks) - 1)))
  counts <- aggregate(cell_id ~ sample_id + end_timepoint_ploidy_group + dose + dose_mg + compartment + neighborhood, combined, length)
  names(counts)[names(counts) == "cell_id"] <- "cells"
  totals <- aggregate(cells ~ sample_id + compartment, counts, sum)
  names(totals)[3] <- "compartment_cells"
  counts <- merge(counts, totals, by = c("sample_id", "compartment"))
  counts$proportion <- counts$cells / counts$compartment_cells
  metadata <- aggregate(pseudotime ~ compartment + neighborhood, combined, mean)
  names(metadata)[3] <- "mean_pseudotime"
  metadata$n_cells <- aggregate(cell_id ~ compartment + neighborhood, combined, length)$cell_id
  write_csv(metadata, file.path(dirs$NeighborhoodDA, "neighborhood_metadata.csv"))
  treat_rows <- lapply(split(counts, paste(counts$compartment, counts$neighborhood, sep = "|")), function(tab) {
    gp <- ifelse(tab$dose_mg == 0, "control", "treated")
    wt <- safe_wilcox(tab$proportion, gp)
    data.frame(compartment = tab$compartment[1], neighborhood = tab$neighborhood[1], n = wt$n, treated_minus_control_mean_proportion = mean(tab$proportion[gp == "treated"], na.rm = TRUE) - mean(tab$proportion[gp == "control"], na.rm = TRUE), p_value = wt$p_value, stringsAsFactors = FALSE)
  })
  treat <- do.call(rbind, treat_rows)
  treat$q_value <- p.adjust(treat$p_value, method = "BH")
  tgi_meta <- rbind(
    results$CellCycle$equal_shift[, c("sample_id", "TGI_percent_auc")],
    results$NonCellCycle$equal_shift[, c("sample_id", "TGI_percent_auc")]
  )
  tgi_meta <- tgi_meta[!duplicated(tgi_meta$sample_id), ]
  counts_tgi <- merge(counts[counts$dose_mg > 0, ], tgi_meta, by = "sample_id", all.x = TRUE)
  tgi_rows <- lapply(split(counts_tgi, paste(counts_tgi$compartment, counts_tgi$neighborhood, sep = "|")), function(tab) {
    cr <- safe_cor_result(tab$proportion, tab$TGI_percent_auc, "spearman")
    data.frame(compartment = tab$compartment[1], neighborhood = tab$neighborhood[1], n = cr$n, spearman_rho = cr$estimate, p_value = cr$p_value, stringsAsFactors = FALSE)
  })
  tgi <- do.call(rbind, tgi_rows)
  tgi$q_value <- p.adjust(tgi$p_value, method = "BH")
  alignment <- merge(metadata, treat, by = c("compartment", "neighborhood"), all.x = TRUE)
  write_csv(treat, file.path(dirs$NeighborhoodDA, "treatment_DA_results.csv"))
  write_csv(tgi, file.path(dirs$NeighborhoodDA, "TGI_DA_results.csv"))
  write_csv(alignment, file.path(dirs$NeighborhoodDA, "neighborhood_pseudotime_alignment.csv"))
  p1 <- ggplot(alignment, aes(x = mean_pseudotime, y = treated_minus_control_mean_proportion, color = compartment)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(aes(size = -log10(pmax(p_value, .Machine$double.xmin)))) +
    labs(title = "Pseudotime-bin neighborhood DA fallback: treatment", x = "mean pseudotime", y = "treated-control proportion", size = "-log10 p") +
    plot_theme
  save_plot(p1, dirs$NeighborhoodDA, "treatment_DA_on_pseudotime_or_umap", 7, 4.8)
  tgi_align <- merge(metadata, tgi, by = c("compartment", "neighborhood"), all.x = TRUE)
  p2 <- ggplot(tgi_align, aes(x = mean_pseudotime, y = spearman_rho, color = compartment)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(aes(size = -log10(pmax(p_value, .Machine$double.xmin)))) +
    labs(title = "Pseudotime-bin neighborhood DA fallback: TGI", x = "mean pseudotime", y = "Spearman rho with TGI", size = "-log10 p") +
    plot_theme
  save_plot(p2, dirs$NeighborhoodDA, "TGI_DA_on_pseudotime_or_umap", 7, 4.8)
  p3 <- ggplot(alignment, aes(x = mean_pseudotime, y = treated_minus_control_mean_proportion, color = compartment)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_line() +
    geom_point() +
    labs(title = "Neighborhood logFC proxy along pseudotime", x = "mean pseudotime", y = "treated-control proportion") +
    plot_theme
  save_plot(p3, dirs$NeighborhoodDA, "neighborhood_logFC_vs_pseudotime", 7, 4.8)
  write_text(c(
    "# Neighborhood DA Summary",
    "",
    "MiloR was not available in the current R environment, so this analysis used pseudotime bins as a fallback neighborhood definition.",
    "Counts were summarized per sample, compartment, and pseudotime bin; tests use sample-level proportions.",
    "This is a robustness analysis for one-dimensional state occupancy, not a replacement for graph-neighborhood DA."
  ), file.path(dirs$NeighborhoodDA, "neighborhood_DA_summary.md"))
  list(treatment = treat, tgi = tgi, metadata = metadata)
}

write_manuscript_outputs <- function(results, robustness, comparison, composition, growth, dose_specific = NULL, pseudotime_ploidy = NULL, auc_tgi = NULL, ploidy_response = NULL, growth_curve_mixed = NULL) {
  cc_dose <- results$CellCycle$dose_tests
  nc_dose <- results$NonCellCycle$dose_tests
  get_dose_p <- function(tab, comparison, strat = FALSE) {
    hit <- tab$comparison == comparison & tab$stratified_by_end_timepoint_ploidy_group == strat
    if (any(hit)) tab$p_ecdf_rmse[which(hit)[1]] else NA_real_
  }
  first_value <- function(tab, col) {
    if (is.null(tab) || nrow(tab) == 0 || !(col %in% names(tab))) return(NA_real_)
    tab[[col]][1]
  }
  cc_primary <- results$CellCycle$tgi_stats[results$CellCycle$tgi_stats$pre_specified_primary, , drop = FALSE]
  nc_primary <- results$NonCellCycle$tgi_stats[results$NonCellCycle$tgi_stats$reference_type == "primary_equal_sample_reference" & results$NonCellCycle$tgi_stats$shift_metric == "ecdf_rmse" & results$NonCellCycle$tgi_stats$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  delta <- comparison$delta
  dose_comp <- if (!is.null(dose_specific)) dose_specific$comparisons else data.frame()
  dose_cor <- if (!is.null(dose_specific)) dose_specific$correlations else data.frame()
  dose_int <- if (!is.null(dose_specific)) dose_specific$interactions else data.frame()
  pseudotime_ploidy_tests <- if (!is.null(pseudotime_ploidy)) pseudotime_ploidy$tests else data.frame()
  cc_auc_dose <- if (nrow(dose_comp) > 0) {
    dose_comp[dose_comp$compartment == "CellCycle" & dose_comp$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  } else data.frame()
  cc_endpoint_dose <- if (nrow(dose_comp) > 0) {
    dose_comp[dose_comp$compartment == "CellCycle" & grepl("^TGI_percent_Day_", dose_comp$tgi_measure), , drop = FALSE]
  } else data.frame()
  nc_auc_dose <- if (nrow(dose_comp) > 0) {
    dose_comp[dose_comp$compartment == "NonCellCycle" & dose_comp$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  } else data.frame()
  nc_endpoint_dose <- if (nrow(dose_comp) > 0) {
    dose_comp[dose_comp$compartment == "NonCellCycle" & grepl("^TGI_percent_Day_", dose_comp$tgi_measure), , drop = FALSE]
  } else data.frame()
  cc_30_cor <- if (nrow(dose_cor) > 0) {
    dose_cor[dose_cor$compartment == "CellCycle" & dose_cor$dose_mg == 30 & dose_cor$pre_specified_primary, , drop = FALSE]
  } else data.frame()
  cc_120_cor <- if (nrow(dose_cor) > 0) {
    dose_cor[dose_cor$compartment == "CellCycle" & dose_cor$dose_mg == 120 & dose_cor$pre_specified_primary, , drop = FALSE]
  } else data.frame()
  nc_30_cor <- if (nrow(dose_cor) > 0) {
    dose_cor[dose_cor$compartment == "NonCellCycle" & dose_cor$dose_mg == 30 & dose_cor$shift_metric == "ecdf_rmse" & dose_cor$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  } else data.frame()
  nc_120_cor <- if (nrow(dose_cor) > 0) {
    dose_cor[dose_cor$compartment == "NonCellCycle" & dose_cor$dose_mg == 120 & dose_cor$shift_metric == "ecdf_rmse" & dose_cor$tgi_measure == "TGI_percent_auc", , drop = FALSE]
  } else data.frame()
  int_hit <- if (nrow(dose_int) > 0) {
    dose_int[dose_int$compartment == "CellCycle" & grepl("scale\\(ecdf_rmse\\).*:factor\\(dose_mg\\)|factor\\(dose_mg\\).*:scale\\(ecdf_rmse\\)", dose_int$term), , drop = FALSE]
  } else data.frame()
  nc_int_hit <- if (nrow(dose_int) > 0) {
    dose_int[dose_int$compartment == "NonCellCycle" & grepl("scale\\(ecdf_rmse\\).*:factor\\(dose_mg\\)|factor\\(dose_mg\\).*:scale\\(ecdf_rmse\\)", dose_int$term), , drop = FALSE]
  } else data.frame()
  cc_pseudotime_ploidy_did <- if (nrow(pseudotime_ploidy_tests) > 0) pseudotime_ploidy_tests[pseudotime_ploidy_tests$compartment == "CellCycle", , drop = FALSE] else data.frame()
  nc_pseudotime_ploidy_did <- if (nrow(pseudotime_ploidy_tests) > 0) pseudotime_ploidy_tests[pseudotime_ploidy_tests$compartment == "NonCellCycle", , drop = FALSE] else data.frame()
  endpoint_min_p <- if (nrow(cc_endpoint_dose) > 0) min(cc_endpoint_dose$permutation_p_two_sided, na.rm = TRUE) else NA_real_
  endpoint_min_q <- if (nrow(cc_endpoint_dose) > 0) min(cc_endpoint_dose$permutation_q_exploratory, na.rm = TRUE) else NA_real_
  nc_endpoint_min_p <- if (nrow(nc_endpoint_dose) > 0) min(nc_endpoint_dose$permutation_p_two_sided, na.rm = TRUE) else NA_real_
  nc_endpoint_min_q <- if (nrow(nc_endpoint_dose) > 0) min(nc_endpoint_dose$permutation_q_exploratory, na.rm = TRUE) else NA_real_
  ploidy_ecdf <- if (!is.null(ploidy_response)) ploidy_response$ecdf else data.frame()
  ploidy_shift <- if (!is.null(ploidy_response)) ploidy_response$shift else data.frame()
  ploidy_tgi <- if (!is.null(ploidy_response)) ploidy_response$tgi else data.frame()
  ploidy_models <- if (!is.null(ploidy_response)) ploidy_response$models else data.frame()
  auc_models <- if (!is.null(auc_tgi)) auc_tgi$models else data.frame()
  growth_lrt <- if (!is.null(growth_curve_mixed)) growth_curve_mixed$lrt else data.frame()
  growth_contrasts <- if (!is.null(growth_curve_mixed)) growth_curve_mixed$contrasts else data.frame()
  cc_ecdf_2n <- if (nrow(ploidy_ecdf) > 0) ploidy_ecdf[ploidy_ecdf$compartment == "CellCycle" & ploidy_ecdf$end_timepoint_ploidy_group == "ETP-lower", , drop = FALSE] else data.frame()
  cc_ecdf_4n <- if (nrow(ploidy_ecdf) > 0) ploidy_ecdf[ploidy_ecdf$compartment == "CellCycle" & ploidy_ecdf$end_timepoint_ploidy_group == "ETP-higher", , drop = FALSE] else data.frame()
  nc_ecdf_2n <- if (nrow(ploidy_ecdf) > 0) ploidy_ecdf[ploidy_ecdf$compartment == "NonCellCycle" & ploidy_ecdf$end_timepoint_ploidy_group == "ETP-lower", , drop = FALSE] else data.frame()
  nc_ecdf_4n <- if (nrow(ploidy_ecdf) > 0) ploidy_ecdf[ploidy_ecdf$compartment == "NonCellCycle" & ploidy_ecdf$end_timepoint_ploidy_group == "ETP-higher", , drop = FALSE] else data.frame()
  cc_shift_ploidy <- if (nrow(ploidy_shift) > 0) ploidy_shift[ploidy_shift$compartment == "CellCycle" & ploidy_shift$response_metric == "ecdf_rmse", , drop = FALSE] else data.frame()
  nc_shift_ploidy <- if (nrow(ploidy_shift) > 0) ploidy_shift[ploidy_shift$compartment == "NonCellCycle" & ploidy_shift$response_metric == "ecdf_rmse", , drop = FALSE] else data.frame()
  tgi_auc_ploidy <- if (nrow(ploidy_tgi) > 0) ploidy_tgi[ploidy_tgi$tgi_measure == "TGI_percent_auc", , drop = FALSE] else data.frame()
  cc_shift_interaction <- if (nrow(ploidy_models) > 0) ploidy_models[ploidy_models$compartment == "CellCycle" & ploidy_models$analysis == "shift_treatment_by_ploidy_interaction" & grepl(":", ploidy_models$term), , drop = FALSE] else data.frame()
  cc_tgi_ploidy <- if (nrow(ploidy_models) > 0) ploidy_models[ploidy_models$compartment == "CellCycle" & ploidy_models$analysis == "treated_TGI_by_dose_and_ploidy" & grepl("end_timepoint_ploidy_group", ploidy_models$term), , drop = FALSE] else data.frame()
  cc_shift_tgi_ploidy <- if (nrow(ploidy_models) > 0) ploidy_models[ploidy_models$compartment == "CellCycle" & ploidy_models$analysis == "treated_TGI_shift_by_ploidy_interaction" & grepl(":", ploidy_models$term), , drop = FALSE] else data.frame()
  cc_auc_adjusted_shift <- if (nrow(auc_models) > 0) auc_models[auc_models$compartment == "CellCycle" & auc_models$analysis == "shift_ploidy_dose_adjusted" & auc_models$term == "scale(ecdf_rmse)", , drop = FALSE] else data.frame()
  cc_auc_adjusted_ploidy <- if (nrow(auc_models) > 0) auc_models[auc_models$compartment == "CellCycle" & auc_models$analysis == "shift_ploidy_dose_adjusted" & grepl("end_timepoint_ploidy_group", auc_models$term), , drop = FALSE] else data.frame()
  cc_auc_shift_ploidy_interaction <- if (nrow(auc_models) > 0) auc_models[auc_models$compartment == "CellCycle" & auc_models$analysis == "shift_by_ploidy_dose_adjusted" & grepl(":", auc_models$term), , drop = FALSE] else data.frame()
  nc_auc_adjusted_shift <- if (nrow(auc_models) > 0) auc_models[auc_models$compartment == "NonCellCycle" & auc_models$analysis == "shift_ploidy_dose_adjusted" & auc_models$term == "scale(ecdf_rmse)", , drop = FALSE] else data.frame()
  nc_auc_shift_ploidy_interaction <- if (nrow(auc_models) > 0) auc_models[auc_models$compartment == "NonCellCycle" & auc_models$analysis == "shift_by_ploidy_dose_adjusted" & grepl(":", auc_models$term), , drop = FALSE] else data.frame()
  growth_primary_lrt <- if (nrow(growth_lrt) > 0) growth_lrt[growth_lrt$analysis == "primary_test_time_dose_ploidy_interaction", , drop = FALSE] else data.frame()
  growth_linear_lrt <- if (nrow(growth_lrt) > 0) growth_lrt[growth_lrt$analysis == "sensitivity_linear_time_dose_ploidy_interaction", , drop = FALSE] else data.frame()
  growth_final_day <- if (nrow(growth_contrasts) > 0) max(growth_contrasts$day, na.rm = TRUE) else NA_real_
  growth_final_did <- if (nrow(growth_contrasts) > 0 && is.finite(growth_final_day)) growth_contrasts[growth_contrasts$contrast_type == "ploidy_difference_in_treatment_effect" & growth_contrasts$day == growth_final_day, , drop = FALSE] else data.frame()
  did30 <- if (nrow(growth_final_did) > 0) growth_final_did[growth_final_did$dose_mg == 30, , drop = FALSE] else data.frame()
  did120 <- if (nrow(growth_final_did) > 0) growth_final_did[growth_final_did$dose_mg == 120, , drop = FALSE] else data.frame()
  tex <- c(
    "\\paragraph{Gemcitabine-associated pseudotime remodeling and in vivo response.}",
    paste0(
      "Using sample-level ECDFs and an equal-sample ploidy-matched untreated reference, gemcitabine-treated tumors showed a redistribution of pseudotime states among cell-cycle-associated tumor cells relative to untreated controls ",
      "(0 vs treated ECDF RMSE permutation $p=", format_p(get_dose_p(cc_dose, "0_vs_30plus120", FALSE)), "$). ",
      "This result remained similar under EndTimePoint-ploidy-stratified sample-label permutation ($p=", format_p(get_dose_p(cc_dose, "0_vs_30plus120", TRUE)), "$). ",
      "The 30 mg/kg and 120 mg/kg treated groups were not clearly separated from each other ($p=", format_p(get_dose_p(cc_dose, "30_vs_120", FALSE)), "$), so the effect is best described as treatment-associated rather than dose-dependent."
    ),
    "",
    paste0(
      "The corresponding non-cell-cycle tumor-cell compartment did not show the same primary treated-vs-control ECDF result ",
      "(0 vs treated $p=", format_p(get_dose_p(nc_dose, "0_vs_30plus120", FALSE)), "$). ",
      "Paired compartment tests and correlation-difference permutations were used as formal sensitivity analyses rather than relying on a significant/non-significant contrast."
    ),
    "",
    paste0(
      "EndTimePoint-ploidy effect modification was then tested on the pseudotime treatment effect itself, using the difference-in-differences curve ",
      "$[\\Delta ECDF_{ETP\\text{-}higher}(t)-\\Delta ECDF_{ETP\\text{-}lower}(t)]$. ",
      "For CellCycle, the ETP-higher-minus-ETP-lower DID ECDF RMSE was ", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse")), " ",
      "(ploidy-label permutation within treatment groups $p=", format_p(first_value(cc_pseudotime_ploidy_did, "p_did_ecdf_rmse")), "$; bootstrap 95\\% CI [", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse_ci025")), ", ", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse_ci975")), "]). ",
      "For NonCellCycle, the analogous DID ECDF RMSE was ", format_num(first_value(nc_pseudotime_ploidy_did, "did_ecdf_rmse")), " ",
      "($p=", format_p(first_value(nc_pseudotime_ploidy_did, "p_did_ecdf_rmse")), "$). ",
      "This analysis belongs to the pseudotime branch: it asks whether gemcitabine-induced pseudotime redistribution differs between ETP-lower and ETP-higher sample groups, not whether TGI differs by ploidy."
    ),
    "",
    paste0(
      "Among treated tumors, the magnitude of the CellCycle pseudotime redistribution showed a positive association with AUC-based TGI ",
      "(Pearson $r=", format_num(cc_primary$pearson_r[1]), "$, exact/permutation $p=", format_p(cc_primary$pearson_p_permutation_two_sided[1]), "$; ",
      "Spearman $\\rho=", format_num(cc_primary$spearman_rho[1]), "$, exact/permutation $p=", format_p(cc_primary$spearman_p_permutation_two_sided[1]), "$). ",
      "After exploratory multiplicity correction and robustness checks, this should be interpreted as a positive nominal or suggestive association, not as a definitive predictive biomarker."
    ),
    "",
    paste0(
      "To address overall tumor-growth inhibition directly, AUC-based TGI was modeled as a sample-level endpoint in treated tumors. ",
      "In CellCycle, the dose-adjusted model `TGI_percent_auc ~ ECDF shift + EndTimePoint ploidy + dose' gave an ECDF-shift term $p=", format_p(first_value(cc_auc_adjusted_shift, "p_value")), "$ and an EndTimePoint-ploidy term $p=", format_p(first_value(cc_auc_adjusted_ploidy, "p_value")), "$. ",
      "The shift-by-ploidy interaction sensitivity model gave $p=", format_p(first_value(cc_auc_shift_ploidy_interaction, "p_value")), "$. ",
      "In NonCellCycle, the analogous ECDF-shift term gave $p=", format_p(first_value(nc_auc_adjusted_shift, "p_value")), "$ and the shift-by-ploidy interaction gave $p=", format_p(first_value(nc_auc_shift_ploidy_interaction, "p_value")), "$. ",
      "Thus AUC-TGI is used as the primary summary of overall inhibition, while ploidy-specific equality is not inferred from these small treated-only regressions."
    ),
    "",
    paste0(
      "Exploratory dose-specific response analyses compared 30 and 120 mg/kg tumors for AUC-based TGI and every available endpoint TGI column. ",
      "For AUC-based TGI, the 120-minus-30 mg/kg mean difference was ", format_num(first_value(cc_auc_dose, "mean_diff_b_minus_a")), " percentage points ",
      "(exact/permutation two-sided $p=", format_p(first_value(cc_auc_dose, "permutation_p_two_sided")), "$). ",
      "Within-dose CellCycle shift--TGI correlations were positive in both strata ",
      "(30 mg/kg Pearson $r=", format_num(first_value(cc_30_cor, "pearson_r")), "$, exact/permutation $p=", format_p(first_value(cc_30_cor, "pearson_p_permutation_two_sided")), "$; ",
      "120 mg/kg Pearson $r=", format_num(first_value(cc_120_cor, "pearson_r")), "$, exact/permutation $p=", format_p(first_value(cc_120_cor, "pearson_p_permutation_two_sided")), "$). ",
      "In NonCellCycle cells, the same sample-level AUC TGI dose comparison gave a 120-minus-30 mg/kg mean difference of ", format_num(first_value(nc_auc_dose, "mean_diff_b_minus_a")), " percentage points ",
      "(exact/permutation two-sided $p=", format_p(first_value(nc_auc_dose, "permutation_p_two_sided")), "$), while within-dose NonCellCycle shift--TGI correlations did not show the same positive direction ",
      "(30 mg/kg Pearson $r=", format_num(first_value(nc_30_cor, "pearson_r")), "$, $p=", format_p(first_value(nc_30_cor, "pearson_p_permutation_two_sided")), "$; ",
      "120 mg/kg Pearson $r=", format_num(first_value(nc_120_cor, "pearson_r")), "$, $p=", format_p(first_value(nc_120_cor, "pearson_p_permutation_two_sided")), "$). ",
      "The treated-only interaction sensitivity model for the AUC TGI--shift slope gave $p=", format_p(first_value(int_hit, "p_value")), "$. ",
      "Because each treated dose stratum contains only four samples, these dose-specific results should be framed as descriptive evidence, including negative evidence for a clear 30 vs 120 difference."
    ),
    "",
    paste0(
      "EndTimePoint-ploidy response analyses did not detect a clear ETP-lower versus ETP-higher response difference. ",
      "Within CellCycle, 0 versus treated ECDF tests were computed separately in ETP-lower and ETP-higher samples ",
      "(ETP-lower $p=", format_p(first_value(cc_ecdf_2n, "p_ecdf_rmse")), "$; ETP-higher $p=", format_p(first_value(cc_ecdf_4n, "p_ecdf_rmse")), "$). ",
      "Among treated tumors, the ETP-higher-minus-ETP-lower CellCycle ECDF RMSE difference was ", format_num(first_value(cc_shift_ploidy, "mean_diff_b_minus_a")), " ",
      "(permutation $p=", format_p(first_value(cc_shift_ploidy, "permutation_p_two_sided")), "$; bootstrap 95\\% CI [", format_num(first_value(cc_shift_ploidy, "bootstrap_ci025")), ", ", format_num(first_value(cc_shift_ploidy, "bootstrap_ci975")), "]). ",
      "The corresponding ETP-higher-minus-ETP-lower AUC TGI difference was ", format_num(first_value(tgi_auc_ploidy, "mean_diff_b_minus_a")), " percentage points ",
      "(permutation $p=", format_p(first_value(tgi_auc_ploidy, "permutation_p_two_sided")), "$; bootstrap 95\\% CI [", format_num(first_value(tgi_auc_ploidy, "bootstrap_ci025")), ", ", format_num(first_value(tgi_auc_ploidy, "bootstrap_ci975")), "]). ",
      "The pre-specified equivalence margins were $\\pm 0.01$ ECDF RMSE and $\\pm 10$ TGI percentage points; the CIs should be used to state whether equivalence is formally established rather than interpreting nonsignificant difference tests as proof of equality."
    ),
    "",
    paste0(
      "For the longitudinal growth-response question, raw tumor-volume measurements were modeled as `log(tumor volume / Day 0 tumor volume)' with a natural-spline mixed model containing sample-level random intercepts. ",
      "The primary likelihood-ratio test compared the full model with `ns(day) \\times dose \\times EndTimePoint ploidy' against a reduced model without that three-way interaction; $p=", format_p(first_value(growth_primary_lrt, "p_value")), "$. ",
      "A linear-time sensitivity model gave $p=", format_p(first_value(growth_linear_lrt, "p_value")), "$. ",
      "At Day ", format_num(growth_final_day, digits = 4), ", the model-estimated ploidy difference-in-differences was ", format_num(first_value(did30, "estimate_log_relative_volume_diff")), " for 30 mg/kg ($p=", format_p(first_value(did30, "p_value")), "$) and ", format_num(first_value(did120, "estimate_log_relative_volume_diff")), " for 120 mg/kg ($p=", format_p(first_value(did120, "p_value")), "$). ",
      "This mixed model describes whether longitudinal tumor-growth response differs by the post-treatment EndTimePoint-ploidy group; it is not a baseline predictive-biomarker analysis."
    ),
    "",
    paste0(
      "The formal CellCycle-minus-NonCellCycle TGI-correlation difference was $\\Delta r=", format_num(delta$observed_delta_r[1]), "$ ",
      "(label-swap permutation two-sided $p=", format_p(delta$p_two_sided[1]), "$), which should be described as directionally informative only if not statistically supported."
    )
  )
  write_text(tex, file.path(dirs$Manuscript, "dose_pseudotime_conclusion_revised.tex"))
  change_log <- c(
    "# dose_pseudotime_conclusion Change Log",
    "",
    "- Replaced pooled-control-cell shift metric with sample-equal ploidy-matched control ECDF as the primary metric.",
    "- Retained pooled-control-cell shift metric as a labeled sensitivity analysis.",
    "- Added exact TGI-label permutation p-values for small-n treated-only TGI associations.",
    "- Added leave-one-out, bootstrap, cell-count downsampling, ploidy/dose residualization, paired compartment tests, composition, neighborhood fallback, pseudobulk/module, and growth-model sensitivity outputs.",
    "- Added dose-specific TGI exploration: 30 vs 120 AUC/endpoint TGI comparisons, within-dose shift-TGI correlations, slope-interaction sensitivity, and dose-specific figures.",
    "- Added explicit AUC-TGI association models for the overall inhibition question: AUC-TGI by pseudotime shift, EndTimePoint ploidy, treated dose, and shift-by-ploidy interaction.",
    "- Added spline mixed-effects growth-curve models for the longitudinal response question: log-relative tumor volume by time, dose, EndTimePoint ploidy, and time-by-dose-by-ploidy interaction.",
    "- Wording was revised to avoid unsupported dose-dependence, global ploidy-independence, predictive-biomarker, or definitive CellCycle-specific TGI-association claims."
  )
  write_text(change_log, file.path(dirs$Manuscript, "dose_pseudotime_conclusion_change_log.md"))
  lit <- c(
    "# Literature and Methodology Note",
    "",
    "A focused literature search was performed during the 2026-07-05 analysis run. The revised analysis follows common single-cell perturbation-analysis principles: use biological samples as the inference unit, avoid cell-level pseudoreplication, report neighborhood or state-abundance changes, and treat small-n response correlations as hypothesis-supporting rather than definitive.",
    "",
    "Key methodological precedents:",
    "",
    "- Milo uses overlapping k-nearest-neighbor graph neighborhoods for single-cell differential abundance testing and is motivated by the limitation that discrete clusters may miss continuous trajectory/state changes: https://www.nature.com/articles/s41587-021-01033-z",
    "- muscat emphasizes multi-sample multi-condition single-cell analysis with sample-level inference and pseudobulk/differential-state strategies rather than treating cells as independent biological replicates: https://www.nature.com/articles/s41467-020-19894-4",
    "- tradeSeq models expression along pseudotime/lineages with generalized additive models, providing the trajectory-DE precedent for pseudotime-resolved gene or module trends: https://www.nature.com/articles/s41467-020-14766-3",
    "- Soneson and Robinson benchmarked single-cell DE methods and highlighted robustness/bias concerns, supporting cautious interpretation and replicate-aware/pseudobulk sensitivity analyses: https://www.nature.com/articles/nmeth.4612",
    "- MELD-style perturbation analysis illustrates continuous perturbation-response modeling across transcriptomic state space rather than only fixed clusters: https://www.nature.com/articles/s41587-020-00803-5",
    "",
    "How those ideas map to this output:",
    "",
    "- Sample-aware ECDF tests and sample-level label permutations are used for pseudotime distribution shifts.",
    "- Pseudobulk-style sample-level expression summaries are used when targeted pseudobulk count tables are available.",
    "- MiloR was unavailable, so pseudotime-bin sample-level abundance is used as a clearly labeled fallback neighborhood analysis.",
    "- Pseudotime-bin module trends are used as a fallback for tradeSeq/GAM-style trajectory-expression analysis because raw cell-level expression is not provided through `--input_root`.",
    "- Leave-one-out, bootstrap, exact label permutation, and a negative/comparison compartment are used to avoid overstating small-n TGI associations."
  )
  write_text(lit, file.path(dirs$Manuscript, "literature_methodology_note.md"))
  report <- c(
    "# Final Analysis Report",
    "",
    "## Executive summary",
    "",
    paste0("CellCycle 0 vs treated sample-aware ECDF p = ", format_p(get_dose_p(cc_dose, "0_vs_30plus120", FALSE)), "; ploidy-stratified p = ", format_p(get_dose_p(cc_dose, "0_vs_30plus120", TRUE)), "."),
    paste0("NonCellCycle 0 vs treated sample-aware ECDF p = ", format_p(get_dose_p(nc_dose, "0_vs_30plus120", FALSE)), "."),
    paste0("Pseudotime A2 CellCycle ETP-higher-minus-ETP-lower treatment-effect DID ECDF RMSE = ", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse")), ", permutation p = ", format_p(first_value(cc_pseudotime_ploidy_did, "p_did_ecdf_rmse")), "."),
    paste0("CellCycle treated-only AUC TGI association: Pearson r = ", format_num(cc_primary$pearson_r[1]), ", permutation p = ", format_p(cc_primary$pearson_p_permutation_two_sided[1]), "; Spearman rho = ", format_num(cc_primary$spearman_rho[1]), ", permutation p = ", format_p(cc_primary$spearman_p_permutation_two_sided[1]), "."),
    paste0("NonCellCycle treated-only AUC TGI association: Pearson r = ", format_num(nc_primary$pearson_r[1]), ", permutation p = ", format_p(nc_primary$pearson_p_permutation_two_sided[1]), "; Spearman rho = ", format_num(nc_primary$spearman_rho[1]), ", permutation p = ", format_p(nc_primary$spearman_p_permutation_two_sided[1]), "."),
    paste0("Question 1 AUC-TGI model: CellCycle dose-adjusted ECDF-shift term p = ", format_p(first_value(cc_auc_adjusted_shift, "p_value")), "; EndTimePoint-ploidy term p = ", format_p(first_value(cc_auc_adjusted_ploidy, "p_value")), "; shift-by-ploidy interaction p = ", format_p(first_value(cc_auc_shift_ploidy_interaction, "p_value")), "."),
    paste0("Dose-specific CellCycle AUC TGI comparison, 120 minus 30 mg/kg mean difference = ", format_num(first_value(cc_auc_dose, "mean_diff_b_minus_a")), ", permutation p = ", format_p(first_value(cc_auc_dose, "permutation_p_two_sided")), "."),
    paste0("Dose-specific NonCellCycle within-dose shift-TGI correlations do not mirror the CellCycle positive pattern: 30 mg/kg Pearson r = ", format_num(first_value(nc_30_cor, "pearson_r")), ", permutation p = ", format_p(first_value(nc_30_cor, "pearson_p_permutation_two_sided")), "; 120 mg/kg Pearson r = ", format_num(first_value(nc_120_cor, "pearson_r")), ", permutation p = ", format_p(first_value(nc_120_cor, "pearson_p_permutation_two_sided")), "."),
    paste0("EndTimePoint-ploidy response analysis: CellCycle treated ECDF RMSE ETP-higher minus ETP-lower mean difference = ", format_num(first_value(cc_shift_ploidy, "mean_diff_b_minus_a")), ", permutation p = ", format_p(first_value(cc_shift_ploidy, "permutation_p_two_sided")), "; AUC TGI ETP-higher minus ETP-lower mean difference = ", format_num(first_value(tgi_auc_ploidy, "mean_diff_b_minus_a")), ", permutation p = ", format_p(first_value(tgi_auc_ploidy, "permutation_p_two_sided")), "."),
    paste0("Question 2 spline mixed model: time-by-dose-by-ploidy LRT p = ", format_p(first_value(growth_primary_lrt, "p_value")), "; Day ", format_num(growth_final_day, digits = 4), " ploidy difference-in-differences p values are 30 mg/kg p = ", format_p(first_value(did30, "p_value")), " and 120 mg/kg p = ", format_p(first_value(did120, "p_value")), "."),
    "",
    "## What changed relative to original analysis",
    "",
    "The primary shift metric now compares each sample with an equal-weighted average of ploidy-matched untreated sample ECDFs. The older pooled-cell control reference is retained only as a sensitivity output.",
    "",
    "## Analysis structure",
    "",
    "A. Pseudotime analysis: A1 tests whether gemcitabine changes CellCycle/NonCellCycle pseudotime distributions; A2 tests whether the pseudotime treatment effect differs between ETP-lower and ETP-higher sample groups.",
    "",
    "B. TGI / tumor-growth response analysis: B1 tests whether pseudotime shift is associated with AUC-based TGI; B2 tests whether tumor-growth response differs between ETP-lower and ETP-higher, using AUC-TGI as the overall inhibition summary and a spline mixed model for longitudinal growth curves.",
    "",
    "## A1: CellCycle pseudotime treatment effect",
    "",
    paste0("The CellCycle treated-vs-untreated ECDF test gives p = ", format_p(get_dose_p(cc_dose, "0_vs_30plus120", FALSE)), "; the EndTimePoint-ploidy-stratified version gives p = ", format_p(get_dose_p(cc_dose, "0_vs_30plus120", TRUE)), "."),
    "",
    "## A1: NonCellCycle pseudotime comparison",
    "",
    paste0("The NonCellCycle treated-vs-untreated ECDF test gives p = ", format_p(get_dose_p(nc_dose, "0_vs_30plus120", FALSE)), ". This is the comparison compartment for specificity. Its treated-only AUC TGI association is not positive like CellCycle: Pearson r = ", format_num(nc_primary$pearson_r[1]), ", permutation p = ", format_p(nc_primary$pearson_p_permutation_two_sided[1]), "."),
    "",
    "## A2: pseudotime ploidy effect modification",
    "",
    paste0("A2 compares the treatment-effect curves themselves: Delta ECDF_ETP-lower(t) = treated ETP-lower minus control ETP-lower, Delta ECDF_ETP-higher(t) = treated ETP-higher minus control ETP-higher, and DID(t) = Delta ECDF_ETP-higher(t) minus Delta ECDF_ETP-lower(t). For CellCycle, DID ECDF RMSE = ", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse")), ", ETP-label permutation p = ", format_p(first_value(cc_pseudotime_ploidy_did, "p_did_ecdf_rmse")), ", bootstrap CI = [", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse_ci025")), ", ", format_num(first_value(cc_pseudotime_ploidy_did, "did_ecdf_rmse_ci975")), "]. For NonCellCycle, DID ECDF RMSE = ", format_num(first_value(nc_pseudotime_ploidy_did, "did_ecdf_rmse")), ", permutation p = ", format_p(first_value(nc_pseudotime_ploidy_did, "p_did_ecdf_rmse")), "."),
    "This A2 test is a pseudotime analysis, not a TGI analysis. It asks whether gemcitabine-induced pseudotime redistribution differs between ETP-lower and ETP-higher sample groups.",
    "",
    "## B1: TGI correlation and robustness",
    "",
    paste0("The primary CellCycle TGI association is positive but should remain nominal/suggestive. Leave-one-out all positive: ", bool_text(robustness$summary$leave_one_out_all_positive[1]), "; bootstrap Pearson r 95% interval = [", format_num(robustness$summary$bootstrap_pearson_r_ci025[1]), ", ", format_num(robustness$summary$bootstrap_pearson_r_ci975[1]), "]."),
    "",
    "## B1: AUC-TGI correlation with pseudotime shift and ploidy",
    "",
    paste0("AUC-based TGI is the primary summary endpoint for overall tumor-growth inhibition. In the treated-only CellCycle model adjusted for dose and EndTimePoint ploidy, the ECDF-shift term has p = ", format_p(first_value(cc_auc_adjusted_shift, "p_value")), "; the EndTimePoint-ploidy term has p = ", format_p(first_value(cc_auc_adjusted_ploidy, "p_value")), "; and the shift-by-ploidy interaction has p = ", format_p(first_value(cc_auc_shift_ploidy_interaction, "p_value")), ". The analogous NonCellCycle adjusted ECDF-shift term has p = ", format_p(first_value(nc_auc_adjusted_shift, "p_value")), ", and the NonCellCycle shift-by-ploidy interaction has p = ", format_p(first_value(nc_auc_shift_ploidy_interaction, "p_value")), "."),
    "",
    "## B2a support: dose-specific and endpoint TGI",
    "",
    paste0("The new dose-specific module asks two questions: whether TGI differs between 30 and 120 mg/kg, and whether the shift-TGI association has the same direction within each dose stratum. For CellCycle AUC-based TGI, 120 minus 30 mg/kg mean difference = ", format_num(first_value(cc_auc_dose, "mean_diff_b_minus_a")), " percentage points; exact/permutation p = ", format_p(first_value(cc_auc_dose, "permutation_p_two_sided")), ". Across endpoint TGI columns, the smallest exploratory endpoint permutation p = ", format_p(endpoint_min_p), " and the smallest BH q = ", format_p(endpoint_min_q), "."),
    paste0("Within-dose CellCycle AUC shift-TGI correlations: 30 mg/kg Pearson r = ", format_num(first_value(cc_30_cor, "pearson_r")), ", permutation p = ", format_p(first_value(cc_30_cor, "pearson_p_permutation_two_sided")), "; 120 mg/kg Pearson r = ", format_num(first_value(cc_120_cor, "pearson_r")), ", permutation p = ", format_p(first_value(cc_120_cor, "pearson_p_permutation_two_sided")), ". The interaction-model slope-difference p = ", format_p(first_value(int_hit, "p_value")), ". Because n = 4 per treated dose, this is descriptive evidence and should be reported as exploratory."),
    paste0("For NonCellCycle, the sample-level AUC TGI 120 minus 30 mg/kg mean difference is the same response contrast, ", format_num(first_value(nc_auc_dose, "mean_diff_b_minus_a")), " percentage points, with exact/permutation p = ", format_p(first_value(nc_auc_dose, "permutation_p_two_sided")), ". Across endpoint TGI columns, the smallest exploratory endpoint permutation p = ", format_p(nc_endpoint_min_p), " and the smallest BH q = ", format_p(nc_endpoint_min_q), ". Within-dose NonCellCycle AUC shift-TGI correlations are negative/weak: 30 mg/kg Pearson r = ", format_num(first_value(nc_30_cor, "pearson_r")), ", permutation p = ", format_p(first_value(nc_30_cor, "pearson_p_permutation_two_sided")), "; 120 mg/kg Pearson r = ", format_num(first_value(nc_120_cor, "pearson_r")), ", permutation p = ", format_p(first_value(nc_120_cor, "pearson_p_permutation_two_sided")), ". The NonCellCycle interaction-model slope-difference p = ", format_p(first_value(nc_int_hit, "p_value")), "."),
    "",
    "## B2a: AUC-TGI response by EndTimePoint ploidy and supporting shift summaries",
    "",
    paste0("This module treats EndTimePoint ploidy as an effect modifier. Within-group CellCycle 0 vs treated ECDF tests gave ETP-lower p = ", format_p(first_value(cc_ecdf_2n, "p_ecdf_rmse")), " and ETP-higher p = ", format_p(first_value(cc_ecdf_4n, "p_ecdf_rmse")), ". Within-group NonCellCycle tests gave ETP-lower p = ", format_p(first_value(nc_ecdf_2n, "p_ecdf_rmse")), " and ETP-higher p = ", format_p(first_value(nc_ecdf_4n, "p_ecdf_rmse")), "."),
    paste0("Among treated tumors, CellCycle ECDF RMSE ETP-higher minus ETP-lower mean difference = ", format_num(first_value(cc_shift_ploidy, "mean_diff_b_minus_a")), "; permutation p = ", format_p(first_value(cc_shift_ploidy, "permutation_p_two_sided")), "; bootstrap 95% CI = [", format_num(first_value(cc_shift_ploidy, "bootstrap_ci025")), ", ", format_num(first_value(cc_shift_ploidy, "bootstrap_ci975")), "]. The pre-specified ECDF RMSE equivalence margin was +/-0.01; CI within margin = ", bool_text(first_value(cc_shift_ploidy, "ci_within_equivalence_margin")), "."),
    paste0("TGI AUC ETP-higher minus ETP-lower mean difference = ", format_num(first_value(tgi_auc_ploidy, "mean_diff_b_minus_a")), " percentage points; permutation p = ", format_p(first_value(tgi_auc_ploidy, "permutation_p_two_sided")), "; bootstrap 95% CI = [", format_num(first_value(tgi_auc_ploidy, "bootstrap_ci025")), ", ", format_num(first_value(tgi_auc_ploidy, "bootstrap_ci975")), "]. The pre-specified AUC TGI equivalence margin was +/-10 percentage points; CI within margin = ", bool_text(first_value(tgi_auc_ploidy, "ci_within_equivalence_margin")), "."),
    paste0("Interaction sensitivities: CellCycle treatment-by-ploidy ECDF RMSE interaction p = ", format_p(first_value(cc_shift_interaction, "p_value")), "; treated TGI ploidy term p = ", format_p(first_value(cc_tgi_ploidy, "p_value")), "; treated shift-by-ploidy TGI interaction p = ", format_p(first_value(cc_shift_tgi_ploidy, "p_value")), ". These small-n models should be interpreted as screening sensitivity rather than definitive proof of equality or difference."),
    "",
    "## B2b: longitudinal growth response by EndTimePoint ploidy",
    "",
    paste0("The primary spline mixed-effects model uses `log(tumor volume / Day 0 tumor volume)` and a sample-level random intercept. The full-vs-reduced likelihood-ratio test for the time-by-dose-by-EndTimePoint-ploidy interaction gives p = ", format_p(first_value(growth_primary_lrt, "p_value")), ". The linear-time sensitivity version gives p = ", format_p(first_value(growth_linear_lrt, "p_value")), "."),
    paste0("At Day ", format_num(growth_final_day, digits = 4), ", the ETP-higher-minus-ETP-lower treatment-effect difference on the log-relative-volume scale is ", format_num(first_value(did30, "estimate_log_relative_volume_diff")), " for 30 mg/kg (p = ", format_p(first_value(did30, "p_value")), ") and ", format_num(first_value(did120, "estimate_log_relative_volume_diff")), " for 120 mg/kg (p = ", format_p(first_value(did120, "p_value")), "). These contrasts compare `(ETP-higher treated - ETP-higher control) - (ETP-lower treated - ETP-lower control)` and directly address whether longitudinal treatment response differs by EndTimePoint ploidy."),
    "Because ETP is determined after treatment, all ETP-stratified growth results are endpoint-defined associations and not baseline predictive effects.",
    "",
    "## Ploidy, dose, and cell-count sensitivity",
    "",
    "Ploidy/dose confounding tables and downsampling summaries are in `C_SupportingAnalyses/C1_Robustness/`. These should be used to say whether the shift is not detectably explained by ploidy or dose; they should not be used to claim strong independence.",
    "",
    "## CellCycle vs NonCellCycle comparison",
    "",
    paste0("The CellCycle-minus-NonCellCycle TGI-correlation delta r is ", format_num(delta$observed_delta_r[1]), " with label-swap two-sided p = ", format_p(delta$p_two_sided[1]), ". If this is not significant, the CellCycle-specific TGI association remains underpowered."),
    "",
    "## Composition and cluster-specific results",
    "",
    "Composition tables test CellCycle fraction and cluster fractions at the sample level. Cluster-specific pseudotime and TGI association outputs identify possible cluster drivers without treating cells as independent samples.",
    "",
    "## Gene program, pseudobulk, neighborhood, and growth-model results",
    "",
    "When available, targeted pseudobulk counts were used for module/pseudobulk summaries. MiloR was unavailable, so neighborhood DA used pseudotime bins as a fallback. The legacy GrowthModel module remains a supportive sensitivity analysis, while the GrowthCurveMixedModel module is the formal longitudinal growth-response analysis for EndTimePoint-ploidy effect modification.",
    "",
    "## Final answer to the scientific question",
    "",
    "The analysis now has two main branches. In the pseudotime branch, A1 tests gemcitabine-associated pseudotime redistribution, and A2 tests whether that effect differs between ETP-lower and ETP-higher samples. In the TGI branch, B1 uses AUC-based TGI for shift-TGI association, while B2 separates overall AUC-TGI response from longitudinal mixed-model growth response. ETP is post-treatment, so these are endpoint-defined associations rather than baseline predictive effects.",
    "",
    "## Supported claims",
    "",
    "- Gemcitabine treatment is associated with a CellCycle pseudotime redistribution.",
    "- The CellCycle treatment effect can be tested with sample-level and ploidy-stratified permutations.",
    "- Pseudotime EndTimePoint-ploidy effect modification is tested directly with Delta ECDF ETP-higher-minus-ETP-lower difference-in-differences curves.",
    "- The treated-only CellCycle TGI association is positive in direction.",
    "- AUC-based TGI is the primary endpoint for the overall tumor-growth inhibition question.",
    "- The spline mixed-effects model is the primary analysis for whether longitudinal tumor-growth response differs by EndTimePoint ploidy.",
    "- Dose-specific TGI checks can be reported as exploratory negative/descriptive evidence.",
    "- EndTimePoint-ploidy response checks compare ETP-lower with ETP-higher, with equivalence judged by the pre-specified CI margins.",
    "",
    "## Suggestive or unsupported claims",
    "",
    "- Predictive biomarker claims are unsupported.",
    "- Dose dependence is unsupported unless 30 vs 120 becomes significant.",
    "- Ploidy-response equality should not be claimed unless equivalence margins are satisfied.",
    "- Longitudinal ploidy-response equality should not be claimed from nonsignificant mixed-model tests without a separate equivalence framework.",
    "- Global ploidy independence should not be claimed; use `not detectably explained by ploidy` if robustness tables support it.",
    "- A statistically proven CellCycle-specific TGI association requires a significant delta-r/interaction test.",
    "",
    "## Recommended manuscript panels",
    "",
    "- CellCycle sample-equal ECDFs by dose/ploidy.",
    "- NonCellCycle sample-equal ECDFs by dose/ploidy.",
    "- CellCycle AUC TGI vs ECDF RMSE.",
    "- Pseudotime ETP-lower vs ETP-higher Delta ECDF effect-modification curves.",
    "- AUC-TGI association panels by EndTimePoint ploidy and dose.",
    "- Dose-specific AUC/endpoint TGI panels and within-dose shift-TGI panel.",
    "- Ploidy-response shift and TGI panels.",
    "- Spline mixed-model predicted growth curves and ploidy difference-in-differences over time.",
    "- Leave-one-out/bootstrap robustness panels.",
    "- Paired CellCycle vs NonCellCycle ECDF RMSE panel.",
    "- Composition and cluster-specific CellCycle panels.",
    "- Legacy growth curves by dose and high/low CellCycle shift."
  )
  write_text(report, file.path(dirs$Manuscript, "final_analysis_report.md"))
  list(cc_primary = cc_primary, nc_primary = nc_primary)
}

sample_exclusion_log <- data.frame(analysis = character(), sample_id = character(), reason = character(), stringsAsFactors = FALSE)

message("Running 04h pseudotime/TGI ", analysis_version)
message("Input root: ", input_root)
message("Results root: ", results_root)
message("Fixed sample-level ETP threshold: ", format_num(endpoint_threshold, 8))

results <- list()
for (name in names(compartments)) {
  message("Compartment: ", name)
  results[[name]] <- run_compartment(name, compartments[[name]], n_perm)
}

message("Running EndTimePoint ploidy continuous-variable exploration")
for (name in names(compartments)) {
  plot_endpoint_continuous_exploration(name, results[[name]])
  plot_sample_level_etp_ploidy_violin(name, results[[name]], compartments[[name]], n_perm)
}

message("Running direct 9-panel sample-level ECDF comparisons")
direct_ecdf <- lapply(names(compartments), function(name) {
  plot_direct_group_ecdf_comparisons(name, results[[name]], compartments[[name]], n_perm)
})
names(direct_ecdf) <- names(compartments)

message("Running robustness analyses")
loo <- leave_one_out_primary(results$CellCycle$equal_shift)
boot <- bootstrap_primary(results$CellCycle$equal_shift, n_boot)
rob_summary <- summarize_bootstrap(boot, loo, results$CellCycle$tgi_stats)
write_csv(loo, file.path(dirs$Robustness, "CellCycle_primary_TGI_leave_one_out.csv"))
write_csv(boot, file.path(dirs$Robustness, "CellCycle_primary_TGI_bootstrap.csv"))
write_csv(rob_summary, file.path(dirs$Robustness, "CellCycle_primary_TGI_robustness_summary.csv"))
plot_robustness(loo, boot)

message("Running cell-count downsampling sensitivity")
down_cc <- downsample_sensitivity("CellCycle", results$CellCycle$df, results$CellCycle$sample_meta, n_downsample, min(n_perm, 1000L))
down_nc <- downsample_sensitivity("NonCellCycle", results$NonCellCycle$df, results$NonCellCycle$sample_meta, n_downsample, min(n_perm, 1000L))
write_csv(down_cc$summary, file.path(dirs$Robustness, "downsampling_CellCycle_summary.csv"))
write_csv(down_cc$iterations, file.path(dirs$Robustness, "downsampling_CellCycle_iterations.csv"))
write_csv(down_nc$summary, file.path(dirs$Robustness, "downsampling_NonCellCycle_summary.csv"))
write_csv(down_nc$iterations, file.path(dirs$Robustness, "downsampling_NonCellCycle_iterations.csv"))
sample_exclusion_log <- rbind(sample_exclusion_log, down_cc$exclusions, down_nc$exclusions)

message("Running ploidy/dose confounding checks")
confounding <- run_confounding(results)

message("Running dose-specific TGI exploration")
dose_specific <- run_dose_specific_tgi_exploration(results, n_perm)

message("Running pseudotime ploidy effect-modification analyses")
pseudotime_ploidy <- run_pseudotime_ploidy_effect_modification(results, n_perm, n_boot)

message("Running AUC-TGI association analyses")
auc_tgi <- run_auc_tgi_association_analysis(results)

message("Running ploidy response analyses")
ploidy_response <- run_ploidy_response_analysis(results, n_perm, n_boot)

message("Running CellCycle vs NonCellCycle comparison")
comparison <- compartment_comparison(results)

message("Running composition analyses")
composition <- composition_analysis(results, n_perm)

message("Running growth model analyses")
growth <- growth_model_analysis(results)

message("Running growth-curve mixed model analyses")
growth_curve_mixed <- run_growth_curve_mixed_model(results)

message("Running neighborhood DA fallback")
neighborhood <- neighborhood_fallback(results)

message("Running optional gene-program and pseudobulk modules")
optional <- run_optional_expression_modules(results)

message("Writing manuscript outputs")
manuscript <- write_manuscript_outputs(
  results,
  robustness = list(loo = loo, boot = boot, summary = rob_summary),
  comparison = comparison,
  composition = composition,
  growth = growth,
  dose_specific = dose_specific,
  pseudotime_ploidy = pseudotime_ploidy,
  auc_tgi = auc_tgi,
  ploidy_response = ploidy_response,
  growth_curve_mixed = growth_curve_mixed
)

sample_exclusion_log <- sample_exclusion_log[!duplicated(sample_exclusion_log), , drop = FALSE]
write_csv(sample_exclusion_log, file.path(dirs$logs, "sample_exclusion_log.csv"))

scalar_or_na <- function(x) {
  if (length(x) > 0) x[1] else NA_real_
}

summary_growth_final_day <- if (nrow(growth_curve_mixed$contrasts) > 0) max(growth_curve_mixed$contrasts$day, na.rm = TRUE) else NA_real_
summary_growth_final_did <- if (nrow(growth_curve_mixed$contrasts) > 0 && is.finite(summary_growth_final_day)) {
  growth_curve_mixed$contrasts[growth_curve_mixed$contrasts$contrast_type == "ploidy_difference_in_treatment_effect" & growth_curve_mixed$contrasts$day == summary_growth_final_day, , drop = FALSE]
} else data.frame()

summary_values <- list(
  results_root = results_root,
  cellcycle_treated_vs_control_p = results$CellCycle$dose_tests$p_ecdf_rmse[results$CellCycle$dose_tests$comparison == "0_vs_30plus120" & !results$CellCycle$dose_tests$stratified_by_end_timepoint_ploidy_group][1],
  cellcycle_treated_vs_control_ploidy_stratified_p = results$CellCycle$dose_tests$p_ecdf_rmse[results$CellCycle$dose_tests$comparison == "0_vs_30plus120" & results$CellCycle$dose_tests$stratified_by_end_timepoint_ploidy_group][1],
  noncellcycle_treated_vs_control_p = results$NonCellCycle$dose_tests$p_ecdf_rmse[results$NonCellCycle$dose_tests$comparison == "0_vs_30plus120" & !results$NonCellCycle$dose_tests$stratified_by_end_timepoint_ploidy_group][1],
  cellcycle_pseudotime_ploidy_did_ecdf_rmse = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse[pseudotime_ploidy$tests$compartment == "CellCycle"]),
  cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p = scalar_or_na(pseudotime_ploidy$tests$p_did_ecdf_rmse[pseudotime_ploidy$tests$compartment == "CellCycle"]),
  cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci025 = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse_ci025[pseudotime_ploidy$tests$compartment == "CellCycle"]),
  cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci975 = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse_ci975[pseudotime_ploidy$tests$compartment == "CellCycle"]),
  noncellcycle_pseudotime_ploidy_did_ecdf_rmse = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse[pseudotime_ploidy$tests$compartment == "NonCellCycle"]),
  noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p = scalar_or_na(pseudotime_ploidy$tests$p_did_ecdf_rmse[pseudotime_ploidy$tests$compartment == "NonCellCycle"]),
  noncellcycle_pseudotime_ploidy_did_ecdf_rmse_ci025 = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse_ci025[pseudotime_ploidy$tests$compartment == "NonCellCycle"]),
  noncellcycle_pseudotime_ploidy_did_ecdf_rmse_ci975 = scalar_or_na(pseudotime_ploidy$tests$did_ecdf_rmse_ci975[pseudotime_ploidy$tests$compartment == "NonCellCycle"]),
  cellcycle_TGI_AUC_pearson_r = manuscript$cc_primary$pearson_r[1],
  cellcycle_TGI_AUC_pearson_perm_p = manuscript$cc_primary$pearson_p_permutation_two_sided[1],
  cellcycle_AUC_TGI_adjusted_shift_p = scalar_or_na(auc_tgi$models$p_value[auc_tgi$models$compartment == "CellCycle" & auc_tgi$models$analysis == "shift_ploidy_dose_adjusted" & auc_tgi$models$term == "scale(ecdf_rmse)"]),
  cellcycle_AUC_TGI_adjusted_end_timepoint_ploidy_group_p = scalar_or_na(auc_tgi$models$p_value[auc_tgi$models$compartment == "CellCycle" & auc_tgi$models$analysis == "shift_ploidy_dose_adjusted" & grepl("end_timepoint_ploidy_group", auc_tgi$models$term)]),
  cellcycle_AUC_TGI_shift_by_ploidy_interaction_p = scalar_or_na(auc_tgi$models$p_value[auc_tgi$models$compartment == "CellCycle" & auc_tgi$models$analysis == "shift_by_ploidy_dose_adjusted" & grepl(":", auc_tgi$models$term)]),
  cellcycle_AUC_TGI_mean_ETP_unadjusted_slope = scalar_or_na(auc_tgi$mean_etp_models$estimate[auc_tgi$mean_etp_models$compartment == "CellCycle" & auc_tgi$mean_etp_models$analysis == "mean_etp_only" & auc_tgi$mean_etp_models$term == "scale(sample_mean_endpoint_ploidy)"]),
  cellcycle_AUC_TGI_mean_ETP_unadjusted_p = scalar_or_na(auc_tgi$mean_etp_models$p_value[auc_tgi$mean_etp_models$compartment == "CellCycle" & auc_tgi$mean_etp_models$analysis == "mean_etp_only" & auc_tgi$mean_etp_models$term == "scale(sample_mean_endpoint_ploidy)"]),
  cellcycle_AUC_TGI_mean_ETP_adjusted_p = scalar_or_na(auc_tgi$mean_etp_models$p_value[auc_tgi$mean_etp_models$compartment == "CellCycle" & auc_tgi$mean_etp_models$analysis == "mean_etp_ploidy_dose_adjusted" & auc_tgi$mean_etp_models$term == "scale(sample_mean_endpoint_ploidy)"]),
  cellcycle_AUC_TGI_mean_ETP_by_ploidy_interaction_p = scalar_or_na(auc_tgi$mean_etp_models$p_value[auc_tgi$mean_etp_models$compartment == "CellCycle" & auc_tgi$mean_etp_models$analysis == "mean_etp_by_ploidy_dose_adjusted" & grepl(":", auc_tgi$mean_etp_models$term)]),
  compartment_delta_r = comparison$delta$observed_delta_r[1],
  compartment_delta_r_perm_p = comparison$delta$p_two_sided[1],
  cellcycle_TGI_AUC_120_minus_30_mean_diff = scalar_or_na(dose_specific$comparisons$mean_diff_b_minus_a[dose_specific$comparisons$compartment == "CellCycle" & dose_specific$comparisons$tgi_measure == "TGI_percent_auc"]),
  cellcycle_TGI_AUC_30_vs_120_perm_p = scalar_or_na(dose_specific$comparisons$permutation_p_two_sided[dose_specific$comparisons$compartment == "CellCycle" & dose_specific$comparisons$tgi_measure == "TGI_percent_auc"]),
  cellcycle_TGI_AUC_shift_30mg_pearson_r = scalar_or_na(dose_specific$correlations$pearson_r[dose_specific$correlations$compartment == "CellCycle" & dose_specific$correlations$dose_mg == 30 & dose_specific$correlations$pre_specified_primary]),
  cellcycle_TGI_AUC_shift_30mg_perm_p = scalar_or_na(dose_specific$correlations$pearson_p_permutation_two_sided[dose_specific$correlations$compartment == "CellCycle" & dose_specific$correlations$dose_mg == 30 & dose_specific$correlations$pre_specified_primary]),
  cellcycle_TGI_AUC_shift_120mg_pearson_r = scalar_or_na(dose_specific$correlations$pearson_r[dose_specific$correlations$compartment == "CellCycle" & dose_specific$correlations$dose_mg == 120 & dose_specific$correlations$pre_specified_primary]),
  cellcycle_TGI_AUC_shift_120mg_perm_p = scalar_or_na(dose_specific$correlations$pearson_p_permutation_two_sided[dose_specific$correlations$compartment == "CellCycle" & dose_specific$correlations$dose_mg == 120 & dose_specific$correlations$pre_specified_primary]),
  cellcycle_shift_TGI_AUC_dose_interaction_p = scalar_or_na(dose_specific$interactions$p_value[dose_specific$interactions$compartment == "CellCycle" & grepl("scale\\(ecdf_rmse\\).*:factor\\(dose_mg\\)|factor\\(dose_mg\\).*:scale\\(ecdf_rmse\\)", dose_specific$interactions$term)]),
  noncellcycle_TGI_AUC_pearson_r = manuscript$nc_primary$pearson_r[1],
  noncellcycle_TGI_AUC_pearson_perm_p = manuscript$nc_primary$pearson_p_permutation_two_sided[1],
  noncellcycle_AUC_TGI_adjusted_shift_p = scalar_or_na(auc_tgi$models$p_value[auc_tgi$models$compartment == "NonCellCycle" & auc_tgi$models$analysis == "shift_ploidy_dose_adjusted" & auc_tgi$models$term == "scale(ecdf_rmse)"]),
  noncellcycle_AUC_TGI_shift_by_ploidy_interaction_p = scalar_or_na(auc_tgi$models$p_value[auc_tgi$models$compartment == "NonCellCycle" & auc_tgi$models$analysis == "shift_by_ploidy_dose_adjusted" & grepl(":", auc_tgi$models$term)]),
  noncellcycle_TGI_AUC_120_minus_30_mean_diff = scalar_or_na(dose_specific$comparisons$mean_diff_b_minus_a[dose_specific$comparisons$compartment == "NonCellCycle" & dose_specific$comparisons$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_TGI_AUC_30_vs_120_perm_p = scalar_or_na(dose_specific$comparisons$permutation_p_two_sided[dose_specific$comparisons$compartment == "NonCellCycle" & dose_specific$comparisons$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_TGI_AUC_shift_30mg_pearson_r = scalar_or_na(dose_specific$correlations$pearson_r[dose_specific$correlations$compartment == "NonCellCycle" & dose_specific$correlations$dose_mg == 30 & dose_specific$correlations$shift_metric == "ecdf_rmse" & dose_specific$correlations$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_TGI_AUC_shift_30mg_perm_p = scalar_or_na(dose_specific$correlations$pearson_p_permutation_two_sided[dose_specific$correlations$compartment == "NonCellCycle" & dose_specific$correlations$dose_mg == 30 & dose_specific$correlations$shift_metric == "ecdf_rmse" & dose_specific$correlations$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_TGI_AUC_shift_120mg_pearson_r = scalar_or_na(dose_specific$correlations$pearson_r[dose_specific$correlations$compartment == "NonCellCycle" & dose_specific$correlations$dose_mg == 120 & dose_specific$correlations$shift_metric == "ecdf_rmse" & dose_specific$correlations$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_TGI_AUC_shift_120mg_perm_p = scalar_or_na(dose_specific$correlations$pearson_p_permutation_two_sided[dose_specific$correlations$compartment == "NonCellCycle" & dose_specific$correlations$dose_mg == 120 & dose_specific$correlations$shift_metric == "ecdf_rmse" & dose_specific$correlations$tgi_measure == "TGI_percent_auc"]),
  noncellcycle_shift_TGI_AUC_dose_interaction_p = scalar_or_na(dose_specific$interactions$p_value[dose_specific$interactions$compartment == "NonCellCycle" & grepl("scale\\(ecdf_rmse\\).*:factor\\(dose_mg\\)|factor\\(dose_mg\\).*:scale\\(ecdf_rmse\\)", dose_specific$interactions$term)]),
  cellcycle_ETP_lower_treated_vs_control_ecdf_p = scalar_or_na(ploidy_response$ecdf$p_ecdf_rmse[ploidy_response$ecdf$compartment == "CellCycle" & ploidy_response$ecdf$end_timepoint_ploidy_group == "ETP-lower"]),
  cellcycle_ETP_higher_treated_vs_control_ecdf_p = scalar_or_na(ploidy_response$ecdf$p_ecdf_rmse[ploidy_response$ecdf$compartment == "CellCycle" & ploidy_response$ecdf$end_timepoint_ploidy_group == "ETP-higher"]),
  noncellcycle_ETP_lower_treated_vs_control_ecdf_p = scalar_or_na(ploidy_response$ecdf$p_ecdf_rmse[ploidy_response$ecdf$compartment == "NonCellCycle" & ploidy_response$ecdf$end_timepoint_ploidy_group == "ETP-lower"]),
  noncellcycle_ETP_higher_treated_vs_control_ecdf_p = scalar_or_na(ploidy_response$ecdf$p_ecdf_rmse[ploidy_response$ecdf$compartment == "NonCellCycle" & ploidy_response$ecdf$end_timepoint_ploidy_group == "ETP-higher"]),
  cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_mean_diff = scalar_or_na(ploidy_response$shift$mean_diff_b_minus_a[ploidy_response$shift$compartment == "CellCycle" & ploidy_response$shift$response_metric == "ecdf_rmse"]),
  cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_perm_p = scalar_or_na(ploidy_response$shift$permutation_p_two_sided[ploidy_response$shift$compartment == "CellCycle" & ploidy_response$shift$response_metric == "ecdf_rmse"]),
  cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci025 = scalar_or_na(ploidy_response$shift$bootstrap_ci025[ploidy_response$shift$compartment == "CellCycle" & ploidy_response$shift$response_metric == "ecdf_rmse"]),
  cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci975 = scalar_or_na(ploidy_response$shift$bootstrap_ci975[ploidy_response$shift$compartment == "CellCycle" & ploidy_response$shift$response_metric == "ecdf_rmse"]),
  cellcycle_ecdf_rmse_equivalence_margin_met = scalar_or_na(ploidy_response$shift$ci_within_equivalence_margin[ploidy_response$shift$compartment == "CellCycle" & ploidy_response$shift$response_metric == "ecdf_rmse"]),
  tgi_auc_ETP_higher_minus_ETP_lower_mean_diff = scalar_or_na(ploidy_response$tgi$mean_diff_b_minus_a[ploidy_response$tgi$tgi_measure == "TGI_percent_auc"]),
  tgi_auc_ETP_higher_minus_ETP_lower_perm_p = scalar_or_na(ploidy_response$tgi$permutation_p_two_sided[ploidy_response$tgi$tgi_measure == "TGI_percent_auc"]),
  tgi_auc_ETP_higher_minus_ETP_lower_ci025 = scalar_or_na(ploidy_response$tgi$bootstrap_ci025[ploidy_response$tgi$tgi_measure == "TGI_percent_auc"]),
  tgi_auc_ETP_higher_minus_ETP_lower_ci975 = scalar_or_na(ploidy_response$tgi$bootstrap_ci975[ploidy_response$tgi$tgi_measure == "TGI_percent_auc"]),
  tgi_auc_equivalence_margin_met = scalar_or_na(ploidy_response$tgi$ci_within_equivalence_margin[ploidy_response$tgi$tgi_measure == "TGI_percent_auc"]),
  cellcycle_treatment_by_ploidy_ecdf_interaction_p = scalar_or_na(ploidy_response$models$p_value[ploidy_response$models$compartment == "CellCycle" & ploidy_response$models$analysis == "shift_treatment_by_ploidy_interaction" & grepl(":", ploidy_response$models$term)]),
  cellcycle_TGI_ploidy_adjusted_term_p = scalar_or_na(ploidy_response$models$p_value[ploidy_response$models$compartment == "CellCycle" & ploidy_response$models$analysis == "treated_TGI_by_dose_and_ploidy" & grepl("end_timepoint_ploidy_group", ploidy_response$models$term)]),
  cellcycle_shift_TGI_ploidy_interaction_p = scalar_or_na(ploidy_response$models$p_value[ploidy_response$models$compartment == "CellCycle" & ploidy_response$models$analysis == "treated_TGI_shift_by_ploidy_interaction" & grepl(":", ploidy_response$models$term)]),
  growth_curve_spline_time_dose_ploidy_lrt_p = scalar_or_na(growth_curve_mixed$lrt$p_value[growth_curve_mixed$lrt$analysis == "primary_test_time_dose_ploidy_interaction"]),
  growth_curve_linear_time_dose_ploidy_lrt_p = scalar_or_na(growth_curve_mixed$lrt$p_value[growth_curve_mixed$lrt$analysis == "sensitivity_linear_time_dose_ploidy_interaction"]),
  growth_curve_final_day = summary_growth_final_day,
  growth_curve_day_final_30mg_ploidy_did_estimate = scalar_or_na(summary_growth_final_did$estimate_log_relative_volume_diff[summary_growth_final_did$dose_mg == 30]),
  growth_curve_day_final_30mg_ploidy_did_p = scalar_or_na(summary_growth_final_did$p_value[summary_growth_final_did$dose_mg == 30]),
  growth_curve_day_final_120mg_ploidy_did_estimate = scalar_or_na(summary_growth_final_did$estimate_log_relative_volume_diff[summary_growth_final_did$dose_mg == 120]),
  growth_curve_day_final_120mg_ploidy_did_p = scalar_or_na(summary_growth_final_did$p_value[summary_growth_final_did$dose_mg == 120])
)
write_json_like(summary_values, file.path(results_root, "analysis_summary.json"))
write_text(c(
  paste0("# 04h Pseudotime/TGI ", analysis_version, " Summary"),
  "",
  paste0("- Results root: `", results_root, "`"),
  paste0("- CellCycle 0 vs treated ECDF p: ", format_p(summary_values$cellcycle_treated_vs_control_p)),
  paste0("- CellCycle ploidy-stratified 0 vs treated ECDF p: ", format_p(summary_values$cellcycle_treated_vs_control_ploidy_stratified_p)),
  paste0("- NonCellCycle 0 vs treated ECDF p: ", format_p(summary_values$noncellcycle_treated_vs_control_p)),
  paste0("- A2 CellCycle pseudotime treatment-effect ETP-higher-minus-ETP-lower DID ECDF RMSE: ", format_num(summary_values$cellcycle_pseudotime_ploidy_did_ecdf_rmse), ", permutation p: ", format_p(summary_values$cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p), ", bootstrap CI: [", format_num(summary_values$cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci025), ", ", format_num(summary_values$cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci975), "]"),
  paste0("- A2 NonCellCycle pseudotime treatment-effect ETP-higher-minus-ETP-lower DID ECDF RMSE: ", format_num(summary_values$noncellcycle_pseudotime_ploidy_did_ecdf_rmse), ", permutation p: ", format_p(summary_values$noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p), ", bootstrap CI: [", format_num(summary_values$noncellcycle_pseudotime_ploidy_did_ecdf_rmse_ci025), ", ", format_num(summary_values$noncellcycle_pseudotime_ploidy_did_ecdf_rmse_ci975), "]"),
  paste0("- CellCycle TGI AUC Pearson r: ", format_num(summary_values$cellcycle_TGI_AUC_pearson_r), ", permutation p: ", format_p(summary_values$cellcycle_TGI_AUC_pearson_perm_p)),
  paste0("- Question 1 CellCycle AUC-TGI adjusted ECDF-shift p: ", format_p(summary_values$cellcycle_AUC_TGI_adjusted_shift_p), ", EndTimePoint-ploidy p: ", format_p(summary_values$cellcycle_AUC_TGI_adjusted_end_timepoint_ploidy_group_p), ", shift-by-ploidy interaction p: ", format_p(summary_values$cellcycle_AUC_TGI_shift_by_ploidy_interaction_p)),
  paste0("- CellCycle AUC-TGI vs sample mean ETP: standardized unadjusted slope ", format_num(summary_values$cellcycle_AUC_TGI_mean_ETP_unadjusted_slope), ", p: ", format_p(summary_values$cellcycle_AUC_TGI_mean_ETP_unadjusted_p), "; dose/ploidy-adjusted p: ", format_p(summary_values$cellcycle_AUC_TGI_mean_ETP_adjusted_p), "; mean-ETP-by-ploidy interaction p: ", format_p(summary_values$cellcycle_AUC_TGI_mean_ETP_by_ploidy_interaction_p)),
  paste0("- CellCycle minus NonCellCycle TGI association delta r: ", format_num(summary_values$compartment_delta_r), ", permutation p: ", format_p(summary_values$compartment_delta_r_perm_p)),
  paste0("- Dose-specific CellCycle AUC TGI 120 minus 30 mg/kg mean difference: ", format_num(summary_values$cellcycle_TGI_AUC_120_minus_30_mean_diff), ", permutation p: ", format_p(summary_values$cellcycle_TGI_AUC_30_vs_120_perm_p)),
  paste0("- Within-dose CellCycle AUC shift-TGI Pearson r, 30 mg/kg: ", format_num(summary_values$cellcycle_TGI_AUC_shift_30mg_pearson_r), ", permutation p: ", format_p(summary_values$cellcycle_TGI_AUC_shift_30mg_perm_p)),
  paste0("- Within-dose CellCycle AUC shift-TGI Pearson r, 120 mg/kg: ", format_num(summary_values$cellcycle_TGI_AUC_shift_120mg_pearson_r), ", permutation p: ", format_p(summary_values$cellcycle_TGI_AUC_shift_120mg_perm_p)),
  paste0("- CellCycle shift-TGI dose-interaction sensitivity p: ", format_p(summary_values$cellcycle_shift_TGI_AUC_dose_interaction_p)),
  paste0("- NonCellCycle TGI AUC Pearson r: ", format_num(summary_values$noncellcycle_TGI_AUC_pearson_r), ", permutation p: ", format_p(summary_values$noncellcycle_TGI_AUC_pearson_perm_p)),
  paste0("- Dose-specific NonCellCycle AUC TGI 120 minus 30 mg/kg mean difference: ", format_num(summary_values$noncellcycle_TGI_AUC_120_minus_30_mean_diff), ", permutation p: ", format_p(summary_values$noncellcycle_TGI_AUC_30_vs_120_perm_p)),
  paste0("- Within-dose NonCellCycle AUC shift-TGI Pearson r, 30 mg/kg: ", format_num(summary_values$noncellcycle_TGI_AUC_shift_30mg_pearson_r), ", permutation p: ", format_p(summary_values$noncellcycle_TGI_AUC_shift_30mg_perm_p)),
  paste0("- Within-dose NonCellCycle AUC shift-TGI Pearson r, 120 mg/kg: ", format_num(summary_values$noncellcycle_TGI_AUC_shift_120mg_pearson_r), ", permutation p: ", format_p(summary_values$noncellcycle_TGI_AUC_shift_120mg_perm_p)),
  paste0("- NonCellCycle shift-TGI dose-interaction sensitivity p: ", format_p(summary_values$noncellcycle_shift_TGI_AUC_dose_interaction_p)),
  paste0("- EndTimePoint-ploidy response CellCycle ETP-lower 0 vs treated ECDF p: ", format_p(summary_values$cellcycle_ETP_lower_treated_vs_control_ecdf_p)),
  paste0("- EndTimePoint-ploidy response CellCycle ETP-higher 0 vs treated ECDF p: ", format_p(summary_values$cellcycle_ETP_higher_treated_vs_control_ecdf_p)),
  paste0("- EndTimePoint-ploidy response NonCellCycle ETP-lower 0 vs treated ECDF p: ", format_p(summary_values$noncellcycle_ETP_lower_treated_vs_control_ecdf_p)),
  paste0("- EndTimePoint-ploidy response NonCellCycle ETP-higher 0 vs treated ECDF p: ", format_p(summary_values$noncellcycle_ETP_higher_treated_vs_control_ecdf_p)),
  paste0("- Treated CellCycle ECDF RMSE ETP-higher minus ETP-lower mean difference: ", format_num(summary_values$cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_mean_diff), ", permutation p: ", format_p(summary_values$cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_perm_p), ", bootstrap CI: [", format_num(summary_values$cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci025), ", ", format_num(summary_values$cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci975), "], equivalence margin met: ", bool_text(summary_values$cellcycle_ecdf_rmse_equivalence_margin_met)),
  paste0("- Treated AUC TGI ETP-higher minus ETP-lower mean difference: ", format_num(summary_values$tgi_auc_ETP_higher_minus_ETP_lower_mean_diff), ", permutation p: ", format_p(summary_values$tgi_auc_ETP_higher_minus_ETP_lower_perm_p), ", bootstrap CI: [", format_num(summary_values$tgi_auc_ETP_higher_minus_ETP_lower_ci025), ", ", format_num(summary_values$tgi_auc_ETP_higher_minus_ETP_lower_ci975), "], equivalence margin met: ", bool_text(summary_values$tgi_auc_equivalence_margin_met)),
  paste0("- CellCycle treatment-by-ploidy ECDF interaction p: ", format_p(summary_values$cellcycle_treatment_by_ploidy_ecdf_interaction_p)),
  paste0("- CellCycle TGI ploidy adjusted term p: ", format_p(summary_values$cellcycle_TGI_ploidy_adjusted_term_p)),
  paste0("- CellCycle shift-TGI ploidy interaction p: ", format_p(summary_values$cellcycle_shift_TGI_ploidy_interaction_p)),
  paste0("- Question 2 spline mixed-model time-by-dose-by-ploidy LRT p: ", format_p(summary_values$growth_curve_spline_time_dose_ploidy_lrt_p), "; linear-time sensitivity LRT p: ", format_p(summary_values$growth_curve_linear_time_dose_ploidy_lrt_p)),
  paste0("- Question 2 Day ", format_num(summary_values$growth_curve_final_day, digits = 4), " ploidy difference-in-differences, 30 mg/kg estimate: ", format_num(summary_values$growth_curve_day_final_30mg_ploidy_did_estimate), ", p: ", format_p(summary_values$growth_curve_day_final_30mg_ploidy_did_p), "; 120 mg/kg estimate: ", format_num(summary_values$growth_curve_day_final_120mg_ploidy_did_estimate), ", p: ", format_p(summary_values$growth_curve_day_final_120mg_ploidy_did_p))
), file.path(results_root, "human_readable_summary.md"))

message("Completed 04h pseudotime/TGI ", analysis_version)
message("Key outputs:")
message("  ", file.path(dirs$CellCycle, "stats", "dose_group_ecdf_tests_v4.csv"))
message("  ", file.path(dirs$CellCycle, "stats", "CellCycle_TGI_associations_v4.csv"))
message("  ", file.path(dirs$Manuscript, "final_analysis_report.md"))
