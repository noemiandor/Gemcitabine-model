#!/usr/bin/env Rscript

# Generate interpretation and cross-version comparison reports for the formal
# 04h pseudotime/TGI initial_ploidy, ETP_fixed_threshold_2_25, ETP_boundary_stress_threshold_2_375, and ETP_reference_balanced_threshold_2_24 result trees.
#
# This script intentionally reads the completed result artifacts instead of
# embedding the main statistical estimates.  Re-running it after a formal
# analysis refresh therefore keeps the reports synchronized with the result
# tables and analysis_summary.json files.

suppressPackageStartupMessages({
  library(jsonlite)
})

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
results_root_arg <- cmd_args[["results_root"]]
results_root <- normalizePath(
  if (!is.null(results_root_arg) && nzchar(results_root_arg)) {
    results_root_arg
  } else {
    Sys.getenv(
      "RESULTS_ROOT",
      unset = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04h_pseudotime_TGI"
    )
  },
  mustWork = FALSE
)
roots <- c(
  initial_ploidy = file.path(results_root, "initial_ploidy"),
  ETP_fixed_threshold_2_25 = file.path(results_root, "ETP_fixed_threshold_2_25"),
  ETP_boundary_stress_threshold_2_375 = file.path(results_root, "ETP_boundary_stress_threshold_2_375"),
  ETP_reference_balanced_threshold_2_24 = file.path(results_root, "ETP_reference_balanced_threshold_2_24")
)
comparison_root <- file.path(results_root, "method_comparison")

required_files <- c(
  file.path(roots, "analysis_summary.json"),
  file.path(roots[c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24")], "sample_end_timepoint_ploidy_assignments.csv"),
  file.path(roots[c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24")], "ETP_group_dose_coverage_audit.csv")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Missing required formal result files:\n", paste(missing_files, collapse = "\n"))
}

dir.create(comparison_root, recursive = TRUE, showWarnings = FALSE)
for (version in c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24")) {
  dir.create(file.path(roots[[version]], "Manuscript"), recursive = TRUE, showWarnings = FALSE)
}

summaries <- lapply(roots, function(root) {
  fromJSON(file.path(root, "analysis_summary.json"), simplifyVector = TRUE)
})

as_number <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_real_)
  out <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(out)) NA_real_ else out
}

value <- function(version, key) {
  as_number(summaries[[version]][[key]])
}

fmt <- function(x, digits = 3L) {
  x <- as_number(x)
  if (is.na(x)) return("NE")
  formatC(x, format = "f", digits = digits)
}

fmt_p <- function(x) {
  x <- as_number(x)
  if (is.na(x)) return("NE")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = 4)
}

fmt_ci <- function(lo, hi, digits = 3L) {
  if (is.na(as_number(lo)) || is.na(as_number(hi))) return("NE")
  sprintf("[%s, %s]", fmt(lo, digits), fmt(hi, digits))
}

yes_no <- function(x) {
  if (is.null(x) || length(x) == 0L || is.na(x[[1]])) return("NE")
  if (isTRUE(x[[1]])) "yes" else "no"
}

read_csv_if_exists <- function(path) {
  if (!file.exists(path)) return(NULL)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

support <- lapply(roots, function(root) {
  robust <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C1_Robustness", "CellCycle_primary_TGI_robustness_summary.csv"
  ))
  ploidy_confounding <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C1_Robustness", "ploidy_confounding_tests.csv"
  ))
  dose_confounding <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C1_Robustness", "dose_confounding_tests.csv"
  ))
  within_dose_tgi <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C1_Robustness", "residualized_TGI_associations.csv"
  ))
  endpoint_tgi_dose <- read_csv_if_exists(file.path(
    root, "B_TGI_GrowthResponse", "B2_DoseSpecific_TGI", "CellCycle", "stats",
    "CellCycle_TGI_AUC_endpoint_30_vs_120_comparisons.csv"
  ))
  fraction <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C3_Composition", "cellcycle_fraction_treatment_tests.csv"
  ))
  cluster_pt <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C3_Composition", "cellcycle_cluster_specific_pseudotime_tests.csv"
  ))
  cluster_tgi <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C3_Composition", "cellcycle_cluster_specific_TGI_associations.csv"
  ))
  module <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C4_GenePrograms", "module_score_pseudotime_trends.csv"
  ))
  de_paths <- c(
    CellCycle_treated_vs_control = file.path(root, "C_SupportingAnalyses", "C5_Pseudobulk", "CellCycle_treated_vs_control_DE.csv"),
    NonCellCycle_treated_vs_control = file.path(root, "C_SupportingAnalyses", "C5_Pseudobulk", "NonCellCycle_treated_vs_control_DE.csv"),
    CellCycle_highTGI_vs_lowTGI = file.path(root, "C_SupportingAnalyses", "C5_Pseudobulk", "CellCycle_highTGI_vs_lowTGI_DE.csv")
  )
  de_counts <- vapply(de_paths, function(path) {
    tab <- read_csv_if_exists(path)
    if (is.null(tab)) return(NA_integer_)
    adjusted_p_col <- intersect(c("padj", "adj.P.Val", "q_value", "FDR"), names(tab))
    if (!length(adjusted_p_col)) return(NA_integer_)
    adjusted_p <- suppressWarnings(as.numeric(tab[[adjusted_p_col[1]]]))
    sum(is.finite(adjusted_p) & adjusted_p < 0.05)
  }, integer(1))
  enrichment <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C5_Pseudobulk", "pathway_enrichment_results.csv"
  ))
  treatment_da <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C6_NeighborhoodDA", "treatment_DA_results.csv"
  ))
  tgi_da <- read_csv_if_exists(file.path(
    root, "C_SupportingAnalyses", "C6_NeighborhoodDA", "TGI_DA_results.csv"
  ))
  list(
    robust = robust,
    ploidy_confounding = ploidy_confounding,
    dose_confounding = dose_confounding,
    within_dose_tgi = within_dose_tgi,
    endpoint_tgi_dose = endpoint_tgi_dose,
    fraction = fraction,
    cluster_pt = cluster_pt,
    cluster_tgi = cluster_tgi,
    module = module,
    de_counts = de_counts,
    enrichment = enrichment,
    treatment_da = treatment_da,
    tgi_da = tgi_da
  )
})

robust_value <- function(version, key) {
  tab <- support[[version]]$robust
  if (is.null(tab) || !key %in% names(tab)) return(NA_real_)
  as_number(tab[[key]][1])
}

endpoint_ploidy_shift_p <- function(version) {
  tab <- support[[version]]$ploidy_confounding
  if (is.null(tab)) return(NA_real_)
  row <- tab[
    tab$compartment == "CellCycle" & tab$sample_set == "treated" &
      tab$shift_metric == "ecdf_rmse" & tab$ploidy_measure == "mean_cell_ploidy" &
      tab$method == "pearson",
    , drop = FALSE
  ]
  if (!nrow(row)) return(NA_real_)
  as_number(row$p_value[1])
}

dose_shift_p <- function(version) {
  tab <- support[[version]]$dose_confounding
  if (is.null(tab)) return(NA_real_)
  row <- tab[
    tab$compartment == "CellCycle" & tab$value == "ecdf_rmse" &
      tab$comparison == "30_vs_120" & tab$test == "t_test",
    , drop = FALSE
  ]
  if (!nrow(row)) return(NA_real_)
  as_number(row$p_value[1])
}

within_dose_centered_value <- function(version, column) {
  tab <- support[[version]]$within_dose_tgi
  if (is.null(tab) || !column %in% names(tab)) return(NA_real_)
  row <- tab[
    tab$compartment == "CellCycle" & tab$analysis == "within_dose_centered" & tab$method == "pearson",
    , drop = FALSE
  ]
  if (!nrow(row)) return(NA_real_)
  as_number(row[[column]][1])
}

endpoint_tgi_day38_value <- function(version, column) {
  tab <- support[[version]]$endpoint_tgi_dose
  if (is.null(tab) || !column %in% names(tab)) return(NA_real_)
  row <- tab[tab$tgi_measure == "TGI_percent_Day_38", , drop = FALSE]
  if (!nrow(row)) return(NA_real_)
  as_number(row[[column]][1])
}

fraction_p <- function(version) {
  tab <- support[[version]]$fraction
  if (is.null(tab)) return(NA_real_)
  row <- tab[tab$fraction == "CellCycle_fraction", , drop = FALSE]
  if (!nrow(row)) return(NA_real_)
  as_number(row$permutation_p[1])
}

cluster10_pt_p <- function(version) {
  tab <- support[[version]]$cluster_pt
  if (is.null(tab)) return(NA_real_)
  row <- tab[tab$cluster == "10" & tab$comparison == "0_vs_30plus120", , drop = FALSE]
  strat_col <- grep("^stratified_by_", names(row), value = TRUE)
  if (length(strat_col)) row <- row[!as.logical(row[[strat_col[1]]]), , drop = FALSE]
  if (!nrow(row)) return(NA_real_)
  as_number(row$p_ecdf_rmse[1])
}

cluster10_tgi <- function(version, key) {
  tab <- support[[version]]$cluster_tgi
  if (is.null(tab)) return(NA_real_)
  row <- tab[tab$cluster == "10", , drop = FALSE]
  if (!nrow(row) || !key %in% names(row)) return(NA_real_)
  as_number(row[[key]][1])
}

module_estimate <- function(version, module, analysis) {
  tab <- support[[version]]$module
  if (is.null(tab)) return(NA_real_)
  row <- tab[tab$module == module & tab$analysis == analysis, , drop = FALSE]
  if (!nrow(row)) return(NA_real_)
  as_number(row$estimate[1])
}

assignments <- lapply(c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24"), function(version) {
  read.csv(
    file.path(roots[[version]], "sample_end_timepoint_ploidy_assignments.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
names(assignments) <- c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24")

thresholds <- vapply(assignments, function(x) unique(x$threshold)[1], numeric(1))

# -----------------------------------------------------------------------------
# Crosswalk and coverage exports
# -----------------------------------------------------------------------------

a41 <- assignments[["ETP_fixed_threshold_2_25"]][, c(
  "sample_id", "original_initial_ploidy", "dose", "dose_mg",
  "n_endpoint_ploidy_cells", "sample_mean_ploidy", "sample_median_ploidy",
  "sample_max_ploidy", "end_timepoint_ploidy_group"
)]
names(a41)[names(a41) == "end_timepoint_ploidy_group"] <- "ETP_fixed_threshold_2_25_ETP_group"
a42 <- assignments[["ETP_boundary_stress_threshold_2_375"]][, c("sample_id", "end_timepoint_ploidy_group")]
names(a42)[2] <- "ETP_boundary_stress_threshold_2_375_ETP_group"
a43 <- assignments[["ETP_reference_balanced_threshold_2_24"]][, c("sample_id", "end_timepoint_ploidy_group")]
names(a43)[2] <- "ETP_reference_balanced_threshold_2_24_ETP_group"
crosswalk <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE, sort = FALSE), list(a41, a42, a43))
crosswalk$initial_ploidy_group <- crosswalk$original_initial_ploidy
crosswalk$ETP_fixed_threshold_2_25_threshold <- thresholds[["ETP_fixed_threshold_2_25"]]
crosswalk$ETP_boundary_stress_threshold_2_375_threshold <- thresholds[["ETP_boundary_stress_threshold_2_375"]]
crosswalk$ETP_reference_balanced_threshold_2_24_threshold <- thresholds[["ETP_reference_balanced_threshold_2_24"]]
crosswalk$ETP_fixed_threshold_2_25_to_ETP_boundary_stress_threshold_2_375_change <- ifelse(
  crosswalk$`ETP_fixed_threshold_2_25_ETP_group` == crosswalk$`ETP_boundary_stress_threshold_2_375_ETP_group`,
  "unchanged",
  paste(crosswalk$`ETP_fixed_threshold_2_25_ETP_group`, "to", crosswalk$`ETP_boundary_stress_threshold_2_375_ETP_group`)
)
crosswalk$ETP_fixed_threshold_2_25_to_ETP_reference_balanced_threshold_2_24_change <- ifelse(
  crosswalk$`ETP_fixed_threshold_2_25_ETP_group` == crosswalk$`ETP_reference_balanced_threshold_2_24_ETP_group`,
  "unchanged",
  paste(crosswalk$`ETP_fixed_threshold_2_25_ETP_group`, "to", crosswalk$`ETP_reference_balanced_threshold_2_24_ETP_group`)
)
crosswalk$ETP_reference_balanced_threshold_2_24_to_ETP_boundary_stress_threshold_2_375_change <- ifelse(
  crosswalk$`ETP_reference_balanced_threshold_2_24_ETP_group` == crosswalk$`ETP_boundary_stress_threshold_2_375_ETP_group`,
  "unchanged",
  paste(crosswalk$`ETP_reference_balanced_threshold_2_24_ETP_group`, "to", crosswalk$`ETP_boundary_stress_threshold_2_375_ETP_group`)
)
crosswalk <- crosswalk[order(crosswalk$dose_mg, crosswalk$sample_id), c(
  "sample_id", "dose", "dose_mg", "initial_ploidy_group",
  "n_endpoint_ploidy_cells", "sample_mean_ploidy", "sample_median_ploidy", "sample_max_ploidy",
  "ETP_fixed_threshold_2_25_threshold", "ETP_fixed_threshold_2_25_ETP_group", "ETP_boundary_stress_threshold_2_375_threshold", "ETP_boundary_stress_threshold_2_375_ETP_group",
  "ETP_reference_balanced_threshold_2_24_threshold", "ETP_reference_balanced_threshold_2_24_ETP_group", "ETP_fixed_threshold_2_25_to_ETP_boundary_stress_threshold_2_375_change",
  "ETP_fixed_threshold_2_25_to_ETP_reference_balanced_threshold_2_24_change", "ETP_reference_balanced_threshold_2_24_to_ETP_boundary_stress_threshold_2_375_change"
)]
write.csv(crosswalk, file.path(comparison_root, "sample_group_crosswalk.csv"), row.names = FALSE, na = "NA")

doses <- sort(unique(crosswalk$dose_mg))
make_coverage <- function(version, group_var, groups, threshold, grouping_variable) {
  rows <- expand.grid(dose_mg = doses, group = groups, stringsAsFactors = FALSE)
  rows$dose <- paste0(rows$dose_mg, "mg/kg")
  rows$n_samples <- vapply(seq_len(nrow(rows)), function(i) {
    sum(crosswalk$dose_mg == rows$dose_mg[i] & crosswalk[[group_var]] == rows$group[i], na.rm = TRUE)
  }, integer(1))
  rows$empty_group <- rows$n_samples == 0
  rows$version <- version
  rows$grouping_variable <- grouping_variable
  rows$threshold <- threshold
  rows$estimability_note <- ifelse(
    rows$n_samples == 0,
    "group absent in this dose; group contrasts/interactions requiring this cell are not estimable",
    ifelse(rows$n_samples == 1, "single sample in this dose-by-group cell; estimates are highly sparse", "represented")
  )
  rows[, c("version", "grouping_variable", "threshold", "dose", "dose_mg", "group", "n_samples", "empty_group", "estimability_note")]
}

coverage <- rbind(
  make_coverage("initial_ploidy", "initial_ploidy_group", c("2N", "4N"), NA_real_, "baseline initial ploidy"),
  make_coverage("ETP_fixed_threshold_2_25", "ETP_fixed_threshold_2_25_ETP_group", c("ETP-lower", "ETP-higher"), thresholds[["ETP_fixed_threshold_2_25"]], "post-treatment endpoint mean ploidy"),
  make_coverage("ETP_boundary_stress_threshold_2_375", "ETP_boundary_stress_threshold_2_375_ETP_group", c("ETP-lower", "ETP-higher"), thresholds[["ETP_boundary_stress_threshold_2_375"]], "post-treatment endpoint mean ploidy"),
  make_coverage("ETP_reference_balanced_threshold_2_24", "ETP_reference_balanced_threshold_2_24_ETP_group", c("ETP-lower", "ETP-higher"), thresholds[["ETP_reference_balanced_threshold_2_24"]], "post-treatment endpoint mean ploidy")
)
write.csv(coverage, file.path(comparison_root, "dose_group_coverage_comparison.csv"), row.names = FALSE, na = "NA")

coverage_sentence <- function(version) {
  tab <- coverage[coverage$version == version, ]
  group_levels <- if (version == "initial_ploidy") c("2N", "4N") else c("ETP-lower", "ETP-higher")
  paste(vapply(doses, function(dose) {
    z <- tab[tab$dose_mg == dose, ]
    counts <- setNames(z$n_samples, z$group)
    sprintf(
      "%d mg/kg: %s n=%d, %s n=%d",
      dose,
      group_levels[1], counts[[group_levels[1]]],
      group_levels[2], counts[[group_levels[2]]]
    )
  }, character(1)), collapse = "; ")
}

untreated_initial_composition <- function(version) {
  z <- assignments[[version]][assignments[[version]]$dose_mg == 0, , drop = FALSE]
  groups <- c("ETP-lower", "ETP-higher")
  paste(vapply(groups, function(group) {
    w <- z[z$end_timepoint_ploidy_group == group, , drop = FALSE]
    n2 <- sum(w$original_initial_ploidy == "2N")
    n4 <- sum(w$original_initial_ploidy == "4N")
    sprintf("%s: initial 2N n=%d, initial 4N n=%d", group, n2, n4)
  }, character(1)), collapse = "; ")
}

# -----------------------------------------------------------------------------
# Continuous endpoint-ploidy exploration (C7)
# -----------------------------------------------------------------------------

continuous_stats <- function(root, compartment) {
  table_root <- file.path(
    root, "C_SupportingAnalyses", "C7_EndTimePointPloidy_Exploration",
    compartment, "tables"
  )
  tgi <- read_csv_if_exists(file.path(
    table_root, paste0(compartment, "_endpoint_ploidy_mean_median_max_vs_TGI_plot_data.csv")
  ))
  pt <- read_csv_if_exists(file.path(
    table_root, paste0(compartment, "_endpoint_ploidy_mean_median_max_vs_pseudotime_plot_data.csv")
  ))
  cells <- read_csv_if_exists(file.path(
    table_root, paste0(compartment, "_single_cell_ploidy_vs_single_cell_pseudotime_plot_data.csv")
  ))
  cor_safe <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 3L || length(unique(x[ok])) < 2L || length(unique(y[ok])) < 2L) return(NA_real_)
    cor(x[ok], y[ok], method = "pearson")
  }
  metric_map <- c(mean = "mean_cell_ploidy", median = "median_cell_ploidy", max = "max_cell_ploidy")
  tgi_all <- tgi_treated <- setNames(rep(NA_real_, length(metric_map)), names(metric_map))
  if (!is.null(tgi)) {
    for (nm in names(metric_map)) {
      z <- tgi[tgi$endpoint_ploidy_metric == metric_map[[nm]] & tgi$tgi_measure == "TGI_percent_auc", , drop = FALSE]
      tgi_all[[nm]] <- cor_safe(z$endpoint_ploidy_value, z$tgi_value)
      z <- z[z$dose_mg > 0, , drop = FALSE]
      tgi_treated[[nm]] <- cor_safe(z$endpoint_ploidy_value, z$tgi_value)
    }
  }
  mean_pt_treated <- NA_real_
  if (!is.null(pt)) {
    z <- pt[
      pt$endpoint_ploidy_metric == "mean_cell_ploidy" &
        pt$pseudotime_metric == "mean_pseudotime" & pt$dose_mg > 0,
      , drop = FALSE
    ]
    mean_pt_treated <- cor_safe(z$endpoint_ploidy_value, z$pseudotime_value)
  }
  cell_pooled <- NA_real_
  cell_by_dose <- setNames(rep(NA_real_, length(doses)), doses)
  if (!is.null(cells)) {
    cell_pooled <- cor_safe(cells$cell_ploidy, cells$pseudotime)
    for (dose in doses) {
      z <- cells[cells$dose_mg == dose, , drop = FALSE]
      cell_by_dose[[as.character(dose)]] <- cor_safe(z$cell_ploidy, z$pseudotime)
    }
  }
  list(
    tgi_all = tgi_all,
    tgi_treated = tgi_treated,
    mean_pt_treated = mean_pt_treated,
    cell_pooled = cell_pooled,
    cell_by_dose = cell_by_dose
  )
}

# C7 values are identical between ETP_fixed_threshold_2_25, ETP_boundary_stress_threshold_2_375, and ETP_reference_balanced_threshold_2_24 because the continuous
# quantities do not depend on the binary threshold. Read from the formal ETP_fixed_threshold_2_25
# tables and state that identity explicitly in the reports.
c7 <- list(
  CellCycle = continuous_stats(roots[["ETP_fixed_threshold_2_25"]], "CellCycle"),
  NonCellCycle = continuous_stats(roots[["ETP_fixed_threshold_2_25"]], "NonCellCycle")
)

# -----------------------------------------------------------------------------
# Master metric comparison table
# -----------------------------------------------------------------------------

metric_spec <- list(
  c("cellcycle_treated_vs_control_p", "CellCycle 0 vs treated permutation p", "A1 treatment redistribution", "cellcycle_treated_vs_control_p", "cellcycle_treated_vs_control_p", "cellcycle_treated_vs_control_p", "Directly comparable; the unstratified test is unchanged."),
  c("cellcycle_stratified_p", "CellCycle stratified 0 vs treated p", "A1 treatment redistribution", "cellcycle_treated_vs_control_ploidy_stratified_p", "cellcycle_treated_vs_control_ploidy_stratified_p", "cellcycle_treated_vs_control_ploidy_stratified_p", "Strata differ: baseline initial ploidy in initial_ploidy versus post-treatment ETP group in ETP."),
  c("noncellcycle_treated_vs_control_p", "NonCellCycle 0 vs treated permutation p", "A1 comparison compartment", "noncellcycle_treated_vs_control_p", "noncellcycle_treated_vs_control_p", "noncellcycle_treated_vs_control_p", "Directly comparable unstratified test."),
  c("cellcycle_group_DID", "CellCycle group treatment-effect DID ECDF RMSE", "A2 effect modification", "cellcycle_pseudotime_ploidy_did_ecdf_rmse", "cellcycle_pseudotime_ploidy_did_ecdf_rmse", "cellcycle_pseudotime_ploidy_did_ecdf_rmse", "Group contrast changes from initial 4N-vs-2N to endpoint higher-vs-lower."),
  c("cellcycle_group_DID_p", "CellCycle group treatment-effect DID p", "A2 effect modification", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "Threshold-sensitive; ETP_boundary_stress_threshold_2_375 is sparse and has no bootstrap CI."),
  c("noncellcycle_group_DID_p", "NonCellCycle group treatment-effect DID p", "A2 effect modification", "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p", "ETP_boundary_stress_threshold_2_375 also yields nominal significance in the comparison compartment."),
  c("cellcycle_shift_TGI_r", "CellCycle shift vs AUC-TGI Pearson r", "B1 response association", "cellcycle_TGI_AUC_pearson_r", "cellcycle_TGI_AUC_pearson_r", "cellcycle_TGI_AUC_pearson_r", "ETP-matched reference changes the sample shift metric across ETP thresholds."),
  c("cellcycle_shift_TGI_p", "CellCycle shift vs AUC-TGI permutation p", "B1 response association", "cellcycle_TGI_AUC_pearson_perm_p", "cellcycle_TGI_AUC_pearson_perm_p", "cellcycle_TGI_AUC_pearson_perm_p", "Association is exploratory and threshold sensitive."),
  c("cellcycle_adjusted_shift_p", "Dose/group-adjusted CellCycle shift term p", "B1 response association", "cellcycle_AUC_TGI_adjusted_shift_p", "cellcycle_AUC_TGI_adjusted_shift_p", "cellcycle_AUC_TGI_adjusted_shift_p", "The group covariate is baseline in initial_ploidy but post-treatment in ETP."),
  c("cellcycle_group_covariate_p", "Adjusted group covariate p", "B1 response association", "cellcycle_AUC_TGI_adjusted_initial_ploidy_p", "cellcycle_AUC_TGI_adjusted_end_timepoint_ploidy_group_p", "cellcycle_AUC_TGI_adjusted_end_timepoint_ploidy_group_p", "Not the same estimand: baseline initial ploidy versus endpoint group."),
  c("cellcycle_mean_ETP_unadjusted_slope", "CellCycle AUC-TGI vs mean ETP standardized slope", "B1 continuous ETP association", "cellcycle_AUC_TGI_mean_ETP_unadjusted_slope", "cellcycle_AUC_TGI_mean_ETP_unadjusted_slope", "cellcycle_AUC_TGI_mean_ETP_unadjusted_slope", "Directly comparable: the same treated samples and continuous mean ETP are used in every method."),
  c("cellcycle_mean_ETP_unadjusted_p", "CellCycle AUC-TGI vs mean ETP unadjusted p", "B1 continuous ETP association", "cellcycle_AUC_TGI_mean_ETP_unadjusted_p", "cellcycle_AUC_TGI_mean_ETP_unadjusted_p", "cellcycle_AUC_TGI_mean_ETP_unadjusted_p", "Directly comparable unadjusted linear association."),
  c("cellcycle_mean_ETP_adjusted_p", "CellCycle AUC-TGI vs mean ETP dose/group-adjusted p", "B1 continuous ETP association", "cellcycle_AUC_TGI_mean_ETP_adjusted_p", "cellcycle_AUC_TGI_mean_ETP_adjusted_p", "cellcycle_AUC_TGI_mean_ETP_adjusted_p", "The continuous ETP predictor is identical, but the adjustment group is baseline initial ploidy versus endpoint ETP group."),
  c("cellcycle_mean_ETP_group_interaction_p", "CellCycle mean-ETP-by-group interaction p", "B1 continuous ETP association", "cellcycle_AUC_TGI_mean_ETP_by_ploidy_interaction_p", "cellcycle_AUC_TGI_mean_ETP_by_ploidy_interaction_p", "cellcycle_AUC_TGI_mean_ETP_by_ploidy_interaction_p", "Exploratory small-n interaction; group definitions differ by method."),
  c("compartment_delta_r", "CellCycle minus NonCellCycle TGI-correlation delta r", "C2 compartment comparison", "compartment_delta_r", "compartment_delta_r", "compartment_delta_r", "Same label-swap question, but ETP reference matching changes each shift metric."),
  c("compartment_delta_r_p", "Compartment delta-r permutation p", "C2 compartment comparison", "compartment_delta_r_perm_p", "compartment_delta_r_perm_p", "compartment_delta_r_perm_p", "None establishes a statistically proven compartment difference in the TGI association."),
  c("dose_TGI_diff", "AUC-TGI 120 minus 30 mg/kg", "B2 dose response", "cellcycle_TGI_AUC_120_minus_30_mean_diff", "cellcycle_TGI_AUC_120_minus_30_mean_diff", "cellcycle_TGI_AUC_120_minus_30_mean_diff", "Response endpoint is identical across versions."),
  c("dose_TGI_diff_p", "AUC-TGI 30 vs 120 permutation p", "B2 dose response", "cellcycle_TGI_AUC_30_vs_120_perm_p", "cellcycle_TGI_AUC_30_vs_120_perm_p", "cellcycle_TGI_AUC_30_vs_120_perm_p", "Directly comparable; no clear dose separation."),
  c("group_lower_treatment_p", "Within lower/2N group CellCycle treatment p", "B3 group response", "cellcycle_2N_treated_vs_control_ecdf_p", "cellcycle_ETP_lower_treated_vs_control_ecdf_p", "cellcycle_ETP_lower_treated_vs_control_ecdf_p", "initial_ploidy and ETP group meanings differ."),
  c("group_higher_treatment_p", "Within higher/4N group CellCycle treatment p", "B3 group response", "cellcycle_4N_treated_vs_control_ecdf_p", "cellcycle_ETP_higher_treated_vs_control_ecdf_p", "cellcycle_ETP_higher_treated_vs_control_ecdf_p", "ETP_boundary_stress_threshold_2_375 is not estimable because endpoint-high contains only one control and one treated sample; 120 mg/kg endpoint-high is also absent."),
  c("group_shift_diff", "Treated group shift difference", "B3 group response", "cellcycle_ecdf_rmse_4N_minus_2N_mean_diff", "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_mean_diff", "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_mean_diff", "Different grouping definitions; ETP_boundary_stress_threshold_2_375 is not estimable."),
  c("group_TGI_diff", "Treated group AUC-TGI difference", "B3 group response", "tgi_auc_4N_minus_2N_mean_diff", "tgi_auc_ETP_higher_minus_ETP_lower_mean_diff", "tgi_auc_ETP_higher_minus_ETP_lower_mean_diff", "Different grouping definitions; ETP_boundary_stress_threshold_2_375 is not estimable."),
  c("growth_spline_interaction_p", "Spline time-by-dose-by-group LRT p", "B4 longitudinal growth", "growth_curve_spline_time_dose_ploidy_lrt_p", "growth_curve_spline_time_dose_ploidy_lrt_p", "growth_curve_spline_time_dose_ploidy_lrt_p", "The group is baseline in initial_ploidy, post-treatment in ETP; ETP_boundary_stress_threshold_2_375 is rank deficient/not estimable."),
  c("growth_linear_interaction_p", "Linear time-by-dose-by-group LRT p", "B4 longitudinal growth", "growth_curve_linear_time_dose_ploidy_lrt_p", "growth_curve_linear_time_dose_ploidy_lrt_p", "growth_curve_linear_time_dose_ploidy_lrt_p", "Sensitivity model; same estimability caveat."),
  c("day38_30_group_DID_p", "Day 38 group DID p at 30 mg/kg", "B4 longitudinal growth", "growth_curve_day_final_30mg_ploidy_did_p", "growth_curve_day_final_30mg_ploidy_did_p", "growth_curve_day_final_30mg_ploidy_did_p", "Different grouping definitions."),
  c("day38_120_group_DID_p", "Day 38 group DID p at 120 mg/kg", "B4 longitudinal growth", "growth_curve_day_final_120mg_ploidy_did_p", "growth_curve_day_final_120mg_ploidy_did_p", "growth_curve_day_final_120mg_ploidy_did_p", "ETP_boundary_stress_threshold_2_375 is not estimable because 120 mg/kg endpoint-high is absent.")
)

master <- do.call(rbind, lapply(metric_spec, function(spec) {
  data.frame(
    metric_id = spec[1], metric_label = spec[2], analysis_family = spec[3],
    initial_ploidy = value("initial_ploidy", spec[4]), `ETP_fixed_threshold_2_25` = value("ETP_fixed_threshold_2_25", spec[5]), `ETP_boundary_stress_threshold_2_375` = value("ETP_boundary_stress_threshold_2_375", spec[6]), `ETP_reference_balanced_threshold_2_24` = value("ETP_reference_balanced_threshold_2_24", spec[5]),
    comparability_note = spec[7], check.names = FALSE, stringsAsFactors = FALSE
  )
}))

extra_master <- data.frame(
  metric_id = c(
    "robust_LOO_min_r", "robust_bootstrap_r_ci025", "robust_bootstrap_r_ci975",
    "cellcycle_fraction_p", "cluster10_shift_p", "cluster10_TGI_r", "cluster10_TGI_p"
  ),
  metric_label = c(
    "Leave-one-out minimum CellCycle shift-TGI r",
    "Bootstrap CellCycle shift-TGI r CI lower",
    "Bootstrap CellCycle shift-TGI r CI upper",
    "CellCycle fraction treatment permutation p",
    "Cluster 10 0-vs-treated ECDF RMSE p",
    "Cluster 10 shift vs AUC-TGI Pearson r",
    "Cluster 10 shift vs AUC-TGI permutation p"
  ),
  analysis_family = c(rep("C1 robustness", 3), "C3 composition", rep("C3 cluster analysis", 3)),
  initial_ploidy = c(
    robust_value("initial_ploidy", "leave_one_out_min_pearson_r"), robust_value("initial_ploidy", "bootstrap_pearson_r_ci025"),
    robust_value("initial_ploidy", "bootstrap_pearson_r_ci975"), fraction_p("initial_ploidy"), cluster10_pt_p("initial_ploidy"),
    cluster10_tgi("initial_ploidy", "pearson_r"), cluster10_tgi("initial_ploidy", "pearson_p_permutation_two_sided")
  ),
  `ETP_fixed_threshold_2_25` = c(
    robust_value("ETP_fixed_threshold_2_25", "leave_one_out_min_pearson_r"), robust_value("ETP_fixed_threshold_2_25", "bootstrap_pearson_r_ci025"),
    robust_value("ETP_fixed_threshold_2_25", "bootstrap_pearson_r_ci975"), fraction_p("ETP_fixed_threshold_2_25"), cluster10_pt_p("ETP_fixed_threshold_2_25"),
    cluster10_tgi("ETP_fixed_threshold_2_25", "pearson_r"), cluster10_tgi("ETP_fixed_threshold_2_25", "pearson_p_permutation_two_sided")
  ),
  `ETP_boundary_stress_threshold_2_375` = c(
    robust_value("ETP_boundary_stress_threshold_2_375", "leave_one_out_min_pearson_r"), robust_value("ETP_boundary_stress_threshold_2_375", "bootstrap_pearson_r_ci025"),
    robust_value("ETP_boundary_stress_threshold_2_375", "bootstrap_pearson_r_ci975"), fraction_p("ETP_boundary_stress_threshold_2_375"), cluster10_pt_p("ETP_boundary_stress_threshold_2_375"),
    cluster10_tgi("ETP_boundary_stress_threshold_2_375", "pearson_r"), cluster10_tgi("ETP_boundary_stress_threshold_2_375", "pearson_p_permutation_two_sided")
  ),
  `ETP_reference_balanced_threshold_2_24` = c(
    robust_value("ETP_reference_balanced_threshold_2_24", "leave_one_out_min_pearson_r"), robust_value("ETP_reference_balanced_threshold_2_24", "bootstrap_pearson_r_ci025"),
    robust_value("ETP_reference_balanced_threshold_2_24", "bootstrap_pearson_r_ci975"), fraction_p("ETP_reference_balanced_threshold_2_24"), cluster10_pt_p("ETP_reference_balanced_threshold_2_24"),
    cluster10_tgi("ETP_reference_balanced_threshold_2_24", "pearson_r"), cluster10_tgi("ETP_reference_balanced_threshold_2_24", "pearson_p_permutation_two_sided")
  ),
  comparability_note = c(
    rep("Reference matching changes the ETP sample-level shift metric.", 3),
    "Composition conclusion is stable.",
    "Unstratified cluster-level treatment contrast; Monte Carlo variation is expected.",
    "Cluster shift metric inherits the version-specific reference.",
    "Cluster shift metric inherits the version-specific reference."
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
master <- rbind(master, extra_master)
write.csv(master, file.path(comparison_root, "comparison_master_metrics.csv"), row.names = FALSE, na = "NA")

markdown_metric_table <- function(rows = c(
  "cellcycle_treated_vs_control_p", "noncellcycle_treated_vs_control_p", "cellcycle_group_DID_p",
  "cellcycle_shift_TGI_r", "cellcycle_shift_TGI_p", "cellcycle_adjusted_shift_p",
  "compartment_delta_r_p", "dose_TGI_diff_p", "growth_spline_interaction_p"
)) {
  tab <- master[match(rows, master$metric_id), ]
  c(
    "| Metric | initial_ploidy | ETP_fixed_threshold_2_25 | ETP_boundary_stress_threshold_2_375 | ETP_reference_balanced_threshold_2_24 |",
    "|---|---:|---:|---:|---:|",
    vapply(seq_len(nrow(tab)), function(i) {
      sprintf(
        "| %s | %s | %s | %s | %s |",
        tab$metric_label[i], fmt(tab$initial_ploidy[i], 4), fmt(tab$`ETP_fixed_threshold_2_25`[i], 4), fmt(tab$`ETP_boundary_stress_threshold_2_375`[i], 4), fmt(tab$`ETP_reference_balanced_threshold_2_24`[i], 4)
      )
    }, character(1))
  )
}

# -----------------------------------------------------------------------------
# Per-version comprehensive reports
# -----------------------------------------------------------------------------

comprehensive_report <- function(version) {
  threshold <- thresholds[[version]]
  is_v42 <- identical(version, "ETP_boundary_stress_threshold_2_375")
  is_v43 <- identical(version, "ETP_reference_balanced_threshold_2_24")
  r <- support[[version]]$robust
  c10_p <- cluster10_pt_p(version)
  c10_r <- cluster10_tgi(version, "pearson_r")
  c10_tgi_p <- cluster10_tgi(version, "pearson_p_permutation_two_sided")
  g2m <- module_estimate(version, "G2M_checkpoint", "score_vs_pseudotime_control")
  e2f <- module_estimate(version, "E2F_S_phase_DNA_replication", "score_vs_pseudotime_control")
  de_counts <- support[[version]]$de_counts
  enrich <- support[[version]]$enrichment
  e2f_q <- if (!is.null(enrich)) as_number(enrich$q_value[enrich$module == "E2F_S_phase_DNA_replication"][1]) else NA_real_
  p53_q <- if (!is.null(enrich)) as_number(enrich$q_value[enrich$module == "p53_apoptosis_senescence"][1]) else NA_real_
  treatment_da_sig <- if (!is.null(support[[version]]$treatment_da)) {
    sum(is.finite(support[[version]]$treatment_da$q_value) & support[[version]]$treatment_da$q_value < 0.05)
  } else NA_integer_
  higher_note <- if (is_v42) {
    "120 mg/kg 没有 ETP-higher 样本，因此任何需要该 dose-by-group 单元的 higher-vs-lower 对比、交互项和纵向混合模型均不可估（NE）。"
  } else if (is_v43) {
    paste0(
      "ETP_reference_balanced_threshold_2_24 的 2.24 阈值专门使两个 untreated ETP reference 都混有 initial 2N 与 initial 4N：",
      untreated_initial_composition(version),
      "。但 Dose 内的 initial 来源与 ETP 仍不完全交叉，因此这只能减少 reference 被单一 initial 来源垄断，不能消除所有初始状态嵌套。"
    )
  } else {
    "每个剂量都同时有 ETP-lower 和 ETP-higher，但 30 mg/kg lower 与 120 mg/kg higher 均只有 1 个样本，分层估计仍高度稀疏。"
  }
  a2_note <- if (is_v42) {
    sprintf(
      paste0(
        "CellCycle DID RMSE=%s、p=%s，看似达到名义显著；但 NonCellCycle 也出现 DID RMSE=%s、p=%s，",
        "两者 bootstrap CI 均不可得。ETP-higher 在 control 与 treated 中都只有 1 个样本，检验仅有 64 个精确排列；",
        "结合极端不平衡和 120 mg/kg 缺组，这更像阈值/可估性敏感信号，不能据此宣称 CellCycle 特异的 ETP effect modification。"
      ),
      fmt(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse"), 4),
      fmt_p(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
      fmt(value(version, "noncellcycle_pseudotime_ploidy_did_ecdf_rmse"), 4),
      fmt_p(value(version, "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
    )
  } else {
    sprintf(
      paste0(
        "CellCycle DID RMSE=%s、p=%s、bootstrap CI=%s；NonCellCycle p=%s。",
        "没有检测到 endpoint-defined group 对治疗相关 pseudotime redistribution 的稳定修饰。"
      ),
      fmt(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse"), 4),
      fmt_p(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
      fmt_ci(
        value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci025"),
        value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_ci975"), 4
      ),
      fmt_p(value(version, "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
    )
  }
  b3_note <- if (is_v42) {
    "ETP-higher 的组内 treated-vs-control、treated higher-vs-lower shift/TGI、equivalence、交互及 spline mixed model 均为 NE；NE 是设计矩阵缺少组别组合，不是“无差异”。"
  } else {
    sprintf(
      paste0(
        "treated ETP-higher minus ETP-lower：CellCycle shift 差=%s、p=%s、CI=%s；",
        "AUC-TGI 差=%s、p=%s、CI=%s。两项等效性界值均未满足，因此既不能说有差异，也不能说等效。"
      ),
      fmt(value(version, "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_mean_diff"), 4),
      fmt_p(value(version, "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_perm_p")),
      fmt_ci(
        value(version, "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci025"),
        value(version, "cellcycle_ecdf_rmse_ETP_higher_minus_ETP_lower_ci975"), 4
      ),
      fmt(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_mean_diff"), 3),
      fmt_p(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_perm_p")),
      fmt_ci(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_ci025"), value(version, "tgi_auc_ETP_higher_minus_ETP_lower_ci975"), 3)
    )
  }
  growth_note <- if (is_v42) {
    sprintf(
      "spline 与 linear time-by-dose-by-ETP-group LRT 均为 NE；Day 38 的 30 mg/kg DID=%s（p=%s），120 mg/kg 因无 endpoint-high 为 NE。",
      fmt(value(version, "growth_curve_day_final_30mg_ploidy_did_estimate"), 4),
      fmt_p(value(version, "growth_curve_day_final_30mg_ploidy_did_p"))
    )
  } else {
    sprintf(
      paste0(
        "spline LRT p=%s，linear-time sensitivity p=%s；Day 38 DID p：30 mg/kg=%s、120 mg/kg=%s。",
        "没有检测到 endpoint group 定义的纵向异质性。"
      ),
      fmt_p(value(version, "growth_curve_spline_time_dose_ploidy_lrt_p")),
      fmt_p(value(version, "growth_curve_linear_time_dose_ploidy_lrt_p")),
      fmt_p(value(version, "growth_curve_day_final_30mg_ploidy_did_p")),
      fmt_p(value(version, "growth_curve_day_final_120mg_ploidy_did_p"))
    )
  }

  c(
    sprintf("# %s 全面结果解读：EndTimePoint ploidy 分组的 pseudotime–TGI 分析", version),
    "",
    "## 1. 分析对象、分组与科学问题",
    "",
    sprintf(
      paste0(
        "%s 将 16 个肿瘤样本的终点细胞 ploidy 先汇总到样本均值，再用固定阈值 %.3f 对**全部样本**分组：",
        "sample mean ploidy > threshold 为 ETP-higher，否则为 ETP-lower。这个分组不受原始 initial 2N/4N 身份限制。"
      ), version, threshold
    ),
    "",
    sprintf("剂量覆盖：%s。", coverage_sentence(version)),
    higher_note,
    "",
    paste0(
      "核心科学问题为：(1) gemcitabine 是否重排 CellCycle pseudotime 分布；",
      "(2) endpoint-defined ploidy group 是否与重排幅度或 TGI 有关；",
      "(3) pseudotime shift 是否与 AUC-TGI 同向；(4)这些结论对剂量、compartment、样本和分析方法是否稳健。"
    ),
    "",
    "**重要 estimand 边界。** ETP 是处理后变量。ETP 不只是把图上的标签从 initial ploidy 换成 ETP group；",
    "它还按 ETP group 构造匹配的 untreated ECDF reference，因此每个样本的 ECDF shift 及其与 TGI 的相关会随阈值改变。",
    "按处理后的 ETP 分层、匹配或调整，可能条件化于治疗下游的 mediator/collider/selection 结构。",
    "所以 ETP 应被表述为 endpoint-defined descriptive association / sensitivity analysis，不是 baseline predictive biomarker 分析，也不支持因果中介结论。",
    "",
    "## 2. 主结论概览",
    "",
    sprintf(
      "- CellCycle 0 vs treated：p=%s；按 ETP group 分层后 p=%s。治疗相关 CellCycle pseudotime redistribution 是最稳定的主结果。",
      fmt_p(value(version, "cellcycle_treated_vs_control_p")),
      fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p"))
    ),
    sprintf(
      "- NonCellCycle 0 vs treated：p=%s，没有相应的总体 redistribution。",
      fmt_p(value(version, "noncellcycle_treated_vs_control_p"))
    ),
    sprintf(
      "- CellCycle shift vs AUC-TGI：Pearson r=%s、permutation p=%s；方向为正，但统计证据受阈值和 n=8 treated 样本限制。",
      fmt(value(version, "cellcycle_TGI_AUC_pearson_r"), 3),
      fmt_p(value(version, "cellcycle_TGI_AUC_pearson_perm_p"))
    ),
    sprintf(
      "- 120 minus 30 mg/kg AUC-TGI=%s percentage points、p=%s，没有清楚的剂量组间 TGI 分离。",
      fmt(value(version, "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 3),
      fmt_p(value(version, "cellcycle_TGI_AUC_30_vs_120_perm_p"))
    ),
    "",
    "## 3. A1：治疗是否改变 pseudotime 分布？",
    "",
    sprintf(
      paste0(
        "CellCycle sample-aware ECDF test 的 0-vs-treated p=%s，ETP-group-stratified p=%s；",
        "NonCellCycle p=%s。由于 permutation 的单位是肿瘤样本而不是细胞，这一主检验保持了正确的生物学重复层级。",
        "结果支持 gemcitabine 与 CellCycle 状态分布重排相关，而不是所有肿瘤细胞 compartment 的普遍 pseudotime 漂移。"
      ),
      fmt_p(value(version, "cellcycle_treated_vs_control_p")),
      fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p")),
      fmt_p(value(version, "noncellcycle_treated_vs_control_p"))
    ),
    "",
    "## 4. A2：ETP group 是否修饰 pseudotime treatment effect？",
    "",
    a2_note,
    "这里的 DID 是 ETP-higher 与 ETP-lower 的 treated-minus-control ECDF 曲线差异；它描述 endpoint-defined strata，不是随机化前 effect modifier。",
    "",
    "## 5. B1/B2：pseudotime shift 与肿瘤生长抑制",
    "",
    sprintf(
      paste0(
        "treated-only CellCycle ECDF RMSE 与 AUC-TGI 的 Pearson r=%s、p=%s。加入 dose 与 ETP group 后，",
        "shift term p=%s、ETP group term p=%s、shift-by-group interaction p=%s。",
        "NonCellCycle 的总体 r=%s、p=%s；compartment delta-r=%s、label-swap p=%s。"
      ),
      fmt(value(version, "cellcycle_TGI_AUC_pearson_r"), 3),
      fmt_p(value(version, "cellcycle_TGI_AUC_pearson_perm_p")),
      fmt_p(value(version, "cellcycle_AUC_TGI_adjusted_shift_p")),
      fmt_p(value(version, "cellcycle_AUC_TGI_adjusted_end_timepoint_ploidy_group_p")),
      fmt_p(value(version, "cellcycle_AUC_TGI_shift_by_ploidy_interaction_p")),
      fmt(value(version, "noncellcycle_TGI_AUC_pearson_r"), 3),
      fmt_p(value(version, "noncellcycle_TGI_AUC_pearson_perm_p")),
      fmt(value(version, "compartment_delta_r"), 3),
      fmt_p(value(version, "compartment_delta_r_perm_p"))
    ),
    paste0(
      "正方向可作为“CellCycle state redistribution 与 response alignment”的假设线索；",
      "但 delta-r 不显著，不能说已统计学证明 CellCycle-specific TGI association。",
      "调整 ETP group 的模型又条件化于处理后信息，因此不能升级为预测或因果结论。"
    ),
    "",
    sprintf(
      paste0(
        "剂量分层：30 mg/kg r=%s、p=%s；120 mg/kg r=%s、p=%s；slope-difference p=%s。",
        "每个 treated dose 只有 4 个样本，因此 within-dose 相关只作探索性方向判断。"
      ),
      fmt(value(version, "cellcycle_TGI_AUC_shift_30mg_pearson_r"), 3),
      fmt_p(value(version, "cellcycle_TGI_AUC_shift_30mg_perm_p")),
      fmt(value(version, "cellcycle_TGI_AUC_shift_120mg_pearson_r"), 3),
      fmt_p(value(version, "cellcycle_TGI_AUC_shift_120mg_perm_p")),
      fmt_p(value(version, "cellcycle_shift_TGI_AUC_dose_interaction_p"))
    ),
    "",
    "## 6. B3：endpoint group 的 shift/TGI response 对比",
    "",
    b3_note,
    "",
    "## 7. B4：纵向肿瘤生长混合模型",
    "",
    growth_note,
    "因为 ETP group 在肿瘤生长过程之后/末端定义，这一模型是按最终状态进行的描述性轨迹分层，不是 baseline response heterogeneity 模型。",
    "",
    "## 8. 稳健性、compartment 和组成",
    "",
    sprintf(
      paste0(
        "CellCycle shift–TGI 的 leave-one-out 最小 r=%s，bootstrap Pearson 95%% interval=%s；",
        "方向没有因删除单一样本而反转，但区间%s，说明小样本不确定性明显。"
      ),
      fmt(robust_value(version, "leave_one_out_min_pearson_r"), 3),
      fmt_ci(robust_value(version, "bootstrap_pearson_r_ci025"), robust_value(version, "bootstrap_pearson_r_ci975"), 3),
      if (robust_value(version, "bootstrap_pearson_r_ci025") <= 0) "跨越 0" else "虽未跨越 0 但很宽"
    ),
    sprintf(
      "CellCycle fraction 的 treatment permutation p=%s，主结果不支持被简单的 compartment abundance 变化解释。",
      fmt_p(fraction_p(version))
    ),
    sprintf(
      paste0(
        "Cluster 10 的 0-vs-treated ECDF RMSE p=%s；其 shift–TGI r=%s、p=%s。",
        "可把 cluster 10 作为潜在贡献亚群，但其 TGI 关联仍是探索性。"
      ), fmt_p(c10_p), fmt(c10_r, 3), fmt_p(c10_tgi_p)
    ),
    "",
    "## 9. C7：连续 endpoint ploidy 探索",
    "",
    sprintf(
      paste0(
        "ETP_fixed_threshold_2_25、ETP_boundary_stress_threshold_2_375 与 ETP_reference_balanced_threshold_2_24 的 C7 连续数据完全相同，因为连续统计不使用二分阈值。CellCycle endpoint mean/median/max ploidy vs AUC-TGI：",
        "all-sample r=%s/%s/%s；treated-only r=%s/%s/%s。NonCellCycle treated-only r=%s/%s/%s。"
      ),
      fmt(c7$CellCycle$tgi_all[["mean"]], 3), fmt(c7$CellCycle$tgi_all[["median"]], 3), fmt(c7$CellCycle$tgi_all[["max"]], 3),
      fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
      fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3)
    ),
    sprintf(
      "mean endpoint ploidy vs mean pseudotime（treated-only）：CellCycle r=%s，NonCellCycle r=%s。",
      fmt(c7$CellCycle$mean_pt_treated, 3), fmt(c7$NonCellCycle$mean_pt_treated, 3)
    ),
    sprintf(
      paste0(
        "single-cell ploidy vs single-cell pseudotime pooled r：CellCycle=%s（dose 0/30/120：%s/%s/%s）；",
        "NonCellCycle=%s（dose 0/30/120：%s/%s/%s）。"
      ),
      fmt(c7$CellCycle$cell_pooled, 3),
      fmt(c7$CellCycle$cell_by_dose[["0"]], 3), fmt(c7$CellCycle$cell_by_dose[["30"]], 3), fmt(c7$CellCycle$cell_by_dose[["120"]], 3),
      fmt(c7$NonCellCycle$cell_pooled, 3),
      fmt(c7$NonCellCycle$cell_by_dose[["0"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["30"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["120"]], 3)
    ),
    paste0(
      "这些是 scatter-plot exploration，不是正式因果或预测检验。尤其 pooled single-cell r 把同一肿瘤内细胞当作点，",
      "细胞不独立，不能用其 p 值或 r 进行样本级推断。连续样本级模型优于任意二分，但仍须保留 treatment timing 与样本数限制。"
    ),
    "",
    "## 10. C4–C6 optional modules",
    "",
    sprintf(
      paste0(
        "C4 的聚合 module-score 结果支持 pseudotime 的细胞周期含义：control G2M rho=%s、E2F/S-phase rho=%s。",
        "这些是按 pseudotime 汇总的 module trend，不应写成 raw single-cell gene dynamics。"
      ), fmt(g2m, 3), fmt(e2f, 3)
    ),
    sprintf(
      paste0(
        "C5 使用 targeted 75-gene pseudobulk，而不是 genome-wide 数据。三个 DE contrast 的 padj<0.05 基因数分别为 %s/%s/%s；",
        "module enrichment 仅 E2F/S-phase/DNA-replication 达 q=%s，p53 module q=%s（未过 0.05）。"
      ),
      ifelse(is.na(de_counts[["CellCycle_treated_vs_control"]]), "NE", de_counts[["CellCycle_treated_vs_control"]]),
      ifelse(is.na(de_counts[["NonCellCycle_treated_vs_control"]]), "NE", de_counts[["NonCellCycle_treated_vs_control"]]),
      ifelse(is.na(de_counts[["CellCycle_highTGI_vs_lowTGI"]]), "NE", de_counts[["CellCycle_highTGI_vs_lowTGI"]]),
      fmt_p(e2f_q), fmt_p(p53_q)
    ),
    sprintf(
      paste0(
        "C6 因 MiloR 不可用，采用 pseudotime-bin fallback；treatment DA 中 q<0.05 的 bin 数=%s。",
        "TGI DA 的 CellCycle bin 4 只有 n=5、rho=-1 并产生极小 q，这是小样本/退化相关，不作为稳健发现或主结论。"
      ), ifelse(is.na(treatment_da_sig), "NE", treatment_da_sig)
    ),
    "",
    "## 11. 可支持与不可支持的表述",
    "",
    "### 可支持",
    "",
    "- Gemcitabine 与 CellCycle tumor-cell pseudotime distribution redistribution 相关，NonCellCycle 未显示同样总体 shift。",
    "- CellCycle shift–AUC-TGI 为正方向的样本级探索性关联；其强度和显著性受阈值与小样本影响。",
    "- 未检测到 30 与 120 mg/kg AUC-TGI 的清晰差异。",
    "- ETP group 结果可作为 endpoint-defined 描述性分层与敏感性分析。",
    "",
    "### 不可支持",
    "",
    "- 不可把 ETP-higher/lower 称为 baseline predictive biomarker。",
    "- 不可把按 ETP 调整后的关联解释为 ploidy-independent causal effect。",
    "- 不可把不显著写成等效；只有预设 equivalence CI 完全落在界值内才可宣称等效。",
    "- 不可用 pooled single-cell 相关代替样本级推断。",
    if (is_v42) {
      "- 不可把 ETP_boundary_stress_threshold_2_375 的 NE 写成“没有差异”；120 mg/kg endpoint-high 缺失导致相关比较不可估。"
    } else if (is_v43) {
      "- ETP_reference_balanced_threshold_2_24 改善了 untreated reference 的 initial-source 构成，但 30 mg/kg ETP-lower 与 120 mg/kg ETP-higher 仍各只有一个、且均来自 initial 4N。"
    } else {
      "- ETP_fixed_threshold_2_25 虽有全剂量两组覆盖，但单样本 strata 使交互与分组结论仍不稳定。"
    },
    "",
    "## 12. 最终科学解释",
    "",
    paste0(
      "最稳健的生物学信息是：gemcitabine 暴露与 cell-cycle-associated tumor states 的重新分布一致；",
      "这些状态变化在部分分析中与更强的整体 tumor-growth inhibition 同向。Endpoint ploidy 本身没有形成稳定、阈值无关的二分响应亚型。",
      if (is_v42) {
        "ETP_boundary_stress_threshold_2_375 的高阈值把 endpoint-high 压缩为 2 个样本并造成 120 mg/kg 缺组，主要价值是暴露模型对分组阈值和可估性的脆弱性。"
      } else if (is_v43) {
        "ETP_reference_balanced_threshold_2_24 使两个 untreated ETP reference 均含 initial 2N 和 initial 4N，因而更适合作为 endpoint-defined reference-sensitivity 版本；但 Dose 内嵌套与处理后分层的因果边界仍存在。"
      } else {
        "ETP_fixed_threshold_2_25 保留了每个剂量的两组覆盖，适合作为主要 endpoint-defined 描述性版本，但不能跨越处理后分层的因果边界。"
      },
      " 下一步优先采用样本级连续 endpoint-ploidy 模型、增加各 dose×ETP 区间样本，并在独立队列中预先固定阈值。"
    ),
    "",
    "## 数据来源",
    "",
    sprintf("- `%s/analysis_summary.json`", roots[[version]]),
    sprintf("- `%s/sample_end_timepoint_ploidy_assignments.csv`", roots[[version]]),
    sprintf("- `%s/ETP_group_dose_coverage_audit.csv`", roots[[version]]),
    sprintf("- `%s/C_SupportingAnalyses/`", roots[[version]]),
    ""
  )
}

# -----------------------------------------------------------------------------
# Slide-by-slide interpretation files
# -----------------------------------------------------------------------------

slide_records <- function(version) {
  is_v42 <- version == "ETP_boundary_stress_threshold_2_375"
  threshold <- thresholds[[version]]
  records <- list(
    list(
      title = "The questions lead the analysis",
      en = sprintf("The deck separates pseudotime redistribution, tumor response, and supporting evidence after all 16 tumors are grouped by endpoint mean ploidy at %.3f.", threshold),
      zh = sprintf("本页将问题分为 pseudotime redistribution、tumor response 与 supporting evidence；全部 16 个肿瘤按终点平均 ploidy 的 %.3f 阈值分组。", threshold),
      caveat_en = "Endpoint ploidy is measured after treatment, so this is descriptive endpoint stratification, not a baseline predictive-biomarker design.",
      caveat_zh = "终点 ploidy 是处理后变量，因此这是描述性终点分层，不是基线预测标志物设计。"
    ),
    list(
      title = "Main result and endpoint-ploidy distribution",
      en = sprintf("CellCycle 0 versus treated ECDF p=%s and ETP-stratified p=%s; the violin plot audits dose distributions inside the global endpoint groups.", fmt_p(value(version, "cellcycle_treated_vs_control_p")), fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p"))),
      zh = sprintf("CellCycle 0-vs-treated ECDF p=%s，ETP-stratified p=%s；violin 图核查全局 endpoint group 内各 dose 的 ploidy 分布。", fmt_p(value(version, "cellcycle_treated_vs_control_p")), fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p"))),
      caveat_en = if (is_v42) "Only two tumors are ETP-higher and 120 mg/kg has none; the plot exposes sparse coverage rather than proving a response subtype." else "All doses contain both groups, but two dose-by-group cells contain one tumor; the plot does not establish a response subtype.",
      caveat_zh = if (is_v42) "ETP-higher 仅 2 个肿瘤且 120 mg/kg 完全缺组；该图揭示稀疏覆盖，而不是证明响应亚型。" else "所有 dose 都有两组，但两个 dose-by-group 单元仅 1 个肿瘤；该图不能建立响应亚型。"
    ),
    list(
      title = "Data and endpoint-group coverage",
      en = sprintf("Sixteen tumors span 0/30/120 mg/kg. Coverage is %s.", coverage_sentence(version)),
      zh = sprintf("16 个肿瘤覆盖 0/30/120 mg/kg。分组覆盖为：%s。", coverage_sentence(version)),
      caveat_en = if (is_v42) "There is no endpoint-high sample at 120 mg/kg; models requiring that cell are not estimable." else "All cells exist, but 30-lower and 120-higher each contain one sample.",
      caveat_zh = if (is_v42) "120 mg/kg 没有 ETP-higher，要求该单元的模型不可估。" else "所有单元都存在，但 30-lower 与 120-higher 各仅 1 个样本。"
    ),
    list(
      title = "TGI definitions remain upstream of the ETP regrouping",
      en = "Endpoint and AUC-based TGI retain their original tumor-growth definitions. ETP does not recompute TGI from ETP groups; it changes the pseudotime grouping and untreated ECDF reference to sample-level ETP matching.",
      zh = "endpoint 与 AUC-based TGI 保留原始 tumor-growth 定义。ETP 不按 ETP group 重新计算 TGI；改变的是 pseudotime 分组以及按样本级 ETP 匹配的 untreated ECDF reference。",
      caveat_en = "The TGI formulas and the ECDF-shift matching rule use different grouping roles and should not be described as the same recalculation.",
      caveat_zh = "TGI 公式与 ECDF-shift matching rule 的分组角色不同，不能写成同一种重新计算。"
    ),
    list(
      title = "Each sample is compared with endpoint-group-matched controls",
      en = "Sample ECDFs are compared with an equal-weighted untreated ECDF reference from the same endpoint group; untreated samples use leave-one-control-out references.",
      zh = "每个样本 ECDF 与同 ETP group 的等权 untreated reference 比较；untreated 样本采用 leave-one-control-out reference。",
      caveat_en = "Changing the threshold changes group membership, reference ECDFs, individual shift metrics, and downstream TGI correlations.",
      caveat_zh = "阈值改变会同时改变组别、reference ECDF、个体 shift metric 与下游 TGI 相关。"
    ),
    list(
      title = "CellCycle shows treatment-associated redistribution",
      en = sprintf("CellCycle 0 versus treated ECDF p=%s; endpoint-group-stratified p=%s.", fmt_p(value(version, "cellcycle_treated_vs_control_p")), fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p"))),
      zh = sprintf("CellCycle 0-vs-treated ECDF p=%s；按 ETP group 分层后 p=%s。", fmt_p(value(version, "cellcycle_treated_vs_control_p")), fmt_p(value(version, "cellcycle_treated_vs_control_ploidy_stratified_p"))),
      caveat_en = "This is the strongest sample-level treatment result; it does not by itself establish a dose gradient.",
      caveat_zh = "这是最强的样本级治疗结果，但其本身不证明剂量梯度。"
    ),
    list(
      title = "NonCellCycle is the comparison compartment",
      en = sprintf("NonCellCycle 0 versus treated p=%s, providing no evidence of the same global redistribution.", fmt_p(value(version, "noncellcycle_treated_vs_control_p"))),
      zh = sprintf("NonCellCycle 0-vs-treated p=%s，没有同样的总体 redistribution 证据。", fmt_p(value(version, "noncellcycle_treated_vs_control_p"))),
      caveat_en = "A null comparison-compartment test supports contrast, but formal specificity also requires the direct compartment test.",
      caveat_zh = "comparison compartment 的阴性结果支持对照，但正式特异性仍需直接 compartment 检验。"
    ),
    list(
      title = "Endpoint-group treatment-effect modification is threshold sensitive",
      en = sprintf("CellCycle DID RMSE=%s (p=%s); NonCellCycle DID p=%s.", fmt(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse"), 4), fmt_p(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value(version, "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))),
      zh = sprintf("CellCycle DID RMSE=%s（p=%s）；NonCellCycle DID p=%s。", fmt(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse"), 4), fmt_p(value(version, "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value(version, "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))),
      caveat_en = if (is_v42) "Both compartments are nominally significant, CIs are unavailable, and endpoint-high is extremely sparse; do not claim CellCycle-specific modification." else "The null result does not prove equivalence, and strata remain small.",
      caveat_zh = if (is_v42) "两个 compartment 均名义显著、CI 不可得且 endpoint-high 极稀疏；不可宣称 CellCycle 特异修饰。" else "阴性结果不证明等效，而且 strata 仍很小。"
    ),
    list(
      title = "AUC-TGI summarizes the whole growth curve",
      en = sprintf("CellCycle shift versus AUC-TGI r=%s, permutation p=%s; adjusted shift p=%s.", fmt(value(version, "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value(version, "cellcycle_TGI_AUC_pearson_perm_p")), fmt_p(value(version, "cellcycle_AUC_TGI_adjusted_shift_p"))),
      zh = sprintf("CellCycle shift vs AUC-TGI r=%s、permutation p=%s；adjusted shift p=%s。", fmt(value(version, "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value(version, "cellcycle_TGI_AUC_pearson_perm_p")), fmt_p(value(version, "cellcycle_AUC_TGI_adjusted_shift_p"))),
      caveat_en = "Only eight treated tumors contribute; ETP adjustment is post-treatment conditioning, not confounder control in a causal sense.",
      caveat_zh = "仅 8 个 treated 肿瘤参与；ETP 调整是处理后条件化，不是因果意义的混杂控制。"
    ),
    list(
      title = "CellCycle dose-specific response remains exploratory",
      en = sprintf("AUC-TGI 120-minus-30=%s (p=%s); within-dose r values are %s at 30 and %s at 120 mg/kg.", fmt(value(version, "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 2), fmt_p(value(version, "cellcycle_TGI_AUC_30_vs_120_perm_p")), fmt(value(version, "cellcycle_TGI_AUC_shift_30mg_pearson_r"), 3), fmt(value(version, "cellcycle_TGI_AUC_shift_120mg_pearson_r"), 3)),
      zh = sprintf("AUC-TGI 120-minus-30=%s（p=%s）；dose 内 r：30 mg/kg=%s、120 mg/kg=%s。", fmt(value(version, "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 2), fmt_p(value(version, "cellcycle_TGI_AUC_30_vs_120_perm_p")), fmt(value(version, "cellcycle_TGI_AUC_shift_30mg_pearson_r"), 3), fmt(value(version, "cellcycle_TGI_AUC_shift_120mg_pearson_r"), 3)),
      caveat_en = "Each treated-dose stratum contains four samples.", caveat_zh = "每个 treated-dose stratum 只有 4 个样本。"
    ),
    list(
      title = "NonCellCycle does not mirror the CellCycle response alignment",
      en = sprintf("NonCellCycle overall shift-TGI r=%s (p=%s); within-dose r values are %s and %s.", fmt(value(version, "noncellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value(version, "noncellcycle_TGI_AUC_pearson_perm_p")), fmt(value(version, "noncellcycle_TGI_AUC_shift_30mg_pearson_r"), 3), fmt(value(version, "noncellcycle_TGI_AUC_shift_120mg_pearson_r"), 3)),
      zh = sprintf("NonCellCycle 总体 shift-TGI r=%s（p=%s）；dose 内 r=%s 与 %s。", fmt(value(version, "noncellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value(version, "noncellcycle_TGI_AUC_pearson_perm_p")), fmt(value(version, "noncellcycle_TGI_AUC_shift_30mg_pearson_r"), 3), fmt(value(version, "noncellcycle_TGI_AUC_shift_120mg_pearson_r"), 3)),
      caveat_en = sprintf("The direct compartment delta-r test remains non-significant (p=%s).", fmt_p(value(version, "compartment_delta_r_perm_p"))),
      caveat_zh = sprintf("直接 compartment delta-r 检验仍不显著（p=%s）。", fmt_p(value(version, "compartment_delta_r_perm_p")))
    ),
    list(
      title = "Endpoint TGI trajectories and endpoint-group response",
      en = if (is_v42) {
        sprintf("Endpoint-high versus endpoint-low AUC-TGI is not estimable. The endpoint-TGI trajectory figure compares 30 and 120 mg/kg across days; Day-38 120-minus-30=%s (dose-label permutation p=%s).", fmt(endpoint_tgi_day38_value(version, "mean_diff_b_minus_a"), 2), fmt_p(endpoint_tgi_day38_value(version, "permutation_p_two_sided")))
      } else {
        sprintf("Endpoint-high minus endpoint-low AUC-TGI=%s (p=%s). The endpoint-TGI trajectory figure compares 30 and 120 mg/kg across days; Day-38 120-minus-30=%s (dose-label permutation p=%s).", fmt(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_mean_diff"), 2), fmt_p(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_perm_p")), fmt(endpoint_tgi_day38_value(version, "mean_diff_b_minus_a"), 2), fmt_p(endpoint_tgi_day38_value(version, "permutation_p_two_sided")))
      },
      zh = if (is_v42) {
        sprintf("endpoint-high vs endpoint-low 的 AUC-TGI 不可估。endpoint-TGI trajectory 图比较 30 与 120 mg/kg 的各日 TGI；Day-38 的 120-minus-30=%s（dose-label permutation p=%s）。", fmt(endpoint_tgi_day38_value(version, "mean_diff_b_minus_a"), 2), fmt_p(endpoint_tgi_day38_value(version, "permutation_p_two_sided")))
      } else {
        sprintf("endpoint-high minus endpoint-low 的 AUC-TGI=%s（p=%s）。endpoint-TGI trajectory 图比较 30 与 120 mg/kg 的各日 TGI；Day-38 的 120-minus-30=%s（dose-label permutation p=%s）。", fmt(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_mean_diff"), 2), fmt_p(value(version, "tgi_auc_ETP_higher_minus_ETP_lower_perm_p")), fmt(endpoint_tgi_day38_value(version, "mean_diff_b_minus_a"), 2), fmt_p(endpoint_tgi_day38_value(version, "permutation_p_two_sided")))
      },
      caveat_en = if (is_v42) "The treated endpoint-high stratum has one tumor and is absent at 120 mg/kg; NE is not evidence of no difference." else "All-day endpoint-TGI trajectories and the AUC group comparison are descriptive; neither establishes a dose or endpoint-group response subtype.",
      caveat_zh = if (is_v42) "treated endpoint-high 仅 1 个肿瘤，且 120 mg/kg 缺组；NE 不代表没有差异。" else "各日 endpoint-TGI trajectory 与 AUC group comparison 均为描述性结果；不能建立 dose 或 endpoint-group 的响应亚型。"
    ),
    list(
      title = "Longitudinal growth models depend on estimable group cells",
      en = sprintf("Spline time-by-dose-by-group p=%s; Day-38 group-DID p values are %s at 30 and %s at 120 mg/kg.", fmt_p(value(version, "growth_curve_spline_time_dose_ploidy_lrt_p")), fmt_p(value(version, "growth_curve_day_final_30mg_ploidy_did_p")), fmt_p(value(version, "growth_curve_day_final_120mg_ploidy_did_p"))),
      zh = sprintf("spline time-by-dose-by-group p=%s；Day-38 group-DID p：30=%s、120=%s。", fmt_p(value(version, "growth_curve_spline_time_dose_ploidy_lrt_p")), fmt_p(value(version, "growth_curve_day_final_30mg_ploidy_did_p")), fmt_p(value(version, "growth_curve_day_final_120mg_ploidy_did_p"))),
      caveat_en = "Endpoint-defined growth curves are descriptive, and ETP_boundary_stress_threshold_2_375 is rank deficient at 120 mg/kg.",
      caveat_zh = "按终点分组的 growth curves 是描述性的；ETP_boundary_stress_threshold_2_375 在 120 mg/kg 处秩亏。"
    ),
    list(
      title = "Robustness preserves direction but not precision",
      en = sprintf("Leave-one-out minimum r=%s; bootstrap 95%% interval=%s.", fmt(robust_value(version, "leave_one_out_min_pearson_r"), 3), fmt_ci(robust_value(version, "bootstrap_pearson_r_ci025"), robust_value(version, "bootstrap_pearson_r_ci975"), 3)),
      zh = sprintf("leave-one-out 最小 r=%s；bootstrap 95%% interval=%s。", fmt(robust_value(version, "leave_one_out_min_pearson_r"), 3), fmt_ci(robust_value(version, "bootstrap_pearson_r_ci025"), robust_value(version, "bootstrap_pearson_r_ci975"), 3)),
      caveat_en = "Wide intervals and one influential flag show small-n fragility.", caveat_zh = "宽区间与一个 influential flag 显示小样本脆弱性。"
    ),
    list(
      title = "Endpoint ploidy and dose do not detectably explain the shift",
      en = sprintf(
        "Treated CellCycle shift versus mean endpoint ploidy p=%s; 30-versus-120 ECDF-shift t-test p=%s. The right-hand within-dose-centered shift–TGI plot has r=%s (p=%s).",
        fmt_p(endpoint_ploidy_shift_p(version)), fmt_p(dose_shift_p(version)),
        fmt(within_dose_centered_value(version, "estimate"), 3), fmt_p(within_dose_centered_value(version, "p_value"))
      ),
      zh = sprintf(
        "treated CellCycle shift 与平均终点 ploidy 的相关 p=%s；30-vs-120 ECDF-shift t-test p=%s。右图的 within-dose-centered shift–TGI 关系为 r=%s（p=%s）。",
        fmt_p(endpoint_ploidy_shift_p(version)), fmt_p(dose_shift_p(version)),
        fmt(within_dose_centered_value(version, "estimate"), 3), fmt_p(within_dose_centered_value(version, "p_value"))
      ),
      caveat_en = "Only eight treated tumors contribute. The figure is an exploratory within-dose-centering sensitivity analysis.",
      caveat_zh = "仅 8 个 treated 肿瘤参与。该图为探索性的 within-dose centering 敏感性分析。"
    ),
    list(
      title = "CellCycle endpoint ploidy is weakly related to sample endpoints",
      en = sprintf("CellCycle mean/median/max endpoint ploidy versus treated AUC-TGI r=%s/%s/%s; mean ploidy versus treated mean pseudotime r=%s.", fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3), fmt(c7$CellCycle$mean_pt_treated, 3)),
      zh = sprintf("CellCycle mean/median/max endpoint ploidy vs treated AUC-TGI r=%s/%s/%s；mean ploidy vs treated mean pseudotime r=%s。", fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3), fmt(c7$CellCycle$mean_pt_treated, 3)),
      caveat_en = "These are exploratory sample-level scatter plots and endpoint ploidy is post-treatment.",
      caveat_zh = "这些是探索性样本级散点图，而且 endpoint ploidy 是处理后变量。"
    ),
    list(
      title = "NonCellCycle endpoint ploidy is weakly related to sample endpoints",
      en = sprintf("NonCellCycle mean/median/max endpoint ploidy versus treated AUC-TGI r=%s/%s/%s; mean ploidy versus treated mean pseudotime r=%s.", fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3), fmt(c7$NonCellCycle$mean_pt_treated, 3)),
      zh = sprintf("NonCellCycle mean/median/max endpoint ploidy vs treated AUC-TGI r=%s/%s/%s；mean ploidy vs treated mean pseudotime r=%s。", fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3), fmt(c7$NonCellCycle$mean_pt_treated, 3)),
      caveat_en = "These are exploratory sample-level scatter plots and do not show a strong monotonic relationship.",
      caveat_zh = "这些是探索性样本级散点图，没有显示强的单调关系。"
    ),
    list(
      title = "Single-cell ploidy-pseudotime correlations are near zero",
      en = sprintf("Pooled single-cell r is %s for CellCycle (dose 0/30/120: %s/%s/%s) and %s for NonCellCycle (dose 0/30/120: %s/%s/%s).", fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$CellCycle$cell_by_dose[["0"]], 3), fmt(c7$CellCycle$cell_by_dose[["30"]], 3), fmt(c7$CellCycle$cell_by_dose[["120"]], 3), fmt(c7$NonCellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_by_dose[["0"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["30"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["120"]], 3)),
      zh = sprintf("pooled single-cell r：CellCycle=%s（dose 0/30/120：%s/%s/%s），NonCellCycle=%s（dose 0/30/120：%s/%s/%s）。", fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$CellCycle$cell_by_dose[["0"]], 3), fmt(c7$CellCycle$cell_by_dose[["30"]], 3), fmt(c7$CellCycle$cell_by_dose[["120"]], 3), fmt(c7$NonCellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_by_dose[["0"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["30"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["120"]], 3)),
      caveat_en = "Cells are nested within tumors; these plots are descriptive and pooled cell-level correlations cannot support sample-level inference.",
      caveat_zh = "细胞嵌套于肿瘤样本；这些图仅作描述，pooled cell-level 相关不能支持样本级推断。"
    ),
    list(
      title = "CellCycle versus NonCellCycle requires a direct test",
      en = sprintf("Compartment delta r=%s with label-swap p=%s.", fmt(value(version, "compartment_delta_r"), 3), fmt_p(value(version, "compartment_delta_r_perm_p"))),
      zh = sprintf("compartment delta r=%s，label-swap p=%s。", fmt(value(version, "compartment_delta_r"), 3), fmt_p(value(version, "compartment_delta_r_perm_p"))),
      caveat_en = "CellCycle specificity is compelling for the treatment distribution shift but not statistically established for the TGI correlation difference.",
      caveat_zh = "CellCycle 对 treatment distribution shift 的特异性较强，但 TGI correlation difference 尚未统计学确立。"
    ),
    list(
      title = "Gene programs support the pseudotime orientation",
      en = sprintf("Aggregated control trends show G2M rho=%s and E2F/S-phase rho=%s; targeted pseudobulk enrichment gives E2F module q=0.0198.", fmt(module_estimate(version, "G2M_checkpoint", "score_vs_pseudotime_control"), 3), fmt(module_estimate(version, "E2F_S_phase_DNA_replication", "score_vs_pseudotime_control"), 3)),
      zh = sprintf("聚合 control trend：G2M rho=%s、E2F/S-phase rho=%s；targeted pseudobulk 的 E2F module q=0.0198。", fmt(module_estimate(version, "G2M_checkpoint", "score_vs_pseudotime_control"), 3), fmt(module_estimate(version, "E2F_S_phase_DNA_replication", "score_vs_pseudotime_control"), 3)),
      caveat_en = "The 75-gene panel and pseudotime aggregation do not establish raw single-cell gene dynamics.",
      caveat_zh = "75-gene panel 与 pseudotime aggregation 不能证明 raw single-cell gene dynamics。"
    )
  )
  stopifnot(length(records) == 20L)
  records
}

slide_interpretation_en <- function(version) {
  recs <- slide_records(version)
  lines <- c(
    sprintf("# Slide-by-slide interpretation: %s endpoint-ploidy pseudotime/TGI analysis", version),
    "",
    sprintf("This document interprets the 20-slide %s deck and is synchronized with the formal result tables.", version),
    "",
    "Global boundary: ETP endpoint ploidy is post-treatment. Endpoint-group matching changes the untreated ECDF reference and therefore the individual shift metrics. The deck describes endpoint-defined associations and sensitivity, not a baseline predictive biomarker.",
    ""
  )
  for (i in seq_along(recs)) {
    r <- recs[[i]]
    lines <- c(
      lines,
      sprintf("## Slide %d. %s", i, r$title),
      "",
      sprintf("**Interpretation.** %s", r$en),
      "",
      sprintf("**Guardrail.** %s", r$caveat_en),
      "",
      sprintf("**Data basis.** `%s/analysis_summary.json` and the corresponding A/B/C module tables and figures.", roots[[version]]),
      ""
    )
  }
  lines
}

slide_interpretation_bilingual <- function(version) {
  recs <- slide_records(version)
  lines <- c(
    sprintf("# %s slide-by-slide interpretation / 逐页解读", version),
    "",
    "This file pairs concise English narration with a Chinese scientific interpretation.",
    "本文件将英文讲述要点与中文科学解读逐页对应。",
    "",
    "Global guardrail / 总体边界: endpoint ploidy is post-treatment; ETP is endpoint-defined descriptive/sensitivity analysis, not baseline predictive-biomarker evidence. / 终点 ploidy 是处理后变量；ETP 是终点定义的描述性/敏感性分析，不是基线预测标志物证据。",
    ""
  )
  for (i in seq_along(recs)) {
    r <- recs[[i]]
    lines <- c(
      lines,
      sprintf("## Slide %d. %s", i, r$title),
      "",
      sprintf("**English interpretation:** %s", r$en),
      "",
      sprintf("**中文解读：** %s", r$zh),
      "",
      sprintf("**English guardrail:** %s", r$caveat_en),
      "",
      sprintf("**中文边界：** %s", r$caveat_zh),
      ""
    )
  }
  lines
}

for (version in c("ETP_fixed_threshold_2_25", "ETP_boundary_stress_threshold_2_375", "ETP_reference_balanced_threshold_2_24")) {
  manuscript <- file.path(roots[[version]], "Manuscript")
  writeLines(comprehensive_report(version), file.path(manuscript, "comprehensive_results_interpretation.md"), useBytes = TRUE)
  writeLines(slide_interpretation_en(version), file.path(manuscript, "slide_by_slide_interpretation.md"), useBytes = TRUE)
  writeLines(slide_interpretation_bilingual(version), file.path(manuscript, "slide_by_slide_interpretation_bilingual.md"), useBytes = TRUE)
}

# -----------------------------------------------------------------------------
# Pairwise and integrated comparison reports
# -----------------------------------------------------------------------------

changed_samples <- crosswalk[crosswalk$`ETP_fixed_threshold_2_25_to_ETP_boundary_stress_threshold_2_375_change` != "unchanged", ]

v41_vs_v42 <- c(
  "# ETP_fixed_threshold_2_25 versus ETP_boundary_stress_threshold_2_375: endpoint-ploidy threshold sensitivity",
  "",
  "## Design difference",
  "",
  sprintf("Both versions group all 16 samples by endpoint mean ploidy. ETP_fixed_threshold_2_25 uses %.3f; ETP_boundary_stress_threshold_2_375 uses %.3f.", thresholds[["ETP_fixed_threshold_2_25"]], thresholds[["ETP_boundary_stress_threshold_2_375"]]),
  sprintf("Raising the threshold reclassifies %d samples; all move from ETP-higher in ETP_fixed_threshold_2_25 to ETP-lower in ETP_boundary_stress_threshold_2_375.", nrow(changed_samples)),
  sprintf("ETP_fixed_threshold_2_25 coverage: %s.", coverage_sentence("ETP_fixed_threshold_2_25")),
  sprintf("ETP_boundary_stress_threshold_2_375 coverage: %s.", coverage_sentence("ETP_boundary_stress_threshold_2_375")),
  "ETP_boundary_stress_threshold_2_375 therefore has only two ETP-higher samples overall and none at 120 mg/kg. This is a structural estimability failure, not a negative biological result.",
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Scientific interpretation",
  "",
  sprintf(
    paste0(
      "The core unstratified treatment result is identical: CellCycle p=%s and NonCellCycle p=%s. ",
      "The shift-TGI association weakens from r=%s (p=%s) in ETP_fixed_threshold_2_25 to r=%s (p=%s) in ETP_boundary_stress_threshold_2_375."
    ),
    fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_treated_vs_control_p")), fmt_p(value("ETP_fixed_threshold_2_25", "noncellcycle_treated_vs_control_p")),
    fmt(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_perm_p"))
  ),
  paste0(
    "ETP_boundary_stress_threshold_2_375 produces nominal A2 DID p values in both CellCycle and NonCellCycle, but the endpoint-high group is extremely sparse, bootstrap CIs are unavailable, ",
    "and group-response/mixed-model outputs become NE. This pattern is evidence of threshold sensitivity and rank deficiency, not a reproducible CellCycle-specific effect modifier."
  ),
  paste0(
    "Crucially, changing the threshold also changes which untreated samples form the matched ECDF reference. Individual ECDF shifts and shift-TGI correlations therefore change by construction. ",
    "The difference between ETP_fixed_threshold_2_25 and ETP_boundary_stress_threshold_2_375 cannot be interpreted as pure biological heterogeneity."
  ),
  "",
  "## Continuous-analysis comparison",
  "",
  sprintf(
    "Continuous C7 outputs are identical across ETP_fixed_threshold_2_25 and ETP_boundary_stress_threshold_2_375 because they do not use the binary threshold. CellCycle treated ploidy-TGI r=%s/%s/%s and NonCellCycle=%s/%s/%s for mean/median/max endpoint ploidy.",
    fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
    fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3)
  ),
  sprintf(
    "Single-cell pooled r is %s for CellCycle and %s for NonCellCycle; these are descriptive because cells are nested within tumors.",
    fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_pooled, 3)
  ),
  "",
  "## Decision",
  "",
  "- Use ETP_fixed_threshold_2_25 as the more interpretable endpoint-defined binary description because every dose has both groups, while stating the single-sample strata.",
  "- Treat ETP_boundary_stress_threshold_2_375 as a stress test showing how the binary analysis collapses at a higher threshold.",
  "- Prefer continuous sample-level endpoint-ploidy models for follow-up; the C7 continuous data are identical across ETP_fixed_threshold_2_25/ETP_boundary_stress_threshold_2_375.",
  "- Neither version supports a baseline predictive-biomarker or causal ploidy claim.",
  "",
  "## Reclassified samples",
  "",
  "See `sample_group_crosswalk.csv`; five samples change from ETP_fixed_threshold_2_25 higher to ETP_boundary_stress_threshold_2_375 lower.",
  ""
)
writeLines(v41_vs_v42, file.path(comparison_root, "ETP_fixed_threshold_2_25_vs_ETP_boundary_stress_threshold_2_375_interpretation.md"), useBytes = TRUE)

v41_vs_v2 <- c(
  "# ETP_fixed_threshold_2_25 versus initial_ploidy: endpoint-defined versus baseline ploidy",
  "",
  "## Estimand change",
  "",
  "initial_ploidy stratifies by initial 2N/4N ploidy measured before treatment and has balanced dose-by-group coverage (4/4 controls, 2/2 at 30, 2/2 at 120).",
  sprintf("ETP_fixed_threshold_2_25 stratifies all samples by post-treatment endpoint mean ploidy at %.3f: %s.", thresholds[["ETP_fixed_threshold_2_25"]], coverage_sentence("ETP_fixed_threshold_2_25")),
  paste0(
    "This is not a relabeling-only comparison. initial_ploidy constructs initial-ploidy-matched untreated ECDF references; ETP_fixed_threshold_2_25 constructs ETP-group-matched references. ",
    "Thus the sample-level shift metric and all shift-TGI analyses change. ETP is downstream of treatment, so ETP_fixed_threshold_2_25 may condition on a mediator/collider/selection variable."
  ),
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Stable findings",
  "",
  sprintf("- CellCycle 0-vs-treated p remains %s; NonCellCycle remains %s.", fmt_p(value("initial_ploidy", "cellcycle_treated_vs_control_p")), fmt_p(value("initial_ploidy", "noncellcycle_treated_vs_control_p"))),
  sprintf("- No clear 30-vs-120 AUC-TGI difference in either version: difference=%s, p=%s.", fmt(value("initial_ploidy", "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_30_vs_120_perm_p"))),
  "- CellCycle fraction, cluster-10 redistribution, and aggregated gene-program orientation remain supportive and are not the basis of the grouping change.",
  "",
  "## Changed or non-comparable findings",
  "",
  sprintf(
    "- Shift-TGI changes from initial_ploidy r=%s, p=%s to ETP_fixed_threshold_2_25 r=%s, p=%s; adjusted shift p changes from %s to %s.",
    fmt(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt_p(value("initial_ploidy", "cellcycle_AUC_TGI_adjusted_shift_p")), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_AUC_TGI_adjusted_shift_p"))
  ),
  paste0(
    "  The p-value movement should not be called biological confirmation: the reference ECDF and adjustment variable changed, and ETP_fixed_threshold_2_25 conditions on a post-treatment endpoint."
  ),
  sprintf(
    "- A2 remains null in both definitions (initial_ploidy initial-ploidy DID p=%s; ETP_fixed_threshold_2_25 ETP-group DID p=%s), but these are different effect-modification questions.",
    fmt_p(value("initial_ploidy", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
  ),
  sprintf(
    "- The spline growth interaction changes from initial_ploidy p=%s to ETP_fixed_threshold_2_25 p=%s. This contrasts baseline-defined heterogeneity with endpoint-defined descriptive trajectories; it is not evidence that one grouping is biologically true and the other false.",
    fmt_p(value("initial_ploidy", "growth_curve_spline_time_dose_ploidy_lrt_p")), fmt_p(value("ETP_fixed_threshold_2_25", "growth_curve_spline_time_dose_ploidy_lrt_p"))
  ),
  "",
  "## Continuous-analysis comparability",
  "",
  sprintf(
    "ETP_fixed_threshold_2_25 CellCycle treated endpoint mean/median/max ploidy versus AUC-TGI r=%s/%s/%s; NonCellCycle=%s/%s/%s. Pooled single-cell r=%s/%s for CellCycle/NonCellCycle.",
    fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
    fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3),
    fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_pooled, 3)
  ),
  "The formal initial_ploidy result tree contains no corresponding C7 continuous-output figures/tables, so a same-definition numerical ETP_fixed_threshold_2_25-versus-initial_ploidy continuous comparison is unavailable.",
  "",
  "## Scientific conclusion",
  "",
  paste0(
    "initial_ploidy is the appropriate framework for baseline response heterogeneity/predictive hypotheses because initial ploidy precedes treatment. ",
    "ETP_fixed_threshold_2_25 answers a different, endpoint-defined descriptive question: tumors ending in different ploidy states can be compared for their observed state redistribution and response, but causal direction is unresolved. ",
    "Agreement on the main CellCycle redistribution strengthens that treatment-associated result; disagreement in secondary correlations and growth interactions defines model/threshold sensitivity rather than pure biology."
  ),
  ""
)
writeLines(v41_vs_v2, file.path(comparison_root, "ETP_fixed_threshold_2_25_vs_initial_ploidy_interpretation.md"), useBytes = TRUE)

v42_vs_v2 <- c(
  "# ETP_boundary_stress_threshold_2_375 versus initial_ploidy: baseline balance versus endpoint-group rank deficiency",
  "",
  "## Estimand and coverage",
  "",
  "initial_ploidy uses balanced pre-treatment initial 2N/4N groups and can estimate all dose-by-group contrasts.",
  sprintf("ETP_boundary_stress_threshold_2_375 uses endpoint mean ploidy > %.3f for ETP-higher: %s.", thresholds[["ETP_boundary_stress_threshold_2_375"]], coverage_sentence("ETP_boundary_stress_threshold_2_375")),
  "At 120 mg/kg all four samples are ETP-lower. Higher-vs-lower response contrasts, key interactions, and the longitudinal spline model are therefore NE.",
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Interpretation",
  "",
  sprintf("The stable core remains CellCycle treatment redistribution (p=%s) with no NonCellCycle redistribution (p=%s).", fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_treated_vs_control_p")), fmt_p(value("ETP_boundary_stress_threshold_2_375", "noncellcycle_treated_vs_control_p"))),
  sprintf(
    "CellCycle shift-TGI weakens from initial_ploidy r=%s (p=%s) to ETP_boundary_stress_threshold_2_375 r=%s (p=%s); the ETP_boundary_stress_threshold_2_375 bootstrap interval %s crosses zero.",
    fmt(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt_ci(robust_value("ETP_boundary_stress_threshold_2_375", "bootstrap_pearson_r_ci025"), robust_value("ETP_boundary_stress_threshold_2_375", "bootstrap_pearson_r_ci975"), 3)
  ),
  sprintf(
    paste0(
      "ETP_boundary_stress_threshold_2_375 A2 gives CellCycle p=%s and NonCellCycle p=%s, unlike the null initial_ploidy initial-ploidy DID. ",
      "Because endpoint-high contains only two samples, both compartments become nominal, CIs are unavailable, and downstream cells are missing, this is not credible evidence of a new CellCycle-specific biological subtype."
    ),
    fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_boundary_stress_threshold_2_375", "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
  ),
  paste0(
    "The reference construction also differs: initial_ploidy matches controls by pre-treatment initial ploidy, whereas ETP_boundary_stress_threshold_2_375 matches by a sparse post-treatment ETP group. ",
    "Consequently, initial_ploidy-to-ETP_boundary_stress_threshold_2_375 changes are a mixture of estimand, reference, threshold, and support differences and cannot be attributed solely to biological heterogeneity."
  ),
  "",
  "## Continuous-analysis comparability",
  "",
  sprintf(
    "ETP_boundary_stress_threshold_2_375 continuous C7 values equal ETP_fixed_threshold_2_25: CellCycle treated mean/median/max endpoint-ploidy versus AUC-TGI r=%s/%s/%s; pooled single-cell r=%s for CellCycle and %s for NonCellCycle.",
    fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
    fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_pooled, 3)
  ),
  "The formal initial_ploidy result tree contains no corresponding C7 continuous-output figures/tables, so a same-definition numerical ETP_boundary_stress_threshold_2_375-versus-initial_ploidy continuous comparison is unavailable.",
  "",
  "## Scientific conclusion",
  "",
  "- Retain initial_ploidy for baseline predictive/effect-modification questions.",
  "- Use ETP_boundary_stress_threshold_2_375 only as a boundary stress test demonstrating that the endpoint binary split becomes non-estimable at 2.375.",
  "- Do not interpret NE as no difference and do not use ETP_boundary_stress_threshold_2_375 nominal A2 p values as a primary discovery.",
  "- Use continuous sample-level endpoint ploidy for descriptive follow-up and collect additional endpoint-high samples, especially at 120 mg/kg.",
  ""
)
writeLines(v42_vs_v2, file.path(comparison_root, "ETP_boundary_stress_threshold_2_375_vs_initial_ploidy_interpretation.md"), useBytes = TRUE)

v43_vs_v41 <- c(
  "# ETP_reference_balanced_threshold_2_24 versus ETP_fixed_threshold_2_25: 2.24 reference-balanced threshold versus 2.25",
  "",
  "## Design and coverage",
  "",
  sprintf("ETP_reference_balanced_threshold_2_24 uses endpoint mean ploidy > %.3f: %s.", thresholds[["ETP_reference_balanced_threshold_2_24"]], coverage_sentence("ETP_reference_balanced_threshold_2_24")),
  sprintf("Its untreated reference composition is deliberately mixed by initial source: %s.", untreated_initial_composition("ETP_reference_balanced_threshold_2_24")),
  sprintf("ETP_fixed_threshold_2_25 uses endpoint mean ploidy > %.3f: %s.", thresholds[["ETP_fixed_threshold_2_25"]], coverage_sentence("ETP_fixed_threshold_2_25")),
  sprintf("Its untreated reference composition is %s; therefore ETP_fixed_threshold_2_25 ETP-higher is initial-4N-only.", untreated_initial_composition("ETP_fixed_threshold_2_25")),
  "ETP_reference_balanced_threshold_2_24 changes one sample assignment relative to ETP_fixed_threshold_2_25 only if its sample mean lies between 2.240 and 2.250, but it changes the group-matched untreated reference composition by design. Individual ECDF shifts can therefore change even when the binary labels are nearly identical.",
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Interpretation",
  "",
  sprintf("The unstratified CellCycle treatment result remains p=%s and the NonCellCycle comparison remains p=%s in both versions.", fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_treated_vs_control_p")), fmt_p(value("ETP_reference_balanced_threshold_2_24", "noncellcycle_treated_vs_control_p"))),
  sprintf("CellCycle treated-sample shift–AUC-TGI is r=%s (p=%s) in ETP_reference_balanced_threshold_2_24 versus r=%s (p=%s) in ETP_fixed_threshold_2_25. This is a sample-level association with Dose pooled in the raw correlation; it is not a group-mean correlation.", fmt(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_perm_p")), fmt(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_perm_p"))),
  sprintf("A2 CellCycle DID p is %s in ETP_reference_balanced_threshold_2_24 versus %s in ETP_fixed_threshold_2_25. Both remain endpoint-defined descriptive contrasts, not baseline effect-modification tests.", fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))),
  "",
  "## Continuous-analysis comparison",
  "",
  "All C7 continuous quantities are identical in ETP_fixed_threshold_2_25 and ETP_reference_balanced_threshold_2_24 because they use the same cells and do not use the binary ETP threshold. The threshold only changes binary labels, matched references, and threshold-drawn plot elements.",
  "",
  "## Decision",
  "",
  "ETP_reference_balanced_threshold_2_24 is preferable to ETP_fixed_threshold_2_25 for an endpoint-defined sensitivity analysis when reference composition is a concern, because neither untreated ETP reference is monopolized by one initial-ploidy source. It does not turn ETP into a baseline biomarker or remove Dose-specific initial-ploidy nesting.",
  ""
)
writeLines(v43_vs_v41, file.path(comparison_root, "ETP_reference_balanced_threshold_2_24_vs_ETP_fixed_threshold_2_25_interpretation.md"), useBytes = TRUE)

v43_vs_v42 <- c(
  "# ETP_reference_balanced_threshold_2_24 versus ETP_boundary_stress_threshold_2_375: reference-balanced 2.24 threshold versus sparse 2.375 threshold",
  "",
  "## Design and coverage",
  "",
  sprintf("ETP_reference_balanced_threshold_2_24: threshold %.3f; %s. Untreated references: %s.", thresholds[["ETP_reference_balanced_threshold_2_24"]], coverage_sentence("ETP_reference_balanced_threshold_2_24"), untreated_initial_composition("ETP_reference_balanced_threshold_2_24")),
  sprintf("ETP_boundary_stress_threshold_2_375: threshold %.3f; %s.", thresholds[["ETP_boundary_stress_threshold_2_375"]], coverage_sentence("ETP_boundary_stress_threshold_2_375")),
  "ETP_boundary_stress_threshold_2_375 has no 120 mg/kg ETP-higher tumor, whereas ETP_reference_balanced_threshold_2_24 has an ETP-higher tumor at every Dose. ETP_boundary_stress_threshold_2_375 ETP-higher untreated reference contains one initial-4N tumor only.",
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Interpretation",
  "",
  sprintf("CellCycle shift–AUC-TGI is r=%s (p=%s) in ETP_reference_balanced_threshold_2_24 and r=%s (p=%s) in ETP_boundary_stress_threshold_2_375. The change is expected to be reference-sensitive because the matched untreated ECDF changes with threshold and composition.", fmt(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_perm_p")), fmt(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_perm_p"))),
  sprintf("A2 CellCycle DID p is %s in ETP_reference_balanced_threshold_2_24 and %s in ETP_boundary_stress_threshold_2_375. The ETP_boundary_stress_threshold_2_375 result remains a sparse/rank-deficient stress test rather than evidence for a distinct biological subtype.", fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))),
  "",
  "## Decision",
  "",
  "Use ETP_reference_balanced_threshold_2_24 rather than ETP_boundary_stress_threshold_2_375 for interpretable endpoint-defined sensitivity analyses. Treat ETP_boundary_stress_threshold_2_375 as an upper-threshold boundary test showing loss of estimability, not as a competing primary result.",
  ""
)
writeLines(v43_vs_v42, file.path(comparison_root, "ETP_reference_balanced_threshold_2_24_vs_ETP_boundary_stress_threshold_2_375_interpretation.md"), useBytes = TRUE)

v43_vs_v2 <- c(
  "# ETP_reference_balanced_threshold_2_24 versus initial_ploidy: endpoint-defined reference-balanced sensitivity versus baseline initial ploidy",
  "",
  "## Estimand boundary",
  "",
  "initial_ploidy uses initial 2N/4N measured before treatment and is the appropriate analysis for baseline predictive/effect-modification questions.",
  sprintf("ETP_reference_balanced_threshold_2_24 uses post-treatment endpoint mean ploidy > %.3f: %s.", thresholds[["ETP_reference_balanced_threshold_2_24"]], coverage_sentence("ETP_reference_balanced_threshold_2_24")),
  sprintf("ETP_reference_balanced_threshold_2_24 untreated ETP references are mixed by initial source: %s.", untreated_initial_composition("ETP_reference_balanced_threshold_2_24")),
  "This improves the endpoint-matched reference composition but does not make ETP_reference_balanced_threshold_2_24 a baseline analysis: ETP remains downstream of gemcitabine treatment and is also used to construct the reference ECDF.",
  "",
  "## Numerical comparison",
  "",
  markdown_metric_table(),
  "",
  "## Interpretation",
  "",
  sprintf("CellCycle 0-vs-treated redistribution is the same (p=%s), while NonCellCycle remains null (p=%s).", fmt_p(value("initial_ploidy", "cellcycle_treated_vs_control_p")), fmt_p(value("initial_ploidy", "noncellcycle_treated_vs_control_p"))),
  sprintf("The treated-sample shift–AUC-TGI association is initial_ploidy r=%s (p=%s) versus ETP_reference_balanced_threshold_2_24 r=%s (p=%s). This numerical movement cannot be called a biological effect of endpoint ploidy because the matched controls, shift values, and adjustment variable differ.", fmt(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_perm_p")), fmt(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_perm_p"))),
  sprintf("The group DID p values are initial_ploidy=%s and ETP_reference_balanced_threshold_2_24=%s; they answer different questions (baseline initial ploidy versus endpoint state).", fmt_p(value("initial_ploidy", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")), fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))),
  "",
  "## Continuous-analysis comparability",
  "",
  "ETP_reference_balanced_threshold_2_24 carries the same C7 continuous endpoint-ploidy figures as the other ETP thresholds. initial_ploidy has no same-definition C7 output in its formal result tree, so direct numerical comparison requires rerunning the continuous module for initial_ploidy.",
  "",
  "## Decision",
  "",
  "Retain initial_ploidy for baseline biological heterogeneity or prediction. Use ETP_reference_balanced_threshold_2_24 only for a more reference-balanced descriptive endpoint-state sensitivity analysis.",
  ""
)
writeLines(v43_vs_v2, file.path(comparison_root, "ETP_reference_balanced_threshold_2_24_vs_initial_ploidy_interpretation.md"), useBytes = TRUE)

integrated <- c(
  "# Integrated scientific interpretation: initial_ploidy, ETP_fixed_threshold_2_25, and ETP_boundary_stress_threshold_2_375",
  "",
  "## Executive answer",
  "",
  paste0(
    "Across all three versions, the most reproducible result is that gemcitabine treatment is associated with a redistribution of CellCycle tumor-cell pseudotime states, ",
    "while NonCellCycle does not show the same overall shift. The positive association between CellCycle shift magnitude and AUC-TGI is directionally recurrent but small-n and sensitive to how controls are matched. ",
    "Binary endpoint-ploidy heterogeneity is not stable across 2.25 and 2.375 thresholds."
  ),
  "",
  markdown_metric_table(),
  "",
  "## What is invariant",
  "",
  sprintf("- CellCycle 0-vs-treated p=%s in initial_ploidy, ETP_fixed_threshold_2_25, and ETP_boundary_stress_threshold_2_375.", fmt_p(value("initial_ploidy", "cellcycle_treated_vs_control_p"))),
  sprintf("- NonCellCycle 0-vs-treated p=%s in all versions.", fmt_p(value("initial_ploidy", "noncellcycle_treated_vs_control_p"))),
  sprintf("- AUC-TGI 120-minus-30=%s with p=%s in all versions.", fmt(value("initial_ploidy", "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_30_vs_120_perm_p"))),
  sprintf("- CellCycle fraction treatment p values remain non-significant (%s/%s/%s).", fmt_p(fraction_p("initial_ploidy")), fmt_p(fraction_p("ETP_fixed_threshold_2_25")), fmt_p(fraction_p("ETP_boundary_stress_threshold_2_375"))),
  sprintf("- Cluster 10 shows treatment redistribution in all runs (p=%s/%s/%s), with suggestive rather than definitive TGI correlation.", fmt_p(cluster10_pt_p("initial_ploidy")), fmt_p(cluster10_pt_p("ETP_fixed_threshold_2_25")), fmt_p(cluster10_pt_p("ETP_boundary_stress_threshold_2_375"))),
  "",
  "## What is threshold/estimand sensitive",
  "",
  sprintf(
    "- CellCycle shift-TGI r/p: initial_ploidy %s/%s; ETP_fixed_threshold_2_25 %s/%s; ETP_boundary_stress_threshold_2_375 %s/%s.",
    fmt(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_perm_p"))
  ),
  sprintf(
    "- A2 CellCycle DID p: initial_ploidy=%s, ETP_fixed_threshold_2_25=%s, ETP_boundary_stress_threshold_2_375=%s; ETP_boundary_stress_threshold_2_375 NonCellCycle p=%s and lacks CIs.",
    fmt_p(value("initial_ploidy", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_boundary_stress_threshold_2_375", "noncellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
  ),
  sprintf(
    "- Spline growth interaction p: initial_ploidy=%s, ETP_fixed_threshold_2_25=%s, ETP_boundary_stress_threshold_2_375=NE.",
    fmt_p(value("initial_ploidy", "growth_curve_spline_time_dose_ploidy_lrt_p")),
    fmt_p(value("ETP_fixed_threshold_2_25", "growth_curve_spline_time_dose_ploidy_lrt_p"))
  ),
  "",
  "## Why initial_ploidy and ETP are not interchangeable",
  "",
  paste0(
    "initial_ploidy asks whether a pre-treatment characteristic (initial 2N/4N) modifies treatment-state and growth response. ",
    "ETP asks how tumors grouped by their post-treatment endpoint state differ descriptively. ETP also uses endpoint-group-matched untreated ECDF references, so changing the grouping changes the measured shift itself. ",
    "Conditioning on ETP can condition on treatment downstream and induce mediator/collider/selection bias. Therefore, ETP cannot validate a baseline predictive biomarker, and initial_ploidy-ETP differences cannot be interpreted as pure biological heterogeneity."
  ),
  "",
  "## Continuous endpoint-ploidy synthesis",
  "",
  sprintf(
    paste0(
      "ETP_fixed_threshold_2_25/ETP_boundary_stress_threshold_2_375 share the same continuous C7 data. CellCycle endpoint mean/median/max ploidy versus treated AUC-TGI r=%s/%s/%s; ",
      "NonCellCycle=%s/%s/%s. Mean endpoint ploidy versus treated mean pseudotime r=%s for CellCycle and %s for NonCellCycle."
    ),
    fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
    fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3),
    fmt(c7$CellCycle$mean_pt_treated, 3), fmt(c7$NonCellCycle$mean_pt_treated, 3)
  ),
  sprintf(
    paste0(
      "Single-cell ploidy-pseudotime pooled r is %s for CellCycle (dose 0/30/120: %s/%s/%s) and %s for NonCellCycle (dose 0/30/120: %s/%s/%s). ",
      "These correlations are descriptive and cannot be used for sample-level inference because cells within tumors are dependent."
    ),
    fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$CellCycle$cell_by_dose[["0"]], 3), fmt(c7$CellCycle$cell_by_dose[["30"]], 3), fmt(c7$CellCycle$cell_by_dose[["120"]], 3),
    fmt(c7$NonCellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_by_dose[["0"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["30"]], 3), fmt(c7$NonCellCycle$cell_by_dose[["120"]], 3)
  ),
  "The formal initial_ploidy result tree has no same-definition C7 continuous outputs; therefore, initial_ploidy cannot be compared numerically for these exploratory plots without rerunning the continuous module.",
  "",
  "## Optional-module synthesis",
  "",
  sprintf("- Aggregated control module trends orient pseudotime biologically: G2M rho=%s and E2F/S-phase rho=%s.", fmt(module_estimate("ETP_fixed_threshold_2_25", "G2M_checkpoint", "score_vs_pseudotime_control"), 3), fmt(module_estimate("ETP_fixed_threshold_2_25", "E2F_S_phase_DNA_replication", "score_vs_pseudotime_control"), 3)),
  "- The targeted 75-gene pseudobulk produces no padj<0.05 gene in any of the three DE contrasts; E2F/S-phase/DNA replication is the only q<0.05 module enrichment (q=0.0198), while p53 is borderline (q=0.0538).",
  "- MiloR was unavailable; pseudotime-bin fallback shows no treatment-DA bin at q<0.05. CellCycle TGI bin 4 (n=5, rho=-1, tiny q) is a degenerate small-n correlation and is excluded from the main conclusion.",
  "",
  "## Final hierarchy of evidence",
  "",
  "1. **Primary supported:** gemcitabine-associated CellCycle pseudotime redistribution under sample-level inference.",
  "2. **Suggestive:** larger CellCycle redistribution aligns with stronger AUC-TGI in direction, but uncertainty and reference sensitivity remain.",
  "3. **Not established:** a threshold-stable endpoint-ploidy response subtype, a statistically proven CellCycle-specific TGI correlation, or 30-vs-120 dose separation.",
  "4. **Not permissible:** baseline predictive-biomarker, causal mediation, or ploidy-independent causal-effect claims from ETP.",
  "",
  "## Recommended next analysis",
  "",
  "- Pre-specify and validate a continuous sample-level endpoint-ploidy model; include nonlinear terms only if supported by sample size.",
  "- Keep baseline initial ploidy for predictive/effect-modification questions and endpoint ploidy for descriptive/mediation-oriented questions with explicit causal assumptions.",
  "- Increase the number of endpoint-high samples within every dose, especially 120 mg/kg, before fitting interactions.",
  "- Validate the CellCycle shift-TGI relation in an independent cohort and report leave-one-out/bootstrap uncertainty.",
  "- If neighborhood DA is a priority, rerun with a prespecified graph method/MiloR environment and adequate per-neighborhood sample support.",
  ""
)
writeLines(integrated, file.path(comparison_root, "integrated_scientific_interpretation.md"), useBytes = TRUE)

integrated_v43 <- c(
  "# Integrated scientific interpretation: initial_ploidy, ETP_fixed_threshold_2_25, ETP_boundary_stress_threshold_2_375, and ETP_reference_balanced_threshold_2_24",
  "",
  "## Executive answer",
  "",
  paste0(
    "Across all four versions, gemcitabine is associated with a redistribution of CellCycle tumor-cell pseudotime states, while NonCellCycle does not show the same overall shift. ",
    "The positive treated-sample shift–AUC-TGI association is directionally recurrent but small-n and sensitive to the group-matched untreated reference. ",
    "ETP_reference_balanced_threshold_2_24 (threshold 2.24) was selected specifically to ensure that each untreated ETP reference contains both initial 2N and initial 4N tumors; it improves reference composition without converting endpoint ploidy into a baseline predictor."
  ),
  "",
  markdown_metric_table(),
  "",
  "## What is invariant",
  "",
  sprintf("- CellCycle 0-vs-treated p=%s in initial_ploidy, ETP_fixed_threshold_2_25, ETP_boundary_stress_threshold_2_375, and ETP_reference_balanced_threshold_2_24.", fmt_p(value("initial_ploidy", "cellcycle_treated_vs_control_p"))),
  sprintf("- NonCellCycle 0-vs-treated p=%s in all versions.", fmt_p(value("initial_ploidy", "noncellcycle_treated_vs_control_p"))),
  sprintf("- AUC-TGI 120-minus-30=%s with p=%s in all versions.", fmt(value("initial_ploidy", "cellcycle_TGI_AUC_120_minus_30_mean_diff"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_30_vs_120_perm_p"))),
  sprintf("- CellCycle fraction treatment p values remain non-significant (initial_ploidy/ETP_fixed_threshold_2_25/ETP_boundary_stress_threshold_2_375/ETP_reference_balanced_threshold_2_24: %s/%s/%s/%s).", fmt_p(fraction_p("initial_ploidy")), fmt_p(fraction_p("ETP_fixed_threshold_2_25")), fmt_p(fraction_p("ETP_boundary_stress_threshold_2_375")), fmt_p(fraction_p("ETP_reference_balanced_threshold_2_24"))),
  "",
  "## ETP_reference_balanced_threshold_2_24 reference-composition gain and remaining limitation",
  "",
  sprintf("- ETP_reference_balanced_threshold_2_24 untreated references: %s.", untreated_initial_composition("ETP_reference_balanced_threshold_2_24")),
  sprintf("- ETP_fixed_threshold_2_25 untreated references: %s; ETP_fixed_threshold_2_25 ETP-higher is therefore initial-4N-only.", untreated_initial_composition("ETP_fixed_threshold_2_25")),
  sprintf("- ETP_boundary_stress_threshold_2_375 untreated references: %s; its ETP-higher reference has one initial-4N sample only.", untreated_initial_composition("ETP_boundary_stress_threshold_2_375")),
  "- No scalar threshold can make every Dose x ETP group contain both initial sources: at 120 mg/kg, all initial-2N sample means are below all initial-4N sample means. ETP_reference_balanced_threshold_2_24 reduces the control-reference confounding but cannot remove that dose-specific nesting.",
  "",
  "## What is threshold- or estimand-sensitive",
  "",
  sprintf(
    "- CellCycle shift–AUC-TGI r/p: initial_ploidy %s/%s; ETP_fixed_threshold_2_25 %s/%s; ETP_boundary_stress_threshold_2_375 %s/%s; ETP_reference_balanced_threshold_2_24 %s/%s.",
    fmt(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("initial_ploidy", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_TGI_AUC_pearson_perm_p")),
    fmt(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_r"), 3), fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_TGI_AUC_pearson_perm_p"))
  ),
  sprintf(
    "- A2 CellCycle DID p: initial_ploidy=%s, ETP_fixed_threshold_2_25=%s, ETP_boundary_stress_threshold_2_375=%s, ETP_reference_balanced_threshold_2_24=%s. ETP_boundary_stress_threshold_2_375 remains sparse/rank-deficient; ETP_reference_balanced_threshold_2_24 is reference-balanced but remains endpoint-defined.",
    fmt_p(value("initial_ploidy", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_fixed_threshold_2_25", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_boundary_stress_threshold_2_375", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p")),
    fmt_p(value("ETP_reference_balanced_threshold_2_24", "cellcycle_pseudotime_ploidy_did_ecdf_rmse_perm_p"))
  ),
  sprintf(
    "- Spline growth interaction p: initial_ploidy=%s, ETP_fixed_threshold_2_25=%s, ETP_boundary_stress_threshold_2_375=NE, ETP_reference_balanced_threshold_2_24=%s.",
    fmt_p(value("initial_ploidy", "growth_curve_spline_time_dose_ploidy_lrt_p")),
    fmt_p(value("ETP_fixed_threshold_2_25", "growth_curve_spline_time_dose_ploidy_lrt_p")),
    fmt_p(value("ETP_reference_balanced_threshold_2_24", "growth_curve_spline_time_dose_ploidy_lrt_p"))
  ),
  "",
  "## Why initial_ploidy and ETP remain different scientific estimands",
  "",
  "initial_ploidy asks whether pre-treatment initial 2N/4N ploidy modifies treatment-state or growth response. ETP asks how tumors grouped by their post-treatment endpoint state differ descriptively. Because ETP uses ETP both to define strata and to select matched untreated ECDF references, it conditions on a treatment-downstream variable; ETP_reference_balanced_threshold_2_24 improves the initial-source mixture of that reference, but cannot establish a baseline predictive biomarker, a causal mediator, or a ploidy-independent causal effect.",
  "",
  "## Continuous endpoint-ploidy synthesis",
  "",
  sprintf(
    "ETP_fixed_threshold_2_25/ETP_boundary_stress_threshold_2_375/ETP_reference_balanced_threshold_2_24 share identical C7 continuous data. CellCycle endpoint mean/median/max ploidy versus treated AUC-TGI r=%s/%s/%s; NonCellCycle=%s/%s/%s. Mean endpoint ploidy versus treated mean pseudotime r=%s for CellCycle and %s for NonCellCycle.",
    fmt(c7$CellCycle$tgi_treated[["mean"]], 3), fmt(c7$CellCycle$tgi_treated[["median"]], 3), fmt(c7$CellCycle$tgi_treated[["max"]], 3),
    fmt(c7$NonCellCycle$tgi_treated[["mean"]], 3), fmt(c7$NonCellCycle$tgi_treated[["median"]], 3), fmt(c7$NonCellCycle$tgi_treated[["max"]], 3),
    fmt(c7$CellCycle$mean_pt_treated, 3), fmt(c7$NonCellCycle$mean_pt_treated, 3)
  ),
  sprintf("Single-cell pooled ploidy–pseudotime r is %s for CellCycle and %s for NonCellCycle. These descriptive cell-level correlations are near zero and cannot be used as sample-level inference because cells are nested within tumors.", fmt(c7$CellCycle$cell_pooled, 3), fmt(c7$NonCellCycle$cell_pooled, 3)),
  "initial_ploidy has no same-definition C7 output in its formal result tree; a direct initial_ploidy continuous comparison requires rerunning that module.",
  "",
  "## Final hierarchy of evidence",
  "",
  "1. **Primary supported:** gemcitabine-associated CellCycle pseudotime redistribution under sample-level inference.",
  "2. **Suggestive:** larger CellCycle redistribution aligns directionally with stronger AUC-TGI, but uncertainty and reference sensitivity remain.",
  "3. **Not established:** a threshold-stable endpoint-ploidy response subtype, a statistically proven CellCycle-specific TGI correlation, or 30-vs-120 dose separation.",
  "4. **Not permissible:** baseline predictive-biomarker, causal mediation, or ploidy-independent causal-effect claims from any ETP threshold.",
  ""
)
writeLines(integrated_v43, file.path(comparison_root, "integrated_scientific_interpretation.md"), useBytes = TRUE)

readme <- c(
  "# 04h cross-version comparison results",
  "",
  "This directory was generated from the formal initial_ploidy, ETP_fixed_threshold_2_25, ETP_boundary_stress_threshold_2_375, and ETP_reference_balanced_threshold_2_24 result artifacts by:",
  "",
  "`Code/in-vivo/04h_v4_results_interpretation_and_comparison.R`",
  "",
  "## Files",
  "",
  "- `sample_group_crosswalk.csv`: one row per tumor, showing baseline initial ploidy and all endpoint-threshold assignments.",
  "- `dose_group_coverage_comparison.csv`: dose-by-group sample counts and estimability flags for initial_ploidy/ETP_fixed_threshold_2_25/ETP_boundary_stress_threshold_2_375/ETP_reference_balanced_threshold_2_24.",
  "- `comparison_master_metrics.csv`: harmonized primary/supporting metrics with comparability notes.",
  "- `ETP_fixed_threshold_2_25_vs_ETP_boundary_stress_threshold_2_375_interpretation.md`: endpoint-threshold sensitivity comparison.",
  "- `ETP_fixed_threshold_2_25_vs_initial_ploidy_interpretation.md`: ETP_fixed_threshold_2_25 endpoint-defined versus initial_ploidy baseline-defined comparison.",
  "- `ETP_boundary_stress_threshold_2_375_vs_initial_ploidy_interpretation.md`: ETP_boundary_stress_threshold_2_375 endpoint-defined/rank-deficient versus initial_ploidy baseline-defined comparison.",
  "- `ETP_reference_balanced_threshold_2_24_vs_ETP_fixed_threshold_2_25_interpretation.md`: 2.24 reference-balanced threshold versus 2.25 comparison.",
  "- `ETP_reference_balanced_threshold_2_24_vs_ETP_boundary_stress_threshold_2_375_interpretation.md`: 2.24 reference-balanced threshold versus sparse 2.375 comparison.",
  "- `ETP_reference_balanced_threshold_2_24_vs_initial_ploidy_interpretation.md`: ETP_reference_balanced_threshold_2_24 endpoint-defined versus initial_ploidy baseline-defined comparison.",
  "- `integrated_scientific_interpretation.md`: evidence hierarchy and integrated scientific conclusion.",
  "",
  "Each ETP `Manuscript` directory also receives:",
  "",
  "- `comprehensive_results_interpretation.md`",
  "- `slide_by_slide_interpretation.md`",
  "- `slide_by_slide_interpretation_bilingual.md`",
  "",
  "## Reproducibility",
  "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("initial_ploidy source: `%s`", roots[["initial_ploidy"]]),
  sprintf("ETP_fixed_threshold_2_25 source: `%s`", roots[["ETP_fixed_threshold_2_25"]]),
  sprintf("ETP_boundary_stress_threshold_2_375 source: `%s`", roots[["ETP_boundary_stress_threshold_2_375"]]),
  sprintf("ETP_reference_balanced_threshold_2_24 source: `%s`", roots[["ETP_reference_balanced_threshold_2_24"]]),
  "",
  "## Interpretive boundary",
  "",
  "initial_ploidy initial ploidy is baseline. ETP endpoint ploidy is post-treatment and is also used to construct matched untreated ECDF references. ETP results are endpoint-defined descriptive/sensitivity analyses, not baseline predictive-biomarker evidence. ETP_reference_balanced_threshold_2_24 uses threshold 2.24 so each untreated ETP reference contains both initial sources, but Dose-specific nesting remains. ETP_boundary_stress_threshold_2_375 has no endpoint-high sample at 120 mg/kg; NE denotes non-estimability rather than no difference.",
  ""
)
writeLines(readme, file.path(comparison_root, "README.md"), useBytes = TRUE)

cat("Generated interpretation reports and comparisons in:\n")
cat("  ", file.path(roots[["ETP_fixed_threshold_2_25"]], "Manuscript"), "\n", sep = "")
cat("  ", file.path(roots[["ETP_boundary_stress_threshold_2_375"]], "Manuscript"), "\n", sep = "")
cat("  ", file.path(roots[["ETP_reference_balanced_threshold_2_24"]], "Manuscript"), "\n", sep = "")
cat("  ", comparison_root, "\n", sep = "")
