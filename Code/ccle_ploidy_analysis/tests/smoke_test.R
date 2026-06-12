#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("--file=", "", hit[1], fixed = TRUE))))
  }
  normalizePath(getwd())
}

tests_dir <- script_dir()
base_dir <- normalizePath(file.path(tests_dir, ".."))
out_dir <- tempfile("ccle_ploidy_smoke_")

source(file.path(base_dir, "src", "run_ccle_ploidy_analysis.R"))
run <- run_ccle_ploidy_analysis(base_dir = base_dir, out_dir = out_dir)

plotted <- run$result$plotted
stopifnot(file.exists(run$plot_file))
stopifnot(file.exists(file.path(out_dir, "tables", "drug_ploidy_correlations_plotted_z_score.tsv")))
stopifnot(!file.exists(file.path(out_dir, "tables", "drug_ploidy_correlations_plotted_ic50.tsv")))
stopifnot(run$result$metric == "Z_SCORE")
stopifnot(run$result$metric_label == "BreastCancerDrugSensitivity Z Score")
stopifnot(run$result$response_column == "Z Score")
stopifnot(isTRUE(run$result$lower_metric_more_sensitive))
stopifnot(nrow(plotted) == 45)
stopifnot("Gemcitabine" %in% plotted$drug)
gem <- plotted[plotted$drug == "Gemcitabine", ]
stopifnot(abs(gem$estimate - 0.4628112903956) < 1e-12)
stopifnot(gem$direction == "Low ploidy is sensitive")

run_config <- read.table(file.path(out_dir, "metadata", "run_config.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
stopifnot(!any(grepl("IC50", run_config$value, fixed = TRUE)))
stopifnot(!grepl("ic50", basename(run$plot_file), ignore.case = TRUE))

alias_out_dir <- tempfile("ccle_ploidy_alias_")
alias_warning <- NULL
alias_run <- withCallingHandlers(
  run_ccle_ploidy_analysis(
    base_dir = base_dir,
    out_dir = alias_out_dir,
    metric = "IC50",
    metric_source = "legacy"
  ),
  warning = function(w) {
    alias_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(!is.null(alias_warning))
stopifnot(grepl("deprecated", alias_warning, fixed = TRUE))
stopifnot(alias_run$result$metric == "Z_SCORE")
stopifnot(file.exists(file.path(alias_out_dir, "tables", "drug_ploidy_correlations_plotted_z_score.tsv")))

regression_report <- file.path(base_dir, "baseline", "regression", "legacy_ic50_labeled_to_z_score_numeric_equivalence.tsv")
stopifnot(file.exists(regression_report))
regression <- read.table(regression_report, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
stopifnot(nrow(regression) == nrow(plotted))
stopifnot(all(regression$numeric_match))

cat("CCLE ploidy smoke tests passed.\n")
