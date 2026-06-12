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
run <- run_ccle_ploidy_analysis(base_dir = base_dir, out_dir = out_dir, metric = "IC50", metric_source = "legacy")

plotted <- run$result$plotted
stopifnot(file.exists(run$plot_file))
stopifnot(file.exists(file.path(out_dir, "tables", "drug_ploidy_correlations_plotted_ic50.tsv")))
stopifnot(nrow(plotted) == 45)
stopifnot("Gemcitabine" %in% plotted$drug)
gem <- plotted[plotted$drug == "Gemcitabine", ]
stopifnot(abs(gem$estimate - 0.4628112903956) < 1e-12)
stopifnot(gem$direction == "Low ploidy is sensitive")

cat("CCLE ploidy smoke tests passed.\n")
