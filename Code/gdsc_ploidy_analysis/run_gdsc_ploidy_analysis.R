#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
hit <- grep("--file=", args, value = TRUE)
if (length(hit) > 0) {
  base_dir <- dirname(normalizePath(sub("--file=", "", hit[1], fixed = TRUE)))
} else {
  base_dir <- normalizePath(getwd())
}

source(file.path(base_dir, "src", "run_gdsc_ploidy_analysis.R"))
