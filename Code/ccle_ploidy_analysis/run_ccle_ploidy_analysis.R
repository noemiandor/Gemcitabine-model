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

base_dir <- script_dir()
source(file.path(base_dir, "src", "run_ccle_ploidy_analysis.R"))
ccle_ploidy_analysis_main(base_dir)
