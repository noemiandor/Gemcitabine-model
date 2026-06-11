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

source(file.path(base_dir, "src", "dependencies.R"))
setup_local_lib(base_dir)
require_packages(c("EnrichIntersect", "xlsx", "plyr"))

source(file.path(tests_dir, "test_error_handling.R"))
source(file.path(tests_dir, "test_pubchem_client.R"))

cat("GDSC ploidy smoke tests passed.\n")
