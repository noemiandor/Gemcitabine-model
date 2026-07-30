#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/figure7/tests/testthat.R"
test_dir <- dirname(normalizePath(path))
Sys.setenv(FIGURE7_MODULE_DIR = dirname(test_dir))
testthat::test_dir(
  test_dir,
  reporter = "summary",
  filter = "contract|gsea-retry|statistics|workflow|seurat-upstream|species"
)
