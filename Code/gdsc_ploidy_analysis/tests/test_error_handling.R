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
source(file.path(base_dir, "analysis_helpers.R"))

expect_error_contains <- function(expr, expected, label) {
  msg <- tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
  if (is.na(msg)) {
    stop(label, ": expected an error but none was thrown", call. = FALSE)
  }
  missing <- expected[!vapply(expected, grepl, logical(1), x = msg, fixed = TRUE)]
  if (length(missing) > 0) {
    stop(label, ": missing expected text in error message: ", paste(missing, collapse = ", "), "\nObserved: ", msg, call. = FALSE)
  }
}

annotations <- data.frame(
  drug = c("A", "B"),
  group = c("SIGNALING", "CYTOTOXIC"),
  stringsAsFactors = FALSE
)
values <- matrix(c(0.1, -0.2), dimnames = list(c("A", "B"), "TEST"))

expect_error_contains(
  run_enrichment_or_stop(
    values,
    annotations,
    cancer = "TEST",
    direction = "low_ploidy_sensitive",
    permute_n = 1,
    pvalue_cutoff = 0.05,
    enrichment_fn = function(...) stop("forced enrichment failure")
  ),
  c("cancer=TEST", "direction=low_ploidy_sensitive", "forced enrichment failure"),
  "enrichment context"
)

expect_error_contains(
  plot_barplot_or_stop(
    c(A = 0.1),
    colors = c("red"),
    cancer = "TEST",
    xlab = "test metric",
    plotting_fn = function(...) stop("forced plot failure")
  ),
  c("cancer=TEST", "plotted_drugs=1", "forced plot failure"),
  "plot context"
)

expect_error_contains(
  validate_required_groups(c("SIGNALING"), required_groups = c("SIGNALING", "CYTOTOXIC")),
  c("Missing required annotation groups", "CYTOTOXIC"),
  "required group validation"
)

cat("Error handling tests passed.\n")
