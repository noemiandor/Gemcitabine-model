#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

baseline_dir <- arg_value("baseline-dir")
candidate_dir <- arg_value("candidate-dir")
report_file <- arg_value("report-file", file.path(candidate_dir %||% ".", "regression_report.md"))
diff_file <- arg_value("diff-file", file.path(candidate_dir %||% ".", "regression_diff_summary.tsv"))
allow_differences <- identical(tolower(arg_value("allow-differences", "false")), "true")

if (is.null(baseline_dir) || is.null(candidate_dir)) {
  stop("Usage: compare_analysis_outputs.R --baseline-dir=<dir> --candidate-dir=<dir> [--allow-differences=true]", call. = FALSE)
}

baseline_dir <- normalizePath(baseline_dir, mustWork = TRUE)
candidate_dir <- normalizePath(candidate_dir, mustWork = TRUE)
dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(diff_file), recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  utils::read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

relative_files <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  rel <- sub(paste0("^", normalizePath(root), "/?"), "", normalizePath(files, mustWork = FALSE))
  rel[
      basename(rel) != "regression_report.md" &
      basename(rel) != "regression_diff_summary.tsv" &
      basename(rel) != "pubchem_drug_annotations.refresh.tsv" &
      basename(rel) != "pubchem_refresh_log.tsv" &
      !grepl("^outputs/", rel)
  ]
}

file_md5 <- function(path) {
  if (!file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

all_rel <- sort(unique(c(relative_files(baseline_dir), relative_files(candidate_dir))))
diffs <- data.frame(
  path = character(),
  status = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

append_diff <- function(path, status, detail) {
  diffs <<- rbind(diffs, data.frame(path = path, status = status, detail = detail, stringsAsFactors = FALSE))
}

for (rel in all_rel) {
  b <- file.path(baseline_dir, rel)
  c <- file.path(candidate_dir, rel)
  if (!file.exists(b)) {
    append_diff(rel, "added", "present only in candidate")
  } else if (!file.exists(c)) {
    append_diff(rel, "removed", "present only in baseline")
  } else {
    b_md5 <- file_md5(b)
    c_md5 <- file_md5(c)
    if (!identical(b_md5, c_md5)) {
      append_diff(rel, "changed", paste("md5", b_md5, "->", c_md5))
    }
  }
}

write_tsv(diffs, diff_file)

changed <- nrow(diffs) > 0
report <- c(
  "# Regression Report",
  "",
  paste0("- Baseline: `", baseline_dir, "`"),
  paste0("- Candidate: `", candidate_dir, "`"),
  paste0("- Final result changed: `", if (changed) "yes" else "no", "`"),
  paste0("- Difference count: `", nrow(diffs), "`"),
  "",
  "## Difference Summary",
  ""
)

if (nrow(diffs) == 0) {
  report <- c(report, "No file-level differences detected.")
} else {
  for (i in seq_len(nrow(diffs))) {
    report <- c(report, paste0("- `", diffs$path[[i]], "`: ", diffs$status[[i]], " (", diffs$detail[[i]], ")"))
  }
}

report <- c(
  report,
  "",
  "## Attribution",
  "",
  if (changed) {
    "Differences require milestone-specific attribution before acceptance."
  } else {
    "No attribution required; candidate matches baseline at the tracked artifact level."
  }
)

writeLines(report, report_file)
cat("Regression report written to: ", report_file, "\n", sep = "")

if (changed && !allow_differences) {
  quit(status = 1)
}
