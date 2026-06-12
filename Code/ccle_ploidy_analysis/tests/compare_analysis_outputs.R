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

arg_value <- function(args, name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[1], fixed = TRUE)
}

tests_dir <- script_dir()
base_dir <- normalizePath(file.path(tests_dir, ".."))
args <- commandArgs(trailingOnly = TRUE)
baseline_dir <- normalizePath(arg_value(args, "baseline-dir", file.path(base_dir, "baseline")), mustWork = TRUE)
candidate_dir <- normalizePath(arg_value(args, "candidate-dir", tempfile("ccle_ploidy_candidate_")), mustWork = FALSE)

source(file.path(base_dir, "src", "run_ccle_ploidy_analysis.R"))
run_ccle_ploidy_analysis(base_dir = base_dir, out_dir = candidate_dir)

baseline <- read.table(file.path(baseline_dir, "checksums.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
candidate_files <- file.path(candidate_dir, baseline$file)
candidate <- data.frame(
  file = baseline$file,
  md5 = unname(tools::md5sum(candidate_files)),
  stringsAsFactors = FALSE
)
comparison <- merge(baseline, candidate, by = "file", suffixes = c("_baseline", "_candidate"), all = TRUE)
comparison$match <- comparison$md5_baseline == comparison$md5_candidate
comparison$match[is.na(comparison$match)] <- FALSE

report_path <- file.path(candidate_dir, "baseline_comparison.tsv")
write.table(comparison, report_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
if (!all(comparison$match)) {
  stop("Candidate outputs differ from baseline. See ", report_path, call. = FALSE)
}
cat("Candidate outputs match baseline.\n")
