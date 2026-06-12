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
baseline_dir <- normalizePath(arg_value(args, "baseline-dir", file.path(base_dir, "baseline")), mustWork = FALSE)
output_dir <- normalizePath(arg_value(args, "output-dir", tempfile("ccle_ploidy_baseline_output_")), mustWork = FALSE)

source(file.path(base_dir, "src", "run_ccle_ploidy_analysis.R"))
run_ccle_ploidy_analysis(base_dir = base_dir, out_dir = output_dir, metric = "IC50", metric_source = "legacy")

dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(baseline_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(baseline_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)

files_to_copy <- c(
  file.path(output_dir, "result_summary.tsv"),
  file.path(output_dir, "tables", "drug_ploidy_correlations_plotted_ic50.tsv"),
  file.path(output_dir, "tables", "drug_ploidy_correlations_all_ic50.tsv"),
  file.path(output_dir, "tables", "drug_coverage_ic50.tsv"),
  file.path(output_dir, "metadata", "input_manifest.tsv"),
  file.path(output_dir, "metadata", "run_config.tsv"),
  file.path(output_dir, "metadata", "session_info.txt")
)
destinations <- file.path(baseline_dir, sub(paste0("^", output_dir, "/?"), "", files_to_copy))
for (i in seq_along(files_to_copy)) {
  dir.create(dirname(destinations[i]), recursive = TRUE, showWarnings = FALSE)
  file.copy(files_to_copy[i], destinations[i], overwrite = TRUE)
}

checksum_files <- destinations[file.exists(destinations)]
checksums <- data.frame(
  file = sub(paste0("^", baseline_dir, "/?"), "", checksum_files),
  md5 = unname(tools::md5sum(checksum_files)),
  stringsAsFactors = FALSE
)
write.table(checksums, file.path(baseline_dir, "checksums.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("Wrote baseline: ", baseline_dir, "\n", sep = "")
