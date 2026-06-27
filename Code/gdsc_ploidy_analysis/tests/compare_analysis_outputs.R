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

stage_for_path <- function(path) {
  if (grepl("drug_class_curation_input|drug_class_curation_evidence", path)) {
    return("curation input")
  }
  if (grepl("drug_class_final_used|drug_class_assignment_diff|drug_class_category_counts|drug_class_low_count|drug_class_reviewed_exclusions", path)) {
    return("final category mapping")
  }
  if (grepl("class_enrichment_selected_drugs", path)) {
    return("selected drug table")
  }
  if (grepl("class_enrichment|drugsVsPloidyCorr_.*\\.xlsx|workbook_", path)) {
    return("enrichment matrix")
  }
  if (grepl("ploidy_enrichment", path)) {
    return("heatmap rendering")
  }
  "other tracked artifact"
}

append_checksum_diffs <- function(b, c) {
  b_sum <- read_tsv(b)
  c_sum <- read_tsv(c)
  b_sum$.key <- paste(b_sum$kind, b_sum$path, sep = "::")
  c_sum$.key <- paste(c_sum$kind, c_sum$path, sep = "::")
  all_keys <- sort(unique(c(b_sum$.key, c_sum$.key)))

  for (key in all_keys) {
    b_row <- b_sum[b_sum$.key == key, , drop = FALSE]
    c_row <- c_sum[c_sum$.key == key, , drop = FALSE]
    out_path <- paste0("checksums:", sub("::", "/", key, fixed = TRUE))
    if (nrow(b_row) == 0) {
      append_diff(out_path, "added", "artifact present only in candidate checksum manifest")
    } else if (nrow(c_row) == 0) {
      append_diff(out_path, "removed", "artifact present only in baseline checksum manifest")
    } else if (!identical(b_row$md5[[1]], c_row$md5[[1]])) {
      append_diff(out_path, "changed", paste("md5", b_row$md5[[1]], "->", c_row$md5[[1]]))
    }
  }
}

for (rel in all_rel) {
  b <- file.path(baseline_dir, rel)
  c <- file.path(candidate_dir, rel)
  if (!file.exists(b)) {
    append_diff(rel, "added", "present only in candidate")
  } else if (!file.exists(c)) {
    append_diff(rel, "removed", "present only in baseline")
  } else if (rel == "checksums.tsv") {
    append_checksum_diffs(b, c)
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
first_changed_stage <- if (changed) stage_for_path(diffs$path[[1]]) else "none"
report <- c(
  "# Regression Report",
  "",
  paste0("- Baseline: `", baseline_dir, "`"),
  paste0("- Candidate: `", candidate_dir, "`"),
  paste0("- Final result changed: `", if (changed) "yes" else "no", "`"),
  paste0("- First changed stage: `", first_changed_stage, "`"),
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
