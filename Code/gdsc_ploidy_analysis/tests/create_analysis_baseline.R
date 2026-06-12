#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]], fixed = TRUE)
}

script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", cmd_args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("--file=", "", hit[[1]], fixed = TRUE))))
  }
  normalizePath(getwd())
}

tests_dir <- script_dir()
base_dir <- normalizePath(file.path(tests_dir, ".."))
baseline_dir <- normalizePath(arg_value("baseline-dir", file.path(base_dir, "baseline")), mustWork = FALSE)
output_dir <- normalizePath(arg_value("output-dir", file.path(base_dir, "output")), mustWork = FALSE)
data_dir <- normalizePath(arg_value("data-dir", file.path(base_dir, "data", "raw")), mustWork = FALSE)

dir.create(baseline_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(baseline_dir, "metadata"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(baseline_dir, "intermediate"), recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

checksum_files <- function(files, root, canonical_paths = NULL) {
  existing <- files[file.exists(files)]
  if (length(existing) == 0) {
    return(data.frame(path = character(), md5 = character()))
  }
  sums <- tools::md5sum(existing)
  if (is.null(canonical_paths)) {
    paths <- sub(paste0("^", normalizePath(root, mustWork = FALSE), "/?"), "", normalizePath(names(sums), mustWork = FALSE))
  } else {
    paths <- canonical_paths[match(names(sums), files)]
  }
  data.frame(
    path = paths,
    md5 = unname(sums),
    stringsAsFactors = FALSE
  )
}

capture.output(sessionInfo(), file = file.path(baseline_dir, "metadata", "session_info.txt"))

input_files <- c(
  file.path(data_dir, "GDSC2_fitted_dose_response_24Jul22.txt"),
  file.path(data_dir, "ploidyAcrossCellLines_V1.txt"),
  file.path(data_dir, "small_molecule_20200407234909.csv"),
  file.path(base_dir, "data", "manual", "custom_set_candidate.tsv"),
  file.path(base_dir, "data", "derived", "pubchem_drug_annotations.tsv")
)

output_files <- c(
  file.path(output_dir, "drugsVsPloidyCorr.RData"),
  file.path(output_dir, "coxIn.RData"),
  file.path(output_dir, "drugsVsPloidyCorr.xlsx"),
  file.path(output_dir, "drugsVsPloidyCorr.pdf"),
  file.path(output_dir, "ploidyVsDrugSensitivity.pdf"),
  file.path(output_dir, "metadata", "run_config.tsv"),
  file.path(output_dir, "metadata", "input_manifest.tsv"),
  file.path(output_dir, "metadata", "session_info.txt"),
  file.path(output_dir, "tables", "drug_category_counts_before_filter.tsv"),
  file.path(output_dir, "tables", "drug_category_counts_after_filter.tsv"),
  file.path(output_dir, "qc", "duplicate_drug_cell_line_records.tsv"),
  file.path(output_dir, "qc", "duplicate_resolution_summary.tsv")
)
generated_detail_files <- unlist(lapply(
  file.path(output_dir, c("metadata", "tables", "qc")),
  function(dir) {
    if (!dir.exists(dir)) {
      return(character())
    }
    list.files(dir, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  }
), use.names = FALSE)
output_files <- sort(unique(c(output_files, generated_detail_files)))
canonical_output_paths <- file.path(
  "output",
  sub(paste0("^", normalizePath(output_dir, mustWork = FALSE), "/?"), "", normalizePath(output_files, mustWork = FALSE))
)

checksums <- rbind(
  transform(checksum_files(input_files, base_dir), kind = "input"),
  transform(
    checksum_files(output_files, base_dir, canonical_paths = canonical_output_paths),
    kind = "output"
  )
)
checksums <- checksums[, c("kind", "path", "md5")]
write_tsv(checksums, file.path(baseline_dir, "checksums.tsv"))

result_summary <- data.frame(
  item = character(),
  value = character(),
  stringsAsFactors = FALSE
)

append_summary <- function(item, value) {
  result_summary <<- rbind(result_summary, data.frame(item = item, value = as.character(value), stringsAsFactors = FALSE))
}

correlation_file <- file.path(output_dir, "drugsVsPloidyCorr.RData")
if (file.exists(correlation_file)) {
  env <- new.env(parent = emptyenv())
  load(correlation_file, envir = env)
  if (exists("R", envir = env)) {
    R <- get("R", envir = env)
    correlation_rows <- do.call(rbind, lapply(names(R), function(cancer) {
      values <- R[[cancer]]
      data.frame(
        cancer = cancer,
        drug = names(values),
        correlation = as.numeric(values),
        stringsAsFactors = FALSE
      )
    }))
    write_tsv(correlation_rows, file.path(baseline_dir, "intermediate", "drug_ploidy_correlations.tsv"))
    append_summary("correlation_cancer_count", length(R))
    append_summary("correlation_row_count", nrow(correlation_rows))
  }
}

annotation_file <- file.path(output_dir, "coxIn.RData")
if (file.exists(annotation_file)) {
  env <- new.env(parent = emptyenv())
  load(annotation_file, envir = env)
  if (exists("coxIn", envir = env)) {
    coxIn <- get("coxIn", envir = env)
    coxIn_out <- coxIn
    coxIn_out$.rowname <- rownames(coxIn_out)
    coxIn_out <- coxIn_out[, c(".rowname", setdiff(colnames(coxIn_out), ".rowname")), drop = FALSE]
    write_tsv(coxIn_out, file.path(baseline_dir, "intermediate", "coxIn_saved.tsv"))
    append_summary("coxIn_row_count", nrow(coxIn))
    if ("group" %in% colnames(coxIn)) {
      category_counts <- as.data.frame(table(coxIn$group, useNA = "ifany"), stringsAsFactors = FALSE)
      colnames(category_counts) <- c("group", "n")
      write_tsv(category_counts, file.path(baseline_dir, "intermediate", "coxIn_group_counts.tsv"))
      append_summary("coxIn_group_count", length(unique(coxIn$group)))
    }
  }
}

workbook_file <- file.path(output_dir, "drugsVsPloidyCorr.xlsx")
if (file.exists(workbook_file) && requireNamespace("openxlsx", quietly = TRUE)) {
  sheets <- openxlsx::getSheetNames(workbook_file)
  append_summary("workbook_sheet_count", length(sheets))
  for (sheet in sheets) {
    sheet_data <- openxlsx::read.xlsx(workbook_file, sheet = sheet, rowNames = TRUE)
    sheet_data$.rowname <- rownames(sheet_data)
    sheet_data <- sheet_data[, c(".rowname", setdiff(colnames(sheet_data), ".rowname")), drop = FALSE]
    write_tsv(sheet_data, file.path(baseline_dir, "intermediate", paste0("workbook_", sheet, ".tsv")))
    append_summary(paste0("workbook_", sheet, "_rows"), nrow(sheet_data))
    append_summary(paste0("workbook_", sheet, "_cols"), ncol(sheet_data) - 1)
  }
} else if (file.exists(workbook_file)) {
  append_summary("workbook_parse_status", "skipped_openxlsx_not_installed")
}

write_tsv(result_summary, file.path(baseline_dir, "result_summary.tsv"))

cat("Baseline written to: ", baseline_dir, "\n", sep = "")
