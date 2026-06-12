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
source(file.path(base_dir, "src", "analysis_helpers.R"))
source(file.path(base_dir, "src", "annotations.R"))

normalized <- normalize_legacy_drug_group(c(
  "Antineoplastic Agents; Signal Transduction Inhibitors",
  "Cytotoxic medicines",
  "(Anti-)Inflammatory"
))
expected <- c("SIGNALING", "CYTOTOXIC", "IMMUNOSUPPRESSIVE AGENTS")
if (!identical(normalized, expected)) {
  stop("Legacy category normalization changed unexpectedly", call. = FALSE)
}

custom_file <- tempfile(fileext = ".tsv")
write.table(
  data.frame(drug = c("DrugA", "DrugB"), group = c("Signaling", "Cytotoxic"), stringsAsFactors = FALSE),
  file = custom_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

annotations <- data.frame(
  drug = c("DrugA", "DrugB"),
  drugName = c("DrugA", "DrugB"),
  drugCategory_Pubchem_cache = c("Signal Transduction Inhibitors", NA),
  drugCategory_Pubchem = c("Signal Transduction Inhibitors", "Cytotoxic"),
  annotation_source = c("cache", "cache+manual_override"),
  annotation_status = c("cached", "manual_override"),
  stringsAsFactors = FALSE
)
audit <- build_drug_annotation_audit(annotations, custom_file)
if (!all(c("legacy_group_used_for_enrichment", "custom_set_candidate_category", "multi_category_flag") %in% colnames(audit))) {
  stop("Annotation audit is missing required columns", call. = FALSE)
}
if (!identical(audit$legacy_group_used_for_enrichment, c("SIGNALING", "CYTOTOXIC"))) {
  stop("Annotation audit did not record normalized legacy enrichment groups", call. = FALSE)
}

template_file <- tempfile(fileext = ".tsv")
template <- write_curated_mapping_template(audit, template_file)
if (!file.exists(template_file) || !all(template$curation_status == "unreviewed")) {
  stop("Curated mapping template was not written with unreviewed status", call. = FALSE)
}

cat("Annotation audit tests passed.\n")
