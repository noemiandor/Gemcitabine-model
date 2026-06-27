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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

expect_error <- function(expr, pattern, label) {
  err <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e$message
  )
  if (is.null(err) || !grepl(pattern, err)) {
    stop("Expected error for ", label, " matching '", pattern, "', got: ", err %||% "<no error>", call. = FALSE)
  }
}

tokens <- parse_pubchem_categories(" Antineoplastic Agents ; ; PARP Inhibitors;  Signal Transduction Inhibitors; ")
expected_tokens <- c("Antineoplastic Agents", "PARP Inhibitors", "Signal Transduction Inhibitors")
if (!identical(tokens, expected_tokens)) {
  stop("PubChem category parsing did not preserve all non-empty tokens", call. = FALSE)
}

schema_file <- file.path(base_dir, "data", "manual", "drug_class_category_schema.tsv")
mapping_file <- file.path(base_dir, "data", "manual", "drug_class_final_curated.tsv")
fixture_file <- file.path(tests_dir, "fixtures", "wrong_assignment_examples.tsv")
schema <- read_category_schema(schema_file)
mapping <- read_curated_drug_mapping(mapping_file, schema, mode = "curated")
fixture <- read_tsv_strings(fixture_file)

anti_candidate <- infer_candidate_categories("(Anti-)Inflammatory", character(), "Other", schema)
if (identical(anti_candidate, "IMMUNOSUPPRESSIVE AGENTS")) {
  stop("Final curated/proposal resolution must not treat (Anti-)Inflammatory as an automatic immunosuppressive category", call. = FALSE)
}

fixture_idx <- match(normalize_drug_key(fixture$drug_key), mapping$drug_key)
if (any(is.na(fixture_idx))) {
  stop("Wrong-assignment fixture contains drugs missing from curated mapping", call. = FALSE)
}
observed <- mapping[fixture_idx, , drop = FALSE]
expected_include <- parse_boolean_column(fixture$expected_include_in_enrichment, "expected_include_in_enrichment")
if (!identical(observed$approved_final_category_id, normalize_drug_key(fixture$expected_category_id))) {
  stop("Curated mapping does not match expected fixture category IDs", call. = FALSE)
}
if (!identical(observed$include_in_enrichment, expected_include)) {
  stop("Curated mapping does not match expected fixture inclusion flags", call. = FALSE)
}
expected_exclusion <- fixture$expected_exclusion_reason_code
expected_exclusion[is.na(expected_exclusion)] <- ""
if (!identical(observed$exclusion_reason_code, expected_exclusion)) {
  stop("Curated mapping does not match expected fixture exclusion reason codes", call. = FALSE)
}

if (mapping$approved_final_category_id[mapping$drug_key == "OLAPARIB"] == "SIGNALING") {
  stop("PARP inhibitors must not resolve to SIGNALING", call. = FALSE)
}
if (mapping$approved_final_category_id[mapping$drug_key == "JQ1"] == "IMMUNOSUPPRESSIVE AGENTS") {
  stop("JQ1 must not resolve to IMMUNOSUPPRESSIVE AGENTS", call. = FALSE)
}
broad_terms <- c("ANTINEOPLASTIC AGENTS", "TARGETED THERAPIES", "ENZYME INHIBITORS", "SIGNAL TRANSDUCTION INHIBITORS")
if (any(mapping$approved_final_category_id %in% broad_terms)) {
  stop("Broad administrative categories are present as final curated category IDs", call. = FALSE)
}

mini_universe <- data.frame(
  drug = c("DrugA", "DrugB"),
  drug_key = c("DRUGA", "DRUGB"),
  stringsAsFactors = FALSE
)
mini_mapping <- mapping[1, , drop = FALSE]
mini_mapping$drug <- "DrugA"
mini_mapping$drug_key <- "DRUGA"
failure_file <- tempfile(fileext = ".tsv")
expect_error(
  validate_curated_mapping_coverage(mini_universe, mini_mapping, failure_file = failure_file),
  "missing correlation-eligible drugs",
  "missing curated key coverage"
)
failure_rows <- read_tsv_strings(failure_file)
if (!identical(failure_rows$drug_key, "DRUGB")) {
  stop("Missing-key validation did not write the expected failure row", call. = FALSE)
}

bad_mapping <- mapping
bad_mapping$curation_status[[1]] <- "unreviewed"
bad_file <- tempfile(fileext = ".tsv")
write_tsv(bad_mapping, bad_file)
expect_error(
  read_curated_drug_mapping(bad_file, schema, mode = "curated"),
  "approved curation_status",
  "unapproved curation status"
)

bad_mapping <- mapping
bad_mapping$approved_final_category_id[[1]] <- "SIGNAL TRANSDUCTION INHIBITORS"
bad_mapping$approved_final_display_label[[1]] <- ""
bad_file <- tempfile(fileext = ".tsv")
write_tsv(bad_mapping, bad_file)
expect_error(
  read_curated_drug_mapping(bad_file, schema, mode = "curated"),
  "outside schema",
  "broad invalid category"
)

bad_schema <- schema
bad_schema$plot_order[[2]] <- bad_schema$plot_order[[1]]
bad_schema_file <- tempfile(fileext = ".tsv")
write_tsv(bad_schema, bad_schema_file)
expect_error(
  read_category_schema(bad_schema_file),
  "duplicated plot_order",
  "duplicate schema plot order"
)

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
  drugCategory_Pubchem_cache = c("Antineoplastic Agents; Signal Transduction Inhibitors", NA),
  drugCategory_Pubchem = c("Antineoplastic Agents; Signal Transduction Inhibitors", "Cytotoxic"),
  annotation_source = c("cache", "cache+manual_override"),
  annotation_status = c("cached", "manual_override"),
  stringsAsFactors = FALSE
)
audit <- build_drug_annotation_audit(annotations, custom_file)
if (!all(c("legacy_group_used_for_enrichment", "custom_set_candidate_category", "parsed_pubchem_category_tokens") %in% colnames(audit))) {
  stop("Annotation audit is missing required columns", call. = FALSE)
}

template_file <- tempfile(fileext = ".tsv")
template <- write_curated_mapping_template(audit, template_file)
if (!file.exists(template_file) || !all(template$curation_status == "unreviewed")) {
  stop("Curated mapping template was not written with unreviewed status", call. = FALSE)
}
if (!all(c("approved_final_category_id", "include_in_enrichment", "exclusion_reason_code") %in% colnames(template))) {
  stop("Curated mapping template does not expose required curated-mapping columns", call. = FALSE)
}

cat("Annotation and curated category tests passed.\n")
