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
source(file.path(base_dir, "src", "pubchem_client.R"))

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Package 'jsonlite' is required for PubChem parser tests.", call. = FALSE)
}

read_fixture <- function(name) {
  jsonlite::fromJSON(file.path(tests_dir, "fixtures", "pubchem", name), simplifyVector = FALSE)
}

assert_equal <- function(observed, expected, label) {
  if (!identical(observed, expected)) {
    stop(sprintf("%s: expected '%s', observed '%s'", label, expected, observed), call. = FALSE)
  }
}

drug_classes <- classify_pubchem_json(read_fixture("drug_classes.json"))
assert_equal(drug_classes$category, "Antineoplastic Agents; Signal Transduction Inhibitors", "drug classes category")
assert_equal(drug_classes$method, "drug_classes", "drug classes method")

mechanism <- classify_pubchem_json(read_fixture("mechanism_only.json"))
assert_equal(mechanism$category, "Signaling", "mechanism fallback category")
assert_equal(mechanism$method, "mechanism_keyword_fallback", "mechanism fallback method")

no_category <- classify_pubchem_json(read_fixture("no_category.json"))
if (!is.na(no_category$category)) {
  stop("no_category fixture should return NA category", call. = FALSE)
}
assert_equal(no_category$method, "unclassified", "no category method")

malformed_path <- file.path(tests_dir, "fixtures", "pubchem", "malformed.json")
malformed_error <- tryCatch(
  {
    jsonlite::fromJSON(malformed_path, simplifyVector = FALSE)
    FALSE
  },
  error = function(e) TRUE
)
if (!malformed_error) {
  stop("Malformed JSON fixture should fail parsing.", call. = FALSE)
}

cache_cols <- pubchem_cache_columns()
cache <- data.frame(
  drug = c("DrugA", "DrugB", "DrugC"),
  drugName = c("DrugA", "DrugB", "DrugC"),
  pubchem_cid = c("1", "2", "3"),
  drugCategory_Pubchem = c("A", "B", "C"),
  annotation_source = "cache",
  annotation_status = "cached",
  retrieved_at = "old",
  pubchem_url = "",
  notes = "",
  stringsAsFactors = FALSE
)
refreshed <- data.frame(
  drug = "DrugB",
  drugName = "DrugB",
  pubchem_cid = "22",
  drugCategory_Pubchem = "B2",
  annotation_source = "pubchem_structured_json",
  annotation_status = "refreshed",
  retrieved_at = "new",
  pubchem_url = "",
  notes = "",
  stringsAsFactors = FALSE
)

merged_full <- merge_pubchem_refresh_results(cache, "DrugB", refreshed, subset_output = FALSE, cache_cols = cache_cols)
if (nrow(merged_full$cache) != 3) {
  stop("Default PubChem subset refresh should preserve unrequested cached rows", call. = FALSE)
}
if (merged_full$cache$drugCategory_Pubchem[merged_full$cache$drug == "DrugB"] != "B2") {
  stop("Default PubChem subset refresh should update the requested drug", call. = FALSE)
}

merged_subset <- merge_pubchem_refresh_results(cache, "DrugB", refreshed, subset_output = TRUE, cache_cols = cache_cols)
if (nrow(merged_subset$cache) != 1 || merged_subset$cache$drug != "DrugB") {
  stop("PubChem subset output should contain only requested drugs when explicitly enabled", call. = FALSE)
}

failed_refresh <- data.frame(
  drug = "DrugD",
  drugName = "DrugD",
  pubchem_cid = NA_character_,
  drugCategory_Pubchem = NA_character_,
  annotation_source = "pubchem_structured_json",
  annotation_status = "failed",
  retrieved_at = "new",
  pubchem_url = NA_character_,
  notes = "not found",
  stringsAsFactors = FALSE
)
missing_error <- tryCatch(
  {
    merge_pubchem_refresh_results(cache, "DrugD", failed_refresh, allow_missing = FALSE, cache_cols = cache_cols)
    FALSE
  },
  error = function(e) grepl("could not resolve", conditionMessage(e), fixed = TRUE)
)
if (!missing_error) {
  stop("PubChem refresh should fail unresolved requested drugs with no cached annotation", call. = FALSE)
}

allowed_missing <- merge_pubchem_refresh_results(cache, "DrugD", failed_refresh, allow_missing = TRUE, cache_cols = cache_cols)
if (!"DrugD" %in% allowed_missing$cache$drug) {
  stop("PubChem refresh with allow_missing should retain the failed annotation row for review", call. = FALSE)
}

cat("PubChem client fixture tests passed.\n")
