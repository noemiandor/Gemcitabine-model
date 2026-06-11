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

cat("PubChem client fixture tests passed.\n")
