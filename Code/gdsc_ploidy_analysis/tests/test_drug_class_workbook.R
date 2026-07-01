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
source(file.path(base_dir, "src", "drug_class_workbook.R"))

expect_error_like <- function(expr, pattern, label) {
  err <- tryCatch(
    {
      force(expr)
      NULL
    },
    error = function(e) e$message
  )
  if (is.null(err) || !grepl(pattern, err)) {
    stop("Expected error for ", label, " matching '", pattern, "', got: ", ifelse(is.null(err), "<no error>", err), call. = FALSE)
  }
}

workbook <- file.path(base_dir, "data", "manual", "drug_class_final_used_with_primary_secondary_corrected.xlsx")
rows <- read_primary_secondary_drug_classes(workbook)

if (nrow(rows) != 286L) {
  stop("Primary/secondary workbook should contain 286 drug rows.", call. = FALSE)
}
if (any(duplicated(rows$drug_key))) {
  stop("Primary/secondary workbook contains duplicate normalized drug keys.", call. = FALSE)
}
if (!all(rows$include_in_enrichment)) {
  stop("Expected all rows in the current workbook to be included in enrichment.", call. = FALSE)
}

numeric_keys <- c("123138", "123829", "150412", "50869", "615590", "630600", "667880",
                  "720427", "729189", "741909", "743380", "765771", "776928")
if (!all(numeric_keys %in% rows$drug_key)) {
  stop("Numeric-looking workbook keys were not preserved exactly after normalization.", call. = FALSE)
}

class_order <- attr(rows, "class_order")
if (length(class_order) != 19L || any(!class_order %in% rows$primary_anticancer_class)) {
  stop("Workbook class order should contain the 19 primary anticancer classes.", call. = FALSE)
}

coverage <- validate_primary_secondary_class_coverage(
  data.frame(drug = rows$drug_key, drug_key = rows$drug_key, stringsAsFactors = FALSE),
  rows
)
if (!all(coverage$status == "ok")) {
  stop("Coverage validation should pass for a matching universe.", call. = FALSE)
}

coxIn <- make_enrichment_class_table(rows)
if (!all(c("drug", "group", "category_id", "secondary_mechanism_note") %in% colnames(coxIn))) {
  stop("Primary/secondary coxIn table is missing required columns.", call. = FALSE)
}
if (!identical(sort(unique(coxIn$group)), sort(unique(rows$primary_anticancer_class)))) {
  stop("coxIn groups must come only from primary_anticancer_class.", call. = FALSE)
}

rows_with_broken_legacy <- rows
for (col in intersect(c("approved_final_category_id", "approved_final_display_label", "category_mode", "legacy_group_used_for_enrichment"), colnames(rows_with_broken_legacy))) {
  rows_with_broken_legacy[[col]] <- paste0("BROKEN_", seq_len(nrow(rows_with_broken_legacy)))
}
coxIn_broken_legacy <- make_enrichment_class_table(rows_with_broken_legacy)
if (!identical(coxIn[, c("drug", "group", "category_id")], coxIn_broken_legacy[, c("drug", "group", "category_id")])) {
  stop("Ignored legacy workbook columns changed primary/secondary enrichment assignments.", call. = FALSE)
}

used <- primary_secondary_used_table(rows)
required_used <- c(
  "drug",
  "drug_key",
  "normalized_drug_key",
  "primary_anticancer_class",
  "primary_anticancer_class_slug",
  "secondary_mechanism_note",
  "include_in_enrichment",
  "source_workbook_sha256"
)
if (!all(required_used %in% colnames(used))) {
  stop("Primary/secondary used table is missing required audit columns.", call. = FALSE)
}

tmp_xlsx <- function(data) {
  path <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, "drug_class_final_used_with_prim")
  openxlsx::writeData(wb, "drug_class_final_used_with_prim", data)
  openxlsx::addWorksheet(wb, "Drug counts")
  counts <- as.data.frame(table(data$primary_anticancer_class), stringsAsFactors = FALSE)
  names(counts) <- c("primary_anticancer_class", "count")
  openxlsx::writeData(wb, "Drug counts", counts)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

bad_missing <- rows[1:3, setdiff(colnames(rows), "primary_anticancer_class"), drop = FALSE]
expect_error_like(
  read_primary_secondary_drug_classes(tmp_xlsx(bad_missing)),
  "missing required columns",
  "missing required primary class column"
)

bad_duplicate <- rows[1:3, , drop = FALSE]
bad_duplicate$drug_key[[2]] <- bad_duplicate$drug_key[[1]]
expect_error_like(
  read_primary_secondary_drug_classes(tmp_xlsx(bad_duplicate)),
  "duplicate normalized drug_key",
  "duplicate normalized workbook key"
)

expect_error_like(
  validate_primary_secondary_class_coverage(
    data.frame(drug = c("A", "B"), drug_key = c("A", "B"), stringsAsFactors = FALSE),
    rows[1, , drop = FALSE]
  ),
  "coverage validation failed",
  "missing and extra coverage"
)

cat("Primary/secondary drug-class workbook tests passed.\n")
