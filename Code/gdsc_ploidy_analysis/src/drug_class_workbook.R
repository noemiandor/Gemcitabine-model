primary_secondary_required_columns <- function() {
  c(
    "drug",
    "drug_key",
    "primary_anticancer_class",
    "secondary_suggested_class",
    "include_in_enrichment"
  )
}

primary_secondary_provenance_columns <- function() {
  c(
    "gdsc_putative_target",
    "gdsc_pathway_name",
    "classification_note",
    "revision_note",
    "source_annotation_used",
    "revision_source_basis",
    "pubchem_query_url",
    "drugbank_query_url"
  )
}

primary_secondary_ignored_legacy_columns <- function() {
  c(
    "approved_final_category_id",
    "approved_final_display_label",
    "category_mode",
    "category_source",
    "legacy_group_used_for_enrichment"
  )
}

normalize_excel_drug_key <- function(x) {
  out <- if (is.numeric(x)) {
    ifelse(is.na(x), NA_character_, format(x, scientific = FALSE, trim = TRUE, digits = 15))
  } else {
    as.character(x)
  }
  out <- trimws(out)
  out <- sub("^([0-9]+)\\.0+$", "\\1", out)
  normalize_drug_key(out)
}

parse_workbook_boolean <- function(x, column_name) {
  if (is.logical(x)) {
    if (any(is.na(x))) {
      stop("Missing boolean values in ", column_name, ".", call. = FALSE)
    }
    return(x)
  }
  value <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(value))
  out[value %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[value %in% c("false", "f", "0", "no", "n")] <- FALSE
  if (any(is.na(out))) {
    stop(
      "Invalid boolean values in ",
      column_name,
      ": ",
      paste(unique(as.character(x[is.na(out)])), collapse = ", "),
      call. = FALSE
    )
  }
  out
}

primary_secondary_class_slug <- function(x) {
  out <- normalize_drug_key(x)
  out <- gsub("[^A-Z0-9]+", "_", out)
  out <- gsub("_+", "_", out)
  gsub("^_|_$", "", out)
}

read_primary_secondary_class_order <- function(path, sheet = "Drug counts") {
  sheets <- openxlsx::getSheetNames(path)
  if (!sheet %in% sheets) {
    return(NULL)
  }
  counts <- openxlsx::read.xlsx(path, sheet = sheet, detectDates = FALSE)
  required <- c("primary_anticancer_class", "count")
  missing_cols <- setdiff(required, colnames(counts))
  if (length(missing_cols) > 0) {
    stop("Drug-count sheet is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  counts <- counts[nonempty(counts$primary_anticancer_class), , drop = FALSE]
  counts$primary_anticancer_class <- trimws(as.character(counts$primary_anticancer_class))
  counts$count <- suppressWarnings(as.integer(counts$count))
  if (any(is.na(counts$count) | counts$count < 0)) {
    stop("Drug-count sheet has invalid count values.", call. = FALSE)
  }
  if (any(duplicated(counts$primary_anticancer_class))) {
    stop(
      "Drug-count sheet has duplicate primary_anticancer_class values: ",
      paste(unique(counts$primary_anticancer_class[duplicated(counts$primary_anticancer_class)]), collapse = ", "),
      call. = FALSE
    )
  }
  counts[, c("primary_anticancer_class", "count"), drop = FALSE]
}

read_primary_secondary_drug_classes <- function(path,
                                                sheet = "drug_class_final_used_with_prim",
                                                count_sheet = "Drug counts") {
  if (!file.exists(path)) {
    stop("Primary/secondary drug-class workbook is missing: ", path, call. = FALSE)
  }
  sheets <- openxlsx::getSheetNames(path)
  if (!sheet %in% sheets) {
    stop("Primary/secondary drug-class workbook is missing sheet: ", sheet, call. = FALSE)
  }
  rows <- openxlsx::read.xlsx(path, sheet = sheet, detectDates = FALSE)
  required_cols <- unique(c(primary_secondary_required_columns(), primary_secondary_provenance_columns()))
  missing_cols <- setdiff(required_cols, colnames(rows))
  if (length(missing_cols) > 0) {
    stop("Primary/secondary workbook is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  for (col in intersect(colnames(rows), c(required_cols, primary_secondary_ignored_legacy_columns()))) {
    rows[[col]] <- trimws(as.character(rows[[col]]))
    rows[[col]][is.na(rows[[col]])] <- ""
  }
  rows$drug_key <- normalize_excel_drug_key(rows$drug_key)
  rows$drug <- trimws(as.character(rows$drug))
  rows$include_in_enrichment <- parse_workbook_boolean(rows$include_in_enrichment, "include_in_enrichment")
  rows$secondary_mechanism_note <- rows$secondary_suggested_class
  rows$primary_anticancer_class_slug <- primary_secondary_class_slug(rows$primary_anticancer_class)
  rows$source_workbook <- normalizePath(path, mustWork = TRUE)
  rows$source_workbook_md5 <- file_checksum(path)
  rows$source_workbook_sha256 <- file_sha256_checksum(path)

  if (any(!nonempty(rows$drug_key))) {
    stop("Primary/secondary workbook has missing drug_key values.", call. = FALSE)
  }
  if (any(duplicated(rows$drug_key))) {
    stop(
      "Primary/secondary workbook has duplicate normalized drug_key rows: ",
      paste(unique(rows$drug_key[duplicated(rows$drug_key)]), collapse = ", "),
      call. = FALSE
    )
  }
  included <- rows$include_in_enrichment
  if (any(included & !nonempty(rows$primary_anticancer_class))) {
    stop(
      "Included workbook rows have missing primary_anticancer_class values: ",
      paste(rows$drug_key[included & !nonempty(rows$primary_anticancer_class)], collapse = ", "),
      call. = FALSE
    )
  }

  count_rows <- read_primary_secondary_class_order(path, sheet = count_sheet)
  if (!is.null(count_rows)) {
    observed <- as.data.frame(table(rows$primary_anticancer_class[included]), stringsAsFactors = FALSE)
    colnames(observed) <- c("primary_anticancer_class", "count")
    observed$count <- as.integer(observed$count)
    observed <- observed[order(observed$primary_anticancer_class), , drop = FALSE]
    expected <- count_rows[order(count_rows$primary_anticancer_class), , drop = FALSE]
    rownames(observed) <- NULL
    rownames(expected) <- NULL
    if (!identical(observed, expected)) {
      stop("Drug-count sheet does not match computed included primary-class counts.", call. = FALSE)
    }
  }

  attr(rows, "class_order") <- if (!is.null(count_rows)) {
    count_rows$primary_anticancer_class
  } else {
    counts <- sort(table(rows$primary_anticancer_class[included]), decreasing = TRUE)
    order <- names(counts)
    c(setdiff(order, "Other"), intersect("Other", order))
  }
  rows
}

validate_primary_secondary_class_coverage <- function(correlation_drugs,
                                                      class_rows,
                                                      allow_extra_class_rows = FALSE) {
  correlation_keys <- normalize_drug_key(correlation_drugs$drug_key)
  included_rows <- class_rows[class_rows$include_in_enrichment, , drop = FALSE]
  class_keys <- normalize_excel_drug_key(included_rows$drug_key)
  missing_keys <- setdiff(correlation_keys, class_keys)
  extra_keys <- setdiff(class_keys, correlation_keys)
  duplicate_keys <- unique(class_keys[duplicated(class_keys)])
  validation_rows <- function(check, keys, status) {
    if (length(keys) == 0) {
      return(data.frame(check = character(), drug_key = character(), status = character(), stringsAsFactors = FALSE))
    }
    data.frame(
      check = rep(check, length(keys)),
      drug_key = keys,
      status = rep(status, length(keys)),
      stringsAsFactors = FALSE
    )
  }

  validation <- rbind(
    validation_rows("missing_from_workbook", missing_keys, "fail"),
    validation_rows("extra_in_workbook", extra_keys, if (allow_extra_class_rows) "warn" else "fail"),
    validation_rows("duplicate_in_workbook", duplicate_keys, "fail")
  )
  if (nrow(validation) == 0) {
    validation <- data.frame(
      check = c("coverage", "duplicates"),
      drug_key = "",
      status = "ok",
      stringsAsFactors = FALSE
    )
  }

  failures <- validation[validation$status == "fail", , drop = FALSE]
  if (nrow(failures) > 0) {
    stop(
      "Primary/secondary workbook coverage validation failed: ",
      paste(unique(failures$check), collapse = ", "),
      ". Review drug_class_workbook_validation.tsv.",
      call. = FALSE
    )
  }
  validation
}

make_enrichment_class_table <- function(class_rows) {
  included <- class_rows[class_rows$include_in_enrichment, , drop = FALSE]
  coxIn <- data.frame(
    drug = normalize_excel_drug_key(included$drug_key),
    group = trimws(included$primary_anticancer_class),
    category_id = primary_secondary_class_slug(included$primary_anticancer_class),
    secondary_mechanism_note = included$secondary_mechanism_note,
    stringsAsFactors = FALSE
  )
  rownames(coxIn) <- coxIn$drug
  coxIn
}

primary_secondary_class_counts <- function(class_rows) {
  included <- class_rows[class_rows$include_in_enrichment, , drop = FALSE]
  counts <- as.data.frame(table(included$primary_anticancer_class), stringsAsFactors = FALSE)
  colnames(counts) <- c("primary_anticancer_class", "count")
  counts$count <- as.integer(counts$count)
  low_count_threshold <- 5L
  counts$low_count_flag <- counts$count <= low_count_threshold
  counts$interpretation_note <- ifelse(
    counts$low_count_flag,
    "exploratory_low_count_class",
    ""
  )
  class_order <- attr(class_rows, "class_order")
  if (!is.null(class_order)) {
    counts$.order <- match(counts$primary_anticancer_class, class_order)
    counts <- counts[order(counts$.order, counts$primary_anticancer_class), , drop = FALSE]
    counts$.order <- NULL
  } else {
    counts <- counts[order(-counts$count, counts$primary_anticancer_class), , drop = FALSE]
  }
  counts
}

primary_secondary_used_table <- function(class_rows) {
  out <- class_rows[class_rows$include_in_enrichment, , drop = FALSE]
  keep <- c(
    "drug",
    "drug_key",
    "primary_anticancer_class",
    "primary_anticancer_class_slug",
    "secondary_mechanism_note",
    "include_in_enrichment",
    primary_secondary_provenance_columns(),
    "source_workbook",
    "source_workbook_sha256",
    "source_workbook_md5"
  )
  keep <- keep[keep %in% colnames(out)]
  names(out)[names(out) == "drug_key"] <- "normalized_drug_key"
  out$drug_key <- out$normalized_drug_key
  out <- out[, unique(c("drug", "drug_key", "normalized_drug_key", setdiff(keep, c("drug", "drug_key")))), drop = FALSE]
  out[order(out$normalized_drug_key), , drop = FALSE]
}

primary_secondary_review_summary <- function(class_counts,
                                             enrichment_long,
                                             low_ploidy_pvalue_cutoff = 0.05,
                                             high_ploidy_pvalue_cutoff = 0.1) {
  classes <- class_counts$primary_anticancer_class
  low <- enrichment_long[enrichment_long$direction == "low_ploidy_sensitive", , drop = FALSE]
  high <- enrichment_long[enrichment_long$direction == "high_ploidy_sensitive", , drop = FALSE]
  summarize_direction <- function(dt, cutoff, prefix) {
    rows <- lapply(classes, function(group) {
      vals <- suppressWarnings(as.numeric(dt$pvalue[dt$group == group]))
      vals <- vals[is.finite(vals)]
      data.frame(
        group = group,
        selected = sum(vals <= cutoff),
        min_pvalue = if (length(vals) == 0) NA_real_ else min(vals),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, rows)
    names(out) <- c("primary_anticancer_class", paste0("n_", prefix, "_selected"), paste0("min_", prefix, "_pvalue"))
    out
  }
  low_summary <- summarize_direction(low, low_ploidy_pvalue_cutoff, "low_ploidy")
  high_summary <- summarize_direction(high, high_ploidy_pvalue_cutoff, "high_ploidy")
  out <- merge(class_counts, low_summary, by = "primary_anticancer_class", all.x = TRUE, sort = FALSE)
  out <- merge(out, high_summary, by = "primary_anticancer_class", all.x = TRUE, sort = FALSE)
  out$.order <- match(out$primary_anticancer_class, classes)
  out <- out[order(out$.order), , drop = FALSE]
  out$.order <- NULL
  names(out)[names(out) == "count"] <- "n_drugs"
  out
}
