load_drug_annotations <- function(drugs, annotation_file, custom_file, allow_unclassified = TRUE) {
  stopifnot(file.exists(annotation_file), file.exists(custom_file))

  required_cols <- c(
    "drug",
    "drugName",
    "drugCategory_Pubchem",
    "annotation_source",
    "annotation_status"
  )
  annotations <- read.table(
    annotation_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  missing_cols <- setdiff(required_cols, colnames(annotations))
  if (length(missing_cols) > 0) {
    stop("Annotation cache is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  annotations$drugCategory_Pubchem[annotations$drugCategory_Pubchem == ""] <- NA_character_

  annotations$.drug_key <- toupper(annotations$drug)
  annotations <- annotations[!duplicated(annotations$.drug_key), , drop = FALSE]

  requested <- data.frame(
    drug = drugs,
    .drug_key = toupper(drugs),
    stringsAsFactors = FALSE
  )
  missing_keys <- setdiff(requested$.drug_key, annotations$.drug_key)
  if (length(missing_keys) > 0) {
    missing_drugs <- requested$drug[requested$.drug_key %in% missing_keys]
    stop(
      "Annotation cache is missing drugs required by analysis: ",
      paste(missing_drugs, collapse = ", "),
      call. = FALSE
    )
  }

  matched <- annotations[match(requested$.drug_key, annotations$.drug_key), , drop = FALSE]
  coxIn <- data.frame(
    drug = requested$drug,
    drugName = matched$drugName,
    drugCategory_Pubchem_cache = matched$drugCategory_Pubchem,
    drugCategory_Pubchem = matched$drugCategory_Pubchem,
    annotation_source = matched$annotation_source,
    annotation_status = matched$annotation_status,
    stringsAsFactors = FALSE
  )

  custom_set <- read.table(
    custom_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!all(c("drug", "group") %in% colnames(custom_set))) {
    stop("Manual override file must contain 'drug' and 'group' columns.", call. = FALSE)
  }

  rownames(custom_set) <- toupper(custom_set$drug)
  rownames(coxIn) <- toupper(coxIn$drug)

  empty_group <- is.na(coxIn$drugCategory_Pubchem) |
    coxIn$drugCategory_Pubchem == "," |
    coxIn$drugCategory_Pubchem == ""
  override_keys <- intersect(rownames(custom_set), rownames(coxIn[empty_group, , drop = FALSE]))
  if (length(override_keys) > 0) {
    coxIn[override_keys, "drugCategory_Pubchem"] <- custom_set[override_keys, "group"]
    coxIn[override_keys, "annotation_source"] <- paste(coxIn[override_keys, "annotation_source"], "manual_override", sep = "+")
    coxIn[override_keys, "annotation_status"] <- "manual_override"
  }

  still_unclassified <- is.na(coxIn$drugCategory_Pubchem) | coxIn$drugCategory_Pubchem == ""
  if (any(still_unclassified) && !allow_unclassified) {
    stop(
      "Unclassified drugs remain after applying annotation cache and overrides: ",
      paste(coxIn$drug[still_unclassified], collapse = ", "),
      call. = FALSE
    )
  }

  rownames(coxIn) <- NULL
  coxIn
}

normalize_legacy_drug_group <- function(group) {
  last_category <- vapply(group, function(value) {
    if (is.na(value) || !nzchar(value)) {
      return(NA_character_)
    }
    parts <- strsplit(value, "; ", fixed = TRUE)[[1]]
    parts[length(parts)]
  }, character(1))

  normalized <- gsub("Cytotoxic medicines", "Cytotoxic", gsub(";", "", gsub(",", "", last_category)))
  normalized[grep("Alkylating", normalized)] <- "Alkylating"
  normalized[grep("Topoisomerase", normalized)] <- "Cytotoxic"
  normalized[grep("Tubulin", normalized)] <- "Cytotoxic"
  normalized[grep("Antimitotic", normalized)] <- "Cytotoxic"
  normalized[grep("Antineoplastic Agents", normalized)] <- "Antineoplastic Agents"
  normalized <- toupper(normalized)
  normalized <- gsub("(ANTI-)INFLAMMATORY", "IMMUNOSUPPRESSIVE AGENTS", normalized, fixed = TRUE)
  normalized[
    normalized %in% c(
      "PARP INHIBITORS",
      "SIGNAL TRANSDUCTION INHIBITORS",
      "JAK INHIBITORS",
      "ENZYME INHIBITORS",
      "TARGETED THERAPIES"
    )
  ] <- "SIGNALING"
  unname(normalized)
}

read_custom_set_categories <- function(custom_file) {
  custom_set <- read.table(
    custom_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  if (!all(c("drug", "group") %in% colnames(custom_set))) {
    stop("Manual override file must contain 'drug' and 'group' columns.", call. = FALSE)
  }
  custom_set$drug_key <- toupper(custom_set$drug)
  custom_set <- custom_set[!duplicated(custom_set$drug_key), , drop = FALSE]
  custom_set
}

build_drug_annotation_audit <- function(annotations, custom_file) {
  required_cols <- c("drug", "drugName", "drugCategory_Pubchem", "annotation_source", "annotation_status")
  missing_cols <- setdiff(required_cols, colnames(annotations))
  if (length(missing_cols) > 0) {
    stop("Annotation audit is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  custom_set <- read_custom_set_categories(custom_file)
  annotations$drug_key <- toupper(annotations$drug)
  custom_idx <- match(annotations$drug_key, custom_set$drug_key)
  pubchem_category <- if ("drugCategory_Pubchem_cache" %in% colnames(annotations)) {
    annotations$drugCategory_Pubchem_cache
  } else {
    annotations$drugCategory_Pubchem
  }
  custom_category <- custom_set$group[custom_idx]
  legacy_group <- normalize_legacy_drug_group(annotations$drugCategory_Pubchem)
  normalized_custom <- normalize_legacy_drug_group(custom_category)
  normalized_pubchem <- normalize_legacy_drug_group(pubchem_category)

  data.frame(
    drug = annotations$drug,
    drug_key = annotations$drug_key,
    pubchem_category = pubchem_category,
    legacy_group_used_for_enrichment = legacy_group,
    custom_set_candidate_category = custom_category,
    would_manual_override_change_legacy_group = !is.na(normalized_custom) &
      (is.na(legacy_group) | normalized_custom != legacy_group),
    would_manual_override_change_pubchem_group = !is.na(normalized_custom) &
      (is.na(normalized_pubchem) | normalized_custom != normalized_pubchem),
    all_pubchem_categories = pubchem_category,
    annotation_source = annotations$annotation_source,
    annotation_status = annotations$annotation_status,
    multi_category_flag = !is.na(pubchem_category) & grepl(";", pubchem_category, fixed = TRUE),
    notes = ifelse(
      annotations$annotation_status == "manual_override",
      "Legacy enrichment category uses manual override because cached PubChem category was missing.",
      ""
    ),
    stringsAsFactors = FALSE
  )
}

write_curated_mapping_template <- function(annotation_audit, path) {
  template <- data.frame(
    drug = annotation_audit$drug,
    drug_key = annotation_audit$drug_key,
    legacy_group_used_for_enrichment = annotation_audit$legacy_group_used_for_enrichment,
    all_pubchem_categories = annotation_audit$all_pubchem_categories,
    custom_set_candidate_category = annotation_audit$custom_set_candidate_category,
    proposed_final_category = "",
    approved_final_category = "",
    curation_status = "unreviewed",
    curator = "",
    curation_notes = "",
    stringsAsFactors = FALSE
  )
  write_tsv(template, path)
  invisible(template)
}
