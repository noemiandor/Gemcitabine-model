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
