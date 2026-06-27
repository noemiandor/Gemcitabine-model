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

normalize_drug_key <- function(x) {
  toupper(trimws(as.character(x)))
}

nonempty <- function(x) {
  !is.na(x) & nzchar(trimws(as.character(x)))
}

collapse_unique <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return("")
  }
  paste(sort(unique(x)), collapse = ";")
}

collapse_value_counts <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return("")
  }
  counts <- sort(table(x), decreasing = TRUE)
  paste(sprintf("%s(n=%d)", names(counts), as.integer(counts)), collapse = ";")
}

read_tsv_strings <- function(path) {
  read.table(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    colClasses = "character",
    stringsAsFactors = FALSE
  )
}

parse_boolean_column <- function(x, column_name) {
  value <- tolower(trimws(as.character(x)))
  bad <- !(value %in% c("true", "false"))
  if (any(bad)) {
    stop(
      "Invalid boolean values in ",
      column_name,
      ": ",
      paste(unique(as.character(x[bad])), collapse = ", "),
      call. = FALSE
    )
  }
  value == "true"
}

file_checksum <- function(path) {
  if (is.null(path) || is.na(path) || !file.exists(path)) {
    return(NA_character_)
  }
  unname(tools::md5sum(path))
}

parse_pubchem_categories <- function(category_string) {
  parse_one <- function(value) {
    if (is.na(value) || !nzchar(trimws(as.character(value)))) {
      return(character())
    }
    tokens <- unlist(strsplit(as.character(value), ";", fixed = TRUE), use.names = FALSE)
    tokens <- normalize_category_token(tokens)
    tokens[nonempty(tokens)]
  }
  if (length(category_string) <= 1) {
    return(parse_one(category_string))
  }
  lapply(category_string, parse_one)
}

normalize_category_token <- function(token) {
  token <- trimws(as.character(token))
  token <- gsub("[,;]+$", "", token)
  trimws(token)
}

# Legacy mode exists only for comparisons with historical outputs. It is not used
# for manuscript category assignment.
normalize_legacy_drug_group <- function(group) {
  last_category <- vapply(group, function(value) {
    tokens <- parse_pubchem_categories(value)
    if (length(tokens) == 0) {
      return(NA_character_)
    }
    tokens[length(tokens)]
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

split_semicolon_values <- function(x) {
  x <- trimws(as.character(x))
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) {
    return(character())
  }
  tokens <- unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE)
  tokens <- trimws(tokens)
  tokens[nonempty(tokens)]
}

read_category_schema <- function(schema_file) {
  if (!file.exists(schema_file)) {
    stop("Category schema file is missing: ", schema_file, call. = FALSE)
  }
  schema <- read_tsv_strings(schema_file)
  required_cols <- c(
    "schema_version",
    "category_id",
    "display_label",
    "parent_family",
    "plot_order",
    "include_in_enrichment_default",
    "manuscript_allowed",
    "min_count_policy",
    "min_count_threshold",
    "collapse_parent_category_id",
    "allowed_aliases",
    "deprecated_or_broad_terms",
    "is_deprecated_final_class",
    "notes"
  )
  missing_cols <- setdiff(required_cols, colnames(schema))
  if (length(missing_cols) > 0) {
    stop("Category schema is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  schema$category_id <- normalize_drug_key(schema$category_id)
  schema$display_label <- trimws(schema$display_label)
  if (any(!nonempty(schema$schema_version))) {
    stop("Category schema has missing schema_version values.", call. = FALSE)
  }
  if (any(!nonempty(schema$category_id))) {
    stop("Category schema has missing category_id values.", call. = FALSE)
  }
  if (any(duplicated(schema$category_id))) {
    stop("Category schema has duplicate category_id values: ", paste(unique(schema$category_id[duplicated(schema$category_id)]), collapse = ", "), call. = FALSE)
  }
  if (any(!nonempty(schema$display_label))) {
    stop("Category schema has missing display_label values.", call. = FALSE)
  }

  schema$include_in_enrichment_default <- parse_boolean_column(schema$include_in_enrichment_default, "include_in_enrichment_default")
  schema$manuscript_allowed <- parse_boolean_column(schema$manuscript_allowed, "manuscript_allowed")
  schema$is_deprecated_final_class <- parse_boolean_column(schema$is_deprecated_final_class, "is_deprecated_final_class")
  schema$plot_order <- suppressWarnings(as.integer(schema$plot_order))
  if (any(is.na(schema$plot_order))) {
    stop("Category schema has missing or non-integer plot_order values.", call. = FALSE)
  }
  manuscript_visible <- schema$manuscript_allowed & !schema$is_deprecated_final_class
  if (any(duplicated(schema$plot_order[manuscript_visible]))) {
    stop("Category schema has duplicated plot_order values among manuscript-visible categories.", call. = FALSE)
  }

  allowed_policies <- c("retain", "collapse_to_parent", "exclude_with_report")
  if (any(!schema$min_count_policy %in% allowed_policies)) {
    stop("Category schema has invalid min_count_policy values.", call. = FALSE)
  }
  schema$min_count_threshold <- suppressWarnings(as.integer(schema$min_count_threshold))
  if (any(is.na(schema$min_count_threshold) | schema$min_count_threshold < 1)) {
    stop("Category schema has invalid min_count_threshold values.", call. = FALSE)
  }
  needs_parent <- schema$min_count_policy == "collapse_to_parent"
  if (any(needs_parent & !schema$collapse_parent_category_id %in% schema$category_id)) {
    stop("Category schema has collapse_to_parent rows without a valid collapse_parent_category_id.", call. = FALSE)
  }

  aliases <- unlist(lapply(schema$allowed_aliases, split_semicolon_values), use.names = FALSE)
  aliases <- normalize_drug_key(aliases)
  aliases <- aliases[nonempty(aliases)]
  duplicated_aliases <- unique(aliases[duplicated(aliases)])
  if (length(duplicated_aliases) > 0) {
    stop("Category schema assigns aliases to multiple categories: ", paste(duplicated_aliases, collapse = ", "), call. = FALSE)
  }

  broad_terms <- normalize_drug_key(c(
    "ANTINEOPLASTIC AGENTS",
    "TARGETED THERAPIES",
    "ENZYME INHIBITORS",
    "SIGNAL TRANSDUCTION INHIBITORS"
  ))
  broad_final <- schema$category_id %in% broad_terms | normalize_drug_key(schema$display_label) %in% broad_terms
  if (any(broad_final & schema$manuscript_allowed & !schema$is_deprecated_final_class)) {
    stop("Broad administrative categories cannot be manuscript final classes.", call. = FALSE)
  }

  schema
}

map_gdsc_evidence_to_category <- function(gdsc_targets, gdsc_pathways) {
  target_text <- toupper(paste(gdsc_targets, collapse = ";"))
  pathway_text <- toupper(paste(gdsc_pathways, collapse = ";"))

  if (grepl("HDAC|BRD|METHYLTRANSFER|DNMT|CHROMATIN", target_text) ||
      grepl("CHROMATIN", pathway_text)) {
    return("EPIGENETICS_TRANSCRIPTION")
  }
  if (grepl("PYRIMIDINE|DHODH|ANTI-OXIDANT|ANTIOXIDANT", target_text)) {
    return("METABOLISM")
  }
  if (grepl("PROTEASOME|HSP90", target_text) ||
      grepl("PROTEIN STABILITY", pathway_text)) {
    return("PROTEIN_DEGRADATION")
  }
  if (grepl("PARP|TOP|DNA|NUCLEOSIDE", target_text) ||
      grepl("DNA REPLICATION|GENOME INTEGRITY", pathway_text)) {
    return("GENOME_INTEGRITY_DNA_DAMAGE")
  }
  if (grepl("MITOSIS|CELL CYCLE|CYTOSKELETON", pathway_text)) {
    return("CELL_CYCLE_MITOSIS")
  }
  if (grepl("APOPTOSIS", pathway_text)) {
    return("APOPTOSIS")
  }
  if (grepl("HORMONE", pathway_text)) {
    return("HORMONE_RELATED")
  }
  if (grepl("P53", pathway_text)) {
    return("P53_PATHWAY")
  }
  if (grepl("SIGNALING|KINASE|MTOR|RTK|WNT", pathway_text)) {
    return("SIGNALING")
  }
  if (grepl("METABOLISM", pathway_text)) {
    return("METABOLISM")
  }
  "OTHER_REVIEWED"
}

infer_candidate_categories <- function(pubchem_tokens = character(),
                                       gdsc_targets = character(),
                                       gdsc_pathways = character(),
                                       schema = NULL) {
  candidate <- map_gdsc_evidence_to_category(gdsc_targets, gdsc_pathways)
  if (!is.null(schema) && !candidate %in% schema$category_id) {
    return("OTHER_REVIEWED")
  }
  candidate
}

derive_correlation_eligible_drugs <- function(R, gdsc_rows = NULL) {
  drugs <- sort(unique(unlist(lapply(R, names), use.names = FALSE)))
  out <- data.frame(
    drug = drugs,
    drug_key = normalize_drug_key(drugs),
    stringsAsFactors = FALSE
  )
  if (!is.null(gdsc_rows)) {
    gdsc_rows$drug_key <- normalize_drug_key(gdsc_rows$DRUG_NAME)
    out$drug_id_list <- vapply(out$drug_key, function(key) {
      collapse_unique(gdsc_rows$DRUG_ID[gdsc_rows$drug_key == key])
    }, character(1))
  }
  out[order(out$drug_key), , drop = FALSE]
}

build_drug_evidence_table <- function(drug_key, gdsc_rows, pubchem_rows, manual_rows) {
  gdsc_sub <- gdsc_rows[normalize_drug_key(gdsc_rows$DRUG_NAME) == normalize_drug_key(drug_key), , drop = FALSE]
  pubchem_sub <- pubchem_rows[normalize_drug_key(pubchem_rows$drug) == normalize_drug_key(drug_key), , drop = FALSE]
  manual_sub <- manual_rows[normalize_drug_key(manual_rows$drug) == normalize_drug_key(drug_key), , drop = FALSE]
  pubchem_category <- if (nrow(pubchem_sub) > 0) pubchem_sub$drugCategory_Pubchem[[1]] else NA_character_
  pubchem_tokens <- parse_pubchem_categories(pubchem_category)
  data.frame(
    drug_key = normalize_drug_key(drug_key),
    drug_id_list = collapse_unique(gdsc_sub$DRUG_ID),
    gdsc_putative_target_values = collapse_unique(gdsc_sub$PUTATIVE_TARGET),
    gdsc_pathway_name_values = collapse_unique(gdsc_sub$PATHWAY_NAME),
    raw_drugCategory_Pubchem = pubchem_category,
    parsed_pubchem_category_tokens = paste(pubchem_tokens, collapse = ";"),
    current_manual_override = if (nrow(manual_sub) > 0) manual_sub$group[[1]] else "",
    proposed_final_category_id = infer_candidate_categories(
      pubchem_tokens,
      gdsc_sub$PUTATIVE_TARGET,
      gdsc_sub$PATHWAY_NAME
    ),
    stringsAsFactors = FALSE
  )
}

build_drug_class_curation_tables <- function(correlation_drugs,
                                             annotations,
                                             custom_file,
                                             gdsc_rows,
                                             wrong_assignment_file = NULL,
                                             schema = NULL) {
  custom_set <- read_custom_set_categories(custom_file)
  gdsc_rows$drug_key <- normalize_drug_key(gdsc_rows$DRUG_NAME)
  annotations$drug_key <- normalize_drug_key(annotations$drug)

  wrong_keys <- character()
  if (!is.null(wrong_assignment_file) && file.exists(wrong_assignment_file)) {
    wrong <- read_tsv_strings(wrong_assignment_file)
    if ("drug_key" %in% colnames(wrong)) {
      wrong_keys <- normalize_drug_key(wrong$drug_key)
    }
  }

  evidence_by_id <- unique(gdsc_rows[
    gdsc_rows$drug_key %in% correlation_drugs$drug_key,
    c("drug_key", "DRUG_ID", "DRUG_NAME", "PUTATIVE_TARGET", "PATHWAY_NAME"),
    drop = FALSE
  ])
  evidence_by_id <- evidence_by_id[order(evidence_by_id$drug_key, evidence_by_id$DRUG_ID), , drop = FALSE]

  rows <- lapply(seq_len(nrow(correlation_drugs)), function(i) {
    key <- correlation_drugs$drug_key[[i]]
    gdsc_sub <- gdsc_rows[gdsc_rows$drug_key == key, , drop = FALSE]
    anno_idx <- match(key, annotations$drug_key)
    anno <- if (!is.na(anno_idx)) annotations[anno_idx, , drop = FALSE] else data.frame()
    custom_idx <- match(key, custom_set$drug_key)
    custom_category <- if (!is.na(custom_idx)) custom_set$group[[custom_idx]] else ""
    pubchem_category <- if (nrow(anno) > 0) anno$drugCategory_Pubchem[[1]] else NA_character_
    pubchem_tokens <- parse_pubchem_categories(pubchem_category)
    target_values <- unique(gdsc_sub$PUTATIVE_TARGET)
    pathway_values <- unique(gdsc_sub$PATHWAY_NAME)
    proposed <- infer_candidate_categories(pubchem_tokens, target_values, pathway_values, schema)
    data.frame(
      drug = correlation_drugs$drug[[i]],
      drug_key = key,
      drug_id_list = collapse_unique(gdsc_sub$DRUG_ID),
      current_legacy_enrichment_group = normalize_legacy_drug_group(pubchem_category),
      raw_drugCategory_Pubchem = pubchem_category,
      parsed_pubchem_category_tokens = paste(pubchem_tokens, collapse = ";"),
      gdsc_putative_target_values = collapse_unique(target_values),
      gdsc_pathway_name_values = collapse_unique(pathway_values),
      gdsc_target_pathway_conflict_flag = length(unique(pathway_values[nonempty(pathway_values)])) > 1,
      rollup_conflict_status = if (length(unique(pathway_values[nonempty(pathway_values)])) > 1) {
        "requires_manual_resolution"
      } else if (length(unique(gdsc_sub$DRUG_ID)) > 1) {
        "resolved_same_category"
      } else {
        "none"
      },
      rollup_resolution_reason = if (length(unique(gdsc_sub$DRUG_ID)) > 1) {
        "Multiple GDSC DRUG_ID values share this drug key but have the same pathway-level category evidence."
      } else {
        ""
      },
      gdsc_putative_target_counts = collapse_value_counts(gdsc_sub$PUTATIVE_TARGET),
      gdsc_pathway_name_counts = collapse_value_counts(gdsc_sub$PATHWAY_NAME),
      current_manual_override = custom_category,
      in_wrong_assignment_fixture = key %in% wrong_keys,
      proposed_final_category_id = proposed,
      evidence_reason = paste0(
        "GDSC pathway: ",
        collapse_unique(pathway_values),
        "; target: ",
        collapse_unique(target_values)
      ),
      stringsAsFactors = FALSE
    )
  })
  list(
    evidence_by_drug_id = evidence_by_id,
    curation_input = do.call(rbind, rows)
  )
}

read_curated_drug_mapping <- function(mapping_file, schema, mode = "curated") {
  if (!file.exists(mapping_file)) {
    stop("Curated drug-class mapping file is missing: ", mapping_file, call. = FALSE)
  }
  mapping <- read_tsv_strings(mapping_file)
  required_cols <- c(
    "drug",
    "drug_key",
    "drug_id_list",
    "approved_final_category_id",
    "approved_final_display_label",
    "include_in_enrichment",
    "exclusion_reason_code",
    "exclusion_reason",
    "rollup_resolution_status",
    "rollup_resolution_reason",
    "curation_status",
    "primary_evidence_source",
    "evidence_summary",
    "curator",
    "curation_date",
    "curation_notes"
  )
  missing_cols <- setdiff(required_cols, colnames(mapping))
  if (length(missing_cols) > 0) {
    stop("Curated mapping is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  mapping$drug_key <- normalize_drug_key(mapping$drug_key)
  mapping$approved_final_category_id <- normalize_drug_key(mapping$approved_final_category_id)
  mapping$include_in_enrichment <- parse_boolean_column(mapping$include_in_enrichment, "include_in_enrichment")
  if (any(duplicated(mapping$drug_key))) {
    stop("Curated mapping has duplicate drug_key rows: ", paste(unique(mapping$drug_key[duplicated(mapping$drug_key)]), collapse = ", "), call. = FALSE)
  }

  approved_status <- tolower(trimws(mapping$curation_status)) == "approved"
  if (mode == "curated" && any(!approved_status)) {
    stop("Curated manuscript mode requires approved curation_status for every row.", call. = FALSE)
  }

  included <- mapping$include_in_enrichment
  missing_included_category <- included & !nonempty(mapping$approved_final_category_id)
  if (any(missing_included_category)) {
    stop("Included curated rows have missing approved_final_category_id: ", paste(mapping$drug_key[missing_included_category], collapse = ", "), call. = FALSE)
  }
  invalid_category <- nonempty(mapping$approved_final_category_id) & !mapping$approved_final_category_id %in% schema$category_id
  if (any(invalid_category)) {
    stop("Curated mapping has category IDs outside schema: ", paste(unique(mapping$approved_final_category_id[invalid_category]), collapse = ", "), call. = FALSE)
  }

  schema_idx <- match(mapping$approved_final_category_id, schema$category_id)
  schema_display <- schema$display_label[schema_idx]
  display_mismatch <- included & nonempty(mapping$approved_final_display_label) &
    !is.na(schema_display) & mapping$approved_final_display_label != schema_display
  if (any(display_mismatch)) {
    stop("Curated mapping display labels do not match schema for: ", paste(mapping$drug_key[display_mismatch], collapse = ", "), call. = FALSE)
  }
  mapping$approved_final_display_label[included] <- schema_display[included]

  disallowed <- included & (!schema$manuscript_allowed[schema_idx] | schema$is_deprecated_final_class[schema_idx])
  disallowed[is.na(disallowed)] <- FALSE
  if (mode == "curated" && any(disallowed)) {
    stop("Curated mapping includes manuscript-disallowed or deprecated categories for: ", paste(mapping$drug_key[disallowed], collapse = ", "), call. = FALSE)
  }

  broad_terms <- normalize_drug_key(c(
    "ANTINEOPLASTIC AGENTS",
    "TARGETED THERAPIES",
    "ENZYME INHIBITORS",
    "SIGNAL TRANSDUCTION INHIBITORS"
  ))
  broad_final <- included & mapping$approved_final_category_id %in% broad_terms
  if (mode == "curated" && any(broad_final)) {
    stop("Broad administrative categories cannot be used as curated manuscript final classes.", call. = FALSE)
  }

  allowed_exclusion_codes <- c(
    "AMBIGUOUS_MECHANISM",
    "NO_SPECIFIC_TARGET_EVIDENCE",
    "BROAD_ONLY_EVIDENCE",
    "DUPLICATE_ALIAS_COLLAPSED",
    "SALT_OR_PRODRUG_COLLAPSED",
    "CONTROL_OR_NON_THERAPEUTIC",
    "LOW_COUNT_CATEGORY_EXCLUDED",
    "CONFLICTING_GDSC_EVIDENCE",
    "OUT_OF_SCOPE_FOR_ENRICHMENT"
  )
  excluded <- !included
  bad_exclusion <- excluded & (!nonempty(mapping$exclusion_reason_code) | !nonempty(mapping$exclusion_reason))
  if (any(bad_exclusion)) {
    stop("Reviewed exclusions require exclusion_reason_code and exclusion_reason: ", paste(mapping$drug_key[bad_exclusion], collapse = ", "), call. = FALSE)
  }
  invalid_exclusion_code <- excluded & !mapping$exclusion_reason_code %in% allowed_exclusion_codes
  if (any(invalid_exclusion_code)) {
    stop("Curated mapping has invalid exclusion_reason_code values.", call. = FALSE)
  }

  mapping
}

validate_curated_mapping_coverage <- function(correlation_drugs, mapping, failure_file = NULL) {
  missing_keys <- setdiff(correlation_drugs$drug_key, mapping$drug_key)
  extra_keys <- setdiff(mapping$drug_key, correlation_drugs$drug_key)
  failures <- rbind(
    data.frame(
      drug_key = missing_keys,
      failure_type = rep("missing_from_curated_mapping", length(missing_keys)),
      stringsAsFactors = FALSE
    ),
    data.frame(
      drug_key = extra_keys,
      failure_type = rep("extra_in_curated_mapping", length(extra_keys)),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(failure_file)) {
    write_tsv(failures, failure_file)
  }
  if (length(missing_keys) > 0) {
    stop("Curated mapping is missing correlation-eligible drugs: ", paste(missing_keys, collapse = ", "), call. = FALSE)
  }
  if (length(extra_keys) > 0) {
    stop("Curated mapping contains drugs outside the correlation-eligible universe: ", paste(extra_keys, collapse = ", "), call. = FALSE)
  }
  invisible(failures)
}

resolve_final_drug_category <- function(drug, evidence, curated_mapping, schema, mode = "curated") {
  key <- normalize_drug_key(drug)
  row <- curated_mapping[curated_mapping$drug_key == key, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Could not resolve exactly one curated category for drug key: ", key, call. = FALSE)
  }
  row
}

resolve_final_drug_categories <- function(correlation_drugs,
                                          annotations,
                                          curated_mapping,
                                          schema,
                                          curation_input = NULL,
                                          mode = "curated",
                                          analysis_mode = "dev",
                                          metric = "Z_SCORE",
                                          schema_file = NULL,
                                          curated_mapping_file = NULL) {
  validate_curated_mapping_coverage(correlation_drugs, curated_mapping)
  annotations$drug_key <- normalize_drug_key(annotations$drug)
  rows <- lapply(seq_len(nrow(correlation_drugs)), function(i) {
    key <- correlation_drugs$drug_key[[i]]
    anno_idx <- match(key, annotations$drug_key)
    anno <- if (!is.na(anno_idx)) annotations[anno_idx, , drop = FALSE] else data.frame()
    map_row <- resolve_final_drug_category(key, NULL, curated_mapping, schema, mode)
    schema_row <- schema[match(map_row$approved_final_category_id, schema$category_id), , drop = FALSE]
    pubchem_category <- if (nrow(anno) > 0) anno$drugCategory_Pubchem[[1]] else NA_character_
    curation_idx <- if (!is.null(curation_input)) match(key, curation_input$drug_key) else NA_integer_
    curation_row <- if (!is.na(curation_idx)) curation_input[curation_idx, , drop = FALSE] else data.frame()
    data.frame(
      drug = correlation_drugs$drug[[i]],
      drug_key = key,
      drug_id_list = map_row$drug_id_list,
      approved_final_category_id = map_row$approved_final_category_id,
      approved_final_display_label = map_row$approved_final_display_label,
      legacy_group_used_for_enrichment = normalize_legacy_drug_group(pubchem_category),
      pubchem_category = pubchem_category,
      gdsc_putative_target = if (nrow(curation_row) > 0) curation_row$gdsc_putative_target_values[[1]] else "",
      gdsc_pathway_name = if (nrow(curation_row) > 0) curation_row$gdsc_pathway_name_values[[1]] else "",
      assignment_source = if (mode == "curated") "curated_mapping" else mode,
      include_in_enrichment = map_row$include_in_enrichment,
      included_in_enrichment = map_row$include_in_enrichment,
      exclusion_reason = map_row$exclusion_reason,
      exclusion_reason_code = map_row$exclusion_reason_code,
      low_count_filter_status = "not_evaluated",
      low_count_filter_reason = "",
      category_source = if (mode == "curated") "curated_mapping" else mode,
      category_mode = mode,
      analysis_metric = metric,
      schema_version = if (nrow(schema_row) > 0) schema_row$schema_version[[1]] else "",
      schema_checksum = file_checksum(schema_file),
      curated_mapping_version = file_checksum(curated_mapping_file),
      curated_mapping_checksum = file_checksum(curated_mapping_file),
      curation_status = map_row$curation_status,
      evidence_summary = map_row$evidence_summary,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$analysis_mode <- analysis_mode
  out
}

apply_category_low_count_policy <- function(final_used, schema) {
  final_used$low_count_filter_status <- ifelse(final_used$included_in_enrichment, "retained", "not_included_by_curation")
  final_used$low_count_filter_reason <- ""
  actions <- data.frame(
    category_id = character(),
    display_label = character(),
    included_count_before = integer(),
    min_count_threshold = integer(),
    min_count_policy = character(),
    action = character(),
    affected_drugs = character(),
    stringsAsFactors = FALSE
  )

  included <- final_used$included_in_enrichment
  counts <- table(final_used$approved_final_category_id[included])
  for (category_id in names(counts)) {
    schema_row <- schema[schema$category_id == category_id, , drop = FALSE]
    if (nrow(schema_row) != 1) {
      next
    }
    count <- as.integer(counts[[category_id]])
    threshold <- schema_row$min_count_threshold[[1]]
    policy <- schema_row$min_count_policy[[1]]
    action <- "retained"
    affected <- final_used$drug_key[included & final_used$approved_final_category_id == category_id]
    if (count < threshold) {
      if (policy == "exclude_with_report") {
        idx <- final_used$included_in_enrichment & final_used$approved_final_category_id == category_id
        final_used$included_in_enrichment[idx] <- FALSE
        final_used$low_count_filter_status[idx] <- "excluded_low_count"
        final_used$low_count_filter_reason[idx] <- paste0("Category count ", count, " below threshold ", threshold)
        action <- "excluded_low_count"
      } else if (policy == "collapse_to_parent") {
        parent_id <- schema_row$collapse_parent_category_id[[1]]
        parent_row <- schema[schema$category_id == parent_id, , drop = FALSE]
        idx <- final_used$included_in_enrichment & final_used$approved_final_category_id == category_id
        final_used$approved_final_category_id[idx] <- parent_id
        final_used$approved_final_display_label[idx] <- parent_row$display_label[[1]]
        final_used$low_count_filter_status[idx] <- "collapsed_to_parent"
        final_used$low_count_filter_reason[idx] <- paste0("Category count ", count, " below threshold ", threshold)
        action <- paste0("collapsed_to_", parent_id)
      } else {
        idx <- final_used$included_in_enrichment & final_used$approved_final_category_id == category_id
        final_used$low_count_filter_status[idx] <- "retained_below_threshold"
        final_used$low_count_filter_reason[idx] <- paste0("Category count ", count, " below threshold ", threshold, " but policy retain")
        action <- "retained_below_threshold"
      }
    }
    actions <- rbind(actions, data.frame(
      category_id = category_id,
      display_label = schema_row$display_label[[1]],
      included_count_before = count,
      min_count_threshold = threshold,
      min_count_policy = policy,
      action = action,
      affected_drugs = paste(affected, collapse = ";"),
      stringsAsFactors = FALSE
    ))
  }
  list(final_used = final_used, actions = actions)
}

build_drug_class_assignment_diff <- function(final_used) {
  data.frame(
    drug = final_used$drug,
    drug_key = final_used$drug_key,
    legacy_group_used_for_enrichment = final_used$legacy_group_used_for_enrichment,
    approved_final_category_id = final_used$approved_final_category_id,
    approved_final_display_label = final_used$approved_final_display_label,
    include_in_enrichment = final_used$include_in_enrichment,
    included_in_enrichment = final_used$included_in_enrichment,
    exclusion_reason_code = final_used$exclusion_reason_code,
    low_count_filter_status = final_used$low_count_filter_status,
    category_changed = final_used$legacy_group_used_for_enrichment != final_used$approved_final_category_id,
    evidence_summary = final_used$evidence_summary,
    stringsAsFactors = FALSE
  )
}

build_drug_class_category_counts <- function(final_used) {
  before <- as.data.frame(table(final_used$approved_final_category_id[final_used$include_in_enrichment]), stringsAsFactors = FALSE)
  colnames(before) <- c("approved_final_category_id", "count")
  before$count_stage <- "before_low_count_filter"
  after <- as.data.frame(table(final_used$approved_final_category_id[final_used$included_in_enrichment]), stringsAsFactors = FALSE)
  colnames(after) <- c("approved_final_category_id", "count")
  after$count_stage <- "after_low_count_filter"
  rbind(before[, c("count_stage", "approved_final_category_id", "count")], after[, c("count_stage", "approved_final_category_id", "count")])
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
  custom_category <- unname(custom_set$group[custom_idx])
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
    parsed_pubchem_category_tokens = unname(vapply(
      pubchem_category,
      function(value) paste(parse_pubchem_categories(value), collapse = ";"),
      character(1)
    )),
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
    drug_id_list = "",
    approved_final_category_id = "",
    approved_final_display_label = "",
    include_in_enrichment = "",
    exclusion_reason_code = "",
    exclusion_reason = "",
    rollup_resolution_status = "",
    rollup_resolution_reason = "",
    curation_status = "unreviewed",
    primary_evidence_source = "",
    evidence_summary = "",
    curator = "",
    curation_date = "",
    curation_notes = "",
    stringsAsFactors = FALSE
  )
  write_tsv(template, path)
  invisible(template)
}
