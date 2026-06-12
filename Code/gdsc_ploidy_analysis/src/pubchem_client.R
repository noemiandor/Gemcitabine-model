pubchem_name_to_cid <- function(drug, cmap_csv) {
  id <- NA_character_
  if (requireNamespace("ChemmineR", quietly = TRUE)) {
    id <- suppressWarnings(ChemmineR::pubchemName2CID(drug)[1])
  }
  if ((length(id) == 0 || is.na(id)) && file.exists(cmap_csv)) {
    drugs <- read.csv(cmap_csv, stringsAsFactors = FALSE)
    if (all(c("Name", "PubChem CID") %in% colnames(drugs))) {
      id <- unique(drugs[["PubChem CID"]][drugs$Name == drug])[1]
    } else if (all(c("Name", "PubChem.CID") %in% colnames(drugs))) {
      id <- unique(drugs[["PubChem.CID"]][drugs$Name == drug])[1]
    }
  }
  if (length(id) == 0 || is.na(id) || !nzchar(as.character(id))) {
    return(NA_character_)
  }
  as.character(id)
}

pubchem_pug_view_url <- function(cid) {
  sprintf("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/JSON", cid)
}

fetch_pubchem_json <- function(cid) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required to refresh PubChem annotations.", call. = FALSE)
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required to refresh PubChem annotations.", call. = FALSE)
  }
  url <- pubchem_pug_view_url(cid)
  response <- httr2::request(url) |>
    httr2::req_timeout(30) |>
    httr2::req_retry(max_tries = 3) |>
    httr2::req_perform()
  text <- httr2::resp_body_string(response)
  jsonlite::fromJSON(text, simplifyVector = FALSE)
}

as_node_list <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (is.list(x) && is.null(names(x))) {
    return(x)
  }
  list(x)
}

collect_pubchem_sections <- function(node) {
  found <- list()
  if (!is.list(node)) {
    return(found)
  }
  if (!is.null(node$TOCHeading)) {
    found <- c(found, list(node))
  }
  for (section in as_node_list(node$Section)) {
    found <- c(found, collect_pubchem_sections(section))
  }
  found
}

extract_pubchem_strings <- function(node) {
  strings <- character()
  if (is.null(node) || !is.list(node)) {
    return(strings)
  }
  if (!is.null(node$String)) {
    strings <- c(strings, unlist(node$String, use.names = FALSE))
  }
  if (!is.null(node$StringWithMarkup)) {
    for (entry in as_node_list(node$StringWithMarkup)) {
      strings <- c(strings, extract_pubchem_strings(entry))
    }
  }
  if (!is.null(node$Value)) {
    strings <- c(strings, extract_pubchem_strings(node$Value))
  }
  if (!is.null(node$Information)) {
    for (entry in as_node_list(node$Information)) {
      strings <- c(strings, extract_pubchem_strings(entry))
    }
  }
  unique(strings[nzchar(strings)])
}

extract_pubchem_section_strings <- function(pubchem_json, heading) {
  root <- pubchem_json$Record %||% pubchem_json
  sections <- collect_pubchem_sections(root)
  matches <- sections[vapply(sections, function(section) identical(section$TOCHeading, heading), logical(1))]
  unique(unlist(lapply(matches, extract_pubchem_strings), use.names = FALSE))
}

classify_pubchem_mechanism_text <- function(text) {
  if (length(text) == 0) {
    return(NA_character_)
  }
  if (any(grepl("inflammatory|immun|lympho|cytokine", text, ignore.case = TRUE))) {
    return("(Anti-)Inflammatory")
  }
  if (any(grepl("platinum|Cytotox|Microtubul|Mitotic|Spindle|Alkylating", text, ignore.case = TRUE))) {
    return("Cytotoxic")
  }
  if (any(grepl("MAPK|ERK kinase|MEK1|MEK2|MTOR|WNT|EGFR|ROS", text)) ||
      any(grepl("kinase inhibitor", text, ignore.case = TRUE))) {
    return("Signaling")
  }
  if (any(grepl("Metabolic Inhibitor|detoxification", text, ignore.case = TRUE))) {
    return("MetabolicInhibitor")
  }
  NA_character_
}

classify_pubchem_json <- function(pubchem_json) {
  drug_classes <- extract_pubchem_section_strings(pubchem_json, "Drug Classes")
  if (length(drug_classes) > 0) {
    return(list(
      category = paste(drug_classes, collapse = "; "),
      method = "drug_classes",
      note = paste("Drug Classes values:", paste(drug_classes, collapse = " | "))
    ))
  }

  mechanism <- extract_pubchem_section_strings(pubchem_json, "Mechanism of Action")
  fallback <- classify_pubchem_mechanism_text(mechanism)
  if (!is.na(fallback)) {
    return(list(
      category = fallback,
      method = "mechanism_keyword_fallback",
      note = paste("Mechanism text matched fallback heuristic:", paste(mechanism, collapse = " | "))
    ))
  }

  list(
    category = NA_character_,
    method = "unclassified",
    note = "No Drug Classes section and no mechanism fallback match"
  )
}

annotate_pubchem_drug_structured <- function(drug, cmap_csv) {
  cid <- pubchem_name_to_cid(drug, cmap_csv)
  if (is.na(cid)) {
    return(list(
      drug = drug,
      drugName = drug,
      pubchem_cid = NA_character_,
      drugCategory_Pubchem = NA_character_,
      annotation_source = "pubchem_structured_json",
      annotation_status = "failed",
      pubchem_url = NA_character_,
      notes = "Could not resolve PubChem CID"
    ))
  }

  pubchem_json <- fetch_pubchem_json(cid)
  classification <- classify_pubchem_json(pubchem_json)
  list(
    drug = drug,
    drugName = drug,
    pubchem_cid = cid,
    drugCategory_Pubchem = classification$category,
    annotation_source = paste("pubchem_structured_json", classification$method, sep = ":"),
    annotation_status = ifelse(is.na(classification$category), "unclassified", "refreshed"),
    pubchem_url = pubchem_pug_view_url(cid),
    notes = classification$note
  )
}

pubchem_cache_columns <- function() {
  c(
    "drug",
    "drugName",
    "pubchem_cid",
    "drugCategory_Pubchem",
    "annotation_source",
    "annotation_status",
    "retrieved_at",
    "pubchem_url",
    "notes"
  )
}

standardize_pubchem_cache <- function(cache, cache_cols = pubchem_cache_columns()) {
  if (is.null(cache) || nrow(cache) == 0) {
    cache <- data.frame(matrix(ncol = length(cache_cols), nrow = 0), stringsAsFactors = FALSE)
    colnames(cache) <- cache_cols
  }
  missing_cols <- setdiff(cache_cols, colnames(cache))
  for (col in missing_cols) {
    cache[[col]] <- NA_character_
  }
  cache <- cache[, cache_cols, drop = FALSE]
  cache$.drug_key <- toupper(cache$drug)
  cache <- cache[!duplicated(cache$.drug_key, fromLast = TRUE), , drop = FALSE]
  cache
}

merge_pubchem_refresh_results <- function(cache,
                                          requested_drugs,
                                          refreshed,
                                          subset_output = FALSE,
                                          allow_missing = FALSE,
                                          cache_cols = pubchem_cache_columns()) {
  cache <- standardize_pubchem_cache(cache, cache_cols)
  refreshed <- standardize_pubchem_cache(refreshed, cache_cols)
  requested_drugs <- unique(requested_drugs[nzchar(requested_drugs)])
  requested_keys <- toupper(requested_drugs)
  refreshed_keys <- refreshed$.drug_key
  cache_keys <- cache$.drug_key

  failed_without_cache <- refreshed$.drug_key[
    refreshed$annotation_status == "failed" &
      !refreshed$.drug_key %in% cache_keys
  ]
  missing_without_cache <- setdiff(requested_keys, union(cache_keys, refreshed_keys))
  unresolved_without_cache <- unique(c(failed_without_cache, missing_without_cache))

  metadata <- data.frame(
    parameter = c(
      "requested_drug_count",
      "refreshed_drug_count",
      "reused_requested_count",
      "preserved_cached_row_count",
      "failed_drug_count",
      "requested_missing_without_cache_count",
      "subset_output",
      "allow_missing"
    ),
    value = c(
      length(requested_keys),
      nrow(refreshed),
      sum(requested_keys %in% cache_keys & !requested_keys %in% refreshed_keys),
      sum(!cache_keys %in% refreshed_keys),
      sum(refreshed$annotation_status == "failed"),
      length(unresolved_without_cache),
      as.character(subset_output),
      as.character(allow_missing)
    ),
    stringsAsFactors = FALSE
  )

  if (length(unresolved_without_cache) > 0 && !allow_missing) {
    stop(
      "PubChem refresh could not resolve requested drugs without cached annotations: ",
      paste(requested_drugs[requested_keys %in% unresolved_without_cache], collapse = ", "),
      call. = FALSE
    )
  }

  updated <- cache
  for (i in seq_len(nrow(refreshed))) {
    key <- refreshed$.drug_key[[i]]
    hit <- match(key, updated$.drug_key)
    if (is.na(hit)) {
      updated <- rbind(updated, refreshed[i, , drop = FALSE])
    } else {
      updated[hit, ] <- refreshed[i, ]
    }
  }

  if (subset_output) {
    keep <- match(requested_keys, updated$.drug_key)
    keep <- keep[!is.na(keep)]
    updated <- updated[keep, , drop = FALSE]
  }

  list(
    cache = updated[, cache_cols, drop = FALSE],
    metadata = metadata,
    unresolved_without_cache = unresolved_without_cache
  )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
