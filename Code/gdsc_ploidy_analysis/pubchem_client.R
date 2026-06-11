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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}
