#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  hit <- grep(needle, args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub(needle, "", hit[1], fixed = TRUE))))
  }
  normalizePath(getwd())
}

arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[1], fixed = TRUE)
}

arg_flag <- function(name) {
  paste0("--", name) %in% commandArgs(trailingOnly = TRUE)
}

read_tsv <- function(path) {
  read.table(path, sep = "\t", header = TRUE, quote = "", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write.table(x, file = path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
}

base_dir <- script_dir()
local_lib <- file.path(base_dir, ".Rlibs")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

existing_cache <- normalizePath(
  arg_value("existing-cache", file.path(base_dir, "data", "derived", "pubchem_drug_annotations.tsv")),
  mustWork = FALSE
)
output_cache <- normalizePath(
  arg_value("output-cache", file.path(base_dir, "data", "derived", "pubchem_drug_annotations.refresh.tsv")),
  mustWork = FALSE
)
cmap_file <- normalizePath(
  arg_value("cmap-file", file.path(base_dir, "data", "small_molecule_20200407234909.csv")),
  mustWork = FALSE
)
drug_file <- arg_value("drug-file", NULL)
log_file <- normalizePath(
  arg_value("log-file", file.path(base_dir, "data", "derived", "pubchem_refresh_log.tsv")),
  mustWork = FALSE
)
limit <- as.integer(arg_value("limit", NA_character_))
force_refresh <- arg_flag("force-refresh")

if (!file.exists(cmap_file)) {
  stop("Missing cmap file: ", cmap_file, call. = FALSE)
}
if (!file.exists(existing_cache) && is.null(drug_file)) {
  stop("Provide --drug-file or an existing cache at: ", existing_cache, call. = FALSE)
}

source(file.path(base_dir, "pubchem_client.R"))

cache_cols <- c(
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

if (file.exists(existing_cache)) {
  cache <- read_tsv(existing_cache)
} else {
  cache <- data.frame(matrix(ncol = length(cache_cols), nrow = 0), stringsAsFactors = FALSE)
  colnames(cache) <- cache_cols
}

missing_cols <- setdiff(cache_cols, colnames(cache))
for (col in missing_cols) {
  cache[[col]] <- NA_character_
}
cache <- cache[, cache_cols, drop = FALSE]

if (!is.null(drug_file)) {
  if (!file.exists(drug_file)) {
    stop("Missing drug file: ", drug_file, call. = FALSE)
  }
  drug_table <- read_tsv(drug_file)
  if ("drug" %in% colnames(drug_table)) {
    requested_drugs <- drug_table$drug
  } else {
    requested_drugs <- drug_table[[1]]
  }
} else {
  requested_drugs <- cache$drug
}
requested_drugs <- unique(requested_drugs[nzchar(requested_drugs)])
if (!is.na(limit)) {
  requested_drugs <- head(requested_drugs, limit)
}

cache$.drug_key <- toupper(cache$drug)
requested_keys <- toupper(requested_drugs)
cached_keys <- cache$.drug_key
refresh_keys <- if (force_refresh) requested_keys else setdiff(requested_keys, cached_keys)
refresh_drugs <- requested_drugs[requested_keys %in% refresh_keys]

log_rows <- data.frame(
  drug = requested_drugs,
  action = ifelse(requested_keys %in% refresh_keys, "refresh", "reuse"),
  status = "pending",
  message = "",
  stringsAsFactors = FALSE
)

refreshed <- data.frame(matrix(ncol = length(cache_cols), nrow = 0), stringsAsFactors = FALSE)
colnames(refreshed) <- cache_cols

if (length(refresh_drugs) > 0) {
  for (drug in refresh_drugs) {
    annotated <- tryCatch(
      annotate_pubchem_drug_structured(drug, cmap_file),
      error = function(e) {
        list(
          drug = drug,
          drugName = drug,
          pubchem_cid = NA_character_,
          drugCategory_Pubchem = NA_character_,
          annotation_source = "pubchem_structured_json",
          annotation_status = "failed",
          pubchem_url = NA_character_,
          notes = conditionMessage(e)
        )
      }
    )
    refreshed <- rbind(refreshed, data.frame(
      drug = annotated$drug,
      drugName = annotated$drugName,
      pubchem_cid = annotated$pubchem_cid,
      drugCategory_Pubchem = annotated$drugCategory_Pubchem,
      annotation_source = annotated$annotation_source,
      annotation_status = annotated$annotation_status,
      retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      pubchem_url = annotated$pubchem_url,
      notes = annotated$notes,
      stringsAsFactors = FALSE
    ))
    log_rows$status[log_rows$drug == drug] <- annotated$annotation_status
    log_rows$message[log_rows$drug == drug] <- annotated$notes
  }
}

reuse <- cache[cache$.drug_key %in% requested_keys & !cache$.drug_key %in% refresh_keys, cache_cols, drop = FALSE]
log_rows$status[log_rows$action == "reuse"] <- "reused"
log_rows$message[log_rows$action == "reuse"] <- "existing cache row reused"

combined <- rbind(reuse, refreshed)
combined$.drug_key <- toupper(combined$drug)
combined <- combined[match(requested_keys[requested_keys %in% combined$.drug_key], combined$.drug_key), cache_cols, drop = FALSE]

write_tsv(combined, output_cache)
write_tsv(log_rows, log_file)

message("Annotation refresh wrote cache to: ", output_cache)
message("Annotation refresh wrote log to: ", log_file)
