#!/usr/bin/env Rscript

required_packages <- c("cloneid", "DBI", "rJava")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(cloneid)
  library(DBI)
  library(rJava)
})

# Population ploidy is computed as a Perspective.size-weighted mean of
# clone-level total ploidy. Total clone ploidy is cloneid::calcPloidy() on
# chromosomes 1-22 plus the residual marker-chromosome ploidy stored in the
# GenomePerspective 9999 row. The Java-decoded root SP_1.0 profile is excluded,
# matching the GenomePerspective plotting code in cloneid_physicell.
#
# If the requested dose-response ID has no exact GenomePerspective record, the
# script resolves the closest same-cell-line Passaging ID with GenomePerspective
# support using the same passaged_from_id1 ancestry definition used by
# cloneid::getPedigreeTree() and cloneid::pedigree_dist().
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_path <- if (length(file_arg) > 0) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
} else {
  normalizePath("Code/in-vitro/drug_response/query_cloneid_fig3h_ploidy.R", mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(script_path), "../../.."), mustWork = TRUE)

default_lineage_map <- file.path(repo_root, "Data/in-vitro/drug_response/fig3h_cloneid_lineage_map.tsv")
default_output <- file.path(repo_root, "Data/in-vitro/drug_response/fig3h_cloneid_ploidy.tsv")
default_details_output <- file.path(repo_root, "Data/in-vitro/drug_response/fig3h_cloneid_ploidy_profile_details.tsv")
default_mapping_output <- file.path(repo_root, "Data/in-vitro/drug_response/fig3h_cloneid_id_mapping.tsv")

parse_args <- function(args) {
  out <- list(
    lineage_map = default_lineage_map,
    output = default_output,
    details_output = default_details_output,
    mapping_output = default_mapping_output,
    check_config = FALSE
  )
  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg == "--lineage-map") {
      i <- i + 1
      out$lineage_map <- args[[i]]
    } else if (arg == "--output") {
      i <- i + 1
      out$output <- args[[i]]
    } else if (arg == "--details-output") {
      i <- i + 1
      out$details_output <- args[[i]]
    } else if (arg == "--mapping-output") {
      i <- i + 1
      out$mapping_output <- args[[i]]
    } else if (arg == "--check-config") {
      out$check_config <- TRUE
    } else {
      stop(sprintf("Unknown argument: %s", arg), call. = FALSE)
    }
    i <- i + 1
  }
  out
}

read_lineage_map <- function(path) {
  lineage_map <- read.delim(
    path,
    sep = "\t",
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("")
  )
  required <- c("fig3h_group", "ploidy_class", "dose_response_label", "cloneid_origin", "mapping_status", "notes")
  missing <- setdiff(required, names(lineage_map))
  if (length(missing) > 0) {
    stop("Lineage map is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  lineage_map$cloneid_origin <- ifelse(is.na(lineage_map$cloneid_origin), "", lineage_map$cloneid_origin)
  lineage_map
}

extract_clone_id <- function(profile_labels) {
  clone_ids <- sub("^.*_ID", "", profile_labels)
  clone_ids[clone_ids == profile_labels] <- NA_character_
  clone_ids
}

extract_profile_size <- function(profile_labels) {
  matched <- regexec("^SP_([0-9.]+)_ID", profile_labels)
  parts <- regmatches(profile_labels, matched)
  vapply(parts, function(part) {
    if (length(part) < 2) {
      return(NA_real_)
    }
    suppressWarnings(as.numeric(part[[2]]))
  }, numeric(1))
}

decode_origin_matrix <- function(origin) {
  perspective_enum <- .jcall(
    "core.utils.Perspectives",
    "Lcore/utils/Perspectives;",
    "valueOf",
    "GenomePerspective"
  )
  profile_map <- .jcall(
    "cloneid.Manager",
    returnSig = "Ljava/util/Map;",
    method = "profiles",
    origin,
    perspective_enum,
    TRUE
  )
  cloneid:::.javamap2Rmatrix(profile_map)
}

parse_genome_perspective_segments <- function(segment_labels) {
  chrom <- suppressWarnings(as.integer(sub(":.*$", "", segment_labels)))
  coords <- sub("^[^:]+:", "", segment_labels)
  start <- suppressWarnings(as.numeric(sub("-.*$", "", coords)))
  end <- suppressWarnings(as.numeric(sub("^.*-", "", coords)))
  width <- end - start + 1
  keep <- !is.na(chrom) & chrom >= 1 & chrom <= 22 & is.finite(width) & width > 0
  data.frame(
    label = segment_labels,
    chrom = chrom,
    start = start,
    end = end,
    width = width,
    keep = keep,
    stringsAsFactors = FALSE
  )
}

clone_ploidy_from_profile_matrix <- function(matrix_24_by_n, origin) {
  segment_info <- parse_genome_perspective_segments(rownames(matrix_24_by_n))
  autosome_keep <- segment_info$keep
  if (!any(autosome_keep)) {
    stop("No autosome segment rows were decoded for origin: ", origin, call. = FALSE)
  }

  keep_columns <- !startsWith(colnames(matrix_24_by_n), "SP_1.0_")
  matrix_no_root <- matrix_24_by_n[, keep_columns, drop = FALSE]
  if (ncol(matrix_no_root) == 0) {
    stop("No non-root GenomePerspective clone profiles were decoded for origin: ", origin, call. = FALSE)
  }

  matrix_22_by_n <- matrix_no_root[autosome_keep, , drop = FALSE]
  autosome_info <- segment_info[autosome_keep, , drop = FALSE]
  whole_chromosome_copy_number <- t(matrix_22_by_n)
  colnames(whole_chromosome_copy_number) <- as.character(autosome_info$chrom)
  chrwhole <- matrix(autosome_info$width, ncol = 1)
  rownames(chrwhole) <- paste0("chr", autosome_info$chrom)

  autosomal_ploidy <- cloneid::calcPloidy(whole_chromosome_copy_number, chrwhole = chrwhole)
  marker_count_label <- "999:1-999"
  marker_ploidy_label <- "9999:1-999"
  marker_chromosome_count <- if (marker_count_label %in% rownames(matrix_no_root)) {
    as.numeric(matrix_no_root[marker_count_label, ])
  } else {
    rep(NA_real_, ncol(matrix_no_root))
  }
  marker_chromosome_ploidy <- if (marker_ploidy_label %in% rownames(matrix_no_root)) {
    as.numeric(matrix_no_root[marker_ploidy_label, ])
  } else {
    rep(0, ncol(matrix_no_root))
  }
  marker_chromosome_ploidy[!is.finite(marker_chromosome_ploidy)] <- 0
  total_ploidy <- as.numeric(autosomal_ploidy) + marker_chromosome_ploidy

  profile_labels <- colnames(matrix_no_root)

  data.frame(
    origin = origin,
    profile_label = profile_labels,
    cloneID = extract_clone_id(profile_labels),
    profile_size = extract_profile_size(profile_labels),
    autosomal_ploidy = as.numeric(autosomal_ploidy),
    marker_chromosome_count = marker_chromosome_count,
    marker_chromosome_ploidy = marker_chromosome_ploidy,
    clone_ploidy = total_ploidy,
    stringsAsFactors = FALSE
  )
}

fetch_perspective_records <- function(con, origin) {
  sql <- paste0(
    "SELECT cloneID, origin, whichPerspective, size, parent, sampleSource, rootID, state, alias, hasChildren ",
    "FROM Perspective ",
    "WHERE whichPerspective = ", DBI::dbQuoteString(con, "GenomePerspective"),
    " AND origin = ", DBI::dbQuoteString(con, origin),
    " ORDER BY hasChildren DESC, cloneID"
  )
  DBI::dbGetQuery(con, sql)
}

collapse_clean_values <- function(x) {
  x <- unique(stats::na.omit(as.character(x)))
  x <- x[nzchar(x) & tolower(x) != "null"]
  if (length(x) == 0) {
    return(NA_character_)
  }
  paste(x, collapse = ";")
}

weighted_quantile <- function(x, w, probs) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  x <- x[ok]
  w <- w[ok]
  if (length(x) == 0) {
    return(rep(NA_real_, length(probs)))
  }
  ord <- order(x)
  x <- x[ord]
  w <- w[ord]
  cumulative_weight <- cumsum(w) / sum(w)
  vapply(probs, function(prob) x[which(cumulative_weight >= prob)[[1]]], numeric(1))
}

fetch_genome_perspective_aliases <- function(con, origin) {
  if (is.na(origin) || !nzchar(origin)) {
    return(NA_character_)
  }
  sql <- paste0(
    "SELECT DISTINCT alias FROM Perspective ",
    "WHERE whichPerspective = ", DBI::dbQuoteString(con, "GenomePerspective"),
    " AND origin = ", DBI::dbQuoteString(con, origin),
    " AND alias IS NOT NULL ",
    "ORDER BY alias"
  )
  aliases <- try(DBI::dbGetQuery(con, sql), silent = TRUE)
  if (inherits(aliases, "try-error") || nrow(aliases) == 0) {
    return(NA_character_)
  }
  collapse_clean_values(aliases$alias)
}

fetch_passaging_record <- function(con, passaging_id) {
  sql <- paste0(
    "SELECT id, cellLine, event, passaged_from_id1, passaged_from_id2, passage, date ",
    "FROM Passaging WHERE id = ", DBI::dbQuoteString(con, passaging_id)
  )
  DBI::dbGetQuery(con, sql)
}

fetch_genome_perspective_origins_for_cell_line <- function(con, cell_line) {
  sql <- paste0(
    "SELECT DISTINCT p.id AS origin, p.event, p.passage, p.date ",
    "FROM Passaging p ",
    "INNER JOIN Perspective pr ON pr.origin = p.id ",
    "WHERE pr.whichPerspective = ", DBI::dbQuoteString(con, "GenomePerspective"),
    " AND p.cellLine = ", DBI::dbQuoteString(con, cell_line),
    " ORDER BY p.date, p.id"
  )
  DBI::dbGetQuery(con, sql)
}

cell_line_passaging_cache <- new.env(parent = emptyenv())

fetch_passaging_records_for_cell_line <- function(con, cell_line) {
  cache_key <- as.character(cell_line)
  if (exists(cache_key, envir = cell_line_passaging_cache, inherits = FALSE)) {
    return(get(cache_key, envir = cell_line_passaging_cache, inherits = FALSE))
  }
  sql <- paste0(
    "SELECT id, cellLine, event, passaged_from_id1, passaged_from_id2, passage, date ",
    "FROM Passaging WHERE cellLine = ", DBI::dbQuoteString(con, cell_line)
  )
  records <- DBI::dbGetQuery(con, sql)
  records$id <- as.character(records$id)
  records$passaged_from_id1 <- as.character(records$passaged_from_id1)
  records$passaged_from_id1[is.na(records$passaged_from_id1) | records$passaged_from_id1 == "NA"] <- NA_character_
  assign(cache_key, records, envir = cell_line_passaging_cache)
  records
}

pedigree_ancestors_from_passaging <- function(passaging_records, start_id) {
  if (!start_id %in% passaging_records$id) {
    return(character())
  }
  ancestors <- character()
  current_id <- start_id
  seen <- character()
  while (!is.na(current_id) && nzchar(current_id) && !current_id %in% seen) {
    ancestors <- c(ancestors, current_id)
    seen <- c(seen, current_id)
    idx <- match(current_id, passaging_records$id)
    if (is.na(idx)) {
      break
    }
    current_id <- passaging_records$passaged_from_id1[[idx]]
  }
  ancestors
}

pedigree_distance_edges <- function(passaging_records, target_id, candidate_id) {
  target_ancestors <- pedigree_ancestors_from_passaging(passaging_records, target_id)
  candidate_ancestors <- pedigree_ancestors_from_passaging(passaging_records, candidate_id)
  common_ancestors <- intersect(target_ancestors, candidate_ancestors)
  mrca <- if (length(common_ancestors) > 0) common_ancestors[[1]] else NA_character_
  if (is.na(mrca) || length(target_ancestors) == 0 || length(candidate_ancestors) == 0) {
    return(list(distance_edges = Inf, mrca = NA_character_, target_depth = NA_integer_, candidate_depth = NA_integer_))
  }
  target_depth <- match(mrca, target_ancestors) - 1L
  candidate_depth <- match(mrca, candidate_ancestors) - 1L
  if (is.na(target_depth) || is.na(candidate_depth)) {
    return(list(distance_edges = Inf, mrca = mrca, target_depth = NA_integer_, candidate_depth = NA_integer_))
  }
  list(
    distance_edges = as.numeric(target_depth + candidate_depth),
    mrca = mrca,
    target_depth = as.integer(target_depth),
    candidate_depth = as.integer(candidate_depth)
  )
}

resolve_analysis_origin <- function(con, lineage_row) {
  target_id <- lineage_row$dose_response_label
  requested_origin <- lineage_row$cloneid_origin
  if (nzchar(requested_origin) && nrow(fetch_perspective_records(con, requested_origin)) > 0) {
    target_record <- fetch_passaging_record(con, requested_origin)
    target_cell_line <- if (nrow(target_record) > 0) as.character(target_record$cellLine[[1]]) else NA_character_
    return(list(
      analysis_origin = requested_origin,
      target_passaging_id = target_id,
      target_cellLine = target_cell_line,
      target_event = if (nrow(target_record) > 0) as.character(target_record$event[[1]]) else NA_character_,
      target_passage = if (nrow(target_record) > 0) suppressWarnings(as.integer(target_record$passage[[1]])) else NA_integer_,
      target_date = if (nrow(target_record) > 0) as.character(target_record$date[[1]]) else NA_character_,
      analysis_event = if (nrow(target_record) > 0) as.character(target_record$event[[1]]) else NA_character_,
      analysis_passage = if (nrow(target_record) > 0) suppressWarnings(as.integer(target_record$passage[[1]])) else NA_integer_,
      analysis_date = if (nrow(target_record) > 0) as.character(target_record$date[[1]]) else NA_character_,
      origin_resolution_status = "exact_mapped_genome_perspective",
      origin_resolution_distance_edges = 0,
      origin_resolution_mrca = requested_origin,
      origin_resolution_target_depth = 0L,
      origin_resolution_candidate_depth = 0L,
      origin_resolution_candidate_count = 1L,
      origin_resolution_note = "Used exact mapped CLONEID GenomePerspective origin."
    ))
  }

  target_record <- fetch_passaging_record(con, target_id)
  if (nrow(target_record) == 0) {
    return(list(
      analysis_origin = "",
      target_passaging_id = target_id,
      target_cellLine = NA_character_,
      target_event = NA_character_,
      target_passage = NA_integer_,
      target_date = NA_character_,
      analysis_event = NA_character_,
      analysis_passage = NA_integer_,
      analysis_date = NA_character_,
      origin_resolution_status = "target_not_in_passaging",
      origin_resolution_distance_edges = NA_real_,
      origin_resolution_mrca = NA_character_,
      origin_resolution_target_depth = NA_integer_,
      origin_resolution_candidate_depth = NA_integer_,
      origin_resolution_candidate_count = 0L,
      origin_resolution_note = "Dose-response label is not a Passaging.id, so no pedigree-nearest GenomePerspective origin could be selected."
    ))
  }

  cell_line <- as.character(target_record$cellLine[[1]])
  candidate_origins <- fetch_genome_perspective_origins_for_cell_line(con, cell_line)
  if (nrow(candidate_origins) == 0) {
    return(list(
      analysis_origin = "",
      target_passaging_id = target_id,
      target_cellLine = cell_line,
      target_event = as.character(target_record$event[[1]]),
      target_passage = suppressWarnings(as.integer(target_record$passage[[1]])),
      target_date = as.character(target_record$date[[1]]),
      analysis_event = NA_character_,
      analysis_passage = NA_integer_,
      analysis_date = NA_character_,
      origin_resolution_status = "no_cell_line_genome_perspective_candidates",
      origin_resolution_distance_edges = NA_real_,
      origin_resolution_mrca = NA_character_,
      origin_resolution_target_depth = NA_integer_,
      origin_resolution_candidate_depth = NA_integer_,
      origin_resolution_candidate_count = 0L,
      origin_resolution_note = "No same-cell-line Passaging IDs have GenomePerspective records."
    ))
  }
  passaging_records <- fetch_passaging_records_for_cell_line(con, cell_line)

  distances <- do.call(rbind, lapply(seq_len(nrow(candidate_origins)), function(i) {
    candidate_id <- as.character(candidate_origins$origin[[i]])
    d <- pedigree_distance_edges(passaging_records, target_id, candidate_id)
    target_date <- suppressWarnings(as.POSIXct(target_record$date[[1]]))
    candidate_date <- suppressWarnings(as.POSIXct(candidate_origins$date[[i]]))
    date_delta_days <- if (is.na(target_date) || is.na(candidate_date)) {
      Inf
    } else {
      abs(as.numeric(difftime(target_date, candidate_date, units = "days")))
    }
    data.frame(
      candidate_origin = candidate_id,
      candidate_event = as.character(candidate_origins$event[[i]]),
      candidate_passage = suppressWarnings(as.integer(candidate_origins$passage[[i]])),
      candidate_date = as.character(candidate_origins$date[[i]]),
      distance_edges = d$distance_edges,
      mrca = d$mrca,
      target_depth = d$target_depth,
      candidate_depth = d$candidate_depth,
      date_delta_days = date_delta_days,
      stringsAsFactors = FALSE
    )
  }))

  finite_distances <- distances[is.finite(distances$distance_edges), , drop = FALSE]
  if (nrow(finite_distances) == 0) {
    return(list(
      analysis_origin = "",
      target_passaging_id = target_id,
      target_cellLine = cell_line,
      target_event = as.character(target_record$event[[1]]),
      target_passage = suppressWarnings(as.integer(target_record$passage[[1]])),
      target_date = as.character(target_record$date[[1]]),
      analysis_event = NA_character_,
      analysis_passage = NA_integer_,
      analysis_date = NA_character_,
      origin_resolution_status = "no_connected_genome_perspective_candidate",
      origin_resolution_distance_edges = NA_real_,
      origin_resolution_mrca = NA_character_,
      origin_resolution_target_depth = NA_integer_,
      origin_resolution_candidate_depth = NA_integer_,
      origin_resolution_candidate_count = nrow(candidate_origins),
      origin_resolution_note = "Same-cell-line GenomePerspective candidates exist, but none were connected by the cloneid pedigree ancestor/MRCA functions."
    ))
  }

  finite_distances <- finite_distances[order(
    finite_distances$distance_edges,
    finite_distances$date_delta_days,
    finite_distances$candidate_origin
  ), , drop = FALSE]
  closest <- finite_distances[1, , drop = FALSE]
  list(
    analysis_origin = closest$candidate_origin[[1]],
    target_passaging_id = target_id,
    target_cellLine = cell_line,
    target_event = as.character(target_record$event[[1]]),
    target_passage = suppressWarnings(as.integer(target_record$passage[[1]])),
    target_date = as.character(target_record$date[[1]]),
    analysis_event = closest$candidate_event[[1]],
    analysis_passage = closest$candidate_passage[[1]],
    analysis_date = closest$candidate_date[[1]],
    origin_resolution_status = "nearest_pedigree_genome_perspective",
    origin_resolution_distance_edges = closest$distance_edges[[1]],
    origin_resolution_mrca = closest$mrca[[1]],
    origin_resolution_target_depth = closest$target_depth[[1]],
    origin_resolution_candidate_depth = closest$candidate_depth[[1]],
    origin_resolution_candidate_count = nrow(candidate_origins),
    origin_resolution_note = sprintf(
      "No exact GenomePerspective origin for requested ID; using closest same-cell-line GenomePerspective origin %s by passaged_from_id1 pedigree distance (%s edge(s), MRCA %s), matching the lineage edge used by cloneid::getPedigreeTree()/pedigree_dist(). Candidate event=%s, passage=%s, date=%s.",
      closest$candidate_origin[[1]],
      closest$distance_edges[[1]],
      closest$mrca[[1]],
      closest$candidate_event[[1]],
      closest$candidate_passage[[1]],
      closest$candidate_date[[1]]
    )
  )
}

summarize_one_origin <- function(con, lineage_row) {
  resolved <- resolve_analysis_origin(con, lineage_row)
  base <- data.frame(
    fig3h_group = lineage_row$fig3h_group,
    ploidy_class = lineage_row$ploidy_class,
    dose_response_label = lineage_row$dose_response_label,
    cloneid_origin = lineage_row$cloneid_origin,
    analysis_cloneid_origin = resolved$analysis_origin,
    target_genome_perspective_alias = NA_character_,
    analysis_genome_perspective_alias = NA_character_,
    target_passaging_id = resolved$target_passaging_id,
    target_cellLine = resolved$target_cellLine,
    target_event = resolved$target_event,
    target_passage = resolved$target_passage,
    target_date = resolved$target_date,
    analysis_event = resolved$analysis_event,
    analysis_passage = resolved$analysis_passage,
    analysis_date = resolved$analysis_date,
    origin_resolution_status = resolved$origin_resolution_status,
    origin_resolution_distance_edges = resolved$origin_resolution_distance_edges,
    origin_resolution_mrca = resolved$origin_resolution_mrca,
    origin_resolution_target_depth = resolved$origin_resolution_target_depth,
    origin_resolution_candidate_depth = resolved$origin_resolution_candidate_depth,
    origin_resolution_candidate_count = resolved$origin_resolution_candidate_count,
    mapping_status = lineage_row$mapping_status,
    mean_population_ploidy = NA_real_,
    median_population_ploidy = NA_real_,
    q90_population_ploidy = NA_real_,
    q10_population_ploidy = NA_real_,
    q75_population_ploidy = NA_real_,
    max_population_ploidy = NA_real_,
    min_population_ploidy = NA_real_,
    unweighted_mean_clone_ploidy = NA_real_,
    n_perspective_records = NA_integer_,
    n_leaf_clones = NA_integer_,
    leaf_size_sum = NA_real_,
    sampleSource = NA_character_,
    rootID = NA_character_,
    ploidy_status = NA_character_,
    notes = paste(lineage_row$notes, resolved$origin_resolution_note, sep = " "),
    stringsAsFactors = FALSE
  )

  if (!nzchar(resolved$analysis_origin)) {
    base$ploidy_status <- resolved$origin_resolution_status
    return(list(summary = base, details = data.frame()))
  }

  perspective <- fetch_perspective_records(con, resolved$analysis_origin)
  if (nrow(perspective) == 0) {
    base$n_perspective_records <- 0L
    base$n_leaf_clones <- 0L
    base$ploidy_status <- "no_exact_genome_perspective"
    return(list(summary = base, details = data.frame()))
  }
  base$analysis_genome_perspective_alias <- collapse_clean_values(perspective$alias)
  if (identical(resolved$target_passaging_id, resolved$analysis_origin)) {
    base$target_genome_perspective_alias <- base$analysis_genome_perspective_alias
  }

  matrix_24_by_n <- try(decode_origin_matrix(resolved$analysis_origin), silent = TRUE)
  if (inherits(matrix_24_by_n, "try-error")) {
    base$n_perspective_records <- nrow(perspective)
    base$ploidy_status <- "profile_decode_failed"
    base$notes <- paste(base$notes, as.character(matrix_24_by_n), sep = " ")
    return(list(summary = base, details = data.frame()))
  }

  details <- clone_ploidy_from_profile_matrix(matrix_24_by_n, resolved$analysis_origin)
  details$profile_order <- seq_len(nrow(details))

  has_children <- suppressWarnings(as.integer(as.character(perspective$hasChildren)))
  perspective$cloneID <- as.character(perspective$cloneID)
  perspective$size <- suppressWarnings(as.numeric(perspective$size))
  leaf_perspective <- perspective[!is.na(has_children) & has_children == 0, , drop = FALSE]
  db_leaf <- leaf_perspective[, c("cloneID", "size", "state", "alias", "parent", "sampleSource", "rootID"), drop = FALSE]
  names(db_leaf)[names(db_leaf) == "size"] <- "db_size"

  details <- merge(details, db_leaf, by = "cloneID", all.x = TRUE, sort = FALSE)
  details <- details[order(details$profile_order), , drop = FALSE]
  details$weight_used <- ifelse(is.finite(details$db_size) & details$db_size > 0, details$db_size, details$profile_size)
  details$fig3h_group <- lineage_row$fig3h_group
  details$ploidy_class <- lineage_row$ploidy_class
  details$dose_response_label <- lineage_row$dose_response_label
  details$analysis_cloneid_origin <- resolved$analysis_origin
  details$origin_resolution_status <- resolved$origin_resolution_status
  details$origin_resolution_distance_edges <- resolved$origin_resolution_distance_edges
  details <- details[, c(
    "fig3h_group", "ploidy_class", "dose_response_label", "analysis_cloneid_origin",
    "origin_resolution_status", "origin_resolution_distance_edges", "origin", "cloneID", "profile_label",
    "profile_size", "db_size", "weight_used", "autosomal_ploidy",
    "marker_chromosome_count", "marker_chromosome_ploidy", "clone_ploidy",
    "state", "alias", "parent",
    "sampleSource", "rootID"
  ), drop = FALSE]

  valid <- is.finite(details$clone_ploidy) & is.finite(details$weight_used) & details$weight_used > 0
  if (!any(valid)) {
    base$ploidy_status <- "no_valid_clone_weights"
  } else {
    base$mean_population_ploidy <- stats::weighted.mean(details$clone_ploidy[valid], details$weight_used[valid])
    ploidy_quantiles <- weighted_quantile(
      details$clone_ploidy[valid],
      details$weight_used[valid],
      probs = c(0.10, 0.50, 0.75, 0.90)
    )
    base$q10_population_ploidy <- ploidy_quantiles[[1]]
    base$median_population_ploidy <- ploidy_quantiles[[2]]
    base$q75_population_ploidy <- ploidy_quantiles[[3]]
    base$q90_population_ploidy <- ploidy_quantiles[[4]]
    base$min_population_ploidy <- min(details$clone_ploidy[valid], na.rm = TRUE)
    base$max_population_ploidy <- max(details$clone_ploidy[valid], na.rm = TRUE)
    base$unweighted_mean_clone_ploidy <- mean(details$clone_ploidy[is.finite(details$clone_ploidy)])
    base$ploidy_status <- "ok"
  }
  base$n_perspective_records <- nrow(perspective)
  base$n_leaf_clones <- nrow(leaf_perspective)
  base$leaf_size_sum <- sum(leaf_perspective$size, na.rm = TRUE)
  base$sampleSource <- paste(unique(stats::na.omit(as.character(perspective$sampleSource))), collapse = ";")
  leaf_root_ids <- unique(stats::na.omit(as.character(leaf_perspective$rootID)))
  leaf_root_ids <- leaf_root_ids[nzchar(leaf_root_ids) & leaf_root_ids != "0"]
  base$rootID <- paste(leaf_root_ids, collapse = ";")

  list(summary = base, details = details)
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  lineage_map <- read_lineage_map(args$lineage_map)

  if (isTRUE(args$check_config)) {
    message("Lineage map rows: ", nrow(lineage_map))
    message("Mapped CLONEID origins: ", sum(nzchar(lineage_map$cloneid_origin)))
    message("Unresolved rows: ", sum(!nzchar(lineage_map$cloneid_origin)))
    message("Required packages available: ", paste(required_packages, collapse = ", "))
    return(invisible(NULL))
  }

  con <- cloneid::connect2DB()
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  results <- lapply(seq_len(nrow(lineage_map)), function(i) summarize_one_origin(con, lineage_map[i, , drop = FALSE]))
  summary <- do.call(rbind, lapply(results, `[[`, "summary"))
  details_list <- lapply(results, `[[`, "details")
  details <- if (length(details_list) == 0 || all(vapply(details_list, nrow, integer(1)) == 0)) {
    data.frame()
  } else {
    do.call(rbind, details_list[vapply(details_list, nrow, integer(1)) > 0])
  }
  id_mapping <- summary[, c(
    "fig3h_group", "ploidy_class", "dose_response_label", "target_passaging_id",
    "target_cellLine", "target_event", "target_passage", "target_date",
    "cloneid_origin", "target_genome_perspective_alias",
    "analysis_cloneid_origin", "analysis_genome_perspective_alias",
    "analysis_event", "analysis_passage", "analysis_date",
    "origin_resolution_status", "origin_resolution_distance_edges",
    "origin_resolution_mrca", "origin_resolution_target_depth",
    "origin_resolution_candidate_depth", "origin_resolution_candidate_count",
    "mean_population_ploidy", "median_population_ploidy",
    "q90_population_ploidy", "q10_population_ploidy",
    "q75_population_ploidy", "max_population_ploidy",
    "min_population_ploidy", "ploidy_status", "notes"
  ), drop = FALSE]
  names(id_mapping)[names(id_mapping) == "cloneid_origin"] <- "requested_cloneid_origin"
  names(id_mapping)[names(id_mapping) == "target_genome_perspective_alias"] <- "target_alias"
  names(id_mapping)[names(id_mapping) == "analysis_cloneid_origin"] <- "closest_genome_perspective_id"
  names(id_mapping)[names(id_mapping) == "analysis_genome_perspective_alias"] <- "closest_alias"
  id_mapping$exact_genome_perspective_match <- id_mapping$origin_resolution_status == "exact_mapped_genome_perspective"
  id_mapping$closest_id_used_for_ploidy <- id_mapping$ploidy_status == "ok"
  id_mapping <- id_mapping[, c(
    "fig3h_group", "ploidy_class", "dose_response_label", "target_passaging_id",
    "target_cellLine", "target_event", "target_passage", "target_date",
    "requested_cloneid_origin", "target_alias", "closest_genome_perspective_id",
    "closest_alias", "exact_genome_perspective_match", "closest_id_used_for_ploidy",
    "analysis_event", "analysis_passage", "analysis_date",
    "origin_resolution_status", "origin_resolution_distance_edges",
    "origin_resolution_mrca", "origin_resolution_target_depth",
    "origin_resolution_candidate_depth", "origin_resolution_candidate_count",
    "mean_population_ploidy", "median_population_ploidy",
    "q90_population_ploidy", "q10_population_ploidy",
    "q75_population_ploidy", "max_population_ploidy",
    "min_population_ploidy", "ploidy_status", "notes"
  ), drop = FALSE]

  dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
  write.table(summary, file = args$output, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(details, file = args$details_output, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(id_mapping, file = args$mapping_output, sep = "\t", quote = FALSE, row.names = FALSE)

  message("Wrote ", args$output)
  message("Wrote ", args$details_output)
  message("Wrote ", args$mapping_output)
}

main()
