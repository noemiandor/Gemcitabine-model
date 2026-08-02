# Input and TGI preparation adapted from the reviewed 04h workflow at 77caec93.

figure7_numeric <- function(x) suppressWarnings(as.numeric(x))

figure7_curated_endpoint_counts <- function() {
  c(
    "2N-A1-0" = 413L,
    "2N-A1-LR" = 369L,
    "2N-A1-R" = 196L,
    "2N-A1-RR" = 1495L,
    "2N-A2-0" = 317L,
    "2N-A2-L" = 505L,
    "2N-A4-R" = 305L,
    "2N-A4-RL" = 1280L,
    "4N-A5-0" = 888L,
    "4N-A5-RR" = 393L,
    "A5-4N-L" = 385L,
    "A5-4N-R" = 358L,
    "A6-4N-O" = 189L,
    "A6-4N-RR" = 660L,
    "4N-A8-RL" = 1832L,
    "4N-A8-RR" = 247L
  )
}

figure7_unique_sample_value <- function(data, column, sample_id, numeric = FALSE) {
  values <- data[data$sample_id == sample_id, column]
  if (numeric) {
    values <- unique(figure7_numeric(values))
    values <- values[is.finite(values)]
  } else {
    values <- unique(trimws(as.character(values)))
    values <- values[!is.na(values) & nzchar(values)]
  }
  if (length(values) != 1L) {
    figure7_stop("Expected one sample-level value for ", column, " in ", sample_id,
                 "; found ", length(values))
  }
  values[[1L]]
}

figure7_read_cell_table <- function(path, compartment, config) {
  delta_measure <- figure7_tgi_delta_measure(config)
  tgi_measure <- figure7_tgi_measure(config)
  required <- c(
    "cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy",
    "growth_curve_harvest",
    "tumor_volume_baseline_day", "tumor_volume_baseline",
    delta_measure, tgi_measure
  )
  if (!file.exists(path)) figure7_stop("Missing ", compartment, " input: ", path)
  data <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(data))
  if (length(missing)) figure7_stop(compartment, " input is missing: ", paste(missing, collapse = ", "))
  numeric_columns <- unique(c(
    "gemcitabine_dose_mg_per_kg", "pseudotime", "cell_ploidy",
    "tumor_volume_baseline", delta_measure, tgi_measure,
    grep("^tumor_volume_Day_[0-9]+$", names(data), value = TRUE)
  ))
  for (column in intersect(numeric_columns, names(data))) data[[column]] <- figure7_numeric(data[[column]])
  if (any(!data$initial_ploidy %in% c("2N", "4N"))) figure7_stop("Unexpected initial-ploidy label in ", compartment)
  data$compartment <- compartment
  data
}

figure7_read_endpoint_ploidy_table <- function(path, config) {
  if (!file.exists(path)) {
    figure7_stop("Missing canonical endpoint-ploidy input: ", path)
  }
  expected_hash <- as.character(
    config$versioned_source_artifacts$panel_l_endpoint_ploidy$sha256
  )
  figure7_verify_checksum(
    path,
    expected_hash,
    "canonical 14,125-cell endpoint-ploidy input"
  )
  cells <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  expected_columns <- c(
    "file", "cell_id", "ploidy", "total_chromosomes", "frac_covered",
    "format"
  )
  if (!identical(names(cells), expected_columns)) {
    figure7_stop(
      "Canonical endpoint-ploidy input must have exactly these columns: ",
      paste(expected_columns, collapse = ", ")
    )
  }
  cells$file <- trimws(as.character(cells$file))
  cells$cell_id <- trimws(as.character(cells$cell_id))
  cells$ploidy <- figure7_numeric(cells$ploidy)
  cells$total_chromosomes <- figure7_numeric(cells$total_chromosomes)
  cells$frac_covered <- figure7_numeric(cells$frac_covered)
  cells$format <- trimws(as.character(cells$format))
  keys <- paste(cells$file, cells$cell_id, sep = "\r")
  if (nrow(cells) != 14125L || length(unique(cells$file)) != 16L ||
      anyNA(cells$file) || any(!nzchar(cells$file)) ||
      any(basename(cells$file) != cells$file) ||
      any(!grepl("[.]sps[.]cbs$", cells$file)) ||
      anyNA(cells$cell_id) || any(!nzchar(cells$cell_id)) ||
      anyDuplicated(keys) || any(!is.finite(cells$ploidy)) ||
      any(cells$ploidy <= 0) ||
      any(!is.finite(cells$total_chromosomes)) ||
      any(cells$total_chromosomes <= 0) ||
      any(!is.finite(cells$frac_covered)) ||
      any(cells$frac_covered <= 0 | cells$frac_covered > 1) ||
      any(cells$format != "wide")) {
    figure7_stop(
      "Canonical endpoint-ploidy input must contain the exact complete ",
      "14,125-cell, 16-file reviewed wide-CBS collection"
    )
  }
  list(
    cells = cells,
    path = normalizePath(path, mustWork = TRUE),
    sha256 = expected_hash,
    n_cells = nrow(cells),
    n_files = length(unique(cells$file)),
    score_policy =
      paste0(
        "arithmetic_mean_of_finite_postprocessed_cell_ploidy_in_exact_",
        "qc_passed_cellcycle_noncellcycle_union_per_sample"
      ),
    mapping_policy =
      paste0(
        "exact_processed_sample_barcode_to_canonical_cbs_file_cell_and_value;",
        "score_universe=reviewed_final_seurat_tumor_cells"
      )
  )
}

figure7_sample_table <- function(
  cellcycle,
  noncellcycle,
  config,
  endpoint_ploidy
) {
  if (is.null(endpoint_ploidy$cells) ||
      !identical(endpoint_ploidy$n_cells, 14125L) ||
      !identical(endpoint_ploidy$n_files, 16L)) {
    figure7_stop(
      "Figure 7 sample preparation requires the validated complete ",
      "14,125-cell endpoint-ploidy input"
    )
  }
  delta_measure <- figure7_tgi_delta_measure(config)
  tgi_measure <- figure7_tgi_measure(config)
  embedded_tgi_measure <- paste0("embedded_", tgi_measure)
  columns <- c("cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose",
               "gemcitabine_dose_mg_per_kg", "growth_curve_harvest",
               "cell_ploidy", "compartment")
  union <- rbind(cellcycle[, columns], noncellcycle[, columns])
  union <- union[is.finite(union$cell_ploidy), , drop = FALSE]
  duplicate_ids <- unique(union$cell_id[duplicated(union$cell_id)])
  if (length(duplicate_ids)) {
    figure7_stop(
      "Figure 7 score-universe tables must be disjoint by cell_id; found ",
      length(duplicate_ids), " duplicated curated cell(s)"
    )
  }

  compartment_counts <- table(union$compartment)
  if (nrow(union) != 9832L ||
      !identical(
        as.integer(compartment_counts[c("CellCycle", "NonCellCycle")]),
        c(2881L, 6951L)
      )) {
    figure7_stop(
      "Figure 7 score universe must be the exact 9,832-cell reviewed ",
      "CellCycle + NonCellCycle tumor union (2,881 + 6,951 cells)"
    )
  }

  # The complete 14,125-cell table remains the immutable upstream CBS
  # inventory.  The Figure 7 analysis universe is narrower by design: it is
  # the exact 9,832 tumor cells retained in the reviewed final Seurat object
  # and represented by the CellCycle + NonCellCycle tables.  Validate every
  # curated key and value against the canonical inventory before aggregating
  # only those curated cells for panel L and its endpoint sensitivities.
  prefix <- paste0(union$sample_id, "_")
  prefix_matches <- startsWith(union$cell_id, prefix)
  if (any(!prefix_matches)) {
    figure7_stop(
      "Processed Figure 7 cell IDs do not preserve the sample_id_barcode contract"
    )
  }
  union$endpoint_cell_id <- substring(union$cell_id, nchar(prefix) + 1L)
  union$endpoint_ploidy_file <- paste0(
    trimws(as.character(union$growth_curve_harvest)),
    ".sps.cbs"
  )
  endpoint_cells <- endpoint_ploidy$cells
  endpoint_keys <- paste(
    endpoint_cells$file,
    endpoint_cells$cell_id,
    sep = "\r"
  )
  union_keys <- paste(
    union$endpoint_ploidy_file,
    union$endpoint_cell_id,
    sep = "\r"
  )
  endpoint_index <- match(union_keys, endpoint_keys)
  if (anyNA(endpoint_index) ||
      any(abs(union$cell_ploidy - endpoint_cells$ploidy[endpoint_index]) >
        1e-12)) {
    figure7_stop(
      "Processed Figure 7 cells do not map exactly to the canonical ",
      "sample-specific endpoint CBS values"
    )
  }

  rows <- lapply(sort(unique(union$sample_id)), function(id) {
    local <- union[union$sample_id == id, , drop = FALSE]
    cc <- cellcycle[cellcycle$sample_id == id, , drop = FALSE]
    if (!nrow(cc)) figure7_stop("Sample has no CellCycle cells: ", id)
    harvest <- figure7_unique_sample_value(
      local,
      "growth_curve_harvest",
      id
    )
    endpoint_file <- paste0(harvest, ".sps.cbs")
    canonical <- endpoint_cells[
      endpoint_cells$file == endpoint_file,
      ,
      drop = FALSE
    ]
    if (!nrow(canonical)) {
      figure7_stop(
        "No canonical endpoint CBS cells mapped to Figure 7 sample ",
        id,
        " via ",
        endpoint_file
      )
    }
    row <- data.frame(
      sample_id = id,
      initial_ploidy = figure7_unique_sample_value(local, "initial_ploidy", id),
      dose = figure7_unique_sample_value(local, "gemcitabine_dose", id),
      dose_mg = figure7_unique_sample_value(local, "gemcitabine_dose_mg_per_kg", id, TRUE),
      growth_curve_harvest = harvest,
      endpoint_ploidy_file = endpoint_file,
      sample_mean_endpoint_ploidy = mean(local$cell_ploidy),
      sample_median_endpoint_ploidy = stats::median(local$cell_ploidy),
      n_endpoint_ploidy_cells = nrow(local),
      endpoint_ploidy_source_total_cells = endpoint_ploidy$n_cells,
      endpoint_ploidy_source_file_count = endpoint_ploidy$n_files,
      endpoint_ploidy_source_sha256 = endpoint_ploidy$sha256,
      endpoint_ploidy_score_universe_total_cells = nrow(union),
      endpoint_ploidy_score_policy = endpoint_ploidy$score_policy,
      endpoint_ploidy_mapping_policy = endpoint_ploidy$mapping_policy,
      n_cellcycle_cells = nrow(cc),
      mean_cellcycle_pseudotime = mean(cc$pseudotime, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    row[[delta_measure]] <- figure7_unique_sample_value(cc, delta_measure, id, TRUE)
    row[[embedded_tgi_measure]] <- figure7_unique_sample_value(cc, tgi_measure, id, TRUE)
    row
  })
  samples <- do.call(rbind, rows)
  file_match <- regexec(
    "^SUM159-(2N|4N)-([0-9]+)-.+_harvest[.]sps[.]cbs$",
    samples$endpoint_ploidy_file
  )
  file_parts <- regmatches(samples$endpoint_ploidy_file, file_match)
  valid_file_parts <- lengths(file_parts) == 3L
  encoded_origin <- rep(NA_character_, nrow(samples))
  encoded_dose <- rep(NA_real_, nrow(samples))
  encoded_origin[valid_file_parts] <- vapply(
    file_parts[valid_file_parts],
    `[[`,
    character(1L),
    2L
  )
  encoded_dose[valid_file_parts] <- figure7_numeric(vapply(
    file_parts[valid_file_parts],
    `[[`,
    character(1L),
    3L
  ))
  expected_curated_counts <- figure7_curated_endpoint_counts()
  observed_curated_counts <- stats::setNames(
    samples$n_endpoint_ploidy_cells,
    samples$sample_id
  )
  if (nrow(samples) != 16L || anyDuplicated(samples$endpoint_ploidy_file) ||
      !setequal(samples$endpoint_ploidy_file, unique(endpoint_cells$file)) ||
      sum(samples$n_endpoint_ploidy_cells) != 9832L ||
      !identical(
        observed_curated_counts[names(expected_curated_counts)],
        expected_curated_counts
      ) ||
      any(!valid_file_parts) ||
      any(encoded_origin != samples$initial_ploidy) ||
      any(encoded_dose != samples$dose_mg)) {
    figure7_stop(
      "Figure 7 sample-to-CBS mapping must validate the complete 14,125-cell ",
      "inventory while selecting the exact reviewed 9,832-cell score ",
      "universe with matching sample, injected origin, and dose"
    )
  }
  samples$matched_control_reference_delta <- NA_real_
  samples$matched_control_n <- NA_integer_
  samples$matched_control_sample_ids <- NA_character_
  for (ploidy in c("2N", "4N")) {
    controls <- samples[samples$initial_ploidy == ploidy & samples$dose_mg == 0, , drop = FALSE]
    if (!nrow(controls)) figure7_stop("No untreated controls for initial ploidy ", ploidy)
    hit <- samples$initial_ploidy == ploidy
    samples$matched_control_reference_delta[hit] <- mean(controls[[delta_measure]])
    samples$matched_control_n[hit] <- nrow(controls)
    samples$matched_control_sample_ids[hit] <- paste(sort(controls$sample_id), collapse = ";")
  }
  samples[[tgi_measure]] <- 100 * (1 - samples[[delta_measure]] /
                                      samples$matched_control_reference_delta)
  treated <- samples$dose_mg %in% as.numeric(unlist(config$tgi$treated_doses_mg_per_kg))
  difference <- abs(samples[[tgi_measure]] - samples[[embedded_tgi_measure]])
  tolerance <- as.numeric(config$tgi$numerical_tolerance)
  if (any(!is.finite(samples[[embedded_tgi_measure]][treated]))) {
    figure7_stop("All treated mice require a finite embedded Day-", figure7_tgi_day(config),
                 " TGI for verification")
  }
  if (any(difference[treated] > tolerance)) {
    figure7_stop("Recomputed Day-", figure7_tgi_day(config), " TGI disagrees with embedded values for: ",
                 paste(samples$sample_id[treated & difference > tolerance], collapse = ", "))
  }
  samples$etp_group <- ifelse(samples$sample_mean_endpoint_ploidy > as.numeric(config$etp$threshold),
                              "ETP-higher", "ETP-lower")
  samples <- figure7_add_tgi_metadata(samples, config)
  samples <- samples[order(samples$dose_mg, samples$initial_ploidy, samples$sample_id), , drop = FALSE]
  treated_doses <- as.numeric(unlist(config$tgi$treated_doses_mg_per_kg))
  treated_rows <- samples$dose_mg %in% treated_doses
  if (sum(treated_rows) != 8L || any(samples$dose_mg > 0 & !treated_rows)) {
    figure7_stop("Panels 7C-7E require exactly eight treated mice at 30 or 120 mg/kg")
  }
  if (any(!is.finite(samples[[tgi_measure]][treated_rows]))) {
    figure7_stop("All eight treated mice require finite recomputed Day-", figure7_tgi_day(config), " TGI")
  }
  rownames(samples) <- NULL
  samples
}

figure7_prepare_cellcycle <- function(cellcycle, samples, config) {
  tgi_measure <- figure7_tgi_measure(config)
  data <- cellcycle[is.finite(cellcycle$pseudotime), , drop = FALSE]
  index <- match(data$sample_id, samples$sample_id)
  if (anyNA(index)) figure7_stop("Missing sample metadata for CellCycle rows")
  data$dose <- samples$dose[index]; data$dose_mg <- samples$dose_mg[index]
  data$sample_mean_endpoint_ploidy <- samples$sample_mean_endpoint_ploidy[index]
  data$etp_group <- samples$etp_group[index]
  data[[tgi_measure]] <- samples[[tgi_measure]][index]
  data
}
