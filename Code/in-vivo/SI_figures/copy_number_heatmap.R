# Helpers for main Figure 7J: downstream NUMBAT-derived CBS states.
#
# The injected 2N- and 4N-origin matrices have different consensus segment
# boundaries and may reflect different coordinate conventions. They must
# therefore never be locus-aligned directly. These helpers validate every
# cell against all_ploidy.tsv and reduce each schema independently to a common
# set of chromosome-level, available-segment length-weighted means.

si_copy_number_truthy <- function(value) {
  tolower(trimws(as.character(value))) %in% c("true", "t", "1", "yes")
}

si_copy_number_sha256 <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("CBS manifest validation requires the digest package", call. = FALSE)
  }
  digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)
}

si_copy_number_validate_manifest <- function(
  cbs_dir,
  manifest_path = file.path(cbs_dir, "cbs_manifest.tsv")
) {
  if (!file.exists(manifest_path)) {
    stop("Missing reviewed CBS manifest: ", manifest_path, call. = FALSE)
  }
  manifest <- utils::read.delim(
    manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  required <- c("filename", "bytes", "sha256", "notes")
  manifest$bytes <- suppressWarnings(as.numeric(manifest$bytes))
  if (!identical(names(manifest), required) || nrow(manifest) != 16L ||
      anyNA(manifest$filename) || any(!nzchar(manifest$filename)) ||
      any(basename(manifest$filename) != manifest$filename) ||
      anyDuplicated(manifest$filename) || any(!is.finite(manifest$bytes)) ||
      any(manifest$bytes <= 0) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256))) {
    stop(
      "Reviewed CBS manifest must contain exactly 16 unique pinned matrices",
      call. = FALSE
    )
  }
  files <- sort(Sys.glob(file.path(cbs_dir, "*.sps.cbs")))
  observed <- sort(basename(files))
  expected <- sort(manifest$filename)
  if (!identical(observed, expected)) {
    stop(
      "CBS manifest inventory mismatch; missing=",
      paste(setdiff(expected, observed), collapse = ","),
      "; unexpected=",
      paste(setdiff(observed, expected), collapse = ","),
      call. = FALSE
    )
  }
  reviewed_manifest <- manifest
  manifest <- manifest[match(basename(files), manifest$filename), , drop = FALSE]
  observed_bytes <- as.numeric(file.info(files)$size)
  observed_hashes <- vapply(files, si_copy_number_sha256, character(1L))
  mismatch <- observed_bytes != manifest$bytes |
    observed_hashes != manifest$sha256
  if (any(mismatch)) {
    stop(
      "CBS matrix differs from its reviewed checksum: ",
      paste(manifest$filename[mismatch], collapse = ", "),
      call. = FALSE
    )
  }
  reviewed_manifest
}

si_copy_number_validate_reference_manifest <- function(
  reference_dir,
  manifest_path = file.path(reference_dir, "reference_manifest.tsv")
) {
  if (!file.exists(manifest_path)) {
    stop("Missing injected-cell reference manifest: ", manifest_path,
         call. = FALSE)
  }
  manifest <- utils::read.delim(
    manifest_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  required <- c(
    "filename", "injected_origin", "reference_label", "n_cells", "bytes",
    "sha256", "source_repository", "source_commit", "designation_basis",
    "chr999_unit", "chr999_interpretation_basis", "ploidy_policy"
  )
  manifest$n_cells <- suppressWarnings(as.integer(manifest$n_cells))
  manifest$bytes <- suppressWarnings(as.numeric(manifest$bytes))
  if (!identical(names(manifest), required) || nrow(manifest) != 2L ||
      !identical(sort(manifest$injected_origin), c("2N", "4N")) ||
      anyNA(manifest$filename) || any(!nzchar(manifest$filename)) ||
      any(basename(manifest$filename) != manifest$filename) ||
      anyDuplicated(manifest$filename) || anyDuplicated(manifest$injected_origin) ||
      anyNA(manifest$reference_label) || any(!nzchar(manifest$reference_label)) ||
      anyNA(manifest$n_cells) || any(manifest$n_cells <= 0L) ||
      any(!is.finite(manifest$bytes)) || any(manifest$bytes <= 0) ||
      any(!grepl("^[0-9a-f]{64}$", manifest$sha256)) ||
      any(manifest$source_repository != "miningcloneid") ||
      any(manifest$source_commit !=
        "c505cd9159fa2a8c0974c7379f6aacd09fe19abc") ||
      anyNA(manifest$designation_basis) ||
      any(!nzchar(manifest$designation_basis)) ||
      any(manifest$chr999_unit !=
        "haploid-genome-equivalent unassigned DNA") ||
      any(manifest$chr999_interpretation_basis !=
        "project-confirmed 2026-08-01") ||
      any(manifest$ploidy_policy != paste(
        "autosomal length-weighted estimate plus chr999",
        "haploid-genome-equivalent unassigned DNA"
      ))) {
    stop(
      "Injected-cell reference manifest must pin one valid 2N and one valid 4N matrix",
      call. = FALSE
    )
  }
  files <- sort(Sys.glob(file.path(reference_dir, "*.sps.cbs")))
  if (!identical(sort(basename(files)), sort(manifest$filename))) {
    stop("Injected-cell reference matrix inventory differs from its manifest",
         call. = FALSE)
  }
  matched <- manifest[match(basename(files), manifest$filename), , drop = FALSE]
  observed_bytes <- as.numeric(file.info(files)$size)
  observed_hashes <- vapply(files, si_copy_number_sha256, character(1L))
  mismatch <- observed_bytes != matched$bytes |
    observed_hashes != matched$sha256
  if (any(mismatch)) {
    stop(
      "Injected-cell reference matrix differs from its checksum: ",
      paste(matched$filename[mismatch], collapse = ", "),
      call. = FALSE
    )
  }
  manifest[match(c("2N", "4N"), manifest$injected_origin), , drop = FALSE]
}

si_copy_number_read_injected_references <- function(reference_dir) {
  manifest <- si_copy_number_validate_reference_manifest(reference_dir)
  rows <- lapply(seq_len(nrow(manifest)), function(index) {
    entry <- manifest[index, , drop = FALSE]
    path <- file.path(reference_dir, entry$filename)
    raw <- utils::read.delim(
      path,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      na.strings = c("NA", "NaN", "")
    )
    required <- c("chr", "seglength")
    if (!all(required %in% names(raw))) {
      stop("Injected-cell reference is not a long-format CBS matrix: ",
           entry$filename, call. = FALSE)
    }
    cell_columns <- grep("^SP_", names(raw), value = TRUE)
    if (length(cell_columns) != entry$n_cells) {
      stop("Injected-cell reference cell count differs from its manifest: ",
           entry$filename, call. = FALSE)
    }
    chromosome <- suppressWarnings(as.integer(raw$chr))
    segment_length <- suppressWarnings(as.numeric(raw$seglength))
    extra_index <- which(chromosome == 999L)
    if (length(extra_index) != 1L) {
      stop(
        "Injected-cell reference must contain exactly one chr999 extra-DNA row: ",
        entry$filename,
        call. = FALSE
      )
    }
    extra_dna_haploid_genome_equivalents <- suppressWarnings(as.numeric(
      raw[extra_index, cell_columns, drop = TRUE]
    ))
    if (length(extra_dna_haploid_genome_equivalents) !=
          length(cell_columns) ||
        any(!is.finite(extra_dna_haploid_genome_equivalents)) ||
        any(extra_dna_haploid_genome_equivalents < 0)) {
      stop("Injected-cell reference has invalid chr999 extra-DNA values: ",
           entry$filename, call. = FALSE)
    }
    keep <- chromosome %in% seq_len(22L)
    if (!any(keep) || any(!is.finite(segment_length[keep])) ||
        any(segment_length[keep] <= 0)) {
      stop("Injected-cell reference has invalid autosomal lengths: ",
           entry$filename, call. = FALSE)
    }
    chromosome <- chromosome[keep]
    segment_length <- segment_length[keep]
    values <- as.matrix(raw[keep, cell_columns, drop = FALSE])
    suppressWarnings(storage.mode(values) <- "double")
    if (any(!is.finite(values) & !is.na(values)) ||
        any(values < 0, na.rm = TRUE)) {
      stop("Injected-cell reference has invalid copy-number values: ",
           entry$filename, call. = FALSE)
    }
    finite <- is.finite(values)
    denominator <- colSums(finite * segment_length)
    assigned_autosomal_ploidy <-
      colSums(values * segment_length, na.rm = TRUE) / denominator
    # chr999 is expressed in haploid-genome-equivalent units of unassigned
    # DNA. It is therefore additive on the ploidy scale, not a multiplicative
    # fraction of the chromosome-assigned karyotype.
    ploidy <- assigned_autosomal_ploidy +
      extra_dna_haploid_genome_equivalents
    chromosome_means <- lapply(sort(unique(chromosome)), function(value) {
      columns <- values[chromosome == value, , drop = FALSE]
      colMeans(columns, na.rm = TRUE)
    })
    assigned_autosomal_chromosomes <- Reduce(`+`, chromosome_means)
    total_chromosomes <- assigned_autosomal_chromosomes +
      22 * extra_dna_haploid_genome_equivalents
    data.frame(
      reference_file = entry$filename,
      reference_cell_id = cell_columns,
      injected_origin = entry$injected_origin,
      reference_label = entry$reference_label,
      designation_basis = entry$designation_basis,
      source_repository = entry$source_repository,
      source_commit = entry$source_commit,
      chr999_unit = entry$chr999_unit,
      chr999_interpretation_basis = entry$chr999_interpretation_basis,
      assigned_autosomal_ploidy = assigned_autosomal_ploidy,
      extra_dna_haploid_genome_equivalents =
        extra_dna_haploid_genome_equivalents,
      ploidy = ploidy,
      assigned_autosomal_chromosomes = assigned_autosomal_chromosomes,
      total_chromosomes = total_chromosomes,
      frac_covered = denominator / sum(segment_length),
      ploidy_policy = entry$ploidy_policy,
      stringsAsFactors = FALSE
    )
  })
  cells <- do.call(rbind, rows)
  rownames(cells) <- NULL
  if (any(!is.finite(cells$ploidy)) ||
      any(!is.finite(cells$total_chromosomes)) ||
      any(!is.finite(cells$frac_covered)) ||
      any(cells$frac_covered <= 0 | cells$frac_covered > 1)) {
    stop("Injected-cell reference calculation produced invalid values",
         call. = FALSE)
  }
  list(manifest = manifest, cells = cells)
}

si_copy_number_parse_segments <- function(segment_names) {
  pattern <- "^([^:]+):([^-]+)-(.+)$"
  matches <- regexec(pattern, as.character(segment_names))
  pieces <- regmatches(as.character(segment_names), matches)
  valid <- lengths(pieces) == 4L
  if (!all(valid)) {
    stop(
      "Malformed CBS segment names: ",
      paste(segment_names[!valid], collapse = ", "),
      call. = FALSE
    )
  }
  chromosome <- suppressWarnings(as.integer(vapply(
    pieces, `[[`, character(1L), 2L
  )))
  start <- suppressWarnings(as.numeric(vapply(
    pieces, `[[`, character(1L), 3L
  )))
  end <- suppressWarnings(as.numeric(vapply(
    pieces, `[[`, character(1L), 4L
  )))
  if (anyNA(chromosome) || any(!is.finite(start)) || any(!is.finite(end)) ||
      any(start < 0) || any(end <= start)) {
    stop("CBS segment coordinates are invalid", call. = FALSE)
  }
  data.frame(
    segment = as.character(segment_names),
    chromosome = chromosome,
    start = start,
    end = end,
    length = end - start,
    stringsAsFactors = FALSE
  )
}

si_copy_number_read_ploidy <- function(path) {
  ploidy <- utils::read.delim(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "",
    comment.char = ""
  )
  required <- c("file", "cell_id", "ploidy", "frac_covered", "format")
  if (!identical(names(ploidy), required)) {
    stop(
      "all_ploidy.tsv must have columns: ",
      paste(required, collapse = ", "),
      call. = FALSE
    )
  }
  ploidy$file <- as.character(ploidy$file)
  ploidy$cell_id <- as.character(ploidy$cell_id)
  ploidy$ploidy <- suppressWarnings(as.numeric(ploidy$ploidy))
  ploidy$frac_covered <- suppressWarnings(as.numeric(ploidy$frac_covered))
  key <- paste(ploidy$file, ploidy$cell_id, sep = "::")
  if (!nrow(ploidy) || anyNA(ploidy$file) || any(!nzchar(ploidy$file)) ||
      anyNA(ploidy$cell_id) || any(!nzchar(ploidy$cell_id)) ||
      anyDuplicated(key) || any(!is.finite(ploidy$ploidy)) ||
      any(!is.finite(ploidy$frac_covered)) ||
      any(ploidy$frac_covered <= 0 | ploidy$frac_covered > 1) ||
      any(ploidy$format != "wide")) {
    stop("all_ploidy.tsv failed its CBS key/value contract", call. = FALSE)
  }
  ploidy
}

si_copy_number_qc_selection <- function(endpoint_audit) {
  required <- c(
    "cell", "context", "initial_ploidy", "endpoint_file",
    "endpoint_cell_id", "endpoint_ploidy", "matched"
  )
  if (!all(required %in% names(endpoint_audit))) {
    stop(
      "Endpoint-ploidy audit lacks the QC-passed cell-universe columns",
      call. = FALSE
    )
  }
  selected <- endpoint_audit[
    si_copy_number_truthy(endpoint_audit$matched),
    required,
    drop = FALSE
  ]
  selected$cell <- as.character(selected$cell)
  selected$context <- as.character(selected$context)
  selected$initial_ploidy <- as.character(selected$initial_ploidy)
  selected$endpoint_file <- as.character(selected$endpoint_file)
  selected$endpoint_cell_id <- as.character(selected$endpoint_cell_id)
  selected$endpoint_ploidy <- suppressWarnings(as.numeric(
    selected$endpoint_ploidy
  ))
  sample_id_length <- nchar(selected$cell) -
    nchar(selected$endpoint_cell_id) - 1L
  selected$sample_id <- ifelse(
    sample_id_length > 0L,
    substring(selected$cell, 1L, sample_id_length),
    ""
  )
  key <- paste(selected$endpoint_file, selected$endpoint_cell_id, sep = "::")
  if (nrow(selected) != 9832L ||
      anyNA(selected[, setdiff(required, "matched"), drop = FALSE]) ||
      any(!nzchar(selected$cell)) || any(!nzchar(selected$endpoint_file)) ||
      any(!nzchar(selected$endpoint_cell_id)) || anyDuplicated(selected$cell) ||
      anyDuplicated(key) || any(selected$context != "Tumor") ||
      any(!nzchar(selected$sample_id)) ||
      any(selected$cell != paste(
        selected$sample_id,
        selected$endpoint_cell_id,
        sep = "_"
      )) ||
      !identical(
        as.integer(table(factor(
          selected$initial_ploidy,
          levels = c("2N", "4N")
        ))),
        c(4880L, 4952L)
      ) ||
      length(unique(selected$endpoint_file)) != 16L ||
      any(!is.finite(selected$endpoint_ploidy))) {
    stop(
      "Endpoint-ploidy audit must identify the exact 9,832-cell ",
      "QC-passed tumor universe",
      call. = FALSE
    )
  }
  selected
}

si_copy_number_sample_lookup <- function(endpoint_audit) {
  required <- c(
    "cell", "initial_ploidy", "endpoint_file", "endpoint_cell_id", "matched"
  )
  if (!all(required %in% names(endpoint_audit))) {
    stop("Endpoint-ploidy audit lacks the sample lookup columns", call. = FALSE)
  }
  matched <- endpoint_audit[
    si_copy_number_truthy(endpoint_audit$matched) &
      !is.na(endpoint_audit$endpoint_file) &
      nzchar(as.character(endpoint_audit$endpoint_file)),
    required,
    drop = FALSE
  ]
  matched$sample_id <- sub("_[^_]+$", "", as.character(matched$cell))
  mapping <- unique(matched[, c(
    "endpoint_file", "sample_id", "initial_ploidy"
  ), drop = FALSE])
  count <- table(mapping$endpoint_file)
  if (!nrow(mapping) || any(count != 1L)) {
    stop("Each CBS file must map to exactly one sequenced mouse", call. = FALSE)
  }
  mapping
}

si_copy_number_filename_metadata <- function(filename) {
  match <- regexec(
    "^SUM159-(2N|4N)-([0-9]+)-(.+)_harvest[.]sps[.]cbs$",
    filename
  )
  pieces <- regmatches(filename, match)[[1L]]
  if (length(pieces) != 4L) {
    stop("Unexpected CBS filename: ", filename, call. = FALSE)
  }
  list(
    initial_ploidy = pieces[[2L]],
    dose_mg_per_kg = as.numeric(pieces[[3L]])
  )
}

si_copy_number_read_collection <- function(
  cbs_dir,
  all_ploidy_path,
  endpoint_audit,
  sample_metadata = NULL,
  tolerance = 1e-10
) {
  reviewed_manifest <- si_copy_number_validate_manifest(cbs_dir)
  # Keep the reviewed row order as part of the derivation contract.  This is
  # the same order used by weighted_ploidy.py --manifest and therefore makes
  # regenerated all_ploidy outputs byte-deterministic rather than glob-order
  # dependent.
  files <- file.path(cbs_dir, reviewed_manifest$filename)
  if (!length(files)) {
    stop("No .sps.cbs matrices found under: ", cbs_dir, call. = FALSE)
  }
  ploidy <- si_copy_number_read_ploidy(all_ploidy_path)
  expected_files <- sort(unique(ploidy$file))
  observed_files <- sort(basename(files))
  if (!identical(observed_files, expected_files)) {
    stop(
      "CBS matrix inventory does not match all_ploidy.tsv; missing=",
      paste(setdiff(expected_files, observed_files), collapse = ","),
      "; unexpected=",
      paste(setdiff(observed_files, expected_files), collapse = ","),
      call. = FALSE
    )
  }
  qc_selection <- si_copy_number_qc_selection(endpoint_audit)
  sample_lookup <- si_copy_number_sample_lookup(endpoint_audit)

  matrices <- vector("list", length(files))
  annotations <- vector("list", length(files))
  names(matrices) <- observed_files
  names(annotations) <- observed_files
  for (index in seq_along(files)) {
    path <- files[[index]]
    filename <- basename(path)
    raw <- utils::read.delim(
      path,
      row.names = 1L,
      check.names = FALSE,
      stringsAsFactors = FALSE,
      quote = "",
      comment.char = "",
      na.strings = c("NA", "NaN", "")
    )
    if (!nrow(raw) || !ncol(raw) || anyDuplicated(rownames(raw))) {
      stop("CBS matrix is empty or has duplicate cell IDs: ", filename,
           call. = FALSE)
    }
    values <- as.matrix(raw)
    suppressWarnings(storage.mode(values) <- "double")
    if (any(!is.finite(values) & !is.na(values)) ||
        any(values < 0, na.rm = TRUE)) {
      stop("CBS matrix has invalid copy-number values: ", filename,
           call. = FALSE)
    }
    segments <- si_copy_number_parse_segments(colnames(values))
    autosomal <- segments$chromosome %in% seq_len(22L)
    if (!any(autosomal)) {
      stop("CBS matrix has no autosomal segments: ", filename,
           call. = FALSE)
    }
    autosomal_values <- values[, autosomal, drop = FALSE]
    lengths <- segments$length[autosomal]
    finite <- is.finite(autosomal_values)
    denominator <- as.numeric(finite %*% lengths)
    numerator <- as.numeric(replace(autosomal_values, !finite, 0) %*% lengths)
    derived_ploidy <- numerator / denominator
    derived_coverage <- denominator / sum(lengths)

    expected <- ploidy[ploidy$file == filename, , drop = FALSE]
    position <- match(rownames(values), expected$cell_id)
    if (nrow(expected) != nrow(values) || anyNA(position) ||
        anyDuplicated(position)) {
      stop("CBS cell keys do not match all_ploidy.tsv: ", filename,
           call. = FALSE)
    }
    expected <- expected[position, , drop = FALSE]
    if (max(abs(derived_ploidy - expected$ploidy)) > tolerance ||
        max(abs(derived_coverage - expected$frac_covered)) > tolerance) {
      stop("CBS-derived ploidy disagrees with all_ploidy.tsv: ", filename,
           call. = FALSE)
    }

    # Validate the complete checksum-pinned CBS source above, but analyze only
    # cells retained in Tao's final QC-curated Seurat object.  The endpoint
    # audit is the exact bridge between that 9,832-cell transcriptional
    # universe and the 14,125-cell downstream CBS source.
    selected <- qc_selection[
      qc_selection$endpoint_file == filename, , drop = FALSE
    ]
    selected_position <- match(expected$cell_id, selected$endpoint_cell_id)
    keep <- !is.na(selected_position)
    if (!nrow(selected) || sum(keep) != nrow(selected) ||
        anyDuplicated(selected_position[keep]) ||
        any(abs(
          expected$ploidy[keep] -
            selected$endpoint_ploidy[selected_position[keep]]
        ) > tolerance)) {
      stop(
        "QC-passed endpoint audit does not match canonical CBS cells: ",
        filename,
        call. = FALSE
      )
    }
    values <- values[keep, , drop = FALSE]
    expected <- expected[keep, , drop = FALSE]

    file_metadata <- si_copy_number_filename_metadata(filename)
    sample_row <- sample_lookup[
      sample_lookup$endpoint_file == filename, , drop = FALSE
    ]
    if (nrow(sample_row) != 1L ||
        any(selected$sample_id != as.character(sample_row$sample_id)) ||
        any(selected$cell != paste(
          as.character(sample_row$sample_id),
          selected$endpoint_cell_id,
          sep = "_"
        )) ||
        !identical(
          as.character(sample_row$initial_ploidy),
          file_metadata$initial_ploidy
        )) {
      stop("CBS filename injected origin disagrees with endpoint audit: ", filename,
           call. = FALSE)
    }
    if (!is.null(sample_metadata)) {
      required_sample_columns <- c("mouse", "initial_ploidy", "dose")
      if (!all(required_sample_columns %in% names(sample_metadata))) {
        stop("Sample metadata lacks mouse, initial_ploidy, or dose",
             call. = FALSE)
      }
      reviewed_sample <- unique(sample_metadata[
        as.character(sample_metadata$mouse) == as.character(sample_row$sample_id),
        required_sample_columns,
        drop = FALSE
      ])
      reviewed_dose <- suppressWarnings(as.numeric(sub(
        "mg/kg$", "", as.character(reviewed_sample$dose)
      )))
      if (nrow(reviewed_sample) != 1L ||
          !identical(
            as.character(reviewed_sample$initial_ploidy),
            file_metadata$initial_ploidy
          ) ||
          length(reviewed_dose) != 1L || !is.finite(reviewed_dose) ||
          !identical(reviewed_dose, file_metadata$dose_mg_per_kg)) {
        stop(
          "CBS filename ploidy/dose disagrees with reviewed sample metadata: ",
          filename,
          call. = FALSE
        )
      }
    }
    heatmap_row_id <- paste(filename, rownames(values), sep = "::")
    rownames(values) <- heatmap_row_id
    matrices[[filename]] <- list(
      values = values,
      segments = segments
    )
    annotations[[filename]] <- data.frame(
      heatmap_row_id = heatmap_row_id,
      file = filename,
      cell_id = expected$cell_id,
      sample_id = as.character(sample_row$sample_id),
      initial_ploidy = file_metadata$initial_ploidy,
      dose_mg_per_kg = file_metadata$dose_mg_per_kg,
      endpoint_ploidy = expected$ploidy,
      frac_covered = expected$frac_covered,
      stringsAsFactors = FALSE
    )
  }
  cell_annotations <- do.call(rbind, annotations)
  rownames(cell_annotations) <- NULL
  treated <- cell_annotations$dose_mg_per_kg > 0
  if (nrow(cell_annotations) != 9832L || sum(treated) != 5335L ||
      sum(!treated) != 4497L ||
      anyDuplicated(cell_annotations$heatmap_row_id)) {
    stop(
      "Combined CBS analysis collection must contain the exact QC-passed ",
      "9,832-cell universe, including 5,335 treated cells",
      call. = FALSE
    )
  }
  list(
    matrices = matrices,
    cell_annotations = cell_annotations,
    ploidy = ploidy,
    qc_selection = qc_selection,
    reviewed_manifest = reviewed_manifest
  )
}

si_copy_number_chromosome_matrix <- function(values, segments) {
  result <- matrix(
    NA_real_,
    nrow = nrow(values),
    ncol = 22L,
    dimnames = list(rownames(values), paste0("chr", seq_len(22L)))
  )
  for (chromosome in seq_len(22L)) {
    selected <- which(segments$chromosome == chromosome)
    if (!length(selected)) next
    chromosome_values <- values[, selected, drop = FALSE]
    weights <- segments$length[selected]
    finite <- is.finite(chromosome_values)
    denominator <- as.numeric(finite %*% weights)
    numerator <- as.numeric(
      replace(chromosome_values, !finite, 0) %*% weights
    )
    result[, chromosome] <- numerator / denominator
    result[denominator <= 0, chromosome] <- NA_real_
  }
  result
}

si_copy_number_schema_audit <- function(collection) {
  annotation <- collection$cell_annotations
  rows <- lapply(names(collection$matrices), function(filename) {
    item <- collection$matrices[[filename]]
    file_annotation <- annotation[annotation$file == filename, , drop = FALSE]
    origin <- unique(file_annotation$initial_ploidy)
    if (length(origin) != 1L) {
      stop("CBS file has an ambiguous injected origin: ", filename,
           call. = FALSE)
    }
    do.call(rbind, lapply(seq_len(22L), function(chromosome) {
      exported <- item$segments$chromosome == chromosome
      finite_by_segment <- colSums(is.finite(item$values)) > 0L
      available <- exported & finite_by_segment
      exported_lengths <- item$segments$length[exported]
      cell_available_fraction <- as.numeric(
        is.finite(item$values[, exported, drop = FALSE]) %*%
          exported_lengths
      ) / sum(exported_lengths)
      data.frame(
        file = filename,
        initial_ploidy = origin,
        chromosome = chromosome,
        n_exported_segments = sum(exported),
        exported_bp = sum(exported_lengths),
        n_available_segments = sum(available),
        n_all_missing_exported_segments = sum(exported & !finite_by_segment),
        represented_bp = sum(item$segments$length[available]),
        represented_bp_fraction_of_exported =
          sum(item$segments$length[available]) / sum(exported_lengths),
        minimum_cell_available_bp_fraction = min(cell_available_fraction),
        mean_cell_available_bp_fraction = mean(cell_available_fraction),
        maximum_cell_available_bp_fraction = max(cell_available_fraction),
        minimum_available_coordinate = if (any(available)) {
          min(item$segments$start[available])
        } else {
          NA_real_
        },
        maximum_available_coordinate = if (any(available)) {
          max(item$segments$end[available])
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
  })
  file_audit <- do.call(rbind, rows)
  schema_signature <- vapply(collection$matrices, function(item) {
    paste(item$segments$segment, collapse = ";")
  }, character(1L))
  signature_id <- match(schema_signature, unique(schema_signature))
  schema_lookup <- data.frame(
    file = names(schema_signature),
    schema_id = paste0("schema_", signature_id),
    stringsAsFactors = FALSE
  )
  file_audit$schema_id <- schema_lookup$schema_id[
    match(file_audit$file, schema_lookup$file)
  ]
  schema_audit <- unique(file_audit[, c(
    "schema_id", "initial_ploidy", "chromosome", "n_exported_segments",
    "exported_bp", "n_available_segments",
    "n_all_missing_exported_segments", "represented_bp",
    "represented_bp_fraction_of_exported",
    "minimum_available_coordinate",
    "maximum_available_coordinate"
  ), drop = FALSE])
  schema_counts <- table(schema_audit$schema_id)
  if (length(unique(schema_signature)) != 2L || any(schema_counts != 22L) ||
      any(schema_audit$n_exported_segments < 1L) ||
      any(!is.finite(schema_audit$exported_bp)) ||
      any(schema_audit$exported_bp <= 0) ||
      any(schema_audit$n_available_segments < 0L) ||
      any(schema_audit$n_all_missing_exported_segments < 0L) ||
      any(!is.finite(schema_audit$represented_bp)) ||
      any(schema_audit$represented_bp < 0) ||
      any(!is.finite(schema_audit$represented_bp_fraction_of_exported)) ||
      any(schema_audit$represented_bp_fraction_of_exported < 0 |
        schema_audit$represented_bp_fraction_of_exported > 1) ||
      any(!is.finite(file_audit$minimum_cell_available_bp_fraction)) ||
      any(!is.finite(file_audit$mean_cell_available_bp_fraction)) ||
      any(!is.finite(file_audit$maximum_cell_available_bp_fraction)) ||
      any(file_audit$minimum_cell_available_bp_fraction < 0 |
        file_audit$maximum_cell_available_bp_fraction > 1)) {
    stop("Expected exactly two complete 22-autosome CBS schemas",
         call. = FALSE)
  }
  schema_audit <- schema_audit[order(
    match(schema_audit$initial_ploidy, c("2N", "4N")),
    schema_audit$chromosome
  ), , drop = FALSE]
  rownames(schema_audit) <- NULL
  list(file = file_audit, schema = schema_audit)
}

si_copy_number_harmonize <- function(collection) {
  chromosome_matrices <- lapply(collection$matrices, function(item) {
    si_copy_number_chromosome_matrix(item$values, item$segments)
  })
  matrix_data <- do.call(rbind, chromosome_matrices)
  annotation <- collection$cell_annotations
  position <- match(rownames(matrix_data), annotation$heatmap_row_id)
  if (anyNA(position) || anyDuplicated(position)) {
    stop("Projected CBS matrix lost its cell annotations", call. = FALSE)
  }
  annotation <- annotation[position, , drop = FALSE]
  lineage_order <- match(annotation$initial_ploidy, c("2N", "4N"))
  sample_table <- unique(annotation[, c(
    "sample_id", "initial_ploidy", "dose_mg_per_kg"
  ), drop = FALSE])
  sample_table <- sample_table[order(
    match(sample_table$initial_ploidy, c("2N", "4N")),
    sample_table$dose_mg_per_kg,
    sample_table$sample_id
  ), , drop = FALSE]
  sample_levels <- sample_table$sample_id
  ordering <- order(
    lineage_order,
    match(annotation$sample_id, sample_levels),
    annotation$endpoint_ploidy,
    annotation$cell_id
  )
  matrix_data <- matrix_data[ordering, , drop = FALSE]
  annotation <- annotation[ordering, , drop = FALSE]
  annotation$display_order <- seq_len(nrow(annotation))
  chromosome_audit <- si_copy_number_schema_audit(collection)
  list(
    matrix = matrix_data,
    cell_annotations = annotation,
    chromosomes = data.frame(
      chromosome = seq_len(22L),
      column = paste0("chr", seq_len(22L)),
      display_order = seq_len(22L),
      fraction_all_cells_available = colMeans(is.finite(matrix_data)),
      stringsAsFactors = FALSE
    ),
    chromosome_schema_audit = chromosome_audit$schema,
    chromosome_file_audit = chromosome_audit$file,
    sample_levels = sample_levels
  )
}

si_copy_number_heatmap <- function(
  harmonized,
  title = paste(
    "NUMBAT-derived copy-number states",
    "(cell-by-chromosome available-segment mean CN)"
  ),
  labels_col = paste0("chr", seq_len(22L)),
  fontsize = 8,
  fontsize_col = 6.5,
  annotation_legend = TRUE,
  angle_col = 0
) {
  annotation <- harmonized$cell_annotations
  matrix_data <- harmonized$matrix
  sample_levels <- harmonized$sample_levels
  visible_labels <- as.character(labels_col)[nzchar(as.character(labels_col))]
  if (length(labels_col) != ncol(matrix_data) || anyNA(labels_col) ||
      !length(visible_labels) || anyDuplicated(visible_labels)) {
    stop("Copy-number heatmap chromosome labels are invalid", call. = FALSE)
  }
  angle_col <- as.character(angle_col)
  if (length(angle_col) != 1L || !angle_col %in% c("0", "45", "90", "270", "315")) {
    stop("Copy-number heatmap column-label angle is invalid", call. = FALSE)
  }
  dose_labels <- c(
    "0" = "Vehicle",
    "30" = "30 mg/kg",
    "120" = "120 mg/kg"
  )
  dose_key <- as.character(annotation$dose_mg_per_kg)
  if (any(!dose_key %in% names(dose_labels))) {
    stop("Copy-number heatmap contains an unexpected gemcitabine dose",
         call. = FALSE)
  }
  row_annotation <- data.frame(
    `Injected origin` = factor(
      annotation$initial_ploidy, levels = c("2N", "4N")
    ),
    `Gemcitabine dose` = factor(
      unname(dose_labels[dose_key]),
      levels = unname(dose_labels)
    ),
    Mouse = factor(annotation$sample_id, levels = sample_levels),
    check.names = FALSE
  )
  rownames(row_annotation) <- annotation$heatmap_row_id
  sample_colors <- stats::setNames(
    grDevices::hcl.colors(length(sample_levels), "Dynamic"),
    sample_levels
  )
  annotation_colors <- list(
    `Injected origin` = c("2N" = "#4C78A8", "4N" = "#E45756"),
    `Gemcitabine dose` = c(
      "Vehicle" = "#666666",
      "30 mg/kg" = "#D95F02",
      "120 mg/kg" = "#1B9E77"
    ),
    Mouse = sample_colors
  )
  sample_counts <- table(factor(annotation$sample_id, levels = sample_levels))
  gaps_row <- head(cumsum(as.integer(sample_counts)), -1L)
  gaps_col <- seq_len(21L)
  labels_col <- as.character(labels_col)
  # Copy number 3 is frequent in these tumors and must remain visible on the
  # white manuscript background.  The gold anchor distinguishes this common
  # one-copy gain from both the blue low-copy states and red high-copy states.
  colors <- grDevices::colorRampPalette(c(
    "#2166AC", "#67A9CF", "#E6C84F", "#F4A582", "#B2182B", "#762A83"
  ))(120L)
  breaks <- seq(0.5, 6.5, length.out = length(colors) + 1L)
  na_color <- "#D9D9D9"
  heatmap <- pheatmap::pheatmap(
    matrix_data,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    gaps_row = gaps_row,
    gaps_col = gaps_col,
    color = colors,
    breaks = breaks,
    legend_breaks = 1:6,
    legend_labels = as.character(1:6),
    na_col = na_color,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = TRUE,
    labels_col = labels_col,
    fontsize = fontsize,
    fontsize_col = fontsize_col,
    angle_col = angle_col,
    annotation_row = row_annotation,
    annotation_colors = annotation_colors,
    annotation_legend = annotation_legend,
    annotation_names_row = FALSE,
    main = if (is.null(title)) NA_character_ else as.character(title),
    silent = TRUE
  )
  list(
    gtable = heatmap$gtable,
    row_annotation = row_annotation,
    annotation_colors = annotation_colors,
    colors = colors,
    breaks = breaks,
    na_color = na_color,
    gaps_row = gaps_row,
    gaps_col = gaps_col,
    labels_col = labels_col,
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
}

si_copy_number_sample_summary <- function(cell_annotations) {
  split_rows <- split(cell_annotations, cell_annotations$sample_id)
  summary <- do.call(rbind, lapply(split_rows, function(data) {
    data.frame(
      sample_id = data$sample_id[[1L]],
      file = data$file[[1L]],
      initial_ploidy = data$initial_ploidy[[1L]],
      dose_mg_per_kg = data$dose_mg_per_kg[[1L]],
      n_cells = nrow(data),
      mean_postprocessed_copy_number_score = mean(data$endpoint_ploidy),
      median_postprocessed_copy_number_score = stats::median(data$endpoint_ploidy),
      postprocessed_copy_number_score_q25 = unname(stats::quantile(
        data$endpoint_ploidy, 0.25, names = FALSE
      )),
      postprocessed_copy_number_score_q75 = unname(stats::quantile(
        data$endpoint_ploidy, 0.75, names = FALSE
      )),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary) <- NULL
  summary[order(
    match(summary$initial_ploidy, c("2N", "4N")),
    summary$dose_mg_per_kg,
    summary$sample_id
  ), , drop = FALSE]
}

si_copy_number_reduction_summary <- function(
  reference_cells,
  sample_summary,
  terminal_cells
) {
  required_reference <- c(
    "injected_origin", "reference_label", "designation_basis",
    "source_repository", "source_commit", "chr999_unit",
    "chr999_interpretation_basis", "assigned_autosomal_ploidy",
    "extra_dna_haploid_genome_equivalents", "ploidy", "ploidy_policy"
  )
  required_sample <- c(
    "sample_id", "initial_ploidy", "mean_postprocessed_copy_number_score"
  )
  required_terminal <- c("initial_ploidy", "endpoint_ploidy")
  if (!all(required_reference %in% names(reference_cells)) ||
      !all(required_sample %in% names(sample_summary)) ||
      !all(required_terminal %in% names(terminal_cells))) {
    stop("Ploidy-reduction summary inputs do not satisfy their schemas",
         call. = FALSE)
  }
  rows <- lapply(c("2N", "4N"), function(origin) {
    reference <- reference_cells[
      reference_cells$injected_origin == origin, , drop = FALSE
    ]
    samples <- sample_summary[
      sample_summary$initial_ploidy == origin, , drop = FALSE
    ]
    cells <- terminal_cells[
      terminal_cells$initial_ploidy == origin, , drop = FALSE
    ]
    if (!nrow(reference) || !nrow(samples) || !nrow(cells)) {
      stop("Missing reference or endpoint observations for origin ", origin,
           call. = FALSE)
    }
    reference_mean <- mean(reference$ploidy)
    endpoint_mean <- mean(samples$mean_postprocessed_copy_number_score)
    reference_policy <- unique(reference$ploidy_policy)
    reference_label <- unique(reference$reference_label)
    designation_basis <- unique(reference$designation_basis)
    source_repository <- unique(reference$source_repository)
    source_commit <- unique(reference$source_commit)
    chr999_unit <- unique(reference$chr999_unit)
    chr999_interpretation_basis <- unique(
      reference$chr999_interpretation_basis
    )
    if (length(reference_policy) != 1L || length(reference_label) != 1L ||
        length(designation_basis) != 1L ||
        length(source_repository) != 1L || length(source_commit) != 1L ||
        length(chr999_unit) != 1L ||
        length(chr999_interpretation_basis) != 1L) {
      stop("Injected-cell reference provenance is ambiguous for origin ",
           origin, call. = FALSE)
    }
    data.frame(
      injected_origin = origin,
      reference_label = reference_label,
      designation_basis = designation_basis,
      reference_source_repository = source_repository,
      reference_source_commit = source_commit,
      reference_chr999_unit = chr999_unit,
      reference_chr999_interpretation_basis = chr999_interpretation_basis,
      reference_ploidy_policy = reference_policy,
      n_reference_cells = nrow(reference),
      reference_mean_assigned_autosomal_ploidy = mean(
        reference$assigned_autosomal_ploidy
      ),
      reference_mean_chr999_extra_dna_haploid_genome_equivalents = mean(
        reference$extra_dna_haploid_genome_equivalents
      ),
      reference_mean_ploidy = reference_mean,
      reference_min_ploidy = min(reference$ploidy),
      reference_max_ploidy = max(reference$ploidy),
      n_endpoint_mice = nrow(samples),
      n_endpoint_cells = nrow(cells),
      endpoint_mouse_balanced_mean_ploidy = endpoint_mean,
      endpoint_min_mouse_mean_ploidy = min(
        samples$mean_postprocessed_copy_number_score
      ),
      endpoint_max_mouse_mean_ploidy = max(
        samples$mean_postprocessed_copy_number_score
      ),
      endpoint_min_cell_ploidy = min(cells$endpoint_ploidy),
      endpoint_max_cell_ploidy = max(cells$endpoint_ploidy),
      absolute_change = endpoint_mean - reference_mean,
      relative_change_percent = 100 * (endpoint_mean / reference_mean - 1),
      all_endpoint_mouse_means_below_reference_min = all(
        samples$mean_postprocessed_copy_number_score < min(reference$ploidy)
      ),
      all_endpoint_cells_below_reference_min = all(
        cells$endpoint_ploidy < min(reference$ploidy)
      ),
      analysis_type = "descriptive_only",
      formal_test_performed = FALSE,
      inference = paste(
        "Descriptive comparison of project-designated injected-cell",
        "karyotype-reference cells with independent endpoint mouse means;",
        "no P value because the reference has one culture-level biological unit."
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

si_copy_number_separation_summary <- function(change_summary) {
  required <- c(
    "injected_origin", "reference_mean_ploidy",
    "endpoint_mouse_balanced_mean_ploidy", "analysis_type",
    "formal_test_performed"
  )
  if (!all(required %in% names(change_summary)) ||
      !identical(as.character(change_summary$injected_origin), c("2N", "4N")) ||
      any(change_summary$analysis_type != "descriptive_only") ||
      any(change_summary$formal_test_performed)) {
    stop("Ploidy-separation summary requires ordered descriptive 2N/4N rows",
         call. = FALSE)
  }
  reference_separation <-
    change_summary$reference_mean_ploidy[[2L]] -
    change_summary$reference_mean_ploidy[[1L]]
  endpoint_separation <-
    change_summary$endpoint_mouse_balanced_mean_ploidy[[2L]] -
    change_summary$endpoint_mouse_balanced_mean_ploidy[[1L]]
  if (!is.finite(reference_separation) || reference_separation <= 0 ||
      !is.finite(endpoint_separation)) {
    stop("Ploidy separation is invalid", call. = FALSE)
  }
  data.frame(
    reference_4n_minus_2n_mean_ploidy = reference_separation,
    endpoint_4n_minus_2n_mouse_balanced_mean_ploidy = endpoint_separation,
    absolute_separation_change = endpoint_separation - reference_separation,
    separation_contraction_percent =
      100 * (1 - endpoint_separation / reference_separation),
    analysis_type = "descriptive_only",
    formal_test_performed = FALSE,
    stringsAsFactors = FALSE
  )
}

si_copy_number_endpoint_summaries <- function(sample_summary) {
  required <- c(
    "sample_id", "file", "initial_ploidy", "dose_mg_per_kg", "n_cells",
    "mean_postprocessed_copy_number_score"
  )
  if (!all(required %in% names(sample_summary)) ||
      anyDuplicated(sample_summary$file) ||
      !setequal(unique(sample_summary$initial_ploidy), c("2N", "4N"))) {
    stop("Copy-number-score sample summary has an invalid mouse-level contract",
         call. = FALSE)
  }
  samples <- sample_summary
  samples$initial_ploidy <- factor(
    samples$initial_ploidy, levels = c("2N", "4N")
  )
  samples$dose_factor <- factor(samples$dose_mg_per_kg)

  origin_split <- split(samples, samples$initial_ploidy, drop = TRUE)
  origin_summary <- do.call(rbind, lapply(origin_split, function(data) {
    data.frame(
      initial_ploidy = as.character(data$initial_ploidy[[1L]]),
      n_mice = nrow(data),
      mean_of_mouse_means = mean(data$mean_postprocessed_copy_number_score),
      min_mouse_mean = min(data$mean_postprocessed_copy_number_score),
      max_mouse_mean = max(data$mean_postprocessed_copy_number_score),
      analysis_type = "descriptive_only",
      summary_unit = "sequenced_mouse_CBS_file",
      formal_test_performed = FALSE,
      stringsAsFactors = FALSE
    )
  }))
  origin_summary <- origin_summary[
    match(c("2N", "4N"), origin_summary$initial_ploidy), , drop = FALSE
  ]
  rownames(origin_summary) <- NULL

  origin_dose_split <- split(
    samples,
    interaction(
      samples$initial_ploidy,
      samples$dose_factor,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  origin_dose_summary <- do.call(rbind, lapply(
    origin_dose_split,
    function(data) {
      data.frame(
        initial_ploidy = as.character(data$initial_ploidy[[1L]]),
        dose_mg_per_kg = data$dose_mg_per_kg[[1L]],
        n_mice = nrow(data),
        mean_of_mouse_means = mean(data$mean_postprocessed_copy_number_score),
        min_mouse_mean = min(data$mean_postprocessed_copy_number_score),
        max_mouse_mean = max(data$mean_postprocessed_copy_number_score),
        analysis_type = "descriptive_only",
        summary_unit = "sequenced_mouse_CBS_file",
        formal_test_performed = FALSE,
        stringsAsFactors = FALSE
      )
    }
  ))
  origin_dose_summary <- origin_dose_summary[order(
    match(origin_dose_summary$initial_ploidy, c("2N", "4N")),
    origin_dose_summary$dose_mg_per_kg
  ), , drop = FALSE]
  rownames(origin_dose_summary) <- NULL
  list(
    samples = samples[, required, drop = FALSE],
    origin_summary = origin_summary,
    origin_dose_summary = origin_dose_summary
  )
}
