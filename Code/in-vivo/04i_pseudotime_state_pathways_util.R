`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

pst_parse_scalar <- function(value, template) {
  if (is.logical(template)) return(tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y"))
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  as.character(value)
}

pst_parse_args <- function(argv, defaults) {
  out <- defaults
  if (length(argv) == 0L) return(out)
  for (arg in argv) {
    if (!grepl("^--", arg)) next
    key_value <- sub("^--", "", arg)
    if (grepl("=", key_value, fixed = TRUE)) {
      key <- sub("=.*$", "", key_value)
      value <- sub("^[^=]*=", "", key_value)
    } else {
      key <- key_value
      value <- "TRUE"
    }
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) {
      stop("Unknown argument: --", key, call. = FALSE)
    }
    out[[key]] <- pst_parse_scalar(value, out[[key]])
  }
  out
}

pst_required_packages <- function() {
  c(
    "Seurat", "Matrix", "yaml", "digest", "edgeR", "limma", "splines",
    "fgsea", "msigdbr", "dplyr", "readr", "tidyr", "ggplot2"
  )
}

pst_check_packages <- function(packages = pst_required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

pst_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

pst_clean_dir <- function(path, overwrite = FALSE) {
  if (dir.exists(path) && !overwrite) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) > 0L) {
      stop("Output directory is not empty. Use --overwrite=TRUE: ", path, call. = FALSE)
    }
  }
  if (dir.exists(path) && overwrite) unlink(path, recursive = TRUE, force = TRUE)
  pst_ensure_dir(path)
}

pst_write_csv <- function(x, path) {
  pst_ensure_dir(dirname(path))
  readr::write_csv(as.data.frame(x, stringsAsFactors = FALSE), path, na = "")
  invisible(path)
}

pst_write_lines <- function(x, path) {
  pst_ensure_dir(dirname(path))
  writeLines(as.character(x), path, useBytes = TRUE)
  invisible(path)
}

pst_file_checksum <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

pst_safe_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

pst_collection_label <- function(collection) {
  out <- c(
    "H" = "hallmark",
    "C2:CP:REACTOME" = "reactome",
    "C5:GO:BP" = "go_bp"
  )
  unname(out[collection] %||% pst_safe_name(collection))
}

pst_normalize_choice <- function(x) {
  x <- tolower(trimws(as.character(x)))
  gsub("-", "_", x)
}

pst_etp_threshold_specs <- function() {
  list(
    ETP_fixed_threshold_2_25 = list(
      method = "ETP_fixed_threshold_2_25",
      threshold_scheme = "fixed_threshold_2_25",
      threshold = 2.25,
      group_column = "ETP_fixed_threshold_2_25_group",
      factor_column = "ETP_fixed_threshold_2_25_factor",
      group_levels = c("ETP-lower", "ETP-higher")
    ),
    ETP_boundary_stress_threshold_2_375 = list(
      method = "ETP_boundary_stress_threshold_2_375",
      threshold_scheme = "boundary_stress_threshold_2_375",
      threshold = 2.375,
      group_column = "ETP_boundary_stress_threshold_2_375_group",
      factor_column = "ETP_boundary_stress_threshold_2_375_factor",
      group_levels = c("ETP-lower", "ETP-higher")
    ),
    ETP_reference_balanced_threshold_2_24 = list(
      method = "ETP_reference_balanced_threshold_2_24",
      threshold_scheme = "reference_balanced_threshold_2_24",
      threshold = 2.24,
      group_column = "ETP_reference_balanced_threshold_2_24_group",
      factor_column = "ETP_reference_balanced_threshold_2_24_factor",
      group_levels = c("ETP-lower", "ETP-higher")
    )
  )
}

pst_resolve_etp_specs <- function(etp_threshold) {
  specs <- pst_etp_threshold_specs()
  aliases <- c(
    all = "all",
    fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
    etp_fixed_threshold_2_25 = "ETP_fixed_threshold_2_25",
    boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
    etp_boundary_stress_threshold_2_375 = "ETP_boundary_stress_threshold_2_375",
    reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24",
    etp_reference_balanced_threshold_2_24 = "ETP_reference_balanced_threshold_2_24"
  )
  requested <- trimws(unlist(strsplit(as.character(etp_threshold), ",", fixed = TRUE)))
  requested <- pst_normalize_choice(requested[nzchar(requested)])
  if (length(requested) == 0L) requested <- "all"
  resolved <- unname(aliases[requested])
  if (anyNA(resolved)) {
    stop(
      "--etp-threshold must be one or more of: all, fixed_threshold_2_25, ",
      "boundary_stress_threshold_2_375, reference_balanced_threshold_2_24",
      call. = FALSE
    )
  }
  if ("all" %in% resolved) return(specs)
  specs[unique(resolved)]
}

pst_model_specs <- function(covariate_mode, etp_specs, include_full_etp_models = FALSE) {
  mode <- pst_normalize_choice(covariate_mode)
  valid <- c("initial_ploidy", "etp_group", "etp_continuous", "initial_plus_etp_group", "all")
  if (!mode %in% valid) {
    stop("--covariate-mode must be one of: ", paste(valid, collapse = ", "), call. = FALSE)
  }
  specs <- list()
  add_primary <- mode %in% c("initial_ploidy", "all")
  add_etp_group <- mode %in% c("etp_group", "all")
  add_etp_continuous <- mode %in% c("etp_continuous", "all")
  add_full <- mode %in% c("initial_plus_etp_group") || (mode == "all" && isTRUE(include_full_etp_models))
  if (add_primary) {
    specs$primary_initial_ploidy <- list(
      model_id = "primary_initial_ploidy",
      model_label = "Primary initial-ploidy-adjusted model",
      covariate_terms = c("initial_ploidy_factor"),
      covariate_mode = "initial_ploidy",
      etp_method = NA_character_,
      etp_threshold = NA_real_,
      is_primary = TRUE
    )
  }
  if (add_etp_continuous) {
    specs$etp_continuous_adjusted <- list(
      model_id = "etp_continuous_adjusted",
      model_label = "Continuous mean-ETP-adjusted model",
      covariate_terms = c("sample_mean_endpoint_ploidy_scaled"),
      covariate_mode = "etp_continuous",
      etp_method = "sample_mean_endpoint_ploidy",
      etp_threshold = NA_real_,
      is_primary = FALSE
    )
  }
  if (add_etp_group || add_full) {
    for (spec in etp_specs) {
      if (add_etp_group) {
        specs[[spec$method]] <- list(
          model_id = spec$method,
          model_label = paste0("ETP group-adjusted model: ", spec$threshold_scheme),
          covariate_terms = c(spec$factor_column),
          covariate_mode = "etp_group",
          etp_method = spec$method,
          etp_threshold = spec$threshold,
          is_primary = FALSE
        )
      }
      if (add_full) {
        full_id <- paste0("initial_ploidy_plus_", spec$method)
        specs[[full_id]] <- list(
          model_id = full_id,
          model_label = paste0("Initial ploidy plus ETP group model: ", spec$threshold_scheme),
          covariate_terms = c("initial_ploidy_factor", spec$factor_column),
          covariate_mode = "initial_plus_etp_group",
          etp_method = spec$method,
          etp_threshold = spec$threshold,
          is_primary = FALSE
        )
      }
    }
  }
  if (length(specs) == 0L) stop("No model specifications selected.", call. = FALSE)
  specs
}

pst_format_pathway_label <- function(pathway, collection = NA_character_) {
  x <- as.character(pathway)
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^GO_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

pst_interval_hit <- function(x, interval) {
  left <- if (isTRUE(interval$include_start)) x >= interval$start else x > interval$start
  right <- if (isTRUE(interval$include_end)) x <= interval$end else x < interval$end
  left & right
}

pst_read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$intervals)) stop("Config is missing intervals.", call. = FALSE)
  intervals <- lapply(names(cfg$intervals), function(name) {
    x <- cfg$intervals[[name]]
    x$name <- name
    x$start <- as.numeric(x$start)
    x$end <- as.numeric(x$end)
    x$include_start <- isTRUE(x$include_start)
    x$include_end <- isTRUE(x$include_end)
    x
  })
  names(intervals) <- names(cfg$intervals)
  cfg$interval_list <- intervals
  cfg
}

pst_intervals_table <- function(cfg) {
  rows <- lapply(cfg$interval_list, function(x) {
    data.frame(
      interval_id = x$name,
      role = x$role %||% NA_character_,
      label = x$label %||% x$name,
      start = x$start,
      end = x$end,
      include_start = x$include_start,
      include_end = x$include_end,
      source = x$source %||% NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pst_default_seurat_candidates <- function(results_root) {
  c(
    file.path(results_root, "03c_cluster_annotation/04_objects/integrated_sct_cca_seurat_cluster_final_annotation.rds"),
    file.path(results_root, "03b_manual_cluster_merge/objects/integrated_sct_cca_seurat_final_manual_merge.rds"),
    file.path(results_root, "03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds"),
    file.path(results_root, "02e_cluster_annotation/04_objects/integrated_sct_cca_seurat_cluster_final_annotation.rds")
  )
}

pst_resolve_seurat_rds <- function(explicit, results_root) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop("Explicit --seurat-rds does not exist: ", explicit, call. = FALSE)
    return(normalizePath(explicit, mustWork = TRUE))
  }
  candidates <- pst_default_seurat_candidates(results_root)
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("No Seurat RDS found. Tried: ", paste(candidates, collapse = "; "), call. = FALSE)
  }
  normalizePath(hit[[1L]], mustWork = TRUE)
}

pst_package_versions <- function(packages = pst_required_packages()) {
  data.frame(
    package = packages,
    version = vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1L)),
    stringsAsFactors = FALSE
  )
}

pst_read_cell_metadata <- function(path) {
  required <- c(
    "cell_id", "sample_id", "cluster", "initial_ploidy",
    "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "pseudotime"
  )
  meta <- readr::read_csv(path, show_col_types = FALSE)
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) stop("Cell metadata missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  meta$pseudotime <- suppressWarnings(as.numeric(meta$pseudotime))
  meta$gemcitabine_dose_mg_per_kg <- suppressWarnings(as.numeric(meta$gemcitabine_dose_mg_per_kg))
  if ("cell_ploidy" %in% names(meta)) meta$cell_ploidy <- suppressWarnings(as.numeric(meta$cell_ploidy))
  meta$cell_id <- as.character(meta$cell_id)
  meta$sample_id <- as.character(meta$sample_id)
  meta$initial_ploidy <- as.character(meta$initial_ploidy)
  meta
}

pst_unique_value <- function(x, label) {
  x <- unique(as.character(x[!is.na(x)]))
  x <- x[nzchar(x)]
  if (length(x) > 1L) stop("Non-unique sample-level value for ", label, ": ", paste(x, collapse = ", "), call. = FALSE)
  if (length(x) == 0L) NA_character_ else x[[1L]]
}

pst_sample_metadata <- function(meta) {
  split_meta <- split(meta, meta$sample_id)
  rows <- lapply(names(split_meta), function(sample_id) {
    x <- split_meta[[sample_id]]
    base <- data.frame(
      sample_id = sample_id,
      dose = unique(as.character(x$gemcitabine_dose))[1L],
      dose_mg = unique(x$gemcitabine_dose_mg_per_kg)[1L],
      initial_ploidy = unique(as.character(x$initial_ploidy))[1L],
      n_cells = nrow(x),
      mean_pseudotime = mean(x$pseudotime, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    optional <- intersect(
      c(
        "n_endpoint_ploidy_cells",
        "sample_mean_endpoint_ploidy",
        "sample_median_endpoint_ploidy",
        "sample_max_endpoint_ploidy",
        "sample_mean_endpoint_ploidy_scaled",
        vapply(pst_etp_threshold_specs(), `[[`, character(1L), "group_column")
      ),
      names(x)
    )
    for (column in optional) {
      if (is.numeric(x[[column]])) {
        values <- unique(x[[column]][is.finite(x[[column]])])
        base[[column]] <- if (length(values) == 0L) NA_real_ else values[[1L]]
      } else {
        base[[column]] <- pst_unique_value(x[[column]], paste(sample_id, column))
      }
    }
    base
  })
  out <- do.call(rbind, rows)
  rownames(out) <- out$sample_id
  out
}

pst_read_etp_compartment <- function(path, compartment) {
  if (!file.exists(path)) stop("Missing ETP metadata input: ", path, call. = FALSE)
  data <- readr::read_csv(path, show_col_types = FALSE)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  required <- c("cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "cell_ploidy")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) stop("ETP metadata input missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  data$cell_id <- as.character(data$cell_id)
  data$sample_id <- as.character(data$sample_id)
  data$initial_ploidy <- as.character(data$initial_ploidy)
  data$gemcitabine_dose <- as.character(data$gemcitabine_dose)
  data$gemcitabine_dose_mg_per_kg <- suppressWarnings(as.numeric(data$gemcitabine_dose_mg_per_kg))
  data$cell_ploidy <- suppressWarnings(as.numeric(data$cell_ploidy))
  data$etp_source_compartment <- compartment
  data[, c("cell_id", "sample_id", "initial_ploidy", "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "cell_ploidy", "etp_source_compartment"), drop = FALSE]
}

pst_build_etp_assignments <- function(cellcycle_path, noncell_path, etp_specs) {
  cellcycle <- pst_read_etp_compartment(cellcycle_path, "CellCycle")
  noncell <- pst_read_etp_compartment(noncell_path, "NonCellCycle")
  union_cells <- rbind(cellcycle, noncell)
  union_cells <- union_cells[is.finite(union_cells$cell_ploidy), , drop = FALSE]
  if (nrow(union_cells) == 0L) stop("No finite cell-level ETP values were found.", call. = FALSE)
  duplicated_ids <- unique(union_cells$cell_id[duplicated(union_cells$cell_id)])
  if (length(duplicated_ids) > 0L) {
    keep <- rep(TRUE, nrow(union_cells))
    for (cell_id in duplicated_ids) {
      idx <- which(union_cells$cell_id == cell_id)
      local <- union_cells[idx, , drop = FALSE]
      same_sample <- length(unique(local$sample_id)) == 1L
      same_ploidy <- length(unique(local$cell_ploidy)) == 1L
      if (!same_sample || !same_ploidy) stop("Conflicting duplicated ETP cell_id: ", cell_id, call. = FALSE)
      keep[idx[-1L]] <- FALSE
    }
    union_cells <- union_cells[keep, , drop = FALSE]
  }
  split_cells <- split(union_cells, union_cells$sample_id)
  rows <- lapply(names(split_cells), function(sample_id) {
    x <- split_cells[[sample_id]]
    data.frame(
      sample_id = sample_id,
      initial_ploidy = pst_unique_value(x$initial_ploidy, paste(sample_id, "initial_ploidy")),
      dose = pst_unique_value(x$gemcitabine_dose, paste(sample_id, "dose")),
      dose_mg = unique(x$gemcitabine_dose_mg_per_kg[is.finite(x$gemcitabine_dose_mg_per_kg)])[1L],
      n_endpoint_ploidy_cells = nrow(x),
      n_cellcycle_endpoint_ploidy_cells = sum(x$etp_source_compartment == "CellCycle"),
      n_noncellcycle_endpoint_ploidy_cells = sum(x$etp_source_compartment == "NonCellCycle"),
      sample_mean_endpoint_ploidy = mean(x$cell_ploidy),
      sample_median_endpoint_ploidy = median(x$cell_ploidy),
      sample_max_endpoint_ploidy = max(x$cell_ploidy),
      stringsAsFactors = FALSE
    )
  })
  assignments <- do.call(rbind, rows)
  assignments$sample_mean_endpoint_ploidy_scaled <- as.numeric(scale(assignments$sample_mean_endpoint_ploidy))
  for (spec in etp_specs) {
    group <- ifelse(assignments$sample_mean_endpoint_ploidy > spec$threshold, "ETP-higher", "ETP-lower")
    assignments[[spec$group_column]] <- factor(group, levels = spec$group_levels)
  }
  assignments[order(assignments$sample_id), , drop = FALSE]
}

pst_attach_etp_assignments <- function(meta, assignments) {
  idx <- match(meta$sample_id, assignments$sample_id)
  if (anyNA(idx)) {
    missing <- unique(meta$sample_id[is.na(idx)])
    stop("Missing ETP assignments for samples: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  add_cols <- setdiff(names(assignments), c("sample_id", "initial_ploidy", "dose", "dose_mg"))
  out <- meta
  for (column in add_cols) out[[column]] <- assignments[[column]][idx]
  out
}

pst_write_etp_audits <- function(assignments, etp_specs, qc_dir) {
  pst_write_csv(assignments, file.path(qc_dir, "etp_covariate_sample_table.csv"))
  rows <- list()
  for (spec in etp_specs) {
    group_col <- spec$group_column
    tab <- as.data.frame(
      xtabs(
        stats::as.formula(paste("~ initial_ploidy + dose +", group_col)),
        assignments
      ),
      stringsAsFactors = FALSE
    )
    names(tab)[names(tab) == group_col] <- "ETP_group"
    tab$etp_method <- spec$method
    tab$etp_threshold <- spec$threshold
    rows[[length(rows) + 1L]] <- tab
  }
  cross_tab <- if (length(rows) > 0L) do.call(rbind, rows) else data.frame()
  pst_write_csv(cross_tab, file.path(qc_dir, "etp_group_cross_tab.csv"))
  invisible(list(assignments = assignments, cross_tab = cross_tab))
}

pst_audit_metadata <- function(meta) {
  duplicate_cells <- meta[duplicated(meta$cell_id) | duplicated(meta$cell_id, fromLast = TRUE), , drop = FALSE]
  per_cell <- stats::aggregate(
    cbind(sample_id = meta$sample_id, dose = meta$gemcitabine_dose, initial_ploidy = meta$initial_ploidy),
    by = list(cell_id = meta$cell_id),
    FUN = function(x) length(unique(as.character(x)))
  )
  inconsistent <- per_cell[per_cell$sample_id > 1L | per_cell$dose > 1L | per_cell$initial_ploidy > 1L, , drop = FALSE]
  list(duplicate_cells = duplicate_cells, inconsistent_cells = inconsistent)
}

pst_check_pseudotime <- function(meta, tolerance = 1e-8) {
  bad <- meta[!is.finite(meta$pseudotime) | meta$pseudotime < -tolerance | meta$pseudotime > 1 + tolerance, , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop("Pseudotime contains non-finite or out-of-[0,1] values. See QC audit outputs.", call. = FALSE)
  }
  invisible(TRUE)
}

pst_load_counts <- function(seurat_rds, assay, counts_layer) {
  obj <- readRDS(seurat_rds)
  counts <- get_assay_matrix(obj, assay = assay, slot_name = counts_layer)
  if (is.null(counts)) {
    stop("Cannot read assay matrix: assay=", assay, ", layer/slot=", counts_layer, call. = FALSE)
  }
  if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  counts
}

pst_audit_counts <- function(counts) {
  values <- counts@x
  non_negative <- length(values) == 0L || all(values >= 0 & is.finite(values))
  integer_like <- length(values) == 0L || max(abs(values - round(values)), na.rm = TRUE) < 1e-8
  data.frame(
    matrix_class = paste(class(counts), collapse = ";"),
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    n_nonzero = length(values),
    sparse = inherits(counts, "sparseMatrix"),
    non_negative = non_negative,
    integer_like = integer_like,
    min_nonzero = if (length(values) > 0L) min(values) else NA_real_,
    max_nonzero = if (length(values) > 0L) max(values) else NA_real_,
    stringsAsFactors = FALSE
  )
}

pst_match_cells <- function(meta, counts, min_match_rate, qc_dir) {
  metadata_ids <- unique(meta$cell_id)
  expression_ids <- colnames(counts)
  meta_match <- metadata_ids %in% expression_ids
  expr_match <- expression_ids %in% metadata_ids
  audit <- data.frame(
    source = c("metadata", "expression"),
    n_ids = c(length(metadata_ids), length(expression_ids)),
    n_matched = c(sum(meta_match), sum(expr_match)),
    n_unmatched = c(sum(!meta_match), sum(!expr_match)),
    match_rate = c(mean(meta_match), mean(expr_match)),
    stringsAsFactors = FALSE
  )
  unmatched_meta_ids <- metadata_ids[!meta_match]
  unmatched_expr_ids <- expression_ids[!expr_match]
  unmatched_meta <- data.frame(
    cell_id = unmatched_meta_ids,
    source = rep("metadata", length(unmatched_meta_ids)),
    stringsAsFactors = FALSE
  )
  unmatched_expr <- data.frame(
    cell_id = unmatched_expr_ids,
    source = rep("expression", length(unmatched_expr_ids)),
    stringsAsFactors = FALSE
  )
  pst_write_csv(audit, file.path(qc_dir, "cell_id_join_audit.csv"))
  pst_write_csv(rbind(unmatched_meta, unmatched_expr), file.path(qc_dir, "cell_id_unmatched_ids.csv"))
  if (audit$match_rate[audit$source == "metadata"] < min_match_rate) {
    stop("Metadata-to-expression cell-ID match rate is below threshold: ", audit$match_rate[audit$source == "metadata"], call. = FALSE)
  }
  keep_meta <- meta$cell_id %in% expression_ids
  meta <- meta[keep_meta, , drop = FALSE]
  counts <- counts[, meta$cell_id, drop = FALSE]
  list(meta = meta, counts = counts, audit = audit)
}

pst_limit_for_smoke <- function(meta, counts, max_cells = 0L, max_genes = 0L, seed = 1L) {
  set.seed(seed)
  if (is.finite(max_cells) && max_cells > 0L && nrow(meta) > max_cells) {
    by_sample <- split(seq_len(nrow(meta)), meta$sample_id)
    n_each <- pmax(1L, floor(max_cells * lengths(by_sample) / nrow(meta)))
    idx <- unlist(Map(function(v, n) sample(v, min(length(v), n)), by_sample, n_each), use.names = FALSE)
    if (length(idx) < max_cells) {
      extra <- setdiff(seq_len(nrow(meta)), idx)
      idx <- c(idx, sample(extra, min(length(extra), max_cells - length(idx))))
    }
    idx <- sort(unique(idx))
    meta <- meta[idx, , drop = FALSE]
    counts <- counts[, meta$cell_id, drop = FALSE]
  }
  if (is.finite(max_genes) && max_genes > 0L && nrow(counts) > max_genes) {
    gene_nnz <- Matrix::rowSums(counts > 0)
    keep <- order(gene_nnz, decreasing = TRUE)[seq_len(max_genes)]
    counts <- counts[keep, , drop = FALSE]
  }
  list(meta = meta, counts = counts)
}

pst_bin_metadata <- function(meta, n_bins) {
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  bin_id <- findInterval(meta$pseudotime, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  bin_id <- pmin(pmax(bin_id, 1L), n_bins)
  meta$bin_id <- bin_id
  meta$sample_bin_id <- paste(meta$sample_id, sprintf("bin%02d", bin_id), sep = "__")
  bin_info <- data.frame(
    bin_id = seq_len(n_bins),
    bin_start = head(breaks, -1L),
    bin_end = tail(breaks, -1L),
    bin_midpoint = (head(breaks, -1L) + tail(breaks, -1L)) / 2,
    stringsAsFactors = FALSE
  )
  list(meta = meta, bin_info = bin_info)
}

pst_construct_pseudobulk <- function(meta, counts, n_bins, min_cells) {
  binned <- pst_bin_metadata(meta, n_bins)
  meta <- binned$meta
  bin_info <- binned$bin_info
  sample_meta <- pst_sample_metadata(meta)
  all_sample_bins <- expand.grid(
    sample_id = rownames(sample_meta),
    bin_id = seq_len(n_bins),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  all_sample_bins$sample_bin_id <- paste(all_sample_bins$sample_id, sprintf("bin%02d", all_sample_bins$bin_id), sep = "__")
  observed <- as.data.frame(table(meta$sample_bin_id), stringsAsFactors = FALSE)
  names(observed) <- c("sample_bin_id", "cell_count")
  sample_bin_meta <- merge(all_sample_bins, observed, by = "sample_bin_id", all.x = TRUE, sort = FALSE)
  sample_bin_meta$cell_count[is.na(sample_bin_meta$cell_count)] <- 0L
  sample_bin_meta <- merge(sample_bin_meta, bin_info, by = "bin_id", all.x = TRUE, sort = FALSE)
  sample_bin_meta <- merge(sample_bin_meta, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  sample_bin_meta <- sample_bin_meta[order(sample_bin_meta$sample_id, sample_bin_meta$bin_id), , drop = FALSE]
  sample_bin_meta$retained_for_model <- sample_bin_meta$cell_count >= min_cells
  sample_bin_meta$exclusion_reason <- ifelse(sample_bin_meta$retained_for_model, "", paste0("cell_count<", min_cells))

  group_factor <- factor(meta$sample_bin_id, levels = sample_bin_meta$sample_bin_id)
  design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)
  pseudobulk_counts <- counts %*% design
  colnames(pseudobulk_counts) <- levels(group_factor)
  rownames(pseudobulk_counts) <- rownames(counts)
  lib_size <- Matrix::colSums(pseudobulk_counts)
  sample_bin_meta$library_size <- as.numeric(lib_size[match(sample_bin_meta$sample_bin_id, names(lib_size))])
  sample_bin_meta$library_size[is.na(sample_bin_meta$library_size)] <- 0
  list(counts = pseudobulk_counts, metadata = sample_bin_meta, cell_metadata = meta)
}

pst_region_coverage <- function(meta, sample_bin_meta, cfg) {
  sample_meta <- unique(meta[, c("sample_id", "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "initial_ploidy"), drop = FALSE])
  rows <- list()
  for (interval in cfg$interval_list) {
    cell_hit <- pst_interval_hit(meta$pseudotime, interval)
    bin_hit <- pst_interval_hit(sample_bin_meta$bin_midpoint, interval)
    for (sample_id in unique(meta$sample_id)) {
      sm <- sample_meta[sample_meta$sample_id == sample_id, , drop = FALSE][1L, , drop = FALSE]
      sb <- sample_bin_meta[sample_bin_meta$sample_id == sample_id & bin_hit, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        sample_id = sample_id,
        interval_id = interval$name,
        role = interval$role %||% NA_character_,
        interval_start = interval$start,
        interval_end = interval$end,
        dose = sm$gemcitabine_dose,
        dose_mg = sm$gemcitabine_dose_mg_per_kg,
        initial_ploidy = sm$initial_ploidy,
        n_cells = sum(cell_hit & meta$sample_id == sample_id),
        n_bins = nrow(sb),
        n_retained_bins = sum(sb$retained_for_model),
        library_size = sum(sb$library_size, na.rm = TRUE),
        contributes = sum(sb$retained_for_model) > 0L,
        exclusion_reason = ifelse(sum(sb$retained_for_model) > 0L, "", "no_retained_bin_in_region"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

pst_check_primary_coverage <- function(coverage) {
  primary <- coverage[coverage$interval_id %in% c("primary_accumulated_state", "left_neighbor", "right_neighbor"), , drop = FALSE]
  usable <- stats::aggregate(contributes ~ sample_id + dose_mg + initial_ploidy, primary, FUN = all)
  usable <- usable[usable$contributes, , drop = FALSE]
  data.frame(
    n_contributing_mice = nrow(usable),
    has_control = any(usable$dose_mg == 0),
    has_treated = any(usable$dose_mg > 0),
    n_initial_ploidy_groups = length(unique(usable$initial_ploidy)),
    passes = nrow(usable) >= 6L && any(usable$dose_mg == 0) && any(usable$dose_mg > 0) && length(unique(usable$initial_ploidy)) >= 2L,
    stringsAsFactors = FALSE
  )
}

pst_prepare_design <- function(meta, spline_df, model_spec, include_dose = TRUE) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta$dose_mg_factor <- factor(meta$dose_mg)
  meta$initial_ploidy_factor <- factor(meta$initial_ploidy)
  for (spec in pst_etp_threshold_specs()) {
    if (spec$group_column %in% names(meta)) {
      meta[[spec$factor_column]] <- factor(as.character(meta[[spec$group_column]]), levels = spec$group_levels)
    }
  }
  if ("sample_mean_endpoint_ploidy" %in% names(meta)) {
    meta$sample_mean_endpoint_ploidy_scaled <- as.numeric(scale(meta$sample_mean_endpoint_ploidy))
  }
  basis <- splines::ns(meta$bin_midpoint, df = spline_df)
  colnames(basis) <- paste0("pt_spline", seq_len(ncol(basis)))
  design_data <- cbind(meta, as.data.frame(basis))
  covariate_terms <- model_spec$covariate_terms %||% character(0L)
  if (!isTRUE(include_dose)) covariate_terms <- setdiff(covariate_terms, "dose_mg_factor")
  model_terms <- c(colnames(basis), if (isTRUE(include_dose)) "dose_mg_factor", covariate_terms)
  missing_terms <- setdiff(model_terms, names(design_data))
  if (length(missing_terms) > 0L) {
    stop("Model ", model_spec$model_id, " missing design terms: ", paste(missing_terms, collapse = ", "), call. = FALSE)
  }
  formula_text <- paste("~", paste(model_terms, collapse = " + "))
  design <- stats::model.matrix(stats::as.formula(formula_text), design_data)
  list(
    design = design,
    design_data = design_data,
    basis = basis,
    formula = stats::as.formula(formula_text),
    model_terms = model_terms,
    original_design_columns = colnames(design),
    model_spec = model_spec
  )
}

pst_new_design_rows <- function(model, pseudotime_grid, reference_meta) {
  new_basis <- predict(model$basis, newx = pseudotime_grid)
  colnames(new_basis) <- colnames(model$basis)
  new_data <- model$design_data[rep(1L, length(pseudotime_grid)), , drop = FALSE]
  rownames(new_data) <- NULL
  new_data$bin_midpoint <- pseudotime_grid
  for (column in colnames(new_basis)) new_data[[column]] <- new_basis[, column]
  design <- stats::model.matrix(model$formula, new_data)
  missing <- setdiff(colnames(model$design), colnames(design))
  for (column in missing) design <- cbind(design, setNames(data.frame(0), column))
  design <- design[, colnames(model$design), drop = FALSE]
  design
}

pst_contrast_vector <- function(model, cfg, contrast_name, grid_size = 501L) {
  if (!is.null(model$design) && !is.null(model$design$basis)) {
    model <- model$design
  }
  grid <- seq(0, 1, length.out = grid_size)
  reference_meta <- model$design_data[1L, , drop = FALSE]
  grid_design <- pst_new_design_rows(model, grid, reference_meta)
  intervals <- cfg$interval_list
  mean_design <- function(mask) {
    if (!any(mask)) stop("Contrast grid has no points for requested interval.", call. = FALSE)
    colMeans(grid_design[mask, , drop = FALSE])
  }
  if (identical(contrast_name, "primary_adjacent_state")) {
    primary <- pst_interval_hit(grid, intervals$primary_accumulated_state)
    left <- pst_interval_hit(grid, intervals$left_neighbor)
    right <- pst_interval_hit(grid, intervals$right_neighbor)
    vec <- mean_design(primary) - 0.5 * mean_design(left) - 0.5 * mean_design(right)
  } else if (identical(contrast_name, "full_trajectory_state_specificity")) {
    primary <- pst_interval_hit(grid, intervals$primary_accumulated_state)
    vec <- mean_design(primary) - mean_design(!primary)
  } else {
    interval <- intervals[[contrast_name]]
    if (is.null(interval)) stop("Unknown contrast interval: ", contrast_name, call. = FALSE)
    inside <- pst_interval_hit(grid, interval)
    vec <- mean_design(inside) - mean_design(!inside)
  }
  vec
}

pst_fit_pseudotime_model <- function(pb_counts, sample_bin_meta, spline_df, min_primary_mice, include_dose = TRUE, model_spec = NULL) {
  model_spec <- model_spec %||% list(
    model_id = "primary_initial_ploidy",
    covariate_terms = c("initial_ploidy_factor"),
    covariate_mode = "initial_ploidy",
    etp_method = NA_character_,
    etp_threshold = NA_real_
  )
  keep_obs <- sample_bin_meta$retained_for_model & sample_bin_meta$library_size > 0
  meta <- sample_bin_meta[keep_obs, , drop = FALSE]
  counts <- pb_counts[, meta$sample_bin_id, drop = FALSE]
  y0 <- edgeR::DGEList(counts = counts)
  cpm0 <- edgeR::cpm(y0)
  min_obs <- max(2L, min_primary_mice)
  keep_gene <- rowSums(cpm0 > 1) >= min_obs
  if (sum(keep_gene) < 10L) stop("Too few genes pass expression filtering.", call. = FALSE)
  counts <- counts[keep_gene, , drop = FALSE]
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design_bundle <- pst_prepare_design(meta, spline_df, model_spec = model_spec, include_dose = include_dose)
  original_design_columns <- colnames(design_bundle$design)
  qr_rank <- qr(design_bundle$design)$rank
  dropped_design_columns <- character(0L)
  if (qr_rank < ncol(design_bundle$design)) {
    keep_cols <- qr(design_bundle$design)$pivot[seq_len(qr_rank)]
    retained <- colnames(design_bundle$design)[sort(keep_cols)]
    dropped_design_columns <- setdiff(colnames(design_bundle$design), retained)
    design_bundle$design <- design_bundle$design[, retained, drop = FALSE]
  }
  v <- limma::voom(y, design_bundle$design, plot = FALSE)
  block <- meta$sample_id
  dup <- limma::duplicateCorrelation(v, design_bundle$design, block = block)
  fit <- limma::lmFit(v, design_bundle$design, block = block, correlation = dup$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  design_bundle$design <- design_bundle$design[, colnames(fit$coefficients), drop = FALSE]
  list(
    fit = fit,
    voom = v,
    dge = y,
    metadata = meta,
    design = design_bundle,
    correlation = dup$consensus.correlation,
    retained_genes = rownames(counts),
    include_dose = include_dose,
    model_spec = model_spec,
    original_design_columns = original_design_columns,
    retained_design_columns = colnames(design_bundle$design),
    dropped_design_columns = dropped_design_columns,
    design_rank = qr_rank,
    design_ncol_original = length(original_design_columns),
    design_n_observations = nrow(design_bundle$design)
  )
}

pst_design_audit <- function(model_fit) {
  data.frame(
    model_id = model_fit$model_spec$model_id %||% NA_character_,
    covariate_mode = model_fit$model_spec$covariate_mode %||% NA_character_,
    etp_method = model_fit$model_spec$etp_method %||% NA_character_,
    etp_threshold = model_fit$model_spec$etp_threshold %||% NA_real_,
    include_dose = isTRUE(model_fit$include_dose),
    n_observations = model_fit$design_n_observations,
    n_design_columns_original = model_fit$design_ncol_original,
    design_rank = model_fit$design_rank,
    rank_deficient = length(model_fit$dropped_design_columns) > 0L,
    retained_design_columns = paste(model_fit$retained_design_columns, collapse = ";"),
    dropped_design_columns = paste(model_fit$dropped_design_columns, collapse = ";"),
    duplicate_correlation = model_fit$correlation,
    stringsAsFactors = FALSE
  )
}

pst_contrast_table <- function(model_fit, contrast_vector, contrast_id) {
  common <- intersect(names(contrast_vector), colnames(model_fit$fit$coefficients))
  vec <- rep(0, ncol(model_fit$fit$coefficients))
  names(vec) <- colnames(model_fit$fit$coefficients)
  vec[common] <- contrast_vector[common]
  fit2 <- limma::contrasts.fit(model_fit$fit, contrasts = vec)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  tt <- limma::topTable(fit2, coef = 1L, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt$gene_symbol <- clean_gene_symbols(tt$gene)
  out <- data.frame(
    gene = tt$gene,
    gene_symbol = tt$gene_symbol,
    contrast_id = contrast_id,
    fitted_log_expr_difference = tt$logFC,
    average_expression = tt$AveExpr,
    t_statistic = tt$t,
    p_value = tt$P.Value,
    fdr = tt$adj.P.Val,
    B = tt$B,
    stringsAsFactors = FALSE
  )
  out[order(out$p_value, -abs(out$t_statistic), out$gene), , drop = FALSE]
}

pst_resolve_gene_symbols <- function(primary_table) {
  df <- primary_table
  df$gene_symbol[is.na(df$gene_symbol) | df$gene_symbol == ""] <- df$gene[is.na(df$gene_symbol) | df$gene_symbol == ""]
  df$abs_t <- abs(df$t_statistic)
  df <- df[order(df$gene_symbol, -df$abs_t, df$p_value), , drop = FALSE]
  df$resolution_rank <- ave(df$abs_t, df$gene_symbol, FUN = function(x) rank(-x, ties.method = "first"))
  df$retained_for_gsea <- df$resolution_rank == 1L
  df[, c("gene", "gene_symbol", "contrast_id", "t_statistic", "p_value", "fdr", "resolution_rank", "retained_for_gsea"), drop = FALSE]
}

pst_ranked_stats <- function(contrast_table, symbol_resolution = NULL) {
  df <- contrast_table
  df$gene_symbol[is.na(df$gene_symbol) | df$gene_symbol == ""] <- df$gene[is.na(df$gene_symbol) | df$gene_symbol == ""]
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" & is.finite(df$t_statistic), , drop = FALSE]
  df$abs_t <- abs(df$t_statistic)
  df <- df[order(df$gene_symbol, -df$abs_t, df$p_value), , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]
  stats <- df$t_statistic
  names(stats) <- df$gene_symbol
  sort(stats, decreasing = TRUE)
}

pst_fetch_gene_sets <- function(collections, species = "Homo sapiens") {
  fetch_one <- function(collection) {
    parts <- strsplit(collection, ":", fixed = TRUE)[[1L]]
    if (identical(collection, "H")) {
      df <- tryCatch(
        msigdbr::msigdbr(species = species, collection = "H"),
        error = function(e) msigdbr::msigdbr(species = species, category = "H")
      )
    } else {
      category <- parts[[1L]]
      subcategory <- paste(parts[-1L], collapse = ":")
      df <- tryCatch(
        msigdbr::msigdbr(species = species, collection = category, subcollection = subcategory),
        error = function(e1) {
          tryCatch(
            msigdbr::msigdbr(species = species, category = category, subcategory = subcategory),
            error = function(e2) {
              all_sets <- tryCatch(
                msigdbr::msigdbr(species = species, collection = category),
                error = function(e3) msigdbr::msigdbr(species = species, category = category)
              )
              sub_col <- if ("gs_subcollection" %in% names(all_sets)) "gs_subcollection" else if ("gs_subcat" %in% names(all_sets)) "gs_subcat" else NA_character_
              if (is.na(sub_col)) stop("Cannot filter msigdbr subcollection for ", collection, call. = FALSE)
              all_sets[as.character(all_sets[[sub_col]]) == subcategory, , drop = FALSE]
            }
          )
        }
      )
    }
    if (nrow(df) == 0L) stop("No gene sets returned for collection: ", collection, call. = FALSE)
    df$gene_symbol <- clean_gene_symbols(df$gene_symbol)
    sets <- split(df$gene_symbol, df$gs_name)
    sets <- lapply(sets, function(x) unique(x[!is.na(x) & x != ""]))
    list(collection = collection, label = pst_collection_label(collection), sets = sets)
  }
  collections <- trimws(unlist(strsplit(collections, ",", fixed = TRUE)))
  collections <- collections[nzchar(collections)]
  out <- lapply(collections, fetch_one)
  names(out) <- vapply(out, `[[`, character(1L), "collection")
  out
}

pst_run_fgsea_collection <- function(stats, collection_obj, ranking_id, min_size, max_size, nperm_simple, seed) {
  set.seed(seed)
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  if (length(stats) == 0L) return(data.frame())
  res <- tryCatch(
    {
      if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
        fgsea::fgseaMultilevel(
          pathways = collection_obj$sets,
          stats = stats,
          minSize = min_size,
          maxSize = max_size,
          nPermSimple = nperm_simple,
          nproc = 1L,
          eps = 0
        )
      } else {
        fgsea::fgsea(
          pathways = collection_obj$sets,
          stats = stats,
          minSize = min_size,
          maxSize = max_size,
          nperm = nperm_simple,
          nproc = 1L
        )
      }
    },
    error = function(e) {
      warning("fgsea failed for ", collection_obj$collection, ": ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(res) || nrow(res) == 0L) return(data.frame())
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  if ("leadingEdge" %in% names(out)) {
    out$leading_edge <- vapply(out$leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1L))
    out$leadingEdge <- NULL
  } else {
    out$leading_edge <- ""
  }
  out$collection <- collection_obj$collection
  out$collection_label <- collection_obj$label
  out$ranking_id <- ranking_id
  out$pathway_label <- pst_format_pathway_label(out$pathway, collection_obj$collection)
  out$direction <- ifelse(out$NES >= 0, "positive", "negative")
  out <- out[order(out$padj, -abs(out$NES), out$pathway), , drop = FALSE]
  rownames(out) <- NULL
  out
}

pst_run_all_gsea <- function(stats, gene_sets, ranking_id, min_size, max_size, nperm_simple, seed) {
  rows <- lapply(seq_along(gene_sets), function(i) {
    pst_run_fgsea_collection(stats, gene_sets[[i]], ranking_id, min_size, max_size, nperm_simple, seed + i)
  })
  do.call(rbind, rows)
}

pst_write_gsea_tables <- function(primary_gsea, secondary_gsea, gsea_dir) {
  for (df in list(primary_gsea, secondary_gsea)) {
    if (is.null(df) || nrow(df) == 0L) next
    ranking_id <- unique(df$ranking_id)
    for (collection_label in unique(df$collection_label)) {
      local <- df[df$collection_label == collection_label, , drop = FALSE]
      file <- paste0(collection_label, "_", ranking_id[[1L]], "_gsea.csv")
      pst_write_csv(local, file.path(gsea_dir, file))
    }
  }
  pst_write_csv(primary_gsea, file.path(gsea_dir, "all_collections_primary_adjacent_state_gsea.csv"))
  pst_write_csv(secondary_gsea, file.path(gsea_dir, "all_collections_full_trajectory_state_specificity_gsea.csv"))
}

pst_leading_edge_table <- function(gsea_df) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0L) return(data.frame())
  rows <- lapply(seq_len(nrow(gsea_df)), function(i) {
    genes <- trimws(unlist(strsplit(gsea_df$leading_edge[[i]], ";", fixed = TRUE)))
    genes <- genes[nzchar(genes)]
    if (length(genes) == 0L) return(NULL)
    data.frame(
      collection = gsea_df$collection[[i]],
      collection_label = gsea_df$collection_label[[i]],
      ranking_id = gsea_df$ranking_id[[i]],
      pathway = gsea_df$pathway[[i]],
      pathway_label = gsea_df$pathway_label[[i]],
      NES = gsea_df$NES[[i]],
      padj = gsea_df$padj[[i]],
      leading_edge_gene = genes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

pst_fitted_grid <- function(model_fit, grid_size = 501L) {
  grid <- seq(0, 1, length.out = grid_size)
  design <- pst_new_design_rows(model_fit$design, grid, model_fit$design$design_data[1L, , drop = FALSE])
  coef <- model_fit$fit$coefficients[, colnames(design), drop = FALSE]
  fitted <- coef %*% t(design)
  colnames(fitted) <- sprintf("%.6f", grid)
  list(grid = grid, fitted = fitted)
}

pst_pathway_activity <- function(fitted_grid, gsea_df, max_pathways = Inf) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0L) return(data.frame())
  fg <- fitted_grid$fitted
  gene_symbols <- clean_gene_symbols(rownames(fg))
  symbol_to_row <- split(seq_along(gene_symbols), gene_symbols)
  source <- gsea_df[order(gsea_df$collection_label, gsea_df$padj, -abs(gsea_df$NES)), , drop = FALSE]
  if (is.finite(max_pathways)) source <- head(source, max_pathways)
  rows <- lapply(seq_len(nrow(source)), function(i) {
    genes <- trimws(unlist(strsplit(source$leading_edge[[i]], ";", fixed = TRUE)))
    genes <- unique(genes[nzchar(genes)])
    row_idx <- unique(unlist(symbol_to_row[intersect(genes, names(symbol_to_row))], use.names = FALSE))
    if (length(row_idx) < 2L) return(NULL)
    mat <- fg[row_idx, , drop = FALSE]
    mat <- t(scale(t(mat)))
    activity <- colMeans(mat, na.rm = TRUE)
    peak <- fitted_grid$grid[which.max(activity)]
    data.frame(
      collection = source$collection[[i]],
      collection_label = source$collection_label[[i]],
      pathway = source$pathway[[i]],
      pathway_label = source$pathway_label[[i]],
      ranking_id = source$ranking_id[[i]],
      NES = source$NES[[i]],
      padj = source$padj[[i]],
      pseudotime = fitted_grid$grid,
      standardized_activity = as.numeric(activity),
      n_leading_edge_genes_used = length(row_idx),
      peak_pseudotime = peak,
      peaks_inside_primary = peak >= 0.30 & peak <= 0.49,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

pst_density_difference_data <- function(meta, cfg, grid_size = 501L) {
  grid <- seq(0, 1, length.out = grid_size)
  sample_ids <- sort(unique(meta$sample_id))
  pooled <- meta$pseudotime[is.finite(meta$pseudotime)]
  bw <- stats::bw.nrd0(pooled)
  if (!is.finite(bw) || bw <= 0) bw <- 0.05
  rows <- lapply(sample_ids, function(sample_id) {
    x <- meta$pseudotime[meta$sample_id == sample_id]
    d <- stats::density(x, bw = bw, from = 0, to = 1, n = grid_size, na.rm = TRUE)
    sm <- meta[match(sample_id, meta$sample_id), , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      dose_mg = sm$gemcitabine_dose_mg_per_kg[[1L]],
      treatment_group = ifelse(sm$gemcitabine_dose_mg_per_kg[[1L]] == 0, "control", "treated"),
      pseudotime = d$x,
      density = d$y,
      common_bandwidth = bw,
      stringsAsFactors = FALSE
    )
  })
  sample_density <- do.call(rbind, rows)
  mean_density <- stats::aggregate(density ~ treatment_group + pseudotime + common_bandwidth, sample_density, mean)
  control <- mean_density[mean_density$treatment_group == "control", c("pseudotime", "density")]
  treated <- mean_density[mean_density$treatment_group == "treated", c("pseudotime", "density")]
  diff <- merge(treated, control, by = "pseudotime", suffixes = c("_treated", "_control"))
  diff$density_difference_treated_minus_control <- diff$density_treated - diff$density_control
  list(sample_density = sample_density, mean_density = mean_density, difference = diff)
}

pst_select_top_pathways <- function(gsea_df, top_each = 5L) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0L) return(data.frame())
  rows <- lapply(split(gsea_df, gsea_df$collection_label), function(df) {
    pos <- df[df$NES > 0 & is.finite(df$padj), , drop = FALSE]
    neg <- df[df$NES < 0 & is.finite(df$padj), , drop = FALSE]
    pos <- head(pos[order(pos$padj, -pos$NES), , drop = FALSE], top_each)
    neg <- head(neg[order(neg$padj, neg$NES), , drop = FALSE], top_each)
    rbind(pos, neg)
  })
  do.call(rbind, rows)
}

pst_plot_density <- function(density_data, intervals, path) {
  rects <- intervals[intervals$interval_id %in% c("primary_accumulated_state", "broad_sensitivity", "strict_support_sensitivity", "depleted_negative_control"), , drop = FALSE]
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = rects,
      ggplot2::aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = interval_id),
      alpha = 0.10,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_line(
      data = density_data$mean_density,
      ggplot2::aes(x = pseudotime, y = density, color = treatment_group),
      linewidth = 0.8
    ) +
    ggplot2::geom_line(
      data = density_data$difference,
      ggplot2::aes(x = pseudotime, y = density_difference_treated_minus_control),
      color = "black",
      linewidth = 0.6
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    ggplot2::labs(
      title = "CellCycle pseudotime density with frozen state intervals",
      subtitle = "Equal-mouse mean density; black line is treated minus control",
      x = "Pseudotime",
      y = "Density / density difference",
      color = "Group",
      fill = "Frozen interval"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(path, p, width = 8, height = 5.5)
  invisible(p)
}

pst_plot_gsea_bar <- function(gsea_df, path, title, top_each = 10L) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0L) return(invisible(NULL))
  df <- rbind(
    head(gsea_df[gsea_df$NES > 0, , drop = FALSE][order(gsea_df$padj[gsea_df$NES > 0], -gsea_df$NES[gsea_df$NES > 0]), , drop = FALSE], top_each),
    head(gsea_df[gsea_df$NES < 0, , drop = FALSE][order(gsea_df$padj[gsea_df$NES < 0], gsea_df$NES[gsea_df$NES < 0]), , drop = FALSE], top_each)
  )
  if (nrow(df) == 0L) return(invisible(NULL))
  df$pathway_label <- factor(df$pathway_label, levels = rev(unique(df$pathway_label[order(df$NES)])))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pathway_label, y = NES, fill = padj)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_viridis_c(option = "C", direction = -1, na.value = "grey80") +
    ggplot2::labs(title = title, x = NULL, y = "NES", fill = "FDR") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(path, p, width = 8, height = 6)
  invisible(p)
}

pst_plot_activity_heatmap <- function(activity, path) {
  if (is.null(activity) || nrow(activity) == 0L) return(invisible(NULL))
  top <- unique(activity[, c("collection_label", "pathway", "pathway_label", "NES", "padj"), drop = FALSE])
  top <- pst_select_top_pathways(top, top_each = 4L)
  df <- activity[activity$pathway %in% top$pathway, , drop = FALSE]
  df$pathway_label <- factor(df$pathway_label, levels = rev(unique(top$pathway_label)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pseudotime, y = pathway_label, fill = standardized_activity)) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(xintercept = c(0.30, 0.49), color = "black", linetype = "22", linewidth = 0.3) +
    ggplot2::facet_grid(collection_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0) +
    ggplot2::labs(title = "Pathway activity over CellCycle pseudotime", x = "Pseudotime", y = NULL, fill = "Activity") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(path, p, width = 9, height = 8)
  invisible(p)
}

pst_plot_activity_curves <- function(activity, path) {
  if (is.null(activity) || nrow(activity) == 0L) return(invisible(NULL))
  top <- unique(activity[, c("collection_label", "pathway", "pathway_label", "NES", "padj"), drop = FALSE])
  top <- pst_select_top_pathways(top, top_each = 2L)
  df <- activity[activity$pathway %in% top$pathway, , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = pseudotime, y = standardized_activity, color = pathway_label)) +
    ggplot2::annotate("rect", xmin = 0.30, xmax = 0.49, ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.08) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::facet_wrap(~ collection_label, scales = "free_y") +
    ggplot2::labs(title = "Top pathway fitted activity curves", x = "Pseudotime", y = "Standardized activity", color = "Pathway") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(path, p, width = 10, height = 7)
  invisible(p)
}

pst_leading_gene_curve_data <- function(fitted_grid, leading_edge, primary_gsea, max_pathways = 3L, max_genes_per_pathway = 8L) {
  if (is.null(leading_edge) || nrow(leading_edge) == 0L) return(data.frame())
  top <- primary_gsea[primary_gsea$NES > 0, , drop = FALSE]
  top <- head(top[order(top$padj, -top$NES), , drop = FALSE], max_pathways)
  if (nrow(top) == 0L) return(data.frame())
  fg <- fitted_grid$fitted
  gene_symbols <- clean_gene_symbols(rownames(fg))
  rows <- lapply(seq_len(nrow(top)), function(i) {
    genes <- leading_edge$leading_edge_gene[leading_edge$pathway == top$pathway[[i]]]
    genes <- head(unique(genes), max_genes_per_pathway)
    idx <- match(genes, gene_symbols)
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0L) return(NULL)
    mat <- fg[idx, , drop = FALSE]
    mat <- t(scale(t(mat)))
    do.call(rbind, lapply(seq_along(idx), function(j) {
      data.frame(
        collection = top$collection[[i]],
        collection_label = top$collection_label[[i]],
        pathway = top$pathway[[i]],
        pathway_label = top$pathway_label[[i]],
        gene_symbol = gene_symbols[[idx[[j]]]],
        pseudotime = fitted_grid$grid,
        standardized_fitted_expression = as.numeric(mat[j, ]),
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}

pst_plot_gene_curves <- function(gene_curves, path) {
  if (is.null(gene_curves) || nrow(gene_curves) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(gene_curves, ggplot2::aes(x = pseudotime, y = standardized_fitted_expression, color = gene_symbol)) +
    ggplot2::annotate("rect", xmin = 0.30, xmax = 0.49, ymin = -Inf, ymax = Inf, fill = "grey50", alpha = 0.08) +
    ggplot2::geom_line(linewidth = 0.6) +
    ggplot2::facet_wrap(~ pathway_label, scales = "free_y") +
    ggplot2::labs(title = "Leading-edge gene fitted curves", x = "Pseudotime", y = "Standardized fitted expression", color = "Gene") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom")
  ggplot2::ggsave(path, p, width = 10, height = 7)
  invisible(p)
}

pst_plot_coverage <- function(coverage, path) {
  p <- ggplot2::ggplot(coverage, ggplot2::aes(x = sample_id, y = n_cells, fill = contributes)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ interval_id, scales = "free_y") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Sample coverage in frozen pseudotime regions", x = "Sample", y = "Cells", fill = "Usable bins") +
    ggplot2::theme_bw(base_size = 9)
  ggplot2::ggsave(path, p, width = 10, height = 8)
  invisible(p)
}

pst_plot_robustness <- function(robustness, path) {
  if (is.null(robustness) || nrow(robustness) == 0L) return(invisible(NULL))
  df <- robustness[is.finite(robustness$NES), , drop = FALSE]
  top_pathways <- head(unique(df$pathway[order(df$primary_padj, -abs(df$primary_NES))]), 30L)
  df <- df[df$pathway %in% top_pathways, , drop = FALSE]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = analysis_id, y = NES, color = collection_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::facet_wrap(~ pathway_label, scales = "free_y") +
    ggplot2::labs(title = "Pathway robustness across sensitivity analyses", x = "Sensitivity analysis", y = "NES", color = "Collection") +
    ggplot2::theme_bw(base_size = 8) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(path, p, width = 11, height = 9)
  invisible(p)
}

pst_sensitivity_one <- function(analysis_id, pb_counts, sample_bin_meta, cfg, gene_sets, args, include_dose = TRUE, n_bins = NULL, spline_df = NULL, base_cell_meta = NULL, base_counts = NULL, model_spec = NULL) {
  if (!is.null(n_bins)) {
    pb <- pst_construct_pseudobulk(base_cell_meta, base_counts, n_bins, args$min_cells_per_sample_bin)
    local_counts <- pb$counts
    local_meta <- pb$metadata
  } else {
    local_counts <- pb_counts
    local_meta <- sample_bin_meta
  }
  spline_df <- spline_df %||% args$spline_df
  coverage <- NULL
  min_mice <- 2L
  fit <- pst_fit_pseudotime_model(local_counts, local_meta, spline_df, min_mice, include_dose = include_dose, model_spec = model_spec)
  contrast_id <- if (analysis_id %in% names(cfg$interval_list)) analysis_id else "primary_adjacent_state"
  contrast <- pst_contrast_vector(fit, cfg, contrast_id, args$grid_size)
  genes <- pst_contrast_table(fit, contrast, contrast_id)
  stats <- pst_ranked_stats(genes)
  gsea <- pst_run_all_gsea(stats, gene_sets, analysis_id, args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  gsea
}

pst_run_sensitivity <- function(pb_counts, sample_bin_meta, cfg, gene_sets, primary_gsea, args, base_cell_meta, base_counts, model_spec) {
  analyses <- list(
    broad_sensitivity = list(interval = "broad_sensitivity"),
    strict_support_sensitivity = list(interval = "strict_support_sensitivity"),
    depleted_negative_control = list(interval = "depleted_negative_control"),
    no_dose_nuisance = list(include_dose = FALSE),
    bins_15 = list(n_bins = 15L),
    bins_25 = list(n_bins = 25L),
    spline_df_4 = list(spline_df = 4L),
    spline_df_6 = list(spline_df = 6L)
  )
  sensitivity_ids <- names(analyses)
  sensitivity_workers <- suppressWarnings(as.integer(args$workers %||% 1L))
  if (!is.finite(sensitivity_workers) || sensitivity_workers < 1L) sensitivity_workers <- 1L
  sensitivity_workers <- max(1L, min(length(sensitivity_ids), sensitivity_workers))
  message("Running sensitivity for ", model_spec$model_id %||% "model", " with workers=", sensitivity_workers)
  worker_fun <- function(id) {
    spec <- analyses[[id]]
    gsea <- tryCatch(
      pst_sensitivity_one(
        id,
        pb_counts,
        sample_bin_meta,
        cfg,
        gene_sets,
        args,
        include_dose = spec$include_dose %||% TRUE,
        n_bins = spec$n_bins %||% NULL,
        spline_df = spec$spline_df %||% NULL,
        base_cell_meta = base_cell_meta,
        base_counts = base_counts,
        model_spec = model_spec
      ),
      error = function(e) {
        warning("Sensitivity failed for ", id, ": ", conditionMessage(e), call. = FALSE)
        data.frame()
      }
    )
    if (nrow(gsea) == 0L) return(NULL)
    gsea$analysis_id <- id
    gsea
  }
  if (.Platform$OS.type == "unix" && sensitivity_workers > 1L) {
    rows <- parallel::mclapply(sensitivity_ids, worker_fun, mc.cores = sensitivity_workers)
  } else {
    rows <- lapply(sensitivity_ids, worker_fun)
  }
  sens <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(sens) || nrow(sens) == 0L) return(data.frame())
  primary <- primary_gsea[, c("collection", "collection_label", "pathway", "pathway_label", "NES", "padj"), drop = FALSE]
  names(primary)[names(primary) == "NES"] <- "primary_NES"
  names(primary)[names(primary) == "padj"] <- "primary_padj"
  out <- merge(sens, primary, by = c("collection", "collection_label", "pathway", "pathway_label"), all.x = TRUE, sort = FALSE)
  out$sign_concordant_with_primary <- sign(out$NES) == sign(out$primary_NES)
  out
}

pst_run_leave_one_out <- function(pb_counts, sample_bin_meta, cfg, gene_sets, args, model_spec) {
  sample_ids <- sort(unique(sample_bin_meta$sample_id))
  loo_workers <- suppressWarnings(as.integer(args$workers %||% 1L))
  if (!is.finite(loo_workers) || loo_workers < 1L) loo_workers <- 1L
  loo_workers <- max(1L, min(length(sample_ids), loo_workers))
  message("Running leave-one-out for ", model_spec$model_id %||% "model", " with workers=", loo_workers)
  worker_fun <- function(sample_id) {
    local_meta <- sample_bin_meta[sample_bin_meta$sample_id != sample_id, , drop = FALSE]
    local_counts <- pb_counts[, local_meta$sample_bin_id, drop = FALSE]
    gsea <- tryCatch(
      pst_sensitivity_one(
        paste0("leave_one_out_", sample_id),
        local_counts,
        local_meta,
        cfg,
        gene_sets,
        args,
        include_dose = TRUE,
        model_spec = model_spec
      ),
      error = function(e) {
        warning("Leave-one-out failed for ", sample_id, ": ", conditionMessage(e), call. = FALSE)
        data.frame()
      }
    )
    if (nrow(gsea) == 0L) return(data.frame())
    gsea$left_out_sample_id <- sample_id
    gsea
  }
  if (.Platform$OS.type == "unix" && loo_workers > 1L) {
    rows <- parallel::mclapply(sample_ids, worker_fun, mc.cores = loo_workers)
  } else {
    rows <- lapply(sample_ids, worker_fun)
  }
  do.call(rbind, rows)
}

pst_leave_one_out_stability <- function(loo, primary_gsea) {
  if (is.null(loo) || nrow(loo) == 0L) return(data.frame())
  primary <- primary_gsea[, c("collection", "collection_label", "pathway", "pathway_label", "NES", "padj"), drop = FALSE]
  names(primary)[names(primary) == "NES"] <- "primary_NES"
  names(primary)[names(primary) == "padj"] <- "primary_padj"
  loo$loo_rank <- ave(
    loo$padj,
    loo$left_out_sample_id,
    loo$collection_label,
    FUN = function(x) rank(x, ties.method = "first", na.last = "keep")
  )
  rows <- lapply(split(loo, paste(loo$collection_label, loo$pathway, sep = "\r")), function(df) {
    data.frame(
      collection = df$collection[[1L]],
      collection_label = df$collection_label[[1L]],
      pathway = df$pathway[[1L]],
      pathway_label = df$pathway_label[[1L]],
      n_leave_one_out_runs = length(unique(df$left_out_sample_id)),
      mean_leave_one_out_NES = mean(df$NES, na.rm = TRUE),
      sd_leave_one_out_NES = stats::sd(df$NES, na.rm = TRUE),
      prop_top20_leave_one_out = mean(df$loo_rank <= 20, na.rm = TRUE),
      prop_sign_concordant_leave_one_out = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- merge(out, primary, by = c("collection", "collection_label", "pathway", "pathway_label"), all.x = TRUE, sort = FALSE)
  out$prop_sign_concordant_leave_one_out <- vapply(seq_len(nrow(out)), function(i) {
    local <- loo[loo$collection == out$collection[[i]] & loo$pathway == out$pathway[[i]], , drop = FALSE]
    mean(sign(local$NES) == sign(out$primary_NES[[i]]), na.rm = TRUE)
  }, numeric(1L))
  out[order(out$collection_label, out$primary_padj, -abs(out$primary_NES), out$pathway), , drop = FALSE]
}

pst_model_dirs <- function(model_root) {
  list(
    manifest = pst_ensure_dir(file.path(model_root, "00_manifest")),
    qc = pst_ensure_dir(file.path(model_root, "01_qc")),
    gene_models = pst_ensure_dir(file.path(model_root, "03_gene_models")),
    gsea = pst_ensure_dir(file.path(model_root, "04_gsea")),
    figures = pst_ensure_dir(file.path(model_root, "05_figures")),
    sensitivity = pst_ensure_dir(file.path(model_root, "06_sensitivity"))
  )
}

pst_model_parameters_table <- function(model_spec) {
  data.frame(
    parameter = c("model_id", "model_label", "covariate_mode", "covariate_terms", "etp_method", "etp_threshold", "is_primary"),
    value = c(
      model_spec$model_id %||% NA_character_,
      model_spec$model_label %||% NA_character_,
      model_spec$covariate_mode %||% NA_character_,
      paste(model_spec$covariate_terms %||% character(0L), collapse = ";"),
      model_spec$etp_method %||% NA_character_,
      as.character(model_spec$etp_threshold %||% NA_real_),
      as.character(isTRUE(model_spec$is_primary))
    ),
    stringsAsFactors = FALSE
  )
}

pst_run_single_model <- function(model_spec, output_root, pb, counts, coverage, coverage_check, cfg, gene_sets, args) {
  model_root <- pst_ensure_dir(file.path(output_root, model_spec$model_id))
  dirs <- pst_model_dirs(model_root)
  pst_write_csv(pst_model_parameters_table(model_spec), file.path(dirs$manifest, "model_parameters.csv"))
  min_primary_mice <- coverage_check$n_contributing_mice[[1L]]

  model_fit <- pst_fit_pseudotime_model(
    pb$counts,
    pb$metadata,
    args$spline_df,
    min_primary_mice,
    include_dose = TRUE,
    model_spec = model_spec
  )
  design_audit <- pst_design_audit(model_fit)
  pst_write_csv(design_audit, file.path(dirs$qc, "model_design_rank_audit.csv"))

  primary_contrast <- pst_contrast_vector(model_fit, cfg, "primary_adjacent_state", args$grid_size)
  secondary_contrast <- pst_contrast_vector(model_fit, cfg, "full_trajectory_state_specificity", args$grid_size)
  primary_genes <- pst_contrast_table(model_fit, primary_contrast, "primary_adjacent_state")
  secondary_genes <- pst_contrast_table(model_fit, secondary_contrast, "full_trajectory_state_specificity")
  primary_genes$model_id <- model_spec$model_id
  secondary_genes$model_id <- model_spec$model_id
  symbol_resolution <- pst_resolve_gene_symbols(primary_genes)
  symbol_resolution$model_id <- model_spec$model_id
  pst_write_csv(primary_genes, file.path(dirs$gene_models, "gene_primary_adjacent_state_contrast.csv"))
  pst_write_csv(secondary_genes, file.path(dirs$gene_models, "gene_full_trajectory_state_specificity.csv"))
  pst_write_csv(symbol_resolution, file.path(dirs$gene_models, "gene_symbol_resolution.csv"))

  primary_stats <- pst_ranked_stats(primary_genes, symbol_resolution)
  secondary_stats <- pst_ranked_stats(secondary_genes)
  primary_gsea <- pst_run_all_gsea(primary_stats, gene_sets, "primary_adjacent_state", args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  secondary_gsea <- pst_run_all_gsea(secondary_stats, gene_sets, "full_trajectory_state_specificity", args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  primary_gsea$model_id <- model_spec$model_id
  secondary_gsea$model_id <- model_spec$model_id
  pst_write_gsea_tables(primary_gsea, secondary_gsea, dirs$gsea)
  leading_edge <- pst_leading_edge_table(rbind(primary_gsea, secondary_gsea))
  if (nrow(leading_edge) > 0L) leading_edge$model_id <- model_spec$model_id
  pst_write_csv(leading_edge, file.path(dirs$gsea, "all_collections_leading_edge_genes.csv"))
  if (nrow(leading_edge) > 0L) {
    for (collection_label in unique(leading_edge$collection_label)) {
      pst_write_csv(
        leading_edge[leading_edge$collection_label == collection_label, , drop = FALSE],
        file.path(dirs$gsea, paste0(collection_label, "_leading_edge_genes.csv"))
      )
    }
  }

  fitted_grid <- pst_fitted_grid(model_fit, args$grid_size)
  activity <- pst_pathway_activity(fitted_grid, primary_gsea)
  if (nrow(activity) > 0L) activity$model_id <- model_spec$model_id
  pst_write_csv(activity, file.path(dirs$gsea, "pathway_activity_over_pseudotime.csv"))

  density_data <- pst_density_difference_data(pb$cell_metadata, cfg, args$grid_size)
  pst_write_csv(density_data$sample_density, file.path(dirs$figures, "pseudotime_density_with_frozen_state_sample_density_plot_data.csv"))
  pst_write_csv(density_data$mean_density, file.path(dirs$figures, "pseudotime_density_with_frozen_state_mean_density_plot_data.csv"))
  pst_write_csv(density_data$difference, file.path(dirs$figures, "pseudotime_density_with_frozen_state_difference_plot_data.csv"))
  pst_plot_density(density_data, pst_intervals_table(cfg), file.path(dirs$figures, "pseudotime_density_with_frozen_state.pdf"))

  for (collection_label in unique(primary_gsea$collection_label)) {
    local <- primary_gsea[primary_gsea$collection_label == collection_label, , drop = FALSE]
    pst_write_csv(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea_plot_data.csv")))
    pst_plot_gsea_bar(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea.pdf")), paste0(model_spec$model_id, ": ", collection_label, " GSEA"))
  }
  pst_write_csv(activity, file.path(dirs$figures, "primary_state_pathway_activity_heatmap_plot_data.csv"))
  pst_plot_activity_heatmap(activity, file.path(dirs$figures, "primary_state_pathway_activity_heatmap.pdf"))
  pst_write_csv(activity, file.path(dirs$figures, "primary_state_top_pathway_curves_plot_data.csv"))
  pst_plot_activity_curves(activity, file.path(dirs$figures, "primary_state_top_pathway_curves.pdf"))
  gene_curves <- pst_leading_gene_curve_data(fitted_grid, leading_edge, primary_gsea)
  if (nrow(gene_curves) > 0L) gene_curves$model_id <- model_spec$model_id
  pst_write_csv(gene_curves, file.path(dirs$figures, "primary_state_leading_edge_gene_curves_plot_data.csv"))
  pst_plot_gene_curves(gene_curves, file.path(dirs$figures, "primary_state_leading_edge_gene_curves.pdf"))
  pst_write_csv(coverage, file.path(dirs$figures, "sample_region_coverage_plot_data.csv"))
  pst_plot_coverage(coverage, file.path(dirs$figures, "sample_region_coverage.pdf"))

  robustness <- data.frame()
  if (isTRUE(args$run_sensitivity)) {
    robustness <- pst_run_sensitivity(pb$counts, pb$metadata, cfg, gene_sets, primary_gsea, args, pb$cell_metadata, counts, model_spec = model_spec)
    if (nrow(robustness) > 0L) robustness$model_id <- model_spec$model_id
    pst_write_csv(robustness, file.path(dirs$sensitivity, "pathway_robustness_summary.csv"))
    pst_write_csv(robustness, file.path(dirs$figures, "pathway_robustness_summary_plot_data.csv"))
    pst_plot_robustness(robustness, file.path(dirs$figures, "pathway_robustness_summary.pdf"))
  }
  loo <- data.frame()
  loo_stability <- data.frame()
  if (isTRUE(args$run_leave_one_out)) {
    loo <- pst_run_leave_one_out(pb$counts, pb$metadata, cfg, gene_sets, args, model_spec = model_spec)
    if (nrow(loo) > 0L) loo$model_id <- model_spec$model_id
    pst_write_csv(loo, file.path(dirs$sensitivity, "leave_one_mouse_out_all_collections_gsea.csv"))
    if (nrow(loo) > 0L) {
      for (collection_label in unique(loo$collection_label)) {
        pst_write_csv(loo[loo$collection_label == collection_label, , drop = FALSE], file.path(dirs$sensitivity, paste0("leave_one_mouse_out_", collection_label, "_gsea.csv")))
      }
      loo_stability <- pst_leave_one_out_stability(loo, primary_gsea)
      if (nrow(loo_stability) > 0L) loo_stability$model_id <- model_spec$model_id
      pst_write_csv(loo_stability, file.path(dirs$sensitivity, "leave_one_mouse_out_pathway_stability.csv"))
    }
  }

  list(
    model_id = model_spec$model_id,
    model_spec = model_spec,
    model_root = model_root,
    design_audit = design_audit,
    primary_genes = primary_genes,
    secondary_genes = secondary_genes,
    primary_gsea = primary_gsea,
    secondary_gsea = secondary_gsea,
    robustness = robustness,
    loo = loo,
    loo_stability = loo_stability
  )
}

pst_safe_cor <- function(x, y, method = "pearson") {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L) return(NA_real_)
  suppressWarnings(stats::cor(x[keep], y[keep], method = method))
}

pst_gene_t_correlations <- function(model_results, reference_id) {
  ref <- model_results[[reference_id]]$primary_genes[, c("gene", "gene_symbol", "t_statistic"), drop = FALSE]
  names(ref)[names(ref) == "t_statistic"] <- "reference_t_statistic"
  rows <- lapply(names(model_results), function(model_id) {
    df <- model_results[[model_id]]$primary_genes[, c("gene", "t_statistic"), drop = FALSE]
    names(df)[names(df) == "t_statistic"] <- "model_t_statistic"
    merged <- merge(ref, df, by = "gene", all = FALSE, sort = FALSE)
    data.frame(
      reference_model_id = reference_id,
      model_id = model_id,
      n_genes = nrow(merged),
      pearson_r = pst_safe_cor(merged$reference_t_statistic, merged$model_t_statistic, "pearson"),
      spearman_rho = pst_safe_cor(merged$reference_t_statistic, merged$model_t_statistic, "spearman"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pst_pathway_model_long <- function(model_results) {
  rows <- lapply(names(model_results), function(model_id) {
    df <- model_results[[model_id]]$primary_gsea
    keep <- intersect(
      c("model_id", "collection", "collection_label", "pathway", "pathway_label", "NES", "pval", "padj", "size", "direction"),
      names(df)
    )
    df[, keep, drop = FALSE]
  })
  do.call(rbind, rows)
}

pst_pathway_model_wide <- function(long) {
  keys <- c("collection", "collection_label", "pathway", "pathway_label")
  models <- unique(long$model_id)
  out <- unique(long[, keys, drop = FALSE])
  for (model_id in models) {
    local <- long[long$model_id == model_id, c(keys, "NES", "padj"), drop = FALSE]
    names(local)[names(local) == "NES"] <- paste0("NES__", model_id)
    names(local)[names(local) == "padj"] <- paste0("padj__", model_id)
    out <- merge(out, local, by = keys, all.x = TRUE, sort = FALSE)
  }
  out
}

pst_gsea_nes_correlations <- function(pathway_long, reference_id) {
  ref <- pathway_long[pathway_long$model_id == reference_id, c("collection", "collection_label", "pathway", "NES"), drop = FALSE]
  names(ref)[names(ref) == "NES"] <- "reference_NES"
  rows <- list()
  for (model_id in unique(pathway_long$model_id)) {
    local <- pathway_long[pathway_long$model_id == model_id, c("collection", "collection_label", "pathway", "NES"), drop = FALSE]
    names(local)[names(local) == "NES"] <- "model_NES"
    merged <- merge(ref, local, by = c("collection", "collection_label", "pathway"), all = FALSE, sort = FALSE)
    for (collection_label in unique(merged$collection_label)) {
      sub <- merged[merged$collection_label == collection_label, , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(
        reference_model_id = reference_id,
        model_id = model_id,
        collection_label = collection_label,
        n_pathways = nrow(sub),
        pearson_r = pst_safe_cor(sub$reference_NES, sub$model_NES, "pearson"),
        spearman_rho = pst_safe_cor(sub$reference_NES, sub$model_NES, "spearman"),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

pst_gsea_significant_overlap <- function(pathway_long, reference_id, alpha = 0.05) {
  ref <- pathway_long[pathway_long$model_id == reference_id, , drop = FALSE]
  rows <- list()
  for (model_id in unique(pathway_long$model_id)) {
    local <- pathway_long[pathway_long$model_id == model_id, , drop = FALSE]
    for (collection_label in unique(pathway_long$collection_label)) {
      ref_set <- ref$pathway[ref$collection_label == collection_label & ref$padj < alpha]
      model_set <- local$pathway[local$collection_label == collection_label & local$padj < alpha]
      overlap <- intersect(ref_set, model_set)
      union <- union(ref_set, model_set)
      rows[[length(rows) + 1L]] <- data.frame(
        reference_model_id = reference_id,
        model_id = model_id,
        collection_label = collection_label,
        alpha = alpha,
        reference_significant_n = length(ref_set),
        model_significant_n = length(model_set),
        overlap_n = length(overlap),
        jaccard = if (length(union) == 0L) NA_real_ else length(overlap) / length(union),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

pst_pathway_rank_shift <- function(pathway_long, reference_id) {
  pathway_long$rank_within_model <- ave(
    pathway_long$padj,
    pathway_long$model_id,
    pathway_long$collection_label,
    FUN = function(x) rank(x, ties.method = "first", na.last = "keep")
  )
  ref <- pathway_long[pathway_long$model_id == reference_id, c("collection", "collection_label", "pathway", "rank_within_model", "NES", "padj"), drop = FALSE]
  names(ref)[names(ref) == "rank_within_model"] <- "reference_rank"
  names(ref)[names(ref) == "NES"] <- "reference_NES"
  names(ref)[names(ref) == "padj"] <- "reference_padj"
  rows <- lapply(setdiff(unique(pathway_long$model_id), character(0L)), function(model_id) {
    local <- pathway_long[pathway_long$model_id == model_id, c("collection", "collection_label", "pathway", "pathway_label", "rank_within_model", "NES", "padj"), drop = FALSE]
    names(local)[names(local) == "rank_within_model"] <- "model_rank"
    names(local)[names(local) == "NES"] <- "model_NES"
    names(local)[names(local) == "padj"] <- "model_padj"
    merged <- merge(ref, local, by = c("collection", "collection_label", "pathway"), all = FALSE, sort = FALSE)
    merged$model_id <- model_id
    merged$reference_model_id <- reference_id
    merged$rank_shift <- merged$model_rank - merged$reference_rank
    merged$NES_delta <- merged$model_NES - merged$reference_NES
    merged
  })
  out <- do.call(rbind, rows)
  out[order(out$model_id, out$collection_label, -abs(out$rank_shift), out$pathway), , drop = FALSE]
}

pst_plot_nes_correlation <- function(pathway_long, reference_id, path) {
  ref <- pathway_long[pathway_long$model_id == reference_id, c("collection", "collection_label", "pathway", "NES"), drop = FALSE]
  names(ref)[names(ref) == "NES"] <- "reference_NES"
  alternative_ids <- setdiff(unique(pathway_long$model_id), reference_id)
  if (length(alternative_ids) == 0L) return(invisible(NULL))
  rows <- lapply(alternative_ids, function(model_id) {
    local <- pathway_long[pathway_long$model_id == model_id, c("collection", "collection_label", "pathway", "NES"), drop = FALSE]
    names(local)[names(local) == "NES"] <- "model_NES"
    merged <- merge(ref, local, by = c("collection", "collection_label", "pathway"), all = FALSE, sort = FALSE)
    merged$model_id <- model_id
    merged
  })
  df <- do.call(rbind, rows)
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = reference_NES, y = model_NES)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.45, size = 0.9) +
    ggplot2::facet_grid(model_id ~ collection_label, scales = "free") +
    ggplot2::labs(title = "Primary-state GSEA NES correlation versus primary model", x = paste0(reference_id, " NES"), y = "Alternative model NES") +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(path, p, width = 12, height = 8)
  invisible(p)
}

pst_plot_rank_shift <- function(rank_shift, path) {
  df <- rank_shift[rank_shift$model_id != rank_shift$reference_model_id, , drop = FALSE]
  if (nrow(df) == 0L) return(invisible(NULL))
  df <- df[order(-abs(df$rank_shift)), , drop = FALSE]
  df <- head(df, 80L)
  df$pathway_label <- factor(df$pathway_label, levels = rev(unique(df$pathway_label)))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = rank_shift, y = pathway_label, color = collection_label)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~ model_id, scales = "free_y") +
    ggplot2::labs(title = "Largest pathway rank shifts versus primary model", x = "Rank shift", y = "Pathway", color = "Collection") +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(path, p, width = 13, height = 9)
  invisible(p)
}

pst_run_model_comparison <- function(model_results, output_root) {
  comparison_root <- pst_ensure_dir(file.path(output_root, "model_comparison"))
  figures_dir <- pst_ensure_dir(file.path(comparison_root, "figures"))
  reference_id <- if ("primary_initial_ploidy" %in% names(model_results)) "primary_initial_ploidy" else names(model_results)[[1L]]
  design_audits <- do.call(rbind, lapply(model_results, `[[`, "design_audit"))
  gene_cor <- pst_gene_t_correlations(model_results, reference_id)
  pathway_long <- pst_pathway_model_long(model_results)
  pathway_wide <- pst_pathway_model_wide(pathway_long)
  nes_cor <- pst_gsea_nes_correlations(pathway_long, reference_id)
  sig_overlap <- pst_gsea_significant_overlap(pathway_long, reference_id)
  rank_shift <- pst_pathway_rank_shift(pathway_long, reference_id)
  pst_write_csv(design_audits, file.path(comparison_root, "model_design_rank_audit_all_models.csv"))
  pst_write_csv(gene_cor, file.path(comparison_root, "gene_t_stat_correlations.csv"))
  pst_write_csv(pathway_long, file.path(comparison_root, "pathway_model_comparison_long.csv"))
  pst_write_csv(pathway_wide, file.path(comparison_root, "pathway_model_comparison_wide.csv"))
  pst_write_csv(nes_cor, file.path(comparison_root, "gsea_NES_correlations.csv"))
  pst_write_csv(sig_overlap, file.path(comparison_root, "gsea_significant_overlap.csv"))
  pst_write_csv(rank_shift, file.path(comparison_root, "pathway_rank_shift.csv"))
  pst_plot_nes_correlation(pathway_long, reference_id, file.path(figures_dir, "gsea_NES_correlation_primary_vs_ETP.pdf"))
  pst_plot_rank_shift(rank_shift, file.path(figures_dir, "top_pathway_rank_shift_primary_vs_ETP.pdf"))
  lines <- c(
    "# 04i model comparison",
    "",
    paste0("Reference model: `", reference_id, "`."),
    "",
    "Tables compare the primary adjacent-state gene statistics and GSEA results across covariate specifications.",
    "",
    "- `gene_t_stat_correlations.csv`: gene-level t-statistic concordance.",
    "- `gsea_NES_correlations.csv`: pathway NES concordance by collection.",
    "- `gsea_significant_overlap.csv`: FDR-significant pathway overlap.",
    "- `pathway_rank_shift.csv`: pathway rank changes versus the reference model.",
    "- `model_design_rank_audit_all_models.csv`: design-rank and dropped-column audit."
  )
  pst_write_lines(lines, file.path(comparison_root, "model_comparison_summary.md"))
  invisible(list(reference_id = reference_id, design_audits = design_audits))
}

pst_materialize_snapshot <- function(output_root, snapshot_root) {
  if (dir.exists(snapshot_root)) unlink(snapshot_root, recursive = TRUE, force = TRUE)
  pst_ensure_dir(snapshot_root)
  figure_files <- list.files(
    output_root,
    pattern = "\\.(pdf|png|svg|tif|tiff)$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(figure_files) == 0L) return(invisible(snapshot_root))
  relative <- substring(figure_files, nchar(output_root) + 2L)
  for (i in seq_along(figure_files)) {
    dst <- file.path(snapshot_root, relative[[i]])
    pst_ensure_dir(dirname(dst))
    file.copy(figure_files[[i]], dst, overwrite = TRUE)
  }
  ds_store <- list.files(snapshot_root, pattern = "^\\.DS_Store$", recursive = TRUE, full.names = TRUE, all.files = TRUE)
  if (length(ds_store) > 0L) unlink(ds_store, force = TRUE)
  invisible(snapshot_root)
}

pst_write_readme_legacy <- function(output_root) {
  lines <- c(
    "# 04i Pseudotime State Pathways",
    "",
    "This workflow annotates the frozen CellCycle pseudotime state in which gemcitabine-treated tumor cells accumulate.",
    "",
    "The treatment-versus-control density comparison is used only to define the frozen pseudotime intervals. Expression modeling uses a common pseudotime smooth and does not test treatment-by-pseudotime expression effects.",
    "",
    "Primary interpretation: positive NES pathways characterize the accumulated pseudotime state relative to adjacent states. They must not be described as directly treatment-upregulated pathways or as evidence that treatment accelerates or slows a pathway.",
    "",
    "Gene-set collections run by default: Hallmark H, Reactome C2:CP:REACTOME, and GO Biological Process C5:GO:BP. FDR is interpreted within each collection and ranking.",
    "",
    "Main directories:",
    "",
    "- `00_manifest`: parameters, checksums, package versions, and frozen interval definitions.",
    "- `01_qc`: cell ID joins, expression source audit, and region coverage.",
    "- `02_pseudobulk`: sample-by-pseudotime-bin metadata.",
    "- `03_gene_models`: primary and secondary gene ranking tables.",
    "- `04_gsea`: complete GSEA results and pathway activity tables.",
    "- `05_figures`: PDFs and their plotting-data CSV files.",
    "- `06_sensitivity`: sensitivity and leave-one-mouse-out pathway results."
  )
  pst_write_lines(lines, file.path(output_root, "README.md"))
}

pst_run_workflow_legacy <- function(args, repo_root, script_dir) {
  pst_check_packages()
  set.seed(args$seed)
  output_root <- pst_clean_dir(args$output_root, overwrite = args$overwrite)
  dirs <- list(
    manifest = pst_ensure_dir(file.path(output_root, "00_manifest")),
    qc = pst_ensure_dir(file.path(output_root, "01_qc")),
    pseudobulk = pst_ensure_dir(file.path(output_root, "02_pseudobulk")),
    gene_models = pst_ensure_dir(file.path(output_root, "03_gene_models")),
    gsea = pst_ensure_dir(file.path(output_root, "04_gsea")),
    figures = pst_ensure_dir(file.path(output_root, "05_figures")),
    sensitivity = pst_ensure_dir(file.path(output_root, "06_sensitivity"))
  )

  cfg <- pst_read_config(args$config)
  seurat_rds <- pst_resolve_seurat_rds(args$seurat_rds, args$results_root)
  params <- data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE)
  pst_write_csv(params, file.path(dirs$manifest, "analysis_parameters.csv"))
  pst_write_csv(pst_intervals_table(cfg), file.path(dirs$manifest, "frozen_interval_definition.csv"))
  pst_write_csv(pst_package_versions(), file.path(dirs$manifest, "package_versions.csv"))

  input_checksums <- data.frame(
    input = c("cell_metadata", "seurat_rds", "config"),
    path = c(args$cell_metadata, seurat_rds, args$config),
    sha256 = c(pst_file_checksum(args$cell_metadata), pst_file_checksum(seurat_rds), pst_file_checksum(args$config)),
    stringsAsFactors = FALSE
  )
  pst_write_csv(input_checksums, file.path(dirs$manifest, "input_checksums.csv"))

  meta <- pst_read_cell_metadata(args$cell_metadata)
  metadata_audit <- pst_audit_metadata(meta)
  pst_write_csv(metadata_audit$duplicate_cells, file.path(dirs$qc, "duplicate_cell_ids.csv"))
  pst_write_csv(metadata_audit$inconsistent_cells, file.path(dirs$qc, "cell_assignment_consistency_audit.csv"))
  pst_check_pseudotime(meta)

  counts <- pst_load_counts(seurat_rds, args$assay, args$counts_layer)
  expression_audit <- pst_audit_counts(counts)
  expression_audit$seurat_rds <- seurat_rds
  expression_audit$assay <- args$assay
  expression_audit$counts_layer <- args$counts_layer
  pst_write_csv(expression_audit, file.path(dirs$qc, "expression_source_audit.csv"))
  if (!isTRUE(expression_audit$non_negative) || !isTRUE(expression_audit$integer_like)) {
    stop("Expression matrix failed raw-count audit. See 01_qc/expression_source_audit.csv.", call. = FALSE)
  }

  matched <- pst_match_cells(meta, counts, args$min_match_rate, dirs$qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- pst_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts

  pb <- pst_construct_pseudobulk(meta, counts, args$n_pseudotime_bins, args$min_cells_per_sample_bin)
  pst_write_csv(pb$metadata, file.path(dirs$pseudobulk, "sample_bin_metadata.csv"))
  coverage <- pst_region_coverage(pb$cell_metadata, pb$metadata, cfg)
  pst_write_csv(coverage, file.path(dirs$qc, "sample_region_coverage.csv"))
  coverage_check <- pst_check_primary_coverage(coverage)
  pst_write_csv(coverage_check, file.path(dirs$qc, "primary_coverage_check.csv"))
  if (!isTRUE(coverage_check$passes)) {
    stop("Primary interval coverage check failed. See 01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  min_primary_mice <- coverage_check$n_contributing_mice[[1L]]
  model_fit <- pst_fit_pseudotime_model(pb$counts, pb$metadata, args$spline_df, min_primary_mice, include_dose = TRUE)
  primary_contrast <- pst_contrast_vector(model_fit, cfg, "primary_adjacent_state", args$grid_size)
  secondary_contrast <- pst_contrast_vector(model_fit, cfg, "full_trajectory_state_specificity", args$grid_size)
  primary_genes <- pst_contrast_table(model_fit, primary_contrast, "primary_adjacent_state")
  secondary_genes <- pst_contrast_table(model_fit, secondary_contrast, "full_trajectory_state_specificity")
  symbol_resolution <- pst_resolve_gene_symbols(primary_genes)
  pst_write_csv(primary_genes, file.path(dirs$gene_models, "gene_primary_adjacent_state_contrast.csv"))
  pst_write_csv(secondary_genes, file.path(dirs$gene_models, "gene_full_trajectory_state_specificity.csv"))
  pst_write_csv(symbol_resolution, file.path(dirs$gene_models, "gene_symbol_resolution.csv"))

  gene_sets <- pst_fetch_gene_sets(args$gene_set_collections)
  primary_stats <- pst_ranked_stats(primary_genes, symbol_resolution)
  secondary_stats <- pst_ranked_stats(secondary_genes)
  primary_gsea <- pst_run_all_gsea(primary_stats, gene_sets, "primary_adjacent_state", args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  secondary_gsea <- pst_run_all_gsea(secondary_stats, gene_sets, "full_trajectory_state_specificity", args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  pst_write_gsea_tables(primary_gsea, secondary_gsea, dirs$gsea)
  leading_edge <- pst_leading_edge_table(rbind(primary_gsea, secondary_gsea))
  pst_write_csv(leading_edge, file.path(dirs$gsea, "all_collections_leading_edge_genes.csv"))
  if (nrow(leading_edge) > 0L) {
    for (collection_label in unique(leading_edge$collection_label)) {
      pst_write_csv(
        leading_edge[leading_edge$collection_label == collection_label, , drop = FALSE],
        file.path(dirs$gsea, paste0(collection_label, "_leading_edge_genes.csv"))
      )
    }
  }

  fitted_grid <- pst_fitted_grid(model_fit, args$grid_size)
  activity <- pst_pathway_activity(fitted_grid, primary_gsea)
  pst_write_csv(activity, file.path(dirs$gsea, "pathway_activity_over_pseudotime.csv"))

  density_data <- pst_density_difference_data(pb$cell_metadata, cfg, args$grid_size)
  pst_write_csv(density_data$sample_density, file.path(dirs$figures, "pseudotime_density_with_frozen_state_sample_density_plot_data.csv"))
  pst_write_csv(density_data$mean_density, file.path(dirs$figures, "pseudotime_density_with_frozen_state_mean_density_plot_data.csv"))
  pst_write_csv(density_data$difference, file.path(dirs$figures, "pseudotime_density_with_frozen_state_difference_plot_data.csv"))
  pst_plot_density(density_data, pst_intervals_table(cfg), file.path(dirs$figures, "pseudotime_density_with_frozen_state.pdf"))

  for (collection_label in unique(primary_gsea$collection_label)) {
    local <- primary_gsea[primary_gsea$collection_label == collection_label, , drop = FALSE]
    pst_write_csv(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea_plot_data.csv")))
    pst_plot_gsea_bar(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea.pdf")), paste0("Primary accumulated state: ", collection_label, " GSEA"))
  }
  pst_write_csv(activity, file.path(dirs$figures, "primary_state_pathway_activity_heatmap_plot_data.csv"))
  pst_plot_activity_heatmap(activity, file.path(dirs$figures, "primary_state_pathway_activity_heatmap.pdf"))
  pst_write_csv(activity, file.path(dirs$figures, "primary_state_top_pathway_curves_plot_data.csv"))
  pst_plot_activity_curves(activity, file.path(dirs$figures, "primary_state_top_pathway_curves.pdf"))
  gene_curves <- pst_leading_gene_curve_data(fitted_grid, leading_edge, primary_gsea)
  pst_write_csv(gene_curves, file.path(dirs$figures, "primary_state_leading_edge_gene_curves_plot_data.csv"))
  pst_plot_gene_curves(gene_curves, file.path(dirs$figures, "primary_state_leading_edge_gene_curves.pdf"))
  pst_write_csv(coverage, file.path(dirs$figures, "sample_region_coverage_plot_data.csv"))
  pst_plot_coverage(coverage, file.path(dirs$figures, "sample_region_coverage.pdf"))

  robustness <- data.frame()
  if (isTRUE(args$run_sensitivity)) {
    robustness <- pst_run_sensitivity(pb$counts, pb$metadata, cfg, gene_sets, primary_gsea, args, pb$cell_metadata, counts)
    pst_write_csv(robustness, file.path(dirs$sensitivity, "pathway_robustness_summary.csv"))
    pst_write_csv(robustness, file.path(dirs$figures, "pathway_robustness_summary_plot_data.csv"))
    pst_plot_robustness(robustness, file.path(dirs$figures, "pathway_robustness_summary.pdf"))
  }
  if (isTRUE(args$run_leave_one_out)) {
    loo <- pst_run_leave_one_out(pb$counts, pb$metadata, cfg, gene_sets, args)
    pst_write_csv(loo, file.path(dirs$sensitivity, "leave_one_mouse_out_all_collections_gsea.csv"))
    if (nrow(loo) > 0L) {
      for (collection_label in unique(loo$collection_label)) {
        pst_write_csv(loo[loo$collection_label == collection_label, , drop = FALSE], file.path(dirs$sensitivity, paste0("leave_one_mouse_out_", collection_label, "_gsea.csv")))
      }
      loo_stability <- pst_leave_one_out_stability(loo, primary_gsea)
      pst_write_csv(loo_stability, file.path(dirs$sensitivity, "leave_one_mouse_out_pathway_stability.csv"))
    }
  }

  pst_write_readme(output_root)
  if (isTRUE(args$make_snapshot)) pst_materialize_snapshot(output_root, args$snapshot_root)
  message("04i workflow complete: ", output_root)
  invisible(output_root)
}

pst_write_readme <- function(output_root) {
  lines <- c(
    "# 04i Pseudotime State Pathways",
    "",
    "This workflow annotates the frozen CellCycle pseudotime state in which gemcitabine-treated tumor cells accumulate.",
    "",
    "The treatment-versus-control density comparison is used only to define the frozen pseudotime intervals. Expression modeling uses a common pseudotime smooth and does not test treatment-by-pseudotime expression effects.",
    "",
    "Primary interpretation: positive NES pathways characterize the accumulated state relative to adjacent states. They must not be described as directly treatment-upregulated pathways or as evidence that treatment accelerates or slows a pathway.",
    "",
    "Gene-set collections run by default: Hallmark H, Reactome C2:CP:REACTOME, and GO Biological Process C5:GO:BP. FDR is interpreted within each collection and ranking.",
    "",
    "Main directories:",
    "",
    "- `00_manifest`: shared parameters, checksums, package versions, and frozen interval definitions.",
    "- `01_qc`: shared cell ID joins, expression source audit, region coverage, and ETP covariate audits.",
    "- `02_pseudobulk`: shared sample-by-pseudotime-bin metadata.",
    "- `primary_initial_ploidy`: primary model specified by the analysis plan.",
    "- `etp_continuous_adjusted`: continuous mean-ETP covariate sensitivity model when selected.",
    "- `ETP_*`: thresholded ETP-group covariate sensitivity models when selected.",
    "- `model_comparison`: gene, pathway, and design-rank comparison across selected models.",
    "",
    "Snapshot policy:",
    "",
    "The repository snapshot copies figure files only. Plot-data tables and statistical CSV files remain in the results directory."
  )
  pst_write_lines(lines, file.path(output_root, "README.md"))
}

pst_run_workflow <- function(args, repo_root, script_dir) {
  pst_check_packages()
  set.seed(args$seed)
  output_root <- pst_clean_dir(args$output_root, overwrite = args$overwrite)
  etp_specs <- pst_resolve_etp_specs(args$etp_threshold)
  model_specs <- pst_model_specs(args$covariate_mode, etp_specs, include_full_etp_models = args$include_full_etp_models)
  dirs <- list(
    manifest = pst_ensure_dir(file.path(output_root, "00_manifest")),
    qc = pst_ensure_dir(file.path(output_root, "01_qc")),
    pseudobulk = pst_ensure_dir(file.path(output_root, "02_pseudobulk"))
  )

  cfg <- pst_read_config(args$config)
  seurat_rds <- pst_resolve_seurat_rds(args$seurat_rds, args$results_root)
  params <- data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE)
  pst_write_csv(params, file.path(dirs$manifest, "analysis_parameters.csv"))
  pst_write_csv(pst_intervals_table(cfg), file.path(dirs$manifest, "frozen_interval_definition.csv"))
  pst_write_csv(pst_package_versions(), file.path(dirs$manifest, "package_versions.csv"))
  selected_model_table <- do.call(rbind, lapply(model_specs, function(x) {
    data.frame(
      model_id = x$model_id,
      model_label = x$model_label,
      covariate_mode = x$covariate_mode,
      covariate_terms = paste(x$covariate_terms, collapse = ";"),
      etp_method = x$etp_method %||% NA_character_,
      etp_threshold = x$etp_threshold %||% NA_real_,
      is_primary = isTRUE(x$is_primary),
      stringsAsFactors = FALSE
    )
  }))
  pst_write_csv(selected_model_table, file.path(dirs$manifest, "selected_model_specifications.csv"))

  input_checksums <- data.frame(
    input = c("cell_metadata", "noncell_metadata", "seurat_rds", "config"),
    path = c(args$cell_metadata, args$noncell_metadata, seurat_rds, args$config),
    sha256 = c(
      pst_file_checksum(args$cell_metadata),
      pst_file_checksum(args$noncell_metadata),
      pst_file_checksum(seurat_rds),
      pst_file_checksum(args$config)
    ),
    stringsAsFactors = FALSE
  )
  pst_write_csv(input_checksums, file.path(dirs$manifest, "input_checksums.csv"))

  meta <- pst_read_cell_metadata(args$cell_metadata)
  etp_assignments <- pst_build_etp_assignments(args$cell_metadata, args$noncell_metadata, etp_specs)
  pst_write_etp_audits(etp_assignments, etp_specs, dirs$qc)
  meta <- pst_attach_etp_assignments(meta, etp_assignments)
  metadata_audit <- pst_audit_metadata(meta)
  pst_write_csv(metadata_audit$duplicate_cells, file.path(dirs$qc, "duplicate_cell_ids.csv"))
  pst_write_csv(metadata_audit$inconsistent_cells, file.path(dirs$qc, "cell_assignment_consistency_audit.csv"))
  pst_check_pseudotime(meta)

  counts <- pst_load_counts(seurat_rds, args$assay, args$counts_layer)
  expression_audit <- pst_audit_counts(counts)
  expression_audit$seurat_rds <- seurat_rds
  expression_audit$assay <- args$assay
  expression_audit$counts_layer <- args$counts_layer
  pst_write_csv(expression_audit, file.path(dirs$qc, "expression_source_audit.csv"))
  if (!isTRUE(expression_audit$non_negative) || !isTRUE(expression_audit$integer_like)) {
    stop("Expression matrix failed raw-count audit. See 01_qc/expression_source_audit.csv.", call. = FALSE)
  }

  matched <- pst_match_cells(meta, counts, args$min_match_rate, dirs$qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- pst_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts

  pb <- pst_construct_pseudobulk(meta, counts, args$n_pseudotime_bins, args$min_cells_per_sample_bin)
  pst_write_csv(pb$metadata, file.path(dirs$pseudobulk, "sample_bin_metadata.csv"))
  coverage <- pst_region_coverage(pb$cell_metadata, pb$metadata, cfg)
  pst_write_csv(coverage, file.path(dirs$qc, "sample_region_coverage.csv"))
  coverage_check <- pst_check_primary_coverage(coverage)
  pst_write_csv(coverage_check, file.path(dirs$qc, "primary_coverage_check.csv"))
  if (!isTRUE(coverage_check$passes)) {
    stop("Primary interval coverage check failed. See 01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  gene_sets <- pst_fetch_gene_sets(args$gene_set_collections)
  model_results <- list()
  for (model_id in names(model_specs)) {
    message("Running 04i model: ", model_id)
    model_results[[model_id]] <- pst_run_single_model(
      model_specs[[model_id]],
      output_root,
      pb,
      counts,
      coverage,
      coverage_check,
      cfg,
      gene_sets,
      args
    )
  }
  pst_run_model_comparison(model_results, output_root)

  pst_write_readme(output_root)
  if (isTRUE(args$make_snapshot)) pst_materialize_snapshot(output_root, args$snapshot_root)
  message("04i workflow complete: ", output_root)
  invisible(output_root)
}

pst_model_specifications_table <- function(model_specs) {
  do.call(rbind, lapply(model_specs, function(x) {
    data.frame(
      model_id = x$model_id,
      model_label = x$model_label,
      covariate_mode = x$covariate_mode,
      covariate_terms = paste(x$covariate_terms, collapse = ";"),
      etp_method = x$etp_method %||% NA_character_,
      etp_threshold = x$etp_threshold %||% NA_real_,
      is_primary = isTRUE(x$is_primary),
      stringsAsFactors = FALSE
    )
  }))
}

pst_effective_workers <- function(args, n_tasks, field = "model_workers") {
  requested <- suppressWarnings(as.integer(args[[field]] %||% 0L))
  total <- suppressWarnings(as.integer(args$workers %||% 1L))
  if (!is.finite(total) || total < 1L) total <- 1L
  if (!is.finite(requested) || requested < 1L) requested <- total
  max(1L, min(as.integer(n_tasks), requested))
}

pst_effective_within_model_workers <- function(args, model_workers) {
  requested <- suppressWarnings(as.integer(args$within_model_workers %||% 1L))
  total <- suppressWarnings(as.integer(args$workers %||% 1L))
  model_workers <- suppressWarnings(as.integer(model_workers %||% 1L))
  if (!is.finite(requested) || requested < 1L) requested <- 1L
  if (!is.finite(total) || total < 1L) total <- requested
  if (!is.finite(model_workers) || model_workers < 1L) model_workers <- 1L
  max_by_total <- max(1L, floor(total / model_workers))
  max(1L, min(requested, max_by_total))
}

pst_parallel_lapply <- function(x, fun, workers = 1L, task_label = "task") {
  if (length(x) == 0L) return(list(results = list(), log = data.frame()))
  workers <- max(1L, min(as.integer(workers), length(x)))
  task_names <- names(x)
  if (is.null(task_names)) task_names <- rep("", length(x))
  runner <- function(i) {
    item <- x[[i]]
    started <- Sys.time()
    name <- task_names[[i]]
    if (!nzchar(name)) name <- as.character(item$model_id %||% item[[1L]] %||% i)
    value <- tryCatch(
      fun(item),
      error = function(e) structure(list(error = conditionMessage(e)), class = "pst_task_error")
    )
    ended <- Sys.time()
    list(
      value = value,
      log = data.frame(
        task_group = task_label,
        task_name = name,
        start_time = format(started, "%Y-%m-%d %H:%M:%S %Z"),
        end_time = format(ended, "%Y-%m-%d %H:%M:%S %Z"),
        elapsed_sec = as.numeric(difftime(ended, started, units = "secs")),
        status = if (inherits(value, "pst_task_error")) "failed" else "success",
        error_message = if (inherits(value, "pst_task_error")) value$error else "",
        stringsAsFactors = FALSE
      )
    )
  }
  idx <- seq_along(x)
  if (.Platform$OS.type == "unix" && workers > 1L) {
    rows <- parallel::mclapply(idx, runner, mc.cores = workers)
  } else {
    rows <- lapply(idx, runner)
  }
  task_log <- do.call(rbind, lapply(rows, `[[`, "log"))
  values <- lapply(rows, `[[`, "value")
  names(values) <- task_names
  failed <- vapply(values, inherits, logical(1L), "pst_task_error")
  if (any(failed)) {
    msg <- paste(task_log$task_name[failed], task_log$error_message[failed], sep = ": ", collapse = "; ")
    stop("Parallel task failure in ", task_label, ": ", msg, call. = FALSE)
  }
  list(results = values, log = task_log)
}

pst_run_model_set <- function(model_specs, workflow_root, pb, counts, coverage, coverage_check, cfg, gene_sets, args, task_label) {
  model_workers <- if (isTRUE(args$parallel)) pst_effective_workers(args, length(model_specs), "model_workers") else 1L
  within_model_workers <- if (isTRUE(args$parallel)) pst_effective_within_model_workers(args, model_workers) else 1L
  worker_args <- args
  worker_args$workers <- within_model_workers
  if (!is.null(worker_args$gsea_workers)) worker_args$gsea_workers <- 1L
  tasks <- model_specs
  names(tasks) <- names(model_specs)
  message("Running ", task_label, " models with model_workers=", model_workers, ", within_model_workers=", within_model_workers)
  out <- pst_parallel_lapply(
    tasks,
    function(model_spec) {
      message("Running ", task_label, " model: ", model_spec$model_id)
      pst_run_single_model(
        model_spec,
        workflow_root,
        pb,
        counts,
        coverage,
        coverage_check,
        cfg,
        gene_sets,
        worker_args
      )
    },
    workers = model_workers,
    task_label = paste0(task_label, "_model")
  )
  pst_write_csv(out$log, file.path(workflow_root, "00_manifest", "parallel_task_log.csv"))
  pst_run_model_comparison(out$results, workflow_root)
  out
}

pst_focus_interval_ids <- function() {
  c("left_neighbor", "primary_accumulated_state", "right_neighbor")
}

pst_construct_interval_pseudobulk <- function(meta, counts, cfg, min_cells) {
  sample_meta <- pst_sample_metadata(meta)
  focus_ids <- pst_focus_interval_ids()
  intervals <- cfg$interval_list[focus_ids]
  all_sample_intervals <- expand.grid(
    sample_id = rownames(sample_meta),
    interval_id = focus_ids,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  all_sample_intervals$sample_interval_id <- paste(all_sample_intervals$sample_id, all_sample_intervals$interval_id, sep = "__")
  all_sample_intervals$interval_start <- vapply(all_sample_intervals$interval_id, function(id) intervals[[id]]$start, numeric(1L))
  all_sample_intervals$interval_end <- vapply(all_sample_intervals$interval_id, function(id) intervals[[id]]$end, numeric(1L))
  all_sample_intervals$interval_midpoint <- (all_sample_intervals$interval_start + all_sample_intervals$interval_end) / 2
  all_sample_intervals$interval_role <- vapply(all_sample_intervals$interval_id, function(id) intervals[[id]]$role %||% id, character(1L))

  assigned <- rep(NA_character_, nrow(meta))
  for (id in focus_ids) {
    hit <- pst_interval_hit(meta$pseudotime, intervals[[id]])
    assigned[hit] <- id
  }
  meta_interval <- meta[!is.na(assigned), , drop = FALSE]
  meta_interval$interval_id <- assigned[!is.na(assigned)]
  meta_interval$sample_interval_id <- paste(meta_interval$sample_id, meta_interval$interval_id, sep = "__")

  observed <- as.data.frame(table(meta_interval$sample_interval_id), stringsAsFactors = FALSE)
  names(observed) <- c("sample_interval_id", "cell_count")
  interval_meta <- merge(all_sample_intervals, observed, by = "sample_interval_id", all.x = TRUE, sort = FALSE)
  interval_meta$cell_count[is.na(interval_meta$cell_count)] <- 0L
  interval_meta <- merge(interval_meta, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  interval_meta$interval_id <- factor(interval_meta$interval_id, levels = focus_ids)
  interval_meta <- interval_meta[order(interval_meta$sample_id, interval_meta$interval_id), , drop = FALSE]
  interval_meta$retained_for_model <- interval_meta$cell_count >= min_cells
  interval_meta$exclusion_reason <- ifelse(interval_meta$retained_for_model, "", paste0("cell_count<", min_cells))

  if (nrow(meta_interval) > 0L) {
    group_factor <- factor(meta_interval$sample_interval_id, levels = interval_meta$sample_interval_id)
    design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
    colnames(design) <- levels(group_factor)
    interval_counts <- counts[, meta_interval$cell_id, drop = FALSE] %*% design
    colnames(interval_counts) <- levels(group_factor)
    rownames(interval_counts) <- rownames(counts)
  } else {
    interval_counts <- counts[, integer(0), drop = FALSE]
    interval_counts <- interval_counts[, interval_meta$sample_interval_id, drop = FALSE]
  }
  missing_cols <- setdiff(interval_meta$sample_interval_id, colnames(interval_counts))
  if (length(missing_cols) > 0L) {
    zero <- Matrix::Matrix(0, nrow = nrow(counts), ncol = length(missing_cols), sparse = TRUE)
    rownames(zero) <- rownames(counts)
    colnames(zero) <- missing_cols
    interval_counts <- cbind(interval_counts, zero)
  }
  interval_counts <- interval_counts[, interval_meta$sample_interval_id, drop = FALSE]
  lib_size <- Matrix::colSums(interval_counts)
  interval_meta$library_size <- as.numeric(lib_size[match(interval_meta$sample_interval_id, names(lib_size))])
  interval_meta$library_size[is.na(interval_meta$library_size)] <- 0
  list(counts = interval_counts, metadata = interval_meta, cell_metadata = meta_interval)
}

pst_interval_region_coverage <- function(interval_meta) {
  rows <- interval_meta[, c(
    "sample_id", "interval_id", "interval_role", "interval_start", "interval_end",
    "dose", "dose_mg", "initial_ploidy", "cell_count", "library_size",
    "retained_for_model", "exclusion_reason"
  ), drop = FALSE]
  names(rows)[names(rows) == "interval_role"] <- "role"
  names(rows)[names(rows) == "cell_count"] <- "n_cells"
  names(rows)[names(rows) == "retained_for_model"] <- "contributes"
  rows$n_bins <- NA_integer_
  rows$n_retained_bins <- ifelse(rows$contributes, 1L, 0L)
  rows$exclusion_reason <- ifelse(rows$contributes, "", rows$exclusion_reason)
  rows[, c(
    "sample_id", "interval_id", "role", "interval_start", "interval_end",
    "dose", "dose_mg", "initial_ploidy", "n_cells", "n_bins",
    "n_retained_bins", "library_size", "contributes", "exclusion_reason"
  ), drop = FALSE]
}

pst_prepare_interval_design <- function(meta, model_spec, include_dose = TRUE) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta$interval_factor <- factor(meta$interval_id, levels = pst_focus_interval_ids())
  meta$dose_mg_factor <- factor(meta$dose_mg)
  meta$initial_ploidy_factor <- factor(meta$initial_ploidy)
  for (spec in pst_etp_threshold_specs()) {
    if (spec$group_column %in% names(meta)) {
      meta[[spec$factor_column]] <- factor(as.character(meta[[spec$group_column]]), levels = spec$group_levels)
    }
  }
  if ("sample_mean_endpoint_ploidy" %in% names(meta)) {
    meta$sample_mean_endpoint_ploidy_scaled <- as.numeric(scale(meta$sample_mean_endpoint_ploidy))
  }
  covariate_terms <- model_spec$covariate_terms %||% character(0L)
  model_terms <- c("interval_factor", if (isTRUE(include_dose)) "dose_mg_factor", covariate_terms)
  missing_terms <- setdiff(model_terms, names(meta))
  if (length(missing_terms) > 0L) {
    stop("Model ", model_spec$model_id, " missing design terms: ", paste(missing_terms, collapse = ", "), call. = FALSE)
  }
  formula_text <- paste("~ 0 +", paste(model_terms, collapse = " + "))
  design <- stats::model.matrix(stats::as.formula(formula_text), meta)
  list(
    design = design,
    design_data = meta,
    formula = stats::as.formula(formula_text),
    model_terms = model_terms,
    original_design_columns = colnames(design),
    model_spec = model_spec
  )
}

pst_fit_interval_model <- function(interval_counts, interval_meta, min_primary_mice, include_dose = TRUE, model_spec) {
  keep_obs <- interval_meta$retained_for_model & interval_meta$library_size > 0
  meta <- interval_meta[keep_obs, , drop = FALSE]
  counts <- interval_counts[, meta$sample_interval_id, drop = FALSE]
  y0 <- edgeR::DGEList(counts = counts)
  cpm0 <- edgeR::cpm(y0)
  min_obs <- max(2L, min_primary_mice)
  keep_gene <- rowSums(cpm0 > 1) >= min_obs
  if (sum(keep_gene) < 10L) stop("Too few genes pass expression filtering.", call. = FALSE)
  counts <- counts[keep_gene, , drop = FALSE]
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  design_bundle <- pst_prepare_interval_design(meta, model_spec = model_spec, include_dose = include_dose)
  original_design_columns <- colnames(design_bundle$design)
  qr_rank <- qr(design_bundle$design)$rank
  dropped_design_columns <- character(0L)
  if (qr_rank < ncol(design_bundle$design)) {
    keep_cols <- qr(design_bundle$design)$pivot[seq_len(qr_rank)]
    retained <- colnames(design_bundle$design)[sort(keep_cols)]
    dropped_design_columns <- setdiff(colnames(design_bundle$design), retained)
    design_bundle$design <- design_bundle$design[, retained, drop = FALSE]
  }
  v <- limma::voom(y, design_bundle$design, plot = FALSE)
  block <- meta$sample_id
  dup <- limma::duplicateCorrelation(v, design_bundle$design, block = block)
  fit <- limma::lmFit(v, design_bundle$design, block = block, correlation = dup$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  design_bundle$design <- design_bundle$design[, colnames(fit$coefficients), drop = FALSE]
  list(
    fit = fit,
    voom = v,
    dge = y,
    metadata = meta,
    design = design_bundle,
    correlation = dup$consensus.correlation,
    retained_genes = rownames(counts),
    include_dose = include_dose,
    model_spec = model_spec,
    original_design_columns = original_design_columns,
    retained_design_columns = colnames(design_bundle$design),
    dropped_design_columns = dropped_design_columns,
    design_rank = qr_rank,
    design_ncol_original = length(original_design_columns),
    design_n_observations = nrow(design_bundle$design)
  )
}

pst_interval_contrast_vector <- function(model_fit) {
  cols <- colnames(model_fit$fit$coefficients)
  vec <- rep(0, length(cols))
  names(vec) <- cols
  map <- c(
    interval_factorprimary_accumulated_state = 1,
    interval_factorleft_neighbor = -0.5,
    interval_factorright_neighbor = -0.5
  )
  common <- intersect(names(map), names(vec))
  vec[common] <- map[common]
  vec
}

pst_run_single_interval_model <- function(model_spec, output_root, interval_pb, coverage, coverage_check, gene_sets, args) {
  model_root <- pst_ensure_dir(file.path(output_root, model_spec$model_id))
  dirs <- pst_model_dirs(model_root)
  pst_write_csv(pst_model_parameters_table(model_spec), file.path(dirs$manifest, "model_parameters.csv"))
  min_primary_mice <- coverage_check$n_contributing_mice[[1L]]

  model_fit <- pst_fit_interval_model(
    interval_pb$counts,
    interval_pb$metadata,
    min_primary_mice,
    include_dose = TRUE,
    model_spec = model_spec
  )
  design_audit <- pst_design_audit(model_fit)
  pst_write_csv(design_audit, file.path(dirs$qc, "model_design_rank_audit.csv"))

  primary_contrast <- pst_interval_contrast_vector(model_fit)
  primary_genes <- pst_contrast_table(model_fit, primary_contrast, "primary_adjacent_state")
  primary_genes$model_id <- model_spec$model_id
  symbol_resolution <- pst_resolve_gene_symbols(primary_genes)
  symbol_resolution$model_id <- model_spec$model_id
  pst_write_csv(primary_genes, file.path(dirs$gene_models, "gene_primary_adjacent_state_contrast.csv"))
  pst_write_csv(symbol_resolution, file.path(dirs$gene_models, "gene_symbol_resolution.csv"))

  primary_stats <- pst_ranked_stats(primary_genes, symbol_resolution)
  primary_gsea <- pst_run_all_gsea(primary_stats, gene_sets, "primary_adjacent_state", args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
  primary_gsea$model_id <- model_spec$model_id
  pst_write_csv(primary_gsea, file.path(dirs$gsea, "all_collections_primary_adjacent_state_gsea.csv"))
  if (nrow(primary_gsea) > 0L) {
    for (collection_label in unique(primary_gsea$collection_label)) {
      local <- primary_gsea[primary_gsea$collection_label == collection_label, , drop = FALSE]
      pst_write_csv(local, file.path(dirs$gsea, paste0(collection_label, "_primary_adjacent_state_gsea.csv")))
      pst_write_csv(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea_plot_data.csv")))
      pst_plot_gsea_bar(local, file.path(dirs$figures, paste0("primary_state_", collection_label, "_gsea.pdf")), paste0(model_spec$model_id, ": ", collection_label, " GSEA"))
    }
  }
  leading_edge <- pst_leading_edge_table(primary_gsea)
  if (nrow(leading_edge) > 0L) leading_edge$model_id <- model_spec$model_id
  pst_write_csv(leading_edge, file.path(dirs$gsea, "all_collections_leading_edge_genes.csv"))
  if (nrow(leading_edge) > 0L) {
    for (collection_label in unique(leading_edge$collection_label)) {
      pst_write_csv(leading_edge[leading_edge$collection_label == collection_label, , drop = FALSE], file.path(dirs$gsea, paste0(collection_label, "_leading_edge_genes.csv")))
    }
  }
  pst_write_csv(coverage, file.path(dirs$figures, "sample_region_coverage_plot_data.csv"))
  pst_plot_coverage(coverage, file.path(dirs$figures, "sample_region_coverage.pdf"))

  list(
    model_id = model_spec$model_id,
    model_spec = model_spec,
    model_root = model_root,
    design_audit = design_audit,
    primary_genes = primary_genes,
    secondary_genes = data.frame(),
    primary_gsea = primary_gsea,
    secondary_gsea = data.frame(),
    robustness = data.frame(),
    loo = data.frame(),
    loo_stability = data.frame()
  )
}

pst_run_interval_model_set <- function(model_specs, workflow_root, interval_pb, coverage, coverage_check, gene_sets, args, task_label) {
  model_workers <- if (isTRUE(args$parallel)) pst_effective_workers(args, length(model_specs), "model_workers") else 1L
  worker_args <- args
  worker_args$workers <- 1L
  if (!is.null(worker_args$gsea_workers)) worker_args$gsea_workers <- 1L
  tasks <- model_specs
  names(tasks) <- names(model_specs)
  message("Running ", task_label, " models with workers=", model_workers)
  out <- pst_parallel_lapply(
    tasks,
    function(model_spec) {
      message("Running ", task_label, " model: ", model_spec$model_id)
      pst_run_single_interval_model(
        model_spec,
        workflow_root,
        interval_pb,
        coverage,
        coverage_check,
        gene_sets,
        worker_args
      )
    },
    workers = model_workers,
    task_label = paste0(task_label, "_model")
  )
  pst_write_csv(out$log, file.path(workflow_root, "00_manifest", "parallel_task_log.csv"))
  pst_run_model_comparison(out$results, workflow_root)
  out
}

pst_cross_workflow_comparison <- function(binning_results, non_binning_results, output_root) {
  comparison_root <- pst_ensure_dir(file.path(output_root, "cross_workflow_comparison"))
  figures_dir <- pst_ensure_dir(file.path(comparison_root, "figures"))
  common_models <- intersect(names(binning_results), names(non_binning_results))
  gene_rows <- lapply(common_models, function(model_id) {
    b <- binning_results[[model_id]]$primary_genes[, c("gene", "gene_symbol", "t_statistic"), drop = FALSE]
    n <- non_binning_results[[model_id]]$primary_genes[, c("gene", "t_statistic"), drop = FALSE]
    names(b)[names(b) == "t_statistic"] <- "binning_t_statistic"
    names(n)[names(n) == "t_statistic"] <- "non_binning_t_statistic"
    merged <- merge(b, n, by = "gene", all = FALSE, sort = FALSE)
    data.frame(
      model_id = model_id,
      n_genes = nrow(merged),
      pearson_r = pst_safe_cor(merged$binning_t_statistic, merged$non_binning_t_statistic, "pearson"),
      spearman_rho = pst_safe_cor(merged$binning_t_statistic, merged$non_binning_t_statistic, "spearman"),
      stringsAsFactors = FALSE
    )
  })
  gene_cor <- do.call(rbind, gene_rows)

  pathway_rows <- lapply(common_models, function(model_id) {
    b <- binning_results[[model_id]]$primary_gsea
    n <- non_binning_results[[model_id]]$primary_gsea
    b$workflow <- "binning"
    n$workflow <- "non_binning"
    rbind(b, n)
  })
  pathway_long <- do.call(rbind, pathway_rows)
  pst_write_csv(pathway_long, file.path(comparison_root, "pathway_workflow_comparison_long.csv"))

  nes_rows <- list()
  overlap_rows <- list()
  rank_rows <- list()
  for (model_id in common_models) {
    b <- binning_results[[model_id]]$primary_gsea
    n <- non_binning_results[[model_id]]$primary_gsea
    merged <- merge(
      b[, c("collection", "collection_label", "pathway", "pathway_label", "NES", "padj"), drop = FALSE],
      n[, c("collection", "collection_label", "pathway", "NES", "padj"), drop = FALSE],
      by = c("collection", "collection_label", "pathway"),
      suffixes = c("_binning", "_non_binning"),
      all = FALSE,
      sort = FALSE
    )
    for (collection_label in unique(merged$collection_label)) {
      sub <- merged[merged$collection_label == collection_label, , drop = FALSE]
      nes_rows[[length(nes_rows) + 1L]] <- data.frame(
        model_id = model_id,
        collection_label = collection_label,
        n_pathways = nrow(sub),
        pearson_r = pst_safe_cor(sub$NES_binning, sub$NES_non_binning, "pearson"),
        spearman_rho = pst_safe_cor(sub$NES_binning, sub$NES_non_binning, "spearman"),
        stringsAsFactors = FALSE
      )
      b_set <- sub$pathway[sub$padj_binning < 0.05]
      n_set <- sub$pathway[sub$padj_non_binning < 0.05]
      overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
        model_id = model_id,
        collection_label = collection_label,
        alpha = 0.05,
        binning_significant_n = length(b_set),
        non_binning_significant_n = length(n_set),
        overlap_n = length(intersect(b_set, n_set)),
        jaccard = if (length(union(b_set, n_set)) == 0L) NA_real_ else length(intersect(b_set, n_set)) / length(union(b_set, n_set)),
        stringsAsFactors = FALSE
      )
      sub$binning_rank <- rank(sub$padj_binning, ties.method = "first", na.last = "keep")
      sub$non_binning_rank <- rank(sub$padj_non_binning, ties.method = "first", na.last = "keep")
      sub$rank_shift <- sub$non_binning_rank - sub$binning_rank
      sub$NES_delta <- sub$NES_non_binning - sub$NES_binning
      sub$model_id <- model_id
      rank_rows[[length(rank_rows) + 1L]] <- sub
    }
  }
  nes_cor <- do.call(rbind, nes_rows)
  sig_overlap <- do.call(rbind, overlap_rows)
  rank_shift <- do.call(rbind, rank_rows)
  rank_shift <- rank_shift[order(rank_shift$model_id, rank_shift$collection_label, -abs(rank_shift$rank_shift), rank_shift$pathway), , drop = FALSE]
  pst_write_csv(gene_cor, file.path(comparison_root, "gene_t_stat_correlations_binning_vs_non_binning.csv"))
  pst_write_csv(nes_cor, file.path(comparison_root, "gsea_NES_correlations_binning_vs_non_binning.csv"))
  pst_write_csv(sig_overlap, file.path(comparison_root, "gsea_significant_overlap_binning_vs_non_binning.csv"))
  pst_write_csv(rank_shift, file.path(comparison_root, "pathway_rank_shift_binning_vs_non_binning.csv"))
  pst_plot_cross_workflow_nes(rank_shift, file.path(figures_dir, "gsea_NES_correlation_binning_vs_non_binning.pdf"))
  pst_plot_cross_workflow_rank_shift(rank_shift, file.path(figures_dir, "top_pathway_rank_shift_binning_vs_non_binning.pdf"))
  invisible(list(gene_cor = gene_cor, nes_cor = nes_cor, sig_overlap = sig_overlap, rank_shift = rank_shift))
}

pst_plot_cross_workflow_nes <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES_binning, y = NES_non_binning)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(alpha = 0.45, size = 0.9) +
    ggplot2::facet_grid(model_id ~ collection_label, scales = "free") +
    ggplot2::labs(title = "Primary-state GSEA NES: binning versus non-binning", x = "Binning NES", y = "Non-binning NES") +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(path, p, width = 12, height = 9)
  invisible(p)
}

pst_plot_cross_workflow_rank_shift <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  top <- head(df[order(-abs(df$rank_shift)), , drop = FALSE], 100L)
  label_col <- if ("pathway_label" %in% names(top)) "pathway_label" else "pathway_label_binning"
  top$pathway_label <- factor(top[[label_col]], levels = rev(unique(top[[label_col]])))
  p <- ggplot2::ggplot(top, ggplot2::aes(x = rank_shift, y = pathway_label, color = collection_label)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~ model_id, scales = "free_y") +
    ggplot2::labs(title = "Largest pathway rank shifts: binning versus non-binning", x = "Non-binning rank minus binning rank", y = "Pathway", color = "Collection") +
    ggplot2::theme_bw(base_size = 8)
  ggplot2::ggsave(path, p, width = 13, height = 10)
  invisible(p)
}

pst_write_readme <- function(output_root) {
  lines <- c(
    "# 04i Pseudotime State Pathways",
    "",
    "This workflow annotates the frozen CellCycle pseudotime state in which gemcitabine-treated tumor cells accumulate.",
    "",
    "The output is organized into two complementary workflows:",
    "",
    "- `binning`: sample-by-pseudotime-bin pseudobulk with a smooth common pseudotime model.",
    "- `non_binning`: interval-level pseudobulk that aggregates cells directly within the frozen left, primary, and right intervals.",
    "- `cross_workflow_comparison`: concordance between the binning and non-binning workflows for matched model specifications.",
    "",
    "Both workflows run the same five covariate specifications: primary initial ploidy, continuous mean ETP, and three thresholded ETP group models.",
    "",
    "Positive NES pathways characterize the accumulated state relative to adjacent states. They must not be described as directly treatment-upregulated pathways."
  )
  pst_write_lines(lines, file.path(output_root, "README.md"))
}

pst_run_workflow <- function(args, repo_root, script_dir) {
  pst_check_packages()
  set.seed(args$seed)
  output_root <- pst_clean_dir(args$output_root, overwrite = args$overwrite)
  etp_specs <- pst_resolve_etp_specs(args$etp_threshold)
  model_specs <- pst_model_specs(args$covariate_mode, etp_specs, include_full_etp_models = args$include_full_etp_models)
  dirs <- list(
    manifest = pst_ensure_dir(file.path(output_root, "00_manifest")),
    shared_qc = pst_ensure_dir(file.path(output_root, "shared_qc")),
    binning = pst_ensure_dir(file.path(output_root, "binning")),
    non_binning = pst_ensure_dir(file.path(output_root, "non_binning"))
  )
  binning_dirs <- list(
    manifest = pst_ensure_dir(file.path(dirs$binning, "00_manifest")),
    qc = pst_ensure_dir(file.path(dirs$binning, "01_qc")),
    pseudobulk = pst_ensure_dir(file.path(dirs$binning, "02_pseudobulk"))
  )
  non_binning_dirs <- list(
    manifest = pst_ensure_dir(file.path(dirs$non_binning, "00_manifest")),
    qc = pst_ensure_dir(file.path(dirs$non_binning, "01_qc")),
    interval_pseudobulk = pst_ensure_dir(file.path(dirs$non_binning, "02_interval_pseudobulk"))
  )

  cfg <- pst_read_config(args$config)
  seurat_rds <- pst_resolve_seurat_rds(args$seurat_rds, args$results_root)
  params <- data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE)
  pst_write_csv(params, file.path(dirs$manifest, "analysis_parameters.csv"))
  pst_write_csv(pst_intervals_table(cfg), file.path(dirs$manifest, "frozen_interval_definition.csv"))
  pst_write_csv(pst_package_versions(), file.path(dirs$manifest, "package_versions.csv"))
  selected_model_table <- pst_model_specifications_table(model_specs)
  pst_write_csv(selected_model_table, file.path(dirs$manifest, "selected_model_specifications.csv"))
  pst_write_csv(params, file.path(binning_dirs$manifest, "analysis_parameters.csv"))
  pst_write_csv(params, file.path(non_binning_dirs$manifest, "analysis_parameters.csv"))

  input_checksums <- data.frame(
    input = c("cell_metadata", "noncell_metadata", "seurat_rds", "config"),
    path = c(args$cell_metadata, args$noncell_metadata, seurat_rds, args$config),
    sha256 = c(
      pst_file_checksum(args$cell_metadata),
      pst_file_checksum(args$noncell_metadata),
      pst_file_checksum(seurat_rds),
      pst_file_checksum(args$config)
    ),
    stringsAsFactors = FALSE
  )
  pst_write_csv(input_checksums, file.path(dirs$manifest, "input_checksums.csv"))

  meta <- pst_read_cell_metadata(args$cell_metadata)
  etp_assignments <- pst_build_etp_assignments(args$cell_metadata, args$noncell_metadata, etp_specs)
  pst_write_etp_audits(etp_assignments, etp_specs, dirs$shared_qc)
  meta <- pst_attach_etp_assignments(meta, etp_assignments)
  metadata_audit <- pst_audit_metadata(meta)
  pst_write_csv(metadata_audit$duplicate_cells, file.path(dirs$shared_qc, "duplicate_cell_ids.csv"))
  pst_write_csv(metadata_audit$inconsistent_cells, file.path(dirs$shared_qc, "cell_assignment_consistency_audit.csv"))
  pst_check_pseudotime(meta)

  counts <- pst_load_counts(seurat_rds, args$assay, args$counts_layer)
  expression_audit <- pst_audit_counts(counts)
  expression_audit$seurat_rds <- seurat_rds
  expression_audit$assay <- args$assay
  expression_audit$counts_layer <- args$counts_layer
  pst_write_csv(expression_audit, file.path(dirs$shared_qc, "expression_source_audit.csv"))
  if (!isTRUE(expression_audit$non_negative) || !isTRUE(expression_audit$integer_like)) {
    stop("Expression matrix failed raw-count audit. See shared_qc/expression_source_audit.csv.", call. = FALSE)
  }

  matched <- pst_match_cells(meta, counts, args$min_match_rate, dirs$shared_qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- pst_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts

  pb <- pst_construct_pseudobulk(meta, counts, args$n_pseudotime_bins, args$min_cells_per_sample_bin)
  pst_write_csv(pb$metadata, file.path(binning_dirs$pseudobulk, "sample_bin_metadata.csv"))
  coverage <- pst_region_coverage(pb$cell_metadata, pb$metadata, cfg)
  pst_write_csv(coverage, file.path(binning_dirs$qc, "sample_region_coverage.csv"))
  coverage_check <- pst_check_primary_coverage(coverage)
  pst_write_csv(coverage_check, file.path(binning_dirs$qc, "primary_coverage_check.csv"))
  if (!isTRUE(coverage_check$passes)) {
    stop("Binning primary interval coverage check failed. See binning/01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  interval_pb <- pst_construct_interval_pseudobulk(meta, counts, cfg, args$min_cells_per_sample_region)
  pst_write_csv(interval_pb$metadata, file.path(non_binning_dirs$interval_pseudobulk, "interval_sample_metadata.csv"))
  interval_coverage <- pst_interval_region_coverage(interval_pb$metadata)
  pst_write_csv(interval_coverage, file.path(non_binning_dirs$qc, "interval_region_coverage.csv"))
  interval_coverage_check <- pst_check_primary_coverage(interval_coverage)
  pst_write_csv(interval_coverage_check, file.path(non_binning_dirs$qc, "primary_coverage_check.csv"))
  if (!isTRUE(interval_coverage_check$passes)) {
    stop("Non-binning primary interval coverage check failed. See non_binning/01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  gene_sets <- pst_fetch_gene_sets(args$gene_set_collections)
  binning_run <- pst_run_model_set(model_specs, dirs$binning, pb, counts, coverage, coverage_check, cfg, gene_sets, args, "binning")
  non_binning_run <- pst_run_interval_model_set(model_specs, dirs$non_binning, interval_pb, interval_coverage, interval_coverage_check, gene_sets, args, "non_binning")
  pst_cross_workflow_comparison(binning_run$results, non_binning_run$results, output_root)

  pst_write_readme(output_root)
  if (isTRUE(args$make_snapshot)) pst_materialize_snapshot(output_root, args$snapshot_root)
  message("04i workflow complete: ", output_root)
  invisible(output_root)
}
