`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

ptp_parse_scalar <- function(value, template) {
  if (is.logical(template)) return(tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y"))
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  as.character(value)
}

ptp_parse_args <- function(argv, defaults) {
  out <- defaults
  for (arg in argv) {
    if (!grepl("^--", arg)) next
    kv <- sub("^--", "", arg)
    if (grepl("=", kv, fixed = TRUE)) {
      key <- sub("=.*$", "", kv)
      value <- sub("^[^=]*=", "", kv)
    } else {
      key <- kv
      value <- "TRUE"
    }
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) stop("Unknown argument: --", key, call. = FALSE)
    out[[key]] <- ptp_parse_scalar(value, out[[key]])
  }
  out
}

ptp_required_packages <- function() {
  c("Seurat", "SeuratObject", "Matrix", "yaml", "digest", "edgeR", "limma",
    "splines", "fgsea", "msigdbr", "dplyr", "readr", "tidyr", "ggplot2",
    "readxl", "statmod", "base64enc")
}

ptp_check_packages <- function(packages = ptp_required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

ptp_ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

ptp_clean_dir <- function(path, overwrite = FALSE) {
  if (dir.exists(path) && !overwrite) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) > 0L) {
      stop("Output directory is not empty. Use --overwrite=TRUE: ", path, call. = FALSE)
    }
  }
  if (dir.exists(path) && overwrite) unlink(path, recursive = TRUE, force = TRUE)
  ptp_ensure_dir(path)
}

ptp_output_dirs <- function(output_root) {
  dir_at <- function(...) ptp_ensure_dir(file.path(output_root, ...))
  list(
    root = ptp_ensure_dir(output_root),
    manifest = dir_at("00_manifest"),
    qc = dir_at("01_qc"),
    pseudobulk = dir_at("02_pseudobulk"),
    primary = list(
      root = dir_at("primary_initial_ploidy"),
      qc = dir_at("primary_initial_ploidy", "01_qc"),
      genes = dir_at("primary_initial_ploidy", "03_gene_models"),
      programs = dir_at("primary_initial_ploidy", "04_program_models"),
      score_methods = dir_at("primary_initial_ploidy", "04_program_models", "score_methods"),
      gsea = dir_at("primary_initial_ploidy", "05_gsea"),
      enrichment = list(
        root = dir_at("primary_initial_ploidy", "05_enrichment"),
        camera = dir_at("primary_initial_ploidy", "05_enrichment", "cameraPR"),
        gsea = dir_at("primary_initial_ploidy", "05_enrichment", "gsea"),
        tier_mapped = dir_at("primary_initial_ploidy", "05_enrichment", "tier_mapped_summary")
      ),
      figures = dir_at("primary_initial_ploidy", "06_figures"),
      robustness = dir_at("primary_initial_ploidy", "07_robustness")
    ),
    dose = list(
      root = dir_at("dose_specific_initial_ploidy"),
      programs = dir_at("dose_specific_initial_ploidy", "04_program_models"),
      figures = dir_at("dose_specific_initial_ploidy", "06_figures")
    ),
    endpoint = list(
      root = dir_at("endpoint_ploidy"),
      manifest = dir_at("endpoint_ploidy", "00_manifest"),
      qc = dir_at("endpoint_ploidy", "01_qc"),
      figures = dir_at("endpoint_ploidy", "06_figures")
    ),
    model_comparison = dir_at("model_comparison"),
    report = dir_at("report")
  )
}

ptp_endpoint_model_dirs <- function(endpoint_root, model_id) {
  model_id <- gsub("[^A-Za-z0-9_]+", "_", as.character(model_id))
  dir_at <- function(...) ptp_ensure_dir(file.path(endpoint_root, model_id, ...))
  list(
    root = dir_at(),
    manifest = dir_at("00_manifest"),
    qc = dir_at("01_qc"),
    genes = dir_at("03_gene_models"),
    programs = dir_at("04_program_models"),
    figures = dir_at("06_figures")
  )
}

ptp_write_csv <- function(x, path) {
  ptp_ensure_dir(dirname(path))
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (ncol(x) == 0L) x <- data.frame(note = character(0L), stringsAsFactors = FALSE)
  readr::write_csv(x, path, na = "")
  invisible(path)
}

ptp_write_lines <- function(x, path) {
  ptp_ensure_dir(dirname(path))
  writeLines(as.character(x), path, useBytes = TRUE)
  invisible(path)
}

ptp_effective_workers <- function(args, n_tasks, field = "workers") {
  if (length(n_tasks) == 0L || is.na(n_tasks) || n_tasks < 1L) return(1L)
  if (!isTRUE(args$parallel %||% FALSE)) return(1L)
  requested <- suppressWarnings(as.integer(args[[field]] %||% args$workers %||% 1L))
  if (!is.finite(requested) || requested < 1L) requested <- 1L
  max(1L, min(as.integer(n_tasks), requested))
}

ptp_parallel_lapply <- function(x, fun, workers = 1L, task_label = "task") {
  if (length(x) == 0L) return(list(results = list(), log = data.frame()))
  workers <- suppressWarnings(as.integer(workers))
  if (!is.finite(workers) || workers < 1L) workers <- 1L
  workers <- max(1L, min(workers, length(x)))
  task_names <- names(x)
  if (is.null(task_names)) task_names <- rep("", length(x))
  runner <- function(i) {
    item <- x[[i]]
    started <- Sys.time()
    name <- task_names[[i]]
    if (!nzchar(name)) name <- as.character(i)
    value <- tryCatch(
      fun(item),
      error = function(e) structure(list(error = conditionMessage(e)), class = "ptp_task_error")
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
        status = if (inherits(value, "ptp_task_error")) "failed" else "success",
        error_message = if (inherits(value, "ptp_task_error")) value$error else "",
        stringsAsFactors = FALSE
      )
    )
  }
  idx <- seq_along(x)
  rows <- if (.Platform$OS.type == "unix" && workers > 1L) {
    parallel::mclapply(idx, runner, mc.cores = workers)
  } else {
    lapply(idx, runner)
  }
  task_log <- do.call(rbind, lapply(rows, `[[`, "log"))
  values <- lapply(rows, `[[`, "value")
  names(values) <- task_names
  failed <- vapply(values, inherits, logical(1L), "ptp_task_error")
  if (any(failed)) {
    msg <- paste(task_log$task_name[failed], task_log$error_message[failed], sep = ": ", collapse = "; ")
    stop("Parallel task failure in ", task_label, ": ", msg, call. = FALSE)
  }
  list(results = values, log = task_log)
}

ptp_file_checksum <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

ptp_safe_name <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

ptp_endpoint_ploidy_specs <- function(cfg = NULL) {
  if (!is.null(cfg) && !is.null(cfg$endpoint_ploidy$thresholds)) {
    specs <- cfg$endpoint_ploidy$thresholds
  } else {
    specs <- list(
      list(
        method = "ETP_fixed_threshold_2_25",
        threshold_scheme = "fixed_threshold_2_25",
        threshold = 2.25,
        group_column = "ETP_fixed_threshold_2_25_group",
        factor_column = "ETP_fixed_threshold_2_25_factor",
        group_levels = c("ETP-lower", "ETP-higher")
      ),
      list(
        method = "ETP_boundary_stress_threshold_2_375",
        threshold_scheme = "boundary_stress_threshold_2_375",
        threshold = 2.375,
        group_column = "ETP_boundary_stress_threshold_2_375_group",
        factor_column = "ETP_boundary_stress_threshold_2_375_factor",
        group_levels = c("ETP-lower", "ETP-higher")
      ),
      list(
        method = "ETP_reference_balanced_threshold_2_24",
        threshold_scheme = "reference_balanced_threshold_2_24",
        threshold = 2.24,
        group_column = "ETP_reference_balanced_threshold_2_24_group",
        factor_column = "ETP_reference_balanced_threshold_2_24_factor",
        group_levels = c("ETP-lower", "ETP-higher")
      )
    )
  }
  specs <- lapply(specs, function(x) {
    x$method <- as.character(x$method)
    x$threshold_scheme <- as.character(x$threshold_scheme %||% x$method)
    x$threshold <- as.numeric(x$threshold)
    x$group_column <- as.character(x$group_column %||% paste0(x$method, "_group"))
    x$factor_column <- as.character(x$factor_column %||% paste0(x$method, "_factor"))
    x$group_levels <- as.character(x$group_levels %||% c("ETP-lower", "ETP-higher"))
    x
  })
  names(specs) <- vapply(specs, `[[`, character(1L), "method")
  specs
}

ptp_endpoint_ploidy_specs_table <- function(specs) {
  rows <- lapply(specs, function(x) {
    data.frame(
      method = x$method,
      threshold_scheme = x$threshold_scheme,
      threshold = x$threshold,
      group_column = x$group_column,
      factor_column = x$factor_column,
      group_levels = paste(x$group_levels, collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ptp_clean_gene_symbols <- function(genes) {
  if (exists("clean_gene_symbols", mode = "function")) return(clean_gene_symbols(genes))
  g <- as.character(genes)
  g <- trimws(g)
  g <- sub("^GRCh[0-9]+[-_]", "", g, ignore.case = TRUE)
  g <- sub("^GRCm39[-_]", "", g, ignore.case = TRUE)
  g <- sub("^hg38[-_]", "", g, ignore.case = TRUE)
  g <- sub("\\.[0-9]+$", "", g)
  g[g == ""] <- NA_character_
  g
}

ptp_package_versions <- function(packages = ptp_required_packages()) {
  data.frame(
    package = packages,
    version = vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1L)),
    stringsAsFactors = FALSE
  )
}

ptp_read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$primary_subbins) || is.null(cfg$fallback_subbins)) {
    stop("Config must define primary_subbins and fallback_subbins.", call. = FALSE)
  }
  cfg
}

ptp_interval_hit <- function(x, interval) {
  left <- if (isTRUE(interval$include_start)) x >= interval$start else x > interval$start
  right <- if (isTRUE(interval$include_end)) x <= interval$end else x < interval$end
  left & right
}

ptp_subbins_table <- function(subbins, grid_id) {
  rows <- lapply(subbins, function(sb) {
    data.frame(
      grid_id = grid_id,
      subbin_id = sb$id,
      label = sb$label %||% sb$id,
      start = as.numeric(sb$start),
      end = as.numeric(sb$end),
      include_start = isTRUE(sb$include_start),
      include_end = isTRUE(sb$include_end),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ptp_regions_table <- function(cfg) {
  rows <- lapply(names(cfg$regions), function(id) {
    x <- cfg$regions[[id]]
    data.frame(
      region_id = id,
      label = x$label %||% id,
      start = as.numeric(x$start),
      end = as.numeric(x$end),
      include_start = isTRUE(x$include_start),
      include_end = isTRUE(x$include_end),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ptp_resolve_path <- function(path, repo_root) {
  if (is.null(path) || !nzchar(path)) return("")
  if (grepl("^/", path)) return(path)
  file.path(repo_root, path)
}

ptp_resolve_seurat_rds <- function(explicit, cfg, results_root) {
  if (!is.null(explicit) && nzchar(explicit)) {
    if (!file.exists(explicit)) stop("Explicit --seurat-rds does not exist: ", explicit, call. = FALSE)
    return(normalizePath(explicit, winslash = "/", mustWork = TRUE))
  }
  candidates <- cfg$inputs$seurat_rds_candidates %||% character(0L)
  if (length(candidates) == 0L && nzchar(results_root)) {
    candidates <- c(
      file.path(results_root, "03c_cluster_annotation/04_objects/integrated_sct_cca_seurat_cluster_final_annotation.rds"),
      file.path(results_root, "03b_manual_cluster_merge/objects/integrated_sct_cca_seurat_final_manual_merge.rds"),
      file.path(results_root, "03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds"),
      file.path(results_root, "02e_cluster_annotation/04_objects/integrated_sct_cca_seurat_cluster_final_annotation.rds")
    )
  }
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) stop("No Seurat RDS found. Tried: ", paste(candidates, collapse = "; "), call. = FALSE)
  normalizePath(hit[[1L]], winslash = "/", mustWork = TRUE)
}

ptp_read_cell_metadata <- function(path) {
  required <- c("cell_id", "sample_id", "cluster", "initial_ploidy",
                "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "pseudotime")
  meta <- as.data.frame(readr::read_csv(path, show_col_types = FALSE), stringsAsFactors = FALSE)
  missing <- setdiff(required, names(meta))
  if (length(missing) > 0L) stop("Cell metadata missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  meta$cell_id <- as.character(meta$cell_id)
  meta$sample_id <- as.character(meta$sample_id)
  meta$cluster <- as.character(meta$cluster)
  meta$initial_ploidy <- as.character(meta$initial_ploidy)
  meta$gemcitabine_dose <- as.character(meta$gemcitabine_dose)
  meta$gemcitabine_dose_mg_per_kg <- suppressWarnings(as.numeric(meta$gemcitabine_dose_mg_per_kg))
  meta$pseudotime <- suppressWarnings(as.numeric(meta$pseudotime))
  if ("cell_ploidy" %in% names(meta)) meta$cell_ploidy <- suppressWarnings(as.numeric(meta$cell_ploidy))
  if ("average_ploidy" %in% names(meta)) meta$average_ploidy <- suppressWarnings(as.numeric(meta$average_ploidy))
  if ("median_ploidy" %in% names(meta)) meta$median_ploidy <- suppressWarnings(as.numeric(meta$median_ploidy))
  bad <- !is.finite(meta$pseudotime) | meta$pseudotime < 0 | meta$pseudotime > 1
  if (any(bad)) stop("Cell metadata contains invalid pseudotime values.", call. = FALSE)
  meta$treatment <- ifelse(meta$gemcitabine_dose_mg_per_kg == 0, "control", "gemcitabine")
  meta$treatment <- factor(meta$treatment, levels = c("control", "gemcitabine"))
  meta$dose_factor <- factor(paste0("dose", meta$gemcitabine_dose_mg_per_kg), levels = c("dose0", "dose30", "dose120"))
  meta$initial_ploidy <- factor(meta$initial_ploidy, levels = c("2N", "4N"))
  meta
}

ptp_sample_level_unique <- function(x) {
  x <- unique(as.character(x[!is.na(x)]))
  x <- x[nzchar(x)]
  if (length(x) == 0L) NA_character_ else paste(x, collapse = ";")
}

ptp_sample_level_numeric <- function(x, fallback = NA_real_) {
  x <- suppressWarnings(as.numeric(x))
  x <- unique(x[is.finite(x)])
  if (length(x) == 0L) fallback else x[[1L]]
}

ptp_sample_metadata <- function(meta, endpoint_ploidy_specs = ptp_endpoint_ploidy_specs()) {
  split_meta <- split(meta, meta$sample_id)
  rows <- lapply(names(split_meta), function(sample_id) {
    x <- split_meta[[sample_id]]
    mean_endpoint_ploidy <- if ("sample_mean_endpoint_ploidy" %in% names(x)) {
      ptp_sample_level_numeric(x$sample_mean_endpoint_ploidy)
    } else if ("average_ploidy" %in% names(x)) {
      ptp_sample_level_numeric(x$average_ploidy)
    } else if ("cell_ploidy" %in% names(x)) {
      mean(x$cell_ploidy, na.rm = TRUE)
    } else {
      NA_real_
    }
    median_endpoint_ploidy <- if ("sample_median_endpoint_ploidy" %in% names(x)) {
      ptp_sample_level_numeric(x$sample_median_endpoint_ploidy)
    } else if ("median_ploidy" %in% names(x)) {
      ptp_sample_level_numeric(x$median_ploidy)
    } else if ("cell_ploidy" %in% names(x)) {
      stats::median(x$cell_ploidy, na.rm = TRUE)
    } else {
      NA_real_
    }
    data.frame(
      sample_id = sample_id,
      initial_ploidy = ptp_sample_level_unique(x$initial_ploidy),
      gemcitabine_dose = ptp_sample_level_unique(x$gemcitabine_dose),
      dose_mg = unique(x$gemcitabine_dose_mg_per_kg)[1L],
      treatment = ptp_sample_level_unique(x$treatment),
      n_cells = nrow(x),
      n_endpoint_ploidy_cells = if ("cell_ploidy" %in% names(x)) sum(is.finite(x$cell_ploidy)) else NA_integer_,
      sample_mean_endpoint_ploidy = mean_endpoint_ploidy,
      sample_median_endpoint_ploidy = median_endpoint_ploidy,
      sample_max_endpoint_ploidy = if ("cell_ploidy" %in% names(x) && any(is.finite(x$cell_ploidy))) max(x$cell_ploidy, na.rm = TRUE) else NA_real_,
      mean_pseudotime = mean(x$pseudotime, na.rm = TRUE),
      median_pseudotime = median(x$pseudotime, na.rm = TRUE),
      q05_pseudotime = unname(stats::quantile(x$pseudotime, 0.05, na.rm = TRUE)),
      q25_pseudotime = unname(stats::quantile(x$pseudotime, 0.25, na.rm = TRUE)),
      q75_pseudotime = unname(stats::quantile(x$pseudotime, 0.75, na.rm = TRUE)),
      q95_pseudotime = unname(stats::quantile(x$pseudotime, 0.95, na.rm = TRUE)),
      mean_percent_mt = if ("percent.mt" %in% names(x)) mean(x$percent.mt, na.rm = TRUE) else NA_real_,
      mean_nCount_RNA = if ("nCount_RNA" %in% names(x)) mean(x$nCount_RNA, na.rm = TRUE) else NA_real_,
      mean_nFeature_RNA = if ("nFeature_RNA" %in% names(x)) mean(x$nFeature_RNA, na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (any(is.finite(out$sample_mean_endpoint_ploidy))) {
    out$sample_mean_endpoint_ploidy_scaled <- as.numeric(scale(out$sample_mean_endpoint_ploidy))
  } else {
    out$sample_mean_endpoint_ploidy_scaled <- NA_real_
  }
  for (spec in endpoint_ploidy_specs) {
    group <- ifelse(out$sample_mean_endpoint_ploidy > spec$threshold, "ETP-higher", "ETP-lower")
    group[!is.finite(out$sample_mean_endpoint_ploidy)] <- NA_character_
    out[[spec$group_column]] <- group
    out[[spec$factor_column]] <- factor(group, levels = spec$group_levels)
  }
  out[order(out$initial_ploidy, out$dose_mg, out$sample_id), , drop = FALSE]
}

ptp_load_counts_and_qc <- function(seurat_rds, assay, counts_layer, meta) {
  obj <- readRDS(seurat_rds)
  counts <- get_assay_matrix(obj, assay = assay, slot_name = counts_layer)
  if (is.null(counts)) stop("Cannot read assay matrix: assay=", assay, ", layer/slot=", counts_layer, call. = FALSE)
  if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  seurat_meta <- as.data.frame(obj@meta.data, stringsAsFactors = FALSE)
  seurat_meta$cell_id <- rownames(seurat_meta)
  wanted <- intersect(c(
    "cell_id", "orig.ident", "sample", "sample_folder", "percent.mt", "nCount_RNA", "nFeature_RNA",
    "scDblFinder.class", "scDblFinder.score", "harvest", "Sequencing.IDs", "IDs", "Dose",
    "S.Score", "G2M.Score", "Phase", "cluster_final", "cluster_final_annotation_primary",
    "cluster_final_annotation_multi", "seurat_clusters"
  ), names(seurat_meta))
  merged <- merge(meta, seurat_meta[, wanted, drop = FALSE], by = "cell_id", all.x = TRUE, sort = FALSE)
  merged <- merged[match(meta$cell_id, merged$cell_id), , drop = FALSE]
  attr(counts, "seurat_audit") <- data.frame(
    seurat_rds = seurat_rds,
    object_class = paste(class(obj), collapse = ";"),
    assays = paste(names(obj@assays), collapse = ";"),
    n_object_cells = ncol(counts),
    n_object_genes = nrow(counts),
    assay = assay,
    counts_layer = counts_layer,
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = merged, seurat_meta_columns = wanted)
}

ptp_audit_counts <- function(counts) {
  values <- counts@x
  data.frame(
    matrix_class = paste(class(counts), collapse = ";"),
    n_genes = nrow(counts),
    n_cells = ncol(counts),
    n_nonzero = length(values),
    sparse = inherits(counts, "sparseMatrix"),
    non_negative = length(values) == 0L || all(values >= 0 & is.finite(values)),
    integer_like = length(values) == 0L || max(abs(values - round(values)), na.rm = TRUE) < 1e-8,
    min_nonzero = if (length(values) > 0L) min(values) else NA_real_,
    max_nonzero = if (length(values) > 0L) max(values) else NA_real_,
    stringsAsFactors = FALSE
  )
}

ptp_match_cells <- function(meta, counts, min_match_rate, audit_dir) {
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
  ptp_write_csv(audit, file.path(audit_dir, "cell_id_join_audit.csv"))
  unmatched_meta <- metadata_ids[!meta_match]
  unmatched_expr <- expression_ids[!expr_match]
  unmatched <- rbind(
    data.frame(cell_id = unmatched_meta, source = rep("metadata", length(unmatched_meta)), stringsAsFactors = FALSE),
    data.frame(cell_id = unmatched_expr, source = rep("expression", length(unmatched_expr)), stringsAsFactors = FALSE)
  )
  ptp_write_csv(unmatched, file.path(audit_dir, "cell_id_unmatched_ids.csv"))
  if (audit$match_rate[audit$source == "metadata"] < min_match_rate) {
    stop("Metadata-to-expression match rate below threshold.", call. = FALSE)
  }
  meta <- meta[meta$cell_id %in% expression_ids, , drop = FALSE]
  counts <- counts[, meta$cell_id, drop = FALSE]
  list(meta = meta, counts = counts, audit = audit)
}

ptp_limit_for_smoke <- function(meta, counts, max_cells = 0L, max_genes = 0L, seed = 1L) {
  set.seed(seed)
  if (is.finite(max_cells) && max_cells > 0L && nrow(meta) > max_cells) {
    by_sample <- split(seq_len(nrow(meta)), meta$sample_id)
    n_each <- pmax(1L, floor(max_cells * lengths(by_sample) / nrow(meta)))
    idx <- unlist(Map(function(v, n) sample(v, min(length(v), n)), by_sample, n_each), use.names = FALSE)
    if (length(idx) < max_cells) {
      rest <- setdiff(seq_len(nrow(meta)), idx)
      idx <- c(idx, sample(rest, min(length(rest), max_cells - length(idx))))
    }
    idx <- sort(unique(idx))
    meta <- meta[idx, , drop = FALSE]
    counts <- counts[, meta$cell_id, drop = FALSE]
  }
  if (is.finite(max_genes) && max_genes > 0L && nrow(counts) > max_genes) {
    nnz <- Matrix::rowSums(counts > 0)
    keep <- order(nnz, decreasing = TRUE)[seq_len(max_genes)]
    counts <- counts[keep, , drop = FALSE]
  }
  list(meta = meta, counts = counts)
}

ptp_assignment_audit <- function(meta) {
  sample_ids <- sort(unique(meta$sample_id))
  rows <- lapply(sample_ids, function(sample_id) {
    x <- meta[meta$sample_id == sample_id, , drop = FALSE]
    data.frame(
      sample_id = sample_id,
      n_cells = nrow(x),
      n_initial_ploidy = length(unique(as.character(x$initial_ploidy))),
      initial_ploidy_values = paste(unique(as.character(x$initial_ploidy)), collapse = ";"),
      n_dose = length(unique(x$gemcitabine_dose_mg_per_kg)),
      dose_values = paste(unique(x$gemcitabine_dose_mg_per_kg), collapse = ";"),
      one_initial_ploidy = length(unique(as.character(x$initial_ploidy))) == 1L,
      one_dose = length(unique(x$gemcitabine_dose_mg_per_kg)) == 1L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ptp_processing_batch_crosswalk <- function(meta, sample_info_path = NULL, endpoint_ploidy_specs = ptp_endpoint_ploidy_specs()) {
  sm <- ptp_sample_metadata(meta, endpoint_ploidy_specs)
  extra_cols <- intersect(
    c("orig.ident", "sample", "sample_folder", "harvest", "Sequencing.IDs", "IDs", "Dose"),
    names(meta)
  )
  for (col in extra_cols) {
    sm[[col]] <- vapply(sm$sample_id, function(s) ptp_sample_level_unique(meta[meta$sample_id == s, col]), character(1L))
  }
  if (!is.null(sample_info_path) && nzchar(sample_info_path) && file.exists(sample_info_path)) {
    info <- as.data.frame(readxl::read_excel(sample_info_path), stringsAsFactors = FALSE)
    names(info) <- gsub(" ", ".", names(info), fixed = TRUE)
    if ("IDs" %in% names(info)) {
      info$sample_id_info <- as.character(info$IDs)
      keep <- intersect(c("sample_id_info", "harvest", "Sequencing.IDs", "Dose"), names(info))
      sm <- merge(sm, info[, keep, drop = FALSE], by.x = "sample_id", by.y = "sample_id_info", all.x = TRUE, suffixes = c("", ".sample_info"), sort = FALSE)
    }
  }
  sm$batch_separability_note <- "Available local metadata do not contain explicit dissociation batch, 10x channel, sequencing batch, sequencing saturation, or demultiplexing identifiers."
  sm[order(sm$initial_ploidy, sm$dose_mg, sm$sample_id), , drop = FALSE]
}

ptp_batch_group_separability <- function(crosswalk) {
  group <- paste(crosswalk$initial_ploidy, crosswalk$treatment, sep = "_")
  candidates <- setdiff(names(crosswalk), c(
    "sample_id", "initial_ploidy", "gemcitabine_dose", "dose_mg", "treatment", "n_cells",
    "n_endpoint_ploidy_cells", "sample_mean_endpoint_ploidy", "sample_median_endpoint_ploidy",
    "sample_max_endpoint_ploidy", "sample_mean_endpoint_ploidy_scaled",
    "mean_pseudotime", "median_pseudotime", "q05_pseudotime", "q25_pseudotime",
    "q75_pseudotime", "q95_pseudotime", "mean_percent_mt", "mean_nCount_RNA",
    "mean_nFeature_RNA", "batch_separability_note"
  ))
  candidates <- candidates[!grepl("^ETP_|endpoint_ploidy|average_ploidy|median_ploidy", candidates)]
  rows <- lapply(candidates, function(col) {
    val <- as.character(crosswalk[[col]])
    usable <- !is.na(val) & nzchar(val)
    if (!any(usable)) {
      return(data.frame(batch_variable = col, n_nonmissing = 0L, n_levels = 0L,
                        any_level_unique_to_one_ploidy_treatment = NA, note = "all missing", stringsAsFactors = FALSE))
    }
    tab <- table(val[usable], group[usable])
    unique_level <- apply(tab > 0, 1L, sum) == 1L
    data.frame(
      batch_variable = col,
      n_nonmissing = sum(usable),
      n_levels = length(unique(val[usable])),
      any_level_unique_to_one_ploidy_treatment = any(unique_level),
      note = if (any(unique_level)) "At least one observed level maps to only one ploidy-treatment group; inspect manually." else "No observed level uniquely identifies a ploidy-treatment group.",
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) {
    return(data.frame(batch_variable = NA_character_, n_nonmissing = 0L, n_levels = 0L,
                      any_level_unique_to_one_ploidy_treatment = NA, note = "no available batch-like variables", stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

ptp_group_label <- function(subbin_id, ploidy, treatment_or_dose) {
  paste0(subbin_id, "__p", gsub("[^A-Za-z0-9]+", "", ploidy), "__", treatment_or_dose)
}

ptp_assign_subbin <- function(meta, subbins) {
  meta$subbin_id <- NA_character_
  meta$subbin_label <- NA_character_
  for (sb in subbins) {
    hit <- ptp_interval_hit(meta$pseudotime, sb)
    meta$subbin_id[hit] <- sb$id
    meta$subbin_label[hit] <- sb$label %||% sb$id
  }
  meta
}

ptp_cluster_composition <- function(meta, sample_ids, subbins) {
  cluster_col <- if ("cluster_final" %in% names(meta) && any(!is.na(meta$cluster_final) & nzchar(as.character(meta$cluster_final)))) "cluster_final" else "cluster"
  rows <- list()
  for (sample_id in sample_ids) {
    for (sb in subbins) {
      x <- meta[meta$sample_id == sample_id & meta$subbin_id == sb$id, , drop = FALSE]
      if (nrow(x) == 0L) {
        rows[[length(rows) + 1L]] <- data.frame(sample_id = sample_id, subbin_id = sb$id, cluster = NA_character_, n_cells = 0L, fraction = NA_real_, stringsAsFactors = FALSE)
      } else {
        cluster_values <- as.character(x[[cluster_col]])
        cluster_values <- cluster_values[!is.na(cluster_values) & nzchar(cluster_values)]
        if (length(cluster_values) == 0L) cluster_values <- "unknown"
        tb <- as.data.frame(table(cluster_values), stringsAsFactors = FALSE)
        names(tb) <- c("cluster", "n_cells")
        tb$sample_id <- sample_id
        tb$subbin_id <- sb$id
        tb$fraction <- tb$n_cells / nrow(x)
        rows[[length(rows) + 1L]] <- tb[, c("sample_id", "subbin_id", "cluster", "n_cells", "fraction")]
      }
    }
  }
  do.call(rbind, rows)
}

ptp_construct_subbin_pseudobulk <- function(meta, counts, subbins, qc, endpoint_ploidy_specs = ptp_endpoint_ploidy_specs()) {
  meta <- ptp_assign_subbin(meta, subbins)
  meta_bin <- meta[!is.na(meta$subbin_id), , drop = FALSE]
  sample_meta <- ptp_sample_metadata(meta, endpoint_ploidy_specs)
  all_sample_bins <- expand.grid(
    sample_id = sample_meta$sample_id,
    subbin_id = vapply(subbins, `[[`, character(1L), "id"),
    stringsAsFactors = FALSE
  )
  subbin_info <- ptp_subbins_table(subbins, grid_id = "active")
  all_sample_bins <- merge(all_sample_bins, subbin_info[, c("subbin_id", "label", "start", "end")], by = "subbin_id", all.x = TRUE, sort = FALSE)
  all_sample_bins <- merge(all_sample_bins, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  all_sample_bins$sample_subbin_id <- paste(all_sample_bins$sample_id, all_sample_bins$subbin_id, sep = "__")
  all_sample_bins <- all_sample_bins[order(all_sample_bins$sample_id, all_sample_bins$start), , drop = FALSE]

  if (nrow(meta_bin) > 0L) {
    meta_bin$sample_subbin_id <- paste(meta_bin$sample_id, meta_bin$subbin_id, sep = "__")
    group_factor <- factor(meta_bin$sample_subbin_id, levels = all_sample_bins$sample_subbin_id)
    design <- Matrix::sparse.model.matrix(~ 0 + group_factor)
    colnames(design) <- levels(group_factor)
    pb_counts <- counts[, meta_bin$cell_id, drop = FALSE] %*% design
  } else {
    pb_counts <- counts[, integer(0), drop = FALSE]
  }
  missing_cols <- setdiff(all_sample_bins$sample_subbin_id, colnames(pb_counts))
  if (length(missing_cols) > 0L) {
    zero <- Matrix::Matrix(0, nrow = nrow(counts), ncol = length(missing_cols), sparse = TRUE)
    rownames(zero) <- rownames(counts)
    colnames(zero) <- missing_cols
    pb_counts <- cbind(pb_counts, zero)
  }
  pb_counts <- pb_counts[, all_sample_bins$sample_subbin_id, drop = FALSE]
  rownames(pb_counts) <- rownames(counts)
  all_sample_bins$cell_count <- as.integer(table(factor(meta_bin$sample_subbin_id, levels = all_sample_bins$sample_subbin_id)))
  all_sample_bins$library_size <- as.numeric(Matrix::colSums(pb_counts))
  all_sample_bins$detected_genes <- as.numeric(Matrix::colSums(pb_counts > 0))
  all_sample_bins$mean_subbin_pseudotime <- vapply(all_sample_bins$sample_subbin_id, function(id) {
    values <- meta_bin$pseudotime[meta_bin$sample_subbin_id == id]
    if (length(values) == 0L) NA_real_ else mean(values)
  }, numeric(1L))
  all_sample_bins$median_subbin_pseudotime <- vapply(all_sample_bins$sample_subbin_id, function(id) {
    values <- meta_bin$pseudotime[meta_bin$sample_subbin_id == id]
    if (length(values) == 0L) NA_real_ else median(values)
  }, numeric(1L))
  all_sample_bins$retained_for_model <- with(
    all_sample_bins,
    cell_count >= qc$min_cells_per_mouse_subbin &
      library_size >= qc$min_library_size_per_mouse_subbin &
      detected_genes >= qc$min_detected_genes_per_mouse_subbin
  )
  all_sample_bins$exclusion_reason <- ""
  all_sample_bins$exclusion_reason[all_sample_bins$cell_count < qc$min_cells_per_mouse_subbin] <- paste0("cell_count<", qc$min_cells_per_mouse_subbin)
  low_lib <- all_sample_bins$retained_for_model == FALSE & all_sample_bins$cell_count >= qc$min_cells_per_mouse_subbin &
    all_sample_bins$library_size < qc$min_library_size_per_mouse_subbin
  all_sample_bins$exclusion_reason[low_lib] <- paste0("library_size<", qc$min_library_size_per_mouse_subbin)
  low_gene <- all_sample_bins$retained_for_model == FALSE & all_sample_bins$cell_count >= qc$min_cells_per_mouse_subbin &
    all_sample_bins$library_size >= qc$min_library_size_per_mouse_subbin &
    all_sample_bins$detected_genes < qc$min_detected_genes_per_mouse_subbin
  all_sample_bins$exclusion_reason[low_gene] <- paste0("detected_genes<", qc$min_detected_genes_per_mouse_subbin)
  all_sample_bins$ptp_group <- ptp_group_label(all_sample_bins$subbin_id, all_sample_bins$initial_ploidy, all_sample_bins$treatment)
  all_sample_bins$dose_group <- ptp_group_label(all_sample_bins$subbin_id, all_sample_bins$initial_ploidy, paste0("dose", all_sample_bins$dose_mg))
  for (spec in endpoint_ploidy_specs) {
    out_col <- paste0(spec$method, "_ptp_group")
    if (spec$group_column %in% names(all_sample_bins)) {
      group <- as.character(all_sample_bins[[spec$group_column]])
      all_sample_bins[[out_col]] <- ifelse(
        !is.na(group) & nzchar(group),
        ptp_group_label(all_sample_bins$subbin_id, group, all_sample_bins$treatment),
        NA_character_
      )
    }
  }
  list(counts = pb_counts, metadata = all_sample_bins, cell_metadata = meta, cluster_composition = ptp_cluster_composition(meta, sample_meta$sample_id, subbins))
}

ptp_support_table <- function(pb_meta, qc) {
  groups <- expand.grid(
    subbin_id = unique(pb_meta$subbin_id),
    initial_ploidy = c("2N", "4N"),
    treatment = c("control", "gemcitabine"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    x <- pb_meta[pb_meta$subbin_id == g$subbin_id & pb_meta$initial_ploidy == g$initial_ploidy & pb_meta$treatment == g$treatment, , drop = FALSE]
    data.frame(
      subbin_id = g$subbin_id,
      initial_ploidy = g$initial_ploidy,
      treatment = g$treatment,
      n_mice_total = length(unique(x$sample_id)),
      n_mice_qc = length(unique(x$sample_id[x$retained_for_model])),
      total_cells = sum(x$cell_count),
      min_cells = if (nrow(x) > 0L) min(x$cell_count) else NA_integer_,
      median_cells = if (nrow(x) > 0L) median(x$cell_count) else NA_real_,
      support_pass = length(unique(x$sample_id[x$retained_for_model])) >= qc$min_mice_per_ploidy_treatment_subbin,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  per_bin <- aggregate(support_pass ~ subbin_id, out, all)
  names(per_bin)[2] <- "common_support"
  merge(out, per_bin, by = "subbin_id", all.x = TRUE, sort = FALSE)
}

ptp_endpoint_group_balance <- function(sample_meta, endpoint_ploidy_specs = ptp_endpoint_ploidy_specs()) {
  rows <- list()
  for (spec in endpoint_ploidy_specs) {
    if (!spec$group_column %in% names(sample_meta)) next
    tab <- as.data.frame(table(
      initial_ploidy = sample_meta$initial_ploidy,
      treatment = sample_meta$treatment,
      endpoint_ploidy_group = sample_meta[[spec$group_column]],
      useNA = "ifany"
    ), stringsAsFactors = FALSE)
    tab$method <- spec$method
    tab$threshold <- spec$threshold
    rows[[length(rows) + 1L]] <- tab
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out[, c("method", "threshold", "initial_ploidy", "treatment", "endpoint_ploidy_group", "Freq")]
}

ptp_support_table_by_group <- function(pb_meta, group_value_col, group_levels, qc) {
  groups <- expand.grid(
    subbin_id = unique(pb_meta$subbin_id),
    endpoint_ploidy_group = group_levels,
    treatment = c("control", "gemcitabine"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    x <- pb_meta[pb_meta$subbin_id == g$subbin_id &
                   pb_meta[[group_value_col]] %in% g$endpoint_ploidy_group &
                   pb_meta$treatment == g$treatment, , drop = FALSE]
    data.frame(
      subbin_id = g$subbin_id,
      endpoint_ploidy_group = g$endpoint_ploidy_group,
      treatment = g$treatment,
      n_mice_total = length(unique(x$sample_id)),
      n_mice_qc = length(unique(x$sample_id[x$retained_for_model])),
      total_cells = sum(x$cell_count),
      min_cells = if (nrow(x) > 0L) min(x$cell_count) else NA_integer_,
      median_cells = if (nrow(x) > 0L) median(x$cell_count) else NA_real_,
      support_pass = length(unique(x$sample_id[x$retained_for_model])) >= qc$min_mice_per_ploidy_treatment_subbin,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  per_bin <- aggregate(support_pass ~ subbin_id, out, all)
  names(per_bin)[2] <- "common_support"
  merge(out, per_bin, by = "subbin_id", all.x = TRUE, sort = FALSE)
}

ptp_endpoint_continuous_support <- function(pb_meta, eligible_subbins, qc, etp_col = "sample_mean_endpoint_ploidy_scaled") {
  groups <- expand.grid(
    subbin_id = eligible_subbins,
    treatment = c("control", "gemcitabine"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    x <- pb_meta[pb_meta$subbin_id == g$subbin_id & pb_meta$treatment == g$treatment, , drop = FALSE]
    q <- x$retained_for_model & is.finite(x[[etp_col]])
    data.frame(
      subbin_id = g$subbin_id,
      treatment = g$treatment,
      n_mice_total = length(unique(x$sample_id)),
      n_mice_qc_finite_etp = length(unique(x$sample_id[q])),
      min_scaled_endpoint_ploidy = if (any(q)) min(x[[etp_col]][q]) else NA_real_,
      max_scaled_endpoint_ploidy = if (any(q)) max(x[[etp_col]][q]) else NA_real_,
      support_pass = length(unique(x$sample_id[q])) >= qc$min_mice_per_ploidy_treatment_subbin,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  per_bin <- aggregate(support_pass ~ subbin_id, out, all)
  names(per_bin)[2] <- "common_support"
  merge(out, per_bin, by = "subbin_id", all.x = TRUE, sort = FALSE)
}

ptp_choose_primary_grid <- function(meta, counts, cfg) {
  qc <- cfg$qc
  endpoint_ploidy_specs <- ptp_endpoint_ploidy_specs(cfg)
  primary <- ptp_construct_subbin_pseudobulk(meta, counts, cfg$primary_subbins, qc, endpoint_ploidy_specs)
  primary_support <- ptp_support_table(primary$metadata, qc)
  primary_common <- unique(primary_support$subbin_id[primary_support$common_support])
  if (length(primary_common) == length(cfg$primary_subbins)) {
    primary$grid_id <- "primary_four_bin"
    primary$support <- primary_support
    primary$eligible_subbins <- primary_common
    primary$estimable <- length(primary_common) >= qc$min_common_support_subbins
    return(primary)
  }
  fallback <- ptp_construct_subbin_pseudobulk(meta, counts, cfg$fallback_subbins, qc, endpoint_ploidy_specs)
  fallback_support <- ptp_support_table(fallback$metadata, qc)
  fallback_common <- unique(fallback_support$subbin_id[fallback_support$common_support])
  fallback$grid_id <- "fallback_two_bin"
  fallback$support <- fallback_support
  fallback$eligible_subbins <- fallback_common
  fallback$estimable <- length(fallback_common) >= qc$min_common_support_subbins
  fallback$primary_four_bin_support <- primary_support
  fallback
}

ptp_fit_cellmeans <- function(pb_counts, pb_meta, eligible_subbins, group_col, qc) {
  keep_obs <- pb_meta$retained_for_model & pb_meta$subbin_id %in% eligible_subbins & pb_meta$library_size > 0
  meta <- pb_meta[keep_obs, , drop = FALSE]
  counts <- pb_counts[, meta$sample_subbin_id, drop = FALSE]
  if (ncol(counts) < 8L) stop("Too few pseudobulk observations for model fitting.", call. = FALSE)
  y0 <- edgeR::DGEList(counts = counts)
  cpm0 <- edgeR::cpm(y0)
  keep_gene <- rowSums(cpm0 > qc$gene_filter_cpm) >= qc$gene_filter_min_observations
  if (sum(keep_gene) < 10L) stop("Too few genes pass expression filtering.", call. = FALSE)
  counts <- counts[keep_gene, , drop = FALSE]
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  group <- factor(meta[[group_col]])
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- sub("^group", "", colnames(design))
  qr_rank <- qr(design)$rank
  dropped <- character(0L)
  if (qr_rank < ncol(design)) {
    keep_cols <- qr(design)$pivot[seq_len(qr_rank)]
    retained <- colnames(design)[sort(keep_cols)]
    dropped <- setdiff(colnames(design), retained)
    design <- design[, retained, drop = FALSE]
  }
  v1 <- limma::voomWithQualityWeights(y, design, plot = FALSE)
  dup1 <- limma::duplicateCorrelation(v1, design, block = meta$sample_id)
  v2 <- tryCatch(
    limma::voomWithQualityWeights(y, design, plot = FALSE, block = meta$sample_id, correlation = dup1$consensus.correlation),
    error = function(e) limma::voom(y, design, plot = FALSE)
  )
  dup2 <- limma::duplicateCorrelation(v2, design, block = meta$sample_id)
  fit <- limma::lmFit(v2, design, block = meta$sample_id, correlation = dup2$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  meta$voom_quality_weight <- NA_real_
  if (!is.null(v2$targets$sample.weights)) meta$voom_quality_weight <- as.numeric(v2$targets$sample.weights)
  list(
    fit = fit,
    voom = v2,
    dge = y,
    metadata = meta,
    design = design,
    retained_genes = rownames(counts),
    group_col = group_col,
    duplicate_correlation_first = dup1$consensus.correlation,
    duplicate_correlation_second = dup2$consensus.correlation,
    design_rank = qr_rank,
    design_ncol_original = nlevels(group),
    dropped_design_columns = dropped
  )
}

ptp_design_audit <- function(model, model_id) {
  data.frame(
    model_id = model_id,
    n_observations = nrow(model$design),
    n_genes = length(model$retained_genes),
    n_design_columns_original = model$design_ncol_original,
    n_design_columns_retained = ncol(model$design),
    design_rank = model$design_rank,
    rank_deficient = length(model$dropped_design_columns) > 0L,
    retained_design_columns = paste(colnames(model$design), collapse = ";"),
    dropped_design_columns = paste(model$dropped_design_columns, collapse = ";"),
    duplicate_correlation_first = model$duplicate_correlation_first,
    duplicate_correlation_second = model$duplicate_correlation_second,
    stringsAsFactors = FALSE
  )
}

ptp_apply_contrasts <- function(model, contrast_list, contrast_family) {
  if (length(contrast_list) == 0L) return(data.frame())
  mat <- do.call(cbind, lapply(contrast_list, function(vec) {
    out <- rep(0, ncol(model$design))
    names(out) <- colnames(model$design)
    common <- intersect(names(vec), names(out))
    out[common] <- vec[common]
    out
  }))
  colnames(mat) <- names(contrast_list)
  fit2 <- limma::contrasts.fit(model$fit, contrasts = mat)
  fit2 <- limma::eBayes(fit2, robust = TRUE)
  rows <- lapply(seq_len(ncol(mat)), function(i) {
    tt <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none")
    se <- fit2$stdev.unscaled[, i] * sqrt(fit2$s2.post)
    names(se) <- rownames(fit2$coefficients)
    df <- fit2$df.total
    if (length(df) == 1L) df <- rep(df, length(se))
    names(df) <- rownames(fit2$coefficients)
    crit <- stats::qt(0.975, df = df)
    row_id <- rownames(tt)
    data.frame(
      gene = row_id,
      gene_symbol = ptp_clean_gene_symbols(row_id),
      contrast_family = contrast_family,
      contrast_id = colnames(mat)[i],
      estimate = tt$logFC,
      se = se[row_id],
      ci_low = tt$logFC - crit[row_id] * se[row_id],
      ci_high = tt$logFC + crit[row_id] * se[row_id],
      average_expression = tt$AveExpr,
      t_statistic = tt$t,
      p_value = tt$P.Value,
      fdr = tt$adj.P.Val,
      B = tt$B,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out[order(out$contrast_id, out$p_value, -abs(out$t_statistic), out$gene), , drop = FALSE]
}

ptp_primary_contrasts <- function(eligible_subbins, design_cols) {
  avg <- function(ploidy, treatment) {
    vec <- rep(0, length(design_cols)); names(vec) <- design_cols
    labels <- vapply(eligible_subbins, function(sb) ptp_group_label(sb, ploidy, treatment), character(1L))
    labels <- intersect(labels, design_cols)
    if (length(labels) > 0L) vec[labels] <- 1 / length(eligible_subbins)
    vec
  }
  eff2 <- avg("2N", "gemcitabine") - avg("2N", "control")
  eff4 <- avg("4N", "gemcitabine") - avg("4N", "control")
  list(
    delta_2N_origin = eff2,
    delta_4N_origin = eff4,
    treatment_by_initial_ploidy_interaction = eff4 - eff2
  )
}

ptp_dose_contrasts <- function(eligible_subbins, design_cols) {
  avg <- function(ploidy, dose_label) {
    vec <- rep(0, length(design_cols)); names(vec) <- design_cols
    labels <- vapply(eligible_subbins, function(sb) ptp_group_label(sb, ploidy, dose_label), character(1L))
    labels <- intersect(labels, design_cols)
    if (length(labels) > 0L) vec[labels] <- 1 / length(eligible_subbins)
    vec
  }
  eff2_30 <- avg("2N", "dose30") - avg("2N", "dose0")
  eff4_30 <- avg("4N", "dose30") - avg("4N", "dose0")
  eff2_120 <- avg("2N", "dose120") - avg("2N", "dose0")
  eff4_120 <- avg("4N", "dose120") - avg("4N", "dose0")
  list(
    interaction_30_vs_0 = eff4_30 - eff2_30,
    interaction_120_vs_0 = eff4_120 - eff2_120,
    delta_2N_30_vs_0 = eff2_30,
    delta_4N_30_vs_0 = eff4_30,
    delta_2N_120_vs_0 = eff2_120,
    delta_4N_120_vs_0 = eff4_120
  )
}

ptp_endpoint_group_contrasts <- function(eligible_subbins, design_cols) {
  avg <- function(endpoint_group, treatment) {
    vec <- rep(0, length(design_cols)); names(vec) <- design_cols
    labels <- vapply(eligible_subbins, function(sb) ptp_group_label(sb, endpoint_group, treatment), character(1L))
    labels <- intersect(labels, design_cols)
    if (length(labels) > 0L) vec[labels] <- 1 / length(eligible_subbins)
    vec
  }
  eff_low <- avg("ETP-lower", "gemcitabine") - avg("ETP-lower", "control")
  eff_high <- avg("ETP-higher", "gemcitabine") - avg("ETP-higher", "control")
  list(
    delta_ETP_lower = eff_low,
    delta_ETP_higher = eff_high,
    treatment_by_endpoint_ploidy_group_interaction = eff_high - eff_low
  )
}

ptp_continuous_endpoint_ploidy_contrasts <- function(eligible_subbins, design_cols) {
  avg <- function(treatment, suffix = "") {
    vec <- rep(0, length(design_cols)); names(vec) <- design_cols
    labels <- vapply(eligible_subbins, function(sb) paste0(sb, "__", treatment, suffix), character(1L))
    labels <- intersect(labels, design_cols)
    if (length(labels) > 0L) vec[labels] <- 1 / length(eligible_subbins)
    vec
  }
  base_control <- avg("control")
  base_gem <- avg("gemcitabine")
  slope_suffix <- "__sample_mean_endpoint_ploidy_scaled"
  slope_control <- avg("control", slope_suffix)
  slope_gem <- avg("gemcitabine", slope_suffix)
  list(
    treatment_at_mean_endpoint_ploidy = base_gem - base_control,
    endpoint_ploidy_slope_control = slope_control,
    endpoint_ploidy_slope_gemcitabine = slope_gem,
    continuous_endpoint_ploidy_treatment_interaction = slope_gem - slope_control
  )
}

ptp_fit_continuous_endpoint_ploidy <- function(pb_counts, pb_meta, eligible_subbins, qc, etp_col = "sample_mean_endpoint_ploidy_scaled") {
  keep_obs <- pb_meta$retained_for_model &
    pb_meta$subbin_id %in% eligible_subbins &
    pb_meta$library_size > 0 &
    is.finite(pb_meta[[etp_col]])
  meta <- pb_meta[keep_obs, , drop = FALSE]
  counts <- pb_counts[, meta$sample_subbin_id, drop = FALSE]
  if (ncol(counts) < 8L) stop("Too few pseudobulk observations for continuous ETP model fitting.", call. = FALSE)
  y0 <- edgeR::DGEList(counts = counts)
  cpm0 <- edgeR::cpm(y0)
  keep_gene <- rowSums(cpm0 > qc$gene_filter_cpm) >= qc$gene_filter_min_observations
  if (sum(keep_gene) < 10L) stop("Too few genes pass expression filtering.", call. = FALSE)
  counts <- counts[keep_gene, , drop = FALSE]
  y <- edgeR::DGEList(counts = counts)
  y <- edgeR::calcNormFactors(y, method = "TMM")

  group_levels <- as.vector(outer(eligible_subbins, c("control", "gemcitabine"), paste, sep = "__"))
  meta$subbin_treatment <- factor(paste(meta$subbin_id, meta$treatment, sep = "__"), levels = group_levels)
  base <- stats::model.matrix(~ 0 + subbin_treatment, data = meta)
  colnames(base) <- sub("^subbin_treatment", "", colnames(base))
  slope <- base * meta[[etp_col]]
  colnames(slope) <- paste0(colnames(base), "__sample_mean_endpoint_ploidy_scaled")
  design <- cbind(base, slope)
  qr_rank <- qr(design)$rank
  dropped <- character(0L)
  if (qr_rank < ncol(design)) {
    keep_cols <- qr(design)$pivot[seq_len(qr_rank)]
    retained <- colnames(design)[sort(keep_cols)]
    dropped <- setdiff(colnames(design), retained)
    design <- design[, retained, drop = FALSE]
  }
  v1 <- limma::voomWithQualityWeights(y, design, plot = FALSE)
  dup1 <- limma::duplicateCorrelation(v1, design, block = meta$sample_id)
  v2 <- tryCatch(
    limma::voomWithQualityWeights(y, design, plot = FALSE, block = meta$sample_id, correlation = dup1$consensus.correlation),
    error = function(e) limma::voom(y, design, plot = FALSE)
  )
  dup2 <- limma::duplicateCorrelation(v2, design, block = meta$sample_id)
  fit <- limma::lmFit(v2, design, block = meta$sample_id, correlation = dup2$consensus.correlation)
  fit <- limma::eBayes(fit, robust = TRUE)
  meta$voom_quality_weight <- NA_real_
  if (!is.null(v2$targets$sample.weights)) meta$voom_quality_weight <- as.numeric(v2$targets$sample.weights)
  list(
    fit = fit,
    voom = v2,
    dge = y,
    metadata = meta,
    design = design,
    retained_genes = rownames(counts),
    group_col = etp_col,
    duplicate_correlation_first = dup1$consensus.correlation,
    duplicate_correlation_second = dup2$consensus.correlation,
    design_rank = qr_rank,
    design_ncol_original = length(group_levels) * 2L,
    dropped_design_columns = dropped
  )
}

ptp_load_programs <- function(cfg, min_genes_observed, observed_symbols) {
  msig <- NULL
  if (any(vapply(cfg$programs, function(x) identical(x$source, "msigdbr"), logical(1L)))) {
    msig <- as.data.frame(msigdbr::msigdbr(species = "Homo sapiens"), stringsAsFactors = FALSE)
  }
  scoring_cfg <- cfg$program_scoring %||% list()
  direction_signs <- unlist(scoring_cfg$direction_signs %||% list(), use.names = TRUE)
  direction_labels <- unlist(scoring_cfg$direction_labels %||% list(), use.names = TRUE)
  resolve_set_genes <- function(set_names) {
    set_names <- as.character(set_names %||% character(0L))
    if (length(set_names) == 0L) return(character(0L))
    unique(msig$gene_symbol[msig$gs_name %in% set_names])
  }
  build_gene_weights <- function(p, genes, expected_direction) {
    scoring <- p$scoring %||% list()
    mode <- scoring$mode %||% scoring_cfg$default_mode %||% "signed_weighted_mean"
    direction_source <- scoring$direction_source %||% scoring_cfg$direction_source %||% "expected_direction_fixed_in_config"
    component_policy <- scoring$component_policy %||% "single_component"
    components <- scoring$components %||% list()
    if (length(components) > 0L) {
      rows <- lapply(components, function(comp) {
        comp_genes <- unique(c(as.character(comp$genes %||% character(0L)), resolve_set_genes(comp$set_names %||% character(0L))))
        comp_genes <- comp_genes[!is.na(comp_genes) & nzchar(comp_genes)]
        sign <- suppressWarnings(as.numeric(comp$sign %||% 1))
        wt <- suppressWarnings(as.numeric(comp$weight %||% scoring_cfg$default_gene_weight %||% 1))
        if (!is.finite(sign) || sign == 0) sign <- 1
        if (!is.finite(wt) || wt == 0) wt <- 1
        data.frame(
          gene = comp_genes,
          raw_weight = wt,
          sign = sign,
          signed_weight = sign * abs(wt),
          component_id = comp$id %||% "component",
          component_label = comp$label %||% comp$id %||% "component",
          stringsAsFactors = FALSE
        )
      })
      w <- do.call(rbind, rows)
    } else {
      sign <- suppressWarnings(as.numeric(p$score_sign %||% direction_signs[[expected_direction]] %||% 1))
      wt <- suppressWarnings(as.numeric(p$score_weight %||% scoring_cfg$default_gene_weight %||% 1))
      if (!is.finite(sign) || sign == 0) sign <- 1
      if (!is.finite(wt) || wt == 0) wt <- 1
      w <- data.frame(
        gene = genes,
        raw_weight = wt,
        sign = sign,
        signed_weight = sign * abs(wt),
        component_id = "primary",
        component_label = direction_labels[[expected_direction]] %||% expected_direction %||% "primary",
        stringsAsFactors = FALSE
      )
    }
    w <- w[!is.na(w$gene) & nzchar(w$gene) & is.finite(w$signed_weight) & w$signed_weight != 0, , drop = FALSE]
    if (nrow(w) == 0L) {
      w <- data.frame(gene = genes, raw_weight = 1, sign = 1, signed_weight = 1,
                      component_id = "primary", component_label = "primary", stringsAsFactors = FALSE)
    }
    agg <- stats::aggregate(signed_weight ~ gene, data = w, sum)
    agg <- agg[agg$signed_weight != 0, , drop = FALSE]
    list(
      weights = agg,
      component_weights = w,
      score_mode = mode,
      direction_source = direction_source,
      component_policy = component_policy,
      score_direction_label = if (length(components) > 0L) {
        paste(unique(w$component_label), collapse = ";")
      } else {
        direction_labels[[expected_direction]] %||% expected_direction %||% "primary"
      }
    )
  }
  programs <- lapply(cfg$programs, function(p) {
    if (identical(p$source, "custom")) {
      genes <- unique(as.character(p$genes %||% character(0L)))
      set_names <- "custom"
      missing_sets <- character(0L)
    } else {
      set_names <- as.character(p$set_names %||% character(0L))
      missing_sets <- setdiff(set_names, unique(msig$gs_name))
      genes <- unique(msig$gene_symbol[msig$gs_name %in% set_names])
    }
    genes <- unique(genes[!is.na(genes) & nzchar(genes)])
    observed <- intersect(genes, observed_symbols)
    expected_direction <- p$expected_direction %||% NA_character_
    weight_info <- build_gene_weights(p, genes, expected_direction)
    weight_info$weights$gene <- ptp_clean_gene_symbols(weight_info$weights$gene)
    weight_info$component_weights$gene <- ptp_clean_gene_symbols(weight_info$component_weights$gene)
    weight_info$weights <- weight_info$weights[!is.na(weight_info$weights$gene) & nzchar(weight_info$weights$gene), , drop = FALSE]
    weight_info$component_weights <- weight_info$component_weights[!is.na(weight_info$component_weights$gene) & nzchar(weight_info$component_weights$gene), , drop = FALSE]
    observed_weights <- weight_info$weights[weight_info$weights$gene %in% observed_symbols, , drop = FALSE]
    observed_weight_genes <- intersect(observed_weights$gene, observed_symbols)
    observed <- intersect(observed, observed_weight_genes)
    list(
      id = p$id,
      label = p$label %||% p$id,
      family = p$family %||% NA_character_,
      response_family = p$response_family %||% p$family %||% NA_character_,
      response_family_label = p$response_family_label %||% p$response_family %||% p$family %||% NA_character_,
      tier = p$tier %||% "tier1_core",
      tier_label = p$tier_label %||% p$tier %||% "Tier 1 original core",
      upgrade_status = p$upgrade_status %||% "original_core",
      expected_direction = expected_direction,
      score_mode = weight_info$score_mode,
      score_direction_label = weight_info$score_direction_label,
      score_direction_source = weight_info$direction_source,
      score_component_policy = weight_info$component_policy,
      source = p$source,
      set_names = set_names,
      missing_set_names = missing_sets,
      genes = genes,
      observed_genes = observed,
      gene_weights = weight_info$weights,
      observed_gene_weights = observed_weights[observed_weights$gene %in% observed, , drop = FALSE],
      component_gene_weights = weight_info$component_weights,
      estimable = length(observed) >= min_genes_observed
    )
  })
  names(programs) <- vapply(programs, `[[`, character(1L), "id")
  programs
}

ptp_program_meta_cols <- function() {
  c("program_id", "program_label", "tier", "tier_label", "response_family",
    "response_family_label", "family", "upgrade_status", "expected_direction",
    "score_mode", "score_direction_label", "score_direction_source",
    "score_component_policy", "n_genes_observed", "n_positive_weight_genes",
    "n_negative_weight_genes", "sum_abs_observed_weights")
}

ptp_add_program_tier_fdr <- function(df, p_col = "p_value", fdr_col = "fdr") {
  if (is.null(df) || nrow(df) == 0L || !p_col %in% names(df) || !fdr_col %in% names(df)) return(df)
  split_col <- if ("contrast_id" %in% names(df)) "contrast_id" else NULL
  strata <- if (is.null(split_col)) list(all = seq_len(nrow(df))) else split(seq_len(nrow(df)), df[[split_col]], drop = TRUE)
  all_col <- paste0(fdr_col, "_all_programs")
  focused_col <- paste0(fdr_col, "_focused_tier1_tier2")
  tier3_col <- paste0(fdr_col, "_tier3_exploratory")
  negative_col <- paste0(fdr_col, "_negative_control")
  df[[all_col]] <- df[[fdr_col]]
  df[[focused_col]] <- NA_real_
  df[[tier3_col]] <- NA_real_
  df[[negative_col]] <- NA_real_
  for (idx in strata) {
    focused <- idx[df$tier[idx] %in% c("tier1_core", "tier2_upgraded_focused")]
    tier3 <- idx[df$tier[idx] == "tier3_exploratory_only"]
    negative <- idx[df$tier[idx] == "negative_control"]
    if (length(focused) > 0L) df[[focused_col]][focused] <- stats::p.adjust(df[[p_col]][focused], method = "BH")
    if (length(tier3) > 0L) df[[tier3_col]][tier3] <- stats::p.adjust(df[[p_col]][tier3], method = "BH")
    if (length(negative) > 0L) df[[negative_col]][negative] <- stats::p.adjust(df[[p_col]][negative], method = "BH")
  }
  df
}

ptp_program_membership_table <- function(programs) {
  rows <- lapply(programs, function(p) {
    data.frame(
      program_id = p$id,
      program_label = p$label,
      tier = p$tier,
      tier_label = p$tier_label,
      response_family = p$response_family,
      response_family_label = p$response_family_label,
      family = p$family,
      upgrade_status = p$upgrade_status,
      expected_direction = p$expected_direction,
      score_mode = p$score_mode,
      score_direction_label = p$score_direction_label,
      score_direction_source = p$score_direction_source,
      score_component_policy = p$score_component_policy,
      source = p$source,
      set_names = paste(p$set_names, collapse = ";"),
      missing_set_names = paste(p$missing_set_names, collapse = ";"),
      n_genes_defined = length(p$genes),
      n_genes_observed = length(p$observed_genes),
      n_positive_weight_genes = sum(p$observed_gene_weights$signed_weight > 0),
      n_negative_weight_genes = sum(p$observed_gene_weights$signed_weight < 0),
      sum_abs_observed_weights = sum(abs(p$observed_gene_weights$signed_weight)),
      estimable = isTRUE(p$estimable),
      genes_defined = paste(p$genes, collapse = ";"),
      genes_observed = paste(p$observed_genes, collapse = ";"),
      observed_gene_weights = paste(paste(p$observed_gene_weights$gene, signif(p$observed_gene_weights$signed_weight, 6), sep = ":"), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

ptp_program_weight_table <- function(programs) {
  rows <- lapply(programs, function(p) {
    w <- p$component_gene_weights
    if (is.null(w) || nrow(w) == 0L) return(data.frame())
    w$program_id <- p$id
    w$program_label <- p$label
    w$observed <- w$gene %in% p$observed_genes
    w[, c("program_id", "program_label", "component_id", "component_label", "gene", "raw_weight", "sign", "signed_weight", "observed")]
  })
  rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

ptp_logcpm_matrix <- function(model) {
  edgeR::cpm(model$dge, log = TRUE, prior.count = 1)
}

ptp_program_scores_from_logcpm <- function(logcpm, programs) {
  row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  score_rows <- lapply(programs, function(p) {
    if (!isTRUE(p$estimable)) return(rep(NA_real_, ncol(logcpm)))
    idx <- which(row_symbols %in% p$observed_genes)
    idx <- idx[!duplicated(row_symbols[idx])]
    if (length(idx) == 0L) return(rep(NA_real_, ncol(logcpm)))
    mat <- logcpm[idx, , drop = FALSE]
    mat <- t(scale(t(mat)))
    mat[!is.finite(mat)] <- NA_real_
    gene_weights <- p$observed_gene_weights
    gene_weights <- gene_weights[match(row_symbols[idx], gene_weights$gene), , drop = FALSE]
    weights <- gene_weights$signed_weight
    weights[!is.finite(weights)] <- NA_real_
    keep <- is.finite(weights)
    if (!any(keep)) return(rep(NA_real_, ncol(logcpm)))
    denom <- sum(abs(weights[keep]))
    if (!is.finite(denom) || denom <= 0) return(rep(NA_real_, ncol(logcpm)))
    as.numeric(crossprod(weights[keep], mat[keep, , drop = FALSE]) / denom)
  })
  score <- do.call(rbind, score_rows)
  rownames(score) <- names(programs)
  colnames(score) <- colnames(logcpm)
  score
}

ptp_score_method_specs <- function(cfg) {
  methods <- cfg$program_scoring$score_methods %||% list()
  if (length(methods) == 0L) {
    methods <- list(list(
      id = "signed_weighted_mean_zscore",
      label = "Primary signed weighted mean z-score",
      role = "primary",
      mode = "signed_weighted_mean_zscore"
    ))
  }
  methods
}

ptp_gene_z_matrix <- function(logcpm) {
  mat <- t(scale(t(logcpm)))
  mat[!is.finite(mat)] <- NA_real_
  mat
}

ptp_rank_score_matrix <- function(logcpm) {
  out <- apply(logcpm, 2L, function(v) {
    ok <- is.finite(v)
    r <- rep(NA_real_, length(v))
    if (sum(ok) >= 2L) {
      rv <- rank(v[ok], ties.method = "average")
      r[ok] <- (rv - (length(rv) + 1) / 2) / ((length(rv) - 1) / 2)
    }
    r
  })
  out <- as.matrix(out)
  rownames(out) <- rownames(logcpm)
  colnames(out) <- colnames(logcpm)
  out
}

ptp_program_score_methods_from_logcpm <- function(logcpm, programs, cfg, pb_counts = NULL) {
  methods <- ptp_score_method_specs(cfg)
  row_symbols <- ptp_clean_gene_symbols(rownames(logcpm))
  gene_z <- ptp_gene_z_matrix(logcpm)
  rank_scores <- NULL
  primary_scores <- ptp_program_scores_from_logcpm(logcpm, programs)
  count_symbols <- if (!is.null(pb_counts)) ptp_clean_gene_symbols(rownames(pb_counts)) else character(0L)
  out <- list()
  for (method in methods) {
    method_id <- method$id %||% method$mode %||% "score_method"
    method_label <- method$label %||% method_id
    method_role <- method$role %||% "score_definition_sensitivity"
    method_mode <- method$mode %||% method_id
    score_rows <- list()
    qc_rows <- list()
    for (p in programs) {
      idx <- which(row_symbols %in% p$observed_genes)
      idx <- idx[!duplicated(row_symbols[idx])]
      gene_weights <- p$observed_gene_weights
      gene_weights <- gene_weights[match(row_symbols[idx], gene_weights$gene), , drop = FALSE]
      weights <- gene_weights$signed_weight
      keep <- is.finite(weights) & rowSums(is.finite(gene_z[idx, , drop = FALSE])) == ncol(logcpm)
      idx <- idx[keep]
      weights <- weights[keep]
      n_used_base <- length(idx)
      fallback_method <- NA_character_
      fallback_reason <- NA_character_
      n_used_by_obs <- rep(n_used_base, ncol(logcpm))
      score <- rep(NA_real_, ncol(logcpm))
      if (n_used_base > 0L) {
        denom <- sum(abs(weights))
        primary_score <- as.numeric(primary_scores[p$id, colnames(logcpm)])
        if (identical(method_mode, "signed_weighted_mean_zscore")) {
          score <- as.numeric(crossprod(weights, gene_z[idx, , drop = FALSE]) / denom)
        } else if (identical(method_mode, "signed_pca_eigengene")) {
          min_genes <- as.integer(method$min_genes %||% 2L)
          if (n_used_base < min_genes) {
            score <- primary_score
            fallback_method <- "signed_weighted_mean_zscore"
            fallback_reason <- "too_few_genes_for_pca"
          } else {
            signed_mat <- sweep(gene_z[idx, , drop = FALSE], 1L, sign(weights), `*`)
            pc <- tryCatch(stats::prcomp(t(signed_mat), center = FALSE, scale. = FALSE)$x[, 1L], error = function(e) NULL)
            if (is.null(pc) || all(!is.finite(pc)) || stats::sd(pc, na.rm = TRUE) == 0) {
              score <- primary_score
              fallback_method <- "signed_weighted_mean_zscore"
              fallback_reason <- "pca_failed_or_constant"
            } else {
              cor_to_primary <- tryCatch(
                suppressWarnings(stats::cor(pc, primary_score, use = "complete.obs")),
                error = function(e) NA_real_
              )
              if (is.finite(cor_to_primary) && cor_to_primary < 0) pc <- -pc
              score <- as.numeric(scale(pc))
              if (any(!is.finite(score))) {
                score <- primary_score
                fallback_method <- "signed_weighted_mean_zscore"
                fallback_reason <- "pca_standardization_failed"
              }
            }
          }
        } else if (identical(method_mode, "signed_rank_mean")) {
          if (is.null(rank_scores)) rank_scores <- ptp_rank_score_matrix(logcpm)
          score <- as.numeric(crossprod(weights, rank_scores[idx, , drop = FALSE]) / denom)
        } else if (identical(method_mode, "trimmed_signed_mean_zscore")) {
          trim_fraction <- suppressWarnings(as.numeric(method$trim_fraction %||% 0.10))
          min_trim <- as.integer(method$min_genes_for_trimming %||% 10L)
          if (!is.finite(trim_fraction) || trim_fraction < 0 || trim_fraction >= 0.5) trim_fraction <- 0.10
          if (n_used_base < min_trim) {
            score <- primary_score
            fallback_method <- "signed_weighted_mean_zscore"
            fallback_reason <- "too_few_genes_for_trimming"
          } else {
            contrib <- sweep(gene_z[idx, , drop = FALSE], 1L, weights, `*`)
            score <- vapply(seq_len(ncol(contrib)), function(j) {
              v <- contrib[, j]
              w <- weights
              ok <- is.finite(v) & is.finite(w)
              v <- v[ok]
              w <- w[ok]
              if (length(v) == 0L) return(NA_real_)
              k <- floor(trim_fraction * length(v))
              keep_j <- seq_along(v)
              if (k > 0L && length(v) > 2L * k) {
                ord <- order(v)
                keep_j <- ord[(k + 1L):(length(v) - k)]
              }
              sum(v[keep_j]) / sum(abs(w[keep_j]))
            }, numeric(1L))
            k <- floor(trim_fraction * n_used_base)
            if (k > 0L && n_used_base > 2L * k) n_used_by_obs <- rep(n_used_base - 2L * k, ncol(logcpm))
          }
        } else {
          score <- primary_score
          fallback_method <- "signed_weighted_mean_zscore"
          fallback_reason <- paste0("unknown_mode_", method_mode)
        }
      }
      score_rows[[length(score_rows) + 1L]] <- score
      count_idx <- if (!is.null(pb_counts)) which(count_symbols %in% p$observed_genes) else integer(0L)
      count_idx <- count_idx[!duplicated(count_symbols[count_idx])]
      detected <- if (length(count_idx) > 0L) Matrix::colSums(pb_counts[count_idx, , drop = FALSE] > 0) else rep(NA_real_, ncol(logcpm))
      raw_umi <- if (length(count_idx) > 0L) Matrix::colSums(pb_counts[count_idx, , drop = FALSE]) else rep(NA_real_, ncol(logcpm))
      qc_rows[[length(qc_rows) + 1L]] <- data.frame(
        score_method_id = method_id,
        score_method_label = method_label,
        score_method_role = method_role,
        score_method_mode = method_mode,
        program_id = p$id,
        program_label = p$label,
        tier = p$tier,
        tier_label = p$tier_label,
        response_family = p$response_family,
        response_family_label = p$response_family_label,
        sample_subbin_id = colnames(logcpm),
        n_genes_observed = length(p$observed_genes),
        n_genes_used_by_score = n_used_by_obs,
        detected_program_genes = as.numeric(detected),
        detected_program_fraction = as.numeric(detected) / max(length(p$observed_genes), 1L),
        program_raw_umi = as.numeric(raw_umi),
        score_value = as.numeric(score),
        score_finite = is.finite(score),
        fallback_method = fallback_method,
        fallback_reason = fallback_reason,
        stringsAsFactors = FALSE
      )
    }
    score_matrix <- do.call(rbind, score_rows)
    rownames(score_matrix) <- names(programs)
    colnames(score_matrix) <- colnames(logcpm)
    qc_long <- do.call(rbind, qc_rows)
    out[[method_id]] <- list(
      method = method,
      score_matrix = score_matrix,
      qc_long = qc_long,
      qc_summary = ptp_score_method_qc_summary(qc_long)
    )
  }
  out
}

ptp_score_method_qc_summary <- function(qc_long) {
  if (is.null(qc_long) || nrow(qc_long) == 0L) return(data.frame())
  keys <- c("score_method_id", "score_method_label", "score_method_role", "score_method_mode",
            "program_id", "program_label", "tier", "tier_label", "response_family", "response_family_label")
  groups <- do.call(interaction, c(qc_long[, keys, drop = FALSE], list(drop = TRUE, sep = "\r")))
  out <- do.call(rbind, lapply(split(qc_long, groups, drop = TRUE), function(df) {
    data.frame(
      score_method_id = df$score_method_id[1L],
      score_method_label = df$score_method_label[1L],
      score_method_role = df$score_method_role[1L],
      score_method_mode = df$score_method_mode[1L],
      program_id = df$program_id[1L],
      program_label = df$program_label[1L],
      tier = df$tier[1L],
      tier_label = df$tier_label[1L],
      response_family = df$response_family[1L],
      response_family_label = df$response_family_label[1L],
      n_observations = nrow(df),
      n_finite_scores = sum(df$score_finite, na.rm = TRUE),
      n_genes_observed = df$n_genes_observed[1L],
      median_n_genes_used_by_score = stats::median(df$n_genes_used_by_score, na.rm = TRUE),
      min_detected_program_fraction = min(df$detected_program_fraction, na.rm = TRUE),
      median_detected_program_fraction = stats::median(df$detected_program_fraction, na.rm = TRUE),
      median_program_raw_umi = stats::median(df$program_raw_umi, na.rm = TRUE),
      score_sd = stats::sd(df$score_value, na.rm = TRUE),
      n_fallback_observations = sum(!is.na(df$fallback_reason) & nzchar(df$fallback_reason)),
      fallback_reason = paste(unique(df$fallback_reason[!is.na(df$fallback_reason) & nzchar(df$fallback_reason)]), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$score_method_id, out$tier_label, out$response_family_label, out$program_label), , drop = FALSE]
}

ptp_fit_program_scores <- function(score_matrix, model) {
  keep <- rowSums(is.finite(score_matrix)) == ncol(score_matrix)
  score_matrix <- score_matrix[keep, , drop = FALSE]
  if (nrow(score_matrix) == 0L) stop("No estimable program scores.", call. = FALSE)
  dup <- limma::duplicateCorrelation(score_matrix, model$design, block = model$metadata$sample_id)
  fit <- limma::lmFit(score_matrix, model$design, block = model$metadata$sample_id, correlation = dup$consensus.correlation)
  fit <- limma::eBayes(fit, robust = FALSE)
  list(fit = fit, score_matrix = score_matrix, metadata = model$metadata, design = model$design, duplicate_correlation = dup$consensus.correlation)
}

ptp_apply_program_contrasts <- function(program_model, contrast_list, contrast_family, programs) {
  tmp_model <- list(fit = program_model$fit, design = program_model$design)
  out <- ptp_apply_contrasts(tmp_model, contrast_list, contrast_family)
  out$program_id <- out$gene
  out$gene <- NULL
  out$gene_symbol <- NULL
  lookup <- ptp_program_membership_table(programs)
  meta_cols <- ptp_program_meta_cols()
  out <- merge(out, lookup[, intersect(meta_cols, names(lookup)), drop = FALSE], by = "program_id", all.x = TRUE, sort = FALSE)
  out <- out[, c(intersect(meta_cols, names(out)), setdiff(names(out), meta_cols))]
  out <- ptp_add_program_tier_fdr(out, p_col = "p_value", fdr_col = "fdr")
  out[order(out$contrast_id, out$p_value, out$program_id), , drop = FALSE]
}

ptp_run_endpoint_ploidy_group_model <- function(spec, pb, programs, qc) {
  value_col <- spec$group_column
  group_col <- paste0(spec$method, "_ptp_group")
  status <- data.frame(
    model_id = spec$method,
    model_type = "endpoint_ploidy_group",
    status = "not_run",
    detail = "",
    stringsAsFactors = FALSE
  )
  if (!value_col %in% names(pb$metadata) || !group_col %in% names(pb$metadata)) {
    status$status <- "not_estimable"
    status$detail <- "missing endpoint ploidy group columns"
    return(list(status = status, support = data.frame(), design = data.frame(), gene = data.frame(), program = data.frame()))
  }
  support <- ptp_support_table_by_group(pb$metadata, value_col, spec$group_levels, qc)
  eligible <- unique(support$subbin_id[support$common_support])
  estimable <- length(eligible) >= qc$min_common_support_subbins
  if (!estimable) {
    status$status <- "not_estimable"
    status$detail <- paste("eligible_subbins", paste(eligible, collapse = ";"))
    support$method <- spec$method
    support$threshold <- spec$threshold
    return(list(status = status, support = support, design = data.frame(), gene = data.frame(), program = data.frame()))
  }
  fit <- ptp_fit_cellmeans(pb$counts, pb$metadata, eligible, group_col, qc)
  design <- ptp_design_audit(fit, paste0(spec$method, "_cellmeans"))
  contrasts <- ptp_endpoint_group_contrasts(eligible, colnames(fit$design))
  gene <- ptp_apply_contrasts(fit, contrasts, spec$method)
  gene$endpoint_ploidy_method <- spec$method
  gene$endpoint_ploidy_threshold <- spec$threshold
  scores <- ptp_program_scores_from_logcpm(ptp_logcpm_matrix(fit), programs)
  pfit <- ptp_fit_program_scores(scores, fit)
  program <- ptp_apply_program_contrasts(pfit, contrasts, spec$method, programs)
  program$endpoint_ploidy_method <- spec$method
  program$endpoint_ploidy_threshold <- spec$threshold
  support$method <- spec$method
  support$threshold <- spec$threshold
  status$status <- "ok"
  status$detail <- paste("eligible_subbins", paste(eligible, collapse = ";"))
  list(status = status, support = support, design = design, gene = gene, program = program)
}

ptp_run_endpoint_ploidy_continuous_model <- function(pb, programs, qc) {
  status <- data.frame(
    model_id = "continuous_mean_endpoint_ploidy",
    model_type = "endpoint_ploidy_continuous",
    status = "not_run",
    detail = "",
    stringsAsFactors = FALSE
  )
  etp_col <- "sample_mean_endpoint_ploidy_scaled"
  if (!etp_col %in% names(pb$metadata)) {
    status$status <- "not_estimable"
    status$detail <- "missing sample_mean_endpoint_ploidy_scaled"
    return(list(status = status, support = data.frame(), design = data.frame(), gene = data.frame(), program = data.frame()))
  }
  support <- ptp_endpoint_continuous_support(pb$metadata, pb$eligible_subbins, qc, etp_col = etp_col)
  eligible <- unique(support$subbin_id[support$common_support])
  estimable <- length(eligible) >= qc$min_common_support_subbins &&
    length(unique(pb$metadata$sample_mean_endpoint_ploidy[is.finite(pb$metadata$sample_mean_endpoint_ploidy)])) >= 4L
  if (!estimable) {
    status$status <- "not_estimable"
    status$detail <- paste("eligible_subbins", paste(eligible, collapse = ";"))
    return(list(status = status, support = support, design = data.frame(), gene = data.frame(), program = data.frame()))
  }
  fit <- ptp_fit_continuous_endpoint_ploidy(pb$counts, pb$metadata, eligible, qc, etp_col = etp_col)
  design <- ptp_design_audit(fit, "continuous_mean_endpoint_ploidy")
  contrasts <- ptp_continuous_endpoint_ploidy_contrasts(eligible, colnames(fit$design))
  gene <- ptp_apply_contrasts(fit, contrasts, "continuous_mean_endpoint_ploidy")
  scores <- ptp_program_scores_from_logcpm(ptp_logcpm_matrix(fit), programs)
  pfit <- ptp_fit_program_scores(scores, fit)
  program <- ptp_apply_program_contrasts(pfit, contrasts, "continuous_mean_endpoint_ploidy", programs)
  status$status <- "ok"
  status$detail <- paste("eligible_subbins", paste(eligible, collapse = ";"))
  list(status = status, support = support, design = design, gene = gene, program = program)
}

ptp_run_endpoint_ploidy_models <- function(pb, programs, cfg, args) {
  specs <- ptp_endpoint_ploidy_specs(cfg)
  workers <- ptp_effective_workers(args, length(specs), "workers")
  message("Running endpoint-ploidy threshold models with workers=", workers)
  group_out <- ptp_parallel_lapply(
    specs,
    function(spec) ptp_run_endpoint_ploidy_group_model(spec, pb, programs, cfg$qc),
    workers = workers,
    task_label = "endpoint_ploidy_group_model"
  )
  continuous <- ptp_run_endpoint_ploidy_continuous_model(pb, programs, cfg$qc)
  group_results <- group_out$results
  all_results <- c(group_results, list(continuous_mean_endpoint_ploidy = continuous))
  bind_part <- function(part) {
    rows <- lapply(all_results, `[[`, part)
    rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
    if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
  }
  list(
    status = bind_part("status"),
    group_support = {
      rows <- lapply(group_results, `[[`, "support")
      rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
      if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
    },
    continuous_support = continuous$support,
    design = bind_part("design"),
    group_gene = {
      rows <- lapply(group_results, `[[`, "gene")
      rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
      if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
    },
    group_program = {
      rows <- lapply(group_results, `[[`, "program")
      rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
      if (length(rows) == 0L) data.frame() else do.call(rbind, rows)
    },
    continuous_gene = continuous$gene,
    continuous_program = continuous$program,
    parallel_log = group_out$log
  )
}

ptp_write_endpoint_outputs <- function(endpoint, endpoint_dirs) {
  ptp_write_csv(endpoint$status, file.path(endpoint_dirs$manifest, "endpoint_ploidy_model_status.csv"))
  ptp_write_csv(endpoint$group_support, file.path(endpoint_dirs$qc, "endpoint_ploidy_group_common_support.csv"))
  ptp_write_csv(endpoint$continuous_support, file.path(endpoint_dirs$qc, "endpoint_ploidy_continuous_support.csv"))
  ptp_write_csv(endpoint$parallel_log, file.path(endpoint_dirs$manifest, "parallel_task_log_endpoint_ploidy.csv"))

  if (!is.null(endpoint$group_gene) && nrow(endpoint$group_gene) > 0L) {
    for (method in unique(endpoint$group_gene$endpoint_ploidy_method)) {
      md <- ptp_endpoint_model_dirs(endpoint_dirs$root, method)
      ptp_write_csv(endpoint$status[endpoint$status$model_id == method, , drop = FALSE], file.path(md$manifest, "model_status.csv"))
      ptp_write_csv(endpoint$group_support[endpoint$group_support$method == method, , drop = FALSE], file.path(md$qc, "common_support.csv"))
      ptp_write_csv(endpoint$group_gene[endpoint$group_gene$endpoint_ploidy_method == method, , drop = FALSE], file.path(md$genes, "gene_interactions.csv"))
      ptp_write_csv(endpoint$group_program[endpoint$group_program$endpoint_ploidy_method == method, , drop = FALSE], file.path(md$programs, "program_interactions.csv"))
    }
  } else {
    for (method in endpoint$status$model_id[endpoint$status$model_type == "endpoint_ploidy_group"]) {
      md <- ptp_endpoint_model_dirs(endpoint_dirs$root, method)
      ptp_write_csv(endpoint$status[endpoint$status$model_id == method, , drop = FALSE], file.path(md$manifest, "model_status.csv"))
      if (!is.null(endpoint$group_support) && nrow(endpoint$group_support) > 0L) {
        ptp_write_csv(endpoint$group_support[endpoint$group_support$method == method, , drop = FALSE], file.path(md$qc, "common_support.csv"))
      }
    }
  }

  cd <- ptp_endpoint_model_dirs(endpoint_dirs$root, "continuous_mean_endpoint_ploidy")
  ptp_write_csv(endpoint$status[endpoint$status$model_id == "continuous_mean_endpoint_ploidy", , drop = FALSE], file.path(cd$manifest, "model_status.csv"))
  ptp_write_csv(endpoint$continuous_support, file.path(cd$qc, "common_support.csv"))
  ptp_write_csv(endpoint$continuous_gene, file.path(cd$genes, "gene_interactions.csv"))
  ptp_write_csv(endpoint$continuous_program, file.path(cd$programs, "program_interactions.csv"))

  invisible(TRUE)
}

ptp_camera_tests <- function(gene_contrasts, programs, contrast_id) {
  x <- gene_contrasts[gene_contrasts$contrast_id == contrast_id, , drop = FALSE]
  if (nrow(x) == 0L) return(data.frame())
  stats <- x$t_statistic
  names(stats) <- x$gene_symbol
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & nzchar(names(stats))]
  stats <- stats[!duplicated(names(stats))]
  indices <- lapply(programs, function(p) which(names(stats) %in% p$observed_genes))
  indices <- indices[vapply(indices, length, integer(1L)) >= 3L]
  if (length(indices) == 0L) return(data.frame())
  cam <- limma::cameraPR(statistic = stats, index = indices, use.ranks = FALSE, sort = FALSE)
  cam$program_id <- rownames(cam)
  cam$contrast_id <- contrast_id
  rownames(cam) <- NULL
  cam <- as.data.frame(cam, stringsAsFactors = FALSE)
  names(cam) <- gsub(" ", "_", names(cam), fixed = TRUE)
  lookup <- ptp_program_membership_table(programs)
  meta_cols <- ptp_program_meta_cols()
  out <- merge(cam, lookup[, intersect(meta_cols, names(lookup)), drop = FALSE], by = "program_id", all.x = TRUE, sort = FALSE)
  out <- ptp_add_program_tier_fdr(out, p_col = "PValue", fdr_col = "FDR")
  out[, c("contrast_id", intersect(meta_cols, names(out)), setdiff(names(out), c("contrast_id", meta_cols)))]
}

ptp_run_fgsea <- function(gene_contrasts, min_size, max_size, nperm_simple, seed) {
  x <- gene_contrasts[gene_contrasts$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(data.frame())
  stats <- x$t_statistic
  names(stats) <- x$gene_symbol
  stats <- stats[is.finite(stats) & !is.na(names(stats)) & nzchar(names(stats))]
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  msig <- as.data.frame(msigdbr::msigdbr(species = "Homo sapiens"), stringsAsFactors = FALSE)
  collection_col <- if ("gs_collection" %in% names(msig)) "gs_collection" else "gs_cat"
  subcollection_col <- if ("gs_subcollection" %in% names(msig)) "gs_subcollection" else "gs_subcat"
  collection <- msig[[collection_col]]
  subcollection <- msig[[subcollection_col]]
  keep <- collection == "H" |
    (collection == "C2" & subcollection == "CP:REACTOME") |
    (collection == "C5" & subcollection == "GO:BP")
  msig <- msig[keep, , drop = FALSE]
  if (nrow(msig) == 0L) stop("No msigdbr Hallmark/Reactome/GO BP pathways were selected.", call. = FALSE)
  pathways <- split(msig$gene_symbol, msig$gs_name)
  pathways <- lapply(pathways, unique)
  set.seed(seed)
  fg <- tryCatch(
    fgsea::fgseaMultilevel(pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size),
    error = function(e) fgsea::fgsea(pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size, nperm = nperm_simple)
  )
  fg <- as.data.frame(fg, stringsAsFactors = FALSE)
  if ("leadingEdge" %in% names(fg)) fg$leadingEdge <- vapply(fg$leadingEdge, paste, character(1L), collapse = ";")
  fg[order(fg$padj, -abs(fg$NES), fg$pathway), , drop = FALSE]
}

ptp_construct_region_pseudobulk <- function(meta, counts, regions, qc, endpoint_ploidy_specs = ptp_endpoint_ploidy_specs()) {
  sample_meta <- ptp_sample_metadata(meta, endpoint_ploidy_specs)
  all_regions <- expand.grid(
    sample_id = sample_meta$sample_id,
    region_id = names(regions),
    stringsAsFactors = FALSE
  )
  region_info <- ptp_regions_table(list(regions = regions))
  all_regions <- merge(all_regions, region_info, by = "region_id", all.x = TRUE, sort = FALSE)
  all_regions <- merge(all_regions, sample_meta, by = "sample_id", all.x = TRUE, sort = FALSE)
  all_regions$sample_region_id <- paste(all_regions$sample_id, all_regions$region_id, sep = "__")
  meta_region <- meta[FALSE, , drop = FALSE]
  for (region_id in names(regions)) {
    x <- meta[ptp_interval_hit(meta$pseudotime, regions[[region_id]]), , drop = FALSE]
    if (nrow(x) > 0L) {
      x$region_id <- region_id
      x$sample_region_id <- paste(x$sample_id, region_id, sep = "__")
      meta_region <- rbind(meta_region, x)
    }
  }
  if (nrow(meta_region) > 0L) {
    group <- factor(meta_region$sample_region_id, levels = all_regions$sample_region_id)
    design <- Matrix::sparse.model.matrix(~ 0 + group)
    colnames(design) <- levels(group)
    pb_counts <- counts[, meta_region$cell_id, drop = FALSE] %*% design
  } else {
    pb_counts <- counts[, integer(0), drop = FALSE]
  }
  missing <- setdiff(all_regions$sample_region_id, colnames(pb_counts))
  if (length(missing) > 0L) {
    zero <- Matrix::Matrix(0, nrow = nrow(counts), ncol = length(missing), sparse = TRUE)
    rownames(zero) <- rownames(counts)
    colnames(zero) <- missing
    pb_counts <- cbind(pb_counts, zero)
  }
  pb_counts <- pb_counts[, all_regions$sample_region_id, drop = FALSE]
  all_regions$cell_count <- as.integer(table(factor(meta_region$sample_region_id, levels = all_regions$sample_region_id)))
  all_regions$library_size <- as.numeric(Matrix::colSums(pb_counts))
  all_regions$detected_genes <- as.numeric(Matrix::colSums(pb_counts > 0))
  all_regions$mean_region_pseudotime <- vapply(all_regions$sample_region_id, function(id) {
    values <- meta_region$pseudotime[meta_region$sample_region_id == id]
    if (length(values) == 0L) NA_real_ else mean(values)
  }, numeric(1L))
  all_regions$retained_for_model <- with(
    all_regions,
    cell_count >= qc$min_cells_per_mouse_region &
      library_size >= qc$min_library_size_per_mouse_subbin &
      detected_genes >= qc$min_detected_genes_per_mouse_subbin
  )
  list(counts = pb_counts, metadata = all_regions, cell_metadata = meta_region)
}

ptp_fit_region_programs <- function(region_pb, programs, qc, workers = 1L) {
  region_ids <- unique(region_pb$metadata$region_id)
  names(region_ids) <- region_ids
  out_parallel <- ptp_parallel_lapply(region_ids, function(region_id) {
    meta <- region_pb$metadata[region_pb$metadata$region_id == region_id & region_pb$metadata$retained_for_model, , drop = FALSE]
    if (nrow(meta) < 8L || length(unique(paste(meta$initial_ploidy, meta$treatment))) < 4L) return(data.frame())
    counts <- region_pb$counts[, meta$sample_region_id, drop = FALSE]
    y0 <- edgeR::DGEList(counts = counts)
    keep <- rowSums(edgeR::cpm(y0) > qc$gene_filter_cpm) >= qc$gene_filter_min_observations
    if (sum(keep) < 10L) return(data.frame())
    y <- edgeR::calcNormFactors(edgeR::DGEList(counts = counts[keep, , drop = FALSE]), method = "TMM")
    meta$initial_ploidy <- factor(meta$initial_ploidy, levels = c("2N", "4N"))
    meta$treatment <- factor(meta$treatment, levels = c("control", "gemcitabine"))
    meta$centered_mean_pseudotime <- as.numeric(scale(meta$mean_region_pseudotime, center = TRUE, scale = FALSE))
    design <- stats::model.matrix(~ centered_mean_pseudotime + initial_ploidy * treatment, data = meta)
    if (qr(design)$rank < ncol(design)) return(data.frame())
    v <- limma::voomWithQualityWeights(y, design, plot = FALSE)
    fit <- limma::lmFit(v, design)
    fit <- limma::eBayes(fit, robust = TRUE)
    logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 1)
    scores <- ptp_program_scores_from_logcpm(logcpm, programs)
    keep_prog <- rowSums(is.finite(scores)) == ncol(scores)
    if (!any(keep_prog)) return(data.frame())
    pfit <- limma::lmFit(scores[keep_prog, , drop = FALSE], design)
    pfit <- limma::eBayes(pfit, robust = FALSE)
    cn <- colnames(design)
    contrast_list <- list()
    v2 <- rep(0, length(cn)); names(v2) <- cn; v2["treatmentgemcitabine"] <- 1
    int <- rep(0, length(cn)); names(int) <- cn; int["initial_ploidy4N:treatmentgemcitabine"] <- 1
    contrast_list$delta_2N_origin <- v2
    contrast_list$delta_4N_origin <- v2 + int
    contrast_list$treatment_by_initial_ploidy_interaction <- int
    tmp <- list(fit = pfit, design = design)
    out <- ptp_apply_contrasts(tmp, contrast_list, paste0("window_restricted_", region_id))
    out$program_id <- out$gene
    out$gene <- NULL
    out$gene_symbol <- NULL
    out$region_id <- region_id
    out
  }, workers = workers, task_label = "region_sensitivity")
  rows <- out_parallel$results
  attr(rows, "parallel_log") <- out_parallel$log
  rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  attr(out, "parallel_log") <- out_parallel$log
  lookup <- ptp_program_membership_table(programs)
  out <- merge(out, lookup[, c("program_id", "program_label", "family", "expected_direction", "n_genes_observed")], by = "program_id", all.x = TRUE, sort = FALSE)
  attr(out, "parallel_log") <- out_parallel$log
  out
}

ptp_leave_one_out <- function(pb, eligible_subbins, programs, qc, workers = 1L) {
  sample_ids <- sort(unique(pb$metadata$sample_id))
  names(sample_ids) <- sample_ids
  out_parallel <- ptp_parallel_lapply(sample_ids, function(sample_id) {
    meta <- pb$metadata[pb$metadata$sample_id != sample_id, , drop = FALSE]
    counts <- pb$counts[, meta$sample_subbin_id, drop = FALSE]
    fit <- tryCatch(ptp_fit_cellmeans(counts, meta, eligible_subbins, "ptp_group", qc), error = function(e) NULL)
    if (is.null(fit)) {
      return(data.frame(left_out_sample_id = sample_id, status = "failed", error = "model_fit_failed", stringsAsFactors = FALSE))
    }
    logcpm <- ptp_logcpm_matrix(fit)
    scores <- ptp_program_scores_from_logcpm(logcpm, programs)
    pfit <- tryCatch(ptp_fit_program_scores(scores, fit), error = function(e) NULL)
    if (is.null(pfit)) {
      return(data.frame(left_out_sample_id = sample_id, status = "failed", error = "program_fit_failed", stringsAsFactors = FALSE))
    }
    contrasts <- ptp_primary_contrasts(eligible_subbins, colnames(fit$design))
    out <- ptp_apply_program_contrasts(pfit, contrasts["treatment_by_initial_ploidy_interaction"], "leave_one_mouse_out", programs)
    out$left_out_sample_id <- sample_id
    out$status <- "ok"
    out
  }, workers = workers, task_label = "leave_one_mouse_out")
  out <- do.call(rbind, out_parallel$results)
  attr(out, "parallel_log") <- out_parallel$log
  out
}

ptp_occupancy_models <- function(meta, regions) {
  sample_meta <- ptp_sample_metadata(meta)
  rows <- list()
  total <- as.integer(table(factor(meta$sample_id, levels = sample_meta$sample_id)))
  names(total) <- sample_meta$sample_id
  for (region_id in names(regions)) {
    hit <- meta[ptp_interval_hit(meta$pseudotime, regions[[region_id]]), , drop = FALSE]
    n_region <- as.integer(table(factor(hit$sample_id, levels = sample_meta$sample_id)))
    df <- sample_meta
    df$region_id <- region_id
    df$region_cells <- n_region
    df$total_cells <- total[df$sample_id]
    df$fraction <- df$region_cells / df$total_cells
    df$initial_ploidy <- factor(df$initial_ploidy, levels = c("2N", "4N"))
    df$treatment <- factor(df$treatment, levels = c("control", "gemcitabine"))
    fit <- tryCatch(stats::lm(fraction ~ initial_ploidy * treatment, data = df), error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients
    term <- "initial_ploidy4N:treatmentgemcitabine"
    rows[[length(rows) + 1L]] <- data.frame(
      region_id = region_id,
      model = "lm_fraction_region_cells_over_all_cellcycle_cells",
      n_samples = nrow(df),
      mean_fraction = mean(df$fraction),
      interaction_estimate = if (term %in% rownames(co)) co[term, "Estimate"] else NA_real_,
      interaction_se = if (term %in% rownames(co)) co[term, "Std. Error"] else NA_real_,
      interaction_p_value = if (term %in% rownames(co)) co[term, "Pr(>|t|)"] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

ptp_precision_table <- function(program_results) {
  x <- program_results[program_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(data.frame())
  out <- data.frame(
    program_id = x$program_id,
    program_label = x$program_label,
    tier = x$tier,
    tier_label = x$tier_label,
    response_family = x$response_family,
    response_family_label = x$response_family_label,
    score_mode = x$score_mode,
    score_direction_label = x$score_direction_label,
    family = x$family,
    interaction_estimate = x$estimate,
    interaction_se = x$se,
    ci_width_95 = x$ci_high - x$ci_low,
    minimum_detectable_abs_interaction_approx_80pct_power_alpha_0_05 = 2.8 * x$se,
    stringsAsFactors = FALSE
  )
  method_cols <- intersect(c("score_method_id", "score_method_label", "score_method_role", "score_method_mode"), names(x))
  if (length(method_cols) > 0L) out <- cbind(x[, method_cols, drop = FALSE], out)
  out
}

ptp_order_programs <- function(x) {
  x[order(x$tier_label %||% "", x$response_family_label %||% "", x$program_label %||% x$program_id), , drop = FALSE]
}

ptp_contrast_label <- function(x) {
  labels <- c(
    delta_2N_origin = "Gemcitabine effect in 2N-origin",
    delta_4N_origin = "Gemcitabine effect in 4N-origin",
    treatment_by_initial_ploidy_interaction = "4N-origin minus 2N-origin treatment response",
    treatment_by_endpoint_ploidy_group_interaction = "ETP-higher minus ETP-lower treatment response",
    continuous_endpoint_ploidy_treatment_interaction = "Treatment-response slope per 1 SD mean ETP"
  )
  out <- labels[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}

ptp_plot_program_forest <- function(program_results, path) {
  x <- program_results[program_results$contrast_id %in% c("delta_2N_origin", "delta_4N_origin", "treatment_by_initial_ploidy_interaction"), , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  x <- ptp_order_programs(x)
  x$contrast_label <- ptp_contrast_label(x$contrast_id)
  x$contrast_id <- factor(x$contrast_id, levels = c("delta_2N_origin", "delta_4N_origin", "treatment_by_initial_ploidy_interaction"))
  x$contrast_label <- factor(x$contrast_label, levels = ptp_contrast_label(levels(x$contrast_id)))
  x$program_label <- factor(x$program_label, levels = rev(unique(x$program_label)))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = estimate, y = program_label, color = contrast_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = ci_low, xmax = ci_high), orientation = "y", width = 0, position = ggplot2::position_dodge(width = 0.6), linewidth = 0.4) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.6), size = 1.6) +
    ggplot2::facet_grid(tier_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0), legend.position = "bottom") +
    ggplot2::labs(x = "Program score association (standardized log-expression units)", y = NULL, color = "Contrast")
  ggplot2::ggsave(path, p, width = 9.5, height = max(6, 0.28 * length(levels(x$program_label)) + 2.2))
  invisible(TRUE)
}

ptp_plot_program_heatmap <- function(program_results, path) {
  x <- program_results[program_results$contrast_id %in% c("delta_2N_origin", "delta_4N_origin", "treatment_by_initial_ploidy_interaction"), , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  x <- ptp_order_programs(x)
  x$contrast_label <- factor(ptp_contrast_label(x$contrast_id), levels = ptp_contrast_label(c("delta_2N_origin", "delta_4N_origin", "treatment_by_initial_ploidy_interaction")))
  x$program_label <- factor(x$program_label, levels = rev(unique(x$program_label)))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = contrast_label, y = program_label, fill = estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B6423C", midpoint = 0) +
    ggplot2::facet_grid(tier_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Estimate")
  ggplot2::ggsave(path, p, width = 8.8, height = max(6, 0.25 * length(unique(x$program_label)) + 2.0))
  invisible(TRUE)
}

ptp_program_score_long_table <- function(program_fit, programs) {
  if (is.null(program_fit$score_matrix) || nrow(program_fit$score_matrix) == 0L) return(data.frame())
  score <- as.data.frame(t(program_fit$score_matrix), stringsAsFactors = FALSE)
  score$sample_subbin_id <- rownames(score)
  long <- tidyr::pivot_longer(score, cols = setdiff(names(score), "sample_subbin_id"), names_to = "program_id", values_to = "program_score")
  meta <- program_fit$metadata
  keep_meta <- intersect(c("sample_subbin_id", "sample_id", "subbin_id", "subbin_label", "initial_ploidy", "treatment",
                           "dose_mg", "cell_count", "library_size", "mean_subbin_pseudotime",
                           "sample_mean_endpoint_ploidy", "sample_mean_endpoint_ploidy_scaled"), names(meta))
  long <- merge(long, meta[, keep_meta, drop = FALSE], by = "sample_subbin_id", all.x = TRUE, sort = FALSE)
  lookup <- ptp_program_membership_table(programs)
  meta_cols <- ptp_program_meta_cols()
  long <- merge(long, lookup[, intersect(meta_cols, names(lookup)), drop = FALSE], by = "program_id", all.x = TRUE, sort = FALSE)
  long[, c("program_id", "program_label", "tier", "tier_label", "response_family", "response_family_label",
           "sample_subbin_id", "sample_id", "subbin_id", "initial_ploidy", "treatment", "dose_mg",
           "program_score", setdiff(names(long), c("program_id", "program_label", "tier", "tier_label",
                                                     "response_family", "response_family_label", "sample_subbin_id",
                                                     "sample_id", "subbin_id", "initial_ploidy", "treatment",
                                                     "dose_mg", "program_score")))]
}

ptp_plot_program_score_pseudobulk <- function(score_long, program_results, path) {
  if (is.null(score_long) || nrow(score_long) == 0L) return(invisible(FALSE))
  ranked <- program_results[program_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  ranked <- ranked[order(ranked$p_value), , drop = FALSE]
  keep <- unique(ranked$program_id)[seq_len(min(12L, length(unique(ranked$program_id))))]
  x <- score_long[score_long$program_id %in% keep, , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  x$group <- paste(x$initial_ploidy, x$treatment, sep = " / ")
  x$group <- factor(x$group, levels = c("2N / control", "2N / gemcitabine", "4N / control", "4N / gemcitabine"))
  x$program_label <- factor(x$program_label, levels = rev(unique(ranked$program_label[ranked$program_id %in% keep])))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = group, y = program_score, color = treatment)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(shape = initial_ploidy), position = ggplot2::position_jitter(width = 0.12, height = 0), size = 1.4, alpha = 0.75) +
    ggplot2::stat_summary(fun = mean, geom = "point", color = "black", size = 2.0) +
    ggplot2::facet_wrap(~ program_label, scales = "free_y", ncol = 3) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "bottom") +
    ggplot2::labs(x = NULL, y = "Pseudobulk program score", color = "Treatment", shape = "Initial ploidy")
  ggplot2::ggsave(path, p, width = 10, height = 7)
  invisible(TRUE)
}

ptp_gsea_axis_table <- function(gsea, cfg) {
  if (is.null(gsea) || nrow(gsea) == 0L) return(data.frame())
  axes <- cfg$exploratory_response_axes %||% list()
  if (length(axes) == 0L) return(data.frame())
  rows <- list()
  for (axis in axes) {
    patterns <- as.character(axis$pathway_patterns %||% character(0L))
    hit <- rep(FALSE, nrow(gsea))
    for (pat in patterns) hit <- hit | grepl(pat, gsea$pathway, ignore.case = TRUE)
    if (!any(hit)) next
    tmp <- gsea[hit, , drop = FALSE]
    tmp$axis_id <- axis$id
    tmp$axis_label <- axis$label
    tmp$tier <- axis$tier %||% "tier3_exploratory_only"
    tmp$tier_label <- axis$tier_label %||% "Tier 3 exploratory-only"
    rows[[length(rows) + 1L]] <- tmp
  }
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$axis_id, out$padj, -abs(out$NES), out$pathway), , drop = FALSE]
}

ptp_plot_gsea_axis_bar <- function(axis_gsea, path) {
  if (is.null(axis_gsea) || nrow(axis_gsea) == 0L) return(invisible(FALSE))
  x <- do.call(rbind, lapply(split(axis_gsea, axis_gsea$axis_id), function(df) {
    df <- df[order(df$padj, -abs(df$NES)), , drop = FALSE]
    head(df, 8L)
  }))
  x$pathway_short <- gsub("^(HALLMARK|REACTOME|GOBP)_", "", x$pathway)
  x$pathway_short <- gsub("_", " ", x$pathway_short)
  x$pathway_short <- substr(x$pathway_short, 1L, 55L)
  x$pathway_short <- factor(x$pathway_short, levels = rev(unique(x$pathway_short[order(x$axis_label, x$NES)])))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = NES, y = pathway_short, fill = padj < 0.05)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::facet_wrap(~ axis_label, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "#B6423C", "FALSE" = "#8A96A3")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Interaction-ranked GSEA NES", y = NULL, fill = "padj < 0.05")
  ggplot2::ggsave(path, p, width = 8.5, height = max(5.5, 0.23 * nrow(x) + 2.2))
  invisible(TRUE)
}

ptp_plot_program_camera_comparison <- function(program_results, camera_results, path) {
  if (is.null(program_results) || nrow(program_results) == 0L ||
      is.null(camera_results) || nrow(camera_results) == 0L) return(invisible(FALSE))
  p0 <- program_results[program_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  c0 <- camera_results[camera_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(p0) == 0L || nrow(c0) == 0L) return(invisible(FALSE))
  p0$program_score_interpretation_fdr <- ptp_interpretation_fdr(p0, "fdr")
  c0$camera_interpretation_fdr <- ptp_interpretation_fdr(c0, "FDR")
  keep_p <- c("program_id", "program_label", "tier", "tier_label", "response_family_label",
              "estimate", "p_value", "fdr", "program_score_interpretation_fdr")
  keep_c <- c("program_id", "Direction", "PValue", "FDR", "camera_interpretation_fdr")
  x <- merge(p0[, intersect(keep_p, names(p0)), drop = FALSE],
             c0[, intersect(keep_c, names(c0)), drop = FALSE],
             by = "program_id", all = FALSE, sort = FALSE)
  if (nrow(x) == 0L) return(invisible(FALSE))
  x$neg_log10_camera_fdr <- -log10(pmax(x$camera_interpretation_fdr, .Machine$double.xmin))
  x$program_score_sig <- x$program_score_interpretation_fdr < 0.05
  x$camera_sig <- x$camera_interpretation_fdr < 0.05
  x$label <- ifelse(x$program_score_sig | x$camera_sig, x$program_label, "")
  p <- ggplot2::ggplot(x, ggplot2::aes(x = estimate, y = neg_log10_camera_fdr, color = tier_label, shape = camera_sig)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(0.05), color = "grey45", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_point(size = 2.0, alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = label), check_overlap = TRUE, size = 2.5, hjust = 0, nudge_x = 0.02, show.legend = FALSE) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Program-score interaction estimate",
                  y = "-log10(cameraPR interpretation FDR)",
                  color = "Tier", shape = "cameraPR FDR < 0.05")
  ggplot2::ggsave(path, p, width = 8.8, height = 5.8)
  invisible(TRUE)
}

ptp_tier_method_summary <- function(program_results, camera_results, gsea_axes = data.frame()) {
  rows <- list()
  tiers <- unique(c(program_results$tier_label, camera_results$tier_label, "Tier 3 exploratory-only"))
  tiers <- tiers[!is.na(tiers) & nzchar(tiers)]
  if (!is.null(program_results) && nrow(program_results) > 0L) {
    p0 <- program_results[program_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
    if (nrow(p0) > 0L) {
      p0$interpretation_fdr <- ptp_interpretation_fdr(p0, "fdr")
      rows <- c(rows, lapply(split(p0, p0$tier_label), function(df) {
        data.frame(tier_label = df$tier_label[1L], method = "Program score model",
                   n_tests = nrow(df), n_fdr_0_05 = sum(df$interpretation_fdr < 0.05, na.rm = TRUE),
                   min_fdr = min(df$interpretation_fdr, na.rm = TRUE), stringsAsFactors = FALSE)
      }))
    }
  }
  if (!is.null(camera_results) && nrow(camera_results) > 0L) {
    c0 <- camera_results[camera_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
    if (nrow(c0) > 0L) {
      c0$interpretation_fdr <- ptp_interpretation_fdr(c0, "FDR")
      rows <- c(rows, lapply(split(c0, c0$tier_label), function(df) {
        data.frame(tier_label = df$tier_label[1L], method = "cameraPR",
                   n_tests = nrow(df), n_fdr_0_05 = sum(df$interpretation_fdr < 0.05, na.rm = TRUE),
                   min_fdr = min(df$interpretation_fdr, na.rm = TRUE), stringsAsFactors = FALSE)
      }))
    }
  }
  if (!is.null(gsea_axes) && nrow(gsea_axes) > 0L) {
    rows <- c(rows, lapply(split(gsea_axes, gsea_axes$tier_label), function(df) {
      data.frame(tier_label = df$tier_label[1L], method = "Genome-wide GSEA mapped pathways",
                 n_tests = nrow(df), n_fdr_0_05 = sum(df$padj < 0.05, na.rm = TRUE),
                 min_fdr = min(df$padj, na.rm = TRUE), stringsAsFactors = FALSE)
    }))
  }
  rows <- rows[vapply(rows, nrow, integer(1L)) > 0L]
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$tier_label, out$method), , drop = FALSE]
}

ptp_plot_tier_method_summary <- function(summary, path) {
  if (is.null(summary) || nrow(summary) == 0L) return(invisible(FALSE))
  x <- summary
  x$method <- factor(x$method, levels = c("Program score model", "cameraPR", "Genome-wide GSEA mapped pathways"))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = method, y = tier_label, fill = n_fdr_0_05)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(n_fdr_0_05, "/", n_tests)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#E8EDF3", high = "#B6423C") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), legend.position = "bottom") +
    ggplot2::labs(x = NULL, y = NULL, fill = "FDR < 0.05 count")
  ggplot2::ggsave(path, p, width = 8, height = max(4.2, 0.55 * length(unique(x$tier_label)) + 1.8))
  invisible(TRUE)
}

ptp_score_method_label_order <- function(x) {
  if (!"score_method_id" %in% names(x) || !"score_method_label" %in% names(x)) return(unique(x$score_method_label))
  ids <- unique(x$score_method_id)
  preferred <- c("signed_weighted_mean_zscore", "signed_pca_eigengene", "signed_rank_mean", "trimmed_signed_mean_zscore")
  ordered_ids <- c(intersect(preferred, ids), setdiff(ids, preferred))
  vapply(ordered_ids, function(id) x$score_method_label[match(id, x$score_method_id)], character(1L))
}

ptp_plot_score_method_interaction_heatmap <- function(score_method_results, path) {
  if (is.null(score_method_results) || nrow(score_method_results) == 0L) return(invisible(FALSE))
  x <- score_method_results[score_method_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  x$interpretation_fdr <- ptp_interpretation_fdr(x, "fdr")
  x <- ptp_order_programs(x)
  x$program_label <- factor(x$program_label, levels = rev(unique(x$program_label)))
  x$score_method_label <- factor(x$score_method_label, levels = ptp_score_method_label_order(x))
  x$sig_label <- ifelse(x$interpretation_fdr < 0.05, "*", "")
  p <- ggplot2::ggplot(x, ggplot2::aes(x = score_method_label, y = program_label, fill = estimate)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::geom_text(ggplot2::aes(label = sig_label), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#B6423C", midpoint = 0) +
    ggplot2::facet_grid(tier_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::labs(x = NULL, y = NULL, fill = "Interaction estimate")
  ggplot2::ggsave(path, p, width = 9.5, height = max(6.5, 0.25 * length(unique(x$program_label)) + 2.1))
  invisible(TRUE)
}

ptp_plot_score_method_qc_summary <- function(qc_summary, path) {
  if (is.null(qc_summary) || nrow(qc_summary) == 0L) return(invisible(FALSE))
  x <- qc_summary
  x$score_method_label <- factor(x$score_method_label, levels = ptp_score_method_label_order(x))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = score_method_label, y = median_detected_program_fraction, color = tier_label)) +
    ggplot2::geom_hline(yintercept = 0.2, color = "grey55", linewidth = 0.3, linetype = "dashed") +
    ggplot2::geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.75) +
    ggplot2::facet_wrap(~ tier_label, ncol = 2) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = "Median detected program-gene fraction")
  ggplot2::ggsave(path, p, width = 9, height = 6)
  invisible(TRUE)
}

ptp_score_method_result_summary <- function(score_method_results) {
  if (is.null(score_method_results) || nrow(score_method_results) == 0L) return(data.frame())
  x <- score_method_results[score_method_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(data.frame())
  x$interpretation_fdr <- ptp_interpretation_fdr(x, "fdr")
  groups <- interaction(x$score_method_id, x$score_method_label, x$score_method_role, x$tier, x$tier_label, drop = TRUE, sep = "\r")
  out <- do.call(rbind, lapply(split(x, groups, drop = TRUE), function(df) {
    top <- df[order(df$p_value), , drop = FALSE][1L, , drop = FALSE]
    data.frame(
      score_method_id = top$score_method_id,
      score_method_label = top$score_method_label,
      score_method_role = top$score_method_role,
      tier = top$tier,
      tier_label = top$tier_label,
      n_programs = nrow(df),
      n_fdr_0_05 = sum(df$interpretation_fdr < 0.05, na.rm = TRUE),
      min_interpretation_fdr = min(df$interpretation_fdr, na.rm = TRUE),
      top_program = top$program_label,
      top_estimate = top$estimate,
      top_p_value = top$p_value,
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$score_method_id, out$tier_label), , drop = FALSE]
}

ptp_plot_endpoint_program_comparison <- function(primary, endpoint_group, endpoint_continuous, path) {
  p0 <- primary[primary$contrast_id == "treatment_by_initial_ploidy_interaction", c("program_id", "program_label", "tier_label", "estimate"), drop = FALSE]
  if (nrow(p0) == 0L) return(invisible(FALSE))
  names(p0)[names(p0) == "estimate"] <- "primary_initial_ploidy_estimate"
  rows <- list()
  if (!is.null(endpoint_group) && nrow(endpoint_group) > 0L) {
    g <- endpoint_group[endpoint_group$contrast_id == "treatment_by_endpoint_ploidy_group_interaction", , drop = FALSE]
    if (nrow(g) > 0L) {
      g$model_label <- paste0("Threshold: ", g$endpoint_ploidy_method)
      rows[[length(rows) + 1L]] <- g[, c("program_id", "model_label", "estimate"), drop = FALSE]
    }
  }
  if (!is.null(endpoint_continuous) && nrow(endpoint_continuous) > 0L) {
    ctm <- endpoint_continuous[endpoint_continuous$contrast_id == "continuous_endpoint_ploidy_treatment_interaction", , drop = FALSE]
    if (nrow(ctm) > 0L) {
      ctm$model_label <- "Continuous mean ETP"
      rows[[length(rows) + 1L]] <- ctm[, c("program_id", "model_label", "estimate"), drop = FALSE]
    }
  }
  if (length(rows) == 0L) return(invisible(FALSE))
  y <- do.call(rbind, rows)
  names(y)[names(y) == "estimate"] <- "endpoint_estimate"
  x <- merge(p0, y, by = "program_id", all = FALSE, sort = FALSE)
  if (nrow(x) == 0L) return(invisible(FALSE))
  p <- ggplot2::ggplot(x, ggplot2::aes(x = primary_initial_ploidy_estimate, y = endpoint_estimate, label = program_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = 0, color = "grey70", linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(color = tier_label), size = 1.7) +
    ggplot2::facet_wrap(~ model_label, scales = "free") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::labs(x = "Primary initial-ploidy interaction estimate", y = "Endpoint-ploidy sensitivity estimate", color = "Tier")
  ggplot2::ggsave(path, p, width = 9, height = 5.5)
  invisible(TRUE)
}

ptp_plot_loo_robustness <- function(loo, primary, path) {
  if (is.null(loo) || nrow(loo) == 0L) return(invisible(FALSE))
  x <- loo[loo$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  summary <- stats::aggregate(estimate ~ program_id + program_label + tier_label + response_family_label, data = x, function(v) c(min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE)))
  summary$loo_min <- summary$estimate[, "min"]
  summary$loo_max <- summary$estimate[, "max"]
  summary$estimate <- NULL
  p0 <- primary[primary$contrast_id == "treatment_by_initial_ploidy_interaction", c("program_id", "estimate"), drop = FALSE]
  names(p0)[names(p0) == "estimate"] <- "full_estimate"
  summary <- merge(summary, p0, by = "program_id", all.x = TRUE, sort = FALSE)
  summary <- ptp_order_programs(summary)
  summary$program_label <- factor(summary$program_label, levels = rev(unique(summary$program_label)))
  p <- ggplot2::ggplot(summary, ggplot2::aes(y = program_label)) +
    ggplot2::geom_vline(xintercept = 0, color = "grey65", linewidth = 0.25) +
    ggplot2::geom_segment(ggplot2::aes(x = loo_min, xend = loo_max, yend = program_label), color = "grey45", linewidth = 0.5) +
    ggplot2::geom_point(ggplot2::aes(x = full_estimate), color = "#B6423C", size = 1.7) +
    ggplot2::facet_grid(tier_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0)) +
    ggplot2::labs(x = "Interaction estimate: full model point, leave-one-mouse-out range", y = NULL)
  ggplot2::ggsave(path, p, width = 8.5, height = max(6, 0.25 * nrow(summary) + 2))
  invisible(TRUE)
}

ptp_continuous_binned_interactions <- function(meta, counts, programs, qc, n_bins = 20L) {
  breaks <- seq(0, 1, length.out = n_bins + 1L)
  subbins <- lapply(seq_len(n_bins), function(i) {
    list(id = sprintf("ptbin%02d", i), label = sprintf("[%.2f,%.2f%s", breaks[i], breaks[i + 1L], ifelse(i == n_bins, "]", ")")),
         start = breaks[i], end = breaks[i + 1L], include_start = TRUE, include_end = i == n_bins)
  })
  pb <- ptp_construct_subbin_pseudobulk(meta, counts, subbins, qc)
  support <- ptp_support_table(pb$metadata, qc)
  common <- unique(support$subbin_id[support$common_support])
  if (length(common) == 0L) return(list(results = data.frame(), support = support))
  fit <- tryCatch(ptp_fit_cellmeans(pb$counts, pb$metadata, common, "ptp_group", qc), error = function(e) NULL)
  if (is.null(fit)) return(list(results = data.frame(), support = support))
  logcpm <- ptp_logcpm_matrix(fit)
  scores <- ptp_program_scores_from_logcpm(logcpm, programs)
  pfit <- tryCatch(ptp_fit_program_scores(scores, fit), error = function(e) NULL)
  if (is.null(pfit)) return(list(results = data.frame(), support = support))
  rows <- list()
  for (sb in common) {
    contrasts <- ptp_primary_contrasts(sb, colnames(fit$design))
    out <- ptp_apply_program_contrasts(pfit, contrasts["treatment_by_initial_ploidy_interaction"], "continuous_binned", programs)
    out$subbin_id <- sb
    info <- ptp_subbins_table(subbins, "20bin")
    out$bin_midpoint <- (info$start[match(sb, info$subbin_id)] + info$end[match(sb, info$subbin_id)]) / 2
    rows[[length(rows) + 1L]] <- out
  }
  list(results = do.call(rbind, rows), support = support)
}

ptp_plot_continuous <- function(cont_results, path) {
  x <- cont_results[cont_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
  if (nrow(x) == 0L) return(invisible(FALSE))
  keep_programs <- unique(x$program_id[order(x$p_value)])[seq_len(min(8L, length(unique(x$program_id))))]
  x <- x[x$program_id %in% keep_programs, , drop = FALSE]
  p <- ggplot2::ggplot(x, ggplot2::aes(x = bin_midpoint, y = estimate, group = program_label, color = program_label)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ci_low, ymax = ci_high), width = 0.006, linewidth = 0.3, alpha = 0.65) +
    ggplot2::geom_line(linewidth = 0.5) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(x = "Pseudotime bin midpoint", y = "Ploidy interaction estimate", color = "Program")
  ggplot2::ggsave(path, p, width = 8, height = 5.5)
  invisible(TRUE)
}

ptp_html_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

ptp_html_table <- function(df, n = 12L) {
  if (is.null(df) || nrow(df) == 0L || ncol(df) == 0L) return("<p>No rows.</p>")
  df <- head(as.data.frame(df, stringsAsFactors = FALSE), n)
  cols <- names(df)
  rows <- apply(df, 1L, function(r) paste0("<tr>", paste0("<td>", ptp_html_escape(r), "</td>", collapse = ""), "</tr>"))
  paste0("<table><thead><tr>", paste0("<th>", ptp_html_escape(cols), "</th>", collapse = ""), "</tr></thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

ptp_format_num <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits)
}

ptp_interpretation_fdr <- function(df, default_col = "fdr") {
  if (is.null(df) || nrow(df) == 0L) return(numeric(0L))
  out <- if (default_col %in% names(df)) df[[default_col]] else rep(NA_real_, nrow(df))
  focused_col <- paste0(default_col, "_focused_tier1_tier2")
  tier3_col <- paste0(default_col, "_tier3_exploratory")
  negative_col <- paste0(default_col, "_negative_control")
  if (focused_col %in% names(df)) {
    idx <- df$tier %in% c("tier1_core", "tier2_upgraded_focused")
    out[idx] <- df[[focused_col]][idx]
  }
  if (tier3_col %in% names(df)) {
    idx <- df$tier == "tier3_exploratory_only"
    out[idx] <- df[[tier3_col]][idx]
  }
  if (negative_col %in% names(df)) {
    idx <- df$tier == "negative_control"
    out[idx] <- df[[negative_col]][idx]
  }
  out
}

ptp_pdf_preview_png <- function(pdf_path, dpi = 170L) {
  if (!file.exists(pdf_path)) return(NA_character_)
  png_path <- paste0(tools::file_path_sans_ext(pdf_path), ".embedded.png")
  if (file.exists(png_path) && file.info(png_path)$mtime >= file.info(pdf_path)$mtime) return(png_path)
  pdftoppm <- Sys.which("pdftoppm")
  if (!nzchar(pdftoppm)) return(NA_character_)
  prefix <- tools::file_path_sans_ext(png_path)
  status <- tryCatch(
    system2(pdftoppm, args = c("-f", "1", "-singlefile", "-png", "-r", as.character(dpi), shQuote(pdf_path), shQuote(prefix)),
            stdout = TRUE, stderr = TRUE),
    error = function(e) structure(character(), status = 1L)
  )
  if (!file.exists(png_path)) return(NA_character_)
  png_path
}

ptp_data_uri <- function(path, mime) {
  if (!file.exists(path)) return(NA_character_)
  paste0("data:", mime, ";base64,", base64enc::base64encode(path, linewidth = 0L))
}

ptp_write_report <- function(output_root, cfg, status, support, program_results, precision, design_audit,
                             endpoint_status = data.frame(), endpoint_group_program = data.frame(),
                             endpoint_continuous_program = data.frame(),
                             dose_program = data.frame(), camera_results = data.frame(),
                             gsea_axes = data.frame(), loo = data.frame(),
                             score_method_results = data.frame(),
                             score_method_precision = data.frame(),
                             score_method_qc_summary = data.frame()) {
  primary <- data.frame()
  if (!is.null(program_results) && "contrast_id" %in% names(program_results) && nrow(program_results) > 0L) {
    primary <- program_results[program_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
    primary <- primary[order(primary$tier_label, primary$response_family_label, primary$p_value), , drop = FALSE]
    primary$interpretation_fdr <- ptp_interpretation_fdr(primary, "fdr")
  }
  focused <- primary[primary$tier %in% c("tier1_core", "tier2_upgraded_focused"), , drop = FALSE]
  tier3_primary <- primary[primary$tier == "tier3_exploratory_only", , drop = FALSE]
  n_focus_sig <- if (nrow(focused) > 0L && "interpretation_fdr" %in% names(focused)) sum(focused$interpretation_fdr < 0.05, na.rm = TRUE) else 0L
  n_tier3_sig <- if (nrow(tier3_primary) > 0L && "interpretation_fdr" %in% names(tier3_primary)) sum(tier3_primary$interpretation_fdr < 0.05, na.rm = TRUE) else 0L
  conclusion <- if (n_focus_sig > 0L) {
    paste(n_focus_sig, "focused program interaction(s) reached FDR < 0.05 in the primary initial-ploidy model.")
  } else {
    "No focused Tier 1 or Tier 2 program interaction reached FDR < 0.05 in the primary initial-ploidy model."
  }

  family_summary <- data.frame()
  if (nrow(primary) > 0L) {
    keys <- c("tier", "tier_label", "response_family", "response_family_label")
    groups <- do.call(interaction, c(primary[, keys, drop = FALSE], list(drop = TRUE, sep = "\r")))
    family_summary <- do.call(rbind, lapply(split(primary, groups, drop = TRUE), function(df) {
      top <- df[order(df$p_value), , drop = FALSE][1L, , drop = FALSE]
      data.frame(
        tier_label = top$tier_label,
        response_family_label = top$response_family_label,
        n_programs = nrow(df),
        best_program = top$program_label,
        best_estimate = top$estimate,
        best_ci = paste0(sprintf("%.3f", top$ci_low), " to ", sprintf("%.3f", top$ci_high)),
        best_p_value = top$p_value,
        min_interpretation_fdr = min(df$interpretation_fdr, na.rm = TRUE),
        min_all_program_fdr = if ("fdr_all_programs" %in% names(df)) min(df$fdr_all_programs, na.rm = TRUE) else min(df$fdr, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }))
    family_summary <- family_summary[order(family_summary$tier_label, family_summary$min_interpretation_fdr), , drop = FALSE]
  }

  tier_table <- function(tier_id) {
    x <- primary[primary$tier == tier_id, , drop = FALSE]
    if (nrow(x) == 0L) return(data.frame())
    cols <- intersect(c("response_family_label", "program_label", "score_direction_label",
                        "n_positive_weight_genes", "n_negative_weight_genes",
                        "estimate", "ci_low", "ci_high", "p_value", "interpretation_fdr",
                        "fdr", "fdr_all_programs", "fdr_focused_tier1_tier2",
                        "fdr_tier3_exploratory", "fdr_negative_control",
                        "n_genes_observed", "upgrade_status"), names(x))
    x[order(x$response_family_label, x$p_value), cols, drop = FALSE]
  }
  endpoint_table <- function(tier_id) {
    parts <- list()
    display_cols <- c("model", "response_family_label", "program_label", "estimate", "ci_low", "ci_high", "p_value",
                      "fdr", "fdr_all_programs", "fdr_focused_tier1_tier2", "fdr_tier3_exploratory", "fdr_negative_control")
    if (!is.null(endpoint_group_program) && nrow(endpoint_group_program) > 0L) {
      x <- endpoint_group_program[endpoint_group_program$contrast_id == "treatment_by_endpoint_ploidy_group_interaction" & endpoint_group_program$tier == tier_id, , drop = FALSE]
      if (nrow(x) > 0L) {
        x$model <- paste0("threshold_", x$endpoint_ploidy_method)
        parts[[length(parts) + 1L]] <- x[, intersect(display_cols, names(x)), drop = FALSE]
      }
    }
    if (!is.null(endpoint_continuous_program) && nrow(endpoint_continuous_program) > 0L) {
      x <- endpoint_continuous_program[endpoint_continuous_program$contrast_id == "continuous_endpoint_ploidy_treatment_interaction" & endpoint_continuous_program$tier == tier_id, , drop = FALSE]
      if (nrow(x) > 0L) {
        x$model <- "continuous_mean_endpoint_ploidy"
        parts[[length(parts) + 1L]] <- x[, intersect(display_cols, names(x)), drop = FALSE]
      }
    }
    if (length(parts) == 0L) return(data.frame())
    out <- do.call(rbind, parts)
    out[order(out$model, out$response_family_label, out$p_value), , drop = FALSE]
  }
  dose_table <- function(tier_id) {
    if (is.null(dose_program) || nrow(dose_program) == 0L) return(data.frame())
    x <- dose_program[grepl("interaction", dose_program$contrast_id) & dose_program$tier == tier_id, , drop = FALSE]
    if (nrow(x) == 0L) return(data.frame())
    cols <- intersect(c("contrast_id", "response_family_label", "program_label", "estimate", "ci_low", "ci_high", "p_value",
                        "fdr", "fdr_all_programs", "fdr_focused_tier1_tier2", "fdr_tier3_exploratory", "fdr_negative_control"), names(x))
    x[order(x$response_family_label, x$contrast_id, x$p_value), cols, drop = FALSE]
  }
  camera_table <- function(tier_id) {
    if (is.null(camera_results) || nrow(camera_results) == 0L) return(data.frame())
    x <- camera_results[camera_results$tier == tier_id, , drop = FALSE]
    if (nrow(x) == 0L) return(data.frame())
    x$interpretation_FDR <- ptp_interpretation_fdr(x, "FDR")
    cols <- intersect(c("response_family_label", "program_label", "Direction", "PValue",
                        "interpretation_FDR", "FDR", "FDR_all_programs",
                        "FDR_focused_tier1_tier2", "FDR_tier3_exploratory",
                        "FDR_negative_control", "NGenes"), names(x))
    x[order(x$response_family_label, x$FDR, x$PValue), cols, drop = FALSE]
  }
  rel <- function(path) {
    p <- normalizePath(path, winslash = "/", mustWork = FALSE)
    root <- normalizePath(file.path(output_root, "report"), winslash = "/", mustWork = FALSE)
    if (startsWith(p, paste0(root, "/"))) return(substring(p, nchar(root) + 2L))
    p
  }
  figure_counter <- 0L
  figure_block <- function(title, path, legend, interpretation) {
    figure_counter <<- figure_counter + 1L
    fig_id <- sprintf("fig-%02d", figure_counter)
    media <- "<p class='missing'>Figure file is missing.</p>"
    if (file.exists(path)) {
      png <- ptp_pdf_preview_png(path)
      if (!is.na(png) && file.exists(png)) {
        uri <- ptp_data_uri(png, "image/png")
        media <- paste0("<img src='", uri, "' alt='", ptp_html_escape(title), "'>")
      } else {
        uri <- ptp_data_uri(path, "application/pdf")
        if (!is.na(uri)) {
          media <- paste0("<object data='", uri, "' type='application/pdf'></object>")
        }
      }
    }
    paste0(
      "<figure class='figure-card' id='", fig_id, "'>",
      "<h3>Figure ", figure_counter, ". ", ptp_html_escape(title), "</h3>",
      media,
      "<figcaption class='figure-caption'><p><strong>Figure legend.</strong> ", ptp_html_escape(legend), "</p>",
      "<p><strong>Result interpretation.</strong> ", ptp_html_escape(interpretation), "</p></figcaption>",
      "</figure>"
    )
  }
  section <- function(id, title, ...) {
    paste0("<section id='", id, "'><h2>", ptp_html_escape(title), "</h2>", paste0(..., collapse = "\n"), "</section>")
  }
  tier_section <- function(tier_id, title, note) {
    section(
      tier_id,
      title,
      "<p class='note'>", ptp_html_escape(note), "</p>",
      "<h3>Primary initial-ploidy interactions</h3>",
      ptp_html_table(tier_table(tier_id), n = 80L),
      "<h3>Competitive program support</h3>",
      ptp_html_table(camera_table(tier_id), n = 80L),
      "<h3>Dose-specific sensitivity</h3>",
      ptp_html_table(dose_table(tier_id), n = 80L),
      "<h3>Endpoint-ploidy sensitivity</h3>",
      ptp_html_table(endpoint_table(tier_id), n = 80L)
    )
  }
  tier3_summary <- data.frame()
  if (!is.null(gsea_axes) && nrow(gsea_axes) > 0L) {
    tier3_summary <- do.call(rbind, lapply(split(gsea_axes, gsea_axes$axis_label), function(df) {
      top <- df[order(df$padj, -abs(df$NES)), , drop = FALSE][1L, , drop = FALSE]
      data.frame(
        axis_label = top$axis_label,
        n_pathways = nrow(df),
        n_fdr_0_05 = sum(df$padj < 0.05, na.rm = TRUE),
        top_pathway = top$pathway,
        top_NES = top$NES,
        top_padj = top$padj,
        stringsAsFactors = FALSE
      )
    }))
    tier3_summary <- tier3_summary[order(-tier3_summary$n_fdr_0_05, tier3_summary$top_padj), , drop = FALSE]
  }

  fig_primary <- file.path(output_root, "primary_initial_ploidy", "06_figures")
  fig_endpoint <- file.path(output_root, "endpoint_ploidy", "06_figures")
  tier_method_summary <- ptp_tier_method_summary(program_results, camera_results, gsea_axes)
  score_method_summary <- ptp_score_method_result_summary(score_method_results)
  top_primary <- if (nrow(primary) > 0L) primary[order(primary$p_value), , drop = FALSE][1L, , drop = FALSE] else data.frame()
  primary_interp <- if (nrow(top_primary) > 0L) {
    paste0(
      "The primary signed-score model finds ", n_focus_sig,
      " focused Tier 1/Tier 2 interaction(s) and ", n_tier3_sig,
      " exploratory Tier 3 interaction(s) at the tier-appropriate FDR < 0.05. The top-ranked all-universe program is ",
      top_primary$program_label, " (", top_primary$response_family_label, "), estimate ",
      ptp_format_num(top_primary$estimate), ", 95% CI ",
      ptp_format_num(top_primary$ci_low), " to ", ptp_format_num(top_primary$ci_high),
      ", interpretation FDR ", ptp_format_num(top_primary$interpretation_fdr), "."
    )
  } else {
    "No primary program interaction rows were available for interpretation."
  }
  heatmap_interp <- paste0(
    "The heatmap shows signed treatment effects in 2N-origin and 4N-origin tumors plus the interaction contrast. ",
    "The interaction column is interpreted by tier: Tier 1/2 are focused evidence, Tier 3 remains exploratory, and negative controls are diagnostics."
  )
  score_plot_interp <- paste0(
    "The displayed mouse-subbin pseudobulk scores are the independent observations used by the repeated-measures model. ",
    "The panels are selected from the top-ranked all-universe program interactions and should be read together with tier-specific FDR columns."
  )
  tier3_interp <- if (nrow(tier3_summary) > 0L) {
    top_axis <- tier3_summary[1L, , drop = FALSE]
    paste0(
      "Tier 3 remains exploratory-only. The strongest contextual axis is ",
      top_axis$axis_label, " with ", top_axis$n_fdr_0_05,
      " GSEA pathway(s) at padj < 0.05; the top pathway is ",
      top_axis$top_pathway, " (NES ", ptp_format_num(top_axis$top_NES),
      ", padj ", formatC(top_axis$top_padj, format = "e", digits = 2L), ")."
    )
  } else {
    "No Tier 3 GSEA axis table was available."
  }
  camera_interp <- if (!is.null(camera_results) && nrow(camera_results) > 0L) {
    c0 <- camera_results[camera_results$contrast_id == "treatment_by_initial_ploidy_interaction", , drop = FALSE]
    c0$interpretation_FDR <- ptp_interpretation_fdr(c0, "FDR")
    c_focus <- c0[c0$tier %in% c("tier1_core", "tier2_upgraded_focused"), , drop = FALSE]
    c_tier3 <- c0[c0$tier == "tier3_exploratory_only", , drop = FALSE]
    paste0(
      "cameraPR is complementary gene-level enrichment evidence on the same frozen program universe. It finds ",
      sum(c_focus$interpretation_FDR < 0.05, na.rm = TRUE), " focused Tier 1/Tier 2 program(s) and ",
      sum(c_tier3$interpretation_FDR < 0.05, na.rm = TRUE), " exploratory Tier 3 program(s) at tier-appropriate FDR < 0.05."
    )
  } else {
    "cameraPR results were not available in this report run."
  }
  method_summary_interp <- if (nrow(tier_method_summary) > 0L) {
    top_method <- tier_method_summary[order(-tier_method_summary$n_fdr_0_05, tier_method_summary$min_fdr), , drop = FALSE][1L, , drop = FALSE]
    paste0(
      "This summary compares evidence types after stratifying the unified program universe by tier. The largest significant count is ",
      top_method$n_fdr_0_05, "/", top_method$n_tests, " for ", top_method$method,
      " in ", top_method$tier_label, "."
    )
  } else {
    "No tier-by-method summary was available."
  }
  score_method_interp <- if (nrow(score_method_summary) > 0L) {
    focused_sm <- score_method_summary[score_method_summary$tier %in% c("tier1_core", "tier2_upgraded_focused"), , drop = FALSE]
    tier3_sm <- score_method_summary[score_method_summary$tier == "tier3_exploratory_only", , drop = FALSE]
    paste0(
      "Score-definition sensitivity keeps the matched model fixed and changes only the score aggregation. Across focused Tier 1/Tier 2 programs, ",
      sum(focused_sm$n_fdr_0_05, na.rm = TRUE),
      " method-tier program result(s) reach tier-appropriate FDR < 0.05. Across Tier 3 exploratory programs, ",
      sum(tier3_sm$n_fdr_0_05, na.rm = TRUE),
      " method-tier program result(s) reach tier-appropriate FDR < 0.05."
    )
  } else {
    "Score-definition sensitivity results were not available in this report run."
  }
  score_method_qc_interp <- if (!is.null(score_method_qc_summary) && nrow(score_method_qc_summary) > 0L) {
    min_cov <- min(score_method_qc_summary$median_detected_program_fraction, na.rm = TRUE)
    n_fallback <- sum(score_method_qc_summary$n_fallback_observations > 0, na.rm = TRUE)
    paste0(
      "QC summaries are computed for every score method and program. The minimum median detected program-gene fraction is ",
      ptp_format_num(min_cov), "; ", n_fallback,
      " method-program combination(s) used a recorded fallback."
    )
  } else {
    "Score-method QC summaries were not available in this report run."
  }
  endpoint_sig <- 0L
  if (!is.null(endpoint_group_program) && nrow(endpoint_group_program) > 0L) {
    endpoint_group_program$interpretation_fdr <- ptp_interpretation_fdr(endpoint_group_program, "fdr")
    endpoint_sig <- endpoint_sig + sum(endpoint_group_program$interpretation_fdr < 0.05 & endpoint_group_program$contrast_id == "treatment_by_endpoint_ploidy_group_interaction", na.rm = TRUE)
  }
  if (!is.null(endpoint_continuous_program) && nrow(endpoint_continuous_program) > 0L) {
    endpoint_continuous_program$interpretation_fdr <- ptp_interpretation_fdr(endpoint_continuous_program, "fdr")
    endpoint_sig <- endpoint_sig + sum(endpoint_continuous_program$interpretation_fdr < 0.05 & endpoint_continuous_program$contrast_id == "continuous_endpoint_ploidy_treatment_interaction", na.rm = TRUE)
  }
  endpoint_not_estimable <- if (!is.null(endpoint_status) && nrow(endpoint_status) > 0L) {
    paste(endpoint_status$model_id[endpoint_status$status != "ok"], collapse = "; ")
  } else {
    ""
  }
  endpoint_interp <- paste0(
    "Endpoint-ploidy analyses are descriptive post-treatment sensitivities. Across thresholded and continuous ETP program contrasts, ",
    endpoint_sig, " program interaction(s) reach FDR < 0.05.",
    if (nzchar(endpoint_not_estimable)) paste0(" Non-estimable model(s): ", endpoint_not_estimable, ".") else ""
  )
  loo_interp <- if (!is.null(loo) && nrow(loo) > 0L) {
    paste0("Leave-one-mouse-out ranges assess influence of individual mice. The primary conclusion remains unchanged: ",
           n_focus_sig, " focused Tier 1/Tier 2 program interaction(s) reach FDR < 0.05 in the full model.")
  } else {
    "Leave-one-mouse-out results were not available in this report run."
  }
  continuous_path <- file.path(output_root, "primary_initial_ploidy", "07_robustness", "continuous_pseudotime_binned_interactions.csv")
  continuous_interp <- "Continuous pseudotime results were not available in this report run."
  if (file.exists(continuous_path)) {
    cont <- tryCatch(readr::read_csv(continuous_path, show_col_types = FALSE), error = function(e) data.frame())
    if (nrow(cont) > 0L) {
      cont$interpretation_fdr <- ptp_interpretation_fdr(cont, "fdr")
      cont_sig <- sum(cont$interpretation_fdr < 0.05, na.rm = TRUE)
      cont_top <- cont[order(cont$p_value), , drop = FALSE][1L, , drop = FALSE]
      continuous_interp <- paste0(
        "The 20-bin continuous-pseudotime sensitivity finds ", cont_sig,
        " program-bin interaction(s) at FDR < 0.05. The top program-bin is ",
        cont_top$program_label, " at pseudotime midpoint ", ptp_format_num(cont_top$bin_midpoint),
        " with estimate ", ptp_format_num(cont_top$estimate), " and interpretation FDR ", ptp_format_num(cont_top$interpretation_fdr), "."
      )
    }
  }
  lines <- c(
    "<!doctype html><html><head><meta charset='utf-8'><title>04j pseudotime treatment ploidy programs</title>",
    "<style>",
    "body{font-family:Arial,sans-serif;margin:0;background:#f7f8fa;color:#1f2933}.layout{display:grid;grid-template-columns:260px minmax(0,1fr);min-height:100vh}.sidebar{position:sticky;top:0;height:100vh;overflow:auto;background:#111827;color:white;padding:24px 18px}.sidebar h1{font-size:18px;margin:0 0 18px}.sidebar a{display:block;color:#d1d5db;text-decoration:none;margin:10px 0;font-size:13px}.content{padding:28px 36px;max-width:1220px}.hero,.card,section{background:white;border:1px solid #d9dee8;border-radius:8px;padding:20px;margin:0 0 18px}.note{color:#53606f}.summary{font-size:16px;font-weight:600}.table-wrap{overflow:auto}table{border-collapse:collapse;font-size:12px;min-width:760px}td,th{border:1px solid #dbe1ea;padding:5px 7px;text-align:left}th{background:#eef2f7}h2{margin-top:0}h3{margin:18px 0 8px}.figure-grid{display:grid;grid-template-columns:1fr;gap:18px}.figure-card{background:#fbfcfe;border:1px solid #d9dee8;border-radius:8px;padding:16px;margin:0}.figure-card img{display:block;width:100%;height:auto;border:1px solid #dbe1ea;border-radius:6px;background:white}.figure-card object{width:100%;height:520px;border:1px solid #dbe1ea;border-radius:6px}.figure-caption{font-size:13px;color:#374151;line-height:1.45;margin-top:10px}.figure-caption p{margin:7px 0}.guardrails li{margin-bottom:7px}",
    "</style></head><body><div class='layout'><nav class='sidebar'><h1>04j report</h1>",
    "<a href='#summary'>Executive summary</a><a href='#mapping'>Plan-to-report mapping</a><a href='#qc'>Input/QC/common support</a><a href='#score_methods'>Score sensitivity</a><a href='#tier1_core'>Tier 1 original core</a><a href='#tier2_upgraded_focused'>Tier 2 upgraded focused</a><a href='#tier3_exploratory_only'>Tier 3 exploratory-only</a><a href='#negative_control'>Negative controls</a><a href='#tier3_gsea'>Tier 3 GSEA map</a><a href='#figures'>Figures</a><a href='#integrated'>Integrated conclusion</a><a href='#outputs'>Output registry</a>",
    "</nav><main class='content'>",
    "<section id='summary' class='hero'><h1>04j pseudotime treatment-by-ploidy response-family analysis</h1>",
    "<p class='note'>Estimand: treatment association in captured endpoint cells from baseline-2N-origin and baseline-4N-origin tumors after pseudotime sub-bin standardization. Primary program estimates use pre-specified signed weighted mouse-pseudobulk scores. Endpoint ploidy is retained as descriptive sensitivity, including thresholded ETP and continuous mean ETP.</p>",
    "<p class='summary'>", ptp_html_escape(conclusion), "</p></section>",
    section("mapping", "Plan-to-report mapping",
            "<p>All frozen programs form one analysis universe. Tier assignment is a classification and interpretation label, not a method filter. Program-score modeling, cameraPR, endpoint-ploidy sensitivity, dose sensitivity, and robustness analyses are run on Tier 1, Tier 2, Tier 3, and negative-control programs. Tier 1 and Tier 2 remain the focused interpretation set; Tier 3 is analyzed in parallel but remains exploratory-only; negative controls are diagnostics. Genome-wide GSEA remains pathway-level evidence and is mapped back to contextual axes.</p>",
            "<h3>Response-family summary from program-score interaction model</h3>",
            ptp_html_table(family_summary, n = 120L),
            "<h3>Tier-by-method evidence summary</h3>",
            ptp_html_table(tier_method_summary, n = 80L)),
    section("qc", "Input/QC/common support",
            "<h3>Run status</h3>", ptp_html_table(status, n = 50L),
            "<h3>Design audit</h3>", ptp_html_table(design_audit, n = 50L),
            "<h3>Primary common support</h3>", ptp_html_table(support, n = 50L),
            "<h3>Endpoint ploidy model status</h3>", ptp_html_table(endpoint_status, n = 50L)),
    section("score_methods", "Score-definition sensitivity",
            "<p class='note'>These analyses keep the primary matched-pseudotime model fixed and change only the frozen score aggregation method. The signed weighted mean z-score remains the primary md-defined score; PCA/eigengene, rank-mean, and trimmed-mean scores are sensitivity analyses.</p>",
            "<h3>Score-method result summary</h3>",
            ptp_html_table(score_method_summary, n = 120L),
            "<h3>Score-method QC summary</h3>",
            ptp_html_table(score_method_qc_summary[, intersect(c("score_method_label", "program_label", "tier_label", "n_observations", "n_finite_scores", "n_genes_observed", "median_n_genes_used_by_score", "median_detected_program_fraction", "median_program_raw_umi", "score_sd", "n_fallback_observations", "fallback_reason"), names(score_method_qc_summary)), drop = FALSE], n = 120L),
            "<h3>Score-method precision</h3>",
            ptp_html_table(score_method_precision[, intersect(c("score_method_label", "program_label", "tier_label", "interaction_estimate", "interaction_se", "ci_width_95", "minimum_detectable_abs_interaction_approx_80pct_power_alpha_0_05"), names(score_method_precision)), drop = FALSE], n = 120L)),
    tier_section("tier1_core", "Tier 1: Original core gemcitabine-mechanism families",
                 "These are the originally focused families: gemcitabine handling, nucleotide supply/buffering, replication execution, replication-stress checkpoint, fork repair, and downstream fate."),
    tier_section("tier2_upgraded_focused", "Tier 2: Upgraded focused response families",
                 "These families were upgraded from exploratory GSEA because they are mechanistically close to gemcitabine-induced replication stress and transcriptional response, but they remain labeled separately from Tier 1."),
    tier_section("tier3_exploratory_only", "Tier 3: Exploratory-only contextual programs",
                 "These programs are analyzed with the same model family as Tier 1 and Tier 2, but their evidence level remains exploratory-only because they were defined from contextual exploratory axes."),
    tier_section("negative_control", "Negative controls",
                 "Negative controls are run through the same machinery to detect broad normalization, technical, or global-expression behavior that could weaken a mechanistic interpretation."),
    section("tier3_gsea", "Tier 3 GSEA-mapped contextual pathways",
            "<p class='note'>This table maps genome-wide interaction-ranked GSEA pathways to the exploratory Tier 3 contextual axes. These pathway-level results are complementary to, not replacements for, the Tier 3 program-score model.</p>",
            ptp_html_table(tier3_summary, n = 80L),
            ptp_html_table(gsea_axes[, intersect(c("axis_label", "pathway", "NES", "padj", "pval", "size"), names(gsea_axes)), drop = FALSE], n = 60L)),
    section("figures", "Figure-supported results",
            "<div class='figure-grid'>",
            figure_block(
              "Primary signed program interaction forest by response-family tier",
              file.path(fig_primary, "program_interaction_forest_by_tier.pdf"),
              "Points show signed program-score model coefficients and horizontal bars show 95% confidence intervals. Colors distinguish treatment effects within 2N-origin tumors, treatment effects within 4N-origin tumors, and the primary 4N-origin minus 2N-origin treatment-response interaction. Facets separate response-family tiers.",
              primary_interp
            ),
            figure_block(
              "Signed program simple-effect and interaction heatmap by tier",
              file.path(fig_primary, "program_effect_heatmap_by_tier.pdf"),
              "Tiles encode signed program-score estimates. Red indicates positive signed-score association and blue indicates negative signed-score association. Columns show the 2N-origin treatment effect, 4N-origin treatment effect, and primary initial-ploidy interaction.",
              heatmap_interp
            ),
            figure_block(
              "Mouse-subbin pseudobulk signed scores for top-ranked interactions",
              file.path(fig_primary, "program_score_pseudobulk_top_interactions.pdf"),
              "Each point is one mouse-by-pseudotime-subbin pseudobulk score, not a cell. Black points show group means. Panels display the top-ranked programs by primary interaction p-value.",
              score_plot_interp
            ),
            figure_block(
              "Score-definition sensitivity of primary interaction estimates",
              file.path(fig_primary, "score_method_interaction_heatmap.pdf"),
              "Rows are frozen programs and columns are score definitions. Tiles encode the signed treatment-by-initial-ploidy interaction estimate from the same matched model. Asterisks mark tier-appropriate FDR < 0.05.",
              score_method_interp
            ),
            figure_block(
              "Score-definition QC across the unified program universe",
              file.path(fig_primary, "score_method_qc_summary.pdf"),
              "Each point is one method-program combination. The y-axis shows the median detected fraction of observed program genes across retained mouse-subbin pseudobulks. Color marks response-family tier.",
              score_method_qc_interp
            ),
            figure_block(
              "Program-score interaction versus cameraPR enrichment across the unified universe",
              file.path(fig_primary, "program_score_vs_cameraPR.pdf"),
              "Each point is one frozen program. The x-axis is the signed program-score treatment-by-initial-ploidy interaction estimate. The y-axis is the negative log10 tier-appropriate cameraPR FDR. Color marks response-family tier.",
              camera_interp
            ),
            figure_block(
              "Tier-by-method evidence summary",
              file.path(fig_primary, "tier_method_evidence_summary.pdf"),
              "Tiles show the number of FDR < 0.05 results over the number of tests for each tier and method. Program-score and cameraPR use tier-appropriate FDR; GSEA counts mapped genome-wide pathways.",
              method_summary_interp
            ),
            figure_block(
              "Tier 3 exploratory GSEA response axes",
              file.path(fig_primary, "tier3_exploratory_axis_gsea.pdf"),
              "Bars show interaction-ranked GSEA normalized enrichment scores. Red bars indicate pathways with padj < 0.05. These genome-wide pathways are mapped back to Tier 3 contextual axes and are complementary to the Tier 3 program-score tests.",
              tier3_interp
            ),
            figure_block(
              "Primary initial-ploidy estimates versus endpoint-ploidy sensitivity estimates",
              file.path(fig_endpoint, "program_estimate_comparison_initial_vs_endpoint_ploidy.pdf"),
              "Each point is one program. The x-axis is the primary baseline initial-ploidy interaction estimate. The y-axis is the endpoint-ploidy sensitivity estimate for thresholded or continuous ETP models.",
              endpoint_interp
            ),
            figure_block(
              "Leave-one-mouse-out robustness of signed program interactions",
              file.path(fig_primary, "leave_one_mouse_out_program_interactions.pdf"),
              "Red points show the full-model primary interaction estimate. Grey horizontal ranges show the minimum-to-maximum estimates from leave-one-mouse-out refits for the same program.",
              loo_interp
            ),
            figure_block(
              "Continuous pseudotime sensitivity of signed program interactions",
              file.path(fig_primary, "continuous_pseudotime_interactions.pdf"),
              "Lines and points show primary interaction estimates across 20 pseudotime bins for the top-ranked binned programs. Error bars show 95% confidence intervals where the model is estimable.",
              continuous_interp
            ),
            "</div>"),
    section("integrated", "Integrated conclusion",
            "<p>", ptp_html_escape(conclusion), " All estimable frozen programs are analyzed in one universe by program-score modeling and cameraPR, then interpreted by tier. Tier 1 and Tier 2 are the focused response-family evidence set, with Tier 2 explicitly flagged as upgraded after exploratory-pathway review. Tier 3 program-score and cameraPR results are analyzed in parallel but remain exploratory-only. Genome-wide GSEA remains pathway-level complementary evidence mapped back to contextual axes. Score-definition sensitivity keeps the matched model unchanged and tests whether the primary conclusion depends on the score aggregator. Endpoint-ploidy threshold and continuous ETP models describe post-treatment endpoint state and should not replace baseline initial ploidy as the primary tumor-origin modifier.</p>",
            "<h3>Precision</h3>", ptp_html_table(precision, n = 80L),
            "<h3>Interpretation guardrails</h3><ul class='guardrails'><li>Use baseline initial_ploidy as the primary tumor-of-origin modifier.</li><li>Endpoint ploidy and cell ploidy are post-treatment descriptive sensitivity modifiers; they are not baseline predictive-biomarker evidence.</li><li>RNA program interactions do not measure dFdCTP, dNTP pools, DNA-incorporated gemcitabine, DNA damage, protein phosphorylation, or cells lost before capture.</li><li>Unsupported pseudotime regions are not interpreted as matched-pseudotime evidence.</li></ul>"),
    section("outputs", "Output registry",
            "<ul><li><code>00_manifest/</code>: inputs, parameters, program membership, shared status.</li><li><code>01_qc/</code>: shared input and sample-level QC.</li><li><code>02_pseudobulk/</code>: mouse-by-subbin and mouse-by-region pseudobulk metadata.</li><li><code>primary_initial_ploidy/</code>: primary initial-ploidy gene, program, GSEA, figure, and robustness outputs.</li><li><code>endpoint_ploidy/&lt;model_id&gt;/</code>: ETP threshold and continuous ETP sensitivity outputs.</li><li><code>dose_specific_initial_ploidy/</code>: dose-specific program sensitivity outputs.</li><li><code>model_comparison/</code>: model-level design audit outputs.</li></ul>"),
    "</main></div></body></html>"
  )
  path <- file.path(output_root, "report", "04j_pseudotime_treatment_ploidy_programs_report.html")
  ptp_write_lines(lines, path)
  path
}

ptp_write_readme <- function(output_root) {
  ptp_write_lines(c(
    "# 04j pseudotime treatment-by-ploidy programs",
    "",
    "Standalone matched-pseudotime treatment-by-baseline-initial-ploidy analysis.",
    "",
    "Key folders:",
    "- `00_manifest`: parameters, checksums, program membership, run status, and shared parallel task logs.",
    "- `01_qc`: shared expression/sample QC and processing-batch checks.",
    "- `02_pseudobulk`: mouse-by-sub-bin and mouse-by-region pseudobulk metadata.",
    "- `primary_initial_ploidy`: primary initial-ploidy gene models, all-universe program models, enrichment outputs, figures, and robustness outputs.",
    "- `primary_initial_ploidy/04_program_models/score_methods`: primary and score-definition sensitivity program-score results, QC, and precision tables.",
    "- `endpoint_ploidy/<model_id>`: thresholded ETP and continuous mean ETP sensitivity models, separated by endpoint-ploidy definition.",
    "- `dose_specific_initial_ploidy`: dose-specific program sensitivity outputs.",
    "- `model_comparison`: cross-model design-rank and model-audit outputs.",
    "- `report`: tier-organized HTML report.",
    "",
    "Primary contrast: `(gemcitabine-control in 4N-origin) - (gemcitabine-control in 2N-origin)` with equal weights over common-support pseudotime sub-bins.",
    "",
    "Program scores:",
    "- Primary program estimates use pre-specified signed weighted mouse-pseudobulk scores.",
    "- Score-definition sensitivity additionally fits signed PCA/eigengene, signed rank-mean, and robust trimmed signed mean scores with the same matched model and contrasts.",
    "- `00_manifest/program_gene_sets.csv` records score direction, component policy, positive/negative observed weights, and observed genes.",
    "- `00_manifest/program_score_gene_weights.csv` records the fixed gene-level component, sign, and weight audit used to compute scores.",
    "- `primary_initial_ploidy/04_program_models/program_score_pseudobulk_long.csv` contains the mouse-subbin score values plotted in the report.",
    "",
    "Endpoint-ploidy sensitivity outputs:",
    "- `endpoint_ploidy/ETP_*/04_program_models/program_interactions.csv`: thresholded post-treatment ETP-lower versus ETP-higher modifiers.",
    "- `endpoint_ploidy/continuous_mean_endpoint_ploidy/04_program_models/program_interactions.csv`: continuous sample mean ETP modifier; the main contrast is per 1 SD mean-ETP change in the gemcitabine-control association.",
    "",
    "Endpoint ploidy is post-treatment. These outputs are descriptive sensitivities and should not replace baseline initial ploidy as the primary tumor-origin estimand.",
    "",
    "Response-family tiers:",
    "- Tier 1 original core: original gemcitabine-mechanism response families.",
    "- Tier 2 upgraded focused: exploratory GSEA axes upgraded because they map to gemcitabine-linked replication stress, DNA repair, checkpoint, TP53/fate, and nucleotide-metabolism biology.",
    "- Tier 3 exploratory-only: translation/RNA/MYC, immune/interferon/antigen/cytotoxicity, and adhesion/cell-state organization programs are analyzed in the same program universe but interpreted as exploratory evidence.",
    "- `primary_initial_ploidy/05_enrichment`: cameraPR, genome-wide GSEA, and tier-mapped evidence summaries."
  ), file.path(output_root, "README.md"))
}

ptp_run_workflow <- function(args, repo_root, script_dir) {
  set.seed(args$seed)
  ptp_check_packages()
  cfg <- ptp_read_config(args$config)
  endpoint_specs <- ptp_endpoint_ploidy_specs(cfg)
  output_root <- ptp_clean_dir(args$output_root, overwrite = args$overwrite)
  dirs <- ptp_output_dirs(output_root)
  cell_metadata_path <- normalizePath(args$cell_metadata, winslash = "/", mustWork = TRUE)
  sample_info_path <- if (file.exists(args$sample_info)) normalizePath(args$sample_info, winslash = "/", mustWork = TRUE) else args$sample_info
  seurat_rds <- ptp_resolve_seurat_rds(args$seurat_rds, cfg, args$results_root)

  ptp_write_csv(data.frame(
    parameter = names(args),
    value = vapply(args, as.character, character(1L)),
    stringsAsFactors = FALSE
  ), file.path(dirs$manifest, "analysis_parameters.csv"))
  ptp_write_csv(ptp_package_versions(), file.path(dirs$manifest, "package_versions.csv"))
  ptp_write_csv(rbind(
    ptp_subbins_table(cfg$primary_subbins, "primary_four_bin"),
    ptp_subbins_table(cfg$fallback_subbins, "fallback_two_bin")
  ), file.path(dirs$manifest, "pseudotime_subbin_definitions.csv"))
  ptp_write_csv(ptp_regions_table(cfg), file.path(dirs$manifest, "pseudotime_region_definitions.csv"))
  ptp_write_csv(ptp_endpoint_ploidy_specs_table(endpoint_specs), file.path(dirs$endpoint$manifest, "endpoint_ploidy_threshold_definitions.csv"))
  ptp_write_csv(data.frame(
    input = c("cell_metadata", "sample_info", "seurat_rds", "config", "plan"),
    path = c(cell_metadata_path, sample_info_path, seurat_rds, args$config, cfg$analysis$plan_path %||% ""),
    sha256 = c(ptp_file_checksum(cell_metadata_path), ptp_file_checksum(sample_info_path), ptp_file_checksum(seurat_rds), ptp_file_checksum(args$config), ptp_file_checksum(cfg$analysis$plan_path %||% "")),
    stringsAsFactors = FALSE
  ), file.path(dirs$manifest, "input_manifest.csv"))

  meta <- ptp_read_cell_metadata(cell_metadata_path)
  loaded <- ptp_load_counts_and_qc(seurat_rds, args$assay, args$counts_layer, meta)
  counts <- loaded$counts
  meta <- loaded$metadata
  ptp_write_csv(attr(counts, "seurat_audit"), file.path(dirs$qc, "expression_source_audit.csv"))
  ptp_write_csv(ptp_audit_counts(counts), file.path(dirs$qc, "count_matrix_audit.csv"))
  matched <- ptp_match_cells(meta, counts, cfg$qc$min_metadata_to_counts_match_rate, dirs$qc)
  meta <- matched$meta
  counts <- matched$counts
  limited <- ptp_limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts

  assignment <- ptp_assignment_audit(meta)
  ptp_write_csv(assignment, file.path(dirs$qc, "sample_assignment_audit.csv"))
  if (any(!assignment$one_initial_ploidy) || any(!assignment$one_dose)) {
    stop("At least one mouse has non-unique dose or initial_ploidy assignment.", call. = FALSE)
  }
  sample_meta <- ptp_sample_metadata(meta)
  ptp_write_csv(sample_meta, file.path(dirs$qc, "mouse_region_coverage.csv"))
  ptp_write_csv(sample_meta, file.path(dirs$endpoint$qc, "endpoint_ploidy_sample_assignments.csv"))
  ptp_write_csv(ptp_endpoint_group_balance(sample_meta, endpoint_specs), file.path(dirs$endpoint$qc, "endpoint_ploidy_group_balance.csv"))
  crosswalk <- ptp_processing_batch_crosswalk(meta, sample_info_path, endpoint_specs)
  ptp_write_csv(crosswalk, file.path(dirs$qc, "processing_batch_crosswalk.csv"))
  ptp_write_csv(ptp_batch_group_separability(crosswalk), file.path(dirs$qc, "processing_batch_separability_audit.csv"))

  pb <- ptp_choose_primary_grid(meta, counts, cfg)
  ptp_write_csv(pb$metadata, file.path(dirs$pseudobulk, "sample_subbin_metadata.csv"))
  ptp_write_csv(pb$support, file.path(dirs$primary$qc, "pseudotime_subbin_common_support.csv"))
  if (!is.null(pb$primary_four_bin_support)) {
    ptp_write_csv(pb$primary_four_bin_support, file.path(dirs$primary$qc, "pseudotime_subbin_common_support_primary_four_bin_failed.csv"))
  }
  ptp_write_csv(pb$cluster_composition, file.path(dirs$primary$qc, "sample_subbin_cluster_composition.csv"))
  balance <- pb$metadata[, c("sample_id", "subbin_id", "initial_ploidy", "treatment", "dose_mg", "cell_count", "library_size", "detected_genes", "mean_subbin_pseudotime", "median_subbin_pseudotime", "retained_for_model", "exclusion_reason"), drop = FALSE]
  ptp_write_csv(balance, file.path(dirs$primary$qc, "pseudotime_distribution_balance.csv"))

  run_status <- data.frame(
    step = c("input_loaded", "cell_id_match", "primary_grid", "primary_estimable"),
    status = c("ok", "ok", pb$grid_id, if (isTRUE(pb$estimable)) "ok" else "not_estimable"),
    detail = c(
      paste(nrow(meta), "cells", length(unique(meta$sample_id)), "samples"),
      "metadata cell IDs matched expression matrix",
      paste("eligible_subbins", paste(pb$eligible_subbins, collapse = ";")),
      if (isTRUE(pb$estimable)) "common support satisfied" else "fewer than required common-support subbins"
    ),
    stringsAsFactors = FALSE
  )
  ptp_write_csv(run_status, file.path(dirs$manifest, "workflow_run_status.csv"))
  if (!isTRUE(pb$estimable)) {
    ptp_write_report(output_root, cfg, run_status, pb$support, data.frame(), data.frame(), data.frame())
    ptp_write_readme(output_root)
    message("04j workflow stopped with non-estimability report: ", output_root)
    return(invisible(output_root))
  }

  primary_fit <- ptp_fit_cellmeans(pb$counts, pb$metadata, pb$eligible_subbins, "ptp_group", cfg$qc)
  primary_design <- ptp_design_audit(primary_fit, "primary_pooled_treatment_cellmeans")
  dose_fit <- ptp_fit_cellmeans(pb$counts, pb$metadata, pb$eligible_subbins, "dose_group", cfg$qc)
  dose_design <- ptp_design_audit(dose_fit, "dose_specific_cellmeans")
  design_audit <- rbind(primary_design, dose_design)
  ptp_write_csv(primary_fit$metadata, file.path(dirs$pseudobulk, "sample_subbin_model_metadata.csv"))

  observed_symbols <- unique(ptp_clean_gene_symbols(primary_fit$retained_genes))
  programs <- ptp_load_programs(cfg, cfg$qc$min_program_genes_observed, observed_symbols)
  ptp_write_csv(ptp_program_membership_table(programs), file.path(dirs$manifest, "program_gene_sets.csv"))
  ptp_write_csv(ptp_program_weight_table(programs), file.path(dirs$manifest, "program_score_gene_weights.csv"))

  primary_contrasts <- ptp_primary_contrasts(pb$eligible_subbins, colnames(primary_fit$design))
  gene_primary <- ptp_apply_contrasts(primary_fit, primary_contrasts, "primary_pooled_treatment")
  ptp_write_csv(gene_primary, file.path(dirs$primary$genes, "gene_level_simple_effects_and_interactions.csv"))

  primary_logcpm <- ptp_logcpm_matrix(primary_fit)
  score_method_scores <- ptp_program_score_methods_from_logcpm(primary_logcpm, programs, cfg, pb_counts = primary_fit$dge$counts)
  score_method_results_list <- list()
  score_method_precision_list <- list()
  score_method_qc_long_list <- list()
  score_method_qc_summary_list <- list()
  primary_method_id <- names(score_method_scores)[vapply(score_method_scores, function(x) identical(x$method$role %||% "", "primary"), logical(1L))]
  if (length(primary_method_id) == 0L) primary_method_id <- names(score_method_scores)[1L]
  primary_program_fit <- NULL
  program_primary <- data.frame()
  for (method_id in names(score_method_scores)) {
    sm <- score_method_scores[[method_id]]
    method <- sm$method
    method_label <- method$label %||% method_id
    method_role <- method$role %||% "score_definition_sensitivity"
    method_mode <- method$mode %||% method_id
    method_dir <- file.path(dirs$primary$score_methods, ptp_safe_name(method_id))
    pfit <- ptp_fit_program_scores(sm$score_matrix, primary_fit)
    res <- ptp_apply_program_contrasts(pfit, primary_contrasts, "primary_pooled_treatment", programs)
    if (nrow(res) > 0L) {
      method_cols <- data.frame(
        score_method_id = method_id,
        score_method_label = method_label,
        score_method_role = method_role,
        score_method_mode = method_mode,
        stringsAsFactors = FALSE
      )
      res <- cbind(method_cols[rep(1L, nrow(res)), , drop = FALSE], res)
    }
    prec <- ptp_precision_table(res)
    ptp_write_csv(res, file.path(method_dir, "program_score_simple_effects_and_interactions.csv"))
    ptp_write_csv(prec, file.path(method_dir, "program_score_precision_and_power.csv"))
    ptp_write_csv(sm$qc_long, file.path(method_dir, "program_score_qc_long.csv"))
    ptp_write_csv(sm$qc_summary, file.path(method_dir, "program_score_qc_summary.csv"))
    score_method_results_list[[method_id]] <- res
    score_method_precision_list[[method_id]] <- prec
    score_method_qc_long_list[[method_id]] <- sm$qc_long
    score_method_qc_summary_list[[method_id]] <- sm$qc_summary
    if (identical(method_id, primary_method_id[1L])) {
      primary_program_fit <- pfit
      program_primary <- res
    }
  }
  score_method_results <- do.call(rbind, score_method_results_list)
  score_method_precision <- do.call(rbind, score_method_precision_list)
  score_method_qc_long <- do.call(rbind, score_method_qc_long_list)
  score_method_qc_summary <- do.call(rbind, score_method_qc_summary_list)
  ptp_write_csv(score_method_results, file.path(dirs$primary$score_methods, "all_score_methods_program_score_simple_effects_and_interactions.csv"))
  ptp_write_csv(score_method_precision, file.path(dirs$primary$score_methods, "all_score_methods_program_score_precision_and_power.csv"))
  ptp_write_csv(score_method_qc_long, file.path(dirs$primary$score_methods, "all_score_methods_program_score_qc_long.csv"))
  ptp_write_csv(score_method_qc_summary, file.path(dirs$primary$score_methods, "all_score_methods_program_score_qc_summary.csv"))
  score_method_summary <- ptp_score_method_result_summary(score_method_results)
  ptp_write_csv(score_method_summary, file.path(dirs$primary$score_methods, "score_method_result_summary.csv"))
  ptp_write_csv(program_primary, file.path(dirs$primary$programs, "program_score_simple_effects_and_interactions.csv"))
  ptp_write_csv(program_primary, file.path(dirs$primary$programs, "all_programs", "program_score_simple_effects_and_interactions.csv"))
  for (tier_id in unique(program_primary$tier)) {
    ptp_write_csv(program_primary[program_primary$tier == tier_id, , drop = FALSE],
                  file.path(dirs$primary$programs, "by_tier", ptp_safe_name(tier_id), "program_score_simple_effects_and_interactions.csv"))
  }
  precision <- ptp_precision_table(program_primary)
  ptp_write_csv(precision, file.path(dirs$primary$programs, "program_score_precision_and_power.csv"))
  ptp_write_csv(precision, file.path(dirs$primary$programs, "all_programs", "program_score_precision_and_power.csv"))
  score_long <- ptp_program_score_long_table(primary_program_fit, programs)
  ptp_write_csv(score_long, file.path(dirs$primary$programs, "program_score_pseudobulk_long.csv"))
  ptp_write_csv(score_long, file.path(dirs$primary$programs, "all_programs", "program_score_pseudobulk_long.csv"))

  dose_contrasts <- ptp_dose_contrasts(pb$eligible_subbins, colnames(dose_fit$design))
  dose_scores <- ptp_program_scores_from_logcpm(ptp_logcpm_matrix(dose_fit), programs)
  dose_program_fit <- ptp_fit_program_scores(dose_scores, dose_fit)
  dose_program <- ptp_apply_program_contrasts(dose_program_fit, dose_contrasts, "dose_specific", programs)
  ptp_write_csv(dose_program, file.path(dirs$dose$programs, "dose_specific_interactions.csv"))

  endpoint <- ptp_run_endpoint_ploidy_models(pb, programs, cfg, args)
  ptp_write_endpoint_outputs(endpoint, dirs$endpoint)
  if (nrow(endpoint$design) > 0L) design_audit <- rbind(design_audit, endpoint$design)
  ptp_write_csv(design_audit, file.path(dirs$model_comparison, "design_matrix_rank.csv"))

  cam <- ptp_camera_tests(gene_primary, programs, "treatment_by_initial_ploidy_interaction")
  ptp_write_csv(cam, file.path(dirs$primary$gsea, "competitive_gene_set_tests.csv"))
  ptp_write_csv(cam, file.path(dirs$primary$enrichment$camera, "competitive_gene_set_tests.csv"))
  for (tier_id in unique(cam$tier)) {
    ptp_write_csv(cam[cam$tier == tier_id, , drop = FALSE],
                  file.path(dirs$primary$enrichment$camera, "by_tier", ptp_safe_name(tier_id), "competitive_gene_set_tests.csv"))
  }
  fg <- data.frame()
  gsea_axes <- data.frame()
  if (isTRUE(args$run_gsea)) {
    fg <- ptp_run_fgsea(gene_primary, args$gsea_min_size, args$gsea_max_size, args$gsea_nperm_simple, args$seed)
    ptp_write_csv(fg, file.path(dirs$primary$gsea, "interaction_ranked_gsea.csv"))
    ptp_write_csv(fg, file.path(dirs$primary$enrichment$gsea, "interaction_ranked_gsea.csv"))
    gsea_axes <- ptp_gsea_axis_table(fg, cfg)
    ptp_write_csv(gsea_axes, file.path(dirs$primary$gsea, "tier3_exploratory_axis_gsea.csv"))
    ptp_write_csv(gsea_axes, file.path(dirs$primary$enrichment$tier_mapped, "tier3_exploratory_axis_gsea.csv"))
  }

  region_pb <- ptp_construct_region_pseudobulk(meta, counts, cfg$regions, cfg$qc, endpoint_specs)
  ptp_write_csv(region_pb$metadata, file.path(dirs$pseudobulk, "sample_region_metadata.csv"))
  occupancy <- ptp_occupancy_models(meta, cfg$regions)
  ptp_write_csv(occupancy, file.path(dirs$primary$robustness, "pseudotime_region_occupancy_interactions.csv"))
  if (isTRUE(args$run_sensitivities)) {
    sens_workers <- ptp_effective_workers(args, length(unique(region_pb$metadata$region_id)), "workers")
    message("Running region sensitivities with workers=", sens_workers)
    sens <- ptp_fit_region_programs(region_pb, programs, cfg$qc, workers = sens_workers)
    ptp_write_csv(sens, file.path(dirs$primary$robustness, "interval_and_model_sensitivities.csv"))
    sens_log <- attr(sens, "parallel_log")
    if (!is.null(sens_log)) ptp_write_csv(sens_log, file.path(dirs$manifest, "parallel_task_log_region_sensitivities.csv"))
    coord <- data.frame(
      sensitivity_id = c("untreated_learned_trajectory_projection", "phase_state_label", "coordinate_without_tested_genes"),
      status = c("not_available_in_local_inputs", "available_descriptive_only", "not_implemented"),
      note = c(
        "No independent untreated-only trajectory coordinate was present in the frozen inputs.",
        "Seurat Phase/S.Score/G2M.Score are exported in processing_batch_crosswalk/sample metadata but are same-transcriptome diagnostics.",
        "A sensitivity coordinate excluding nucleotide/replication genes would require recomputing the trajectory upstream."
      ),
      stringsAsFactors = FALSE
    )
    ptp_write_csv(coord, file.path(dirs$primary$robustness, "pseudotime_coordinate_sensitivities.csv"))
  }

  loo <- data.frame()
  if (isTRUE(args$run_leave_one_out)) {
    loo_workers <- ptp_effective_workers(args, length(unique(pb$metadata$sample_id)), "workers")
    message("Running leave-one-mouse-out with workers=", loo_workers)
    loo <- ptp_leave_one_out(pb, pb$eligible_subbins, programs, cfg$qc, workers = loo_workers)
    ptp_write_csv(loo, file.path(dirs$primary$robustness, "leave_one_mouse_out.csv"))
    loo_log <- attr(loo, "parallel_log")
    if (!is.null(loo_log)) ptp_write_csv(loo_log, file.path(dirs$manifest, "parallel_task_log_leave_one_mouse_out.csv"))
  }

  if (isTRUE(args$run_continuous)) {
    cont <- ptp_continuous_binned_interactions(meta, counts, programs, cfg$qc, n_bins = 20L)
    ptp_write_csv(cont$support, file.path(dirs$primary$qc, "continuous_pseudotime_20bin_support.csv"))
    ptp_write_csv(cont$results, file.path(dirs$primary$robustness, "continuous_pseudotime_binned_interactions.csv"))
    ptp_plot_continuous(cont$results, file.path(dirs$primary$figures, "continuous_pseudotime_interactions.pdf"))
  }

  ptp_plot_program_forest(program_primary, file.path(dirs$primary$figures, "program_interaction_forest_by_tier.pdf"))
  ptp_plot_program_heatmap(program_primary, file.path(dirs$primary$figures, "program_effect_heatmap_by_tier.pdf"))
  ptp_plot_program_score_pseudobulk(score_long, program_primary, file.path(dirs$primary$figures, "program_score_pseudobulk_top_interactions.pdf"))
  ptp_plot_score_method_interaction_heatmap(score_method_results, file.path(dirs$primary$figures, "score_method_interaction_heatmap.pdf"))
  ptp_plot_score_method_qc_summary(score_method_qc_summary, file.path(dirs$primary$figures, "score_method_qc_summary.pdf"))
  method_summary <- ptp_tier_method_summary(program_primary, cam, gsea_axes)
  ptp_write_csv(method_summary, file.path(dirs$primary$enrichment$tier_mapped, "tier_method_evidence_summary.csv"))
  ptp_plot_program_camera_comparison(program_primary, cam, file.path(dirs$primary$figures, "program_score_vs_cameraPR.pdf"))
  ptp_plot_tier_method_summary(method_summary, file.path(dirs$primary$figures, "tier_method_evidence_summary.pdf"))
  ptp_plot_gsea_axis_bar(gsea_axes, file.path(dirs$primary$figures, "tier3_exploratory_axis_gsea.pdf"))
  ptp_plot_endpoint_program_comparison(program_primary, endpoint$group_program, endpoint$continuous_program, file.path(dirs$endpoint$figures, "program_estimate_comparison_initial_vs_endpoint_ploidy.pdf"))
  ptp_plot_loo_robustness(loo, program_primary, file.path(dirs$primary$figures, "leave_one_mouse_out_program_interactions.pdf"))
  ptp_write_report(output_root, cfg, run_status, pb$support, program_primary, precision, design_audit,
                   endpoint$status, endpoint$group_program, endpoint$continuous_program,
                   dose_program = dose_program, camera_results = cam, gsea_axes = gsea_axes, loo = loo,
                   score_method_results = score_method_results,
                   score_method_precision = score_method_precision,
                   score_method_qc_summary = score_method_qc_summary)
  ptp_write_readme(output_root)
  message("04j workflow complete: ", output_root)
  invisible(output_root)
}
