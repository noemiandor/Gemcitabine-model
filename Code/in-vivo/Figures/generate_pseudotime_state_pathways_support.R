#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

parse_cli_args <- function(argv, defaults) {
  out <- defaults
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!grepl("^--", token)) {
      i <- i + 1L
      next
    }
    key_value <- sub("^--", "", token)
    if (grepl("=", key_value, fixed = TRUE)) {
      key <- sub("=.*$", "", key_value)
      value <- sub("^[^=]*=", "", key_value)
      i <- i + 1L
    } else {
      key <- key_value
      if (i < length(argv) && !grepl("^--", argv[[i + 1L]])) {
        value <- argv[[i + 1L]]
        i <- i + 2L
      } else {
        value <- "TRUE"
        i <- i + 1L
      }
    }
    key <- gsub("-", "_", key)
    if (!key %in% names(out)) stop("Unknown argument: --", key, call. = FALSE)
    out[[key]] <- parse_scalar(value, out[[key]])
  }
  out
}

parse_scalar <- function(value, template) {
  if (is.logical(template)) return(tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y"))
  if (is.integer(template)) return(as.integer(value))
  if (is.numeric(template)) return(as.numeric(value))
  as.character(value)
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript Code/in-vivo/Figures/generate_pseudotime_state_pathways_support.R \\",
      "    --output_root Data/in-vivo/pseudotime_state_pathways \\",
      "    --seurat_rds /path/to/integrated_sct_cca_seurat_final_reclustered.rds",
      "",
      "Default input files:",
      "  Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "  Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "  Data/in-vivo/scRNA_Seq_Data/seurat_obj_annotated/integrated_sct_cca_seurat_final_reclustered.rds",
      "  Code/in-vivo/04i_pseudotime_state_pathways_config.yaml",
      "",
      "This standalone script writes only the Figure 7 support files requested under:",
      "  <output_root>/00_manifest",
      "  <output_root>/binning",
      "  <output_root>/binning/ETP_reference_balanced_threshold_2_24",
      sep = "\n"
    ),
    "\n"
  )
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0L) return(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE))
  NA_character_
}

is_abs_path <- function(path) {
  grepl("^/", path)
}

resolve_path <- function(path, repo_root, must_work = FALSE) {
  if (is.null(path) || length(path) == 0L || is.na(path) || !nzchar(path)) return("")
  out <- if (is_abs_path(path)) path else file.path(repo_root, path)
  normalizePath(out, mustWork = must_work)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, mustWork = TRUE)
}

clean_dir <- function(path, overwrite = FALSE) {
  if (dir.exists(path) && !overwrite) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) > 0L) {
      stop("Output directory is not empty. Use --overwrite=TRUE: ", path, call. = FALSE)
    }
  }
  if (dir.exists(path) && overwrite) unlink(path, recursive = TRUE, force = TRUE)
  ensure_dir(path)
}

write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  readr::write_csv(as.data.frame(x, stringsAsFactors = FALSE), path, na = "")
  invisible(path)
}

required_packages <- function() {
  c(
    "Seurat", "Matrix", "yaml", "digest", "edgeR", "limma", "splines",
    "fgsea", "msigdbr", "dplyr", "readr", "tidyr", "ggplot2"
  )
}

check_packages <- function(packages = required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

package_versions <- function(packages = required_packages()) {
  data.frame(
    package = packages,
    version = vapply(packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1L)),
    stringsAsFactors = FALSE
  )
}

file_checksum <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

safe_name <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

collection_label <- function(collection) {
  out <- c("H" = "hallmark", "C2:CP:REACTOME" = "reactome", "C5:GO:BP" = "go_bp")
  unname(out[collection] %||% safe_name(collection))
}

normalize_choice <- function(x) {
  gsub("-", "_", tolower(trimws(as.character(x))))
}

format_pathway_label <- function(pathway, collection = NA_character_) {
  x <- as.character(pathway)
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^GO_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

clean_gene_symbols <- function(genes) {
  g <- as.character(genes)
  g <- trimws(g)
  g <- sub("^GRCh[0-9]+[-_]", "", g, ignore.case = TRUE)
  g <- sub("^GRCm39[-_]", "", g, ignore.case = TRUE)
  g <- sub("^hg38[-_]", "", g, ignore.case = TRUE)
  g <- sub("\\.[0-9]+$", "", g)
  g[g == ""] <- NA_character_
  g
}

get_assay_data_slot <- function(obj, assay = "RNA", slot_name = "counts") {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = slot_name),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay, layer = slot_name),
        error = function(e2) {
          tryCatch(
            SeuratObject::LayerData(obj[[assay]], layer = slot_name),
            error = function(e3) NULL
          )
        }
      )
    }
  )
}

get_assay_matrix <- function(obj, assay = "RNA", slot_name = "counts") {
  get_assay_data_slot(obj, assay = assay, slot_name = slot_name)
}

etp_threshold_specs <- function() {
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

model_spec_etp_reference_balanced <- function() {
  spec <- etp_threshold_specs()[["ETP_reference_balanced_threshold_2_24"]]
  list(
    model_id = spec$method,
    model_label = paste0("ETP group-adjusted model: ", spec$threshold_scheme),
    covariate_terms = c(spec$factor_column),
    covariate_mode = "etp_group",
    etp_method = spec$method,
    etp_threshold = spec$threshold,
    is_primary = FALSE
  )
}

model_parameters_table <- function(model_spec) {
  data.frame(
    parameter = c("model_id", "model_label", "covariate_mode", "covariate_terms", "etp_method", "etp_threshold", "is_primary"),
    value = c(
      model_spec$model_id %||% NA_character_,
      model_spec$model_label %||% NA_character_,
      model_spec$covariate_mode %||% NA_character_,
      paste(model_spec$covariate_terms %||% character(), collapse = ";"),
      model_spec$etp_method %||% NA_character_,
      as.character(model_spec$etp_threshold %||% NA_real_),
      as.character(isTRUE(model_spec$is_primary))
    ),
    stringsAsFactors = FALSE
  )
}

read_config <- function(path) {
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

intervals_table <- function(cfg) {
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

interval_hit <- function(x, interval) {
  left <- if (isTRUE(interval$include_start)) x >= interval$start else x > interval$start
  right <- if (isTRUE(interval$include_end)) x <= interval$end else x < interval$end
  left & right
}

read_cell_metadata <- function(path) {
  required <- c(
    "cell_id", "sample_id", "cluster", "initial_ploidy",
    "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "pseudotime"
  )
  meta <- as.data.frame(readr::read_csv(path, show_col_types = FALSE), stringsAsFactors = FALSE)
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

unique_value <- function(x, label) {
  x <- unique(as.character(x[!is.na(x)]))
  x <- x[nzchar(x)]
  if (length(x) > 1L) stop("Non-unique sample-level value for ", label, ": ", paste(x, collapse = ", "), call. = FALSE)
  if (length(x) == 0L) NA_character_ else x[[1L]]
}

read_etp_compartment <- function(path, compartment) {
  if (!file.exists(path)) stop("Missing ETP metadata input: ", path, call. = FALSE)
  data <- as.data.frame(readr::read_csv(path, show_col_types = FALSE), stringsAsFactors = FALSE)
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

build_etp_assignments <- function(cellcycle_path, noncell_path, etp_specs) {
  union_cells <- rbind(
    read_etp_compartment(cellcycle_path, "CellCycle"),
    read_etp_compartment(noncell_path, "NonCellCycle")
  )
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
      initial_ploidy = unique_value(x$initial_ploidy, paste(sample_id, "initial_ploidy")),
      dose = unique_value(x$gemcitabine_dose, paste(sample_id, "dose")),
      dose_mg = unique(x$gemcitabine_dose_mg_per_kg[is.finite(x$gemcitabine_dose_mg_per_kg)])[1L],
      n_endpoint_ploidy_cells = nrow(x),
      n_cellcycle_endpoint_ploidy_cells = sum(x$etp_source_compartment == "CellCycle"),
      n_noncellcycle_endpoint_ploidy_cells = sum(x$etp_source_compartment == "NonCellCycle"),
      sample_mean_endpoint_ploidy = mean(x$cell_ploidy),
      sample_median_endpoint_ploidy = stats::median(x$cell_ploidy),
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

attach_etp_assignments <- function(meta, assignments) {
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

audit_metadata <- function(meta) {
  duplicate_cells <- meta[duplicated(meta$cell_id) | duplicated(meta$cell_id, fromLast = TRUE), , drop = FALSE]
  per_cell <- stats::aggregate(
    cbind(sample_id = meta$sample_id, dose = meta$gemcitabine_dose, initial_ploidy = meta$initial_ploidy),
    by = list(cell_id = meta$cell_id),
    FUN = function(x) length(unique(as.character(x)))
  )
  inconsistent <- per_cell[per_cell$sample_id > 1L | per_cell$dose > 1L | per_cell$initial_ploidy > 1L, , drop = FALSE]
  list(duplicate_cells = duplicate_cells, inconsistent_cells = inconsistent)
}

check_pseudotime <- function(meta, tolerance = 1e-8) {
  bad <- meta[!is.finite(meta$pseudotime) | meta$pseudotime < -tolerance | meta$pseudotime > 1 + tolerance, , drop = FALSE]
  if (nrow(bad) > 0L) stop("Pseudotime contains non-finite or out-of-[0,1] values.", call. = FALSE)
  invisible(TRUE)
}

load_counts <- function(seurat_rds, assay, counts_layer) {
  obj <- readRDS(seurat_rds)
  counts <- get_assay_matrix(obj, assay = assay, slot_name = counts_layer)
  if (is.null(counts)) stop("Cannot read assay matrix: assay=", assay, ", layer/slot=", counts_layer, call. = FALSE)
  if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  counts
}

match_cells <- function(meta, counts, min_match_rate) {
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
  if (audit$match_rate[audit$source == "metadata"] < min_match_rate) {
    stop("Metadata-to-expression cell-ID match rate is below threshold: ", audit$match_rate[audit$source == "metadata"], call. = FALSE)
  }
  keep_meta <- meta$cell_id %in% expression_ids
  meta <- meta[keep_meta, , drop = FALSE]
  counts <- counts[, meta$cell_id, drop = FALSE]
  list(meta = meta, counts = counts, audit = audit)
}

limit_for_smoke <- function(meta, counts, max_cells = 0L, max_genes = 0L, seed = 1L) {
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

sample_metadata <- function(meta) {
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
        vapply(etp_threshold_specs(), `[[`, character(1L), "group_column")
      ),
      names(x)
    )
    for (column in optional) {
      if (is.numeric(x[[column]])) {
        values <- unique(x[[column]][is.finite(x[[column]])])
        base[[column]] <- if (length(values) == 0L) NA_real_ else values[[1L]]
      } else {
        base[[column]] <- unique_value(x[[column]], paste(sample_id, column))
      }
    }
    base
  })
  out <- do.call(rbind, rows)
  rownames(out) <- out$sample_id
  out
}

bin_metadata <- function(meta, n_bins) {
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

construct_pseudobulk <- function(meta, counts, n_bins, min_cells) {
  binned <- bin_metadata(meta, n_bins)
  meta <- binned$meta
  bin_info <- binned$bin_info
  sample_meta <- sample_metadata(meta)
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

region_coverage <- function(meta, sample_bin_meta, cfg) {
  sample_meta <- unique(meta[, c("sample_id", "gemcitabine_dose", "gemcitabine_dose_mg_per_kg", "initial_ploidy"), drop = FALSE])
  rows <- list()
  for (interval in cfg$interval_list) {
    cell_hit <- interval_hit(meta$pseudotime, interval)
    bin_hit <- interval_hit(sample_bin_meta$bin_midpoint, interval)
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

check_primary_coverage <- function(coverage) {
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

prepare_design <- function(meta, spline_df, model_spec, include_dose = TRUE) {
  meta <- as.data.frame(meta, stringsAsFactors = FALSE)
  meta$dose_mg_factor <- factor(meta$dose_mg)
  meta$initial_ploidy_factor <- factor(meta$initial_ploidy)
  for (spec in etp_threshold_specs()) {
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
  covariate_terms <- model_spec$covariate_terms %||% character()
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

new_design_rows <- function(model, pseudotime_grid) {
  new_basis <- predict(model$basis, newx = pseudotime_grid)
  colnames(new_basis) <- colnames(model$basis)
  new_data <- model$design_data[rep(1L, length(pseudotime_grid)), , drop = FALSE]
  rownames(new_data) <- NULL
  new_data$bin_midpoint <- pseudotime_grid
  for (column in colnames(new_basis)) new_data[[column]] <- new_basis[, column]
  design <- stats::model.matrix(model$formula, new_data)
  missing <- setdiff(colnames(model$design), colnames(design))
  for (column in missing) design <- cbind(design, setNames(data.frame(0), column))
  design[, colnames(model$design), drop = FALSE]
}

contrast_vector <- function(model, cfg, contrast_name, grid_size = 501L) {
  if (!is.null(model$design) && !is.null(model$design$basis)) model <- model$design
  grid <- seq(0, 1, length.out = grid_size)
  grid_design <- new_design_rows(model, grid)
  intervals <- cfg$interval_list
  mean_design <- function(mask) {
    if (!any(mask)) stop("Contrast grid has no points for requested interval.", call. = FALSE)
    colMeans(grid_design[mask, , drop = FALSE])
  }
  if (identical(contrast_name, "primary_adjacent_state")) {
    primary <- interval_hit(grid, intervals$primary_accumulated_state)
    left <- interval_hit(grid, intervals$left_neighbor)
    right <- interval_hit(grid, intervals$right_neighbor)
    mean_design(primary) - 0.5 * mean_design(left) - 0.5 * mean_design(right)
  } else {
    stop("Unsupported contrast in standalone script: ", contrast_name, call. = FALSE)
  }
}

fit_pseudotime_model <- function(pb_counts, sample_bin_meta, spline_df, min_primary_mice, model_spec) {
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
  design_bundle <- prepare_design(meta, spline_df, model_spec = model_spec, include_dose = TRUE)
  original_design_columns <- colnames(design_bundle$design)
  qr_rank <- qr(design_bundle$design)$rank
  dropped_design_columns <- character()
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
    include_dose = TRUE,
    model_spec = model_spec,
    original_design_columns = original_design_columns,
    retained_design_columns = colnames(design_bundle$design),
    dropped_design_columns = dropped_design_columns,
    design_rank = qr_rank,
    design_ncol_original = length(original_design_columns),
    design_n_observations = nrow(design_bundle$design)
  )
}

design_audit <- function(model_fit) {
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

contrast_table <- function(model_fit, contrast_vector, contrast_id) {
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

resolve_gene_symbols <- function(primary_table) {
  df <- primary_table
  missing_symbol <- is.na(df$gene_symbol) | df$gene_symbol == ""
  df$gene_symbol[missing_symbol] <- df$gene[missing_symbol]
  df$abs_t <- abs(df$t_statistic)
  df <- df[order(df$gene_symbol, -df$abs_t, df$p_value), , drop = FALSE]
  df$resolution_rank <- ave(df$abs_t, df$gene_symbol, FUN = function(x) rank(-x, ties.method = "first"))
  df$retained_for_gsea <- df$resolution_rank == 1L
  df[, c("gene", "gene_symbol", "contrast_id", "t_statistic", "p_value", "fdr", "resolution_rank", "retained_for_gsea"), drop = FALSE]
}

ranked_stats <- function(contrast_table, symbol_resolution = NULL) {
  df <- contrast_table
  missing_symbol <- is.na(df$gene_symbol) | df$gene_symbol == ""
  df$gene_symbol[missing_symbol] <- df$gene[missing_symbol]
  df <- df[!is.na(df$gene_symbol) & df$gene_symbol != "" & is.finite(df$t_statistic), , drop = FALSE]
  df$abs_t <- abs(df$t_statistic)
  df <- df[order(df$gene_symbol, -df$abs_t, df$p_value), , drop = FALSE]
  df <- df[!duplicated(df$gene_symbol), , drop = FALSE]
  stats <- df$t_statistic
  names(stats) <- df$gene_symbol
  sort(stats, decreasing = TRUE)
}

fetch_gene_sets <- function(collections, species = "Homo sapiens") {
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
    list(collection = collection, label = collection_label(collection), sets = sets)
  }
  collections <- trimws(unlist(strsplit(collections, ",", fixed = TRUE)))
  collections <- collections[nzchar(collections)]
  out <- lapply(collections, fetch_one)
  names(out) <- vapply(out, `[[`, character(1L), "collection")
  out
}

run_fgsea_collection <- function(stats, collection_obj, ranking_id, min_size, max_size, nperm_simple, seed) {
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
  out$pathway_label <- format_pathway_label(out$pathway, collection_obj$collection)
  out$direction <- ifelse(out$NES >= 0, "positive", "negative")
  out <- out[order(out$padj, -abs(out$NES), out$pathway), , drop = FALSE]
  rownames(out) <- NULL
  out
}

run_all_gsea <- function(stats, gene_sets, ranking_id, min_size, max_size, nperm_simple, seed) {
  rows <- lapply(seq_along(gene_sets), function(i) {
    run_fgsea_collection(stats, gene_sets[[i]], ranking_id, min_size, max_size, nperm_simple, seed + i)
  })
  do.call(rbind, rows)
}

leading_edge_table <- function(gsea_df) {
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
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) data.frame() else out
}

fitted_grid <- function(model_fit, grid_size = 501L) {
  grid <- seq(0, 1, length.out = grid_size)
  design <- new_design_rows(model_fit$design, grid)
  coef <- model_fit$fit$coefficients[, colnames(design), drop = FALSE]
  fitted <- coef %*% t(design)
  colnames(fitted) <- sprintf("%.6f", grid)
  list(grid = grid, fitted = fitted)
}

pathway_activity <- function(fitted_grid, gsea_df, max_pathways = Inf) {
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
  out <- do.call(rbind, Filter(Negate(is.null), rows))
  if (is.null(out)) data.frame() else out
}

resolve_seurat_rds <- function(seurat_rds) {
  if (is.null(seurat_rds) || !nzchar(seurat_rds)) {
    stop("--seurat_rds is required when the default path is intentionally blank.", call. = FALSE)
  }
  if (!file.exists(seurat_rds)) {
    stop("--seurat_rds does not exist: ", seurat_rds, call. = FALSE)
  }
  normalizePath(seurat_rds, mustWork = TRUE)
}

run_support_workflow <- function(args, repo_root) {
  check_packages()
  set.seed(args$seed)

  output_root <- clean_dir(args$output_root, overwrite = args$overwrite)
  manifest_dir <- ensure_dir(file.path(output_root, "00_manifest"))
  binning_root <- ensure_dir(file.path(output_root, "binning"))
  binning_manifest <- ensure_dir(file.path(binning_root, "00_manifest"))
  binning_qc <- ensure_dir(file.path(binning_root, "01_qc"))
  binning_pseudobulk <- ensure_dir(file.path(binning_root, "02_pseudobulk"))
  model_root <- ensure_dir(file.path(binning_root, "ETP_reference_balanced_threshold_2_24"))
  model_manifest <- ensure_dir(file.path(model_root, "00_manifest"))
  model_qc <- ensure_dir(file.path(model_root, "01_qc"))
  model_gene_models <- ensure_dir(file.path(model_root, "03_gene_models"))
  model_gsea <- ensure_dir(file.path(model_root, "04_gsea"))

  cfg <- read_config(args$config)
  seurat_rds <- resolve_seurat_rds(args$seurat_rds)
  params <- data.frame(parameter = names(args), value = vapply(args, as.character, character(1L)), stringsAsFactors = FALSE)

  write_csv(intervals_table(cfg), file.path(manifest_dir, "frozen_interval_definition.csv"))
  write_csv(package_versions(), file.path(manifest_dir, "package_versions.csv"))
  input_checksums <- data.frame(
    input = c("cell_metadata", "noncell_metadata", "seurat_rds", "config"),
    path = c(args$cell_metadata, args$noncell_metadata, seurat_rds, args$config),
    sha256 = c(
      file_checksum(args$cell_metadata),
      file_checksum(args$noncell_metadata),
      file_checksum(seurat_rds),
      file_checksum(args$config)
    ),
    stringsAsFactors = FALSE
  )
  write_csv(input_checksums, file.path(manifest_dir, "input_checksums.csv"))
  write_csv(params, file.path(binning_manifest, "analysis_parameters.csv"))

  message("Reading cell metadata")
  meta <- read_cell_metadata(args$cell_metadata)
  etp_assignments <- build_etp_assignments(args$cell_metadata, args$noncell_metadata, etp_threshold_specs())
  meta <- attach_etp_assignments(meta, etp_assignments)
  metadata_audit <- audit_metadata(meta)
  if (nrow(metadata_audit$inconsistent_cells) > 0L) {
    stop("Cell metadata has inconsistent sample/dose/ploidy assignments.", call. = FALSE)
  }
  check_pseudotime(meta)

  message("Reading Seurat counts: ", seurat_rds)
  counts <- load_counts(seurat_rds, args$assay, args$counts_layer)
  matched <- match_cells(meta, counts, args$min_match_rate)
  meta <- matched$meta
  counts <- matched$counts
  limited <- limit_for_smoke(meta, counts, args$max_cells, args$max_genes, args$seed)
  meta <- limited$meta
  counts <- limited$counts

  message("Constructing sample-bin pseudobulk")
  pb <- construct_pseudobulk(meta, counts, args$n_pseudotime_bins, args$min_cells_per_sample_bin)
  write_csv(pb$metadata, file.path(binning_pseudobulk, "sample_bin_metadata.csv"))
  coverage <- region_coverage(pb$cell_metadata, pb$metadata, cfg)
  coverage_check <- check_primary_coverage(coverage)
  write_csv(coverage_check, file.path(binning_qc, "primary_coverage_check.csv"))
  if (!isTRUE(coverage_check$passes)) {
    stop("Binning primary interval coverage check failed. See binning/01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  model_spec <- model_spec_etp_reference_balanced()
  write_csv(model_parameters_table(model_spec), file.path(model_manifest, "model_parameters.csv"))

  message("Fitting model: ", model_spec$model_id)
  min_primary_mice <- coverage_check$n_contributing_mice[[1L]]
  model_fit <- fit_pseudotime_model(pb$counts, pb$metadata, args$spline_df, min_primary_mice, model_spec)
  write_csv(design_audit(model_fit), file.path(model_qc, "model_design_rank_audit.csv"))

  primary_contrast <- contrast_vector(model_fit, cfg, "primary_adjacent_state", args$grid_size)
  primary_genes <- contrast_table(model_fit, primary_contrast, "primary_adjacent_state")
  primary_genes$model_id <- model_spec$model_id
  symbol_resolution <- resolve_gene_symbols(primary_genes)
  symbol_resolution$model_id <- model_spec$model_id
  write_csv(primary_genes, file.path(model_gene_models, "gene_primary_adjacent_state_contrast.csv"))
  write_csv(symbol_resolution, file.path(model_gene_models, "gene_symbol_resolution.csv"))

  message("Running GSEA")
  gene_sets <- fetch_gene_sets(args$gene_set_collections)
  primary_stats <- ranked_stats(primary_genes, symbol_resolution)
  primary_gsea <- run_all_gsea(
    primary_stats,
    gene_sets,
    "primary_adjacent_state",
    args$gsea_min_size,
    args$gsea_max_size,
    args$gsea_nperm_simple,
    args$seed
  )
  primary_gsea$model_id <- model_spec$model_id
  write_csv(primary_gsea, file.path(model_gsea, "all_collections_primary_adjacent_state_gsea.csv"))

  leading_edge <- leading_edge_table(primary_gsea)
  if (nrow(leading_edge) > 0L) leading_edge$model_id <- model_spec$model_id
  write_csv(leading_edge, file.path(model_gsea, "all_collections_leading_edge_genes.csv"))

  activity <- pathway_activity(fitted_grid(model_fit, args$grid_size), primary_gsea)
  if (nrow(activity) > 0L) activity$model_id <- model_spec$model_id
  write_csv(activity, file.path(model_gsea, "pathway_activity_over_pseudotime.csv"))

  message("Completed Figure 7 support folder: ", output_root)
  invisible(output_root)
}

main <- function() {
  path <- script_path()
  script_dir <- if (!is.na(path)) dirname(path) else getwd()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  defaults <- list(
    cell_metadata = file.path(repo_root, "Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    noncell_metadata = file.path(repo_root, "Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
    seurat_rds = file.path(repo_root, "Data/in-vivo/scRNA_Seq_Data/seurat_obj_annotated/integrated_sct_cca_seurat_final_reclustered.rds"),
    config = file.path(repo_root, "Code/in-vivo/04i_pseudotime_state_pathways_config.yaml"),
    output_root = file.path(repo_root, "Data/in-vivo/pseudotime_state_pathways"),
    snapshot_root = file.path(repo_root, "Figs/pseudotime_state_pathways"),
    assay = "RNA",
    counts_layer = "counts",
    n_pseudotime_bins = 20L,
    min_cells_per_sample_bin = 5L,
    min_cells_per_sample_region = 10L,
    spline_df = 5L,
    gene_set_collections = "H,C2:CP:REACTOME,C5:GO:BP",
    seed = 1L,
    workers = 4L,
    parallel = FALSE,
    model_workers = 0L,
    within_model_workers = 1L,
    gsea_workers = 0L,
    overwrite = FALSE,
    min_match_rate = 0.99,
    gsea_min_size = 15L,
    gsea_max_size = 500L,
    gsea_nperm_simple = 10000L,
    grid_size = 501L,
    max_cells = 0L,
    max_genes = 0L,
    run_sensitivity = FALSE,
    run_leave_one_out = FALSE,
    covariate_mode = "etp_group",
    etp_threshold = "reference_balanced_threshold_2_24",
    include_full_etp_models = FALSE,
    make_snapshot = FALSE
  )

  raw_args <- commandArgs(trailingOnly = TRUE)
  if (any(raw_args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0L)
  }
  args <- parse_cli_args(raw_args, defaults)
  args$cell_metadata <- resolve_path(args$cell_metadata, repo_root, must_work = TRUE)
  args$noncell_metadata <- resolve_path(args$noncell_metadata, repo_root, must_work = TRUE)
  args$config <- resolve_path(args$config, repo_root, must_work = TRUE)
  args$output_root <- resolve_path(args$output_root, repo_root, must_work = FALSE)
  args$seurat_rds <- resolve_path(args$seurat_rds, repo_root, must_work = FALSE)

  run_support_workflow(args, repo_root)
}

if (identical(environment(), globalenv())) {
  main()
}
