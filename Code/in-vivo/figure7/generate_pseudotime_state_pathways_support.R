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
      "  Rscript Code/in-vivo/figure7/generate_pseudotime_state_pathways_support.R \\",
      "    --output_root Results/in-vivo/figure7/intermediates/state_pathways \\",
      "    --cell_metadata /path/to/CellCycleCells.csv \\",
      "    --noncell_metadata /path/to/NonCellCycleCells.csv \\",
      "    --seurat_rds /path/to/integrated_sct_cca_seurat_final_reclustered.rds",
      "",
      "Default input files:",
      "  Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "  Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
      "  Data/in-vivo/scRNA_Seq_Data/seurat_obj_annotated/integrated_sct_cca_seurat_final_reclustered.rds",
      "  Code/in-vivo/figure7/figure7_config.yaml",
      "",
      "This standalone script writes only the Figure 7 support files requested under:",
      "  <output_root>/00_manifest",
      "  <output_root>/binning",
      "  <output_root>/binning/initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4",
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
    "Seurat", "SeuratObject", "Matrix", "yaml", "digest", "edgeR",
    "limma", "splines", "fgsea", "BiocParallel", "msigdbr", "readr"
  )
}

recorded_packages <- function() {
  required_packages()
}

check_packages <- function(packages = required_packages()) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

package_versions <- function(packages = recorded_packages()) {
  data.frame(
    package = packages,
    version = vapply(packages, function(pkg) {
      if (!requireNamespace(pkg, quietly = TRUE) && !nzchar(system.file(package = pkg))) return(NA_character_)
      as.character(utils::packageVersion(pkg))
    }, character(1L)),
    stringsAsFactors = FALSE
  )
}

file_checksum <- function(path) {
  if (is.null(path) || is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  digest::digest(file = path, algo = "sha256")
}

portable_locator <- function(path, repo_root) {
  normalized <- normalizePath(path, mustWork = TRUE)
  root <- normalizePath(repo_root, mustWork = TRUE)
  prefix <- paste0(root, .Platform$file.sep)
  if (startsWith(normalized, prefix)) {
    return(substring(normalized, nchar(prefix) + 1L))
  }
  paste0("external:", basename(normalized))
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

format_pathway_label <- function(pathway, collection = NA_character_) {
  x <- as.character(pathway)
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^REACTOME_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- sub("^GO_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

clean_gene_symbols <- function(
  genes,
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  g <- as.character(genes)
  g <- trimws(g)
  if (any(startsWith(g, mouse_prefix))) {
    stop(
      "Mouse-prefixed features reached human-only symbol cleanup",
      call. = FALSE
    )
  }
  human <- startsWith(g, human_prefix)
  g[human] <- substring(g[human], nchar(human_prefix) + 1L)
  g <- sub("\\.[0-9]+$", "", g)
  g[g == ""] <- NA_character_
  g
}

get_assay_data_slot <- function(obj, assay = "RNA", slot_name = "counts") {
  assay_obj <- tryCatch({
    if (methods::is(obj, "Seurat") && assay %in% names(obj@assays)) {
      obj@assays[[assay]]
    } else {
      obj[[assay]]
    }
  }, error = function(e) NULL)
  if (is.null(assay_obj)) return(NULL)

  direct_slot <- tryCatch({
    if (methods::is(assay_obj, "Assay") && slot_name %in% slotNames(assay_obj)) {
      methods::slot(assay_obj, slot_name)
    } else {
      NULL
    }
  }, error = function(e) NULL)
  if (!is.null(direct_slot)) return(direct_slot)

  tryCatch(
    SeuratObject::LayerData(assay_obj, layer = slot_name),
    error = function(e1) {
      tryCatch(
        Seurat::GetAssayData(obj, assay = assay, slot = slot_name),
        error = function(e2) {
          tryCatch(
            Seurat::GetAssayData(obj, assay = assay, layer = slot_name),
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

join_assay_layers_for_counts <- function(
  obj,
  assay = "RNA",
  counts_layer = "counts"
) {
  if (!methods::is(obj, "Seurat") ||
      !assay %in% names(obj@assays) ||
      !"Layers" %in% getNamespaceExports("SeuratObject")) {
    return(obj)
  }
  count_layers <- SeuratObject::Layers(obj[[assay]])
  count_layers <- count_layers[
    count_layers == counts_layer |
      startsWith(count_layers, paste0(counts_layer, "."))
  ]
  if (!length(count_layers)) {
    stop(
      "Cannot find assay counts layer: assay=",
      assay,
      ", layer=",
      counts_layer,
      call. = FALSE
    )
  }
  if (length(count_layers) > 1L ||
      !identical(count_layers[[1L]], counts_layer)) {
    if (!"JoinLayers" %in% getNamespaceExports("SeuratObject")) {
      stop(
        "Multiple/split RNA counts layers require SeuratObject::JoinLayers",
        call. = FALSE
      )
    }
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = assay),
      error = function(error) {
        stop(
          "Cannot join split assay layers before panel-7F counts extraction: ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
  }
  joined_count_layers <- SeuratObject::Layers(obj[[assay]])
  joined_count_layers <- joined_count_layers[
    joined_count_layers == counts_layer |
      startsWith(joined_count_layers, paste0(counts_layer, "."))
  ]
  if (!identical(joined_count_layers, counts_layer)) {
    stop(
      "Panel 7F requires exactly one joined assay counts layer named ",
      counts_layer,
      "; observed ",
      paste(joined_count_layers, collapse = ","),
      call. = FALSE
    )
  }
  obj
}

model_spec_initial_ploidy <- function(
  model_id = "initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4"
) {
  list(
    model_id = model_id,
    model_label = "Injected-initial-ploidy-adjusted GRCh-only model",
    covariate_terms = "initial_ploidy_factor",
    covariate_mode = "initial_ploidy",
    initial_ploidy_levels = c("2N", "4N"),
    is_primary = FALSE
  )
}

model_parameters_table <- function(model_spec) {
  data.frame(
    parameter = c(
      "model_id", "model_label", "covariate_mode", "covariate_terms",
      "initial_ploidy_levels", "is_primary"
    ),
    value = c(
      model_spec$model_id %||% NA_character_,
      model_spec$model_label %||% NA_character_,
      model_spec$covariate_mode %||% NA_character_,
      paste(model_spec$covariate_terms %||% character(), collapse = ";"),
      paste(model_spec$initial_ploidy_levels %||% character(), collapse = ";"),
      as.character(isTRUE(model_spec$is_primary))
    ),
    stringsAsFactors = FALSE
  )
}

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  if (is.null(cfg$intervals)) stop("Config is missing intervals.", call. = FALSE)
  species <- cfg$feature_species
  expected_species <- c(
    policy_id = "grch_human_tumor_only_v2",
    human_prefix = "GRCh38-",
    mouse_prefix = "GRCm39-",
    unknown_feature_policy = "reject"
  )
  observed_species <- c(
    policy_id = as.character(species$policy_id),
    human_prefix = as.character(species$human_prefix),
    mouse_prefix = as.character(species$mouse_prefix),
    unknown_feature_policy =
      as.character(species$unknown_feature_policy)
  )
  if (!identical(observed_species, expected_species)) {
    stop(
      "Config must use the exact fail-closed GRCh38-/GRCm39- policy",
      call. = FALSE
    )
  }
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

load_counts <- function(
  seurat_rds,
  assay,
  counts_layer,
  human_prefix = "GRCh38-",
  mouse_prefix = "GRCm39-"
) {
  obj <- readRDS(seurat_rds)
  obj <- join_assay_layers_for_counts(
    obj,
    assay = assay,
    counts_layer = counts_layer
  )
  counts <- get_assay_matrix(obj, assay = assay, slot_name = counts_layer)
  if (is.null(counts)) stop("Cannot read assay matrix: assay=", assay, ", layer/slot=", counts_layer, call. = FALSE)
  if (!inherits(counts, "Matrix")) counts <- Matrix::Matrix(as.matrix(counts), sparse = TRUE)
  filtered <- figure7_filter_human_feature_matrix(
    counts,
    analysis = "Panel 7F RNA counts",
    human_prefix = human_prefix,
    mouse_prefix = mouse_prefix
  )
  list(counts = filtered$matrix, species_audit = filtered$audit)
}

match_cells <- function(meta, counts, min_match_rate) {
  metadata_ids <- unique(meta$cell_id)
  expression_ids <- colnames(counts)
  if (!isTRUE(all.equal(as.numeric(min_match_rate), 1, tolerance = 0))) {
    stop(
      "The pointwise-interval state model requires --min_match_rate=1 so its expression universe exactly matches the 2,881 cells used for localization.",
      call. = FALSE
    )
  }
  if (anyDuplicated(meta$cell_id) || anyDuplicated(expression_ids)) {
    stop(
      "Cell metadata and expression matrices require unique cell IDs.",
      call. = FALSE
    )
  }
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
  if (!all(meta_match)) {
    stop(
      "The exact 2,881-cell localization universe is not fully represented in the expression matrix: matched ",
      sum(meta_match), " of ", length(meta_match),
      call. = FALSE
    )
  }
  keep_meta <- meta$cell_id %in% expression_ids
  meta <- meta[keep_meta, , drop = FALSE]
  counts <- counts[, meta$cell_id, drop = FALSE]
  list(meta = meta, counts = counts, audit = audit)
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
  meta$initial_ploidy_factor <- factor(
    meta$initial_ploidy,
    levels = c("2N", "4N")
  )
  if (anyNA(meta$initial_ploidy_factor) ||
      !identical(levels(droplevels(meta$initial_ploidy_factor)), c("2N", "4N"))) {
    stop(
      "State-pathway model requires both injected initial-ploidy levels 2N and 4N",
      call. = FALSE
    )
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
    initial_ploidy_levels = paste(
      model_fit$model_spec$initial_ploidy_levels %||% character(),
      collapse = ";"
    ),
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
    if (!"db_version" %in% names(df)) {
      stop(
        "msigdbr output does not expose db_version; cannot verify the ",
        "MSigDB release",
        call. = FALSE
      )
    }
    database_release <- unique(as.character(df$db_version))
    database_release <- database_release[nzchar(database_release)]
    if (length(database_release) != 1L) {
      stop("Collection has an ambiguous MSigDB release: ", collection, call. = FALSE)
    }
    df$gene_symbol <- clean_gene_symbols(df$gene_symbol)
    sets <- split(df$gene_symbol, df$gs_name)
    sets <- lapply(sets, function(x) unique(x[!is.na(x) & x != ""]))
    list(
      collection = collection,
      label = collection_label(collection),
      database_release = database_release,
      sets = sets
    )
  }
  collections <- trimws(unlist(strsplit(collections, ",", fixed = TRUE)))
  collections <- collections[nzchar(collections)]
  out <- lapply(collections, fetch_one)
  names(out) <- vapply(out, `[[`, character(1L), "collection")
  out
}

gene_set_membership_table <- function(gene_sets) {
  rows <- unlist(lapply(gene_sets, function(collection) {
    lapply(names(collection$sets), function(pathway) {
      data.frame(
        collection_id = collection$collection,
        pathway_id = pathway,
        gene_symbol = sort(unique(as.character(collection$sets[[pathway]]))),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  out <- do.call(rbind, rows)
  out <- out[
    order(out$collection_id, out$pathway_id, out$gene_symbol),
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

run_fgsea_once <- function(
  pathways,
  stats,
  min_size,
  max_size,
  nperm_simple,
  seed
) {
  set.seed(seed)
  if (!requireNamespace("fgsea", quietly = TRUE) ||
      !"fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
    stop(
      "Figure 7 GSEA requires fgsea::fgseaMultilevel",
      call. = FALSE
    )
  }
  fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    nPermSimple = nperm_simple,
    nproc = 1L,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
    eps = 0
  )
}

fgsea_retry_budgets <- function(initial, maximum, multiplier) {
  values <- suppressWarnings(as.numeric(c(initial, maximum, multiplier)))
  if (length(values) != 3L || any(!is.finite(values)) ||
      any(values != floor(values)) || values[[1L]] < 1L ||
      values[[2L]] < values[[1L]] || values[[3L]] < 2L) {
    stop(
      "Adaptive GSEA requires positive integer initial/max nPermSimple ",
      "values with max >= initial and multiplier >= 2",
      call. = FALSE
    )
  }
  budgets <- as.integer(values[[1L]])
  maximum <- as.integer(values[[2L]])
  multiplier <- as.integer(values[[3L]])
  while (tail(budgets, 1L) < maximum) {
    next_budget <- min(
      as.double(tail(budgets, 1L)) * multiplier,
      maximum
    )
    if (!is.finite(next_budget) ||
        next_budget <= tail(budgets, 1L)) {
      stop("Adaptive GSEA retry schedule cannot advance", call. = FALSE)
    }
    budgets <- c(budgets, as.integer(next_budget))
  }
  budgets
}

fgsea_eligible_pathways <- function(pathways, stats, min_size, max_size) {
  stats_genes <- names(stats)
  if (is.null(stats_genes) || anyNA(stats_genes) ||
      any(!nzchar(stats_genes)) || anyDuplicated(stats_genes)) {
    stop("GSEA ranking requires unique nonempty gene names", call. = FALSE)
  }
  sizes <- vapply(
    pathways,
    function(genes) {
      length(intersect(unique(as.character(genes)), stats_genes))
    },
    integer(1L)
  )
  pathways[sizes >= min_size & sizes <= max_size]
}

fgsea_core_finite <- function(result) {
  core <- c("pval", "padj", "ES", "NES", "size")
  missing <- setdiff(core, names(result))
  if (length(missing)) {
    stop(
      "fgsea result is missing core field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  Reduce(
    `&`,
    lapply(
      result[core],
      function(column) is.finite(suppressWarnings(as.numeric(column)))
    )
  )
}

fgsea_retry_resolved <- function(result) {
  required <- c("pval", "ES", "NES")
  missing <- setdiff(required, names(result))
  if (length(missing)) {
    stop(
      "fgsea result is missing retry field(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  Reduce(
    `&`,
    lapply(
      result[required],
      function(column) is.finite(suppressWarnings(as.numeric(column)))
    )
  )
}

validate_fgsea_attempt <- function(result, requested_pathways, collection) {
  out <- as.data.frame(result, stringsAsFactors = FALSE)
  if (!"pathway" %in% names(out)) {
    stop("fgsea result is missing pathway identifiers", call. = FALSE)
  }
  out$pathway <- as.character(out$pathway)
  requested <- names(requested_pathways)
  if (anyNA(out$pathway) || any(!nzchar(out$pathway)) ||
      anyDuplicated(out$pathway) ||
      !identical(sort(out$pathway), sort(requested))) {
    stop(
      "fgsea returned incomplete or duplicate pathway keys for ",
      collection,
      call. = FALSE
    )
  }
  if (!"size" %in% names(out)) {
    stop("fgsea result is missing pathway sizes", call. = FALSE)
  }
  sizes <- suppressWarnings(as.numeric(out$size))
  if (any(!is.finite(sizes)) || any(sizes != floor(sizes)) ||
      any(sizes < 1L)) {
    stop(
      "fgsea returned invalid pathway sizes for ",
      collection,
      call. = FALSE
    )
  }
  out
}

run_fgsea_collection <- function(
  stats,
  collection_obj,
  ranking_id,
  min_size,
  max_size,
  nperm_simple,
  seed,
  nperm_simple_max = nperm_simple,
  nperm_simple_multiplier = 10L,
  fgsea_runner = run_fgsea_once
) {
  stats <- sort(stats[is.finite(stats)], decreasing = TRUE)
  if (length(stats) == 0L) return(data.frame())
  eligible <- fgsea_eligible_pathways(
    collection_obj$sets,
    stats,
    min_size,
    max_size
  )
  if (!length(eligible)) return(data.frame())
  budgets <- fgsea_retry_budgets(
    nperm_simple,
    nperm_simple_max,
    nperm_simple_multiplier
  )
  merged <- NULL
  unresolved <- names(eligible)
  for (retry_index in seq_along(budgets)) {
    if (!length(unresolved)) break
    requested <- eligible[unresolved]
    budget <- budgets[[retry_index]]
    message(
      "GSEA ",
      collection_obj$collection,
      ": nPermSimple=", budget,
      "; pathways=", length(requested),
      if (retry_index == 1L) "" else " (adaptive retry)"
    )
    attempt <- tryCatch(
      fgsea_runner(
        pathways = requested,
        stats = stats,
        min_size = min_size,
        max_size = max_size,
        nperm_simple = budget,
        seed = seed
      ),
      error = function(error) {
        stop(
          "fgsea failed for ",
          collection_obj$collection,
          " at nPermSimple=",
          budget,
          ": ",
          conditionMessage(error),
          call. = FALSE
        )
      }
    )
    attempt <- validate_fgsea_attempt(
      attempt,
      requested,
      collection_obj$collection
    )
    attempt$nPermSimple <- as.integer(budget)
    attempt$retry_round <- as.integer(retry_index - 1L)
    if (is.null(merged)) {
      merged <- attempt
    } else {
      replace_index <- match(attempt$pathway, merged$pathway)
      if (anyNA(replace_index)) {
        stop(
          "Adaptive GSEA retry returned an unknown pathway for ",
          collection_obj$collection,
          call. = FALSE
        )
      }
      merged[replace_index, names(attempt)] <- attempt
    }
    unresolved <- attempt$pathway[!fgsea_retry_resolved(attempt)]
  }
  if (length(unresolved)) {
    stop(
      "Adaptive GSEA exhausted nPermSimple=",
      tail(budgets, 1L),
      " with unresolved pathways in ",
      collection_obj$collection,
      ": ",
      paste(sort(unresolved), collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(merged) != length(eligible) ||
      anyDuplicated(merged$pathway) ||
      !identical(sort(merged$pathway), sort(names(eligible)))) {
    stop(
      "Adaptive GSEA merge is incomplete for ",
      collection_obj$collection,
      call. = FALSE
    )
  }
  merged$padj <- stats::p.adjust(
    as.numeric(merged$pval),
    method = "BH"
  )
  if (any(!fgsea_core_finite(merged))) {
    stop(
      "Adaptive GSEA merge has nonfinite final statistics for ",
      collection_obj$collection,
      call. = FALSE
    )
  }
  merged$direction <- ifelse(
    as.numeric(merged$NES) >= 0,
    "positive",
    "negative"
  )
  if ("leadingEdge" %in% names(merged)) {
    merged$leading_edge <- vapply(
      merged$leadingEdge,
      function(genes) paste(as.character(genes), collapse = ";"),
      character(1L)
    )
    merged$leadingEdge <- NULL
  } else {
    merged$leading_edge <- ""
  }
  merged$collection <- collection_obj$collection
  merged$collection_label <- collection_obj$label
  merged$ranking_id <- ranking_id
  merged$pathway_label <- format_pathway_label(
    merged$pathway,
    collection_obj$collection
  )
  merged <- merged[
    order(merged$padj, -abs(merged$NES), merged$pathway),
    ,
    drop = FALSE
  ]
  rownames(merged) <- NULL
  merged
}

run_all_gsea <- function(
  stats,
  gene_sets,
  ranking_id,
  min_size,
  max_size,
  nperm_simple,
  seed,
  nperm_simple_max = nperm_simple,
  nperm_simple_multiplier = 10L,
  fgsea_runner = run_fgsea_once
) {
  rows <- lapply(seq_along(gene_sets), function(i) {
    run_fgsea_collection(
      stats,
      gene_sets[[i]],
      ranking_id,
      min_size,
      max_size,
      nperm_simple,
      seed + i,
      nperm_simple_max,
      nperm_simple_multiplier,
      fgsea_runner
    )
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

pathway_activity <- function(
  fitted_grid,
  gsea_df,
  primary_interval,
  max_pathways = Inf
) {
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
      peaks_inside_primary = interval_hit(peak, primary_interval),
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

  final_output_root <- args$output_root
  if (!nzchar(final_output_root)) {
    stop("--output_root is required and must be run-scoped.", call. = FALSE)
  }
  if (dir.exists(final_output_root)) {
    existing <- list.files(
      final_output_root,
      all.files = TRUE,
      no.. = TRUE
    )
    if (length(existing) && !isTRUE(args$overwrite)) {
      stop(
        "Output directory is not empty. Use a new run-scoped --output_root: ",
        final_output_root,
        call. = FALSE
      )
    }
    if (length(existing)) unlink(final_output_root, recursive = TRUE, force = TRUE)
    if (dir.exists(final_output_root)) unlink(final_output_root, recursive = TRUE)
  }
  staging_root <- paste0(final_output_root, ".tmp-", Sys.getpid())
  if (dir.exists(staging_root)) unlink(staging_root, recursive = TRUE, force = TRUE)
  output_root <- clean_dir(staging_root, overwrite = FALSE)
  completed <- FALSE
  on.exit(
    if (!completed && dir.exists(staging_root)) {
      unlink(staging_root, recursive = TRUE, force = TRUE)
    },
    add = TRUE
  )
  manifest_dir <- ensure_dir(file.path(output_root, "00_manifest"))
  binning_root <- ensure_dir(file.path(output_root, "binning"))
  binning_manifest <- ensure_dir(file.path(binning_root, "00_manifest"))
  binning_qc <- ensure_dir(file.path(binning_root, "01_qc"))
  binning_pseudobulk <- ensure_dir(file.path(binning_root, "02_pseudobulk"))
  cfg <- read_config(args$config)
  cfg <- figure7_attach_density_localization_config(
    cfg,
    args$density_config
  )
  model_id <- as.character(cfg$state_pathways$model)
  if (!identical(
        model_id,
        "initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4"
      )) {
    stop(
      paste(
        "State-pathway config must request",
        "initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4"
      ),
      call. = FALSE
    )
  }
  model_root <- ensure_dir(file.path(binning_root, model_id))
  model_manifest <- ensure_dir(file.path(model_root, "00_manifest"))
  model_qc <- ensure_dir(file.path(model_root, "01_qc"))
  model_gene_models <- ensure_dir(file.path(model_root, "03_gene_models"))
  model_gsea <- ensure_dir(file.path(model_root, "04_gsea"))

  species <- cfg$feature_species
  configured_gsea_retry <- as.integer(c(
    cfg$state_pathways$gsea_nperm_simple,
    cfg$state_pathways$gsea_nperm_simple_max,
    cfg$state_pathways$gsea_nperm_simple_multiplier
  ))
  requested_gsea_retry <- as.integer(c(
    args$gsea_nperm_simple,
    args$gsea_nperm_simple_max,
    args$gsea_nperm_simple_multiplier
  ))
  if (!identical(requested_gsea_retry, configured_gsea_retry)) {
    stop(
      "Adaptive GSEA CLI values must exactly match state_pathways config: ",
      paste(configured_gsea_retry, collapse = ","),
      call. = FALSE
    )
  }
  species_policy_path <- file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src",
    "feature_species_policy.R"
  )
  support_script_path <- file.path(
    repo_root,
    "Code", "in-vivo", "figure7",
    "generate_pseudotime_state_pathways_support.R"
  )
  density_localization_code_path <- file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src",
    "tgi_statistics.R"
  )
  common_io_code_path <- file.path(
    repo_root,
    "Code", "in-vivo", "figure7", "src",
    "common_io.R"
  )
  if (!file.exists(species_policy_path)) {
    stop("Missing shared feature-species policy helper", call. = FALSE)
  }
  observed_msigdbr <- as.character(utils::packageVersion("msigdbr"))
  expected_msigdbr <- as.character(cfg$gene_sets$package_version)
  if (!identical(observed_msigdbr, expected_msigdbr)) {
    stop(
      "Generated human-only panel-7F rebuild requires msigdbr ",
      expected_msigdbr,
      "; observed ",
      observed_msigdbr,
      call. = FALSE
    )
  }
  seurat_rds <- resolve_seurat_rds(args$seurat_rds)
  params <- data.frame(
    parameter = c(
      names(args), "feature_species_policy_id",
      "human_feature_prefix", "mouse_feature_prefix",
      "unknown_feature_policy"
    ),
    value = c(
      vapply(args, as.character, character(1L)),
      as.character(species$policy_id),
      as.character(species$human_prefix),
      as.character(species$mouse_prefix),
      as.character(species$unknown_feature_policy)
    ),
    stringsAsFactors = FALSE
  )

  write_csv(package_versions(), file.path(manifest_dir, "package_versions.csv"))
  input_checksums <- data.frame(
    input = c(
      "cell_metadata", "noncell_metadata", "seurat_rds", "config",
      "support_script", "density_localization_config",
      "feature_species_policy_code",
      "density_localization_code", "common_io_code"
    ),
    locator = vapply(
      c(
        args$cell_metadata, args$noncell_metadata, seurat_rds,
        args$config, support_script_path, args$density_config,
        species_policy_path,
        density_localization_code_path,
        common_io_code_path
      ),
      portable_locator,
      character(1L),
      repo_root = repo_root
    ),
    sha256 = c(
      file_checksum(args$cell_metadata),
      file_checksum(args$noncell_metadata),
      file_checksum(seurat_rds),
      file_checksum(args$config),
      file_checksum(support_script_path),
      file_checksum(args$density_config),
      file_checksum(species_policy_path),
      file_checksum(density_localization_code_path),
      file_checksum(common_io_code_path)
    ),
    stringsAsFactors = FALSE
  )
  write_csv(input_checksums, file.path(manifest_dir, "input_checksums.csv"))
  write_csv(params, file.path(binning_manifest, "analysis_parameters.csv"))

  message("Reading cell metadata")
  meta <- read_cell_metadata(args$cell_metadata)
  metadata_audit <- audit_metadata(meta)
  if (nrow(metadata_audit$inconsistent_cells) > 0L) {
    stop("Cell metadata has inconsistent sample/dose/ploidy assignments.", call. = FALSE)
  }
  check_pseudotime(meta)

  message("Computing the state interval from exact density-localization support")
  density_localization <- figure7_density_localization(
    meta,
    sample_metadata(meta),
    cfg
  )
  cfg <- figure7_apply_density_supported_state_intervals(
    cfg,
    density_localization$intervals
  )
  state_interval_definition <- intervals_table(cfg)
  state_interval_definition$derivation_analysis_id <-
    density_localization$test$analysis_id[[1L]]
  state_interval_definition$derivation_support_type <-
    "positive_pointwise_two_sided"
  state_interval_definition$derivation_pointwise_alpha <-
    density_localization$test$pointwise_alpha[[1L]]
  state_interval_definition$derivation_n_permutations <-
    density_localization$test$n_permutations[[1L]]
  write_csv(
    state_interval_definition,
    file.path(manifest_dir, "state_interval_definition.csv")
  )
  write_csv(
    density_localization$grid,
    file.path(manifest_dir, "state_interval_localization_grid.csv")
  )
  write_csv(
    density_localization$intervals,
    file.path(manifest_dir, "state_interval_localization_support.csv")
  )
  write_csv(
    density_localization$test,
    file.path(manifest_dir, "state_interval_localization_test.csv")
  )

  message("Reading Seurat counts: ", seurat_rds)
  count_source <- load_counts(
    seurat_rds,
    args$assay,
    args$counts_layer,
    human_prefix = as.character(species$human_prefix),
    mouse_prefix = as.character(species$mouse_prefix)
  )
  counts <- count_source$counts
  write_csv(
    count_source$species_audit,
    file.path(manifest_dir, "feature_species_audit.csv")
  )
  matched <- match_cells(meta, counts, args$min_match_rate)
  write_csv(
    matched$audit,
    file.path(manifest_dir, "cell_expression_match_audit.csv")
  )
  meta <- matched$meta
  counts <- matched$counts
  if (nrow(meta) != 2881L || ncol(counts) != 2881L) {
    stop(
      "State-pathway expression modeling must retain exactly the same 2,881 cells used for density localization.",
      call. = FALSE
    )
  }
  message("Constructing sample-bin pseudobulk")
  pb <- construct_pseudobulk(meta, counts, args$n_pseudotime_bins, args$min_cells_per_sample_bin)
  write_csv(pb$metadata, file.path(binning_pseudobulk, "sample_bin_metadata.csv"))
  coverage <- region_coverage(pb$cell_metadata, pb$metadata, cfg)
  coverage_check <- check_primary_coverage(coverage)
  write_csv(coverage_check, file.path(binning_qc, "primary_coverage_check.csv"))
  if (!isTRUE(coverage_check$passes)) {
    stop("Binning primary interval coverage check failed. See binning/01_qc/primary_coverage_check.csv.", call. = FALSE)
  }

  model_spec <- model_spec_initial_ploidy(model_id)
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
  observed_releases <- unique(vapply(
    gene_sets,
    `[[`,
    character(1L),
    "database_release"
  ))
  expected_release <- as.character(cfg$gene_sets$database_release)
  if (!identical(observed_releases, expected_release)) {
    stop(
      "Generated human-only panel-7F rebuild requires MSigDB ",
      expected_release,
      "; observed ",
      paste(observed_releases, collapse = ","),
      call. = FALSE
    )
  }
  gene_set_membership <- gene_set_membership_table(gene_sets)
  write_csv(
    gene_set_membership,
    file.path(manifest_dir, "gene_set_membership.csv")
  )
  write_csv(
    data.frame(
      key = c(
        "provider", "msigdbr_package_version", "database_release",
        "species", "collections", "membership_sha256"
      ),
      value = c(
        "msigdbr",
        observed_msigdbr,
        observed_releases,
        as.character(cfg$gene_sets$species),
        paste(names(gene_sets), collapse = ","),
        file_checksum(file.path(manifest_dir, "gene_set_membership.csv"))
      ),
      stringsAsFactors = FALSE
    ),
    file.path(manifest_dir, "gene_set_contract.csv")
  )
  primary_stats <- ranked_stats(primary_genes, symbol_resolution)
  primary_gsea <- run_all_gsea(
    primary_stats,
    gene_sets,
    "primary_adjacent_state",
    args$gsea_min_size,
    args$gsea_max_size,
    args$gsea_nperm_simple,
    args$seed,
    args$gsea_nperm_simple_max,
    args$gsea_nperm_simple_multiplier
  )
  primary_gsea$model_id <- model_spec$model_id
  write_csv(primary_gsea, file.path(model_gsea, "all_collections_primary_adjacent_state_gsea.csv"))

  leading_edge <- leading_edge_table(primary_gsea)
  if (nrow(leading_edge) > 0L) leading_edge$model_id <- model_spec$model_id
  write_csv(leading_edge, file.path(model_gsea, "all_collections_leading_edge_genes.csv"))

  activity <- pathway_activity(
    fitted_grid(model_fit, args$grid_size),
    primary_gsea,
    cfg$interval_list$primary_accumulated_state
  )
  if (nrow(activity) > 0L) activity$model_id <- model_spec$model_id
  write_csv(activity, file.path(model_gsea, "pathway_activity_over_pseudotime.csv"))

  dir.create(dirname(final_output_root), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(staging_root, final_output_root)) {
    stop(
      "Could not atomically finalize Figure 7 support folder: ",
      final_output_root,
      call. = FALSE
    )
  }
  completed <- TRUE
  message("Completed Figure 7 support folder: ", final_output_root)
  invisible(final_output_root)
}

main <- function() {
  path <- script_path()
  script_dir <- if (!is.na(path)) dirname(path) else getwd()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
  sys.source(
    file.path(script_dir, "src", "common_io.R"),
    envir = .GlobalEnv
  )
  sys.source(
    file.path(script_dir, "src", "feature_species_policy.R"),
    envir = .GlobalEnv
  )
  sys.source(
    file.path(script_dir, "src", "tgi_statistics.R"),
    envir = .GlobalEnv
  )
  defaults <- list(
    cell_metadata = "",
    noncell_metadata = "",
    seurat_rds = "",
    config = file.path(repo_root, "Code/in-vivo/figure7/figure7_config.yaml"),
    density_config = file.path(
      repo_root,
      "Code/in-vivo/figure7/density_localization_config.yaml"
    ),
    output_root = "",
    assay = "RNA",
    counts_layer = "counts",
    n_pseudotime_bins = 20L,
    min_cells_per_sample_bin = 5L,
    spline_df = 5L,
    gene_set_collections = "H,C2:CP:REACTOME,C5:GO:BP",
    seed = 1L,
    overwrite = FALSE,
    min_match_rate = 1,
    gsea_min_size = 15L,
    gsea_max_size = 500L,
    gsea_nperm_simple = 10000L,
    gsea_nperm_simple_max = 1000000L,
    gsea_nperm_simple_multiplier = 10L,
    grid_size = 501L
  )

  raw_args <- commandArgs(trailingOnly = TRUE)
  if (any(raw_args %in% c("--help", "-h"))) {
    usage()
    quit(save = "no", status = 0L)
  }
  args <- parse_cli_args(raw_args, defaults)
  for (required in c(
    "cell_metadata", "noncell_metadata", "seurat_rds", "output_root"
  )) {
    if (!nzchar(args[[required]])) {
      stop("--", required, " is required.", call. = FALSE)
    }
  }
  args$cell_metadata <- resolve_path(args$cell_metadata, repo_root, must_work = TRUE)
  args$noncell_metadata <- resolve_path(args$noncell_metadata, repo_root, must_work = TRUE)
  args$config <- resolve_path(args$config, repo_root, must_work = TRUE)
  args$density_config <- resolve_path(
    args$density_config,
    repo_root,
    must_work = TRUE
  )
  args$output_root <- resolve_path(args$output_root, repo_root, must_work = FALSE)
  args$seurat_rds <- resolve_path(args$seurat_rds, repo_root, must_work = FALSE)

  run_support_workflow(args, repo_root)
}

if (identical(environment(), globalenv())) {
  main()
}
