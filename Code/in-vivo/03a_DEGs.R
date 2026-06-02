#!/usr/bin/env Rscript

resolve_in_vivo_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  candidate_files <- character(0)
  if (length(file_match) > 0) {
    candidate_files <- c(candidate_files, sub("--file=", "", file_match[1]))
  }

  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  candidate_files <- c(candidate_files, frame_files[!is.na(frame_files)])

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
    if (nzchar(active_path)) candidate_files <- c(candidate_files, active_path)
  }

  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  cwd_parts <- strsplit(cwd, .Platform$file.sep, fixed = TRUE)[[1]]
  parent_dirs <- vapply(
    seq_along(cwd_parts),
    function(i) {
      paste(c(cwd_parts[seq_len(length(cwd_parts) - i + 1)]), collapse = .Platform$file.sep)
    },
    character(1)
  )
  parent_dirs <- parent_dirs[nzchar(parent_dirs)]
  parent_dirs <- if (grepl("^/", cwd)) paste0("/", sub("^/+", "", parent_dirs)) else parent_dirs
  candidate_dirs <- unique(c(
    candidate_dirs,
    cwd,
    file.path(cwd, "Code", "in-vivo"),
    parent_dirs,
    file.path(parent_dirs, "Code", "in-vivo")
  ))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

required_packages <- c(
  "Seurat",
  "dplyr",
  "future",
  "readr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(future)
  library(readr)
  library(tibble)
})

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  KMP_DUPLICATE_LIB_OK = "TRUE",
  KMP_INIT_AT_FORK = "FALSE"
)

set.seed(1234)

input_rds <- file.path(
  results_root,
  "03_final_cluster",
  "03_objects",
  "integrated_sct_cca_seurat_final_reclustered.rds"
)
output_root <- file.path(results_root, "03a_DEGs")

assay_use <- "RNA"
cluster_col <- "clusters"
ploidy_col <- "Ploidy"
tn_col <- "TN"
dose_col <- "Dose"
dose_deg_col <- "Dose_DEG"
ids_col <- "IDs"

min_cells_per_group <- 3L
expected_dose_levels <- c("0mg/kg", "30mg/kg", "120mg/kg")
deg_parallel_workers <- suppressWarnings(as.integer(Sys.getenv("DEG_PARALLEL_WORKERS", unset = "10")))
if (is.na(deg_parallel_workers) || deg_parallel_workers < 1L) {
  deg_parallel_workers <- 10L
}
future_globals_maxsize_gb <- suppressWarnings(as.numeric(Sys.getenv("DEG_FUTURE_GLOBALS_MAXSIZE_GB", unset = "80")))
if (is.na(future_globals_maxsize_gb) || future_globals_maxsize_gb <= 0) {
  future_globals_maxsize_gb <- 80
}

safe_file_component <- function(x) {
  x <- sanitize_path_component(x)
  x <- gsub("[.]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "value")
}

metadata_levels <- function(meta, col_name, preferred_levels = NULL) {
  values <- as.character(meta[[col_name]])
  values <- values[!is.na(values) & values != ""]
  observed <- unique(values)

  if (!is.null(preferred_levels)) {
    return(c(
      intersect(preferred_levels, observed),
      setdiff(sort_maybe_numeric(observed), preferred_levels)
    ))
  }

  if (is.factor(meta[[col_name]])) {
    factor_levels <- as.character(levels(meta[[col_name]]))
    return(c(
      intersect(factor_levels, observed),
      setdiff(sort_maybe_numeric(observed), factor_levels)
    ))
  }

  sort_maybe_numeric(observed)
}

derive_ploidy_tn_if_needed <- function(obj) {
  if (!(ploidy_col %in% colnames(obj@meta.data))) {
    assert_required_column(obj@meta.data, ids_col)
    ids_values <- as.character(obj@meta.data[[ids_col]])
    ids_values[is.na(ids_values)] <- ""
    ploidy_values <- ifelse(
      grepl("2N", ids_values, fixed = TRUE),
      "2N",
      ifelse(grepl("4N", ids_values, fixed = TRUE), "4N", NA_character_)
    )
    obj@meta.data[[ploidy_col]] <- ploidy_values
  }

  if (!(tn_col %in% colnames(obj@meta.data))) {
    assert_required_column(obj@meta.data, ids_col)
    ids_values <- as.character(obj@meta.data[[ids_col]])
    ids_values[is.na(ids_values)] <- ""
    tn_values <- ifelse(grepl("Cell-Culture", ids_values, fixed = TRUE), "CellLine", "Tumor")
    obj@meta.data[[tn_col]] <- factor(tn_values, levels = c("CellLine", "Tumor"))
  }

  ploidy_values <- as.character(obj@meta.data[[ploidy_col]])
  invalid_ploidy <- is.na(ploidy_values) | !(ploidy_values %in% c("2N", "4N"))
  if (any(invalid_ploidy)) {
    invalid_values <- sort(unique(ploidy_values[invalid_ploidy]))
    invalid_values[is.na(invalid_values)] <- "NA"
    stop(
      "Ploidy must contain only 2N and 4N. Invalid values: ",
      paste(utils::head(invalid_values, 20), collapse = ", "),
      call. = FALSE
    )
  }
  obj@meta.data[[ploidy_col]] <- factor(ploidy_values, levels = c("2N", "4N"))

  obj
}

format_deg_table <- function(de, comparison_meta) {
  de <- as.data.frame(de, stringsAsFactors = FALSE)
  de$gene <- rownames(de)
  de$gene_symbol <- clean_gene_symbols(de$gene)

  if (nrow(de) > 0) {
    lfc_col <- tryCatch(resolve_lfc_col(de), error = function(e) NA_character_)
    if ("p_val_adj" %in% colnames(de) && !is.na(lfc_col)) {
      de <- de[order(de$p_val_adj, -abs(de[[lfc_col]]), de$gene), , drop = FALSE]
    } else if ("p_val" %in% colnames(de)) {
      de <- de[order(de$p_val, de$gene), , drop = FALSE]
    }
  }

  rownames(de) <- NULL
  meta_df <- as.data.frame(comparison_meta, stringsAsFactors = FALSE)
  meta_df <- meta_df[rep(1, max(nrow(de), 1L)), , drop = FALSE]
  if (nrow(de) == 0) {
    de <- data.frame(gene = character(0), gene_symbol = character(0), stringsAsFactors = FALSE)
    meta_df <- meta_df[0, , drop = FALSE]
  }
  out <- cbind(meta_df, de)

  front_cols <- c(
    names(comparison_meta),
    "gene",
    "gene_symbol"
  )
  out[, c(front_cols, setdiff(colnames(out), front_cols)), drop = FALSE]
}

make_summary_row <- function(
  scope,
  comparison,
  group_col,
  ident_1,
  ident_2,
  n_ident_1,
  n_ident_2,
  status,
  n_markers = NA_integer_,
  output_file = NA_character_,
  message = NA_character_,
  cluster = NA_character_,
  subset_group = NA_character_,
  dose = NA_character_
) {
  data.frame(
    scope = scope,
    cluster = cluster,
    subset_group = subset_group,
    dose = dose,
    comparison = comparison,
    group_col = group_col,
    ident_1 = paste(as.character(ident_1), collapse = ","),
    ident_2 = if (is.null(ident_2)) "rest" else paste(as.character(ident_2), collapse = ","),
    n_ident_1 = as.integer(n_ident_1),
    n_ident_2 = as.integer(n_ident_2),
    status = status,
    n_markers = as.integer(n_markers),
    output_file = output_file,
    message = message,
    stringsAsFactors = FALSE
  )
}

run_findmarkers_default <- function(
  obj,
  group_col,
  ident_1,
  ident_2 = NULL,
  comparison,
  output_file,
  assay,
  scope,
  cluster = NA_character_,
  subset_group = NA_character_,
  dose = NA_character_,
  min_cells = 3L
) {
  if (ncol(obj) == 0) {
    return(make_summary_row(
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose,
      comparison = comparison,
      group_col = group_col,
      ident_1 = ident_1,
      ident_2 = ident_2,
      n_ident_1 = 0L,
      n_ident_2 = 0L,
      status = "skipped",
      message = "No cells in subset."
    ))
  }

  assert_required_column(obj@meta.data, group_col)
  group_values <- as.character(obj@meta.data[[group_col]])
  valid_flag <- !is.na(group_values) & group_values != ""
  valid_cells <- rownames(obj@meta.data)[valid_flag]

  if (length(valid_cells) == 0) {
    return(make_summary_row(
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose,
      comparison = comparison,
      group_col = group_col,
      ident_1 = ident_1,
      ident_2 = ident_2,
      n_ident_1 = 0L,
      n_ident_2 = 0L,
      status = "skipped",
      message = paste0("No non-missing values in ", group_col, ".")
    ))
  }

  obj_use <- subset(obj, cells = valid_cells)
  group_values <- as.character(obj_use@meta.data[[group_col]])
  ident_1 <- as.character(ident_1)
  ident_2_chr <- if (is.null(ident_2)) NULL else as.character(ident_2)

  cells_1 <- rownames(obj_use@meta.data)[group_values %in% ident_1]
  cells_2 <- if (is.null(ident_2_chr)) {
    rownames(obj_use@meta.data)[!(group_values %in% ident_1)]
  } else {
    rownames(obj_use@meta.data)[group_values %in% ident_2_chr]
  }

  n_ident_1 <- length(cells_1)
  n_ident_2 <- length(cells_2)

  if (n_ident_1 < min_cells || n_ident_2 < min_cells) {
    return(make_summary_row(
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose,
      comparison = comparison,
      group_col = group_col,
      ident_1 = ident_1,
      ident_2 = ident_2_chr,
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      status = "skipped",
      message = paste0("Fewer than ", min_cells, " cells in one or both groups.")
    ))
  }

  if (file.exists(output_file) && isTRUE(file.info(output_file)$size > 0)) {
    return(make_summary_row(
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose,
      comparison = comparison,
      group_col = group_col,
      ident_1 = ident_1,
      ident_2 = ident_2_chr,
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      status = "cached",
      output_file = output_file,
      message = "Existing DEG file found; skipped FindMarkers."
    ))
  }

  group_levels <- unique(c(
    ident_1,
    ident_2_chr,
    metadata_levels(obj_use@meta.data, group_col)
  ))
  obj_use@meta.data[[".deg_group"]] <- factor(group_values, levels = group_levels)
  Idents(obj_use) <- obj_use@meta.data[[".deg_group"]]

  .ensure_dir(dirname(output_file))
  comparison_meta <- list(
    scope = scope,
    cluster = cluster,
    subset_group = subset_group,
    dose = dose,
    comparison = comparison,
    group_col = group_col,
    ident_1 = paste(ident_1, collapse = ","),
    ident_2 = if (is.null(ident_2_chr)) "rest" else paste(ident_2_chr, collapse = ","),
    n_ident_1 = n_ident_1,
    n_ident_2 = n_ident_2
  )

  result <- tryCatch(
    {
      de <- Seurat::FindMarkers(
        object = obj_use,
        ident.1 = ident_1,
        ident.2 = ident_2_chr,
        assay = assay,
        verbose = FALSE
      )
      de_out <- format_deg_table(de, comparison_meta)
      write_table_csv(de_out, output_file)
      make_summary_row(
        scope = scope,
        cluster = cluster,
        subset_group = subset_group,
        dose = dose,
        comparison = comparison,
        group_col = group_col,
        ident_1 = ident_1,
        ident_2 = ident_2_chr,
        n_ident_1 = n_ident_1,
        n_ident_2 = n_ident_2,
        status = "success",
        n_markers = nrow(de_out),
        output_file = output_file,
        message = NA_character_
      )
    },
    error = function(e) {
      error_file <- file.path(
        dirname(output_file),
        paste0(tools::file_path_sans_ext(basename(output_file)), "_FindMarkers_error.txt")
      )
      writeLines(conditionMessage(e), con = error_file)
      make_summary_row(
        scope = scope,
        cluster = cluster,
        subset_group = subset_group,
        dose = dose,
        comparison = comparison,
        group_col = group_col,
        ident_1 = ident_1,
        ident_2 = ident_2_chr,
        n_ident_1 = n_ident_1,
        n_ident_2 = n_ident_2,
        status = "error",
        output_file = error_file,
        message = conditionMessage(e)
      )
    }
  )

  result
}

subset_by_metadata <- function(obj, col_name, value) {
  assert_required_column(obj@meta.data, col_name)
  vals <- as.character(obj@meta.data[[col_name]])
  cells <- rownames(obj@meta.data)[!is.na(vals) & vals == value]
  if (length(cells) == 0) return(NULL)
  subset(obj, cells = cells)
}

run_dose_vs_rest <- function(
  obj,
  output_dir,
  assay,
  scope,
  cluster = NA_character_,
  subset_group = NA_character_,
  min_cells = 3L
) {
  summaries <- list()
  dose_values <- as.character(obj@meta.data[[dose_deg_col]])
  present_doses <- unique(dose_values[!is.na(dose_values) & dose_values != ""])

  for (dose_level in expected_dose_levels) {
    rest_levels <- setdiff(present_doses, dose_level)
    dose_stub <- safe_file_component(dose_level)
    summaries[[dose_level]] <- run_findmarkers_default(
      obj = obj,
      group_col = dose_deg_col,
      ident_1 = dose_level,
      ident_2 = rest_levels,
      comparison = paste0("Dose_", dose_level, "_vs_rest"),
      output_file = file.path(output_dir, paste0("Dose_", dose_stub, "_vs_rest_markers.csv")),
      assay = assay,
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose_level,
      min_cells = min_cells
    )
  }

  dplyr::bind_rows(summaries)
}

make_findmarkers_task <- function(
  cells = NULL,
  group_col,
  ident_1,
  ident_2 = NULL,
  comparison,
  output_file,
  assay,
  scope,
  cluster = NA_character_,
  subset_group = NA_character_,
  dose = NA_character_,
  min_cells = 3L
) {
  list(
    type = "findmarkers",
    cells = cells,
    group_col = group_col,
    ident_1 = ident_1,
    ident_2 = ident_2,
    comparison = comparison,
    output_file = output_file,
    assay = assay,
    scope = scope,
    cluster = cluster,
    subset_group = subset_group,
    dose = dose,
    min_cells = min_cells
  )
}

make_skip_task <- function(
  scope,
  comparison,
  group_col,
  ident_1,
  ident_2,
  n_ident_1 = 0L,
  n_ident_2 = 0L,
  message,
  cluster = NA_character_,
  subset_group = NA_character_,
  dose = NA_character_
) {
  list(
    type = "skip",
    scope = scope,
    comparison = comparison,
    group_col = group_col,
    ident_1 = ident_1,
    ident_2 = ident_2,
    n_ident_1 = n_ident_1,
    n_ident_2 = n_ident_2,
    message = message,
    cluster = cluster,
    subset_group = subset_group,
    dose = dose
  )
}

run_deg_task <- function(task) {
  if (identical(task$type, "skip")) {
    return(make_summary_row(
      scope = task$scope,
      cluster = task$cluster,
      subset_group = task$subset_group,
      dose = task$dose,
      comparison = task$comparison,
      group_col = task$group_col,
      ident_1 = task$ident_1,
      ident_2 = task$ident_2,
      n_ident_1 = task$n_ident_1,
      n_ident_2 = task$n_ident_2,
      status = "skipped",
      message = task$message
    ))
  }

  task_cells <- if (is.null(task$cells)) {
    colnames(obj)
  } else {
    as.character(task$cells)
  }
  task_cells <- task_cells[!is.na(task_cells) & task_cells != ""]
  task_cells <- intersect(task_cells, colnames(obj))

  if (length(task_cells) == 0) {
    return(make_summary_row(
      scope = task$scope,
      cluster = task$cluster,
      subset_group = task$subset_group,
      dose = task$dose,
      comparison = task$comparison,
      group_col = task$group_col,
      ident_1 = task$ident_1,
      ident_2 = task$ident_2,
      n_ident_1 = 0L,
      n_ident_2 = 0L,
      status = "skipped",
      message = "No valid cells in task after removing NA or unmatched cell names."
    ))
  }

  obj_task <- if (is.null(task$cells)) {
    obj
  } else {
    subset(obj, cells = task_cells)
  }

  run_findmarkers_default(
    obj = obj_task,
    group_col = task$group_col,
    ident_1 = task$ident_1,
    ident_2 = task$ident_2,
    comparison = task$comparison,
    output_file = task$output_file,
    assay = task$assay,
    scope = task$scope,
    cluster = task$cluster,
    subset_group = task$subset_group,
    dose = task$dose,
    min_cells = task$min_cells
  )
}

make_dose_vs_rest_tasks <- function(
  cells,
  output_dir,
  assay,
  scope,
  cluster = NA_character_,
  subset_group = NA_character_,
  min_cells = 3L
) {
  cells <- as.character(cells)
  cells <- cells[!is.na(cells) & cells != ""]
  cells <- intersect(cells, rownames(obj@meta.data))
  dose_values <- as.character(obj@meta.data[cells, dose_deg_col, drop = TRUE])
  present_doses <- unique(dose_values[!is.na(dose_values) & dose_values != ""])
  tasks <- list()

  for (dose_level in expected_dose_levels) {
    rest_levels <- setdiff(present_doses, dose_level)
    dose_stub <- safe_file_component(dose_level)
    tasks[[length(tasks) + 1L]] <- make_findmarkers_task(
      cells = cells,
      group_col = dose_deg_col,
      ident_1 = dose_level,
      ident_2 = rest_levels,
      comparison = paste0("Dose_", dose_level, "_vs_rest"),
      output_file = file.path(output_dir, paste0("Dose_", dose_stub, "_vs_rest_markers.csv")),
      assay = assay,
      scope = scope,
      cluster = cluster,
      subset_group = subset_group,
      dose = dose_level,
      min_cells = min_cells
    )
  }

  tasks
}

run_deg_tasks_parallel <- function(tasks, workers, future_max_size_gb) {
  if (length(tasks) == 0) return(list())

  workers <- max(1L, min(as.integer(workers), length(tasks)))
  if (workers <= 1L) {
    configure_future_for_seurat(
      max_size_gb = future_max_size_gb,
      strategy = "sequential"
    )
    return(lapply(tasks, run_deg_task))
  }

  configure_future_for_seurat(
    max_size_gb = future_max_size_gb,
    strategy = "multisession",
    workers = workers
  )
  on.exit(future::plan(future::sequential), add = TRUE)

  task_chunks <- split(tasks, rep(seq_len(workers), length.out = length(tasks)))
  message("Running ", length(tasks), " DEG tasks with ", workers, " multisession workers.")
  futures <- lapply(task_chunks, function(chunk) {
    future::future(
      {
        lapply(chunk, run_deg_task)
      },
      packages = c("Seurat", "dplyr", "readr", "tibble"),
      seed = TRUE
    )
  })
  unlist(lapply(futures, future::value), recursive = FALSE)
}

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_clusters <- .ensure_dir(file.path(output_root, "01_clusters_vs_rest"))
out_ploidy <- .ensure_dir(file.path(output_root, "02_Ploidy_4N_vs_2N"))
out_tn <- .ensure_dir(file.path(output_root, "03_TN_Tumor_vs_CellLine"))
out_dose <- .ensure_dir(file.path(output_root, "04_Dose_vs_rest"))
out_within_clusters <- .ensure_dir(file.path(output_root, "05_within_clusters"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

if (!(assay_use %in% names(obj@assays))) {
  stop("Assay is missing in Seurat object: ", assay_use, call. = FALSE)
}

assert_required_column(obj@meta.data, cluster_col)
assert_required_column(obj@meta.data, dose_col)

obj <- derive_ploidy_tn_if_needed(obj)
assert_required_column(obj@meta.data, ploidy_col)
assert_required_column(obj@meta.data, tn_col)

message("Preparing RNA assay for Seurat::FindMarkers.")
obj <- prepare_obj_for_markers(obj, assay = assay_use)
DefaultAssay(obj) <- assay_use

cluster_levels <- metadata_levels(obj@meta.data, cluster_col)
obj@meta.data[[cluster_col]] <- factor(
  as.character(obj@meta.data[[cluster_col]]),
  levels = cluster_levels
)
obj@meta.data[[ploidy_col]] <- factor(
  as.character(obj@meta.data[[ploidy_col]]),
  levels = c("2N", "4N")
)
obj@meta.data[[tn_col]] <- factor(
  as.character(obj@meta.data[[tn_col]]),
  levels = c("CellLine", "Tumor")
)

dose_raw <- obj@meta.data[[dose_col]]
dose_standardized <- standardize_in_vivo_dose(dose_raw)
missing_dose <- is.na(dose_standardized)
if (any(missing_dose)) {
  missing_dose_df <- data.frame(
    cell = rownames(obj@meta.data)[missing_dose],
    Dose_raw = as.character(dose_raw)[missing_dose],
    stringsAsFactors = FALSE
  )
  write_table_csv(missing_dose_df, file.path(out_summary, "cells_with_missing_Dose.csv"))
  message("Dose has ", sum(missing_dose), " NA/blank cell(s); Dose DEG analyses will exclude them.")
}
obj@meta.data[[dose_deg_col]] <- factor(dose_standardized, levels = expected_dose_levels)

metadata_count_list <- list(
  clusters = data.frame(
    metadata = cluster_col,
    value = names(table(obj@meta.data[[cluster_col]])),
    n_cells = as.integer(table(obj@meta.data[[cluster_col]])),
    stringsAsFactors = FALSE
  ),
  Ploidy = data.frame(
    metadata = ploidy_col,
    value = names(table(obj@meta.data[[ploidy_col]])),
    n_cells = as.integer(table(obj@meta.data[[ploidy_col]])),
    stringsAsFactors = FALSE
  ),
  TN = data.frame(
    metadata = tn_col,
    value = names(table(obj@meta.data[[tn_col]])),
    n_cells = as.integer(table(obj@meta.data[[tn_col]])),
    stringsAsFactors = FALSE
  ),
  Dose = data.frame(
    metadata = dose_deg_col,
    value = names(table(obj@meta.data[[dose_deg_col]], useNA = "no")),
    n_cells = as.integer(table(obj@meta.data[[dose_deg_col]], useNA = "no")),
    stringsAsFactors = FALSE
  )
)
write_table_csv(dplyr::bind_rows(metadata_count_list), file.path(out_summary, "metadata_cell_counts.csv"))

deg_tasks <- list()
add_task <- function(task) {
  task$task_id <- length(deg_tasks) + 1L
  deg_tasks[[task$task_id]] <<- task
  invisible(NULL)
}
add_tasks <- function(tasks) {
  for (task in tasks) add_task(task)
  invisible(NULL)
}

message("Preparing cluster vs rest DEG tasks.")
all_cells <- rownames(obj@meta.data)
for (cluster_id in cluster_levels) {
  cluster_stub <- safe_file_component(cluster_id)
  add_task(make_findmarkers_task(
    cells = all_cells,
    group_col = cluster_col,
    ident_1 = cluster_id,
    ident_2 = NULL,
    comparison = paste0("cluster_", cluster_id, "_vs_rest"),
    output_file = file.path(out_clusters, paste0("cluster_", cluster_stub, "_vs_rest_markers.csv")),
    assay = assay_use,
    scope = "clusters_vs_rest",
    cluster = cluster_id,
    min_cells = min_cells_per_group
  ))
}

message("Preparing global Ploidy DEG task: 4N vs 2N.")
add_task(make_findmarkers_task(
  cells = all_cells,
  group_col = ploidy_col,
  ident_1 = "4N",
  ident_2 = "2N",
  comparison = "Ploidy_4N_vs_2N",
  output_file = file.path(out_ploidy, "Ploidy_4N_vs_2N_markers.csv"),
  assay = assay_use,
  scope = "global_Ploidy_4N_vs_2N",
  min_cells = min_cells_per_group
))

message("Preparing global TN DEG task: Tumor vs CellLine.")
add_task(make_findmarkers_task(
  cells = all_cells,
  group_col = tn_col,
  ident_1 = "Tumor",
  ident_2 = "CellLine",
  comparison = "TN_Tumor_vs_CellLine",
  output_file = file.path(out_tn, "TN_Tumor_vs_CellLine_markers.csv"),
  assay = assay_use,
  scope = "global_TN_Tumor_vs_CellLine",
  min_cells = min_cells_per_group
))

message("Preparing global Dose vs rest DEG tasks.")
add_tasks(make_dose_vs_rest_tasks(
  cells = all_cells,
  output_dir = out_dose,
  assay = assay_use,
  scope = "global_Dose_vs_rest",
  min_cells = min_cells_per_group
))

message("Preparing within-cluster DEG tasks.")
for (cluster_id in cluster_levels) {
  cluster_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[cluster_col]]) == cluster_id]
  if (length(cluster_cells) == 0) next

  cluster_stub <- safe_file_component(cluster_id)
  cluster_dir <- .ensure_dir(file.path(out_within_clusters, paste0("cluster_", cluster_stub)))
  cluster_ploidy_dir <- .ensure_dir(file.path(cluster_dir, "01_Ploidy_4N_vs_2N"))
  cluster_tn_dir <- .ensure_dir(file.path(cluster_dir, "02_TN_Tumor_vs_CellLine"))
  cluster_dose_dir <- .ensure_dir(file.path(cluster_dir, "03_Dose_vs_rest"))
  cluster_tumor_ploidy_dir <- .ensure_dir(file.path(cluster_dir, "04_Tumor_only_Ploidy_4N_vs_2N"))
  cluster_cellline_ploidy_dir <- .ensure_dir(file.path(cluster_dir, "05_CellLine_only_Ploidy_4N_vs_2N"))
  cluster_dose_ploidy_dir <- .ensure_dir(file.path(cluster_dir, "06_Ploidy_within_Dose_4N_vs_2N"))

  add_task(make_findmarkers_task(
    cells = cluster_cells,
    group_col = ploidy_col,
    ident_1 = "4N",
    ident_2 = "2N",
    comparison = paste0("cluster_", cluster_id, "_Ploidy_4N_vs_2N"),
    output_file = file.path(cluster_ploidy_dir, "Ploidy_4N_vs_2N_markers.csv"),
    assay = assay_use,
    scope = "within_cluster_Ploidy_4N_vs_2N",
    cluster = cluster_id,
    min_cells = min_cells_per_group
  ))

  add_task(make_findmarkers_task(
    cells = cluster_cells,
    group_col = tn_col,
    ident_1 = "Tumor",
    ident_2 = "CellLine",
    comparison = paste0("cluster_", cluster_id, "_TN_Tumor_vs_CellLine"),
    output_file = file.path(cluster_tn_dir, "TN_Tumor_vs_CellLine_markers.csv"),
    assay = assay_use,
    scope = "within_cluster_TN_Tumor_vs_CellLine",
    cluster = cluster_id,
    min_cells = min_cells_per_group
  ))

  add_tasks(make_dose_vs_rest_tasks(
    cells = cluster_cells,
    output_dir = cluster_dose_dir,
    assay = assay_use,
    scope = "within_cluster_Dose_vs_rest",
    cluster = cluster_id,
    min_cells = min_cells_per_group
  ))

  tumor_cells <- cluster_cells[as.character(obj@meta.data[cluster_cells, tn_col, drop = TRUE]) %in% "Tumor"]
  if (length(tumor_cells) == 0) {
    add_task(make_skip_task(
      scope = "within_cluster_Tumor_only_Ploidy_4N_vs_2N",
      cluster = cluster_id,
      subset_group = "Tumor",
      comparison = paste0("cluster_", cluster_id, "_Tumor_only_Ploidy_4N_vs_2N"),
      group_col = ploidy_col,
      ident_1 = "4N",
      ident_2 = "2N",
      message = "No Tumor cells in cluster."
    ))
  } else {
    add_task(make_findmarkers_task(
      cells = tumor_cells,
      group_col = ploidy_col,
      ident_1 = "4N",
      ident_2 = "2N",
      comparison = paste0("cluster_", cluster_id, "_Tumor_only_Ploidy_4N_vs_2N"),
      output_file = file.path(cluster_tumor_ploidy_dir, "Tumor_only_Ploidy_4N_vs_2N_markers.csv"),
      assay = assay_use,
      scope = "within_cluster_Tumor_only_Ploidy_4N_vs_2N",
      cluster = cluster_id,
      subset_group = "Tumor",
      min_cells = min_cells_per_group
    ))
  }

  cellline_cells <- cluster_cells[as.character(obj@meta.data[cluster_cells, tn_col, drop = TRUE]) %in% "CellLine"]
  if (length(cellline_cells) == 0) {
    add_task(make_skip_task(
      scope = "within_cluster_CellLine_only_Ploidy_4N_vs_2N",
      cluster = cluster_id,
      subset_group = "CellLine",
      comparison = paste0("cluster_", cluster_id, "_CellLine_only_Ploidy_4N_vs_2N"),
      group_col = ploidy_col,
      ident_1 = "4N",
      ident_2 = "2N",
      message = "No CellLine cells in cluster."
    ))
  } else {
    add_task(make_findmarkers_task(
      cells = cellline_cells,
      group_col = ploidy_col,
      ident_1 = "4N",
      ident_2 = "2N",
      comparison = paste0("cluster_", cluster_id, "_CellLine_only_Ploidy_4N_vs_2N"),
      output_file = file.path(cluster_cellline_ploidy_dir, "CellLine_only_Ploidy_4N_vs_2N_markers.csv"),
      assay = assay_use,
      scope = "within_cluster_CellLine_only_Ploidy_4N_vs_2N",
      cluster = cluster_id,
      subset_group = "CellLine",
      min_cells = min_cells_per_group
    ))
  }

  for (dose_level in expected_dose_levels) {
    dose_cells <- cluster_cells[as.character(obj@meta.data[cluster_cells, dose_deg_col, drop = TRUE]) %in% dose_level]
    dose_stub <- safe_file_component(dose_level)
    if (length(dose_cells) == 0) {
      add_task(make_skip_task(
        scope = "within_cluster_Dose_Ploidy_4N_vs_2N",
        cluster = cluster_id,
        dose = dose_level,
        comparison = paste0("cluster_", cluster_id, "_Dose_", dose_level, "_Ploidy_4N_vs_2N"),
        group_col = ploidy_col,
        ident_1 = "4N",
        ident_2 = "2N",
        message = "No cells for this Dose in cluster."
      ))
    } else {
      add_task(make_findmarkers_task(
        cells = dose_cells,
        group_col = ploidy_col,
        ident_1 = "4N",
        ident_2 = "2N",
        comparison = paste0("cluster_", cluster_id, "_Dose_", dose_level, "_Ploidy_4N_vs_2N"),
        output_file = file.path(
          cluster_dose_ploidy_dir,
          paste0("Dose_", dose_stub, "_Ploidy_4N_vs_2N_markers.csv")
        ),
        assay = assay_use,
        scope = "within_cluster_Dose_Ploidy_4N_vs_2N",
        cluster = cluster_id,
        dose = dose_level,
        min_cells = min_cells_per_group
      ))
    }
  }
}

deg_parallel_workers_used <- if (length(deg_tasks) == 0) {
  0L
} else {
  max(1L, min(as.integer(deg_parallel_workers), length(deg_tasks)))
}
message(
  "Prepared ", length(deg_tasks), " DEG tasks. ",
  "Requested workers: ", deg_parallel_workers,
  "; effective workers: ", deg_parallel_workers_used, "."
)
comparison_summary <- run_deg_tasks_parallel(
  tasks = deg_tasks,
  workers = deg_parallel_workers_used,
  future_max_size_gb = future_globals_maxsize_gb
)

comparison_summary_df <- dplyr::bind_rows(comparison_summary)
write_table_csv(comparison_summary_df, file.path(out_summary, "DEG_comparison_summary.csv"))

run_summary <- c(
  "03a DEGs completed.",
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Assay: ", assay_use),
  "DEG method: Seurat::FindMarkers with default thresholds and default test settings.",
  paste0("DEG parallel backend: future::multisession when workers > 1, otherwise sequential."),
  paste0("DEG parallel workers requested: ", deg_parallel_workers),
  paste0("DEG parallel workers used: ", deg_parallel_workers_used),
  paste0("Future globals max size GB: ", future_globals_maxsize_gb),
  paste0("Minimum cells per group guard: ", min_cells_per_group),
  paste0("Cluster column: ", cluster_col),
  paste0("Ploidy comparison: 4N vs 2N"),
  paste0("TN comparison: Tumor vs CellLine"),
  paste0("Dose comparisons: ", paste(expected_dose_levels, collapse = ", "), " vs rest after filtering missing Dose values."),
  paste0("Total comparisons recorded: ", nrow(comparison_summary_df)),
  paste0("Successful comparisons: ", sum(comparison_summary_df$status == "success", na.rm = TRUE)),
  paste0("Cached comparisons skipped: ", sum(comparison_summary_df$status == "cached", na.rm = TRUE)),
  paste0("Skipped comparisons: ", sum(comparison_summary_df$status == "skipped", na.rm = TRUE)),
  paste0("Errored comparisons: ", sum(comparison_summary_df$status == "error", na.rm = TRUE))
)
writeLines(run_summary, con = file.path(out_summary, "run_summary.txt"))

message("03a_DEGs finished.")
message("Output root: ", output_root)
