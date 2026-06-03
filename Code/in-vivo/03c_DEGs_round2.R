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

cfg_value <- function(config, name, env_name = NULL, default = NULL) {
  if (!is.null(env_name)) {
    env_value <- trimws(Sys.getenv(env_name, unset = ""))
    if (nzchar(env_value)) return(env_value)
  }
  value <- config[[name]]
  if (!is.null(value) && length(value) > 0 && nzchar(trimws(as.character(value[1])))) {
    return(as.character(value[1]))
  }
  default
}

resolve_first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) normalizePath(hit[1], mustWork = TRUE) else paths[1]
}

input_rds <- cfg_value(
  config,
  "Trajectory_input_rds",
  "TRAJECTORY_SEURAT_RDS",
  default = resolve_first_existing(c(
    file.path(results_root, "03c_cluster_annotation", "04_objects", "integrated_sct_cca_seurat_cluster_final_annotation.rds"),
    file.path(results_root, "03b_manual_cluster_merge", "objects", "integrated_sct_cca_seurat_final_manual_merge.rds"),
    file.path(results_root, "03_final_cluster", "03_objects", "integrated_sct_cca_seurat_final_reclustered.rds")
  ))
)
output_root <- file.path(results_root, "03c_DEGs_round2")

assay_use <- "RNA"
cluster_col <- "cluster_final"
cluster_fallback_col <- "clusters"
ploidy_col <- "Ploidy"
tn_col <- "TN"
id_col_candidates <- unique(c("IDs", cfg_value(config, "Trajectory_id_col", "TRAJECTORY_ID_COL", default = "ID"), "ID"))
sample_folder_col_candidates <- unique(c(
  cfg_value(config, "Cellranger_sample_folder_col", "CELLRANGER_SAMPLE_FOLDER_COL", default = "sample_folder"),
  "sample_folder",
  "Sequencing.IDs",
  "sample",
  "orig.ident",
  "IDs"
))

min_cells_per_group <- 3L
deg_parallel_workers_setting <- cfg_value(config, "DEG_parallel_workers", "DEG_PARALLEL_WORKERS", default = "auto")
deg_worker_memory_factor <- suppressWarnings(as.numeric(cfg_value(config, "DEG_worker_memory_factor", "DEG_WORKER_MEMORY_FACTOR", default = "2.5")))
if (is.na(deg_worker_memory_factor) || !is.finite(deg_worker_memory_factor) || deg_worker_memory_factor <= 0) {
  deg_worker_memory_factor <- 2.5
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

parse_worker_setting <- function(value) {
  value <- trimws(as.character(value)[1])
  if (!nzchar(value) || tolower(value) %in% c("auto", "system")) return(NA_integer_)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) {
    stop("DEG_parallel_workers/DEG_PARALLEL_WORKERS must be a positive integer or 'auto'.", call. = FALSE)
  }
  out
}

detect_available_cores <- function() {
  cores <- suppressWarnings(tryCatch(future::availableCores(), error = function(e) NA_integer_))
  cores <- as.integer(cores[1])
  if (is.na(cores) || cores < 1L) {
    cores <- suppressWarnings(as.integer(parallel::detectCores(logical = TRUE)))
  }
  if (is.na(cores) || cores < 1L) cores <- 1L
  cores
}

detect_available_memory_gb <- function(available_cores) {
  env_gb <- trimws(Sys.getenv("DEG_AVAILABLE_MEMORY_GB", unset = ""))
  if (nzchar(env_gb)) {
    out <- suppressWarnings(as.numeric(env_gb))
    if (!is.na(out) && is.finite(out) && out > 0) return(out)
  }

  slurm_mem_node_mb <- trimws(Sys.getenv("SLURM_MEM_PER_NODE", unset = ""))
  if (nzchar(slurm_mem_node_mb)) {
    out <- suppressWarnings(as.numeric(slurm_mem_node_mb))
    if (!is.na(out) && is.finite(out) && out > 0) return(out / 1024)
  }

  slurm_mem_cpu_mb <- trimws(Sys.getenv("SLURM_MEM_PER_CPU", unset = ""))
  if (nzchar(slurm_mem_cpu_mb)) {
    out <- suppressWarnings(as.numeric(slurm_mem_cpu_mb))
    if (!is.na(out) && is.finite(out) && out > 0) return(out * available_cores / 1024)
  }

  meminfo <- "/proc/meminfo"
  if (file.exists(meminfo)) {
    lines <- readLines(meminfo, warn = FALSE)
    hit <- grep("^MemAvailable:", lines, value = TRUE)
    if (length(hit) > 0) {
      kb <- suppressWarnings(as.numeric(sub("^MemAvailable:\\s*([0-9]+)\\s+kB.*$", "\\1", hit[1])))
      if (!is.na(kb) && is.finite(kb) && kb > 0) return(kb / 1024^2)
    }
  }

  NA_real_
}

format_numeric_or_na <- function(x, digits = 2) {
  if (length(x) == 0 || is.na(x) || !is.finite(x)) return("NA")
  format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE)
}

plan_deg_workers <- function(worker_setting, tasks, obj, worker_memory_factor) {
  n_tasks <- length(tasks)
  available_cores <- detect_available_cores()
  available_memory_gb <- detect_available_memory_gb(available_cores)
  obj_size_gb <- as.numeric(utils::object.size(obj)) / 1024^3
  manual_workers <- parse_worker_setting(worker_setting)

  if (n_tasks == 0) {
    return(list(
      workers = 0L,
      mode = "none",
      setting = as.character(worker_setting)[1],
      requested_workers = 0L,
      available_cores = available_cores,
      cpu_limited_workers = 0L,
      available_memory_gb = available_memory_gb,
      object_size_gb = obj_size_gb,
      memory_limited_workers = 0L,
      worker_memory_factor = worker_memory_factor
    ))
  }

  if (!is.na(manual_workers)) {
    return(list(
      workers = max(1L, min(manual_workers, n_tasks)),
      mode = "manual",
      setting = as.character(worker_setting)[1],
      requested_workers = manual_workers,
      available_cores = available_cores,
      cpu_limited_workers = max(1L, min(available_cores, n_tasks)),
      available_memory_gb = available_memory_gb,
      object_size_gb = obj_size_gb,
      memory_limited_workers = NA_integer_,
      worker_memory_factor = worker_memory_factor
    ))
  }

  cpu_limited_workers <- max(1L, available_cores - 1L)
  memory_limited_workers <- cpu_limited_workers
  if (!is.na(available_memory_gb) && is.finite(available_memory_gb) && available_memory_gb > 0) {
    reserved_memory_gb <- max(2, min(16, available_memory_gb * 0.10))
    usable_memory_gb <- max(1, available_memory_gb - reserved_memory_gb)
    per_worker_memory_gb <- max(1, obj_size_gb * worker_memory_factor)
    memory_limited_workers <- max(1L, floor(usable_memory_gb / per_worker_memory_gb))
  }

  list(
    workers = max(1L, min(n_tasks, cpu_limited_workers, memory_limited_workers)),
    mode = "auto",
    setting = as.character(worker_setting)[1],
    requested_workers = NA_integer_,
    available_cores = available_cores,
    cpu_limited_workers = cpu_limited_workers,
    available_memory_gb = available_memory_gb,
    object_size_gb = obj_size_gb,
    memory_limited_workers = memory_limited_workers,
    worker_memory_factor = worker_memory_factor
  )
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

infer_ploidy_from_text <- function(values) {
  values <- trimws(as.character(values))
  values_upper <- toupper(values)
  out <- rep(NA_character_, length(values))
  out[!is.na(values_upper) & grepl("2N", values_upper, fixed = TRUE)] <- "2N"
  out[is.na(out) & !is.na(values_upper) & grepl("4N", values_upper, fixed = TRUE)] <- "4N"
  out
}

infer_ploidy_from_sample_label <- function(values) {
  values <- trimws(as.character(values))
  out <- infer_ploidy_from_text(values)
  out[is.na(out) & !is.na(values) & grepl("^A5", values, ignore.case = TRUE)] <- "4N"
  out[is.na(out) & !is.na(values) & grepl("^A6", values, ignore.case = TRUE)] <- "4N"
  out
}

infer_trajectory_context_from_ids <- function(values) {
  values <- trimws(as.character(values))
  values_lower <- tolower(values)
  ifelse(!is.na(values_lower) & grepl("cell-culture", values_lower, fixed = TRUE), "CellLine", "Tumor")
}

add_trajectory_subset_metadata <- function(obj) {
  meta <- obj@meta.data
  id_col <- resolve_col_case_insensitive(meta, id_col_candidates)
  if (is.na(id_col)) {
    stop(
      "A metadata column is required to derive CellLine/Tumor context. Tried: ",
      paste(unique(id_col_candidates), collapse = ", "),
      call. = FALSE
    )
  }

  sample_folder_col <- resolve_col_case_insensitive(meta, sample_folder_col_candidates)
  n_cells <- nrow(meta)
  ploidy_from_id <- infer_ploidy_from_text(meta[[id_col]])
  ploidy_from_sample_folder <- rep(NA_character_, n_cells)
  ploidy_from_sample <- rep(NA_character_, n_cells)

  if (!is.na(sample_folder_col)) {
    ploidy_from_sample_folder <- infer_ploidy_from_sample_label(meta[[sample_folder_col]])
  }
  if ("sample" %in% colnames(meta)) {
    ploidy_from_sample <- infer_ploidy_from_sample_label(meta$sample)
  }

  ploidy_values <- ploidy_from_id
  ploidy_source <- ifelse(!is.na(ploidy_from_id), id_col, NA_character_)
  missing_ploidy <- is.na(ploidy_values)
  ploidy_values[missing_ploidy] <- ploidy_from_sample_folder[missing_ploidy]
  ploidy_source[missing_ploidy & !is.na(ploidy_from_sample_folder)] <- sample_folder_col
  missing_ploidy <- is.na(ploidy_values)
  ploidy_values[missing_ploidy] <- ploidy_from_sample[missing_ploidy]
  ploidy_source[missing_ploidy & !is.na(ploidy_from_sample)] <- "sample"

  if (any(is.na(ploidy_values))) {
    tried_cols <- c(if (!is.na(id_col)) id_col, if (!is.na(sample_folder_col)) sample_folder_col, "sample")
    tried_cols <- unique(tried_cols[tried_cols %in% colnames(meta)])
    bad_cells <- unique(rownames(meta)[is.na(ploidy_values)])
    stop(
      "Could not infer ploidy for 03c_DEGs_round2.R. Tried metadata column(s): ",
      paste(tried_cols, collapse = ", "),
      ". Example unresolved cell(s): ",
      paste(head(bad_cells, 20), collapse = ", "),
      call. = FALSE
    )
  }

  meta[[ploidy_col]] <- factor(ploidy_values, levels = c("2N", "4N"))
  meta[["ploidy_source"]] <- ploidy_source
  meta[[tn_col]] <- factor(infer_trajectory_context_from_ids(meta[[id_col]]), levels = c("CellLine", "Tumor"))
  meta[["trajectory_context"]] <- meta[[tn_col]]

  obj@meta.data <- meta
  obj
}

get_analysis_group_cells <- function(meta, tn_scope, ploidy_scope) {
  idx <- rep(TRUE, nrow(meta))
  if (!identical(tn_scope, "All")) {
    idx <- idx & as.character(meta[[tn_col]]) == tn_scope
  }
  if (ploidy_scope %in% c("2N", "4N")) {
    idx <- idx & as.character(meta[[ploidy_col]]) == ploidy_scope
  }
  rownames(meta)[idx %in% TRUE]
}

build_analysis_group_df <- function(meta) {
  context_df <- data.frame(
    context_label = c("All_cells", "Tumor", "CellLine"),
    tn_scope = c("All", "Tumor", "CellLine"),
    stringsAsFactors = FALSE
  )
  ploidy_scopes <- c("ploidy_all", "2N", "4N")
  rows <- list()

  for (i in seq_len(nrow(context_df))) {
    for (ploidy_scope in ploidy_scopes) {
      context_label <- context_df$context_label[i]
      tn_scope <- context_df$tn_scope[i]
      group_name <- if (identical(context_label, "All_cells") && identical(ploidy_scope, "ploidy_all")) {
        "All_cells"
      } else {
        paste(context_label, ploidy_scope, sep = "_")
      }
      cells <- get_analysis_group_cells(meta, tn_scope, ploidy_scope)
      rows[[length(rows) + 1L]] <- data.frame(
        analysis_group = group_name,
        tn_scope = tn_scope,
        ploidy_scope = ploidy_scope,
        context_label = context_label,
        group_safe = file.path(context_label, ploidy_scope),
        n_cells = length(cells),
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(rows)
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
  analysis_group,
  tn_scope,
  ploidy_scope,
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
    analysis_group = analysis_group,
    tn_scope = tn_scope,
    ploidy_scope = ploidy_scope,
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
  analysis_group,
  tn_scope,
  ploidy_scope,
  cluster = NA_character_,
  subset_group = NA_character_,
  dose = NA_character_,
  min_cells = 3L
) {
  if (ncol(obj) == 0) {
    return(make_summary_row(
      scope = scope,
      analysis_group = analysis_group,
      tn_scope = tn_scope,
      ploidy_scope = ploidy_scope,
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
      analysis_group = analysis_group,
      tn_scope = tn_scope,
      ploidy_scope = ploidy_scope,
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
      analysis_group = analysis_group,
      tn_scope = tn_scope,
      ploidy_scope = ploidy_scope,
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
      analysis_group = analysis_group,
      tn_scope = tn_scope,
      ploidy_scope = ploidy_scope,
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
    analysis_group = analysis_group,
    tn_scope = tn_scope,
    ploidy_scope = ploidy_scope,
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
        analysis_group = analysis_group,
        tn_scope = tn_scope,
        ploidy_scope = ploidy_scope,
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
        analysis_group = analysis_group,
        tn_scope = tn_scope,
        ploidy_scope = ploidy_scope,
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

make_findmarkers_task <- function(
  cells = NULL,
  group_col,
  ident_1,
  ident_2 = NULL,
  comparison,
  output_file,
  assay,
  scope,
  analysis_group,
  tn_scope,
  ploidy_scope,
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
    analysis_group = analysis_group,
    tn_scope = tn_scope,
    ploidy_scope = ploidy_scope,
    cluster = cluster,
    subset_group = subset_group,
    dose = dose,
    min_cells = min_cells
  )
}

run_deg_task <- function(task) {
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
      analysis_group = task$analysis_group,
      tn_scope = task$tn_scope,
      ploidy_scope = task$ploidy_scope,
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
    analysis_group = task$analysis_group,
    tn_scope = task$tn_scope,
    ploidy_scope = task$ploidy_scope,
    cluster = task$cluster,
    subset_group = task$subset_group,
    dose = task$dose,
    min_cells = task$min_cells
  )
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

write_cluster_counts_by_analysis_group <- function(meta, analysis_group_df, analysis_group_cells, cluster_levels, output_file) {
  rows <- lapply(seq_len(nrow(analysis_group_df)), function(i) {
    group_name <- analysis_group_df$analysis_group[i]
    cells <- analysis_group_cells[[group_name]]
    counts <- table(factor(as.character(meta[cells, cluster_col, drop = TRUE]), levels = cluster_levels))
    data.frame(
      analysis_group = group_name,
      tn_scope = analysis_group_df$tn_scope[i],
      ploidy_scope = analysis_group_df$ploidy_scope[i],
      cluster = names(counts),
      n_cells = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
  write_table_csv(dplyr::bind_rows(rows), output_file)
}

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_cluster_vs_rest <- .ensure_dir(file.path(output_root, "01_cluster_final_vs_rest"))

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

cluster_source_col <- cluster_col
if (!(cluster_col %in% colnames(obj@meta.data))) {
  if (cluster_fallback_col %in% colnames(obj@meta.data)) {
    message("Metadata column '", cluster_col, "' is missing; deriving it from '", cluster_fallback_col, "'.")
    obj@meta.data[[cluster_col]] <- as.character(obj@meta.data[[cluster_fallback_col]])
    cluster_source_col <- cluster_fallback_col
  } else {
    stop(
      "Cannot find required cluster metadata column '",
      cluster_col,
      "' or fallback column '",
      cluster_fallback_col,
      "'.",
      call. = FALSE
    )
  }
}
assert_required_column(obj@meta.data, cluster_col)

message("Preparing trajectory-consistent subset metadata.")
obj <- add_trajectory_subset_metadata(obj)
assert_required_column(obj@meta.data, ploidy_col)
assert_required_column(obj@meta.data, tn_col)

message("Preparing RNA assay for Seurat::FindMarkers.")
obj <- prepare_obj_for_markers(obj, assay = assay_use)
DefaultAssay(obj) <- assay_use

cluster_levels <- metadata_levels(obj@meta.data, cluster_col)
if (length(cluster_levels) == 0) {
  stop("No non-missing cluster_final values were found.", call. = FALSE)
}
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
obj@meta.data[["trajectory_context"]] <- obj@meta.data[[tn_col]]

analysis_group_df_all <- build_analysis_group_df(obj@meta.data)
write_table_csv(analysis_group_df_all, file.path(out_summary, "analysis_group_cell_counts.csv"))
analysis_group_df <- analysis_group_df_all %>%
  dplyr::filter(.data$n_cells > 0)
if (nrow(analysis_group_df) == 0) {
  stop("No cells were found for the requested DEG analysis groups.", call. = FALSE)
}

analysis_group_cells <- stats::setNames(vector("list", nrow(analysis_group_df)), analysis_group_df$analysis_group)
for (i in seq_len(nrow(analysis_group_df))) {
  analysis_group_cells[[analysis_group_df$analysis_group[i]]] <- get_analysis_group_cells(
    obj@meta.data,
    analysis_group_df$tn_scope[i],
    analysis_group_df$ploidy_scope[i]
  )
}

metadata_count_list <- list(
  cluster_final = data.frame(
    metadata = cluster_col,
    value = names(table(obj@meta.data[[cluster_col]], useNA = "no")),
    n_cells = as.integer(table(obj@meta.data[[cluster_col]], useNA = "no")),
    stringsAsFactors = FALSE
  ),
  Ploidy = data.frame(
    metadata = ploidy_col,
    value = names(table(obj@meta.data[[ploidy_col]], useNA = "no")),
    n_cells = as.integer(table(obj@meta.data[[ploidy_col]], useNA = "no")),
    stringsAsFactors = FALSE
  ),
  TN = data.frame(
    metadata = tn_col,
    value = names(table(obj@meta.data[[tn_col]], useNA = "no")),
    n_cells = as.integer(table(obj@meta.data[[tn_col]], useNA = "no")),
    stringsAsFactors = FALSE
  )
)
write_table_csv(dplyr::bind_rows(metadata_count_list), file.path(out_summary, "metadata_cell_counts.csv"))
write_cluster_counts_by_analysis_group(
  meta = obj@meta.data,
  analysis_group_df = analysis_group_df,
  analysis_group_cells = analysis_group_cells,
  cluster_levels = cluster_levels,
  output_file = file.path(out_summary, "cluster_final_cell_counts_by_analysis_group.csv")
)

deg_tasks <- list()
add_task <- function(task) {
  task$task_id <- length(deg_tasks) + 1L
  deg_tasks[[task$task_id]] <<- task
  invisible(NULL)
}

message("Preparing cluster_final vs rest DEG tasks by trajectory-style subsets.")
for (i in seq_len(nrow(analysis_group_df))) {
  group_name <- analysis_group_df$analysis_group[i]
  group_cells <- analysis_group_cells[[group_name]]
  group_dir <- .ensure_dir(file.path(out_cluster_vs_rest, analysis_group_df$group_safe[i]))
  subset_group <- if (identical(analysis_group_df$tn_scope[i], "All")) "All_cells" else analysis_group_df$tn_scope[i]

  for (cluster_id in cluster_levels) {
    cluster_stub <- safe_file_component(cluster_id)
    add_task(make_findmarkers_task(
      cells = group_cells,
      group_col = cluster_col,
      ident_1 = cluster_id,
      ident_2 = NULL,
      comparison = paste0(group_name, "_cluster_final_", cluster_id, "_vs_rest"),
      output_file = file.path(group_dir, paste0("cluster_final_", cluster_stub, "_vs_rest_markers.csv")),
      assay = assay_use,
      scope = "cluster_final_vs_rest_by_trajectory_subset",
      analysis_group = group_name,
      tn_scope = analysis_group_df$tn_scope[i],
      ploidy_scope = analysis_group_df$ploidy_scope[i],
      cluster = cluster_id,
      subset_group = subset_group,
      min_cells = min_cells_per_group
    ))
  }
}

deg_worker_plan <- plan_deg_workers(
  worker_setting = deg_parallel_workers_setting,
  tasks = deg_tasks,
  obj = obj,
  worker_memory_factor = deg_worker_memory_factor
)
deg_parallel_workers_used <- deg_worker_plan$workers
message(
  "Prepared ", length(deg_tasks), " DEG tasks. ",
  "Worker mode: ", deg_worker_plan$mode,
  "; available cores: ", deg_worker_plan$available_cores,
  "; available memory GB: ", format_numeric_or_na(deg_worker_plan$available_memory_gb),
  "; Seurat object size GB: ", format_numeric_or_na(deg_worker_plan$object_size_gb),
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
  "03c DEGs round2 completed.",
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Assay: ", assay_use),
  "DEG method: Seurat::FindMarkers with default thresholds and default test settings.",
  paste0("DEG parallel backend: future::multisession when workers > 1, otherwise sequential."),
  paste0("DEG parallel worker setting: ", deg_worker_plan$setting),
  paste0("DEG parallel worker mode: ", deg_worker_plan$mode),
  paste0("DEG parallel workers used: ", deg_parallel_workers_used),
  paste0("Available CPU cores detected: ", deg_worker_plan$available_cores),
  paste0("CPU-limited workers: ", deg_worker_plan$cpu_limited_workers),
  paste0("Available memory GB detected: ", format_numeric_or_na(deg_worker_plan$available_memory_gb)),
  paste0("Seurat object size GB: ", format_numeric_or_na(deg_worker_plan$object_size_gb)),
  paste0("Memory-limited workers: ", ifelse(is.na(deg_worker_plan$memory_limited_workers), "NA", deg_worker_plan$memory_limited_workers)),
  paste0("DEG worker memory factor: ", deg_worker_plan$worker_memory_factor),
  paste0("Future globals max size GB: ", future_globals_maxsize_gb),
  paste0("Minimum cells per group guard: ", min_cells_per_group),
  paste0("Cluster column: ", cluster_col),
  paste0("Cluster source column: ", cluster_source_col),
  paste0("Trajectory-style TN scopes: All, Tumor, CellLine."),
  paste0("Trajectory-style ploidy scopes in each TN scope: ploidy_all, 2N, 4N."),
  paste0("Analysis groups with cells: ", paste(analysis_group_df$analysis_group, collapse = ", ")),
  paste0("Total comparisons recorded: ", nrow(comparison_summary_df)),
  paste0("Successful comparisons: ", sum(comparison_summary_df$status == "success", na.rm = TRUE)),
  paste0("Cached comparisons skipped: ", sum(comparison_summary_df$status == "cached", na.rm = TRUE)),
  paste0("Skipped comparisons: ", sum(comparison_summary_df$status == "skipped", na.rm = TRUE)),
  paste0("Errored comparisons: ", sum(comparison_summary_df$status == "error", na.rm = TRUE))
)
writeLines(run_summary, con = file.path(out_summary, "run_summary.txt"))

message("03c_DEGs_round2 finished.")
message("Output root: ", output_root)
