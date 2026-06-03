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
  "ggplot2",
  "Matrix",
  "readr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})

trajectory_character_cols <- c(
  "cell",
  "group_1",
  "group_2",
  "Dose",
  "ploidy",
  "Ploidy",
  "TN",
  "Dose_DEG",
  "trajectory_context",
  "trajectory_group",
  "trajectory_analysis_group",
  "trajectory_shape_group",
  "trajectory_tn_scope",
  "trajectory_ploidy_scope",
  "trajectory_branch",
  "cluster_final",
  "clusters",
  "sample",
  "sample_type",
  "cluster_final_annotation_primary"
)

coerce_trajectory_character_cols <- function(df) {
  for (col in intersect(trajectory_character_cols, colnames(df))) {
    df[[col]] <- as.character(df[[col]])
  }
  df
}

read_trajectory_csv <- function(path, ...) {
  coerce_trajectory_character_cols(readr::read_csv(path, show_col_types = FALSE, ...))
}

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

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- tolower(trimws(as.character(x)[1]))
  if (!nzchar(x)) return(default)
  if (x %in% c("1", "true", "t", "yes", "y")) return(TRUE)
  if (x %in% c("0", "false", "f", "no", "n")) return(FALSE)
  stop("Cannot parse logical value: ", x, call. = FALSE)
}

parse_int_cfg <- function(value, default, name) {
  if (is.null(value) || !nzchar(trimws(as.character(value)))) return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  out
}

parse_numeric_cfg <- function(value, default, name) {
  if (is.null(value) || !nzchar(trimws(as.character(value)))) return(default)
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out) || !is.finite(out) || out <= 0) {
    stop(name, " must be a positive number.", call. = FALSE)
  }
  out
}

collapse_unique_nonmissing <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x) & x != "NA"]
  x <- sort(unique(x))
  if (length(x) == 0) NA_character_ else paste(x, collapse = ";")
}

semicolon_contains <- function(x, value) {
  value <- as.character(value)[1]
  vapply(
    strsplit(as.character(x), ";", fixed = TRUE),
    function(parts) value %in% trimws(parts),
    logical(1)
  )
}

resolve_first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) normalizePath(hit[1], mustWork = TRUE) else paths[1]
}

resolve_first_existing_col <- function(df, candidates) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) hit[1] else NA_character_
}

make_pretty_umap_theme <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 3),
      plot.subtitle = element_text(color = "grey35"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
}

save_both <- function(plot_obj, file_stub, width = 9, height = 7, dpi = 300) {
  ggsave(paste0(file_stub, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(file_stub, ".png"), plot_obj, width = width, height = height, dpi = dpi)
}

find_velocity_files <- function(root_dir) {
  if (is.null(root_dir) || !nzchar(root_dir) || !dir.exists(root_dir)) {
    return(character(0))
  }
  list.files(
    root_dir,
    pattern = "\\.(loom|h5ad)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
}

find_first_existing_file <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) normalizePath(hit[1], mustWork = TRUE) else NA_character_
}

find_cellranger_bam <- function(sample_dir) {
  candidates <- c(
    file.path(sample_dir, "outs", "possorted_genome_bam.bam"),
    file.path(sample_dir, "possorted_genome_bam.bam")
  )
  bam <- find_first_existing_file(candidates)
  if (!is.na(bam)) return(bam)
  if (!dir.exists(sample_dir)) return(NA_character_)
  recursive_hits <- list.files(
    sample_dir,
    pattern = "^possorted_genome_bam\\.bam$",
    full.names = TRUE,
    recursive = TRUE
  )
  find_first_existing_file(recursive_hits)
}

find_cellranger_barcodes <- function(sample_dir) {
  candidates <- c(
    file.path(sample_dir, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz"),
    file.path(sample_dir, "outs", "filtered_feature_bc_matrix", "barcodes.tsv"),
    file.path(sample_dir, "outs", "barcodes.tsv.gz"),
    file.path(sample_dir, "outs", "barcodes.tsv")
  )
  barcodes <- find_first_existing_file(candidates)
  if (!is.na(barcodes)) return(barcodes)
  if (!dir.exists(sample_dir)) return(NA_character_)
  recursive_hits <- list.files(
    sample_dir,
    pattern = "^barcodes\\.tsv(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE
  )
  find_first_existing_file(recursive_hits)
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

find_cellranger_loom <- function(sample_dir, sample_folder, extra_dirs = character(0)) {
  sample_name <- basename(sample_folder)
  candidates <- c(
    file.path(extra_dirs, paste0(sample_name, ".loom")),
    file.path(extra_dirs, sample_name, paste0(sample_name, ".loom")),
    file.path(sample_dir, "velocyto", paste0(sample_name, ".loom")),
    file.path(sample_dir, "outs", "velocyto", paste0(sample_name, ".loom"))
  )
  loom <- find_first_existing_file(candidates)
  if (!is.na(loom)) return(loom)

  search_dirs <- c(file.path(sample_dir, "velocyto"), file.path(sample_dir, "outs", "velocyto"), extra_dirs)
  search_dirs <- search_dirs[dir.exists(search_dirs)]
  if (length(search_dirs) == 0) return(NA_character_)
  loom_hits <- unlist(lapply(search_dirs, function(path) {
    list.files(path, pattern = "\\.loom$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  }), use.names = FALSE)
  find_first_existing_file(loom_hits)
}

build_cellranger_sample_manifest <- function(
  meta_all,
  cellranger_root,
  sample_folder_col,
  dose_col,
  ploidy_col,
  velocyto_work_root,
  velocyto_output_root,
  velocyto_run_tag
) {
  if (!(sample_folder_col %in% colnames(meta_all))) {
    return(data.frame())
  }

  sample_df <- meta_all %>%
    dplyr::filter(!is.na(.data[[sample_folder_col]]), .data[[sample_folder_col]] != "") %>%
    dplyr::group_by(.data[[sample_folder_col]]) %>%
    dplyr::summarise(
      Dose = collapse_unique_nonmissing(.data[[dose_col]]),
      ploidy = if (ploidy_col %in% colnames(meta_all)) collapse_unique_nonmissing(.data[[ploidy_col]]) else NA_character_,
      trajectory_context = if ("trajectory_context" %in% colnames(meta_all)) collapse_unique_nonmissing(.data[["trajectory_context"]]) else NA_character_,
      trajectory_group = if ("trajectory_group" %in% colnames(meta_all)) collapse_unique_nonmissing(.data[["trajectory_group"]]) else NA_character_,
      sample = if ("sample" %in% colnames(meta_all)) collapse_unique_nonmissing(.data[["sample"]]) else NA_character_,
      sample_type = if ("sample_type" %in% colnames(meta_all)) collapse_unique_nonmissing(.data[["sample_type"]]) else NA_character_,
      n_cells_in_seurat = dplyr::n(),
      .groups = "drop"
    )
  colnames(sample_df)[colnames(sample_df) == sample_folder_col] <- "sample_folder"

  if (nrow(sample_df) == 0) return(data.frame())

  sample_df$cellranger_dir <- file.path(cellranger_root, sample_df$sample_folder)
  sample_df$velocyto_output_dir <- file.path(velocyto_output_root, sanitize_path_component(sample_df$sample_folder))
  sample_df$velocyto_work_sample_root <- file.path(
    velocyto_work_root,
    sanitize_path_component(sample_df$sample_folder)
  )
  sample_df$velocyto_run_tag <- sanitize_path_component(velocyto_run_tag)[1]
  sample_df$velocyto_writable_sample_dir <- file.path(
    sample_df$velocyto_work_sample_root,
    sample_df$velocyto_run_tag,
    basename(sample_df$sample_folder)
  )
  sample_df$cellranger_dir_exists <- dir.exists(sample_df$cellranger_dir)
  sample_df$bam_file <- vapply(sample_df$cellranger_dir, find_cellranger_bam, character(1))
  sample_df$bam_exists <- !is.na(sample_df$bam_file) & file.exists(sample_df$bam_file)
  sample_df$barcodes_file <- vapply(sample_df$cellranger_dir, find_cellranger_barcodes, character(1))
  sample_df$barcodes_exists <- !is.na(sample_df$barcodes_file) & file.exists(sample_df$barcodes_file)
  sample_df$velocyto_loom_file <- vapply(
    seq_len(nrow(sample_df)),
    function(i) find_cellranger_loom(
      sample_df$cellranger_dir[i],
      sample_df$sample_folder[i],
      extra_dirs = c(
        sample_df$velocyto_output_dir[i],
        sample_df$velocyto_work_sample_root[i],
        file.path(sample_df$velocyto_writable_sample_dir[i], "velocyto")
      )
    ),
    character(1)
  )
  sample_df$velocyto_loom_exists <- !is.na(sample_df$velocyto_loom_file) & file.exists(sample_df$velocyto_loom_file)
  sample_df
}

write_velocyto_runner <- function(
  run_script,
  velocyto_bin,
  samtools_bin,
  samtools_threads,
  source_sample_dir,
  bam_file,
  barcodes_file,
  writable_sample_dir,
  gtf_file,
  output_dir,
  log_file
) {
  .ensure_dir(dirname(run_script))
  if (
    !is.null(gtf_file) && nzchar(gtf_file) && file.exists(gtf_file) &&
      !is.na(bam_file) && nzchar(bam_file) && file.exists(bam_file) &&
      !is.na(barcodes_file) && nzchar(barcodes_file) && file.exists(barcodes_file) &&
      !is.na(samtools_bin) && nzchar(samtools_bin) && file.exists(samtools_bin)
  ) {
    command_lines <- c(
      paste(
        "export PATH=",
        shQuote(dirname(samtools_bin)),
        ":$PATH",
        sep = ""
      ),
      paste(
        "mkdir -p",
        shQuote(output_dir),
        shQuote(dirname(writable_sample_dir))
      ),
      paste(
        "mkdir -p",
        shQuote(writable_sample_dir),
        shQuote(file.path(writable_sample_dir, "outs"))
      ),
      paste(
        "find",
        shQuote(file.path(source_sample_dir, "outs")),
        "-mindepth 1 -maxdepth 1 -exec ln -s {}",
        shQuote(file.path(writable_sample_dir, "outs")),
        "\\;"
      ),
      paste(
        "if [ ! -e",
        shQuote(file.path(writable_sample_dir, "outs", "filtered_feature_bc_matrix", "barcodes.tsv")),
        "] && [ ! -e",
        shQuote(file.path(writable_sample_dir, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")),
        "]; then mkdir -p",
        shQuote(file.path(writable_sample_dir, "outs", "filtered_feature_bc_matrix")),
        "; ln -sf",
        shQuote(barcodes_file),
        shQuote(file.path(writable_sample_dir, "outs", "filtered_feature_bc_matrix")),
        "; fi"
      ),
      paste(
        "ln -sf",
        shQuote(bam_file),
        shQuote(file.path(writable_sample_dir, "outs", "possorted_genome_bam.bam"))
      ),
      paste(
        "chmod u+rwx",
        shQuote(writable_sample_dir)
      ),
      paste(
        shQuote(velocyto_bin),
        "run10x",
        "--samtools-threads",
        as.character(samtools_threads),
        shQuote(writable_sample_dir),
        shQuote(gtf_file),
        ">",
        shQuote(log_file),
        "2>&1"
      ),
      paste(
        "for loom in",
        paste0(shQuote(file.path(writable_sample_dir, "velocyto")), "/*.loom;"),
        "do [ -e \"$loom\" ] && cp -f \"$loom\"",
        shQuote(output_dir),
        "|| true; done"
      )
    )
  } else {
    reason <- if (is.null(gtf_file) || !nzchar(gtf_file) || !file.exists(gtf_file)) {
      "VELOCYTO_GTF is not set or does not exist. Set VELOCYTO_GTF to an uncompressed gene annotation .gtf file before generating loom files."
    } else if (is.na(samtools_bin) || !nzchar(samtools_bin) || !file.exists(samtools_bin)) {
      "samtools was not found. Set SAMTOOLS_BIN to an existing samtools executable."
    } else if (is.na(bam_file) || !nzchar(bam_file) || !file.exists(bam_file)) {
      paste0("Cell Ranger BAM was not found for sample folder: ", source_sample_dir)
    } else if (is.na(barcodes_file) || !nzchar(barcodes_file) || !file.exists(barcodes_file)) {
      paste0("Cell Ranger barcode file was not found for sample folder: ", source_sample_dir)
    } else {
      paste0("Cell Ranger inputs were not found for sample folder: ", source_sample_dir)
    }
    command_lines <- paste("echo", shQuote(reason), ">&2; exit 2")
  }

  writeLines(c("#!/usr/bin/env bash", "set -euo pipefail", command_lines), con = run_script)
  Sys.chmod(run_script, mode = "0755")
}

prepare_velocyto_runners <- function(sample_manifest, output_dir, velocyto_bin, samtools_bin, samtools_threads, gtf_file) {
  if (is.null(sample_manifest) || nrow(sample_manifest) == 0) return(sample_manifest)

  sample_manifest$velocyto_samtools_threads <- samtools_threads
  runner_dir <- .ensure_dir(file.path(output_dir, "velocyto_run10x"))
  sample_manifest$velocyto_run_script <- file.path(
    runner_dir,
    paste0(sanitize_path_component(sample_manifest$sample_folder), "_run_velocyto.sh")
  )
  sample_manifest$velocyto_log_file <- file.path(
    runner_dir,
    paste0(sanitize_path_component(sample_manifest$sample_folder), "_velocyto.log")
  )

  for (i in seq_len(nrow(sample_manifest))) {
    write_velocyto_runner(
      run_script = sample_manifest$velocyto_run_script[i],
      velocyto_bin = velocyto_bin,
      samtools_bin = samtools_bin,
      samtools_threads = samtools_threads,
      source_sample_dir = sample_manifest$cellranger_dir[i],
      bam_file = sample_manifest$bam_file[i],
      barcodes_file = sample_manifest$barcodes_file[i],
      writable_sample_dir = sample_manifest$velocyto_writable_sample_dir[i],
      gtf_file = gtf_file,
      output_dir = sample_manifest$velocyto_output_dir[i],
      log_file = sample_manifest$velocyto_log_file[i]
    )
  }
  sample_manifest
}

resolve_sample_loom_inputs_by_dose <- function(sample_manifest, dose_levels) {
  out <- stats::setNames(vector("list", length(dose_levels)), dose_levels)
  for (dose in dose_levels) {
    rows <- sample_manifest[semicolon_contains(sample_manifest$Dose, dose), , drop = FALSE]
    if (nrow(rows) > 0 && all(rows$velocyto_loom_exists)) {
      out[[dose]] <- normalizePath(rows$velocyto_loom_file, mustWork = FALSE)
    } else {
      out[[dose]] <- character(0)
    }
  }
  out
}

resolve_sample_loom_inputs_by_group <- function(sample_manifest, group_df) {
  out <- stats::setNames(vector("list", nrow(group_df)), group_df$trajectory_group)
  for (i in seq_len(nrow(group_df))) {
    if ("trajectory_group" %in% colnames(sample_manifest)) {
      rows <- sample_manifest[
        semicolon_contains(sample_manifest$trajectory_group, group_df$trajectory_group[i]),
        ,
        drop = FALSE
      ]
    } else {
      rows <- sample_manifest[
        semicolon_contains(sample_manifest$Dose, group_df$Dose[i]) &
          semicolon_contains(sample_manifest$ploidy, group_df$ploidy[i]),
        ,
        drop = FALSE
      ]
    }
    if (nrow(rows) > 0 && all(rows$velocyto_loom_exists)) {
      out[[group_df$trajectory_group[i]]] <- normalizePath(rows$velocyto_loom_file, mustWork = FALSE)
    } else {
      out[[group_df$trajectory_group[i]]] <- character(0)
    }
  }
  out
}

resolve_velocity_inputs <- function(dose_levels, manifest_path, velocity_input_root) {
  dose_levels <- as.character(dose_levels)
  empty <- stats::setNames(rep(NA_character_, length(dose_levels)), dose_levels)

  if (!is.null(manifest_path) && nzchar(manifest_path) && file.exists(manifest_path)) {
    manifest <- readr::read_csv(manifest_path, show_col_types = FALSE)
    required_cols <- c("dose", "input_file")
    missing_cols <- setdiff(required_cols, colnames(manifest))
    if (length(missing_cols) > 0) {
      stop("Velocity manifest is missing required column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    manifest$dose <- standardize_in_vivo_dose(manifest$dose)
    base_dir <- dirname(normalizePath(manifest_path, mustWork = TRUE))
    for (dose in dose_levels) {
      rows <- manifest[manifest$dose == dose, , drop = FALSE]
      if (nrow(rows) >= 1) {
        input_files <- rows$input_file
        input_files <- ifelse(grepl("^/", input_files), input_files, file.path(base_dir, input_files))
        empty[[dose]] <- paste(normalizePath(input_files, mustWork = FALSE), collapse = ",")
      }
    }
    return(empty)
  }

  files <- find_velocity_files(velocity_input_root)
  if (length(files) == 0) {
    return(empty)
  }
  files <- normalizePath(files, mustWork = FALSE)
  if (length(files) == 1) {
    empty[] <- files[1]
    return(empty)
  }

  file_base <- tolower(basename(files))
  for (dose in dose_levels) {
    dose_safe <- tolower(sanitize_path_component(gsub("/", "", dose)))
    dose_num <- sub("mg/kg$", "", dose)
    idx <- grepl(dose_safe, file_base, fixed = TRUE) |
      grepl(paste0(dose_num, "mg"), file_base, fixed = TRUE) |
      grepl(paste0("dose", dose_num), file_base, fixed = TRUE)
    if (sum(idx) >= 1) {
      empty[[dose]] <- paste(files[idx], collapse = ",")
    }
  }
  empty
}

as_velocity_input_list <- function(velocity_inputs) {
  stats::setNames(
    lapply(velocity_inputs, function(x) {
      x <- as.character(x)
      x <- x[!is.na(x) & nzchar(x)]
      if (length(x) == 0) return(character(0))
      x <- unlist(strsplit(x, "[,;]", perl = TRUE), use.names = FALSE)
      x <- trimws(x)
      x <- x[nzchar(x)]
      x
    }),
    names(velocity_inputs)
  )
}

format_input_files <- function(files) {
  files <- as.character(files)
  files <- files[!is.na(files) & nzchar(files)]
  paste(files, collapse = ";")
}

all_input_files_exist <- function(files) {
  files <- as.character(files)
  files <- files[!is.na(files) & nzchar(files)]
  length(files) > 0 && all(file.exists(files))
}

prefer_plain_gtf <- function(gtf_path) {
  gtf_path <- path.expand(as.character(gtf_path)[1])
  if (grepl("\\.gz$", gtf_path, ignore.case = TRUE)) {
    plain_gtf <- sub("\\.gz$", "", gtf_path, ignore.case = TRUE)
    if (file.exists(plain_gtf)) return(plain_gtf)
  }
  gtf_path
}

make_scvelo_input_arg <- function(files) {
  files <- as.character(files)
  files <- files[!is.na(files) & nzchar(files)]
  paste(files, collapse = ",")
}

get_named_input_files <- function(input_list, name) {
  out <- input_list[[as.character(name)]]
  if (is.null(out)) return(character(0))
  out <- as.character(out)
  out[!is.na(out) & nzchar(out)]
}

resolve_python_from_env_name <- function(env_name) {
  virtual_env <- Sys.getenv("VIRTUAL_ENV", unset = "")
  virtual_env_python <- if (nzchar(virtual_env) && basename(virtual_env) == env_name) {
    file.path(virtual_env, "bin", "python")
  } else {
    ""
  }
  pyenv_root <- Sys.getenv("PYENV_ROOT", unset = file.path(Sys.getenv("HOME", unset = ""), ".pyenv"))
  candidates <- c(
    Sys.getenv(paste0(toupper(env_name), "_PYTHON"), unset = ""),
    virtual_env_python,
    file.path(Sys.getenv("WORKON_HOME", unset = ""), env_name, "bin", "python"),
    file.path(Sys.getenv("HOME", unset = ""), "venvs", env_name, "bin", "python"),
    file.path("/home/4482173/venvs", env_name, "bin", "python"),
    file.path(pyenv_root, "versions", env_name, "bin", "python"),
    file.path(pyenv_root, "versions", env_name, "bin", "python3"),
    file.path(Sys.getenv("CONDA_PREFIX", unset = ""), "envs", env_name, "bin", "python"),
    file.path(Sys.getenv("HOME", unset = ""), "miniconda3", "envs", env_name, "bin", "python"),
    file.path(Sys.getenv("HOME", unset = ""), "anaconda3", "envs", env_name, "bin", "python"),
    file.path(Sys.getenv("HOME", unset = ""), "mambaforge", "envs", env_name, "bin", "python"),
    file.path(Sys.getenv("HOME", unset = ""), "micromamba", "envs", env_name, "bin", "python"),
    file.path(Sys.getenv("MAMBA_ROOT_PREFIX", unset = ""), "envs", env_name, "bin", "python")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) > 0) {
    return(path.expand(hit[1]))
  }

  pyenv_bin <- Sys.which("pyenv")
  if (nzchar(pyenv_bin)) {
    pyenv_prefix <- tryCatch(
      system2(pyenv_bin, c("prefix", env_name), stdout = TRUE, stderr = TRUE),
      error = function(e) character(0)
    )
    pyenv_prefix <- pyenv_prefix[file.exists(pyenv_prefix)]
    if (length(pyenv_prefix) > 0) {
      python_path <- file.path(pyenv_prefix[1], "bin", "python")
      if (file.exists(python_path)) {
        return(path.expand(python_path))
      }
    }
  }

  conda_bin <- Sys.which("conda")
  if (nzchar(conda_bin)) {
    env_list <- tryCatch(
      system2(conda_bin, c("env", "list"), stdout = TRUE, stderr = TRUE),
      error = function(e) character(0)
    )
    env_line <- env_list[grepl(paste0("(^|[[:space:]])", env_name, "[[:space:]]"), env_list)]
    if (length(env_line) > 0) {
      env_path <- sub("^.*[[:space:]](/.*)$", "\\1", env_line[1])
      python_path <- file.path(env_path, "bin", "python")
      if (file.exists(python_path)) {
        return(path.expand(python_path))
      }
    }
  }

  NA_character_
}

resolve_python_executable <- function(python_bin) {
  python_bin <- path.expand(as.character(python_bin)[1])
  if (file.exists(python_bin)) return(python_bin)
  Sys.which(python_bin)
}

resolve_executable_path <- function(x) {
  x <- path.expand(as.character(x)[1])
  if (!is.na(x) && nzchar(x) && file.exists(x)) return(x)
  hit <- Sys.which(x)
  if (nzchar(hit)) hit else NA_character_
}

python_matches_required_env <- function(python_path, env_name) {
  if (!nzchar(python_path)) return(FALSE)
  path_parts <- strsplit(path.expand(python_path), .Platform$file.sep, fixed = TRUE)[[1]]
  if (env_name %in% path_parts) return(TRUE)

  prefix <- tryCatch(
    system2(
      python_path,
      c("-c", "import sys; print(sys.prefix)"),
      stdout = TRUE,
      stderr = FALSE
    ),
    error = function(e) character(0)
  )
  prefix <- prefix[!is.na(prefix) & nzchar(prefix)]
  if (length(prefix) == 0) return(FALSE)
  prefix_parts <- strsplit(path.expand(prefix[1]), .Platform$file.sep, fixed = TRUE)[[1]]
  env_name %in% prefix_parts || basename(prefix[1]) == env_name
}

assert_required_python_env <- function(python_bin, env_name) {
  python_path <- resolve_python_executable(python_bin)
  if (!nzchar(python_path)) {
    stop(
      "Python binary for required environment '",
      env_name,
      "' was not found. Set SCVELO_PYTHON to .../",
      env_name,
      "/bin/python.",
      call. = FALSE
    )
  }

  if (!python_matches_required_env(python_path, env_name)) {
    stop(
      "04_trajectory.R must use Python environment '",
      env_name,
      "', but resolved Python is: ",
      python_path,
      ". Set SCVELO_PYTHON to the python inside ",
      env_name,
      ". For pyenv virtualenv, use: Sys.setenv(SCVELO_PYTHON = path.expand('~/.pyenv/versions/",
      env_name,
      "/bin/python')).",
      call. = FALSE
    )
  }

  python_path
}

select_required_python_env <- function(configured_python_bin, env_name) {
  python_path <- resolve_python_executable(configured_python_bin)

  if (nzchar(python_path) && python_matches_required_env(python_path, env_name)) {
    return(python_path)
  }

  env_python <- resolve_python_from_env_name(env_name)
  if (!is.na(env_python) && python_matches_required_env(env_python, env_name)) {
    if (nzchar(python_path)) {
      warning(
        "Ignoring configured Python because it is not from required environment '",
        env_name,
        "': ",
        python_path,
        ". Using: ",
        env_python,
        call. = FALSE
      )
    }
    return(env_python)
  }

  assert_required_python_env(configured_python_bin, env_name)
}

quote_args <- function(args) {
  paste(vapply(args, shQuote, character(1)), collapse = " ")
}

write_scvelo_runner <- function(
  run_script,
  python_bin,
  python_script,
  args,
  log_file
) {
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste(
      shQuote(python_bin),
      shQuote(python_script),
      quote_args(args),
      ">",
      shQuote(log_file),
      "2>&1"
    )
  )
  writeLines(lines, con = run_script)
  Sys.chmod(run_script, mode = "0755")
}

write_paga_runner <- function(
  run_script,
  python_bin,
  python_script,
  args,
  log_file
) {
  lines <- c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste(
      shQuote(python_bin),
      shQuote(python_script),
      quote_args(args),
      ">",
      shQuote(log_file),
      "2>&1"
    )
  )
  writeLines(lines, con = run_script)
  Sys.chmod(run_script, mode = "0755")
}

is_nonempty_assay_matrix <- function(mat) {
  !is.null(mat) && length(dim(mat)) == 2L && nrow(mat) > 0 && ncol(mat) > 0
}

get_paga_expression_matrix <- function(obj, cells, assay = "RNA", slot_name = "data") {
  if (!(assay %in% names(obj@assays))) {
    stop("PAGA expression assay not found in Seurat object: ", assay, call. = FALSE)
  }

  slot_candidates <- unique(c(slot_name, if (!identical(slot_name, "counts")) "counts" else "data"))
  for (slot_candidate in slot_candidates) {
    mat <- get_assay_matrix(obj, assay = assay, slot_name = slot_candidate)
    if (!is_nonempty_assay_matrix(mat)) next

    missing_cells <- setdiff(cells, colnames(mat))
    if (length(missing_cells) > 0) {
      stop(
        "PAGA expression matrix is missing ",
        length(missing_cells),
        " cell(s), e.g. ",
        paste(utils::head(missing_cells, 5), collapse = ", "),
        call. = FALSE
      )
    }

    mat <- mat[, cells, drop = FALSE]
    return(list(matrix = mat, assay = assay, slot = slot_candidate))
  }

  stop(
    "Cannot find a non-empty PAGA expression matrix for assay ",
    assay,
    " in slot/layer candidates: ",
    paste(slot_candidates, collapse = ", "),
    call. = FALSE
  )
}

write_paga_expression_export <- function(
  obj,
  cells,
  output_dir,
  assay = "RNA",
  slot_name = "data"
) {
  cells <- as.character(cells)
  cells <- cells[!is.na(cells) & nzchar(cells)]
  if (length(cells) == 0) {
    return(list(args = character(0), matrix_file = "", genes_file = "", cells_file = "", source = ""))
  }

  expr_info <- get_paga_expression_matrix(obj, cells = cells, assay = assay, slot_name = slot_name)
  expr_mat <- expr_info$matrix
  if (!inherits(expr_mat, "sparseMatrix")) {
    expr_mat <- Matrix::Matrix(expr_mat, sparse = TRUE)
  }

  matrix_file <- file.path(output_dir, "expression_matrix.mtx")
  genes_file <- file.path(output_dir, "expression_genes.csv")
  cells_file <- file.path(output_dir, "expression_cells.csv")
  Matrix::writeMM(expr_mat, matrix_file)
  write_table_csv(data.frame(gene = rownames(expr_mat), stringsAsFactors = FALSE), genes_file)
  write_table_csv(data.frame(cell = colnames(expr_mat), stringsAsFactors = FALSE), cells_file)

  source_label <- paste(expr_info$assay, expr_info$slot, sep = ":")
  list(
    args = c(
      "--expression-matrix", matrix_file,
      "--expression-genes", genes_file,
      "--expression-cells", cells_file,
      "--expression-source", source_label
    ),
    matrix_file = matrix_file,
    genes_file = genes_file,
    cells_file = cells_file,
    source = source_label
  )
}

plot_metric_umap <- function(df, metric, file_stub, title) {
  if (!(metric %in% colnames(df))) return(invisible(NULL))
  values <- suppressWarnings(as.numeric(df[[metric]]))
  if (!any(is.finite(values))) return(invisible(NULL))

  plot_df <- df
  plot_df[[metric]] <- values
  facet_col <- if ("trajectory_analysis_group" %in% colnames(plot_df)) {
    "trajectory_analysis_group"
  } else if ("trajectory_group" %in% colnames(plot_df)) {
    "trajectory_group"
  } else {
    "Dose"
  }
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[metric]])) +
    geom_point(size = 0.25, alpha = 0.85, stroke = 0) +
    coord_equal() +
    facet_wrap(stats::as.formula(paste("~", facet_col))) +
    scale_color_gradientn(colors = c("#17324D", "#2D8C88", "#F5C542", "#C43B3B"), na.value = "grey88") +
    labs(title = title, color = metric) +
    make_pretty_umap_theme()
  save_both(p, file_stub, width = 11, height = 7)
  invisible(p)
}

prepare_or_run_paga_analysis <- function(
  group_name,
  group_safe,
  analysis_type,
  analysis_root,
  metadata_df,
  seurat_obj = NULL,
  export_expression = TRUE,
  expression_assay = "RNA",
  expression_slot = "data",
  color_cols,
  dose = NA_character_,
  ploidy = NA_character_,
  trajectory_context = NA_character_,
  base_groups = character(0),
  tn_scope = NA_character_,
  ploidy_scope = "ploidy_all"
) {
  analysis_dir <- .ensure_dir(file.path(analysis_root, group_safe))
  metadata_df <- metadata_df %>%
    dplyr::mutate(
      trajectory_analysis_group = group_name,
      trajectory_tn_scope = tn_scope,
      trajectory_ploidy_scope = ploidy_scope,
      trajectory_branch = trajectory_branch
    ) %>%
    dplyr::arrange(.data$cell)

  metadata_file <- file.path(analysis_dir, "cells_metadata_umap.csv")
  cells_file <- file.path(analysis_dir, "cells.txt")
  write_table_csv(metadata_df, metadata_file)
  writeLines(metadata_df$cell, con = cells_file)

  expression_export <- list(
    args = character(0),
    matrix_file = "",
    genes_file = "",
    cells_file = "",
    source = ""
  )
  if (isTRUE(export_expression)) {
    if (is.null(seurat_obj)) {
      stop("PAGA expression export requested but seurat_obj is NULL.", call. = FALSE)
    }
    expression_export <- write_paga_expression_export(
      obj = seurat_obj,
      cells = metadata_df$cell,
      output_dir = analysis_dir,
      assay = expression_assay,
      slot_name = expression_slot
    )
  }

  log_file <- file.path(analysis_dir, "paga.log")
  run_script <- file.path(analysis_dir, "run_paga.sh")
  paga_args <- c(
    "--metadata", metadata_file,
    "--output-dir", analysis_dir,
    "--analysis-label", group_name,
    "--basis", "umap",
    "--cluster-col", cluster_col,
    "--color-cols", paste(color_cols, collapse = ","),
    "--pca-prefix", paga_pca_prefix,
    "--n-pcs", as.character(paga_n_pcs),
    "--n-neighbors", as.character(paga_n_neighbors),
    "--min-cells", as.character(paga_min_cells),
    "--paga-threshold", as.character(paga_threshold)
  )
  paga_args <- c(paga_args, expression_export$args)
  if (nzchar(paga_root_clusters)) {
    paga_args <- c(paga_args, "--root-clusters", paga_root_clusters)
  }

  write_paga_runner(
    run_script = run_script,
    python_bin = python_bin,
    python_script = paga_python_script,
    args = paga_args,
    log_file = log_file
  )

  status <- "prepared"
  exit_status <- NA_integer_
  missing_required_cols <- setdiff(c("cell", "UMAP_1", "UMAP_2", cluster_col), colnames(metadata_df))
  if (length(missing_required_cols) > 0) {
    status <- paste0("missing_metadata_columns:", paste(missing_required_cols, collapse = ";"))
  } else if (nrow(metadata_df) < paga_min_cells) {
    status <- "skipped_too_few_cells"
  } else if (run_paga) {
    message("Running PAGA for ", group_name)
    if (!nzchar(Sys.which(python_bin)) && !file.exists(python_bin)) {
      stop("Python binary not found: ", python_bin, call. = FALSE)
    }
    exit_status <- system2(
      command = python_bin,
      args = c(paga_python_script, paga_args),
      stdout = log_file,
      stderr = log_file
    )
    status <- if (identical(exit_status, 0L)) "completed" else "failed"
    if (!identical(exit_status, 0L)) {
      stop("PAGA failed for ", group_name, ". See log: ", log_file, call. = FALSE)
    }
  } else {
    status <- "prepared_only"
  }

  data.frame(
    analysis_group = group_name,
    analysis_type = analysis_type,
    trajectory_group = group_name,
    trajectory_branch = trajectory_branch,
    tn_scope = as.character(tn_scope),
    ploidy_scope = as.character(ploidy_scope),
    base_groups = format_input_files(base_groups),
    Dose = as.character(dose),
    ploidy = as.character(ploidy),
    trajectory_context = as.character(trajectory_context),
    paga_dir = analysis_dir,
    n_cells = nrow(metadata_df),
    cluster_col = cluster_col,
    run_paga = run_paga,
    status = status,
    exit_status = exit_status,
    n_pcs = paga_n_pcs,
    n_neighbors = paga_n_neighbors,
    pca_prefix = paga_pca_prefix,
    paga_threshold = paga_threshold,
    root_clusters = paga_root_clusters,
    expression_export = export_expression,
    expression_source = expression_export$source,
    expression_matrix = expression_export$matrix_file,
    expression_genes = expression_export$genes_file,
    expression_cells = expression_export$cells_file,
    run_script = run_script,
    log_file = log_file,
    stringsAsFactors = FALSE
  )
}

prepare_or_run_scvelo_analysis <- function(
  group_name,
  group_safe,
  analysis_type,
  analysis_root,
  metadata_df,
  input_files,
  color_cols,
  dose = NA_character_,
  ploidy = NA_character_,
  trajectory_context = NA_character_,
  base_groups = character(0),
  tn_scope = NA_character_,
  ploidy_scope = "ploidy_all",
  shape_col = "",
  shape_order = character(0),
  shape_markers = character(0)
) {
  analysis_dir <- .ensure_dir(file.path(analysis_root, group_safe))
  metadata_df <- metadata_df %>%
    dplyr::mutate(
      trajectory_analysis_group = group_name,
      trajectory_tn_scope = tn_scope,
      trajectory_ploidy_scope = ploidy_scope,
      trajectory_branch = trajectory_branch
    ) %>%
    dplyr::arrange(.data$cell)

  metadata_file <- file.path(analysis_dir, "cells_metadata_umap.csv")
  cells_file <- file.path(analysis_dir, "cells.txt")
  write_table_csv(metadata_df, metadata_file)
  writeLines(metadata_df$cell, con = cells_file)

  input_files <- unique(as.character(input_files))
  input_files <- input_files[!is.na(input_files) & nzchar(input_files)]
  input_exists <- all_input_files_exist(input_files)
  input_file_arg <- if (input_exists) {
    make_scvelo_input_arg(input_files)
  } else {
    file.path(analysis_dir, "MISSING_velocity_input.loom_or_h5ad")
  }
  log_file <- file.path(analysis_dir, "scvelo.log")
  run_script <- file.path(analysis_dir, "run_scvelo.sh")

  scvelo_args <- c(
    "--input", input_file_arg,
    "--metadata", metadata_file,
    "--output-dir", analysis_dir,
    "--dose-label", group_name,
    "--basis", "umap",
    "--cluster-col", cluster_col,
    "--color-cols", paste(color_cols, collapse = ","),
    "--mode", scvelo_mode,
    "--min-shared-counts", as.character(scvelo_min_shared_counts),
    "--n-top-genes", as.character(scvelo_n_top_genes),
    "--n-pcs", as.character(scvelo_n_pcs),
    "--n-neighbors", as.character(scvelo_n_neighbors),
    "--n-jobs", as.character(scvelo_n_jobs),
    "--min-matched-cells", as.character(scvelo_min_matched_cells),
    "--stream-density", as.character(scvelo_stream_density)
  )
  if (isTRUE(scvelo_use_metadata_pca)) {
    scvelo_args <- c(scvelo_args, "--use-metadata-pca", "--pca-prefix", scvelo_pca_prefix)
  }
  if (nzchar(scvelo_root_clusters)) {
    scvelo_args <- c(scvelo_args, "--root-clusters", scvelo_root_clusters)
  }
  if (nzchar(scvelo_end_clusters)) {
    scvelo_args <- c(scvelo_args, "--end-clusters", scvelo_end_clusters)
  }
  if (nzchar(shape_col)) {
    scvelo_args <- c(scvelo_args, "--shape-col", shape_col)
    if (length(shape_order) > 0) {
      scvelo_args <- c(scvelo_args, "--shape-order", paste(shape_order, collapse = ","))
    }
    if (length(shape_markers) > 0) {
      scvelo_args <- c(scvelo_args, "--shape-markers", paste(shape_markers, collapse = ","))
    }
  }
  if (!is.null(cell_map_file) && nzchar(cell_map_file) && file.exists(cell_map_file)) {
    scvelo_args <- c(scvelo_args, "--cell-map", normalizePath(cell_map_file, mustWork = TRUE))
  }

  write_scvelo_runner(
    run_script = run_script,
    python_bin = python_bin,
    python_script = python_script,
    args = scvelo_args,
    log_file = log_file
  )

  status <- "prepared"
  exit_status <- NA_integer_
  if (nrow(metadata_df) == 0) {
    status <- "missing_cells"
  } else if (run_scvelo && input_exists) {
    message("Running scVelo for ", group_name)
    if (!nzchar(Sys.which(python_bin)) && !file.exists(python_bin)) {
      stop("Python binary not found: ", python_bin, call. = FALSE)
    }
    exit_status <- system2(
      command = python_bin,
      args = c(python_script, scvelo_args),
      stdout = log_file,
      stderr = log_file
    )
    status <- if (identical(exit_status, 0L)) "completed" else "failed"
    if (!identical(exit_status, 0L)) {
      stop("scVelo failed for ", group_name, ". See log: ", log_file, call. = FALSE)
    }
  } else if (run_scvelo && !input_exists) {
    status <- "missing_velocity_input"
  } else {
    status <- "prepared_only"
  }

  data.frame(
    analysis_group = group_name,
    analysis_type = analysis_type,
    trajectory_group = group_name,
    trajectory_branch = trajectory_branch,
    tn_scope = as.character(tn_scope),
    ploidy_scope = as.character(ploidy_scope),
    base_groups = format_input_files(base_groups),
    Dose = as.character(dose),
    ploidy = as.character(ploidy),
    trajectory_context = as.character(trajectory_context),
    dose_dir = analysis_dir,
    n_cells = nrow(metadata_df),
    velocity_input_file = format_input_files(input_files),
    n_velocity_input_files = length(input_files),
    input_exists = input_exists,
    run_scvelo = run_scvelo,
    status = status,
    exit_status = exit_status,
    use_metadata_pca = scvelo_use_metadata_pca,
    pca_prefix = scvelo_pca_prefix,
    root_clusters = scvelo_root_clusters,
    end_clusters = scvelo_end_clusters,
    neighbor_representation = ifelse(scvelo_use_metadata_pca, "Seurat_PCA", "scVelo_auto"),
    shape_col = shape_col,
    shape_order = paste(shape_order, collapse = ";"),
    run_script = run_script,
    log_file = log_file,
    stringsAsFactors = FALSE
  )
}

is_missing_task_value <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  x <- as.character(x)[1]
  is.na(x) || !nzchar(trimws(x)) || x %in% c("NA", "NaN", "NULL")
}

split_task_values <- function(x) {
  if (is_missing_task_value(x)) return(character(0))
  vals <- trimws(unlist(strsplit(as.character(x)[1], ",", fixed = TRUE)))
  vals <- vals[!is.na(vals) & nzchar(vals) & !(vals %in% c("NA", "NaN", "NULL"))]
  vals
}

comparison_ident_label <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(NA_character_)
  paste(values, collapse = "+")
}

comparison_rest_label <- function(task_row) {
  ident_2_values <- split_task_values(task_row$ident_2)
  if (length(ident_2_values) == 0 || identical(tolower(ident_2_values[1]), "rest")) {
    return("rest")
  }
  if (length(ident_2_values) > 1 || grepl("_vs_rest$", as.character(task_row$comparison), ignore.case = TRUE)) {
    return("rest")
  }
  ident_2_values[1]
}

comparison_task_stub <- function(task_row) {
  pieces <- c(
    sprintf("task_%03d", as.integer(task_row$trajectory_task_id)),
    as.character(task_row$scope),
    if (!is_missing_task_value(task_row$cluster)) paste0("cluster_", task_row$cluster) else NA_character_,
    if (!is_missing_task_value(task_row$subset_group)) as.character(task_row$subset_group) else NA_character_,
    if (!is_missing_task_value(task_row$dose)) as.character(task_row$dose) else NA_character_,
    as.character(task_row$comparison)
  )
  sanitize_path_component(paste(pieces[!is.na(pieces) & nzchar(pieces)], collapse = "__"))[1]
}

make_velocity_comparison_skip <- function(
  task_row,
  status,
  message,
  n_ident_1 = 0L,
  n_ident_2 = 0L,
  n_cells = 0L,
  cells = character(0),
  group_by_cell = character(0),
  original_group_by_cell = character(0),
  selected_sample_folders = character(0),
  missing_sample_folders = character(0),
  available_input_files = character(0),
  input_files = character(0)
) {
  ident_1_values <- split_task_values(task_row$ident_1)
  list(
    status = status,
    message = message,
    task = task_row,
    cells = cells,
    group_by_cell = group_by_cell,
    original_group_by_cell = original_group_by_cell,
    ident_1_label = comparison_ident_label(ident_1_values),
    ident_2_label = comparison_rest_label(task_row),
    n_ident_1 = as.integer(n_ident_1),
    n_ident_2 = as.integer(n_ident_2),
    n_cells = as.integer(n_cells),
    selected_sample_folders = selected_sample_folders,
    missing_sample_folders = missing_sample_folders,
    available_input_files = available_input_files,
    input_files = input_files,
    comparison_stub = comparison_task_stub(task_row)
  )
}

build_velocity_comparison_task <- function(
  meta_all,
  task_row,
  cluster_match_col,
  sample_folder_col,
  sample_manifest,
  min_cells_per_group,
  min_total_cells
) {
  meta <- as.data.frame(meta_all, stringsAsFactors = FALSE)
  if (!("cell" %in% colnames(meta))) {
    return(make_velocity_comparison_skip(task_row, "skipped", "Input metadata is missing cell column."))
  }
  rownames(meta) <- meta$cell

  group_col <- as.character(task_row$group_col)[1]
  if (!(group_col %in% colnames(meta))) {
    return(make_velocity_comparison_skip(task_row, "skipped", paste0("Missing group_col in metadata: ", group_col)))
  }

  base_cells <- rownames(meta)
  scope <- as.character(task_row$scope)[1]

  if (grepl("^within_cluster_", scope) && !is_missing_task_value(task_row$cluster)) {
    if (!(cluster_match_col %in% colnames(meta))) {
      return(make_velocity_comparison_skip(task_row, "skipped", paste0("Missing cluster match column in metadata: ", cluster_match_col)))
    }
    cluster_values <- as.character(meta[[cluster_match_col]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(cluster_values) & cluster_values == as.character(task_row$cluster)])
  }

  if (!is_missing_task_value(task_row$subset_group)) {
    if (!("TN" %in% colnames(meta))) {
      return(make_velocity_comparison_skip(task_row, "skipped", "Missing TN metadata for subset_group filter."))
    }
    tn_values <- as.character(meta[["TN"]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(tn_values) & tn_values == as.character(task_row$subset_group)])
  }

  if (!is_missing_task_value(task_row$dose) && !identical(group_col, "Dose_DEG")) {
    if (!("Dose_DEG" %in% colnames(meta))) {
      return(make_velocity_comparison_skip(task_row, "skipped", "Missing Dose_DEG metadata for dose filter."))
    }
    dose_values <- as.character(meta[["Dose_DEG"]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(dose_values) & dose_values == as.character(task_row$dose)])
  }

  base_cells <- intersect(base_cells, rownames(meta))
  if (length(base_cells) == 0) {
    return(make_velocity_comparison_skip(task_row, "skipped", "No cells in comparison universe after scope filters."))
  }

  group_values <- as.character(meta[base_cells, group_col, drop = TRUE])
  valid_group <- !is.na(group_values) & nzchar(group_values)
  base_cells <- base_cells[valid_group]
  group_values <- group_values[valid_group]
  if (length(base_cells) == 0) {
    return(make_velocity_comparison_skip(task_row, "skipped", paste0("No non-missing values in ", group_col, ".")))
  }

  ident_1_values <- split_task_values(task_row$ident_1)
  if (length(ident_1_values) == 0) {
    return(make_velocity_comparison_skip(task_row, "skipped", "Missing ident_1."))
  }
  ident_2_values <- split_task_values(task_row$ident_2)
  ident_2_is_rest <- length(ident_2_values) == 0 || identical(tolower(ident_2_values[1]), "rest")

  cells_1 <- base_cells[group_values %in% ident_1_values]
  cells_2 <- if (ident_2_is_rest) {
    base_cells[!(group_values %in% ident_1_values)]
  } else {
    base_cells[group_values %in% ident_2_values]
  }

  n_ident_1 <- length(cells_1)
  n_ident_2 <- length(cells_2)
  selected_cells <- unique(c(cells_1, cells_2))
  n_cells <- length(selected_cells)

  if (n_ident_1 < min_cells_per_group || n_ident_2 < min_cells_per_group) {
    return(make_velocity_comparison_skip(
      task_row,
      "skipped",
      paste0("Fewer than ", min_cells_per_group, " cells in one or both trajectory comparison groups."),
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells
    ))
  }
  if (n_cells < min_total_cells) {
    return(make_velocity_comparison_skip(
      task_row,
      "skipped",
      paste0("Fewer than ", min_total_cells, " total cells in trajectory comparison universe."),
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells
    ))
  }

  ident_1_label <- comparison_ident_label(ident_1_values)
  ident_2_label <- comparison_rest_label(task_row)
  group_by_cell <- stats::setNames(rep(NA_character_, n_cells), selected_cells)
  group_by_cell[cells_1] <- ident_1_label
  group_by_cell[cells_2] <- ident_2_label
  original_group_by_cell <- stats::setNames(as.character(meta[selected_cells, group_col, drop = TRUE]), selected_cells)

  if (!(sample_folder_col %in% colnames(meta))) {
    return(make_velocity_comparison_skip(
      task_row,
      "missing_velocity_input",
      paste0("Missing sample folder metadata column for loom lookup: ", sample_folder_col),
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells,
      cells = selected_cells,
      group_by_cell = group_by_cell,
      original_group_by_cell = original_group_by_cell
    ))
  }

  selected_samples <- sort(unique(as.character(meta[selected_cells, sample_folder_col, drop = TRUE])))
  selected_samples <- selected_samples[!is.na(selected_samples) & nzchar(selected_samples)]
  if (length(selected_samples) == 0) {
    return(make_velocity_comparison_skip(
      task_row,
      "missing_velocity_input",
      paste0("No sample folder values found in selected cells: ", sample_folder_col),
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells,
      cells = selected_cells,
      group_by_cell = group_by_cell,
      original_group_by_cell = original_group_by_cell
    ))
  }
  if (is.null(sample_manifest) || nrow(sample_manifest) == 0) {
    return(make_velocity_comparison_skip(
      task_row,
      "missing_velocity_input",
      "No sample-level loom manifest was available.",
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells,
      cells = selected_cells,
      group_by_cell = group_by_cell,
      original_group_by_cell = original_group_by_cell,
      selected_sample_folders = selected_samples
    ))
  }

  manifest_rows <- sample_manifest[sample_manifest$sample_folder %in% selected_samples, , drop = FALSE]
  missing_manifest_samples <- setdiff(selected_samples, manifest_rows$sample_folder)
  missing_loom_idx <- is.na(manifest_rows$velocyto_loom_exists) | !manifest_rows$velocyto_loom_exists
  missing_loom_samples <- manifest_rows$sample_folder[missing_loom_idx]
  missing_samples <- sort(unique(c(missing_manifest_samples, missing_loom_samples)))
  available_input_files <- unique(as.character(manifest_rows$velocyto_loom_file[manifest_rows$velocyto_loom_exists]))
  available_input_files <- available_input_files[!is.na(available_input_files) & nzchar(available_input_files)]
  input_files <- if (length(missing_samples) == 0) available_input_files else character(0)
  if (length(input_files) == 0 && length(missing_samples) == 0) {
    missing_samples <- selected_samples
  }
  input_status <- if (length(missing_samples) == 0 && length(input_files) > 0) "ready" else "missing_velocity_input"
  input_message <- if (identical(input_status, "ready")) {
    NA_character_
  } else {
    paste0("Missing existing loom for sample folder(s): ", paste(missing_samples, collapse = ";"))
  }

  make_velocity_comparison_skip(
    task_row,
    input_status,
    input_message,
    n_ident_1 = n_ident_1,
    n_ident_2 = n_ident_2,
    n_cells = n_cells,
    cells = selected_cells,
    group_by_cell = group_by_cell,
    original_group_by_cell = original_group_by_cell,
    selected_sample_folders = selected_samples,
    missing_sample_folders = missing_samples,
    available_input_files = available_input_files,
    input_files = input_files
  )
}

attach_velocity_comparison_metadata <- function(df, task_info) {
  task <- task_info$task
  df$TrajectoryComparisonID <- rep(as.character(task$trajectory_task_id), nrow(df))
  df$TrajectoryScope <- rep(as.character(task$scope), nrow(df))
  df$TrajectoryComparison <- rep(as.character(task$comparison), nrow(df))
  df$TrajectoryIdent1 <- rep(as.character(task$ident_1), nrow(df))
  df$TrajectoryIdent2 <- rep(as.character(task$ident_2), nrow(df))
  df$TrajectoryGroupCol <- rep(as.character(task$group_col), nrow(df))
  df$TrajectoryCluster <- rep(as.character(task$cluster), nrow(df))
  df$TrajectorySubsetGroup <- rep(as.character(task$subset_group), nrow(df))
  df$TrajectoryDose <- rep(as.character(task$dose), nrow(df))
  df
}

make_velocity_task_summary_row <- function(task_info, status = task_info$status, output_dir = NA_character_) {
  task <- task_info$task
  data.frame(
    trajectory_task_id = as.integer(task$trajectory_task_id),
    scope = as.character(task$scope),
    cluster = as.character(task$cluster),
    subset_group = as.character(task$subset_group),
    dose = as.character(task$dose),
    comparison = as.character(task$comparison),
    group_col = as.character(task$group_col),
    ident_1 = as.character(task$ident_1),
    ident_2 = as.character(task$ident_2),
    n_ident_1 = as.integer(task_info$n_ident_1),
    n_ident_2 = as.integer(task_info$n_ident_2),
    n_cells = as.integer(task_info$n_cells),
    n_selected_sample_folders = length(task_info$selected_sample_folders),
    selected_sample_folders = paste(task_info$selected_sample_folders, collapse = ";"),
    missing_sample_folders = paste(task_info$missing_sample_folders, collapse = ";"),
    available_velocity_input_file = format_input_files(task_info$available_input_files),
    velocity_input_file = format_input_files(task_info$input_files),
    status = status,
    output_dir = output_dir,
    message = as.character(task_info$message),
    stringsAsFactors = FALSE
  )
}

prepare_or_run_deg_comparison_scvelo <- function(
  meta_all,
  comparison_summary_file,
  output_dir,
  summary_dir,
  cluster_match_col,
  sample_folder_col,
  sample_manifest,
  base_color_cols,
  min_cells_per_group,
  min_total_cells
) {
  if (!file.exists(comparison_summary_file)) {
    warning("Skipping 03a DEG-comparison velocity analyses because comparison summary is missing: ", comparison_summary_file)
    return(list(run_rows = list(), task_manifest = data.frame()))
  }

  comparison_summary <- readr::read_csv(comparison_summary_file, show_col_types = FALSE)
  comparison_summary <- as.data.frame(comparison_summary, stringsAsFactors = FALSE)
  required_cols <- c("scope", "comparison", "group_col", "ident_1", "ident_2", "status")
  missing_cols <- setdiff(required_cols, colnames(comparison_summary))
  if (length(missing_cols) > 0) {
    stop("DEG comparison summary is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  for (optional_col in c("cluster", "subset_group", "dose", "n_ident_1", "n_ident_2")) {
    if (!(optional_col %in% colnames(comparison_summary))) comparison_summary[[optional_col]] <- NA
  }

  comparison_tasks <- comparison_summary %>%
    dplyr::filter(.data$status %in% c("success", "cached")) %>%
    dplyr::mutate(trajectory_task_id = dplyr::row_number()) %>%
    as.data.frame(stringsAsFactors = FALSE)
  if (nrow(comparison_tasks) == 0) {
    warning("No successful or cached DEG comparisons found for velocity analysis.")
    return(list(run_rows = list(), task_manifest = data.frame()))
  }

  task_infos <- lapply(seq_len(nrow(comparison_tasks)), function(i) {
    build_velocity_comparison_task(
      meta_all = meta_all,
      task_row = comparison_tasks[i, , drop = FALSE],
      cluster_match_col = cluster_match_col,
      sample_folder_col = sample_folder_col,
      sample_manifest = sample_manifest,
      min_cells_per_group = min_cells_per_group,
      min_total_cells = min_total_cells
    )
  })

  run_rows <- list()
  task_summary_rows <- list()

  for (task_info in task_infos) {
    task <- task_info$task
    group_safe <- file.path(sanitize_path_component(as.character(task$scope)), task_info$comparison_stub)
    analysis_dir <- file.path(output_dir, group_safe)
    .ensure_dir(analysis_dir)

    if (length(task_info$cells) > 0 && task_info$status %in% c("ready", "missing_velocity_input")) {
      comparison_meta <- meta_all[match(task_info$cells, meta_all$cell), , drop = FALSE]
      comparison_meta <- comparison_meta[!is.na(comparison_meta$cell), , drop = FALSE]
      comparison_meta$TrajectoryComparisonGroup <- unname(task_info$group_by_cell[comparison_meta$cell])
      comparison_meta$TrajectoryOriginalGroup <- unname(task_info$original_group_by_cell[comparison_meta$cell])
      comparison_meta$TrajectoryComparisonGroup <- factor(
        as.character(comparison_meta$TrajectoryComparisonGroup),
        levels = c(task_info$ident_1_label, task_info$ident_2_label)
      )
      comparison_meta <- attach_velocity_comparison_metadata(comparison_meta, task_info)
      comparison_meta$trajectory_shape_group <- comparison_meta$TrajectoryComparisonGroup

      comparison_color_cols <- unique(c(
        base_color_cols,
        "TrajectoryComparisonGroup",
        "TrajectoryOriginalGroup",
        "TrajectoryComparison",
        "TrajectoryScope",
        "TrajectoryGroupCol"
      ))
      run_row <- prepare_or_run_scvelo_analysis(
        group_name = paste0("task_", task$trajectory_task_id, " ", task$comparison),
        group_safe = group_safe,
        analysis_type = "deg_comparison",
        analysis_root = output_dir,
        metadata_df = comparison_meta,
        input_files = task_info$input_files,
        color_cols = comparison_color_cols,
        dose = if (!is_missing_task_value(task$dose)) as.character(task$dose) else NA_character_,
        ploidy = NA_character_,
        trajectory_context = if (!is_missing_task_value(task$subset_group)) as.character(task$subset_group) else NA_character_,
        base_groups = c(task_info$ident_1_label, task_info$ident_2_label),
        shape_col = "TrajectoryComparisonGroup",
        shape_order = c(task_info$ident_1_label, task_info$ident_2_label),
        shape_markers = c("o", "^")
      )
      if (identical(task_info$status, "missing_velocity_input")) {
        run_row$status <- "missing_velocity_input"
      }
      run_row <- cbind(
        run_row,
        make_velocity_task_summary_row(task_info, status = task_info$status, output_dir = analysis_dir)[
          ,
          c(
            "trajectory_task_id",
            "scope",
            "cluster",
            "subset_group",
            "comparison",
            "group_col",
            "ident_1",
            "ident_2",
            "n_ident_1",
            "n_ident_2",
            "n_selected_sample_folders",
            "selected_sample_folders",
            "missing_sample_folders",
            "available_velocity_input_file",
            "message"
          ),
          drop = FALSE
        ]
      )
      run_rows[[task_info$comparison_stub]] <- run_row
    }

    task_summary_rows[[task_info$comparison_stub]] <- make_velocity_task_summary_row(
      task_info,
      status = task_info$status,
      output_dir = analysis_dir
    )
  }

  task_manifest <- dplyr::bind_rows(task_summary_rows)
  run_summary <- dplyr::bind_rows(run_rows)
  write_table_csv(task_manifest, file.path(summary_dir, "trajectory_DEG_comparison_velocity_task_manifest.csv"))
  write_table_csv(task_manifest, file.path(summary_dir, "resolved_velocity_inputs_by_DEG_comparison.csv"))
  write_table_csv(run_summary, file.path(summary_dir, "scvelo_run_summary_DEG_comparisons.csv"))

  list(run_rows = run_rows, task_manifest = task_manifest)
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
output_root <- file.path(results_root, "04_trajectory")
default_velocyto_output_root <- resolve_first_existing(c(
  file.path(getwd(), "Results", "04_trajectory", "00_inputs", "velocyto_loom"),
  file.path(output_root, "00_inputs", "velocyto_loom")
))
velocyto_output_root <- cfg_value(
  config,
  "Velocyto_output_root",
  "VELOCYTO_OUTPUT_ROOT",
  default = default_velocyto_output_root
)
velocyto_work_root <- cfg_value(
  config,
  "Velocyto_work_root",
  "VELOCYTO_WORK_ROOT",
  default = velocyto_output_root
)
velocyto_run_tag <- cfg_value(
  config,
  "Velocyto_run_tag",
  "VELOCYTO_RUN_TAG",
  default = paste0(
    "run_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_pid",
    Sys.getpid(),
    "_",
    sample.int(.Machine$integer.max, 1)
  )
)
velocity_input_root <- cfg_value(
  config,
  "Velocity_input_root",
  "VELOCITY_INPUT_ROOT",
  default = file.path(results_root, "00_velocity_inputs")
)
velocity_manifest <- cfg_value(config, "Velocity_input_manifest", "VELOCITY_INPUT_MANIFEST", default = "")
cell_map_file <- cfg_value(config, "Velocity_cell_map", "VELOCITY_CELL_MAP", default = "")
cellranger_root <- cfg_value(
  config,
  "Cellranger_root",
  "CELLRANGER_ROOT",
  default = "/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/SUM-159/A02_cellRanger"
)
cellranger_sample_folder_col <- cfg_value(
  config,
  "Cellranger_sample_folder_col",
  "CELLRANGER_SAMPLE_FOLDER_COL",
  default = "sample_folder"
)
velocyto_gtf <- cfg_value(
  config,
  "Velocyto_gtf",
  "VELOCYTO_GTF",
  default = "/share/lab_crd/lab_crd/taoli/Project/BreastCancerOrthotopicModels/gene/refdata-gex-GRCh38_and_GRCm39-2024-A/genes/genes.gtf"
)
velocyto_gtf <- prefer_plain_gtf(velocyto_gtf)
run_velocyto <- parse_bool(cfg_value(config, "Run_velocyto", "RUN_VELOCYTO", default = "FALSE"), default = FALSE)
required_python_env <- cfg_value(config, "Scvelo_python_env", "SCVELO_PYTHON_ENV", default = "rna_velocity_py310")
default_python_bin <- resolve_python_from_env_name(required_python_env)
if (is.na(default_python_bin)) {
  default_python_bin <- file.path(Sys.getenv("HOME", unset = "~"), "venvs", required_python_env, "bin", "python")
}
python_bin <- cfg_value(config, "Scvelo_python", "SCVELO_PYTHON", default = default_python_bin)
python_bin <- select_required_python_env(python_bin, required_python_env)
default_velocyto_bin <- file.path(dirname(python_bin), "velocyto")
if (!file.exists(default_velocyto_bin)) {
  default_velocyto_bin <- "velocyto"
}
velocyto_bin <- cfg_value(config, "Velocyto_bin", "VELOCYTO_BIN", default = default_velocyto_bin)
velocyto_bin <- resolve_executable_path(velocyto_bin)
default_samtools_candidates <- c(
  file.path(dirname(python_bin), "samtools"),
  path.expand("~/samtools-1.23.1/install/bin/samtools"),
  "/home/4482173/samtools-1.23.1/install/bin/samtools",
  "samtools"
)
default_samtools_hit <- default_samtools_candidates[file.exists(default_samtools_candidates)]
default_samtools_bin <- if (length(default_samtools_hit) > 0) default_samtools_hit[1] else "samtools"
samtools_bin <- cfg_value(config, "Samtools_bin", "SAMTOOLS_BIN", default = default_samtools_bin)
samtools_bin <- resolve_executable_path(samtools_bin)
velocyto_samtools_threads <- parse_int_cfg(
  cfg_value(config, "Velocyto_samtools_threads", "VELOCYTO_SAMTOOLS_THREADS", default = "62"),
  62L,
  "Velocyto_samtools_threads"
)
python_script <- file.path(script_dir, "04_trajectory_scvelo.py")
run_scvelo <- parse_bool(cfg_value(config, "Run_scvelo", "RUN_SCVELO", default = "TRUE"), default = TRUE)

scvelo_mode <- cfg_value(config, "Scvelo_mode", "SCVELO_MODE", default = "stochastic")
if (!(scvelo_mode %in% c("deterministic", "stochastic", "dynamical"))) {
  stop("Scvelo_mode must be one of deterministic, stochastic, dynamical.", call. = FALSE)
}
scvelo_n_jobs <- parse_int_cfg(cfg_value(config, "Scvelo_n_jobs", "SCVELO_N_JOBS", default = "4"), 4L, "Scvelo_n_jobs")
scvelo_min_shared_counts <- parse_int_cfg(cfg_value(config, "Scvelo_min_shared_counts", "SCVELO_MIN_SHARED_COUNTS", default = "20"), 20L, "Scvelo_min_shared_counts")
scvelo_n_top_genes <- parse_int_cfg(cfg_value(config, "Scvelo_n_top_genes", "SCVELO_N_TOP_GENES", default = "2000"), 2000L, "Scvelo_n_top_genes")
scvelo_n_pcs <- parse_int_cfg(cfg_value(config, "Scvelo_n_pcs", "SCVELO_N_PCS", default = "30"), 30L, "Scvelo_n_pcs")
scvelo_n_neighbors <- parse_int_cfg(cfg_value(config, "Scvelo_n_neighbors", "SCVELO_N_NEIGHBORS", default = "30"), 30L, "Scvelo_n_neighbors")
scvelo_min_matched_cells <- parse_int_cfg(cfg_value(config, "Scvelo_min_matched_cells", "SCVELO_MIN_MATCHED_CELLS", default = "50"), 50L, "Scvelo_min_matched_cells")
scvelo_stream_density <- parse_numeric_cfg(cfg_value(config, "Scvelo_stream_density", "SCVELO_STREAM_DENSITY", default = "1.2"), 1.2, "Scvelo_stream_density")
scvelo_use_metadata_pca <- TRUE
scvelo_pca_prefix <- "PCA_"
scvelo_root_clusters <- "14"
scvelo_end_clusters <- "13"
trajectory_branch <- "root14_end13"
paga_python_script <- file.path(script_dir, "04_trajectory_paga.py")
run_paga <- parse_bool(cfg_value(config, "Run_paga", "RUN_PAGA", default = "TRUE"), default = TRUE)
paga_n_neighbors <- parse_int_cfg(cfg_value(config, "Paga_n_neighbors", "PAGA_N_NEIGHBORS", default = "15"), 15L, "Paga_n_neighbors")
paga_n_pcs <- parse_int_cfg(cfg_value(config, "Paga_n_pcs", "PAGA_N_PCS", default = as.character(scvelo_n_pcs)), scvelo_n_pcs, "Paga_n_pcs")
paga_min_cells <- parse_int_cfg(cfg_value(config, "Paga_min_cells", "PAGA_MIN_CELLS", default = "20"), 20L, "Paga_min_cells")
paga_threshold <- parse_numeric_cfg(cfg_value(config, "Paga_threshold", "PAGA_THRESHOLD", default = "0.03"), 0.03, "Paga_threshold")
paga_pca_prefix <- cfg_value(config, "Paga_pca_prefix", "PAGA_PCA_PREFIX", default = scvelo_pca_prefix)
paga_root_clusters <- cfg_value(config, "Paga_root_clusters", "PAGA_ROOT_CLUSTERS", default = scvelo_root_clusters)
paga_export_expression <- parse_bool(cfg_value(config, "Paga_export_expression", "PAGA_EXPORT_EXPRESSION", default = "TRUE"), default = TRUE)
paga_expression_assay <- cfg_value(config, "Paga_expression_assay", "PAGA_EXPRESSION_ASSAY", default = "RNA")
paga_expression_slot <- cfg_value(config, "Paga_expression_slot", "PAGA_EXPRESSION_SLOT", default = "data")

dose_col <- "Dose"
id_col_candidates <- unique(c("IDs", cfg_value(config, "Trajectory_id_col", "TRAJECTORY_ID_COL", default = "ID"), "ID"))
ploidy_col <- "ploidy"
expected_doses <- c("0mg/kg", "30mg/kg", "120mg/kg")
cluster_col <- "cluster_final"

message("Preparing output directories.")
.ensure_dir(output_root)
out_inputs <- .ensure_dir(file.path(output_root, "00_inputs"))
out_velocity_groups <- .ensure_dir(file.path(output_root, "01_velocity_groups"))
out_summary <- .ensure_dir(file.path(output_root, "02_summary"))
out_plots <- .ensure_dir(file.path(output_root, "03_plots"))
out_paga_groups <- .ensure_dir(file.path(output_root, "04_paga_groups"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}
if (!file.exists(python_script)) {
  stop("Python worker script does not exist: ", python_script, call. = FALSE)
}
if (run_paga && !file.exists(paga_python_script)) {
  stop("PAGA Python worker script does not exist: ", paga_python_script, call. = FALSE)
}

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}
assert_required_column(obj@meta.data, dose_col)
if (!("umap" %in% names(obj@reductions))) {
  stop("UMAP reduction not found in Seurat object.", call. = FALSE)
}
if (!("pca" %in% names(obj@reductions))) {
  stop("PCA reduction not found in Seurat object. scVelo now requires Seurat PCA coordinates for neighbor construction.", call. = FALSE)
}

configured_sample_folder_col <- cellranger_sample_folder_col
resolved_sample_folder_col <- resolve_first_existing_col(
  obj@meta.data,
  c(configured_sample_folder_col, "sample_folder", "Sequencing.IDs", "sample", "orig.ident", "IDs")
)
if (!is.na(resolved_sample_folder_col)) {
  if (!identical(resolved_sample_folder_col, configured_sample_folder_col)) {
    message(
      "Using metadata column '",
      resolved_sample_folder_col,
      "' for sample-level loom discovery because configured column '",
      configured_sample_folder_col,
      "' was not found."
    )
  }
  cellranger_sample_folder_col <- resolved_sample_folder_col
}

if (!(cluster_col %in% colnames(obj@meta.data))) {
  cluster_col <- if ("clusters" %in% colnames(obj@meta.data)) {
    "clusters"
  } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    "seurat_clusters"
  } else {
    NA_character_
  }
}
if (is.na(cluster_col)) {
  stop("Cannot find cluster_final, clusters, or seurat_clusters metadata column.", call. = FALSE)
}

obj@meta.data[[dose_col]] <- standardize_in_vivo_dose(obj@meta.data[[dose_col]])
dose_values <- as.character(obj@meta.data[[dose_col]])
observed_dose_values <- unique(dose_values[!is.na(dose_values)])
missing_expected_doses <- setdiff(expected_doses, observed_dose_values)
if (length(missing_expected_doses) > 0) {
  stop("Missing expected Dose group(s): ", paste(missing_expected_doses, collapse = ", "), call. = FALSE)
}
dose_levels <- c(expected_doses, sort(setdiff(observed_dose_values, expected_doses)))
dose_levels <- dose_levels[dose_levels %in% observed_dose_values]
obj@meta.data[[dose_col]] <- factor(dose_values, levels = dose_levels)

umap <- as.data.frame(Embeddings(obj, "umap"))
if (ncol(umap) < 2) {
  stop("UMAP embedding has fewer than 2 dimensions.", call. = FALSE)
}
umap <- umap[, seq_len(2), drop = FALSE]
colnames(umap) <- c("UMAP_1", "UMAP_2")

pca <- as.data.frame(Embeddings(obj, "pca"))
if (ncol(pca) < 1) {
  stop("PCA embedding has no usable dimensions.", call. = FALSE)
}
pca <- pca[, seq_len(min(ncol(pca), scvelo_n_pcs)), drop = FALSE]
colnames(pca) <- paste0(scvelo_pca_prefix, seq_len(ncol(pca)))

meta_cols <- c(
  dose_col,
  id_col_candidates,
  cluster_col,
  "clusters",
  "sample",
  cellranger_sample_folder_col,
  "Sequencing.IDs",
  "orig.ident",
  "sample_type",
  "cluster_final_annotation_primary",
  "cluster_final_annotation_multi",
  "seurat_clusters",
  "nCount_RNA",
  "nFeature_RNA",
  "percent.mt"
)
meta_cols <- unique(meta_cols[meta_cols %in% colnames(obj@meta.data)])
meta_all <- obj@meta.data[, meta_cols, drop = FALSE] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::bind_cols(umap[.$cell, , drop = FALSE]) %>%
  dplyr::bind_cols(pca[.$cell, , drop = FALSE]) %>%
  dplyr::mutate(
    Dose = factor(as.character(.data[[dose_col]]), levels = dose_levels),
    dose_safe = sanitize_path_component(as.character(.data[[dose_col]]))
  )

if (!("sample_type" %in% colnames(meta_all)) && "sample" %in% colnames(meta_all)) {
  meta_all$sample_type <- infer_in_vivo_sample_type(meta_all$sample)
}

id_col <- resolve_col_case_insensitive(meta_all, id_col_candidates)
if (is.na(id_col)) {
  stop(
    "A metadata column is required to derive CellLine/Tumor context for 04_trajectory.R. Tried: ",
    paste(unique(id_col_candidates), collapse = ", "),
    call. = FALSE
  )
}
n_cells <- nrow(meta_all)
ploidy_from_id <- rep(NA_character_, n_cells)
ploidy_from_sample_folder <- rep(NA_character_, n_cells)
ploidy_from_sample <- rep(NA_character_, n_cells)

ploidy_from_id <- infer_ploidy_from_text(meta_all[[id_col]])
if (cellranger_sample_folder_col %in% colnames(meta_all)) {
  ploidy_from_sample_folder <- infer_ploidy_from_sample_label(meta_all[[cellranger_sample_folder_col]])
}
if ("sample" %in% colnames(meta_all)) {
  ploidy_from_sample <- infer_ploidy_from_sample_label(meta_all$sample)
}

meta_all[[ploidy_col]] <- ploidy_from_id
meta_all$ploidy_source <- ifelse(!is.na(ploidy_from_id), id_col, NA_character_)
missing_ploidy <- is.na(meta_all[[ploidy_col]])
meta_all[[ploidy_col]][missing_ploidy] <- ploidy_from_sample_folder[missing_ploidy]
meta_all$ploidy_source[missing_ploidy & !is.na(ploidy_from_sample_folder)] <- cellranger_sample_folder_col
missing_ploidy <- is.na(meta_all[[ploidy_col]])
meta_all[[ploidy_col]][missing_ploidy] <- ploidy_from_sample[missing_ploidy]
meta_all$ploidy_source[missing_ploidy & !is.na(ploidy_from_sample)] <- "sample"

if (any(is.na(meta_all[[ploidy_col]]))) {
  tried_cols <- c(if (!is.na(id_col)) id_col, cellranger_sample_folder_col, "sample")
  tried_cols <- unique(tried_cols[tried_cols %in% colnames(meta_all)])
  bad_cells <- unique(meta_all$cell[is.na(meta_all[[ploidy_col]])])
  stop(
    "Could not infer ploidy for 04_trajectory.R. Tried metadata column(s): ",
    paste(tried_cols, collapse = ", "),
    ". Example unresolved cell(s): ",
    paste(head(bad_cells, 20), collapse = ", "),
    call. = FALSE
  )
}
ploidy_levels <- c("2N", "4N")
ploidy_levels <- ploidy_levels[ploidy_levels %in% unique(meta_all[[ploidy_col]])]
meta_all[[ploidy_col]] <- factor(meta_all[[ploidy_col]], levels = ploidy_levels)

meta_all$trajectory_context <- factor(
  infer_trajectory_context_from_ids(meta_all[[id_col]]),
  levels = c("CellLine", "Tumor")
)
meta_all$Ploidy <- factor(as.character(meta_all[[ploidy_col]]), levels = c("2N", "4N"))
meta_all$TN <- factor(as.character(meta_all$trajectory_context), levels = c("CellLine", "Tumor"))
meta_all$Dose_DEG <- factor(as.character(meta_all[[dose_col]]), levels = expected_doses)
if (!("clusters" %in% colnames(meta_all))) {
  meta_all$clusters <- as.character(meta_all[[cluster_col]])
}
meta_all$trajectory_group <- factor(as.character(meta_all$TN), levels = c("CellLine", "Tumor"))

tn_ploidy_tasks <- expand.grid(
  tn_scope = c("CellLine", "Tumor"),
  ploidy_scope = c("ploidy_all", "2N", "4N"),
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    trajectory_group = paste(.data$tn_scope, .data$ploidy_scope, sep = "_"),
    trajectory_context = .data$tn_scope,
    Dose = NA_character_,
    ploidy = ifelse(.data$ploidy_scope %in% c("2N", "4N"), .data$ploidy_scope, NA_character_),
    group_safe = file.path(.data$tn_scope, .data$ploidy_scope),
    analysis_type = paste0(tolower(.data$tn_scope), "_", .data$ploidy_scope)
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    n_cells = {
      idx <- as.character(meta_all$TN) == .data$tn_scope
      if (.data$ploidy_scope %in% c("2N", "4N")) {
        idx <- idx & as.character(meta_all$Ploidy) == .data$ploidy_scope
      }
      sum(idx, na.rm = TRUE)
    }
  ) %>%
  dplyr::ungroup()
trajectory_group_df <- dplyr::bind_rows(
  data.frame(
    trajectory_group = "All_cells",
    trajectory_context = "All",
    Dose = NA_character_,
    ploidy = NA_character_,
    group_safe = "All_cells",
    analysis_type = "all_cells",
    tn_scope = "All",
    ploidy_scope = "ploidy_all",
    n_cells = nrow(meta_all),
    stringsAsFactors = FALSE
  ),
  tn_ploidy_tasks
)
trajectory_group_df <- trajectory_group_df %>%
  dplyr::filter(.data$n_cells > 0)
if (nrow(trajectory_group_df) == 0) {
  stop("No cells were found for the requested velocity analyses.", call. = FALSE)
}
trajectory_group_levels <- unique(trajectory_group_df$trajectory_group)

cellranger_sample_manifest <- data.frame()
if (cellranger_sample_folder_col %in% colnames(meta_all)) {
  message("Resolving Cell Ranger BAM and velocyto loom paths from metadata column: ", cellranger_sample_folder_col)
  cellranger_sample_manifest <- build_cellranger_sample_manifest(
    meta_all = meta_all,
    cellranger_root = cellranger_root,
    sample_folder_col = cellranger_sample_folder_col,
    dose_col = dose_col,
    ploidy_col = ploidy_col,
    velocyto_work_root = velocyto_work_root,
    velocyto_output_root = velocyto_output_root,
    velocyto_run_tag = velocyto_run_tag
  )
  cellranger_sample_manifest <- prepare_velocyto_runners(
    sample_manifest = cellranger_sample_manifest,
    output_dir = out_inputs,
    velocyto_bin = velocyto_bin,
    samtools_bin = samtools_bin,
    samtools_threads = velocyto_samtools_threads,
    gtf_file = velocyto_gtf
  )

  if (run_velocyto && nrow(cellranger_sample_manifest) > 0) {
    if (!nzchar(velocyto_gtf) || !file.exists(velocyto_gtf)) {
      stop("Run_velocyto is TRUE, but VELOCYTO_GTF/Velocyto_gtf is not set to an existing uncompressed .gtf file.", call. = FALSE)
    }
    if (grepl("\\.gz$", velocyto_gtf, ignore.case = TRUE)) {
      stop(
        "This velocyto version cannot read compressed .gtf.gz files. Decompress genes.gtf.gz and set VELOCYTO_GTF to the uncompressed genes.gtf.",
        call. = FALSE
      )
    }
    if (is.na(velocyto_bin) || !nzchar(velocyto_bin) || !file.exists(velocyto_bin)) {
      stop("velocyto binary not found: ", velocyto_bin, call. = FALSE)
    }
    if (is.na(samtools_bin) || !nzchar(samtools_bin) || !file.exists(samtools_bin)) {
      stop(
        "samtools binary not found. Install samtools in the velocity environment, load a samtools module, or set SAMTOOLS_BIN.",
        call. = FALSE
      )
    }

    n_existing_loom <- sum(cellranger_sample_manifest$velocyto_loom_exists, na.rm = TRUE)
    if (n_existing_loom > 0) {
      message("Found existing velocyto loom for ", n_existing_loom, " sample folder(s); skipping generation for those samples.")
    }
    missing_loom <- cellranger_sample_manifest %>%
      dplyr::filter(.data$cellranger_dir_exists, !.data$velocyto_loom_exists)
    if (nrow(missing_loom) > 0) {
      message("Running velocyto run10x for ", nrow(missing_loom), " sample folder(s).")
      for (i in seq_len(nrow(missing_loom))) {
        exit_status <- system2("bash", missing_loom$velocyto_run_script[i])
        if (!identical(exit_status, 0L)) {
          stop(
            "velocyto run10x failed for sample_folder ",
            missing_loom$sample_folder[i],
            ". See log: ",
            missing_loom$velocyto_log_file[i],
            call. = FALSE
          )
        }
      }
      cellranger_sample_manifest$velocyto_loom_file <- vapply(
        seq_len(nrow(cellranger_sample_manifest)),
        function(i) find_cellranger_loom(
          cellranger_sample_manifest$cellranger_dir[i],
          cellranger_sample_manifest$sample_folder[i],
          extra_dirs = c(
            cellranger_sample_manifest$velocyto_output_dir[i],
            cellranger_sample_manifest$velocyto_work_sample_root[i],
            file.path(cellranger_sample_manifest$velocyto_writable_sample_dir[i], "velocyto")
          )
        ),
        character(1)
      )
      cellranger_sample_manifest$velocyto_loom_exists <- !is.na(cellranger_sample_manifest$velocyto_loom_file) &
        file.exists(cellranger_sample_manifest$velocyto_loom_file)
    }
  }

  write_table_csv(cellranger_sample_manifest, file.path(out_inputs, "cellranger_sample_bam_manifest.csv"))
} else {
  warning(
    "Metadata column '",
    cellranger_sample_folder_col,
    "' was not found. Cell Ranger BAM/velocyto loom discovery will be skipped."
  )
}

write_table_csv(meta_all, file.path(out_inputs, "seurat_cells_metadata_umap.csv"))
write_table_csv(trajectory_group_df, file.path(out_summary, "trajectory_group_levels.csv"))

cell_count_by_dose <- meta_all %>%
  dplyr::count(Dose, name = "n_cells") %>%
  dplyr::arrange(match(as.character(Dose), dose_levels))
write_table_csv(cell_count_by_dose, file.path(out_summary, "dose_cell_counts.csv"))

cell_count_by_group <- trajectory_group_df %>%
  dplyr::select(
    trajectory_group,
    tn_scope,
    ploidy_scope,
    trajectory_context,
    Dose,
    ploidy,
    n_cells,
    analysis_type,
    group_safe
  )
write_table_csv(cell_count_by_group, file.path(out_summary, "trajectory_group_cell_counts.csv"))

message("Writing Seurat UMAP overview plots.")
dose_colors <- c(
  "0mg/kg" = "#607D8B",
  "30mg/kg" = "#D9902F",
  "120mg/kg" = "#B33434"
)
observed_dose_colors <- dose_colors[names(dose_colors) %in% as.character(dose_levels)]

p_dose <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2, color = Dose)) +
  geom_point(size = 0.25, alpha = 0.85, stroke = 0) +
  coord_equal() +
  scale_color_manual(values = observed_dose_colors, na.value = "grey80") +
  labs(title = "UMAP by Dose", subtitle = "Input cells for trajectory velocity analysis") +
  make_pretty_umap_theme()
save_both(p_dose, file.path(out_plots, "umap_by_Dose"), width = 9, height = 7)

cluster_levels <- sort_maybe_numeric(as.character(meta_all[[cluster_col]]))
cluster_palette <- setNames(grDevices::hcl.colors(length(cluster_levels), palette = "Dark 3"), cluster_levels)
p_cluster <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2, color = as.character(.data[[cluster_col]]))) +
  geom_point(size = 0.25, alpha = 0.85, stroke = 0) +
  coord_equal() +
  scale_color_manual(values = cluster_palette, na.value = "grey80") +
  labs(title = "UMAP by cluster", color = cluster_col) +
  make_pretty_umap_theme()
save_both(p_cluster, file.path(out_plots, "umap_by_cluster"), width = 9, height = 7)

p_cluster_by_dose <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2, color = as.character(.data[[cluster_col]]))) +
  geom_point(size = 0.18, alpha = 0.85, stroke = 0) +
  coord_equal() +
  facet_wrap(~Dose) +
  scale_color_manual(values = cluster_palette, na.value = "grey80") +
  labs(title = "UMAP by cluster within each Dose", color = cluster_col) +
  make_pretty_umap_theme()
save_both(p_cluster_by_dose, file.path(out_plots, "umap_by_cluster_facet_Dose"), width = 11, height = 7)

p_cluster_by_tn <- ggplot(meta_all, aes(x = UMAP_1, y = UMAP_2, color = as.character(.data[[cluster_col]]))) +
  geom_point(size = 0.18, alpha = 0.85, stroke = 0) +
  coord_equal() +
  facet_wrap(~TN) +
  scale_color_manual(values = cluster_palette, na.value = "grey80") +
  labs(title = "UMAP by cluster within Tumor and CellLine cells", color = cluster_col) +
  make_pretty_umap_theme()
save_both(p_cluster_by_tn, file.path(out_plots, "umap_by_cluster_facet_TN"), width = 10, height = 7)

velocity_inputs <- as_velocity_input_list(resolve_velocity_inputs(
  dose_levels = dose_levels,
  manifest_path = velocity_manifest,
  velocity_input_root = velocity_input_root
))

sample_loom_inputs <- if (nrow(cellranger_sample_manifest) > 0) {
  resolve_sample_loom_inputs_by_dose(cellranger_sample_manifest, dose_levels)
} else {
  stats::setNames(vector("list", length(dose_levels)), dose_levels)
}
for (dose in dose_levels) {
  if (length(velocity_inputs[[dose]]) == 0 && length(sample_loom_inputs[[dose]]) > 0) {
    velocity_inputs[[dose]] <- sample_loom_inputs[[dose]]
  }
}
global_velocity_input_files <- unique(unlist(velocity_inputs, use.names = FALSE))
global_velocity_input_files <- global_velocity_input_files[!is.na(global_velocity_input_files) & nzchar(global_velocity_input_files)]
has_single_global_velocity_input <- length(global_velocity_input_files) == 1 &&
  all(vapply(velocity_inputs, function(files) identical(files, global_velocity_input_files), logical(1)))

velocity_manifest_df <- data.frame(
  Dose = names(velocity_inputs),
  velocity_input_file = vapply(velocity_inputs, format_input_files, character(1)),
  input_exists = vapply(velocity_inputs, all_input_files_exist, logical(1)),
  n_velocity_input_files = vapply(velocity_inputs, length, integer(1)),
  source = vapply(names(velocity_inputs), function(dose) {
    if (length(velocity_inputs[[dose]]) == 0) return("missing")
    if (identical(velocity_inputs[[dose]], sample_loom_inputs[[dose]])) return("cellranger_sample_loom")
    "velocity_input_root_or_manifest"
  }, character(1)),
  stringsAsFactors = FALSE
)
write_table_csv(velocity_manifest_df, file.path(out_summary, "resolved_velocity_inputs.csv"))

color_cols <- unique(c(cluster_col, "clusters", dose_col, ploidy_col, "Ploidy", "TN", "Dose_DEG", "trajectory_context", "trajectory_group", "trajectory_analysis_group", "trajectory_shape_group", "trajectory_tn_scope", "trajectory_ploidy_scope", "trajectory_branch", "sample_type", "cluster_final_annotation_primary"))
color_cols <- color_cols[color_cols %in% c(colnames(meta_all), "trajectory_analysis_group", "trajectory_shape_group", "trajectory_tn_scope", "trajectory_ploidy_scope", "trajectory_branch")]
get_analysis_group_meta <- function(task_row) {
  group_name <- task_row$trajectory_group[1]
  if (identical(group_name, "All_cells")) {
    meta_all %>%
      dplyr::mutate(
        trajectory_shape_group = factor(as.character(.data$TN), levels = c("CellLine", "Tumor"))
      )
  } else {
    tn_scope <- task_row$tn_scope[1]
    ploidy_scope <- task_row$ploidy_scope[1]
    out <- meta_all %>%
      dplyr::filter(as.character(.data$TN) == tn_scope)
    if (ploidy_scope %in% c("2N", "4N")) {
      out <- out %>% dplyr::filter(as.character(.data$Ploidy) == ploidy_scope)
    }
    out %>%
      dplyr::mutate(
        trajectory_shape_group = ifelse(
          ploidy_scope == "ploidy_all",
          as.character(.data$Ploidy),
          ploidy_scope
        )
      )
  }
}

resolve_analysis_group_inputs <- function(group_meta) {
  selected_samples <- character(0)
  missing_samples <- character(0)
  available_sample_inputs <- character(0)
  if (
    nrow(group_meta) > 0 &&
      nrow(cellranger_sample_manifest) > 0 &&
      cellranger_sample_folder_col %in% colnames(group_meta) &&
      "sample_folder" %in% colnames(cellranger_sample_manifest)
  ) {
    selected_samples <- sort(unique(as.character(group_meta[[cellranger_sample_folder_col]])))
    selected_samples <- selected_samples[!is.na(selected_samples) & nzchar(selected_samples)]
    manifest_rows <- cellranger_sample_manifest[cellranger_sample_manifest$sample_folder %in% selected_samples, , drop = FALSE]
    missing_manifest_samples <- setdiff(selected_samples, manifest_rows$sample_folder)
    missing_loom_idx <- is.na(manifest_rows$velocyto_loom_exists) | !manifest_rows$velocyto_loom_exists
    missing_loom_samples <- manifest_rows$sample_folder[missing_loom_idx]
    missing_samples <- sort(unique(c(missing_manifest_samples, missing_loom_samples)))
    available_sample_inputs <- unique(as.character(manifest_rows$velocyto_loom_file[manifest_rows$velocyto_loom_exists]))
    available_sample_inputs <- available_sample_inputs[!is.na(available_sample_inputs) & nzchar(available_sample_inputs)]
    if (length(selected_samples) > 0 && length(missing_samples) == 0 && length(available_sample_inputs) > 0) {
      return(list(
        input_files = available_sample_inputs,
        source = "cellranger_sample_loom",
        selected_sample_folders = selected_samples,
        missing_sample_folders = character(0),
        available_sample_inputs = available_sample_inputs
      ))
    }
  }

  fallback_files <- if (isTRUE(has_single_global_velocity_input)) {
    global_velocity_input_files
  } else if (length(global_velocity_input_files) > 0) {
    global_velocity_input_files
  } else {
    character(0)
  }
  fallback_source <- if (length(fallback_files) == 0) {
    "missing"
  } else if (isTRUE(has_single_global_velocity_input)) {
    "global_velocity_input"
  } else {
    "combined_dose_velocity_input"
  }
  if (length(selected_samples) > 0 && length(missing_samples) > 0 && length(fallback_files) == 0) {
    fallback_source <- "cellranger_sample_loom_incomplete"
  }
  list(
    input_files = fallback_files,
    source = fallback_source,
    selected_sample_folders = selected_samples,
    missing_sample_folders = missing_samples,
    available_sample_inputs = available_sample_inputs
  )
}

message("Preparing PAGA and scVelo analyses: All_cells plus CellLine/Tumor ploidy_all, 2N, and 4N.")
run_rows <- list()
paga_rows <- list()
analysis_group_meta <- stats::setNames(vector("list", length(trajectory_group_levels)), trajectory_group_levels)
analysis_group_inputs <- stats::setNames(vector("list", length(trajectory_group_levels)), trajectory_group_levels)
analysis_group_input_info <- stats::setNames(vector("list", length(trajectory_group_levels)), trajectory_group_levels)

for (i in seq_len(nrow(trajectory_group_df))) {
  group_name <- trajectory_group_df$trajectory_group[i]
  group_meta <- get_analysis_group_meta(trajectory_group_df[i, , drop = FALSE])
  input_info <- resolve_analysis_group_inputs(group_meta)
  analysis_group_meta[[group_name]] <- group_meta
  analysis_group_inputs[[group_name]] <- input_info$input_files
  analysis_group_input_info[[group_name]] <- input_info

  paga_rows[[group_name]] <- prepare_or_run_paga_analysis(
    group_name = group_name,
    group_safe = trajectory_group_df$group_safe[i],
    analysis_type = trajectory_group_df$analysis_type[i],
    analysis_root = out_paga_groups,
    metadata_df = group_meta,
    seurat_obj = obj,
    export_expression = paga_export_expression,
    expression_assay = paga_expression_assay,
    expression_slot = paga_expression_slot,
    color_cols = color_cols,
    dose = NA_character_,
    ploidy = trajectory_group_df$ploidy[i],
    trajectory_context = trajectory_group_df$trajectory_context[i],
    base_groups = group_name,
    tn_scope = trajectory_group_df$tn_scope[i],
    ploidy_scope = trajectory_group_df$ploidy_scope[i]
  )

  run_rows[[group_name]] <- prepare_or_run_scvelo_analysis(
    group_name = group_name,
    group_safe = trajectory_group_df$group_safe[i],
    analysis_type = trajectory_group_df$analysis_type[i],
    analysis_root = out_velocity_groups,
    metadata_df = group_meta,
    input_files = input_info$input_files,
    color_cols = color_cols,
    dose = NA_character_,
    ploidy = trajectory_group_df$ploidy[i],
    trajectory_context = trajectory_group_df$trajectory_context[i],
    base_groups = group_name,
    tn_scope = trajectory_group_df$tn_scope[i],
    ploidy_scope = trajectory_group_df$ploidy_scope[i],
    shape_col = if (identical(group_name, "All_cells") || identical(trajectory_group_df$ploidy_scope[i], "ploidy_all")) "trajectory_shape_group" else "",
    shape_order = if (identical(group_name, "All_cells")) c("CellLine", "Tumor") else if (identical(trajectory_group_df$ploidy_scope[i], "ploidy_all")) c("2N", "4N") else character(0),
    shape_markers = if (identical(group_name, "All_cells")) c("o", "^") else if (identical(trajectory_group_df$ploidy_scope[i], "ploidy_all")) c("o", "^") else character(0)
  )
}

paga_run_summary_df <- dplyr::bind_rows(paga_rows)
write_table_csv(paga_run_summary_df, file.path(out_summary, "paga_run_summary.csv"))
paga_edge_files <- file.path(paga_run_summary_df$paga_dir, "paga_edges.csv")
paga_edge_files <- paga_edge_files[file.exists(paga_edge_files)]
if (length(paga_edge_files) > 0) {
  paga_edges_all <- dplyr::bind_rows(lapply(paga_edge_files, function(path) {
    group_name <- paga_run_summary_df$analysis_group[match(dirname(path), paga_run_summary_df$paga_dir)]
    read_trajectory_csv(
      path,
      col_types = readr::cols(
        group_1 = readr::col_character(),
        group_2 = readr::col_character(),
        connectivity = readr::col_double(),
        above_threshold = readr::col_logical()
      )
    ) %>%
      dplyr::mutate(analysis_group = group_name, paga_dir = dirname(path), .before = 1)
  }))
  write_table_csv(paga_edges_all, file.path(out_summary, "paga_edges_all_trajectory_analyses.csv"))
}
paga_metric_files <- file.path(paga_run_summary_df$paga_dir, "paga_cell_metrics.csv")
paga_metric_files <- paga_metric_files[file.exists(paga_metric_files)]
if (length(paga_metric_files) > 0) {
  paga_metrics_all <- dplyr::bind_rows(lapply(paga_metric_files, function(path) {
    group_name <- paga_run_summary_df$analysis_group[match(dirname(path), paga_run_summary_df$paga_dir)]
    read_trajectory_csv(path) %>%
      dplyr::mutate(analysis_group = group_name, paga_dir = dirname(path), .before = 1)
  }))
  write_table_csv(paga_metrics_all, file.path(out_summary, "paga_cell_metrics_all_trajectory_analyses.csv"))
}

group_velocity_manifest_df <- data.frame(
  trajectory_group = trajectory_group_df$trajectory_group,
  tn_scope = trajectory_group_df$tn_scope,
  ploidy_scope = trajectory_group_df$ploidy_scope,
  trajectory_context = as.character(trajectory_group_df$trajectory_context),
  Dose = as.character(trajectory_group_df$Dose),
  ploidy = as.character(trajectory_group_df$ploidy),
  velocity_input_file = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    format_input_files(get_named_input_files(analysis_group_inputs, group_name))
  }, character(1)),
  input_exists = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    all_input_files_exist(get_named_input_files(analysis_group_inputs, group_name))
  }, logical(1)),
  n_velocity_input_files = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    length(get_named_input_files(analysis_group_inputs, group_name))
  }, integer(1)),
  source = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    analysis_group_input_info[[group_name]]$source
  }, character(1)),
  selected_sample_folders = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    paste(analysis_group_input_info[[group_name]]$selected_sample_folders, collapse = ";")
  }, character(1)),
  missing_sample_folders = vapply(trajectory_group_df$trajectory_group, function(group_name) {
    paste(analysis_group_input_info[[group_name]]$missing_sample_folders, collapse = ";")
  }, character(1)),
  stringsAsFactors = FALSE
)
write_table_csv(group_velocity_manifest_df, file.path(out_summary, "resolved_velocity_inputs_by_trajectory_group.csv"))

run_summary_df <- dplyr::bind_rows(run_rows)
write_table_csv(run_summary_df, file.path(out_summary, "scvelo_run_summary.csv"))

missing_inputs <- run_summary_df %>%
  dplyr::filter(status == "missing_velocity_input")
if (nrow(missing_inputs) > 0) {
  requirement_lines <- c(
    "RNA velocity input is missing for one or more trajectory analyses.",
    "",
    "Required input:",
    "  A .loom or .h5ad containing spliced and unspliced layers.",
    "",
    "Supported ways to provide input:",
    "  1. Set VELOCITY_INPUT_ROOT to a directory containing one combined .loom/.h5ad, or per-Dose files named with 0mgkg, 30mgkg, 120mgkg.",
    "  2. Set VELOCITY_INPUT_MANIFEST to a CSV with columns: dose,input_file.",
    "  3. If cell IDs differ between Seurat and velocity input, set VELOCITY_CELL_MAP to a CSV with columns: seurat_cell,velocity_cell.",
    "  4. Put existing per-sample loom files under 00_inputs/velocyto_loom/<sample_folder>/<sample_folder>.loom; this is checked before any velocyto generation.",
    "  5. If Cell Ranger folders are available, set CELLRANGER_ROOT and CELLRANGER_SAMPLE_FOLDER_COL; this script writes velocyto run10x scripts from sample_folder and can use existing per-sample velocyto .loom files.",
    "  6. To generate missing loom files from Cell Ranger BAMs, set VELOCYTO_GTF to a gene annotation GTF and RUN_VELOCYTO=TRUE.",
    "",
    paste0("Default searched root: ", velocity_input_root),
    paste0("Per-sample loom root: ", velocyto_output_root),
    paste0("Cell Ranger root: ", cellranger_root),
    paste0("Cell Ranger sample folder metadata column: ", cellranger_sample_folder_col),
    paste0("Velocyto writable work root: ", velocyto_work_root),
    paste0("Velocyto output root: ", velocyto_output_root),
    paste0("samtools binary: ", ifelse(is.na(samtools_bin), "not found", samtools_bin)),
    paste0("Velocyto samtools threads: ", velocyto_samtools_threads),
    "",
    "Per-analysis Seurat metadata and run_scvelo.sh files have already been written under 01_velocity_groups/."
  )
  writeLines(requirement_lines, con = file.path(out_summary, "velocity_input_requirements.txt"))
  stop(
    "Missing velocity input for trajectory analysis group(s): ",
    paste(missing_inputs$analysis_group, collapse = ", "),
    ". See ",
    file.path(out_summary, "velocity_input_requirements.txt"),
    call. = FALSE
  )
}

metric_files <- file.path(run_summary_df$dose_dir, "scvelo_cell_metrics.csv")
metric_files <- metric_files[file.exists(metric_files)]
if (length(metric_files) > 0) {
  metric_df <- dplyr::bind_rows(lapply(metric_files, function(path) {
    read_trajectory_csv(path)
  }))
  metric_join_cols <- unique(c(
    "cell",
    "UMAP_1",
    "UMAP_2",
    "Dose",
    "ploidy",
    "Ploidy",
    "TN",
    "Dose_DEG",
    "trajectory_context",
    "trajectory_group",
    "trajectory_analysis_group",
    "trajectory_shape_group",
    "trajectory_tn_scope",
    "trajectory_ploidy_scope",
    "trajectory_branch",
    "clusters"
  ))
  metric_join_cols <- metric_join_cols[metric_join_cols %in% colnames(meta_all)]
  metric_df <- metric_df %>%
    dplyr::left_join(
      meta_all %>% dplyr::select(dplyr::all_of(metric_join_cols)),
      by = "cell",
      suffix = c("", "_seurat")
    )
  for (col in c("Dose", "ploidy", "Ploidy", "TN", "Dose_DEG", "trajectory_context", "trajectory_group", "trajectory_analysis_group", "trajectory_shape_group", "trajectory_tn_scope", "trajectory_ploidy_scope", "trajectory_branch", "clusters")) {
    if (!(col %in% colnames(metric_df))) {
      metric_df[[col]] <- NA_character_
    }
    seurat_col <- paste0(col, "_seurat")
    if (!(seurat_col %in% colnames(metric_df))) {
      metric_df[[seurat_col]] <- metric_df[[col]]
    }
  }
  metric_df <- metric_df %>%
    dplyr::mutate(
      Dose = dplyr::coalesce(.data$Dose, .data$Dose_seurat),
      Dose = factor(as.character(Dose), levels = dose_levels),
      ploidy = dplyr::coalesce(.data$ploidy, .data$ploidy_seurat),
      ploidy = factor(as.character(ploidy), levels = ploidy_levels),
      Ploidy = dplyr::coalesce(.data$Ploidy, .data$Ploidy_seurat),
      Ploidy = factor(as.character(Ploidy), levels = c("2N", "4N")),
      TN = dplyr::coalesce(.data$TN, .data$TN_seurat),
      TN = factor(as.character(TN), levels = c("CellLine", "Tumor")),
      Dose_DEG = dplyr::coalesce(.data$Dose_DEG, .data$Dose_DEG_seurat),
      Dose_DEG = factor(as.character(Dose_DEG), levels = expected_doses),
      trajectory_context = dplyr::coalesce(.data$trajectory_context, .data$trajectory_context_seurat),
      trajectory_group = dplyr::coalesce(.data$trajectory_group, .data$trajectory_group_seurat),
      trajectory_group = factor(as.character(trajectory_group), levels = c("CellLine", "Tumor")),
      trajectory_analysis_group = dplyr::coalesce(.data$trajectory_analysis_group, .data$trajectory_analysis_group_seurat),
      trajectory_analysis_group = factor(as.character(trajectory_analysis_group), levels = trajectory_group_levels),
      trajectory_shape_group = dplyr::coalesce(.data$trajectory_shape_group, .data$trajectory_shape_group_seurat),
      trajectory_tn_scope = dplyr::coalesce(.data$trajectory_tn_scope, .data$trajectory_tn_scope_seurat),
      trajectory_ploidy_scope = dplyr::coalesce(.data$trajectory_ploidy_scope, .data$trajectory_ploidy_scope_seurat),
      trajectory_branch = dplyr::coalesce(.data$trajectory_branch, .data$trajectory_branch_seurat)
    )
  write_table_csv(metric_df, file.path(out_summary, "scvelo_cell_metrics_all_trajectory_analyses.csv"))

  metric_single_df <- metric_df %>%
    dplyr::filter(is.na(.data$trajectory_analysis_group) | .data$trajectory_analysis_group %in% trajectory_group_levels)

  plot_metric_umap(
    metric_single_df,
    metric = "velocity_pseudotime",
    file_stub = file.path(out_plots, "umap_velocity_pseudotime_by_velocity_analysis_group"),
    title = "Velocity pseudotime by velocity analysis group"
  )
  plot_metric_umap(
    metric_single_df,
    metric = "latent_time",
    file_stub = file.path(out_plots, "umap_latent_time_by_velocity_analysis_group"),
    title = "Latent time by velocity analysis group"
  )
  plot_metric_umap(
    metric_single_df,
    metric = "velocity_confidence",
    file_stub = file.path(out_plots, "umap_velocity_confidence_by_velocity_analysis_group"),
    title = "Velocity confidence by velocity analysis group"
  )
}

summary_lines <- c(
  "04 trajectory velocity workflow completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Trajectory branch: ", trajectory_branch),
  paste0("scVelo neighbors: ", ifelse(scvelo_use_metadata_pca, paste0("Seurat PCA columns with prefix ", scvelo_pca_prefix), "scVelo automatic PCA/neighbors")),
  paste0("scVelo root cluster(s): ", scvelo_root_clusters),
  paste0("scVelo end cluster(s): ", ifelse(nzchar(scvelo_end_clusters), scvelo_end_clusters, "not specified")),
  paste0("Run PAGA immediately: ", run_paga),
  paste0("PAGA neighbors: ", paga_n_neighbors),
  paste0("PAGA PCs: ", paga_n_pcs),
  paste0("PAGA threshold: ", paga_threshold),
  paste0("PAGA root cluster(s) for DPT: ", ifelse(nzchar(paga_root_clusters), paga_root_clusters, "not specified")),
  paste0("PAGA h5ad expression export: ", paga_export_expression),
  paste0("PAGA h5ad expression source requested: ", paga_expression_assay, ":", paga_expression_slot),
  paste0("Dose groups: ", paste(dose_levels, collapse = ", ")),
  paste0("Ploidy groups: ", paste(ploidy_levels, collapse = ", ")),
  paste0("Velocity analysis groups: ", paste(trajectory_group_levels, collapse = ", ")),
  paste0("Cluster column: ", cluster_col),
  paste0("Velocity input root: ", velocity_input_root),
  paste0("Velocity input manifest: ", ifelse(nzchar(velocity_manifest), velocity_manifest, "not set")),
  paste0("Per-sample loom root: ", velocyto_output_root),
  "Per-sample loom layout: 00_inputs/velocyto_loom/<sample_folder>/<sample_folder>.loom",
  paste0("Cell Ranger root: ", cellranger_root),
  paste0("Cell Ranger sample folder column: ", cellranger_sample_folder_col),
  paste0("Velocyto binary: ", velocyto_bin),
  paste0("samtools binary: ", ifelse(is.na(samtools_bin), "not found", samtools_bin)),
  paste0("Velocyto samtools threads: ", velocyto_samtools_threads),
  paste0("Velocyto GTF: ", ifelse(nzchar(velocyto_gtf), velocyto_gtf, "not set")),
  paste0("Velocyto writable work root: ", velocyto_work_root),
  paste0("Velocyto output root: ", velocyto_output_root),
  paste0("Run velocyto immediately: ", run_velocyto),
  paste0("Cell map file: ", ifelse(nzchar(cell_map_file), cell_map_file, "not set")),
  paste0("Required Python environment: ", required_python_env),
  paste0("Python binary: ", python_bin),
  paste0("scVelo mode: ", scvelo_mode),
  paste0("Run scVelo immediately: ", run_scvelo),
  "",
  "Velocity outputs:",
  "  00_inputs/cellranger_sample_bam_manifest.csv",
  "  00_inputs/velocyto_run10x/*_run_velocyto.sh",
  "  01_velocity_groups/All_cells/cells_metadata_umap.csv",
  "  01_velocity_groups/CellLine/ploidy_all/cells_metadata_umap.csv",
  "  01_velocity_groups/CellLine/2N/cells_metadata_umap.csv",
  "  01_velocity_groups/CellLine/4N/cells_metadata_umap.csv",
  "  01_velocity_groups/Tumor/ploidy_all/cells_metadata_umap.csv",
  "  01_velocity_groups/Tumor/2N/cells_metadata_umap.csv",
  "  01_velocity_groups/Tumor/4N/cells_metadata_umap.csv",
  "  01_velocity_groups/**/run_scvelo.sh",
  "  01_velocity_groups/**/plots/velocity_stream.pdf(.png), after scVelo runs",
  "",
  "PAGA outputs:",
  "  04_paga_groups/All_cells/cells_metadata_umap.csv",
  "  04_paga_groups/**/expression_matrix.mtx, expression_genes.csv, expression_cells.csv",
  "  04_paga_groups/**/paga_result.h5ad",
  "  04_paga_groups/CellLine/ploidy_all/paga_edges.csv",
  "  04_paga_groups/Tumor/ploidy_all/paga_connectivities.csv",
  "  04_paga_groups/**/plots/paga_graph.pdf(.png)",
  "  04_paga_groups/**/plots/paga_umap_overlay.pdf(.png)",
  "  02_summary/paga_run_summary.csv",
  "  02_summary/paga_edges_all_trajectory_analyses.csv",
  "  02_summary/paga_cell_metrics_all_trajectory_analyses.csv",
  "",
  "R overview plots:",
  "  03_plots/umap_by_Dose.pdf(.png)",
  "  03_plots/umap_by_cluster.pdf(.png)",
  "  03_plots/umap_by_cluster_facet_Dose.pdf(.png)",
  "  03_plots/umap_by_cluster_facet_TN.pdf(.png)"
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("04 trajectory workflow finished.")
message("Output root: ", normalizePath(output_root, mustWork = FALSE))
