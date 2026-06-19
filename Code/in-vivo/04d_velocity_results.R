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
  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  candidate_dirs <- unique(c(candidate_dirs, cwd, file.path(cwd, "Code", "in-vivo")))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04c_extra_helpers.R"))

required_packages <- c("dplyr", "ggplot2", "readr", "readxl", "tibble", "tidyr")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(readxl)
  library(tibble)
  library(tidyr)
})

set.seed(1234)

root_label <- "6"
scenario_id <- "ROOT_6_END_NULL"
monocle3_root_id <- "ROOT_6"
analysis_objects <- c("All_cells", "CellLine/ploidy_all", "Tumor/ploidy_all")

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

is_truthy_env <- function(name, default = FALSE) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  tolower(trimws(value)) %in% c("1", "true", "t", "yes", "y")
}

safe_num <- function(x) suppressWarnings(as.numeric(x))

safe_file_component <- function(x) {
  if (exists("sanitize_path_component", mode = "function")) {
    return(sanitize_path_component(x))
  }
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "value")
}

standardize_output_text <- function(x) {
  if (is.null(x)) return(x)
  if (!is.character(x)) return(x)
  out <- x
  out <- gsub("continuous[- ]cell[- ]ploidy", "Cell Ploidy", out, ignore.case = TRUE)
  out <- gsub("continuous[- ]ploidy", "Cell Ploidy", out, ignore.case = TRUE)
  out <- gsub("cell[- ]ploidy", "Cell Ploidy", out, ignore.case = TRUE)
  out <- gsub("annotation[- ]ploidy", "Initial Ploidy", out, ignore.case = TRUE)
  out <- gsub("Within-annotation", "Within-Initial-Ploidy", out, fixed = TRUE)
  out <- gsub("within annotation", "within Initial Ploidy", out, ignore.case = TRUE)
  out[out == "Ploidy"] <- "Initial Ploidy"
  out
}

standardize_output_path <- function(path) {
  if (is.null(path)) return(path)
  out <- as.character(path)
  pre_replacements <- c(
    "continuous_cell_ploidy" = "cell_ploidy",
    "continuous_ploidy" = "cell_ploidy",
    "continuous_high" = "cell_ploidy_high",
    "continuous_mean" = "cell_ploidy_mean",
    "continuous_subset" = "cell_ploidy_subset",
    "tumor_continuous" = "tumor_cell_ploidy",
    "cell_level_continuous" = "cell_level_cell_ploidy",
    "sample_level_continuous" = "sample_level_cell_ploidy",
    "n_continuous" = "n_cell_ploidy",
    "continuous_" = "cell_ploidy_",
    "continuous" = "cell_ploidy"
  )
  for (old in names(pre_replacements)) {
    out <- gsub(old, pre_replacements[[old]], out, fixed = TRUE)
  }
  out <- gsub("initial_ploidy", "__INITIAL_PLOIDY__", out, fixed = TRUE)
  out <- gsub("initial_4N", "__INITIAL_4N__", out, fixed = TRUE)
  out <- gsub("cell_ploidy", "__CELL_PLOIDY__", out, fixed = TRUE)
  replacements <- c(
    "08_scvelo_ploidy_dose_cluster_deep_dive" = "08_scvelo_initial_ploidy_cell_ploidy_dose_cluster_deep_dive",
    "07_tumor_ploidy_dose_response" = "07_tumor_initial_ploidy_dose_response",
    "04_tumor_ploidy_dose_pseudotime" = "04_tumor_initial_ploidy_dose_pseudotime",
    "01_ploidy_pseudotime" = "01_initial_ploidy_pseudotime",
    "annotation_ploidy" = "initial_ploidy",
    "annotation_4N" = "initial_4N",
    "dose_within_ploidy" = "dose_within_initial_ploidy",
    "within_ploidy" = "within_initial_ploidy",
    "ploidy_pseudotime" = "initial_ploidy_pseudotime",
    "tumor_ploidy_dose" = "tumor_initial_ploidy_dose"
  )
  for (old in names(replacements)) {
    out <- gsub(old, replacements[[old]], out, fixed = TRUE)
  }
  out <- gsub("__INITIAL_PLOIDY__", "initial_ploidy", out, fixed = TRUE)
  out <- gsub("__INITIAL_4N__", "initial_4N", out, fixed = TRUE)
  out <- gsub("__CELL_PLOIDY__", "cell_ploidy", out, fixed = TRUE)
  out
}

standardize_output_columns <- function(df) {
  if (!is.data.frame(df)) return(df)
  names(df) <- standardize_output_path(names(df))
  char_cols <- vapply(df, is.character, logical(1))
  df[char_cols] <- lapply(df[char_cols], function(x) {
    x <- standardize_output_path(x)
    standardize_output_text(x)
  })
  df
}

standardize_ggplot_output_text <- function(plot_obj) {
  if (!inherits(plot_obj, "ggplot")) return(plot_obj)
  if (!is.null(plot_obj$labels) && length(plot_obj$labels) > 0) {
    for (label_name in names(plot_obj$labels)) {
      plot_obj$labels[[label_name]] <- standardize_output_text(plot_obj$labels[[label_name]])
    }
  }
  if (!is.null(plot_obj$scales) && length(plot_obj$scales$scales) > 0) {
    for (i in seq_along(plot_obj$scales$scales)) {
      scale_obj <- plot_obj$scales$scales[[i]]
      if (!is.null(scale_obj$name) && is.character(scale_obj$name)) {
        plot_obj$scales$scales[[i]]$name <- standardize_output_text(scale_obj$name)
      }
    }
  }
  plot_obj
}

base_ensure_dir <- .ensure_dir
.ensure_dir <- function(path) {
  base_ensure_dir(standardize_output_path(path))
}

save_pdf <- function(plot_obj, file, width = 8, height = 6) {
  file <- standardize_output_path(file)
  plot_obj <- standardize_ggplot_output_text(plot_obj)
  .ensure_dir(dirname(file))
  ggplot2::ggsave(file, plot_obj, width = width, height = height, limitsize = FALSE)
  invisible(file)
}

write_csv_safe <- function(df, file) {
  file <- standardize_output_path(file)
  df <- standardize_output_columns(df)
  .ensure_dir(dirname(file))
  write_table_csv(df, file)
  invisible(file)
}

first_existing_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) > 0) return(hit[1])
  if (isTRUE(required)) {
    stop("Cannot find ", label, ". Tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NA_character_
}

read_sample_info_map <- function(sample_info_file) {
  if (!file.exists(sample_info_file)) {
    stop("Missing sample info file: ", sample_info_file, call. = FALSE)
  }
  sample_info <- readxl::read_excel(sample_info_file)
  sample_info <- as.data.frame(sample_info, stringsAsFactors = FALSE)
  harvest_col <- first_existing_col(sample_info, c("harvest"), label = "harvest column")
  sample_col <- first_existing_col(sample_info, c("IDs", "sample", "sample_id", "sampleID"), label = "sample ID column")
  seq_col <- first_existing_col(sample_info, c("Sequencing IDs", "Sequencing.IDs", "sequencing_ids"), required = FALSE)
  dose_col <- first_existing_col(sample_info, c("Dose", "dose"), required = FALSE)

  out <- data.frame(
    sample = as_clean_chr(sample_info[[sample_col]]),
    harvest = as_clean_chr(sample_info[[harvest_col]]),
    sequencing_id = if (!is.na(seq_col)) as_clean_chr(sample_info[[seq_col]]) else NA_character_,
    sample_info_dose = if (!is.na(dose_col)) as_clean_chr(sample_info[[dose_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$sample) & nzchar(out$sample) & !is.na(out$harvest) & nzchar(out$harvest), , drop = FALSE]
  out <- out[!duplicated(out$sample), , drop = FALSE]
  out$sample_info_dose <- standardize_in_vivo_dose(out$sample_info_dose)
  out
}

read_cell_ploidy_map <- function(ploidy_file) {
  if (!file.exists(ploidy_file)) {
    stop("Missing cell ploidy file: ", ploidy_file, call. = FALSE)
  }
  ploidy_raw <- readr::read_tsv(
    ploidy_file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  required <- c("file", "cell_id", "ploidy")
  missing <- setdiff(required, names(ploidy_raw))
  if (length(missing) > 0) {
    stop("Cell ploidy file is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ploidy_df <- ploidy_raw %>%
    dplyr::transmute(
      harvest = sub("\\.sps\\.cbs$", "", as_clean_chr(.data$file)),
      barcode_raw = as_clean_chr(.data$cell_id),
      cell_ploidy = safe_num(.data$ploidy),
      ploidy_frac_covered = if ("frac_covered" %in% names(ploidy_raw)) safe_num(.data$frac_covered) else NA_real_
    ) %>%
    dplyr::filter(!is.na(.data$harvest), nzchar(.data$harvest), !is.na(.data$barcode_raw), nzchar(.data$barcode_raw))

  ploidy_df %>%
    dplyr::group_by(.data$harvest, .data$barcode_raw) %>%
    dplyr::summarise(
      cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      ploidy_frac_covered = mean(.data$ploidy_frac_covered, na.rm = TRUE),
      n_ploidy_rows = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cell_ploidy = ifelse(is.nan(.data$cell_ploidy), NA_real_, .data$cell_ploidy),
      ploidy_frac_covered = ifelse(is.nan(.data$ploidy_frac_covered), NA_real_, .data$ploidy_frac_covered)
    )
}

standardize_annotation_ploidy <- function(x) {
  x <- as_clean_chr(x)
  x <- toupper(x)
  x <- gsub("\\s+", "", x)
  x[x %in% c("2", "2.0", "2N")] <- "2N"
  x[x %in% c("4", "4.0", "4N")] <- "4N"
  x[!(x %in% c("2N", "4N"))] <- NA_character_
  x
}

standardize_dose_label <- function(x) {
  out <- standardize_in_vivo_dose(x)
  out <- gsub("mg/kg", "mg/kg", out, fixed = TRUE)
  out
}

extract_sample_from_cell <- function(cell) {
  cell <- as_clean_chr(cell)
  sample <- sub("_.*$", "", cell)
  sample[is.na(cell) | !grepl("_", cell)] <- NA_character_
  sample
}

extract_barcode_from_cell <- function(cell, sample) {
  cell <- as_clean_chr(cell)
  sample <- as_clean_chr(sample)
  barcode <- cell
  has_sample <- !is.na(cell) & !is.na(sample) & nzchar(sample) & startsWith(cell, paste0(sample, "_"))
  barcode[has_sample] <- substring(cell[has_sample], nchar(sample[has_sample]) + 2L)
  fallback <- !has_sample & !is.na(cell) & grepl("_", cell)
  barcode[fallback] <- sub("^[^_]+_", "", cell[fallback])
  barcode
}

cluster_levels <- function(x) {
  x <- as_clean_chr(x)
  ux <- unique(x[!is.na(x) & nzchar(x)])
  nums <- suppressWarnings(as.numeric(ux))
  if (length(ux) == 0) return(character(0))
  if (all(is.finite(nums))) ux[order(nums)] else sort(ux)
}

load_method_metrics <- function(file, method_label, analysis_object, time_col) {
  if (!file.exists(file)) {
    warning("Missing ", method_label, " input for ", analysis_object, ": ", file, call. = FALSE)
    return(data.frame())
  }
  df <- readr::read_csv(
    file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!(time_col %in% names(df))) {
    stop("Input file lacks time column '", time_col, "': ", file, call. = FALSE)
  }
  if (!("cell" %in% names(df))) {
    stop("Input file lacks required cell column: ", file, call. = FALSE)
  }

  ploidy_source <- rep(NA_character_, nrow(df))
  for (candidate in c("Ploidy", "ploidy", "Karyotype")) {
    if (candidate %in% names(df)) {
      missing <- is.na(ploidy_source) | !nzchar(ploidy_source)
      ploidy_source[missing] <- as_clean_chr(df[[candidate]])[missing]
    }
  }

  dose_source <- rep(NA_character_, nrow(df))
  for (candidate in c("Dose", "Dose_DEG", "sample_info_dose")) {
    if (candidate %in% names(df)) {
      missing <- is.na(dose_source) | !nzchar(dose_source)
      dose_source[missing] <- as_clean_chr(df[[candidate]])[missing]
    }
  }

  sample_source <- rep(NA_character_, nrow(df))
  for (candidate in c("sample", "sampleID", "sample_id", "IDs", "orig.ident", "Sequencing.IDs")) {
    if (candidate %in% names(df)) {
      missing <- is.na(sample_source) | !nzchar(sample_source)
      sample_source[missing] <- as_clean_chr(df[[candidate]])[missing]
    }
  }
  missing_sample <- is.na(sample_source) | !nzchar(sample_source)
  sample_source[missing_sample] <- extract_sample_from_cell(df$cell)[missing_sample]

  cluster_source <- rep(NA_character_, nrow(df))
  for (candidate in c("clusters", "cluster", "cluster_final")) {
    if (candidate %in% names(df)) {
      missing <- is.na(cluster_source) | !nzchar(cluster_source)
      cluster_source[missing] <- as_clean_chr(df[[candidate]])[missing]
    }
  }

  out <- df %>%
    dplyr::mutate(
      method = method_label,
      analysis_object = analysis_object,
      analysis_object_safe = safe_file_component(analysis_object),
      source_file = normalizePath(file, mustWork = FALSE),
      pseudotime = safe_num(.data[[time_col]]),
      annotation_ploidy = standardize_annotation_ploidy(ploidy_source),
      Dose = standardize_dose_label(dose_source),
      sample = as_clean_chr(sample_source),
      clusters = as_clean_chr(cluster_source),
      barcode_raw = extract_barcode_from_cell(.data$cell, sample_source)
    )

  object_parts <- strsplit(analysis_object, "/", fixed = TRUE)[[1]]
  out$analysis_tn_scope <- object_parts[1]
  out$analysis_ploidy_scope <- if (length(object_parts) > 1) object_parts[2] else "all"
  if (!("UMAP_1" %in% names(out))) out$UMAP_1 <- NA_real_
  if (!("UMAP_2" %in% names(out))) out$UMAP_2 <- NA_real_
  out$UMAP_1 <- safe_num(out$UMAP_1)
  out$UMAP_2 <- safe_num(out$UMAP_2)
  out
}

metric_manifest <- function(results_root) {
  data.frame(
    method = rep(c("scVelo", "PAGA", "Monocle3"), each = length(analysis_objects)),
    analysis_object = rep(analysis_objects, times = 3L),
    time_col = rep(c("velocity_pseudotime", "dpt_pseudotime", "pseudotime"), each = length(analysis_objects)),
    input_file = c(
      file.path(results_root, "04_trajectory", scenario_id, "scvelo", "01_velocity_groups", analysis_objects, "scvelo_cell_metrics.csv"),
      file.path(results_root, "04_trajectory", scenario_id, "PAGA", "04_paga_groups", analysis_objects, "paga_cell_metrics.csv"),
      file.path(results_root, "04a_psudo_rajectory", monocle3_root_id, "01_pseudotrajectory_groups", analysis_objects, "cells_pseudotime.csv")
    ),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      scenario_id = scenario_id,
      root_cluster = root_label,
      file_exists = file.exists(.data$input_file)
    )
}

load_all_metrics <- function(manifest) {
  rows <- lapply(seq_len(nrow(manifest)), function(i) {
    load_method_metrics(
      file = manifest$input_file[i],
      method_label = manifest$method[i],
      analysis_object = manifest$analysis_object[i],
      time_col = manifest$time_col[i]
    )
  })
  dplyr::bind_rows(rows)
}

join_cell_ploidy <- function(metrics, sample_info_map, ploidy_map) {
  metrics %>%
    dplyr::left_join(sample_info_map, by = "sample") %>%
    dplyr::mutate(
      Dose = dplyr::coalesce(.data$Dose, .data$sample_info_dose)
    ) %>%
    dplyr::left_join(ploidy_map, by = c("harvest", "barcode_raw"))
}

fill_umap_from_reference <- function(metrics) {
  if (!("UMAP_1" %in% names(metrics))) metrics$UMAP_1 <- NA_real_
  if (!("UMAP_2" %in% names(metrics))) metrics$UMAP_2 <- NA_real_
  metrics$UMAP_1 <- safe_num(metrics$UMAP_1)
  metrics$UMAP_2 <- safe_num(metrics$UMAP_2)
  ref <- metrics %>%
    dplyr::filter(!is.na(.data$cell), nzchar(.data$cell), is.finite(.data$UMAP_1), is.finite(.data$UMAP_2)) %>%
    dplyr::select(cell, UMAP_1_ref = UMAP_1, UMAP_2_ref = UMAP_2) %>%
    dplyr::distinct(.data$cell, .keep_all = TRUE)
  if (nrow(ref) == 0) return(metrics)
  metrics %>%
    dplyr::left_join(ref, by = "cell") %>%
    dplyr::mutate(
      UMAP_1 = ifelse(is.finite(.data$UMAP_1), .data$UMAP_1, .data$UMAP_1_ref),
      UMAP_2 = ifelse(is.finite(.data$UMAP_2), .data$UMAP_2, .data$UMAP_2_ref)
    ) %>%
    dplyr::select(-dplyr::any_of(c("UMAP_1_ref", "UMAP_2_ref")))
}

summarise_mapping <- function(metrics) {
  metrics %>%
    dplyr::group_by(.data$method, .data$analysis_object) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_finite_pseudotime = sum(is.finite(.data$pseudotime), na.rm = TRUE),
      n_with_sample = sum(!is.na(.data$sample) & nzchar(.data$sample), na.rm = TRUE),
      n_with_harvest = sum(!is.na(.data$harvest) & nzchar(.data$harvest), na.rm = TRUE),
      n_with_barcode_raw = sum(!is.na(.data$barcode_raw) & nzchar(.data$barcode_raw), na.rm = TRUE),
      n_with_cell_ploidy = sum(is.finite(.data$cell_ploidy), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      cell_ploidy_mapping_rate = ifelse(.data$n_cells > 0, .data$n_with_cell_ploidy / .data$n_cells, NA_real_)
    )
}

choose_numeric_binwidth <- function(x, min_bins = 12L, max_bins = 60L, fallback_bins = 30L) {
  x <- safe_num(x)
  x <- x[is.finite(x)]
  if (length(x) < 2) return(NA_real_)
  xr <- range(x, na.rm = TRUE)
  span <- diff(xr)
  if (!is.finite(span) || span <= 0) return(NA_real_)
  fd <- 2 * stats::IQR(x, na.rm = TRUE) / (length(x)^(1 / 3))
  if (!is.finite(fd) || fd <= 0) fd <- span / fallback_bins
  n_bins <- ceiling(span / fd)
  n_bins <- max(min_bins, min(max_bins, n_bins))
  span / n_bins
}

make_frequency_bins <- function(df, value_col, group_cols, binwidth = NULL, x_range = NULL) {
  x <- safe_num(df[[value_col]])
  keep <- is.finite(x)
  if (!any(keep)) return(data.frame())
  df <- df[keep, , drop = FALSE]
  x <- x[keep]
  if (is.null(x_range)) x_range <- range(x, na.rm = TRUE)
  if (!is.finite(diff(x_range)) || diff(x_range) <= 0) return(data.frame())
  if (is.null(binwidth) || !is.finite(binwidth) || binwidth <= 0) {
    binwidth <- choose_numeric_binwidth(x)
  }
  if (!is.finite(binwidth) || binwidth <= 0) return(data.frame())
  breaks <- seq(x_range[1], x_range[2] + binwidth, by = binwidth)
  if (length(breaks) < 2) breaks <- c(x_range[1], x_range[2])
  df$.value <- x
  df$.bin <- cut(df$.value, breaks = breaks, include.lowest = TRUE, right = FALSE)
  bin_table <- data.frame(
    .bin = levels(df$.bin),
    bin_left = head(breaks, -1),
    bin_right = tail(breaks, -1),
    bin_mid = head(breaks, -1) + diff(breaks) / 2,
    stringsAsFactors = FALSE
  )
  df %>%
    dplyr::filter(!is.na(.data$.bin)) %>%
    dplyr::count(dplyr::across(dplyr::all_of(c(group_cols, ".bin"))), name = "n") %>%
    dplyr::left_join(bin_table, by = ".bin") %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::mutate(frequency = .data$n / sum(.data$n)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(binwidth = binwidth)
}

plot_frequency_histogram <- function(hist_df, group_col, title, x_label) {
  if (nrow(hist_df) == 0) return(NULL)
  ggplot(hist_df, aes(x = .data$bin_mid, y = .data$frequency, fill = .data[[group_col]], color = .data[[group_col]])) +
    geom_col(width = unique(hist_df$binwidth)[1], alpha = 0.35, position = "identity") +
    geom_line(linewidth = 0.7) +
    labs(title = title, x = x_label, y = "Frequency", fill = group_col, color = group_col) +
    theme_extra(base_size = 11)
}

save_sample_frequency_pages <- function(df, group_col, output_pdf, title_prefix, x_label, bins = 60L) {
  output_pdf <- standardize_output_path(output_pdf)
  title_prefix <- standardize_output_text(title_prefix)
  x_label <- standardize_output_text(x_label)
  df <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), !is.na(.data[[group_col]]), nzchar(as.character(.data[[group_col]])), !is.na(.data$sample), nzchar(as.character(.data$sample)))
  if (nrow(df) == 0) {
    return(data.frame(output_pdf = output_pdf, pages = 0L, samples = 0L, stringsAsFactors = FALSE))
  }
  samples <- sort_maybe_numeric(unique(df$sample))
  pages <- split(samples, ceiling(seq_along(samples) / 9L))
  x_range <- range(df$pseudotime, na.rm = TRUE)
  binwidth <- diff(x_range) / bins
  .ensure_dir(dirname(output_pdf))
  grDevices::pdf(output_pdf, width = 12, height = 9)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (page_samples in pages) {
    page_df <- df %>% dplyr::filter(.data$sample %in% page_samples)
    hist_df <- make_frequency_bins(page_df, "pseudotime", c("sample", group_col), binwidth = binwidth, x_range = x_range)
    if (nrow(hist_df) == 0) next
    p <- ggplot(hist_df, aes(x = .data$bin_mid, y = .data$frequency, fill = .data[[group_col]], color = .data[[group_col]])) +
      geom_col(width = binwidth, alpha = 0.35, position = "identity") +
      geom_line(linewidth = 0.55) +
      facet_wrap(~sample, ncol = 3, scales = "fixed") +
      labs(title = title_prefix, x = x_label, y = "Frequency", fill = group_col, color = group_col) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    print(standardize_ggplot_output_text(p))
  }
  data.frame(output_pdf = output_pdf, pages = length(pages), samples = length(samples), stringsAsFactors = FALSE)
}

summary_by_group <- function(df, group_cols, value_col = "pseudotime") {
  df %>%
    dplyr::filter(is.finite(.data[[value_col]])) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.data[[value_col]], na.rm = TRUE),
      median = stats::median(.data[[value_col]], na.rm = TRUE),
      q25 = as.numeric(stats::quantile(.data[[value_col]], 0.25, na.rm = TRUE, names = FALSE)),
      q75 = as.numeric(stats::quantile(.data[[value_col]], 0.75, na.rm = TRUE, names = FALSE)),
      .groups = "drop"
    )
}

group_test <- function(df, value_col, group_col, test_label) {
  dat <- df %>%
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(.data[[group_col]]), nzchar(as.character(.data[[group_col]])))
  groups <- unique(dat[[group_col]])
  if (nrow(dat) == 0 || length(groups) < 2) {
    return(data.frame(
      test = test_label,
      value_col = value_col,
      group_col = group_col,
      statistical_test = NA_character_,
      n = nrow(dat),
      n_groups = length(groups),
      statistic = NA_real_,
      p_value = NA_real_,
      note = "skipped: fewer than two groups",
      stringsAsFactors = FALSE
    ))
  }
  dat$.group <- factor(dat[[group_col]], levels = sort_maybe_numeric(groups))
  if (length(groups) == 2) {
    fit <- tryCatch(stats::wilcox.test(dat[[value_col]] ~ dat$.group), error = function(e) NULL)
    statistical_test <- "Wilcoxon rank-sum"
  } else {
    fit <- tryCatch(stats::kruskal.test(dat[[value_col]] ~ dat$.group), error = function(e) NULL)
    statistical_test <- "Kruskal-Wallis"
  }
  if (is.null(fit)) {
    return(data.frame(
      test = test_label,
      value_col = value_col,
      group_col = group_col,
      statistical_test = statistical_test,
      n = nrow(dat),
      n_groups = length(groups),
      statistic = NA_real_,
      p_value = NA_real_,
      note = "test failed",
      stringsAsFactors = FALSE
    ))
  }
  data.frame(
    test = test_label,
    value_col = value_col,
    group_col = group_col,
    statistical_test = statistical_test,
    n = nrow(dat),
    n_groups = length(groups),
    statistic = unname(fit$statistic),
    p_value = fit$p.value,
    note = NA_character_,
    stringsAsFactors = FALSE
  )
}

pairwise_group_tests <- function(df, value_col, group_col, test_label) {
  dat <- df %>%
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(.data[[group_col]]), nzchar(as.character(.data[[group_col]])))
  groups <- unique(dat[[group_col]])
  if (length(groups) < 3) return(data.frame())
  dat$.group <- factor(dat[[group_col]], levels = sort_maybe_numeric(groups))
  fit <- tryCatch(
    stats::pairwise.wilcox.test(dat[[value_col]], dat$.group, p.adjust.method = "BH"),
    error = function(e) NULL
  )
  if (is.null(fit)) return(data.frame())
  mat <- as.data.frame(as.table(fit$p.value), stringsAsFactors = FALSE)
  names(mat) <- c("group_1", "group_2", "p_adj")
  mat <- mat[is.finite(mat$p_adj), , drop = FALSE]
  mat$test <- test_label
  mat$value_col <- value_col
  mat$group_col <- group_col
  mat
}

tidy_lm_coefficients <- function(model, model_label) {
  if (is.null(model)) return(data.frame())
  coef_df <- as.data.frame(summary(model)$coefficients, stringsAsFactors = FALSE)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  names(coef_df)[seq_len(4)] <- c("estimate", "std_error", "statistic", "p_value")
  coef_df$model <- model_label
  coef_df[, c("model", "term", "estimate", "std_error", "statistic", "p_value"), drop = FALSE]
}

tidy_lm_anova <- function(model, model_label) {
  if (is.null(model)) return(data.frame())
  anova_df <- as.data.frame(stats::anova(model), stringsAsFactors = FALSE)
  anova_df$term <- rownames(anova_df)
  rownames(anova_df) <- NULL
  names(anova_df) <- gsub(" ", "_", names(anova_df), fixed = TRUE)
  anova_df$model <- model_label
  anova_df
}

fit_lm_safe <- function(formula, df) {
  tryCatch(stats::lm(formula, data = df), error = function(e) NULL)
}

cor_test_one <- function(x, y, method) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3 || length(unique(x[keep])) < 2 || length(unique(y[keep])) < 2) {
    return(c(estimate = NA_real_, p_value = NA_real_))
  }
  fit <- tryCatch(stats::cor.test(x[keep], y[keep], method = method, exact = FALSE), error = function(e) NULL)
  if (is.null(fit)) return(c(estimate = NA_real_, p_value = NA_real_))
  c(estimate = unname(fit$estimate), p_value = fit$p.value)
}

cluster_ploidy_correlations <- function(df) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), is.finite(.data$cell_ploidy), !is.na(.data$clusters), nzchar(.data$clusters))
  if (nrow(dat) == 0) return(data.frame())
  rows <- lapply(cluster_levels(dat$clusters), function(cl) {
    sub <- dat[dat$clusters == cl, , drop = FALSE]
    sp <- cor_test_one(sub$cell_ploidy, sub$pseudotime, "spearman")
    pe <- cor_test_one(sub$cell_ploidy, sub$pseudotime, "pearson")
    data.frame(
      clusters = cl,
      n = nrow(sub),
      spearman_rho = sp[["estimate"]],
      spearman_p = sp[["p_value"]],
      pearson_r = pe[["estimate"]],
      pearson_p = pe[["p_value"]],
      stringsAsFactors = FALSE
    )
  })
  out <- dplyr::bind_rows(rows)
  out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
  out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
  out
}

plot_cluster_ploidy_correlation <- function(df, title) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), is.finite(.data$cell_ploidy), !is.na(.data$clusters), nzchar(.data$clusters))
  if (nrow(dat) == 0) return(NULL)
  dat$clusters <- factor(dat$clusters, levels = cluster_levels(dat$clusters))
  ggplot(dat, aes(x = .data$cell_ploidy, y = .data$pseudotime)) +
    geom_point(alpha = 0.25, size = 0.45, color = "#3A4A5B") +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, color = "#C0392B") +
    facet_wrap(~clusters, ncol = 3, scales = "fixed") +
    labs(title = title, x = "Cell ploidy", y = "Pseudotime") +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0))
}

plot_cluster_cell_ploidy_distribution <- function(df, subset_label, title) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$cell_ploidy), !is.na(.data$clusters), nzchar(.data$clusters))
  if (subset_label %in% c("2N", "4N")) {
    dat <- dat %>% dplyr::filter(.data$annotation_ploidy == subset_label)
  }
  if (nrow(dat) == 0) return(list(plot = NULL, hist = data.frame(), binwidth = NA_real_))
  dat$clusters <- factor(dat$clusters, levels = cluster_levels(dat$clusters))
  x_range <- range(dat$cell_ploidy, na.rm = TRUE)
  binwidth <- choose_numeric_binwidth(dat$cell_ploidy, min_bins = 12L, max_bins = 45L)
  hist_df <- make_frequency_bins(dat, "cell_ploidy", c("clusters"), binwidth = binwidth, x_range = x_range)
  density_rows <- lapply(levels(dat$clusters), function(cl) {
    x <- dat$cell_ploidy[dat$clusters == cl]
    x <- x[is.finite(x)]
    if (length(x) < 3 || length(unique(x)) < 2) return(data.frame())
    den <- stats::density(x, from = x_range[1], to = x_range[2], na.rm = TRUE)
    data.frame(clusters = cl, x = den$x, y = den$y * binwidth, stringsAsFactors = FALSE)
  })
  density_df <- dplyr::bind_rows(density_rows)
  p <- ggplot(hist_df, aes(x = .data$bin_mid, y = .data$frequency)) +
    geom_col(width = binwidth, fill = "#4B7895", alpha = 0.55) +
    facet_wrap(~clusters, ncol = 3, scales = "fixed") +
    labs(title = title, x = "Cell ploidy", y = "Frequency") +
    coord_cartesian(xlim = x_range) +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0))
  if (nrow(density_df) > 0) {
    density_df$clusters <- factor(density_df$clusters, levels = levels(dat$clusters))
    p <- p + geom_line(data = density_df, aes(x = .data$x, y = .data$y), inherit.aes = FALSE, color = "#B03A2E", linewidth = 0.65)
  }
  list(plot = p, hist = hist_df, binwidth = binwidth)
}

cell_ploidy_time_correlations <- function(df, group_cols, coverage_label) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$cell_ploidy), is.finite(.data$pseudotime)) %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(group_cols), ~ !is.na(.x) & nzchar(as.character(.x))))
  if (nrow(dat) == 0) return(data.frame())
  grouped <- dat %>% dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)))
  keys <- dplyr::group_keys(grouped)
  splits <- dplyr::group_split(grouped)
  rows <- lapply(seq_along(splits), function(i) {
    sub <- splits[[i]]
    sp <- cor_test_one(sub$cell_ploidy, sub$pseudotime, "spearman")
    pe <- cor_test_one(sub$cell_ploidy, sub$pseudotime, "pearson")
    dplyr::bind_cols(
      keys[i, , drop = FALSE],
      data.frame(
        coverage_filter = coverage_label,
        n = nrow(sub),
        spearman_rho = sp[["estimate"]],
        spearman_p = sp[["p_value"]],
        pearson_r = pe[["estimate"]],
        pearson_p = pe[["p_value"]],
        stringsAsFactors = FALSE
      )
    )
  })
  out <- dplyr::bind_rows(rows)
  out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
  out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
  out
}

plot_cell_ploidy_time_by_dose <- function(df, split_by_ploidy = FALSE, title) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$cell_ploidy), is.finite(.data$pseudotime), !is.na(.data$Dose), nzchar(as.character(.data$Dose)))
  if (isTRUE(split_by_ploidy)) {
    dat <- dat %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N"))
  }
  if (nrow(dat) == 0) return(NULL)
  dat$Dose <- factor(dat$Dose, levels = sort_maybe_numeric(unique(dat$Dose)))
  if (isTRUE(split_by_ploidy)) {
    dat$annotation_ploidy <- factor(dat$annotation_ploidy, levels = c("2N", "4N"))
    ggplot(dat, aes(x = .data$cell_ploidy, y = .data$pseudotime)) +
      geom_point(alpha = 0.25, size = 0.45, color = "#3A4A5B") +
      geom_smooth(method = "lm", se = FALSE, color = "#B03A2E", linewidth = 0.7) +
      facet_grid(annotation_ploidy ~ Dose, scales = "fixed") +
      labs(title = title, x = "Cell ploidy", y = "Pseudotime") +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
  } else {
    ggplot(dat, aes(x = .data$cell_ploidy, y = .data$pseudotime, color = .data$annotation_ploidy)) +
      geom_point(alpha = 0.25, size = 0.45) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "#B03A2E", linewidth = 0.7) +
      facet_wrap(~Dose, nrow = 1, scales = "fixed") +
      scale_color_manual(values = c("2N" = "#2C7FB8", "4N" = "#D95F0E"), na.value = "grey55", name = "Ploidy") +
      labs(title = title, x = "Cell ploidy", y = "Pseudotime") +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
  }
}

save_two_panel_pdf <- function(left_plot, right_plot, output_pdf, width = 11.5, height = 5.6) {
  output_pdf <- standardize_output_path(output_pdf)
  left_plot <- standardize_ggplot_output_text(left_plot)
  right_plot <- standardize_ggplot_output_text(right_plot)
  .ensure_dir(dirname(output_pdf))
  grDevices::pdf(output_pdf, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(left_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(right_plot, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  invisible(output_pdf)
}

plot_tumor_dose_ploidy_umap_pair <- function(df, title, output_pdf) {
  dat <- df %>%
    dplyr::filter(is.finite(.data$UMAP_1), is.finite(.data$UMAP_2))
  if (nrow(dat) == 0) {
    return(data.frame(output_pdf = output_pdf, n_cells = 0L, n_with_cell_ploidy = 0L, n_with_dose = 0L, status = "skipped_no_umap", stringsAsFactors = FALSE))
  }
  dat$Dose <- factor(dat$Dose, levels = sort_maybe_numeric(unique(dat$Dose[!is.na(dat$Dose) & nzchar(dat$Dose)])))
  umap_coord <- coord_square_umap(dat$UMAP_1, dat$UMAP_2)
  base_theme <- theme_extra(base_size = 9) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  dose_df <- dat %>% dplyr::filter(!is.na(.data$Dose), nzchar(as.character(.data$Dose)))
  p_dose <- ggplot() +
    geom_point(data = dat, aes(x = .data$UMAP_1, y = .data$UMAP_2), color = "grey88", size = 0.10, alpha = 0.35) +
    geom_point(data = dose_df, aes(x = .data$UMAP_1, y = .data$UMAP_2, color = .data$Dose), size = 0.16, alpha = 0.75) +
    umap_coord +
    labs(title = "Dose", color = "Dose") +
    base_theme

  ploidy_df <- dat %>% dplyr::filter(is.finite(.data$cell_ploidy))
  p_ploidy <- ggplot() +
    geom_point(data = dat, aes(x = .data$UMAP_1, y = .data$UMAP_2), color = "grey88", size = 0.10, alpha = 0.35) +
    geom_point(data = ploidy_df, aes(x = .data$UMAP_1, y = .data$UMAP_2, color = .data$cell_ploidy), size = 0.16, alpha = 0.78) +
    scale_color_gradientn(
      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
      na.value = "grey88",
      name = "Cell\nploidy"
    ) +
    umap_coord +
    labs(title = "Cell ploidy") +
    base_theme

  save_two_panel_pdf(p_dose, p_ploidy + labs(subtitle = title), output_pdf)
  data.frame(
    output_pdf = output_pdf,
    n_cells = nrow(dat),
    n_with_cell_ploidy = nrow(ploidy_df),
    n_with_dose = nrow(dose_df),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

run_tumor_cell_ploidy_time_added_outputs <- function(df, out_dir, method_label) {
  added_dir <- .ensure_dir(file.path(out_dir, "cell_ploidy_time_by_dose"))
  filters <- list(
    all_cells = df %>% dplyr::filter(is.finite(.data$cell_ploidy)),
    coverage90 = df %>% dplyr::filter(is.finite(.data$cell_ploidy), is.finite(.data$ploidy_frac_covered), .data$ploidy_frac_covered >= 0.9)
  )

  status_rows <- list()
  for (filter_name in names(filters)) {
    dat <- filters[[filter_name]]
    by_dose <- cell_ploidy_time_correlations(dat, "Dose", filter_name)
    by_dose_ploidy <- cell_ploidy_time_correlations(dat %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")), c("Dose", "annotation_ploidy"), filter_name)
    write_csv_safe(by_dose, file.path(added_dir, paste0(filter_name, "_cell_ploidy_time_correlations_by_dose.csv")))
    write_csv_safe(by_dose_ploidy, file.path(added_dir, paste0(filter_name, "_cell_ploidy_time_correlations_by_dose_ploidy.csv")))

    p1 <- plot_cell_ploidy_time_by_dose(
      dat,
      split_by_ploidy = FALSE,
      title = paste(method_label, "Tumor cell ploidy vs pseudotime by Dose:", filter_name)
    )
    p2 <- plot_cell_ploidy_time_by_dose(
      dat,
      split_by_ploidy = TRUE,
      title = paste(method_label, "Tumor cell ploidy vs pseudotime by Dose and 2N/4N:", filter_name)
    )
    scatter_by_dose <- file.path(added_dir, paste0(filter_name, "_cell_ploidy_time_by_dose_scatter.pdf"))
    scatter_by_dose_ploidy <- file.path(added_dir, paste0(filter_name, "_cell_ploidy_time_by_dose_ploidy_scatter.pdf"))
    if (!is.null(p1)) save_pdf(p1, scatter_by_dose, width = 11, height = 4.8)
    if (!is.null(p2)) save_pdf(p2, scatter_by_dose_ploidy, width = 11, height = 6.8)

    umap_status <- plot_tumor_dose_ploidy_umap_pair(
      dat,
      title = paste(method_label, filter_name),
      output_pdf = file.path(added_dir, paste0(filter_name, "_tumor_umap_dose_cell_ploidy.pdf"))
    )
    status_rows[[filter_name]] <- data.frame(
      coverage_filter = filter_name,
      n_cells = nrow(dat),
      n_by_dose_tests = nrow(by_dose),
      n_by_dose_ploidy_tests = nrow(by_dose_ploidy),
      scatter_by_dose_pdf = if (!is.null(p1)) scatter_by_dose else NA_character_,
      scatter_by_dose_ploidy_pdf = if (!is.null(p2)) scatter_by_dose_ploidy else NA_character_,
      umap_pdf = umap_status$output_pdf,
      umap_status = umap_status$status,
      stringsAsFactors = FALSE
    )
  }
  status <- dplyr::bind_rows(status_rows)
  write_csv_safe(status, file.path(added_dir, "cell_ploidy_time_by_dose_added_outputs_status.csv"))
  status
}

method_order <- c("scVelo", "PAGA", "Monocle3")
method_colors <- c("scVelo" = "#1B9E77", "PAGA" = "#7570B3", "Monocle3" = "#D95F02")

filter_complete_group_values <- function(df, cols) {
  if (length(cols) == 0) return(df)
  df %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(cols), ~ !is.na(.x) & nzchar(as.character(.x))))
}

paired_method_tests <- function(df, split_cols, test_label) {
  dat <- df %>%
    dplyr::filter(.data$method %in% method_order, is.finite(.data$pseudotime), !is.na(.data$cell), nzchar(.data$cell)) %>%
    filter_complete_group_values(split_cols)
  if (nrow(dat) == 0) {
    return(list(global = data.frame(), pairwise = data.frame()))
  }

  grouped <- dat %>% dplyr::group_by(dplyr::across(dplyr::all_of(split_cols)))
  keys <- dplyr::group_keys(grouped)
  splits <- dplyr::group_split(grouped)

  global_rows <- list()
  pairwise_rows <- list()
  for (i in seq_along(splits)) {
    sub <- splits[[i]]
    wide <- sub %>%
      dplyr::group_by(.data$cell, .data$method) %>%
      dplyr::summarise(pseudotime = mean(.data$pseudotime, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = method, values_from = pseudotime)
    for (method in method_order) {
      if (!(method %in% names(wide))) wide[[method]] <- NA_real_
    }
    complete <- wide %>% dplyr::filter(dplyr::if_all(dplyr::all_of(method_order), is.finite))
    n_complete <- nrow(complete)
    key <- keys[i, , drop = FALSE]
    if (n_complete >= 3) {
      mat <- as.matrix(complete[, method_order, drop = FALSE])
      fit <- tryCatch(stats::friedman.test(mat), error = function(e) NULL)
      global_rows[[i]] <- dplyr::bind_cols(
        key,
        data.frame(
          test = test_label,
          method_test = "Friedman paired by cell",
          n_complete_cells = n_complete,
          statistic = if (!is.null(fit)) unname(fit$statistic) else NA_real_,
          p_value = if (!is.null(fit)) fit$p.value else NA_real_,
          status = if (!is.null(fit)) "ok" else "failed",
          stringsAsFactors = FALSE
        )
      )

      pairs <- utils::combn(method_order, 2, simplify = FALSE)
      pairwise_rows[[i]] <- dplyr::bind_rows(lapply(pairs, function(pair) {
        wt <- tryCatch(stats::wilcox.test(complete[[pair[1]]], complete[[pair[2]]], paired = TRUE, exact = FALSE), error = function(e) NULL)
        dplyr::bind_cols(
          key,
          data.frame(
            test = test_label,
            method_1 = pair[1],
            method_2 = pair[2],
            n_complete_cells = n_complete,
            statistic = if (!is.null(wt)) unname(wt$statistic) else NA_real_,
            p_value = if (!is.null(wt)) wt$p.value else NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }))
    } else {
      global_rows[[i]] <- dplyr::bind_cols(
        key,
        data.frame(
          test = test_label,
          method_test = "Friedman paired by cell",
          n_complete_cells = n_complete,
          statistic = NA_real_,
          p_value = NA_real_,
          status = "skipped_fewer_than_three_complete_cells",
          stringsAsFactors = FALSE
        )
      )
    }
  }

  global <- dplyr::bind_rows(global_rows)
  pairwise <- dplyr::bind_rows(pairwise_rows)
  if (nrow(global) > 0) global$p_adj <- stats::p.adjust(global$p_value, method = "BH")
  if (nrow(pairwise) > 0) pairwise$p_adj <- stats::p.adjust(pairwise$p_value, method = "BH")
  list(global = global, pairwise = pairwise)
}

pairwise_method_correlations <- function(df, split_cols, test_label) {
  dat <- df %>%
    dplyr::filter(.data$method %in% method_order, is.finite(.data$pseudotime), !is.na(.data$cell), nzchar(as.character(.data$cell))) %>%
    filter_complete_group_values(split_cols)
  if (nrow(dat) == 0) return(data.frame())
  grouped <- dat %>% dplyr::group_by(dplyr::across(dplyr::all_of(split_cols)))
  keys <- dplyr::group_keys(grouped)
  splits <- dplyr::group_split(grouped)
  pairs <- utils::combn(method_order, 2, simplify = FALSE)
  rows <- list()
  idx <- 1L
  for (i in seq_along(splits)) {
    sub <- splits[[i]]
    wide <- sub %>%
      dplyr::group_by(.data$cell, .data$method) %>%
      dplyr::summarise(pseudotime = mean(.data$pseudotime, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = method, values_from = pseudotime)
    key <- keys[i, , drop = FALSE]
    for (pair in pairs) {
      if (!all(pair %in% names(wide))) next
      keep <- is.finite(wide[[pair[1]]]) & is.finite(wide[[pair[2]]])
      n <- sum(keep)
      sp <- if (n >= 3) cor_test_one(wide[[pair[1]]][keep], wide[[pair[2]]][keep], "spearman") else c(estimate = NA_real_, p_value = NA_real_)
      pe <- if (n >= 3) cor_test_one(wide[[pair[1]]][keep], wide[[pair[2]]][keep], "pearson") else c(estimate = NA_real_, p_value = NA_real_)
      rows[[idx]] <- dplyr::bind_cols(
        key,
        data.frame(
          test = test_label,
          method_1 = pair[1],
          method_2 = pair[2],
          n_complete_cells = n,
          spearman_rho = sp[["estimate"]],
          spearman_p = sp[["p_value"]],
          pearson_r = pe[["estimate"]],
          pearson_p = pe[["p_value"]],
          stringsAsFactors = FALSE
        )
      )
      idx <- idx + 1L
    }
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) > 0) {
    out$spearman_p_adj <- stats::p.adjust(out$spearman_p, method = "BH")
    out$pearson_p_adj <- stats::p.adjust(out$pearson_p, method = "BH")
  }
  out
}

compare_cluster_order_pair <- function(x, method_1, method_2) {
  wide <- x %>%
    dplyr::select(clusters, method, median_pseudotime) %>%
    tidyr::pivot_wider(names_from = method, values_from = median_pseudotime)
  if (!all(c(method_1, method_2) %in% names(wide))) {
    return(data.frame(n_clusters = 0L, n_cluster_pairs = 0L, concordant_pairs = 0L, discordant_pairs = 0L, tied_pairs = 0L, concordance_rate = NA_real_))
  }
  wide <- wide %>% dplyr::filter(is.finite(.data[[method_1]]), is.finite(.data[[method_2]]))
  if (nrow(wide) < 2) {
    return(data.frame(n_clusters = nrow(wide), n_cluster_pairs = 0L, concordant_pairs = 0L, discordant_pairs = 0L, tied_pairs = 0L, concordance_rate = NA_real_))
  }
  pairs <- utils::combn(seq_len(nrow(wide)), 2, simplify = FALSE)
  signs <- lapply(pairs, function(pair_idx) {
    s1 <- sign(wide[[method_1]][pair_idx[1]] - wide[[method_1]][pair_idx[2]])
    s2 <- sign(wide[[method_2]][pair_idx[1]] - wide[[method_2]][pair_idx[2]])
    c(s1 = s1, s2 = s2)
  })
  signs <- do.call(rbind, signs)
  tied <- signs[, "s1"] == 0 | signs[, "s2"] == 0
  comparable <- !tied
  concordant <- comparable & signs[, "s1"] == signs[, "s2"]
  discordant <- comparable & signs[, "s1"] != signs[, "s2"]
  data.frame(
    n_clusters = nrow(wide),
    n_cluster_pairs = length(pairs),
    concordant_pairs = sum(concordant),
    discordant_pairs = sum(discordant),
    tied_pairs = sum(tied),
    concordance_rate = ifelse(sum(comparable) > 0, sum(concordant) / sum(comparable), NA_real_)
  )
}

compare_root_relative_pair <- function(x, method_1, method_2, root_cluster = "6") {
  wide <- x %>%
    dplyr::select(clusters, method, median_pseudotime) %>%
    tidyr::pivot_wider(names_from = method, values_from = median_pseudotime)
  if (!all(c(method_1, method_2) %in% names(wide)) || !(root_cluster %in% wide$clusters)) {
    return(data.frame(n_compared_clusters = 0L, concordant_clusters = 0L, discordant_clusters = 0L, tied_clusters = 0L, root_relative_concordance_rate = NA_real_))
  }
  root_row <- wide[wide$clusters == root_cluster, , drop = FALSE]
  comp <- wide %>%
    dplyr::filter(.data$clusters != root_cluster, is.finite(.data[[method_1]]), is.finite(.data[[method_2]]))
  if (nrow(comp) == 0 || !is.finite(root_row[[method_1]][1]) || !is.finite(root_row[[method_2]][1])) {
    return(data.frame(n_compared_clusters = nrow(comp), concordant_clusters = 0L, discordant_clusters = 0L, tied_clusters = 0L, root_relative_concordance_rate = NA_real_))
  }
  s1 <- sign(comp[[method_1]] - root_row[[method_1]][1])
  s2 <- sign(comp[[method_2]] - root_row[[method_2]][1])
  tied <- s1 == 0 | s2 == 0
  comparable <- !tied
  concordant <- comparable & s1 == s2
  discordant <- comparable & s1 != s2
  data.frame(
    n_compared_clusters = nrow(comp),
    concordant_clusters = sum(concordant),
    discordant_clusters = sum(discordant),
    tied_clusters = sum(tied),
    root_relative_concordance_rate = ifelse(sum(comparable) > 0, sum(concordant) / sum(comparable), NA_real_)
  )
}

run_trajectory_direction_consistency_outputs <- function(dat, out_dir, root_cluster = "6") {
  direction_dir <- .ensure_dir(file.path(out_dir, "trajectory_direction_consistency"))
  trend_df <- dat %>%
    dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all", .data$method %in% method_order, is.finite(.data$pseudotime), !is.na(.data$clusters), nzchar(as.character(.data$clusters))) %>%
    dplyr::group_by(.data$analysis_object, .data$method, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$analysis_object, .data$method) %>%
    dplyr::mutate(
      cluster_rank = rank(.data$median_pseudotime, ties.method = "average"),
      n_clusters = dplyr::n(),
      rank_percentile = ifelse(.data$n_clusters > 1, (.data$cluster_rank - 1) / (.data$n_clusters - 1), NA_real_)
    ) %>%
    dplyr::ungroup()
  write_csv_safe(trend_df, file.path(direction_dir, "cluster_median_pseudotime_rank_by_method.csv"))

  pairs <- utils::combn(method_order, 2, simplify = FALSE)
  cluster_cor_rows <- list()
  order_rows <- list()
  root_rows <- list()
  idx <- 1L
  for (obj in unique(trend_df$analysis_object)) {
    obj_df <- trend_df %>% dplyr::filter(.data$analysis_object == obj)
    wide <- obj_df %>%
      dplyr::select(clusters, method, median_pseudotime) %>%
      tidyr::pivot_wider(names_from = method, values_from = median_pseudotime)
    for (pair in pairs) {
      if (!all(pair %in% names(wide))) next
      keep <- is.finite(wide[[pair[1]]]) & is.finite(wide[[pair[2]]])
      n <- sum(keep)
      sp <- if (n >= 3) cor_test_one(wide[[pair[1]]][keep], wide[[pair[2]]][keep], "spearman") else c(estimate = NA_real_, p_value = NA_real_)
      kd <- if (n >= 3) cor_test_one(wide[[pair[1]]][keep], wide[[pair[2]]][keep], "kendall") else c(estimate = NA_real_, p_value = NA_real_)
      cluster_cor_rows[[idx]] <- data.frame(
        analysis_object = obj,
        method_1 = pair[1],
        method_2 = pair[2],
        n_common_clusters = n,
        spearman_rho = sp[["estimate"]],
        spearman_p = sp[["p_value"]],
        kendall_tau = kd[["estimate"]],
        kendall_p = kd[["p_value"]],
        stringsAsFactors = FALSE
      )
      order_stats <- compare_cluster_order_pair(obj_df, pair[1], pair[2])
      order_rows[[idx]] <- cbind(
        data.frame(analysis_object = obj, method_1 = pair[1], method_2 = pair[2], stringsAsFactors = FALSE),
        order_stats
      )
      root_stats <- compare_root_relative_pair(obj_df, pair[1], pair[2], root_cluster = root_cluster)
      root_rows[[idx]] <- cbind(
        data.frame(analysis_object = obj, method_1 = pair[1], method_2 = pair[2], root_cluster = root_cluster, stringsAsFactors = FALSE),
        root_stats
      )
      idx <- idx + 1L
    }
  }
  cluster_cor <- dplyr::bind_rows(cluster_cor_rows)
  if (nrow(cluster_cor) > 0) {
    cluster_cor$spearman_p_adj <- stats::p.adjust(cluster_cor$spearman_p, method = "BH")
    cluster_cor$kendall_p_adj <- stats::p.adjust(cluster_cor$kendall_p, method = "BH")
  }
  order_concordance <- dplyr::bind_rows(order_rows)
  root_concordance <- dplyr::bind_rows(root_rows)
  write_csv_safe(cluster_cor, file.path(direction_dir, "cluster_median_pairwise_method_correlations.csv"))
  write_csv_safe(order_concordance, file.path(direction_dir, "cluster_pair_order_concordance_by_method_pair.csv"))
  write_csv_safe(root_concordance, file.path(direction_dir, "root6_relative_cluster_direction_concordance.csv"))

  root_summary <- trend_df %>%
    dplyr::filter(.data$clusters == root_cluster) %>%
    dplyr::select(analysis_object, method, clusters, n_cells, median_pseudotime, cluster_rank, n_clusters, rank_percentile)
  write_csv_safe(root_summary, file.path(direction_dir, "root6_cluster_rank_by_method.csv"))

  cell_cor_object <- pairwise_method_correlations(dat %>% dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all"), "analysis_object", "cell_level_method_pseudotime_correlation_by_object")
  cell_cor_ploidy <- pairwise_method_correlations(dat %>% dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all", .data$annotation_ploidy %in% c("2N", "4N")), c("analysis_object", "annotation_ploidy"), "cell_level_method_pseudotime_correlation_by_object_ploidy")
  cell_cor_dose <- pairwise_method_correlations(dat %>% dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all", !is.na(.data$Dose), nzchar(as.character(.data$Dose))), c("analysis_object", "Dose"), "cell_level_method_pseudotime_correlation_by_object_dose")
  cell_cor_ploidy_dose <- pairwise_method_correlations(dat %>% dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all", .data$annotation_ploidy %in% c("2N", "4N"), !is.na(.data$Dose), nzchar(as.character(.data$Dose))), c("analysis_object", "annotation_ploidy", "Dose"), "cell_level_method_pseudotime_correlation_by_object_ploidy_dose")
  write_csv_safe(cell_cor_object, file.path(direction_dir, "cell_level_pairwise_method_correlations_by_object.csv"))
  write_csv_safe(cell_cor_ploidy, file.path(direction_dir, "cell_level_pairwise_method_correlations_by_object_ploidy.csv"))
  write_csv_safe(cell_cor_dose, file.path(direction_dir, "cell_level_pairwise_method_correlations_by_object_dose.csv"))
  write_csv_safe(cell_cor_ploidy_dose, file.path(direction_dir, "cell_level_pairwise_method_correlations_by_object_ploidy_dose.csv"))

  plot_df <- trend_df
  plot_df$method <- factor(plot_df$method, levels = method_order)
  plot_df$clusters <- factor(plot_df$clusters, levels = cluster_levels(plot_df$clusters))
  p_rank <- ggplot(plot_df, aes(x = .data$method, y = .data$rank_percentile, group = .data$clusters, color = .data$clusters)) +
    geom_line(alpha = 0.75, linewidth = 0.55) +
    geom_point(size = 1.7, alpha = 0.9) +
    facet_wrap(~analysis_object, nrow = 1) +
    labs(
      title = paste0("Cluster pseudotime rank consistency across methods, root cluster ", root_cluster),
      x = "Method",
      y = "Cluster median pseudotime rank percentile\n(0 = earliest, 1 = latest)",
      color = "Cluster"
    ) +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0))
  save_pdf(p_rank, file.path(direction_dir, "cluster_pseudotime_rank_percentile_by_method.pdf"), width = 11, height = 5.8)

  data.frame(
    n_cluster_rank_rows = nrow(trend_df),
    n_cluster_correlation_rows = nrow(cluster_cor),
    n_order_concordance_rows = nrow(order_concordance),
    n_root_concordance_rows = nrow(root_concordance),
    n_cell_correlation_rows = nrow(cell_cor_object) + nrow(cell_cor_ploidy) + nrow(cell_cor_dose) + nrow(cell_cor_ploidy_dose),
    output_dir = direction_dir,
    stringsAsFactors = FALSE
  )
}

plot_method_violin_box <- function(df, x_col, title, x_label = NULL, facet_formula = NULL) {
  dat <- df %>%
    dplyr::filter(.data$method %in% method_order, is.finite(.data$pseudotime)) %>%
    filter_complete_group_values(x_col)
  if (nrow(dat) == 0) return(NULL)
  dat$method <- factor(dat$method, levels = method_order)
  dat[[x_col]] <- factor(dat[[x_col]], levels = sort_maybe_numeric(unique(dat[[x_col]])))
  p <- ggplot(dat, aes(x = .data[[x_col]], y = .data$pseudotime, fill = .data$method)) +
    geom_violin(position = position_dodge(width = 0.82), scale = "width", trim = TRUE, alpha = 0.55, linewidth = 0.2) +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.size = 0.15, alpha = 0.82, linewidth = 0.2) +
    scale_fill_manual(values = method_colors, drop = FALSE, name = "Method") +
    labs(title = title, x = x_label %||% x_col, y = "Pseudotime") +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(facet_formula)) p <- p + facet_wrap(facet_formula, scales = "free_x")
  p
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || !nzchar(x)) y else x
}

run_cross_method_comparison_outputs <- function(metrics, output_base) {
  out_dir <- .ensure_dir(file.path(output_base, "06_cross_method_comparison"))
  dat <- metrics %>%
    dplyr::filter(.data$analysis_ploidy_scope == "ploidy_all", .data$method %in% method_order, is.finite(.data$pseudotime)) %>%
    dplyr::mutate(
      method = factor(.data$method, levels = method_order),
      annotation_ploidy = factor(.data$annotation_ploidy, levels = c("2N", "4N")),
      Dose = factor(.data$Dose, levels = sort_maybe_numeric(unique(.data$Dose[!is.na(.data$Dose) & nzchar(as.character(.data$Dose))])))
    )
  if (nrow(dat) == 0) {
    warning("No ploidy_all metrics available for cross-method comparison.", call. = FALSE)
    return(data.frame())
  }

  object_status <- list()
  for (obj in unique(dat$analysis_object)) {
    obj_df <- dat %>% dplyr::filter(.data$analysis_object == obj)
    obj_dir <- .ensure_dir(file.path(out_dir, safe_file_component(obj)))

    cluster_plot <- plot_method_violin_box(
      obj_df,
      "clusters",
      paste(obj, "cluster pseudotime by method"),
      x_label = "Cluster"
    )
    ploidy_plot <- plot_method_violin_box(
      obj_df %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")),
      "annotation_ploidy",
      paste(obj, "2N/4N pseudotime by method"),
      x_label = "Ploidy"
    )
    dose_plot <- plot_method_violin_box(
      obj_df %>% dplyr::filter(!is.na(.data$Dose), nzchar(as.character(.data$Dose))),
      "Dose",
      paste(obj, "Dose pseudotime by method"),
      x_label = "Dose"
    )
    dose_ploidy_plot <- plot_method_violin_box(
      obj_df %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N"), !is.na(.data$Dose), nzchar(as.character(.data$Dose))),
      "Dose",
      paste(obj, "Dose pseudotime by method within 2N/4N"),
      x_label = "Dose",
      facet_formula = stats::as.formula("~annotation_ploidy")
    )
    if (!is.null(cluster_plot)) save_pdf(cluster_plot, file.path(obj_dir, "cluster_pseudotime_violin_box_by_method.pdf"), width = 12, height = 6)
    if (!is.null(ploidy_plot)) save_pdf(ploidy_plot, file.path(obj_dir, "ploidy_pseudotime_violin_box_by_method.pdf"), width = 7, height = 5.5)
    if (!is.null(dose_plot)) save_pdf(dose_plot, file.path(obj_dir, "dose_pseudotime_violin_box_by_method.pdf"), width = 8, height = 5.5)
    if (!is.null(dose_ploidy_plot)) save_pdf(dose_ploidy_plot, file.path(obj_dir, "dose_within_ploidy_pseudotime_violin_box_by_method.pdf"), width = 11, height = 5.5)
    object_status[[obj]] <- data.frame(
      analysis_object = obj,
      n_rows = nrow(obj_df),
      n_methods = dplyr::n_distinct(obj_df$method),
      cluster_plot = !is.null(cluster_plot),
      ploidy_plot = !is.null(ploidy_plot),
      dose_plot = !is.null(dose_plot),
      dose_within_ploidy_plot = !is.null(dose_ploidy_plot),
      stringsAsFactors = FALSE
    )
  }

  cluster_tests <- paired_method_tests(dat, c("analysis_object", "clusters"), "method_difference_within_cluster")
  ploidy_tests <- paired_method_tests(dat %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")), c("analysis_object", "annotation_ploidy"), "method_difference_within_ploidy")
  dose_tests <- paired_method_tests(dat %>% dplyr::filter(!is.na(.data$Dose), nzchar(as.character(.data$Dose))), c("analysis_object", "Dose"), "method_difference_within_dose")
  dose_ploidy_tests <- paired_method_tests(dat %>% dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N"), !is.na(.data$Dose), nzchar(as.character(.data$Dose))), c("analysis_object", "annotation_ploidy", "Dose"), "method_difference_within_ploidy_dose")

  write_csv_safe(cluster_tests$global, file.path(out_dir, "method_difference_within_cluster_friedman_tests.csv"))
  write_csv_safe(cluster_tests$pairwise, file.path(out_dir, "method_difference_within_cluster_pairwise_paired_wilcox.csv"))
  write_csv_safe(ploidy_tests$global, file.path(out_dir, "method_difference_within_ploidy_friedman_tests.csv"))
  write_csv_safe(ploidy_tests$pairwise, file.path(out_dir, "method_difference_within_ploidy_pairwise_paired_wilcox.csv"))
  write_csv_safe(dose_tests$global, file.path(out_dir, "method_difference_within_dose_friedman_tests.csv"))
  write_csv_safe(dose_tests$pairwise, file.path(out_dir, "method_difference_within_dose_pairwise_paired_wilcox.csv"))
  write_csv_safe(dose_ploidy_tests$global, file.path(out_dir, "method_difference_within_ploidy_dose_friedman_tests.csv"))
  write_csv_safe(dose_ploidy_tests$pairwise, file.path(out_dir, "method_difference_within_ploidy_dose_pairwise_paired_wilcox.csv"))
  direction_status <- run_trajectory_direction_consistency_outputs(dat, out_dir, root_cluster = root_label)

  status <- dplyr::bind_rows(object_status)
  status$direction_consistency_output_dir <- direction_status$output_dir[1]
  write_csv_safe(status, file.path(out_dir, "cross_method_comparison_plot_status.csv"))
  status
}

dose_effect_by_group_tests <- function(df, value_col, group_cols, test_label) {
  dat <- df %>%
    dplyr::filter(is.finite(.data[[value_col]]), !is.na(.data$Dose), nzchar(as.character(.data$Dose))) %>%
    filter_complete_group_values(group_cols)
  if (nrow(dat) == 0) return(list(global = data.frame(), pairwise = data.frame()))
  grouped <- dat %>% dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)))
  keys <- dplyr::group_keys(grouped)
  splits <- dplyr::group_split(grouped)
  global_rows <- list()
  pairwise_rows <- list()
  for (i in seq_along(splits)) {
    sub <- splits[[i]]
    sub$Dose <- factor(sub$Dose, levels = sort_maybe_numeric(unique(sub$Dose)))
    key <- keys[i, , drop = FALSE]
    global_rows[[i]] <- dplyr::bind_cols(key, group_test(sub, value_col, "Dose", test_label))
    pairwise <- pairwise_group_tests(sub, value_col, "Dose", paste0(test_label, "_pairwise"))
    if (nrow(pairwise) > 0) pairwise_rows[[i]] <- dplyr::bind_cols(key, pairwise)
  }
  global <- dplyr::bind_rows(global_rows)
  pairwise <- dplyr::bind_rows(pairwise_rows)
  if (nrow(global) > 0) global$p_adj <- stats::p.adjust(global$p_value, method = "BH")
  if (nrow(pairwise) > 0 && "p_adj" %in% names(pairwise)) pairwise$p_adj_global <- stats::p.adjust(pairwise$p_adj, method = "BH")
  list(global = global, pairwise = pairwise)
}

run_tumor_ploidy_dose_response_outputs <- function(metrics, output_base) {
  out_dir <- .ensure_dir(file.path(output_base, "07_tumor_ploidy_dose_response"))
  dat <- metrics %>%
    dplyr::filter(
      .data$analysis_object == "Tumor/ploidy_all",
      .data$analysis_ploidy_scope == "ploidy_all",
      .data$method %in% method_order,
      .data$annotation_ploidy %in% c("2N", "4N"),
      !is.na(.data$Dose),
      nzchar(as.character(.data$Dose)),
      is.finite(.data$pseudotime)
    ) %>%
    dplyr::mutate(
      method = factor(.data$method, levels = method_order),
      annotation_ploidy = factor(.data$annotation_ploidy, levels = c("2N", "4N")),
      Dose = factor(.data$Dose, levels = sort_maybe_numeric(unique(.data$Dose)))
    )
  if (nrow(dat) == 0) {
    warning("No Tumor/ploidy_all data available for ploidy-dose response analysis.", call. = FALSE)
    return(data.frame())
  }

  cell_summary <- summary_by_group(dat, c("method", "annotation_ploidy", "Dose"), "pseudotime")
  sample_summary <- dat %>%
    dplyr::filter(!is.na(.data$sample), nzchar(.data$sample)) %>%
    dplyr::group_by(.data$method, .data$annotation_ploidy, .data$Dose, .data$sample) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      sample_median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      sample_mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(cell_summary, file.path(out_dir, "tumor_ploidy_dose_cell_pseudotime_summary.csv"))
  write_csv_safe(sample_summary, file.path(out_dir, "tumor_ploidy_dose_sample_median_pseudotime_summary.csv"))

  cell_dose_tests <- dose_effect_by_group_tests(dat, "pseudotime", c("method", "annotation_ploidy"), "dose_effect_within_method_ploidy_cell_level")
  sample_dose_tests <- dose_effect_by_group_tests(sample_summary, "sample_median_pseudotime", c("method", "annotation_ploidy"), "dose_effect_within_method_ploidy_sample_level")
  write_csv_safe(cell_dose_tests$global, file.path(out_dir, "dose_effect_within_method_ploidy_cell_level_tests.csv"))
  write_csv_safe(cell_dose_tests$pairwise, file.path(out_dir, "dose_effect_within_method_ploidy_cell_level_pairwise.csv"))
  write_csv_safe(sample_dose_tests$global, file.path(out_dir, "dose_effect_within_method_ploidy_sample_level_tests.csv"))
  write_csv_safe(sample_dose_tests$pairwise, file.path(out_dir, "dose_effect_within_method_ploidy_sample_level_pairwise.csv"))

  model_rows <- list()
  anova_rows <- list()
  sample_model_rows <- list()
  sample_anova_rows <- list()
  for (method in method_order) {
    sub <- dat %>% dplyr::filter(.data$method == .env$method)
    model <- if (nrow(sub) > 0 && dplyr::n_distinct(sub$annotation_ploidy) > 1 && dplyr::n_distinct(sub$Dose) > 1) {
      fit_lm_safe(pseudotime ~ annotation_ploidy * Dose, sub)
    } else {
      NULL
    }
    model_rows[[method]] <- tidy_lm_coefficients(model, paste0(method, "_cell_pseudotime_ploidy_dose"))
    anova_rows[[method]] <- tidy_lm_anova(model, paste0(method, "_cell_pseudotime_ploidy_dose"))

    sample_sub <- sample_summary %>% dplyr::filter(.data$method == .env$method)
    sample_model <- if (nrow(sample_sub) > 0 && dplyr::n_distinct(sample_sub$annotation_ploidy) > 1 && dplyr::n_distinct(sample_sub$Dose) > 1) {
      fit_lm_safe(sample_median_pseudotime ~ annotation_ploidy * Dose, sample_sub)
    } else {
      NULL
    }
    sample_model_rows[[method]] <- tidy_lm_coefficients(sample_model, paste0(method, "_sample_median_pseudotime_ploidy_dose"))
    sample_anova_rows[[method]] <- tidy_lm_anova(sample_model, paste0(method, "_sample_median_pseudotime_ploidy_dose"))
  }
  write_csv_safe(dplyr::bind_rows(model_rows), file.path(out_dir, "cell_level_ploidy_dose_lm_coefficients_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(anova_rows), file.path(out_dir, "cell_level_ploidy_dose_lm_anova_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(sample_model_rows), file.path(out_dir, "sample_level_ploidy_dose_lm_coefficients_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(sample_anova_rows), file.path(out_dir, "sample_level_ploidy_dose_lm_anova_by_method.csv"))

  continuous_dat <- dat %>%
    dplyr::filter(is.finite(.data$cell_ploidy))
  continuous_cell_summary <- continuous_dat %>%
    dplyr::group_by(.data$method, .data$Dose) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      .groups = "drop"
    )
  continuous_cell_ploidy_summary <- continuous_dat %>%
    dplyr::group_by(.data$method, .data$annotation_ploidy, .data$Dose) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      .groups = "drop"
    )
  continuous_sample_summary <- continuous_dat %>%
    dplyr::filter(!is.na(.data$sample), nzchar(as.character(.data$sample))) %>%
    dplyr::group_by(.data$method, .data$Dose, .data$sample) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      sample_median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      sample_mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      sample_median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      .groups = "drop"
    )
  continuous_corr_by_dose <- cell_ploidy_time_correlations(continuous_dat, c("method", "Dose"), "all_cells")
  continuous_corr_by_dose_ploidy <- cell_ploidy_time_correlations(continuous_dat, c("method", "annotation_ploidy", "Dose"), "all_cells")
  write_csv_safe(continuous_cell_summary, file.path(out_dir, "continuous_cell_ploidy_dose_summary_by_method.csv"))
  write_csv_safe(continuous_cell_ploidy_summary, file.path(out_dir, "continuous_cell_ploidy_dose_summary_by_method_ploidy.csv"))
  write_csv_safe(continuous_sample_summary, file.path(out_dir, "continuous_cell_ploidy_dose_sample_summary_by_method.csv"))
  write_csv_safe(continuous_corr_by_dose, file.path(out_dir, "continuous_cell_ploidy_time_correlations_by_method_dose.csv"))
  write_csv_safe(continuous_corr_by_dose_ploidy, file.path(out_dir, "continuous_cell_ploidy_time_correlations_by_method_ploidy_dose.csv"))

  continuous_cell_model_rows <- list()
  continuous_cell_anova_rows <- list()
  continuous_sample_model_rows <- list()
  continuous_sample_anova_rows <- list()
  for (method in method_order) {
    sub <- continuous_dat %>% dplyr::filter(.data$method == .env$method)
    continuous_model <- if (nrow(sub) > 0 && dplyr::n_distinct(sub$Dose) > 1 && dplyr::n_distinct(sub$cell_ploidy) > 1) {
      fit_lm_safe(pseudotime ~ cell_ploidy * Dose, sub)
    } else {
      NULL
    }
    continuous_cell_model_rows[[method]] <- tidy_lm_coefficients(continuous_model, paste0(method, "_cell_pseudotime_continuous_cell_ploidy_dose"))
    continuous_cell_anova_rows[[method]] <- tidy_lm_anova(continuous_model, paste0(method, "_cell_pseudotime_continuous_cell_ploidy_dose"))

    sample_sub <- continuous_sample_summary %>% dplyr::filter(.data$method == .env$method)
    continuous_sample_model <- if (nrow(sample_sub) > 0 && dplyr::n_distinct(sample_sub$Dose) > 1 && dplyr::n_distinct(sample_sub$sample_mean_cell_ploidy) > 1) {
      fit_lm_safe(sample_median_pseudotime ~ sample_mean_cell_ploidy * Dose, sample_sub)
    } else {
      NULL
    }
    continuous_sample_model_rows[[method]] <- tidy_lm_coefficients(continuous_sample_model, paste0(method, "_sample_median_pseudotime_mean_cell_ploidy_dose"))
    continuous_sample_anova_rows[[method]] <- tidy_lm_anova(continuous_sample_model, paste0(method, "_sample_median_pseudotime_mean_cell_ploidy_dose"))
  }
  write_csv_safe(dplyr::bind_rows(continuous_cell_model_rows), file.path(out_dir, "cell_level_continuous_cell_ploidy_dose_lm_coefficients_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(continuous_cell_anova_rows), file.path(out_dir, "cell_level_continuous_cell_ploidy_dose_lm_anova_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(continuous_sample_model_rows), file.path(out_dir, "sample_level_continuous_cell_ploidy_dose_lm_coefficients_by_method.csv"))
  write_csv_safe(dplyr::bind_rows(continuous_sample_anova_rows), file.path(out_dir, "sample_level_continuous_cell_ploidy_dose_lm_anova_by_method.csv"))

  cell_plot <- ggplot(dat, aes(x = .data$Dose, y = .data$pseudotime, fill = .data$method)) +
    geom_violin(position = position_dodge(width = 0.82), scale = "width", trim = TRUE, alpha = 0.55, linewidth = 0.2) +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.16, outlier.size = 0.15, alpha = 0.82, linewidth = 0.2) +
    facet_wrap(~annotation_ploidy, nrow = 1) +
    scale_fill_manual(values = method_colors, drop = FALSE, name = "Method") +
    labs(title = "Tumor/ploidy_all pseudotime response to Dose within 2N and 4N", x = "Dose", y = "Pseudotime") +
    theme_extra(base_size = 10)
  save_pdf(cell_plot, file.path(out_dir, "tumor_ploidy_dose_response_three_methods_violin_box.pdf"), width = 11, height = 5.5)

  sample_plot <- ggplot(sample_summary, aes(x = .data$Dose, y = .data$sample_median_pseudotime, color = .data$method, group = .data$method)) +
    geom_point(position = position_dodge(width = 0.45), size = 2.1, alpha = 0.85) +
    geom_line(position = position_dodge(width = 0.45), alpha = 0.55) +
    facet_wrap(~annotation_ploidy, nrow = 1) +
    scale_color_manual(values = method_colors, drop = FALSE, name = "Method") +
    labs(title = "Tumor/ploidy_all sample-median pseudotime response to Dose within 2N and 4N", x = "Dose", y = "Sample median pseudotime") +
    theme_extra(base_size = 10)
  save_pdf(sample_plot, file.path(out_dir, "tumor_ploidy_dose_response_sample_median_three_methods.pdf"), width = 11, height = 5.5)

  if (nrow(continuous_dat) > 0) {
    continuous_scatter <- ggplot(continuous_dat, aes(x = .data$cell_ploidy, y = .data$pseudotime, color = .data$annotation_ploidy)) +
      geom_point(alpha = 0.18, size = 0.35) +
      geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 0.55) +
      facet_grid(method ~ Dose) +
      scale_color_manual(values = c("2N" = "#2C7FB8", "4N" = "#D95F0E"), drop = FALSE, name = "Ploidy") +
      labs(title = "Tumor/ploidy_all continuous cell ploidy vs pseudotime by method and Dose", x = "Continuous cell ploidy", y = "Pseudotime") +
      theme_extra(base_size = 9) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(continuous_scatter, file.path(out_dir, "tumor_continuous_cell_ploidy_pseudotime_by_method_dose_scatter.pdf"), width = 12, height = 8.5)
  }

  status <- data.frame(
    analysis_object = "Tumor/ploidy_all",
    n_cell_rows = nrow(dat),
    n_sample_rows = nrow(sample_summary),
    n_continuous_cell_rows = nrow(continuous_dat),
    n_continuous_sample_rows = nrow(continuous_sample_summary),
    n_cell_dose_tests = nrow(cell_dose_tests$global),
    n_sample_dose_tests = nrow(sample_dose_tests$global),
    n_continuous_dose_correlations = nrow(continuous_corr_by_dose),
    n_continuous_ploidy_dose_correlations = nrow(continuous_corr_by_dose_ploidy),
    stringsAsFactors = FALSE
  )
  write_csv_safe(status, file.path(out_dir, "tumor_ploidy_dose_response_status.csv"))
  status
}

extract_dose_numeric <- function(x) {
  out <- safe_num(gsub("[^0-9.]+", "", as.character(x)))
  out[!is.finite(out)] <- NA_real_
  out
}

dose_levels_numeric <- function(x) {
  ux <- unique(as_clean_chr(x))
  ux <- ux[!is.na(ux) & nzchar(ux)]
  if (length(ux) == 0) return(character(0))
  nums <- extract_dose_numeric(ux)
  ux[order(ifelse(is.finite(nums), nums, Inf), ux)]
}

read_cluster_annotation_summary <- function(results_root) {
  file <- file.path(results_root, "03b_cluster_annotation_and_GSEA", "00_summary", "cluster_annotation_summary.csv")
  if (!file.exists(file)) {
    warning("Missing cluster annotation summary: ", file, call. = FALSE)
    return(data.frame())
  }
  readr::read_csv(file, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE, progress = FALSE) %>%
    dplyr::mutate(cluster = as_clean_chr(.data$cluster))
}

read_stratified_hallmark_gsea <- function(results_root) {
  file <- file.path(results_root, "03b_cluster_annotation_and_GSEA", "00_summary", "stratified_Hallmark_GSEA_all.csv")
  if (!file.exists(file)) {
    warning("Missing stratified Hallmark GSEA summary: ", file, call. = FALSE)
    return(data.frame())
  }
  readr::read_csv(file, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE, progress = FALSE) %>%
    dplyr::mutate(
      cluster = as_clean_chr(.data$cluster),
      padj = safe_num(.data$padj),
      pval = safe_num(.data$pval),
      NES = safe_num(.data$NES),
      ES = safe_num(.data$ES),
      n_ident_1 = safe_num(.data$n_ident_1),
      n_ident_2 = safe_num(.data$n_ident_2),
      hallmark_label = as_clean_chr(.data$hallmark_label),
      gsea_enrichment = dplyr::case_when(
        is.finite(.data$NES) & .data$NES > 0 ~ "4N_enriched",
        is.finite(.data$NES) & .data$NES < 0 ~ "2N_enriched",
        TRUE ~ NA_character_
      )
    )
}

summarise_cluster_ploidy_shift <- function(dat) {
  summary_wide <- dat %>%
    dplyr::filter(
      .data$annotation_ploidy %in% c("2N", "4N"),
      is.finite(.data$pseudotime),
      !is.na(.data$clusters),
      nzchar(as.character(.data$clusters))
    ) %>%
    dplyr::group_by(.data$clusters, .data$annotation_ploidy) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = stats::median(.data$pseudotime, na.rm = TRUE),
      mean = mean(.data$pseudotime, na.rm = TRUE),
      q25 = as.numeric(stats::quantile(.data$pseudotime, 0.25, na.rm = TRUE, names = FALSE)),
      q75 = as.numeric(stats::quantile(.data$pseudotime, 0.75, na.rm = TRUE, names = FALSE)),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = annotation_ploidy,
      values_from = c(n, median, mean, q25, q75),
      names_glue = "{.value}_{annotation_ploidy}"
    )

  test_rows <- lapply(cluster_levels(dat$clusters), function(cl) {
    sub <- dat %>%
      dplyr::filter(.data$clusters == cl, .data$annotation_ploidy %in% c("2N", "4N"), is.finite(.data$pseudotime))
    test <- group_test(sub, "pseudotime", "annotation_ploidy", "scvelo_cluster_4N_vs_2N_pseudotime")
    dplyr::bind_cols(data.frame(clusters = cl, stringsAsFactors = FALSE), test)
  })
  tests <- dplyr::bind_rows(test_rows) %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, method = "BH"))

  summary_wide %>%
    dplyr::left_join(tests, by = "clusters") %>%
    dplyr::mutate(
      delta_median_4N_minus_2N = .data$median_4N - .data$median_2N,
      delta_mean_4N_minus_2N = .data$mean_4N - .data$mean_2N,
      time_shift_direction = dplyr::case_when(
        !is.finite(.data$delta_median_4N_minus_2N) ~ "insufficient",
        .data$delta_median_4N_minus_2N > 0 ~ "4N_delayed",
        .data$delta_median_4N_minus_2N < 0 ~ "4N_advanced",
        TRUE ~ "similar"
      )
    )
}

summarise_top_ploidy_gsea <- function(gsea_df, n_top = 5L) {
  if (nrow(gsea_df) == 0) return(data.frame())
  gsea_sig <- gsea_df %>%
    dplyr::filter(
      .data$scope == "within_cluster_Tumor_only_Ploidy_4N_vs_2N",
      !is.na(.data$cluster),
      nzchar(.data$cluster),
      is.finite(.data$NES),
      is.finite(.data$padj),
      .data$padj < 0.05
    )
  positive <- gsea_sig %>%
    dplyr::filter(.data$NES > 0) %>%
    dplyr::arrange(.data$cluster, .data$padj, dplyr::desc(.data$NES)) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gsea_enrichment = "4N_enriched")
  negative <- gsea_sig %>%
    dplyr::filter(.data$NES < 0) %>%
    dplyr::arrange(.data$cluster, .data$padj, .data$NES) %>%
    dplyr::group_by(.data$cluster) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gsea_enrichment = "2N_enriched")
  dplyr::bind_rows(positive, negative)
}

add_cluster_time_bins <- function(dat) {
  cluster_time <- dat %>%
    dplyr::filter(is.finite(.data$pseudotime), !is.na(.data$clusters), nzchar(as.character(.data$clusters))) %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_cluster_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      mean_cluster_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$median_cluster_pseudotime, .data$clusters) %>%
    dplyr::mutate(
      cluster_time_rank = dplyr::row_number(),
      n_clusters = dplyr::n(),
      rank_percentile = ifelse(.data$n_clusters > 1, (.data$cluster_time_rank - 1) / (.data$n_clusters - 1), NA_real_),
      time_bin = dplyr::case_when(
        .data$cluster_time_rank <= ceiling(.data$n_clusters / 3) ~ "early",
        .data$cluster_time_rank <= ceiling(2 * .data$n_clusters / 3) ~ "middle",
        TRUE ~ "late"
      )
    )
  list(
    data = dat %>% dplyr::left_join(cluster_time, by = "clusters"),
    cluster_time = cluster_time
  )
}

build_high_ploidy_long <- function(dat) {
  x <- dat$cell_ploidy[is.finite(dat$cell_ploidy)]
  if (length(x) == 0) return(list(data = data.frame(), thresholds = data.frame()))
  thresholds <- data.frame(
    high_ploidy_definition = c("top_quartile", "top_decile"),
    cell_ploidy_threshold = as.numeric(stats::quantile(x, c(0.75, 0.90), na.rm = TRUE, names = FALSE)),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(nrow(thresholds)), function(i) {
    dat %>%
      dplyr::mutate(
        high_ploidy_definition = thresholds$high_ploidy_definition[i],
        cell_ploidy_threshold = thresholds$cell_ploidy_threshold[i],
        high_cell_ploidy = .data$cell_ploidy >= thresholds$cell_ploidy_threshold[i]
      )
  })
  list(data = dplyr::bind_rows(rows), thresholds = thresholds)
}

tidy_model_coefficients <- function(model, model_label) {
  tidy_lm_coefficients(model, model_label)
}

fit_glm_safe <- function(formula, df, family) {
  tryCatch(stats::glm(formula, data = df, family = family), error = function(e) NULL)
}

run_scvelo_deep_dive_outputs <- function(metrics, output_base, results_root) {
  out_dir <- .ensure_dir(file.path(output_base, "08_scvelo_ploidy_dose_cluster_deep_dive"))
  dat <- metrics %>%
    dplyr::filter(
      .data$method == "scVelo",
      .data$analysis_object == "Tumor/ploidy_all",
      .data$analysis_ploidy_scope == "ploidy_all",
      is.finite(.data$pseudotime),
      !is.na(.data$clusters),
      nzchar(as.character(.data$clusters))
    ) %>%
    dplyr::mutate(
      clusters = as_clean_chr(.data$clusters),
      annotation_ploidy = factor(.data$annotation_ploidy, levels = c("2N", "4N")),
      Dose = factor(.data$Dose, levels = dose_levels_numeric(.data$Dose)),
      dose_numeric = extract_dose_numeric(.data$Dose)
    )
  if (nrow(dat) == 0) {
    warning("No scVelo Tumor/ploidy_all rows available for deep-dive analysis.", call. = FALSE)
    return(data.frame(status = "skipped_no_scvelo_tumor_ploidy_all", output_dir = out_dir, stringsAsFactors = FALSE))
  }

  cluster_bins <- add_cluster_time_bins(dat)
  dat <- cluster_bins$data
  cluster_time <- cluster_bins$cluster_time
  cluster_order <- cluster_time$clusters
  write_csv_safe(cluster_time, file.path(out_dir, "scvelo_cluster_median_pseudotime_time_bins.csv"))

  annotations <- read_cluster_annotation_summary(results_root)
  gsea_all <- read_stratified_hallmark_gsea(results_root)
  gsea_ploidy <- if (nrow(gsea_all) > 0) {
    gsea_all %>%
      dplyr::filter(.data$scope == "within_cluster_Tumor_only_Ploidy_4N_vs_2N", !is.na(.data$cluster), nzchar(.data$cluster))
  } else {
    data.frame()
  }
  top_gsea <- summarise_top_ploidy_gsea(gsea_all, n_top = 5L)
  write_csv_safe(gsea_ploidy, file.path(out_dir, "tumor_only_4N_vs_2N_hallmark_gsea_all_clusters.csv"))
  write_csv_safe(top_gsea, file.path(out_dir, "tumor_only_4N_vs_2N_top_hallmark_gsea_by_cluster.csv"))

  top_gsea_text <- if (nrow(top_gsea) > 0) {
    top_gsea %>%
      dplyr::group_by(.data$cluster, .data$gsea_enrichment) %>%
      dplyr::summarise(
        top_hallmark_pathways = paste(.data$hallmark_label, collapse = "; "),
        top_hallmark_NES = paste(round(.data$NES, 3), collapse = "; "),
        .groups = "drop"
      ) %>%
      tidyr::pivot_wider(
        names_from = gsea_enrichment,
        values_from = c(top_hallmark_pathways, top_hallmark_NES),
        names_glue = "{gsea_enrichment}_{.value}"
      )
  } else {
    data.frame()
  }

  shift_summary <- summarise_cluster_ploidy_shift(dat) %>%
    dplyr::left_join(cluster_time, by = "clusters") %>%
    dplyr::left_join(annotations, by = c("clusters" = "cluster")) %>%
    dplyr::left_join(top_gsea_text, by = c("clusters" = "cluster")) %>%
    dplyr::arrange(.data$median_cluster_pseudotime, .data$clusters)
  write_csv_safe(shift_summary, file.path(out_dir, "cluster_4N_vs_2N_scvelo_pseudotime_summary_with_gsea.csv"))

  sample_shift <- dat %>%
    dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N"), !is.na(.data$sample), nzchar(as.character(.data$sample))) %>%
    dplyr::group_by(.data$clusters, .data$sample, .data$annotation_ploidy) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      sample_median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = annotation_ploidy, values_from = c(n_cells, sample_median_pseudotime), names_glue = "{.value}_{annotation_ploidy}") %>%
    dplyr::mutate(sample_delta_4N_minus_2N = .data$sample_median_pseudotime_4N - .data$sample_median_pseudotime_2N)
  sample_shift_tests <- sample_shift %>%
    dplyr::filter(is.finite(.data$sample_delta_4N_minus_2N)) %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(
      n_samples = dplyr::n(),
      median_sample_delta_4N_minus_2N = stats::median(.data$sample_delta_4N_minus_2N, na.rm = TRUE),
      p_value = if (dplyr::n() >= 3) {
        tryCatch(stats::wilcox.test(.data$sample_delta_4N_minus_2N, mu = 0, exact = FALSE)$p.value, error = function(e) NA_real_)
      } else {
        NA_real_
      },
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = stats::p.adjust(.data$p_value, method = "BH"))
  write_csv_safe(sample_shift, file.path(out_dir, "cluster_4N_vs_2N_sample_median_pseudotime_delta.csv"))
  write_csv_safe(sample_shift_tests, file.path(out_dir, "cluster_4N_vs_2N_sample_median_delta_tests.csv"))

  gsea_by_shift <- if (nrow(gsea_ploidy) > 0) {
    gsea_ploidy %>%
      dplyr::filter(is.finite(.data$padj), .data$padj < 0.05, is.finite(.data$NES)) %>%
      dplyr::inner_join(
        shift_summary %>%
          dplyr::select(
            cluster = clusters,
            time_shift_direction,
            delta_median_4N_minus_2N,
            shift_p_adj = p_adj,
            annotation_primary,
            annotation_secondary,
            annotation_tertiary
          ),
        by = "cluster"
      ) %>%
      dplyr::filter(.data$time_shift_direction %in% c("4N_advanced", "4N_delayed")) %>%
      dplyr::group_by(.data$time_shift_direction, .data$gsea_enrichment, .data$hallmark_label) %>%
      dplyr::summarise(
        n_clusters = dplyr::n_distinct(.data$cluster),
        clusters = paste(sort_maybe_numeric(unique(.data$cluster)), collapse = ";"),
        min_padj = min(.data$padj, na.rm = TRUE),
        median_NES = stats::median(.data$NES, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$time_shift_direction, .data$gsea_enrichment, dplyr::desc(.data$n_clusters), .data$min_padj)
  } else {
    data.frame()
  }
  write_csv_safe(gsea_by_shift, file.path(out_dir, "pathway_enrichment_by_4N_time_shift_direction.csv"))

  gsea_by_significant_shift <- if (nrow(gsea_ploidy) > 0) {
    gsea_ploidy %>%
      dplyr::filter(is.finite(.data$padj), .data$padj < 0.05, is.finite(.data$NES)) %>%
      dplyr::inner_join(
        shift_summary %>%
          dplyr::filter(is.finite(.data$p_adj), .data$p_adj < 0.05) %>%
          dplyr::select(
            cluster = clusters,
            time_shift_direction,
            delta_median_4N_minus_2N,
            shift_p_adj = p_adj,
            annotation_primary,
            annotation_secondary,
            annotation_tertiary
          ),
        by = "cluster"
      ) %>%
      dplyr::filter(.data$time_shift_direction %in% c("4N_advanced", "4N_delayed")) %>%
      dplyr::group_by(.data$time_shift_direction, .data$gsea_enrichment, .data$hallmark_label) %>%
      dplyr::summarise(
        n_clusters = dplyr::n_distinct(.data$cluster),
        clusters = paste(sort_maybe_numeric(unique(.data$cluster)), collapse = ";"),
        min_padj = min(.data$padj, na.rm = TRUE),
        median_NES = stats::median(.data$NES, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(.data$time_shift_direction, .data$gsea_enrichment, dplyr::desc(.data$n_clusters), .data$min_padj)
  } else {
    data.frame()
  }
  write_csv_safe(gsea_by_significant_shift, file.path(out_dir, "pathway_enrichment_by_significant_4N_time_shift_direction.csv"))

  plot_shift <- shift_summary %>%
    dplyr::filter(is.finite(.data$delta_median_4N_minus_2N)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_order),
      significance_label = dplyr::case_when(
        is.finite(.data$p_adj) & .data$p_adj < 0.001 ~ "***",
        is.finite(.data$p_adj) & .data$p_adj < 0.01 ~ "**",
        is.finite(.data$p_adj) & .data$p_adj < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
  p_delta <- ggplot(plot_shift, aes(x = .data$clusters, y = .data$delta_median_4N_minus_2N, color = .data$time_shift_direction)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey45") +
    geom_segment(aes(xend = .data$clusters, y = 0, yend = .data$delta_median_4N_minus_2N), linewidth = 0.85) +
    geom_point(size = 2.5) +
    geom_text(aes(label = .data$significance_label), vjust = ifelse(plot_shift$delta_median_4N_minus_2N >= 0, -0.8, 1.4), size = 3.2, show.legend = FALSE) +
    scale_color_manual(values = c("4N_advanced" = "#2C7FB8", "4N_delayed" = "#D95F0E", "similar" = "grey50", "insufficient" = "grey70"), drop = FALSE) +
    labs(
      title = "scVelo Tumor/ploidy_all: 4N versus 2N median pseudotime shift by cluster",
      x = "Cluster ordered by scVelo median pseudotime",
      y = "Median pseudotime difference: 4N - 2N",
      color = "Shift"
    ) +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0))
  save_pdf(p_delta, file.path(out_dir, "cluster_4N_vs_2N_pseudotime_delta_lollipop.pdf"), width = 9, height = 5.5)

  plot_dat <- dat %>%
    dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")) %>%
    dplyr::mutate(clusters = factor(.data$clusters, levels = cluster_order))
  p_violin <- ggplot(plot_dat, aes(x = .data$clusters, y = .data$pseudotime, fill = .data$annotation_ploidy)) +
    geom_violin(position = position_dodge(width = 0.82), scale = "width", trim = TRUE, alpha = 0.55, linewidth = 0.2) +
    geom_boxplot(position = position_dodge(width = 0.82), width = 0.15, outlier.size = 0.12, alpha = 0.85, linewidth = 0.2) +
    scale_fill_manual(values = c("2N" = "#2C7FB8", "4N" = "#D95F0E"), drop = FALSE, name = "Ploidy") +
    labs(
      title = "scVelo Tumor/ploidy_all: 2N and 4N pseudotime distributions by cluster",
      x = "Cluster ordered by scVelo median pseudotime",
      y = "Pseudotime"
    ) +
    theme_extra(base_size = 10) +
    theme(axis.text.x = element_text(angle = 0))
  save_pdf(p_violin, file.path(out_dir, "cluster_4N_vs_2N_pseudotime_violin_box.pdf"), width = 11, height = 5.8)

  if (nrow(gsea_ploidy) > 0) {
    selected_pathways <- gsea_ploidy %>%
      dplyr::filter(is.finite(.data$padj), .data$padj < 0.05, is.finite(.data$NES)) %>%
      dplyr::group_by(.data$hallmark_label) %>%
      dplyr::summarise(n_clusters = dplyr::n_distinct(.data$cluster), min_padj = min(.data$padj, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(.data$n_clusters), .data$min_padj) %>%
      dplyr::slice_head(n = 22) %>%
      dplyr::pull(.data$hallmark_label)
    heat_df <- gsea_ploidy %>%
      dplyr::filter(.data$hallmark_label %in% selected_pathways, .data$cluster %in% cluster_order) %>%
      dplyr::mutate(
        cluster = factor(.data$cluster, levels = cluster_order),
        hallmark_label = factor(.data$hallmark_label, levels = rev(selected_pathways))
      )
    if (nrow(heat_df) > 0) {
      p_gsea <- ggplot(heat_df, aes(x = .data$cluster, y = .data$hallmark_label, fill = .data$NES)) +
        geom_tile(color = "white", linewidth = 0.25) +
        scale_fill_gradient2(low = "#2C7FB8", mid = "white", high = "#D95F0E", midpoint = 0, na.value = "grey90", name = "NES\n4N vs 2N") +
        labs(
          title = "Tumor-only 4N vs 2N Hallmark GSEA across scVelo-ordered clusters",
          x = "Cluster ordered by scVelo median pseudotime",
          y = "Hallmark pathway"
        ) +
        theme_extra(base_size = 9) +
        theme(axis.text.x = element_text(angle = 0), axis.text.y = element_text(size = 7.5))
      save_pdf(p_gsea, file.path(out_dir, "cluster_4N_vs_2N_top_hallmark_gsea_heatmap.pdf"), width = 10, height = 7.5)
    }
  }

  annotation_plot_df <- shift_summary %>%
    dplyr::filter(is.finite(.data$delta_median_4N_minus_2N)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_order),
      annotation_primary = ifelse(is.na(.data$annotation_primary) | !nzchar(.data$annotation_primary), "NA", .data$annotation_primary)
    )
  p_annotation <- ggplot(annotation_plot_df, aes(x = .data$rank_percentile, y = .data$delta_median_4N_minus_2N, color = .data$time_shift_direction)) +
    geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    geom_point(aes(size = -log10(.data$p_adj)), alpha = 0.85) +
    geom_text(aes(label = .data$clusters), vjust = -0.8, size = 3, show.legend = FALSE) +
    facet_wrap(~annotation_primary, scales = "free_x") +
    scale_color_manual(values = c("4N_advanced" = "#2C7FB8", "4N_delayed" = "#D95F0E", "similar" = "grey50", "insufficient" = "grey70"), drop = FALSE) +
    labs(
      title = "scVelo 4N time shift over cluster trajectory with 03b primary annotations",
      x = "Cluster pseudotime rank percentile",
      y = "Median pseudotime difference: 4N - 2N",
      color = "Shift",
      size = "-log10(adj. p)"
    ) +
    theme_extra(base_size = 9)
  save_pdf(p_annotation, file.path(out_dir, "cluster_4N_vs_2N_time_shift_with_03b_annotations.pdf"), width = 12, height = 7.2)

  continuous_dat <- dat %>%
    dplyr::filter(is.finite(.data$cell_ploidy), !is.na(.data$Dose), nzchar(as.character(.data$Dose))) %>%
    dplyr::mutate(
      Dose = factor(.data$Dose, levels = dose_levels_numeric(.data$Dose)),
      time_bin = factor(.data$time_bin, levels = c("early", "middle", "late")),
      clusters = factor(.data$clusters, levels = cluster_order)
    )
  hp <- build_high_ploidy_long(continuous_dat)
  high_df <- hp$data
  write_csv_safe(hp$thresholds, file.path(out_dir, "continuous_cell_ploidy_high_ploidy_thresholds.csv"))

  cluster_dose_summary <- continuous_dat %>%
    dplyr::group_by(.data$clusters, .data$cluster_time_rank, .data$time_bin, .data$Dose) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      q25_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      q75_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      annotation_4N_fraction = mean(.data$annotation_ploidy == "4N", na.rm = TRUE),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(cluster_dose_summary, file.path(out_dir, "continuous_ploidy_cluster_dose_summary.csv"))

  time_bin_dose_summary <- if (nrow(high_df) > 0) {
    high_df %>%
      dplyr::group_by(.data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$time_bin, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        annotation_4N_fraction = mean(.data$annotation_ploidy == "4N", na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(time_bin_dose_summary, file.path(out_dir, "continuous_ploidy_time_bin_summary_by_dose.csv"))

  late_summary <- if (nrow(time_bin_dose_summary) > 0) {
    time_bin_dose_summary %>%
      dplyr::filter(.data$time_bin == "late") %>%
      dplyr::arrange(.data$high_ploidy_definition, .data$Dose)
  } else {
    data.frame()
  }
  write_csv_safe(late_summary, file.path(out_dir, "continuous_ploidy_late_cluster_summary_by_dose.csv"))

  high_distribution <- if (nrow(high_df) > 0) {
    high_df %>%
      dplyr::group_by(.data$high_ploidy_definition, .data$Dose, .data$high_cell_ploidy, .data$time_bin) %>%
      dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
      dplyr::group_by(.data$high_ploidy_definition, .data$Dose, .data$high_cell_ploidy) %>%
      dplyr::mutate(fraction_in_time_bin = .data$n_cells / sum(.data$n_cells)) %>%
      dplyr::ungroup()
  } else {
    data.frame()
  }
  write_csv_safe(high_distribution, file.path(out_dir, "continuous_ploidy_time_bin_composition_by_dose_and_high_ploidy.csv"))

  cluster_high_fraction <- if (nrow(high_df) > 0) {
    high_df %>%
      dplyr::group_by(.data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$clusters, .data$cluster_time_rank, .data$time_bin, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(cluster_high_fraction, file.path(out_dir, "continuous_ploidy_high_fraction_by_cluster_dose.csv"))

  model_rows <- list()
  if (nrow(high_df) > 0) {
    for (def in unique(high_df$high_ploidy_definition)) {
      sub <- high_df %>%
        dplyr::filter(.data$high_ploidy_definition == def, is.finite(.data$dose_numeric), !is.na(.data$time_bin)) %>%
        dplyr::mutate(
          high_cell_ploidy_int = as.integer(.data$high_cell_ploidy),
          late_cluster_int = as.integer(.data$time_bin == "late")
        )
      if (nrow(sub) > 0 && dplyr::n_distinct(sub$high_cell_ploidy_int) > 1 && dplyr::n_distinct(sub$dose_numeric) > 1) {
        high_model <- fit_glm_safe(high_cell_ploidy_int ~ dose_numeric * time_bin, sub, family = stats::binomial())
        model_rows[[paste0(def, "_high")]] <- tidy_model_coefficients(high_model, paste0("high_cell_ploidy_", def, "_dose_time_bin_logistic"))
      }
      if (nrow(sub) > 0 && dplyr::n_distinct(sub$late_cluster_int) > 1 && dplyr::n_distinct(sub$dose_numeric) > 1 && dplyr::n_distinct(sub$cell_ploidy) > 1) {
        late_model <- fit_glm_safe(late_cluster_int ~ cell_ploidy * dose_numeric, sub %>% dplyr::distinct(.data$cell, .keep_all = TRUE), family = stats::binomial())
        model_rows[[paste0(def, "_late")]] <- tidy_model_coefficients(late_model, paste0("late_cluster_", def, "_cell_ploidy_dose_logistic"))
      }
    }
  }
  model_coefficients <- dplyr::bind_rows(model_rows)
  write_csv_safe(model_coefficients, file.path(out_dir, "continuous_ploidy_late_cluster_dose_logistic_coefficients.csv"))

  sample_late_summary <- continuous_dat %>%
    dplyr::filter(!is.na(.data$sample), nzchar(as.character(.data$sample))) %>%
    dplyr::group_by(.data$sample, .data$Dose) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      late_cluster_fraction = mean(.data$time_bin == "late", na.rm = TRUE),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      annotation_4N_fraction = mean(.data$annotation_ploidy == "4N", na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(sample_late_summary, file.path(out_dir, "continuous_ploidy_sample_late_cluster_summary_by_dose.csv"))

  if (nrow(continuous_dat) > 0) {
    p_ploidy_cluster_dose <- ggplot(continuous_dat, aes(x = .data$clusters, y = .data$cell_ploidy, fill = .data$Dose)) +
      geom_boxplot(position = position_dodge(width = 0.78), width = 0.62, outlier.size = 0.1, alpha = 0.82, linewidth = 0.2) +
      labs(
        title = "scVelo Tumor/ploidy_all: continuous cell ploidy by cluster and Dose",
        x = "Cluster ordered by scVelo median pseudotime",
        y = "Continuous cell ploidy",
        fill = "Dose"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_ploidy_cluster_dose, file.path(out_dir, "continuous_cell_ploidy_by_cluster_and_dose_boxplot.pdf"), width = 11, height = 5.8)
  }
  if (nrow(time_bin_dose_summary) > 0) {
    p_late_fraction <- ggplot(time_bin_dose_summary, aes(x = .data$Dose, y = .data$high_cell_ploidy_fraction, color = .data$time_bin, group = .data$time_bin)) +
      geom_point(size = 2.2) +
      geom_line(linewidth = 0.75) +
      facet_wrap(~high_ploidy_definition, nrow = 1) +
      scale_color_manual(values = c("early" = "#1B9E77", "middle" = "#7570B3", "late" = "#D95F02"), drop = FALSE) +
      labs(
        title = "High continuous-ploidy fraction across scVelo early/middle/late clusters by Dose",
        x = "Dose",
        y = "High cell-ploidy fraction",
        color = "Cluster time bin"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_late_fraction, file.path(out_dir, "continuous_high_ploidy_fraction_by_time_bin_and_dose.pdf"), width = 9.5, height = 5.2)

    p_mean_timebin <- ggplot(time_bin_dose_summary, aes(x = .data$Dose, y = .data$mean_cell_ploidy, color = .data$time_bin, group = .data$time_bin)) +
      geom_point(size = 2.2) +
      geom_line(linewidth = 0.75) +
      scale_color_manual(values = c("early" = "#1B9E77", "middle" = "#7570B3", "late" = "#D95F02"), drop = FALSE) +
      labs(
        title = "Mean continuous cell ploidy across scVelo early/middle/late clusters by Dose",
        x = "Dose",
        y = "Mean continuous cell ploidy",
        color = "Cluster time bin"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_mean_timebin, file.path(out_dir, "continuous_mean_cell_ploidy_by_time_bin_and_dose.pdf"), width = 8.5, height = 5.2)
  }
  if (nrow(high_distribution) > 0) {
    p_composition <- high_distribution %>%
      dplyr::filter(.data$high_cell_ploidy) %>%
      ggplot(aes(x = .data$Dose, y = .data$fraction_in_time_bin, fill = .data$time_bin)) +
      geom_col(width = 0.75, color = "white", linewidth = 0.2) +
      facet_wrap(~high_ploidy_definition, nrow = 1) +
      scale_fill_manual(values = c("early" = "#1B9E77", "middle" = "#7570B3", "late" = "#D95F02"), drop = FALSE) +
      labs(
        title = "Time-bin composition of high continuous-ploidy cells by Dose",
        x = "Dose",
        y = "Fraction of high-ploidy cells",
        fill = "Cluster time bin"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_composition, file.path(out_dir, "continuous_high_ploidy_time_bin_composition_by_dose.pdf"), width = 9.5, height = 5.2)
  }
  if (nrow(cluster_high_fraction) > 0) {
    p_heat <- cluster_high_fraction %>%
      dplyr::mutate(clusters = factor(.data$clusters, levels = rev(cluster_order))) %>%
      ggplot(aes(x = .data$Dose, y = .data$clusters, fill = .data$high_cell_ploidy_fraction)) +
      geom_tile(color = "white", linewidth = 0.25) +
      facet_wrap(~high_ploidy_definition, nrow = 1) +
      scale_fill_gradient(low = "white", high = "#B03A2E", na.value = "grey90", name = "High\nfraction") +
      labs(
        title = "High continuous-ploidy fraction by Dose and scVelo-ordered cluster",
        x = "Dose",
        y = "Cluster"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_heat, file.path(out_dir, "continuous_high_ploidy_fraction_cluster_dose_heatmap.pdf"), width = 9.5, height = 6.2)
  }

  annotation_continuous_dat <- continuous_dat %>%
    dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")) %>%
    dplyr::mutate(annotation_ploidy = factor(.data$annotation_ploidy, levels = c("2N", "4N")))
  annotation_high_global <- if (nrow(high_df) > 0) {
    high_df %>%
      dplyr::filter(.data$annotation_ploidy %in% c("2N", "4N")) %>%
      dplyr::mutate(annotation_ploidy = factor(.data$annotation_ploidy, levels = c("2N", "4N")))
  } else {
    data.frame()
  }

  annotation_dose_summary <- if (nrow(annotation_high_global) > 0) {
    annotation_high_global %>%
      dplyr::group_by(.data$annotation_ploidy, .data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
        q25_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
        q75_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(annotation_dose_summary, file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_summary_global_thresholds.csv"))

  annotation_time_bin_dose_summary <- if (nrow(annotation_high_global) > 0) {
    annotation_high_global %>%
      dplyr::group_by(.data$annotation_ploidy, .data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$time_bin, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(annotation_time_bin_dose_summary, file.path(out_dir, "annotation_ploidy_continuous_ploidy_time_bin_dose_summary_global_thresholds.csv"))
  write_csv_safe(
    annotation_time_bin_dose_summary %>% dplyr::filter(.data$time_bin == "late"),
    file.path(out_dir, "annotation_ploidy_continuous_ploidy_late_cluster_summary_by_dose_global_thresholds.csv")
  )

  annotation_cluster_dose_summary <- annotation_continuous_dat %>%
    dplyr::group_by(.data$annotation_ploidy, .data$clusters, .data$cluster_time_rank, .data$time_bin, .data$Dose) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      q25_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.25, na.rm = TRUE, names = FALSE)),
      q75_cell_ploidy = as.numeric(stats::quantile(.data$cell_ploidy, 0.75, na.rm = TRUE, names = FALSE)),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    )
  write_csv_safe(annotation_cluster_dose_summary, file.path(out_dir, "annotation_ploidy_continuous_ploidy_cluster_dose_summary.csv"))

  annotation_dose_tests <- dose_effect_by_group_tests(
    annotation_continuous_dat,
    "cell_ploidy",
    c("annotation_ploidy"),
    "dose_effect_on_continuous_cell_ploidy_within_annotation_ploidy"
  )
  write_csv_safe(annotation_dose_tests$global, file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_global_tests.csv"))
  write_csv_safe(annotation_dose_tests$pairwise, file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_pairwise_tests.csv"))

  annotation_model_rows <- list()
  annotation_model_anova_rows <- list()
  annotation_cor_rows <- list()
  for (ploidy_label in c("2N", "4N")) {
    sub <- annotation_continuous_dat %>%
      dplyr::filter(.data$annotation_ploidy == .env$ploidy_label, is.finite(.data$dose_numeric), is.finite(.data$cell_ploidy))
    model <- if (nrow(sub) > 0 && dplyr::n_distinct(sub$Dose) > 1) {
      fit_lm_safe(cell_ploidy ~ Dose, sub)
    } else {
      NULL
    }
    annotation_model_rows[[ploidy_label]] <- tidy_lm_coefficients(model, paste0("continuous_cell_ploidy_by_Dose_within_", ploidy_label))
    annotation_model_anova_rows[[ploidy_label]] <- tidy_lm_anova(model, paste0("continuous_cell_ploidy_by_Dose_within_", ploidy_label))
    sp <- cor_test_one(sub$dose_numeric, sub$cell_ploidy, "spearman")
    pe <- cor_test_one(sub$dose_numeric, sub$cell_ploidy, "pearson")
    annotation_cor_rows[[ploidy_label]] <- data.frame(
      annotation_ploidy = ploidy_label,
      n_cells = nrow(sub),
      spearman_rho = sp[["estimate"]],
      spearman_p = sp[["p_value"]],
      pearson_r = pe[["estimate"]],
      pearson_p = pe[["p_value"]],
      stringsAsFactors = FALSE
    )
  }
  annotation_cor <- dplyr::bind_rows(annotation_cor_rows)
  annotation_cor$spearman_p_adj <- stats::p.adjust(annotation_cor$spearman_p, method = "BH")
  annotation_cor$pearson_p_adj <- stats::p.adjust(annotation_cor$pearson_p, method = "BH")
  write_csv_safe(dplyr::bind_rows(annotation_model_rows), file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_lm_coefficients.csv"))
  write_csv_safe(dplyr::bind_rows(annotation_model_anova_rows), file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_lm_anova.csv"))
  write_csv_safe(annotation_cor, file.path(out_dir, "annotation_ploidy_continuous_ploidy_dose_numeric_correlations.csv"))

  annotation_threshold_rows <- lapply(c("2N", "4N"), function(ploidy_label) {
    sub <- annotation_continuous_dat %>% dplyr::filter(.data$annotation_ploidy == .env$ploidy_label)
    x <- sub$cell_ploidy[is.finite(sub$cell_ploidy)]
    if (length(x) == 0) return(data.frame())
    data.frame(
      annotation_ploidy = ploidy_label,
      high_ploidy_definition = c("within_top_quartile", "within_top_decile"),
      cell_ploidy_threshold = as.numeric(stats::quantile(x, c(0.75, 0.90), na.rm = TRUE, names = FALSE)),
      stringsAsFactors = FALSE
    )
  })
  annotation_thresholds <- dplyr::bind_rows(annotation_threshold_rows)
  write_csv_safe(annotation_thresholds, file.path(out_dir, "annotation_ploidy_specific_high_ploidy_thresholds.csv"))

  annotation_high_within <- if (nrow(annotation_thresholds) > 0) {
    threshold_splits <- split(annotation_thresholds, seq_len(nrow(annotation_thresholds)))
    dplyr::bind_rows(lapply(threshold_splits, function(thr) {
      annotation_continuous_dat %>%
        dplyr::filter(.data$annotation_ploidy == thr$annotation_ploidy[1]) %>%
        dplyr::mutate(
          high_ploidy_definition = thr$high_ploidy_definition[1],
          cell_ploidy_threshold = thr$cell_ploidy_threshold[1],
          high_cell_ploidy = .data$cell_ploidy >= thr$cell_ploidy_threshold[1]
        )
    }))
  } else {
    data.frame()
  }
  annotation_high_within_time_bin_summary <- if (nrow(annotation_high_within) > 0) {
    annotation_high_within %>%
      dplyr::group_by(.data$annotation_ploidy, .data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$time_bin, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(annotation_high_within_time_bin_summary, file.path(out_dir, "annotation_ploidy_specific_high_ploidy_time_bin_dose_summary.csv"))
  write_csv_safe(
    annotation_high_within_time_bin_summary %>% dplyr::filter(.data$time_bin == "late"),
    file.path(out_dir, "annotation_ploidy_specific_high_ploidy_late_cluster_summary_by_dose.csv")
  )

  annotation_high_within_cluster_summary <- if (nrow(annotation_high_within) > 0) {
    annotation_high_within %>%
      dplyr::group_by(.data$annotation_ploidy, .data$high_ploidy_definition, .data$cell_ploidy_threshold, .data$clusters, .data$cluster_time_rank, .data$time_bin, .data$Dose) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        high_cell_ploidy_n = sum(.data$high_cell_ploidy, na.rm = TRUE),
        high_cell_ploidy_fraction = mean(.data$high_cell_ploidy, na.rm = TRUE),
        mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    data.frame()
  }
  write_csv_safe(annotation_high_within_cluster_summary, file.path(out_dir, "annotation_ploidy_specific_high_ploidy_cluster_dose_summary.csv"))

  if (nrow(annotation_continuous_dat) > 0) {
    p_annotation_dose <- ggplot(annotation_continuous_dat, aes(x = .data$Dose, y = .data$cell_ploidy, fill = .data$Dose)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.5, linewidth = 0.2) +
      geom_boxplot(width = 0.18, outlier.size = 0.1, alpha = 0.82, linewidth = 0.2) +
      facet_wrap(~annotation_ploidy, nrow = 1) +
      labs(
        title = "scVelo Tumor/ploidy_all: continuous cell ploidy response to Dose within 2N and 4N",
        x = "Dose",
        y = "Continuous cell ploidy",
        fill = "Dose"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_annotation_dose, file.path(out_dir, "annotation_ploidy_continuous_cell_ploidy_by_dose_violin_box.pdf"), width = 9.5, height = 5.2)
  }
  if (nrow(annotation_time_bin_dose_summary) > 0) {
    p_annotation_timebin <- annotation_time_bin_dose_summary %>%
      dplyr::filter(.data$high_ploidy_definition == "top_quartile") %>%
      ggplot(aes(x = .data$Dose, y = .data$mean_cell_ploidy, color = .data$time_bin, group = .data$time_bin)) +
      geom_point(size = 2.1) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~annotation_ploidy, nrow = 1) +
      scale_color_manual(values = c("early" = "#1B9E77", "middle" = "#7570B3", "late" = "#D95F02"), drop = FALSE) +
      labs(
        title = "Mean continuous cell ploidy by Dose within annotation ploidy and scVelo time bin",
        x = "Dose",
        y = "Mean continuous cell ploidy",
        color = "Cluster time bin"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_annotation_timebin, file.path(out_dir, "annotation_ploidy_mean_cell_ploidy_by_time_bin_and_dose.pdf"), width = 10.5, height = 5.2)
  }
  if (nrow(annotation_high_within_time_bin_summary) > 0) {
    p_annotation_high <- annotation_high_within_time_bin_summary %>%
      dplyr::filter(.data$high_ploidy_definition == "within_top_quartile") %>%
      ggplot(aes(x = .data$Dose, y = .data$high_cell_ploidy_fraction, color = .data$time_bin, group = .data$time_bin)) +
      geom_point(size = 2.1) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~annotation_ploidy, nrow = 1) +
      scale_color_manual(values = c("early" = "#1B9E77", "middle" = "#7570B3", "late" = "#D95F02"), drop = FALSE) +
      labs(
        title = "Within-annotation high continuous-ploidy fraction by Dose and scVelo time bin",
        x = "Dose",
        y = "Within-annotation top-quartile fraction",
        color = "Cluster time bin"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_annotation_high, file.path(out_dir, "annotation_ploidy_specific_high_ploidy_fraction_by_time_bin_and_dose.pdf"), width = 10.5, height = 5.2)
  }
  if (nrow(annotation_cluster_dose_summary) > 0) {
    p_annotation_cluster_mean <- annotation_cluster_dose_summary %>%
      dplyr::mutate(clusters = factor(.data$clusters, levels = rev(cluster_order))) %>%
      ggplot(aes(x = .data$Dose, y = .data$clusters, fill = .data$mean_cell_ploidy)) +
      geom_tile(color = "white", linewidth = 0.25) +
      facet_wrap(~annotation_ploidy, nrow = 1) +
      scale_fill_gradientn(
        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
        na.value = "grey90",
        name = "Mean\nploidy"
      ) +
      labs(
        title = "Mean continuous cell ploidy by Dose and scVelo-ordered cluster within 2N and 4N",
        x = "Dose",
        y = "Cluster"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_annotation_cluster_mean, file.path(out_dir, "annotation_ploidy_mean_cell_ploidy_cluster_dose_heatmap.pdf"), width = 10, height = 6.2)
  }
  if (nrow(annotation_high_within_cluster_summary) > 0) {
    p_annotation_cluster_high <- annotation_high_within_cluster_summary %>%
      dplyr::filter(.data$high_ploidy_definition == "within_top_quartile") %>%
      dplyr::mutate(clusters = factor(.data$clusters, levels = rev(cluster_order))) %>%
      ggplot(aes(x = .data$Dose, y = .data$clusters, fill = .data$high_cell_ploidy_fraction)) +
      geom_tile(color = "white", linewidth = 0.25) +
      facet_wrap(~annotation_ploidy, nrow = 1) +
      scale_fill_gradient(low = "white", high = "#B03A2E", na.value = "grey90", name = "Within top\nquartile") +
      labs(
        title = "Within-annotation top-quartile continuous-ploidy fraction by Dose and cluster",
        x = "Dose",
        y = "Cluster"
      ) +
      theme_extra(base_size = 10) +
      theme(axis.text.x = element_text(angle = 0))
    save_pdf(p_annotation_cluster_high, file.path(out_dir, "annotation_ploidy_specific_high_ploidy_fraction_cluster_dose_heatmap.pdf"), width = 10, height = 6.2)
  }

  status <- data.frame(
    status = "ok",
    output_dir = out_dir,
    n_scvelo_tumor_cells = nrow(dat),
    n_clusters = dplyr::n_distinct(dat$clusters),
    n_shift_summary_rows = nrow(shift_summary),
    n_top_gsea_rows = nrow(top_gsea),
    n_gsea_by_shift_rows = nrow(gsea_by_shift),
    n_gsea_by_significant_shift_rows = nrow(gsea_by_significant_shift),
    n_continuous_cells = nrow(continuous_dat),
    n_late_summary_rows = nrow(late_summary),
    n_annotation_continuous_cells = nrow(annotation_continuous_dat),
    n_annotation_dose_summary_rows = nrow(annotation_dose_summary),
    n_annotation_specific_threshold_rows = nrow(annotation_thresholds),
    stringsAsFactors = FALSE
  )
  write_csv_safe(status, file.path(out_dir, "scvelo_deep_dive_status.csv"))
  status
}

run_ploidy_pseudotime_outputs <- function(df, method_dir, method_label, object_label, object_safe) {
  out_dir <- .ensure_dir(file.path(method_dir, object_safe, "01_ploidy_pseudotime"))
  dat <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), .data$annotation_ploidy %in% c("2N", "4N"))
  summary_df <- summary_by_group(dat, "annotation_ploidy", "pseudotime")
  tests <- group_test(dat, "pseudotime", "annotation_ploidy", "pseudotime_by_annotation_ploidy")
  write_csv_safe(summary_df, file.path(out_dir, "pseudotime_by_annotation_ploidy_summary.csv"))
  write_csv_safe(tests, file.path(out_dir, "pseudotime_by_annotation_ploidy_test.csv"))

  hist_df <- make_frequency_bins(dat, "pseudotime", "annotation_ploidy", binwidth = diff(range(dat$pseudotime, na.rm = TRUE)) / 80)
  write_csv_safe(hist_df, file.path(out_dir, "ploidy_pseudotime_frequency_histogram_data.csv"))
  p <- plot_frequency_histogram(
    hist_df,
    "annotation_ploidy",
    paste(method_label, object_label, "2N vs 4N pseudotime frequency"),
    "Pseudotime"
  )
  if (!is.null(p)) {
    save_pdf(p, file.path(out_dir, "ploidy_pseudotime_frequency_histogram.pdf"), width = 8.5, height = 6)
  }

  sample_status <- save_sample_frequency_pages(
    dat,
    "annotation_ploidy",
    file.path(out_dir, "sample_pages", "all_samples_ploidy_pseudotime_frequency_histograms.pdf"),
    paste(method_label, object_label, "sample-level 2N vs 4N pseudotime frequency"),
    "Pseudotime"
  )
  write_csv_safe(sample_status, file.path(out_dir, "sample_pages", "sample_histogram_status.csv"))
  list(summary = summary_df, tests = tests)
}

run_dose_pseudotime_outputs <- function(df, method_dir, method_label, object_label, object_safe) {
  out_dir <- .ensure_dir(file.path(method_dir, object_safe, "02_dose_pseudotime"))
  dat <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), !is.na(.data$Dose), nzchar(as.character(.data$Dose)))
  summary_df <- summary_by_group(dat, "Dose", "pseudotime")
  tests <- group_test(dat, "pseudotime", "Dose", "pseudotime_by_dose")
  pairwise <- pairwise_group_tests(dat, "pseudotime", "Dose", "pseudotime_by_dose_pairwise")
  write_csv_safe(summary_df, file.path(out_dir, "pseudotime_by_dose_summary.csv"))
  write_csv_safe(tests, file.path(out_dir, "pseudotime_by_dose_test.csv"))
  write_csv_safe(pairwise, file.path(out_dir, "pseudotime_by_dose_pairwise_tests.csv"))

  hist_df <- make_frequency_bins(dat, "pseudotime", "Dose", binwidth = diff(range(dat$pseudotime, na.rm = TRUE)) / 80)
  write_csv_safe(hist_df, file.path(out_dir, "dose_pseudotime_frequency_histogram_data.csv"))
  p <- plot_frequency_histogram(
    hist_df,
    "Dose",
    paste(method_label, object_label, "Dose pseudotime frequency"),
    "Pseudotime"
  )
  if (!is.null(p)) {
    save_pdf(p, file.path(out_dir, "dose_pseudotime_frequency_histogram.pdf"), width = 8.5, height = 6)
  }

  sample_status <- save_sample_frequency_pages(
    dat,
    "Dose",
    file.path(out_dir, "sample_pages", "all_samples_dose_pseudotime_frequency_histograms.pdf"),
    paste(method_label, object_label, "sample-level Dose pseudotime frequency"),
    "Pseudotime"
  )
  write_csv_safe(sample_status, file.path(out_dir, "sample_pages", "sample_histogram_status.csv"))
  list(summary = summary_df, tests = tests, pairwise = pairwise)
}

run_cluster_correlation_outputs <- function(df, method_dir, method_label, object_label, object_safe) {
  out_dir <- .ensure_dir(file.path(method_dir, object_safe, "03_cluster_cell_ploidy_pseudotime_correlation"))
  corr_df <- cluster_ploidy_correlations(df)
  write_csv_safe(corr_df, file.path(out_dir, "cluster_cell_ploidy_pseudotime_correlations.csv"))
  p <- plot_cluster_ploidy_correlation(
    df,
    paste(method_label, object_label, "cell ploidy and pseudotime by cluster")
  )
  if (!is.null(p)) {
    save_pdf(p, file.path(out_dir, "cluster_cell_ploidy_pseudotime_correlation.pdf"), width = 10, height = 8)
  }
  corr_df
}

run_tumor_ploidy_dose_outputs <- function(df, method_dir, method_label, object_safe) {
  out_dir <- .ensure_dir(file.path(method_dir, object_safe, "04_tumor_ploidy_dose_pseudotime"))
  dat <- df %>%
    dplyr::filter(is.finite(.data$pseudotime), .data$annotation_ploidy %in% c("2N", "4N"), !is.na(.data$Dose), nzchar(as.character(.data$Dose)))

  annotation_tests <- dplyr::bind_rows(
    group_test(dat, "pseudotime", "annotation_ploidy", "tumor_pseudotime_by_annotation_ploidy"),
    group_test(dat, "pseudotime", "Dose", "tumor_pseudotime_by_dose")
  )
  annotation_pairwise <- pairwise_group_tests(dat, "pseudotime", "Dose", "tumor_pseudotime_by_dose_pairwise")
  annotation_summary <- summary_by_group(dat, c("annotation_ploidy", "Dose"), "pseudotime")
  annotation_model_df <- dat %>%
    dplyr::filter(!is.na(.data$annotation_ploidy), !is.na(.data$Dose)) %>%
    dplyr::mutate(annotation_ploidy = factor(.data$annotation_ploidy), Dose = factor(.data$Dose, levels = sort_maybe_numeric(unique(.data$Dose))))
  annotation_model <- if (nrow(annotation_model_df) > 0 && length(unique(annotation_model_df$annotation_ploidy)) > 1 && length(unique(annotation_model_df$Dose)) > 1) {
    fit_lm_safe(pseudotime ~ annotation_ploidy * Dose, annotation_model_df)
  } else {
    NULL
  }

  write_csv_safe(annotation_summary, file.path(out_dir, "annotation_ploidy_dose_pseudotime_summary.csv"))
  write_csv_safe(annotation_tests, file.path(out_dir, "annotation_ploidy_dose_pseudotime_tests.csv"))
  write_csv_safe(annotation_pairwise, file.path(out_dir, "annotation_ploidy_dose_pseudotime_pairwise_tests.csv"))
  write_csv_safe(tidy_lm_coefficients(annotation_model, "pseudotime_annotation_ploidy_dose"), file.path(out_dir, "annotation_ploidy_dose_lm_coefficients.csv"))
  write_csv_safe(tidy_lm_anova(annotation_model, "pseudotime_annotation_ploidy_dose"), file.path(out_dir, "annotation_ploidy_dose_lm_anova.csv"))

  continuous_dat <- dat %>% dplyr::filter(is.finite(.data$cell_ploidy))
  continuous_tests <- dplyr::bind_rows(
    group_test(continuous_dat, "cell_ploidy", "annotation_ploidy", "cell_ploidy_by_annotation_ploidy"),
    group_test(continuous_dat, "cell_ploidy", "Dose", "cell_ploidy_by_dose"),
    group_test(continuous_dat, "pseudotime", "Dose", "continuous_subset_pseudotime_by_dose")
  )
  continuous_model_df <- continuous_dat %>%
    dplyr::filter(is.finite(.data$cell_ploidy), !is.na(.data$Dose)) %>%
    dplyr::mutate(Dose = factor(.data$Dose, levels = sort_maybe_numeric(unique(.data$Dose))))
  continuous_model <- if (nrow(continuous_model_df) > 0 && length(unique(continuous_model_df$Dose)) > 1 && length(unique(continuous_model_df$cell_ploidy)) > 1) {
    fit_lm_safe(pseudotime ~ cell_ploidy * Dose, continuous_model_df)
  } else {
    NULL
  }
  cluster_continuous_summary <- continuous_dat %>%
    dplyr::filter(!is.na(.data$clusters), nzchar(.data$clusters)) %>%
    dplyr::group_by(.data$clusters, .data$annotation_ploidy, .data$Dose) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE),
      median_cell_ploidy = stats::median(.data$cell_ploidy, na.rm = TRUE),
      mean_pseudotime = mean(.data$pseudotime, na.rm = TRUE),
      median_pseudotime = stats::median(.data$pseudotime, na.rm = TRUE),
      .groups = "drop"
    )
  cluster_model_df <- cluster_continuous_summary %>%
    dplyr::filter(is.finite(.data$mean_cell_ploidy), is.finite(.data$median_pseudotime), !is.na(.data$Dose)) %>%
    dplyr::mutate(Dose = factor(.data$Dose, levels = sort_maybe_numeric(unique(.data$Dose))))
  cluster_model <- if (nrow(cluster_model_df) > 0 && length(unique(cluster_model_df$Dose)) > 1 && length(unique(cluster_model_df$mean_cell_ploidy)) > 1) {
    fit_lm_safe(median_pseudotime ~ mean_cell_ploidy * Dose, cluster_model_df)
  } else {
    NULL
  }

  write_csv_safe(continuous_tests, file.path(out_dir, "continuous_cell_ploidy_tests.csv"))
  write_csv_safe(tidy_lm_coefficients(continuous_model, "pseudotime_cell_ploidy_dose"), file.path(out_dir, "continuous_cell_ploidy_dose_lm_coefficients.csv"))
  write_csv_safe(tidy_lm_anova(continuous_model, "pseudotime_cell_ploidy_dose"), file.path(out_dir, "continuous_cell_ploidy_dose_lm_anova.csv"))
  write_csv_safe(cluster_continuous_summary, file.path(out_dir, "cluster_mean_cell_ploidy_pseudotime_summary.csv"))
  write_csv_safe(tidy_lm_coefficients(cluster_model, "cluster_median_pseudotime_mean_cell_ploidy_dose"), file.path(out_dir, "cluster_mean_cell_ploidy_dose_lm_coefficients.csv"))
  write_csv_safe(tidy_lm_anova(cluster_model, "cluster_median_pseudotime_mean_cell_ploidy_dose"), file.path(out_dir, "cluster_mean_cell_ploidy_dose_lm_anova.csv"))
  added_outputs <- run_tumor_cell_ploidy_time_added_outputs(dat, out_dir, method_label)

  data.frame(
    annotation_rows = nrow(dat),
    continuous_rows = nrow(continuous_dat),
    cluster_summary_rows = nrow(cluster_continuous_summary),
    added_output_rows = nrow(added_outputs),
    stringsAsFactors = FALSE
  )
}

run_cluster_ploidy_distribution_outputs <- function(df, method_dir, method_label, object_label, object_safe) {
  out_dir <- .ensure_dir(file.path(method_dir, object_safe, "05_cluster_cell_ploidy_distribution"))
  status_rows <- list()
  for (subset_label in c("all", "2N", "4N")) {
    res <- plot_cluster_cell_ploidy_distribution(
      df,
      subset_label,
      paste(method_label, object_label, subset_label, "cell ploidy distribution by cluster")
    )
    hist_file <- file.path(out_dir, paste0("cluster_cell_ploidy_distribution_", subset_label, "_histogram_data.csv"))
    write_csv_safe(res$hist, hist_file)
    pdf_file <- file.path(out_dir, paste0("cluster_cell_ploidy_distribution_", subset_label, ".pdf"))
    if (!is.null(res$plot)) {
      save_pdf(res$plot, pdf_file, width = 10, height = 8)
    }
    status_rows[[subset_label]] <- data.frame(
      subset = subset_label,
      n_hist_rows = nrow(res$hist),
      binwidth = res$binwidth,
      output_pdf = if (!is.null(res$plot)) pdf_file else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  status <- dplyr::bind_rows(status_rows)
  write_csv_safe(status, file.path(out_dir, "cluster_cell_ploidy_distribution_status.csv"))
  status
}

run_method_object_analysis <- function(df, output_base) {
  method_label <- unique(df$method)
  object_label <- unique(df$analysis_object)
  if (length(method_label) != 1 || length(object_label) != 1) {
    stop("run_method_object_analysis expects one method/object subset.", call. = FALSE)
  }
  method_dir <- .ensure_dir(file.path(output_base, method_label))
  object_safe <- safe_file_component(object_label)
  object_dir <- .ensure_dir(file.path(method_dir, object_safe))
  write_csv_safe(df, file.path(object_dir, "cell_metrics_with_cell_ploidy.csv"))

  ploidy_res <- run_ploidy_pseudotime_outputs(df, method_dir, method_label, object_label, object_safe)
  dose_res <- NULL
  if (identical(object_label, "Tumor/ploidy_all")) {
    dose_res <- run_dose_pseudotime_outputs(df, method_dir, method_label, object_label, object_safe)
  }
  corr_df <- run_cluster_correlation_outputs(df, method_dir, method_label, object_label, object_safe)
  tumor_res <- NULL
  if (identical(object_label, "Tumor/ploidy_all")) {
    tumor_res <- run_tumor_ploidy_dose_outputs(df, method_dir, method_label, object_safe)
  }
  dist_status <- run_cluster_ploidy_distribution_outputs(df, method_dir, method_label, object_label, object_safe)

  data.frame(
    method = method_label,
    analysis_object = object_label,
    output_dir = object_dir,
    n_cells = nrow(df),
    n_finite_pseudotime = sum(is.finite(df$pseudotime), na.rm = TRUE),
    n_with_annotation_ploidy = sum(df$annotation_ploidy %in% c("2N", "4N"), na.rm = TRUE),
    n_with_dose = sum(!is.na(df$Dose) & nzchar(df$Dose), na.rm = TRUE),
    n_with_cell_ploidy = sum(is.finite(df$cell_ploidy), na.rm = TRUE),
    n_ploidy_test_rows = nrow(ploidy_res$tests),
    n_dose_test_rows = if (is.null(dose_res)) 0L else nrow(dose_res$tests),
    n_cluster_correlation_rows = nrow(corr_df),
    n_tumor_ploidy_dose_rows = if (is.null(tumor_res)) 0L else tumor_res$continuous_rows,
    n_distribution_status_rows = nrow(dist_status),
    stringsAsFactors = FALSE
  )
}

config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
repo_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
results_root <- get_env_scalar("EXTRA_RESULTS_ROOT", get_results_root(config))
output_base <- get_env_scalar("EXTRA_04D_OUTPUT_ROOT", file.path(results_root, "04d_velocity_results", scenario_id))
ploidy_file <- get_env_scalar("EXTRA_PLOIDY_FILE", file.path(repo_root, "Data", "in-vivo", "SUM-159", "all_ploidy.tsv"))
sample_info_file <- get_env_scalar("EXTRA_SAMPLE_INFO_FILE", file.path(repo_root, "Data", "in-vivo", "SUM-159", "sample_info.xlsx"))
dry_run <- is_truthy_env("VELOCITY_RESULTS_DRY_RUN", default = FALSE)

message("04d velocity results analysis")
message("  scenario: ", scenario_id)
message("  results_root: ", normalizePath(results_root, mustWork = FALSE))
message("  output_base: ", normalizePath(output_base, mustWork = FALSE))
message("  ploidy_file: ", normalizePath(ploidy_file, mustWork = FALSE))
message("  sample_info_file: ", normalizePath(sample_info_file, mustWork = FALSE))

.ensure_dir(output_base)
out_inputs <- .ensure_dir(file.path(output_base, "00_inputs"))

manifest <- metric_manifest(results_root)
write_csv_safe(manifest, file.path(out_inputs, "input_metric_manifest.csv"))
if (!any(manifest$file_exists)) {
  stop("No input metric files were found for ", scenario_id, " under: ", results_root, call. = FALSE)
}

sample_info_map <- read_sample_info_map(sample_info_file)
cell_ploidy_map <- read_cell_ploidy_map(ploidy_file)
write_csv_safe(sample_info_map, file.path(out_inputs, "sample_info_harvest_map.csv"))
write_csv_safe(
  cell_ploidy_map %>%
    dplyr::group_by(.data$harvest) %>%
    dplyr::summarise(n_barcodes = dplyr::n(), mean_cell_ploidy = mean(.data$cell_ploidy, na.rm = TRUE), .groups = "drop"),
  file.path(out_inputs, "cell_ploidy_source_summary_by_harvest.csv")
)

metrics <- load_all_metrics(manifest) %>%
  dplyr::filter(is.finite(.data$pseudotime))
if (nrow(metrics) == 0) {
  stop("No finite pseudotime rows were loaded for ", scenario_id, ".", call. = FALSE)
}
metrics <- fill_umap_from_reference(metrics)
metrics <- join_cell_ploidy(metrics, sample_info_map, cell_ploidy_map)
mapping_summary <- summarise_mapping(metrics)
write_csv_safe(mapping_summary, file.path(out_inputs, "cell_ploidy_mapping_summary.csv"))
write_csv_safe(
  metrics %>%
    dplyr::filter(is.na(.data$cell_ploidy)) %>%
    dplyr::select(dplyr::any_of(c("method", "analysis_object", "cell", "sample", "harvest", "barcode_raw", "annotation_ploidy", "Dose", "clusters"))) %>%
    utils::head(5000),
  file.path(out_inputs, "unmapped_cell_ploidy_examples.csv")
)

if (isTRUE(dry_run)) {
  message("VELOCITY_RESULTS_DRY_RUN is true; stopping after input discovery and mapping summaries.")
  quit(save = "no", status = 0)
}

summary_rows <- metrics %>%
  dplyr::group_by(.data$method, .data$analysis_object) %>%
  dplyr::group_split() %>%
  lapply(run_method_object_analysis, output_base = output_base)
summary_df <- dplyr::bind_rows(summary_rows)
write_csv_safe(summary_df, file.path(output_base, "04d_velocity_results_summary.csv"))
cross_method_status <- run_cross_method_comparison_outputs(metrics, output_base)
tumor_ploidy_dose_status <- run_tumor_ploidy_dose_response_outputs(metrics, output_base)
scvelo_deep_dive_status <- run_scvelo_deep_dive_outputs(metrics, output_base, results_root)

run_summary_lines <- c(
  "04d velocity results analysis completed.",
  paste0("Scenario: ", scenario_id),
  paste0("Results root: ", normalizePath(results_root, mustWork = FALSE)),
  paste0("Output root: ", normalizePath(output_base, mustWork = FALSE)),
  paste0("Input method/object tables: ", nrow(manifest)),
  paste0("Loaded finite pseudotime rows: ", nrow(metrics)),
  paste0("Rows with mapped Cell Ploidy: ", sum(is.finite(metrics$cell_ploidy), na.rm = TRUE)),
  "",
  "Main outputs:",
  "  00_inputs/input_metric_manifest.csv",
  "  00_inputs/cell_ploidy_mapping_summary.csv",
  "  <method>/<object>/01_initial_ploidy_pseudotime/*",
  "  <method>/<object>/02_dose_pseudotime/* for Tumor/ploidy_all",
  "  <method>/<object>/03_cluster_cell_ploidy_pseudotime_correlation/*",
  "  <method>/<object>/04_tumor_initial_ploidy_dose_pseudotime/* for Tumor/ploidy_all",
  "  <method>/<object>/05_cluster_cell_ploidy_distribution/*",
  "  06_cross_method_comparison/*",
  "  07_tumor_initial_ploidy_dose_response/*",
  "  08_scvelo_initial_ploidy_cell_ploidy_dose_cluster_deep_dive/*",
  paste0("Cross-method comparison object rows: ", nrow(cross_method_status)),
  paste0("Tumor ploidy-dose response status rows: ", nrow(tumor_ploidy_dose_status)),
  paste0("scVelo deep-dive status rows: ", nrow(scvelo_deep_dive_status))
)
writeLines(run_summary_lines, file.path(output_base, "run_summary.txt"))
message(paste(run_summary_lines, collapse = "\n"))
