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
source(file.path(script_dir, "04c_extra_helpers.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

required_packages <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr")
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
  library(tibble)
  library(tidyr)
})

set.seed(1234)

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

trajectory_root <- get_env_scalar("EXTRA_TRAJECTORY_ROOT", file.path(results_root, "04_trajectory"))
pseudotrajectory_root <- get_env_scalar("EXTRA_PSEUDOTRAJECTORY_ROOT", file.path(results_root, "04a_psudo_rajectory"))
deg_root <- get_env_scalar("EXTRA_DEG_ROOT", file.path(results_root, "03a_DEGs"))
deg_round2_root <- get_env_scalar("EXTRA_DEG_ROUND2_ROOT", file.path(results_root, "03c_DEGs_round2"))
gsea_root <- get_env_scalar("EXTRA_GSEA_ROOT", file.path(results_root, "03b_cluster_annotation_and_GSEA"))
output_root <- get_env_scalar("EXTRA_TRAJECTORY_OUTPUT_ROOT", file.path(results_root, "04b_extra_trajectory"))

method_scvelo <- "scVelo"
method_monocle3 <- "Monocle3"
method_dir_scvelo <- "scVelo"
method_dir_monocle3 <- "Monocle3"

out_inputs <- .ensure_dir(file.path(output_root, "00_inputs"))
out_tables <- .ensure_dir(file.path(output_root, "01_tables"))
out_stats <- .ensure_dir(file.path(output_root, "02_stats"))
out_plots <- .ensure_dir(file.path(output_root, "03_plots"))
out_summary <- .ensure_dir(file.path(output_root, "04_question_summary"))
out_separate <- .ensure_dir(file.path(output_root, "05_method_separate"))
out_scvelo <- .ensure_dir(file.path(out_separate, method_dir_scvelo))
out_monocle3 <- .ensure_dir(file.path(out_separate, method_dir_monocle3))
out_biology <- .ensure_dir(file.path(output_root, "06_deg_gsea_pseudotime"))
out_umap <- .ensure_dir(file.path(output_root, "07_umap_overlays"))
out_umap_scvelo <- .ensure_dir(file.path(out_umap, method_dir_scvelo))
out_umap_monocle3 <- .ensure_dir(file.path(out_umap, method_dir_monocle3))
out_flow <- .ensure_dir(file.path(output_root, "08_flow_overlays"))
out_flow_scvelo <- .ensure_dir(file.path(out_flow, method_dir_scvelo))
out_flow_monocle3 <- .ensure_dir(file.path(out_flow, method_dir_monocle3))
out_tn_ploidy <- .ensure_dir(file.path(output_root, "09_tn_ploidy_cluster_response"))
out_tn_ploidy_scvelo <- .ensure_dir(file.path(out_tn_ploidy, method_dir_scvelo))
out_tn_ploidy_monocle3 <- .ensure_dir(file.path(out_tn_ploidy, method_dir_monocle3))
out_tn_ploidy_flow <- .ensure_dir(file.path(out_tn_ploidy, "scVelo_flow_umap"))
out_cluster_timing <- .ensure_dir(file.path(output_root, "10_cluster_timing"))
out_cluster_timing_scvelo <- .ensure_dir(file.path(out_cluster_timing, method_dir_scvelo))
out_cluster_timing_monocle3 <- .ensure_dir(file.path(out_cluster_timing, method_dir_monocle3))
out_paga_extra <- .ensure_dir(file.path(output_root, "11_paga_extra"))
out_final_summary <- .ensure_dir(file.path(output_root, "summary"))

save_both <- function(plot_obj, file_stub, width = 9, height = 7, dpi = 300) {
  ggplot2::ggsave(paste0(file_stub, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(file_stub, ".png"), plot_obj, width = width, height = height, dpi = dpi)
}

save_square_panel_both <- function(plot_obj, file_stub, panel_size = 4.8, dpi = 300) {
  grob <- ggplot2::ggplotGrob(plot_obj)
  panel_idx <- grepl("^panel", grob$layout$name)
  if (!any(panel_idx)) {
    save_both(plot_obj, file_stub, width = 8, height = 7, dpi = dpi)
    return(invisible(plot_obj))
  }
  panel_cols <- sort(unique(grob$layout$l[panel_idx]))
  panel_rows <- sort(unique(grob$layout$t[panel_idx]))
  grob$widths[panel_cols] <- grid::unit(panel_size, "in")
  grob$heights[panel_rows] <- grid::unit(panel_size, "in")
  width <- tryCatch(
    grid::convertWidth(sum(grob$widths), "in", valueOnly = TRUE),
    error = function(e) NA_real_
  )
  height <- tryCatch(
    grid::convertHeight(sum(grob$heights), "in", valueOnly = TRUE),
    error = function(e) NA_real_
  )
  if (!is.finite(width) || width <= 0) width <- length(panel_cols) * panel_size + 3
  if (!is.finite(height) || height <= 0) height <- length(panel_rows) * panel_size + 2
  ggplot2::ggsave(paste0(file_stub, ".pdf"), grob, width = width, height = height, limitsize = FALSE)
  ggplot2::ggsave(paste0(file_stub, ".png"), grob, width = width, height = height, dpi = dpi, limitsize = FALSE)
  invisible(grob)
}

theme_extra <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey35"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      strip.background = element_rect(fill = "grey94", color = NA),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

square_umap_limits <- function(x, y, pad_fraction = 0.03) {
  x <- safe_numeric(x)
  y <- safe_numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) == 0 || length(y) == 0) return(list(xlim = NULL, ylim = NULL))
  x_range <- range(x, na.rm = TRUE)
  y_range <- range(y, na.rm = TRUE)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  if (!is.finite(x_span) || x_span <= 0) x_span <- 1
  if (!is.finite(y_span) || y_span <= 0) y_span <- 1
  span <- max(x_span, y_span) * (1 + 2 * pad_fraction)
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  list(
    xlim = x_mid + c(-0.5, 0.5) * span,
    ylim = y_mid + c(-0.5, 0.5) * span
  )
}

coord_square_umap <- function(x, y, pad_fraction = 0.03) {
  lim <- square_umap_limits(x, y, pad_fraction = pad_fraction)
  ggplot2::coord_fixed(ratio = 1, xlim = lim$xlim, ylim = lim$ylim, expand = FALSE)
}

theme_square_umap_panel <- function() {
  theme(aspect.ratio = 1)
}

as_clean_chr <- function(x) {
  out <- trimws(as.character(x))
  out[out %in% c("", "NA", "NaN", "NULL", "<NA>")] <- NA_character_
  out
}

coalesce_chr <- function(...) {
  vals <- list(...)
  if (length(vals) == 0) return(character(0))
  out <- as_clean_chr(vals[[1]])
  if (length(vals) > 1) {
    for (i in seq(2, length(vals))) {
      y <- as_clean_chr(vals[[i]])
      if (length(y) == 1 && length(out) != 1) y <- rep(y, length(out))
      replace_idx <- is.na(out) & !is.na(y)
      out[replace_idx] <- y[replace_idx]
    }
  }
  out
}

ensure_columns <- function(df, cols, default = NA_character_) {
  for (col in cols) {
    if (!(col %in% names(df))) df[[col]] <- default
  }
  df
}

standardize_group_labels <- function(x) {
  x <- as_clean_chr(x)
  x <- gsub("mg/kg", "mg/kg", x, fixed = TRUE)
  x
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

median_finite <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
}

quantile_finite <- function(x, prob) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  as.numeric(stats::quantile(x, probs = prob, na.rm = TRUE, names = FALSE))
}

detect_velocity_context <- function(path) {
  path <- normalizePath(path, mustWork = FALSE)
  if (grepl("/01_velocity_groups/", path, fixed = TRUE)) return("velocity_groups")
  if (grepl("/07_by_03a_DEG_comparison/", path, fixed = TRUE)) return("DEG_comparison")
  if (grepl("/01_merged_groups/", path, fixed = TRUE)) return("merged_groups")
  if (grepl("/01_by_dose/", path, fixed = TRUE)) return("by_dose_or_dose_ploidy")
  "other"
}

velocity_source_label <- function(file, search_root) {
  source_dir <- normalizePath(dirname(file), mustWork = FALSE)
  root_dir <- normalizePath(search_root, mustWork = FALSE)
  prefix <- paste0(root_dir, .Platform$file.sep)
  if (startsWith(source_dir, prefix)) {
    rel <- substring(source_dir, nchar(prefix) + 1L)
    rel <- gsub(.Platform$file.sep, "/", rel, fixed = TRUE)
    if (nzchar(rel)) return(rel)
  }
  basename(source_dir)
}

read_csv_selected <- function(file, keep_cols) {
  df <- readr::read_csv(
    file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- ensure_columns(df, keep_cols)
  df[, keep_cols, drop = FALSE]
}

metadata_cols <- c(
  "cell",
  "UMAP_1",
  "UMAP_2",
  "Dose",
  "ploidy",
  "Ploidy",
  "Karyotype",
  "TN",
  "Dose_DEG",
  "trajectory_context",
  "trajectory_group",
  "trajectory_analysis_group",
  "trajectory_shape_group",
  "trajectory_tn_scope",
  "trajectory_ploidy_scope",
  "trajectory_branch",
  "AnalysisGroup",
  "AnalysisGroupColumn",
  "Analysis",
  "TrajectoryComparisonID",
  "TrajectoryScope",
  "TrajectoryComparison",
  "TrajectoryIdent1",
  "TrajectoryIdent2",
  "TrajectoryGroupCol",
  "TrajectoryCluster",
  "TrajectorySubsetGroup",
  "TrajectoryDose",
  "TrajectoryComparisonGroup",
  "TrajectoryOriginalGroup",
  "TrajectoryUniverseKey",
  "clusters",
  "cluster",
  "cluster_final",
  "sample",
  "sampleID",
  "sample_id",
  "IDs",
  "orig.ident",
  "Sequencing.IDs",
  "sample_type",
  "velocity_batch"
)

velocity_metric_cols <- c(
  "velocity_pseudotime",
  "velocity_length",
  "velocity_confidence",
  "velocity_confidence_transition",
  "root_cells_manual",
  "end_points_manual"
)
pseudo_metric_cols <- c("pseudotime")

harmonize_common_metadata <- function(df) {
  n <- nrow(df)
  if (!("Ploidy" %in% names(df))) df$Ploidy <- NA_character_
  if (!("ploidy" %in% names(df))) df$ploidy <- NA_character_
  if (!("Karyotype" %in% names(df))) df$Karyotype <- NA_character_
  df$Ploidy <- coalesce_chr(df$Ploidy, df$ploidy, df$Karyotype)
  df$ploidy <- coalesce_chr(df$ploidy, df$Ploidy, df$Karyotype)
  if (!("sample" %in% names(df))) df$sample <- NA_character_
  if (!("sampleID" %in% names(df))) df$sampleID <- NA_character_
  if (!("sample_id" %in% names(df))) df$sample_id <- NA_character_
  if (!("IDs" %in% names(df))) df$IDs <- NA_character_
  if (!("orig.ident" %in% names(df))) df[["orig.ident"]] <- NA_character_
  if (!("Sequencing.IDs" %in% names(df))) df[["Sequencing.IDs"]] <- NA_character_
  df$sampleID <- coalesce_chr(df$sampleID, df$sample_id, df$sample, df$IDs, df[["orig.ident"]], df[["Sequencing.IDs"]])
  df$sample <- coalesce_chr(df$sample, df$sampleID, df$sample_id, df$IDs, df[["orig.ident"]], df[["Sequencing.IDs"]])

  if (!("clusters" %in% names(df))) df$clusters <- NA_character_
  if (!("cluster" %in% names(df))) df$cluster <- NA_character_
  if (!("cluster_final" %in% names(df))) df$cluster_final <- NA_character_
  df$clusters <- coalesce_chr(df$clusters, df$cluster, df$cluster_final)

  for (col in c("UMAP_1", "UMAP_2")) {
    if (!(col %in% names(df))) df[[col]] <- NA_real_
    df[[col]] <- safe_numeric(df[[col]])
  }

  for (col in c(
    "Dose", "Dose_DEG", "TN", "trajectory_context", "trajectory_group",
    "trajectory_analysis_group", "trajectory_shape_group", "trajectory_tn_scope",
    "trajectory_ploidy_scope", "trajectory_branch", "AnalysisGroup",
    "AnalysisGroupColumn", "Analysis", "TrajectoryComparisonID",
    "TrajectoryScope", "TrajectoryComparison", "TrajectoryIdent1",
    "TrajectoryIdent2", "TrajectoryGroupCol", "TrajectoryCluster",
    "TrajectorySubsetGroup", "TrajectoryDose", "TrajectoryComparisonGroup",
    "TrajectoryOriginalGroup", "TrajectoryUniverseKey", "sample",
    "sampleID", "sample_id", "IDs", "orig.ident", "Sequencing.IDs",
    "sample_type", "velocity_batch"
  )) {
    if (!(col %in% names(df))) df[[col]] <- rep(NA_character_, n)
    df[[col]] <- as_clean_chr(df[[col]])
  }

  df$Dose <- factor(as.character(df$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
  df$Dose_DEG <- factor(as.character(df$Dose_DEG), levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
  df$Ploidy <- factor(as.character(df$Ploidy), levels = c("2N", "4N"))
  df$ploidy <- factor(as.character(df$ploidy), levels = c("2N", "4N"))
  df$clusters <- factor(as.character(df$clusters), levels = sort_maybe_numeric(as.character(df$clusters)))
  df$TrajectoryComparisonGroup <- standardize_group_labels(df$TrajectoryComparisonGroup)
  df
}

merge_same_dir_umap <- function(df, source_dir) {
  umap_file <- file.path(source_dir, "cells_metadata_umap.csv")
  if (!file.exists(umap_file) || file.info(umap_file)$size <= 0) return(df)

  umap_cols <- unique(c(
    "cell", "UMAP_1", "UMAP_2", "Dose", "ploidy", "Ploidy", "TN",
    "Dose_DEG", "trajectory_context", "trajectory_group",
    "trajectory_analysis_group", "trajectory_shape_group",
    "trajectory_tn_scope", "trajectory_ploidy_scope", "trajectory_branch",
    "clusters", "sample", "sample_type"
  ))
  umap_df <- read_csv_selected(umap_file, umap_cols) %>%
    dplyr::mutate(cell = as_clean_chr(.data$cell)) %>%
    dplyr::filter(!is.na(.data$cell)) %>%
    dplyr::distinct(.data$cell, .keep_all = TRUE)
  if (nrow(umap_df) == 0) return(df)

  joined <- dplyr::left_join(df, umap_df, by = "cell", suffix = c("", ".umap"))
  for (col in setdiff(umap_cols, "cell")) {
    umap_col <- paste0(col, ".umap")
    if (umap_col %in% names(joined)) {
      joined[[col]] <- coalesce_chr(joined[[col]], joined[[umap_col]])
      joined[[umap_col]] <- NULL
    }
  }
  joined
}

read_velocity_metrics <- function(root_dir) {
  velocity_group_root <- file.path(root_dir, "01_velocity_groups")
  search_root <- if (dir.exists(velocity_group_root)) velocity_group_root else root_dir
  files <- list.files(
    search_root,
    pattern = "^scvelo_cell_metrics\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  files <- files[file.info(files)$size > 0]
  if (length(files) == 0) {
    warning("No scvelo_cell_metrics.csv files found under ", search_root, call. = FALSE)
    return(data.frame())
  }

  keep_cols <- unique(c(metadata_cols, velocity_metric_cols))
  out <- vector("list", length(files))
  manifest <- data.frame(
    source_file = files,
    source_dir = dirname(files),
    source_context = vapply(files, detect_velocity_context, character(1)),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(files)) {
    message("Reading velocity metrics: ", i, "/", length(files), " ", files[i])
    df <- read_csv_selected(files[i], keep_cols)
    df <- merge_same_dir_umap(df, dirname(files[i]))
    for (col in intersect(metadata_cols, names(df))) df[[col]] <- as.character(df[[col]])
    for (col in intersect(velocity_metric_cols, names(df))) df[[col]] <- safe_numeric(df[[col]])
    df$source_file <- files[i]
    df$source_dir <- dirname(files[i])
    df$source_context <- manifest$source_context[i]
    df$source_label <- velocity_source_label(files[i], search_root)
    df$method <- method_scvelo
    df$metric_family <- "velocity"
    out[[i]] <- df
  }

  write_table_csv(manifest, file.path(out_inputs, "scVelo_metric_files.csv"))
  df_all <- dplyr::bind_rows(out)
  df_all <- harmonize_common_metadata(df_all)
  for (metric in velocity_metric_cols) df_all[[metric]] <- safe_numeric(df_all[[metric]])
  df_all
}

read_pseudotrajectory_metrics <- function(root_dir) {
  summary_dir <- file.path(root_dir, "02_summary")
  files <- c(
    pseudotrajectory_groups = file.path(summary_dir, "pseudotime_cells_all_pseudotrajectory_groups.csv"),
    by_dose = file.path(summary_dir, "pseudotime_cells_all_doses.csv"),
    by_dose_karyotype = file.path(summary_dir, "pseudotime_cells_all_dose_karyotype.csv"),
    cellline_0mg_tumor_independent = file.path(summary_dir, "pseudotime_cells_all_cellline_0mg_tumor_independent.csv"),
    cellline_0mg_tumor_combined = file.path(summary_dir, "pseudotime_cells_all_cellline_0mg_tumor_combined.csv"),
    DEG_comparison = file.path(summary_dir, "pseudotime_cells_all_DEG_comparisons.csv")
  )
  files <- files[file.exists(files) & file.info(files)$size > 0]
  if (length(files) == 0) {
    warning("No pseudotime summary files found under ", summary_dir, call. = FALSE)
    return(data.frame())
  }

  keep_cols <- unique(c(metadata_cols, pseudo_metric_cols, "pseudotime_raw"))
  out <- vector("list", length(files))
  manifest <- data.frame(
    source_dataset = names(files),
    source_file = unname(files),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(files)) {
    message("Reading Monocle3 metrics: ", i, "/", length(files), " ", files[i])
    df <- read_csv_selected(files[i], keep_cols)
    for (col in intersect(metadata_cols, names(df))) df[[col]] <- as.character(df[[col]])
    for (col in intersect(c(pseudo_metric_cols, "pseudotime_raw"), names(df))) df[[col]] <- safe_numeric(df[[col]])
    df$source_file <- unname(files[i])
    df$source_dir <- dirname(unname(files[i]))
    source_context <- dplyr::recode(names(files)[i], pseudotrajectory_groups = "monocle3_groups", .default = names(files)[i])
    df$source_context <- source_context
    df$source_label <- source_context
    df$method <- method_monocle3
    df$metric_family <- "pseudotime"
    out[[i]] <- df
  }

  write_table_csv(manifest, file.path(out_inputs, "Monocle3_metric_files.csv"))
  df_all <- dplyr::bind_rows(out)
  df_all <- harmonize_common_metadata(df_all)
  df_all$pseudotime <- safe_numeric(df_all$pseudotime)
  df_all$pseudotime_raw <- safe_numeric(df_all$pseudotime_raw)
  df_all
}

group_summary_table <- function(df, split_cols, group_col, metric_col, method, analysis_family, level) {
  needed <- unique(c(split_cols, group_col, metric_col))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .metric = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(.group), is.finite(.metric))
  if (nrow(df) == 0) return(data.frame())

  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_cols)), .group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean = mean(.metric, na.rm = TRUE),
      median = stats::median(.metric, na.rm = TRUE),
      q25 = stats::quantile(.metric, probs = 0.25, na.rm = TRUE, names = FALSE),
      q75 = stats::quantile(.metric, probs = 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    dplyr::rename(group = .group) %>%
    dplyr::mutate(
      method = method,
      metric = metric_col,
      analysis_family = analysis_family,
      analysis_level = level,
      group_col = group_col
    )
}

run_one_test <- function(df, group_col, metric_col, min_group_n = 3L, preferred_groups = NULL) {
  df <- df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .metric = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(.group), is.finite(.metric))

  group_order <- unique(c(as_clean_chr(preferred_groups), df$.group))
  group_order <- group_order[!is.na(group_order)]
  group_counts <- table(factor(df$.group, levels = group_order))
  keep_groups <- names(group_counts)[group_counts >= min_group_n]
  df <- df[df$.group %in% keep_groups, , drop = FALSE]
  group_order <- group_order[group_order %in% keep_groups]
  group_counts <- table(factor(df$.group, levels = group_order))
  groups <- names(group_counts)

  if (length(groups) < 2) {
    return(data.frame(
      n_total = nrow(df),
      n_groups = length(groups),
      tested_groups = paste(groups, collapse = ";"),
      primary_test = "not_tested",
      primary_p = NA_real_,
      ks_p = NA_real_,
      group_1 = ifelse(length(groups) >= 1, groups[1], NA_character_),
      group_2 = ifelse(length(groups) >= 2, groups[2], NA_character_),
      n_group_1 = ifelse(length(groups) >= 1, as.integer(group_counts[[groups[1]]]), NA_integer_),
      n_group_2 = ifelse(length(groups) >= 2, as.integer(group_counts[[groups[2]]]), NA_integer_),
      median_group_1 = NA_real_,
      median_group_2 = NA_real_,
      median_diff_group_1_minus_group_2 = NA_real_,
      median_range = NA_real_,
      top_median_group = NA_character_,
      bottom_median_group = NA_character_,
      rank_biserial = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  medians <- tapply(df$.metric, df$.group, stats::median, na.rm = TRUE)
  medians <- medians[groups]
  top_group <- names(medians)[which.max(medians)]
  bottom_group <- names(medians)[which.min(medians)]

  if (length(groups) == 2) {
    x <- df$.metric[df$.group == groups[1]]
    y <- df$.metric[df$.group == groups[2]]
    wt <- tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
      error = function(e) NULL
    )
    kt <- tryCatch(
      suppressWarnings(stats::ks.test(x, y, exact = FALSE)),
      error = function(e) NULL
    )
    u_stat <- if (!is.null(wt) && length(wt$statistic) > 0) as.numeric(wt$statistic) else NA_real_
    rb <- if (is.finite(u_stat) && length(x) > 0 && length(y) > 0) {
      (2 * u_stat / (length(x) * length(y))) - 1
    } else {
      NA_real_
    }
    return(data.frame(
      n_total = nrow(df),
      n_groups = 2L,
      tested_groups = paste(groups, collapse = ";"),
      primary_test = "Wilcoxon_rank_sum",
      primary_p = if (!is.null(wt)) wt$p.value else NA_real_,
      ks_p = if (!is.null(kt)) kt$p.value else NA_real_,
      group_1 = groups[1],
      group_2 = groups[2],
      n_group_1 = as.integer(group_counts[[groups[1]]]),
      n_group_2 = as.integer(group_counts[[groups[2]]]),
      median_group_1 = unname(medians[[groups[1]]]),
      median_group_2 = unname(medians[[groups[2]]]),
      median_diff_group_1_minus_group_2 = unname(medians[[groups[1]]] - medians[[groups[2]]]),
      median_range = unname(max(medians, na.rm = TRUE) - min(medians, na.rm = TRUE)),
      top_median_group = top_group,
      bottom_median_group = bottom_group,
      rank_biserial = rb,
      stringsAsFactors = FALSE
    ))
  }

  kt <- tryCatch(
    stats::kruskal.test(.metric ~ .group, data = df),
    error = function(e) NULL
  )
  data.frame(
    n_total = nrow(df),
    n_groups = length(groups),
    tested_groups = paste(groups, collapse = ";"),
    primary_test = "Kruskal_Wallis",
    primary_p = if (!is.null(kt)) kt$p.value else NA_real_,
    ks_p = NA_real_,
    group_1 = NA_character_,
    group_2 = NA_character_,
    n_group_1 = NA_integer_,
    n_group_2 = NA_integer_,
    median_group_1 = NA_real_,
    median_group_2 = NA_real_,
    median_diff_group_1_minus_group_2 = NA_real_,
    median_range = unname(max(medians, na.rm = TRUE) - min(medians, na.rm = TRUE)),
    top_median_group = top_group,
    bottom_median_group = bottom_group,
    rank_biserial = NA_real_,
    stringsAsFactors = FALSE
  )
}

run_tests_by_split <- function(df, split_cols, group_col, metric_col, method, analysis_family, level) {
  needed <- unique(c(split_cols, group_col, metric_col))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .metric = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(.group), is.finite(.metric))
  if (nrow(df) == 0) return(data.frame())

  if (length(split_cols) == 0) {
    min_group_n <- if (identical(level, "sample")) 2L else 3L
    test_df <- run_one_test(df, group_col, metric_col, min_group_n = min_group_n)
    return(test_df %>%
      dplyr::mutate(
        method = method,
        metric = metric_col,
        analysis_family = analysis_family,
        analysis_level = level,
        group_col = group_col
      ))
  }

  split_df <- df %>%
    dplyr::mutate(.split_key = apply(dplyr::across(dplyr::all_of(split_cols)), 1, paste, collapse = "\r"))
  split_keys <- unique(split_df$.split_key)
  out <- vector("list", length(split_keys))
  for (i in seq_along(split_keys)) {
    sub_df <- split_df[split_df$.split_key == split_keys[i], , drop = FALSE]
    min_group_n <- if (identical(level, "sample")) 2L else 3L
    preferred_groups <- NULL
    if (
      identical(group_col, "TrajectoryComparisonGroup") &&
        all(c("TrajectoryIdent1", "TrajectoryIdent2") %in% names(sub_df))
    ) {
      preferred_groups <- c(sub_df$TrajectoryIdent1[1], sub_df$TrajectoryIdent2[1])
    }
    test_df <- run_one_test(
      sub_df,
      group_col,
      metric_col,
      min_group_n = min_group_n,
      preferred_groups = preferred_groups
    )
    split_values <- sub_df[1, split_cols, drop = FALSE]
    out[[i]] <- cbind(split_values, test_df, stringsAsFactors = FALSE)
  }

  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      method = method,
      metric = metric_col,
      analysis_family = analysis_family,
      analysis_level = level,
      group_col = group_col
    )
}

make_sample_level <- function(df, split_cols, group_col, metric_col) {
  sample_col <- if ("sample" %in% names(df)) "sample" else NA_character_
  if (is.na(sample_col)) return(data.frame())
  needed <- unique(c(split_cols, group_col, sample_col, metric_col))
  df <- ensure_columns(df, needed)
  df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .sample = as_clean_chr(.data[[sample_col]]),
      .metric = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(.group), !is.na(.sample), is.finite(.metric)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_cols)), .group, .sample) %>%
    dplyr::summarise(
      sample_median = stats::median(.metric, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(!!group_col := .group, sample = .sample) %>%
    dplyr::rename(!!metric_col := sample_median)
}

add_adjusted_p <- function(df) {
  if (nrow(df) == 0 || !("primary_p" %in% names(df))) return(df)
  df %>%
    dplyr::group_by(.data$method, .data$metric, .data$analysis_level) %>%
    dplyr::mutate(primary_p_adj_global = stats::p.adjust(primary_p, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(.data$method, .data$metric, .data$analysis_level, .data$analysis_family) %>%
    dplyr::mutate(primary_p_adj_family = stats::p.adjust(primary_p, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      signed_log10_fdr = dplyr::case_when(
        is.finite(primary_p_adj_family) & is.finite(median_diff_group_1_minus_group_2) ~
          sign(median_diff_group_1_minus_group_2) * -log10(pmax(primary_p_adj_family, .Machine$double.xmin)),
        is.finite(primary_p_adj_family) & is.finite(median_range) ~
          -log10(pmax(primary_p_adj_family, .Machine$double.xmin)),
        TRUE ~ NA_real_
      ),
      significant_family_fdr_0_05 = is.finite(primary_p_adj_family) & primary_p_adj_family < 0.05,
      significant_global_fdr_0_05 = is.finite(primary_p_adj_global) & primary_p_adj_global < 0.05
    )
}

run_comparison_suite <- function(df, metric_cols, method) {
  tests <- list()
  summaries <- list()

  deg_df <- df %>%
    dplyr::filter(.data$source_context == "DEG_comparison", !is.na(.data$TrajectoryComparisonGroup))
  deg_split_cols <- c(
    "TrajectoryComparisonID",
    "TrajectoryScope",
    "TrajectoryComparison",
    "TrajectoryIdent1",
    "TrajectoryIdent2",
    "TrajectoryGroupCol",
    "TrajectoryCluster",
    "TrajectorySubsetGroup",
    "TrajectoryDose"
  )

  for (metric in metric_cols) {
    if (!(metric %in% names(df))) next
    tests[[paste("deg", metric, sep = "_")]] <- run_tests_by_split(
      deg_df,
      split_cols = deg_split_cols,
      group_col = "TrajectoryComparisonGroup",
      metric_col = metric,
      method = method,
      analysis_family = "03a_DEG_comparison_aligned",
      level = "cell"
    )
    summaries[[paste("deg", metric, sep = "_")]] <- group_summary_table(
      deg_df,
      split_cols = deg_split_cols,
      group_col = "TrajectoryComparisonGroup",
      metric_col = metric,
      method = method,
      analysis_family = "03a_DEG_comparison_aligned",
      level = "cell"
    )

    sample_deg <- make_sample_level(deg_df, deg_split_cols, "TrajectoryComparisonGroup", metric)
    if (nrow(sample_deg) > 0) {
      tests[[paste("deg_sample", metric, sep = "_")]] <- run_tests_by_split(
        sample_deg,
        split_cols = deg_split_cols,
        group_col = "TrajectoryComparisonGroup",
        metric_col = metric,
        method = method,
        analysis_family = "03a_DEG_comparison_aligned",
        level = "sample"
      )
      summaries[[paste("deg_sample", metric, sep = "_")]] <- group_summary_table(
        sample_deg,
        split_cols = deg_split_cols,
        group_col = "TrajectoryComparisonGroup",
        metric_col = metric,
        method = method,
        analysis_family = "03a_DEG_comparison_aligned",
        level = "sample"
      )
    }
  }

  standard_df <- df %>%
    dplyr::filter(!(.data$source_context == "DEG_comparison"))

  for (metric in metric_cols) {
    if (!(metric %in% names(df))) next

    cluster_split <- c("source_context", "trajectory_analysis_group", "AnalysisGroup", "source_label")
    cluster_df <- standard_df %>%
      dplyr::mutate(
        trajectory_analysis_group = coalesce_chr(.data$trajectory_analysis_group, .data$AnalysisGroup, .data$source_label),
        AnalysisGroup = coalesce_chr(.data$AnalysisGroup, .data$trajectory_analysis_group, .data$source_label)
      ) %>%
      dplyr::filter(!is.na(.data$clusters))

    tests[[paste("cluster_omnibus", metric, sep = "_")]] <- run_tests_by_split(
      cluster_df,
      split_cols = cluster_split,
      group_col = "clusters",
      metric_col = metric,
      method = method,
      analysis_family = "cluster_omnibus_by_method_run",
      level = "cell"
    )
    summaries[[paste("cluster_omnibus", metric, sep = "_")]] <- group_summary_table(
      cluster_df,
      split_cols = cluster_split,
      group_col = "clusters",
      metric_col = metric,
      method = method,
      analysis_family = "cluster_omnibus_by_method_run",
      level = "cell"
    )

    if (method == method_scvelo) {
      merged_df <- df %>%
        dplyr::filter(.data$source_context == "merged_groups") %>%
        dplyr::mutate(
          trajectory_analysis_group = coalesce_chr(.data$trajectory_analysis_group, .data$source_label)
        )
      merged_split <- c("source_context", "trajectory_analysis_group", "source_label")

      tests[[paste("merged_ploidy", metric, sep = "_")]] <- run_tests_by_split(
        merged_df,
        split_cols = merged_split,
        group_col = "Ploidy",
        metric_col = metric,
        method = method,
        analysis_family = "scVelo_merged_Ploidy",
        level = "cell"
      )
      summaries[[paste("merged_ploidy", metric, sep = "_")]] <- group_summary_table(
        merged_df,
        split_cols = merged_split,
        group_col = "Ploidy",
        metric_col = metric,
        method = method,
        analysis_family = "scVelo_merged_Ploidy",
        level = "cell"
      )

      tests[[paste("merged_dose", metric, sep = "_")]] <- run_tests_by_split(
        merged_df,
        split_cols = merged_split,
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "scVelo_merged_Dose",
        level = "cell"
      )
      summaries[[paste("merged_dose", metric, sep = "_")]] <- group_summary_table(
        merged_df,
        split_cols = merged_split,
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "scVelo_merged_Dose",
        level = "cell"
      )
    }

    if (method == method_monocle3) {
      global_dose_universe <- df %>%
        dplyr::filter(
          .data$source_context == "DEG_comparison",
          .data$TrajectoryScope == "global_Dose_vs_rest",
          !is.na(.data$TrajectoryUniverseKey)
        ) %>%
        dplyr::arrange(.data$TrajectoryComparisonID) %>%
        dplyr::distinct(.data$TrajectoryUniverseKey, .data$cell, .keep_all = TRUE)
      tests[[paste("pseudo_dose_universe", metric, sep = "_")]] <- run_tests_by_split(
        global_dose_universe,
        split_cols = c("TrajectoryUniverseKey"),
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "Monocle3_global_Dose_universe",
        level = "cell"
      )
      summaries[[paste("pseudo_dose_universe", metric, sep = "_")]] <- group_summary_table(
        global_dose_universe,
        split_cols = c("TrajectoryUniverseKey"),
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "Monocle3_global_Dose_universe",
        level = "cell"
      )

      tests[[paste("pseudo_dose_by_ploidy_universe", metric, sep = "_")]] <- run_tests_by_split(
        global_dose_universe,
        split_cols = c("TrajectoryUniverseKey", "Ploidy"),
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "Monocle3_Dose_within_Ploidy_universe",
        level = "cell"
      )
      summaries[[paste("pseudo_dose_by_ploidy_universe", metric, sep = "_")]] <- group_summary_table(
        global_dose_universe,
        split_cols = c("TrajectoryUniverseKey", "Ploidy"),
        group_col = "Dose",
        metric_col = metric,
        method = method,
        analysis_family = "Monocle3_Dose_within_Ploidy_universe",
        level = "cell"
      )
    }
  }

  list(
    tests = add_adjusted_p(dplyr::bind_rows(tests)),
    summaries = dplyr::bind_rows(summaries)
  )
}

scope_question_map <- data.frame(
  question_id = c(
    "cluster_overall",
    "ploidy_overall",
    "dose_overall",
    "tn_overall",
    "within_cluster_ploidy",
    "within_cluster_dose",
    "within_cluster_dose_ploidy",
    "within_cluster_cellline_ploidy",
    "within_cluster_tumor_ploidy"
  ),
  question = c(
    "Different clusters vs rest differ by scVelo or Monocle3 pseudotime.",
    "Ploidy 4N vs 2N differs overall.",
    "Dose levels differ overall by dose-vs-rest contrasts.",
    "Tumor vs CellLine differs overall.",
    "Within each cluster, Ploidy 4N vs 2N differs.",
    "Within each cluster, Dose levels differ by dose-vs-rest contrasts.",
    "Within each cluster and Dose stratum, Ploidy 4N vs 2N differs.",
    "Within each CellLine cluster, Ploidy 4N vs 2N differs.",
    "Within each Tumor cluster, Ploidy 4N vs 2N differs."
  ),
  TrajectoryScope = c(
    "clusters_vs_rest",
    "global_Ploidy_4N_vs_2N",
    "global_Dose_vs_rest",
    "global_TN_Tumor_vs_CellLine",
    "within_cluster_Ploidy_4N_vs_2N",
    "within_cluster_Dose_vs_rest",
    "within_cluster_Dose_Ploidy_4N_vs_2N",
    "within_cluster_CellLine_only_Ploidy_4N_vs_2N",
    "within_cluster_Tumor_only_Ploidy_4N_vs_2N"
  ),
  stringsAsFactors = FALSE
)

summarize_questions <- function(stat_df) {
  time_metric_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$metric %in% c("velocity_pseudotime", "pseudotime")
    )
  if (nrow(time_metric_df) == 0) return(data.frame())

  out <- list()
  idx <- 1L
  for (i in seq_len(nrow(scope_question_map))) {
    scope_i <- scope_question_map$TrajectoryScope[i]
    sub <- time_metric_df %>% dplyr::filter(.data$TrajectoryScope == scope_i)
    if (nrow(sub) == 0) next
    keys <- unique(sub[, c("method", "metric", "analysis_level"), drop = FALSE])
    for (j in seq_len(nrow(keys))) {
      sub2 <- sub %>%
        dplyr::filter(
          .data$method == keys$method[j],
          .data$metric == keys$metric[j],
          .data$analysis_level == keys$analysis_level[j]
        )
      sig <- sub2 %>% dplyr::filter(.data$significant_family_fdr_0_05)
      strongest <- sub2 %>%
        dplyr::arrange(.data$primary_p_adj_family) %>%
        dplyr::slice(1)
      finite_fdr <- sub2$primary_p_adj_family[is.finite(sub2$primary_p_adj_family)]
      out[[idx]] <- data.frame(
        question_id = scope_question_map$question_id[i],
        question = scope_question_map$question[i],
        TrajectoryScope = scope_i,
        method = keys$method[j],
        metric = keys$metric[j],
        analysis_level = keys$analysis_level[j],
        n_tests = nrow(sub2),
        n_significant_fdr_0_05 = nrow(sig),
        fraction_significant = nrow(sig) / nrow(sub2),
        min_fdr = if (length(finite_fdr) > 0) min(finite_fdr, na.rm = TRUE) else NA_real_,
        strongest_comparison = if (nrow(strongest) > 0) strongest$TrajectoryComparison[1] else NA_character_,
        strongest_group_1 = if (nrow(strongest) > 0) strongest$group_1[1] else NA_character_,
        strongest_group_2 = if (nrow(strongest) > 0) strongest$group_2[1] else NA_character_,
        strongest_median_diff_group_1_minus_group_2 = if (nrow(strongest) > 0) strongest$median_diff_group_1_minus_group_2[1] else NA_real_,
        answer = if (nrow(sig) > 0) "yes_at_FDR_0.05" else "not_detected_at_FDR_0.05",
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

plot_scope_summary <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      !is.na(.data$TrajectoryScope)
    ) %>%
    dplyr::group_by(.data$method, .data$metric, .data$TrajectoryScope) %>%
    dplyr::summarise(
      n_tests = dplyr::n(),
      n_significant = sum(.data$significant_family_fdr_0_05, na.rm = TRUE),
      fraction_significant = n_significant / n_tests,
      .groups = "drop"
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = TrajectoryScope, y = fraction_significant, fill = method)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    geom_text(
      aes(label = paste0(n_significant, "/", n_tests)),
      position = position_dodge(width = 0.75),
      vjust = -0.25,
      size = 3
    ) +
    facet_wrap(~metric, nrow = 1) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
    labs(
      title = "Significant pseudotime differences by 03a/03b scope",
      x = NULL,
      y = "Fraction significant at within-family FDR < 0.05",
      fill = "Method"
    ) +
    theme_extra(10)
  save_both(p, file.path(out_plots, "scope_significance_fraction_time_metrics"), width = 13, height = 6)
}

plot_deg_heatmap <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$primary_p_adj_family)
    ) %>%
    dplyr::mutate(
      comparison_label = paste0(
        .data$TrajectoryComparison,
        ifelse(!is.na(.data$TrajectoryCluster), paste0(" | cluster ", .data$TrajectoryCluster), ""),
        ifelse(!is.na(.data$TrajectoryDose), paste0(" | ", .data$TrajectoryDose), "")
      ),
      neg_log10_fdr = -log10(pmax(.data$primary_p_adj_family, .Machine$double.xmin)),
      direction = dplyr::case_when(
        is.finite(.data$median_diff_group_1_minus_group_2) & .data$median_diff_group_1_minus_group_2 > 0 ~ paste0(.data$group_1, " higher"),
        is.finite(.data$median_diff_group_1_minus_group_2) & .data$median_diff_group_1_minus_group_2 < 0 ~ paste0(.data$group_2, " higher"),
        TRUE ~ "multi-group/range"
      )
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df <- plot_df %>%
    dplyr::group_by(.data$TrajectoryScope) %>%
    dplyr::mutate(comparison_label = factor(.data$comparison_label, levels = unique(.data$comparison_label[order(.data$neg_log10_fdr)]))) %>%
    dplyr::ungroup()

  p <- ggplot(plot_df, aes(x = method, y = comparison_label, fill = signed_log10_fdr)) +
    geom_tile(color = "white", linewidth = 0.2, width = 0.95, height = 0.95) +
    facet_grid(TrajectoryScope ~ metric, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Signed -log10 FDR"
    ) +
    labs(
      title = "03a/03b-aligned scVelo and Monocle3 pseudotime differences",
      subtitle = "Positive values mean group_1 median is higher; negative values mean group_2 median is higher.",
      x = NULL,
      y = NULL
    ) +
    theme_extra(8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_both(p, file.path(out_plots, "deg_aligned_time_metric_signed_fdr_heatmap"), width = 13, height = 24)
}

plot_cluster_vs_rest <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$TrajectoryScope == "clusters_vs_rest",
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$median_diff_group_1_minus_group_2)
    ) %>%
    dplyr::mutate(
      cluster = factor(as.character(.data$TrajectoryCluster), levels = sort_maybe_numeric(as.character(.data$TrajectoryCluster))),
      significant = .data$significant_family_fdr_0_05
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = cluster, y = median_diff_group_1_minus_group_2, color = method, shape = significant)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
    geom_point(size = 2.6, position = position_dodge(width = 0.45)) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Cluster vs rest time-metric differences",
      x = "Cluster",
      y = "Median difference: cluster minus rest",
      color = "Method",
      shape = "FDR < 0.05"
    ) +
    theme_extra(11)
  save_both(p, file.path(out_plots, "cluster_vs_rest_time_metric_effects"), width = 10, height = 5.8)
}

plot_global_effects <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$TrajectoryScope %in% c("global_Ploidy_4N_vs_2N", "global_Dose_vs_rest", "global_TN_Tumor_vs_CellLine"),
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$median_diff_group_1_minus_group_2)
    ) %>%
    dplyr::mutate(
      comparison_label = paste0(.data$TrajectoryScope, "\n", .data$TrajectoryComparison),
      significant = .data$significant_family_fdr_0_05
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = comparison_label, y = median_diff_group_1_minus_group_2, color = method, shape = significant)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey55") +
    geom_point(size = 2.7, position = position_dodge(width = 0.45)) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Global Ploidy, Dose, and TN pseudotime effects",
      x = NULL,
      y = "Median difference: group_1 minus group_2",
      color = "Method",
      shape = "FDR < 0.05"
    ) +
    theme_extra(10)
  save_both(p, file.path(out_plots, "global_ploidy_dose_tn_time_metric_effects"), width = 12, height = 6.5)
}

plot_within_cluster_ploidy <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$TrajectoryScope %in% c(
        "within_cluster_Ploidy_4N_vs_2N",
        "within_cluster_Dose_Ploidy_4N_vs_2N",
        "within_cluster_CellLine_only_Ploidy_4N_vs_2N",
        "within_cluster_Tumor_only_Ploidy_4N_vs_2N"
      ),
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$median_diff_group_1_minus_group_2)
    ) %>%
    dplyr::mutate(
      cluster = factor(as.character(.data$TrajectoryCluster), levels = sort_maybe_numeric(as.character(.data$TrajectoryCluster))),
      stratum = dplyr::case_when(
        !is.na(.data$TrajectoryDose) ~ paste0("Dose ", .data$TrajectoryDose),
        !is.na(.data$TrajectorySubsetGroup) ~ .data$TrajectorySubsetGroup,
        TRUE ~ "All"
      )
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = cluster, y = stratum, fill = signed_log10_fdr)) +
    geom_tile(color = "white", linewidth = 0.2, width = 0.95, height = 0.95) +
    facet_grid(TrajectoryScope + metric ~ method, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Signed -log10 FDR"
    ) +
    labs(
      title = "Within-cluster Ploidy differences",
      subtitle = "Positive values mean group_1 median is higher; negative values mean group_2 median is higher.",
      x = "Cluster",
      y = NULL
    ) +
    theme_extra(8)
  save_both(p, file.path(out_plots, "within_cluster_ploidy_time_metric_signed_fdr_heatmap"), width = 13, height = 14)
}

plot_omnibus_cluster <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "cluster_omnibus_by_method_run",
      .data$analysis_level == "cell",
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$primary_p_adj_family)
    ) %>%
    dplyr::mutate(
      analysis_label = coalesce_chr(.data$trajectory_analysis_group, .data$AnalysisGroup, .data$source_label),
      neg_log10_fdr = -log10(pmax(.data$primary_p_adj_family, .Machine$double.xmin))
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = method, y = analysis_label, fill = neg_log10_fdr)) +
    geom_tile(color = "white", linewidth = 0.2, width = 0.95, height = 0.95) +
    facet_grid(source_context + metric ~ ., scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "white", high = "#7B3294", name = "-log10 FDR") +
    labs(
      title = "Cluster omnibus differences within method runs",
      x = NULL,
      y = NULL
    ) +
    theme_extra(8)
  save_both(p, file.path(out_plots, "cluster_omnibus_existing_runs_fdr_heatmap"), width = 9, height = 16)
}

plot_merged_dose_ploidy <- function(stat_df) {
  plot_df <- stat_df %>%
    dplyr::filter(
      .data$analysis_family %in% c("scVelo_merged_Ploidy", "scVelo_merged_Dose", "Monocle3_global_Dose_universe", "Monocle3_Dose_within_Ploidy_universe"),
      .data$analysis_level == "cell",
      .data$metric %in% c("velocity_pseudotime", "pseudotime"),
      is.finite(.data$primary_p_adj_family)
    ) %>%
    dplyr::mutate(
      analysis_label = coalesce_chr(.data$trajectory_analysis_group, .data$TrajectoryUniverseKey, .data$source_label),
      neg_log10_fdr = -log10(pmax(.data$primary_p_adj_family, .Machine$double.xmin))
    )
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = group_col, y = analysis_label, fill = neg_log10_fdr)) +
    geom_tile(color = "white", linewidth = 0.2, width = 0.95, height = 0.95) +
    facet_grid(analysis_family + metric ~ method, scales = "free_y", space = "free_y") +
    scale_fill_gradient(low = "white", high = "#008837", name = "-log10 FDR") +
    labs(
      title = "Dose and Ploidy tests from merged/universe method runs",
      x = "Tested grouping",
      y = NULL
    ) +
    theme_extra(8)
  save_both(p, file.path(out_plots, "merged_and_universe_dose_ploidy_fdr_heatmap"), width = 10, height = 10)
}

safe_slug <- function(x) {
  x <- as_clean_chr(x)
  x[is.na(x)] <- "NA"
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

coalesce_plot_label <- function(...) {
  out <- coalesce_chr(...)
  out[is.na(out)] <- "Unknown"
  out
}

read_cluster_annotations <- function(root_dir) {
  file <- file.path(root_dir, "00_summary", "cluster_annotation_summary.csv")
  if (!file.exists(file) || file.info(file)$size <= 0) return(data.frame())
  df <- readr::read_csv(
    file,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  )
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df <- ensure_columns(df, c(
    "cluster", "annotation_primary", "annotation_secondary", "annotation_tertiary",
    "annotation_multi", "n_significant_hallmarks"
  ))
  df %>%
    dplyr::mutate(
      cluster = as_clean_chr(.data$cluster),
      n_significant_hallmarks = safe_numeric(.data$n_significant_hallmarks)
    )
}

select_descriptive_cells <- function(df, metric_col, method_label) {
  needed <- c(
    "cell", "source_context", "source_label", "TrajectoryScope", "TrajectoryComparisonID",
    "UMAP_1", "UMAP_2", "clusters", "Dose", "Ploidy", "sample", metric_col
  )
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(.metric = safe_numeric(.data[[metric_col]])) %>%
    dplyr::filter(!is.na(.data$cell), is.finite(.data$.metric))

  candidates <- list(
    velocity_all_cells = df %>%
      dplyr::filter(
        .data$source_context == "velocity_groups",
        coalesce_chr(.data$trajectory_analysis_group, .data$source_label) == "All_cells"
      ),
    velocity_three_group_runs = df %>%
      dplyr::filter(.data$source_context == "velocity_groups"),
    monocle3_all_cells = df %>%
      dplyr::filter(
        .data$source_context == "monocle3_groups",
        coalesce_chr(.data$trajectory_analysis_group, .data$AnalysisGroup, .data$source_label) == "All_cells"
      ),
    monocle3_three_group_runs = df %>%
      dplyr::filter(.data$source_context == "monocle3_groups"),
    DEG_global_ploidy = df %>%
      dplyr::filter(.data$source_context == "DEG_comparison", .data$TrajectoryScope == "global_Ploidy_4N_vs_2N"),
    DEG_global_tn = df %>%
      dplyr::filter(.data$source_context == "DEG_comparison", .data$TrajectoryScope == "global_TN_Tumor_vs_CellLine"),
    DEG_global_dose = df %>%
      dplyr::filter(.data$source_context == "DEG_comparison", .data$TrajectoryScope == "global_Dose_vs_rest"),
    pseudo_dose_karyotype = df %>%
      dplyr::filter(.data$source_context == "by_dose_karyotype"),
    existing_by_dose = df %>%
      dplyr::filter(.data$source_context == "by_dose_or_dose_ploidy")
  )

  chosen_name <- NA_character_
  chosen <- data.frame()
  for (nm in names(candidates)) {
    candidate <- candidates[[nm]]
    if (nrow(candidate) > 0) {
      chosen_name <- nm
      chosen <- candidate
      break
    }
  }
  if (nrow(chosen) == 0) {
    chosen_name <- "all_available_distinct_cells"
    chosen <- df
  }

  chosen %>%
    dplyr::arrange(.data$TrajectoryComparisonID, .data$source_label, .data$cell) %>%
    dplyr::distinct(.data$cell, .keep_all = TRUE) %>%
    dplyr::mutate(
      descriptor_source = chosen_name,
      method = method_label,
      metric = metric_col,
      metric_value = safe_numeric(.data[[metric_col]]),
      Dose_plot = coalesce_plot_label(.data$Dose),
      Ploidy_plot = coalesce_plot_label(.data$Ploidy),
      cluster_plot = coalesce_plot_label(.data$clusters)
    )
}

add_stratum_label <- function(df, cols) {
  if (nrow(df) == 0) return(df)
  label_df <- as.data.frame(lapply(df[, cols, drop = FALSE], as.character), stringsAsFactors = FALSE)
  df$stratum_label <- apply(label_df, 1, paste, collapse = " | ")
  df
}

summarize_method_strata <- function(desc_df, metric_col) {
  strata <- list(
    cluster = c("clusters"),
    ploidy = c("Ploidy"),
    dose = c("Dose"),
    ploidy_dose = c("Ploidy", "Dose"),
    cluster_ploidy = c("clusters", "Ploidy"),
    cluster_dose = c("clusters", "Dose"),
    cluster_ploidy_dose = c("clusters", "Ploidy", "Dose")
  )
  out <- list()
  idx <- 1L
  for (nm in names(strata)) {
    cols <- strata[[nm]]
    tmp <- desc_df %>%
      dplyr::filter(is.finite(.data$metric_value)) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(cols), ~ !is.na(.x))) %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(cols))) %>%
      dplyr::summarise(
        n_cells = dplyr::n(),
        n_samples = dplyr::n_distinct(.data$sample[!is.na(.data$sample)]),
        mean = mean(.data$metric_value, na.rm = TRUE),
        median = stats::median(.data$metric_value, na.rm = TRUE),
        q25 = stats::quantile(.data$metric_value, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(.data$metric_value, 0.75, na.rm = TRUE, names = FALSE),
        .groups = "drop"
      )
    if (nrow(tmp) > 0) {
      tmp <- add_stratum_label(tmp, cols)
      tmp$stratum_type <- nm
      tmp$metric <- metric_col
      out[[idx]] <- tmp
      idx <- idx + 1L
    }
  }
  dplyr::bind_rows(out)
}

dose_to_numeric <- function(x) {
  x <- as_clean_chr(x)
  dplyr::case_when(
    x == "0mg/kg" ~ 0,
    x == "30mg/kg" ~ 30,
    x == "120mg/kg" ~ 120,
    TRUE ~ NA_real_
  )
}

safe_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(unique(x)) < 2 || length(y) < 2) return(NA_real_)
  fit <- tryCatch(stats::lm(y ~ x), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  unname(stats::coef(fit)[["x"]])
}

make_cluster_commonality <- function(desc_df, cluster_annotations = data.frame()) {
  base <- desc_df %>%
    dplyr::filter(!is.na(.data$clusters), is.finite(.data$metric_value))
  if (nrow(base) == 0) return(data.frame())

  cluster_overall <- base %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_time = stats::median(.data$metric_value, na.rm = TRUE),
      mean_time = mean(.data$metric_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      time_rank_early_to_late = rank(.data$median_time, ties.method = "average"),
      relative_timing = dplyr::case_when(
        .data$median_time <= stats::quantile(.data$median_time, 0.25, na.rm = TRUE, names = FALSE) ~ "early_or_low_time",
        .data$median_time >= stats::quantile(.data$median_time, 0.75, na.rm = TRUE, names = FALSE) ~ "late_or_high_time",
        TRUE ~ "middle"
      )
    )

  ploidy_wide <- base %>%
    dplyr::filter(!is.na(.data$Ploidy)) %>%
    dplyr::group_by(.data$clusters, .data$Ploidy) %>%
    dplyr::summarise(median_time = stats::median(.data$metric_value, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = Ploidy, values_from = median_time, names_prefix = "ploidy_")
  ploidy_wide <- ensure_columns(ploidy_wide, c("ploidy_4N", "ploidy_2N"), default = NA_real_)
  ploidy_wide$ploidy_4N <- safe_numeric(ploidy_wide$ploidy_4N)
  ploidy_wide$ploidy_2N <- safe_numeric(ploidy_wide$ploidy_2N)
  ploidy_wide <- ploidy_wide %>%
    dplyr::mutate(
      ploidy_4N_minus_2N = .data$ploidy_4N - .data$ploidy_2N,
      ploidy_pattern = dplyr::case_when(
        is.finite(.data$ploidy_4N_minus_2N) & .data$ploidy_4N_minus_2N > 0 ~ "4N_later_or_higher",
        is.finite(.data$ploidy_4N_minus_2N) & .data$ploidy_4N_minus_2N < 0 ~ "2N_later_or_higher",
        is.finite(.data$ploidy_4N_minus_2N) ~ "similar",
        TRUE ~ "not_enough_groups"
      )
    )

  dose_slope <- base %>%
    dplyr::filter(!is.na(.data$Ploidy), !is.na(.data$Dose)) %>%
    dplyr::mutate(dose_num = dose_to_numeric(.data$Dose)) %>%
    dplyr::filter(is.finite(.data$dose_num)) %>%
    dplyr::group_by(.data$clusters, .data$Ploidy, .data$Dose, .data$dose_num) %>%
    dplyr::summarise(median_time = stats::median(.data$metric_value, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(.data$clusters, .data$Ploidy) %>%
    dplyr::summarise(
      dose_slope_per_mgkg = safe_slope(.data$dose_num, .data$median_time),
      n_dose_levels = dplyr::n_distinct(.data$Dose),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Ploidy,
      values_from = c(dose_slope_per_mgkg, n_dose_levels),
      names_sep = "_"
    )
  dose_slope <- ensure_columns(
    dose_slope,
    c("dose_slope_per_mgkg_2N", "dose_slope_per_mgkg_4N", "n_dose_levels_2N", "n_dose_levels_4N"),
    default = NA_real_
  )
  for (col in c("dose_slope_per_mgkg_2N", "dose_slope_per_mgkg_4N", "n_dose_levels_2N", "n_dose_levels_4N")) {
    dose_slope[[col]] <- safe_numeric(dose_slope[[col]])
  }
  dose_slope <- dose_slope %>%
    dplyr::mutate(
      dose_pattern_2N = dplyr::case_when(
        is.finite(.data$dose_slope_per_mgkg_2N) & .data$dose_slope_per_mgkg_2N > 0 ~ "higher_dose_later_or_higher",
        is.finite(.data$dose_slope_per_mgkg_2N) & .data$dose_slope_per_mgkg_2N < 0 ~ "higher_dose_earlier_or_lower",
        is.finite(.data$dose_slope_per_mgkg_2N) ~ "flat",
        TRUE ~ "not_enough_doses"
      ),
      dose_pattern_4N = dplyr::case_when(
        is.finite(.data$dose_slope_per_mgkg_4N) & .data$dose_slope_per_mgkg_4N > 0 ~ "higher_dose_later_or_higher",
        is.finite(.data$dose_slope_per_mgkg_4N) & .data$dose_slope_per_mgkg_4N < 0 ~ "higher_dose_earlier_or_lower",
        is.finite(.data$dose_slope_per_mgkg_4N) ~ "flat",
        TRUE ~ "not_enough_doses"
      ),
      dose_commonality = dplyr::case_when(
        .data$dose_pattern_2N == .data$dose_pattern_4N &
          !.data$dose_pattern_2N %in% c("not_enough_doses") ~ .data$dose_pattern_2N,
        TRUE ~ "ploidy_specific_or_incomplete"
      )
    )

  out <- cluster_overall %>%
    dplyr::left_join(ploidy_wide, by = "clusters") %>%
    dplyr::left_join(dose_slope, by = "clusters") %>%
    dplyr::arrange(.data$time_rank_early_to_late)

  if (nrow(cluster_annotations) > 0) {
    out <- out %>%
      dplyr::left_join(cluster_annotations, by = c("clusters" = "cluster"))
  }
  out
}

format_p_value_label <- function(p) {
  p <- safe_numeric(p)[1]
  if (!is.finite(p)) return("p=NA")
  if (p < 1e-4) return("p<1e-4")
  paste0("p=", signif(p, 2))
}

p_value_stars <- function(p) {
  p <- safe_numeric(p)[1]
  dplyr::case_when(
    !is.finite(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

run_violin_panel_test <- function(df, group_col, value_col, min_group_n = 3L) {
  test_df <- df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .value = safe_numeric(.data[[value_col]])
    ) %>%
    dplyr::filter(!is.na(.data$.group), is.finite(.data$.value))
  if (nrow(test_df) == 0) {
    return(data.frame(test = "not_tested", p = NA_real_, n_groups = 0L, label = "not tested", stringsAsFactors = FALSE))
  }

  group_counts <- table(test_df$.group)
  keep_groups <- names(group_counts)[group_counts >= min_group_n]
  test_df <- test_df[test_df$.group %in% keep_groups, , drop = FALSE]
  groups <- sort_maybe_numeric(unique(test_df$.group))
  if (length(groups) < 2L) {
    return(data.frame(test = "not_tested", p = NA_real_, n_groups = length(groups), label = "not tested", stringsAsFactors = FALSE))
  }

  if (length(groups) == 2L) {
    x <- test_df$.value[test_df$.group == groups[1]]
    y <- test_df$.value[test_df$.group == groups[2]]
    wt <- tryCatch(
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)),
      error = function(e) NULL
    )
    p <- if (!is.null(wt)) wt$p.value else NA_real_
    test_name <- "Wilcoxon"
  } else {
    kt <- tryCatch(
      stats::kruskal.test(.value ~ .group, data = test_df),
      error = function(e) NULL
    )
    p <- if (!is.null(kt)) kt$p.value else NA_real_
    test_name <- "Kruskal-Wallis"
  }

  data.frame(
    test = test_name,
    p = p,
    n_groups = length(groups),
    label = paste0(test_name, "\n", format_p_value_label(p), " ", p_value_stars(p)),
    stringsAsFactors = FALSE
  )
}

annotation_y_position <- function(x, pad_mult = 0.08) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) span <- max(abs(rng), na.rm = TRUE) * 0.05
  if (!is.finite(span) || span <= 0) span <- 0.05
  rng[2] + span * pad_mult
}

make_violin_stat_annotations <- function(
  df,
  split_cols,
  group_col,
  value_col = "metric_value",
  x_col = NULL,
  x_value = NULL,
  min_group_n = 3L,
  label_prefix = NULL
) {
  needed <- unique(c(split_cols, group_col, value_col, x_col))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::filter(is.finite(safe_numeric(.data[[value_col]])))
  if (nrow(df) == 0) return(data.frame())

  if (length(split_cols) == 0) {
    split_keys <- "all"
    df$.split_key <- "all"
  } else {
    split_label_df <- as.data.frame(lapply(df[, split_cols, drop = FALSE], as.character), stringsAsFactors = FALSE)
    df$.split_key <- apply(split_label_df, 1, paste, collapse = "\r")
    split_keys <- unique(df$.split_key)
  }

  rows <- vector("list", length(split_keys))
  for (i in seq_along(split_keys)) {
    sub_df <- df[df$.split_key == split_keys[i], , drop = FALSE]
    if (nrow(sub_df) == 0) next
    test_row <- run_violin_panel_test(sub_df, group_col = group_col, value_col = value_col, min_group_n = min_group_n)
    if (identical(test_row$test[1], "not_tested")) next
    split_values <- if (length(split_cols) == 0) data.frame(.dummy_split = "all") else sub_df[1, split_cols, drop = FALSE]
    label <- test_row$label[1]
    if (!is.null(label_prefix) && nzchar(label_prefix)) {
      label <- paste0(label_prefix, "\n", label)
    }
    x_plot <- if (!is.null(x_col) && x_col %in% names(sub_df)) {
      as.character(sub_df[[x_col]][1])
    } else if (!is.null(x_value)) {
      as.character(x_value)
    } else {
      group_levels <- sort_maybe_numeric(unique(as_clean_chr(sub_df[[group_col]])))
      group_levels[ceiling(length(group_levels) / 2)]
    }
    rows[[i]] <- cbind(
      split_values,
      data.frame(
        x_plot = x_plot,
        y_plot = annotation_y_position(sub_df[[value_col]]),
        stat_test = test_row$test[1],
        stat_p = test_row$p[1],
        stat_n_groups = test_row$n_groups[1],
        label = label,
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) {
    split_empty <- as.data.frame(
      stats::setNames(rep(list(character()), length(split_cols)), split_cols),
      stringsAsFactors = FALSE
    )
    return(cbind(
      split_empty,
      data.frame(
        x_plot = character(),
        y_plot = numeric(),
        stat_test = character(),
        stat_p = numeric(),
        stat_n_groups = integer(),
        label = character(),
        stringsAsFactors = FALSE
      ),
      stringsAsFactors = FALSE
    ))
  }
  if (".dummy_split" %in% names(out)) out$.dummy_split <- NULL
  out
}

plot_method_distributions <- function(desc_df, method_label, metric_col, out_dir, plot_prefix) {
  unlink(paste0(file.path(out_dir, paste0(plot_prefix, "_cluster_ploidy_dose_median_time_heatmap")), c(".pdf", ".png")))
  unlink(paste0(file.path(out_dir, paste0(plot_prefix, "_cluster_median_time_order")), c(".pdf", ".png")))

  plot_df <- desc_df %>%
    dplyr::filter(
      !is.na(.data$clusters),
      !is.na(.data$Dose),
      !is.na(.data$Ploidy),
      is.finite(.data$metric_value)
    ) %>%
    dplyr::mutate(
      clusters = factor(as.character(.data$clusters), levels = sort_maybe_numeric(as.character(.data$clusters))),
      Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N"))
    ) %>%
    dplyr::add_count(.data$clusters, .data$Ploidy, .data$Dose, name = ".stratum_n") %>%
    dplyr::filter(.data$.stratum_n >= 2)
  if (nrow(plot_df) > 0) {
    set.seed(1234)
    point_df <- plot_df %>%
      dplyr::group_by(.data$clusters, .data$Ploidy, .data$Dose) %>%
      dplyr::mutate(.rand = stats::runif(dplyr::n())) %>%
      dplyr::arrange(.data$.rand, .by_group = TRUE) %>%
      dplyr::slice_head(n = 180) %>%
      dplyr::ungroup()
    stat_labels <- make_violin_stat_annotations(
      plot_df,
      split_cols = c("clusters", "Ploidy"),
      group_col = "Dose",
      value_col = "metric_value",
      x_value = "30mg/kg",
      min_group_n = 3L,
      label_prefix = "Dose"
    )

    p <- ggplot(plot_df, aes(x = Dose, y = metric_value, fill = Dose)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.75, linewidth = 0.15, color = "grey30") +
      geom_point(
        data = point_df,
        aes(x = Dose, y = metric_value),
        inherit.aes = FALSE,
        position = position_jitter(width = 0.13, height = 0),
        size = 0.18,
        alpha = 0.18,
        color = "black"
      ) +
      facet_grid(clusters ~ Ploidy) +
      scale_fill_manual(values = c("0mg/kg" = "#4B8BBE", "30mg/kg" = "#F0A202", "120mg/kg" = "#C44536"), drop = FALSE) +
      labs(
        title = paste0(method_label, " cluster pseudotime distributions by Ploidy and Dose"),
        subtitle = "Violin shapes show per-cell distributions; panel labels report Kruskal-Wallis tests across Dose.",
        x = "Dose",
        y = metric_col,
        fill = "Dose"
      ) +
      geom_text(
        data = stat_labels,
        aes(x = x_plot, y = y_plot, label = label),
        inherit.aes = FALSE,
        size = 2.25,
        lineheight = 0.88,
        fontface = "bold",
        color = "black"
      ) +
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.22))) +
      theme_extra(8) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_cluster_ploidy_dose_time_distribution")), width = 10, height = 13)
  }

  cluster_df <- desc_df %>%
    dplyr::filter(!is.na(.data$clusters), is.finite(.data$metric_value)) %>%
    dplyr::mutate(
      clusters = factor(as.character(.data$clusters), levels = sort_maybe_numeric(as.character(.data$clusters)))
    ) %>%
    dplyr::add_count(.data$clusters, name = ".cluster_n") %>%
    dplyr::filter(.data$.cluster_n >= 2)
  if (nrow(cluster_df) > 0) {
    set.seed(1234)
    point_df <- cluster_df %>%
      dplyr::group_by(.data$clusters) %>%
      dplyr::mutate(.rand = stats::runif(dplyr::n())) %>%
      dplyr::arrange(.data$.rand, .by_group = TRUE) %>%
      dplyr::slice_head(n = 250) %>%
      dplyr::ungroup()
    cluster_levels_for_plot <- levels(cluster_df$clusters)
    stat_labels <- make_violin_stat_annotations(
      cluster_df,
      split_cols = character(0),
      group_col = "clusters",
      value_col = "metric_value",
      x_value = cluster_levels_for_plot[ceiling(length(cluster_levels_for_plot) / 2)],
      min_group_n = 3L,
      label_prefix = "Cluster"
    )
    p <- ggplot(cluster_df, aes(x = clusters, y = metric_value, fill = clusters)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.8, linewidth = 0.15, color = "grey30") +
      geom_point(
        data = point_df,
        aes(x = clusters, y = metric_value),
        inherit.aes = FALSE,
        position = position_jitter(width = 0.12, height = 0),
        size = 0.18,
        alpha = 0.18,
        color = "black"
      ) +
      labs(
        title = paste0(method_label, " cluster-level pseudotime distributions"),
        subtitle = "Label reports the Kruskal-Wallis test across clusters.",
        x = "Cluster",
        y = metric_col
      ) +
      geom_label(
        data = stat_labels,
        aes(x = x_plot, y = y_plot, label = label),
        inherit.aes = FALSE,
        size = 3,
        lineheight = 0.9,
        fontface = "bold",
        fill = "white",
        color = "black",
        alpha = 0.86,
        linewidth = 0.18
      ) +
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.18))) +
      guides(fill = "none") +
      theme_extra(10)
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_cluster_time_distribution_by_cluster")), width = 9.5, height = 5.5)
  }
}

plot_umap_overlays <- function(desc_df, method_label, metric_col, out_dir, plot_prefix, out_umap_dir) {
  plot_df <- desc_df %>%
    dplyr::filter(is.finite(.data$UMAP_1), is.finite(.data$UMAP_2), is.finite(.data$metric_value))
  if (nrow(plot_df) == 0) return(invisible(NULL))
  if (nrow(plot_df) > 60000) {
    set.seed(1234)
    plot_df <- plot_df[sort(sample(seq_len(nrow(plot_df)), 60000)), , drop = FALSE]
  }

  p_time <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = metric_value)) +
    geom_point(size = 0.22, alpha = 0.75) +
    facet_grid(Ploidy_plot ~ Dose_plot) +
	    scale_color_gradientn(
	      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
	      name = metric_col
	    ) +
	    coord_square_umap(plot_df$UMAP_1, plot_df$UMAP_2) +
	    labs(
      title = paste0(method_label, " UMAP overlay by Ploidy and Dose"),
      subtitle = paste0("Cells colored by ", metric_col),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()
	  save_square_panel_both(p_time, file.path(out_dir, paste0(plot_prefix, "_umap_", metric_col, "_by_ploidy_dose")), panel_size = 3.2)
	  save_square_panel_both(p_time, file.path(out_umap_dir, paste0(plot_prefix, "_umap_", metric_col, "_by_ploidy_dose")), panel_size = 3.2)

  cluster_median <- desc_df %>%
    dplyr::filter(!is.na(.data$clusters), is.finite(.data$metric_value)) %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(cluster_median_time = stats::median(.data$metric_value, na.rm = TRUE), .groups = "drop")
  cluster_plot <- plot_df %>%
    dplyr::left_join(cluster_median, by = "clusters")
  if (nrow(cluster_plot) > 0 && any(is.finite(cluster_plot$cluster_median_time))) {
    centers <- cluster_plot %>%
      dplyr::filter(!is.na(.data$clusters)) %>%
      dplyr::group_by(.data$clusters) %>%
      dplyr::summarise(
        UMAP_1 = stats::median(.data$UMAP_1, na.rm = TRUE),
        UMAP_2 = stats::median(.data$UMAP_2, na.rm = TRUE),
        .groups = "drop"
      )
    p_cluster <- ggplot(cluster_plot, aes(x = UMAP_1, y = UMAP_2, color = cluster_median_time)) +
      geom_point(size = 0.22, alpha = 0.65) +
      geom_text(
        data = centers,
        aes(x = UMAP_1, y = UMAP_2, label = clusters),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold",
        color = "black"
      ) +
	      scale_color_gradientn(
	        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
	        name = "Cluster median time"
	      ) +
	      coord_square_umap(cluster_plot$UMAP_1, cluster_plot$UMAP_2) +
	      labs(
        title = paste0(method_label, " UMAP overlay with cluster-level pseudotime order"),
        x = "UMAP 1",
        y = "UMAP 2"
      ) +
	      theme_extra(9) +
	      theme(axis.text.x = element_text(angle = 0)) +
	      theme_square_umap_panel()
	    save_square_panel_both(p_cluster, file.path(out_dir, paste0(plot_prefix, "_umap_cluster_median_time")), panel_size = 5.2)
	    save_square_panel_both(p_cluster, file.path(out_umap_dir, paste0(plot_prefix, "_umap_cluster_median_time")), panel_size = 5.2)
  }
}

write_method_separate_outputs <- function(df, metric_col, method_label, out_dir, plot_prefix, out_umap_dir, cluster_annotations = data.frame()) {
  desc <- select_descriptive_cells(df, metric_col, method_label)
  manifest <- data.frame(
    method = method_label,
    metric = metric_col,
    descriptor_source = unique(desc$descriptor_source)[1],
    n_cells = nrow(desc),
    n_cells_with_umap = sum(is.finite(desc$UMAP_1) & is.finite(desc$UMAP_2)),
    n_clusters = dplyr::n_distinct(desc$clusters[!is.na(desc$clusters)]),
    stringsAsFactors = FALSE
  )
  write_table_csv(manifest, file.path(out_dir, paste0(plot_prefix, "_descriptive_dataset_manifest.csv")))

  strata_summary <- summarize_method_strata(desc, metric_col)
  commonality <- make_cluster_commonality(desc, cluster_annotations)
  write_table_csv(strata_summary, file.path(out_dir, paste0(plot_prefix, "_stratified_time_summary.csv")))
  write_table_csv(commonality, file.path(out_dir, paste0(plot_prefix, "_cluster_ploidy_dose_commonality.csv")))

  plot_method_distributions(desc, method_label, metric_col, out_dir, plot_prefix)
  plot_umap_overlays(desc, method_label, metric_col, out_dir, plot_prefix, out_umap_dir)

  list(desc = desc, strata_summary = strata_summary, commonality = commonality, manifest = manifest)
}

make_tn_ploidy_base <- function(desc_df, metric_col) {
  needed <- c("TN", "Ploidy", "Dose", "clusters", "sample", metric_col, "metric_value")
  desc_df <- ensure_columns(desc_df, needed)
  desc_df %>%
    dplyr::mutate(
      TN = as_clean_chr(.data$TN),
      Ploidy = as_clean_chr(.data$Ploidy),
      Dose = as_clean_chr(.data$Dose),
      clusters = as_clean_chr(.data$clusters),
      sample = as_clean_chr(.data$sample),
      metric_value = safe_numeric(.data$metric_value),
      TN_Ploidy = dplyr::case_when(
        .data$TN == "Tumor" & .data$Ploidy == "2N" ~ "Tumor_2N",
        .data$TN == "Tumor" & .data$Ploidy == "4N" ~ "Tumor_4N",
        .data$TN == "CellLine" & .data$Ploidy == "2N" ~ "CellLine_2N",
        .data$TN == "CellLine" & .data$Ploidy == "4N" ~ "CellLine_4N",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(
      .data$TN %in% c("Tumor", "CellLine"),
      .data$Ploidy %in% c("2N", "4N"),
      !is.na(.data$clusters),
      is.finite(.data$metric_value)
    )
}

summarize_tn_ploidy_clusters <- function(base_df, metric_col) {
  if (nrow(base_df) == 0) return(data.frame())
  base_df %>%
    dplyr::group_by(.data$TN, .data$Ploidy, .data$TN_Ploidy, .data$Dose, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(.data$sample[!is.na(.data$sample)]),
      mean_time = mean(.data$metric_value, na.rm = TRUE),
      median_time = stats::median(.data$metric_value, na.rm = TRUE),
      q25_time = stats::quantile(.data$metric_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_time = stats::quantile(.data$metric_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(metric = metric_col)
}

comparison_direction <- function(delta, tolerance = 0.02) {
  dplyr::case_when(
    !is.finite(delta) ~ "not_enough_groups",
    abs(delta) < tolerance ~ "similar_within_0.02",
    delta > 0 ~ "4N_later_or_higher_time",
    delta < 0 ~ "2N_later_or_higher_time",
    TRUE ~ "not_enough_groups"
  )
}

build_ploidy_comparison_from_summary <- function(summary_df, split_cols, comparison_scope) {
  if (nrow(summary_df) == 0) return(data.frame())
  needed <- unique(c(split_cols, "Ploidy", "n_cells", "n_samples", "median_time", "mean_time", "q25_time", "q75_time"))
  summary_df <- ensure_columns(summary_df, needed)
  wide <- summary_df %>%
    dplyr::filter(.data$Ploidy %in% c("2N", "4N")) %>%
    dplyr::select(dplyr::all_of(needed)) %>%
    tidyr::pivot_wider(
      names_from = Ploidy,
      values_from = c(n_cells, n_samples, median_time, mean_time, q25_time, q75_time),
      names_sep = "_"
    )
  wide <- ensure_columns(
    wide,
    c(
      "n_cells_2N", "n_cells_4N", "n_samples_2N", "n_samples_4N",
      "median_time_2N", "median_time_4N", "mean_time_2N", "mean_time_4N",
      "q25_time_2N", "q25_time_4N", "q75_time_2N", "q75_time_4N"
    ),
    default = NA_real_
  )
  for (col in setdiff(names(wide), split_cols)) {
    if (grepl("^(n_cells|n_samples|median_time|mean_time|q25_time|q75_time)_", col)) {
      wide[[col]] <- safe_numeric(wide[[col]])
    }
  }
  wide %>%
    dplyr::mutate(
      comparison_scope = comparison_scope,
      median_4N_minus_2N = .data$median_time_4N - .data$median_time_2N,
      mean_4N_minus_2N = .data$mean_time_4N - .data$mean_time_2N,
      abs_median_4N_minus_2N = abs(.data$median_4N_minus_2N),
      ploidy_timing_direction = comparison_direction(.data$median_4N_minus_2N)
    )
}

add_ploidy_test_results <- function(comparison_df, test_df, by_cols) {
  if (nrow(comparison_df) == 0 || nrow(test_df) == 0) return(comparison_df)
  test_keep <- test_df %>%
    dplyr::select(
      dplyr::all_of(by_cols),
      ploidy_test = primary_test,
      ploidy_p = primary_p,
      ploidy_fdr_global = primary_p_adj_global,
      ploidy_fdr_family = primary_p_adj_family,
      ploidy_significant_family_fdr_0_05 = significant_family_fdr_0_05,
      rank_biserial
    )
  dplyr::left_join(comparison_df, test_keep, by = by_cols)
}

make_tn_ploidy_cluster_comparisons <- function(base_df, summary_df, metric_col, method_label) {
  out <- list()
  idx <- 1L

  tumor_summary_by_dose <- summary_df %>%
    dplyr::filter(.data$TN == "Tumor", !is.na(.data$Dose))
  tumor_by_dose <- build_ploidy_comparison_from_summary(
    tumor_summary_by_dose,
    split_cols = c("TN", "Dose", "clusters"),
    comparison_scope = "Tumor_2N_vs_4N_by_Dose"
  )
  tumor_test_by_dose <- run_tests_by_split(
    base_df %>% dplyr::filter(.data$TN == "Tumor", !is.na(.data$Dose)),
    split_cols = c("TN", "Dose", "clusters"),
    group_col = "Ploidy",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tn_ploidy_cluster_Tumor_2N_vs_4N_by_Dose",
    level = "cell"
  ) %>%
    add_adjusted_p()
  tumor_by_dose <- add_ploidy_test_results(tumor_by_dose, tumor_test_by_dose, c("TN", "Dose", "clusters"))
  out[[idx]] <- tumor_by_dose
  idx <- idx + 1L

  tumor_summary_pooled <- base_df %>%
    dplyr::filter(.data$TN == "Tumor") %>%
    dplyr::group_by(.data$TN, .data$Ploidy, .data$TN_Ploidy, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(.data$sample[!is.na(.data$sample)]),
      mean_time = mean(.data$metric_value, na.rm = TRUE),
      median_time = stats::median(.data$metric_value, na.rm = TRUE),
      q25_time = stats::quantile(.data$metric_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_time = stats::quantile(.data$metric_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
  tumor_pooled <- build_ploidy_comparison_from_summary(
    tumor_summary_pooled,
    split_cols = c("TN", "clusters"),
    comparison_scope = "Tumor_2N_vs_4N_pooled_Dose"
  )
  tumor_test_pooled <- run_tests_by_split(
    base_df %>% dplyr::filter(.data$TN == "Tumor"),
    split_cols = c("TN", "clusters"),
    group_col = "Ploidy",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tn_ploidy_cluster_Tumor_2N_vs_4N_pooled_Dose",
    level = "cell"
  ) %>%
    add_adjusted_p()
  tumor_pooled <- add_ploidy_test_results(tumor_pooled, tumor_test_pooled, c("TN", "clusters"))
  out[[idx]] <- tumor_pooled
  idx <- idx + 1L

  cellline_summary <- base_df %>%
    dplyr::filter(.data$TN == "CellLine") %>%
    dplyr::group_by(.data$TN, .data$Ploidy, .data$TN_Ploidy, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(.data$sample[!is.na(.data$sample)]),
      mean_time = mean(.data$metric_value, na.rm = TRUE),
      median_time = stats::median(.data$metric_value, na.rm = TRUE),
      q25_time = stats::quantile(.data$metric_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_time = stats::quantile(.data$metric_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )
  cellline_comp <- build_ploidy_comparison_from_summary(
    cellline_summary,
    split_cols = c("TN", "clusters"),
    comparison_scope = "CellLine_2N_vs_4N"
  )
  cellline_test <- run_tests_by_split(
    base_df %>% dplyr::filter(.data$TN == "CellLine"),
    split_cols = c("TN", "clusters"),
    group_col = "Ploidy",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tn_ploidy_cluster_CellLine_2N_vs_4N",
    level = "cell"
  ) %>%
    add_adjusted_p()
  cellline_comp <- add_ploidy_test_results(cellline_comp, cellline_test, c("TN", "clusters"))
  out[[idx]] <- cellline_comp

  dplyr::bind_rows(out) %>%
    dplyr::mutate(method = method_label, metric = metric_col) %>%
    dplyr::arrange(.data$comparison_scope, .data$TN, .data$Dose, factor(as.character(.data$clusters), levels = sort_maybe_numeric(as.character(.data$clusters))))
}

make_tumor_dose_response_deltas <- function(summary_df, metric_col, method_label) {
  tumor <- summary_df %>%
    dplyr::filter(.data$TN == "Tumor", .data$Ploidy %in% c("2N", "4N"), !is.na(.data$Dose))
  if (nrow(tumor) == 0) return(list(dose_delta = data.frame(), interaction = data.frame()))

  baseline <- tumor %>%
    dplyr::filter(.data$Dose == "0mg/kg") %>%
    dplyr::select(TN, Ploidy, clusters, baseline_median_time = median_time, baseline_n_cells = n_cells)
  dose_delta <- tumor %>%
    dplyr::filter(.data$Dose %in% c("30mg/kg", "120mg/kg")) %>%
    dplyr::left_join(baseline, by = c("TN", "Ploidy", "clusters")) %>%
    dplyr::mutate(
      dose_vs_baseline = paste0(.data$Dose, "_vs_0mg/kg"),
      median_delta_vs_0 = .data$median_time - .data$baseline_median_time,
      dose_response_direction = dplyr::case_when(
        !is.finite(.data$median_delta_vs_0) ~ "no_0mg_baseline",
        abs(.data$median_delta_vs_0) < 0.02 ~ "flat_within_0.02",
        .data$median_delta_vs_0 > 0 ~ "treatment_later_or_higher_time",
        .data$median_delta_vs_0 < 0 ~ "treatment_earlier_or_lower_time",
        TRUE ~ "no_0mg_baseline"
      ),
      method = method_label,
      metric = metric_col
    )

  interaction <- dose_delta %>%
    dplyr::select(clusters, Dose, dose_vs_baseline, Ploidy, median_delta_vs_0) %>%
    tidyr::pivot_wider(names_from = Ploidy, values_from = median_delta_vs_0, names_prefix = "delta_")
  interaction <- ensure_columns(interaction, c("delta_2N", "delta_4N"), default = NA_real_)
  interaction$delta_2N <- safe_numeric(interaction$delta_2N)
  interaction$delta_4N <- safe_numeric(interaction$delta_4N)
  interaction <- interaction %>%
    dplyr::mutate(
      interaction_delta_4N_minus_2N = .data$delta_4N - .data$delta_2N,
      abs_interaction_delta = abs(.data$interaction_delta_4N_minus_2N),
      interaction_pattern = dplyr::case_when(
        !is.finite(.data$interaction_delta_4N_minus_2N) ~ "not_enough_groups",
        abs(.data$interaction_delta_4N_minus_2N) < 0.02 ~ "similar_dose_response_within_0.02",
        .data$interaction_delta_4N_minus_2N > 0 ~ "2N_shifts_earlier_more_than_4N",
        .data$interaction_delta_4N_minus_2N < 0 ~ "4N_shifts_earlier_more_than_2N",
        TRUE ~ "not_enough_groups"
      ),
      method = method_label,
      metric = metric_col
    )

  list(dose_delta = dose_delta, interaction = interaction)
}

make_pairwise_tests_by_split <- function(df, split_cols, group_col, metric_col, method, analysis_family, level, pairs) {
  needed <- unique(c(split_cols, group_col, metric_col))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(
      .group = as_clean_chr(.data[[group_col]]),
      .metric = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(!is.na(.data$.group), is.finite(.data$.metric))
  if (nrow(df) == 0) return(data.frame())

  if (length(split_cols) == 0) {
    df$.split_key <- "all"
    split_keys <- "all"
  } else {
    split_df <- as.data.frame(lapply(df[, split_cols, drop = FALSE], as.character), stringsAsFactors = FALSE)
    df$.split_key <- apply(split_df, 1, paste, collapse = "\r")
    split_keys <- unique(df$.split_key)
  }

  rows <- list()
  idx <- 1L
  min_group_n <- if (identical(level, "sample")) 2L else 3L
  for (split_key in split_keys) {
    sub_df <- df[df$.split_key == split_key, , drop = FALSE]
    split_values <- if (length(split_cols) == 0) {
      data.frame(.dummy_split = "all", stringsAsFactors = FALSE)
    } else {
      sub_df[1, split_cols, drop = FALSE]
    }
    for (pair in pairs) {
      pair <- as_clean_chr(pair)
      pair <- pair[!is.na(pair)]
      if (length(pair) != 2L) next
      pair_df <- sub_df[sub_df$.group %in% pair, , drop = FALSE]
      test_df <- run_one_test(
        pair_df,
        group_col = group_col,
        metric_col = metric_col,
        min_group_n = min_group_n,
        preferred_groups = pair
      )
      test_df <- test_df %>%
        dplyr::mutate(
          contrast = paste0(pair[2], "_vs_", pair[1]),
          contrast_group_1 = pair[1],
          contrast_group_2 = pair[2],
          median_diff_group_2_minus_group_1 = .data$median_group_2 - .data$median_group_1
        )
      rows[[idx]] <- cbind(split_values, test_df, stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }

  out <- dplyr::bind_rows(rows)
  if (nrow(out) == 0) return(data.frame())
  if (".dummy_split" %in% names(out)) out$.dummy_split <- NULL
  out %>%
    dplyr::mutate(
      method = method,
      metric = metric_col,
      analysis_family = analysis_family,
      analysis_level = level,
      group_col = group_col
    )
}

compare_lm_models_p <- function(reduced_formula, full_formula, df) {
  reduced <- tryCatch(stats::lm(reduced_formula, data = df), error = function(e) NULL)
  full <- tryCatch(stats::lm(full_formula, data = df), error = function(e) NULL)
  if (is.null(reduced) || is.null(full)) return(NA_real_)
  an <- tryCatch(stats::anova(reduced, full), error = function(e) NULL)
  if (is.null(an) || nrow(an) < 2 || !("Pr(>F)" %in% names(an))) return(NA_real_)
  p <- as.numeric(an[2, "Pr(>F)"])
  if (length(p) == 0) NA_real_ else p
}

make_tumor_ploidy_dose_feature_summary <- function(tumor_df, metric_col, method_label, effect_tolerance = 0.02) {
  needed <- c("TN", "Ploidy", "Dose", "clusters", "sample", metric_col)
  tumor_df <- ensure_columns(tumor_df, needed)
  sample_df <- tumor_df %>%
    dplyr::mutate(
      TN = as_clean_chr(.data$TN),
      Ploidy = as_clean_chr(.data$Ploidy),
      Dose = as_clean_chr(.data$Dose),
      clusters = as_clean_chr(.data$clusters),
      sample = as_clean_chr(.data$sample),
      metric_value = safe_numeric(.data[[metric_col]])
    ) %>%
    dplyr::filter(
      .data$TN == "Tumor",
      .data$Ploidy %in% c("2N", "4N"),
      .data$Dose %in% c("0mg/kg", "30mg/kg", "120mg/kg"),
      !is.na(.data$clusters),
      !is.na(.data$sample),
      is.finite(.data$metric_value)
    ) %>%
    dplyr::group_by(.data$clusters, .data$sample, .data$Ploidy, .data$Dose) %>%
    dplyr::summarise(
      sample_median_time = stats::median(.data$metric_value, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(.data$sample_median_time))
  if (nrow(sample_df) == 0) return(data.frame())

  cluster_levels <- sort_maybe_numeric(unique(sample_df$clusters))
  rows <- vector("list", length(cluster_levels))
  for (i in seq_along(cluster_levels)) {
    cl <- cluster_levels[i]
    cl_df <- sample_df %>%
      dplyr::filter(.data$clusters == cl) %>%
      dplyr::mutate(
        Ploidy = factor(.data$Ploidy, levels = c("2N", "4N")),
        Dose = factor(.data$Dose, levels = c("0mg/kg", "30mg/kg", "120mg/kg"))
      )
    n_ploidy <- dplyr::n_distinct(as.character(cl_df$Ploidy[!is.na(cl_df$Ploidy)]))
    n_dose <- dplyr::n_distinct(as.character(cl_df$Dose[!is.na(cl_df$Dose)]))
    n_sample_strata <- nrow(cl_df)
    n_samples <- dplyr::n_distinct(cl_df$sample)

    ploidy_medians <- tapply(cl_df$sample_median_time, cl_df$Ploidy, stats::median, na.rm = TRUE)
    ploidy_effect <- if (all(c("2N", "4N") %in% names(ploidy_medians))) {
      unname(ploidy_medians[["4N"]] - ploidy_medians[["2N"]])
    } else {
      NA_real_
    }
    dose_medians <- tapply(cl_df$sample_median_time, cl_df$Dose, stats::median, na.rm = TRUE)
    dose_effect_range <- if (sum(is.finite(dose_medians)) >= 2L) {
      unname(max(dose_medians, na.rm = TRUE) - min(dose_medians, na.rm = TRUE))
    } else {
      NA_real_
    }
    finite_dose_medians <- dose_medians[is.finite(dose_medians)]
    dose_top_group <- if (length(finite_dose_medians) >= 1L) names(finite_dose_medians)[which.max(finite_dose_medians)] else NA_character_
    dose_bottom_group <- if (length(finite_dose_medians) >= 1L) names(finite_dose_medians)[which.min(finite_dose_medians)] else NA_character_

    diff_by_dose <- cl_df %>%
      dplyr::group_by(.data$Dose, .data$Ploidy) %>%
      dplyr::summarise(median_time = stats::median(.data$sample_median_time, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Ploidy, values_from = median_time)
    diff_by_dose <- ensure_columns(diff_by_dose, c("2N", "4N"), default = NA_real_)
    diff_by_dose$ploidy_diff_4N_minus_2N <- safe_numeric(diff_by_dose[["4N"]]) - safe_numeric(diff_by_dose[["2N"]])
    finite_diffs <- diff_by_dose$ploidy_diff_4N_minus_2N[is.finite(diff_by_dose$ploidy_diff_4N_minus_2N)]
    interaction_effect_max_delta <- if (length(finite_diffs) >= 2L) {
      max(finite_diffs, na.rm = TRUE) - min(finite_diffs, na.rm = TRUE)
    } else {
      NA_real_
    }
    get_dose_diff <- function(dose_label) {
      vals <- diff_by_dose$ploidy_diff_4N_minus_2N[as.character(diff_by_dose$Dose) == dose_label]
      if (length(vals) == 0) NA_real_ else vals[1]
    }

    p_ploidy <- p_dose <- p_interaction <- NA_real_
    model_status <- "not_fit"
    if (n_sample_strata >= 6L && n_ploidy >= 2L && n_dose >= 2L) {
      p_ploidy <- compare_lm_models_p(sample_median_time ~ Dose, sample_median_time ~ Ploidy + Dose, cl_df)
      p_dose <- compare_lm_models_p(sample_median_time ~ Ploidy, sample_median_time ~ Ploidy + Dose, cl_df)
      p_interaction <- compare_lm_models_p(sample_median_time ~ Ploidy + Dose, sample_median_time ~ Ploidy * Dose, cl_df)
      model_status <- if (any(is.finite(c(p_ploidy, p_dose, p_interaction)))) "ok" else "singular_or_insufficient_df"
    } else {
      model_status <- "insufficient_sample_strata"
    }

    rows[[i]] <- data.frame(
      method = method_label,
      metric = metric_col,
      TN = "Tumor",
      clusters = cl,
      n_sample_strata = n_sample_strata,
      n_samples = n_samples,
      n_ploidy_levels = n_ploidy,
      n_dose_levels = n_dose,
      model_status = model_status,
      p_ploidy_additive = p_ploidy,
      p_dose_additive = p_dose,
      p_ploidy_dose_interaction = p_interaction,
      ploidy_effect_4N_minus_2N = ploidy_effect,
      abs_ploidy_effect_4N_minus_2N = abs(ploidy_effect),
      dose_effect_range = dose_effect_range,
      abs_dose_effect_range = abs(dose_effect_range),
      dose_top_median_group = dose_top_group,
      dose_bottom_median_group = dose_bottom_group,
      ploidy_diff_4N_minus_2N_0mg_kg = get_dose_diff("0mg/kg"),
      ploidy_diff_4N_minus_2N_30mg_kg = get_dose_diff("30mg/kg"),
      ploidy_diff_4N_minus_2N_120mg_kg = get_dose_diff("120mg/kg"),
      interaction_effect_max_delta = interaction_effect_max_delta,
      abs_interaction_effect_max_delta = abs(interaction_effect_max_delta),
      stringsAsFactors = FALSE
    )
  }

  dplyr::bind_rows(rows) %>%
    dplyr::mutate(
      ploidy_fdr = stats::p.adjust(.data$p_ploidy_additive, method = "BH"),
      dose_fdr = stats::p.adjust(.data$p_dose_additive, method = "BH"),
      ploidy_dose_interaction_fdr = stats::p.adjust(.data$p_ploidy_dose_interaction, method = "BH"),
      ploidy_significant_fdr_0_05 = is.finite(.data$ploidy_fdr) & .data$ploidy_fdr < 0.05,
      dose_significant_fdr_0_05 = is.finite(.data$dose_fdr) & .data$dose_fdr < 0.05,
      interaction_significant_fdr_0_05 = is.finite(.data$ploidy_dose_interaction_fdr) & .data$ploidy_dose_interaction_fdr < 0.05,
      primary_timing_feature = dplyr::case_when(
        .data$model_status != "ok" ~ "Not_tested",
        .data$interaction_significant_fdr_0_05 & .data$abs_interaction_effect_max_delta >= effect_tolerance ~ "Ploidy_Dose_interaction",
        .data$ploidy_significant_fdr_0_05 & !.data$dose_significant_fdr_0_05 ~ "Ploidy_feature",
        .data$dose_significant_fdr_0_05 & !.data$ploidy_significant_fdr_0_05 ~ "Dose_feature",
        .data$ploidy_significant_fdr_0_05 & .data$dose_significant_fdr_0_05 &
          .data$abs_ploidy_effect_4N_minus_2N >= .data$abs_dose_effect_range ~ "Ploidy_feature_with_Dose_effect",
        .data$ploidy_significant_fdr_0_05 & .data$dose_significant_fdr_0_05 &
          .data$abs_dose_effect_range > .data$abs_ploidy_effect_4N_minus_2N ~ "Dose_feature_with_Ploidy_effect",
        !.data$ploidy_significant_fdr_0_05 & !.data$dose_significant_fdr_0_05 &
          .data$abs_ploidy_effect_4N_minus_2N >= effect_tolerance &
          .data$abs_ploidy_effect_4N_minus_2N >= .data$abs_dose_effect_range ~ "Ploidy_effect_size_only",
        !.data$ploidy_significant_fdr_0_05 & !.data$dose_significant_fdr_0_05 &
          .data$abs_dose_effect_range >= effect_tolerance &
          .data$abs_dose_effect_range > .data$abs_ploidy_effect_4N_minus_2N ~ "Dose_effect_size_only",
        TRUE ~ "No_clear_feature"
      ),
      feature_interpretation = dplyr::case_when(
        .data$primary_timing_feature == "Ploidy_Dose_interaction" ~ "Ploidy effect changes by Dose; interpret Ploidy and Dose jointly.",
        .data$primary_timing_feature %in% c("Ploidy_feature", "Ploidy_feature_with_Dose_effect") ~ "Cluster timing is more Ploidy-like.",
        .data$primary_timing_feature %in% c("Dose_feature", "Dose_feature_with_Ploidy_effect") ~ "Cluster timing is more Dose-like.",
        .data$primary_timing_feature == "Ploidy_effect_size_only" ~ "Largest detected effect size is Ploidy-like, but FDR is not significant.",
        .data$primary_timing_feature == "Dose_effect_size_only" ~ "Largest detected effect size is Dose-like, but FDR is not significant.",
        .data$primary_timing_feature == "Not_tested" ~ "Not enough sample-level strata for the Ploidy x Dose model.",
        TRUE ~ "No clear Ploidy-like or Dose-like timing feature detected."
      )
    ) %>%
    dplyr::arrange(factor(as.character(.data$clusters), levels = cluster_levels))
}

make_tumor_ploidy_dose_question_outputs <- function(base_df, metric_col, method_label) {
  tumor <- base_df %>%
    dplyr::filter(
      .data$TN == "Tumor",
      .data$Ploidy %in% c("2N", "4N"),
      .data$Dose %in% c("0mg/kg", "30mg/kg", "120mg/kg"),
      !is.na(.data$clusters),
      is.finite(.data$metric_value)
    )
  empty <- list(
    stratum_summary = data.frame(),
    ploidy_within_dose_cell_tests = data.frame(),
    ploidy_within_dose_sample_tests = data.frame(),
    dose_within_ploidy_cell_tests = data.frame(),
    dose_within_ploidy_sample_tests = data.frame(),
    dose_within_ploidy_pairwise_cell_tests = data.frame(),
    dose_within_ploidy_pairwise_sample_tests = data.frame(),
    feature_summary = data.frame()
  )
  if (nrow(tumor) == 0) return(empty)

  stratum_summary <- tumor %>%
    dplyr::group_by(.data$TN, .data$Ploidy, .data$Dose, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(.data$sample[!is.na(.data$sample)]),
      median_time = stats::median(.data$metric_value, na.rm = TRUE),
      mean_time = mean(.data$metric_value, na.rm = TRUE),
      q25_time = stats::quantile(.data$metric_value, 0.25, na.rm = TRUE, names = FALSE),
      q75_time = stats::quantile(.data$metric_value, 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(method = method_label, metric = metric_col)

  ploidy_pairs <- list(c("2N", "4N"))
  dose_pairs <- list(c("0mg/kg", "30mg/kg"), c("0mg/kg", "120mg/kg"), c("30mg/kg", "120mg/kg"))

  q1_cell <- make_pairwise_tests_by_split(
    tumor,
    split_cols = c("TN", "Dose", "clusters"),
    group_col = "Ploidy",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q1_ploidy_within_dose",
    level = "cell",
    pairs = ploidy_pairs
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(
      question = "Q1_same_Dose_2N_vs_4N",
      median_4N_minus_2N = .data$median_diff_group_2_minus_group_1
    )
  q1_sample_input <- make_sample_level(tumor, c("TN", "Dose", "clusters"), "Ploidy", metric_col)
  q1_sample <- make_pairwise_tests_by_split(
    q1_sample_input,
    split_cols = c("TN", "Dose", "clusters"),
    group_col = "Ploidy",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q1_ploidy_within_dose",
    level = "sample",
    pairs = ploidy_pairs
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(
      question = "Q1_same_Dose_2N_vs_4N",
      median_4N_minus_2N = .data$median_diff_group_2_minus_group_1
    )

  q2_cell <- run_tests_by_split(
    tumor,
    split_cols = c("TN", "Ploidy", "clusters"),
    group_col = "Dose",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q2_dose_within_ploidy",
    level = "cell"
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(question = "Q2_same_Ploidy_Dose_difference")
  q2_sample_input <- make_sample_level(tumor, c("TN", "Ploidy", "clusters"), "Dose", metric_col)
  q2_sample <- run_tests_by_split(
    q2_sample_input,
    split_cols = c("TN", "Ploidy", "clusters"),
    group_col = "Dose",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q2_dose_within_ploidy",
    level = "sample"
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(question = "Q2_same_Ploidy_Dose_difference")

  q2_pairwise_cell <- make_pairwise_tests_by_split(
    tumor,
    split_cols = c("TN", "Ploidy", "clusters"),
    group_col = "Dose",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q2_dose_pairwise_within_ploidy",
    level = "cell",
    pairs = dose_pairs
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(
      question = "Q2_same_Ploidy_pairwise_Dose_difference",
      dose_median_diff_group_2_minus_group_1 = .data$median_diff_group_2_minus_group_1
    )
  q2_pairwise_sample <- make_pairwise_tests_by_split(
    q2_sample_input,
    split_cols = c("TN", "Ploidy", "clusters"),
    group_col = "Dose",
    metric_col = metric_col,
    method = method_label,
    analysis_family = "tumor_Q2_dose_pairwise_within_ploidy",
    level = "sample",
    pairs = dose_pairs
  ) %>%
    add_adjusted_p() %>%
    dplyr::mutate(
      question = "Q2_same_Ploidy_pairwise_Dose_difference",
      dose_median_diff_group_2_minus_group_1 = .data$median_diff_group_2_minus_group_1
    )

  feature_summary <- make_tumor_ploidy_dose_feature_summary(tumor, metric_col, method_label)

  list(
    stratum_summary = stratum_summary,
    ploidy_within_dose_cell_tests = q1_cell,
    ploidy_within_dose_sample_tests = q1_sample,
    dose_within_ploidy_cell_tests = q2_cell,
    dose_within_ploidy_sample_tests = q2_sample,
    dose_within_ploidy_pairwise_cell_tests = q2_pairwise_cell,
    dose_within_ploidy_pairwise_sample_tests = q2_pairwise_sample,
    feature_summary = feature_summary
  )
}

plot_tumor_ploidy_dose_question_outputs <- function(question_outputs, method_label, metric_col, out_dir, plot_prefix) {
  cluster_source <- dplyr::bind_rows(
    question_outputs$stratum_summary %>% dplyr::select(dplyr::any_of("clusters")),
    question_outputs$feature_summary %>% dplyr::select(dplyr::any_of("clusters"))
  )
  if (nrow(cluster_source) == 0) return(invisible(NULL))
  cluster_levels <- sort_maybe_numeric(unique(as.character(cluster_source$clusters)))

  q1_plot <- question_outputs$ploidy_within_dose_sample_tests
  q1_level <- "sample"
  if (nrow(q1_plot) == 0) {
    q1_plot <- question_outputs$ploidy_within_dose_cell_tests
    q1_level <- "cell"
  }
  if (nrow(q1_plot) > 0) {
    q1_plot <- q1_plot %>%
      dplyr::filter(.data$contrast == "4N_vs_2N") %>%
      dplyr::mutate(
        clusters = factor(.data$clusters, levels = cluster_levels),
        Dose = factor(.data$Dose, levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
        label = ifelse(
          dplyr::coalesce(.data$significant_family_fdr_0_05, FALSE),
          paste0(sprintf("%.3f", .data$median_4N_minus_2N), "*"),
          sprintf("%.3f", .data$median_4N_minus_2N)
        )
      )
    p <- ggplot(q1_plot, aes(x = Dose, y = clusters, fill = median_4N_minus_2N)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.5) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "4N - 2N\nmedian time") +
      labs(
        title = paste0(method_label, " tumor Q1: Ploidy difference within each Dose"),
        subtitle = paste0("Labels use ", q1_level, "-level tests; * family FDR < 0.05."),
        x = "Dose",
        y = "Cluster"
      ) +
      theme_extra(9) +
      theme(legend.position = "right")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_Q1_ploidy_within_dose_heatmap")), width = 8.8, height = 5.8)
  }

  q2_plot <- question_outputs$dose_within_ploidy_sample_tests
  q2_level <- "sample"
  if (nrow(q2_plot) == 0) {
    q2_plot <- question_outputs$dose_within_ploidy_cell_tests
    q2_level <- "cell"
  }
  if (nrow(q2_plot) > 0) {
    q2_plot <- q2_plot %>%
      dplyr::mutate(
        clusters = factor(.data$clusters, levels = cluster_levels),
        Ploidy = factor(.data$Ploidy, levels = c("2N", "4N")),
        label = ifelse(
          dplyr::coalesce(.data$significant_family_fdr_0_05, FALSE),
          paste0(sprintf("%.3f", .data$median_range), "*"),
          sprintf("%.3f", .data$median_range)
        )
      )
    p <- ggplot(q2_plot, aes(x = Ploidy, y = clusters, fill = median_range)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.5) +
      scale_fill_gradient(low = "white", high = "#2A7F62", name = "Dose median\nrange") +
      labs(
        title = paste0(method_label, " tumor Q2: Dose difference within each Ploidy"),
        subtitle = paste0("Kruskal-Wallis by cluster and Ploidy using ", q2_level, "-level values; * family FDR < 0.05."),
        x = "Ploidy",
        y = "Cluster"
      ) +
      theme_extra(9) +
      theme(legend.position = "right")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_dose_within_ploidy_heatmap")), width = 7.4, height = 5.8)
  }

  q2_pairwise <- question_outputs$dose_within_ploidy_pairwise_sample_tests
  q2_pairwise_level <- "sample"
  if (nrow(q2_pairwise) == 0) {
    q2_pairwise <- question_outputs$dose_within_ploidy_pairwise_cell_tests
    q2_pairwise_level <- "cell"
  }
  if (nrow(q2_pairwise) > 0) {
    contrast_levels <- c("30mg/kg_vs_0mg/kg", "120mg/kg_vs_0mg/kg", "120mg/kg_vs_30mg/kg")
    q2_pairwise <- q2_pairwise %>%
      dplyr::mutate(
        clusters = factor(.data$clusters, levels = cluster_levels),
        Ploidy = factor(.data$Ploidy, levels = c("2N", "4N")),
        contrast = factor(.data$contrast, levels = contrast_levels),
        label = ifelse(
          dplyr::coalesce(.data$significant_family_fdr_0_05, FALSE),
          paste0(sprintf("%.3f", .data$dose_median_diff_group_2_minus_group_1), "*"),
          sprintf("%.3f", .data$dose_median_diff_group_2_minus_group_1)
        )
      )
    p <- ggplot(q2_pairwise, aes(x = contrast, y = clusters, fill = dose_median_diff_group_2_minus_group_1)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.3) +
      facet_wrap(~Ploidy, nrow = 1) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Dose contrast\nmedian time") +
      labs(
        title = paste0(method_label, " tumor Q2: pairwise Dose contrasts within each Ploidy"),
        subtitle = paste0("Labels use ", q2_pairwise_level, "-level Wilcoxon tests; * family FDR < 0.05."),
        x = "Dose contrast",
        y = "Cluster"
      ) +
      theme_extra(8) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_pairwise_dose_within_ploidy_heatmap")), width = 11, height = 5.8)
  }

  q3 <- question_outputs$feature_summary
  if (nrow(q3) > 0) {
    fdr_long <- q3 %>%
      dplyr::select(clusters, ploidy_fdr, dose_fdr, ploidy_dose_interaction_fdr) %>%
      tidyr::pivot_longer(
        cols = c("ploidy_fdr", "dose_fdr", "ploidy_dose_interaction_fdr"),
        names_to = "effect_type",
        values_to = "fdr"
      ) %>%
      dplyr::mutate(
        clusters = factor(.data$clusters, levels = cluster_levels),
        effect_type = dplyr::recode(
          .data$effect_type,
          ploidy_fdr = "Ploidy",
          dose_fdr = "Dose",
          ploidy_dose_interaction_fdr = "Ploidy:Dose"
        ),
        effect_type = factor(.data$effect_type, levels = c("Ploidy", "Dose", "Ploidy:Dose")),
        neg_log10_fdr = dplyr::if_else(is.finite(.data$fdr), -log10(pmax(.data$fdr, .Machine$double.xmin)), NA_real_),
        label = dplyr::case_when(
          !is.finite(.data$fdr) ~ "NA",
          .data$fdr < 0.001 ~ paste0(signif(.data$fdr, 2), "***"),
          .data$fdr < 0.01 ~ paste0(signif(.data$fdr, 2), "**"),
          .data$fdr < 0.05 ~ paste0(signif(.data$fdr, 2), "*"),
          TRUE ~ as.character(signif(.data$fdr, 2))
        )
      )
    p_fdr <- ggplot(fdr_long, aes(x = effect_type, y = clusters, fill = neg_log10_fdr)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.4) +
      scale_fill_gradient(low = "white", high = "#5A4A8F", na.value = "grey92", name = "-log10(FDR)") +
      labs(
        title = paste0(method_label, " tumor Q3: Ploidy, Dose, and interaction model terms"),
        subtitle = "Sample-level median model: time ~ Ploidy + Dose plus Ploidy:Dose interaction test.",
        x = "Model term",
        y = "Cluster"
      ) +
      theme_extra(9) +
      theme(legend.position = "right")
    save_both(p_fdr, file.path(out_dir, paste0(plot_prefix, "_tumor_Q3_ploidy_dose_model_fdr_heatmap")), width = 8.5, height = 5.8)

    feature_levels <- c(
      "Ploidy_feature",
      "Dose_feature",
      "Ploidy_Dose_interaction",
      "Ploidy_feature_with_Dose_effect",
      "Dose_feature_with_Ploidy_effect",
      "Ploidy_effect_size_only",
      "Dose_effect_size_only",
      "No_clear_feature",
      "Not_tested"
    )
    q3_plot <- q3 %>%
      dplyr::mutate(
        clusters = factor(.data$clusters, levels = cluster_levels),
        primary_timing_feature = factor(.data$primary_timing_feature, levels = feature_levels)
      )
    p_feature <- ggplot(q3_plot, aes(x = primary_timing_feature, y = clusters, fill = primary_timing_feature)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_manual(
        values = c(
          "Ploidy_feature" = "#3B73B9",
          "Dose_feature" = "#2A7F62",
          "Ploidy_Dose_interaction" = "#B23A48",
          "Ploidy_feature_with_Dose_effect" = "#7EA7D8",
          "Dose_feature_with_Ploidy_effect" = "#7DBE9E",
          "Ploidy_effect_size_only" = "#B7CAE8",
          "Dose_effect_size_only" = "#B7DEC9",
          "No_clear_feature" = "grey78",
          "Not_tested" = "grey92"
        ),
        drop = FALSE,
        name = "Primary feature"
      ) +
      labs(
        title = paste0(method_label, " tumor Q3: inferred timing feature"),
        x = "Primary timing feature",
        y = "Cluster"
      ) +
      theme_extra(8) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom")
    save_both(p_feature, file.path(out_dir, paste0(plot_prefix, "_tumor_Q3_primary_timing_feature")), width = 12, height = 5.8)
  }

  invisible(NULL)
}

make_tn_ploidy_key_conclusions <- function(comparison_df, dose_interaction_df, method_label, metric_col) {
  key_ploidy <- comparison_df %>%
    dplyr::filter(
      is.finite(.data$median_4N_minus_2N),
      (.data$n_cells_2N >= 20 | is.na(.data$n_cells_2N)),
      (.data$n_cells_4N >= 20 | is.na(.data$n_cells_4N))
    ) %>%
    dplyr::mutate(
      evidence_type = "within_TN_cluster_2N_vs_4N",
      key_effect = .data$median_4N_minus_2N,
      key_effect_abs = abs(.data$key_effect),
      interpretation = dplyr::case_when(
        .data$ploidy_timing_direction == "4N_later_or_higher_time" ~ "4N has later/higher cluster time than 2N in this stratum.",
        .data$ploidy_timing_direction == "2N_later_or_higher_time" ~ "2N has later/higher cluster time than 4N in this stratum.",
        .data$ploidy_timing_direction == "similar_within_0.02" ~ "2N and 4N cluster times are similar within 0.02.",
        TRUE ~ "Not enough cells to interpret this stratum."
      )
    ) %>%
    dplyr::select(
      method, metric, evidence_type, comparison_scope, TN,
      Dose, clusters, key_effect, key_effect_abs,
      ploidy_timing_direction, ploidy_fdr_family,
      ploidy_significant_family_fdr_0_05, interpretation
    )

  key_interaction <- dose_interaction_df %>%
    dplyr::filter(is.finite(.data$interaction_delta_4N_minus_2N)) %>%
    dplyr::mutate(
      evidence_type = "Tumor_dose_response_interaction",
      comparison_scope = "Tumor_Dose_delta_4N_minus_2N",
      TN = "Tumor",
      key_effect = .data$interaction_delta_4N_minus_2N,
      key_effect_abs = abs(.data$key_effect),
      ploidy_timing_direction = .data$interaction_pattern,
      ploidy_fdr_family = NA_real_,
      ploidy_significant_family_fdr_0_05 = NA,
      interpretation = dplyr::case_when(
        .data$interaction_pattern == "2N_shifts_earlier_more_than_4N" ~ "Treatment moves 2N cluster time earlier/lower more strongly than 4N.",
        .data$interaction_pattern == "4N_shifts_earlier_more_than_2N" ~ "Treatment moves 4N cluster time earlier/lower more strongly than 2N.",
        .data$interaction_pattern == "similar_dose_response_within_0.02" ~ "2N and 4N have similar treatment-associated timing shifts within 0.02.",
        TRUE ~ "Not enough groups to interpret this dose response contrast."
      )
    ) %>%
    dplyr::mutate(Dose = .data$dose_vs_baseline) %>%
    dplyr::select(
      method, metric, evidence_type, comparison_scope, TN,
      Dose, clusters, key_effect, key_effect_abs,
      ploidy_timing_direction, ploidy_fdr_family,
      ploidy_significant_family_fdr_0_05, interpretation
    )

  dplyr::bind_rows(key_ploidy, key_interaction) %>%
    dplyr::mutate(method = method_label, metric = metric_col) %>%
    dplyr::arrange(dplyr::desc(.data$key_effect_abs), .data$comparison_scope, .data$clusters)
}

plot_tn_ploidy_cluster_response <- function(base_df, comparison_df, dose_delta_df, dose_interaction_df, method_label, metric_col, out_dir, plot_prefix) {
  cluster_levels <- sort_maybe_numeric(base_df$clusters)

  tumor_plot <- base_df %>%
    dplyr::filter(.data$TN == "Tumor", !is.na(.data$Dose)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      Dose = factor(.data$Dose, levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      Ploidy = factor(.data$Ploidy, levels = c("2N", "4N"))
    ) %>%
    dplyr::add_count(.data$clusters, .data$Ploidy, .data$Dose, name = ".stratum_n") %>%
    dplyr::filter(.data$.stratum_n >= 2)
  if (nrow(tumor_plot) > 0) {
    stat_labels <- make_violin_stat_annotations(
      tumor_plot,
      split_cols = c("clusters", "Dose"),
      group_col = "Ploidy",
      value_col = "metric_value",
      x_col = "Dose",
      min_group_n = 3L,
      label_prefix = "2N vs 4N"
    )
    p <- ggplot(tumor_plot, aes(x = Dose, y = metric_value, fill = Ploidy)) +
      geom_violin(position = position_dodge(width = 0.82), scale = "width", trim = TRUE, alpha = 0.78, linewidth = 0.15, color = "grey30") +
      geom_boxplot(position = position_dodge(width = 0.82), width = 0.15, outlier.shape = NA, alpha = 0.8, linewidth = 0.15) +
      facet_wrap(~clusters, ncol = 3, scales = "free_y") +
      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
      labs(
        title = paste0(method_label, " tumor cluster time distributions by Ploidy and Dose"),
        subtitle = "Each panel is one cluster; labels report Wilcoxon tests comparing 2N vs 4N within each Dose.",
        x = "Dose",
        y = metric_col,
        fill = "Ploidy"
      ) +
      geom_text(
        data = stat_labels,
        aes(x = x_plot, y = y_plot, label = label),
        inherit.aes = FALSE,
        size = 2.15,
        lineheight = 0.86,
        fontface = "bold",
        color = "black"
      ) +
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.28))) +
      theme_extra(9) +
      theme(legend.position = "bottom")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_cluster_ploidy_dose_distribution")), width = 11, height = 10.5)
  }

  cellline_plot <- base_df %>%
    dplyr::filter(.data$TN == "CellLine") %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      Ploidy = factor(.data$Ploidy, levels = c("2N", "4N"))
    ) %>%
    dplyr::add_count(.data$clusters, .data$Ploidy, name = ".stratum_n") %>%
    dplyr::filter(.data$.stratum_n >= 2)
  if (nrow(cellline_plot) > 0) {
    stat_labels <- make_violin_stat_annotations(
      cellline_plot,
      split_cols = c("clusters"),
      group_col = "Ploidy",
      value_col = "metric_value",
      x_value = "2N",
      min_group_n = 3L,
      label_prefix = "2N vs 4N"
    )
    p <- ggplot(cellline_plot, aes(x = Ploidy, y = metric_value, fill = Ploidy)) +
      geom_violin(scale = "width", trim = TRUE, alpha = 0.78, linewidth = 0.15, color = "grey30") +
      geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.8, linewidth = 0.15) +
      facet_wrap(~clusters, ncol = 3, scales = "free_y") +
      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
      labs(
        title = paste0(method_label, " cell-line cluster time distributions by Ploidy"),
        subtitle = "Labels report Wilcoxon tests comparing 2N vs 4N within each cluster.",
        x = "Ploidy",
        y = metric_col,
        fill = "Ploidy"
      ) +
      geom_text(
        data = stat_labels,
        aes(x = x_plot, y = y_plot, label = label),
        inherit.aes = FALSE,
        size = 2.35,
        lineheight = 0.88,
        hjust = -0.05,
        fontface = "bold",
        color = "black"
      ) +
      scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.24))) +
      theme_extra(9) +
      theme(legend.position = "bottom")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_cellline_cluster_ploidy_distribution")), width = 10.5, height = 10)
  }

  tumor_delta <- comparison_df %>%
    dplyr::filter(.data$comparison_scope == "Tumor_2N_vs_4N_by_Dose", is.finite(.data$median_4N_minus_2N)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      Dose = factor(.data$Dose, levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      label = ifelse(
        dplyr::coalesce(.data$ploidy_significant_family_fdr_0_05, FALSE),
        paste0(sprintf("%.3f", .data$median_4N_minus_2N), "*"),
        sprintf("%.3f", .data$median_4N_minus_2N)
      )
    )
  if (nrow(tumor_delta) > 0) {
    p <- ggplot(tumor_delta, aes(x = Dose, y = clusters, fill = median_4N_minus_2N)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.6) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "4N - 2N\nmedian time") +
      labs(
        title = paste0(method_label, " tumor within-cluster Ploidy timing difference"),
        subtitle = "Positive values mean 4N is later/higher; negative values mean 2N is later/higher. * family FDR < 0.05.",
        x = "Dose",
        y = "Cluster"
      ) +
      theme_extra(9) +
      theme(legend.position = "right")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_cluster_4N_minus_2N_by_dose_heatmap")), width = 8.8, height = 5.8)
  }

  cellline_delta <- comparison_df %>%
    dplyr::filter(.data$comparison_scope == "CellLine_2N_vs_4N", is.finite(.data$median_4N_minus_2N)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      significant = ifelse(dplyr::coalesce(.data$ploidy_significant_family_fdr_0_05, FALSE), "FDR < 0.05", "not FDR < 0.05")
    )
  if (nrow(cellline_delta) > 0) {
    p <- ggplot(cellline_delta, aes(x = clusters, y = median_4N_minus_2N, fill = significant)) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.25) +
      geom_col(width = 0.72) +
      scale_fill_manual(values = c("FDR < 0.05" = "#4B8BBE", "not FDR < 0.05" = "grey75"), drop = FALSE) +
      labs(
        title = paste0(method_label, " cell-line within-cluster Ploidy timing difference"),
        subtitle = "Positive values mean 4N is later/higher; negative values mean 2N is later/higher.",
        x = "Cluster",
        y = "4N - 2N median time",
        fill = "Ploidy test"
      ) +
      theme_extra(10) +
      theme(legend.position = "bottom")
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_cellline_cluster_4N_minus_2N_barplot")), width = 8.5, height = 5.4)
  }

  dose_delta_plot <- dose_delta_df %>%
    dplyr::filter(is.finite(.data$median_delta_vs_0)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      Ploidy = factor(.data$Ploidy, levels = c("2N", "4N")),
      dose_vs_baseline = factor(.data$dose_vs_baseline, levels = c("30mg/kg_vs_0mg/kg", "120mg/kg_vs_0mg/kg")),
      label = sprintf("%.3f", .data$median_delta_vs_0)
    )
  if (nrow(dose_delta_plot) > 0) {
    p <- ggplot(dose_delta_plot, aes(x = Ploidy, y = clusters, fill = median_delta_vs_0)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.5) +
      facet_wrap(~dose_vs_baseline, nrow = 1) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Dose - 0mg/kg\nmedian time") +
      labs(
        title = paste0(method_label, " tumor dose-associated cluster timing shift"),
        subtitle = "Negative values mean treatment shifts cells earlier/lower relative to 0mg/kg within the same cluster and Ploidy.",
        x = "Ploidy",
        y = "Cluster"
      ) +
      theme_extra(9)
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_dose_delta_from_0_heatmap")), width = 8.4, height = 5.8)
  }

  interaction_plot <- dose_interaction_df %>%
    dplyr::filter(is.finite(.data$interaction_delta_4N_minus_2N)) %>%
    dplyr::mutate(
      clusters = factor(.data$clusters, levels = cluster_levels),
      dose_vs_baseline = factor(.data$dose_vs_baseline, levels = c("30mg/kg_vs_0mg/kg", "120mg/kg_vs_0mg/kg")),
      label = sprintf("%.3f", .data$interaction_delta_4N_minus_2N)
    )
  if (nrow(interaction_plot) > 0) {
    p <- ggplot(interaction_plot, aes(x = dose_vs_baseline, y = clusters, fill = interaction_delta_4N_minus_2N)) +
      geom_tile(color = "white", linewidth = 0.25) +
      geom_text(aes(label = label), size = 2.5) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "(4N shift) -\n(2N shift)") +
      labs(
        title = paste0(method_label, " tumor Dose x Ploidy pseudotime interaction"),
        subtitle = "Positive values mean 2N shifts earlier/lower more than 4N; negative values mean 4N shifts earlier/lower more than 2N.",
        x = "Dose contrast",
        y = "Cluster"
      ) +
      theme_extra(9)
    save_both(p, file.path(out_dir, paste0(plot_prefix, "_tumor_dose_ploidy_interaction_heatmap")), width = 8.2, height = 5.8)
  }
}

write_tn_ploidy_cluster_outputs <- function(desc_df, metric_col, method_label, out_dir, plot_prefix, cluster_annotations = data.frame()) {
  base <- make_tn_ploidy_base(desc_df, metric_col)
  summary_df <- summarize_tn_ploidy_clusters(base, metric_col)
  summary_for_calc <- summary_df
  if (nrow(cluster_annotations) > 0 && nrow(summary_df) > 0) {
    summary_df <- summary_df %>%
      dplyr::left_join(
        cluster_annotations %>%
          dplyr::select(-dplyr::any_of(c("n_cells", "n_samples", "mean_time", "median_time", "q25_time", "q75_time", "metric"))),
        by = c("clusters" = "cluster")
      )
  }

  comparisons <- make_tn_ploidy_cluster_comparisons(base, summary_for_calc, metric_col, method_label)
  dose_response <- make_tumor_dose_response_deltas(summary_for_calc, metric_col, method_label)
  tumor_ploidy_dose <- make_tumor_ploidy_dose_question_outputs(base, metric_col, method_label)
  key_conclusions <- make_tn_ploidy_key_conclusions(comparisons, dose_response$interaction, method_label, metric_col)

  write_table_csv(summary_df, file.path(out_dir, paste0(plot_prefix, "_tn_ploidy_cluster_time_summary.csv")))
  write_table_csv(comparisons, file.path(out_dir, paste0(plot_prefix, "_tn_cluster_2N_vs_4N_comparisons.csv")))
  write_table_csv(dose_response$dose_delta, file.path(out_dir, paste0(plot_prefix, "_tumor_cluster_dose_delta_from_0.csv")))
  write_table_csv(dose_response$interaction, file.path(out_dir, paste0(plot_prefix, "_tumor_cluster_dose_ploidy_interaction.csv")))
  write_table_csv(tumor_ploidy_dose$stratum_summary, file.path(out_dir, paste0(plot_prefix, "_tumor_Q0_ploidy_dose_cluster_stratum_summary.csv")))
  write_table_csv(tumor_ploidy_dose$ploidy_within_dose_cell_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q1_ploidy_within_dose_cell_tests.csv")))
  write_table_csv(tumor_ploidy_dose$ploidy_within_dose_sample_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q1_ploidy_within_dose_sample_tests.csv")))
  write_table_csv(tumor_ploidy_dose$dose_within_ploidy_cell_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_dose_within_ploidy_cell_tests.csv")))
  write_table_csv(tumor_ploidy_dose$dose_within_ploidy_sample_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_dose_within_ploidy_sample_tests.csv")))
  write_table_csv(tumor_ploidy_dose$dose_within_ploidy_pairwise_cell_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_dose_pairwise_within_ploidy_cell_tests.csv")))
  write_table_csv(tumor_ploidy_dose$dose_within_ploidy_pairwise_sample_tests, file.path(out_dir, paste0(plot_prefix, "_tumor_Q2_dose_pairwise_within_ploidy_sample_tests.csv")))
  write_table_csv(tumor_ploidy_dose$feature_summary, file.path(out_dir, paste0(plot_prefix, "_tumor_Q3_ploidy_dose_feature_summary.csv")))
  write_table_csv(key_conclusions, file.path(out_dir, paste0(plot_prefix, "_tn_ploidy_cluster_key_conclusions.csv")))

  plot_tn_ploidy_cluster_response(
    base,
    comparisons,
    dose_response$dose_delta,
    dose_response$interaction,
	    method_label,
	    metric_col,
	    out_dir,
	    plot_prefix
	  )
  plot_tumor_ploidy_dose_question_outputs(
    tumor_ploidy_dose,
    method_label,
    metric_col,
    out_dir,
    plot_prefix
  )

	  list(
	    base = base,
	    summary = summary_df,
	    comparisons = comparisons,
	    dose_delta = dose_response$dose_delta,
	    dose_interaction = dose_response$interaction,
	    tumor_ploidy_dose = tumor_ploidy_dose,
	    key_conclusions = key_conclusions
	  )
	}

min_finite_value <- function(x) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  min(x)
}

collapse_cluster_values <- function(x) {
  x <- sort_maybe_numeric(unique(as_clean_chr(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0) return("none")
  paste(x, collapse = ", ")
}

max_abs_effect_row <- function(df, effect_col, signed_effect_col = effect_col) {
  df <- ensure_columns(df, c("clusters", effect_col, signed_effect_col), default = NA_real_)
  df[[effect_col]] <- safe_numeric(df[[effect_col]])
  df[[signed_effect_col]] <- safe_numeric(df[[signed_effect_col]])
  df <- df[is.finite(df[[effect_col]]), , drop = FALSE]
  if (nrow(df) == 0) {
    return(data.frame(cluster = NA_character_, abs_effect = NA_real_, signed_effect = NA_real_, stringsAsFactors = FALSE))
  }
  df <- df[order(df[[effect_col]], decreasing = TRUE), , drop = FALSE]
  data.frame(
    cluster = as.character(df$clusters[1]),
    abs_effect = as.numeric(df[[effect_col]][1]),
    signed_effect = as.numeric(df[[signed_effect_col]][1]),
    stringsAsFactors = FALSE
  )
}

make_tumor_method_conclusion_overview <- function(method_label, result_obj) {
  q1 <- ensure_columns(
    result_obj$tumor_ploidy_dose$ploidy_within_dose_sample_tests,
    c("primary_p_adj_family", "significant_family_fdr_0_05")
  )
  q2 <- ensure_columns(
    result_obj$tumor_ploidy_dose$dose_within_ploidy_sample_tests,
    c("primary_p_adj_family", "significant_family_fdr_0_05")
  )
  q2_pair <- ensure_columns(
    result_obj$tumor_ploidy_dose$dose_within_ploidy_pairwise_sample_tests,
    c("primary_p_adj_family", "significant_family_fdr_0_05")
  )
  q3 <- ensure_columns(
    result_obj$tumor_ploidy_dose$feature_summary,
    c(
      "clusters", "model_status", "primary_timing_feature", "ploidy_fdr", "dose_fdr",
      "ploidy_dose_interaction_fdr", "abs_ploidy_effect_4N_minus_2N",
      "ploidy_effect_4N_minus_2N", "abs_dose_effect_range", "dose_effect_range",
      "abs_interaction_effect_max_delta", "interaction_effect_max_delta"
    )
  )

  ploidy_features <- c("Ploidy_feature", "Ploidy_feature_with_Dose_effect", "Ploidy_effect_size_only")
  dose_features <- c("Dose_feature", "Dose_feature_with_Ploidy_effect", "Dose_effect_size_only")
  interaction_features <- c("Ploidy_Dose_interaction")
  max_ploidy <- max_abs_effect_row(q3, "abs_ploidy_effect_4N_minus_2N", "ploidy_effect_4N_minus_2N")
  max_dose <- max_abs_effect_row(q3, "abs_dose_effect_range", "dose_effect_range")
  max_interaction <- max_abs_effect_row(q3, "abs_interaction_effect_max_delta", "interaction_effect_max_delta")

  data.frame(
    method = method_label,
    q1_question = "Same Dose: 2N vs 4N",
    q1_sample_tests = nrow(q1),
    q1_sample_fdr_significant_n = sum(q1$significant_family_fdr_0_05 %in% TRUE, na.rm = TRUE),
    q1_sample_min_fdr = min_finite_value(q1$primary_p_adj_family),
    q2_question = "Same Ploidy: Dose difference",
    q2_sample_tests = nrow(q2),
    q2_sample_fdr_significant_n = sum(q2$significant_family_fdr_0_05 %in% TRUE, na.rm = TRUE),
    q2_sample_min_fdr = min_finite_value(q2$primary_p_adj_family),
    q2_pairwise_sample_tests = nrow(q2_pair),
    q2_pairwise_sample_fdr_significant_n = sum(q2_pair$significant_family_fdr_0_05 %in% TRUE, na.rm = TRUE),
    q2_pairwise_sample_min_fdr = min_finite_value(q2_pair$primary_p_adj_family),
    q3_question = "Ploidy-like, Dose-like, or interaction-like timing",
    q3_model_ok_clusters = sum(q3$model_status == "ok", na.rm = TRUE),
    q3_not_tested_clusters = collapse_cluster_values(q3$clusters[q3$primary_timing_feature == "Not_tested"]),
    q3_ploidy_like_clusters = collapse_cluster_values(q3$clusters[q3$primary_timing_feature %in% ploidy_features]),
    q3_dose_like_clusters = collapse_cluster_values(q3$clusters[q3$primary_timing_feature %in% dose_features]),
    q3_interaction_like_clusters = collapse_cluster_values(q3$clusters[q3$primary_timing_feature %in% interaction_features]),
    q3_no_clear_clusters = collapse_cluster_values(q3$clusters[q3$primary_timing_feature == "No_clear_feature"]),
    q3_min_ploidy_fdr = min_finite_value(q3$ploidy_fdr),
    q3_min_dose_fdr = min_finite_value(q3$dose_fdr),
    q3_min_interaction_fdr = min_finite_value(q3$ploidy_dose_interaction_fdr),
    max_abs_ploidy_effect_cluster = max_ploidy$cluster,
    max_abs_ploidy_effect = max_ploidy$abs_effect,
    max_signed_ploidy_effect_4N_minus_2N = max_ploidy$signed_effect,
    max_abs_dose_effect_cluster = max_dose$cluster,
    max_abs_dose_effect = max_dose$abs_effect,
    max_signed_dose_effect_range = max_dose$signed_effect,
    max_abs_interaction_effect_cluster = max_interaction$cluster,
    max_abs_interaction_effect = max_interaction$abs_effect,
    max_signed_interaction_effect = max_interaction$signed_effect,
    stringsAsFactors = FALSE
  )
}

make_tumor_conclusion_markdown <- function(method_overview, q3_all) {
  scvelo <- method_overview %>% dplyr::filter(.data$method == method_scvelo)
  monocle3 <- method_overview %>% dplyr::filter(.data$method == method_monocle3)
  lines <- c(
    "# Tumor Ploidy x Dose Timing Summary",
    "",
    "## Interpretation rule",
    "- Primary evidence uses sample-level tests, where each sample is summarized first.",
    "- Cell-level tests are more powered but cells are not biological replicates; they are treated as supporting evidence only.",
    "- Feature labels in Q3 are effect-size classifications unless the corresponding FDR is below 0.05.",
    "",
    "## Question definitions",
    "- Q1, same Dose: within Tumor cells, compare 2N versus 4N pseudotime within each Dose and cluster.",
    "- Q2, same Ploidy: within Tumor cells, compare Dose groups within each Ploidy and cluster. This includes an omnibus Dose test and pairwise Dose contrasts.",
    "- Q3, timing feature: within each Tumor cluster, classify whether the timing pattern is more Ploidy-like, Dose-like, Ploidy:Dose interaction-like, or not clearly attributable to either factor, using sample-level Ploidy, Dose, and Ploidy:Dose model terms plus effect sizes.",
    "",
    "## scVelo conclusion",
    paste0("- Q1 same-Dose 2N vs 4N sample-level FDR-significant tests: ", scvelo$q1_sample_fdr_significant_n, "."),
    paste0("- Q2 same-Ploidy Dose sample-level FDR-significant tests: ", scvelo$q2_sample_fdr_significant_n, "."),
    paste0("- Q3 Dose-like clusters: ", scvelo$q3_dose_like_clusters, "."),
    paste0("- Q3 Ploidy-like clusters: ", scvelo$q3_ploidy_like_clusters, "."),
    paste0("- Q3 interaction-like clusters: ", scvelo$q3_interaction_like_clusters, "."),
    paste0("- Q3 not-tested clusters: ", scvelo$q3_not_tested_clusters, "."),
    "- Conservative conclusion: scVelo does not show sample-level FDR-significant Ploidy, Dose, or interaction effects, but the effect-size pattern is more Dose-like than Ploidy-like.",
    "",
    "## Monocle3 conclusion",
    paste0("- Q1 same-Dose 2N vs 4N sample-level FDR-significant tests: ", monocle3$q1_sample_fdr_significant_n, "."),
    paste0("- Q2 same-Ploidy Dose sample-level FDR-significant tests: ", monocle3$q2_sample_fdr_significant_n, "."),
    paste0("- Q3 Dose-like clusters: ", monocle3$q3_dose_like_clusters, "."),
    paste0("- Q3 Ploidy-like clusters: ", monocle3$q3_ploidy_like_clusters, "."),
    paste0("- Q3 interaction-like clusters: ", monocle3$q3_interaction_like_clusters, "."),
    paste0("- Q3 no-clear clusters: ", monocle3$q3_no_clear_clusters, "."),
    paste0("- Q3 not-tested clusters: ", monocle3$q3_not_tested_clusters, "."),
    "- Conservative conclusion: Monocle3 does not show sample-level FDR-significant evidence for Ploidy, Dose, or Ploidy:Dose interaction. The largest effect-size signal is Ploidy-like in cluster 10, while cluster 6 is Dose-like; the overall Monocle3 pattern is weaker and less broadly Dose-like than scVelo.",
    "",
    "## Cross-method readout",
    "- Both methods agree that sample-level evidence is not FDR-significant for the three Tumor questions.",
    "- scVelo suggests a broader Dose-associated timing pattern; Monocle3 suggests a narrower pattern with cluster 10 more Ploidy-like and cluster 6 more Dose-like.",
    "- Neither method supports a robust Ploidy:Dose interaction after sample-level modeling.",
    "",
    "## Files in this folder",
    "- tumor_ploidy_dose_method_overview.csv: one-row-per-method summary for Q1/Q2/Q3.",
    "- tumor_ploidy_dose_feature_summary_all_methods.csv: cluster-level Q3 model and feature calls.",
    "- supporting_tables/: copied CSV tables that support the conclusions.",
    "- supporting_figures/: copied PDF/PNG figures that support the conclusions."
  )
  lines
}

make_tumor_conclusion_markdown_zh <- function(method_overview, q3_all) {
  scvelo <- method_overview %>% dplyr::filter(.data$method == method_scvelo)
  monocle3 <- method_overview %>% dplyr::filter(.data$method == method_monocle3)
  c(
    "# Tumor Ploidy x Dose 时序总结",
    "",
    "## 解读原则",
    "- 主要结论以样本层级检验为准：先在每个样本内取中位数，再进行统计检验。",
    "- cell-level 检验只作为辅助证据，因为细胞不是独立生物学重复，细胞数大时很容易显著。",
    "- Q3 的 Ploidy-like / Dose-like / interaction-like 标签主要是效应量分类；只有对应 FDR < 0.05 时才作为显著统计证据。",
    "",
    "## scVelo 结论",
    paste0("- Q1：同一 Dose 下 2N vs 4N，样本层级 FDR 显著数 = ", scvelo$q1_sample_fdr_significant_n, "。"),
    paste0("- Q2：同一 Ploidy 下不同 Dose，样本层级 FDR 显著数 = ", scvelo$q2_sample_fdr_significant_n, "。"),
    paste0("- Q3：Dose-like clusters = ", scvelo$q3_dose_like_clusters, "。"),
    paste0("- Q3：Ploidy-like clusters = ", scvelo$q3_ploidy_like_clusters, "。"),
    paste0("- Q3：interaction-like clusters = ", scvelo$q3_interaction_like_clusters, "。"),
    paste0("- Q3：未能检验的 clusters = ", scvelo$q3_not_tested_clusters, "。"),
    "- 保守结论：scVelo 没有样本层级 FDR 显著的 Ploidy、Dose 或 Ploidy:Dose interaction 证据；但效应量模式整体更偏 Dose-associated timing，而不是稳定的 Ploidy-specific timing。",
    "",
    "## Monocle3 结论",
    paste0("- Q1：同一 Dose 下 2N vs 4N，样本层级 FDR 显著数 = ", monocle3$q1_sample_fdr_significant_n, "。"),
    paste0("- Q2：同一 Ploidy 下不同 Dose，样本层级 FDR 显著数 = ", monocle3$q2_sample_fdr_significant_n, "。"),
    paste0("- Q3：Dose-like clusters = ", monocle3$q3_dose_like_clusters, "。"),
    paste0("- Q3：Ploidy-like clusters = ", monocle3$q3_ploidy_like_clusters, "。"),
    paste0("- Q3：interaction-like clusters = ", monocle3$q3_interaction_like_clusters, "。"),
    paste0("- Q3：no-clear clusters = ", monocle3$q3_no_clear_clusters, "。"),
    paste0("- Q3：未能检验的 clusters = ", monocle3$q3_not_tested_clusters, "。"),
    "- 保守结论：Monocle3 没有样本层级 FDR 显著的 Ploidy、Dose 或 Ploidy:Dose interaction 证据。效应量上，cluster 10 更偏 Ploidy-like，cluster 6 更偏 Dose-like；整体信号比 scVelo 更弱，也没有 scVelo 那种较广泛的 Dose-like 模式。",
    "",
    "## 两种方法的共同结论",
    "- 两种方法都不支持样本层级显著的 2N/4N 差异、Dose 差异或 Ploidy:Dose interaction。",
    "- scVelo 更偏向提示广泛的 Dose-associated timing；Monocle3 的提示更局限，主要集中在 cluster 10 和 cluster 6。",
    "- 因此最终表述应当是：有 Dose/Ploidy 相关的效应量趋势，但缺乏样本层级 FDR 显著统计证据。",
    "",
    "## 本文件夹内容",
    "- tumor_ploidy_dose_method_overview.csv：每个方法一行的 Q1/Q2/Q3 总结。",
    "- tumor_ploidy_dose_feature_summary_all_methods.csv：每个 cluster 的 Q3 模型与特征分类。",
    "- supporting_tables/：支撑结论的 CSV 表格拷贝。",
    "- supporting_figures/：支撑结论的 PDF/PNG 图拷贝。"
  )
}

copy_summary_supporting_files <- function(output_root, summary_dir) {
  table_root <- .ensure_dir(file.path(summary_dir, "supporting_tables"))
  figure_root <- .ensure_dir(file.path(summary_dir, "supporting_figures"))
  manifest <- list()
  idx <- 1L
  copy_one <- function(src, dest_dir, category, method = "all_methods") {
    .ensure_dir(dest_dir)
    dest <- file.path(dest_dir, basename(src))
    ok <- file.exists(src) && file.copy(src, dest, overwrite = TRUE)
    manifest[[idx]] <<- data.frame(
      category = category,
      method = method,
      source = src,
      destination = dest,
      copied = isTRUE(ok),
      stringsAsFactors = FALSE
    )
    idx <<- idx + 1L
  }

  all_method_tables <- file.path(
    output_root,
    "09_tn_ploidy_cluster_response",
    c(
      "tn_ploidy_cluster_key_conclusions_all_methods.csv",
      "tumor_Q1_ploidy_within_dose_sample_tests_all_methods.csv",
      "tumor_Q2_dose_within_ploidy_sample_tests_all_methods.csv",
      "tumor_Q2_dose_pairwise_within_ploidy_sample_tests_all_methods.csv",
      "tumor_Q3_ploidy_dose_feature_summary_all_methods.csv"
    )
  )
  for (src in all_method_tables) {
    copy_one(src, .ensure_dir(file.path(table_root, "all_methods")), "table", "all_methods")
  }

  method_specs <- data.frame(
    method = c(method_dir_scvelo, method_dir_monocle3),
    prefix = c("scVelo", "Monocle3"),
    stringsAsFactors = FALSE
  )
  table_suffixes <- c(
    "_tumor_Q0_ploidy_dose_cluster_stratum_summary.csv",
    "_tumor_Q1_ploidy_within_dose_sample_tests.csv",
    "_tumor_Q1_ploidy_within_dose_cell_tests.csv",
    "_tumor_Q2_dose_within_ploidy_sample_tests.csv",
    "_tumor_Q2_dose_pairwise_within_ploidy_sample_tests.csv",
    "_tumor_Q3_ploidy_dose_feature_summary.csv",
    "_tn_cluster_2N_vs_4N_comparisons.csv",
    "_tumor_cluster_dose_delta_from_0.csv",
    "_tumor_cluster_dose_ploidy_interaction.csv"
  )
  figure_stubs <- c(
    "_tumor_Q1_ploidy_within_dose_heatmap",
    "_tumor_Q2_dose_within_ploidy_heatmap",
    "_tumor_Q2_pairwise_dose_within_ploidy_heatmap",
    "_tumor_Q3_ploidy_dose_model_fdr_heatmap",
    "_tumor_Q3_primary_timing_feature",
    "_tumor_cluster_ploidy_dose_distribution",
    "_tumor_cluster_4N_minus_2N_by_dose_heatmap",
    "_tumor_dose_delta_from_0_heatmap",
    "_tumor_dose_ploidy_interaction_heatmap"
  )
  for (i in seq_len(nrow(method_specs))) {
    method <- method_specs$method[i]
    prefix <- method_specs$prefix[i]
    method_table_dir <- .ensure_dir(file.path(table_root, method))
    method_figure_dir <- .ensure_dir(file.path(figure_root, method))
    source_method_dir <- file.path(output_root, "09_tn_ploidy_cluster_response", method)
    for (suffix in table_suffixes) {
      copy_one(file.path(source_method_dir, paste0(prefix, suffix)), method_table_dir, "table", method)
    }
    for (stub in figure_stubs) {
      for (ext in c(".pdf", ".png")) {
        copy_one(file.path(source_method_dir, paste0(prefix, stub, ext)), method_figure_dir, "figure", method)
      }
    }
    hist_dir <- file.path(output_root, "10_cluster_timing", method, "ploidy_pseudotime_histograms", "tumor_sample_pages")
    for (ext in c(".pdf", ".png")) {
      copy_one(file.path(hist_dir, paste0("all_samples_ploidy_pseudotime_frequency_histograms_16x1", ext)), method_figure_dir, "figure", method)
    }
  }
  dplyr::bind_rows(manifest)
}

write_final_tumor_summary_outputs <- function(output_root, summary_dir, scvelo_result, monocle3_result) {
  .ensure_dir(summary_dir)
  method_overview <- dplyr::bind_rows(
    make_tumor_method_conclusion_overview(method_scvelo, scvelo_result),
    make_tumor_method_conclusion_overview(method_monocle3, monocle3_result)
  )
  q3_all <- dplyr::bind_rows(
    scvelo_result$tumor_ploidy_dose$feature_summary,
    monocle3_result$tumor_ploidy_dose$feature_summary
  )
  write_table_csv(method_overview, file.path(summary_dir, "tumor_ploidy_dose_method_overview.csv"))
  write_table_csv(q3_all, file.path(summary_dir, "tumor_ploidy_dose_feature_summary_all_methods.csv"))
  writeLines(make_tumor_conclusion_markdown(method_overview, q3_all), file.path(summary_dir, "tumor_ploidy_dose_summary.md"))
  writeLines(make_tumor_conclusion_markdown_zh(method_overview, q3_all), file.path(summary_dir, "tumor_ploidy_dose_summary_zh.md"), useBytes = TRUE)
  copy_manifest <- copy_summary_supporting_files(output_root, summary_dir)
  write_table_csv(copy_manifest, file.path(summary_dir, "supporting_file_manifest.csv"))
  list(method_overview = method_overview, q3_all = q3_all, copy_manifest = copy_manifest)
}

read_gsea_inputs <- function(root_dir) {
  files <- c(
    cluster = file.path(root_dir, "00_summary", "cluster_Hallmark_GSEA_all.csv"),
    stratified = file.path(root_dir, "00_summary", "stratified_Hallmark_GSEA_all.csv")
  )
  out <- list()
  idx <- 1L
  for (nm in names(files)) {
    file <- files[[nm]]
    if (!file.exists(file) || file.info(file)$size <= 0) next
    df <- readr::read_csv(
      file,
      col_types = readr::cols(.default = readr::col_character()),
      show_col_types = FALSE,
      progress = FALSE
    )
    df <- as.data.frame(df, stringsAsFactors = FALSE)
    df <- ensure_columns(df, c(
      "scope", "cluster", "subset_group", "dose", "comparison", "group_col",
      "ident_1", "ident_2", "plot_label", "marker_file", "pathway", "pval",
      "padj", "ES", "NES", "size", "leadingEdge", "hallmark_label", "direction"
    ))
    df$gsea_source <- nm
    out[[idx]] <- df
    idx <- idx + 1L
  }
  if (length(out) == 0) return(data.frame())
  dplyr::bind_rows(out) %>%
    dplyr::mutate(
      scope = as_clean_chr(.data$scope),
      cluster = as_clean_chr(.data$cluster),
      subset_group = as_clean_chr(.data$subset_group),
      dose = as_clean_chr(.data$dose),
      comparison = as_clean_chr(.data$comparison),
      pathway = as_clean_chr(.data$pathway),
      hallmark_label = as_clean_chr(.data$hallmark_label),
      marker_file = as_clean_chr(.data$marker_file),
      marker_basename = basename(.data$marker_file),
      pval = safe_numeric(.data$pval),
      padj = safe_numeric(.data$padj),
      ES = safe_numeric(.data$ES),
      NES = safe_numeric(.data$NES),
      size = safe_numeric(.data$size)
    )
}

classify_pathway <- function(pathway, label) {
  x <- toupper(coalesce_chr(pathway, label))
  dplyr::case_when(
    grepl("HYPOXIA|HIF", x) ~ "Hypoxia",
    grepl("GLYCOLYSIS|OXIDATIVE|FATTY|ADIPOGENESIS|CHOLESTEROL|HEME|BILE|XENOBIOTIC", x) ~ "Metabolism",
    grepl("E2F|G2M|MITOTIC|MYC|DNA_REPAIR|TELOMERE", x) ~ "Proliferation",
    grepl("TNFA|INTERFERON|INFLAMMATORY|COMPLEMENT|IL6|IL2|ALLOGRAFT", x) ~ "Immune_inflammation",
    grepl("EPITHELIAL_MESENCHYMAL|ANGIOGENESIS|COAGULATION|KRAS|APICAL|MYOGENESIS", x) ~ "State_transition",
    grepl("P53|APOPTOSIS|UV|UNFOLDED|REACTIVE_OXYGEN|PEROXISOME", x) ~ "Stress_death",
    TRUE ~ "Other"
  )
}

read_deg_marker_summaries <- function(root_dir) {
  files <- list.files(
    root_dir,
    pattern = "_markers\\.csv$",
    full.names = TRUE,
    recursive = TRUE
  )
  files <- files[file.info(files)$size > 0]
  if (length(files) == 0) return(data.frame())

  keep_cols <- c(
    "scope", "cluster", "subset_group", "dose", "comparison", "group_col",
    "ident_1", "ident_2", "gene", "gene_symbol", "p_val", "avg_log2FC",
    "avg_logFC", "pct.1", "pct.2", "p_val_adj"
  )
  out <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- read_csv_selected(files[i], keep_cols)
    lfc <- if ("avg_log2FC" %in% names(df) && any(!is.na(df$avg_log2FC))) df$avg_log2FC else df$avg_logFC
    df$lfc <- safe_numeric(lfc)
    df$p_val_adj_num <- safe_numeric(df$p_val_adj)
    df$gene_label <- coalesce_chr(df$gene_symbol, df$gene)
    ranked <- df %>%
      dplyr::filter(!is.na(.data$gene_label), is.finite(.data$lfc))
    sig <- ranked %>%
      dplyr::filter(is.finite(.data$p_val_adj_num), .data$p_val_adj_num < 0.05)
    rank_source <- if (nrow(sig) > 0) sig else ranked
    top_up <- rank_source %>%
      dplyr::filter(.data$lfc > 0) %>%
      dplyr::arrange(dplyr::desc(.data$lfc)) %>%
      dplyr::slice_head(n = 15)
    top_down <- rank_source %>%
      dplyr::filter(.data$lfc < 0) %>%
      dplyr::arrange(.data$lfc) %>%
      dplyr::slice_head(n = 15)
    first_row <- df[1, , drop = FALSE]
    out[[i]] <- data.frame(
      marker_file_local = files[i],
      marker_basename = basename(files[i]),
      scope = as_clean_chr(first_row$scope),
      cluster = as_clean_chr(first_row$cluster),
      subset_group = as_clean_chr(first_row$subset_group),
      dose = as_clean_chr(first_row$dose),
      comparison = as_clean_chr(first_row$comparison),
      group_col = as_clean_chr(first_row$group_col),
      ident_1 = as_clean_chr(first_row$ident_1),
      ident_2 = as_clean_chr(first_row$ident_2),
      n_genes_tested = nrow(ranked),
      n_deg_fdr_0_05 = nrow(sig),
      n_up_deg_fdr_0_05 = sum(sig$lfc > 0, na.rm = TRUE),
      n_down_deg_fdr_0_05 = sum(sig$lfc < 0, na.rm = TRUE),
      top_up_genes = paste(top_up$gene_label, collapse = ";"),
      top_down_genes = paste(top_down$gene_label, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(out)
}

join_trajectory_with_gsea <- function(stat_df, gsea_df, deg_marker_summary) {
  if (nrow(stat_df) == 0 || nrow(gsea_df) == 0) return(data.frame())
  stats_key <- stat_df %>%
    dplyr::filter(
      .data$analysis_family == "03a_DEG_comparison_aligned",
      .data$analysis_level == "cell",
      .data$metric %in% c("velocity_pseudotime", "pseudotime")
    ) %>%
    dplyr::mutate(
      scope = as_clean_chr(.data$TrajectoryScope),
      comparison = as_clean_chr(.data$TrajectoryComparison),
      cluster = as_clean_chr(.data$TrajectoryCluster),
      subset_group = as_clean_chr(.data$TrajectorySubsetGroup),
      dose = as_clean_chr(.data$TrajectoryDose),
      group_1 = as_clean_chr(.data$group_1),
      group_2 = as_clean_chr(.data$group_2)
    )
  if (nrow(stats_key) == 0) return(data.frame())

  gsea_key <- gsea_df %>%
    dplyr::mutate(
      pathway_category = classify_pathway(.data$pathway, .data$hallmark_label),
      hypoxia_pathway = .data$pathway_category == "Hypoxia"
    )

  joined <- dplyr::left_join(
    stats_key,
    gsea_key,
    by = c("scope", "comparison", "cluster", "subset_group", "dose"),
    relationship = "many-to-many"
  )

  if (nrow(deg_marker_summary) > 0 && "marker_basename" %in% names(joined)) {
    joined <- joined %>%
      dplyr::left_join(
        deg_marker_summary %>%
          dplyr::select(
            marker_basename, marker_file_local, n_genes_tested,
            n_deg_fdr_0_05, n_up_deg_fdr_0_05, n_down_deg_fdr_0_05,
            top_up_genes, top_down_genes
          ) %>%
          dplyr::distinct(.data$marker_basename, .keep_all = TRUE),
        by = "marker_basename"
      )
  }

  joined %>%
    dplyr::mutate(
      gsea_significant_fdr_0_05 = is.finite(.data$padj) & .data$padj < 0.05,
      trajectory_significant_fdr_0_05 = is.finite(.data$primary_p_adj_family) & .data$primary_p_adj_family < 0.05,
      enriched_group = dplyr::case_when(
        is.finite(.data$NES) & .data$NES > 0 ~ .data$group_1,
        is.finite(.data$NES) & .data$NES < 0 ~ .data$group_2,
        TRUE ~ NA_character_
      ),
      enriched_time_delta = dplyr::case_when(
        is.finite(.data$NES) & .data$NES > 0 ~ .data$median_diff_group_1_minus_group_2,
        is.finite(.data$NES) & .data$NES < 0 ~ -.data$median_diff_group_1_minus_group_2,
        TRUE ~ NA_real_
      ),
      pathway_timing = dplyr::case_when(
        is.finite(.data$enriched_time_delta) & .data$enriched_time_delta > 0 ~ "late_or_high_time",
        is.finite(.data$enriched_time_delta) & .data$enriched_time_delta < 0 ~ "early_or_low_time",
        is.finite(.data$enriched_time_delta) ~ "neutral",
        TRUE ~ NA_character_
      ),
      interpretation_strength = dplyr::case_when(
        .data$gsea_significant_fdr_0_05 & .data$trajectory_significant_fdr_0_05 ~ "GSEA_and_pseudotime_FDR_0.05",
        .data$gsea_significant_fdr_0_05 ~ "GSEA_FDR_0.05_only",
        .data$trajectory_significant_fdr_0_05 ~ "pseudotime_FDR_0.05_only",
        TRUE ~ "descriptive"
      )
    )
}

summarize_pathway_timing <- function(joined_df) {
  if (nrow(joined_df) == 0 || !("pathway" %in% names(joined_df))) return(data.frame())
  joined_df %>%
    dplyr::filter(!is.na(.data$pathway), is.finite(.data$enriched_time_delta)) %>%
    dplyr::group_by(.data$method, .data$metric, .data$pathway_category, .data$pathway, .data$hallmark_label) %>%
    dplyr::summarise(
      n_comparisons = dplyr::n(),
      n_gsea_fdr_0_05 = sum(.data$gsea_significant_fdr_0_05, na.rm = TRUE),
      n_trajectory_fdr_0_05 = sum(.data$trajectory_significant_fdr_0_05, na.rm = TRUE),
      n_joint_fdr_0_05 = sum(.data$gsea_significant_fdr_0_05 & .data$trajectory_significant_fdr_0_05, na.rm = TRUE),
      n_late_or_high_time = sum(.data$pathway_timing == "late_or_high_time", na.rm = TRUE),
      n_early_or_low_time = sum(.data$pathway_timing == "early_or_low_time", na.rm = TRUE),
      median_enriched_time_delta = stats::median(.data$enriched_time_delta, na.rm = TRUE),
      median_abs_NES = stats::median(abs(.data$NES), na.rm = TRUE),
      min_gsea_fdr = min(.data$padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      dominant_timing = dplyr::case_when(
        .data$n_late_or_high_time > .data$n_early_or_low_time ~ "mostly_late_or_high_time",
        .data$n_early_or_low_time > .data$n_late_or_high_time ~ "mostly_early_or_low_time",
        TRUE ~ "mixed_or_neutral"
      )
    ) %>%
    dplyr::arrange(.data$method, .data$metric, dplyr::desc(.data$n_joint_fdr_0_05), .data$min_gsea_fdr)
}

plot_gsea_trajectory_association <- function(joined_df) {
  needed <- c(
    "pathway", "NES", "median_diff_group_1_minus_group_2",
    "gsea_significant_fdr_0_05", "hypoxia_pathway", "trajectory_significant_fdr_0_05",
    "pathway_category", "method", "metric"
  )
  if (nrow(joined_df) == 0 || !all(needed %in% names(joined_df))) return(invisible(NULL))
  plot_df <- joined_df %>%
    dplyr::filter(!is.na(.data$pathway), is.finite(.data$NES), is.finite(.data$median_diff_group_1_minus_group_2)) %>%
    dplyr::filter(.data$gsea_significant_fdr_0_05 | .data$hypoxia_pathway | .data$trajectory_significant_fdr_0_05)
  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = median_diff_group_1_minus_group_2, y = NES)) +
    geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_point(aes(color = pathway_category, shape = hypoxia_pathway), alpha = 0.72, size = 2) +
    facet_grid(method ~ metric, scales = "free") +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), name = "Hypoxia") +
    labs(
      title = "GSEA enrichment versus pseudotime effect",
      subtitle = "X axis is group_1 minus group_2 median time; NES sign follows the 03a DEG contrast.",
      x = "Pseudotime median difference",
      y = "Hallmark NES",
      color = "Pathway class"
    ) +
    theme_extra(9)
  save_both(p, file.path(out_biology, "gsea_nes_vs_pseudotime_effect"), width = 11, height = 7.5)
}

plot_hypoxia_timing <- function(hypoxia_df) {
  needed <- c("enriched_time_delta", "TrajectoryComparison", "TrajectoryCluster", "TrajectoryDose", "method", "metric", "TrajectoryScope")
  if (nrow(hypoxia_df) == 0 || !all(needed %in% names(hypoxia_df))) return(invisible(NULL))
  plot_df <- hypoxia_df %>%
    dplyr::filter(is.finite(.data$enriched_time_delta)) %>%
    dplyr::mutate(
      comparison_label = paste0(
        .data$TrajectoryComparison,
        ifelse(!is.na(.data$TrajectoryCluster), paste0(" | cluster ", .data$TrajectoryCluster), ""),
        ifelse(!is.na(.data$TrajectoryDose), paste0(" | ", .data$TrajectoryDose), "")
      )
    ) %>%
    dplyr::group_by(.data$method, .data$metric) %>%
    dplyr::arrange(dplyr::desc(abs(.data$enriched_time_delta)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 80) %>%
    dplyr::ungroup()
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df <- plot_df %>%
    dplyr::mutate(
      comparison_label = factor(.data$comparison_label, levels = unique(.data$comparison_label[order(.data$enriched_time_delta)]))
    )
  p <- ggplot(plot_df, aes(x = comparison_label, y = method, fill = enriched_time_delta)) +
    geom_tile(color = "white", linewidth = 0.25, width = 0.95, height = 0.95) +
    facet_grid(metric ~ TrajectoryScope, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Hypoxia-enriched\ngroup time delta"
    ) +
    labs(
      title = "HALLMARK_HYPOXIA timing inferred from pseudotime and GSEA",
      subtitle = "Positive means hypoxia-enriched group is later or higher along the time metric; negative means earlier or lower.",
      x = NULL,
      y = NULL
    ) +
    theme_extra(7) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1))
  save_both(p, file.path(out_biology, "hypoxia_timing_by_method_scope"), width = 18, height = 8)
}

plot_pathway_timing_summary <- function(pathway_summary) {
  needed <- c(
    "median_enriched_time_delta", "n_gsea_fdr_0_05", "pathway_category",
    "method", "metric", "pathway", "n_joint_fdr_0_05"
  )
  if (nrow(pathway_summary) == 0 || !all(needed %in% names(pathway_summary))) return(invisible(NULL))
  plot_df <- pathway_summary %>%
    dplyr::filter(is.finite(.data$median_enriched_time_delta)) %>%
    dplyr::filter(.data$n_gsea_fdr_0_05 > 0 | .data$pathway_category %in% c("Hypoxia", "Proliferation", "State_transition", "Metabolism")) %>%
    dplyr::group_by(.data$method, .data$metric) %>%
    dplyr::arrange(dplyr::desc(.data$n_joint_fdr_0_05), dplyr::desc(abs(.data$median_enriched_time_delta)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 30) %>%
    dplyr::ungroup()
  if (nrow(plot_df) == 0) return(invisible(NULL))

  plot_df <- plot_df %>%
    dplyr::mutate(pathway_label = factor(.data$pathway, levels = unique(.data$pathway[order(.data$pathway_category, .data$median_enriched_time_delta)])))
  p <- ggplot(plot_df, aes(x = method, y = pathway_label, fill = median_enriched_time_delta)) +
    geom_tile(color = "white", linewidth = 0.25, width = 0.95, height = 0.95) +
    geom_point(aes(size = n_joint_fdr_0_05), shape = 21, color = "black", fill = "white", alpha = 0.75) +
    facet_grid(pathway_category ~ metric, scales = "free_y", space = "free_y") +
    scale_fill_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      name = "Median pathway\ntime delta"
    ) +
    scale_size_continuous(name = "Joint FDR<0.05", range = c(1.2, 4)) +
    labs(
      title = "Pathway timing summary from DEG/GSEA and pseudotime effects",
      x = NULL,
      y = NULL
    ) +
    theme_extra(8) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  save_both(p, file.path(out_biology, "pathway_timing_summary_by_method"), width = 10, height = 12)
}

plot_cluster_hypoxia_umap <- function(desc_df, gsea_df, method_label, metric_col, out_dir, plot_prefix, out_umap_dir) {
  hypoxia_cluster <- gsea_df %>%
    dplyr::filter(.data$scope == "clusters_vs_rest", .data$pathway == "HALLMARK_HYPOXIA") %>%
    dplyr::select(cluster = cluster, hypoxia_NES = NES, hypoxia_fdr = padj) %>%
    dplyr::distinct(.data$cluster, .keep_all = TRUE)
  if (nrow(hypoxia_cluster) == 0) return(invisible(NULL))
  plot_df <- desc_df %>%
    dplyr::filter(is.finite(.data$UMAP_1), is.finite(.data$UMAP_2), !is.na(.data$clusters)) %>%
    dplyr::left_join(hypoxia_cluster, by = c("clusters" = "cluster")) %>%
    dplyr::filter(is.finite(.data$hypoxia_NES))
  if (nrow(plot_df) == 0) return(invisible(NULL))
  if (nrow(plot_df) > 60000) {
    set.seed(1234)
    plot_df <- plot_df[sort(sample(seq_len(nrow(plot_df)), 60000)), , drop = FALSE]
  }

  centers <- plot_df %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(
      UMAP_1 = stats::median(.data$UMAP_1, na.rm = TRUE),
      UMAP_2 = stats::median(.data$UMAP_2, na.rm = TRUE),
      .groups = "drop"
    )
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = hypoxia_NES)) +
    geom_point(size = 0.22, alpha = 0.7) +
    geom_text(
      data = centers,
      aes(x = UMAP_1, y = UMAP_2, label = clusters),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold",
      color = "black"
    ) +
	    scale_color_gradient2(
	      low = "#2C7BB6",
	      mid = "white",
	      high = "#D7191C",
	      midpoint = 0,
	      name = "Cluster hypoxia NES"
	    ) +
	    coord_square_umap(plot_df$UMAP_1, plot_df$UMAP_2) +
	    labs(
      title = paste0(method_label, " UMAP overlay with cluster HALLMARK_HYPOXIA NES"),
      subtitle = paste0("UMAP cells are from the ", metric_col, " descriptive pseudotime dataset"),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()
	  save_square_panel_both(p, file.path(out_dir, paste0(plot_prefix, "_umap_cluster_hypoxia_NES")), panel_size = 5.2)
	  save_square_panel_both(p, file.path(out_umap_dir, paste0(plot_prefix, "_umap_cluster_hypoxia_NES")), panel_size = 5.2)
}

format_cluster_timing_list <- function(df, n = 5L) {
  if (nrow(df) == 0) return(NA_character_)
  df <- df %>% dplyr::slice_head(n = n)
  label <- paste0(df$clusters, ":", signif(df$median_time, 3))
  if ("annotation_primary" %in% names(df)) {
    annotation <- as_clean_chr(df$annotation_primary)
    label <- ifelse(!is.na(annotation), paste0(label, "(", annotation, ")"), label)
  }
  paste(label, collapse = ";")
}

make_method_answer_summary <- function(method_result, method_label) {
  commonality <- method_result$commonality
  manifest <- method_result$manifest
  if (nrow(commonality) == 0) {
    return(data.frame(
      method = method_label,
      topic = "pseudotime_features",
      summary = "No cluster-level descriptive pseudotime table was available.",
      stringsAsFactors = FALSE
    ))
  }

  early <- commonality %>% dplyr::arrange(.data$median_time)
  late <- commonality %>% dplyr::arrange(dplyr::desc(.data$median_time))
  ploidy_counts <- commonality %>%
    dplyr::count(.data$ploidy_pattern, name = "n") %>%
    dplyr::mutate(label = paste0(.data$ploidy_pattern, "=", .data$n))
  dose_counts <- commonality %>%
    dplyr::count(.data$dose_commonality, name = "n") %>%
    dplyr::mutate(label = paste0(.data$dose_commonality, "=", .data$n))

  data.frame(
    method = method_label,
    topic = c(
      "descriptive_dataset",
      "early_clusters",
      "late_clusters",
      "ploidy_pattern_across_clusters",
      "dose_pattern_commonality"
    ),
    summary = c(
      paste0(
        "source=", manifest$descriptor_source[1],
        "; n_cells=", manifest$n_cells[1],
        "; n_clusters=", manifest$n_clusters[1]
      ),
      format_cluster_timing_list(early, 5L),
      format_cluster_timing_list(late, 5L),
      paste(ploidy_counts$label, collapse = ";"),
      paste(dose_counts$label, collapse = ";")
    ),
    stringsAsFactors = FALSE
  )
}

get_config_scalar <- function(name, default = "") {
  value <- config[[name]]
  if (!is.null(value) && length(value) > 0 && nzchar(trimws(as.character(value[1])))) {
    return(as.character(value[1]))
  }
  default
}

resolve_scvelo_python_for_extra <- function() {
  env_name <- Sys.getenv("SCVELO_PYTHON_ENV", unset = get_config_scalar("Scvelo_python_env", "rna_velocity_py310"))
  configured <- Sys.getenv("SCVELO_PYTHON", unset = get_config_scalar("Scvelo_python", ""))
  pyenv_root <- Sys.getenv("PYENV_ROOT", unset = file.path(Sys.getenv("HOME", unset = ""), ".pyenv"))
  candidates <- c(
    configured,
    file.path(pyenv_root, "versions", env_name, "bin", "python"),
    file.path(pyenv_root, "versions", env_name, "bin", "python3"),
    Sys.which("python3"),
    Sys.which("python")
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0) return(NA_character_)
  normalizePath(hit[1], mustWork = TRUE)
}

build_scvelo_flow_manifest <- function(root_dir) {
  h5ads <- list.files(root_dir, pattern = "^scvelo_result\\.h5ad$", full.names = TRUE, recursive = TRUE)
  h5ads <- h5ads[file.exists(h5ads) & file.info(h5ads)$size > 0]
  if (length(h5ads) == 0) return(data.frame())

  root_norm <- normalizePath(root_dir, mustWork = TRUE)
  rows <- lapply(sort(normalizePath(h5ads, mustWork = TRUE)), function(h5ad) {
    rel_file <- substring(h5ad, nchar(root_norm) + 2L)
    rel_dir <- dirname(rel_file)
    rel_parts <- strsplit(rel_dir, .Platform$file.sep, fixed = TRUE)[[1]]
    rel_parts <- rel_parts[nzchar(rel_parts)]
    top_level <- if (length(rel_parts) >= 1L) rel_parts[1] else "scVelo"
    analysis_scope <- if (length(rel_parts) >= 2L) rel_parts[2] else top_level
    task_name <- if (length(rel_parts) >= 1L) rel_parts[length(rel_parts)] else "scVelo"
    label_parts <- rel_parts
    if (length(label_parts) > 3L) label_parts <- c(label_parts[1:2], label_parts[length(label_parts)])
    label <- gsub("_", " ", paste(label_parts, collapse = " / "), fixed = TRUE)
    if (nchar(label) > 150L) label <- paste0(substr(label, 1L, 147L), "...")
    data.frame(
      h5ad = h5ad,
      source_relative_dir = rel_dir,
      source_top_level = top_level,
      source_analysis_scope = analysis_scope,
      source_task_name = task_name,
      output_subdir = rel_dir,
      output_prefix = "scVelo",
      label = label,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(rows)
}

read_h5ad_obs_vector <- function(h5_file, obs_name, n_obs) {
  if (!(obs_name %in% names(h5_file[["obs"]]))) return(rep(NA_character_, n_obs))
  obj <- h5_file[["obs"]][[obs_name]]
  out <- rep(NA_character_, n_obs)
  if (inherits(obj, "H5Group") && all(c("codes", "categories") %in% names(obj))) {
    codes <- as.integer(obj[["codes"]][])
    categories <- as.character(obj[["categories"]][])
    ok <- !is.na(codes) & codes >= 0 & codes < length(categories)
    out[seq_along(codes)[ok]] <- categories[codes[ok] + 1L]
  } else {
    values <- as.character(obj[])
    out[seq_along(values)] <- values
  }
  out
}

read_scvelo_h5ad_flow <- function(h5ad_file, max_points = 60000L, max_vector_cells = 50000L) {
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop("Package hdf5r is required to read scvelo_result.h5ad.", call. = FALSE)
  }
  f <- hdf5r::H5File$new(h5ad_file, "r")
  on.exit(f$close_all(), add = TRUE)
  if (!("obsm" %in% names(f)) || !("uns" %in% names(f))) {
    stop("Missing obsm or uns in h5ad: ", h5ad_file, call. = FALSE)
  }
  coords_raw <- f[["obsm/X_umap"]][,]
  n_obs <- length(f[["obs/_index"]][])
  coords <- if (nrow(coords_raw) == 2 && ncol(coords_raw) == n_obs) {
    t(coords_raw)
  } else {
    coords_raw
  }
  if (ncol(coords) != 2) stop("X_umap is not two-dimensional in: ", h5ad_file, call. = FALSE)

  time <- safe_numeric(f[["obs/velocity_pseudotime"]][])
  cells <- as.character(f[["obs/_index"]][])
  clusters <- read_h5ad_obs_vector(f, "clusters", n_obs)
  if (all(is.na(clusters))) clusters <- read_h5ad_obs_vector(f, "seurat_clusters", n_obs)
  point_df <- data.frame(
    cell = cells,
    UMAP_1 = coords[, 1],
    UMAP_2 = coords[, 2],
    velocity_pseudotime = time,
    clusters = clusters,
    stringsAsFactors = FALSE
  )
  point_df <- point_df[
    is.finite(point_df$UMAP_1) &
      is.finite(point_df$UMAP_2) &
      is.finite(point_df$velocity_pseudotime),
    ,
    drop = FALSE
  ]
  if (nrow(point_df) > max_points) {
    set.seed(1234)
    point_df <- point_df[sort(sample(seq_len(nrow(point_df)), max_points)), , drop = FALSE]
  }

  if (!("velocity_graph" %in% names(f[["uns"]]))) {
    stop("Missing uns/velocity_graph in h5ad: ", h5ad_file, call. = FALSE)
  }
  indptr <- as.integer(f[["uns/velocity_graph/indptr"]][])
  indices <- as.integer(f[["uns/velocity_graph/indices"]][]) + 1L
  weights <- safe_numeric(f[["uns/velocity_graph/data"]][])
  n_graph <- length(indptr) - 1L
  vector_candidates <- which(
    seq_along(time) <= n_graph &
      is.finite(time) &
      is.finite(coords[, 1]) &
      is.finite(coords[, 2])
  )
  if (length(vector_candidates) > max_vector_cells) {
    set.seed(1234)
    vector_candidates <- sort(sample(vector_candidates, max_vector_cells))
  }

  vector_rows <- vector("list", length(vector_candidates))
  out_idx <- 1L
  for (i in vector_candidates) {
    start <- indptr[i] + 1L
    end <- indptr[i + 1L]
    if (end < start) next
    nbr <- indices[start:end]
    w <- weights[start:end]
    ok <- is.finite(w) &
      w > 0 &
      nbr >= 1L &
      nbr <= nrow(coords) &
      is.finite(coords[nbr, 1]) &
      is.finite(coords[nbr, 2])
    if (!any(ok)) next
    nbr <- nbr[ok]
    w <- w[ok]
    target <- c(
      stats::weighted.mean(coords[nbr, 1], w, na.rm = TRUE),
      stats::weighted.mean(coords[nbr, 2], w, na.rm = TRUE)
    )
    dx <- target[1] - coords[i, 1]
    dy <- target[2] - coords[i, 2]
    len <- sqrt(dx^2 + dy^2)
    if (!is.finite(len) || len <= 0) next
    vector_rows[[out_idx]] <- data.frame(
      UMAP_1 = coords[i, 1],
      UMAP_2 = coords[i, 2],
      velocity_dx = dx,
      velocity_dy = dy,
      velocity_pseudotime = time[i],
      clusters = clusters[i],
      vector_length = len,
      stringsAsFactors = FALSE
    )
    out_idx <- out_idx + 1L
  }
  vector_df <- dplyr::bind_rows(vector_rows)
  list(points = point_df, vectors = vector_df, n_obs = n_obs)
}

make_scvelo_cluster_layers <- function(point_df, nx = 120L, ny = 120L, core_fraction = 0.90) {
  point_df <- point_df %>%
    dplyr::mutate(
      clusters = as_clean_chr(.data$clusters),
      velocity_pseudotime = safe_numeric(.data$velocity_pseudotime)
    ) %>%
    dplyr::filter(
      !is.na(.data$clusters),
      is.finite(.data$UMAP_1),
      is.finite(.data$UMAP_2),
      is.finite(.data$velocity_pseudotime)
    )
  if (nrow(point_df) == 0) {
    return(list(labels = data.frame(), hulls = data.frame(), boundaries = data.frame(), tiles = data.frame()))
  }

  x_range <- range(point_df$UMAP_1, finite = TRUE)
  y_range <- range(point_df$UMAP_2, finite = TRUE)
  if (diff(x_range) <= 0) x_range <- x_range + c(-0.5, 0.5)
  if (diff(y_range) <= 0) y_range <- y_range + c(-0.5, 0.5)
  x_pad <- diff(x_range) * 0.02
  y_pad <- diff(y_range) * 0.02
  x_range <- x_range + c(-x_pad, x_pad)
  y_range <- y_range + c(-y_pad, y_pad)
  x_breaks <- seq(x_range[1], x_range[2], length.out = nx + 1L)
  y_breaks <- seq(y_range[1], y_range[2], length.out = ny + 1L)
  tile_width <- diff(x_range) / nx
  tile_height <- diff(y_range) / ny

  point_df$grid_x <- pmin(pmax(findInterval(point_df$UMAP_1, x_breaks, all.inside = TRUE), 1L), nx)
  point_df$grid_y <- pmin(pmax(findInterval(point_df$UMAP_2, y_breaks, all.inside = TRUE), 1L), ny)
  cluster_grid <- point_df %>%
    dplyr::group_by(.data$grid_x, .data$grid_y, .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_time = mean(.data$velocity_pseudotime, na.rm = TRUE),
      .groups = "drop"
    )
  core_tiles <- cluster_grid %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::arrange(dplyr::desc(.data$n_cells), .data$grid_x, .data$grid_y, .by_group = TRUE) %>%
    dplyr::mutate(
      total_cluster_cells = sum(.data$n_cells, na.rm = TRUE),
      cumulative_cells = cumsum(.data$n_cells),
      previous_cumulative_cells = dplyr::lag(.data$cumulative_cells, default = 0L),
      core_fraction_covered = .data$cumulative_cells / .data$total_cluster_cells,
      keep_core_tile = .data$previous_cumulative_cells < core_fraction * .data$total_cluster_cells
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$keep_core_tile)
  tiles <- core_tiles %>%
    dplyr::group_by(.data$grid_x, .data$grid_y) %>%
    dplyr::slice_max(.data$n_cells, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      UMAP_1 = (x_breaks[.data$grid_x] + x_breaks[.data$grid_x + 1L]) / 2,
      UMAP_2 = (y_breaks[.data$grid_y] + y_breaks[.data$grid_y + 1L]) / 2
    )
  attr(tiles, "tile_width") <- tile_width
  attr(tiles, "tile_height") <- tile_height
  attr(tiles, "core_fraction") <- core_fraction

  labels <- tiles %>%
    dplyr::group_by(.data$clusters) %>%
    dplyr::summarise(
      UMAP_1 = stats::weighted.mean(.data$UMAP_1, .data$n_cells, na.rm = TRUE),
      UMAP_2 = stats::weighted.mean(.data$UMAP_2, .data$n_cells, na.rm = TRUE),
      median_time = stats::weighted.mean(.data$mean_time, .data$n_cells, na.rm = TRUE),
      n_cells = sum(.data$n_cells, na.rm = TRUE),
      .groups = "drop"
    )
  cluster_order <- sort_maybe_numeric(labels$clusters)
  labels <- labels %>%
    dplyr::mutate(clusters = factor(.data$clusters, levels = cluster_order)) %>%
    dplyr::arrange(.data$clusters) %>%
    dplyr::mutate(clusters = as.character(.data$clusters))

  hull_rows <- lapply(split(tiles, tiles$clusters), function(df) {
    pts <- unique(df[, c("UMAP_1", "UMAP_2"), drop = FALSE])
    if (nrow(pts) < 3) return(data.frame())
    idx <- chull(pts$UMAP_1, pts$UMAP_2)
    idx <- c(idx, idx[1])
    data.frame(
      clusters = df$clusters[1],
      UMAP_1 = pts$UMAP_1[idx],
      UMAP_2 = pts$UMAP_2[idx],
      stringsAsFactors = FALSE
    )
  })
  hulls <- dplyr::bind_rows(hull_rows)

  tile_boundary_segments <- function(tiles_df) {
    if (nrow(tiles_df) == 0) return(data.frame())
    key <- paste(tiles_df$clusters, tiles_df$grid_x, tiles_df$grid_y, sep = "\r")
    make_edges <- function(mask, side) {
      if (!any(mask)) return(data.frame())
      df <- tiles_df[mask, , drop = FALSE]
      half_w <- tile_width / 2
      half_h <- tile_height / 2
      if (side == "left") {
        x <- df$UMAP_1 - half_w
        y <- df$UMAP_2 - half_h
        xend <- df$UMAP_1 - half_w
        yend <- df$UMAP_2 + half_h
      } else if (side == "right") {
        x <- df$UMAP_1 + half_w
        y <- df$UMAP_2 - half_h
        xend <- df$UMAP_1 + half_w
        yend <- df$UMAP_2 + half_h
      } else if (side == "bottom") {
        x <- df$UMAP_1 - half_w
        y <- df$UMAP_2 - half_h
        xend <- df$UMAP_1 + half_w
        yend <- df$UMAP_2 - half_h
      } else {
        x <- df$UMAP_1 - half_w
        y <- df$UMAP_2 + half_h
        xend <- df$UMAP_1 + half_w
        yend <- df$UMAP_2 + half_h
      }
      data.frame(
        clusters = df$clusters,
        UMAP_1 = x,
        UMAP_2 = y,
        UMAP_1_to = xend,
        UMAP_2_to = yend,
        stringsAsFactors = FALSE
      )
    }
    dplyr::bind_rows(
      make_edges(!(paste(tiles_df$clusters, tiles_df$grid_x - 1L, tiles_df$grid_y, sep = "\r") %in% key), "left"),
      make_edges(!(paste(tiles_df$clusters, tiles_df$grid_x + 1L, tiles_df$grid_y, sep = "\r") %in% key), "right"),
      make_edges(!(paste(tiles_df$clusters, tiles_df$grid_x, tiles_df$grid_y - 1L, sep = "\r") %in% key), "bottom"),
      make_edges(!(paste(tiles_df$clusters, tiles_df$grid_x, tiles_df$grid_y + 1L, sep = "\r") %in% key), "top")
    )
  }
  boundaries <- tile_boundary_segments(tiles)

  list(labels = labels, hulls = hulls, boundaries = boundaries, tiles = tiles)
}

make_scvelo_velocity_grid <- function(vector_df, nx = 32L, ny = 32L, min_bin_n = 5L) {
  if (nrow(vector_df) == 0) return(data.frame())
  vector_df <- vector_df %>%
    dplyr::filter(
      is.finite(.data$UMAP_1),
      is.finite(.data$UMAP_2),
      is.finite(.data$velocity_dx),
      is.finite(.data$velocity_dy),
      is.finite(.data$vector_length),
      .data$vector_length > 0
    )
  if (nrow(vector_df) == 0) return(data.frame())

  len_cap <- stats::quantile(vector_df$vector_length, 0.995, na.rm = TRUE, names = FALSE)
  if (is.finite(len_cap) && len_cap > 0) {
    vector_df <- vector_df %>% dplyr::filter(.data$vector_length <= len_cap)
  }
  if (nrow(vector_df) == 0) return(data.frame())

  x_range <- range(vector_df$UMAP_1, finite = TRUE)
  y_range <- range(vector_df$UMAP_2, finite = TRUE)
  if (!all(is.finite(x_range)) || !all(is.finite(y_range))) return(data.frame())
  if (diff(x_range) <= 0) x_range <- x_range + c(-0.5, 0.5)
  if (diff(y_range) <= 0) y_range <- y_range + c(-0.5, 0.5)
  x_pad <- diff(x_range) * 0.02
  y_pad <- diff(y_range) * 0.02
  x_range <- x_range + c(-x_pad, x_pad)
  y_range <- y_range + c(-y_pad, y_pad)

  x_breaks <- seq(x_range[1], x_range[2], length.out = nx + 1L)
  y_breaks <- seq(y_range[1], y_range[2], length.out = ny + 1L)
  bin_w <- diff(x_range) / nx
  bin_h <- diff(y_range) / ny
  step_size <- min(bin_w, bin_h)
  if (!is.finite(step_size) || step_size <= 0) return(data.frame())

  vector_df$grid_x <- pmin(pmax(findInterval(vector_df$UMAP_1, x_breaks, all.inside = TRUE), 1L), nx)
  vector_df$grid_y <- pmin(pmax(findInterval(vector_df$UMAP_2, y_breaks, all.inside = TRUE), 1L), ny)

  grid_df <- vector_df %>%
    dplyr::group_by(.data$grid_x, .data$grid_y) %>%
    dplyr::summarise(
      UMAP_1 = mean(.data$UMAP_1, na.rm = TRUE),
      UMAP_2 = mean(.data$UMAP_2, na.rm = TRUE),
      velocity_dx = mean(.data$velocity_dx, na.rm = TRUE),
      velocity_dy = mean(.data$velocity_dy, na.rm = TRUE),
      velocity_pseudotime = mean(.data$velocity_pseudotime, na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(vector_length = sqrt(.data$velocity_dx^2 + .data$velocity_dy^2)) %>%
    dplyr::filter(
      .data$n_cells >= min_bin_n,
      is.finite(.data$vector_length),
      .data$vector_length > 0
    )
  if (nrow(grid_df) == 0) return(data.frame())

  mag_ref <- stats::quantile(grid_df$vector_length, 0.90, na.rm = TRUE, names = FALSE)
  if (!is.finite(mag_ref) || mag_ref <= 0) mag_ref <- max(grid_df$vector_length, na.rm = TRUE)
  if (!is.finite(mag_ref) || mag_ref <= 0) return(data.frame())
  arrow_scale <- 0.72 * step_size / mag_ref
  grid_df <- grid_df %>%
    dplyr::mutate(
      UMAP_1_to = .data$UMAP_1 + .data$velocity_dx * arrow_scale,
      UMAP_2_to = .data$UMAP_2 + .data$velocity_dy * arrow_scale,
      scaled_length = .data$vector_length * arrow_scale
    )
  attr(grid_df, "x_range") <- x_range
  attr(grid_df, "y_range") <- y_range
  attr(grid_df, "step_size") <- step_size
  grid_df
}

make_scvelo_velocity_streamlines <- function(grid_df, max_streamlines = 140L, forward_steps = 34L, backward_steps = 10L) {
  if (nrow(grid_df) == 0) return(data.frame())
  step_size <- attr(grid_df, "step_size")
  x_range <- attr(grid_df, "x_range")
  y_range <- attr(grid_df, "y_range")
  if (!is.finite(step_size) || step_size <= 0 || !all(is.finite(x_range)) || !all(is.finite(y_range))) {
    return(data.frame())
  }
  field_df <- grid_df %>%
    dplyr::filter(
      is.finite(.data$UMAP_1),
      is.finite(.data$UMAP_2),
      is.finite(.data$velocity_dx),
      is.finite(.data$velocity_dy),
      is.finite(.data$vector_length),
      .data$vector_length > 0
    )
  if (nrow(field_df) == 0) return(data.frame())

  seed_df <- field_df %>%
    dplyr::filter(.data$vector_length >= stats::quantile(.data$vector_length, 0.35, na.rm = TRUE, names = FALSE))
  if (nrow(seed_df) == 0) seed_df <- field_df
  if (nrow(seed_df) > max_streamlines) {
    set.seed(1234)
    sample_prob <- pmax(seed_df$vector_length, 0)
    if (!any(is.finite(sample_prob)) || sum(sample_prob, na.rm = TRUE) <= 0) sample_prob <- rep(1, nrow(seed_df))
    seed_idx <- sort(sample(seq_len(nrow(seed_df)), max_streamlines, prob = sample_prob))
    seed_df <- seed_df[seed_idx, , drop = FALSE]
  }

  nearest_vector <- function(x, y) {
    d2 <- (field_df$UMAP_1 - x)^2 + (field_df$UMAP_2 - y)^2
    k <- which.min(d2)
    if (length(k) == 0 || !is.finite(d2[k]) || sqrt(d2[k]) > step_size * 2.2) return(NULL)
    vlen <- field_df$vector_length[k]
    if (!is.finite(vlen) || vlen <= 0) return(NULL)
    c(dx = field_df$velocity_dx[k] / vlen, dy = field_df$velocity_dy[k] / vlen)
  }

  trace_direction <- function(x, y, direction, n_steps) {
    pts <- vector("list", n_steps)
    out_idx <- 1L
    for (step in seq_len(n_steps)) {
      unit_v <- nearest_vector(x, y)
      if (is.null(unit_v)) break
      x <- x + direction * unit_v[["dx"]] * step_size * 0.55
      y <- y + direction * unit_v[["dy"]] * step_size * 0.55
      if (!is.finite(x) || !is.finite(y) || x < x_range[1] || x > x_range[2] || y < y_range[1] || y > y_range[2]) break
      pts[[out_idx]] <- c(UMAP_1 = x, UMAP_2 = y)
      out_idx <- out_idx + 1L
    }
    pts <- pts[seq_len(max(out_idx - 1L, 0L))]
    if (length(pts) == 0) return(matrix(numeric(0), ncol = 2))
    do.call(rbind, pts)
  }

  stream_rows <- vector("list", nrow(seed_df))
  out_idx <- 1L
  for (i in seq_len(nrow(seed_df))) {
    x0 <- seed_df$UMAP_1[i]
    y0 <- seed_df$UMAP_2[i]
    back <- trace_direction(x0, y0, direction = -1, n_steps = backward_steps)
    fwd <- trace_direction(x0, y0, direction = 1, n_steps = forward_steps)
    seed_mat <- matrix(c(x0, y0), nrow = 1, dimnames = list(NULL, c("UMAP_1", "UMAP_2")))
    path_mat <- rbind(if (nrow(back) > 0) back[nrow(back):1, , drop = FALSE] else NULL, seed_mat, fwd)
    if (nrow(path_mat) < 4) next
    stream_rows[[out_idx]] <- data.frame(
      stream_id = out_idx,
      path_order = seq_len(nrow(path_mat)),
      UMAP_1 = path_mat[, "UMAP_1"],
      UMAP_2 = path_mat[, "UMAP_2"],
      stringsAsFactors = FALSE
    )
    out_idx <- out_idx + 1L
  }
  dplyr::bind_rows(stream_rows)
}

cleanup_scvelo_point_to_point_outputs <- function(out_dir) {
  old_files <- c(
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_flow_vectors.pdf")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_flow_vectors.png"))
  )
  if (length(old_files) > 0) unlink(old_files)
  invisible(old_files)
}

cleanup_scvelo_flat_overlay_outputs <- function(out_dir) {
  old_files <- c(
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_grid.pdf")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_grid.png")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_stream.pdf")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_stream.png")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_cluster_grid.pdf")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_cluster_grid.png")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_cluster_stream.pdf")),
    Sys.glob(file.path(out_dir, "*_velocity_pseudotime_cluster_stream.png"))
  )
  if (length(old_files) > 0) unlink(old_files)
  invisible(old_files)
}

plot_scvelo_flow_one <- function(h5ad_file, label, output_prefix, out_dir, out_umap_dir) {
  flow <- read_scvelo_h5ad_flow(h5ad_file)
  if (nrow(flow$points) == 0 || nrow(flow$vectors) == 0) {
    stop("No usable points or velocity vectors for: ", h5ad_file, call. = FALSE)
  }
  grid_df <- data.frame()
  grid_min_bin_n <- NA_integer_
  for (candidate_min_bin_n in c(5L, 3L, 2L, 1L)) {
    grid_df <- make_scvelo_velocity_grid(flow$vectors, min_bin_n = candidate_min_bin_n)
    if (nrow(grid_df) > 0) {
      grid_min_bin_n <- candidate_min_bin_n
      attr(grid_df, "min_bin_n") <- candidate_min_bin_n
      break
    }
  }
  if (nrow(grid_df) == 0) {
    stop("No usable grid arrows after smoothing velocity vectors for: ", h5ad_file, call. = FALSE)
  }
	  stream_df <- make_scvelo_velocity_streamlines(grid_df)
	  cluster_layers <- make_scvelo_cluster_layers(flow$points, core_fraction = 0.90)
	  cluster_labels <- cluster_layers$labels
	  flow_umap_x <- c(flow$points$UMAP_1, grid_df$UMAP_1, grid_df$UMAP_1_to, stream_df$UMAP_1)
	  flow_umap_y <- c(flow$points$UMAP_2, grid_df$UMAP_2, grid_df$UMAP_2_to, stream_df$UMAP_2)

	  base_plot <- ggplot(flow$points, aes(x = UMAP_1, y = UMAP_2, color = velocity_pseudotime)) +
	    geom_point(size = 0.18, alpha = 0.62) +
    scale_color_gradientn(
	      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
	      name = "scVelo velocity\npseudotime"
	    ) +
	    coord_square_umap(flow_umap_x, flow_umap_y) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()

  p_grid <- base_plot +
    geom_segment(
      data = grid_df,
      aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      linewidth = 0.28,
      color = "black",
      alpha = 0.72,
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
    ) +
    labs(
      title = paste0("scVelo ", label, " RNA velocity grid over velocity pseudotime"),
      subtitle = "Cells are colored by scVelo velocity_pseudotime; arrows are binned and smoothed UMAP velocity directions from scVelo velocity_graph.",
      x = "UMAP 1",
      y = "UMAP 2"
    )
  if (nrow(cluster_labels) > 0) {
    p_grid <- p_grid +
      geom_label(
        data = cluster_labels,
        aes(x = UMAP_1, y = UMAP_2, label = clusters),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold",
        linewidth = 0.18,
        fill = "white",
        color = "black",
        alpha = 0.82
      )
  }
	  save_square_panel_both(p_grid, file.path(out_dir, paste0(output_prefix, "_velocity_pseudotime_grid")), panel_size = 5.2)
	  save_square_panel_both(p_grid, file.path(out_umap_dir, paste0(output_prefix, "_velocity_pseudotime_grid")), panel_size = 5.2)

  p_stream <- base_plot +
    labs(
      title = paste0("scVelo ", label, " RNA velocity stream over velocity pseudotime"),
      subtitle = "Cells are colored by scVelo velocity_pseudotime; stream lines trace the same smoothed UMAP velocity field as the grid plot.",
      x = "UMAP 1",
      y = "UMAP 2"
    )
  if (nrow(stream_df) > 0) {
    p_stream <- p_stream +
      geom_path(
        data = stream_df,
        aes(x = UMAP_1, y = UMAP_2, group = stream_id),
        inherit.aes = FALSE,
        linewidth = 0.28,
        color = "black",
        alpha = 0.62,
        lineend = "round",
        arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
      )
  }
  if (nrow(cluster_labels) > 0) {
    p_stream <- p_stream +
      geom_label(
        data = cluster_labels,
        aes(x = UMAP_1, y = UMAP_2, label = clusters),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold",
        linewidth = 0.18,
        fill = "white",
        color = "black",
        alpha = 0.82
      )
  }
	  save_square_panel_both(p_stream, file.path(out_dir, paste0(output_prefix, "_velocity_pseudotime_stream")), panel_size = 5.2)
	  save_square_panel_both(p_stream, file.path(out_umap_dir, paste0(output_prefix, "_velocity_pseudotime_stream")), panel_size = 5.2)

  tile_df <- cluster_layers$tiles
  hull_df <- cluster_layers$hulls
  boundary_df <- cluster_layers$boundaries
  tile_width <- attr(tile_df, "tile_width")
  tile_height <- attr(tile_df, "tile_height")
  if (
    nrow(tile_df) > 0 &&
      length(tile_width) == 1L &&
      length(tile_height) == 1L &&
      is.finite(tile_width) &&
      is.finite(tile_height)
  ) {
    cluster_base <- ggplot() +
      geom_tile(
        data = tile_df,
        aes(x = UMAP_1, y = UMAP_2, fill = mean_time),
        width = tile_width,
        height = tile_height,
        alpha = 0.78
      ) +
	      scale_fill_gradientn(
	        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
	        name = "Mean scVelo\nvelocity pseudotime"
	      ) +
	      coord_square_umap(flow_umap_x, flow_umap_y) +
	      theme_extra(9) +
	      theme(axis.text.x = element_text(angle = 0)) +
	      theme_square_umap_panel()
    if (nrow(boundary_df) > 0) {
      cluster_base <- cluster_base +
        geom_segment(
          data = boundary_df,
          aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
          inherit.aes = FALSE,
          color = "grey12",
          linewidth = 0.12,
          alpha = 0.42
        )
    }

    p_cluster_grid <- cluster_base +
      geom_segment(
        data = grid_df,
        aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
        inherit.aes = FALSE,
        linewidth = 0.30,
        color = "black",
        alpha = 0.82,
        arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
      ) +
      labs(
        title = paste0("scVelo ", label, " RNA velocity grid over cluster-level velocity pseudotime"),
        subtitle = "Cluster regions show the densest grid cells covering 90% of each cluster; flow is drawn above the cluster time layer.",
        x = "UMAP 1",
        y = "UMAP 2"
      )
    if (nrow(cluster_labels) > 0) {
      p_cluster_grid <- p_cluster_grid +
        geom_label(
          data = cluster_labels,
          aes(x = UMAP_1, y = UMAP_2, label = clusters),
          inherit.aes = FALSE,
          size = 3,
          fontface = "bold",
          linewidth = 0.18,
          fill = "white",
          color = "black",
          alpha = 0.86
        )
    }
	    save_square_panel_both(p_cluster_grid, file.path(out_dir, paste0(output_prefix, "_velocity_pseudotime_cluster_grid")), panel_size = 5.2)
	    save_square_panel_both(p_cluster_grid, file.path(out_umap_dir, paste0(output_prefix, "_velocity_pseudotime_cluster_grid")), panel_size = 5.2)

    p_cluster_stream <- cluster_base +
      labs(
        title = paste0("scVelo ", label, " RNA velocity stream over cluster-level velocity pseudotime"),
        subtitle = "Cluster regions show the densest grid cells covering 90% of each cluster; stream lines are drawn above the cluster time layer.",
        x = "UMAP 1",
        y = "UMAP 2"
      )
    if (nrow(stream_df) > 0) {
      p_cluster_stream <- p_cluster_stream +
        geom_path(
          data = stream_df,
          aes(x = UMAP_1, y = UMAP_2, group = stream_id),
          inherit.aes = FALSE,
          linewidth = 0.30,
          color = "black",
          alpha = 0.70,
          lineend = "round",
          arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
        )
    }
    if (nrow(cluster_labels) > 0) {
      p_cluster_stream <- p_cluster_stream +
        geom_label(
          data = cluster_labels,
          aes(x = UMAP_1, y = UMAP_2, label = clusters),
          inherit.aes = FALSE,
          size = 3,
          fontface = "bold",
          linewidth = 0.18,
          fill = "white",
          color = "black",
          alpha = 0.86
        )
    }
	    save_square_panel_both(p_cluster_stream, file.path(out_dir, paste0(output_prefix, "_velocity_pseudotime_cluster_stream")), panel_size = 5.2)
	    save_square_panel_both(p_cluster_stream, file.path(out_umap_dir, paste0(output_prefix, "_velocity_pseudotime_cluster_stream")), panel_size = 5.2)
  }

  data.frame(
    h5ad = h5ad_file,
    label = label,
    output_prefix = output_prefix,
    n_obs = flow$n_obs,
    n_points_plotted = nrow(flow$points),
    n_vectors_used = nrow(flow$vectors),
    n_grid_arrows_plotted = nrow(grid_df),
    grid_min_bin_n = grid_min_bin_n,
    n_streamlines_plotted = length(unique(stream_df$stream_id)),
    n_clusters_labeled = nrow(cluster_labels),
    n_cluster_tiles_plotted = nrow(tile_df),
    n_cluster_hull_points = nrow(hull_df),
    n_cluster_boundary_segments = nrow(boundary_df),
    cluster_core_fraction = attr(tile_df, "core_fraction"),
    status = "ok",
    message = "",
    stringsAsFactors = FALSE
  )
}

safe_file_stub <- function(x) {
  x <- as_clean_chr(x)
  x[is.na(x)] <- "NA"
  x <- gsub("mg/kg", "mg_kg", x, fixed = TRUE)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

build_tn_ploidy_flow_manifest <- function(root_dir) {
  group_rows <- expand.grid(
    TN = c("CellLine", "Tumor"),
    Ploidy = c("2N", "4N"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      plot_scope = "group_cluster_time",
      Dose = NA_character_,
      source_group = paste(.data$TN, .data$Ploidy),
      h5ad = file.path(root_dir, "01_velocity_groups", .data$TN, .data$Ploidy, "scvelo_result.h5ad"),
      output_subdir = "group_cluster_time",
      output_prefix = paste0(tolower(.data$TN), "_", .data$Ploidy),
      label = paste(.data$TN, .data$Ploidy),
      fill_source = "median_time"
    )

  compare_rows <- data.frame(TN = c("CellLine", "Tumor"), stringsAsFactors = FALSE) %>%
    dplyr::mutate(
      plot_scope = "ploidy_comparison",
      Ploidy = NA_character_,
      Dose = NA_character_,
      source_group = paste(.data$TN, "4N-2N"),
      h5ad = file.path(root_dir, "01_velocity_groups", .data$TN, "ploidy_all", "scvelo_result.h5ad"),
      output_subdir = "ploidy_comparison",
      output_prefix = paste0(tolower(.data$TN), "_4N_minus_2N"),
      label = paste(.data$TN, "4N minus 2N"),
      fill_source = "median_4N_minus_2N"
    )

  dplyr::bind_rows(group_rows, compare_rows) %>%
    dplyr::mutate(h5ad_exists = file.exists(.data$h5ad))
}

get_tn_ploidy_cluster_fill <- function(manifest_row, summary_df, comparison_df) {
  plot_scope <- manifest_row$plot_scope[1]
  if (identical(plot_scope, "group_cluster_time")) {
    df <- summary_df %>%
      dplyr::mutate(
        TN = as_clean_chr(.data$TN),
        Ploidy = as_clean_chr(.data$Ploidy),
        Dose = as_clean_chr(.data$Dose),
        clusters = as_clean_chr(.data$clusters),
        median_time = safe_numeric(.data$median_time)
	      ) %>%
      dplyr::filter(.data$TN == manifest_row$TN[1], .data$Ploidy == manifest_row$Ploidy[1])
    if (!is.na(manifest_row$Dose[1]) && nzchar(as.character(manifest_row$Dose[1]))) {
      df <- df %>% dplyr::filter(.data$Dose == manifest_row$Dose[1])
    } else if (manifest_row$TN[1] == "CellLine") {
      df <- df %>% dplyr::filter(is.na(.data$Dose))
    }
    return(df %>%
      dplyr::transmute(
        clusters = .data$clusters,
        plot_value = .data$median_time,
        n_cells = .data$n_cells,
        value_label = sprintf("%.3f", .data$median_time)
      ))
  }

  df <- comparison_df %>%
    dplyr::mutate(
      TN = as_clean_chr(.data$TN),
      Dose = as_clean_chr(.data$Dose),
      clusters = as_clean_chr(.data$clusters),
      median_4N_minus_2N = safe_numeric(.data$median_4N_minus_2N)
    )
  if (manifest_row$TN[1] == "Tumor" && !is.na(manifest_row$Dose[1]) && nzchar(as.character(manifest_row$Dose[1]))) {
    df <- df %>%
      dplyr::filter(
        .data$comparison_scope == "Tumor_2N_vs_4N_by_Dose",
        .data$Dose == manifest_row$Dose[1]
      )
  } else if (manifest_row$TN[1] == "Tumor") {
    df <- df %>%
      dplyr::filter(.data$comparison_scope == "Tumor_2N_vs_4N_pooled_Dose")
  } else {
    df <- df %>%
      dplyr::filter(.data$comparison_scope == "CellLine_2N_vs_4N")
  }
  df %>%
    dplyr::transmute(
      clusters = .data$clusters,
      plot_value = .data$median_4N_minus_2N,
      n_cells = pmin(.data$n_cells_2N, .data$n_cells_4N, na.rm = TRUE),
      value_label = sprintf("%.3f", .data$median_4N_minus_2N)
    )
}

plot_scvelo_cluster_value_flow_one <- function(h5ad_file, label, output_prefix, out_dir, cluster_value_df, fill_label, fill_mode = c("time", "difference")) {
  fill_mode <- match.arg(fill_mode)
  flow <- read_scvelo_h5ad_flow(h5ad_file)
  if (nrow(flow$points) == 0 || nrow(flow$vectors) == 0) {
    stop("No usable points or velocity vectors for: ", h5ad_file, call. = FALSE)
  }

  grid_df <- data.frame()
  grid_min_bin_n <- NA_integer_
  for (candidate_min_bin_n in c(5L, 3L, 2L, 1L)) {
    grid_df <- make_scvelo_velocity_grid(flow$vectors, min_bin_n = candidate_min_bin_n)
    if (nrow(grid_df) > 0) {
      grid_min_bin_n <- candidate_min_bin_n
      attr(grid_df, "min_bin_n") <- candidate_min_bin_n
      break
    }
  }
  if (nrow(grid_df) == 0) {
    stop("No usable grid arrows after smoothing velocity vectors for: ", h5ad_file, call. = FALSE)
  }

  stream_df <- make_scvelo_velocity_streamlines(grid_df)
  cluster_layers <- make_scvelo_cluster_layers(flow$points, core_fraction = 0.90)
	  tile_df <- cluster_layers$tiles
	  boundary_df <- cluster_layers$boundaries
	  cluster_labels <- cluster_layers$labels
	  tile_width <- attr(tile_df, "tile_width")
	  tile_height <- attr(tile_df, "tile_height")
	  flow_umap_x <- c(flow$points$UMAP_1, grid_df$UMAP_1, grid_df$UMAP_1_to, stream_df$UMAP_1)
	  flow_umap_y <- c(flow$points$UMAP_2, grid_df$UMAP_2, grid_df$UMAP_2_to, stream_df$UMAP_2)

  if (nrow(tile_df) == 0 || !is.finite(tile_width) || !is.finite(tile_height)) {
    stop("No usable cluster core tiles for: ", h5ad_file, call. = FALSE)
  }

  cluster_value_df <- cluster_value_df %>%
    dplyr::mutate(
      clusters = as_clean_chr(.data$clusters),
      plot_value = safe_numeric(.data$plot_value)
    ) %>%
    dplyr::filter(!is.na(.data$clusters))

  tile_df <- tile_df %>%
    dplyr::left_join(cluster_value_df, by = "clusters")
  cluster_labels <- cluster_labels %>%
    dplyr::left_join(cluster_value_df, by = "clusters")

  cluster_base <- ggplot() +
	    geom_tile(
	      data = tile_df,
	      aes(x = UMAP_1, y = UMAP_2, fill = plot_value),
	      width = tile_width,
	      height = tile_height,
	      alpha = 0.82
	    ) +
	    coord_square_umap(flow_umap_x, flow_umap_y) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()

  if (fill_mode == "difference") {
    cluster_base <- cluster_base +
      scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        na.value = "grey92",
        name = fill_label
      )
  } else {
    cluster_base <- cluster_base +
      scale_fill_gradientn(
        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
        na.value = "grey92",
        name = fill_label
      )
  }

  if (nrow(boundary_df) > 0) {
    cluster_base <- cluster_base +
      geom_segment(
        data = boundary_df,
        aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
        inherit.aes = FALSE,
        color = "grey12",
        linewidth = 0.12,
        alpha = 0.45
      )
  }

  p_grid <- cluster_base +
    geom_segment(
      data = grid_df,
      aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      linewidth = 0.30,
      color = "black",
      alpha = 0.82,
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
    ) +
    labs(
      title = paste0("scVelo ", label, " cluster result UMAP with RNA velocity grid flow"),
      subtitle = "Cluster regions cover the densest 90% of cells per cluster; cluster labels are shown on top of the flow layer.",
      x = "UMAP 1",
      y = "UMAP 2"
    )

  p_stream <- cluster_base +
    labs(
      title = paste0("scVelo ", label, " cluster result UMAP with RNA velocity stream flow"),
      subtitle = "Cluster regions cover the densest 90% of cells per cluster; cluster labels are shown on top of the flow layer.",
      x = "UMAP 1",
      y = "UMAP 2"
    )
  if (nrow(stream_df) > 0) {
    p_stream <- p_stream +
      geom_path(
        data = stream_df,
        aes(x = UMAP_1, y = UMAP_2, group = stream_id),
        inherit.aes = FALSE,
        linewidth = 0.30,
        color = "black",
        alpha = 0.70,
        lineend = "round",
        arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
      )
  }

  if (nrow(cluster_labels) > 0) {
    label_layer <- geom_label(
      data = cluster_labels,
      aes(x = UMAP_1, y = UMAP_2, label = clusters),
      inherit.aes = FALSE,
      size = 3,
      fontface = "bold",
      linewidth = 0.18,
      fill = "white",
      color = "black",
      alpha = 0.88
    )
    p_grid <- p_grid + label_layer
    p_stream <- p_stream + label_layer
  }

	  save_square_panel_both(p_grid, file.path(out_dir, paste0(output_prefix, "_cluster_result_flow_grid")), panel_size = 5.2)
	  save_square_panel_both(p_stream, file.path(out_dir, paste0(output_prefix, "_cluster_result_flow_stream")), panel_size = 5.2)

  data.frame(
    h5ad = h5ad_file,
    label = label,
    output_prefix = output_prefix,
    n_obs = flow$n_obs,
    n_points_plotted = nrow(flow$points),
    n_vectors_used = nrow(flow$vectors),
    n_grid_arrows_plotted = nrow(grid_df),
    grid_min_bin_n = grid_min_bin_n,
    n_streamlines_plotted = length(unique(stream_df$stream_id)),
    n_clusters_labeled = nrow(cluster_labels),
    n_clusters_with_fill_value = sum(is.finite(cluster_labels$plot_value)),
    status = "ok",
    message = "",
    stringsAsFactors = FALSE
  )
}

plot_tn_ploidy_flow_umaps <- function(root_dir, summary_df, comparison_df, out_dir) {
  manifest <- build_tn_ploidy_flow_manifest(root_dir)
  write_table_csv(manifest, file.path(out_dir, "tn_ploidy_flow_umap_manifest.csv"))
  if (nrow(manifest) == 0) return(invisible(data.frame()))

  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    message("Drawing TN/Ploidy cluster flow UMAP: ", i, "/", nrow(manifest), " ", manifest$label[i])
    plot_out_dir <- .ensure_dir(file.path(out_dir, manifest$output_subdir[i]))
    status_rows[[i]] <- tryCatch({
      if (!isTRUE(manifest$h5ad_exists[i])) {
        stop("Missing h5ad: ", manifest$h5ad[i], call. = FALSE)
      }
      cluster_fill <- get_tn_ploidy_cluster_fill(manifest[i, , drop = FALSE], summary_df, comparison_df)
      fill_mode <- if (identical(manifest$fill_source[i], "median_4N_minus_2N")) "difference" else "time"
      fill_label <- if (identical(fill_mode, "difference")) "4N - 2N\nmedian time" else "Cluster\nmedian time"
      result <- plot_scvelo_cluster_value_flow_one(
        h5ad_file = manifest$h5ad[i],
        label = manifest$label[i],
        output_prefix = manifest$output_prefix[i],
        out_dir = plot_out_dir,
        cluster_value_df = cluster_fill,
        fill_label = fill_label,
        fill_mode = fill_mode
      )
      result$plot_scope <- manifest$plot_scope[i]
      result$TN <- manifest$TN[i]
      result$Ploidy <- manifest$Ploidy[i]
      result$Dose <- manifest$Dose[i]
      result$output_subdir <- manifest$output_subdir[i]
      result
    }, error = function(e) {
      data.frame(
        h5ad = manifest$h5ad[i],
        label = manifest$label[i],
        output_prefix = manifest$output_prefix[i],
        n_obs = NA_integer_,
        n_points_plotted = NA_integer_,
        n_vectors_used = NA_integer_,
        n_grid_arrows_plotted = NA_integer_,
        grid_min_bin_n = NA_integer_,
        n_streamlines_plotted = NA_integer_,
        n_clusters_labeled = NA_integer_,
        n_clusters_with_fill_value = NA_integer_,
        status = "failed",
        message = conditionMessage(e),
        plot_scope = manifest$plot_scope[i],
        TN = manifest$TN[i],
        Ploidy = manifest$Ploidy[i],
        Dose = manifest$Dose[i],
        output_subdir = manifest$output_subdir[i],
        stringsAsFactors = FALSE
      )
    })
  }
  status_df <- dplyr::bind_rows(status_rows)
  write_table_csv(status_df, file.path(out_dir, "tn_ploidy_flow_umap_status.csv"))
  if (any(status_df$status != "ok")) {
    warning("Some TN/Ploidy cluster flow UMAPs failed; see tn_ploidy_flow_umap_status.csv.", call. = FALSE)
  }
  invisible(status_df)
}

add_method_group_label <- function(df) {
  df <- ensure_columns(df, c(
    "source_context", "source_label", "trajectory_analysis_group", "AnalysisGroup",
    "trajectory_group", "trajectory_context", "trajectory_tn_scope", "trajectory_ploidy_scope",
    "trajectory_branch", "TN", "Ploidy", "sampleID", "sample"
  ))
  df %>%
    dplyr::mutate(
      analysis_group = coalesce_chr(.data$trajectory_analysis_group, .data$AnalysisGroup, .data$source_label, .data$trajectory_group, .data$trajectory_context),
      analysis_group = dplyr::case_when(
        .data$analysis_group %in% c("All_cells", "Tumor", "CellLine", "CellLine_ploidy_all", "CellLine_2N", "CellLine_4N", "Tumor_ploidy_all", "Tumor_2N", "Tumor_4N") ~ .data$analysis_group,
        .data$source_label %in% c("All_cells", "Tumor", "CellLine", "CellLine/ploidy_all", "CellLine/2N", "CellLine/4N", "Tumor/ploidy_all", "Tumor/2N", "Tumor/4N") ~ gsub("/", "_", .data$source_label, fixed = TRUE),
        TRUE ~ .data$analysis_group
      ),
      analysis_tn_scope = coalesce_chr(
        .data$trajectory_tn_scope,
        dplyr::case_when(
          .data$analysis_group == "All_cells" ~ "All",
          grepl("^CellLine(_|$)", .data$analysis_group) ~ "CellLine",
          grepl("^Tumor(_|$)", .data$analysis_group) ~ "Tumor",
          TRUE ~ as_clean_chr(.data$TN)
        )
      ),
      analysis_ploidy_scope = coalesce_chr(
        .data$trajectory_ploidy_scope,
        dplyr::case_when(
          grepl("_(2N|4N)$", .data$analysis_group) ~ sub("^.*_", "", .data$analysis_group),
          TRUE ~ "ploidy_all"
        )
      ),
      TN = as_clean_chr(.data$TN),
      Ploidy = as_clean_chr(.data$Ploidy),
      sampleID = coalesce_chr(.data$sampleID, .data$sample)
    )
}

summarize_group_clusters <- function(df, split_cols = character(), time_col = "velocity_pseudotime") {
  needed <- unique(c(split_cols, "clusters", time_col, "sampleID"))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(
      clusters = as_clean_chr(.data$clusters),
      .time = safe_numeric(.data[[time_col]]),
      sampleID = as_clean_chr(.data$sampleID)
    ) %>%
    dplyr::filter(!is.na(.data$clusters), is.finite(.data$.time))
  if (nrow(df) == 0) return(data.frame())
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_cols)), .data$clusters) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_samples = dplyr::n_distinct(.data$sampleID[!is.na(.data$sampleID)]),
      mean_time = mean(.data$.time, na.rm = TRUE),
      median_time = stats::median(.data$.time, na.rm = TRUE),
      q05_time = stats::quantile(.data$.time, 0.05, na.rm = TRUE, names = FALSE),
      q25_time = stats::quantile(.data$.time, 0.25, na.rm = TRUE, names = FALSE),
      q75_time = stats::quantile(.data$.time, 0.75, na.rm = TRUE, names = FALSE),
      q95_time = stats::quantile(.data$.time, 0.95, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    ) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(split_cols))) %>%
    dplyr::mutate(time_rank_early_to_late = rank(.data$median_time, ties.method = "average")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::across(dplyr::all_of(split_cols)), .data$time_rank_early_to_late)
}

build_cluster_time_layers_for_split <- function(df, time_col = "velocity_pseudotime", split_cols = character(), core_fraction = 0.90) {
  needed <- unique(c("UMAP_1", "UMAP_2", "clusters", time_col, split_cols))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(
      UMAP_1 = safe_numeric(.data$UMAP_1),
      UMAP_2 = safe_numeric(.data$UMAP_2),
      clusters = as_clean_chr(.data$clusters),
      velocity_pseudotime = safe_numeric(.data[[time_col]])
    ) %>%
    dplyr::filter(
      is.finite(.data$UMAP_1),
      is.finite(.data$UMAP_2),
      is.finite(.data$velocity_pseudotime),
      !is.na(.data$clusters)
    )
  if (nrow(df) == 0) {
    return(list(points = data.frame(), tiles = data.frame(), boundaries = data.frame(), labels = data.frame()))
  }
  if (length(split_cols) == 0) {
    df$.split_key <- "all"
  } else {
    split_label_df <- as.data.frame(lapply(df[, split_cols, drop = FALSE], as.character), stringsAsFactors = FALSE)
    df$.split_key <- apply(split_label_df, 1, paste, collapse = "\r")
  }

  split_keys <- unique(df$.split_key)
  tile_rows <- list()
  boundary_rows <- list()
  label_rows <- list()
  for (i in seq_along(split_keys)) {
    sub_df <- df[df$.split_key == split_keys[i], , drop = FALSE]
    layers <- make_scvelo_cluster_layers(sub_df, core_fraction = core_fraction)
    split_values <- if (length(split_cols) == 0) data.frame(.dummy_split = "all") else sub_df[1, split_cols, drop = FALSE]

    tile_df <- layers$tiles
    if (nrow(tile_df) > 0) {
      tile_width <- attr(tile_df, "tile_width")
      tile_height <- attr(tile_df, "tile_height")
      tile_df$tile_width <- tile_width
      tile_df$tile_height <- tile_height
      tile_df$xmin <- tile_df$UMAP_1 - tile_width / 2
      tile_df$xmax <- tile_df$UMAP_1 + tile_width / 2
      tile_df$ymin <- tile_df$UMAP_2 - tile_height / 2
      tile_df$ymax <- tile_df$UMAP_2 + tile_height / 2
      tile_rows[[i]] <- cbind(split_values[rep(1, nrow(tile_df)), , drop = FALSE], tile_df, stringsAsFactors = FALSE)
    }

    boundary_df <- layers$boundaries
    if (nrow(boundary_df) > 0) {
      boundary_rows[[i]] <- cbind(split_values[rep(1, nrow(boundary_df)), , drop = FALSE], boundary_df, stringsAsFactors = FALSE)
    }

    label_df <- layers$labels
    if (nrow(label_df) > 0) {
      label_rows[[i]] <- cbind(split_values[rep(1, nrow(label_df)), , drop = FALSE], label_df, stringsAsFactors = FALSE)
    }
  }

  points <- df
  if (nrow(points) > 80000) {
    set.seed(1234)
    points <- points[sort(sample(seq_len(nrow(points)), 80000)), , drop = FALSE]
  }
  list(
    points = points,
    tiles = dplyr::bind_rows(tile_rows),
    boundaries = dplyr::bind_rows(boundary_rows),
    labels = dplyr::bind_rows(label_rows)
  )
}

plot_cluster_time_core_umap <- function(df, label, output_stub, time_col = "velocity_pseudotime", split_cols = character(), fill_label = "Mean time") {
  layers <- build_cluster_time_layers_for_split(df, time_col = time_col, split_cols = split_cols, core_fraction = 0.90)
  if (nrow(layers$points) == 0) return(invisible(data.frame()))
  umap_coord <- coord_square_umap(layers$points$UMAP_1, layers$points$UMAP_2)

  p <- ggplot() +
    geom_point(
      data = layers$points,
      aes(x = UMAP_1, y = UMAP_2),
      inherit.aes = FALSE,
      color = "grey88",
      size = 0.08,
      alpha = 0.25
    )
  if (nrow(layers$tiles) > 0) {
    p <- p +
      geom_rect(
        data = layers$tiles,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = mean_time),
        inherit.aes = FALSE,
        alpha = 0.82
      ) +
      scale_fill_gradientn(
        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
        na.value = "grey92",
        name = fill_label
      )
  } else {
    p <- p +
      geom_point(
        data = layers$points,
        aes(x = UMAP_1, y = UMAP_2, color = velocity_pseudotime),
        inherit.aes = FALSE,
        size = 0.18,
        alpha = 0.7
      ) +
      scale_color_gradientn(
        colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
        name = fill_label
      )
  }
  if (nrow(layers$boundaries) > 0) {
    p <- p +
      geom_segment(
        data = layers$boundaries,
        aes(x = UMAP_1, y = UMAP_2, xend = UMAP_1_to, yend = UMAP_2_to),
        inherit.aes = FALSE,
        color = "grey12",
        linewidth = 0.12,
        alpha = 0.45
      )
  }
  if (nrow(layers$labels) > 0) {
    p <- p +
      geom_label(
        data = layers$labels,
        aes(x = UMAP_1, y = UMAP_2, label = clusters),
        inherit.aes = FALSE,
        size = 3,
        fontface = "bold",
        linewidth = 0.18,
        fill = "white",
        color = "black",
        alpha = 0.88
      )
  }
  if (length(split_cols) > 0) {
    p <- p + facet_wrap(stats::as.formula(paste("~", paste(split_cols, collapse = "+"))))
  }
		  p <- p +
		    umap_coord +
		    labs(
	      title = label,
	      x = "UMAP 1",
	      y = "UMAP 2"
	    ) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()
	  save_square_panel_both(p, output_stub, panel_size = ifelse(length(split_cols) > 0, 4.8, 5.2))

  data.frame(
    output_stub = output_stub,
    label = label,
    n_points = nrow(layers$points),
    n_cluster_tiles = nrow(layers$tiles),
    n_cluster_labels = nrow(layers$labels),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

choose_time_bin_width <- function(x, n_bins = 200L) {
  x <- safe_numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  rng <- range(x, na.rm = TRUE)
  span <- diff(rng)
  if (!is.finite(span) || span <= 0) return(0.005)
  span / n_bins
}

make_frequency_histogram_data <- function(df, time_col, group_cols, bin_width = NULL, n_bins = 200L) {
  needed <- unique(c(time_col, group_cols))
  df <- ensure_columns(df, needed)
  df <- df %>%
    dplyr::mutate(.time = safe_numeric(.data[[time_col]])) %>%
    dplyr::filter(is.finite(.data$.time), dplyr::if_all(dplyr::all_of(group_cols), ~ !is.na(as_clean_chr(.x))))
  if (nrow(df) == 0) return(data.frame())
  if (is.null(bin_width) || !is.finite(bin_width) || bin_width <= 0) {
    bin_width <- choose_time_bin_width(df$.time, n_bins = n_bins)
  }
  rng <- range(df$.time, na.rm = TRUE)
  breaks <- seq(rng[1], rng[2] + bin_width, by = bin_width)
  if (length(breaks) < 3L) breaks <- seq(rng[1] - bin_width, rng[2] + bin_width, by = bin_width)
  df$.bin <- pmin(pmax(findInterval(df$.time, breaks, all.inside = TRUE), 1L), length(breaks) - 1L)
  df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_cols)), .data$.bin) %>%
    dplyr::summarise(n_cells = dplyr::n(), .groups = "drop_last") %>%
    dplyr::mutate(total_cells_in_group = sum(.data$n_cells, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      time_bin_left = breaks[.data$.bin],
      time_bin_right = breaks[.data$.bin + 1L],
      time_bin_mid = (.data$time_bin_left + .data$time_bin_right) / 2,
      bin_width = bin_width,
      cell_percent = 100 * .data$n_cells / .data$total_cells_in_group
    ) %>%
    dplyr::select(dplyr::all_of(group_cols), time_bin_left, time_bin_right, time_bin_mid, bin_width, n_cells, total_cells_in_group, cell_percent)
}

plot_ploidy_frequency_histogram <- function(hist_df, label, output_stub, x_label, facet_col = NULL) {
  if (nrow(hist_df) == 0) return(invisible(NULL))
  hist_df <- hist_df %>%
    dplyr::mutate(Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N")))
  p <- ggplot(hist_df, aes(x = time_bin_mid, y = cell_percent, fill = Ploidy, color = Ploidy)) +
    geom_col(aes(width = bin_width), position = "identity", alpha = 0.36, linewidth = 0.08) +
    geom_step(linewidth = 0.45, alpha = 0.95) +
    scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
    scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
    labs(
      title = label,
      subtitle = "Y axis is percent of cells within each Ploidy group; bins use a shared small time width.",
      x = x_label,
	        y = "Cell percent",
	        fill = "Ploidy",
	        color = "Ploidy"
	      ) +
	      theme_extra(9) +
	      theme(axis.text.x = element_text(angle = 0), legend.position = "bottom")
  if (!is.null(facet_col) && facet_col %in% names(hist_df)) {
    p <- p + facet_wrap(stats::as.formula(paste("~", facet_col)), scales = "free_y")
  }
  save_both(p, output_stub, width = ifelse(is.null(facet_col), 8.8, 12), height = ifelse(is.null(facet_col), 5.4, 8))
  invisible(p)
}

plot_sample_frequency_histogram_pages <- function(df, tn_label, output_pdf, output_png_dir, time_col = "velocity_pseudotime", method_label, time_label, x_label) {
  df <- ensure_columns(df, c(time_col, "Ploidy", "sampleID", "Dose"))
  first_present <- function(x) {
    x <- as_clean_chr(x)
    x <- x[!is.na(x)]
    if (length(x) == 0) NA_character_ else x[1]
  }
  ploidy_levels <- c("2N", "4N")
  dose_levels <- c("0mg/kg", "30mg/kg", "120mg/kg")
  df <- df %>%
    dplyr::mutate(
      Ploidy = as_clean_chr(.data$Ploidy),
      Dose = as_clean_chr(.data$Dose),
      sampleID = as_clean_chr(.data$sampleID),
      .time = safe_numeric(.data[[time_col]])
    ) %>%
    dplyr::filter(is.finite(.data$.time), !is.na(.data$Ploidy), !is.na(.data$sampleID))
  if (nrow(df) == 0) return(data.frame())
  sample_info <- df %>%
    dplyr::group_by(.data$sampleID) %>%
    dplyr::summarise(
      sample_ploidy = first_present(.data$Ploidy),
      sample_dose = first_present(.data$Dose),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      .ploidy_order = match(.data$sample_ploidy, ploidy_levels),
      .dose_order = match(.data$sample_dose, dose_levels),
      .ploidy_order = ifelse(is.na(.data$.ploidy_order), length(ploidy_levels) + 1L, .data$.ploidy_order),
      .dose_order = ifelse(is.na(.data$.dose_order), length(dose_levels) + 1L, .data$.dose_order)
    ) %>%
    dplyr::arrange(.data$.ploidy_order, .data$.dose_order, .data$sampleID) %>%
    dplyr::mutate(
      sample_ploidy_label = ifelse(is.na(.data$sample_ploidy), "Ploidy NA", .data$sample_ploidy),
      sample_dose_label = ifelse(is.na(.data$sample_dose), "Dose NA", .data$sample_dose),
      sample_facet_label = ifelse(
        is.na(.data$sample_dose),
        paste0(.data$sample_ploidy_label, "\n", .data$sampleID),
        paste0(.data$sample_ploidy_label, " | ", .data$sample_dose_label, "\n", .data$sampleID)
      )
    )
  sample_ids <- sample_info$sampleID
  sample_facet_levels <- sample_info$sample_facet_label
  sample_info_for_join <- sample_info %>%
    dplyr::select(
      "sampleID",
      SamplePloidyLabel = "sample_ploidy_label",
      SampleDoseLabel = "sample_dose_label",
      sampleFacetLabel = "sample_facet_label"
    )
  bin_width <- choose_time_bin_width(df$.time, n_bins = 200L)
  hist_all <- make_frequency_histogram_data(df, ".time", c("sampleID", "Ploidy"), bin_width = bin_width)
  hist_all <- hist_all %>%
    dplyr::mutate(sampleID = as_clean_chr(.data$sampleID)) %>%
    dplyr::left_join(sample_info_for_join, by = "sampleID") %>%
    dplyr::mutate(
      sampleID = factor(as.character(.data$sampleID), levels = sample_ids),
      Ploidy = factor(as.character(.data$Ploidy), levels = ploidy_levels),
      SamplePloidyLabel = factor(as.character(.data$SamplePloidyLabel), levels = c(ploidy_levels, "Ploidy NA")),
      SampleDoseLabel = factor(as.character(.data$SampleDoseLabel), levels = c(dose_levels, "Dose NA")),
      sampleFacetLabel = factor(as.character(.data$sampleFacetLabel), levels = sample_facet_levels)
    )
  write_table_csv(hist_all, sub("\\.pdf$", "_data.csv", output_pdf))

  grDevices::pdf(output_pdf, width = 8.8, height = 5.4, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)
  .ensure_dir(output_png_dir)
  for (sid in sample_ids) {
    hist_df <- hist_all %>% dplyr::filter(.data$sampleID == sid)
    p <- ggplot(hist_df, aes(x = time_bin_mid, y = cell_percent, fill = Ploidy, color = Ploidy)) +
      geom_col(aes(width = bin_width), position = "identity", alpha = 0.38, linewidth = 0.08) +
      geom_step(linewidth = 0.45) +
      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
      scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
      labs(
        title = paste0(method_label, " ", tn_label, " sample ", sid, " ", time_label, " distribution"),
        subtitle = "Y axis is percent of cells within each Ploidy group for this sample.",
        x = x_label,
        y = "Cell percent",
        fill = "Ploidy",
        color = "Ploidy"
	      ) +
	      theme_extra(9) +
	      theme(axis.text.x = element_text(angle = 0), legend.position = "bottom")
    print(p)
    ggsave(file.path(output_png_dir, paste0(safe_file_stub(sid), "_ploidy_time_frequency_histogram.png")), p, width = 8.8, height = 5.4, dpi = 300)
  }

	  combined_stub <- file.path(output_png_dir, "all_samples_ploidy_pseudotime_frequency_histograms")
	  combined_16x1_stub <- file.path(output_png_dir, "all_samples_ploidy_pseudotime_frequency_histograms_16x1")
	  if (nrow(hist_all) > 0) {
	    x_limits <- range(c(hist_all$time_bin_left, hist_all$time_bin_right), na.rm = TRUE)
	    y_limit <- max(hist_all$cell_percent, na.rm = TRUE)
    if (!is.finite(y_limit) || y_limit <= 0) y_limit <- 1
    p_combined <- ggplot(hist_all, aes(x = time_bin_mid, y = cell_percent, fill = Ploidy, color = Ploidy)) +
      geom_col(aes(width = bin_width), position = "identity", alpha = 0.38, linewidth = 0.08) +
      geom_step(linewidth = 0.45) +
      facet_wrap(~sampleFacetLabel, ncol = 4, strip.position = "right", scales = "fixed") +
      scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0))) +
      scale_y_continuous(limits = c(0, y_limit * 1.05), expand = ggplot2::expansion(mult = c(0, 0.02))) +
      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
      scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
      labs(
        title = paste0(method_label, " ", tn_label, " sample-level ", time_label, " distributions"),
        x = x_label,
        y = "Cell percent",
        fill = "Ploidy",
        color = "Ploidy"
      ) +
      guides(fill = guide_legend(nrow = 1, byrow = TRUE), color = "none") +
      theme_extra(7) +
      theme(
        axis.text.x = element_text(angle = 0),
        legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.right = element_text(angle = 0, hjust = 0.5),
        panel.spacing.x = grid::unit(0.45, "lines"),
        panel.spacing.y = grid::unit(0.45, "lines")
	      )
	    save_both(p_combined, combined_stub, width = 18, height = 6)

	    p_combined_16x1 <- ggplot(hist_all, aes(x = time_bin_mid, y = cell_percent, fill = Ploidy, color = Ploidy)) +
	      geom_col(aes(width = bin_width), position = "identity", alpha = 0.38, linewidth = 0.08) +
	      geom_step(linewidth = 0.45) +
	      facet_grid(rows = vars(SamplePloidyLabel, SampleDoseLabel, sampleID), scales = "fixed", drop = TRUE) +
	      scale_x_continuous(limits = x_limits, expand = ggplot2::expansion(mult = c(0, 0))) +
	      scale_y_continuous(limits = c(0, y_limit * 1.05), expand = ggplot2::expansion(mult = c(0, 0.02))) +
	      scale_fill_manual(values = c("2N" = "#3B73B9", "4N" = "#D95F02"), drop = FALSE) +
	      scale_color_manual(values = c("2N" = "#1D4E89", "4N" = "#A84600"), drop = FALSE) +
	      labs(
	        title = paste0(method_label, " ", tn_label, " sample-level ", time_label, " distributions"),
	        x = x_label,
	        y = "Cell percent",
	        fill = "Ploidy",
	        color = "Ploidy"
	      ) +
	      guides(fill = guide_legend(nrow = 1, byrow = TRUE), color = "none") +
	      theme_extra(6) +
	      theme(
	        axis.text.x = element_text(angle = 0),
	        legend.position = "bottom",
	        strip.placement = "outside",
	        strip.text.y.right = element_text(angle = 0, hjust = 0.5, size = 5.5),
	        panel.spacing.x = grid::unit(0.35, "lines"),
	        panel.spacing.y = grid::unit(0.28, "lines")
	      )
	    save_both(p_combined_16x1, combined_16x1_stub, width = 7.2, height = max(18, length(sample_ids) * 1.15))
	  }

	  data.frame(
	    output_pdf = output_pdf,
	    combined_output_stub = combined_stub,
	    combined_16x1_output_stub = combined_16x1_stub,
	    n_samples = length(sample_ids),
	    bin_width = bin_width,
    n_histogram_rows = nrow(hist_all),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

run_group_cluster_timing_outputs <- function(df, out_dir, method_label, time_col, source_context_filter, time_label, x_label) {
  tables_dir <- .ensure_dir(file.path(out_dir, "tables"))
  umap_dir <- .ensure_dir(file.path(out_dir, "umap_cluster_pseudotime"))
  hist_dir <- .ensure_dir(file.path(out_dir, "ploidy_pseudotime_histograms"))

	  df <- add_method_group_label(df) %>%
    dplyr::filter(
      .data$source_context %in% source_context_filter,
      .data$analysis_group %in% c(
        "All_cells", "Tumor", "CellLine",
        "CellLine_ploidy_all", "CellLine_2N", "CellLine_4N",
        "Tumor_ploidy_all", "Tumor_2N", "Tumor_4N"
      ),
      is.finite(safe_numeric(.data[[time_col]]))
    ) %>%
    dplyr::mutate(
      .time_value = safe_numeric(.data[[time_col]]),
      UMAP_1 = safe_numeric(.data$UMAP_1),
      UMAP_2 = safe_numeric(.data$UMAP_2),
      clusters = as_clean_chr(.data$clusters),
	      Ploidy = as_clean_chr(.data$Ploidy),
	      TN = as_clean_chr(.data$TN),
      analysis_group = as_clean_chr(.data$analysis_group),
      analysis_tn_scope = as_clean_chr(.data$analysis_tn_scope),
      analysis_ploidy_scope = as_clean_chr(.data$analysis_ploidy_scope),
      sampleID = coalesce_chr(.data$sampleID, .data$sample)
    )
  if (nrow(df) == 0) {
    warning("No ", method_label, " cells were available for group cluster timing outputs.", call. = FALSE)
    return(list(status = data.frame(), cluster_summary = data.frame(), histogram_summary = data.frame()))
  }

	  cluster_summary <- summarize_group_clusters(
    df,
    split_cols = c("analysis_group", "analysis_tn_scope", "analysis_ploidy_scope"),
    time_col = ".time_value"
  )
  cluster_ploidy_summary <- summarize_group_clusters(
    df %>% dplyr::filter(.data$analysis_tn_scope %in% c("Tumor", "CellLine"), .data$analysis_ploidy_scope == "ploidy_all", .data$Ploidy %in% c("2N", "4N")),
    split_cols = c("analysis_group", "analysis_tn_scope", "analysis_ploidy_scope", "Ploidy"),
    time_col = ".time_value"
  )
  write_table_csv(cluster_summary, file.path(tables_dir, "group_cluster_pseudotime_summary.csv"))
  write_table_csv(cluster_ploidy_summary, file.path(tables_dir, "group_ploidy_cluster_pseudotime_summary.csv"))

	  status_rows <- list()
  analysis_order <- c(
    "All_cells",
    "CellLine_ploidy_all", "CellLine_2N", "CellLine_4N",
    "Tumor_ploidy_all", "Tumor_2N", "Tumor_4N",
    "CellLine", "Tumor"
  )
  analysis_order <- analysis_order[analysis_order %in% unique(df$analysis_group)]
  analysis_order <- unique(c(analysis_order, setdiff(unique(df$analysis_group), analysis_order)))
  for (grp in analysis_order) {
    group_df <- df %>% dplyr::filter(.data$analysis_group == grp)
    if (nrow(group_df) == 0) next
    status_rows[[paste0("umap_", safe_file_stub(grp))]] <- plot_cluster_time_core_umap(
      group_df,
      label = paste0(method_label, " ", grp, " cluster ", time_label, " UMAP"),
      output_stub = file.path(umap_dir, paste0(safe_file_stub(grp), "_cluster_pseudotime_umap")),
      time_col = ".time_value",
      split_cols = character(0),
      fill_label = paste0("Mean ", method_label, "\n", time_label)
    )
  }

	  histogram_rows <- list()
  for (tn_label in c("Tumor", "CellLine")) {
    tn_df <- df %>%
      dplyr::filter(.data$analysis_tn_scope == tn_label, .data$analysis_ploidy_scope == "ploidy_all", .data$Ploidy %in% c("2N", "4N"))
    if (nrow(tn_df) == 0) next

    status_rows[[paste0("umap_", tn_label, "_ploidy")]] <- plot_cluster_time_core_umap(
      tn_df,
      label = paste0(method_label, " ", tn_label, " 2N and 4N cluster ", time_label, " UMAP"),
      output_stub = file.path(umap_dir, paste0(tolower(tn_label), "_2N_4N_cluster_pseudotime_umap")),
      time_col = ".time_value",
      split_cols = c("Ploidy"),
      fill_label = paste0("Mean ", method_label, "\n", time_label)
    )

    bin_width <- choose_time_bin_width(tn_df$.time_value, n_bins = 200L)
    hist_df <- make_frequency_histogram_data(tn_df, ".time_value", c("Ploidy"), bin_width = bin_width)
    write_table_csv(hist_df, file.path(tables_dir, paste0(tolower(tn_label), "_ploidy_pseudotime_frequency_histogram.csv")))
    plot_ploidy_frequency_histogram(
      hist_df,
      label = paste0(method_label, " ", tn_label, " 2N vs 4N ", time_label, " frequency histogram"),
      output_stub = file.path(hist_dir, paste0(tolower(tn_label), "_ploidy_pseudotime_frequency_histogram")),
      x_label = x_label
    )

    histogram_rows[[paste0(tn_label, "_sample_pages")]] <- plot_sample_frequency_histogram_pages(
      tn_df,
      tn_label = tn_label,
      output_pdf = file.path(hist_dir, paste0(tolower(tn_label), "_sample_ploidy_pseudotime_frequency_histograms.pdf")),
      output_png_dir = file.path(hist_dir, paste0(tolower(tn_label), "_sample_pages")),
      time_col = ".time_value",
      method_label = method_label,
      time_label = time_label,
      x_label = x_label
    )
  }

  status_df <- dplyr::bind_rows(status_rows)
  histogram_summary <- dplyr::bind_rows(histogram_rows)
  write_table_csv(status_df, file.path(tables_dir, "group_cluster_pseudotime_umap_status.csv"))
  write_table_csv(histogram_summary, file.path(tables_dir, "group_sample_histogram_status.csv"))
  list(status = status_df, cluster_summary = cluster_summary, histogram_summary = histogram_summary)
}

run_scvelo_flow_overlays <- function(root_dir, out_dir, out_umap_dir) {
  manifest <- build_scvelo_flow_manifest(root_dir)
  write_table_csv(manifest, file.path(out_dir, "scvelo_flow_overlay_manifest.csv"))
  if (nrow(manifest) == 0) {
    warning("No scvelo_result.h5ad files found for flow overlays.", call. = FALSE)
    return(invisible(NULL))
  }
  cleanup_scvelo_point_to_point_outputs(out_dir)
  cleanup_scvelo_point_to_point_outputs(out_umap_dir)
  cleanup_scvelo_flat_overlay_outputs(out_dir)
  cleanup_scvelo_flat_overlay_outputs(out_umap_dir)

  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    message("Drawing scVelo flow overlay: ", i, "/", nrow(manifest), " ", manifest$label[i])
    plot_out_dir <- .ensure_dir(file.path(out_dir, manifest$output_subdir[i]))
    plot_umap_dir <- .ensure_dir(file.path(out_umap_dir, manifest$output_subdir[i]))
    status_rows[[i]] <- tryCatch(
      {
        result <- plot_scvelo_flow_one(
          h5ad_file = manifest$h5ad[i],
          label = manifest$label[i],
          output_prefix = manifest$output_prefix[i],
          out_dir = plot_out_dir,
          out_umap_dir = plot_umap_dir
        )
        result$source_relative_dir <- manifest$source_relative_dir[i]
        result$source_top_level <- manifest$source_top_level[i]
        result$source_analysis_scope <- manifest$source_analysis_scope[i]
        result$source_task_name <- manifest$source_task_name[i]
        result$output_subdir <- manifest$output_subdir[i]
        result
      },
      error = function(e) {
        data.frame(
          h5ad = manifest$h5ad[i],
          label = manifest$label[i],
          source_relative_dir = manifest$source_relative_dir[i],
          source_top_level = manifest$source_top_level[i],
          source_analysis_scope = manifest$source_analysis_scope[i],
          source_task_name = manifest$source_task_name[i],
          output_subdir = manifest$output_subdir[i],
          output_prefix = manifest$output_prefix[i],
          n_obs = NA_integer_,
          n_points_plotted = NA_integer_,
          n_vectors_used = NA_integer_,
          n_grid_arrows_plotted = NA_integer_,
          grid_min_bin_n = NA_integer_,
          n_streamlines_plotted = NA_integer_,
          n_clusters_labeled = NA_integer_,
          n_cluster_tiles_plotted = NA_integer_,
          n_cluster_hull_points = NA_integer_,
          n_cluster_boundary_segments = NA_integer_,
          cluster_core_fraction = NA_real_,
          status = "failed",
          message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }
  status_df <- dplyr::bind_rows(status_rows) %>%
    dplyr::left_join(
      manifest %>%
        dplyr::select(dplyr::all_of(c(
          "h5ad",
          "source_relative_dir",
          "source_top_level",
          "source_analysis_scope",
          "source_task_name",
          "output_subdir"
        ))),
      by = "h5ad",
      suffix = c("", ".manifest")
    ) %>%
    dplyr::mutate(
      source_relative_dir = dplyr::coalesce(.data$source_relative_dir, .data$source_relative_dir.manifest),
      source_top_level = dplyr::coalesce(.data$source_top_level, .data$source_top_level.manifest),
      source_analysis_scope = dplyr::coalesce(.data$source_analysis_scope, .data$source_analysis_scope.manifest),
      source_task_name = dplyr::coalesce(.data$source_task_name, .data$source_task_name.manifest),
      output_subdir = dplyr::coalesce(.data$output_subdir, .data$output_subdir.manifest)
    ) %>%
    dplyr::select(-dplyr::ends_with(".manifest"))
  write_table_csv(status_df, file.path(out_dir, "scvelo_flow_overlay_status.csv"))
  if (any(status_df$status != "ok")) {
    warning("Some scVelo flow overlays failed; see scvelo_flow_overlay_status.csv.", call. = FALSE)
  }
  invisible(NULL)
}

read_csv_if_exists <- function(file) {
  if (!file.exists(file) || file.info(file)$size <= 0) return(data.frame())
  readr::read_csv(file, show_col_types = FALSE, progress = FALSE)
}

plot_monocle3_graph_overlays <- function(root_dir, out_dir, out_umap_dir) {
  group_dir <- file.path(root_dir, "01_pseudotrajectory_groups")
  legacy_dir <- file.path(root_dir, "04_by_dose_karyotype")
  search_dir <- if (dir.exists(group_dir)) group_dir else legacy_dir
  dirs <- list.dirs(search_dir, recursive = FALSE, full.names = TRUE)
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs) == 0) return(invisible(NULL))

  cells_list <- list()
  edges_list <- list()
  nodes_list <- list()
  manifest <- list()
  idx <- 1L
  for (dir_i in dirs) {
    cell_file <- file.path(dir_i, "cells_pseudotime.csv")
    edge_file <- file.path(dir_i, "principal_graph_edges.csv")
    node_file <- file.path(dir_i, "principal_graph_nodes.csv")
    if (!file.exists(cell_file) || !file.exists(edge_file) || !file.exists(node_file)) next
    cells <- as.data.frame(read_csv_if_exists(cell_file), stringsAsFactors = FALSE)
    edges <- as.data.frame(read_csv_if_exists(edge_file), stringsAsFactors = FALSE)
    nodes <- as.data.frame(read_csv_if_exists(node_file), stringsAsFactors = FALSE)
    if (nrow(cells) == 0 || nrow(edges) == 0) next
    cells[] <- lapply(cells, as.character)
    edges[] <- lapply(edges, as.character)
    nodes[] <- lapply(nodes, as.character)
    label <- basename(dir_i)
    cells$source_label <- label
    edges$source_label <- label
    nodes$source_label <- label
    cells_list[[idx]] <- cells
    edges_list[[idx]] <- edges
    nodes_list[[idx]] <- nodes
    manifest[[idx]] <- data.frame(
      source_label = label,
      source_dir = dir_i,
      n_cells = nrow(cells),
      n_edges = nrow(edges),
      stringsAsFactors = FALSE
    )
    idx <- idx + 1L
  }

  cell_all <- dplyr::bind_rows(cells_list)
  edge_all <- dplyr::bind_rows(edges_list)
  node_all <- dplyr::bind_rows(nodes_list)
  manifest_df <- dplyr::bind_rows(manifest)
  write_table_csv(manifest_df, file.path(out_dir, "Monocle3_graph_overlay_manifest.csv"))
  if (nrow(cell_all) == 0 || nrow(edge_all) == 0) return(invisible(NULL))

  cell_all <- ensure_columns(cell_all, c("UMAP_1", "UMAP_2", "pseudotime", "Dose", "Ploidy", "clusters", "source_label"))
  edge_all <- ensure_columns(edge_all, c(
    "UMAP_1_from", "UMAP_2_from", "UMAP_1_to", "UMAP_2_to", "Dose", "Ploidy",
    "from_root_distance", "to_root_distance", "source_label"
  ))
  node_all <- ensure_columns(node_all, c("graph_node", "root_graph_node", "UMAP_1", "UMAP_2", "Dose", "Ploidy", "source_label"))

  cell_all <- cell_all %>%
    dplyr::mutate(
      UMAP_1 = safe_numeric(.data$UMAP_1),
      UMAP_2 = safe_numeric(.data$UMAP_2),
      pseudotime = safe_numeric(.data$pseudotime),
      Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N"))
    ) %>%
    dplyr::filter(is.finite(.data$UMAP_1), is.finite(.data$UMAP_2), is.finite(.data$pseudotime))
  edge_all <- edge_all %>%
    dplyr::mutate(
      UMAP_1_from = safe_numeric(.data$UMAP_1_from),
      UMAP_2_from = safe_numeric(.data$UMAP_2_from),
      UMAP_1_to = safe_numeric(.data$UMAP_1_to),
      UMAP_2_to = safe_numeric(.data$UMAP_2_to),
      from_root_distance = safe_numeric(.data$from_root_distance),
      to_root_distance = safe_numeric(.data$to_root_distance),
      Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N"))
    ) %>%
    dplyr::filter(
      is.finite(.data$UMAP_1_from), is.finite(.data$UMAP_2_from),
      is.finite(.data$UMAP_1_to), is.finite(.data$UMAP_2_to)
    )
	  root_nodes <- node_all %>%
	    dplyr::mutate(
	      UMAP_1 = safe_numeric(.data$UMAP_1),
      UMAP_2 = safe_numeric(.data$UMAP_2),
      Dose = factor(as.character(.data$Dose), levels = c("0mg/kg", "30mg/kg", "120mg/kg")),
      Ploidy = factor(as.character(.data$Ploidy), levels = c("2N", "4N"))
	    ) %>%
	    dplyr::filter(.data$graph_node == .data$root_graph_node, is.finite(.data$UMAP_1), is.finite(.data$UMAP_2))
	  graph_umap_x <- c(cell_all$UMAP_1, edge_all$UMAP_1_from, edge_all$UMAP_1_to, root_nodes$UMAP_1)
	  graph_umap_y <- c(cell_all$UMAP_2, edge_all$UMAP_2_from, edge_all$UMAP_2_to, root_nodes$UMAP_2)

	  p <- ggplot(cell_all, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 0.18, alpha = 0.62) +
    geom_segment(
      data = edge_all,
      aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "black",
      alpha = 0.72,
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
    ) +
    geom_point(
      data = root_nodes,
      aes(x = UMAP_1, y = UMAP_2),
      inherit.aes = FALSE,
      shape = 8,
      size = 2.5,
      stroke = 0.8,
      color = "black"
    ) +
    facet_grid(Ploidy ~ Dose) +
    scale_color_gradientn(
	      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
	      name = "Monocle3\npseudotime"
	    ) +
	    coord_square_umap(graph_umap_x, graph_umap_y) +
	    labs(
      title = "Monocle3 pseudotime UMAP with principal graph direction",
      subtitle = "Cells are colored by Monocle3 pseudotime; arrows follow graph distance away from the root node; stars mark root nodes.",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
	    theme_extra(9) +
	    theme(axis.text.x = element_text(angle = 0)) +
	    theme_square_umap_panel()
	  save_square_panel_both(p, file.path(out_dir, "Monocle3_graph_pseudotime_by_ploidy_dose"), panel_size = 3.2)
	  save_square_panel_both(p, file.path(out_umap_dir, "Monocle3_graph_pseudotime_by_ploidy_dose"), panel_size = 3.2)
	}

message("Input roots:")
message("  trajectory_root: ", trajectory_root)
message("  pseudotrajectory_root: ", pseudotrajectory_root)
message("  deg_root: ", deg_root)
message("  deg_round2_root: ", deg_round2_root)
message("  gsea_root: ", gsea_root)
message("  output_root: ", output_root)

if (!dir.exists(trajectory_root)) stop("Missing trajectory root: ", trajectory_root, call. = FALSE)
if (!dir.exists(pseudotrajectory_root)) stop("Missing pseudotrajectory root: ", pseudotrajectory_root, call. = FALSE)
if (!dir.exists(deg_round2_root)) warning("Missing 03c DEG round2 root: ", deg_round2_root, call. = FALSE)
if (!dir.exists(gsea_root)) warning("Missing 03b GSEA root: ", gsea_root, call. = FALSE)

comparison_summary_file <- file.path(deg_root, "00_summary", "DEG_comparison_summary.csv")
if (file.exists(comparison_summary_file)) {
  comparison_summary <- readr::read_csv(comparison_summary_file, show_col_types = FALSE, progress = FALSE)
  write_table_csv(comparison_summary, file.path(out_inputs, "DEG_comparison_summary_used.csv"))
}

cluster_annotations <- if (dir.exists(gsea_root)) read_cluster_annotations(gsea_root) else data.frame()
gsea_inputs <- if (dir.exists(gsea_root)) read_gsea_inputs(gsea_root) else data.frame()
deg_marker_summary <- if (dir.exists(deg_root)) read_deg_marker_summaries(deg_root) else data.frame()
write_table_csv(cluster_annotations, file.path(out_inputs, "cluster_annotation_summary_used.csv"))
write_table_csv(gsea_inputs, file.path(out_inputs, "hallmark_gsea_tables_used.csv"))
write_table_csv(deg_marker_summary, file.path(out_inputs, "deg_marker_top_gene_summary_used.csv"))

velocity_cells <- read_velocity_metrics(trajectory_root)
monocle3_cells <- read_pseudotrajectory_metrics(pseudotrajectory_root)

write_table_csv(velocity_cells, file.path(out_tables, "scVelo_cell_metrics_long.csv"))
write_table_csv(monocle3_cells, file.path(out_tables, "Monocle3_cell_metrics_long.csv"))

velocity_suite <- run_comparison_suite(velocity_cells, velocity_metric_cols, method_scvelo)
monocle3_suite <- run_comparison_suite(monocle3_cells, pseudo_metric_cols, method_monocle3)

stat_tests <- add_adjusted_p(dplyr::bind_rows(velocity_suite$tests, monocle3_suite$tests))
group_summaries <- dplyr::bind_rows(velocity_suite$summaries, monocle3_suite$summaries)
trajectory_split_cols <- c(
  "TrajectoryComparisonID", "TrajectoryScope", "TrajectoryComparison",
  "TrajectoryIdent1", "TrajectoryIdent2", "TrajectoryGroupCol",
  "TrajectoryCluster", "TrajectorySubsetGroup", "TrajectoryDose",
  "TrajectoryComparisonGroup", "TrajectoryOriginalGroup", "TrajectoryUniverseKey"
)
stat_tests <- ensure_columns(stat_tests, trajectory_split_cols)
group_summaries <- ensure_columns(group_summaries, trajectory_split_cols)

write_table_csv(stat_tests, file.path(out_stats, "scVelo_Monocle3_stat_tests.csv"))
write_table_csv(group_summaries, file.path(out_stats, "scVelo_Monocle3_group_summaries.csv"))

question_summary <- summarize_questions(stat_tests)
write_table_csv(question_summary, file.path(out_summary, "question_level_summary.csv"))

time_metric_tests <- stat_tests %>%
  dplyr::filter(.data$metric %in% c("velocity_pseudotime", "pseudotime"))
write_table_csv(time_metric_tests, file.path(out_summary, "time_metric_tests_for_interpretation.csv"))

scvelo_separate <- write_method_separate_outputs(
  velocity_cells,
  metric_col = "velocity_pseudotime",
  method_label = method_scvelo,
  out_dir = out_scvelo,
  plot_prefix = "scVelo",
  out_umap_dir = out_umap_scvelo,
  cluster_annotations = cluster_annotations
)
monocle3_separate <- write_method_separate_outputs(
  monocle3_cells,
  metric_col = "pseudotime",
  method_label = method_monocle3,
  out_dir = out_monocle3,
  plot_prefix = "Monocle3",
  out_umap_dir = out_umap_monocle3,
  cluster_annotations = cluster_annotations
)
method_answer_summary <- dplyr::bind_rows(
  make_method_answer_summary(scvelo_separate, method_scvelo),
  make_method_answer_summary(monocle3_separate, method_monocle3)
)
write_table_csv(method_answer_summary, file.path(out_summary, "method_separate_answer_summary.csv"))

scvelo_tn_ploidy <- write_tn_ploidy_cluster_outputs(
  scvelo_separate$desc,
  metric_col = "velocity_pseudotime",
  method_label = method_scvelo,
  out_dir = out_tn_ploidy_scvelo,
  plot_prefix = "scVelo",
  cluster_annotations = cluster_annotations
)
monocle3_tn_ploidy <- write_tn_ploidy_cluster_outputs(
  monocle3_separate$desc,
  metric_col = "pseudotime",
  method_label = method_monocle3,
  out_dir = out_tn_ploidy_monocle3,
  plot_prefix = "Monocle3",
  cluster_annotations = cluster_annotations
)
write_table_csv(
  dplyr::bind_rows(scvelo_tn_ploidy$key_conclusions, monocle3_tn_ploidy$key_conclusions),
  file.path(out_tn_ploidy, "tn_ploidy_cluster_key_conclusions_all_methods.csv")
)
write_table_csv(
  dplyr::bind_rows(
    scvelo_tn_ploidy$tumor_ploidy_dose$ploidy_within_dose_sample_tests,
    monocle3_tn_ploidy$tumor_ploidy_dose$ploidy_within_dose_sample_tests
  ),
  file.path(out_tn_ploidy, "tumor_Q1_ploidy_within_dose_sample_tests_all_methods.csv")
)
write_table_csv(
  dplyr::bind_rows(
    scvelo_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_sample_tests,
    monocle3_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_sample_tests
  ),
  file.path(out_tn_ploidy, "tumor_Q2_dose_within_ploidy_sample_tests_all_methods.csv")
)
write_table_csv(
  dplyr::bind_rows(
    scvelo_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_pairwise_sample_tests,
    monocle3_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_pairwise_sample_tests
  ),
  file.path(out_tn_ploidy, "tumor_Q2_dose_pairwise_within_ploidy_sample_tests_all_methods.csv")
)
write_table_csv(
  dplyr::bind_rows(
    scvelo_tn_ploidy$tumor_ploidy_dose$feature_summary,
    monocle3_tn_ploidy$tumor_ploidy_dose$feature_summary
  ),
  file.path(out_tn_ploidy, "tumor_Q3_ploidy_dose_feature_summary_all_methods.csv")
)
scvelo_group_timing <- run_group_cluster_timing_outputs(
  velocity_cells,
  out_dir = out_cluster_timing_scvelo,
  method_label = method_scvelo,
  time_col = "velocity_pseudotime",
  source_context_filter = c("velocity_groups"),
  time_label = "velocity pseudotime",
  x_label = "scVelo velocity pseudotime"
)
monocle3_group_timing <- run_group_cluster_timing_outputs(
  monocle3_cells,
  out_dir = out_cluster_timing_monocle3,
  method_label = method_monocle3,
  time_col = "pseudotime",
  source_context_filter = c("monocle3_groups"),
  time_label = "pseudotime",
  x_label = "Monocle3 pseudotime"
)

pathway_trajectory <- join_trajectory_with_gsea(stat_tests, gsea_inputs, deg_marker_summary)
pathway_timing_summary <- summarize_pathway_timing(pathway_trajectory)
hypoxia_timing <- if (nrow(pathway_trajectory) > 0 && all(c("hypoxia_pathway", "pathway") %in% names(pathway_trajectory))) {
  pathway_trajectory %>%
    dplyr::filter(.data$hypoxia_pathway, !is.na(.data$pathway))
} else {
  data.frame()
}
hypoxia_method_summary <- if (nrow(hypoxia_timing) > 0 && "enriched_time_delta" %in% names(hypoxia_timing)) {
  hypoxia_timing %>%
    dplyr::filter(is.finite(.data$enriched_time_delta)) %>%
    dplyr::group_by(.data$method, .data$metric) %>%
    dplyr::summarise(
      n_hypoxia_comparisons = dplyr::n(),
      n_late_or_high_time = sum(.data$pathway_timing == "late_or_high_time", na.rm = TRUE),
      n_early_or_low_time = sum(.data$pathway_timing == "early_or_low_time", na.rm = TRUE),
      n_joint_gsea_pseudotime_fdr_0_05 = sum(.data$gsea_significant_fdr_0_05 & .data$trajectory_significant_fdr_0_05, na.rm = TRUE),
      median_hypoxia_time_delta = stats::median(.data$enriched_time_delta, na.rm = TRUE),
      dominant_hypoxia_timing = dplyr::case_when(
        .data$n_late_or_high_time > .data$n_early_or_low_time ~ "mostly_late_or_high_time",
        .data$n_early_or_low_time > .data$n_late_or_high_time ~ "mostly_early_or_low_time",
        TRUE ~ "mixed_or_neutral"
      ),
      .groups = "drop"
    )
} else {
  data.frame()
}
write_table_csv(pathway_trajectory, file.path(out_biology, "pseudotime_gsea_deg_joined.csv"))
write_table_csv(pathway_timing_summary, file.path(out_biology, "pathway_timing_summary_by_method.csv"))
write_table_csv(hypoxia_timing, file.path(out_biology, "hypoxia_pseudotime_timing.csv"))
write_table_csv(hypoxia_method_summary, file.path(out_biology, "hypoxia_timing_method_summary.csv"))

plot_scope_summary(stat_tests)
plot_deg_heatmap(stat_tests)
plot_cluster_vs_rest(stat_tests)
plot_global_effects(stat_tests)
plot_within_cluster_ploidy(stat_tests)
plot_omnibus_cluster(stat_tests)
plot_merged_dose_ploidy(stat_tests)
plot_gsea_trajectory_association(pathway_trajectory)
plot_hypoxia_timing(hypoxia_timing)
plot_pathway_timing_summary(pathway_timing_summary)
run_scvelo_flow_overlays(trajectory_root, out_flow_scvelo, out_umap_scvelo)
plot_monocle3_graph_overlays(pseudotrajectory_root, out_flow_monocle3, out_umap_monocle3)
if (nrow(gsea_inputs) > 0) {
  plot_cluster_hypoxia_umap(
    scvelo_separate$desc,
    gsea_inputs,
    method_label = method_scvelo,
    metric_col = "velocity_pseudotime",
    out_dir = out_scvelo,
    plot_prefix = "scVelo",
    out_umap_dir = out_umap_scvelo
  )
  plot_cluster_hypoxia_umap(
    monocle3_separate$desc,
    gsea_inputs,
    method_label = method_monocle3,
    metric_col = "pseudotime",
    out_dir = out_monocle3,
    plot_prefix = "Monocle3",
    out_umap_dir = out_umap_monocle3
  )
}

paga_extra <- run_paga_extra_analysis(
  trajectory_root = trajectory_root,
  deg_round2_root = deg_round2_root,
  output_root = out_paga_extra,
  analysis_label = "04b PAGA",
  top_n = 10L,
  max_cascade_genes = 30L
)

final_tumor_summary <- write_final_tumor_summary_outputs(
  output_root,
  out_final_summary,
  scvelo_tn_ploidy,
  monocle3_tn_ploidy
)

run_summary_lines <- c(
  "04b extra scVelo/Monocle3 pseudotime analysis completed.",
  paste0("Output root: ", normalizePath(output_root, mustWork = FALSE)),
  paste0("scVelo cell rows: ", nrow(velocity_cells)),
  paste0("Monocle3 cell rows: ", nrow(monocle3_cells)),
  paste0("Stat test rows: ", nrow(stat_tests)),
  paste0("Group summary rows: ", nrow(group_summaries)),
  paste0("GSEA pseudotime rows: ", nrow(pathway_trajectory)),
  paste0("Hypoxia timing rows: ", nrow(hypoxia_timing)),
  paste0("TN/Ploidy scVelo key conclusion rows: ", nrow(scvelo_tn_ploidy$key_conclusions)),
  paste0("TN/Ploidy Monocle3 key conclusion rows: ", nrow(monocle3_tn_ploidy$key_conclusions)),
  paste0("scVelo group cluster timing UMAP rows: ", nrow(scvelo_group_timing$status)),
  paste0("scVelo group cluster timing summary rows: ", nrow(scvelo_group_timing$cluster_summary)),
  paste0("scVelo group sample histogram rows: ", nrow(scvelo_group_timing$histogram_summary)),
  paste0("Monocle3 group cluster timing UMAP rows: ", nrow(monocle3_group_timing$status)),
  paste0("Monocle3 group cluster timing summary rows: ", nrow(monocle3_group_timing$cluster_summary)),
  paste0("Monocle3 group sample histogram rows: ", nrow(monocle3_group_timing$histogram_summary)),
  paste0("PAGA extra groups processed: ", if (is.data.frame(paga_extra$status)) nrow(paga_extra$status) else 0L),
  paste0("PAGA extra completed groups: ", if (is.data.frame(paga_extra$status) && nrow(paga_extra$status) > 0) sum(paga_extra$status$status == "completed", na.rm = TRUE) else 0L),
  paste0("PAGA extra partial groups: ", if (is.data.frame(paga_extra$status) && nrow(paga_extra$status) > 0) sum(paga_extra$status$status == "partial", na.rm = TRUE) else 0L),
  paste0("Final summary copied supporting files: ", sum(final_tumor_summary$copy_manifest$copied, na.rm = TRUE)),
  "",
  "Main time metrics:",
  "  scVelo: velocity_pseudotime",
  "  Monocle3: pseudotime",
  "",
  "Interpretation note:",
  "  Cell-level tests compare cell distributions and can be highly powered because cells are not biological replicates.",
  "  Sample-level tests summarize each sample by median first and are included as a more conservative companion table.",
  "",
  "Key outputs:",
  "  01_tables/scVelo_cell_metrics_long.csv",
  "  01_tables/Monocle3_cell_metrics_long.csv",
  "  02_stats/scVelo_Monocle3_stat_tests.csv",
  "  02_stats/scVelo_Monocle3_group_summaries.csv",
  "  04_question_summary/question_level_summary.csv",
  "  04_question_summary/method_separate_answer_summary.csv",
  "  05_method_separate/scVelo/*",
  "  05_method_separate/Monocle3/*",
  "  06_deg_gsea_pseudotime/*",
  "  07_umap_overlays/*",
  "  08_flow_overlays/*",
  "  09_tn_ploidy_cluster_response/*",
  "  10_cluster_timing/scVelo/*",
  "  10_cluster_timing/Monocle3/*",
  "  11_paga_extra/*",
  "  summary/*",
  "  03_plots/*.pdf and *.png"
)
writeLines(run_summary_lines, file.path(output_root, "run_summary.txt"))

message(paste(run_summary_lines, collapse = "\n"))
