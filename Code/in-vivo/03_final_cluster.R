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
  library(ggplot2)
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
  "02d_manual_cluster_merge",
  "objects",
  "integrated_sct_cca_seurat_cluster_refine_manual_merge_test.rds"
)
output_root <- file.path(results_root, "03_final_cluster")
final_output_rds <- file.path(
  output_root,
  "03_objects",
  "integrated_sct_cca_seurat_final_reclustered.rds"
)
plot_only_mode <- tolower(Sys.getenv("FINAL_CLUSTER_PLOT_ONLY", unset = "false")) %in%
  c("1", "true", "yes", "y")
plot_only_input_rds <- Sys.getenv("FINAL_CLUSTER_PLOT_INPUT_RDS", unset = final_output_rds)
plot_only_output_root <- Sys.getenv(
  "FINAL_CLUSTER_PLOT_OUTPUT_ROOT",
  unset = dirname(dirname(normalizePath(plot_only_input_rds, mustWork = FALSE)))
)

filter_cluster_col <- "manual_merge_test"
cluster_col <- "clusters"
sample_col <- "sample"
remove_cells_from_clusters <- c("9", "4", "3", "9c")

preferred_assays_for_pca <- c("integrated", "SCT", "RNA")
npcs <- 30L
pca_npcs_requested <- 50L
fallback_variable_features <- 3000L
future_globals_maxsize_gb <- 60

safe_file_component <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "value")
}

save_plot_pdf_png_tiff <- function(plot_obj, file_stub, width = 9, height = 7, dpi = 300) {
  ggplot2::ggsave(paste0(file_stub, ".pdf"), plot_obj, width = width, height = height)
  ggplot2::ggsave(paste0(file_stub, ".png"), plot_obj, width = width, height = height, dpi = dpi)
  ggplot2::ggsave(
    paste0(file_stub, ".tiff"),
    plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    device = "tiff",
    compression = "lzw"
  )
}

write_numeric_score_umap <- function(
  obj,
  output_dir,
  score_values,
  score_label,
  file_stub,
  title,
  width = 9,
  height = 7
) {
  if (!("umap" %in% names(obj@reductions))) {
    warning("Skipping ", score_label, " UMAP because the Seurat object has no 'umap' reduction.", call. = FALSE)
    return(FALSE)
  }

  umap_mat <- Embeddings(obj, reduction = "umap")
  if (ncol(umap_mat) < 2L) {
    warning("Skipping ", score_label, " UMAP because the 'umap' reduction has fewer than two dimensions.", call. = FALSE)
    return(FALSE)
  }
  if (length(score_values) != nrow(obj@meta.data)) {
    warning("Skipping ", score_label, " UMAP because score length does not match metadata rows.", call. = FALSE)
    return(FALSE)
  }

  umap_df <- as.data.frame(umap_mat[, seq_len(2), drop = FALSE])
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  names(score_values) <- rownames(obj@meta.data)
  umap_df[["score_value"]] <- suppressWarnings(as.numeric(score_values[rownames(umap_df)]))
  if (all(is.na(umap_df[["score_value"]]))) {
    warning("Skipping ", score_label, " UMAP because all score values are NA or non-numeric.", call. = FALSE)
    return(FALSE)
  }

  umap_df <- umap_df[order(is.na(umap_df[["score_value"]]), umap_df[["score_value"]]), , drop = FALSE]
  p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = score_value)) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
      na.value = "grey85",
      name = score_label
    ) +
    coord_equal() +
    labs(
      title = title,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))

  save_plot_pdf_png_tiff(
    plot_obj = p,
    file_stub = file.path(output_dir, file_stub),
    width = width,
    height = height
  )
  TRUE
}

write_s_score_umap <- function(obj, output_dir, file_stub = "umap_by_S_Score", width = 9, height = 7) {
  if (!("S.Score" %in% colnames(obj@meta.data))) {
    warning("Skipping S.Score UMAP because metadata column 'S.Score' is missing.", call. = FALSE)
    return(FALSE)
  }
  write_numeric_score_umap(
    obj = obj,
    output_dir = output_dir,
    score_values = suppressWarnings(as.numeric(obj@meta.data[["S.Score"]])),
    score_label = "S.Score",
    file_stub = file_stub,
    title = "UMAP colored by S phase score",
    width = width,
    height = height
  )
}

write_g1_score_umap <- function(obj, output_dir, file_stub = "umap_by_G1_Score", width = 9, height = 7) {
  if ("G1.Score" %in% colnames(obj@meta.data)) {
    g1_score <- suppressWarnings(as.numeric(obj@meta.data[["G1.Score"]]))
    score_label <- "G1.Score"
    title <- "UMAP colored by G1 score"
  } else if (all(c("S.Score", "G2M.Score") %in% colnames(obj@meta.data))) {
    s_score <- suppressWarnings(as.numeric(obj@meta.data[["S.Score"]]))
    g2m_score <- suppressWarnings(as.numeric(obj@meta.data[["G2M.Score"]]))
    cycling_score <- pmax(s_score, g2m_score, na.rm = TRUE)
    cycling_score[is.infinite(cycling_score)] <- NA_real_
    g1_score <- -cycling_score
    score_label <- "G1.Score (derived)"
    title <- "UMAP colored by derived G1 score"
  } else {
    warning(
      "Skipping G1.Score UMAP because metadata column 'G1.Score' is missing and S.Score/G2M.Score are not both available.",
      call. = FALSE
    )
    return(FALSE)
  }

  write_numeric_score_umap(
    obj = obj,
    output_dir = output_dir,
    score_values = g1_score,
    score_label = score_label,
    file_stub = file_stub,
    title = title,
    width = width,
    height = height
  )
}

choose_pca_assay <- function(obj, preferred_assays) {
  hit <- preferred_assays[preferred_assays %in% names(obj@assays)]
  if (length(hit) == 0) {
    stop(
      "None of the preferred PCA assays are present: ",
      paste(preferred_assays, collapse = ", "),
      call. = FALSE
    )
  }
  hit[1]
}

prepare_assay_for_pca <- function(obj, assay, fallback_nfeatures) {
  DefaultAssay(obj) <- assay
  if (assay == "RNA") {
    obj <- maybe_join_layers(obj, assay = assay)
    data_mat <- get_assay_data_slot(obj, assay = assay, slot_name = "data")
    if (is.null(data_mat) || nrow(data_mat) == 0 || ncol(data_mat) == 0) {
      obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
    }
    if (length(VariableFeatures(obj)) < 2L) {
      obj <- FindVariableFeatures(
        obj,
        assay = assay,
        nfeatures = fallback_nfeatures,
        verbose = FALSE
      )
    }
    obj <- ScaleData(
      obj,
      assay = assay,
      features = VariableFeatures(obj),
      verbose = FALSE
    )
    return(obj)
  }

  if (length(VariableFeatures(obj)) < 2L) {
    scale_mat <- get_assay_data_slot(obj, assay = assay, slot_name = "scale.data")
    data_mat <- get_assay_data_slot(obj, assay = assay, slot_name = "data")
    candidate_features <- if (!is.null(scale_mat) && nrow(scale_mat) > 1L) {
      rownames(scale_mat)
    } else if (!is.null(data_mat) && nrow(data_mat) > 1L) {
      rownames(data_mat)
    } else {
      character(0)
    }

    if (length(candidate_features) < 2L) {
      stop("Cannot find enough features in assay for PCA: ", assay, call. = FALSE)
    }
    VariableFeatures(obj) <- head(candidate_features, min(length(candidate_features), fallback_nfeatures))
  }

  scale_mat <- get_assay_data_slot(obj, assay = assay, slot_name = "scale.data")
  if (is.null(scale_mat) || nrow(scale_mat) == 0 || ncol(scale_mat) == 0) {
    obj <- ScaleData(
      obj,
      assay = assay,
      features = VariableFeatures(obj),
      verbose = FALSE
    )
  }
  obj
}

make_observed_levels <- function(values, source_col = NULL, preferred_levels = NULL) {
  values <- as.character(values)
  values[is.na(values) | values == ""] <- "Unknown"
  observed <- unique(values)

  candidate_levels <- if (!is.null(preferred_levels)) {
    as.character(preferred_levels)
  } else if (!is.null(source_col) && is.factor(source_col)) {
    as.character(levels(source_col))
  } else {
    sort_maybe_numeric(values)
  }

  c(
    intersect(candidate_levels, observed),
    setdiff(sort_maybe_numeric(observed), candidate_levels)
  )
}

add_derived_tn_ploidy_metadata <- function(meta, overwrite = FALSE) {
  needs_ploidy <- overwrite || !"Ploidy" %in% colnames(meta)
  needs_tn <- overwrite || !"TN" %in% colnames(meta)
  if (!needs_ploidy && !needs_tn) return(meta)

  assert_required_column(meta, "IDs")
  ids_values <- as.character(meta[["IDs"]])
  ids_values[is.na(ids_values)] <- ""

  if (needs_ploidy) {
    ploidy_values <- ifelse(
      grepl("2N", ids_values, fixed = TRUE),
      "2N",
      ifelse(grepl("4N", ids_values, fixed = TRUE), "4N", NA_character_)
    )
    if (any(is.na(ploidy_values))) {
      unresolved_ploidy <- sort(unique(ids_values[is.na(ploidy_values)]))
      stop(
        "Ploidy must be either 2N or 4N, but some IDs cannot be parsed: ",
        paste(utils::head(unresolved_ploidy, 20), collapse = ", "),
        call. = FALSE
      )
    }
    meta[["Ploidy"]] <- factor(ploidy_values, levels = c("2N", "4N"))
  }

  if (needs_tn) {
    tn_values <- ifelse(
      grepl("Cell-Culture", ids_values, fixed = TRUE),
      "CellLine",
      "Tumor"
    )
    meta[["TN"]] <- factor(tn_values, levels = c("CellLine", "Tumor"))
  }

  meta
}

write_cluster_stack_plots <- function(
  meta,
  group_col,
  cluster_col,
  output_dir,
  summary_dir,
  file_stub,
  group_label = group_col,
  cluster_label = cluster_col,
  preferred_group_levels = NULL,
  width = 10,
  height = 6,
  dpi = 300
) {
  assert_required_column(meta, group_col)
  assert_required_column(meta, cluster_col)

  group_vals <- as.character(meta[[group_col]])
  cluster_vals <- as.character(meta[[cluster_col]])
  group_vals[is.na(group_vals) | group_vals == ""] <- "Unknown"
  cluster_vals[is.na(cluster_vals) | cluster_vals == ""] <- "Unknown"

  group_levels <- make_observed_levels(
    group_vals,
    source_col = meta[[group_col]],
    preferred_levels = preferred_group_levels
  )
  cluster_levels <- make_observed_levels(
    cluster_vals,
    source_col = meta[[cluster_col]]
  )

  plot_df <- data.frame(
    group = factor(group_vals, levels = group_levels),
    cluster = factor(cluster_vals, levels = cluster_levels),
    stringsAsFactors = FALSE
  )

  count_df <- as.data.frame(
    table(group = plot_df$group, cluster = plot_df$cluster),
    stringsAsFactors = FALSE
  )
  names(count_df)[names(count_df) == "Freq"] <- "count"
  count_df$group <- factor(count_df$group, levels = group_levels)
  count_df$cluster <- factor(count_df$cluster, levels = cluster_levels)
  count_df <- count_df %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(
      total = sum(count),
      proportion = ifelse(total > 0, count / total, 0),
      percent = 100 * proportion
    ) %>%
    dplyr::ungroup()

  count_export <- count_df
  names(count_export)[names(count_export) == "group"] <- group_col
  names(count_export)[names(count_export) == "cluster"] <- cluster_col
  write_table_csv(count_export, file.path(summary_dir, paste0(file_stub, "_counts.csv")))

  percent_export <- count_export
  write_table_csv(percent_export, file.path(summary_dir, paste0(file_stub, "_percent.csv")))

  p_count <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = group, y = count, fill = cluster)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = paste0(cluster_label, " cell counts by ", group_label),
      x = group_label,
      y = "Cell count",
      fill = cluster_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p_percent <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = group, y = proportion, fill = cluster)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0(cluster_label, " percentage by ", group_label),
      x = group_label,
      y = "Percentage",
      fill = cluster_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  save_plot_pdf_png_tiff(
    plot_obj = p_count,
    file_stub = file.path(output_dir, paste0(file_stub, "_count")),
    width = width,
    height = height,
    dpi = dpi
  )
  save_plot_pdf_png_tiff(
    plot_obj = p_percent,
    file_stub = file.path(output_dir, paste0(file_stub, "_percent")),
    width = width,
    height = height,
    dpi = dpi
  )

  invisible(count_df)
}

write_cluster_composition_stack_plots <- function(
  meta,
  group_col,
  cluster_col,
  output_dir,
  summary_dir,
  file_stub,
  group_label = group_col,
  cluster_label = cluster_col,
  preferred_group_levels = NULL,
  fill_values = NULL,
  width = 10,
  height = 6,
  dpi = 300
) {
  assert_required_column(meta, group_col)
  assert_required_column(meta, cluster_col)

  group_vals <- as.character(meta[[group_col]])
  cluster_vals <- as.character(meta[[cluster_col]])
  group_vals[is.na(group_vals) | group_vals == ""] <- "Unknown"
  cluster_vals[is.na(cluster_vals) | cluster_vals == ""] <- "Unknown"

  group_levels <- make_observed_levels(
    group_vals,
    source_col = meta[[group_col]],
    preferred_levels = preferred_group_levels
  )
  cluster_levels <- make_observed_levels(
    cluster_vals,
    source_col = meta[[cluster_col]]
  )

  plot_df <- data.frame(
    cluster = factor(cluster_vals, levels = cluster_levels),
    group = factor(group_vals, levels = group_levels),
    stringsAsFactors = FALSE
  )

  count_df <- as.data.frame(
    table(cluster = plot_df$cluster, group = plot_df$group),
    stringsAsFactors = FALSE
  )
  names(count_df)[names(count_df) == "Freq"] <- "count"
  count_df$cluster <- factor(count_df$cluster, levels = cluster_levels)
  count_df$group <- factor(count_df$group, levels = group_levels)
  count_df <- count_df %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      total = sum(count),
      proportion = ifelse(total > 0, count / total, 0),
      percent = 100 * proportion
    ) %>%
    dplyr::ungroup()

  count_export <- count_df
  names(count_export)[names(count_export) == "cluster"] <- cluster_col
  names(count_export)[names(count_export) == "group"] <- group_col
  write_table_csv(count_export, file.path(summary_dir, paste0(file_stub, "_counts.csv")))
  percent_export <- count_export
  write_table_csv(percent_export, file.path(summary_dir, paste0(file_stub, "_percent.csv")))

  p_count <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = cluster, y = count, fill = group)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(
      title = paste0(group_label, " cell counts by ", cluster_label),
      x = cluster_label,
      y = "Cell count",
      fill = group_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  p_percent <- ggplot2::ggplot(
    count_df,
    ggplot2::aes(x = cluster, y = proportion, fill = group)
  ) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(round(100 * x), "%"),
      breaks = seq(0, 1, 0.25),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0(group_label, " percentage by ", cluster_label),
      x = cluster_label,
      y = "Percentage",
      fill = group_label
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(fill_values)) {
    fill_values <- fill_values[names(fill_values) %in% group_levels]
    if (length(fill_values) > 0) {
      p_count <- p_count + ggplot2::scale_fill_manual(values = fill_values, drop = FALSE)
      p_percent <- p_percent + ggplot2::scale_fill_manual(values = fill_values, drop = FALSE)
    }
  }

  save_plot_pdf_png_tiff(
    plot_obj = p_count,
    file_stub = file.path(output_dir, paste0(file_stub, "_count")),
    width = width,
    height = height,
    dpi = dpi
  )
  save_plot_pdf_png_tiff(
    plot_obj = p_percent,
    file_stub = file.path(output_dir, paste0(file_stub, "_percent")),
    width = width,
    height = height,
    dpi = dpi
  )

  invisible(count_df)
}

write_requested_cluster_stack_plots <- function(meta, output_dir, summary_dir, cluster_col, sample_col) {
  meta <- add_derived_tn_ploidy_metadata(meta)
  assert_required_column(meta, cluster_col)
  assert_required_column(meta, sample_col)
  assert_required_column(meta, "TN")
  assert_required_column(meta, "Ploidy")

  cluster_n <- length(make_observed_levels(as.character(meta[[cluster_col]]), source_col = meta[[cluster_col]]))
  sample_n <- length(make_observed_levels(as.character(meta[[sample_col]]), source_col = meta[[sample_col]]))
  cluster_stack_width <- max(9, 0.45 * cluster_n + 5)
  sample_stack_height <- max(6, 0.18 * sample_n + 5)

  write_cluster_composition_stack_plots(
    meta = meta,
    group_col = "TN",
    cluster_col = cluster_col,
    output_dir = output_dir,
    summary_dir = summary_dir,
    file_stub = "stack_TN_by_cluster",
    group_label = "CellLine/Tumor",
    cluster_label = "cluster",
    preferred_group_levels = c("CellLine", "Tumor"),
    fill_values = c(CellLine = "#4E79A7", Tumor = "#E15759"),
    width = cluster_stack_width,
    height = 6
  )
  write_cluster_composition_stack_plots(
    meta = meta,
    group_col = sample_col,
    cluster_col = cluster_col,
    output_dir = output_dir,
    summary_dir = summary_dir,
    file_stub = "stack_sample_by_cluster",
    group_label = "sample",
    cluster_label = "cluster",
    width = cluster_stack_width,
    height = sample_stack_height
  )
  write_cluster_composition_stack_plots(
    meta = meta,
    group_col = "Ploidy",
    cluster_col = cluster_col,
    output_dir = output_dir,
    summary_dir = summary_dir,
    file_stub = "stack_Ploidy_by_cluster",
    group_label = "Ploidy",
    cluster_label = "cluster",
    preferred_group_levels = c("2N", "4N"),
    fill_values = c(`2N` = "#59A14F", `4N` = "#B07AA1"),
    width = cluster_stack_width,
    height = 6
  )

  invisible(meta)
}

run_final_cluster_plot_only <- function(input_path, output_root, cluster_col, sample_col) {
  message("Plot-only mode: reading final Seurat object: ", input_path)
  if (!file.exists(input_path)) {
    stop("Plot-only input Seurat object does not exist: ", input_path, call. = FALSE)
  }

  .ensure_dir(output_root)
  out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
  out_plots <- .ensure_dir(file.path(output_root, "02_plots"))

  obj <- readRDS(input_path)
  if (!inherits(obj, "Seurat")) {
    stop("Plot-only input file is not a Seurat object.", call. = FALSE)
  }

  meta <- write_requested_cluster_stack_plots(
    meta = obj@meta.data,
    output_dir = out_plots,
    summary_dir = out_summary,
    cluster_col = cluster_col,
    sample_col = sample_col
  )
  s_score_umap_written <- write_s_score_umap(obj, out_plots)
  g1_score_umap_written <- write_g1_score_umap(obj, out_plots)

  plot_only_summary <- c(
    "03 final cluster plot-only completed.",
    paste0("Input object: ", input_path),
    paste0("Output root: ", output_root),
    paste0("Output plots: ", out_plots),
    paste0("Cluster column: ", cluster_col),
    paste0("Sample column: ", sample_col),
    paste0("Cells plotted: ", nrow(meta)),
    paste0("S.Score UMAP written: ", s_score_umap_written),
    paste0("G1.Score UMAP written: ", g1_score_umap_written),
    "New cluster-x stacked bar plots written for TN, sample, and Ploidy as count and percent views.",
    "Integration workflow: skipped",
    "PCA/UMAP workflow: skipped",
    "FindNeighbors: skipped",
    "FindClusters: skipped"
  )
  writeLines(plot_only_summary, con = file.path(out_summary, "plot_only_cluster_stack_summary.txt"))

  message("Plot-only cluster stack plots finished.")
  message("Output plots: ", out_plots)
}

if (plot_only_mode) {
  run_final_cluster_plot_only(
    input_path = plot_only_input_rds,
    output_root = plot_only_output_root,
    cluster_col = cluster_col,
    sample_col = sample_col
  )
  quit(save = "no", status = 0, runLast = FALSE)
}

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_plots <- .ensure_dir(file.path(output_root, "02_plots"))
out_objects <- .ensure_dir(file.path(output_root, "03_objects"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}

configure_future_for_seurat(max_size_gb = future_globals_maxsize_gb)

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}

assert_required_column(obj@meta.data, filter_cluster_col)
assert_required_column(obj@meta.data, sample_col)
assert_required_column(obj@meta.data, "IDs")

filter_clusters <- as.character(obj@meta.data[[filter_cluster_col]])
if (any(is.na(filter_clusters) | filter_clusters == "")) {
  stop("Missing values detected in ", filter_cluster_col, ".", call. = FALSE)
}

filter_cluster_levels <- if (is.factor(obj@meta.data[[filter_cluster_col]])) {
  as.character(levels(obj@meta.data[[filter_cluster_col]]))
} else {
  sort_maybe_numeric(filter_clusters)
}

missing_remove_clusters <- setdiff(remove_cells_from_clusters, unique(filter_clusters))
if (length(missing_remove_clusters) > 0) {
  stop(
    "Configured clusters to remove were not found in ",
    filter_cluster_col,
    ": ",
    paste(missing_remove_clusters, collapse = ", "),
    call. = FALSE
  )
}

keep_flag <- !(filter_clusters %in% remove_cells_from_clusters)
keep_cells <- rownames(obj@meta.data)[keep_flag]
removed_cells <- rownames(obj@meta.data)[!keep_flag]
retained_filter_levels <- setdiff(filter_cluster_levels, remove_cells_from_clusters)

filter_summary_df <- data.frame(
  metric = c("input_cells", "removed_cells", "retained_cells", "removed_cluster_n", "retained_cluster_n"),
  value = c(
    ncol(obj),
    length(removed_cells),
    length(keep_cells),
    length(remove_cells_from_clusters),
    length(retained_filter_levels)
  ),
  stringsAsFactors = FALSE
)
write_table_csv(filter_summary_df, file.path(out_summary, "filter_summary.csv"))

cluster_counts_before <- data.frame(
  cluster = names(table(filter_clusters)),
  n_cells = as.integer(table(filter_clusters)),
  status = ifelse(names(table(filter_clusters)) %in% remove_cells_from_clusters, "removed", "retained"),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, filter_cluster_levels))
write_table_csv(cluster_counts_before, file.path(out_summary, "manual_merge_test_counts_before_filter.csv"))

removed_cell_df <- obj@meta.data[!keep_flag, , drop = FALSE] %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::mutate(
    !!filter_cluster_col := as.character(.data[[filter_cluster_col]])
  ) %>%
  dplyr::select(dplyr::any_of(c("cell", filter_cluster_col, "seurat_clusters", sample_col, "Dose", "IDs")))
write_table_csv(removed_cell_df, file.path(out_summary, "removed_cells.csv"))

message("Subsetting retained cells and renaming manual_merge_test to clusters.")
obj_filtered <- subset(obj, cells = keep_cells)
if (cluster_col %in% colnames(obj_filtered@meta.data) && cluster_col != filter_cluster_col) {
  obj_filtered@meta.data[[paste0(cluster_col, "_orig")]] <- obj_filtered@meta.data[[cluster_col]]
}
obj_filtered@meta.data[[cluster_col]] <- factor(
  as.character(obj_filtered@meta.data[[filter_cluster_col]]),
  levels = retained_filter_levels
)
if (filter_cluster_col != cluster_col) {
  obj_filtered@meta.data[[filter_cluster_col]] <- NULL
}
Idents(obj_filtered) <- obj_filtered@meta.data[[cluster_col]]

cluster_counts_after_filter <- data.frame(
  cluster = names(table(obj_filtered@meta.data[[cluster_col]])),
  n_cells = as.integer(table(obj_filtered@meta.data[[cluster_col]])),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(match(cluster, retained_filter_levels))
write_table_csv(cluster_counts_after_filter, file.path(out_summary, "clusters_counts_after_filter.csv"))

message("Adding Ploidy and TN metadata from IDs.")
obj_filtered@meta.data <- add_derived_tn_ploidy_metadata(obj_filtered@meta.data, overwrite = TRUE)

metadata_summary_df <- data.frame(
  metadata = c("Ploidy", "Ploidy", "TN", "TN"),
  value = c("2N", "4N", "CellLine", "Tumor"),
  n_cells = c(
    sum(obj_filtered@meta.data[["Ploidy"]] == "2N", na.rm = TRUE),
    sum(obj_filtered@meta.data[["Ploidy"]] == "4N", na.rm = TRUE),
    sum(obj_filtered@meta.data[["TN"]] == "CellLine", na.rm = TRUE),
    sum(obj_filtered@meta.data[["TN"]] == "Tumor", na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
write_table_csv(metadata_summary_df, file.path(out_summary, "derived_metadata_counts.csv"))

pca_assay <- choose_pca_assay(obj_filtered, preferred_assays_for_pca)
message("Preparing assay for PCA: ", pca_assay)
obj_filtered <- prepare_assay_for_pca(
  obj_filtered,
  assay = pca_assay,
  fallback_nfeatures = fallback_variable_features
)
DefaultAssay(obj_filtered) <- pca_assay

pca_features <- VariableFeatures(obj_filtered)
if (length(pca_features) < 2L) {
  stop("Fewer than 2 variable features are available for PCA.", call. = FALSE)
}

pca_npcs_compute <- min(
  as.integer(pca_npcs_requested),
  as.integer(length(pca_features) - 1L),
  as.integer(ncol(obj_filtered) - 1L)
)
if (pca_npcs_compute < 2L) {
  stop("Too few cells/features remain after filtering to run PCA.", call. = FALSE)
}

message("Running PCA and UMAP only; integration and clustering are skipped.")
obj_filtered <- RunPCA(
  obj_filtered,
  assay = pca_assay,
  features = pca_features,
  npcs = pca_npcs_compute,
  verbose = FALSE
)
pca_available <- ncol(Embeddings(obj_filtered, reduction = "pca"))
dims_use <- seq_len(min(as.integer(npcs), pca_available))
obj_filtered <- RunUMAP(
  obj_filtered,
  reduction = "pca",
  dims = dims_use,
  reduction.name = "umap",
  reduction.key = "UMAP_",
  verbose = FALSE
)

pca_umap_parameters_df <- data.frame(
  parameter = c(
    "input_rds",
    "output_root",
    "filter_cluster_col",
    "cluster_col",
    "remove_cells_from_clusters",
    "preferred_assays_for_pca",
    "pca_assay_used",
    "variable_features_used",
    "pca_npcs_requested",
    "pca_npcs_computed",
    "umap_dims_used",
    "reintegration_run",
    "find_neighbors_run",
    "find_clusters_run",
    "cluster_labels_changed"
  ),
  value = c(
    input_rds,
    output_root,
    filter_cluster_col,
    cluster_col,
    paste(remove_cells_from_clusters, collapse = ","),
    paste(preferred_assays_for_pca, collapse = ","),
    pca_assay,
    as.character(length(pca_features)),
    as.character(pca_npcs_requested),
    as.character(pca_npcs_compute),
    paste(dims_use, collapse = ","),
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE"
  ),
  stringsAsFactors = FALSE
)
write_table_csv(pca_umap_parameters_df, file.path(out_summary, "pca_umap_parameters.csv"))

plot_group_cols <- unique(c(
  cluster_col,
  sample_col,
  "Dose",
  "IDs"
))
plot_group_cols <- plot_group_cols[plot_group_cols %in% colnames(obj_filtered@meta.data)]

for (group_col in plot_group_cols) {
  label_plot <- group_col == cluster_col
  p <- DimPlot(
    obj_filtered,
    reduction = "umap",
    group.by = group_col,
    label = label_plot,
    repel = TRUE,
    pt.size = 0.30,
    raster = FALSE
  ) +
    labs(title = paste0("UMAP by ", group_col))

  save_plot_pdf_png_tiff(
    plot_obj = p,
    file_stub = file.path(out_plots, paste0("umap_by_", safe_file_component(group_col))),
    width = 9,
    height = 7
  )
}
s_score_umap_written <- write_s_score_umap(obj_filtered, out_plots)
g1_score_umap_written <- write_g1_score_umap(obj_filtered, out_plots)

message("Writing cluster stacked bar plots.")
sample_stack_width <- max(10, 0.45 * length(unique(as.character(obj_filtered@meta.data[[sample_col]]))) + 5)
write_cluster_stack_plots(
  meta = obj_filtered@meta.data,
  group_col = sample_col,
  cluster_col = cluster_col,
  output_dir = out_plots,
  summary_dir = out_summary,
  file_stub = "stack_clusters_by_sample",
  group_label = "sample",
  cluster_label = "clusters",
  width = sample_stack_width,
  height = 6
)
write_cluster_stack_plots(
  meta = obj_filtered@meta.data,
  group_col = "Ploidy",
  cluster_col = cluster_col,
  output_dir = out_plots,
  summary_dir = out_summary,
  file_stub = "stack_clusters_by_Ploidy",
  group_label = "Ploidy",
  cluster_label = "clusters",
  preferred_group_levels = c("2N", "4N"),
  width = 7,
  height = 6
)
write_cluster_stack_plots(
  meta = obj_filtered@meta.data,
  group_col = "TN",
  cluster_col = cluster_col,
  output_dir = out_plots,
  summary_dir = out_summary,
  file_stub = "stack_clusters_by_TN",
  group_label = "TN",
  cluster_label = "clusters",
  preferred_group_levels = c("CellLine", "Tumor"),
  width = 7,
  height = 6
)
write_requested_cluster_stack_plots(
  meta = obj_filtered@meta.data,
  output_dir = out_plots,
  summary_dir = out_summary,
  cluster_col = cluster_col,
  sample_col = sample_col
)

message("Saving final Seurat object.")
output_rds <- final_output_rds
saveRDS(obj_filtered, output_rds)

summary_lines <- c(
  "03 final cluster completed.",
  paste0("Input object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Output object: ", output_rds),
  paste0("Filter cluster column: ", filter_cluster_col),
  paste0("Output cluster column: ", cluster_col),
  paste0("Removed cells from clusters: ", paste(remove_cells_from_clusters, collapse = ", ")),
  paste0("Cells removed: ", length(removed_cells)),
  paste0("Cells retained: ", length(keep_cells)),
  paste0("Retained ", cluster_col, " labels: ", paste(retained_filter_levels, collapse = ", ")),
  paste0("PCA assay: ", pca_assay),
  paste0("PCA components computed: ", pca_npcs_compute),
  paste0("UMAP dims used: ", paste(dims_use, collapse = ",")),
  paste0("S.Score UMAP written: ", s_score_umap_written),
  paste0("G1.Score UMAP written: ", g1_score_umap_written),
  "Derived metadata columns added from IDs: Ploidy and TN.",
  "Cluster stacked bar plots written for sample, Ploidy, and TN as count and percent views.",
  "Cluster-x stacked bar plots written for TN, sample, and Ploidy as count and percent views.",
  "Plots are written as PDF, PNG, and TIFF.",
  "Integration workflow: skipped",
  "FindNeighbors: skipped",
  "FindClusters: skipped",
  "Cluster labels were kept unchanged; manual_merge_test was renamed to clusters after filtering.",
  "The output object filename is kept for downstream compatibility, although no reclustering was performed."
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("03_final_cluster finished.")
message("Output root: ", output_root)
