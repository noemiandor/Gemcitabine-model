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
  "dplyr",
  "fgsea",
  "ggplot2",
  "msigdbr",
  "pheatmap",
  "readr",
  "tibble",
  "tidyr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(ggplot2)
  library(msigdbr)
  library(pheatmap)
  library(readr)
  library(tibble)
  library(tidyr)
})

set.seed(12345)

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

get_env_positive_integer <- function(name, default) {
  value <- get_env_scalar(name, as.character(default))
  value <- suppressWarnings(as.integer(value))
  if (is.na(value) || value < 1L) {
    stop("Environment variable ", name, " must be a positive integer when set.", call. = FALSE)
  }
  value
}

input_deg_root <- get_env_scalar("CLUSTER_GSEA_INPUT_DEG_ROOT", file.path(results_root, "03a_DEGs"))
output_root <- get_env_scalar("CLUSTER_GSEA_OUTPUT_ROOT", file.path(results_root, "03b_cluster_annotation_and_GSEA"))

comparison_summary_file <- file.path(input_deg_root, "00_summary", "DEG_comparison_summary.csv")
metadata_counts_file <- file.path(input_deg_root, "00_summary", "metadata_cell_counts.csv")

cluster_scope <- "clusters_vs_rest"
analysis_cluster_label <- "clusters"

deg_padj_cutoff <- 0.05
deg_lfc_cutoff_for_ora <- 0.25
deg_abs_delta_pct_cutoff_for_ora <- 0.05
deg_top_n_for_ora <- 100L
ora_annotation_fdr_cutoff <- 0.05
annotation_weight_expression <- 1 / 3
annotation_weight_detection <- 1 / 3
annotation_weight_overlap <- 1 / 3
min_hallmark_size <- 15L
max_hallmark_size <- 500L
ora_min_overlap <- 3L
top_terms_per_cluster <- 3L
top_plot_terms <- 20L
gsea_min_size <- 15L
gsea_max_size <- 500L
gsea_nperm_simple <- get_env_positive_integer("CLUSTER_GSEA_NPERM_SIMPLE", 10000L)

safe_file_component <- function(x) {
  x <- sanitize_path_component(x)
  x <- gsub("[.]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  ifelse(nzchar(x), x, "value")
}

as_chr1 <- function(x, default = NA_character_) {
  if (is.null(x) || length(x) == 0) return(default)
  x <- as.character(x)[1]
  if (is.na(x) || !nzchar(x)) return(default)
  x
}

not_missing <- function(x) {
  !is.na(x) & nzchar(as.character(x))
}

path_from_03a_suffix <- function(path) {
  path <- as_chr1(path)
  if (is.na(path) || !grepl("03a_DEGs", path, fixed = TRUE)) return(character(0))
  suffix <- sub("^.*03a_DEGs[/\\\\]?", "", path)
  if (!nzchar(suffix) || identical(suffix, path)) return(character(0))
  file.path(input_deg_root, suffix)
}

locate_marker_file <- function(task_row) {
  candidates <- character(0)
  output_file <- as_chr1(task_row$output_file)
  if (!is.na(output_file)) {
    candidates <- c(candidates, output_file, path_from_03a_suffix(output_file))
  }

  scope <- as_chr1(task_row$scope)
  cluster <- as_chr1(task_row$cluster)
  dose <- as_chr1(task_row$dose)

  if (identical(scope, "clusters_vs_rest") && !is.na(cluster)) {
    cluster_stub <- safe_file_component(cluster)
    candidates <- c(
      candidates,
      file.path(input_deg_root, "01_clusters_vs_rest", paste0("cluster_", cluster_stub, "_vs_rest_markers.csv"))
    )
  } else if (identical(scope, "global_Ploidy_4N_vs_2N")) {
    candidates <- c(
      candidates,
      file.path(input_deg_root, "02_Ploidy_4N_vs_2N", "Ploidy_4N_vs_2N_markers.csv")
    )
  } else if (identical(scope, "global_TN_Tumor_vs_CellLine")) {
    candidates <- c(
      candidates,
      file.path(input_deg_root, "03_TN_Tumor_vs_CellLine", "TN_Tumor_vs_CellLine_markers.csv")
    )
  } else if (identical(scope, "global_Dose_vs_rest") && !is.na(dose)) {
    dose_stub <- safe_file_component(dose)
    candidates <- c(
      candidates,
      file.path(input_deg_root, "04_Dose_vs_rest", paste0("Dose_", dose_stub, "_vs_rest_markers.csv"))
    )
  } else if (grepl("^within_cluster_", scope) && !is.na(cluster)) {
    cluster_stub <- safe_file_component(cluster)
    cluster_dir <- file.path(input_deg_root, "05_within_clusters", paste0("cluster_", cluster_stub))

    if (identical(scope, "within_cluster_Ploidy_4N_vs_2N")) {
      candidates <- c(candidates, file.path(cluster_dir, "01_Ploidy_4N_vs_2N", "Ploidy_4N_vs_2N_markers.csv"))
    } else if (identical(scope, "within_cluster_TN_Tumor_vs_CellLine")) {
      candidates <- c(candidates, file.path(cluster_dir, "02_TN_Tumor_vs_CellLine", "TN_Tumor_vs_CellLine_markers.csv"))
    } else if (identical(scope, "within_cluster_Dose_vs_rest") && !is.na(dose)) {
      dose_stub <- safe_file_component(dose)
      candidates <- c(candidates, file.path(cluster_dir, "03_Dose_vs_rest", paste0("Dose_", dose_stub, "_vs_rest_markers.csv")))
    } else if (identical(scope, "within_cluster_Tumor_only_Ploidy_4N_vs_2N")) {
      candidates <- c(candidates, file.path(cluster_dir, "04_Tumor_only_Ploidy_4N_vs_2N", "Tumor_only_Ploidy_4N_vs_2N_markers.csv"))
    } else if (identical(scope, "within_cluster_CellLine_only_Ploidy_4N_vs_2N")) {
      candidates <- c(candidates, file.path(cluster_dir, "05_CellLine_only_Ploidy_4N_vs_2N", "CellLine_only_Ploidy_4N_vs_2N_markers.csv"))
    } else if (identical(scope, "within_cluster_Dose_Ploidy_4N_vs_2N") && !is.na(dose)) {
      dose_stub <- safe_file_component(dose)
      candidates <- c(
        candidates,
        file.path(cluster_dir, "06_Ploidy_within_Dose_4N_vs_2N", paste0("Dose_", dose_stub, "_Ploidy_4N_vs_2N_markers.csv"))
      )
    }
  }

  candidates <- unique(candidates[not_missing(candidates)])
  hits <- candidates[file.exists(candidates) & file.info(candidates)$size > 0]
  if (length(hits) > 0) return(normalizePath(hits[1], mustWork = TRUE))
  NA_character_
}

read_deg_table <- function(marker_file) {
  de <- readr::read_csv(marker_file, show_col_types = FALSE)
  de <- as.data.frame(de, stringsAsFactors = FALSE)
  if (!("gene" %in% colnames(de))) {
    de$gene <- rownames(de)
  }
  if (!("gene_symbol" %in% colnames(de))) {
    de$gene_symbol <- clean_gene_symbols(de$gene)
  } else {
    de$gene_symbol <- clean_gene_symbols(de$gene_symbol)
    missing_symbol <- is.na(de$gene_symbol) | de$gene_symbol == ""
    de$gene_symbol[missing_symbol] <- clean_gene_symbols(de$gene[missing_symbol])
  }
  if (!("p_val_adj" %in% colnames(de)) && "p_val" %in% colnames(de)) {
    de$p_val_adj <- de$p_val
  }
  de
}

prepare_ranked_stats <- function(de) {
  de <- as.data.frame(de, stringsAsFactors = FALSE)
  lfc_col <- resolve_lfc_col(de)
  p_col <- if ("p_val_adj" %in% colnames(de)) "p_val_adj" else if ("p_val" %in% colnames(de)) "p_val" else NA_character_
  if (is.na(p_col)) {
    stop("DEG table is missing p_val_adj and p_val.", call. = FALSE)
  }

  gene_symbol <- if ("gene_symbol" %in% colnames(de)) clean_gene_symbols(de$gene_symbol) else rep(NA_character_, nrow(de))
  missing_symbol <- is.na(gene_symbol) | gene_symbol == ""
  if ("gene" %in% colnames(de)) {
    gene_symbol[missing_symbol] <- clean_gene_symbols(de$gene[missing_symbol])
  }

  lfc <- suppressWarnings(as.numeric(de[[lfc_col]]))
  p_value <- suppressWarnings(as.numeric(de[[p_col]]))
  finite_positive_p <- p_value[is.finite(p_value) & p_value > 0]
  p_floor <- if (length(finite_positive_p) > 0) {
    max(min(finite_positive_p, na.rm = TRUE) * 0.1, 1e-300)
  } else {
    1e-300
  }
  p_value[!is.finite(p_value) | p_value <= 0] <- p_floor
  p_value[p_value > 1] <- 1

  rank_metric <- sign(lfc) * (-log10(p_value) + abs(lfc) * 1e-6)
  rank_df <- data.frame(
    gene_symbol = gene_symbol,
    lfc = lfc,
    p_value = p_value,
    rank_metric = rank_metric,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      !is.na(gene_symbol),
      gene_symbol != "",
      is.finite(rank_metric),
      rank_metric != 0
    ) |>
    dplyr::arrange(dplyr::desc(abs(rank_metric)), dplyr::desc(abs(lfc)), gene_symbol) |>
    dplyr::group_by(gene_symbol) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(rank_metric), gene_symbol)

  stats <- rank_df$rank_metric
  names(stats) <- rank_df$gene_symbol
  stats
}

run_hallmark_gsea <- function(stats, pathways, min_size = 15L, max_size = 500L) {
  stats <- stats[is.finite(stats)]
  stats <- sort(stats, decreasing = TRUE)
  if (length(stats) == 0) return(data.frame())

  gsea_res <- tryCatch(
    {
      if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
        fgsea::fgseaMultilevel(
          pathways = pathways,
          stats = stats,
          minSize = min_size,
          maxSize = max_size,
          nPermSimple = gsea_nperm_simple,
          eps = 0
        )
      } else {
        fgsea::fgsea(
          pathways = pathways,
          stats = stats,
          minSize = min_size,
          maxSize = max_size,
          nperm = 10000
        )
      }
    },
    error = function(e) {
      warning("fgsea failed: ", conditionMessage(e), call. = FALSE)
      NULL
    }
  )
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(data.frame())

  out <- as.data.frame(gsea_res, stringsAsFactors = FALSE)
  if ("leadingEdge" %in% colnames(out)) {
    out$leadingEdge <- vapply(out$leadingEdge, function(x) paste(as.character(x), collapse = ";"), character(1))
  }
  out$hallmark_label <- beautify_hallmark_name(out$pathway)
  out$direction <- ifelse(out$NES >= 0, "positive", "negative")
  out <- out[order(out$padj, -abs(out$NES), out$pathway), , drop = FALSE]
  rownames(out) <- NULL
  out
}

select_top_gsea_terms <- function(gsea_df, top_n = 20L) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) return(data.frame())
  df <- as.data.frame(gsea_df, stringsAsFactors = FALSE)
  df <- df[is.finite(df$NES) & !is.na(df$pathway), , drop = FALSE]
  if (nrow(df) == 0) return(data.frame())

  df$padj_sort <- ifelse(is.na(df$padj), Inf, df$padj)
  df$sig_rank <- ifelse(!is.na(df$padj) & df$padj < 0.05, 0L, 1L)
  df <- df[order(df$sig_rank, -abs(df$NES), df$padj_sort, df$pathway), , drop = FALSE]
  head(df, top_n)
}

write_empty_plot <- function(out_pdf, out_png = NULL, message, width = 9, height = 5) {
  .ensure_dir(dirname(out_pdf))
  grDevices::pdf(out_pdf, width = width, height = height)
  graphics::plot.new()
  graphics::text(0.5, 0.5, message)
  grDevices::dev.off()
  if (!is.null(out_png)) {
    .ensure_dir(dirname(out_png))
    grDevices::png(out_png, width = width * 300, height = height * 300, res = 300)
    graphics::plot.new()
    graphics::text(0.5, 0.5, message)
    grDevices::dev.off()
  }
  invisible(NULL)
}

plot_gsea_top_bar <- function(gsea_df, out_pdf, out_png, title, top_n = 20L) {
  df <- select_top_gsea_terms(gsea_df, top_n = top_n)
  if (nrow(df) == 0) {
    return(write_empty_plot(out_pdf, out_png, paste0(title, "\nNo GSEA result")))
  }

  padj_plot <- suppressWarnings(as.numeric(df$padj))
  padj_plot[!is.finite(padj_plot) | is.na(padj_plot)] <- 1
  df$neg_log10_fdr <- -log10(pmax(padj_plot, 1e-300))
  df$hallmark_label <- factor(df$hallmark_label, levels = rev(df$hallmark_label))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = NES, y = hallmark_label)) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
    ggplot2::geom_point(
      ggplot2::aes(size = neg_log10_fdr, color = NES),
      alpha = 0.9
    ) +
    ggplot2::scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c", midpoint = 0) +
    ggplot2::scale_size_continuous(range = c(2, 6)) +
    ggplot2::labs(
      title = title,
      x = "Normalized enrichment score",
      y = NULL,
      color = "NES",
      size = "-log10(FDR)"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  plot_height <- max(5, 0.28 * nrow(df) + 1.5)
  ggplot2::ggsave(out_pdf, p, width = 9.5, height = plot_height)
  ggplot2::ggsave(out_png, p, width = 9.5, height = plot_height, dpi = 300)
  invisible(df)
}

plot_gsea_running_top <- function(gsea_df, stats, pathways, out_pdf, title, top_n = 20L) {
  df <- select_top_gsea_terms(gsea_df, top_n = top_n)
  if (nrow(df) == 0) {
    return(write_empty_plot(out_pdf, NULL, paste0(title, "\nNo GSEA result")))
  }

  .ensure_dir(dirname(out_pdf))
  grDevices::pdf(out_pdf, width = 9, height = 6, onefile = TRUE)
  on.exit(grDevices::dev.off(), add = TRUE)

  plotted <- 0L
  for (i in seq_len(nrow(df))) {
    pathway_name <- df$pathway[i]
    pathway_genes <- intersect(as.character(pathways[[pathway_name]]), names(stats))
    if (length(pathway_genes) < gsea_min_size) next

    p <- tryCatch(
      fgsea::plotEnrichment(pathway_genes, stats) +
        ggplot2::labs(
          title = paste0(title, " | ", df$hallmark_label[i]),
          subtitle = paste0(
            "NES=", signif(df$NES[i], 3),
            "; FDR=", signif(df$padj[i], 3),
            "; ES=", signif(df$ES[i], 3)
          )
        ) +
        ggplot2::theme_classic(base_size = 11),
      error = function(e) {
        warning("Failed to plot GSEA enrichment for ", pathway_name, ": ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (!is.null(p)) {
      print(p)
      plotted <- plotted + 1L
    }
  }

  if (plotted == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, paste0(title, "\nNo plottable top GSEA pathways"))
  }
  invisible(df)
}

make_heatmap_colors <- function(n = 101L) {
  grDevices::colorRampPalette(c("#2c7bb6", "white", "#d7191c"))(n)
}

plot_nes_heatmap <- function(
  gsea_df,
  column_col,
  out_pdf,
  title,
  top_n = 20L,
  column_levels = NULL,
  cluster_cols = FALSE
) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0) {
    return(write_empty_plot(out_pdf, NULL, paste0(title, "\nNo GSEA result"), width = 9, height = 5))
  }

  df <- as.data.frame(gsea_df, stringsAsFactors = FALSE)
  if (!(column_col %in% colnames(df))) {
    stop("Missing heatmap column: ", column_col, call. = FALSE)
  }
  df$column_value <- as.character(df[[column_col]])
  df <- df |>
    dplyr::filter(
      !is.na(column_value),
      column_value != "",
      !is.na(hallmark_label),
      is.finite(NES)
    )
  if (nrow(df) == 0) {
    return(write_empty_plot(out_pdf, NULL, paste0(title, "\nNo finite NES values"), width = 9, height = 5))
  }

  top_pathways <- df |>
    dplyr::group_by(hallmark_label) |>
    dplyr::summarise(
      best_abs_nes = max(abs(NES), na.rm = TRUE),
      best_fdr = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(best_abs_nes), best_fdr, hallmark_label) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(hallmark_label)

  column_values <- if (is.null(column_levels)) unique(df$column_value) else as.character(column_levels)
  column_values <- column_values[column_values %in% unique(df$column_value)]
  if (length(column_values) == 0 || length(top_pathways) == 0) {
    return(write_empty_plot(out_pdf, NULL, paste0(title, "\nNo heatmap entries"), width = 9, height = 5))
  }

  heatmap_use <- df |>
    dplyr::filter(hallmark_label %in% top_pathways, column_value %in% column_values) |>
    dplyr::arrange(dplyr::desc(abs(NES)), padj) |>
    dplyr::group_by(hallmark_label, column_value) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  heatmap_mat <- matrix(
    0,
    nrow = length(top_pathways),
    ncol = length(column_values),
    dimnames = list(top_pathways, column_values)
  )
  row_idx <- match(heatmap_use$hallmark_label, rownames(heatmap_mat))
  col_idx <- match(heatmap_use$column_value, colnames(heatmap_mat))
  keep_idx <- !is.na(row_idx) & !is.na(col_idx)
  heatmap_mat[cbind(row_idx[keep_idx], col_idx[keep_idx])] <- heatmap_use$NES[keep_idx]

  max_abs <- max(abs(heatmap_mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  breaks <- seq(-max_abs, max_abs, length.out = 102)
  width <- max(9, min(34, 0.28 * ncol(heatmap_mat) + 5))
  height <- max(6, min(16, 0.28 * nrow(heatmap_mat) + 3))

  .ensure_dir(dirname(out_pdf))
  grDevices::pdf(out_pdf, width = width, height = height)
  pheatmap::pheatmap(
    heatmap_mat,
    color = make_heatmap_colors(101),
    breaks = breaks,
    cluster_rows = nrow(heatmap_mat) > 2,
    cluster_cols = isTRUE(cluster_cols) && ncol(heatmap_mat) > 2,
    border_color = NA,
    fontsize_row = 8,
    fontsize_col = if (ncol(heatmap_mat) > 40) 5 else 7,
    angle_col = 45,
    main = title
  )
  grDevices::dev.off()
  invisible(heatmap_mat)
}

add_task_plot_label <- function(tasks) {
  tasks <- as.data.frame(tasks, stringsAsFactors = FALSE)
  label <- tasks$comparison
  has_cluster <- not_missing(tasks$cluster)
  label[has_cluster] <- paste0("cluster ", tasks$cluster[has_cluster], " | ", label[has_cluster])
  label <- gsub("_", " ", label, fixed = TRUE)
  tasks$plot_label <- label
  tasks
}

message("Preparing output directories.")
.ensure_dir(output_root)
out_summary <- .ensure_dir(file.path(output_root, "00_summary"))
out_ora <- .ensure_dir(file.path(output_root, "01_cluster_Hallmark_ORA_annotation"))
out_cluster_gsea <- .ensure_dir(file.path(output_root, "02_cluster_Hallmark_GSEA"))
out_stratified_gsea <- .ensure_dir(file.path(output_root, "03_stratified_Hallmark_GSEA"))
out_plots <- .ensure_dir(file.path(output_root, "04_plots"))
out_heatmaps <- .ensure_dir(file.path(output_root, "05_heatmaps"))

if (!dir.exists(input_deg_root)) {
  stop("Input 03a_DEGs directory does not exist: ", input_deg_root, call. = FALSE)
}
if (!file.exists(comparison_summary_file)) {
  stop("Missing DEG comparison summary: ", comparison_summary_file, call. = FALSE)
}

message("Reading DEG comparison summary.")
comparison_summary <- readr::read_csv(comparison_summary_file, show_col_types = FALSE)
comparison_summary <- as.data.frame(comparison_summary, stringsAsFactors = FALSE)
required_summary_cols <- c("scope", "comparison", "group_col", "ident_1", "ident_2", "status", "output_file")
missing_summary_cols <- setdiff(required_summary_cols, colnames(comparison_summary))
if (length(missing_summary_cols) > 0) {
  stop("DEG comparison summary is missing required columns: ", paste(missing_summary_cols, collapse = ", "), call. = FALSE)
}
for (optional_col in c("cluster", "subset_group", "dose", "n_ident_1", "n_ident_2")) {
  if (!(optional_col %in% colnames(comparison_summary))) comparison_summary[[optional_col]] <- NA
}

marker_tasks <- comparison_summary |>
  dplyr::filter(status %in% c("success", "cached")) |>
  as.data.frame(stringsAsFactors = FALSE)
marker_tasks$marker_file <- vapply(
  seq_len(nrow(marker_tasks)),
  function(i) locate_marker_file(marker_tasks[i, , drop = FALSE]),
  character(1)
)
marker_tasks$has_marker_file <- !is.na(marker_tasks$marker_file) & file.exists(marker_tasks$marker_file)

missing_marker_tasks <- marker_tasks[!marker_tasks$has_marker_file, , drop = FALSE]
if (nrow(missing_marker_tasks) > 0) {
  write_table_csv(
    missing_marker_tasks,
    file.path(out_summary, "DEG_tasks_missing_marker_files.csv")
  )
  warning("Missing marker files for ", nrow(missing_marker_tasks), " successful/cached DEG comparison(s). See 00_summary/DEG_tasks_missing_marker_files.csv.", call. = FALSE)
}
marker_tasks <- marker_tasks[marker_tasks$has_marker_file, , drop = FALSE]
marker_tasks <- add_task_plot_label(marker_tasks)
write_table_csv(marker_tasks, file.path(out_summary, "DEG_marker_files_used.csv"))

message("Preparing Hallmark gene sets.")
hallmark_sets <- get_hallmark_sets(species = "Homo sapiens")

metadata_counts <- if (file.exists(metadata_counts_file)) {
  as.data.frame(readr::read_csv(metadata_counts_file, show_col_types = FALSE), stringsAsFactors = FALSE)
} else {
  data.frame()
}

cluster_tasks <- marker_tasks |>
  dplyr::filter(scope == cluster_scope) |>
  as.data.frame(stringsAsFactors = FALSE)
if (nrow(cluster_tasks) == 0) {
  stop("No cluster vs rest marker files were found under 03a_DEGs.", call. = FALSE)
}
cluster_levels_from_counts <- if (nrow(metadata_counts) > 0 && all(c("metadata", "value") %in% colnames(metadata_counts))) {
  as.character(metadata_counts$value[metadata_counts$metadata == analysis_cluster_label])
} else {
  character(0)
}
cluster_levels <- c(
  intersect(cluster_levels_from_counts, cluster_tasks$cluster),
  setdiff(as.character(cluster_tasks$cluster), cluster_levels_from_counts)
)
cluster_levels <- unique(cluster_levels[not_missing(cluster_levels)])
cluster_tasks <- cluster_tasks[match(cluster_levels, cluster_tasks$cluster), , drop = FALSE]

message("Building ORA universe from cluster vs rest DEG tables.")
cluster_deg_cache <- list()
universe_genes <- character(0)
for (i in seq_len(nrow(cluster_tasks))) {
  cluster_id <- as.character(cluster_tasks$cluster[i])
  de <- read_deg_table(cluster_tasks$marker_file[i])
  cluster_deg_cache[[cluster_id]] <- de
  universe_genes <- union(universe_genes, unique(stats::na.omit(de$gene_symbol)))
}
universe_genes <- sort(unique(universe_genes[!is.na(universe_genes) & universe_genes != ""]))
write_table_csv(
  data.frame(gene_symbol = universe_genes, stringsAsFactors = FALSE),
  file.path(out_summary, "hallmark_ora_universe_genes_from_cluster_DEGs.csv")
)

message("Running cluster Hallmark ORA annotation and cluster Hallmark GSEA.")
deg_summary_list <- list()
ora_all_list <- list()
annotation_summary_list <- list()
cluster_gsea_list <- list()
gsea_summary_list <- list()

for (i in seq_len(nrow(cluster_tasks))) {
  task <- cluster_tasks[i, , drop = FALSE]
  cluster_id <- as.character(task$cluster)
  cluster_stub <- safe_file_component(cluster_id)
  de <- cluster_deg_cache[[cluster_id]]
  lfc_col <- resolve_lfc_col(de)

  n_cells_cluster <- if (nrow(metadata_counts) > 0 && all(c("metadata", "value", "n_cells") %in% colnames(metadata_counts))) {
    hit <- metadata_counts$n_cells[metadata_counts$metadata == analysis_cluster_label & metadata_counts$value == cluster_id]
    if (length(hit) > 0) as.integer(hit[1]) else as.integer(task$n_ident_1)
  } else {
    as.integer(task$n_ident_1)
  }

  cluster_ora_dir <- .ensure_dir(file.path(out_ora, paste0("cluster_", cluster_stub)))
  cluster_gsea_dir <- .ensure_dir(file.path(out_cluster_gsea, paste0("cluster_", cluster_stub)))

  ora_deg <- select_deg_top_up_for_ora(
    df = de,
    lfc_col = lfc_col,
    padj_max = deg_padj_cutoff,
    abs_logfc_min = deg_lfc_cutoff_for_ora,
    abs_delta_pct_min = deg_abs_delta_pct_cutoff_for_ora,
    top_n = deg_top_n_for_ora
  )
  write_table_csv(
    ora_deg,
    file.path(cluster_ora_dir, paste0("cluster_", cluster_stub, "_vs_rest_top100_up_for_Hallmark_ORA.csv"))
  )

  ora_res <- run_ora_hypergeom(
    query_genes = ora_deg$gene_symbol,
    universe_genes = universe_genes,
    pathways = hallmark_sets,
    min_size = min_hallmark_size,
    max_size = max_hallmark_size,
    min_overlap = ora_min_overlap
  )
  ora_res <- score_ora_for_annotation(
    ora_df = ora_res,
    deg_df = ora_deg,
    weight_expression = annotation_weight_expression,
    weight_detection = annotation_weight_detection,
    weight_overlap = annotation_weight_overlap,
    fdr_cutoff = ora_annotation_fdr_cutoff
  )
  if (nrow(ora_res) > 0) {
    ora_res$cluster <- cluster_id
    ora_res <- ora_res |>
      dplyr::select(
        cluster,
        annotation_rank,
        pathway,
        hallmark_label,
        set_size,
        query_size,
        overlap,
        overlap_query_fraction,
        overlap_set_fraction,
        overlap_score_raw,
        marker_mean_abs_logfc,
        marker_mean_pct1,
        marker_mean_abs_delta_pct,
        n_overlap_up,
        n_overlap_down,
        expression_score_raw,
        detection_score_raw,
        expression_score_scaled,
        detection_score_scaled,
        overlap_score_scaled,
        annotation_score,
        gene_ratio,
        bg_ratio,
        odds_ratio,
        p_value,
        p_adj,
        overlap_genes
      )
    ora_all_list[[cluster_id]] <- ora_res
  }
  write_table_csv(ora_res, file.path(cluster_ora_dir, paste0("cluster_", cluster_stub, "_Hallmark_ORA.csv")))
  plot_ora_top(
    ora_df = ora_res,
    out_pdf = file.path(cluster_ora_dir, paste0("cluster_", cluster_stub, "_Hallmark_ORA_top20.pdf")),
    out_png = file.path(cluster_ora_dir, paste0("cluster_", cluster_stub, "_Hallmark_ORA_top20.png")),
    title = paste0("cluster ", cluster_id, " vs rest | Hallmark ORA"),
    top_n = top_plot_terms
  )

  annotation_summary_list[[cluster_id]] <- summarize_cluster_annotation(
    cluster_id = cluster_id,
    n_cells = n_cells_cluster,
    deg_df = ora_deg,
    ora_df = ora_res,
    top_n = top_terms_per_cluster,
    fdr_cutoff = ora_annotation_fdr_cutoff
  )
  deg_summary_list[[cluster_id]] <- data.frame(
    cluster = cluster_id,
    n_cells = n_cells_cluster,
    tested_genes = nrow(de),
    sig_genes = sum(!is.na(de$p_val_adj) & de$p_val_adj < deg_padj_cutoff),
    deg_for_ora = nrow(ora_deg),
    up_genes_for_ora = nrow(ora_deg),
    marker_file = task$marker_file,
    stringsAsFactors = FALSE
  )

  stats <- prepare_ranked_stats(de)
  gsea_res <- run_hallmark_gsea(
    stats = stats,
    pathways = hallmark_sets,
    min_size = gsea_min_size,
    max_size = gsea_max_size
  )
  if (nrow(gsea_res) > 0) {
    gsea_res$cluster <- cluster_id
    gsea_res$scope <- cluster_scope
    gsea_res$comparison <- as.character(task$comparison)
    gsea_res$marker_file <- task$marker_file
    gsea_res <- gsea_res |>
      dplyr::select(scope, cluster, comparison, marker_file, dplyr::everything())
    cluster_gsea_list[[cluster_id]] <- gsea_res
  }
  write_table_csv(gsea_res, file.path(cluster_gsea_dir, paste0("cluster_", cluster_stub, "_Hallmark_GSEA.csv")))

  plot_gsea_top_bar(
    gsea_df = gsea_res,
    out_pdf = file.path(cluster_gsea_dir, paste0("cluster_", cluster_stub, "_Hallmark_GSEA_top20_NES.pdf")),
    out_png = file.path(cluster_gsea_dir, paste0("cluster_", cluster_stub, "_Hallmark_GSEA_top20_NES.png")),
    title = paste0("cluster ", cluster_id, " vs rest | Hallmark GSEA"),
    top_n = top_plot_terms
  )
  plot_gsea_running_top(
    gsea_df = gsea_res,
    stats = stats,
    pathways = hallmark_sets,
    out_pdf = file.path(cluster_gsea_dir, paste0("cluster_", cluster_stub, "_Hallmark_GSEA_running_top20.pdf")),
    title = paste0("cluster ", cluster_id, " vs rest"),
    top_n = top_plot_terms
  )

  gsea_summary_list[[paste0("cluster_", cluster_id)]] <- data.frame(
    scope = cluster_scope,
    cluster = cluster_id,
    comparison = as.character(task$comparison),
    marker_file = task$marker_file,
    ranked_genes = length(stats),
    tested_pathways = nrow(gsea_res),
    significant_pathways = if (nrow(gsea_res) > 0) sum(!is.na(gsea_res$padj) & gsea_res$padj < 0.05) else 0L,
    stringsAsFactors = FALSE
  )
}

deg_summary_df <- dplyr::bind_rows(deg_summary_list)
write_table_csv(deg_summary_df, file.path(out_summary, "cluster_deg_for_annotation_summary.csv"))

ora_all_df <- if (length(ora_all_list) > 0) dplyr::bind_rows(ora_all_list) else data.frame()
write_table_csv(ora_all_df, file.path(out_summary, "cluster_Hallmark_ORA_all.csv"))

annotation_summary_df <- dplyr::bind_rows(annotation_summary_list) |>
  dplyr::arrange(match(cluster, cluster_levels))
write_table_csv(annotation_summary_df, file.path(out_summary, "cluster_annotation_summary.csv"))

cluster_gsea_all <- if (length(cluster_gsea_list) > 0) dplyr::bind_rows(cluster_gsea_list) else data.frame()
write_table_csv(cluster_gsea_all, file.path(out_summary, "cluster_Hallmark_GSEA_all.csv"))

if (nrow(cluster_gsea_all) > 0) {
  plot_nes_heatmap(
    gsea_df = cluster_gsea_all,
    column_col = "cluster",
    out_pdf = file.path(out_heatmaps, "cluster_Hallmark_GSEA_NES_heatmap_top20.pdf"),
    title = "Cluster Hallmark GSEA NES",
    top_n = top_plot_terms,
    column_levels = cluster_levels,
    cluster_cols = FALSE
  )
}

message("Running Hallmark GSEA for non-cluster stratified DEG comparisons.")
stratified_tasks <- marker_tasks |>
  dplyr::filter(scope != cluster_scope) |>
  dplyr::arrange(scope, cluster, subset_group, dose, comparison) |>
  as.data.frame(stringsAsFactors = FALSE)

stratified_gsea_list <- list()
if (nrow(stratified_tasks) > 0) {
  for (i in seq_len(nrow(stratified_tasks))) {
    task <- stratified_tasks[i, , drop = FALSE]
    scope <- as.character(task$scope)
    comparison <- as.character(task$comparison)
    scope_stub <- safe_file_component(scope)
    comparison_stub <- safe_file_component(comparison)
    if (!is.na(task$cluster) && nzchar(as.character(task$cluster))) {
      comparison_stub <- paste0("cluster_", safe_file_component(task$cluster), "_", comparison_stub)
    }
    if (!is.na(task$dose) && nzchar(as.character(task$dose))) {
      comparison_stub <- paste0(comparison_stub, "_", safe_file_component(task$dose))
    }

    comparison_dir <- .ensure_dir(file.path(out_stratified_gsea, scope_stub, comparison_stub))
    de <- read_deg_table(task$marker_file)
    stats <- prepare_ranked_stats(de)
    gsea_res <- run_hallmark_gsea(
      stats = stats,
      pathways = hallmark_sets,
      min_size = gsea_min_size,
      max_size = gsea_max_size
    )

    if (nrow(gsea_res) > 0) {
      gsea_res$scope <- scope
      gsea_res$cluster <- task$cluster
      gsea_res$subset_group <- task$subset_group
      gsea_res$dose <- task$dose
      gsea_res$comparison <- comparison
      gsea_res$group_col <- task$group_col
      gsea_res$ident_1 <- task$ident_1
      gsea_res$ident_2 <- task$ident_2
      gsea_res$n_ident_1 <- task$n_ident_1
      gsea_res$n_ident_2 <- task$n_ident_2
      gsea_res$plot_label <- task$plot_label
      gsea_res$marker_file <- task$marker_file
      gsea_res <- gsea_res |>
        dplyr::select(
          scope,
          cluster,
          subset_group,
          dose,
          comparison,
          group_col,
          ident_1,
          ident_2,
          n_ident_1,
          n_ident_2,
          plot_label,
          marker_file,
          dplyr::everything()
        )
      stratified_gsea_list[[length(stratified_gsea_list) + 1L]] <- gsea_res
    }

    write_table_csv(gsea_res, file.path(comparison_dir, paste0(comparison_stub, "_Hallmark_GSEA.csv")))
    plot_gsea_top_bar(
      gsea_df = gsea_res,
      out_pdf = file.path(comparison_dir, paste0(comparison_stub, "_Hallmark_GSEA_top20_NES.pdf")),
      out_png = file.path(comparison_dir, paste0(comparison_stub, "_Hallmark_GSEA_top20_NES.png")),
      title = paste0(gsub("_", " ", scope, fixed = TRUE), " | ", gsub("_", " ", comparison, fixed = TRUE)),
      top_n = top_plot_terms
    )
    plot_gsea_running_top(
      gsea_df = gsea_res,
      stats = stats,
      pathways = hallmark_sets,
      out_pdf = file.path(comparison_dir, paste0(comparison_stub, "_Hallmark_GSEA_running_top20.pdf")),
      title = paste0(gsub("_", " ", scope, fixed = TRUE), " | ", gsub("_", " ", comparison, fixed = TRUE)),
      top_n = top_plot_terms
    )

    gsea_summary_list[[paste0("stratified_", i)]] <- data.frame(
      scope = scope,
      cluster = task$cluster,
      comparison = comparison,
      marker_file = task$marker_file,
      ranked_genes = length(stats),
      tested_pathways = nrow(gsea_res),
      significant_pathways = if (nrow(gsea_res) > 0) sum(!is.na(gsea_res$padj) & gsea_res$padj < 0.05) else 0L,
      stringsAsFactors = FALSE
    )
  }
}

stratified_gsea_all <- if (length(stratified_gsea_list) > 0) dplyr::bind_rows(stratified_gsea_list) else data.frame()
write_table_csv(stratified_gsea_all, file.path(out_summary, "stratified_Hallmark_GSEA_all.csv"))

gsea_summary_df <- dplyr::bind_rows(gsea_summary_list)
write_table_csv(gsea_summary_df, file.path(out_summary, "GSEA_comparison_summary.csv"))

if (nrow(stratified_gsea_all) > 0) {
  plot_nes_heatmap(
    gsea_df = stratified_gsea_all,
    column_col = "plot_label",
    out_pdf = file.path(out_heatmaps, "stratified_Hallmark_GSEA_NES_heatmap_top20_all.pdf"),
    title = "Stratified Hallmark GSEA NES",
    top_n = top_plot_terms,
    cluster_cols = FALSE
  )

  for (scope_i in unique(stratified_gsea_all$scope)) {
    scope_df <- stratified_gsea_all[stratified_gsea_all$scope == scope_i, , drop = FALSE]
    scope_stub <- safe_file_component(scope_i)
    plot_nes_heatmap(
      gsea_df = scope_df,
      column_col = "plot_label",
      out_pdf = file.path(out_heatmaps, paste0(scope_stub, "_Hallmark_GSEA_NES_heatmap_top20.pdf")),
      title = paste0(gsub("_", " ", scope_i, fixed = TRUE), " Hallmark GSEA NES"),
      top_n = top_plot_terms,
      cluster_cols = FALSE
    )
  }
}

if (nrow(ora_all_df) > 0) {
  heatmap_df <- ora_all_df |>
    dplyr::mutate(score = annotation_score)
  top_pathways <- heatmap_df |>
    dplyr::group_by(hallmark_label) |>
    dplyr::summarise(best_score = max(score, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(best_score), hallmark_label) |>
    dplyr::slice_head(n = top_plot_terms) |>
    dplyr::pull(hallmark_label)
  heatmap_use <- heatmap_df |>
    dplyr::filter(hallmark_label %in% top_pathways) |>
    dplyr::group_by(hallmark_label, cluster) |>
    dplyr::summarise(score = max(score, na.rm = TRUE), .groups = "drop")
  heatmap_mat <- matrix(
    0,
    nrow = length(top_pathways),
    ncol = length(cluster_levels),
    dimnames = list(top_pathways, cluster_levels)
  )
  row_idx <- match(heatmap_use$hallmark_label, rownames(heatmap_mat))
  col_idx <- match(heatmap_use$cluster, colnames(heatmap_mat))
  keep_idx <- !is.na(row_idx) & !is.na(col_idx)
  heatmap_mat[cbind(row_idx[keep_idx], col_idx[keep_idx])] <- heatmap_use$score[keep_idx]

  grDevices::pdf(file.path(out_heatmaps, "cluster_Hallmark_ORA_annotation_score_heatmap_top20.pdf"), width = 10, height = 8)
  pheatmap::pheatmap(
    heatmap_mat,
    cluster_rows = nrow(heatmap_mat) > 2,
    cluster_cols = FALSE,
    border_color = NA,
    main = "Cluster Hallmark ORA annotation score"
  )
  grDevices::dev.off()
}

summary_lines <- c(
  "03b cluster annotation and GSEA completed.",
  paste0("Input DEG root: ", input_deg_root),
  paste0("Output root: ", output_root),
  paste0("Cluster DEG source: ", cluster_scope),
  paste0("Cluster labels: ", paste(cluster_levels, collapse = ", ")),
  paste0("Marker files used: ", nrow(marker_tasks)),
  paste0("Cluster ORA annotations generated for: ", nrow(annotation_summary_df), " cluster(s)."),
  paste0("Cluster GSEA comparisons: ", length(cluster_gsea_list)),
  paste0("Stratified GSEA comparisons: ", length(stratified_gsea_list)),
  paste0(
    "Hallmark ORA DEG rule: top ",
    deg_top_n_for_ora,
    " up genes after p_adj < ",
    deg_padj_cutoff,
    ", abs(logFC) >= ",
    deg_lfc_cutoff_for_ora,
    ", abs(pct.1 - pct.2) >= ",
    deg_abs_delta_pct_cutoff_for_ora,
    "."
  ),
  paste0(
    "Cluster annotation score: ",
    annotation_weight_expression,
    " * marker_strength + ",
    annotation_weight_detection,
    " * marker_detection_fraction + ",
    annotation_weight_overlap,
    " * overlap_degree, after Hallmark ORA FDR < ",
    ora_annotation_fdr_cutoff,
    "."
  ),
  paste0("GSEA method: fgsea preranked Hallmark analysis with minSize=", gsea_min_size, ", maxSize=", gsea_max_size, "."),
  paste0("fgseaMultilevel nPermSimple: ", gsea_nperm_simple),
  "GSEA ranking metric: sign(logFC) * (-log10(adjusted P) + abs(logFC) * 1e-6).",
  paste0("Top pathways plotted per comparison: ", top_plot_terms)
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("03b_cluster_annotation_and_GSEA finished.")
message("Output root: ", output_root)
