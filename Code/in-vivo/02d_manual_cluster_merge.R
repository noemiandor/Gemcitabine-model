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

set.seed(12345)

input_rds <- file.path(results_root, "02b_cluster_refine", "objects", "integrated_sct_cca_seurat_cluster_refine.rds")
output_root <- file.path(results_root, "02d_manual_cluster_merge")
base_cluster_col <- "seurat_cluster_refine"

manual_merge_map <- list(
  `0` = c("0", "1", "7"),
  `10` = c("10", "11", "12")
)

marker_min_pct <- 0.10
marker_logfc_threshold <- 0.25
marker_test_use <- "wilcox"
top_markers_n <- 20L
top_dotplot_n <- 5L
cosine_top_n_each <- 50L

prepare_obj_for_markers <- function(obj) {
  if (!("RNA" %in% names(obj@assays))) {
    stop("RNA assay not found in Seurat object.", call. = FALSE)
  }
  DefaultAssay(obj) <- "RNA"
  obj <- maybe_join_layers(obj, assay = "RNA")
  rna_data <- get_assay_data_slot(obj, assay = "RNA", slot_name = "data")
  if (is.null(rna_data) || nrow(rna_data) == 0 || ncol(rna_data) == 0) {
    message("RNA data slot is empty. Running NormalizeData on RNA assay.")
    obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  obj
}

build_manual_merge_labels <- function(cluster_vec, merge_map, cluster_order = NULL) {
  cluster_vec <- as.character(cluster_vec)
  out <- cluster_vec

  for (merge_name in names(merge_map)) {
    members <- as.character(merge_map[[merge_name]])
    idx <- cluster_vec %in% members
    out[idx] <- merge_name
  }

  original_order <- if (!is.null(cluster_order)) as.character(cluster_order) else sort_maybe_numeric(cluster_vec)
  merged_level_order <- character(0)
  seen_labels <- character(0)
  for (cl in original_order) {
    target_label <- if (cl %in% unlist(merge_map, use.names = FALSE)) {
      names(merge_map)[vapply(merge_map, function(x) cl %in% x, logical(1))][1]
    } else {
      cl
    }
    if (!(target_label %in% seen_labels)) {
      merged_level_order <- c(merged_level_order, target_label)
      seen_labels <- c(seen_labels, target_label)
    }
  }

  factor(out, levels = merged_level_order)
}

write_table_csv <- function(df, file_path) {
  readr::write_csv(df, file_path)
}

run_markers_vs_rest <- function(obj, ident_col, groups, out_dir, top_n = 20L) {
  .ensure_dir(out_dir)
  stopifnot(length(groups) > 0)

  Idents(obj) <- obj@meta.data[[ident_col]]

  all_markers <- list()
  top_markers <- list()
  marker_tables <- list()

  for (grp in groups) {
    marker_file <- file.path(out_dir, paste0("markers_", grp, "_vs_rest.csv"))
    top_file <- file.path(out_dir, paste0("top_positive_markers_", grp, "_vs_rest.csv"))

    if (file.exists(marker_file)) {
      message("[markers vs rest] Reusing existing DEG file: ", marker_file)
      de <- readr::read_csv(marker_file, show_col_types = FALSE)
      de <- as.data.frame(de, stringsAsFactors = FALSE)
    } else {
      message("[markers vs rest] ", ident_col, " = ", grp)
      de <- FindMarkers(
        object = obj,
        ident.1 = grp,
        min.pct = marker_min_pct,
        logfc.threshold = marker_logfc_threshold,
        test.use = marker_test_use,
        verbose = FALSE
      )
      de <- de %>%
        tibble::rownames_to_column("gene")
    }
    lfc_col <- resolve_lfc_col(de)

    de$group <- grp
    de$comparison <- paste0(grp, "_vs_rest")

    if (!file.exists(marker_file)) {
      write_table_csv(de, marker_file)
    }

    if (file.exists(top_file)) {
      top_pos <- readr::read_csv(top_file, show_col_types = FALSE)
      top_pos <- as.data.frame(top_pos, stringsAsFactors = FALSE)
    } else {
      top_pos <- de %>%
        filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) %>%
        arrange(desc(.data[[lfc_col]]), p_val_adj, gene) %>%
        slice_head(n = top_n) %>%
        mutate(rank_within_group = row_number())

      write_table_csv(top_pos, top_file)
    }

    all_markers[[grp]] <- de
    top_markers[[grp]] <- top_pos
    marker_tables[[grp]] <- list(full = de, top_positive = top_pos, lfc_col = lfc_col)
  }

  list(
    full = bind_rows(all_markers),
    top_positive = bind_rows(top_markers),
    by_group = marker_tables
  )
}

run_pairwise_markers <- function(obj, ident_col, pair_df, out_dir) {
  .ensure_dir(out_dir)
  Idents(obj) <- obj@meta.data[[ident_col]]

  summary_rows <- list()
  pair_tables <- list()

  for (i in seq_len(nrow(pair_df))) {
    g1 <- as.character(pair_df$group_1[i])
    g2 <- as.character(pair_df$group_2[i])
    pair_name <- paste0(g1, "_vs_", g2)
    marker_file <- file.path(out_dir, paste0("markers_", pair_name, ".csv"))
    if (file.exists(marker_file)) {
      message("[pairwise markers] Reusing existing DEG file: ", marker_file)
      de <- readr::read_csv(marker_file, show_col_types = FALSE)
      de <- as.data.frame(de, stringsAsFactors = FALSE)
    } else {
      message("[pairwise markers] ", ident_col, " : ", pair_name)
      de <- FindMarkers(
        object = obj,
        ident.1 = g1,
        ident.2 = g2,
        min.pct = marker_min_pct,
        logfc.threshold = marker_logfc_threshold,
        test.use = marker_test_use,
        verbose = FALSE
      ) %>%
        tibble::rownames_to_column("gene")
    }

    lfc_col <- resolve_lfc_col(de)
    if (!file.exists(marker_file)) {
      write_table_csv(de, marker_file)
    }

    n_sig <- sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05)
    n_sig_abs025 <- sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05 & abs(de[[lfc_col]]) >= 0.25)
    top_up_gene <- de %>%
      filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) %>%
      arrange(desc(.data[[lfc_col]]), p_val_adj, gene) %>%
      slice_head(n = 1) %>%
      pull(gene)
    top_down_gene <- de %>%
      filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] < 0) %>%
      arrange(.data[[lfc_col]], p_val_adj, gene) %>%
      slice_head(n = 1) %>%
      pull(gene)

    summary_rows[[pair_name]] <- data.frame(
      ident_col = ident_col,
      group_1 = g1,
      group_2 = g2,
      n_sig = n_sig,
      n_sig_abs_log2fc_0.25 = n_sig_abs025,
      median_abs_log2fc = median(abs(de[[lfc_col]]), na.rm = TRUE),
      top_up_gene = ifelse(length(top_up_gene) == 0, NA_character_, top_up_gene[1]),
      top_down_gene = ifelse(length(top_down_gene) == 0, NA_character_, top_down_gene[1]),
      stringsAsFactors = FALSE
    )
    pair_tables[[pair_name]] <- de
  }

  summary_df <- bind_rows(summary_rows)
  write_table_csv(summary_df, file.path(out_dir, "pairwise_marker_summary.csv"))
  invisible(list(summary = summary_df, by_pair = pair_tables))
}

extract_top_positive_genes <- function(marker_tbl, n_top = 50L) {
  lfc_col <- resolve_lfc_col(marker_tbl)
  marker_tbl %>%
    filter(!is.na(p_val_adj), p_val_adj < 0.05, .data[[lfc_col]] > 0) %>%
    arrange(desc(.data[[lfc_col]]), p_val_adj, gene) %>%
    slice_head(n = n_top) %>%
    pull(gene) %>%
    unique()
}

build_marker_overlap_summary <- function(original_markers, merged_markers, merge_map, n_top = 50L) {
  out <- list()
  idx <- 1L

  for (merge_name in names(merge_map)) {
    merged_tbl <- merged_markers$by_group[[merge_name]]$full
    merged_genes <- extract_top_positive_genes(merged_tbl, n_top = n_top)

    source_clusters <- as.character(merge_map[[merge_name]])
    source_gene_sets <- lapply(source_clusters, function(cl) {
      extract_top_positive_genes(original_markers$by_group[[cl]]$full, n_top = n_top)
    })
    names(source_gene_sets) <- source_clusters

    union_genes <- unique(unlist(source_gene_sets, use.names = FALSE))
    intersect_genes <- Reduce(intersect, source_gene_sets)

    compare_sets <- c(source_gene_sets, list(union_of_sources = union_genes, intersect_of_sources = intersect_genes))
    for (ref_name in names(compare_sets)) {
      ref_genes <- compare_sets[[ref_name]]
      out[[idx]] <- data.frame(
        merged_group = merge_name,
        reference_set = ref_name,
        merged_top_n = length(merged_genes),
        reference_top_n = length(ref_genes),
        overlap_n = length(intersect(merged_genes, ref_genes)),
        jaccard = safe_jaccard(merged_genes, ref_genes),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  bind_rows(out)
}

normalize_similarity_gene_key <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^GRCh[0-9]+[-_]", "", x, ignore.case = TRUE)
  x <- sub("^GRCm39[-_]", "", x, ignore.case = TRUE)
  x <- sub("^hg38[-_]", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  x <- toupper(x)
  x[x == ""] <- NA_character_
  x
}

resolve_similarity_gene_label <- function(df) {
  gene_symbol <- if ("gene_symbol" %in% colnames(df)) as.character(df$gene_symbol) else rep(NA_character_, nrow(df))
  gene <- if ("gene" %in% colnames(df)) as.character(df$gene) else rownames(df)
  out <- ifelse(!is.na(gene_symbol) & gene_symbol != "", gene_symbol, gene)
  out[is.na(out) | out == ""] <- gene[is.na(out) | out == ""]
  out
}

prepare_marker_table_for_similarity <- function(df) {
  lfc_col <- resolve_lfc_col(df)

  df$gene_label <- resolve_similarity_gene_label(df)
  df$gene_key <- normalize_similarity_gene_key(df$gene_label)
  df$lfc_value <- as.numeric(df[[lfc_col]])
  df$p_val_adj_num <- as.numeric(df$p_val_adj)
  df$delta_pct <- as.numeric(df$pct.1) - as.numeric(df$pct.2)
  df$abs_log2fc <- abs(df$lfc_value)
  df$abs_delta_pct <- abs(df$delta_pct)

  df %>%
    filter(
      !is.na(gene_key),
      !is.na(lfc_value),
      !is.na(p_val_adj_num),
      !is.na(delta_pct)
    ) %>%
    arrange(desc(abs_log2fc), p_val_adj_num) %>%
    group_by(gene_key) %>%
    slice(1) %>%
    ungroup()
}

filter_similarity_feature_table <- function(df, spec) {
  keep <- df$p_val_adj_num < spec$padj_max &
    df$abs_log2fc >= spec$abs_log2fc_min &
    df$abs_delta_pct >= spec$abs_delta_pct_min

  filtered <- df[keep, , drop = FALSE]

  up_tbl <- filtered %>%
    filter(lfc_value > 0) %>%
    arrange(desc(lfc_value), p_val_adj_num) %>%
    slice_head(n = spec$top_n_up)

  down_tbl <- filtered %>%
    filter(lfc_value < 0) %>%
    arrange(lfc_value, p_val_adj_num) %>%
    slice_head(n = spec$top_n_down)

  bind_rows(up_tbl, down_tbl) %>%
    arrange(desc(abs_log2fc), p_val_adj_num) %>%
    distinct(gene_key, .keep_all = TRUE)
}

make_similarity_feature_object <- function(marker_tbl, spec) {
  feature_tbl <- prepare_marker_table_for_similarity(marker_tbl)
  feature_tbl <- filter_similarity_feature_table(feature_tbl, spec)

  signed_vec <- feature_tbl$lfc_value
  names(signed_vec) <- feature_tbl$gene_key

  list(
    table = feature_tbl,
    signed_vec = signed_vec
  )
}

make_union_vectors <- function(v1, v2) {
  genes <- union(names(v1), names(v2))
  if (length(genes) == 0) {
    return(list(x = numeric(0), y = numeric(0)))
  }

  x <- setNames(rep(0, length(genes)), genes)
  y <- setNames(rep(0, length(genes)), genes)
  x[names(v1)] <- v1
  y[names(v2)] <- v2

  list(x = unname(x), y = unname(y))
}

safe_signed_lfc_cosine <- function(v1, v2) {
  vv <- make_union_vectors(v1, v2)
  x <- vv$x
  y <- vv$y
  if (length(x) == 0) return(NA_real_)
  denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(x * y) / denom
}

build_manual_merge_cosine_summary <- function(
  original_markers,
  merged_markers,
  original_groups,
  merged_groups,
  merge_map,
  spec
) {
  original_feature_map <- setNames(
    lapply(original_groups, function(grp) make_similarity_feature_object(original_markers$by_group[[grp]]$full, spec)),
    original_groups
  )
  merged_feature_map <- setNames(
    lapply(merged_groups, function(grp) make_similarity_feature_object(merged_markers$by_group[[grp]]$full, spec)),
    merged_groups
  )

  rows <- list()
  idx <- 1L
  for (merge_name in merged_groups) {
    for (original_name in original_groups) {
      merged_feature <- merged_feature_map[[merge_name]]
      original_feature <- original_feature_map[[original_name]]

      rows[[idx]] <- data.frame(
        merged_group = merge_name,
        original_cluster = original_name,
        is_source_cluster = original_name %in% as.character(merge_map[[merge_name]]),
        merged_n_features = nrow(merged_feature$table),
        original_n_features = nrow(original_feature$table),
        signed_lfc_cosine = safe_signed_lfc_cosine(merged_feature$signed_vec, original_feature$signed_vec),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  bind_rows(rows)
}

cosine_summary_to_matrix <- function(cosine_df, merged_groups, original_groups) {
  mat <- matrix(
    NA_real_,
    nrow = length(merged_groups),
    ncol = length(original_groups),
    dimnames = list(merged_groups, original_groups)
  )

  for (i in seq_len(nrow(cosine_df))) {
    mat[
      as.character(cosine_df$merged_group[i]),
      as.character(cosine_df$original_cluster[i])
    ] <- as.numeric(cosine_df$signed_lfc_cosine[i])
  }

  mat
}

write_matrix_csv <- function(mat, file_path, row_id = "row_id") {
  df <- as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column(row_id)
  write_table_csv(df, file_path)
}

plot_cosine_heatmap <- function(mat, title_text, file_stub) {
  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("merged_group", "original_cluster", "signed_lfc_cosine")
  plot_df$merged_group <- factor(plot_df$merged_group, levels = rev(rownames(mat)))
  plot_df$original_cluster <- factor(plot_df$original_cluster, levels = colnames(mat))

  p <- ggplot(plot_df, aes(x = original_cluster, y = merged_group, fill = signed_lfc_cosine)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(is.na(signed_lfc_cosine), "NA", sprintf("%.2f", signed_lfc_cosine))),
      size = 3
    ) +
    scale_fill_gradient2(
      low = "#1f4e79",
      mid = "white",
      high = "#b03a2e",
      midpoint = 0,
      limits = c(-1, 1),
      na.value = "grey90"
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = title_text,
      x = "Original cluster",
      y = "Merged group",
      fill = "Cosine"
    )

  save_plot_pdf_png(
    p,
    file_stub,
    width = max(8, 0.8 * ncol(mat) + 3),
    height = max(4, 0.8 * nrow(mat) + 2)
  )
}

extract_features_for_dotplot <- function(top_marker_df, group_col = "group", n_per_group = 5L) {
  top_marker_df %>%
    group_by(.data[[group_col]]) %>%
    arrange(rank_within_group, .by_group = TRUE) %>%
    slice_head(n = n_per_group) %>%
    ungroup() %>%
    pull(gene) %>%
    unique()
}

plot_dotplot_safe <- function(obj, features, group.by, title_text, file_stub) {
  features <- unique(features[!is.na(features) & features != ""])
  if (length(features) == 0) {
    message("Skipping dotplot for ", group.by, " because no features were selected.")
    return(invisible(NULL))
  }
  p <- DotPlot(obj, features = features, group.by = group.by, assay = "RNA") +
    RotatedAxis() +
    labs(title = title_text)
  save_plot_pdf_png(p, file_stub, width = max(8, 0.32 * length(features) + 3), height = 6)
  invisible(p)
}

.ensure_dir(output_root)
out_plots <- .ensure_dir(file.path(output_root, "plots"))
out_markers_original <- .ensure_dir(file.path(output_root, "markers_original"))
out_markers_merged <- .ensure_dir(file.path(output_root, "markers_merged"))
out_markers_internal <- .ensure_dir(file.path(output_root, "markers_internal_pairs"))
out_summary <- .ensure_dir(file.path(output_root, "summaries"))
out_objects <- .ensure_dir(file.path(output_root, "objects"))

  if (!file.exists(input_rds)) {
    stop("Input RDS does not exist: ", input_rds, call. = FALSE)
  }

  configure_future_for_seurat(max_size_gb = 30)

  message("Reading Seurat object: ", input_rds)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input file is not a Seurat object.", call. = FALSE)
  }

  assert_required_column(obj@meta.data, base_cluster_col)
  if (!("umap" %in% names(obj@reductions))) {
    stop("UMAP reduction not found in Seurat object.", call. = FALSE)
  }

  obj <- prepare_obj_for_markers(obj)

  cluster_vec <- as.character(obj@meta.data[[base_cluster_col]])
  cluster_levels <- if (is.factor(obj@meta.data[[base_cluster_col]])) {
    as.character(levels(obj@meta.data[[base_cluster_col]]))
  } else {
    sort_maybe_numeric(cluster_vec)
  }
  obj@meta.data[[base_cluster_col]] <- factor(cluster_vec, levels = cluster_levels)
  obj@meta.data[["manual_merge_test"]] <- build_manual_merge_labels(cluster_vec, manual_merge_map, cluster_order = cluster_levels)

  merge_map_df <- bind_rows(lapply(names(manual_merge_map), function(x) {
    data.frame(
      merged_group = x,
      source_cluster = as.character(manual_merge_map[[x]]),
      stringsAsFactors = FALSE
    )
  }))
  write_table_csv(merge_map_df, file.path(out_summary, "manual_merge_mapping.csv"))

  cluster_count_df <- data.frame(
    base_cluster = names(table(obj@meta.data[[base_cluster_col]])),
    n_cells = as.integer(table(obj@meta.data[[base_cluster_col]])),
    stringsAsFactors = FALSE
  )
  merged_count_df <- data.frame(
    manual_merge_test = names(table(obj@meta.data[["manual_merge_test"]])),
    n_cells = as.integer(table(obj@meta.data[["manual_merge_test"]])),
    stringsAsFactors = FALSE
  )
  write_table_csv(cluster_count_df, file.path(out_summary, "original_cluster_cell_counts.csv"))
  write_table_csv(merged_count_df, file.path(out_summary, "manual_merge_group_cell_counts.csv"))

  message("Writing UMAP plots.")
  p_umap_original <- DimPlot(obj, reduction = "umap", group.by = base_cluster_col, label = FALSE, repel = FALSE, pt.size = 0.30, raster = FALSE) +
    labs(title = paste("Base clusters:", base_cluster_col))
  save_plot_pdf_png(p_umap_original, file.path(out_plots, "umap_base_clusters"), width = 9, height = 7)

  p_umap_merged <- DimPlot(obj, reduction = "umap", group.by = "manual_merge_test", label = FALSE, repel = FALSE, pt.size = 0.30, raster = FALSE) +
    labs(title = "Manual Merge Test: 0/1/7 -> 0 and 10/11/12 -> 10")
  save_plot_pdf_png(p_umap_merged, file.path(out_plots, "umap_manual_merge_test"), width = 9, height = 7)

  focus_clusters <- unique(unlist(manual_merge_map, use.names = FALSE))
  focus_cells <- rownames(obj@meta.data)[as.character(obj@meta.data[[base_cluster_col]]) %in% focus_clusters]
  focus_obj <- subset(obj, cells = focus_cells)

  p_focus_original <- DimPlot(focus_obj, reduction = "umap", group.by = base_cluster_col, label = TRUE, repel = TRUE, pt.size = 0.35, raster = FALSE) +
    labs(title = "Focus UMAP: refined clusters 0, 1, 7, 10, 11, 12")
  save_plot_pdf_png(p_focus_original, file.path(out_plots, "umap_focus_original_clusters"), width = 8, height = 6)

  p_focus_merged <- DimPlot(focus_obj, reduction = "umap", group.by = "manual_merge_test", label = TRUE, repel = TRUE, pt.size = 0.35, raster = FALSE) +
    labs(title = "Focus UMAP: manual merged groups")
  save_plot_pdf_png(p_focus_merged, file.path(out_plots, "umap_focus_manual_merge_test"), width = 8, height = 6)

  message("Running original-cluster markers.")
  original_marker_groups <- c("0", "1", "7", "10", "11", "12")
  original_markers <- run_markers_vs_rest(
    obj = obj,
    ident_col = base_cluster_col,
    groups = original_marker_groups,
    out_dir = out_markers_original,
    top_n = top_markers_n
  )
  write_table_csv(original_markers$full, file.path(out_markers_original, "markers_original_focus_clusters_all.csv"))
  write_table_csv(original_markers$top_positive, file.path(out_markers_original, "markers_original_focus_clusters_top_positive.csv"))

  message("Running merged-group markers.")
  merged_marker_groups <- names(manual_merge_map)
  merged_markers <- run_markers_vs_rest(
    obj = obj,
    ident_col = "manual_merge_test",
    groups = merged_marker_groups,
    out_dir = out_markers_merged,
    top_n = top_markers_n
  )
  write_table_csv(merged_markers$full, file.path(out_markers_merged, "markers_manual_merged_groups_all.csv"))
  write_table_csv(merged_markers$top_positive, file.path(out_markers_merged, "markers_manual_merged_groups_top_positive.csv"))

  message("Running internal pairwise marker checks.")
  pair_df <- data.frame(
    group_1 = c("0", "0", "1", "10", "10", "11"),
    group_2 = c("1", "7", "7", "11", "12", "12"),
    stringsAsFactors = FALSE
  )
  internal_pairwise <- run_pairwise_markers(
    obj = focus_obj,
    ident_col = base_cluster_col,
    pair_df = pair_df,
    out_dir = out_markers_internal
  )

  message("Computing top50_each signed-LFC cosine similarities.")
  cosine_feature_spec <- list(
    padj_max = 0.05,
    abs_log2fc_min = 0.25,
    abs_delta_pct_min = 0.05,
    top_n_up = cosine_top_n_each,
    top_n_down = cosine_top_n_each
  )
  cosine_summary <- build_manual_merge_cosine_summary(
    original_markers = original_markers,
    merged_markers = merged_markers,
    original_groups = original_marker_groups,
    merged_groups = merged_marker_groups,
    merge_map = manual_merge_map,
    spec = cosine_feature_spec
  )
  write_table_csv(cosine_summary, file.path(out_summary, "manual_merge_vs_original_top50_each_cosine_similarity.csv"))

  cosine_mat <- cosine_summary_to_matrix(
    cosine_df = cosine_summary,
    merged_groups = merged_marker_groups,
    original_groups = original_marker_groups
  )
  write_matrix_csv(
    cosine_mat,
    file.path(out_summary, "manual_merge_vs_original_top50_each_cosine_matrix.csv"),
    row_id = "merged_group"
  )
  plot_cosine_heatmap(
    mat = cosine_mat,
    title_text = "Merged groups vs original clusters | top50_each signed-LFC cosine",
    file_stub = file.path(out_plots, "heatmap_manual_merge_vs_original_top50_each_cosine")
  )

  message("Writing dotplots.")
  original_dotplot_features <- extract_features_for_dotplot(
    top_marker_df = original_markers$top_positive,
    group_col = "group",
    n_per_group = top_dotplot_n
  )
  plot_dotplot_safe(
    obj = focus_obj,
    features = original_dotplot_features,
    group.by = base_cluster_col,
    title_text = "Top refined-cluster markers for clusters 0, 1, 7, 10, 11, 12",
    file_stub = file.path(out_plots, "dotplot_original_focus_clusters")
  )

  merged_dotplot_features <- extract_features_for_dotplot(
    top_marker_df = merged_markers$top_positive,
    group_col = "group",
    n_per_group = top_dotplot_n
  )
  focus_merged_obj <- subset(obj, subset = manual_merge_test %in% merged_marker_groups)
  plot_dotplot_safe(
    obj = focus_merged_obj,
    features = merged_dotplot_features,
    group.by = "manual_merge_test",
    title_text = "Top merged-group markers after manual merge test",
    file_stub = file.path(out_plots, "dotplot_manual_merged_groups")
  )

  message("Saving annotated Seurat object.")
  saveRDS(obj, file.path(out_objects, "integrated_sct_cca_seurat_cluster_refine_manual_merge_test.rds"))

  summary_lines <- c(
    "Manual cluster merge test completed.",
    paste0("Input object: ", input_rds),
    paste0("Output root: ", output_root),
    paste0("Base cluster column: ", base_cluster_col),
    paste0("Merged groups tested: ", paste(names(manual_merge_map), collapse = ", ")),
    paste0("Marker settings: min.pct=", marker_min_pct, ", logfc.threshold=", marker_logfc_threshold, ", test.use=", marker_test_use),
    paste0(
      "Cosine feature rule: padj < ", cosine_feature_spec$padj_max,
      ", abs(log2FC) >= ", cosine_feature_spec$abs_log2fc_min,
      ", abs(pct.1 - pct.2) >= ", cosine_feature_spec$abs_delta_pct_min,
      ", top ", cosine_feature_spec$top_n_up, " up + top ", cosine_feature_spec$top_n_down, " down genes"
    ),
    "",
    "Merged cluster membership:",
    paste0("  ", names(manual_merge_map), ": ", vapply(manual_merge_map, function(x) paste(x, collapse = ", "), character(1))),
    "",
    "Cosine outputs:",
    "  summaries/manual_merge_vs_original_top50_each_cosine_similarity.csv",
    "  summaries/manual_merge_vs_original_top50_each_cosine_matrix.csv",
    "  plots/heatmap_manual_merge_vs_original_top50_each_cosine.pdf(.png)"
  )
  writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("Manual merge test finished.")
message("Outputs written to: ", normalizePath(output_root))
