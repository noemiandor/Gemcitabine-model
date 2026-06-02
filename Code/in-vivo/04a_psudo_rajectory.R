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
  "monocle3",
  "SingleCellExperiment",
  "dplyr",
  "digest",
  "ggplot2",
  "readr",
  "tibble",
  "Matrix",
  "igraph"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SingleCellExperiment)
  library(dplyr)
  library(digest)
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

parse_int_cfg <- function(value, default, name) {
  if (is.null(value) || !nzchar(trimws(as.character(value)))) return(default)
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) {
    stop(name, " must be a positive integer.", call. = FALSE)
  }
  out
}

resolve_first_existing <- function(paths) {
  paths <- paths[!is.na(paths) & nzchar(paths)]
  hit <- paths[file.exists(paths)]
  if (length(hit) > 0) normalizePath(hit[1], mustWork = TRUE) else paths[1]
}

median_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  stats::median(x)
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
      strip.text = element_text(face = "bold"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

save_both <- function(plot_obj, file_stub, width = 9, height = 7, dpi = 300) {
  ggsave(paste0(file_stub, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(file_stub, ".png"), plot_obj, width = width, height = height, dpi = dpi)
}

add_analysis_columns <- function(df, group_col, group_label, analysis_name = NULL) {
  group_col <- as.character(group_col)[1]
  if (is.null(group_col) || !nzchar(group_col)) group_col <- "AnalysisGroup"
  group_label <- as.character(group_label)[1]
  n <- nrow(df)
  df[[group_col]] <- rep(group_label, n)
  df$AnalysisGroup <- rep(group_label, n)
  df$AnalysisGroupColumn <- rep(group_col, n)
  if (!is.null(analysis_name) && nzchar(as.character(analysis_name)[1])) {
    df$Analysis <- rep(as.character(analysis_name)[1], n)
  }
  df
}

constant_metadata_values <- function(df, metadata_cols) {
  metadata_cols <- metadata_cols[metadata_cols %in% colnames(df)]
  out <- stats::setNames(vector("list", length(metadata_cols)), metadata_cols)
  for (col in metadata_cols) {
    vals <- unique(as.character(df[[col]]))
    vals <- vals[!is.na(vals) & nzchar(vals)]
    out[[col]] <- if (length(vals) == 1) vals else if (length(vals) == 0) NA_character_ else "mixed"
  }
  out
}

append_constant_metadata <- function(df, metadata_values, skip_cols = character(0)) {
  for (col in setdiff(names(metadata_values), skip_cols)) {
    if (!(col %in% colnames(df))) {
      df[[col]] <- rep(metadata_values[[col]], nrow(df))
    }
  }
  df
}

write_group_counts <- function(values, col_name, out_file, order_levels = NULL) {
  values_chr <- as.character(values)
  values_chr[is.na(values_chr) | !nzchar(values_chr)] <- NA_character_
  tab <- table(values_chr, useNA = "ifany")
  count_df <- data.frame(
    group = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
  names(count_df)[1] <- col_name
  if (!is.null(order_levels)) {
    order_levels <- as.character(order_levels)
    count_df <- count_df %>%
      dplyr::mutate(.order = match(.data[[col_name]], order_levels)) %>%
      dplyr::arrange(is.na(.order), .order, .data[[col_name]]) %>%
      dplyr::select(-.order)
  }
  write_table_csv(count_df, out_file)
  invisible(count_df)
}

get_expression_matrix_for_monocle3 <- function(obj, assay = "RNA") {
  mat <- get_assay_matrix(obj, assay = assay, slot_name = "counts")
  matrix_source <- "counts"
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    mat <- get_assay_matrix(obj, assay = assay, slot_name = "data")
    matrix_source <- "data"
  }
  if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Cannot find counts or data matrix for assay: ", assay, call. = FALSE)
  }
  list(matrix = mat, source = matrix_source)
}

choose_root_cells <- function(cds, cluster_col, root_cluster = NULL) {
  meta <- as.data.frame(colData(cds))
  clusters <- as.character(meta[[cluster_col]])
  cluster_levels <- sort_maybe_numeric(unique(clusters))
  root_cluster_use <- if (!is.null(root_cluster) && root_cluster %in% clusters) {
    root_cluster
  } else {
    cluster_levels[1]
  }
  root_cells <- rownames(meta)[clusters == root_cluster_use]
  if (length(root_cells) == 0) {
    root_cells <- colnames(cds)[1]
  }
  list(root_cluster = root_cluster_use, root_cells = root_cells)
}

normalize_graph_vertex_ids <- function(vertex_ids, graph_nodes) {
  out <- as.character(vertex_ids)
  out[!nzchar(out) | out %in% c("NA", "NaN")] <- NA_character_

  missing <- which(!is.na(out) & !(out %in% graph_nodes))
  if (length(missing) > 0) {
    idx <- suppressWarnings(as.integer(out[missing]))
    valid_idx <- !is.na(idx) & idx >= 1L & idx <= length(graph_nodes)
    out[missing[valid_idx]] <- graph_nodes[idx[valid_idx]]
  }

  missing <- which(!is.na(out) & !(out %in% graph_nodes))
  if (length(missing) > 0) {
    prefixed <- paste0("Y_", out[missing])
    valid_prefixed <- prefixed %in% graph_nodes
    out[missing[valid_prefixed]] <- prefixed[valid_prefixed]
  }

  out[!(out %in% graph_nodes)] <- NA_character_
  out
}

extract_closest_graph_vertex <- function(closest, cell_names, graph_nodes) {
  out <- stats::setNames(rep(NA_character_, length(cell_names)), cell_names)
  if (is.null(closest)) return(out)

  closest_vec <- if (is.matrix(closest) || is.data.frame(closest)) {
    vec <- as.character(closest[, 1])
    row_ids <- rownames(closest)
    if (!is.null(row_ids) && length(row_ids) == length(vec)) {
      names(vec) <- row_ids
    }
    vec
  } else {
    as.character(closest)
  }

  if (is.null(names(closest_vec)) && length(closest_vec) == length(cell_names)) {
    names(closest_vec) <- cell_names
  }
  closest_vec <- normalize_graph_vertex_ids(closest_vec, graph_nodes)

  common_cells <- intersect(cell_names, names(closest_vec))
  if (length(common_cells) > 0) {
    out[common_cells] <- closest_vec[common_cells]
  }
  out
}

choose_root_graph_node <- function(closest_vec, root_cells, graph_nodes) {
  root_nodes <- closest_vec[intersect(root_cells, names(closest_vec))]
  root_nodes <- root_nodes[!is.na(root_nodes) & root_nodes %in% graph_nodes]
  if (length(root_nodes) == 0) return(graph_nodes[1])

  root_counts <- sort(table(root_nodes), decreasing = TRUE)
  names(root_counts)[1]
}

compute_root_graph_distances <- function(graph, graph_nodes, root_graph_node) {
  out <- stats::setNames(rep(Inf, length(graph_nodes)), graph_nodes)
  if (!(root_graph_node %in% graph_nodes)) {
    root_graph_node <- graph_nodes[1]
  }
  dist_mat <- igraph::distances(graph, v = root_graph_node, to = graph_nodes, weights = NA)
  out[colnames(dist_mat)] <- as.numeric(dist_mat[1, ])
  out
}

extract_monocle_graph_tables <- function(cds, group_label, root_cells, root_cluster, group_col = "Dose", analysis_name = NULL, metadata_cols = character(0)) {
  empty_nodes <- data.frame(
    graph_node = character(0),
    UMAP_1 = numeric(0),
    UMAP_2 = numeric(0),
    graph_node_pseudotime_raw = numeric(0),
    graph_node_root_distance = numeric(0),
    root_graph_node = character(0),
    root_cluster = character(0),
    n_cells_projected = integer(0),
    stringsAsFactors = FALSE
  )
  empty_edges <- data.frame(
    from = character(0),
    to = character(0),
    UMAP_1_from = numeric(0),
    UMAP_2_from = numeric(0),
    UMAP_1_to = numeric(0),
    UMAP_2_to = numeric(0),
    from_pseudotime_raw = numeric(0),
    to_pseudotime_raw = numeric(0),
    from_root_distance = numeric(0),
    to_root_distance = numeric(0),
    root_graph_node = character(0),
    root_cluster = character(0),
    stringsAsFactors = FALSE
  )
  meta_constants <- constant_metadata_values(as.data.frame(colData(cds)), metadata_cols)
  empty_nodes <- append_constant_metadata(empty_nodes, meta_constants, skip_cols = group_col)
  empty_edges <- append_constant_metadata(empty_edges, meta_constants, skip_cols = group_col)
  empty_nodes <- add_analysis_columns(empty_nodes, group_col, group_label, analysis_name)
  empty_edges <- add_analysis_columns(empty_edges, group_col, group_label, analysis_name)

  graph_list <- tryCatch(monocle3::principal_graph(cds), error = function(e) NULL)
  graph_aux <- tryCatch(monocle3::principal_graph_aux(cds), error = function(e) NULL)
  if (is.null(graph_list) || is.null(graph_list[["UMAP"]]) || is.null(graph_aux) || is.null(graph_aux[["UMAP"]])) {
    return(list(nodes = empty_nodes, edges = empty_edges))
  }

  graph <- graph_list[["UMAP"]]
  coords <- graph_aux[["UMAP"]]$dp_mst
  if (is.null(coords) || ncol(coords) == 0) {
    return(list(nodes = empty_nodes, edges = empty_edges))
  }

  coords <- t(as.matrix(coords))
  coords <- coords[, seq_len(min(2, ncol(coords))), drop = FALSE]
  if (ncol(coords) < 2) {
    return(list(nodes = empty_nodes, edges = empty_edges))
  }
  colnames(coords) <- c("UMAP_1", "UMAP_2")
  graph_nodes <- rownames(coords)
  node_df <- data.frame(
    graph_node = graph_nodes,
    UMAP_1 = coords[, "UMAP_1"],
    UMAP_2 = coords[, "UMAP_2"],
    stringsAsFactors = FALSE
  )

  closest <- graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vec <- extract_closest_graph_vertex(closest, colnames(cds), graph_nodes)
  root_graph_node <- choose_root_graph_node(closest_vec, root_cells, graph_nodes)
  root_dist <- compute_root_graph_distances(graph, graph_nodes, root_graph_node)

  pt <- monocle3::pseudotime(cds)
  if (is.null(names(pt))) {
    names(pt) <- colnames(cds)
  }
  cell_node_df <- data.frame(
    cell = names(pt),
    graph_node = closest_vec[names(pt)],
    pseudotime_raw = as.numeric(pt),
    stringsAsFactors = FALSE
  )
  node_pt <- cell_node_df %>%
    dplyr::filter(!is.na(graph_node)) %>%
    dplyr::group_by(graph_node) %>%
    dplyr::summarise(
      graph_node_pseudotime_raw = median_finite(pseudotime_raw),
      n_cells_projected = dplyr::n(),
      .groups = "drop"
    )
  node_df <- node_df %>%
    dplyr::left_join(node_pt, by = "graph_node")
  node_df$graph_node_root_distance <- unname(root_dist[node_df$graph_node])
  node_df$root_graph_node <- root_graph_node
  node_df$root_cluster <- root_cluster
  node_df <- append_constant_metadata(node_df, meta_constants, skip_cols = group_col)
  node_df <- add_analysis_columns(node_df, group_col, group_label, analysis_name)

  edge_mat <- igraph::as_edgelist(graph, names = TRUE)
  if (nrow(edge_mat) == 0) {
    return(list(nodes = node_df, edges = empty_edges))
  }

  pt_lookup <- stats::setNames(node_df$graph_node_pseudotime_raw, node_df$graph_node)
  dist_lookup <- stats::setNames(node_df$graph_node_root_distance, node_df$graph_node)
  raw_edges <- data.frame(
    from_raw = edge_mat[, 1],
    to_raw = edge_mat[, 2],
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      from_pt = unname(pt_lookup[from_raw]),
      to_pt = unname(pt_lookup[to_raw]),
      from_dist = unname(dist_lookup[from_raw]),
      to_dist = unname(dist_lookup[to_raw]),
      reverse_edge = is.finite(from_dist) & is.finite(to_dist) & from_dist > to_dist,
      from = ifelse(reverse_edge, to_raw, from_raw),
      to = ifelse(reverse_edge, from_raw, to_raw),
      from_pseudotime_raw = ifelse(reverse_edge, to_pt, from_pt),
      to_pseudotime_raw = ifelse(reverse_edge, from_pt, to_pt),
      from_root_distance = ifelse(reverse_edge, to_dist, from_dist),
      to_root_distance = ifelse(reverse_edge, from_dist, to_dist)
    )

  coord_lookup <- node_df %>% dplyr::select(graph_node, UMAP_1, UMAP_2)
  edge_df <- raw_edges %>%
    dplyr::select(from, to, from_pseudotime_raw, to_pseudotime_raw, from_root_distance, to_root_distance) %>%
    dplyr::left_join(coord_lookup, by = c("from" = "graph_node")) %>%
    dplyr::rename(UMAP_1_from = UMAP_1, UMAP_2_from = UMAP_2) %>%
    dplyr::left_join(coord_lookup, by = c("to" = "graph_node")) %>%
    dplyr::rename(UMAP_1_to = UMAP_1, UMAP_2_to = UMAP_2) %>%
    dplyr::mutate(
      root_graph_node = root_graph_node,
      root_cluster = root_cluster
    )
  edge_df <- append_constant_metadata(edge_df, meta_constants, skip_cols = group_col)
  edge_df <- add_analysis_columns(edge_df, group_col, group_label, analysis_name)

  list(nodes = node_df, edges = edge_df)
}

plot_dose_trajectory <- function(cell_df, cluster_df, edge_df, cluster_col, dose_label, out_dir, cds = NULL) {
  p_time <- ggplot(cell_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 1.5, alpha = 0.92, stroke = 0) +
    geom_segment(
      data = edge_df,
      aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed"),
      linewidth = 0.5,
      color = "black",
      alpha = 1
    ) +
    geom_label(
      data = cluster_df,
      aes(x = UMAP_1, y = UMAP_2, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      linewidth = 0.15,
      fill = "white",
      alpha = 0.9
    ) +
    coord_equal() +
    scale_color_gradientn(colors = c("#1B2A41", "#216869", "#F2C14E", "#C43B3B"), na.value = "grey85") +
    labs(
      title = paste0(dose_label, " Monocle3 pseudotime"),
      subtitle = "Trajectory graph is learned by monocle3 from expression data",
      color = "Pseudotime"
    ) +
    make_pretty_umap_theme()
  save_both(p_time, file.path(out_dir, "umap_monocle3_pseudotime"), width = 8.5, height = 7)

  cluster_levels <- sort_maybe_numeric(as.character(cell_df[[cluster_col]]))
  cluster_palette <- stats::setNames(grDevices::hcl.colors(length(cluster_levels), palette = "Dark 3"), cluster_levels)
  p_cluster <- ggplot(cell_df, aes(x = UMAP_1, y = UMAP_2, color = as.character(.data[[cluster_col]]))) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    geom_segment(
      data = edge_df,
      aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed"),
      linewidth = 0.5,
      color = "black",
      alpha = 1
    ) +
    geom_label(
      data = cluster_df,
      aes(x = UMAP_1, y = UMAP_2, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      linewidth = 0.15,
      fill = "white",
      alpha = 0.9
    ) +
    coord_equal() +
    scale_color_manual(values = cluster_palette, na.value = "grey80") +
    labs(title = paste0(dose_label, " Monocle3 trajectory by cluster"), color = cluster_col) +
    make_pretty_umap_theme()
  save_both(p_cluster, file.path(out_dir, "umap_monocle3_cluster"), width = 8.5, height = 7)

  p_dist <- ggplot(cell_df, aes(x = factor(.data[[cluster_col]], levels = cluster_levels), y = pseudotime, fill = as.character(.data[[cluster_col]]))) +
    geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, na.rm = TRUE) +
    geom_boxplot(width = 0.13, outlier.size = 0.15, alpha = 0.7, na.rm = TRUE) +
    scale_fill_manual(values = cluster_palette, guide = "none") +
    labs(title = paste0(dose_label, " pseudotime distribution"), x = cluster_col, y = "Pseudotime") +
    theme_classic(base_size = 11)
  save_both(p_dist, file.path(out_dir, "pseudotime_by_cluster"), width = max(7, 0.45 * length(cluster_levels) + 4), height = 5.5)

  if (!is.null(cds)) {
    p_monocle <- tryCatch(
      monocle3::plot_cells(
        cds,
        color_cells_by = "pseudotime",
        label_cell_groups = FALSE,
        label_leaves = TRUE,
        label_branch_points = TRUE,
        graph_label_size = 2,
        trajectory_graph_color = "black",
        trajectory_graph_segment_size = 0.5,
        cell_size = 1.65,
        cell_stroke = 0,
        alpha = 0.95
      ),
      error = function(e) NULL
    )
    if (!is.null(p_monocle)) {
      save_both(p_monocle, file.path(out_dir, "monocle3_plot_cells_pseudotime"), width = 8.5, height = 7)
    }
  }

  invisible(list(pseudotime = p_time, cluster = p_cluster, distribution = p_dist))
}

plot_metadata_trajectory <- function(cell_df, cluster_df, edge_df, cluster_col, color_col, title_label, out_dir, file_stub) {
  if (!(color_col %in% colnames(cell_df))) return(invisible(NULL))
  color_values <- as.character(cell_df[[color_col]])
  color_levels <- if (is.factor(cell_df[[color_col]])) {
    as.character(levels(cell_df[[color_col]]))
  } else {
    unique(color_values)
  }
  color_levels <- color_levels[!is.na(color_levels) & nzchar(color_levels) & color_levels %in% unique(color_values)]
  if (length(color_levels) == 0) return(invisible(NULL))
  color_palette <- stats::setNames(grDevices::hcl.colors(length(color_levels), palette = "Dark 3"), color_levels)

  p_umap <- ggplot(cell_df, aes(x = UMAP_1, y = UMAP_2, color = factor(.data[[color_col]], levels = color_levels))) +
    geom_point(size = 1.5, alpha = 0.9, stroke = 0) +
    geom_segment(
      data = edge_df,
      aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.11, "inches"), type = "closed"),
      linewidth = 0.5,
      color = "black",
      alpha = 1
    ) +
    geom_label(
      data = cluster_df,
      aes(x = UMAP_1, y = UMAP_2, label = cluster),
      inherit.aes = FALSE,
      size = 3,
      linewidth = 0.15,
      fill = "white",
      alpha = 0.9
    ) +
    coord_equal() +
    scale_color_manual(values = color_palette, na.value = "grey80") +
    labs(title = paste0(title_label, " Monocle3 trajectory by ", color_col), color = color_col) +
    make_pretty_umap_theme()
  save_both(p_umap, file.path(out_dir, file_stub), width = 8.5, height = 7)

  p_dist <- ggplot(cell_df, aes(x = factor(.data[[color_col]], levels = color_levels), y = pseudotime, fill = factor(.data[[color_col]], levels = color_levels))) +
    geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, na.rm = TRUE) +
    geom_boxplot(width = 0.14, outlier.size = 0.12, alpha = 0.75, na.rm = TRUE) +
    scale_fill_manual(values = color_palette, guide = "none") +
    labs(title = paste0(title_label, " pseudotime distribution by ", color_col), x = color_col, y = "Pseudotime") +
    theme_classic(base_size = 11)
  save_both(p_dist, file.path(out_dir, paste0(file_stub, "_pseudotime_distribution")), width = max(6.5, 1.2 * length(color_levels) + 4), height = 5)

  p_density <- ggplot(cell_df, aes(x = pseudotime, fill = factor(.data[[color_col]], levels = color_levels))) +
    geom_density(alpha = 0.45, linewidth = 0.4, na.rm = TRUE) +
    scale_fill_manual(values = color_palette, name = color_col) +
    labs(title = paste0(title_label, " pseudotime density by ", color_col), x = "Pseudotime", y = "Density") +
    theme_classic(base_size = 11)
  save_both(p_density, file.path(out_dir, paste0(file_stub, "_pseudotime_density")), width = 7.5, height = 5)

  invisible(list(umap = p_umap, distribution = p_dist, density = p_density))
}

plot_group_summary <- function(cell_all, graph_edge_all, group_col, group_title, out_dir, file_suffix) {
  if (nrow(cell_all) == 0 || !(group_col %in% colnames(cell_all))) return(invisible(NULL))
  group_values <- as.character(cell_all[[group_col]])
  group_levels <- if (is.factor(cell_all[[group_col]])) {
    as.character(levels(cell_all[[group_col]]))
  } else {
    unique(group_values)
  }
  group_levels <- group_levels[!is.na(group_levels) & nzchar(group_levels) & group_levels %in% unique(group_values)]
  if (length(group_levels) == 0) return(invisible(NULL))
  cell_all[[group_col]] <- factor(as.character(cell_all[[group_col]]), levels = group_levels)
  if (!(group_col %in% colnames(graph_edge_all))) {
    graph_edge_all[[group_col]] <- character(nrow(graph_edge_all))
  }
  graph_edge_all[[group_col]] <- factor(as.character(graph_edge_all[[group_col]]), levels = group_levels)
  use_dose_karyotype_grid <- identical(group_col, "Dose_Karyotype") &&
    all(c("Dose", "Karyotype") %in% colnames(cell_all)) &&
    all(c("Dose", "Karyotype") %in% colnames(graph_edge_all))
  if (use_dose_karyotype_grid) {
    dose_order <- c("0mg/kg", "30mg/kg", "120mg/kg")
    dose_extra <- sort(setdiff(unique(as.character(cell_all$Dose)), dose_order))
    dose_levels <- c(dose_order, dose_extra)
    dose_levels <- dose_levels[dose_levels %in% unique(as.character(cell_all$Dose))]
    karyotype_order <- c("2N", "4N")
    karyotype_extra <- sort(setdiff(unique(as.character(cell_all$Karyotype)), karyotype_order))
    karyotype_levels <- c(karyotype_order, karyotype_extra)
    karyotype_levels <- karyotype_levels[karyotype_levels %in% unique(as.character(cell_all$Karyotype))]
    cell_all$Dose <- factor(as.character(cell_all$Dose), levels = dose_levels)
    cell_all$Karyotype <- factor(as.character(cell_all$Karyotype), levels = karyotype_levels)
    graph_edge_all$Dose <- factor(as.character(graph_edge_all$Dose), levels = dose_levels)
    graph_edge_all$Karyotype <- factor(as.character(graph_edge_all$Karyotype), levels = karyotype_levels)
    facet_layer <- facet_grid(Karyotype ~ Dose)
    plot_width <- 12
    plot_height <- 8
  } else {
    facet_formula <- stats::as.formula(paste("~", group_col))
    facet_layer <- facet_wrap(facet_formula)
    plot_width <- max(8, min(22, 3.8 * length(group_levels)))
    plot_height <- if (length(group_levels) > 4) 10 else 7.5
  }

  p_all <- ggplot(cell_all, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
    geom_point(size = 1.26, alpha = 0.94, stroke = 0) +
    geom_segment(
      data = graph_edge_all,
      aes(x = UMAP_1_from, y = UMAP_2_from, xend = UMAP_1_to, yend = UMAP_2_to),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.09, "inches"), type = "closed"),
      linewidth = 0.5,
      color = "black",
      alpha = 1
    ) +
    facet_layer +
    coord_equal() +
    scale_color_gradientn(colors = c("#1B2A41", "#216869", "#F2C14E", "#C43B3B"), na.value = "grey85") +
    labs(
      title = paste0("Monocle3 pseudotime by ", group_title),
      subtitle = paste0("Each ", group_title, " group is analyzed independently from the expression matrix"),
      color = "Pseudotime"
    ) +
    make_pretty_umap_theme()
  save_both(p_all, file.path(out_dir, paste0("umap_monocle3_pseudotime_by_", file_suffix)), width = plot_width, height = plot_height)

  p_density <- ggplot(cell_all, aes(x = pseudotime, fill = .data[[group_col]])) +
    geom_density(alpha = 0.45, linewidth = 0.4, na.rm = TRUE) +
    labs(title = paste0("Monocle3 pseudotime density by ", group_title), x = "Pseudotime", y = "Density", fill = group_title) +
    theme_classic(base_size = 11)
  save_both(p_density, file.path(out_dir, paste0("pseudotime_density_by_", file_suffix)), width = 8.5, height = 5)

  p_box <- ggplot(cell_all, aes(x = .data[[group_col]], y = pseudotime, fill = .data[[group_col]])) +
    geom_violin(scale = "width", linewidth = 0.2, alpha = 0.75, na.rm = TRUE) +
    geom_boxplot(width = 0.14, outlier.size = 0.12, alpha = 0.75, na.rm = TRUE) +
    labs(title = paste0("Monocle3 pseudotime distribution by ", group_title), x = NULL, y = "Pseudotime") +
    theme_classic(base_size = 11) +
    theme(legend.position = "none")
  save_both(p_box, file.path(out_dir, paste0("pseudotime_distribution_by_", file_suffix)), width = max(7, 0.8 * length(group_levels) + 4), height = 5)

  invisible(list(umap = p_all, density = p_density, distribution = p_box))
}

run_dose_monocle3 <- function(
  obj,
  dose_label,
  dose_cells,
  output_dir,
  assay,
  cluster_col,
  nfeatures,
  npcs,
  root_cluster = NULL,
  group_col = "Dose",
  analysis_name = "by_dose",
  tn_scope = NA_character_,
  ploidy_scope = "ploidy_all",
  branch_label = NA_character_
) {
  .ensure_dir(output_dir)
  sub_obj <- subset(obj, cells = dose_cells)
  sub_obj <- prepare_obj_for_markers(sub_obj, assay = assay)

  if (!(cluster_col %in% colnames(sub_obj@meta.data))) {
    stop("Cluster column is missing in group subset: ", cluster_col, call. = FALSE)
  }
  if (!("umap" %in% names(sub_obj@reductions))) {
    stop("UMAP reduction is missing in group subset: ", dose_label, call. = FALSE)
  }

  sub_obj <- FindVariableFeatures(sub_obj, assay = assay, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
  features <- VariableFeatures(sub_obj)
  expr_info <- get_expression_matrix_for_monocle3(sub_obj, assay = assay)
  expr_mat <- expr_info$matrix
  expr_mat <- expr_mat[, colnames(sub_obj), drop = FALSE]
  features <- unique(features[features %in% rownames(expr_mat)])
  if (length(features) < 20) {
    stop("Too few variable features for Monocle3 in group ", dose_label, ".", call. = FALSE)
  }
  expr_mat <- expr_mat[features, , drop = FALSE]

  cell_meta <- sub_obj@meta.data[colnames(expr_mat), , drop = FALSE]
  cell_meta[[cluster_col]] <- as.character(cell_meta[[cluster_col]])
  cell_meta$trajectory_analysis_group <- as.character(dose_label)
  cell_meta$trajectory_tn_scope <- as.character(tn_scope)
  cell_meta$trajectory_ploidy_scope <- as.character(ploidy_scope)
  cell_meta$trajectory_branch <- as.character(branch_label)
  cell_meta$trajectory_shape_group <- if (identical(as.character(dose_label), "All_cells")) {
    as.character(cell_meta$TN)
  } else if (identical(as.character(ploidy_scope), "ploidy_all")) {
    as.character(cell_meta$Ploidy)
  } else {
    as.character(ploidy_scope)
  }
  gene_meta <- data.frame(
    gene_short_name = rownames(expr_mat),
    row.names = rownames(expr_mat),
    stringsAsFactors = FALSE
  )

  cds <- monocle3::new_cell_data_set(
    expr_mat,
    cell_metadata = cell_meta,
    gene_metadata = gene_meta
  )

  npcs_use <- min(npcs, length(features) - 1L, ncol(cds) - 1L)
  if (npcs_use < 2L) {
    stop("Too few cells/features for Monocle3 PCA in group ", dose_label, ".", call. = FALSE)
  }
  cds <- monocle3::preprocess_cds(
    cds,
    method = "PCA",
    num_dim = npcs_use,
    norm_method = "log",
    use_genes = features
  )

  umap_mat <- as.matrix(Embeddings(sub_obj, "umap")[colnames(cds), seq_len(2), drop = FALSE])
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
  reducedDim(cds, "UMAP") <- umap_mat

  cds <- monocle3::cluster_cells(cds, reduction_method = "UMAP", verbose = FALSE)
  cds <- monocle3::learn_graph(cds, use_partition = FALSE, verbose = FALSE)

  root_info <- choose_root_cells(cds, cluster_col = cluster_col, root_cluster = root_cluster)
  cds <- monocle3::order_cells(
    cds,
    reduction_method = "UMAP",
    root_cells = root_info$root_cells
  )

  pt_raw <- monocle3::pseudotime(cds)
  if (is.null(names(pt_raw))) {
    names(pt_raw) <- colnames(cds)
  }
  pt_scaled <- safe_scale01(pt_raw)
  names(pt_scaled) <- names(pt_raw)
  pt_scaled[!is.finite(pt_raw)] <- NA_real_

  meta_cols <- c(
    "Dose",
    "Karyotype",
    "SampleOrigin",
    "Dose_Karyotype",
    "OriginDoseGroup",
    "TrajectorySet",
    "TrajectoryUniverse",
    "TrajectoryComparison",
    "TrajectoryComparisonID",
    "TrajectoryComparisonGroup",
    "TrajectoryOriginalGroup",
    "TrajectoryScope",
    "TrajectoryIdent1",
  "TrajectoryIdent2",
  "trajectory_context",
  "trajectory_group",
  "trajectory_analysis_group",
  "trajectory_shape_group",
  "trajectory_tn_scope",
  "trajectory_ploidy_scope",
  "trajectory_branch",
  "IDs",
  "clusters",
  "Ploidy",
    "TN",
    "Dose_DEG",
    cluster_col,
    "sample",
    "sample_type",
    "cluster_final_annotation_primary",
    "cluster_final_annotation_multi",
    "seurat_clusters",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  )
  meta_df <- as.data.frame(colData(cds))
  meta_cols <- unique(meta_cols[meta_cols %in% colnames(meta_df)])
  umap_df <- as.data.frame(reducedDim(cds, "UMAP"))
  colnames(umap_df) <- c("UMAP_1", "UMAP_2")
  umap_df$cell <- rownames(umap_df)

  cell_df <- meta_df[, meta_cols, drop = FALSE] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::left_join(umap_df, by = "cell") %>%
    dplyr::mutate(
      pseudotime_raw = as.numeric(pt_raw[cell]),
      pseudotime = as.numeric(pt_scaled[cell]),
      root_cluster = root_info$root_cluster,
      monocle3_matrix_source = expr_info$source
    )
  cell_df <- add_analysis_columns(cell_df, group_col, dose_label, analysis_name)

  cluster_values <- as.character(cell_df[[cluster_col]])
  cluster_levels <- if (is.factor(sub_obj@meta.data[[cluster_col]])) {
    as.character(levels(sub_obj@meta.data[[cluster_col]]))
  } else {
    sort_maybe_numeric(cluster_values)
  }
  cluster_levels <- cluster_levels[cluster_levels %in% unique(cluster_values)]
  cell_df[[cluster_col]] <- factor(cluster_values, levels = cluster_levels)

  cluster_df <- cell_df %>%
    dplyr::group_by(cluster = as.character(.data[[cluster_col]])) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      UMAP_1 = mean(UMAP_1, na.rm = TRUE),
      UMAP_2 = mean(UMAP_2, na.rm = TRUE),
      pseudotime_raw_median = median_finite(pseudotime_raw),
      pseudotime_median = median_finite(pseudotime),
      pseudotime_missing_cells = sum(!is.finite(pseudotime_raw) | is.na(pseudotime)),
      root_cluster = root_info$root_cluster,
      .groups = "drop"
    ) %>%
    dplyr::arrange(match(cluster, cluster_levels))
  metadata_cols_for_constants <- c(
    "Dose",
    "Karyotype",
    "SampleOrigin",
    "Dose_Karyotype",
    "OriginDoseGroup",
    "TrajectorySet",
    "TrajectoryUniverse",
    "TrajectoryComparison",
    "TrajectoryComparisonID",
    "TrajectoryScope",
    "TrajectoryIdent1",
    "TrajectoryIdent2",
    "trajectory_context",
    "trajectory_group",
    "trajectory_analysis_group",
    "trajectory_shape_group",
    "trajectory_tn_scope",
    "trajectory_ploidy_scope",
    "trajectory_branch"
  )
  meta_constants <- constant_metadata_values(cell_df, metadata_cols_for_constants)
  cluster_df <- append_constant_metadata(cluster_df, meta_constants, skip_cols = group_col)
  cluster_df <- add_analysis_columns(cluster_df, group_col, dose_label, analysis_name)

  graph_tables <- extract_monocle_graph_tables(
    cds = cds,
    group_label = dose_label,
    root_cells = root_info$root_cells,
    root_cluster = root_info$root_cluster,
    group_col = group_col,
    analysis_name = analysis_name,
    metadata_cols = metadata_cols_for_constants
  )
  graph_node_df <- graph_tables$nodes
  graph_edge_df <- graph_tables$edges

  readr::write_csv(cell_df, file.path(output_dir, "cells_pseudotime.csv"))
  readr::write_csv(cluster_df, file.path(output_dir, "cluster_pseudotime_summary.csv"))
  readr::write_csv(graph_node_df, file.path(output_dir, "principal_graph_nodes.csv"))
  readr::write_csv(graph_edge_df, file.path(output_dir, "principal_graph_edges.csv"))
  readr::write_csv(data.frame(feature = features, stringsAsFactors = FALSE), file.path(output_dir, "monocle3_variable_features.csv"))

  plot_dose_trajectory(
    cell_df = cell_df,
    cluster_df = cluster_df,
    edge_df = graph_edge_df,
    cluster_col = cluster_col,
    dose_label = dose_label,
    out_dir = output_dir,
    cds = cds
  )

  list(
    cells = cell_df,
    clusters = cluster_df,
    graph_nodes = graph_node_df,
    graph_edges = graph_edge_df,
    root_cluster = root_info$root_cluster,
    matrix_source = expr_info$source
  )
}

bind_result_table <- function(results, table_name) {
  dplyr::bind_rows(lapply(results, function(x) {
    if (!is.null(x[[table_name]])) x[[table_name]] else NULL
  }))
}

write_result_collection_outputs <- function(results, group_col, analysis_name, group_title, table_suffix, plot_suffix, out_summary, out_plots) {
  cell_all <- bind_result_table(results, "cells")
  cluster_all <- bind_result_table(results, "clusters")
  graph_node_all <- bind_result_table(results, "graph_nodes")
  graph_edge_all <- bind_result_table(results, "graph_edges")

  write_table_csv(cell_all, file.path(out_summary, paste0("pseudotime_cells_all_", table_suffix, ".csv")))
  write_table_csv(cluster_all, file.path(out_summary, paste0("cluster_pseudotime_summary_all_", table_suffix, ".csv")))
  write_table_csv(graph_node_all, file.path(out_summary, paste0("principal_graph_nodes_all_", table_suffix, ".csv")))
  write_table_csv(graph_edge_all, file.path(out_summary, paste0("principal_graph_edges_all_", table_suffix, ".csv")))

  plot_group_summary(
    cell_all = cell_all,
    graph_edge_all = graph_edge_all,
    group_col = group_col,
    group_title = group_title,
    out_dir = out_plots,
    file_suffix = plot_suffix
  )

  root_summary <- data.frame(
    root_cluster = vapply(results, `[[`, character(1), "root_cluster"),
    matrix_source = vapply(results, `[[`, character(1), "matrix_source"),
    stringsAsFactors = FALSE
  )
  root_summary[[group_col]] <- names(results)
  root_summary$AnalysisGroup <- names(results)
  root_summary$AnalysisGroupColumn <- group_col
  root_summary$Analysis <- analysis_name
  root_summary <- root_summary %>%
    dplyr::select(dplyr::all_of(group_col), AnalysisGroup, AnalysisGroupColumn, Analysis, root_cluster, matrix_source)
  write_table_csv(root_summary, file.path(out_summary, paste0("root_cluster_by_", plot_suffix, ".csv")))

  invisible(list(
    cells = cell_all,
    clusters = cluster_all,
    graph_nodes = graph_node_all,
    graph_edges = graph_edge_all,
    roots = root_summary
  ))
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

comparison_ident_label <- function(values) {
  values <- as.character(values)
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) == 0) return(NA_character_)
  paste(values, collapse = "+")
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
  safe <- sanitize_path_component(paste(pieces[!is.na(pieces) & nzchar(pieces)], collapse = "__"))
  safe[1]
}

make_trajectory_skip <- function(task_row, status, message, n_ident_1 = 0L, n_ident_2 = 0L, n_cells = 0L) {
  list(
    status = status,
    message = message,
    task = task_row,
    cells = character(0),
    group_by_cell = character(0),
    original_group_by_cell = character(0),
    ident_1_label = as.character(task_row$ident_1),
    ident_2_label = as.character(task_row$ident_2),
    n_ident_1 = as.integer(n_ident_1),
    n_ident_2 = as.integer(n_ident_2),
    n_cells = as.integer(n_cells),
    universe_key = NA_character_,
    comparison_stub = comparison_task_stub(task_row)
  )
}

build_trajectory_comparison_task <- function(obj, task_row, cluster_match_col, min_cells_per_group, min_total_cells) {
  meta <- obj@meta.data
  group_col <- as.character(task_row$group_col)[1]
  if (!(group_col %in% colnames(meta))) {
    return(make_trajectory_skip(task_row, "skipped", paste0("Missing group_col in metadata: ", group_col)))
  }

  base_cells <- rownames(meta)
  scope <- as.character(task_row$scope)[1]

  if (grepl("^within_cluster_", scope) && !is_missing_task_value(task_row$cluster)) {
    if (!(cluster_match_col %in% colnames(meta))) {
      return(make_trajectory_skip(task_row, "skipped", paste0("Missing cluster match column in metadata: ", cluster_match_col)))
    }
    cluster_values <- as.character(meta[[cluster_match_col]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(cluster_values) & cluster_values == as.character(task_row$cluster)])
  }

  if (!is_missing_task_value(task_row$subset_group)) {
    if (!("TN" %in% colnames(meta))) {
      return(make_trajectory_skip(task_row, "skipped", "Missing TN metadata for subset_group filter."))
    }
    tn_values <- as.character(meta[[ "TN" ]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(tn_values) & tn_values == as.character(task_row$subset_group)])
  }

  if (!is_missing_task_value(task_row$dose) && !identical(group_col, "Dose_DEG")) {
    if (!("Dose_DEG" %in% colnames(meta))) {
      return(make_trajectory_skip(task_row, "skipped", "Missing Dose_DEG metadata for dose filter."))
    }
    dose_values <- as.character(meta[[ "Dose_DEG" ]])
    base_cells <- intersect(base_cells, rownames(meta)[!is.na(dose_values) & dose_values == as.character(task_row$dose)])
  }

  base_cells <- intersect(base_cells, rownames(meta))
  if (length(base_cells) == 0) {
    return(make_trajectory_skip(task_row, "skipped", "No cells in comparison universe after scope filters."))
  }

  group_values <- as.character(meta[base_cells, group_col, drop = TRUE])
  valid_group <- !is.na(group_values) & nzchar(group_values)
  base_cells <- base_cells[valid_group]
  group_values <- group_values[valid_group]
  if (length(base_cells) == 0) {
    return(make_trajectory_skip(task_row, "skipped", paste0("No non-missing values in ", group_col, ".")))
  }

  ident_1_values <- split_task_values(task_row$ident_1)
  if (length(ident_1_values) == 0) {
    return(make_trajectory_skip(task_row, "skipped", "Missing ident_1."))
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
    return(make_trajectory_skip(
      task_row,
      "skipped",
      paste0("Fewer than ", min_cells_per_group, " cells in one or both trajectory comparison groups."),
      n_ident_1 = n_ident_1,
      n_ident_2 = n_ident_2,
      n_cells = n_cells
    ))
  }
  if (n_cells < min_total_cells) {
    return(make_trajectory_skip(
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

  original_group_by_cell <- stats::setNames(
    as.character(meta[selected_cells, group_col, drop = TRUE]),
    selected_cells
  )

  universe_key <- digest::digest(sort(selected_cells), algo = "xxhash64")
  list(
    status = "ready",
    message = NA_character_,
    task = task_row,
    cells = selected_cells,
    group_by_cell = group_by_cell,
    original_group_by_cell = original_group_by_cell,
    ident_1_label = ident_1_label,
    ident_2_label = ident_2_label,
    n_ident_1 = as.integer(n_ident_1),
    n_ident_2 = as.integer(n_ident_2),
    n_cells = as.integer(n_cells),
    universe_key = universe_key,
    comparison_stub = comparison_task_stub(task_row)
  )
}

attach_comparison_metadata <- function(df, task_info) {
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
  df$TrajectoryUniverseKey <- rep(task_info$universe_key, nrow(df))
  df
}

summarize_comparison_pseudotime <- function(cell_df) {
  out <- cell_df %>%
    dplyr::filter(!is.na(TrajectoryComparisonGroup), is.finite(pseudotime)) %>%
    dplyr::group_by(TrajectoryComparisonGroup) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      pseudotime_mean = mean(pseudotime, na.rm = TRUE),
      pseudotime_median = median_finite(pseudotime),
      pseudotime_q25 = stats::quantile(pseudotime, probs = 0.25, na.rm = TRUE, names = FALSE),
      pseudotime_q75 = stats::quantile(pseudotime, probs = 0.75, na.rm = TRUE, names = FALSE),
      .groups = "drop"
    )

  groups <- unique(as.character(cell_df$TrajectoryComparisonGroup))
  groups <- groups[!is.na(groups) & nzchar(groups)]
  if (length(groups) == 2) {
    x <- cell_df$pseudotime[as.character(cell_df$TrajectoryComparisonGroup) == groups[1]]
    y <- cell_df$pseudotime[as.character(cell_df$TrajectoryComparisonGroup) == groups[2]]
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    wilcox_p <- if (length(x) > 0 && length(y) > 0) {
      tryCatch(stats::wilcox.test(x, y)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    ks_p <- if (length(x) > 0 && length(y) > 0) {
      tryCatch(stats::ks.test(x, y)$p.value, error = function(e) NA_real_)
    } else {
      NA_real_
    }
    out$wilcox_p_two_group <- wilcox_p
    out$ks_p_two_group <- ks_p
  } else {
    out$wilcox_p_two_group <- NA_real_
    out$ks_p_two_group <- NA_real_
  }
  out
}

make_trajectory_task_summary_row <- function(task_info, status, message = NA_character_, universe_reused = NA) {
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
    status = status,
    universe_key = as.character(task_info$universe_key),
    universe_reused = as.logical(universe_reused),
    output_dir = NA_character_,
    message = message,
    stringsAsFactors = FALSE
  )
}

run_deg_comparison_trajectories <- function(
  obj,
  comparison_summary_file,
  output_dir,
  universe_dir,
  summary_dir,
  assay,
  cluster_col,
  cluster_match_col,
  nfeatures,
  npcs,
  root_cluster = NULL,
  min_cells_per_group = 10L,
  min_total_cells = 30L
) {
  if (!file.exists(comparison_summary_file)) {
    warning("Skipping DEG-comparison trajectories because comparison summary is missing: ", comparison_summary_file)
    return(invisible(NULL))
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
    dplyr::filter(status %in% c("success", "cached")) %>%
    dplyr::mutate(trajectory_task_id = dplyr::row_number()) %>%
    as.data.frame(stringsAsFactors = FALSE)
  if (nrow(comparison_tasks) == 0) {
    warning("No successful or cached DEG comparisons found for trajectory analysis.")
    return(invisible(NULL))
  }

  task_infos <- lapply(seq_len(nrow(comparison_tasks)), function(i) {
    build_trajectory_comparison_task(
      obj = obj,
      task_row = comparison_tasks[i, , drop = FALSE],
      cluster_match_col = cluster_match_col,
      min_cells_per_group = min_cells_per_group,
      min_total_cells = min_total_cells
    )
  })

  universe_cache <- new.env(parent = emptyenv())
  task_summary_list <- list()
  group_summary_list <- list()
  all_cells_list <- list()
  all_clusters_list <- list()
  all_nodes_list <- list()
  all_edges_list <- list()

  for (task_info in task_infos) {
    task <- task_info$task
    comparison_label <- as.character(task$comparison)
    comparison_dir <- .ensure_dir(file.path(output_dir, sanitize_path_component(as.character(task$scope)), task_info$comparison_stub))

    if (!identical(task_info$status, "ready")) {
      row <- make_trajectory_task_summary_row(task_info, task_info$status, task_info$message, universe_reused = NA)
      row$output_dir <- comparison_dir
      task_summary_list[[length(task_summary_list) + 1L]] <- row
      next
    }

    universe_key <- task_info$universe_key
    universe_reused <- exists(universe_key, envir = universe_cache, inherits = FALSE)
    if (universe_reused) {
      trajectory_result <- get(universe_key, envir = universe_cache)
    } else {
      universe_output_dir <- .ensure_dir(file.path(universe_dir, paste0("universe_", universe_key)))
      obj@meta.data[task_info$cells, "TrajectoryUniverse"] <- paste0("universe_", universe_key)

      trajectory_result <- tryCatch(
        run_dose_monocle3(
          obj = obj,
          dose_label = paste0("universe_", universe_key),
          dose_cells = task_info$cells,
          output_dir = universe_output_dir,
          assay = assay,
          cluster_col = cluster_col,
          nfeatures = nfeatures,
          npcs = npcs,
          root_cluster = root_cluster,
          group_col = "TrajectoryUniverse",
          analysis_name = "by_03a_DEG_comparison_universe"
        ),
        error = function(e) e
      )

      if (inherits(trajectory_result, "error")) {
        row <- make_trajectory_task_summary_row(task_info, "error", conditionMessage(trajectory_result), universe_reused = FALSE)
        row$output_dir <- comparison_dir
        task_summary_list[[length(task_summary_list) + 1L]] <- row
        next
      }
      assign(universe_key, trajectory_result, envir = universe_cache)
    }

    cell_df <- trajectory_result$cells
    cell_df$TrajectoryComparisonGroup <- unname(task_info$group_by_cell[cell_df$cell])
    cell_df$TrajectoryOriginalGroup <- unname(task_info$original_group_by_cell[cell_df$cell])
    cell_df$TrajectoryComparisonGroup <- factor(
      as.character(cell_df$TrajectoryComparisonGroup),
      levels = c(task_info$ident_1_label, task_info$ident_2_label)
    )
    cell_df <- attach_comparison_metadata(cell_df, task_info)

    cluster_df <- attach_comparison_metadata(trajectory_result$clusters, task_info)
    graph_node_df <- attach_comparison_metadata(trajectory_result$graph_nodes, task_info)
    graph_edge_df <- attach_comparison_metadata(trajectory_result$graph_edges, task_info)

    write_table_csv(cell_df, file.path(comparison_dir, "cells_pseudotime.csv"))
    write_table_csv(cluster_df, file.path(comparison_dir, "cluster_pseudotime_summary.csv"))
    write_table_csv(graph_node_df, file.path(comparison_dir, "principal_graph_nodes.csv"))
    write_table_csv(graph_edge_df, file.path(comparison_dir, "principal_graph_edges.csv"))

    group_summary <- summarize_comparison_pseudotime(cell_df)
    group_summary <- attach_comparison_metadata(group_summary, task_info)
    write_table_csv(group_summary, file.path(comparison_dir, "comparison_group_pseudotime_summary.csv"))

    plot_metadata_trajectory(
      cell_df = cell_df,
      cluster_df = cluster_df,
      edge_df = graph_edge_df,
      cluster_col = cluster_col,
      color_col = "TrajectoryComparisonGroup",
      title_label = comparison_label,
      out_dir = comparison_dir,
      file_stub = "umap_monocle3_by_comparison_group"
    )

    original_values <- unique(as.character(cell_df$TrajectoryOriginalGroup))
    original_values <- original_values[!is.na(original_values) & nzchar(original_values)]
    comparison_values <- unique(as.character(cell_df$TrajectoryComparisonGroup))
    comparison_values <- comparison_values[!is.na(comparison_values) & nzchar(comparison_values)]
    if (length(original_values) > length(comparison_values)) {
      plot_metadata_trajectory(
        cell_df = cell_df,
        cluster_df = cluster_df,
        edge_df = graph_edge_df,
        cluster_col = cluster_col,
        color_col = "TrajectoryOriginalGroup",
        title_label = paste0(comparison_label, " original groups"),
        out_dir = comparison_dir,
        file_stub = "umap_monocle3_by_original_group"
      )
    }

    row <- make_trajectory_task_summary_row(task_info, "success", NA_character_, universe_reused = universe_reused)
    row$output_dir <- comparison_dir
    task_summary_list[[length(task_summary_list) + 1L]] <- row
    group_summary_list[[length(group_summary_list) + 1L]] <- group_summary
    all_cells_list[[length(all_cells_list) + 1L]] <- cell_df
    all_clusters_list[[length(all_clusters_list) + 1L]] <- cluster_df
    all_nodes_list[[length(all_nodes_list) + 1L]] <- graph_node_df
    all_edges_list[[length(all_edges_list) + 1L]] <- graph_edge_df
  }

  task_summary <- dplyr::bind_rows(task_summary_list)
  group_summary_all <- dplyr::bind_rows(group_summary_list)
  write_table_csv(task_summary, file.path(summary_dir, "trajectory_DEG_comparison_task_summary.csv"))
  write_table_csv(group_summary_all, file.path(summary_dir, "trajectory_DEG_comparison_group_pseudotime_summary.csv"))
  write_table_csv(dplyr::bind_rows(all_cells_list), file.path(summary_dir, "pseudotime_cells_all_DEG_comparisons.csv"))
  write_table_csv(dplyr::bind_rows(all_clusters_list), file.path(summary_dir, "cluster_pseudotime_summary_all_DEG_comparisons.csv"))
  write_table_csv(dplyr::bind_rows(all_nodes_list), file.path(summary_dir, "principal_graph_nodes_all_DEG_comparisons.csv"))
  write_table_csv(dplyr::bind_rows(all_edges_list), file.path(summary_dir, "principal_graph_edges_all_DEG_comparisons.csv"))

  invisible(list(
    task_summary = task_summary,
    group_summary = group_summary_all
  ))
}

input_rds <- cfg_value(
  config,
  "PseudoTrajectory_input_rds",
  "PSEUDOTRAJ_SEURAT_RDS",
  default = resolve_first_existing(c(
    file.path(results_root, "03c_cluster_annotation", "04_objects", "integrated_sct_cca_seurat_cluster_final_annotation.rds"),
    file.path(results_root, "03b_manual_cluster_merge", "objects", "integrated_sct_cca_seurat_final_manual_merge.rds"),
    file.path(results_root, "03_final_cluster", "03_objects", "integrated_sct_cca_seurat_final_reclustered.rds")
  ))
)
output_root <- cfg_value(
  config,
  "PseudoTrajectory_output_root",
  "PSEUDOTRAJ_OUTPUT_ROOT",
  default = file.path(results_root, "04a_psudo_rajectory")
)
assay_use <- cfg_value(config, "PseudoTrajectory_assay", "PSEUDOTRAJ_ASSAY", default = "RNA")
nfeatures_use <- parse_int_cfg(cfg_value(config, "PseudoTrajectory_nfeatures", "PSEUDOTRAJ_NFEATURES", default = "2000"), 2000L, "PseudoTrajectory_nfeatures")
npcs_use <- parse_int_cfg(cfg_value(config, "PseudoTrajectory_npcs", "PSEUDOTRAJ_NPCS", default = "30"), 30L, "PseudoTrajectory_npcs")
root_cluster <- cfg_value(config, "PseudoTrajectory_root_cluster", "PSEUDOTRAJ_ROOT_CLUSTER", default = "")
root_cluster <- if (is.null(root_cluster) || !nzchar(root_cluster)) NULL else root_cluster
trajectory_branch <- "root14_end13"
missing_dose_action <- tolower(cfg_value(
  config,
  "PseudoTrajectory_missing_dose_action",
  "PSEUDOTRAJ_MISSING_DOSE_ACTION",
  default = "exclude"
))
if (!(missing_dose_action %in% c("exclude", "stop"))) {
  stop("PseudoTrajectory_missing_dose_action must be 'exclude' or 'stop'.", call. = FALSE)
}

dose_col <- "Dose"
ids_col <- "IDs"
expected_doses <- c("0mg/kg", "30mg/kg", "120mg/kg")
cluster_col <- "cluster_final"

message("Preparing output directories.")
.ensure_dir(output_root)
out_pseudotrajectory_groups <- .ensure_dir(file.path(output_root, "01_pseudotrajectory_groups"))
out_summary <- .ensure_dir(file.path(output_root, "02_summary"))
out_plots <- .ensure_dir(file.path(output_root, "03_plots"))

if (!file.exists(input_rds)) {
  stop("Input Seurat object does not exist: ", input_rds, call. = FALSE)
}

message("Reading Seurat object: ", input_rds)
obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) {
  stop("Input file is not a Seurat object.", call. = FALSE)
}
assert_required_column(obj@meta.data, dose_col)
assert_required_column(obj@meta.data, ids_col)
if (!(assay_use %in% names(obj@assays))) {
  stop("Assay is missing in Seurat object: ", assay_use, call. = FALSE)
}
if (!("umap" %in% names(obj@reductions))) {
  stop("UMAP reduction is missing in Seurat object.", call. = FALSE)
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
  stop("Cannot find cluster_final or seurat_clusters metadata column.", call. = FALSE)
}

dose_raw <- obj@meta.data[[dose_col]]
dose_standardized <- standardize_in_vivo_dose(dose_raw)
missing_dose <- is.na(dose_standardized)
if (any(missing_dose)) {
  missing_dose_file <- file.path(
    out_summary,
    paste0("cells_with_missing_Dose_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  )
  utils::write.csv(
    data.frame(
      cell = rownames(obj@meta.data)[missing_dose],
      Dose_raw = as.character(dose_raw)[missing_dose],
      stringsAsFactors = FALSE
    ),
    missing_dose_file,
    row.names = FALSE
  )
  missing_msg <- paste0(
    "Dose contains NA or blank values in ", sum(missing_dose),
    " cell(s). Details were written to: ", missing_dose_file
  )
  if (missing_dose_action == "stop") {
    stop(missing_msg, call. = FALSE)
  }
  message(missing_msg)
  message("Keeping missing-Dose cells in the Seurat object; excluding them only from analyses that require Dose.")
}
obj@meta.data[[dose_col]] <- dose_standardized
dose_values <- as.character(obj@meta.data[[dose_col]])
dose_values_present <- dose_values[!is.na(dose_values)]
missing_expected_doses <- setdiff(expected_doses, unique(dose_values_present))
if (length(missing_expected_doses) > 0) {
  stop("Missing expected Dose group(s): ", paste(missing_expected_doses, collapse = ", "), call. = FALSE)
}
dose_levels <- c(expected_doses, sort(setdiff(unique(dose_values_present), expected_doses)))
dose_levels <- dose_levels[dose_levels %in% unique(dose_values_present)]
obj@meta.data[[dose_col]] <- factor(dose_values, levels = dose_levels)

ids_values <- as.character(obj@meta.data[[ids_col]])
ids_values[is.na(ids_values)] <- ""
has_2n <- grepl("2N", ids_values, fixed = TRUE)
has_4n <- grepl("4N", ids_values, fixed = TRUE)
karyotype <- rep(NA_character_, length(ids_values))
karyotype[has_2n & !has_4n] <- "2N"
karyotype[has_4n & !has_2n] <- "4N"
obj@meta.data$Karyotype <- factor(karyotype, levels = c("2N", "4N"))
obj@meta.data$SampleOrigin <- factor(
  ifelse(grepl("Cell-Culture", ids_values, fixed = TRUE), "CellLine", "Tumor"),
  levels = c("CellLine", "Tumor")
)
obj@meta.data$Ploidy <- factor(as.character(obj@meta.data$Karyotype), levels = c("2N", "4N"))
obj@meta.data$TN <- factor(as.character(obj@meta.data$SampleOrigin), levels = c("CellLine", "Tumor"))
obj@meta.data$Dose_DEG <- factor(as.character(obj@meta.data[[dose_col]]), levels = expected_doses)
obj@meta.data$trajectory_context <- factor(as.character(obj@meta.data$SampleOrigin), levels = c("CellLine", "Tumor"))
obj@meta.data$trajectory_group <- factor(as.character(obj@meta.data$TN), levels = c("CellLine", "Tumor"))
if (!("clusters" %in% colnames(obj@meta.data))) {
  obj@meta.data$clusters <- as.character(obj@meta.data[[cluster_col]])
}

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
      idx <- as.character(obj@meta.data$TN) == .data$tn_scope
      if (.data$ploidy_scope %in% c("2N", "4N")) {
        idx <- idx & as.character(obj@meta.data$Ploidy) == .data$ploidy_scope
      }
      sum(idx, na.rm = TRUE)
    }
  ) %>%
  dplyr::ungroup()
pseudotrajectory_group_df <- dplyr::bind_rows(
  data.frame(
    trajectory_group = "All_cells",
    trajectory_context = "All",
    Dose = NA_character_,
    ploidy = NA_character_,
    group_safe = "All_cells",
    analysis_type = "all_cells",
    tn_scope = "All",
    ploidy_scope = "ploidy_all",
    n_cells = nrow(obj@meta.data),
    stringsAsFactors = FALSE
  ),
  tn_ploidy_tasks
)
pseudotrajectory_group_df <- pseudotrajectory_group_df %>%
  dplyr::filter(.data$n_cells > 0)
if (nrow(pseudotrajectory_group_df) == 0) {
  stop("No cells were found for the requested pseudotrajectory analyses.", call. = FALSE)
}
pseudotrajectory_group_levels <- unique(pseudotrajectory_group_df$trajectory_group)

unresolved_karyotype <- is.na(karyotype) | (has_2n & has_4n)
if (any(unresolved_karyotype)) {
  unresolved_karyotype_file <- file.path(
    out_summary,
    paste0("cells_with_unresolved_Karyotype_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
  )
  utils::write.csv(
    data.frame(
      cell = rownames(obj@meta.data)[unresolved_karyotype],
      IDs = ids_values[unresolved_karyotype],
      has_2N = has_2n[unresolved_karyotype],
      has_4N = has_4n[unresolved_karyotype],
      stringsAsFactors = FALSE
    ),
    unresolved_karyotype_file,
    row.names = FALSE
  )
  message(
    "Karyotype could not be uniquely resolved from IDs in ", sum(unresolved_karyotype),
    " cell(s). Details were written to: ", unresolved_karyotype_file
  )
}

dose_counts <- write_group_counts(obj@meta.data[[dose_col]], "Dose", file.path(out_summary, "dose_cell_counts.csv"), order_levels = dose_levels)
write_group_counts(obj@meta.data$Karyotype, "Karyotype", file.path(out_summary, "karyotype_cell_counts.csv"), order_levels = c("2N", "4N"))
write_group_counts(obj@meta.data$Ploidy, "Ploidy", file.path(out_summary, "ploidy_cell_counts.csv"), order_levels = c("2N", "4N"))
write_group_counts(obj@meta.data$SampleOrigin, "SampleOrigin", file.path(out_summary, "sample_origin_cell_counts.csv"), order_levels = c("CellLine", "Tumor"))
write_group_counts(obj@meta.data$TN, "TN", file.path(out_summary, "tn_cell_counts.csv"), order_levels = c("CellLine", "Tumor"))
write_group_counts(obj@meta.data$Dose_DEG, "Dose_DEG", file.path(out_summary, "dose_deg_cell_counts.csv"), order_levels = expected_doses)
write_group_counts(obj@meta.data$trajectory_context, "trajectory_context", file.path(out_summary, "trajectory_context_cell_counts.csv"), order_levels = c("CellLine", "Tumor"))
write_group_counts(obj@meta.data$trajectory_group, "trajectory_group", file.path(out_summary, "trajectory_group_cell_counts_by_cell_metadata.csv"), order_levels = c("CellLine", "Tumor"))
write_table_csv(pseudotrajectory_group_df, file.path(out_summary, "pseudotrajectory_group_levels.csv"))
write_table_csv(
  pseudotrajectory_group_df %>%
    dplyr::select(trajectory_group, tn_scope, ploidy_scope, trajectory_context, Dose, ploidy, n_cells, analysis_type, group_safe),
  file.path(out_summary, "pseudotrajectory_group_cell_counts.csv")
)

get_pseudotrajectory_group_cells <- function(task_row) {
  group_name <- task_row$trajectory_group[1]
  if (identical(group_name, "All_cells")) {
    return(rownames(obj@meta.data))
  }
  idx <- as.character(obj@meta.data$TN) == task_row$tn_scope[1]
  if (task_row$ploidy_scope[1] %in% c("2N", "4N")) {
    idx <- idx & as.character(obj@meta.data$Ploidy) == task_row$ploidy_scope[1]
  }
  rownames(obj@meta.data)[idx]
}

message("Running Monocle3 pseudotrajectory analyses: All_cells plus CellLine/Tumor ploidy_all, 2N, and 4N.")
pseudotrajectory_results <- list()
for (i in seq_len(nrow(pseudotrajectory_group_df))) {
  group_name <- pseudotrajectory_group_df$trajectory_group[i]
  group_cells <- get_pseudotrajectory_group_cells(pseudotrajectory_group_df[i, , drop = FALSE])
  group_dir <- .ensure_dir(file.path(out_pseudotrajectory_groups, pseudotrajectory_group_df$group_safe[i]))
  message("Running Monocle3 pseudotime for pseudotrajectory group ", group_name)
  pseudotrajectory_results[[group_name]] <- run_dose_monocle3(
    obj = obj,
    dose_label = group_name,
    dose_cells = group_cells,
    output_dir = group_dir,
    assay = assay_use,
    cluster_col = cluster_col,
    nfeatures = nfeatures_use,
    npcs = npcs_use,
    root_cluster = root_cluster,
    group_col = "trajectory_analysis_group",
    analysis_name = pseudotrajectory_group_df$analysis_type[i],
    tn_scope = pseudotrajectory_group_df$tn_scope[i],
    ploidy_scope = pseudotrajectory_group_df$ploidy_scope[i],
    branch_label = trajectory_branch
  )
}

pseudotrajectory_outputs <- write_result_collection_outputs(
  results = pseudotrajectory_results,
  group_col = "trajectory_analysis_group",
  analysis_name = "pseudotrajectory_groups",
  group_title = "Pseudo trajectory analysis group",
  table_suffix = "pseudotrajectory_groups",
  plot_suffix = "pseudotrajectory_groups",
  out_summary = out_summary,
  out_plots = out_plots
)

summary_lines <- c(
  "04a Monocle3 expression-only pseudotime workflow completed.",
  paste0("Input Seurat object: ", input_rds),
  paste0("Output root: ", output_root),
  paste0("Assay used for Monocle3 matrix: ", assay_use),
  paste0("Requested pseudotrajectory groups: ", paste(pseudotrajectory_group_levels, collapse = ", ")),
  paste0("Observed Dose groups retained as metadata: ", paste(dose_levels, collapse = ", ")),
  paste0("Cluster column: ", cluster_col),
  paste0("Variable features per group: ", nfeatures_use),
  paste0("PCA dimensions requested: ", npcs_use),
  paste0("Root cluster override: ", ifelse(is.null(root_cluster), "not set; smallest cluster label per group used", root_cluster)),
  paste0("Aligned trajectory branch label: ", trajectory_branch),
  paste0("Missing Dose action: ", missing_dose_action),
  "",
  "Metadata derivation:",
  "  Karyotype/Ploidy is parsed from IDs: contains 2N -> 2N; contains 4N -> 4N.",
  "  SampleOrigin is parsed from IDs: contains Cell-Culture -> CellLine; otherwise Tumor.",
  "  Active pseudotrajectory groups are aligned to 04_trajectory.R: All_cells plus CellLine/Tumor ploidy_all, 2N, and 4N.",
  "",
  "Important interpretation:",
  "  This is Monocle3 pseudotime from the expression matrix, not true RNA velocity.",
  "  True RNA velocity requires spliced/unspliced counts from BAM/loom/h5ad.",
  "  The arrows show Monocle3 principal graph direction oriented away from the default root graph node.",
  "",
  "Main outputs:",
  "  01_pseudotrajectory_groups/All_cells/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/CellLine/ploidy_all/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/CellLine/2N/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/CellLine/4N/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/Tumor/ploidy_all/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/Tumor/2N/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/Tumor/4N/cells_pseudotime.csv",
  "  01_pseudotrajectory_groups/*/cluster_pseudotime_summary.csv",
  "  01_pseudotrajectory_groups/*/principal_graph_nodes.csv",
  "  01_pseudotrajectory_groups/*/principal_graph_edges.csv",
  "  02_summary/pseudotime_cells_all_pseudotrajectory_groups.csv",
  "  02_summary/cluster_pseudotime_summary_all_pseudotrajectory_groups.csv",
  "  02_summary/principal_graph_nodes_all_pseudotrajectory_groups.csv",
  "  02_summary/principal_graph_edges_all_pseudotrajectory_groups.csv",
  "  02_summary/pseudotrajectory_group_levels.csv",
  "  02_summary/pseudotrajectory_group_cell_counts.csv",
  "  03_plots/umap_monocle3_pseudotime_by_pseudotrajectory_groups.pdf(.png)",
  "  03_plots/pseudotime_density_by_pseudotrajectory_groups.pdf(.png)",
  "  03_plots/pseudotime_distribution_by_pseudotrajectory_groups.pdf(.png)"
)
writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

saveRDS(obj, file.path(output_root, "seurat_object_with_pseudotrajectory_metadata.rds"))

message("04a Monocle3 pseudotime workflow finished.")
message("Output root: ", normalizePath(output_root, mustWork = FALSE))
