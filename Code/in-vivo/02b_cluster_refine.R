#!/usr/bin/env Rscript

script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_match <- grep(file_arg, cmd_args, value = TRUE)
if (length(file_match) > 0) {
  script_path <- normalizePath(sub(file_arg, "", file_match[1]), mustWork = FALSE)
}
if (is.null(script_path) || !nzchar(script_path)) {
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = FALSE) else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    script_path <- frame_files[length(frame_files)]
  }
}
script_dir <- if (!is.null(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
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

input_rds <- file.path(results_root, "01a_cell_cycle", "integrated_sct_cca_seurat_cell_cycle.rds")
output_root <- file.path(results_root, "02b_cluster_refine")

umap_reduction <- "umap"
core_cell_cycle_clusters <- c("6", "10", "11", "12")
neighbor_clusters <- c("0", "3")
candidate_clusters <- c("4", "9")
region_clusters <- c(core_cell_cycle_clusters, neighbor_clusters)
candidate_label_map <- c("4" = "4c", "9" = "9c")

knn_k <- 50L
min_region_neighbor_frac <- 0.60
min_core_neighbor_frac <- 0.20
require_inside_reference_hull <- TRUE
component_expand_k <- 5L
component_expand_min_region_neighbor_frac <- 0.00
component_expand_min_core_neighbor_frac <- 0.00
component_expand_require_inside_reference_hull <- FALSE
component_expand_use_mutual_knn <- FALSE

build_refined_levels <- function(original_levels, label_map) {
  out <- character(0)
  for (cl in original_levels) {
    out <- c(out, cl)
    if (cl %in% names(label_map)) {
      out <- c(out, unname(label_map[[cl]]))
    }
  }
  unique(out)
}

safe_fraction <- function(x) {
  if (length(x) == 0) return(NA_real_)
  mean(x, na.rm = TRUE)
}

point_in_polygon <- function(x, y, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    xi <- poly_x[i]
    yi <- poly_y[i]
    xj <- poly_x[j]
    yj <- poly_y[j]
    crosses <- ((yi > y) != (yj > y)) &
      (x < ((xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi))
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

compute_knn_indices <- function(coords, k) {
  if (!requireNamespace("RANN", quietly = TRUE)) {
    stop("Package 'RANN' is required for UMAP kNN refinement but is not installed.", call. = FALSE)
  }
  nn <- RANN::nn2(data = coords, query = coords, k = k + 1L)
  nn$nn.idx[, -1, drop = FALSE]
}

compute_knn_graph_components <- function(coords, k, mutual = FALSE) {
  n <- nrow(coords)
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)

  k_use <- min(as.integer(k), n - 1L)
  if (!is.finite(k_use) || k_use < 1L) {
    return(seq_len(n))
  }

  nn_idx <- compute_knn_indices(coords, k = k_use)
  if (!is.matrix(nn_idx)) {
    nn_idx <- matrix(nn_idx, nrow = n, ncol = k_use)
  }

  adj <- vector("list", n)
  for (i in seq_len(n)) adj[[i]] <- integer(0)

  if (mutual) {
    nn_sets <- lapply(seq_len(n), function(i) unique(as.integer(nn_idx[i, ])))
    for (i in seq_len(n)) {
      for (j in nn_sets[[i]]) {
        if (length(j) == 0L || is.na(j) || j < 1L || j > n) next
        if (i %in% nn_sets[[j]]) {
          adj[[i]] <- c(adj[[i]], j)
          adj[[j]] <- c(adj[[j]], i)
        }
      }
    }
  } else {
    for (i in seq_len(n)) {
      nbrs <- unique(as.integer(nn_idx[i, ]))
      nbrs <- nbrs[!is.na(nbrs) & nbrs >= 1L & nbrs <= n]
      if (length(nbrs) == 0L) next
      adj[[i]] <- c(adj[[i]], nbrs)
      for (j in nbrs) {
        adj[[j]] <- c(adj[[j]], i)
      }
    }
  }

  adj <- lapply(adj, unique)
  component <- integer(n)
  component_id <- 0L

  for (i in seq_len(n)) {
    if (component[i] != 0L) next
    component_id <- component_id + 1L
    queue <- i
    component[i] <- component_id
    head <- 1L

    while (head <= length(queue)) {
      node <- queue[head]
      head <- head + 1L
      nbrs <- adj[[node]]
      if (length(nbrs) == 0L) next
      new_nodes <- nbrs[component[nbrs] == 0L]
      if (length(new_nodes) == 0L) next
      component[new_nodes] <- component_id
      queue <- c(queue, new_nodes)
    }
  }

  component
}

expand_candidate_component <- function(
  candidate_idx,
  seed_idx,
  umap_mat,
  frac_region,
  frac_core,
  inside_reference_hull,
  min_region_frac,
  min_core_frac,
  require_inside_hull = TRUE,
  graph_k = 5L,
  mutual_knn = TRUE
) {
  if (length(candidate_idx) == 0L || length(seed_idx) == 0L) {
    return(list(selected_idx = integer(0), eligible_idx = integer(0)))
  }

  eligible_mask <- frac_region[candidate_idx] >= min_region_frac &
    frac_core[candidate_idx] >= min_core_frac &
    (!require_inside_hull | inside_reference_hull[candidate_idx])
  eligible_idx <- candidate_idx[eligible_mask]
  if (length(eligible_idx) == 0L) {
    return(list(selected_idx = integer(0), eligible_idx = integer(0)))
  }

  seed_in_eligible <- eligible_idx %in% seed_idx
  if (!any(seed_in_eligible)) {
    return(list(selected_idx = integer(0), eligible_idx = eligible_idx))
  }

  component <- compute_knn_graph_components(
    coords = umap_mat[eligible_idx, , drop = FALSE],
    k = graph_k,
    mutual = mutual_knn
  )
  selected_components <- unique(component[seed_in_eligible])
  selected_idx <- eligible_idx[component %in% selected_components]

  list(selected_idx = selected_idx, eligible_idx = eligible_idx)
}

build_focus_status <- function(cluster_vec, refined_vec, region_clusters) {
  out <- rep("other", length(cluster_vec))
  out[cluster_vec %in% region_clusters] <- "reference_region"
  out[cluster_vec == "4"] <- "cluster4_other"
  out[cluster_vec == "9"] <- "cluster9_other"
  out[refined_vec == "4c"] <- "cluster4c_selected"
  out[refined_vec == "9c"] <- "cluster9c_selected"
  factor(
    out,
    levels = c("reference_region", "cluster4_other", "cluster4c_selected", "cluster9_other", "cluster9c_selected", "other")
  )
}

.ensure_dir(output_root)
out_plots <- .ensure_dir(file.path(output_root, "plots"))
out_summary <- .ensure_dir(file.path(output_root, "summaries"))
out_objects <- .ensure_dir(file.path(output_root, "objects"))
out_cells <- .ensure_dir(file.path(output_root, "selected_cells"))

  if (!file.exists(input_rds)) {
    stop("Input RDS does not exist: ", input_rds, call. = FALSE)
  }

  message("Reading Seurat object: ", input_rds)
  obj <- readRDS(input_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input file is not a Seurat object.", call. = FALSE)
  }

  assert_required_column(obj@meta.data, "seurat_clusters")
  if (!(umap_reduction %in% names(obj@reductions))) {
    stop("UMAP reduction '", umap_reduction, "' not found in Seurat object.", call. = FALSE)
  }

  cluster_vec <- as.character(obj@meta.data[["seurat_clusters"]])
  cluster_levels <- sort_maybe_numeric(cluster_vec)
  obj@meta.data[["seurat_clusters"]] <- factor(cluster_vec, levels = cluster_levels)

  missing_config_clusters <- setdiff(unique(c(region_clusters, candidate_clusters)), unique(cluster_vec))
  if (length(missing_config_clusters) > 0) {
    stop(
      "Configured clusters were not found in seurat_clusters: ",
      paste(missing_config_clusters, collapse = ", "),
      call. = FALSE
    )
  }

  umap_mat <- Seurat::Embeddings(obj, reduction = umap_reduction)
  if (ncol(umap_mat) < 2) {
    stop("UMAP embedding must contain at least 2 dimensions.", call. = FALSE)
  }
  umap_mat <- umap_mat[, 1:2, drop = FALSE]
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

  message("Computing UMAP neighborhood statistics.")
  nn_idx <- compute_knn_indices(umap_mat, k = knn_k)
  neighbor_labels <- matrix(cluster_vec[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  neighbor_is_region <- matrix(
    neighbor_labels %in% region_clusters,
    nrow = nrow(neighbor_labels),
    ncol = ncol(neighbor_labels)
  )
  neighbor_is_core <- matrix(
    neighbor_labels %in% core_cell_cycle_clusters,
    nrow = nrow(neighbor_labels),
    ncol = ncol(neighbor_labels)
  )
  frac_region <- rowMeans(neighbor_is_region, na.rm = TRUE)
  frac_core <- rowMeans(neighbor_is_core, na.rm = TRUE)
  n_region_neighbors <- rowSums(neighbor_is_region, na.rm = TRUE)
  n_core_neighbors <- rowSums(neighbor_is_core, na.rm = TRUE)

  reference_idx <- which(cluster_vec %in% region_clusters)
  reference_umap <- umap_mat[reference_idx, , drop = FALSE]
  hull_order <- grDevices::chull(reference_umap[, 1], reference_umap[, 2])
  hull_xy <- reference_umap[hull_order, , drop = FALSE]
  hull_xy <- rbind(hull_xy, hull_xy[1, , drop = FALSE])

  inside_reference_hull <- point_in_polygon(
    x = umap_mat[, "UMAP_1"],
    y = umap_mat[, "UMAP_2"],
    poly_x = hull_xy[, "UMAP_1"],
    poly_y = hull_xy[, "UMAP_2"]
  )

  candidate_idx <- which(cluster_vec %in% candidate_clusters)
  seed_idx <- candidate_idx[
    frac_region[candidate_idx] >= min_region_neighbor_frac &
      frac_core[candidate_idx] >= min_core_neighbor_frac &
      (!require_inside_reference_hull | inside_reference_hull[candidate_idx])
  ]

  expanded_idx <- integer(0)
  expansion_eligible_idx <- integer(0)
  for (candidate_cluster in candidate_clusters) {
    cluster_candidate_idx <- which(cluster_vec == candidate_cluster)
    cluster_seed_idx <- seed_idx[cluster_vec[seed_idx] == candidate_cluster]
    expanded_res <- expand_candidate_component(
      candidate_idx = cluster_candidate_idx,
      seed_idx = cluster_seed_idx,
      umap_mat = umap_mat,
      frac_region = frac_region,
      frac_core = frac_core,
      inside_reference_hull = inside_reference_hull,
      min_region_frac = component_expand_min_region_neighbor_frac,
      min_core_frac = component_expand_min_core_neighbor_frac,
      require_inside_hull = component_expand_require_inside_reference_hull,
      graph_k = component_expand_k,
      mutual_knn = component_expand_use_mutual_knn
    )
    expanded_idx <- union(expanded_idx, expanded_res$selected_idx)
    expansion_eligible_idx <- union(expansion_eligible_idx, expanded_res$eligible_idx)
  }

  selected_idx <- sort(unique(c(seed_idx, expanded_idx)))
  selected_seed_cells <- rownames(obj@meta.data)[seed_idx]
  expanded_only_cells <- setdiff(rownames(obj@meta.data)[selected_idx], selected_seed_cells)
  selected_cells <- rownames(obj@meta.data)[selected_idx]

  refined_vec <- cluster_vec
  selected_clusters <- cluster_vec[selected_idx]
  refined_vec[selected_idx] <- unname(candidate_label_map[selected_clusters])
  refined_levels <- build_refined_levels(cluster_levels, candidate_label_map)
  obj@meta.data[["seurat_cluster_refine"]] <- factor(refined_vec, levels = refined_levels)

  focus_status <- build_focus_status(cluster_vec, refined_vec, region_clusters)
  obj@meta.data[["cluster_refine_focus_status"]] <- focus_status

  cell_df <- data.frame(
    cell = rownames(obj@meta.data),
    original_cluster = cluster_vec,
    seurat_cluster_refine = refined_vec,
    UMAP_1 = umap_mat[, "UMAP_1"],
    UMAP_2 = umap_mat[, "UMAP_2"],
    inside_reference_hull = inside_reference_hull,
    frac_region_neighbors = frac_region,
    frac_core_neighbors = frac_core,
    n_region_neighbors = n_region_neighbors,
    n_core_neighbors = n_core_neighbors,
    selected_as_seed = rownames(obj@meta.data) %in% selected_seed_cells,
    selected_by_component_expansion = rownames(obj@meta.data) %in% expanded_only_cells,
    component_expand_eligible = rownames(obj@meta.data) %in% rownames(obj@meta.data)[expansion_eligible_idx],
    selected_for_refine = rownames(obj@meta.data) %in% selected_cells,
    stringsAsFactors = FALSE
  )

  if ("Phase" %in% colnames(obj@meta.data)) {
    cell_df$Phase <- as.character(obj@meta.data[["Phase"]])
  }
  if ("S.Score" %in% colnames(obj@meta.data)) {
    cell_df$S.Score <- as.numeric(obj@meta.data[["S.Score"]])
  }
  if ("G2M.Score" %in% colnames(obj@meta.data)) {
    cell_df$G2M.Score <- as.numeric(obj@meta.data[["G2M.Score"]])
  }

  config_df <- bind_rows(
    data.frame(parameter = "input_rds", value = input_rds, stringsAsFactors = FALSE),
    data.frame(parameter = "output_root", value = output_root, stringsAsFactors = FALSE),
    data.frame(parameter = "umap_reduction", value = umap_reduction, stringsAsFactors = FALSE),
    data.frame(parameter = "core_cell_cycle_clusters", value = paste(core_cell_cycle_clusters, collapse = ","), stringsAsFactors = FALSE),
    data.frame(parameter = "neighbor_clusters", value = paste(neighbor_clusters, collapse = ","), stringsAsFactors = FALSE),
    data.frame(parameter = "candidate_clusters", value = paste(candidate_clusters, collapse = ","), stringsAsFactors = FALSE),
    data.frame(parameter = "knn_k", value = as.character(knn_k), stringsAsFactors = FALSE),
    data.frame(parameter = "min_region_neighbor_frac", value = as.character(min_region_neighbor_frac), stringsAsFactors = FALSE),
    data.frame(parameter = "min_core_neighbor_frac", value = as.character(min_core_neighbor_frac), stringsAsFactors = FALSE),
    data.frame(parameter = "require_inside_reference_hull", value = as.character(require_inside_reference_hull), stringsAsFactors = FALSE),
    data.frame(parameter = "component_expand_k", value = as.character(component_expand_k), stringsAsFactors = FALSE),
    data.frame(parameter = "component_expand_min_region_neighbor_frac", value = as.character(component_expand_min_region_neighbor_frac), stringsAsFactors = FALSE),
    data.frame(parameter = "component_expand_min_core_neighbor_frac", value = as.character(component_expand_min_core_neighbor_frac), stringsAsFactors = FALSE),
    data.frame(parameter = "component_expand_require_inside_reference_hull", value = as.character(component_expand_require_inside_reference_hull), stringsAsFactors = FALSE),
    data.frame(parameter = "component_expand_use_mutual_knn", value = as.character(component_expand_use_mutual_knn), stringsAsFactors = FALSE)
  )
  readr::write_csv(config_df, file.path(out_summary, "refine_parameters.csv"))

  cluster_counts_before <- data.frame(
    cluster = names(table(cluster_vec)),
    n_cells = as.integer(table(cluster_vec)),
    stringsAsFactors = FALSE
  ) %>%
    arrange(match(cluster, cluster_levels))
  readr::write_csv(cluster_counts_before, file.path(out_summary, "cluster_counts_before_refine.csv"))

  cluster_counts_after <- data.frame(
    seurat_cluster_refine = names(table(obj@meta.data[["seurat_cluster_refine"]])),
    n_cells = as.integer(table(obj@meta.data[["seurat_cluster_refine"]])),
    stringsAsFactors = FALSE
  )
  readr::write_csv(cluster_counts_after, file.path(out_summary, "cluster_counts_after_refine.csv"))

  candidate_score_df <- cell_df %>%
    filter(original_cluster %in% candidate_clusters) %>%
    arrange(
      original_cluster,
      desc(selected_for_refine),
      desc(selected_as_seed),
      desc(frac_region_neighbors),
      desc(frac_core_neighbors)
    )
  readr::write_csv(candidate_score_df, file.path(out_summary, "candidate_cluster4_9_region_scores.csv"))

  selected_4c_df <- candidate_score_df %>%
    filter(seurat_cluster_refine == "4c")
  readr::write_csv(selected_4c_df, file.path(out_cells, "cells_selected_as_4c.csv"))

  selected_9c_df <- candidate_score_df %>%
    filter(seurat_cluster_refine == "9c")
  readr::write_csv(selected_9c_df, file.path(out_cells, "cells_selected_as_9c.csv"))

  if ("Phase" %in% colnames(candidate_score_df)) {
    selected_phase_summary <- candidate_score_df %>%
      filter(selected_for_refine) %>%
      count(seurat_cluster_refine, Phase, name = "n_cells") %>%
      arrange(seurat_cluster_refine, Phase)
    readr::write_csv(selected_phase_summary, file.path(out_summary, "selected_refined_cells_phase_counts.csv"))
  }

  saveRDS(obj, file.path(out_objects, "integrated_sct_cca_seurat_cluster_refine.rds"))

  message("Writing diagnostic plots.")
  p_umap_original <- DimPlot(
    obj,
    reduction = umap_reduction,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    pt.size = 0.30
  ) + labs(title = "Original seurat_clusters")
  save_plot_pdf_png(p_umap_original, file.path(out_plots, "umap_original_seurat_clusters"), width = 9, height = 7)

  p_umap_refined <- DimPlot(
    obj,
    reduction = umap_reduction,
    group.by = "seurat_cluster_refine",
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    pt.size = 0.30
  ) + labs(title = "Refined clusters: seurat_cluster_refine")
  save_plot_pdf_png(p_umap_refined, file.path(out_plots, "umap_seurat_cluster_refine"), width = 9, height = 7)

  focus_obj <- subset(obj, subset = seurat_clusters %in% unique(c(region_clusters, candidate_clusters)))

  p_focus_original <- DimPlot(
    focus_obj,
    reduction = umap_reduction,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.35
  ) + labs(title = "Focus UMAP: region clusters 0,3,6,10,11,12 and candidates 4,9")
  save_plot_pdf_png(p_focus_original, file.path(out_plots, "umap_focus_original_clusters"), width = 8.5, height = 6.5)

  p_focus_refined <- DimPlot(
    focus_obj,
    reduction = umap_reduction,
    group.by = "seurat_cluster_refine",
    label = TRUE,
    repel = TRUE,
    pt.size = 0.35
  ) + labs(title = "Focus UMAP: refined 4c and 9c")
  save_plot_pdf_png(p_focus_refined, file.path(out_plots, "umap_focus_seurat_cluster_refine"), width = 8.5, height = 6.5)

  focus_cells <- colnames(focus_obj)
  focus_plot_df <- cell_df %>%
    filter(cell %in% focus_cells) %>%
    mutate(
      focus_status = factor(
        obj@meta.data[cell, "cluster_refine_focus_status", drop = TRUE],
        levels = levels(obj@meta.data[["cluster_refine_focus_status"]])
      )
    )

  p_focus_status <- ggplot(focus_plot_df, aes(x = UMAP_1, y = UMAP_2, color = focus_status)) +
    geom_point(size = 0.45, alpha = 0.9) +
    geom_path(
      data = as.data.frame(hull_xy),
      aes(x = UMAP_1, y = UMAP_2),
      inherit.aes = FALSE,
      linewidth = 0.8,
      linetype = 2,
      color = "black"
    ) +
    scale_color_manual(
      values = c(
        reference_region = "#1f77b4",
        cluster4_other = "#bdbdbd",
        cluster4c_selected = "#d95f02",
        cluster9_other = "#969696",
        cluster9c_selected = "#7570b3",
        other = "#d9d9d9"
      ),
      drop = FALSE
    ) +
    theme_bw(base_size = 11) +
    labs(
      title = "UMAP focus: reference region and selected 4c / 9c cells",
      subtitle = "Dashed line shows convex hull of reference region clusters (0,3,6,10,11,12)",
      color = "Cell status"
    )
  save_plot_pdf_png(p_focus_status, file.path(out_plots, "umap_focus_refine_status"), width = 8.5, height = 6.5)

  summary_lines <- c(
    "Cluster refinement completed.",
    paste0("Input object: ", input_rds),
    paste0("Output root: ", output_root),
    paste0("Reference region clusters: ", paste(region_clusters, collapse = ", ")),
    paste0("Core cell-cycle clusters: ", paste(core_cell_cycle_clusters, collapse = ", ")),
    paste0("Candidate clusters for splitting: ", paste(candidate_clusters, collapse = ", ")),
    paste0("kNN settings on UMAP: k = ", knn_k),
    paste0("Thresholds: frac_region_neighbors >= ", min_region_neighbor_frac, ", frac_core_neighbors >= ", min_core_neighbor_frac),
    paste0("Require inside reference hull: ", require_inside_reference_hull),
    paste0(
      "Component expansion: kNN connected components within full candidate cluster, k = ", component_expand_k,
      ", min frac_region_neighbors >= ", component_expand_min_region_neighbor_frac,
      ", min frac_core_neighbors >= ", component_expand_min_core_neighbor_frac,
      ", require inside hull = ", component_expand_require_inside_reference_hull,
      ", mutual graph = ", component_expand_use_mutual_knn
    ),
    paste0("Selected 4c seed cells: ", sum(cluster_vec[seed_idx] == "4")),
    paste0("Selected 4c cells: ", sum(refined_vec == "4c")),
    paste0("Selected 9c seed cells: ", sum(cluster_vec[seed_idx] == "9")),
    paste0("Selected 9c cells: ", sum(refined_vec == "9c")),
    "",
    "Interpretation:",
    "Seed cells from original cluster 4 or 9 were first identified by embedding in the UMAP neighborhood of the cell-cycle-associated region defined by clusters 0, 3, 6, 10, 11, and 12.",
    "Final 4c and 9c labels were then expanded within each original candidate cluster by taking the kNN-connected component(s) containing those seed cells."
  )
  writeLines(summary_lines, con = file.path(out_summary, "run_summary.txt"))

message("Cluster refinement finished.")
message("Outputs written to: ", normalizePath(output_root))
