extra_all_cells_display_label <- function(analysis_object) {
  out <- as.character(analysis_object)
  out[out == "All_cells"] <- "all_cells"
  out[out %in% c("CellLine/ploidy_all", "CellLine_ploidy_all")] <- "all_cells, CellLine"
  out[out %in% c("Tumor/ploidy_all", "Tumor_ploidy_all")] <- "all_cells, Tumor"
  other <- !(out %in% c("all_cells", "all_cells, CellLine", "all_cells, Tumor"))
  out[other] <- gsub("/", ", ", out[other], fixed = TRUE)
  out
}

extra_all_cells_group_name <- function(analysis_object) {
  out <- gsub("/", "_", as.character(analysis_object), fixed = TRUE)
  out[out == ""] <- "All_cells"
  out
}

extra_all_cells_tn_scope <- function(analysis_object) {
  x <- as.character(analysis_object)
  dplyr::case_when(
    x == "All_cells" ~ "All",
    grepl("^CellLine", x) ~ "CellLine",
    grepl("^Tumor", x) ~ "Tumor",
    TRUE ~ NA_character_
  )
}

extra_all_cells_ploidy_scope <- function(analysis_object) {
  x <- as.character(analysis_object)
  dplyr::case_when(
    grepl("2N$", x) ~ "2N",
    grepl("4N$", x) ~ "4N",
    TRUE ~ "ploidy_all"
  )
}

extra_all_cells_standardize_tn_label <- function(x) {
  out <- as_clean_chr(x)
  lower <- tolower(out)
  out[lower %in% c("cellline", "cell_line", "cell line", "cell-culture", "cell culture", "culture")] <- "CellLine"
  out[lower %in% c("tumor", "tumour", "in-vivo", "invivo", "xenograft")] <- "Tumor"
  out[!(out %in% c("CellLine", "Tumor"))] <- NA_character_
  out
}

extra_all_cells_infer_tn_from_id <- function(x, default_non_cell_culture_to_tumor = TRUE) {
  x <- as_clean_chr(x)
  out <- rep(NA_character_, length(x))
  has_value <- !is.na(x) & nzchar(x)
  out[has_value & grepl("Cell-Culture|CellLine|Cell.Line|Cell_Culture", x, ignore.case = TRUE)] <- "CellLine"
  if (isTRUE(default_non_cell_culture_to_tumor)) {
    out[has_value & is.na(out)] <- "Tumor"
  } else {
    out[has_value & grepl("Tumor|Tumour|Xenograft|In[-_ ]?Vivo", x, ignore.case = TRUE)] <- "Tumor"
  }
  out
}

extra_all_cells_derive_tn_split <- function(df) {
  tn <- rep(NA_character_, nrow(df))
  source <- rep(NA_character_, nrow(df))

  explicit_cols <- c("TN", "SampleOrigin", "trajectory_tn_scope", "sample_type")
  for (candidate in explicit_cols) {
    if (candidate %in% names(df)) {
      values <- extra_all_cells_standardize_tn_label(df[[candidate]])
      fill <- (is.na(tn) | !nzchar(tn)) & !is.na(values) & nzchar(values)
      tn[fill] <- values[fill]
      source[fill] <- candidate
    }
  }

  if ("IDs" %in% names(df)) {
    values <- extra_all_cells_infer_tn_from_id(df$IDs, default_non_cell_culture_to_tumor = TRUE)
    fill <- (is.na(tn) | !nzchar(tn)) & !is.na(values) & nzchar(values)
    tn[fill] <- values[fill]
    source[fill] <- "IDs"
  }

  fallback_cols <- c("sample", "sampleID", "sample_id", "orig.ident", "Sequencing.IDs", "cell")
  for (candidate in fallback_cols) {
    if (candidate %in% names(df)) {
      values <- extra_all_cells_infer_tn_from_id(df[[candidate]], default_non_cell_culture_to_tumor = FALSE)
      fill <- (is.na(tn) | !nzchar(tn)) & !is.na(values) & nzchar(values)
      tn[fill] <- values[fill]
      source[fill] <- candidate
    }
  }

  data.frame(TN = tn, tn_split_source = source, stringsAsFactors = FALSE)
}

extra_all_cells_assign_tn <- function(df) {
  if (nrow(df) == 0) return(df)
  split <- extra_all_cells_derive_tn_split(df)
  df$TN <- split$TN
  df$all_cells_tn_split_source <- split$tn_split_source
  df
}

extra_all_cells_update_analysis_fields <- function(df, analysis_object) {
  group_name <- extra_all_cells_group_name(analysis_object)
  display_label <- extra_all_cells_display_label(analysis_object)
  df$source_label <- analysis_object
  df$trajectory_analysis_group <- group_name
  df$trajectory_group <- group_name
  df$AnalysisGroup <- group_name
  df$trajectory_tn_scope <- extra_all_cells_tn_scope(analysis_object)
  df$trajectory_ploidy_scope <- extra_all_cells_ploidy_scope(analysis_object)
  df$analysis_object <- analysis_object
  df$analysis_group <- group_name
  df$analysis_display_label <- display_label
  df$all_cells_input_source_label <- "All_cells"
  df$all_cells_derived <- analysis_object != "All_cells"
  df
}

extra_all_cells_expand_metrics <- function(metrics) {
  if (nrow(metrics) == 0) return(metrics)
  metrics <- extra_all_cells_assign_tn(metrics)
  resolved_tn <- metrics$TN %in% c("CellLine", "Tumor")
  if (!any(resolved_tn)) {
    stop("Cannot derive CellLine/Tumor outputs because no All_cells rows have TN = CellLine or Tumor.", call. = FALSE)
  }
  if (any(!resolved_tn)) {
    warning(sum(!resolved_tn), " All_cells rows have unresolved TN and will remain only in All_cells.", call. = FALSE)
  }

  rows <- list(extra_all_cells_update_analysis_fields(metrics, "All_cells"))
  for (tn_label in c("CellLine", "Tumor")) {
    sub <- metrics %>% dplyr::filter(.data$TN == .env$tn_label)
    if (nrow(sub) == 0) {
      warning("No All_cells rows were assigned to ", tn_label, ".", call. = FALSE)
      next
    }
    rows[[length(rows) + 1L]] <- extra_all_cells_update_analysis_fields(sub, paste0(tn_label, "/ploidy_all"))
  }
  dplyr::bind_rows(rows)
}

extra_all_cells_tn_split_summary <- function(metrics) {
  metrics <- extra_all_cells_assign_tn(metrics)
  metrics %>%
    dplyr::mutate(
      TN = dplyr::if_else(.data$TN %in% c("CellLine", "Tumor"), .data$TN, "Unresolved"),
      all_cells_tn_split_source = dplyr::if_else(
        !is.na(.data$all_cells_tn_split_source) & nzchar(.data$all_cells_tn_split_source),
        .data$all_cells_tn_split_source,
        "unresolved"
      )
    ) %>%
    dplyr::count(.data$method, .data$source_context, .data$source_label, .data$TN, .data$all_cells_tn_split_source, name = "n_cells") %>%
    dplyr::arrange(.data$method, .data$source_context, .data$source_label, .data$TN, .data$all_cells_tn_split_source)
}

write_method_separate_outputs_all_cells_derived <- function(df, metric_col, method_label, out_dir, plot_prefix, out_umap_dir, cluster_annotations = data.frame()) {
  analysis_objects <- c("All_cells", "CellLine/ploidy_all", "Tumor/ploidy_all")
  analysis_objects <- analysis_objects[analysis_objects %in% unique(as.character(df$source_label))]
  if (length(analysis_objects) == 0) {
    stop("No All_cells-derived analysis objects were available for method separate outputs.", call. = FALSE)
  }
  cleanup_method_separate_top_level_all_cells_outputs(out_dir, out_umap_dir, plot_prefix)

  object_results <- list()
  status_rows <- list()
  for (obj in analysis_objects) {
    obj_df <- df %>% dplyr::filter(.data$source_label == .env$obj)
    if (nrow(obj_df) == 0) next
    obj_safe <- safe_file_stub(gsub("/", "_", obj, fixed = TRUE))
    obj_display <- extra_all_cells_display_label(obj)
    obj_out_dir <- .ensure_dir(file.path(out_dir, obj_safe))
    obj_umap_dir <- .ensure_dir(file.path(out_umap_dir, obj_safe))
    obj_label <- paste(method_label, obj_display)
    result <- write_method_separate_outputs(
      obj_df,
      metric_col = metric_col,
      method_label = obj_label,
      out_dir = obj_out_dir,
      plot_prefix = plot_prefix,
      out_umap_dir = obj_umap_dir,
      cluster_annotations = cluster_annotations
    )
    result$analysis_object <- obj
    result$analysis_display_label <- obj_display
    object_results[[obj]] <- result
    status_rows[[obj]] <- data.frame(
      analysis_object = obj,
      analysis_display_label = obj_display,
      output_dir = obj_out_dir,
      umap_output_dir = obj_umap_dir,
      n_input_rows = nrow(obj_df),
      n_descriptive_rows = nrow(result$desc),
      stringsAsFactors = FALSE
    )
  }

  status <- dplyr::bind_rows(status_rows)
  write_table_csv(status, file.path(out_dir, paste0(plot_prefix, "_all_cells_derived_method_separate_status.csv")))
  all_cells_result <- object_results[["All_cells"]]
  if (is.null(all_cells_result)) all_cells_result <- object_results[[1]]
  all_cells_result$all_cells_derived_objects <- object_results
  all_cells_result$all_cells_derived_status <- status
  all_cells_result
}

make_method_answer_summary_all_cells_derived <- function(method_result, method_label) {
  object_results <- method_result$all_cells_derived_objects
  if (is.null(object_results) || length(object_results) == 0) {
    return(make_method_answer_summary(method_result, method_label))
  }
  rows <- lapply(names(object_results), function(obj) {
    obj_display <- extra_all_cells_display_label(obj)
    out <- make_method_answer_summary(object_results[[obj]], paste(method_label, obj_display))
    out$analysis_object <- obj
    out$analysis_display_label <- obj_display
    out
  })
  dplyr::bind_rows(rows)
}

cleanup_method_separate_top_level_all_cells_outputs <- function(out_dir, out_umap_dir, plot_prefix) {
  status_file <- file.path(out_dir, paste0(plot_prefix, "_all_cells_derived_method_separate_status.csv"))
  candidates <- c(
    list.files(out_dir, pattern = paste0("^", plot_prefix, "_"), full.names = TRUE, recursive = FALSE),
    list.files(out_umap_dir, pattern = paste0("^", plot_prefix, "_"), full.names = TRUE, recursive = FALSE)
  )
  candidates <- setdiff(candidates, status_file)
  if (length(candidates) > 0) unlink(candidates, recursive = FALSE, force = TRUE)
  invisible(candidates)
}

plot_cluster_hypoxia_umap_all_cells_derived <- function(method_result, gsea_inputs, method_label, metric_col, out_dir, plot_prefix, out_umap_dir) {
  object_results <- method_result$all_cells_derived_objects
  if (is.null(object_results) || length(object_results) == 0) {
    return(plot_cluster_hypoxia_umap(
      method_result$desc,
      gsea_inputs,
      method_label = method_label,
      metric_col = metric_col,
      out_dir = out_dir,
      plot_prefix = plot_prefix,
      out_umap_dir = out_umap_dir
    ))
  }
  for (obj in names(object_results)) {
    obj_safe <- safe_file_stub(gsub("/", "_", obj, fixed = TRUE))
    obj_display <- extra_all_cells_display_label(obj)
    plot_cluster_hypoxia_umap(
      object_results[[obj]]$desc,
      gsea_inputs,
      method_label = paste(method_label, obj_display),
      metric_col = metric_col,
      out_dir = .ensure_dir(file.path(out_dir, obj_safe)),
      plot_prefix = plot_prefix,
      out_umap_dir = .ensure_dir(file.path(out_umap_dir, obj_safe))
    )
  }
  invisible(NULL)
}

extra_all_cells_read_velocity_metrics <- function(root_dir) {
  file <- file.path(root_dir, "01_velocity_groups", "All_cells", "scvelo_cell_metrics.csv")
  if (!file.exists(file) || file.info(file)$size <= 0) {
    stop("Missing All_cells scVelo metric file: ", file, call. = FALSE)
  }
  keep_cols <- unique(c(metadata_cols, velocity_metric_cols))
  message("Reading All_cells scVelo metrics: ", file)
  df <- read_csv_selected(file, keep_cols)
  df <- merge_same_dir_umap(df, dirname(file))
  for (col in intersect(metadata_cols, names(df))) df[[col]] <- as.character(df[[col]])
  for (col in intersect(velocity_metric_cols, names(df))) df[[col]] <- safe_numeric(df[[col]])
  df$source_file <- normalizePath(file, mustWork = FALSE)
  df$source_dir <- normalizePath(dirname(file), mustWork = FALSE)
  df$source_context <- "velocity_groups"
  df$source_label <- "All_cells"
  df$method <- method_scvelo
  df$metric_family <- "velocity"
  manifest <- data.frame(
    source_file = df$source_file[1],
    source_dir = df$source_dir[1],
    source_context = "velocity_groups",
    source_label = "All_cells",
    read_mode = "All_cells_only",
    stringsAsFactors = FALSE
  )
  write_table_csv(manifest, file.path(out_inputs, "scVelo_all_cells_metric_files.csv"))
  df <- harmonize_common_metadata(df)
  for (metric in velocity_metric_cols) df[[metric]] <- safe_numeric(df[[metric]])
  df
}

extra_all_cells_read_monocle3_metrics <- function(root_dir) {
  file <- file.path(root_dir, "01_pseudotrajectory_groups", "All_cells", "cells_pseudotime.csv")
  if (!file.exists(file) || file.info(file)$size <= 0) {
    stop("Missing All_cells Monocle3 metric file: ", file, call. = FALSE)
  }
  keep_cols <- unique(c(metadata_cols, pseudo_metric_cols, "pseudotime_raw", "root_cluster", "monocle3_matrix_source"))
  message("Reading All_cells Monocle3 metrics: ", file)
  df <- read_csv_selected(file, keep_cols)
  for (col in intersect(metadata_cols, names(df))) df[[col]] <- as.character(df[[col]])
  for (col in intersect(c(pseudo_metric_cols, "pseudotime_raw"), names(df))) df[[col]] <- safe_numeric(df[[col]])
  df$source_file <- normalizePath(file, mustWork = FALSE)
  df$source_dir <- normalizePath(dirname(file), mustWork = FALSE)
  df$source_context <- "monocle3_groups"
  df$source_label <- "All_cells"
  df$method <- method_monocle3
  df$metric_family <- "pseudotime"
  manifest <- data.frame(
    source_file = df$source_file[1],
    source_dir = df$source_dir[1],
    source_context = "monocle3_groups",
    source_label = "All_cells",
    read_mode = "All_cells_only",
    stringsAsFactors = FALSE
  )
  write_table_csv(manifest, file.path(out_inputs, "Monocle3_all_cells_metric_files.csv"))
  df <- harmonize_common_metadata(df)
  df$pseudotime <- safe_numeric(df$pseudotime)
  df$pseudotime_raw <- safe_numeric(df$pseudotime_raw)
  df
}

run_scvelo_flow_overlays_all_cells <- function(root_dir, out_dir, out_umap_dir) {
  manifest <- build_scvelo_flow_manifest(root_dir)
  if (nrow(manifest) > 0) {
    rel <- gsub("\\\\", "/", manifest$source_relative_dir)
    manifest <- manifest[rel == "01_velocity_groups/All_cells", , drop = FALSE]
  }
  if (nrow(manifest) == 0) {
    direct_h5ad <- file.path(root_dir, "01_velocity_groups", "All_cells", "scvelo_result.h5ad")
    if (file.exists(direct_h5ad) && file.info(direct_h5ad)$size > 0) {
      manifest <- data.frame(
        h5ad = normalizePath(direct_h5ad, mustWork = TRUE),
        source_relative_dir = "01_velocity_groups/All_cells",
        source_top_level = "01_velocity_groups",
        source_analysis_scope = "All_cells",
        source_task_name = "All_cells",
        output_subdir = "All_cells",
        output_prefix = "scVelo",
        label = "all_cells",
        stringsAsFactors = FALSE
      )
    }
  }
  if (nrow(manifest) > 0) {
    manifest$output_subdir <- "All_cells"
    manifest$output_prefix <- "scVelo"
    manifest$label <- "all_cells"
    manifest$read_mode <- "All_cells_only"
  }
  write_table_csv(manifest, file.path(out_dir, "scvelo_flow_overlay_manifest.csv"))
  if (nrow(manifest) == 0) {
    warning("No All_cells scvelo_result.h5ad found for flow overlays.", call. = FALSE)
    return(invisible(NULL))
  }
  cleanup_scvelo_point_to_point_outputs(out_dir)
  cleanup_scvelo_point_to_point_outputs(out_umap_dir)
  cleanup_scvelo_flat_overlay_outputs(out_dir)
  cleanup_scvelo_flat_overlay_outputs(out_umap_dir)

  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    message("Drawing All_cells scVelo flow overlay: ", i, "/", nrow(manifest), " ", manifest$label[i])
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
        result$read_mode <- "All_cells_only"
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
          read_mode = "All_cells_only",
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
  status_df <- dplyr::bind_rows(status_rows)
  write_table_csv(status_df, file.path(out_dir, "scvelo_flow_overlay_status.csv"))
  if (any(status_df$status != "ok")) {
    warning("Some All_cells scVelo flow overlays failed; see scvelo_flow_overlay_status.csv.", call. = FALSE)
  }
  invisible(NULL)
}

plot_monocle3_graph_overlays_all_cells <- function(root_dir, out_dir, out_umap_dir) {
  group_dir <- file.path(root_dir, "01_pseudotrajectory_groups")
  search_dir <- if (dir.exists(group_dir)) group_dir else file.path(root_dir, "04_by_dose_karyotype")
  dirs <- file.path(search_dir, "All_cells")
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs) == 0) {
    warning("No All_cells Monocle3 graph directory found under ", search_dir, call. = FALSE)
    return(invisible(NULL))
  }

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
    cells$source_label <- "all_cells"
    edges$source_label <- "all_cells"
    nodes$source_label <- "all_cells"
    cells_list[[idx]] <- cells
    edges_list[[idx]] <- edges
    nodes_list[[idx]] <- nodes
    manifest[[idx]] <- data.frame(
      source_label = "all_cells",
      source_dir = dir_i,
      n_cells = nrow(cells),
      n_edges = nrow(edges),
      read_mode = "All_cells_only",
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
  p <- ggplot2::ggplot(cell_all, ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2, color = .data$pseudotime)) +
    ggplot2::geom_point(size = 0.18, alpha = 0.62) +
    ggplot2::geom_segment(
      data = edge_all,
      ggplot2::aes(x = .data$UMAP_1_from, y = .data$UMAP_2_from, xend = .data$UMAP_1_to, yend = .data$UMAP_2_to),
      inherit.aes = FALSE,
      linewidth = 0.32,
      color = "black",
      alpha = 0.72,
      arrow = grid::arrow(length = grid::unit(0.055, "inches"), type = "closed")
    ) +
    ggplot2::geom_point(
      data = root_nodes,
      ggplot2::aes(x = .data$UMAP_1, y = .data$UMAP_2),
      inherit.aes = FALSE,
      shape = 8,
      size = 2.5,
      stroke = 0.8,
      color = "black"
    ) +
    ggplot2::facet_grid(Ploidy ~ Dose) +
    ggplot2::scale_color_gradientn(
      colors = c("#313695", "#74ADD1", "#FFFFBF", "#F46D43", "#A50026"),
      name = "Monocle3\npseudotime"
    ) +
    coord_square_umap(graph_umap_x, graph_umap_y) +
    ggplot2::labs(
      title = "Monocle3 all_cells pseudotime UMAP with principal graph direction",
      subtitle = "Cells are colored by Monocle3 pseudotime; arrows follow graph distance away from the root node; stars mark root nodes.",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_extra(9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0)) +
    theme_square_umap_panel()
  save_square_panel_both(p, file.path(out_dir, "Monocle3_graph_pseudotime_by_ploidy_dose"), panel_size = 3.2)
  save_square_panel_both(p, file.path(out_umap_dir, "Monocle3_graph_pseudotime_by_ploidy_dose"), panel_size = 3.2)
  invisible(NULL)
}

extra_all_cells_paga_manifest <- function(trajectory_root) {
  manifest <- paga_extra_discover_runs(trajectory_root)
  if (nrow(manifest) == 0) return(manifest)
  rel <- gsub("\\\\", "/", manifest$source_relative_dir)
  all_cells <- manifest[rel == "All_cells", , drop = FALSE]
  if (nrow(all_cells) == 0) {
    warning("No All_cells PAGA run found under ", file.path(trajectory_root, "04_paga_groups"), call. = FALSE)
    return(all_cells)
  }
  analysis_objects <- c("All_cells", "CellLine/ploidy_all", "Tumor/ploidy_all")
  rows <- lapply(analysis_objects, function(obj) {
    row <- all_cells[1, , drop = FALSE]
    row$analysis_object <- obj
    row$analysis_group <- extra_all_cells_group_name(obj)
    row$analysis_display_label <- extra_all_cells_display_label(obj)
    row$tn_filter <- if (obj == "All_cells") NA_character_ else sub("/.*$", "", obj)
    row$source_relative_dir <- "All_cells"
    row$deg_group_relative_dir <- if (obj == "All_cells") "All_cells/ploidy_all" else obj
    row$output_subdir <- obj
    row$output_prefix <- paga_extra_safe_stub(gsub("/", "_", obj, fixed = TRUE))
    row$label <- extra_all_cells_display_label(obj)
    row$read_mode <- "All_cells_only_derived"
    row
  })
  dplyr::bind_rows(rows)
}

extra_all_cells_paga_read_points <- function(run_row) {
  points <- paga_extra_read_points(run_row)
  points <- extra_all_cells_assign_tn(points)
  tn_filter <- as.character(run_row$tn_filter[1])
  if (!is.na(tn_filter) && nzchar(tn_filter)) {
    points <- points %>% dplyr::filter(.data$TN == .env$tn_filter)
  }
  if (nrow(points) == 0) {
    stop("No PAGA points remained for derived analysis object: ", run_row$analysis_object[1], call. = FALSE)
  }
  points <- extra_all_cells_update_analysis_fields(points, run_row$analysis_object[1])
  points$clusters <- factor(as.character(points$clusters), levels = paga_extra_sort_levels(points$clusters))
  points
}

extra_all_cells_paga_run_one <- function(run_row, output_root, max_points = 60000L, top_n = 10L, max_cascade_genes = 30L, max_neighbor_edges = 50000L) {
  out_dir <- file.path(output_root, run_row$output_subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  points <- extra_all_cells_paga_read_points(run_row)
  edges <- paga_extra_read_edges(run_row$paga_edges)
  valid_clusters <- unique(as.character(points$clusters))
  if (nrow(edges) > 0) {
    edges <- edges[
      as.character(edges$group_1) %in% valid_clusters &
        as.character(edges$group_2) %in% valid_clusters,
      ,
      drop = FALSE
    ]
  }
  cluster_summary <- paga_extra_cluster_summary(points, edges)
  state_summary <- paga_extra_macrostates(cluster_summary, edges)
  paga_extra_write_csv(points, file.path(out_dir, "paga_cells_used.csv"))
  paga_extra_write_csv(edges, file.path(out_dir, "paga_edges_used.csv"))
  paga_extra_write_csv(cluster_summary, file.path(out_dir, "paga_cluster_order.csv"))
  paga_extra_write_csv(state_summary, file.path(out_dir, "paga_cluster_state_summary.csv"))
  paga_extra_plot_umap_edges(
    points,
    edges,
    cluster_summary,
    run_row$label,
    file.path(out_dir, "paga_umap_edges_on_data"),
    max_points = max_points
  )
  n_neighbor_edges_plotted <- paga_extra_plot_umap_neighbor_edges(
    points,
    run_row$h5ad,
    cluster_summary,
    run_row$label,
    file.path(out_dir, "paga_umap_neighbor_edges_on_data"),
    max_points = max_points,
    max_edges = max_neighbor_edges
  )
  paga_extra_plot_dpt_umap(points, run_row$label, file.path(out_dir, "paga_dpt_pseudotime_umap"), max_points = max_points)
  paga_extra_plot_macrostates(points, state_summary, run_row$label, file.path(out_dir, "paga_macrostate_umap"), max_points = max_points)
  paga_extra_plot_terminal_states(points, state_summary, run_row$label, file.path(out_dir, "paga_terminal_candidate_umap"), max_points = max_points)
  paga_extra_plot_graph(edges, state_summary, run_row$label, file.path(out_dir, "paga_cluster_graph_ordered_by_dpt"))

  expr_status <- tryCatch(
    paga_extra_run_expression_outputs(
      run_row = run_row,
      points = points,
      cluster_order = cluster_summary,
      out_dir = out_dir,
      top_n = top_n,
      max_cascade_genes = max_cascade_genes
    ),
    error = function(e) {
      list(
        status = "failed",
        message = conditionMessage(e),
        n_top_deg = NA_integer_,
        n_expr_genes = NA_integer_,
        heatmap_rows = NA_integer_,
        cascade_rows = NA_integer_,
        n_top_deg_clusters = NA_integer_,
        n_heatmap_source_clusters = NA_integer_,
        missing_deg_clusters = NA_character_
      )
    }
  )
  data.frame(
    h5ad = run_row$h5ad,
    source_relative_dir = run_row$source_relative_dir,
    analysis_object = run_row$analysis_object,
    analysis_group = run_row$analysis_group,
    analysis_display_label = run_row$analysis_display_label,
    deg_group_relative_dir = run_row$deg_group_relative_dir,
    output_dir = out_dir,
    read_mode = run_row$read_mode,
    status = if (identical(expr_status$status, "completed")) "completed" else "partial",
    structural_status = "completed",
    expression_status = expr_status$status,
    expression_message = expr_status$message,
    n_cells = nrow(points),
    n_clusters = nrow(cluster_summary),
    n_edges = nrow(edges),
    n_neighbor_edges_plotted = n_neighbor_edges_plotted,
    n_terminal_candidates = sum(state_summary$is_terminal_candidate, na.rm = TRUE),
    n_top_deg = expr_status$n_top_deg,
    n_expr_genes = expr_status$n_expr_genes,
    heatmap_rows = expr_status$heatmap_rows,
    cascade_rows = expr_status$cascade_rows,
    n_top_deg_clusters = expr_status$n_top_deg_clusters,
    n_heatmap_source_clusters = expr_status$n_heatmap_source_clusters,
    missing_deg_clusters = expr_status$missing_deg_clusters,
    stringsAsFactors = FALSE
  )
}

run_paga_extra_analysis_all_cells_derived <- function(
  trajectory_root,
  deg_round2_root,
  output_root,
  analysis_label = "PAGA",
  top_n = 10L,
  max_cascade_genes = 30L,
  max_points = 60000L,
  max_neighbor_edges = 50000L
) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  manifest <- extra_all_cells_paga_manifest(trajectory_root)
  if (nrow(manifest) == 0) {
    empty_status <- data.frame()
    paga_extra_write_csv(manifest, file.path(output_root, "paga_run_manifest.csv"))
    paga_extra_write_csv(empty_status, file.path(output_root, "paga_extra_status.csv"))
    return(list(manifest = manifest, status = empty_status))
  }
  manifest$deg_round2_root <- deg_round2_root
  manifest$analysis_label <- analysis_label
  paga_extra_write_csv(manifest, file.path(output_root, "paga_run_manifest.csv"))

  status_rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    message("Drawing All_cells-derived PAGA outputs: ", i, "/", nrow(manifest), " ", manifest$label[i])
    run_row <- manifest[i, , drop = FALSE]
    status_rows[[i]] <- tryCatch(
      extra_all_cells_paga_run_one(
        run_row = run_row,
        output_root = output_root,
        max_points = max_points,
        top_n = top_n,
        max_cascade_genes = max_cascade_genes,
        max_neighbor_edges = max_neighbor_edges
      ),
      error = function(e) {
        out_dir <- file.path(output_root, run_row$output_subdir)
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
        data.frame(
          h5ad = run_row$h5ad,
          source_relative_dir = run_row$source_relative_dir,
          analysis_object = run_row$analysis_object,
          analysis_group = run_row$analysis_group,
          analysis_display_label = run_row$analysis_display_label,
          deg_group_relative_dir = run_row$deg_group_relative_dir,
          output_dir = out_dir,
          read_mode = run_row$read_mode,
          status = "failed",
          structural_status = "failed",
          expression_status = "not_attempted",
          expression_message = conditionMessage(e),
          n_cells = NA_integer_,
          n_clusters = NA_integer_,
          n_edges = NA_integer_,
          n_neighbor_edges_plotted = NA_integer_,
          n_terminal_candidates = NA_integer_,
          n_top_deg = NA_integer_,
          n_expr_genes = NA_integer_,
          heatmap_rows = NA_integer_,
          cascade_rows = NA_integer_,
          n_top_deg_clusters = NA_integer_,
          n_heatmap_source_clusters = NA_integer_,
          missing_deg_clusters = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    )
  }
  status <- dplyr::bind_rows(status_rows)
  paga_extra_write_csv(status, file.path(output_root, "paga_extra_status.csv"))
  if (any(status$status == "failed" | status$expression_status == "failed", na.rm = TRUE)) {
    warning("Some All_cells-derived PAGA extra outputs failed; see paga_extra_status.csv.", call. = FALSE)
  }
  list(manifest = manifest, status = status)
}
