# Shared production plots for the integrated Tumor/CellLine context.
#
# This file is the single plotting implementation for Supplementary Figure
# 4A-C/E and Supplementary Figure 7A/B.  Callers must source
# normalized_composition.R into the same environment before calling the SI4
# builder.  Display tags are configurable so the same reviewed panels can be
# reused in another composite without changing their deterministic point
# ordering.

shared_context_figure_theme <- function(base_size = 10) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", color = "#222222", size = base_size + 1
      ),
      plot.subtitle = ggplot2::element_text(
        color = "#444444", size = base_size - 1
      ),
      plot.tag = ggplot2::element_text(face = "bold", size = base_size + 3),
      legend.title = ggplot2::element_text(face = "bold"),
      strip.background = ggplot2::element_rect(
        fill = "grey94", color = "grey75", linewidth = 0.35
      ),
      strip.text = ggplot2::element_text(face = "bold", color = "#333333"),
      axis.title = ggplot2::element_text(color = "#333333")
    )
}

shared_context_umap_theme <- function(base_size = 10) {
  shared_context_figure_theme(base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(color = "grey35", linewidth = 0.35),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold", size = base_size + 2),
      legend.text = ggplot2::element_text(size = base_size + 1),
      legend.key.height = grid::unit(0.85, "lines"),
      legend.key.width = grid::unit(0.90, "lines")
    )
}

shared_context_tag_is_empty <- function(tag) {
  is.null(tag) || !length(tag) || is.na(tag[[1L]]) || !nzchar(tag[[1L]])
}

shared_context_add_tag <- function(plot, tag) {
  if (shared_context_tag_is_empty(tag)) {
    return(plot)
  }
  plot + ggplot2::labs(tag = tag)
}

shared_context_shuffle_cells <- function(data, seed) {
  set.seed(seed)
  data[sample.int(nrow(data)), , drop = FALSE]
}

shared_context_shuffle_offset <- function(key) {
  code_points <- utf8ToInt(as.character(key))
  if (length(code_points) != 1L) {
    stop("UMAP shuffle keys must be one character", call. = FALSE)
  }
  code_points[[1L]]
}

shared_context_make_umap_discrete <- function(
  data,
  field,
  colors,
  title,
  legend_title,
  tag,
  point_size,
  plot_seed,
  subtitle = NULL,
  labels = FALSE,
  repel_labels = FALSE,
  shuffle_key = tag
) {
  data <- shared_context_shuffle_cells(
    data,
    as.integer(plot_seed) + shared_context_shuffle_offset(shuffle_key)
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.76, stroke = 0) +
    ggplot2::scale_color_manual(
      values = colors,
      drop = FALSE,
      name = legend_title
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(size = 3, alpha = 1, stroke = 0)
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    shared_context_umap_theme(10)
  if (labels) {
    centers <- aggregate(cbind(UMAP_1, UMAP_2) ~ cluster, data = data, FUN = median)
    if (repel_labels) {
      if (!requireNamespace("ggrepel", quietly = TRUE)) {
        stop("Repelled UMAP labels require package: ggrepel", call. = FALSE)
      }
      plot <- plot + ggrepel::geom_label_repel(
        data = centers,
        ggplot2::aes(UMAP_1, UMAP_2, label = cluster),
        inherit.aes = FALSE,
        size = 2.5,
        label.size = 0.2,
        fill = "white",
        color = "#222222",
        alpha = 0.86,
        label.padding = grid::unit(0.10, "lines"),
        seed = as.integer(plot_seed) + shared_context_shuffle_offset(shuffle_key),
        box.padding = 0.20,
        point.padding = 0.05,
        force = 2,
        force_pull = 0.25,
        min.segment.length = 0,
        segment.color = "grey45",
        segment.size = 0.20,
        max.time = Inf,
        max.iter = 50000,
        max.overlaps = Inf
      )
    } else {
      plot <- plot + ggplot2::geom_label(
        data = centers,
        ggplot2::aes(UMAP_1, UMAP_2, label = cluster),
        inherit.aes = FALSE,
        size = 2.5,
        linewidth = 0.2,
        fill = "white",
        color = "#222222",
        alpha = 0.86,
        label.padding = grid::unit(0.10, "lines")
      )
    }
  }
  shared_context_add_tag(plot, tag)
}

shared_context_make_umap_continuous <- function(
  data,
  field,
  title,
  legend_title,
  tag,
  point_size,
  plot_seed,
  limits = NULL,
  subtitle = NULL,
  diverging = FALSE,
  shuffle_key = tag
) {
  data <- shared_context_shuffle_cells(
    data,
    as.integer(plot_seed) + 1000L +
      shared_context_shuffle_offset(shuffle_key)
  )
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(UMAP_1, UMAP_2, color = .data[[field]])
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.80, stroke = 0) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = title,
      subtitle = subtitle,
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    shared_context_umap_theme(10)
  if (diverging) {
    plot <- plot + ggplot2::scale_color_gradient2(
      low = "#2C7BB6",
      mid = "white",
      high = "#D7191C",
      midpoint = 0,
      limits = limits,
      name = legend_title
    )
  } else {
    plot <- plot + ggplot2::scale_color_gradient(
      low = "#2C7BB6",
      high = "#D7191C",
      limits = limits,
      name = legend_title
    )
  }
  shared_context_add_tag(plot, tag)
}

shared_context_validate_tags <- function(tags, required_names, label) {
  if (is.null(names(tags)) ||
      anyDuplicated(names(tags)) ||
      !setequal(names(tags), required_names) ||
      length(tags) != length(required_names)) {
    stop(
      label,
      " tags must be a named vector containing exactly: ",
      paste(required_names, collapse = ", "),
      call. = FALSE
    )
  }
  stats::setNames(as.character(tags[required_names]), required_names)
}

shared_context_build_si4_panels <- function(
  data,
  cluster_levels,
  cluster_colors,
  ploidy_colors,
  context_colors,
  point_size,
  plot_seed,
  tags = c(
    cluster = "A",
    initial_ploidy = "B",
    context = "C",
    composition = "E"
  ),
  repel_cluster_labels = FALSE,
  composition_builder = make_normalized_composition_plot
) {
  tags <- shared_context_validate_tags(
    tags,
    c("cluster", "initial_ploidy", "context", "composition"),
    "SI4 context-panel"
  )
  if (!is.function(composition_builder)) {
    stop(
      paste(
        "SI4 context-panel construction requires",
        "make_normalized_composition_plot()"
      ),
      call. = FALSE
    )
  }

  plots <- list(
    cluster = shared_context_make_umap_discrete(
      data = data,
      field = "cluster",
      colors = cluster_colors,
      title = "UMAP by cluster",
      legend_title = "Cluster",
      tag = tags[["cluster"]],
      point_size = point_size,
      plot_seed = plot_seed,
      subtitle = sprintf(
        "Tumor and CellLine cells; n = %s",
        format(nrow(data), big.mark = ",")
      ),
      labels = TRUE,
      repel_labels = repel_cluster_labels,
      shuffle_key = "A"
    ),
    initial_ploidy = shared_context_make_umap_discrete(
      data = data,
      field = "initial_ploidy",
      colors = ploidy_colors,
      title = "UMAP by initial ploidy",
      legend_title = "Initial ploidy",
      tag = tags[["initial_ploidy"]],
      point_size = point_size,
      plot_seed = plot_seed,
      subtitle = "Tumor and CellLine cells",
      shuffle_key = "B"
    ),
    context = shared_context_make_umap_discrete(
      data = data,
      field = "context",
      colors = context_colors,
      title = "UMAP by Tumor/CellLine context",
      legend_title = "Context",
      tag = tags[["context"]],
      point_size = point_size,
      plot_seed = plot_seed,
      subtitle = "Tumor and CellLine cells",
      shuffle_key = "C"
    )
  )
  composition_result <- composition_builder(
    data = data,
    unit_col = "mouse",
    cluster_col = "cluster",
    group_col = "context",
    cluster_levels = cluster_levels,
    group_levels = c("Tumor", "CellLine"),
    fill_colors = context_colors,
    bar_axis = "cluster",
    strata_cols = "initial_ploidy",
    title = "Equal-sample cluster composition by context",
    subtitle = paste(
      "Each tumor or CellLine sample contributes equally within context;",
      "stars mark context enrichment within ploidy at BH FDR <= 0.05"
    ),
    x_title = "Cluster",
    y_title = "Mean within-sample cluster proportion",
    legend_title = "Context",
    tag = tags[["composition"]],
    theme_function = shared_context_figure_theme,
    x_text_angle = 0
  )
  if (shared_context_tag_is_empty(tags[["composition"]])) {
    composition_result$plot <- composition_result$plot +
      ggplot2::labs(tag = NULL)
  }
  plots$composition <- composition_result$plot

  list(
    plots = plots,
    composition_result = composition_result
  )
}

shared_context_reader_pathway_labels <- function(labels) {
  replacements <- c(
    "E2f Targets" = "E2F targets",
    "G2m Checkpoint" = "G2/M checkpoint",
    "Myc Targets V1" = "MYC targets V1",
    "Oxidative Phosphorylation" = "Oxidative phosphorylation",
    "Mitotic Spindle" = "Mitotic spindle",
    "Tnfa Signaling via Nfkb" = "TNF-alpha signaling via NF-kB",
    "Unfolded Protein Response" = "Unfolded protein response",
    "Interferon Gamma Response" = "Interferon-gamma response",
    "Uv Response Dn" = "UV response down",
    "Epithelial Mesenchymal Transition" = "Epithelial-mesenchymal transition",
    "Mtorc1 Signaling" = "mTORC1 signaling",
    "Interferon Alpha Response" = "Interferon-alpha response",
    "Xenobiotic Metabolism" = "Xenobiotic metabolism"
  )
  output <- as.character(labels)
  matched <- match(output, names(replacements))
  output[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  output
}

shared_context_build_heatmap <- function(
  matrix_data,
  title,
  diverging,
  tag,
  fontsize_row = 8,
  fontsize_col = 7,
  reader_labels = FALSE,
  mark_zero_missing = FALSE,
  treeheight_row = 50,
  treeheight_col = 50
) {
  tag_prefix <- if (shared_context_tag_is_empty(tag)) {
    ""
  } else {
    paste0(tag[[1L]], "  ")
  }
  arguments <- list(
    mat = matrix_data,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    border_color = NA,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    angle_col = 45,
    treeheight_row = treeheight_row,
    treeheight_col = treeheight_col,
    labels_row = if (reader_labels) {
      shared_context_reader_pathway_labels(rownames(matrix_data))
    } else {
      rownames(matrix_data)
    },
    silent = TRUE
  )
  display_title <- paste0(tag_prefix, title)
  if (nzchar(display_title)) arguments$main <- display_title
  if (mark_zero_missing) {
    arguments$display_numbers <- ifelse(matrix_data == 0, "x", "")
    arguments$number_color <- "#777777"
    arguments$fontsize_number <- 7
  }
  if (diverging) {
    max_abs <- max(abs(matrix_data))
    arguments$color <- grDevices::colorRampPalette(
      c("#2C7BB6", "white", "#D7191C")
    )(101)
    arguments$breaks <- seq(
      -max_abs,
      max_abs,
      length.out = length(arguments$color) + 1L
    )
  }
  do.call(pheatmap::pheatmap, arguments)
}

shared_context_heatmap_plot <- function(
  matrix_data,
  title,
  diverging,
  tag,
  fontsize_row = 8,
  fontsize_col = 7,
  reader_labels = FALSE,
  mark_zero_missing = FALSE,
  treeheight_row = 50,
  treeheight_col = 50
) {
  heatmap <- shared_context_build_heatmap(
    matrix_data,
    title,
    diverging,
    tag,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    reader_labels = reader_labels,
    mark_zero_missing = mark_zero_missing,
    treeheight_row = treeheight_row,
    treeheight_col = treeheight_col
  )
  patchwork::wrap_elements(full = heatmap$gtable)
}

shared_context_build_si7_heatmap_panels <- function(
  ora_matrix,
  gsea_matrix,
  tags = c(ora = "A", gsea = "B")
) {
  tags <- shared_context_validate_tags(
    tags,
    c("ora", "gsea"),
    "SI7 heatmap-panel"
  )
  list(
    plots = list(
      ora = shared_context_heatmap_plot(
        ora_matrix,
        "Cluster Hallmark ORA annotation score",
        FALSE,
        tags[["ora"]]
      ),
      gsea = shared_context_heatmap_plot(
        gsea_matrix,
        "Cluster Hallmark GSEA NES",
        TRUE,
        tags[["gsea"]]
      )
    )
  )
}
