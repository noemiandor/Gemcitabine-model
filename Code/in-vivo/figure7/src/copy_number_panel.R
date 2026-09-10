# Main-figure wrapper for the exact final-QC NUMBAT copy-number heatmap.
# The shared SI helper owns cell selection, CBS validation, chromosome
# reduction, row ordering, and annotation colors. This wrapper only binds the
# reviewed repository inputs and returns a live grob for Figure 7 assembly.

figure7_copy_number_annotation_key <- function(
  annotation_colors,
  font_scale = 1
) {
  expected_groups <- c("Injected origin", "Gemcitabine dose", "Mouse")
  if (!identical(names(annotation_colors), expected_groups) ||
      !identical(
        names(annotation_colors$`Injected origin`), c("2N", "4N")
      ) ||
      !identical(
        names(annotation_colors$`Gemcitabine dose`),
        c("Vehicle", "30 mg/kg", "120 mg/kg")
      ) ||
      length(annotation_colors$Mouse) != 16L ||
      any(!nzchar(names(annotation_colors$Mouse)))) {
    figure7_stop("Figure 7J annotation-color contract is incomplete")
  }
  if (length(font_scale) != 1L || !is.numeric(font_scale) ||
      !is.finite(font_scale) || font_scale <= 0) {
    figure7_stop("Figure 7J annotation-key font scale must be positive")
  }

  origin_display_labels <- figure7_sum159_origin_labels()
  key_data <- rbind(
    data.frame(
      group = "Injected origin",
      label = names(annotation_colors$`Injected origin`),
      color = unname(annotation_colors$`Injected origin`),
      stringsAsFactors = FALSE
    ),
    data.frame(
      group = "Gemcitabine dose",
      label = names(annotation_colors$`Gemcitabine dose`),
      color = unname(annotation_colors$`Gemcitabine dose`),
      stringsAsFactors = FALSE
    ),
    data.frame(
      group = "Mouse",
      label = names(annotation_colors$Mouse),
      color = unname(annotation_colors$Mouse),
      stringsAsFactors = FALSE
    )
  )

  title_grob <- function(label, y) {
    grid::textGrob(
      label,
      x = grid::unit(0.02, "npc"), y = grid::unit(y, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(
        fontfamily = "sans", fontface = "bold", fontsize = 5.7 * font_scale,
        col = "#222222"
      )
    )
  }
  entry_grob <- function(label, color, x, y) {
    grid::grobTree(
      grid::rectGrob(
        x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
        width = grid::unit(0.060, "npc"),
        height = grid::unit(0.038, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(fill = color, col = NA)
      ),
      grid::textGrob(
        label,
        x = grid::unit(x + 0.075, "npc"), y = grid::unit(y, "npc"),
        just = c("left", "center"),
        gp = grid::gpar(
          fontfamily = "sans", fontsize = 5.25 * font_scale, col = "#222222"
        )
      )
    )
  }

  children <- list(
    grid::rectGrob(
      width = grid::unit(0.99, "npc"),
      height = grid::unit(0.99, "npc"),
      gp = grid::gpar(fill = "white", col = "black", lwd = 0.8)
    ),
    title_grob("Injected origin", 0.975),
    entry_grob(
      unname(origin_display_labels[[
        names(annotation_colors$`Injected origin`)[[1L]]
      ]]),
      unname(annotation_colors$`Injected origin`[[1L]]), 0.02, 0.925
    ),
    entry_grob(
      unname(origin_display_labels[[
        names(annotation_colors$`Injected origin`)[[2L]]
      ]]),
      unname(annotation_colors$`Injected origin`[[2L]]), 0.02, 0.875
    ),
    title_grob("Dose (mg/kg)", 0.800)
  )
  dose_x <- c(0.02, 0.35, 0.68)
  dose_display_labels <- c("Vehicle", "30", "120")
  for (index in seq_along(annotation_colors$`Gemcitabine dose`)) {
    children[[length(children) + 1L]] <- entry_grob(
      dose_display_labels[[index]],
      unname(annotation_colors$`Gemcitabine dose`[[index]]),
      dose_x[[index]], 0.740
    )
  }
  children[[length(children) + 1L]] <- title_grob("Mouse", 0.660)
  mouse_names <- names(figure7_mouse_display_label_map())
  if (!setequal(mouse_names, names(annotation_colors$Mouse))) {
    figure7_stop("Figure 7J mouse display-label inventory is incomplete")
  }
  mouse_colors <- annotation_colors$Mouse[mouse_names]
  mouse_display_names <- figure7_mouse_display_labels(mouse_names)
  for (index in seq_along(mouse_colors)) {
    column <- (index - 1L) %/% 8L
    row <- (index - 1L) %% 8L
    children[[length(children) + 1L]] <- entry_grob(
      mouse_display_names[[index]],
      unname(mouse_colors[[index]]),
      0.02 + 0.49 * column,
      0.60 - 0.075 * row
    )
  }
  key <- do.call(
    grid::grobTree,
    c(children, list(name = "figure7_j_annotation_key"))
  )
  attr(key, "figure7_annotation_key_data") <- key_data
  attr(key, "figure7_annotation_key_origin_display_labels") <-
    origin_display_labels
  attr(key, "figure7_annotation_key_mouse_display_labels") <-
    stats::setNames(mouse_display_names, mouse_names)
  key
}

figure7_copy_number_standalone_plot <- function(
  heatmap,
  annotation_key,
  cell_annotations
) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    figure7_stop(
      "R package 'patchwork' is required for standalone Figure 7J"
    )
  }
  if (!requireNamespace("gtable", quietly = TRUE)) {
    figure7_stop("R package 'gtable' is required for standalone Figure 7J")
  }
  if (!is.list(heatmap) || is.null(heatmap$gtable) ||
      !inherits(annotation_key, "grob") || !is.data.frame(cell_annotations) ||
      !all(c("sample_id", "heatmap_row_id") %in% names(cell_annotations)) ||
      !nrow(cell_annotations)) {
    figure7_stop("Standalone Figure 7J received an incomplete heatmap contract")
  }

  sample_runs <- rle(as.character(cell_annotations$sample_id))
  sample_ids <- sample_runs$values
  sample_counts <- sample_runs$lengths
  if (anyDuplicated(sample_ids)) {
    figure7_stop(
      "Standalone Figure 7J mouse rows must form contiguous display blocks"
    )
  }
  sample_centers <- 1 -
    (cumsum(sample_counts) - sample_counts / 2) / sum(sample_counts)
  mouse_labels <- figure7_mouse_display_labels(sample_ids)
  mouse_label_children <- lapply(seq_along(mouse_labels), function(index) {
    grid::textGrob(
      mouse_labels[[index]],
      x = grid::unit(0.98, "npc"),
      y = grid::unit(sample_centers[[index]], "npc"),
      just = c("right", "center"),
      gp = grid::gpar(
        fontfamily = "sans", fontsize = 5.7, col = "#333333"
      )
    )
  })
  mouse_labels_grob <- do.call(
    grid::grobTree,
    c(
      mouse_label_children,
      list(name = "figure7_j_mouse_block_labels")
    )
  )

  heatmap_gtable <- heatmap$gtable
  matrix_layout <- heatmap_gtable$layout[
    heatmap_gtable$layout$name == "matrix", , drop = FALSE
  ]
  annotation_layout <- heatmap_gtable$layout[
    heatmap_gtable$layout$name == "row_annotation", , drop = FALSE
  ]
  legend_index <- which(heatmap_gtable$layout$name == "legend")
  if (nrow(matrix_layout) != 1L || nrow(annotation_layout) != 1L ||
      length(legend_index) != 1L) {
    figure7_stop("Standalone Figure 7J gtable layout is incomplete")
  }

  legend_colors <- heatmap$colors
  legend_n <- length(legend_colors)
  legend_breaks <- 1:6
  legend_range <- range(heatmap$breaks)
  legend_y <- (legend_breaks - legend_range[[1L]]) / diff(legend_range)
  heatmap_gtable$grobs[[legend_index]] <- grid::grobTree(
    grid::rectGrob(
      x = grid::unit(0, "npc"),
      y = grid::unit((seq_len(legend_n) - 0.5) / legend_n, "npc"),
      width = grid::unit(10, "bigpts"),
      height = grid::unit(1 / legend_n, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fill = legend_colors, col = NA)
    ),
    grid::textGrob(
      as.character(legend_breaks),
      x = grid::unit(14, "bigpts"),
      y = grid::unit(legend_y, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(fontfamily = "sans", fontsize = 8, col = "#222222")
    ),
    name = "figure7_j_full_height_copy_number_legend"
  )
  heatmap_gtable$layout$clip[[legend_index]] <- "off"
  legend_column <- heatmap_gtable$layout$l[[legend_index]]
  legend_gap_column <- legend_column - 1L
  heatmap_gtable$widths[[legend_gap_column]] <-
    heatmap_gtable$widths[[legend_gap_column]] + grid::unit(2, "mm")
  heatmap_gtable$widths[[legend_column]] <-
    heatmap_gtable$widths[[legend_column]] + grid::unit(4, "mm")

  heatmap_gtable <- gtable::gtable_add_cols(
    heatmap_gtable,
    grid::unit(21, "mm"),
    pos = annotation_layout$l - 1L
  )
  matrix_layout <- heatmap_gtable$layout[
    heatmap_gtable$layout$name == "matrix", , drop = FALSE
  ]
  annotation_layout <- heatmap_gtable$layout[
    heatmap_gtable$layout$name == "row_annotation", , drop = FALSE
  ]
  mouse_axis_column <- annotation_layout$l - 1L
  heatmap_gtable$widths[[annotation_layout$l]] <- grid::unit(44, "bigpts")
  heatmap_gtable <- gtable::gtable_add_grob(
    heatmap_gtable,
    mouse_labels_grob,
    t = matrix_layout$t,
    b = matrix_layout$b,
    l = mouse_axis_column,
    r = mouse_axis_column,
    clip = "off",
    name = "figure7_j_mouse_axis"
  )
  col_names_layout <- heatmap_gtable$layout[
    heatmap_gtable$layout$name == "col_names", , drop = FALSE
  ]
  if (nrow(col_names_layout) != 1L) {
    figure7_stop("Standalone Figure 7J chromosome-label row is unavailable")
  }
  annotation_axis <- grid::grobTree(
    grid::textGrob(
      c("Mouse", "Dose", "Original\ncell line"),
      x = grid::unit(c(1, 3, 5) / 6, "npc"),
      y = grid::unit(0.5, "npc"),
      gp = grid::gpar(
        fontfamily = "sans", fontface = "bold", fontsize = 4.6,
        lineheight = 0.82, col = "#222222"
      )
    ),
    name = "figure7_j_annotation_axis_labels"
  )
  heatmap_gtable <- gtable::gtable_add_grob(
    heatmap_gtable,
    annotation_axis,
    t = col_names_layout$t, b = col_names_layout$b,
    l = annotation_layout$l, r = annotation_layout$r,
    clip = "off",
    name = "figure7_j_annotation_axes"
  )
  heatmap_gtable <- gtable::gtable_add_rows(
    heatmap_gtable,
    grid::unit(7, "mm"),
    pos = length(heatmap_gtable$heights)
  )
  axis_row <- length(heatmap_gtable$heights)
  heatmap_gtable <- gtable::gtable_add_grob(
    heatmap_gtable,
    grid::textGrob(
      "Chromosome",
      y = grid::unit(0.92, "npc"),
      just = c("center", "top"),
      gp = grid::gpar(
        fontfamily = "sans", fontface = "bold", fontsize = 9,
        col = "#222222"
      )
    ),
    t = axis_row, b = axis_row,
    l = matrix_layout$l, r = matrix_layout$r,
    clip = "off",
    name = "figure7_j_chromosome_axis"
  )

  title_grob <- grid::grobTree(
    grid::textGrob(
      "A",
      x = grid::unit(0.01, "npc"), y = grid::unit(0.5, "npc"),
      just = c("left", "center"),
      gp = grid::gpar(
        fontfamily = "sans", fontface = "bold", fontsize = 13,
        col = "#111111"
      )
    ),
    grid::textGrob(
      "Tumor-cell copy-number landscape",
      x = grid::unit(0.5, "npc"), y = grid::unit(0.5, "npc"),
      just = c("center", "center"),
      gp = grid::gpar(
        fontfamily = "sans", fontface = "bold", fontsize = 11,
        col = "#222222"
      )
    )
  )
  y_axis_title <- grid::textGrob(
    "Tumor cells",
    rot = 90,
    gp = grid::gpar(
      fontfamily = "sans", fontface = "bold", fontsize = 9,
      col = "#222222"
    )
  )
  body <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = y_axis_title),
    patchwork::wrap_elements(full = heatmap_gtable),
    patchwork::plot_spacer(),
    patchwork::wrap_elements(full = annotation_key),
    nrow = 1,
    widths = grid::unit(
      c(8, 4.65, 7, 1.35),
      c("mm", "null", "mm", "null")
    )
  )
  standalone <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = title_grob),
    body,
    ncol = 1,
    heights = grid::unit(c(10, 1), c("mm", "null"))
  )
  patchwork::wrap_elements(full = patchwork::patchworkGrob(standalone))
}

figure7_copy_number_reorder_standalone_blocks <- function(harmonized) {
  if (!is.list(harmonized) || is.null(harmonized$matrix) ||
      !is.data.frame(harmonized$cell_annotations) ||
      !all(c("sample_id", "heatmap_row_id") %in%
        names(harmonized$cell_annotations))) {
    figure7_stop("Standalone Figure 7J block ordering received invalid data")
  }
  annotations <- harmonized$cell_annotations
  sample_ids <- unique(as.character(annotations$sample_id))
  display_labels <- figure7_mouse_display_labels(sample_ids)
  parsed <- regexec(
    "^([24]N)-(0|30|120)-M([1-4])$",
    display_labels
  )
  fields <- regmatches(display_labels, parsed)
  if (any(lengths(fields) != 4L)) {
    figure7_stop("Standalone Figure 7J mouse display labels are not sortable")
  }
  order_table <- data.frame(
    sample_id = sample_ids,
    display_label = display_labels,
    origin = vapply(fields, `[[`, character(1L), 2L),
    dose = as.integer(vapply(fields, `[[`, character(1L), 3L)),
    mouse_number = as.integer(vapply(fields, `[[`, character(1L), 4L)),
    stringsAsFactors = FALSE
  )
  order_table <- order_table[order(
    match(order_table$origin, c("2N", "4N")),
    order_table$dose,
    order_table$mouse_number
  ), , drop = FALSE]
  row_index <- unlist(lapply(order_table$sample_id, function(sample_id) {
    which(as.character(annotations$sample_id) == sample_id)
  }), use.names = FALSE)
  if (!identical(sort(row_index), seq_len(nrow(annotations)))) {
    figure7_stop("Standalone Figure 7J block ordering lost tumor-cell rows")
  }
  standalone <- harmonized
  standalone$matrix <- harmonized$matrix[row_index, , drop = FALSE]
  standalone$cell_annotations <- annotations[row_index, , drop = FALSE]
  standalone$sample_levels <- order_table$sample_id
  if (!is.null(harmonized$row_order_audit)) {
    audit_match <- match(
      standalone$cell_annotations$heatmap_row_id,
      harmonized$row_order_audit$heatmap_row_id
    )
    standalone$row_order_audit <- harmonized$row_order_audit[
      audit_match, , drop = FALSE
    ]
  }
  standalone$standalone_sample_order <- order_table
  standalone
}

figure7_build_copy_number_panel <- function(repo_root) {
  required_functions <- c(
    "si_copy_number_read_collection",
    "si_copy_number_harmonize",
    "si_copy_number_cluster_rows_within_samples",
    "si_copy_number_heatmap"
  )
  missing_functions <- required_functions[!vapply(
    required_functions, exists, logical(1L), mode = "function"
  )]
  if (length(missing_functions)) {
    figure7_stop(
      "Figure 7 copy-number helper functions are unavailable: ",
      paste(missing_functions, collapse = ", ")
    )
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    figure7_stop("R package 'patchwork' is required for Figure 7 copy-number panel")
  }

  cache_root <- file.path(repo_root, "Data", "in-vivo", "SIfigures")
  cell_metadata_path <- file.path(
    cache_root, "si_figures_cell_metadata.csv"
  )
  endpoint_audit_path <- file.path(
    cache_root, "si_figure6_endpoint_ploidy_join_audit.csv"
  )
  all_ploidy_path <- file.path(repo_root, "Data", "in-vivo", "all_ploidy.tsv")
  cbs_root <- file.path(repo_root, "Data", "in-vivo", "scRNAseq_Numbat")
  required_paths <- c(
    cell_metadata_path,
    endpoint_audit_path,
    all_ploidy_path,
    file.path(cbs_root, "cbs_manifest.tsv")
  )
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths)) {
    figure7_stop(
      "Figure 7 copy-number panel is missing input(s): ",
      paste(missing_paths, collapse = ", ")
    )
  }

  cells <- utils::read.csv(
    cell_metadata_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  endpoint_audit <- utils::read.csv(
    endpoint_audit_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  required_cell_columns <- c(
    "sample_id", "initial_ploidy", "dose", "context",
    "included_in_si_figures"
  )
  if (!all(required_cell_columns %in% names(cells))) {
    figure7_stop("Figure 7 copy-number metadata contract is incomplete")
  }
  included <- tolower(as.character(cells$included_in_si_figures)) %in%
    c("true", "t", "1")
  tumor <- cells[included & as.character(cells$context) == "Tumor", , drop = FALSE]
  sample_metadata <- unique(data.frame(
    mouse = as.character(tumor$sample_id),
    initial_ploidy = as.character(tumor$initial_ploidy),
    dose = as.character(tumor$dose),
    stringsAsFactors = FALSE
  ))
  sample_metadata <- sample_metadata[order(
    match(sample_metadata$initial_ploidy, c("2N", "4N")),
    suppressWarnings(as.numeric(sub("mg/kg$", "", sample_metadata$dose))),
    sample_metadata$mouse
  ), , drop = FALSE]
  if (nrow(tumor) != 9832L || nrow(sample_metadata) != 16L ||
      anyDuplicated(sample_metadata$mouse)) {
    figure7_stop(
      "Figure 7 copy-number panel requires the exact 9,832-cell, 16-mouse ",
      "final-QC tumor universe"
    )
  }

  collection <- si_copy_number_read_collection(
    cbs_root,
    all_ploidy_path,
    endpoint_audit,
    sample_metadata
  )
  harmonized <- si_copy_number_harmonize(collection)
  harmonized <- si_copy_number_cluster_rows_within_samples(harmonized)
  chromosome_labels <- rep("", 22L)
  chromosome_labels[c(1L, 5L, 9L, 13L, 17L, 22L)] <-
    as.character(c(1L, 5L, 9L, 13L, 17L, 22L))
  heatmap <- si_copy_number_heatmap(
    harmonized,
    title = NULL,
    labels_col = chromosome_labels,
    fontsize = 6,
    fontsize_col = 5.5,
    annotation_legend = FALSE,
    angle_col = 0,
    column_width_multiplier = 2
  )
  standalone_harmonized <-
    figure7_copy_number_reorder_standalone_blocks(harmonized)
  standalone_heatmap <- si_copy_number_heatmap(
    standalone_harmonized,
    title = NULL,
    labels_col = as.character(seq_len(22L)),
    fontsize = 8,
    fontsize_col = 7.5,
    annotation_legend = FALSE,
    angle_col = 0,
    column_width_multiplier = 1.15
  )
  annotation_key <- figure7_copy_number_annotation_key(
    heatmap$annotation_colors
  )
  standalone_annotation_key <- figure7_copy_number_annotation_key(
    standalone_heatmap$annotation_colors,
    font_scale = 1.25
  )
  annotation_key_gap_mm <- 2
  panel_plot <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = heatmap$gtable),
    patchwork::plot_spacer(),
    patchwork::wrap_elements(full = annotation_key),
    nrow = 1,
    widths = grid::unit(
      c(3.10, annotation_key_gap_mm, 1.15),
      c("null", "mm", "null")
    )
  )
  panel_grob <- patchwork::patchworkGrob(panel_plot)
  standalone_plot <- figure7_copy_number_standalone_plot(
    standalone_heatmap,
    standalone_annotation_key,
    standalone_harmonized$cell_annotations
  )
  treated <- harmonized$cell_annotations$dose_mg_per_kg > 0
  if (!identical(dim(harmonized$matrix), c(9832L, 22L)) ||
      sum(treated) != 5335L ||
      !identical(
        names(heatmap$row_annotation),
        c("Injected origin", "Gemcitabine dose", "Mouse")
      )) {
    figure7_stop(
      "Figure 7 copy-number panel failed its 9,832-cell/dose annotation contract"
    )
  }

  list(
    plot = patchwork::wrap_elements(full = panel_grob),
    heatmap = heatmap,
    standalone_plot = standalone_plot,
    standalone_heatmap = standalone_heatmap,
    standalone_annotation_key = standalone_annotation_key,
    standalone_cell_annotations = standalone_harmonized$cell_annotations,
    standalone_sample_order = standalone_harmonized$standalone_sample_order,
    standalone_column_width_multiplier =
      standalone_heatmap$column_width_multiplier,
    annotation_key = annotation_key,
    annotation_key_data = attr(
      annotation_key, "figure7_annotation_key_data", exact = TRUE
    ),
    annotation_key_origin_display_labels = attr(
      annotation_key,
      "figure7_annotation_key_origin_display_labels",
      exact = TRUE
    ),
    annotation_key_mouse_display_labels = attr(
      annotation_key,
      "figure7_annotation_key_mouse_display_labels",
      exact = TRUE
    ),
    annotation_key_gap_mm = annotation_key_gap_mm,
    cell_annotations = harmonized$cell_annotations,
    n_cells = nrow(harmonized$matrix),
    n_treated_cells = sum(treated),
    n_mice = length(unique(harmonized$cell_annotations$sample_id)),
    n_chromosomes = ncol(harmonized$matrix),
    column_width_multiplier = heatmap$column_width_multiplier,
    row_order_audit = harmonized$row_order_audit,
    row_ordering_policy = harmonized$row_ordering_policy,
    row_distance_method = harmonized$row_distance_method,
    row_linkage_method = harmonized$row_linkage_method,
    row_tie_break_method = harmonized$row_tie_break_method,
    input_paths = c(
      cell_metadata_path,
      endpoint_audit_path,
      all_ploidy_path,
      file.path(cbs_root, "cbs_manifest.tsv")
    )
  )
}
