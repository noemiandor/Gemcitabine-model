# Main-figure wrapper for the exact final-QC NUMBAT copy-number heatmap.
# The shared SI helper owns cell selection, CBS validation, chromosome
# reduction, row ordering, and annotation colors. This wrapper only binds the
# reviewed repository inputs and returns a live grob for Figure 7 assembly.

figure7_copy_number_annotation_key <- function(annotation_colors) {
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
        fontfamily = "sans", fontface = "bold", fontsize = 5.7,
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
          fontfamily = "sans", fontsize = 5.25, col = "#222222"
        )
      )
    )
  }

  children <- list(
    grid::rectGrob(gp = grid::gpar(fill = "white", col = NA)),
    title_grob("Injected origin", 0.975),
    entry_grob(
      names(annotation_colors$`Injected origin`)[[1L]],
      unname(annotation_colors$`Injected origin`[[1L]]), 0.02, 0.915
    ),
    entry_grob(
      names(annotation_colors$`Injected origin`)[[2L]],
      unname(annotation_colors$`Injected origin`[[2L]]), 0.50, 0.915
    ),
    title_grob("Dose (mg/kg)", 0.835)
  )
  dose_x <- c(0.02, 0.35, 0.68)
  dose_display_labels <- c("Vehicle", "30", "120")
  for (index in seq_along(annotation_colors$`Gemcitabine dose`)) {
    children[[length(children) + 1L]] <- entry_grob(
      dose_display_labels[[index]],
      unname(annotation_colors$`Gemcitabine dose`[[index]]),
      dose_x[[index]], 0.775
    )
  }
  children[[length(children) + 1L]] <- title_grob("Mouse", 0.695)
  mouse_names <- names(annotation_colors$Mouse)
  for (index in seq_along(annotation_colors$Mouse)) {
    column <- (index - 1L) %/% 8L
    row <- (index - 1L) %% 8L
    children[[length(children) + 1L]] <- entry_grob(
      mouse_names[[index]],
      unname(annotation_colors$Mouse[[index]]),
      0.02 + 0.49 * column,
      0.63 - 0.078 * row
    )
  }
  key <- do.call(
    grid::grobTree,
    c(children, list(name = "figure7_j_annotation_key"))
  )
  attr(key, "figure7_annotation_key_data") <- key_data
  key
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
  annotation_key <- figure7_copy_number_annotation_key(
    heatmap$annotation_colors
  )
  panel_plot <- patchwork::wrap_plots(
    patchwork::wrap_elements(full = heatmap$gtable),
    patchwork::wrap_elements(full = annotation_key),
    nrow = 1,
    widths = c(3.10, 1.15)
  )
  panel_grob <- patchwork::patchworkGrob(panel_plot)
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
    annotation_key = annotation_key,
    annotation_key_data = attr(
      annotation_key, "figure7_annotation_key_data", exact = TRUE
    ),
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
