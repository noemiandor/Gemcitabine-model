# Main-figure wrapper for the exact final-QC NUMBAT copy-number heatmap.
# The shared SI helper owns cell selection, CBS validation, chromosome
# reduction, row ordering, and annotation colors. This wrapper only binds the
# reviewed repository inputs and returns a live grob for Figure 7 assembly.

figure7_build_copy_number_panel <- function(repo_root) {
  required_functions <- c(
    "si_copy_number_read_collection",
    "si_copy_number_harmonize",
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
    angle_col = 0
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
    plot = patchwork::wrap_elements(full = heatmap$gtable),
    heatmap = heatmap,
    cell_annotations = harmonized$cell_annotations,
    n_cells = nrow(harmonized$matrix),
    n_treated_cells = sum(treated),
    n_mice = length(unique(harmonized$cell_annotations$sample_id)),
    n_chromosomes = ncol(harmonized$matrix),
    input_paths = c(
      cell_metadata_path,
      endpoint_audit_path,
      all_ploidy_path,
      file.path(cbs_root, "cbs_manifest.tsv")
    )
  )
}
