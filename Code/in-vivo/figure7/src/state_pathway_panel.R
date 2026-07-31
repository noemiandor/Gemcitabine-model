# Frozen plotting-table contract and renderer for panel 7F.

figure7_state_required_files <- function() c(
  "panel_7F_pathway_activity_plot_data.tsv",
  "panel_7F_selected_pathway_gsea.tsv",
  "panel_7F_leading_edge_genes.tsv",
  "state_pathway_gene_ranking_complete.tsv",
  "state_pathway_gsea_complete.tsv",
  "state_pathway_sample_bin_coverage.tsv",
  "state_pathway_design_qc.tsv",
  "state_pathway_provenance.tsv"
)

figure7_validate_historical_state_reference <- function(
  path,
  config,
  verify_checksums = TRUE
) {
  expected_id <- as.character(config$state_pathways$reference_id)
  if (!dir.exists(path)) {
    figure7_stop("Missing historical panel-7F saved-state directory: ", path,
                 ". Export the byte-pinned 04i compact audit tables; the report PDF/HTML is not a substitute.")
  }
  if (!identical(basename(normalizePath(path)), expected_id)) {
    figure7_stop("Saved-state directory must use frozen reference ID ", expected_id, ": ", path)
  }
  files <- figure7_state_required_files()
  missing <- files[!file.exists(file.path(path, files))]
  if (length(missing)) figure7_stop("Panel-7F saved state is incomplete; missing: ", paste(missing, collapse = ", "))
  expected_hashes <- config$state_pathways$expected_files
  if (isTRUE(verify_checksums)) {
    for (file in files) figure7_verify_checksum(file.path(path, file), expected_hashes[[file]], file)
  }
  activity <- figure7_read_tsv(file.path(path, files[[1L]]), c(
    "collection_id", "collection_label", "collection_display_order",
    "pathway_id", "pathway_label", "pathway_display_order",
    "selected_direction", "selected_rank_within_direction",
    "pseudotime", "standardized_activity"
  ), "panel-7F activity plotting table")
  activity$pseudotime <- figure7_numeric(activity$pseudotime)
  activity$standardized_activity <- figure7_numeric(activity$standardized_activity)
  if (any(!is.finite(activity$pseudotime)) || any(!is.finite(activity$standardized_activity))) {
    figure7_stop("Panel-7F activity table contains non-finite grid/activity values")
  }
  pathways <- unique(activity[, c(
    "collection_id", "collection_label", "collection_display_order", "pathway_id", "pathway_label",
    "pathway_display_order", "selected_direction", "selected_rank_within_direction"
  )])
  if (any(duplicated(pathways[, c("collection_id", "pathway_id")]))) {
    figure7_stop("Panel-7F pathway metadata are not stable within collection/pathway keys")
  }
  expected_collections <- as.character(unlist(config$state_pathways$collections))
  observed_collections <- pathways$collection_id[order(pathways$collection_display_order)]
  observed_collections <- observed_collections[!duplicated(observed_collections)]
  if (!identical(observed_collections, expected_collections)) {
    figure7_stop("Panel-7F collection IDs/order disagree with frozen config")
  }
  counts <- table(pathways$collection_id)
  if (any(counts[expected_collections] != 8L)) {
    figure7_stop("Panel 7F requires exactly eight selected pathways per configured collection")
  }
  activity_key <- interaction(activity$collection_id, activity$pathway_id, drop = TRUE)
  reference_grid <- NULL
  for (key in levels(activity_key)) {
    local <- activity[activity_key == key, , drop = FALSE]
    if (nrow(local) != as.integer(config$state_pathways$grid_size)) {
      figure7_stop("Panel-7F pathway/grid row count must equal ", config$state_pathways$grid_size)
    }
    if (is.unsorted(local$pseudotime, strictly = TRUE)) figure7_stop("Panel-7F pseudotime grid must be strictly increasing")
    if (is.null(reference_grid)) reference_grid <- local$pseudotime
    if (!identical(local$pseudotime, reference_grid)) {
      figure7_stop("Every panel-7F pathway must use the identical pseudotime grid")
    }
  }
  selected <- figure7_read_tsv(file.path(path, "panel_7F_selected_pathway_gsea.tsv"), c(
    "collection_id", "collection_label", "collection_display_order", "pathway_id", "pathway_label",
    "pathway_display_order", "NES", "padj", "selected_direction", "selected_rank_within_direction"
  ))
  if (!setequal(paste(pathways$collection_id, pathways$pathway_id), paste(selected$collection_id, selected$pathway_id))) {
    figure7_stop("Panel-7F activity and selected-GSEA pathway keys disagree")
  }
  metadata_columns <- c("collection_id", "collection_label", "collection_display_order", "pathway_id",
                        "pathway_label", "pathway_display_order", "selected_direction",
                        "selected_rank_within_direction")
  activity_metadata <- pathways[do.call(order, pathways[c("collection_display_order", "pathway_display_order")]), metadata_columns]
  selected_metadata <- selected[do.call(order, selected[c("collection_display_order", "pathway_display_order")]), metadata_columns]
  rownames(activity_metadata) <- NULL; rownames(selected_metadata) <- NULL
  if (!identical(activity_metadata, selected_metadata)) {
    figure7_stop("Panel-7F activity and selected-GSEA metadata/order disagree")
  }
  leading <- figure7_read_tsv(file.path(path, "panel_7F_leading_edge_genes.tsv"),
    c("collection_id", "pathway_id", "leading_edge_gene_id"))
  leading$leading_edge_gene_id <- trimws(as.character(leading$leading_edge_gene_id))
  leading <- leading[!is.na(leading$leading_edge_gene_id) & nzchar(leading$leading_edge_gene_id), , drop = FALSE]
  if (any(!paste(leading$collection_id, leading$pathway_id) %in%
          paste(pathways$collection_id, pathways$pathway_id))) {
    figure7_stop("Panel-7F leading-edge table contains an unselected pathway")
  }
  leading_counts <- table(paste(leading$collection_id, leading$pathway_id))
  selected_keys <- paste(pathways$collection_id, pathways$pathway_id)
  if (any(is.na(leading_counts[selected_keys])) || any(leading_counts[selected_keys] < 1L)) {
    figure7_stop("Every selected panel-7F pathway requires nonempty leading-edge membership")
  }
  for (file in files[4:7]) {
    table <- figure7_read_tsv(file.path(path, file))
    if (!nrow(table)) figure7_stop("Panel-7F audit table is empty: ", file)
  }
  complete_gsea <- figure7_read_tsv(file.path(path, "state_pathway_gsea_complete.tsv"),
    c("collection_id", "pathway_id", "NES", "padj"))
  if (any(duplicated(complete_gsea[, c("collection_id", "pathway_id")]))) {
    figure7_stop("Complete panel-7F GSEA table has duplicated collection/pathway keys")
  }
  expected_rows <- list()
  for (collection in expected_collections) {
    local <- complete_gsea[complete_gsea$collection_id == collection &
                             is.finite(figure7_numeric(complete_gsea$padj)) &
                             is.finite(figure7_numeric(complete_gsea$NES)), , drop = FALSE]
    positive <- local[figure7_numeric(local$NES) > 0, , drop = FALSE]
    positive <- positive[order(figure7_numeric(positive$padj), -figure7_numeric(positive$NES), positive$pathway_id), , drop = FALSE]
    positive <- head(positive, as.integer(config$state_pathways$top_positive_per_collection))
    positive$selected_direction <- "positive"; positive$selected_rank_within_direction <- seq_len(nrow(positive))
    negative <- local[figure7_numeric(local$NES) < 0, , drop = FALSE]
    negative <- negative[order(figure7_numeric(negative$padj), figure7_numeric(negative$NES), negative$pathway_id), , drop = FALSE]
    negative <- head(negative, as.integer(config$state_pathways$top_negative_per_collection))
    negative$selected_direction <- "negative"; negative$selected_rank_within_direction <- seq_len(nrow(negative))
    expected_rows[[collection]] <- rbind(positive, negative)
  }
  expected_selection <- do.call(rbind, expected_rows)
  selector_columns <- c("collection_id", "pathway_id", "selected_direction", "selected_rank_within_direction")
  actual_selection <- selected[, selector_columns]
  expected_selection <- expected_selection[, selector_columns]
  actual_selection <- actual_selection[order(actual_selection$collection_id, actual_selection$selected_direction,
                                               actual_selection$selected_rank_within_direction), ]
  expected_selection <- expected_selection[order(expected_selection$collection_id, expected_selection$selected_direction,
                                                   expected_selection$selected_rank_within_direction), ]
  rownames(actual_selection) <- NULL; rownames(expected_selection) <- NULL
  if (!identical(actual_selection, expected_selection)) {
    figure7_stop("Selected panel-7F pathways do not reproduce the frozen top-four-per-sign selector")
  }
  provenance <- figure7_read_tsv(file.path(path, "state_pathway_provenance.tsv"), c("key", "value"))
  required_provenance <- c(
    "full_analysis_run_dir", "report_identifier", "code_revision_04i", "seurat_rds_sha256",
    "cellcycle_metadata_sha256", "noncellcycle_metadata_sha256", "interval_config_sha256",
    "assay", "counts_layer", "etp_method", "etp_threshold", "spline_df", "pseudotime_bins",
    "minimum_cells_per_sample_bin", "grid_size", "seed", "gene_set_source", "gene_set_release",
    "gene_set_species", "gene_set_collections", "gene_set_min_size", "gene_set_max_size",
    "expression_filter", "normalization", "observation_model", "mouse_block", "nuisance_terms",
    "treatment_by_pseudotime_interaction", "empirical_bayes", "contrast", "gsea_rank_statistic",
    "pathway_activity", "feature_species_policy", "gsea_ranking_rule", "pathway_selection_rule",
    "activity_table_sha256", "canonical_reference_id", "workflow_id", "model_id",
    "accumulated_interval", "left_neighbor_interval", "right_neighbor_interval",
    "report_html_sha256", "report_html_relative_path", "source_results_id",
    "source_activity_sha256", "source_primary_gsea_sha256", "source_leading_edge_sha256",
    "source_gene_contrast_sha256", "source_gene_resolution_sha256",
    "source_sample_bin_metadata_sha256", "source_design_audit_sha256",
    "source_primary_coverage_sha256"
  )
  if (!all(required_provenance %in% provenance$key)) {
    figure7_stop("Panel-7F provenance is missing: ", paste(setdiff(required_provenance, provenance$key), collapse = ", "))
  }
  provenance_value <- function(key) as.character(provenance$value[match(key, provenance$key)])
  runtime_path_keys <- c("report_html_path", "export_source_results_root")
  if (any(runtime_path_keys %in% provenance$key)) {
    figure7_stop(
      "Canonical panel-7F provenance must not contain runtime filesystem paths: ",
      paste(intersect(runtime_path_keys, provenance$key), collapse = ", ")
    )
  }
  character_expectations <- c(
    assay = as.character(config$state_pathways$assay),
    counts_layer = as.character(config$state_pathways$counts_layer),
    etp_method = as.character(config$etp$method),
    expression_filter = as.character(config$state_pathways$expression_filter),
    normalization = as.character(config$state_pathways$normalization),
    observation_model = as.character(config$state_pathways$observation_model),
    mouse_block = as.character(config$state_pathways$mouse_block),
    empirical_bayes = as.character(config$state_pathways$empirical_bayes),
    contrast = as.character(config$state_pathways$contrast),
    gsea_rank_statistic = as.character(config$state_pathways$gsea_rank_statistic),
    pathway_activity = as.character(config$state_pathways$pathway_activity),
    pathway_selection_rule = as.character(config$state_pathways$pathway_selector),
    canonical_reference_id = as.character(config$state_pathways$reference_id),
    workflow_id = "binning",
    model_id = as.character(config$state_pathways$model),
    accumulated_interval = "[0.30,0.49]",
    left_neighbor_interval = "[0.11,0.30)",
    right_neighbor_interval = "(0.49,0.68]",
    report_html_relative_path = "report/04i_pseudotime_state_pathways_report.html",
    source_results_id = "04i_pseudotime_state_pathways"
  )
  for (key in names(character_expectations)) {
    if (!identical(provenance_value(key), character_expectations[[key]])) {
      figure7_stop("Panel-7F provenance/config mismatch for ", key)
    }
  }
  numeric_expectations <- c(
    etp_threshold = as.numeric(config$etp$threshold), spline_df = as.numeric(config$state_pathways$spline_df),
    pseudotime_bins = as.numeric(config$state_pathways$pseudotime_bins),
    minimum_cells_per_sample_bin = as.numeric(config$state_pathways$minimum_cells_per_sample_bin),
    grid_size = as.numeric(config$state_pathways$grid_size), seed = as.numeric(config$statistics$seed)
  )
  for (key in names(numeric_expectations)) {
    if (!isTRUE(all.equal(figure7_numeric(provenance_value(key)), numeric_expectations[[key]], tolerance = 0))) {
      figure7_stop("Panel-7F provenance/config mismatch for ", key)
    }
  }
  normalize_list <- function(x) paste(trimws(unlist(strsplit(as.character(x), "[,;]"))), collapse = ",")
  if (!identical(normalize_list(provenance_value("nuisance_terms")),
                 paste(as.character(unlist(config$state_pathways$nuisance_terms)), collapse = ",")) ||
      !identical(normalize_list(provenance_value("gene_set_collections")),
                 paste(as.character(unlist(config$state_pathways$collections)), collapse = ",")) ||
      !identical(tolower(provenance_value("treatment_by_pseudotime_interaction")),
                 tolower(as.character(config$state_pathways$treatment_by_pseudotime_interaction)))) {
    figure7_stop("Panel-7F provenance/config mismatch for list or interaction settings")
  }
  declared_activity_hash <- provenance$value[match("activity_table_sha256", provenance$key)]
  observed_activity_hash <- figure7_sha256(file.path(path, "panel_7F_pathway_activity_plot_data.tsv"))
  if (!identical(declared_activity_hash, observed_activity_hash)) {
    figure7_stop("Panel-7F provenance activity_table_sha256 does not match the plotting table")
  }
  source_hash_keys <- grep("^source_.*_sha256$", provenance$key, value = TRUE)
  if (length(source_hash_keys) != 8L ||
      any(!grepl("^[0-9a-f]{64}$", vapply(source_hash_keys, provenance_value, character(1L))))) {
    figure7_stop("Panel-7F provenance must contain eight valid source-table SHA-256 values")
  }
  reference_paths <- stats::setNames(file.path(path, files), files)
  list(
    activity = activity,
    pathways = pathways,
    selected = selected,
    leading = leading,
    provenance = provenance,
    files = reference_paths,
    reference_id = expected_id,
    reference_kind = as.character(config$state_pathways$reference_kind),
    canonical_publication_allowed = FALSE,
    observed_hashes = vapply(
      reference_paths,
      figure7_sha256,
      character(1L)
    )
  )
}

figure7_validate_reviewed_state_reference <- function(
  path,
  config,
  verify_checksums = TRUE
) {
  expected_id <- as.character(config$state_pathways$reviewed_reference_id)
  expected_kind <- as.character(
    config$state_pathways$reviewed_reference_kind
  )
  if (!dir.exists(path)) {
    figure7_stop("Missing reviewed human-only panel-7F reference: ", path)
  }
  if (!identical(basename(normalizePath(path)), expected_id)) {
    figure7_stop(
      "Reviewed panel-7F reference must use ID ",
      expected_id,
      ": ",
      path
    )
  }
  files <- figure7_state_required_files()
  observed_files <- sort(list.files(path, all.files = FALSE))
  if (!identical(observed_files, sort(files))) {
    figure7_stop(
      "Reviewed panel-7F reference must contain exactly eight files; missing=",
      paste(setdiff(files, observed_files), collapse = ","),
      "; unexpected=",
      paste(setdiff(observed_files, files), collapse = ",")
    )
  }
  paths <- stats::setNames(file.path(path, files), files)
  if (any(file.info(paths)$size <= 0)) {
    figure7_stop("Reviewed panel-7F reference contains an empty file")
  }
  expected_hashes <- unlist(
    config$state_pathways$reviewed_expected_files,
    use.names = TRUE
  )
  if (!identical(sort(names(expected_hashes)), sort(files))) {
    figure7_stop("Reviewed panel-7F checksum inventory is incomplete")
  }
  observed_hashes <- vapply(paths, figure7_sha256, character(1L))
  if (isTRUE(verify_checksums)) {
    for (file in files) {
      figure7_verify_checksum(
        paths[[file]],
        expected_hashes[[file]],
        paste("reviewed panel-7F", file)
      )
    }
  }

  provenance <- figure7_read_tsv(
    paths[["state_pathway_provenance.tsv"]],
    c("key", "value")
  )
  if (anyDuplicated(provenance$key)) {
    figure7_stop("Reviewed panel-7F provenance contains duplicate keys")
  }
  value <- stats::setNames(
    as.character(provenance$value),
    provenance$key
  )
  required_review <- c(
    "reference_kind", "canonical_publication_allowed",
    "canonical_reference_id", "reviewed_source_reference_id",
    "reviewed_source_run_id", "reviewed_source_provenance_sha256",
    "reviewed_source_input_manifest_sha256",
    "reviewed_source_output_manifest_sha256",
    "reviewed_source_run_config_sha256",
    "reviewed_source_panel_7F_pdf_sha256",
    "reviewed_source_panel_7F_png_sha256",
    "reviewed_on", "reviewed_decision", "generated_reference_id"
  )
  missing_review <- setdiff(required_review, names(value))
  if (length(missing_review)) {
    figure7_stop(
      "Reviewed panel-7F provenance is missing: ",
      paste(missing_review, collapse = ", ")
    )
  }
  expected_review <- c(
    reference_kind = expected_kind,
    canonical_publication_allowed = "true",
    canonical_reference_id = expected_id,
    reviewed_source_reference_id =
      as.character(config$state_pathways$generated_reference_id),
    reviewed_source_run_id =
      "grch_human_only_v2_20260730_fdr_filtered_retry7_figure7",
    reviewed_source_provenance_sha256 =
      "c325359b1d0fe67c3ad1f13523e993cb80e5168d164d06e1b2fa1b1b584be888",
    reviewed_source_input_manifest_sha256 =
      "c7f26ef6ea0e7206e47468017bde905ac367d22dce65d5c851339f14135667d3",
    reviewed_source_output_manifest_sha256 =
      "f9053044381faefa5213a92ebe77e7edd94c285d4d668f01da5639a7da238fdf",
    reviewed_source_run_config_sha256 =
      "14306ba1b4a3f7c9c3ce4389b2c701a0db664d9110f669bcfd72834dadf7b4d6",
    reviewed_source_panel_7F_pdf_sha256 =
      "d49f6e1e408b0ebdf51f87bbaf66168ba785470e64f307ed1f3760bce41b9b8c",
    reviewed_source_panel_7F_png_sha256 =
      "bc49a6a6d44cdd1681cf0b17eab1a757b1e1a28fbb6444ca9725c9f617ce93f9",
    reviewed_on = "2026-07-30",
    generated_reference_id =
      as.character(config$state_pathways$generated_reference_id)
  )
  if (any(value[names(expected_review)] != expected_review)) {
    figure7_stop(
      "Reviewed panel-7F provenance does not identify the exact approved ",
      "retry7 candidate"
    )
  }
  expected_decision <- paste(
    "Approved exact GRCh-only retry7 result after requiring collection-wide",
    "BH-adjusted P <= 0.05 and prohibiting nonsignificant pathway backfill."
  )
  if (!identical(value[["reviewed_decision"]], expected_decision)) {
    figure7_stop("Reviewed panel-7F approval decision is invalid")
  }
  review_hash_keys <- c(
    "reviewed_source_provenance_sha256",
    "reviewed_source_input_manifest_sha256",
    "reviewed_source_output_manifest_sha256",
    "reviewed_source_run_config_sha256",
    "reviewed_source_panel_7F_pdf_sha256",
    "reviewed_source_panel_7F_png_sha256"
  )
  if (any(!grepl("^[0-9a-f]{64}$", value[review_hash_keys]))) {
    figure7_stop("Reviewed panel-7F source attestation contains an invalid hash")
  }
  provenance_values <- as.character(provenance$value)
  if (any(grepl("^(/|[A-Za-z]:[\\\\/])", provenance_values)) ||
      any(grepl("/(private|share|Users)/", provenance_values))) {
    figure7_stop(
      "Reviewed panel-7F provenance must not contain runtime filesystem paths"
    )
  }

  # Reuse the strict generated scientific contract after changing only the
  # publication identity in a temporary copy. This recomputes the BH selector,
  # validates the GRCh-only ranking/audit chain, and leaves the frozen bytes
  # untouched.
  temporary_root <- tempfile("figure7_reviewed_reference_validation_")
  temporary_reference <- file.path(
    temporary_root,
    as.character(config$state_pathways$generated_reference_id)
  )
  dir.create(temporary_reference, recursive = TRUE)
  on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
  copied <- file.copy(paths, file.path(temporary_reference, files))
  if (any(!copied)) {
    figure7_stop("Could not stage the reviewed panel-7F validation copy")
  }
  generated_provenance_path <- file.path(
    temporary_reference,
    "state_pathway_provenance.tsv"
  )
  generated_provenance <- provenance
  generated_provenance$value[
    match("reference_kind", generated_provenance$key)
  ] <- as.character(config$state_pathways$generated_reference_kind)
  generated_provenance$value[
    match("canonical_publication_allowed", generated_provenance$key)
  ] <- "false"
  figure7_write_tsv(generated_provenance, generated_provenance_path)
  scientific <- figure7_validate_generated_state_reference(
    temporary_reference,
    config,
    expected_inputs = NULL,
    config_path = NULL
  )

  collection_counts <- table(factor(
    scientific$pathways$collection_id,
    levels = as.character(unlist(config$state_pathways$collections))
  ))
  expected_counts <- c(5L, 8L, 8L)
  if (nrow(scientific$selected) != 21L ||
      nrow(scientific$activity) != 21L *
        as.integer(config$state_pathways$grid_size) ||
      !identical(as.integer(collection_counts), expected_counts) ||
      any(figure7_numeric(scientific$selected$padj) >
        figure7_generated_pathway_fdr_threshold())) {
    figure7_stop(
      "Reviewed panel-7F must contain the exact approved 21-pathway ",
      "5/8/8 FDR-significant selection"
    )
  }

  scientific$provenance <- provenance
  scientific$files <- paths
  scientific$reference_id <- expected_id
  scientific$reference_kind <- expected_kind
  scientific$canonical_publication_allowed <- TRUE
  scientific$observed_hashes <- observed_hashes
  scientific
}

figure7_validate_state_reference <- function(
  path,
  config,
  verify_checksums = TRUE
) {
  figure7_validate_reviewed_state_reference(
    path,
    config,
    verify_checksums = verify_checksums
  )
}

figure7_panel_f_plot <- function(activity, config) {
  path_order <- unique(activity[, c(
    "collection_id", "collection_display_order", "pathway_id",
    "pathway_label", "pathway_display_order"
  )])
  path_order <- path_order[order(path_order$collection_display_order, path_order$pathway_display_order), , drop = FALSE]
  collection_order <- unique(activity[, c("collection_label", "collection_display_order")])
  collection_order <- collection_order[order(collection_order$collection_display_order), , drop = FALSE]
  path_order$pathway_plot_key <- paste(
    path_order$collection_id,
    path_order$pathway_id,
    sep = "\r"
  )
  activity$pathway_plot_key <- paste(
    activity$collection_id,
    activity$pathway_id,
    sep = "\r"
  )
  if (anyDuplicated(path_order$pathway_plot_key)) {
    figure7_stop("Panel-7F collection/pathway keys must be unique")
  }
  activity$pathway_plot_key <- factor(
    activity$pathway_plot_key,
    levels = rev(path_order$pathway_plot_key)
  )
  pathway_axis_labels <- stats::setNames(
    path_order$pathway_label,
    path_order$pathway_plot_key
  )
  activity$collection_label <- factor(activity$collection_label, levels = collection_order$collection_label)
  bounds <- c(as.numeric(config$state_pathways$accumulated_interval$start),
              as.numeric(config$state_pathways$accumulated_interval$end))
  ggplot2::ggplot(activity, ggplot2::aes(pseudotime, pathway_plot_key, fill = standardized_activity)) +
    ggplot2::geom_tile() +
    ggplot2::geom_vline(xintercept = bounds, color = "black", linetype = "22", linewidth = 0.3) +
    ggplot2::facet_grid(collection_label ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", midpoint = 0) +
    ggplot2::scale_y_discrete(labels = pathway_axis_labels) +
    ggplot2::labs(title = "Pathway activity over CellCycle pseudotime", x = "Pseudotime", y = NULL, fill = "Activity") +
    ggplot2::theme_bw(base_size = 9)
}

figure7_build_f <- function(reference, output_dir, config) {
  for (file in setdiff(figure7_state_required_files(), "state_pathway_provenance.tsv")) {
    figure7_copy_file(reference$files[[file]], file.path(output_dir, "tables", file))
  }
  figure7_copy_file(reference$files[["state_pathway_provenance.tsv"]],
                    file.path(output_dir, "metadata", "state_pathway_provenance.tsv"))
  comparison <- data.frame(
    file = names(reference$files), sha256 = vapply(reference$files, figure7_sha256, character(1L)),
    checksum_match = TRUE,
    reference_id = if (!is.null(reference$reference_id)) {
      reference$reference_id
    } else {
      config$state_pathways$reviewed_reference_id
    },
    stringsAsFactors = FALSE
  )
  figure7_write_tsv(comparison, file.path(output_dir, "tables", "state_pathway_frozen_reference_comparison.tsv"))
  figure7_save_panel(figure7_panel_f_plot(reference$activity, config),
                   file.path(output_dir, "figures", config$panels$filenames[["7F"]]), 9, 8)
  invisible(reference)
}
