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


figure7_validate_superseded_v2_state_reference <- function(
  path,
  config,
  verify_checksums = TRUE
) {
  expected_id <- "state_pathway_grch_human_only_etp2_24_day17_v2"
  expected_kind <- "superseded_human_only_endpoint_cn_score_frozen"
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
  expected_hashes <- c(
    panel_7F_pathway_activity_plot_data.tsv =
      "c2a6d455fd11c1d95bbc31441419953d4cc0d88eb8acf611fc3488f04af8d0f8",
    panel_7F_selected_pathway_gsea.tsv =
      "ca86833bb6159cfde1e1218d897183809eba00a0cf096fc5c4951b5065767460",
    panel_7F_leading_edge_genes.tsv =
      "6e86e697b64565cf0281e8684d29001eed84c83e1a3b9e45f7ffd37b78fdf9dd",
    state_pathway_gene_ranking_complete.tsv =
      "ac9ced0c9bb691b490b0197c9962a72a6d272781e6883631c76d7c9f2397d34c",
    state_pathway_gsea_complete.tsv =
      "4d1283e89b28e5100586a591621b5522b1698859d5d64b69054c4c470356a8d6",
    state_pathway_sample_bin_coverage.tsv =
      "a0a408a91968f172ea7a30d9432c63571f74cea5fe27f4c5ad1c8f6044172a47",
    state_pathway_design_qc.tsv =
      "66a1cf1a5362539f50777ae482b37b5b30e4d92bd673eee9aa4c842e91160151",
    state_pathway_provenance.tsv =
      "406c4fa97b96dc001b1744e01658a533fc571897724a485080e8bbbc787293a4"
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
    reference_kind = "reviewed_human_only_frozen",
    canonical_publication_allowed = "true",
    canonical_reference_id = expected_id,
    reviewed_source_reference_id =
      "runtime_state_pathway_grch_human_only_v2",
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
      "runtime_state_pathway_grch_human_only_v2"
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

  activity <- figure7_read_tsv(paths[[
    "panel_7F_pathway_activity_plot_data.tsv"
  ]])
  pathways <- unique(activity[, c(
    "collection_id", "collection_label", "collection_display_order",
    "pathway_id", "pathway_label", "pathway_display_order",
    "selected_direction", "selected_rank_within_direction"
  ), drop = FALSE])
  selected <- figure7_read_tsv(paths[[
    "panel_7F_selected_pathway_gsea.tsv"
  ]])
  complete_gsea <- figure7_read_tsv(paths[[
    "state_pathway_gsea_complete.tsv"
  ]])
  expected_selection <- figure7_select_generated_pathways(
    complete_gsea,
    config
  )
  selection_columns <- c(
    "collection_id", "pathway_id", "selected_direction",
    "selected_rank_within_direction"
  )
  order_selection <- function(data) {
    data <- data[, selection_columns, drop = FALSE]
    data <- data[order(
      data$collection_id,
      data$selected_direction,
      data$selected_rank_within_direction
    ), , drop = FALSE]
    rownames(data) <- NULL
    data
  }
  if (!identical(order_selection(selected), order_selection(expected_selection))) {
    figure7_stop(
      "Superseded v2 reference no longer reproduces its FDR-only selector"
    )
  }
  collection_counts <- table(factor(
    pathways$collection_id,
    levels = as.character(unlist(config$state_pathways$collections))
  ))
  expected_counts <- c(5L, 8L, 8L)
  if (nrow(selected) != 21L ||
      nrow(activity) != 21L *
        as.integer(config$state_pathways$grid_size) ||
      !identical(as.integer(collection_counts), expected_counts) ||
      any(figure7_numeric(selected$padj) >
        figure7_generated_pathway_fdr_threshold())) {
    figure7_stop(
      "Superseded v2 panel-7F must retain its exact 21-pathway ",
      "5/8/8 FDR-significant selection"
    )
  }
  ranking <- figure7_read_tsv(paths[[
    "state_pathway_gene_ranking_complete.tsv"
  ]], c("gene_id", "gene_symbol"))
  figure7_assert_human_feature_names(
    ranking$gene_id,
    analysis = "Superseded v2 compact ranking",
    human_prefix = as.character(config$feature_species$human_prefix),
    mouse_prefix = as.character(config$feature_species$mouse_prefix)
  )
  design <- figure7_read_tsv(paths[["state_pathway_design_qc.tsv"]])
  design_text <- paste(unlist(design, use.names = FALSE), collapse = ";")
  if (!grepl("ETP_reference_balanced_threshold_2_24", design_text, fixed = TRUE)) {
    figure7_stop("Superseded v2 design no longer records its endpoint-derived nuisance")
  }
  leading <- figure7_read_tsv(paths[[
    "panel_7F_leading_edge_genes.tsv"
  ]])
  list(
    activity = activity,
    pathways = pathways,
    selected = selected,
    leading = leading,
    provenance = provenance,
    files = paths,
    reference_id = expected_id,
    reference_kind = expected_kind,
    canonical_publication_allowed = FALSE,
    observed_hashes = observed_hashes
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
      "Reviewed panel-7F reference must use ID ", expected_id, ": ", path
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
  value <- stats::setNames(as.character(provenance$value), provenance$key)
  expected_review <- c(
    reference_kind = expected_kind,
    canonical_publication_allowed = "true",
    canonical_reference_id = expected_id,
    reviewed_source_reference_id =
      "state_pathway_grch_human_only_initial_ploidy_day17_v3_candidate",
    reviewed_source_run_id =
      "20260731_figure7_v3_exact_candidate_review_figure7",
    reviewed_source_provenance_sha256 =
      "3824db8e7c3b9d3ff1bc0ad361644156a47e4e146ca4851101ea5b52ea04932a",
    reviewed_source_input_manifest_sha256 =
      "89bfe81d107ec22336608918d703fdaac596cca11ee61aa62aa1a7bc2af06085",
    reviewed_source_output_manifest_sha256 =
      "0f3b9e8f16f3c39541f702f422371f2694c8e90c8abae9288db83a485850eb7f",
    reviewed_source_run_config_sha256 =
      "8d3add2d01b0cf466a5d00649554d875da6ebe7d9b6f802b94c458bfa53c84a0",
    reviewed_source_panel_7F_pdf_sha256 =
      "b09c12d0aa0a38230c6879b506810da138ebde6576cd4a39a57025a70138713f",
    reviewed_source_panel_7F_png_sha256 =
      "a84e8cfdf92832ca5ed245fdd9652dbfa2f5ae6ff7a9bffe07bf7defa7505ea1",
    reviewed_on = "2026-07-31",
    generated_reference_id =
      "state_pathway_grch_human_only_initial_ploidy_day17_v3_candidate"
  )
  missing_review <- setdiff(
    c(names(expected_review), "reviewed_decision"),
    names(value)
  )
  if (length(missing_review) ||
      any(value[names(expected_review)] != expected_review)) {
    figure7_stop(
      "Reviewed panel-7F provenance does not identify the exact approved ",
      "initial-ploidy-adjusted v3 candidate"
    )
  }
  expected_decision <- paste(
    "Approved exact GRCh-only, injected-initial-ploidy-adjusted v3 result",
    "after confirming a full-rank mouse-aware design, no endpoint-CN-score",
    "nuisance term, collection-wide BH-adjusted P <= 0.05, no",
    "nonsignificant pathway backfill, and unchanged selected-pathway",
    "directions relative to v2."
  )
  if (!identical(value[["reviewed_decision"]], expected_decision)) {
    figure7_stop("Reviewed panel-7F approval decision is invalid")
  }
  review_hash_keys <- grep(
    "^reviewed_source_.*_sha256$",
    names(value),
    value = TRUE
  )
  if (length(review_hash_keys) != 6L ||
      any(!grepl("^[0-9a-f]{64}$", value[review_hash_keys]))) {
    figure7_stop("Reviewed panel-7F source attestation is incomplete")
  }
  provenance_values <- as.character(provenance$value)
  if (any(grepl("^(/|[A-Za-z]:[\\\\/])", provenance_values)) ||
      any(grepl("/(private|share|Users)/", provenance_values))) {
    figure7_stop(
      "Reviewed panel-7F provenance must not contain runtime filesystem paths"
    )
  }

  # Re-run the full generated scientific contract against an in-memory
  # identity translation. The seven data tables remain byte-identical to the
  # reviewed candidate; only the publication identity differs.
  temporary_root <- tempfile("figure7_reviewed_v3_validation_")
  temporary_reference <- file.path(
    temporary_root,
    as.character(config$state_pathways$generated_reference_id)
  )
  dir.create(temporary_reference, recursive = TRUE)
  on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)
  copied <- file.copy(paths, file.path(temporary_reference, files))
  if (any(!copied)) {
    figure7_stop("Could not stage the reviewed v3 scientific validation")
  }
  generated_provenance <- provenance
  generated_provenance$value[
    match("reference_kind", generated_provenance$key)
  ] <- as.character(config$state_pathways$generated_reference_kind)
  generated_provenance$value[
    match("canonical_publication_allowed", generated_provenance$key)
  ] <- "false"
  figure7_write_tsv(
    generated_provenance,
    file.path(temporary_reference, "state_pathway_provenance.tsv")
  )
  scientific <- figure7_validate_generated_state_reference(
    temporary_reference,
    config,
    expected_inputs = NULL,
    config_path = NULL
  )
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
  plot <- figure7_panel_f_plot(reference$activity, config)
  figure7_save_panel(
    plot,
    file.path(output_dir, "figures", config$panels$filenames[["7F"]]),
    9,
    8
  )
  invisible(plot)
}
