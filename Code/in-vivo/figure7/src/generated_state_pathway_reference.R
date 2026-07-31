# Runtime-generated panel-7F reference contract.
#
# Generated GRCh-only references are run-scoped, explicitly noncanonical, and
# validated here. The reviewed-v2 validator reuses this scientific contract
# only after independently enforcing its exact frozen bytes and approval
# identity; the historical mixed-v1 validator remains separate.

figure7_generated_state_required_files <- function() {
  figure7_state_required_files()
}

figure7_validate_generated_state_reference <- function(
  path,
  config,
  expected_inputs = NULL,
  config_path = NULL
) {
  expected_id <- as.character(config$state_pathways$generated_reference_id)
  expected_kind <- as.character(config$state_pathways$generated_reference_kind)
  if (!dir.exists(path)) {
    figure7_stop("Missing generated panel-7F reference: ", path)
  }
  if (!identical(basename(normalizePath(path)), expected_id)) {
    figure7_stop(
      "Generated panel-7F reference must use ID ",
      expected_id,
      ": ",
      path
    )
  }
  files <- figure7_generated_state_required_files()
  observed_files <- sort(list.files(path, all.files = FALSE))
  if (!identical(observed_files, sort(files))) {
    figure7_stop(
      "Generated panel-7F reference must contain exactly eight files; missing=",
      paste(setdiff(files, observed_files), collapse = ","),
      "; unexpected=",
      paste(setdiff(observed_files, files), collapse = ",")
    )
  }
  paths <- stats::setNames(file.path(path, files), files)
  if (any(file.info(paths)$size <= 0)) {
    figure7_stop("Generated panel-7F reference contains an empty file")
  }
  reference_files <- figure7_state_required_files()
  reference_paths <- paths[reference_files]
  observed_hashes <- vapply(reference_paths, figure7_sha256, character(1L))

  activity <- figure7_read_tsv(
    paths[["panel_7F_pathway_activity_plot_data.tsv"]],
    c(
      "collection_id", "collection_label", "collection_display_order",
      "pathway_id", "pathway_label", "pathway_display_order",
      "selected_direction", "selected_rank_within_direction",
      "pseudotime", "standardized_activity"
    )
  )
  activity$pseudotime <- figure7_numeric(activity$pseudotime)
  activity$standardized_activity <- figure7_numeric(
    activity$standardized_activity
  )
  if (any(!is.finite(activity$pseudotime)) ||
      any(!is.finite(activity$standardized_activity))) {
    figure7_stop("Generated panel-7F activity contains nonfinite values")
  }
  pathway_columns <- c(
    "collection_id", "collection_label", "collection_display_order",
    "pathway_id", "pathway_label", "pathway_display_order",
    "selected_direction", "selected_rank_within_direction"
  )
  pathways <- unique(activity[, pathway_columns, drop = FALSE])
  pathway_keys <- paste(pathways$collection_id, pathways$pathway_id, sep = "\r")
  collections <- as.character(unlist(config$state_pathways$collections))
  max_per_collection <-
    as.integer(config$state_pathways$top_positive_per_collection) +
    as.integer(config$state_pathways$top_negative_per_collection)
  collection_counts <- table(factor(
    pathways$collection_id,
    levels = collections
  ))
  observed_collections <- pathways$collection_id[
    order(pathways$collection_display_order)
  ]
  observed_collections <- observed_collections[!duplicated(observed_collections)]
  collection_order <- suppressWarnings(as.integer(
    pathways$collection_display_order
  ))
  pathway_order <- suppressWarnings(as.integer(pathways$pathway_display_order))
  if (!nrow(pathways) || anyDuplicated(pathway_keys) ||
      !identical(observed_collections, collections) ||
      any(collection_counts < 1L) ||
      any(collection_counts > max_per_collection) ||
      anyNA(collection_order) ||
      any(collection_order != match(pathways$collection_id, collections)) ||
      anyNA(pathway_order) ||
      !identical(sort(pathway_order), seq_len(nrow(pathways)))) {
    figure7_stop("Generated panel-7F collection contract is invalid")
  }
  activity_ordered <- activity[
    order(activity$pseudotime),
    ,
    drop = FALSE
  ]
  activity_groups <- split(
    activity_ordered,
    paste(
      activity_ordered$collection_id,
      activity_ordered$pathway_id,
      sep = "\r"
    )
  )
  grids <- lapply(activity_groups, function(data) data$pseudotime)
  if (any(lengths(grids) != as.integer(config$state_pathways$grid_size)) ||
      !all(vapply(grids, identical, logical(1L), grids[[1L]])) ||
      any(diff(grids[[1L]]) <= 0) ||
      !identical(grids[[1L]][[1L]], 0) ||
      !identical(grids[[1L]][[length(grids[[1L]])]], 1)) {
    figure7_stop(
      "Generated panel-7F pathways must share a strictly increasing ",
      "configured 0-to-1 pseudotime grid"
    )
  }

  selected <- figure7_read_tsv(
    paths[["panel_7F_selected_pathway_gsea.tsv"]],
    c(
      pathway_columns,
      "NES",
      "padj"
    )
  )
  selected$NES <- figure7_numeric(selected$NES)
  selected$padj <- figure7_numeric(selected$padj)
  selected_keys <- paste(
    selected$collection_id,
    selected$pathway_id,
    sep = "\r"
  )
  if (nrow(selected) != nrow(pathways) ||
      anyDuplicated(selected_keys) ||
      !setequal(pathway_keys, selected_keys)) {
    figure7_stop("Generated panel-7F activity/selection pathway keys disagree")
  }
  if (any(!is.finite(selected$NES)) ||
      any(selected$NES == 0) ||
      any(!is.finite(selected$padj)) ||
      any(selected$padj < 0) ||
      any(selected$padj >
        figure7_generated_pathway_fdr_threshold()) ||
      any(
        (selected$NES > 0 & selected$selected_direction != "positive") |
          (selected$NES < 0 & selected$selected_direction != "negative")
      )) {
    figure7_stop(
      "Generated panel-7F selected pathways require finite nonzero NES, ",
      "sign-consistent direction, and BH-adjusted P <= 0.05"
    )
  }
  ordered_activity <- pathways[
    order(pathways$collection_display_order, pathways$pathway_display_order),
    pathway_columns,
    drop = FALSE
  ]
  ordered_selected <- selected[
    order(selected$collection_display_order, selected$pathway_display_order),
    pathway_columns,
    drop = FALSE
  ]
  rownames(ordered_activity) <- NULL
  rownames(ordered_selected) <- NULL
  if (!identical(ordered_activity, ordered_selected)) {
    figure7_stop("Generated panel-7F pathway metadata/order disagree")
  }

  leading <- figure7_read_tsv(
    paths[["panel_7F_leading_edge_genes.tsv"]],
    c("collection_id", "pathway_id", "leading_edge_gene_id")
  )
  leading$leading_edge_gene_id <- trimws(
    as.character(leading$leading_edge_gene_id)
  )
  leading <- leading[nzchar(leading$leading_edge_gene_id), , drop = FALSE]
  leading_keys <- paste(
    leading$collection_id,
    leading$pathway_id,
    sep = "\r"
  )
  if (any(!leading_keys %in% pathway_keys) ||
      any(!pathway_keys %in% leading_keys)) {
    figure7_stop(
      "Generated panel-7F leading edges must cover only/all selected pathways"
    )
  }

  ranking <- figure7_read_tsv(
    paths[["state_pathway_gene_ranking_complete.tsv"]],
    c("gene_id", "gene_symbol")
  )
  complete_gsea <- figure7_read_tsv(
    paths[["state_pathway_gsea_complete.tsv"]],
    c(
      "collection_id", "pathway_id", "pval", "padj", "ES", "NES",
      "size", "nPermSimple", "retry_round"
    )
  )
  coverage <- figure7_read_tsv(
    paths[["state_pathway_sample_bin_coverage.tsv"]]
  )
  design <- figure7_read_tsv(paths[["state_pathway_design_qc.tsv"]])
  if (!nrow(ranking) || !nrow(complete_gsea) || !nrow(coverage) || !nrow(design)) {
    figure7_stop("Generated panel-7F audit chain contains an empty table")
  }
  if (nrow(ranking) < 10000L || anyDuplicated(ranking$gene_symbol)) {
    figure7_stop(
      "Generated panel-7F human-only ranking is implausibly small or ",
      "contains duplicate gene symbols"
    )
  }
  figure7_assert_human_feature_names(
    as.character(ranking$gene_id),
    analysis = "Generated panel-7F compact ranking",
    human_prefix = as.character(config$feature_species$human_prefix),
    mouse_prefix = as.character(config$feature_species$mouse_prefix)
  )
  if (anyDuplicated(complete_gsea[, c("collection_id", "pathway_id")])) {
    figure7_stop("Generated complete GSEA contains duplicate pathway keys")
  }
  finite_columns <- c("pval", "padj", "ES", "NES", "size")
  if (any(!vapply(
        complete_gsea[finite_columns],
        function(column) all(is.finite(figure7_numeric(column))),
        logical(1L)
      ))) {
    figure7_stop("Generated complete GSEA contains unresolved nonfinite values")
  }
  initial_nperm <- as.integer(config$state_pathways$gsea_nperm_simple)
  maximum_nperm <- as.integer(
    config$state_pathways$gsea_nperm_simple_max
  )
  retry_multiplier <- as.integer(
    config$state_pathways$gsea_nperm_simple_multiplier
  )
  allowed_nperm <- initial_nperm
  while (tail(allowed_nperm, 1L) < maximum_nperm) {
    allowed_nperm <- c(
      allowed_nperm,
      min(tail(allowed_nperm, 1L) * retry_multiplier, maximum_nperm)
    )
  }
  observed_nperm <- as.integer(figure7_numeric(
    complete_gsea$nPermSimple
  ))
  observed_round <- as.integer(figure7_numeric(
    complete_gsea$retry_round
  ))
  expected_round <- match(observed_nperm, allowed_nperm) - 1L
  if (anyNA(observed_nperm) || anyNA(observed_round) ||
      anyNA(expected_round) ||
      any(observed_round != expected_round)) {
    figure7_stop("Generated complete GSEA has invalid adaptive-retry metadata")
  }
  for (collection in collections) {
    local <- complete_gsea[
      complete_gsea$collection_id == collection,
      ,
      drop = FALSE
    ]
    if (!nrow(local)) {
      figure7_stop("Generated complete GSEA is missing collection ", collection)
    }
    expected_padj <- stats::p.adjust(
      figure7_numeric(local$pval),
      method = "BH"
    )
    if (!isTRUE(all.equal(
      figure7_numeric(local$padj),
      expected_padj,
      tolerance = 1e-12,
      check.attributes = FALSE
    ))) {
      figure7_stop(
        "Generated complete GSEA adjusted P values are not ",
        "collection-wide BH values for ",
        collection
      )
    }
  }

  expected_selection <- figure7_select_generated_pathways(
    complete_gsea,
    config
  )
  selector_columns <- c(
    "collection_id", "pathway_id", "selected_direction",
    "selected_rank_within_direction"
  )
  selector_order <- function(data) {
    data <- data[, selector_columns, drop = FALSE]
    data <- data[
      order(
        data$collection_id,
        data$selected_direction,
        data$selected_rank_within_direction
      ),
      ,
      drop = FALSE
    ]
    rownames(data) <- NULL
    data
  }
  if (!identical(selector_order(selected), selector_order(expected_selection))) {
    figure7_stop(
      "Generated panel-7F selection does not reproduce the configured selector"
    )
  }
  selected_display_order <- order(figure7_numeric(
    selected$pathway_display_order
  ))
  selected_display_keys <- selected_keys[selected_display_order]
  expected_display_keys <- paste(
    expected_selection$collection_id,
    expected_selection$pathway_id,
    sep = "\r"
  )
  if (!identical(selected_display_keys, expected_display_keys)) {
    figure7_stop(
      "Generated panel-7F pathway display order disagrees with the selector"
    )
  }
  selected_match <- match(selected_keys, paste(
    expected_selection$collection_id,
    expected_selection$pathway_id,
    sep = "\r"
  ))
  if (anyNA(selected_match) ||
      !isTRUE(all.equal(
        selected$NES,
        figure7_numeric(expected_selection$NES[selected_match]),
        tolerance = 1e-12,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        selected$padj,
        figure7_numeric(expected_selection$padj[selected_match]),
        tolerance = 1e-12,
        check.attributes = FALSE
      ))) {
    figure7_stop(
      "Generated panel-7F selected statistics disagree with complete GSEA"
    )
  }

  provenance <- figure7_read_tsv(
    paths[["state_pathway_provenance.tsv"]],
    c("key", "value")
  )
  if (anyDuplicated(provenance$key)) {
    figure7_stop("Generated panel-7F provenance contains duplicate keys")
  }
  value <- stats::setNames(as.character(provenance$value), provenance$key)
  required <- c(
    "reference_kind", "canonical_publication_allowed",
    "generated_reference_id", "source_code_revision",
    "exporter_script_sha256", "exporter_common_io_sha256",
    "seurat_rds_sha256", "cellcycle_metadata_sha256",
    "noncellcycle_metadata_sha256", "assay", "counts_layer", "etp_method",
    "etp_threshold", "spline_df", "pseudotime_bins",
    "minimum_cells_per_sample_bin", "grid_size", "seed", "gene_set_source",
    "gene_set_species", "gene_set_collections", "gene_set_release",
    "msigdbr_package_version", "gene_set_membership_sha256",
    "figure7_config_sha256", "figure7_config_contract_sha256",
    "feature_species_policy_id", "feature_species_policy",
    "human_feature_prefix", "mouse_feature_prefix",
    "unknown_feature_policy", "n_input_features",
    "n_human_features_retained", "n_mouse_features_excluded",
    "n_ambiguous_features", "feature_species_audit_sha256",
    "feature_species_policy_code_sha256",
    "gsea_nperm_simple", "gsea_nperm_simple_max",
    "gsea_nperm_simple_multiplier", "gsea_nperm_simple_usage",
    "gsea_adaptive_retry_rule",
    "pathway_selection_fdr_threshold", "pathway_selection_rule",
    "activity_table_sha256"
  )
  missing <- setdiff(required, names(value))
  if (length(missing)) {
    figure7_stop(
      "Generated panel-7F provenance is missing: ",
      paste(missing, collapse = ", ")
    )
  }
  if (!identical(value[["reference_kind"]], expected_kind) ||
      !identical(value[["canonical_publication_allowed"]], "false") ||
      !identical(value[["generated_reference_id"]], expected_id)) {
    figure7_stop("Generated panel-7F reference identity/publishing guard is invalid")
  }
  hash_keys <- c(
    "seurat_rds_sha256", "cellcycle_metadata_sha256",
    "noncellcycle_metadata_sha256", "gene_set_membership_sha256",
    "figure7_config_sha256", "activity_table_sha256",
    "feature_species_audit_sha256",
    "feature_species_policy_code_sha256",
    "exporter_script_sha256", "exporter_common_io_sha256"
  )
  if (any(!grepl("^[0-9a-f]{64}$", value[hash_keys])) ||
      !grepl("^sha256:[0-9a-f]{64}$", value[["source_code_revision"]])) {
    figure7_stop("Generated panel-7F provenance contains an invalid lineage hash")
  }
  expected_values <- c(
    assay = as.character(config$state_pathways$assay),
    counts_layer = as.character(config$state_pathways$counts_layer),
    etp_method = as.character(config$etp$method),
    gene_set_release = as.character(config$gene_sets$database_release),
    msigdbr_package_version = as.character(config$gene_sets$package_version),
    gene_set_species = "Homo sapiens",
    gene_set_collections = paste(collections, collapse = ","),
    feature_species_policy_id =
      as.character(config$feature_species$policy_id),
    feature_species_policy =
      as.character(config$state_pathways$generated_feature_species_policy),
    human_feature_prefix =
      as.character(config$feature_species$human_prefix),
    mouse_feature_prefix =
      as.character(config$feature_species$mouse_prefix),
    unknown_feature_policy =
      as.character(config$feature_species$unknown_feature_policy),
    pathway_selection_rule = figure7_generated_pathway_selection_rule(),
    gsea_nperm_simple_usage = {
      usage <- table(observed_nperm)
      paste(
        names(usage),
        as.integer(usage),
        sep = ":",
        collapse = ","
      )
    },
    gsea_adaptive_retry_rule =
      paste(
        "retry only unresolved pathways at geometric nPermSimple",
        "increments; merge by pathway; recompute collection-wide BH;",
        "fail closed at cap"
      )
  )
  if (any(value[names(expected_values)] != expected_values)) {
    figure7_stop("Generated panel-7F provenance/config values disagree")
  }
  species_counts <- suppressWarnings(as.integer(value[c(
    "n_input_features", "n_human_features_retained",
    "n_mouse_features_excluded", "n_ambiguous_features"
  )]))
  names(species_counts) <- c("input", "human", "mouse", "ambiguous")
  if (anyNA(species_counts) ||
      species_counts[["human"]] < 10000L ||
      species_counts[["mouse"]] < 0L ||
      species_counts[["ambiguous"]] != 0L ||
      species_counts[["input"]] !=
        species_counts[["human"]] + species_counts[["mouse"]]) {
    figure7_stop("Generated panel-7F species-audit counts are invalid")
  }
  expected_numeric <- c(
    etp_threshold = as.numeric(config$etp$threshold),
    spline_df = as.numeric(config$state_pathways$spline_df),
    pseudotime_bins = as.numeric(config$state_pathways$pseudotime_bins),
    minimum_cells_per_sample_bin = as.numeric(
      config$state_pathways$minimum_cells_per_sample_bin
    ),
    grid_size = as.numeric(config$state_pathways$grid_size),
    seed = as.numeric(config$statistics$seed),
    gsea_nperm_simple =
      as.numeric(config$state_pathways$gsea_nperm_simple),
    gsea_nperm_simple_max =
      as.numeric(config$state_pathways$gsea_nperm_simple_max),
    gsea_nperm_simple_multiplier =
      as.numeric(config$state_pathways$gsea_nperm_simple_multiplier),
    pathway_selection_fdr_threshold =
      figure7_generated_pathway_fdr_threshold()
  )
  for (key in names(expected_numeric)) {
    if (!isTRUE(all.equal(
      figure7_numeric(value[[key]]),
      expected_numeric[[key]],
      tolerance = 0
    ))) {
      figure7_stop("Generated panel-7F numeric provenance mismatch for ", key)
    }
  }
  if (!identical(
    value[["activity_table_sha256"]],
    figure7_sha256(paths[["panel_7F_pathway_activity_plot_data.tsv"]])
  )) {
    figure7_stop("Generated panel-7F activity hash is invalid")
  }
  compact_hash_keys <- paste0(
    sub(
      "[.]tsv$",
      "",
      setdiff(
        figure7_state_required_files(),
        "state_pathway_provenance.tsv"
      )
    ),
    "_sha256"
  )
  compact_files <- setdiff(
    figure7_state_required_files(),
    "state_pathway_provenance.tsv"
  )
  if (!all(compact_hash_keys %in% names(value)) ||
      !identical(
        unname(value[compact_hash_keys]),
        unname(vapply(paths[compact_files], figure7_sha256, character(1L)))
      )) {
    figure7_stop(
      "Generated panel-7F provenance does not authenticate all seven ",
      "compact data tables"
    )
  }
  if (!identical(
    value[["figure7_config_contract_sha256"]],
    figure7_state_config_contract_sha256(config)
  )) {
    figure7_stop(
      "Generated panel-7F consumed-config contract does not match this run"
    )
  }
  if (!is.null(config_path)) {
    policy_path <- file.path(
      dirname(config_path), "src", "feature_species_policy.R"
    )
    if (!file.exists(policy_path) ||
        !identical(
          value[["feature_species_policy_code_sha256"]],
          figure7_sha256(policy_path)
        )) {
      figure7_stop(
        "Generated panel-7F feature-species implementation changed"
      )
    }
    code_paths <- c(
      exporter_script_sha256 =
        file.path(dirname(config_path), "export_state_pathway_reference.R"),
      exporter_common_io_sha256 =
        file.path(dirname(config_path), "src", "common_io.R")
    )
    if (any(!file.exists(code_paths)) ||
        !identical(
          unname(value[names(code_paths)]),
          unname(vapply(code_paths, figure7_sha256, character(1L)))
        )) {
      figure7_stop(
        "Generated panel-7F exporter/selector implementation changed"
      )
    }
  }
  if (!is.null(expected_inputs)) {
    required_input_names <- c("cellcycle", "noncellcycle", "seurat_rds")
    if (!all(required_input_names %in% names(expected_inputs))) {
      figure7_stop(
        "Generated panel-7F validation requires cellcycle, noncellcycle, ",
        "and seurat_rds lineage inputs"
      )
    }
    available <- vapply(
      expected_inputs[required_input_names],
      function(input) {
        input <- as.character(input)
        length(input) == 1L && nzchar(input) && file.exists(input)
      },
      logical(1L)
    )
    lineage_keys <- c(
      cellcycle = "cellcycle_metadata_sha256",
      noncellcycle = "noncellcycle_metadata_sha256",
      seurat_rds = "seurat_rds_sha256"
    )
    expected_hashes <- vapply(
      expected_inputs[required_input_names][available],
      figure7_sha256,
      character(1L)
    )
    names(expected_hashes) <- unname(
      lineage_keys[names(expected_hashes)]
    )
    if (any(value[names(expected_hashes)] != expected_hashes)) {
      figure7_stop(
        "Generated panel-7F reference was built from different upstream inputs"
      )
    }
  }

  list(
    activity = activity,
    pathways = pathways,
    selected = selected,
    leading = leading,
    provenance = provenance,
    files = reference_paths,
    reference_id = expected_id,
    reference_kind = expected_kind,
    canonical_publication_allowed = FALSE,
    observed_hashes = observed_hashes
  )
}

figure7_build_generated_f <- function(reference, output_dir, config) {
  for (file in setdiff(figure7_state_required_files(), "state_pathway_provenance.tsv")) {
    figure7_copy_file(
      reference$files[[file]],
      file.path(output_dir, "tables", file)
    )
  }
  figure7_copy_file(
    reference$files[["state_pathway_provenance.tsv"]],
    file.path(output_dir, "metadata", "state_pathway_provenance.tsv")
  )
  comparison <- data.frame(
    file = names(reference$files),
    sha256 = unname(reference$observed_hashes[names(reference$files)]),
    reference_id = reference$reference_id,
    reference_kind = reference$reference_kind,
    canonical_publication_allowed = FALSE,
    stringsAsFactors = FALSE
  )
  figure7_write_tsv(
    comparison,
    file.path(
      output_dir,
      "tables",
      "state_pathway_generated_reference_comparison.tsv"
    )
  )
  plot <- figure7_panel_f_plot(reference$activity, config)
  figure7_save_panel(
    plot,
    file.path(output_dir, "figures", config$panels$filenames[["7F"]]),
    9,
    8
  )
  invisible(plot)
}
