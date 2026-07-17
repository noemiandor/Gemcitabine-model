module_dir <- normalizePath(Sys.getenv("FIGURE7_MODULE_DIR"), mustWork = TRUE)
repo_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
for (file in c("common_io.R", "tgi_data.R", "tgi_statistics.R", "tgi_panels.R",
               "state_pathway_panel.R", "state_pathway_analysis.R")) {
  source(file.path(module_dir, "src", file), local = FALSE)
}

figure7_test_inputs <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
    cellcycle <- figure7_read_cell_table(file.path(repo_root, "Data/in-vivo/figure7/processed",
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"), "CellCycle")
    noncellcycle <- figure7_read_cell_table(file.path(repo_root, "Data/in-vivo/figure7/processed",
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"), "NonCellCycle")
    samples <- figure7_sample_table(cellcycle, noncellcycle, config)
    cache <<- list(config = config, cellcycle = cellcycle, noncellcycle = noncellcycle,
                   samples = samples, data = figure7_prepare_cellcycle(cellcycle, samples))
    cache
  }
})

figure7_test_state_reference <- function() {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  parent <- tempfile("figure7_state_fixture_"); dir.create(parent)
  path <- file.path(parent, config$state_pathways$reference_id); dir.create(path)
  collections <- as.character(unlist(config$state_pathways$collections))
  collection_labels <- c("Hallmark", "Reactome", "GO biological process")
  grid <- seq(0, 1, length.out = 501)
  selected_rows <- list(); activity_rows <- list(); leading_rows <- list(); complete_rows <- list()
  display <- 0L
  for (i in seq_along(collections)) {
    collection <- collections[[i]]
    for (rank in 1:4) {
      for (direction in c("positive", "negative")) {
        display <- display + 1L
        pathway <- paste(collection, direction, rank, sep = "_")
        nes <- if (direction == "positive") 5 - rank else -(5 - rank)
        padj <- rank / 100
        meta <- data.frame(
          collection_id = collection, collection_label = collection_labels[[i]], collection_display_order = i,
          pathway_id = pathway, pathway_label = pathway, pathway_display_order = display,
          NES = nes, padj = padj, selected_direction = direction,
          selected_rank_within_direction = rank, stringsAsFactors = FALSE
        )
        selected_rows[[length(selected_rows) + 1L]] <- meta
        leading_rows[[length(leading_rows) + 1L]] <- data.frame(
          collection_id = collection, pathway_id = pathway, leading_edge_gene_id = paste0("GENE", display))
        activity_rows[[length(activity_rows) + 1L]] <- cbind(meta[, setdiff(names(meta), c("NES", "padj"))],
          pseudotime = grid, standardized_activity = sin(2 * pi * grid) + display / 100)
        complete_rows[[length(complete_rows) + 1L]] <- meta[, c("collection_id", "pathway_id", "NES", "padj")]
      }
    }
    complete_rows[[length(complete_rows) + 1L]] <- data.frame(
      collection_id = collection, pathway_id = paste0(collection, "_extra"), NES = 0.1, padj = 0.9)
  }
  selected <- do.call(rbind, selected_rows); activity <- do.call(rbind, activity_rows)
  figure7_write_tsv(activity, file.path(path, "panel_7F_pathway_activity_plot_data.tsv"))
  figure7_write_tsv(selected, file.path(path, "panel_7F_selected_pathway_gsea.tsv"))
  figure7_write_tsv(do.call(rbind, leading_rows), file.path(path, "panel_7F_leading_edge_genes.tsv"))
  figure7_write_tsv(data.frame(gene_id = c("G1", "G2"), moderated_t = c(2, -2)),
                    file.path(path, "state_pathway_gene_ranking_complete.tsv"))
  figure7_write_tsv(do.call(rbind, complete_rows), file.path(path, "state_pathway_gsea_complete.tsv"))
  figure7_write_tsv(data.frame(sample_id = "mouse1", pseudotime_bin = 1, n_cells = 20),
                    file.path(path, "state_pathway_sample_bin_coverage.tsv"))
  figure7_write_tsv(data.frame(metric = "design_rank", value = 8),
                    file.path(path, "state_pathway_design_qc.tsv"))
  activity_hash <- figure7_sha256(file.path(path, "panel_7F_pathway_activity_plot_data.tsv"))
  provenance_keys <- c(
    "full_analysis_run_dir", "report_identifier", "code_revision_04i", "seurat_rds_sha256",
    "cellcycle_metadata_sha256", "noncellcycle_metadata_sha256", "interval_config_sha256",
    "assay", "counts_layer", "etp_method", "etp_threshold", "spline_df", "pseudotime_bins",
    "minimum_cells_per_sample_bin", "grid_size", "seed", "gene_set_source", "gene_set_release",
    "gene_set_species", "gene_set_collections", "gene_set_min_size", "gene_set_max_size",
    "expression_filter", "normalization", "observation_model", "mouse_block", "nuisance_terms",
    "treatment_by_pseudotime_interaction", "empirical_bayes", "contrast", "gsea_rank_statistic",
    "pathway_activity", "feature_species_policy", "gsea_ranking_rule", "pathway_selection_rule",
    "activity_table_sha256"
  )
  values <- rep("fixture", length(provenance_keys)); names(values) <- provenance_keys
  values[c("assay", "counts_layer", "expression_filter", "normalization", "observation_model", "mouse_block",
           "empirical_bayes", "contrast", "gsea_rank_statistic", "pathway_activity")] <- vapply(
    c("assay", "counts_layer", "expression_filter", "normalization", "observation_model", "mouse_block",
      "empirical_bayes", "contrast", "gsea_rank_statistic", "pathway_activity"),
    function(key) as.character(config$state_pathways[[key]]), character(1L))
  values["etp_method"] <- as.character(config$etp$method); values["etp_threshold"] <- as.character(config$etp$threshold)
  values[c("spline_df", "pseudotime_bins", "minimum_cells_per_sample_bin", "grid_size")] <- vapply(
    c("spline_df", "pseudotime_bins", "minimum_cells_per_sample_bin", "grid_size"),
    function(key) as.character(config$state_pathways[[key]]), character(1L))
  values["seed"] <- as.character(config$statistics$seed)
  values["nuisance_terms"] <- paste(unlist(config$state_pathways$nuisance_terms), collapse = ",")
  values["gene_set_collections"] <- paste(unlist(config$state_pathways$collections), collapse = ",")
  values["treatment_by_pseudotime_interaction"] <- as.character(config$state_pathways$treatment_by_pseudotime_interaction)
  values["pathway_selection_rule"] <- as.character(config$state_pathways$pathway_selector)
  values["activity_table_sha256"] <- activity_hash
  figure7_write_tsv(data.frame(key = provenance_keys, value = values),
                    file.path(path, "state_pathway_provenance.tsv"))
  for (file in figure7_state_required_files()) {
    config$state_pathways$expected_files[[file]] <- figure7_sha256(file.path(path, file))
  }
  list(path = path, config = config)
}
