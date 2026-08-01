module_dir <- normalizePath(Sys.getenv("FIGURE7_MODULE_DIR"), mustWork = TRUE)
repo_root <- normalizePath(file.path(module_dir, "..", "..", ".."), mustWork = TRUE)
for (file in c("common_io.R", "feature_species_policy.R", "tgi_data.R", "tgi_statistics.R", "tgi_panels.R",
               "context_panels.R", "state_pathway_panel.R", "generated_state_pathway_reference.R",
               "state_pathway_analysis.R")) {
  source(file.path(module_dir, "src", file), local = FALSE)
}
source(
  file.path(repo_root, "Code/in-vivo/SI_figures/normalized_composition.R"),
  local = FALSE
)
source(
  file.path(repo_root, "Code/in-vivo/SI_figures/shared_context_panels.R"),
  local = FALSE
)

figure7_test_inputs <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
    cellcycle <- figure7_read_cell_table(file.path(repo_root, "Data/in-vivo/figure7/processed",
      "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"), "CellCycle", config)
    noncellcycle <- figure7_read_cell_table(file.path(repo_root, "Data/in-vivo/figure7/processed",
      "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"), "NonCellCycle", config)
    endpoint_ploidy <- figure7_read_endpoint_ploidy_table(
      file.path(
        repo_root,
        "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
      ),
      config
    )
    samples <- figure7_sample_table(
      cellcycle,
      noncellcycle,
      config,
      endpoint_ploidy
    )
    cache <<- list(config = config, cellcycle = cellcycle, noncellcycle = noncellcycle,
                   endpoint_ploidy = endpoint_ploidy, samples = samples,
                   data = figure7_prepare_cellcycle(cellcycle, samples, config))
    cache
  }
})

figure7_test_reviewed_state_reference <- function() {
  config <- figure7_read_config(file.path(module_dir, "figure7_config.yaml"))
  source <- file.path(
    repo_root,
    config$state_pathways$reviewed_reference_root,
    config$state_pathways$reviewed_reference_id
  )
  parent <- tempfile("figure7_reviewed_state_fixture_")
  path <- file.path(parent, config$state_pathways$reviewed_reference_id)
  dir.create(path, recursive = TRUE)
  files <- figure7_state_required_files()
  copied <- file.copy(file.path(source, files), file.path(path, files))
  if (!all(copied)) {
    stop("Could not stage the reviewed panel-7F test fixture", call. = FALSE)
  }
  list(path = path, config = config)
}
