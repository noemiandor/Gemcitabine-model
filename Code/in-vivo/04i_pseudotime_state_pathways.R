#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/04i_pseudotime_state_pathways.R"
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)

source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04i_pseudotime_state_pathways_util.R"))

defaults <- list(
  cell_metadata = file.path(repo_root, "Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
  noncell_metadata = file.path(repo_root, "Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
  seurat_rds = "",
  config = file.path(script_dir, "04i_pseudotime_state_pathways_config.yaml"),
  output_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04i_pseudotime_state_pathways",
  snapshot_root = file.path(repo_root, "Figs/pseudotime_state_pathways"),
  results_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results",
  assay = "RNA",
  counts_layer = "counts",
  n_pseudotime_bins = 20L,
  min_cells_per_sample_bin = 10L,
  spline_df = 5L,
  gene_set_collections = "H,C2:CP:REACTOME,C5:GO:BP",
  seed = 1L,
  workers = 4L,
  overwrite = FALSE,
  min_match_rate = 0.99,
  gsea_min_size = 15L,
  gsea_max_size = 500L,
  gsea_nperm_simple = 10000L,
  grid_size = 501L,
  max_cells = 0L,
  max_genes = 0L,
  run_sensitivity = TRUE,
  run_leave_one_out = TRUE,
  covariate_mode = "all",
  etp_threshold = "all",
  include_full_etp_models = FALSE,
  make_snapshot = TRUE
)

args <- pst_parse_args(commandArgs(trailingOnly = TRUE), defaults)
pst_run_workflow(args, repo_root = repo_root, script_dir = script_dir)
