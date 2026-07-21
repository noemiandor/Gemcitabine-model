#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v3.R"
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)

source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v3_util.R"))

defaults <- list(
  cell_metadata = file.path(repo_root, "Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
  noncell_metadata = file.path(repo_root, "Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
  sample_info = file.path(repo_root, "Data/in-vivo/sample_info.xlsx"),
  seurat_rds = "",
  config = file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v3_config.yaml"),
  output_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs_v3",
  results_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results",
  assay = "RNA",
  counts_layer = "counts",
  overwrite = FALSE,
  parallel = TRUE,
  workers = 4L,
  rotation_test_workers = 1L,
  permutation_workers = 1L,
  simulation_workers = 1L,
  crossfit_workers = 1L,
  multidim_workers = 1L,
  trajectory_workers = 1L,
  seed = 73157L,
  max_cells = 0L,
  max_genes = 0L,
  run_gsea = TRUE,
  run_leave_one_out = TRUE,
  run_sensitivities = TRUE,
  run_continuous = TRUE,
  run_simulation = TRUE,
  run_crossfit_04i = TRUE,
  run_exact_permutation = TRUE,
  run_control_reference = TRUE,
  run_covariance_whitened = TRUE,
  run_multidimensional = TRUE,
  run_trajectory_global = TRUE,
  gsea_min_size = 15L,
  gsea_max_size = 500L,
  gsea_nperm_simple = 10000L
)

args <- ptp_parse_args(commandArgs(trailingOnly = TRUE), defaults)
ptpv3_run_workflow(args, repo_root = repo_root, script_dir = script_dir)
