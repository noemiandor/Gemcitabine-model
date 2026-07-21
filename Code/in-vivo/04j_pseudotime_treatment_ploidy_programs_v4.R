#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v4.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)

source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_v4_util.R"))

defaults <- list(
  cell_metadata = file.path(
    repo_root,
    "Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ),
  noncell_metadata = file.path(
    repo_root,
    "Data/in-vivo/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
  ),
  sample_info = file.path(repo_root, "Data/in-vivo/sample_info.xlsx"),
  seurat_rds = "",
  config = file.path(
    script_dir,
    "04j_pseudotime_treatment_ploidy_programs_v4_config.yaml"
  ),
  output_root = paste0(
    "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/",
    "04j_pseudotime_treatment_ploidy_programs_v4"
  ),
  results_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results",
  assay = "RNA",
  counts_layer = "counts",
  overwrite = FALSE,
  parallel = TRUE,
  workers = 4L,
  permutation_workers = 1L,
  simulation_workers = 1L,
  crossfit_workers = 1L,
  adaptive_workers = 1L,
  multiscale_workers = 1L,
  trajectory_workers = 1L,
  dose_workers = 1L,
  coordinate_workers = 1L,
  rotation_test_workers = 1L,
  model_workers = 1L,
  within_model_workers = 1L,
  seed = 74157L,
  max_cells = 0L,
  max_genes = 0L,
  smoke = FALSE,
  smoke_null_replicates = 8L,
  smoke_power_replicates = 4L,
  smoke_permutations = 40L,
  smoke_dose_assignments = 80L
)

args <- ptp_parse_args(commandArgs(trailingOnly = TRUE), defaults)
ptpv4_run_workflow(args, repo_root = repo_root, script_dir = script_dir)
