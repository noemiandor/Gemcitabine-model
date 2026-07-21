#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg) > 0L) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs.R"
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)

source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_util.R"))

defaults <- list(
  cell_metadata = file.path(repo_root, "Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"),
  sample_info = file.path(repo_root, "Data/in-vivo/sample_info.xlsx"),
  seurat_rds = "",
  config = file.path(script_dir, "04j_pseudotime_treatment_ploidy_programs_config.yaml"),
  output_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs",
  results_root = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results",
  assay = "RNA",
  counts_layer = "counts",
  overwrite = FALSE,
  parallel = TRUE,
  workers = 4L,
  seed = 1L,
  max_cells = 0L,
  max_genes = 0L,
  run_gsea = TRUE,
  run_leave_one_out = TRUE,
  run_sensitivities = TRUE,
  run_continuous = TRUE,
  gsea_min_size = 15L,
  gsea_max_size = 500L,
  gsea_nperm_simple = 10000L
)

args <- ptp_parse_args(commandArgs(trailingOnly = TRUE), defaults)
ptp_run_workflow(args, repo_root = repo_root, script_dir = script_dir)
