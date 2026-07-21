#!/usr/bin/env Rscript

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_v5.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
repo_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
module_dir <- file.path(script_dir, "04j_v5_non_score")

source(file.path(script_dir, "Utils.R"))
source(file.path(module_dir, "common.R"))
source(file.path(module_dir, "01_exact_ranked_enrichment.R"))
source(file.path(module_dir, "02_state_decomposition.R"))
source(file.path(module_dir, "03_bayesian_gene_vector.R"))
source(file.path(module_dir, "04_exact_global_kernel.R"))
source(file.path(module_dir, "05_gene_wise_trajectory.R"))
source(file.path(module_dir, "06_simulation_and_synthesis.R"))
source(file.path(module_dir, "07_report_and_audit.R"))

defaults <- list(
  config = file.path(
    script_dir,
    "04j_pseudotime_treatment_ploidy_programs_v5_config.yaml"
  ),
  output_root = paste0(
    "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/",
    "04j_pseudotime_treatment_ploidy_programs_v5_non_score"
  ),
  overwrite = FALSE,
  workers = 4L,
  smoke = FALSE,
  smoke_assignments = 40L,
  smoke_memberships = 8L,
  smoke_families = 3L,
  smoke_genes = 512L,
  seed = 75157L
)

args <- ptp_parse_args(commandArgs(trailingOnly = TRUE), defaults)
ptpv5_run_workflow(
  args = args,
  repo_root = repo_root,
  script_dir = script_dir,
  module_dir = module_dir
)
