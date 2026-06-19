#!/usr/bin/env Rscript

resolve_in_vivo_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  candidate_files <- character(0)
  if (length(file_match) > 0) {
    candidate_files <- c(candidate_files, sub("--file=", "", file_match[1]))
  }
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) x$ofile else NA_character_
    },
    character(1)
  )
  candidate_files <- c(candidate_files, frame_files[!is.na(frame_files)])
  candidate_files <- candidate_files[!is.na(candidate_files) & nzchar(candidate_files)]
  candidate_dirs <- unique(dirname(normalizePath(candidate_files, mustWork = FALSE)))
  cwd <- normalizePath(getwd(), mustWork = FALSE)
  candidate_dirs <- unique(c(candidate_dirs, cwd, file.path(cwd, "Code", "in-vivo")))

  utils_paths <- file.path(candidate_dirs, "Utils.R")
  hit <- candidate_dirs[file.exists(utils_paths)]
  if (length(hit) > 0) return(normalizePath(hit[1], mustWork = TRUE))
  stop("Cannot locate Code/in-vivo/Utils.R from script path or working directory: ", cwd, call. = FALSE)
}

script_dir <- resolve_in_vivo_script_dir()
source(file.path(script_dir, "Utils.R"))
source(file.path(script_dir, "04c_extra_helpers.R"))
source(file.path(script_dir, "04e_all_cells_extra_helpers.R"))

required_packages <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr", "hdf5r", "Matrix", "igraph")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(tidyr)
})

set.seed(1234)

get_env_scalar <- function(name, default) {
  value <- Sys.getenv(name, unset = NA_character_)
  if (is.na(value) || !nzchar(trimws(value))) return(default)
  value
}

get_env_int <- function(name, default) {
  value <- suppressWarnings(as.integer(get_env_scalar(name, as.character(default))))
  if (!is.finite(value) || is.na(value)) return(default)
  value
}

select_scenarios <- function(specs) {
  requested <- get_env_scalar("EXTRA_SCENARIO", "")
  if (!nzchar(requested)) return(specs)
  wanted <- trimws(strsplit(requested, ",", fixed = TRUE)[[1]])
  selected <- specs[specs$scenario_id %in% wanted, , drop = FALSE]
  if (nrow(selected) == 0) {
    stop("EXTRA_SCENARIO did not match available scenarios: ", paste(specs$scenario_id, collapse = ", "), call. = FALSE)
  }
  selected
}

run_paga_scenario <- function(spec, results_root, output_base, deg_round2_root, top_n, max_cascade_genes, max_neighbor_edges) {
  trajectory_root <- if ("trajectory_root" %in% names(spec) && nzchar(as.character(spec$trajectory_root[1]))) {
    as.character(spec$trajectory_root[1])
  } else {
    file.path(results_root, spec$trajectory_subdir[1])
  }
  output_root <- file.path(output_base, spec$scenario_id[1])

  message("Input roots:")
  message("  scenario: ", spec$scenario_id)
  message("  trajectory_branch: ", spec$trajectory_branch)
  if ("root_cluster" %in% names(spec)) message("  root_cluster: ", spec$root_cluster)
  if ("end_cluster" %in% names(spec)) message("  end_cluster: ", spec$end_cluster)
  message("  trajectory_root: ", trajectory_root)
  message("  deg_round2_root: ", deg_round2_root)
  message("  output_root: ", output_root)

  if (!dir.exists(trajectory_root)) stop("Missing trajectory root: ", trajectory_root, call. = FALSE)
  if (!dir.exists(file.path(trajectory_root, "04_paga_groups"))) {
    stop("Missing PAGA result root: ", file.path(trajectory_root, "04_paga_groups"), call. = FALSE)
  }
  if (!dir.exists(deg_round2_root)) warning("Missing 03c DEG round2 root: ", deg_round2_root, call. = FALSE)

  paga_extra <- run_paga_extra_analysis_all_cells_derived(
    trajectory_root = trajectory_root,
    deg_round2_root = deg_round2_root,
    output_root = output_root,
    analysis_label = paste("04e2 PAGA All_cells-derived", spec$scenario_id),
    top_n = top_n,
    max_cascade_genes = max_cascade_genes,
    max_neighbor_edges = max_neighbor_edges
  )

  status <- if (is.data.frame(paga_extra$status)) paga_extra$status else data.frame()
  manifest <- if (is.data.frame(paga_extra$manifest)) paga_extra$manifest else data.frame()
  run_summary_lines <- c(
    "04e2 extra PAGA All_cells-derived analysis completed.",
    paste0("Scenario: ", spec$scenario_id),
    paste0("Trajectory branch: ", spec$trajectory_branch),
    paste0("Trajectory root: ", normalizePath(trajectory_root, mustWork = FALSE)),
    paste0("Output root: ", normalizePath(output_root, mustWork = FALSE)),
    paste0("All_cells-derived PAGA groups processed: ", nrow(status)),
    paste0("All_cells-derived PAGA groups completed: ", if (nrow(status) > 0) sum(status$status == "completed", na.rm = TRUE) else 0L),
    paste0("All_cells-derived PAGA groups partial: ", if (nrow(status) > 0) sum(status$status == "partial", na.rm = TRUE) else 0L),
    paste0("All_cells-derived PAGA groups failed: ", if (nrow(status) > 0) sum(status$status == "failed", na.rm = TRUE) else 0L),
    "",
    "Main inputs:",
    "  04_paga_groups/All_cells/paga_result.h5ad",
    "  04_paga_groups/All_cells/paga_edges.csv",
    "  04_paga_groups/All_cells/paga_cell_metrics.csv",
    "  03c_DEGs_round2 cluster DEG tables",
    "",
    "Key outputs:",
    "  paga_run_manifest.csv",
    "  paga_extra_status.csv",
    "  */paga_umap_edges_on_data.pdf",
    "  */paga_umap_neighbor_edges_on_data.pdf",
    "  */paga_dpt_pseudotime_umap.pdf",
    "  */paga_top10_deg_ordered_cluster_heatmap.pdf",
    "  */paga_top_deg_expression_cascades.pdf"
  )
  writeLines(run_summary_lines, file.path(output_root, "run_summary.txt"))
  message(paste(run_summary_lines, collapse = "\n"))

  data.frame(
    scenario_id = spec$scenario_id,
    trajectory_branch = spec$trajectory_branch,
    trajectory_root = trajectory_root,
    output_root = output_root,
    paga_manifest_rows = nrow(manifest),
    paga_status_rows = nrow(status),
    paga_completed = if (nrow(status) > 0) sum(status$status == "completed", na.rm = TRUE) else 0L,
    paga_partial = if (nrow(status) > 0) sum(status$status == "partial", na.rm = TRUE) else 0L,
    paga_failed = if (nrow(status) > 0) sum(status$status == "failed", na.rm = TRUE) else 0L,
    stringsAsFactors = FALSE
  )
}

config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_env_scalar("EXTRA_RESULTS_ROOT", get_results_root(config))
deg_round2_root <- get_env_scalar("EXTRA_DEG_ROUND2_ROOT", file.path(results_root, "03c_DEGs_round2"))
output_base <- get_env_scalar("EXTRA_04E2_PAGA_OUTPUT_ROOT", file.path(results_root, "04e2_extra_PAGA"))
trajectory_result_root <- get_env_scalar("EXTRA_TRAJECTORY_ROOT", file.path(results_root, "04_trajectory"))
top_n <- get_env_int("EXTRA_PAGA_TOP_N", 10L)
max_cascade_genes <- get_env_int("EXTRA_PAGA_MAX_CASCADE_GENES", 30L)
max_neighbor_edges <- get_env_int("EXTRA_PAGA_MAX_NEIGHBOR_EDGES", 50000L)

scenario_specs <- extra_discover_paga_scenarios(trajectory_result_root, config)
if (nrow(scenario_specs) == 0) {
  stop("No PAGA scenarios found or configured under: ", trajectory_result_root, call. = FALSE)
}
scenario_specs <- select_scenarios(scenario_specs)

summary_rows <- lapply(seq_len(nrow(scenario_specs)), function(i) {
  run_paga_scenario(scenario_specs[i, , drop = FALSE], results_root, output_base, deg_round2_root, top_n, max_cascade_genes, max_neighbor_edges)
})
summary_df <- dplyr::bind_rows(summary_rows)
.ensure_dir(output_base)
write_table_csv(summary_df, file.path(output_base, "04e2_PAGA_scenario_summary.csv"))
message("04e2 extra PAGA All_cells-derived scenarios completed: ", paste(summary_df$scenario_id, collapse = ", "))
