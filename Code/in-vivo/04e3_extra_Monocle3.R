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

required_packages <- c("dplyr", "ggplot2", "readr", "tibble", "tidyr")
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

setup_monocle3_outputs <- function(output_root) {
  method_scvelo <<- "scVelo"
  method_monocle3 <<- "Monocle3"
  method_dir_scvelo <<- "scVelo"
  method_dir_monocle3 <<- "Monocle3"
  out_inputs <<- .ensure_dir(file.path(output_root, "00_inputs"))
  out_tables <<- .ensure_dir(file.path(output_root, "01_tables"))
  out_stats <<- .ensure_dir(file.path(output_root, "02_stats"))
  out_plots <<- .ensure_dir(file.path(output_root, "03_plots"))
  out_summary <<- .ensure_dir(file.path(output_root, "04_question_summary"))
  out_separate <<- .ensure_dir(file.path(output_root, "05_method_separate"))
  out_monocle3 <<- .ensure_dir(file.path(out_separate, method_dir_monocle3))
  out_biology <<- .ensure_dir(file.path(output_root, "06_deg_gsea_pseudotime"))
  out_umap <<- .ensure_dir(file.path(output_root, "07_umap_overlays"))
  out_umap_monocle3 <<- .ensure_dir(file.path(out_umap, method_dir_monocle3))
  out_flow <<- .ensure_dir(file.path(output_root, "08_flow_overlays"))
  out_flow_monocle3 <<- .ensure_dir(file.path(out_flow, method_dir_monocle3))
  out_tn_ploidy <<- .ensure_dir(file.path(output_root, "09_tn_ploidy_cluster_response"))
  out_tn_ploidy_monocle3 <<- .ensure_dir(file.path(out_tn_ploidy, method_dir_monocle3))
  out_cluster_timing <<- .ensure_dir(file.path(output_root, "10_cluster_timing"))
  out_cluster_timing_monocle3 <<- .ensure_dir(file.path(out_cluster_timing, method_dir_monocle3))
  invisible(output_root)
}

ensure_trajectory_split_cols <- function(df) {
  ensure_columns(df, c(
    "TrajectoryComparisonID", "TrajectoryScope", "TrajectoryComparison",
    "TrajectoryIdent1", "TrajectoryIdent2", "TrajectoryGroupCol",
    "TrajectoryCluster", "TrajectorySubsetGroup", "TrajectoryDose",
    "TrajectoryComparisonGroup", "TrajectoryOriginalGroup", "TrajectoryUniverseKey"
  ))
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

run_monocle3_scenario <- function(spec, output_base, deg_root, gsea_root) {
  pseudotrajectory_root <- if ("trajectory_root" %in% names(spec) && nzchar(as.character(spec$trajectory_root[1]))) {
    as.character(spec$trajectory_root[1])
  } else {
    as.character(spec$branch_root[1])
  }
  output_root <- file.path(output_base, spec$scenario_id[1])
  setup_monocle3_outputs(output_root)

  message("Input roots:")
  message("  scenario: ", spec$scenario_id)
  message("  trajectory_branch: ", spec$trajectory_branch)
  if ("root_cluster" %in% names(spec)) message("  root_cluster: ", spec$root_cluster)
  message("  pseudotrajectory_root: ", pseudotrajectory_root)
  message("  deg_root: ", deg_root)
  message("  gsea_root: ", gsea_root)
  message("  output_root: ", output_root)

  if (!dir.exists(pseudotrajectory_root)) stop("Missing pseudotrajectory root: ", pseudotrajectory_root, call. = FALSE)
  if (!dir.exists(gsea_root)) warning("Missing 03b GSEA root: ", gsea_root, call. = FALSE)

  comparison_summary_file <- file.path(deg_root, "00_summary", "DEG_comparison_summary.csv")
  if (file.exists(comparison_summary_file)) {
    comparison_summary <- readr::read_csv(comparison_summary_file, show_col_types = FALSE, progress = FALSE)
    write_table_csv(comparison_summary, file.path(out_inputs, "DEG_comparison_summary_used.csv"))
  }

  cluster_annotations <- if (dir.exists(gsea_root)) read_cluster_annotations(gsea_root) else data.frame()
  gsea_inputs <- if (dir.exists(gsea_root)) read_gsea_inputs(gsea_root) else data.frame()
  deg_marker_summary <- if (dir.exists(deg_root)) read_deg_marker_summaries(deg_root) else data.frame()
  write_table_csv(cluster_annotations, file.path(out_inputs, "cluster_annotation_summary_used.csv"))
  write_table_csv(gsea_inputs, file.path(out_inputs, "hallmark_gsea_tables_used.csv"))
  write_table_csv(deg_marker_summary, file.path(out_inputs, "deg_marker_top_gene_summary_used.csv"))

  raw_monocle3_cells <- extra_all_cells_read_monocle3_metrics(pseudotrajectory_root)
  if (nrow(raw_monocle3_cells) == 0) stop("No All_cells Monocle3 metric rows found for ", spec$scenario_id, call. = FALSE)
  raw_monocle3_cells$extra_scenario <- spec$scenario_id
  raw_monocle3_cells$extra_trajectory_branch <- spec$trajectory_branch
  write_table_csv(
    extra_all_cells_tn_split_summary(raw_monocle3_cells),
    file.path(out_inputs, "all_cells_tn_split_summary.csv")
  )
  monocle3_cells <- extra_all_cells_expand_metrics(raw_monocle3_cells)
  write_table_csv(monocle3_cells, file.path(out_tables, "Monocle3_cell_metrics_long.csv"))

  monocle3_suite <- run_comparison_suite(monocle3_cells, pseudo_metric_cols, method_monocle3)
  stat_tests <- ensure_trajectory_split_cols(add_adjusted_p(monocle3_suite$tests))
  group_summaries <- ensure_trajectory_split_cols(monocle3_suite$summaries)
  write_table_csv(stat_tests, file.path(out_stats, "Monocle3_stat_tests.csv"))
  write_table_csv(group_summaries, file.path(out_stats, "Monocle3_group_summaries.csv"))

  question_summary <- summarize_questions(stat_tests)
  write_table_csv(question_summary, file.path(out_summary, "question_level_summary.csv"))
  time_metric_tests <- stat_tests %>%
    dplyr::filter(.data$metric == "pseudotime")
  write_table_csv(time_metric_tests, file.path(out_summary, "time_metric_tests_for_interpretation.csv"))

  monocle3_separate <- write_method_separate_outputs_all_cells_derived(
    monocle3_cells,
    metric_col = "pseudotime",
    method_label = method_monocle3,
    out_dir = out_monocle3,
    plot_prefix = "Monocle3",
    out_umap_dir = out_umap_monocle3,
    cluster_annotations = cluster_annotations
  )
  method_answer_summary <- make_method_answer_summary_all_cells_derived(monocle3_separate, method_monocle3)
  write_table_csv(method_answer_summary, file.path(out_summary, "method_separate_answer_summary.csv"))

  monocle3_tn_ploidy <- write_tn_ploidy_cluster_outputs(
    monocle3_separate$desc,
    metric_col = "pseudotime",
    method_label = method_monocle3,
    out_dir = out_tn_ploidy_monocle3,
    plot_prefix = "Monocle3",
    cluster_annotations = cluster_annotations
  )
  write_table_csv(monocle3_tn_ploidy$key_conclusions, file.path(out_tn_ploidy, "Monocle3_tn_ploidy_cluster_key_conclusions.csv"))
  write_table_csv(monocle3_tn_ploidy$tumor_ploidy_dose$ploidy_within_dose_sample_tests, file.path(out_tn_ploidy, "Monocle3_tumor_Q1_ploidy_within_dose_sample_tests.csv"))
  write_table_csv(monocle3_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_sample_tests, file.path(out_tn_ploidy, "Monocle3_tumor_Q2_dose_within_ploidy_sample_tests.csv"))
  write_table_csv(monocle3_tn_ploidy$tumor_ploidy_dose$dose_within_ploidy_pairwise_sample_tests, file.path(out_tn_ploidy, "Monocle3_tumor_Q2_dose_pairwise_within_ploidy_sample_tests.csv"))
  write_table_csv(monocle3_tn_ploidy$tumor_ploidy_dose$feature_summary, file.path(out_tn_ploidy, "Monocle3_tumor_Q3_ploidy_dose_feature_summary.csv"))

  monocle3_group_timing <- run_group_cluster_timing_outputs(
    monocle3_cells,
    out_dir = out_cluster_timing_monocle3,
    method_label = method_monocle3,
    time_col = "pseudotime",
    source_context_filter = c("monocle3_groups"),
    time_label = "pseudotime",
    x_label = "Monocle3 pseudotime"
  )

  pathway_trajectory <- join_trajectory_with_gsea(stat_tests, gsea_inputs, deg_marker_summary)
  pathway_timing_summary <- summarize_pathway_timing(pathway_trajectory)
  hypoxia_timing <- if (nrow(pathway_trajectory) > 0 && all(c("hypoxia_pathway", "pathway") %in% names(pathway_trajectory))) {
    pathway_trajectory %>%
      dplyr::filter(.data$hypoxia_pathway, !is.na(.data$pathway))
  } else {
    data.frame()
  }
  write_table_csv(pathway_trajectory, file.path(out_biology, "Monocle3_pseudotime_gsea_deg_joined.csv"))
  write_table_csv(pathway_timing_summary, file.path(out_biology, "Monocle3_pathway_timing_summary.csv"))
  write_table_csv(hypoxia_timing, file.path(out_biology, "Monocle3_hypoxia_pseudotime_timing.csv"))

  plot_scope_summary(stat_tests)
  plot_deg_heatmap(stat_tests)
  plot_cluster_vs_rest(stat_tests)
  plot_global_effects(stat_tests)
  plot_within_cluster_ploidy(stat_tests)
  plot_omnibus_cluster(stat_tests)
  plot_merged_dose_ploidy(stat_tests)
  plot_gsea_trajectory_association(pathway_trajectory)
  plot_hypoxia_timing(hypoxia_timing)
  plot_pathway_timing_summary(pathway_timing_summary)
  plot_monocle3_graph_overlays_all_cells(pseudotrajectory_root, out_flow_monocle3, out_umap_monocle3)
  if (nrow(gsea_inputs) > 0) {
    plot_cluster_hypoxia_umap_all_cells_derived(
      monocle3_separate,
      gsea_inputs,
      method_label = method_monocle3,
      metric_col = "pseudotime",
      out_dir = out_monocle3,
      plot_prefix = "Monocle3",
      out_umap_dir = out_umap_monocle3
    )
  }

  run_summary_lines <- c(
    "04e3 extra Monocle3 All_cells-derived analysis completed.",
    paste0("Scenario: ", spec$scenario_id),
    paste0("Trajectory branch: ", spec$trajectory_branch),
    paste0("Pseudotrajectory root: ", normalizePath(pseudotrajectory_root, mustWork = FALSE)),
    paste0("Output root: ", normalizePath(output_root, mustWork = FALSE)),
    paste0("All_cells Monocle3 input rows: ", nrow(raw_monocle3_cells)),
    paste0("All_cells-derived Monocle3 analysis rows: ", nrow(monocle3_cells)),
    paste0("Stat test rows: ", nrow(stat_tests)),
    paste0("Group summary rows: ", nrow(group_summaries)),
    paste0("GSEA pseudotime rows: ", nrow(pathway_trajectory)),
    paste0("TN/Ploidy Monocle3 key conclusion rows: ", nrow(monocle3_tn_ploidy$key_conclusions)),
    paste0("Monocle3 group cluster timing UMAP rows: ", nrow(monocle3_group_timing$status)),
    "",
    "Main time metric:",
    "  Monocle3: pseudotime",
    "",
    "Key outputs:",
    "  01_tables/Monocle3_cell_metrics_long.csv",
    "  02_stats/Monocle3_stat_tests.csv",
    "  02_stats/Monocle3_group_summaries.csv",
    "  05_method_separate/Monocle3/*",
    "  07_umap_overlays/Monocle3/*",
    "  08_flow_overlays/Monocle3/*",
    "  09_tn_ploidy_cluster_response/Monocle3/*",
    "  10_cluster_timing/Monocle3/*"
  )
  writeLines(run_summary_lines, file.path(output_root, "run_summary.txt"))
  message(paste(run_summary_lines, collapse = "\n"))

  data.frame(
    scenario_id = spec$scenario_id,
    trajectory_branch = spec$trajectory_branch,
    root_cluster = if ("root_cluster" %in% names(spec)) spec$root_cluster else NA_character_,
    pseudotrajectory_root = pseudotrajectory_root,
    output_root = output_root,
    monocle3_cell_rows = nrow(monocle3_cells),
    stat_test_rows = nrow(stat_tests),
    stringsAsFactors = FALSE
  )
}

config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_env_scalar("EXTRA_RESULTS_ROOT", get_results_root(config))
pseudotrajectory_result_root <- get_env_scalar("EXTRA_PSEUDOTRAJECTORY_ROOT", file.path(results_root, "04a_psudo_rajectory"))
deg_root <- get_env_scalar("EXTRA_DEG_ROOT", file.path(results_root, "03a_DEGs"))
gsea_root <- get_env_scalar("EXTRA_GSEA_ROOT", file.path(results_root, "03b_cluster_annotation_and_GSEA"))
output_base <- get_env_scalar("EXTRA_04E3_MONOCLE3_OUTPUT_ROOT", file.path(results_root, "04e3_extra_Monocle3"))

scenario_specs <- extra_discover_monocle3_scenarios(pseudotrajectory_result_root, config)
if (nrow(scenario_specs) == 0) {
  stop("No Monocle3 scenarios found or configured under: ", pseudotrajectory_result_root, call. = FALSE)
}
scenario_specs <- select_scenarios(scenario_specs)

summary_rows <- lapply(seq_len(nrow(scenario_specs)), function(i) {
  run_monocle3_scenario(scenario_specs[i, , drop = FALSE], output_base, deg_root, gsea_root)
})
summary_df <- dplyr::bind_rows(summary_rows)
.ensure_dir(output_base)
write_table_csv(summary_df, file.path(output_base, "04e3_Monocle3_scenario_summary.csv"))
message("04e3 extra Monocle3 All_cells-derived scenarios completed: ", paste(summary_df$scenario_id, collapse = ", "))
