#!/usr/bin/env Rscript

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1L]]) else "Code/in-vivo/figure7/run_figure7.R"
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = FALSE)
for (file in c("common_io.R", "tgi_data.R", "tgi_statistics.R", "tgi_panels.R",
               "state_pathway_panel.R", "state_pathway_analysis.R")) {
  sys.source(file.path(script_dir, "src", file), envir = .GlobalEnv)
}

if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("ggrepel", quietly = TRUE)) {
  figure7_stop("R packages 'ggplot2' and 'ggrepel' are required")
}

args <- figure7_parse_args(commandArgs(trailingOnly = TRUE))
mode <- figure7_arg(args, "mode", "standard")
allowed_modes <- c("standard", "full-analysis", "render-only")
if (!mode %in% allowed_modes) figure7_stop("Unknown Figure 7 mode: ", mode)
config_path <- normalizePath(figure7_arg(args, "config", file.path(script_dir, "figure7_config.yaml")), mustWork = FALSE)
config <- figure7_read_config(config_path)
output_dir <- normalizePath(figure7_arg(args, "output-dir", required = TRUE), mustWork = FALSE)

write_metadata <- function(output_dir, mode, config, config_path) {
  run_config <- data.frame(
    key = c("module", "mode", "tgi_outcome", "tgi_day", "tgi_measure",
            "matched_control_summary", "matched_control_group", "etp_method", "etp_threshold",
            "state_pathway_reference_id", "state_interval_start", "state_interval_end", "config_sha256"),
    value = c("in_vivo_figure7", mode, "day", "17", "TGI_percent_Day_17", "mean", "initial_ploidy",
              config$etp$method, as.character(config$etp$threshold), config$state_pathways$reference_id,
              as.character(config$state_pathways$accumulated_interval$start),
              as.character(config$state_pathways$accumulated_interval$end), figure7_sha256(config_path)),
    stringsAsFactors = FALSE
  )
  contract <- data.frame(
    panel_id = names(config$panels$filenames),
    filename = unname(unlist(config$panels$filenames)),
    tgi_outcome = "day", tgi_day = 17L, tgi_measure = "TGI_percent_Day_17",
    matched_control_summary = "mean", matched_control_group = "initial_ploidy",
    stringsAsFactors = FALSE
  )
  figure7_write_tsv(run_config, file.path(output_dir, "metadata", "run_config.tsv"))
  figure7_write_tsv(contract, file.path(output_dir, "metadata", "panel_contract.tsv"))
  utils::capture.output(sessionInfo(), file = file.path(output_dir, "metadata", "session_info.txt"))
}

render_from_run <- function(source_dir, output_dir, config, config_path) {
  source_dir <- normalizePath(source_dir, mustWork = TRUE)
  if (identical(source_dir, normalizePath(output_dir, mustWork = FALSE))) {
    figure7_stop("render-only output must differ from the immutable source run")
  }
  source_run_config <- figure7_read_tsv(file.path(source_dir, "metadata", "run_config.tsv"), c("key", "value"))
  source_contract <- figure7_read_tsv(file.path(source_dir, "metadata", "panel_contract.tsv"),
                                      c("panel_id", "filename", "tgi_outcome", "tgi_day", "tgi_measure",
                                        "matched_control_summary", "matched_control_group"))
  config_hash <- source_run_config$value[match("config_sha256", source_run_config$key)]
  if (is.na(config_hash) || !identical(config_hash, figure7_sha256(config_path))) {
    figure7_stop("render-only source run was produced with a different Figure 7 config")
  }
  expected_contract <- data.frame(panel_id = names(config$panels$filenames),
                                  filename = unname(unlist(config$panels$filenames)), stringsAsFactors = FALSE)
  observed_contract <- source_contract[, c("panel_id", "filename")]
  rownames(observed_contract) <- NULL
  if (!identical(observed_contract, expected_contract)) {
    figure7_stop("render-only source panel contract disagrees with the current six-panel contract")
  }
  if (any(source_contract$tgi_outcome != "day") || any(figure7_numeric(source_contract$tgi_day) != 17) ||
      any(source_contract$tgi_measure != "TGI_percent_Day_17") ||
      any(source_contract$matched_control_summary != "mean") ||
      any(source_contract$matched_control_group != "initial_ploidy")) {
    figure7_stop("render-only source panel contract has incompatible TGI metadata")
  }
  table_path <- function(name) file.path(source_dir, "tables", name)
  a <- figure7_read_tsv(table_path("panel_7A_plot_data.tsv"),
    c("sample_id", "initial_ploidy", "dose", "day", "tumor_volume_change", "series_type", "is_highlight_day"))
  a$is_highlight_day <- as.character(a$is_highlight_day) %in% c("TRUE", "T", "1")
  b <- figure7_read_tsv(table_path("panel_7B_plot_data.tsv"),
    c("comparison_id", "panel", "pseudotime", "mean_ecdf", "color_group", "curve_label", "line_group"))
  bt <- figure7_read_tsv(table_path("panel_7B_tests.tsv"), c("comparison_id", "panel", "annotation"))
  cdata <- figure7_read_tsv(table_path("panel_7C_plot_data.tsv"), c("sample_id", "initial_ploidy", "dose", "TGI_percent_Day_17"))
  ct <- figure7_read_tsv(table_path("panel_7C_test.tsv"),
    c("dose_adjusted_difference_high_minus_low", "permutation_p_two_sided", "n_group_low", "n_group_high"))
  d <- figure7_read_tsv(table_path("panel_7D_plot_data.tsv"), c("sample_id", "shift_centered", "tgi_centered", "dose", "etp_group"))
  dt <- figure7_read_tsv(table_path("panel_7D_test.tsv"), c("n", "estimate", "permutation_p_two_sided"))
  e <- figure7_read_tsv(table_path("panel_7E_plot_data.tsv"), c("sample_id", "sample_mean_endpoint_ploidy", "TGI_percent_Day_17", "dose", "etp_group"))
  et <- figure7_read_tsv(table_path("panel_7E_test.tsv"), c("n", "estimate", "permutation_p_two_sided"))
  f <- figure7_read_tsv(table_path("panel_7F_pathway_activity_plot_data.tsv"),
    c("collection_id", "collection_label", "collection_display_order", "pathway_id", "pathway_label",
      "pathway_display_order", "selected_direction", "selected_rank_within_direction",
      "pseudotime", "standardized_activity"))
  for (tab in list(a, b, bt, cdata, ct, d, dt, e, et)) {
    required_metadata <- c("tgi_outcome", "tgi_day", "tgi_measure", "matched_control_summary", "matched_control_group")
    missing_metadata <- setdiff(required_metadata, names(tab))
    if (length(missing_metadata) || any(tab$tgi_outcome != "day") || any(figure7_numeric(tab$tgi_day) != 17) ||
        any(tab$tgi_measure != "TGI_percent_Day_17") || any(tab$matched_control_summary != "mean") ||
        any(tab$matched_control_group != "initial_ploidy")) {
      figure7_stop("render-only source A-E tables have missing or incompatible frozen TGI metadata")
    }
  }
  pathway_meta <- unique(f[, c("collection_id", "pathway_id")])
  if (nrow(pathway_meta) != 24L || any(base::table(pathway_meta$collection_id) != 8L)) {
    figure7_stop("render-only panel-7F table must contain eight pathways in each of three collections")
  }
  f_key <- interaction(f$collection_id, f$pathway_id, drop = TRUE)
  grids <- lapply(levels(f_key), function(key) sort(figure7_numeric(f$pseudotime[f_key == key])))
  if (any(lengths(grids) != as.integer(config$state_pathways$grid_size)) ||
      !all(vapply(grids, identical, logical(1L), grids[[1L]]))) {
    figure7_stop("render-only panel-7F pathways must share the identical frozen 501-point grid")
  }
  provenance <- file.path(source_dir, "metadata", "state_pathway_provenance.tsv")
  if (!file.exists(provenance)) figure7_stop("render-only source is missing state_pathway_provenance.tsv")
  figure7_prepare_output(output_dir)
  for (file in list.files(file.path(source_dir, "tables"), full.names = TRUE)) {
    figure7_copy_file(file, file.path(output_dir, "tables", basename(file)))
  }
  figure7_copy_file(provenance, file.path(output_dir, "metadata", basename(provenance)))
  filenames <- config$panels$filenames
  figure7_save_pdf(figure7_panel_a_plot(a), file.path(output_dir, "figures", filenames[["7A"]]), 10, 6.5)
  figure7_save_pdf(figure7_panel_b_plot(b, bt), file.path(output_dir, "figures", filenames[["7B"]]), 15, 5.5)
  figure7_save_pdf(figure7_panel_c_plot(cdata, ct), file.path(output_dir, "figures", filenames[["7C"]]), 6.8, 6.4)
  figure7_save_pdf(figure7_scatter_plot(d, "shift_centered", "tgi_centered", dt,
    "CellCycle TGI association after within-dose centering", "Dose-centered ECDF RMSE", "Dose-centered Day 17 TGI (%)") +
      ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.35),
    file.path(output_dir, "figures", filenames[["7D"]]), 6.6, 6.6)
  figure7_save_pdf(figure7_scatter_plot(e, "sample_mean_endpoint_ploidy", "TGI_percent_Day_17", et,
    "Cell-cycle-associated tumor cells: Day 17 TGI vs sample mean ETP", "Sample mean ETP", "Day 17 TGI (%)"),
    file.path(output_dir, "figures", filenames[["7E"]]), 6.8, 6.8)
  figure7_save_pdf(figure7_panel_f_plot(f, config), file.path(output_dir, "figures", filenames[["7F"]]), 9, 8)
  write_metadata(output_dir, mode, config, config_path)
  figure7_validate_figure_inventory(output_dir, config)
}

if (identical(mode, "render-only")) {
  render_from_run(figure7_arg(args, "source-run-dir", required = TRUE), output_dir, config, config_path)
  message("Rendered six Figure 7 panels from immutable plotting tables: ", output_dir)
  quit(save = "no", status = 0L)
}

cellcycle_path <- normalizePath(figure7_arg(args, "cellcycle-input", required = TRUE), mustWork = FALSE)
noncellcycle_path <- normalizePath(figure7_arg(args, "non-cellcycle-input", required = TRUE), mustWork = FALSE)
saved_dir <- normalizePath(figure7_arg(args, "saved-state-pathway-dir", required = TRUE), mustWork = FALSE)
figure7_assert_empty_output(output_dir)
figure7_verify_checksum(cellcycle_path, config$inputs$cellcycle_sha256, "CellCycle processed input")
figure7_verify_checksum(noncellcycle_path, config$inputs$noncellcycle_sha256, "NonCellCycle processed input")
reference <- figure7_validate_state_reference(saved_dir, config, verify_checksums = TRUE)

if (identical(mode, "full-analysis")) {
  figure7_full_analysis(
    figure7_arg(args, "seurat-rds", required = TRUE),
    figure7_arg(args, "gene-set-artifact", required = TRUE),
    reference, output_dir, config
  )
}

cellcycle <- figure7_read_cell_table(cellcycle_path, "CellCycle")
noncellcycle <- figure7_read_cell_table(noncellcycle_path, "NonCellCycle")
samples <- figure7_sample_table(cellcycle, noncellcycle, config)
data <- figure7_prepare_cellcycle(cellcycle, samples)
figure7_prepare_output(output_dir)
figure7_build_ae(cellcycle, data, samples, output_dir, config)
figure7_build_f(reference, output_dir, config)
write_metadata(output_dir, mode, config, config_path)
figure7_validate_figure_inventory(output_dir, config)
message("Generated exactly six Figure 7 source panels: ", output_dir)
