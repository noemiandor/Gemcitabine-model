#!/usr/bin/env Rscript

# Export the compact, immutable panel-7F reference from the exact completed
# 04i result tree. This script does not fit a model, query gene sets, or modify
# the source result tree.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "Code/in-vivo/figure7/export_04i_state_pathway_reference.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
sys.source(file.path(script_dir, "src", "common_io.R"), envir = .GlobalEnv)

required_packages <- c("readr", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) {
  figure7_stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

cli <- figure7_parse_args(commandArgs(trailingOnly = TRUE))
allowed_args <- c("results_root", "report_html", "output_dir")
unknown_args <- setdiff(names(cli), allowed_args)
if (length(unknown_args)) figure7_stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))

results_root <- normalizePath(figure7_arg(cli, "results-root", required = TRUE), mustWork = TRUE)
report_html_relative_path <- "report/04i_pseudotime_state_pathways_report.html"
# Accepted for backward-compatible commands, but the HTML file is not read or checksummed.
report_html <- figure7_arg(cli, "report-html", file.path(results_root, report_html_relative_path))
report_html_basename <- basename(report_html)
if (!nzchar(report_html_basename)) report_html_basename <- basename(report_html_relative_path)
output_dir <- normalizePath(
  figure7_arg(
    cli,
    "output-dir",
    file.path(repo_root, "Data", "in-vivo", "figure7", "saved_state_pathway", "taoli_04i_etp2_24_day17_v1")
  ),
  mustWork = FALSE
)

reference_id <- "taoli_04i_etp2_24_day17_v1"
workflow_id <- "binning"
model_id <- "ETP_reference_balanced_threshold_2_24"
code_revision_04i <- "dc751eab928bc40f3edb063baec447fe32a69d73"
expected_report_sha256 <- "b9644b1da0399043a6aba28178a2b61375b661780c4fb08a148c7724da00bfa1"
expected_input_sha256 <- c(
  cell_metadata = "bd0c6fcf3f6691114445a7c8a9fb8de7e55ef6580a9737284cd38f5bba216186",
  noncell_metadata = "1e44794fb2b2c402b49cf9d17abcc54cd7690dc5d1d611e38a6aa3b996093d3c",
  seurat_rds = "727b8a5e5868da911c3b0873838fb1b0023377ed21ea5498dc6493acbbef6d98",
  config = "814e2d0d2d1eb92c62dbd6f00eeac23277bc602cdb5b836da72671c1e97ff49c"
)

if (!identical(basename(output_dir), reference_id)) {
  figure7_stop("Output directory must use canonical reference ID ", reference_id, ": ", output_dir)
}
if (dir.exists(output_dir)) {
  figure7_stop("Canonical output directory already exists; refusing to overwrite: ", output_dir)
}

read_csv <- function(path) {
  if (!file.exists(path)) figure7_stop("Missing 04i source table: ", path)
  as.data.frame(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE, name_repair = "minimal"),
    stringsAsFactors = FALSE
  )
}

assert_columns <- function(data, required, label) {
  missing <- setdiff(required, names(data))
  if (length(missing)) figure7_stop(label, " is missing column(s): ", paste(missing, collapse = ", "))
  invisible(data)
}

parameter_value <- function(data, key, label) {
  hit <- which(as.character(data$parameter) == key)
  if (length(hit) != 1L) figure7_stop(label, " must contain exactly one parameter: ", key)
  as.character(data$value[[hit]])
}

manifest_root <- file.path(results_root, "00_manifest")
workflow_root <- file.path(results_root, workflow_id)
model_root <- file.path(workflow_root, model_id)
source_paths <- c(
  primary_gsea = file.path(model_root, "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"),
  activity = file.path(model_root, "04_gsea", "pathway_activity_over_pseudotime.csv"),
  leading_edge = file.path(model_root, "04_gsea", "all_collections_leading_edge_genes.csv"),
  gene_contrast = file.path(model_root, "03_gene_models", "gene_primary_adjacent_state_contrast.csv"),
  gene_resolution = file.path(model_root, "03_gene_models", "gene_symbol_resolution.csv"),
  sample_bin_metadata = file.path(workflow_root, "02_pseudobulk", "sample_bin_metadata.csv"),
  design_audit = file.path(model_root, "01_qc", "model_design_rank_audit.csv"),
  primary_coverage = file.path(workflow_root, "01_qc", "primary_coverage_check.csv")
)
expected_source_sha256 <- c(
  primary_gsea = "aaca8d189c1288c78502882a15b80a9ca35c98da42bafd89d140d6f4e68faee1",
  activity = "0f545a39cd9ce25e0335c8db2aed69b44834e3c86d11164e492419dada58db76",
  leading_edge = "916a8d3db4b3e936595714c08d4cc666f2e77663f641aaad7f6581b387ebd2aa",
  gene_contrast = "b371b2931dabafa1850ab55313d724a4862b8723c5c9e1f65e5d8b0eab47be5f",
  gene_resolution = "b5973dab047ebf92b645007972f46e6a9bb7629d97eabd8cb3d62d75d2c32100",
  sample_bin_metadata = "89d757f074932eb81951f22adc4f93d1c4e2e62de59f742467c5765f6df011d5",
  design_audit = "1dd5dc5b7e3453dec06780f8092150a7f7c1c31c679596e3eda85019a910146e",
  primary_coverage = "9fecc5339a86cf2e072dccfb05358580fc0db25fd3a8e9f71af301a8ea423141"
)
figure7_verify_named_checksums(source_paths, expected_source_sha256)

analysis_parameters <- read_csv(file.path(workflow_root, "00_manifest", "analysis_parameters.csv"))
intervals <- read_csv(file.path(manifest_root, "frozen_interval_definition.csv"))
input_checksums <- read_csv(file.path(manifest_root, "input_checksums.csv"))
package_versions <- read_csv(file.path(manifest_root, "package_versions.csv"))
model_parameters <- read_csv(file.path(model_root, "00_manifest", "model_parameters.csv"))
assert_columns(analysis_parameters, c("parameter", "value"), "analysis parameters")
assert_columns(intervals, c("interval_id", "start", "end", "include_start", "include_end"), "frozen intervals")
assert_columns(input_checksums, c("input", "path", "sha256"), "input checksums")
assert_columns(package_versions, c("package", "version"), "package versions")
assert_columns(model_parameters, c("parameter", "value"), "model parameters")

observed_inputs <- setNames(as.character(input_checksums$sha256), as.character(input_checksums$input))
if (!identical(unname(observed_inputs[names(expected_input_sha256)]), unname(expected_input_sha256))) {
  figure7_stop("04i input checksums do not match the reviewed result tree")
}
expected_parameters <- c(
  assay = "RNA",
  counts_layer = "counts",
  n_pseudotime_bins = "20",
  min_cells_per_sample_bin = "5",
  spline_df = "5",
  gene_set_collections = "H,C2:CP:REACTOME,C5:GO:BP",
  seed = "1",
  gsea_min_size = "15",
  gsea_max_size = "500",
  grid_size = "501"
)
for (key in names(expected_parameters)) {
  if (!identical(parameter_value(analysis_parameters, key, "analysis parameters"), expected_parameters[[key]])) {
    figure7_stop("Unexpected 04i analysis parameter ", key)
  }
}
if (!identical(parameter_value(model_parameters, "model_id", "model parameters"), model_id) ||
    !identical(parameter_value(model_parameters, "covariate_mode", "model parameters"), "etp_group") ||
    !isTRUE(all.equal(as.numeric(parameter_value(model_parameters, "etp_threshold", "model parameters")), 2.24, tolerance = 0))) {
  figure7_stop("Unexpected 04i ETP model parameters")
}

expected_intervals <- data.frame(
  interval_id = c("primary_accumulated_state", "left_neighbor", "right_neighbor"),
  start = c(0.30, 0.11, 0.49), end = c(0.49, 0.30, 0.68),
  include_start = c(TRUE, TRUE, FALSE), include_end = c(TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(expected_intervals))) {
  local <- intervals[intervals$interval_id == expected_intervals$interval_id[[i]], , drop = FALSE]
  if (nrow(local) != 1L ||
      !isTRUE(all.equal(as.numeric(local$start), expected_intervals$start[[i]], tolerance = 0)) ||
      !isTRUE(all.equal(as.numeric(local$end), expected_intervals$end[[i]], tolerance = 0)) ||
      !identical(as.logical(local$include_start), expected_intervals$include_start[[i]]) ||
      !identical(as.logical(local$include_end), expected_intervals$include_end[[i]])) {
    figure7_stop("Frozen interval mismatch for ", expected_intervals$interval_id[[i]])
  }
}

gsea_path <- source_paths[["primary_gsea"]]
gsea <- read_csv(gsea_path)
assert_columns(
  gsea,
  c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leading_edge", "collection",
    "collection_label", "ranking_id", "pathway_label", "direction", "model_id"),
  "primary-adjacent GSEA"
)
if (nrow(gsea) != 4870L || any(gsea$ranking_id != "primary_adjacent_state") || any(gsea$model_id != model_id)) {
  figure7_stop("Primary-adjacent GSEA does not match the reviewed ETP model")
}

collection_ids <- c("H", "C2:CP:REACTOME", "C5:GO:BP")
collection_labels <- c("H" = "Hallmark", "C2:CP:REACTOME" = "Reactome", "C5:GO:BP" = "GO biological process")
selected_rows <- list()
for (collection_id in collection_ids) {
  local <- gsea[gsea$collection == collection_id & is.finite(gsea$padj) & is.finite(gsea$NES), , drop = FALSE]
  positive <- local[local$NES > 0, , drop = FALSE]
  positive <- positive[order(positive$padj, -positive$NES, positive$pathway), , drop = FALSE]
  positive <- head(positive, 4L)
  positive$selected_direction <- "positive"
  positive$selected_rank_within_direction <- seq_len(nrow(positive))
  negative <- local[local$NES < 0, , drop = FALSE]
  negative <- negative[order(negative$padj, negative$NES, negative$pathway), , drop = FALSE]
  negative <- head(negative, 4L)
  negative$selected_direction <- "negative"
  negative$selected_rank_within_direction <- seq_len(nrow(negative))
  if (nrow(positive) != 4L || nrow(negative) != 4L) {
    figure7_stop("Expected four positive and four negative pathways for ", collection_id)
  }
  selected_rows[[collection_id]] <- rbind(positive, negative)
}
selected <- do.call(rbind, selected_rows)
selected$collection_display_order <- match(selected$collection, collection_ids)
selected$pathway_display_order <- seq_len(nrow(selected))
selected_export <- data.frame(
  collection_id = selected$collection,
  collection_label = unname(collection_labels[selected$collection]),
  collection_display_order = selected$collection_display_order,
  pathway_id = selected$pathway,
  pathway_label = selected$pathway_label,
  pathway_display_order = selected$pathway_display_order,
  NES = selected$NES,
  padj = selected$padj,
  selected_direction = selected$selected_direction,
  selected_rank_within_direction = selected$selected_rank_within_direction,
  stringsAsFactors = FALSE
)
selected_keys <- paste(selected_export$collection_id, selected_export$pathway_id, sep = "\r")
if (nrow(selected_export) != 24L || anyDuplicated(selected_keys)) {
  figure7_stop("Canonical pathway selection must contain 24 unique collection/pathway keys")
}

activity_path <- source_paths[["activity"]]
if (!file.exists(activity_path)) figure7_stop("Missing 04i activity source: ", activity_path)
activity_pieces <- list()
callback <- readr::SideEffectChunkCallback$new(function(chunk, position) {
  assert_columns(chunk, c("collection", "pathway", "ranking_id", "model_id", "pseudotime", "standardized_activity"),
                 "pathway activity chunk")
  key <- paste(chunk$collection, chunk$pathway, sep = "\r")
  keep <- key %in% selected_keys
  if (any(keep)) {
    local <- as.data.frame(chunk[keep, c("collection", "pathway", "ranking_id", "model_id", "pseudotime", "standardized_activity")])
    activity_pieces[[length(activity_pieces) + 1L]] <<- local
  }
})
invisible(readr::read_csv_chunked(
  activity_path,
  callback = callback,
  chunk_size = 50000L,
  show_col_types = FALSE,
  progress = FALSE
))
if (!length(activity_pieces)) figure7_stop("No selected pathways were found in the activity source")
activity <- do.call(rbind, activity_pieces)
activity_key <- paste(activity$collection, activity$pathway, sep = "\r")
metadata_index <- match(activity_key, selected_keys)
if (anyNA(metadata_index) || any(activity$ranking_id != "primary_adjacent_state") || any(activity$model_id != model_id)) {
  figure7_stop("Activity rows do not match the selected primary-adjacent ETP model")
}
activity_export <- cbind(
  selected_export[metadata_index, c(
    "collection_id", "collection_label", "collection_display_order", "pathway_id", "pathway_label",
    "pathway_display_order", "selected_direction", "selected_rank_within_direction"
  ), drop = FALSE],
  pseudotime = as.numeric(activity$pseudotime),
  standardized_activity = as.numeric(activity$standardized_activity)
)
activity_export <- activity_export[order(activity_export$pathway_display_order, activity_export$pseudotime), , drop = FALSE]
rownames(activity_export) <- NULL
activity_groups <- split(activity_export$pseudotime, activity_export$pathway_display_order)
if (nrow(activity_export) != 24L * 501L || any(lengths(activity_groups) != 501L) ||
    !all(vapply(activity_groups, identical, logical(1L), activity_groups[[1L]]))) {
  figure7_stop("Canonical activity table must contain 24 pathways on an identical 501-point grid")
}

leading_path <- source_paths[["leading_edge"]]
leading <- read_csv(leading_path)
assert_columns(leading, c("collection", "pathway", "ranking_id", "leading_edge_gene", "model_id"), "leading-edge genes")
leading$.source_row <- seq_len(nrow(leading))
leading_key <- paste(leading$collection, leading$pathway, sep = "\r")
leading <- leading[
  leading$ranking_id == "primary_adjacent_state" & leading$model_id == model_id & leading_key %in% selected_keys,
  , drop = FALSE
]
leading$pathway_display_order <- match(paste(leading$collection, leading$pathway, sep = "\r"), selected_keys)
leading <- leading[order(leading$pathway_display_order, leading$.source_row), , drop = FALSE]
leading$leading_edge_rank <- ave(seq_len(nrow(leading)), leading$pathway_display_order, FUN = seq_along)
leading_export <- data.frame(
  collection_id = leading$collection,
  pathway_id = leading$pathway,
  leading_edge_gene_id = trimws(as.character(leading$leading_edge_gene)),
  leading_edge_rank = as.integer(leading$leading_edge_rank),
  stringsAsFactors = FALSE
)
leading_counts <- table(paste(leading_export$collection_id, leading_export$pathway_id, sep = "\r"))
if (any(!selected_keys %in% names(leading_counts)) || any(leading_counts[selected_keys] < 1L)) {
  figure7_stop("Every selected pathway must have at least one leading-edge gene")
}

contrast_path <- source_paths[["gene_contrast"]]
resolution_path <- source_paths[["gene_resolution"]]
contrast <- read_csv(contrast_path)
resolution <- read_csv(resolution_path)
assert_columns(
  contrast,
  c("gene", "gene_symbol", "contrast_id", "fitted_log_expr_difference", "average_expression", "t_statistic",
    "p_value", "fdr", "B", "model_id"),
  "gene contrast"
)
assert_columns(
  resolution,
  c("gene", "gene_symbol", "contrast_id", "t_statistic", "p_value", "fdr", "resolution_rank",
    "retained_for_gsea", "model_id"),
  "gene-symbol resolution"
)
if (nrow(contrast) != 28922L || nrow(resolution) != 28922L) figure7_stop("Unexpected complete gene-table row count")
retained <- resolution[as.logical(resolution$retained_for_gsea), , drop = FALSE]
contrast_index <- match(retained$gene, contrast$gene)
if (anyNA(contrast_index)) figure7_stop("Retained GSEA genes are missing from the contrast table")
ranked <- contrast[contrast_index, , drop = FALSE]
ranked$resolution_rank <- retained$resolution_rank
ranked <- ranked[order(-ranked$t_statistic, ranked$gene_symbol), , drop = FALSE]
if (nrow(ranked) != 28908L || anyDuplicated(ranked$gene_symbol)) {
  figure7_stop("Canonical GSEA ranking must contain 28,908 unique gene symbols")
}
gene_ranking_export <- data.frame(
  rank = seq_len(nrow(ranked)),
  gene_id = ranked$gene,
  gene_symbol = ranked$gene_symbol,
  moderated_t = ranked$t_statistic,
  p_value = ranked$p_value,
  fdr = ranked$fdr,
  fitted_log_expr_difference = ranked$fitted_log_expr_difference,
  average_expression = ranked$average_expression,
  B = ranked$B,
  contrast_id = ranked$contrast_id,
  model_id = ranked$model_id,
  symbol_resolution_rank = ranked$resolution_rank,
  retained_for_gsea = TRUE,
  stringsAsFactors = FALSE
)

gsea_complete_export <- data.frame(
  collection_id = gsea$collection,
  collection_label = unname(collection_labels[gsea$collection]),
  pathway_id = gsea$pathway,
  pathway_label = gsea$pathway_label,
  ranking_id = gsea$ranking_id,
  pval = gsea$pval,
  padj = gsea$padj,
  log2err = gsea$log2err,
  ES = gsea$ES,
  NES = gsea$NES,
  size = gsea$size,
  direction = gsea$direction,
  leading_edge_genes = gsea$leading_edge,
  model_id = gsea$model_id,
  stringsAsFactors = FALSE
)
if (anyNA(gsea_complete_export$collection_label) || anyDuplicated(gsea_complete_export[c("collection_id", "pathway_id")])) {
  figure7_stop("Complete GSEA table contains unexpected collections or duplicated pathway keys")
}

sample_bins_path <- source_paths[["sample_bin_metadata"]]
sample_bins <- read_csv(sample_bins_path)
assert_columns(
  sample_bins,
  c("sample_id", "bin_id", "sample_bin_id", "cell_count", "bin_start", "bin_end", "bin_midpoint",
    "dose", "dose_mg", "initial_ploidy", "retained_for_model", "exclusion_reason", "library_size",
    "ETP_reference_balanced_threshold_2_24_group"),
  "sample-bin metadata"
)
if (nrow(sample_bins) != 320L) figure7_stop("Unexpected sample-bin metadata row count")
names(sample_bins)[names(sample_bins) == "bin_id"] <- "pseudotime_bin"
sample_bin_export <- sample_bins

design_path <- source_paths[["design_audit"]]
coverage_path <- source_paths[["primary_coverage"]]
design <- read_csv(design_path)
coverage <- read_csv(coverage_path)
assert_columns(
  design,
  c("model_id", "covariate_mode", "etp_method", "etp_threshold", "include_dose", "n_observations",
    "n_design_columns_original", "design_rank", "rank_deficient", "retained_design_columns",
    "dropped_design_columns", "duplicate_correlation"),
  "design audit"
)
assert_columns(coverage, c("n_contributing_mice", "has_control", "has_treated", "n_initial_ploidy_groups", "passes"),
               "primary coverage check")
if (nrow(design) != 1L || nrow(coverage) != 1L || design$n_observations[[1L]] != 211L ||
    coverage$n_contributing_mice[[1L]] != 14L || !isTRUE(as.logical(coverage$passes[[1L]]))) {
  figure7_stop("Design or primary-coverage audit does not match the reviewed run")
}
design_qc_export <- cbind(
  design,
  n_sample_bins_total = nrow(sample_bins),
  n_sample_bins_retained = sum(as.logical(sample_bins$retained_for_model)),
  n_contributing_mice_primary_interval = coverage$n_contributing_mice,
  has_control_primary_interval = coverage$has_control,
  has_treated_primary_interval = coverage$has_treated,
  n_initial_ploidy_groups_primary_interval = coverage$n_initial_ploidy_groups,
  primary_coverage_passes = coverage$passes,
  minimum_expression_observations = coverage$n_contributing_mice,
  n_genes_passing_expression_filter = nrow(contrast),
  n_unique_gene_symbols_ranked_for_gsea = nrow(gene_ranking_export),
  n_pathways_tested = nrow(gsea_complete_export)
)

staging_dir <- paste0(output_dir, ".tmp-", Sys.getpid())
if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(staging_dir)) figure7_stop("Could not create staging directory: ", staging_dir)
completed <- FALSE
on.exit(if (!completed && dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)

write_reference <- function(data, filename) {
  figure7_write_tsv(as.data.frame(data, stringsAsFactors = FALSE), file.path(staging_dir, filename))
}
write_reference(activity_export, "panel_7F_pathway_activity_plot_data.tsv")
write_reference(selected_export, "panel_7F_selected_pathway_gsea.tsv")
write_reference(leading_export, "panel_7F_leading_edge_genes.tsv")
write_reference(gene_ranking_export, "state_pathway_gene_ranking_complete.tsv")
write_reference(gsea_complete_export, "state_pathway_gsea_complete.tsv")
write_reference(sample_bin_export, "state_pathway_sample_bin_coverage.tsv")
write_reference(design_qc_export, "state_pathway_design_qc.tsv")

package_value <- function(package) {
  hit <- which(as.character(package_versions$package) == package)
  if (length(hit) != 1L) figure7_stop("Package manifest must contain exactly one row for ", package)
  as.character(package_versions$version[[hit]])
}
input_value <- function(input, column) {
  hit <- which(as.character(input_checksums$input) == input)
  if (length(hit) != 1L) figure7_stop("Input manifest must contain exactly one row for ", input)
  as.character(input_checksums[[column]][[hit]])
}

activity_output_path <- file.path(staging_dir, "panel_7F_pathway_activity_plot_data.tsv")
provenance <- c(
  full_analysis_run_dir = parameter_value(analysis_parameters, "output_root", "analysis parameters"),
  report_identifier = paste0(report_html_basename, "@sha256:", expected_report_sha256),
  code_revision_04i = code_revision_04i,
  seurat_rds_sha256 = input_value("seurat_rds", "sha256"),
  cellcycle_metadata_sha256 = input_value("cell_metadata", "sha256"),
  noncellcycle_metadata_sha256 = input_value("noncell_metadata", "sha256"),
  interval_config_sha256 = input_value("config", "sha256"),
  assay = "RNA",
  counts_layer = "counts",
  etp_method = "reference_balanced",
  etp_threshold = "2.24",
  spline_df = "5",
  pseudotime_bins = "20",
  minimum_cells_per_sample_bin = "5",
  grid_size = "501",
  seed = "1",
  gene_set_source = paste0("msigdbr runtime query; package version ", package_value("msigdbr")),
  gene_set_release = "not_recorded_in_04i_manifest",
  gene_set_species = "Homo sapiens",
  gene_set_collections = "H,C2:CP:REACTOME,C5:GO:BP",
  gene_set_min_size = "15",
  gene_set_max_size = "500",
  expression_filter = "CPM > 1 in at least 14 retained sample-bin observations",
  normalization = "TMM",
  observation_model = "voom",
  mouse_block = "duplicateCorrelation",
  nuisance_terms = "dose_mg_factor,ETP_reference_balanced_threshold_2_24_factor",
  treatment_by_pseudotime_interaction = "false",
  empirical_bayes = "robust",
  contrast = "mean(primary grid) - 0.5 * mean(left grid) - 0.5 * mean(right grid)",
  gsea_rank_statistic = "moderated_t",
  pathway_activity = "row-standardize fitted expression per leading-edge gene, then average genes",
  feature_species_policy = paste(
    "mixed GRCh38 and GRCm39 features; strip GRCh/GRCm/hg38 prefixes and version suffixes;",
    "preserve symbol case; resolve duplicate symbols by largest absolute moderated t then p value;",
    "query Homo sapiens MSigDB gene symbols"
  ),
  gsea_ranking_rule = "one row per cleaned gene symbol; descending moderated t; ties retain cleaned-symbol order",
  pathway_selection_rule = "finite adjusted P; no FDR cutoff; top four per sign and collection",
  activity_table_sha256 = figure7_sha256(activity_output_path),
  canonical_reference_id = reference_id,
  workflow_id = workflow_id,
  model_id = model_id,
  accumulated_interval = "[0.30,0.49]",
  left_neighbor_interval = "[0.11,0.30)",
  right_neighbor_interval = "(0.49,0.68]",
  report_html_sha256 = expected_report_sha256,
  report_html_relative_path = report_html_relative_path,
  source_results_id = "04i_pseudotime_state_pathways",
  recorded_cell_metadata_path = input_value("cell_metadata", "path"),
  recorded_noncell_metadata_path = input_value("noncell_metadata", "path"),
  recorded_seurat_rds_path = input_value("seurat_rds", "path"),
  msigdbr_package_version = package_value("msigdbr"),
  fgsea_package_version = package_value("fgsea"),
  source_activity_sha256 = expected_source_sha256[["activity"]],
  source_primary_gsea_sha256 = expected_source_sha256[["primary_gsea"]],
  source_leading_edge_sha256 = expected_source_sha256[["leading_edge"]],
  source_gene_contrast_sha256 = expected_source_sha256[["gene_contrast"]],
  source_gene_resolution_sha256 = expected_source_sha256[["gene_resolution"]],
  source_sample_bin_metadata_sha256 = expected_source_sha256[["sample_bin_metadata"]],
  source_design_audit_sha256 = expected_source_sha256[["design_audit"]],
  source_primary_coverage_sha256 = expected_source_sha256[["primary_coverage"]]
)
provenance_export <- data.frame(key = names(provenance), value = unname(provenance), stringsAsFactors = FALSE)
write_reference(provenance_export, "state_pathway_provenance.tsv")

expected_files <- c(
  "panel_7F_pathway_activity_plot_data.tsv",
  "panel_7F_selected_pathway_gsea.tsv",
  "panel_7F_leading_edge_genes.tsv",
  "state_pathway_gene_ranking_complete.tsv",
  "state_pathway_gsea_complete.tsv",
  "state_pathway_sample_bin_coverage.tsv",
  "state_pathway_design_qc.tsv",
  "state_pathway_provenance.tsv"
)
observed_files <- sort(list.files(staging_dir, all.files = FALSE, recursive = FALSE))
if (!identical(observed_files, sort(expected_files)) || any(file.info(file.path(staging_dir, expected_files))$size <= 0)) {
  figure7_stop("Canonical staging directory does not contain exactly eight nonempty reference tables")
}

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
if (!file.rename(staging_dir, output_dir)) figure7_stop("Could not atomically materialize canonical reference: ", output_dir)
completed <- TRUE
message("Exported canonical 04i panel-7F reference: ", output_dir)
message("Recorded reviewed report SHA-256: ", expected_report_sha256)
message("Activity rows: ", nrow(activity_export), "; selected pathways: ", nrow(selected_export))
