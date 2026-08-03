#!/usr/bin/env Rscript

# Export a compact, run-scoped panel-7F reference from this repository's
# generated state-pathway support tree. The historical byte-pinned 04i audit
# export uses a separate exporter and provenance contract.

args_full <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_path <- if (length(file_arg)) {
  sub("^--file=", "", file_arg[[1L]])
} else {
  "Code/in-vivo/figure7/export_state_pathway_reference.R"
}
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)
sys.source(file.path(script_dir, "src", "common_io.R"), envir = .GlobalEnv)
sys.source(
  file.path(script_dir, "src", "feature_species_policy.R"),
  envir = .GlobalEnv
)
sys.source(
  file.path(script_dir, "src", "tgi_statistics.R"),
  envir = .GlobalEnv
)

required_packages <- c("readr", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) {
  figure7_stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "))
}

cli <- figure7_parse_args(commandArgs(trailingOnly = TRUE))
allowed_args <- c("results_root", "output_dir", "config")
unknown_args <- setdiff(names(cli), allowed_args)
if (length(unknown_args)) figure7_stop("Unknown argument(s): ", paste(unknown_args, collapse = ", "))

results_root <- normalizePath(figure7_arg(cli, "results-root", required = TRUE), mustWork = TRUE)
config_path <- normalizePath(
  figure7_arg(cli, "config", file.path(script_dir, "figure7_config.yaml")),
  mustWork = TRUE
)
config <- figure7_read_config(config_path)
config <- figure7_attach_density_localization_config(
  config,
  file.path(script_dir, "density_localization_config.yaml")
)
output_dir <- normalizePath(
  figure7_arg(
    cli,
    "output-dir",
    file.path(
      repo_root, "Results", "in-vivo", "figure7", "intermediates",
      as.character(config$state_pathways$generated_reference_id)
    )
  ),
  mustWork = FALSE
)

reference_id <- as.character(config$state_pathways$generated_reference_id)
workflow_id <- "binning"
model_id <- as.character(config$state_pathways$model)
recorded_source_code_revision <- paste0(
  "sha256:",
  figure7_sha256(
    file.path(script_dir, "generate_pseudotime_state_pathways_support.R")
  )
)

if (!identical(basename(output_dir), reference_id)) {
  figure7_stop("Output directory must use generated reference ID ", reference_id, ": ", output_dir)
}
if (dir.exists(output_dir)) {
  figure7_stop("Generated output directory already exists; refusing to overwrite: ", output_dir)
}

read_csv <- function(path) {
  if (!file.exists(path)) figure7_stop("Missing state-pathway source table: ", path)
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
  feature_species_audit =
    file.path(manifest_root, "feature_species_audit.csv"),
  interval_definition =
    file.path(manifest_root, "state_interval_definition.csv"),
  interval_localization_grid =
    file.path(manifest_root, "state_interval_localization_grid.csv"),
  interval_localization_support =
    file.path(manifest_root, "state_interval_localization_support.csv"),
  interval_localization_test =
    file.path(manifest_root, "state_interval_localization_test.csv"),
  cell_expression_match_audit =
    file.path(manifest_root, "cell_expression_match_audit.csv"),
  primary_gsea = file.path(model_root, "04_gsea", "all_collections_primary_adjacent_state_gsea.csv"),
  activity = file.path(model_root, "04_gsea", "pathway_activity_over_pseudotime.csv"),
  leading_edge = file.path(model_root, "04_gsea", "all_collections_leading_edge_genes.csv"),
  gene_contrast = file.path(model_root, "03_gene_models", "gene_primary_adjacent_state_contrast.csv"),
  gene_resolution = file.path(model_root, "03_gene_models", "gene_symbol_resolution.csv"),
  sample_bin_metadata = file.path(workflow_root, "02_pseudobulk", "sample_bin_metadata.csv"),
  design_audit = file.path(model_root, "01_qc", "model_design_rank_audit.csv"),
  primary_coverage = file.path(workflow_root, "01_qc", "primary_coverage_check.csv")
)
missing_source_paths <- source_paths[!file.exists(source_paths)]
if (length(missing_source_paths)) {
  figure7_stop("Missing state-pathway source table(s): ", paste(missing_source_paths, collapse = ", "))
}
observed_source_sha256 <- vapply(source_paths, figure7_sha256, character(1L))

analysis_parameters <- read_csv(file.path(workflow_root, "00_manifest", "analysis_parameters.csv"))
intervals <- read_csv(source_paths[["interval_definition"]])
interval_localization_grid <- read_csv(
  source_paths[["interval_localization_grid"]]
)
interval_localization_support <- read_csv(
  source_paths[["interval_localization_support"]]
)
interval_localization_test <- read_csv(
  source_paths[["interval_localization_test"]]
)
cell_expression_match_audit <- read_csv(
  source_paths[["cell_expression_match_audit"]]
)
input_checksums <- read_csv(file.path(manifest_root, "input_checksums.csv"))
package_versions <- read_csv(file.path(manifest_root, "package_versions.csv"))
gene_set_contract <- read_csv(file.path(manifest_root, "gene_set_contract.csv"))
feature_species_audit <- read_csv(
  source_paths[["feature_species_audit"]]
)
gene_set_membership_path <- file.path(manifest_root, "gene_set_membership.csv")
if (!file.exists(gene_set_membership_path)) {
  figure7_stop("Missing generated gene-set membership manifest")
}
model_parameters <- read_csv(file.path(model_root, "00_manifest", "model_parameters.csv"))
assert_columns(analysis_parameters, c("parameter", "value"), "analysis parameters")
assert_columns(
  intervals,
  c(
    "interval_id", "start", "end", "include_start", "include_end",
    "derivation_analysis_id", "derivation_support_type",
    "derivation_pointwise_alpha", "derivation_n_permutations"
  ),
  "computed state intervals"
)
assert_columns(
  interval_localization_grid,
  c(
    "pseudotime", "treated_minus_vehicle_density",
    "pointwise_p_two_sided", "pointwise_positive_supported",
    "modeled_state_interval"
  ),
  "state-interval localization grid"
)
assert_columns(
  interval_localization_support,
  c("support_type", "start", "end", "alpha"),
  "state-interval localization support"
)
assert_columns(
  interval_localization_test,
  c(
    "analysis_id", "n_cells", "n_samples", "n_permutations",
    "pointwise_test", "pointwise_alpha", "pointwise_start",
    "pointwise_end"
  ),
  "state-interval localization test"
)
assert_columns(
  cell_expression_match_audit,
  c("source", "n_ids", "n_matched", "n_unmatched", "match_rate"),
  "cell-expression match audit"
)
assert_columns(input_checksums, c("input", "locator", "sha256"), "input checksums")
assert_columns(package_versions, c("package", "version"), "package versions")
assert_columns(gene_set_contract, c("key", "value"), "gene-set contract")
assert_columns(
  feature_species_audit,
  c(
    "policy_id", "policy", "human_prefix", "mouse_prefix",
    "n_input_features", "n_human_features_retained",
    "n_mouse_features_excluded", "n_ambiguous_features"
  ),
  "feature-species audit"
)
assert_columns(model_parameters, c("parameter", "value"), "model parameters")
if (anyDuplicated(gene_set_contract$key)) {
  figure7_stop("Gene-set contract contains duplicate keys")
}
gene_set_value <- setNames(
  as.character(gene_set_contract$value),
  as.character(gene_set_contract$key)
)
required_gene_set_keys <- c(
  "provider", "msigdbr_package_version", "database_release", "species",
  "collections", "membership_sha256"
)
if (!all(required_gene_set_keys %in% names(gene_set_value)) ||
    !identical(gene_set_value[["provider"]], "msigdbr") ||
    !identical(
      gene_set_value[["msigdbr_package_version"]],
      as.character(config$gene_sets$package_version)
    ) ||
    !identical(
      gene_set_value[["database_release"]],
      as.character(config$gene_sets$database_release)
    ) ||
    !identical(
      gene_set_value[["species"]],
      as.character(config$gene_sets$species)
    ) ||
    !identical(
      gene_set_value[["collections"]],
      paste(as.character(unlist(config$state_pathways$collections)), collapse = ",")
    ) ||
    !identical(
      gene_set_value[["membership_sha256"]],
      figure7_sha256(gene_set_membership_path)
    )) {
  figure7_stop("Generated gene-set contract does not match config/membership")
}

observed_inputs <- setNames(as.character(input_checksums$sha256), as.character(input_checksums$input))
if (!all(c(
  "cell_metadata", "noncell_metadata", "seurat_rds", "config",
  "support_script", "density_localization_config",
  "feature_species_policy_code",
  "density_localization_code", "common_io_code"
) %in%
         names(observed_inputs)) ||
    any(!grepl("^[0-9a-f]{64}$", observed_inputs))) {
  figure7_stop("State-pathway input checksum lineage is incomplete")
}
current_source_inputs <- c(
  config = config_path,
  support_script = file.path(
    script_dir,
    "generate_pseudotime_state_pathways_support.R"
  ),
  density_localization_config = file.path(
    script_dir,
    "density_localization_config.yaml"
  ),
  feature_species_policy_code = file.path(
    script_dir,
    "src",
    "feature_species_policy.R"
  ),
  density_localization_code = file.path(
    script_dir,
    "src",
    "tgi_statistics.R"
  ),
  common_io_code = file.path(script_dir, "src", "common_io.R")
)
current_source_hashes <- vapply(
  current_source_inputs,
  figure7_sha256,
  character(1L)
)
if (!identical(
      unname(observed_inputs[names(current_source_hashes)]),
      unname(current_source_hashes)
    )) {
  figure7_stop(
    "State-pathway support results were generated with different code or config and cannot be re-attested by this exporter"
  )
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
  gsea_nperm_simple =
    as.character(config$state_pathways$gsea_nperm_simple),
  gsea_nperm_simple_max =
    as.character(config$state_pathways$gsea_nperm_simple_max),
  gsea_nperm_simple_multiplier =
    as.character(config$state_pathways$gsea_nperm_simple_multiplier),
  grid_size = "501",
  min_match_rate = "1",
  feature_species_policy_id =
    as.character(config$feature_species$policy_id),
  human_feature_prefix =
    as.character(config$feature_species$human_prefix),
  mouse_feature_prefix =
    as.character(config$feature_species$mouse_prefix),
  unknown_feature_policy =
    as.character(config$feature_species$unknown_feature_policy)
)
for (key in names(expected_parameters)) {
  if (!identical(parameter_value(analysis_parameters, key, "analysis parameters"), expected_parameters[[key]])) {
    figure7_stop("Unexpected state-pathway analysis parameter ", key)
  }
}
if (nrow(feature_species_audit) != 1L ||
    !identical(
      as.character(feature_species_audit$policy_id),
      as.character(config$feature_species$policy_id)
    ) ||
    !identical(
      as.character(feature_species_audit$human_prefix),
      as.character(config$feature_species$human_prefix)
    ) ||
    !identical(
      as.character(feature_species_audit$mouse_prefix),
      as.character(config$feature_species$mouse_prefix)
    )) {
  figure7_stop("Generated feature-species audit/config identity is invalid")
}
species_counts <- c(
  input = as.integer(feature_species_audit$n_input_features),
  human = as.integer(feature_species_audit$n_human_features_retained),
  mouse = as.integer(feature_species_audit$n_mouse_features_excluded),
  ambiguous = as.integer(feature_species_audit$n_ambiguous_features)
)
if (anyNA(species_counts) ||
    species_counts[["human"]] < 10000L ||
    species_counts[["mouse"]] < 0L ||
    species_counts[["ambiguous"]] != 0L ||
    species_counts[["input"]] !=
      species_counts[["human"]] + species_counts[["mouse"]]) {
  figure7_stop("Generated feature-species audit counts are invalid")
}
if (!identical(
      model_id,
      "initial_ploidy_adjusted_grch_human_only_pointwise_interval_v4"
    ) ||
    !identical(parameter_value(model_parameters, "model_id", "model parameters"), model_id) ||
    !identical(parameter_value(model_parameters, "covariate_mode", "model parameters"), "initial_ploidy") ||
    !identical(parameter_value(model_parameters, "covariate_terms", "model parameters"), "initial_ploidy_factor") ||
    !identical(parameter_value(model_parameters, "initial_ploidy_levels", "model parameters"), "2N;4N")) {
  figure7_stop("Unexpected state-pathway injected-initial-ploidy model parameters")
}

derived_intervals <- figure7_state_intervals_from_density_support(
  interval_localization_support
)
expected_intervals <- do.call(rbind, lapply(derived_intervals, function(x) {
  data.frame(
    interval_id = x$name,
    start = x$start,
    end = x$end,
    include_start = x$include_start,
    include_end = x$include_end,
    stringsAsFactors = FALSE
  )
}))
pointwise_grid <- as.logical(
  interval_localization_grid$pointwise_positive_supported
)
modeled_grid <- as.logical(interval_localization_grid$modeled_state_interval)
pointwise_p <- as.numeric(interval_localization_grid$pointwise_p_two_sided)
localization_pseudotime <- as.numeric(
  interval_localization_grid$pseudotime
)
localization_density <- as.numeric(
  interval_localization_grid$treated_minus_vehicle_density
)
recomputed_pointwise_grid <-
  localization_density > 0 & pointwise_p <= 0.05
recomputed_pointwise_intervals <- figure7_supported_intervals(
  localization_pseudotime,
  recomputed_pointwise_grid,
  "positive_pointwise_two_sided",
  0.05
)
reported_pointwise_intervals <- interval_localization_support[
  interval_localization_support$support_type ==
    "positive_pointwise_two_sided",
  c("support_type", "start", "end", "width", "alpha"),
  drop = FALSE
]
if (nrow(interval_localization_test) != 1L ||
    as.integer(interval_localization_test$n_cells[[1L]]) != 2881L ||
    as.integer(interval_localization_test$n_samples[[1L]]) != 16L ||
    as.integer(interval_localization_test$n_permutations[[1L]]) != 4900L ||
    !identical(
      as.character(interval_localization_test$pointwise_test[[1L]]),
      "two-sided exact permutation at each grid point"
    ) ||
    !isTRUE(all.equal(
      as.numeric(interval_localization_test$pointwise_alpha[[1L]]),
      0.05,
      tolerance = 0
    )) ||
    !identical(pointwise_grid, recomputed_pointwise_grid) ||
    !identical(pointwise_grid, modeled_grid) ||
    !isTRUE(all.equal(
      reported_pointwise_intervals,
      recomputed_pointwise_intervals,
      tolerance = 1e-12,
      check.attributes = FALSE
    )) ||
    !isTRUE(all.equal(
      as.numeric(interval_localization_test$pointwise_start[[1L]]),
      min(localization_pseudotime[pointwise_grid]),
      tolerance = 1e-12
    )) ||
    !isTRUE(all.equal(
      as.numeric(interval_localization_test$pointwise_end[[1L]]),
      max(localization_pseudotime[pointwise_grid]),
      tolerance = 1e-12
    )) ||
    any(!is.finite(pointwise_p[modeled_grid])) ||
    any(pointwise_p[modeled_grid] > 0.05)) {
  figure7_stop(
    "Computed state interval does not reproduce its exact pointwise density-support contract"
  )
}
metadata_match <- cell_expression_match_audit[
  cell_expression_match_audit$source == "metadata",
  ,
  drop = FALSE
]
if (nrow(metadata_match) != 1L ||
    as.integer(metadata_match$n_ids[[1L]]) != 2881L ||
    as.integer(metadata_match$n_matched[[1L]]) != 2881L ||
    as.integer(metadata_match$n_unmatched[[1L]]) != 0L ||
    !isTRUE(all.equal(
      as.numeric(metadata_match$match_rate[[1L]]),
      1,
      tolerance = 0
    ))) {
  figure7_stop(
    "The expression model does not use the exact 2,881-cell localization universe"
  )
}
for (i in seq_len(nrow(expected_intervals))) {
  local <- intervals[intervals$interval_id == expected_intervals$interval_id[[i]], , drop = FALSE]
  if (nrow(local) != 1L ||
      !isTRUE(all.equal(as.numeric(local$start), expected_intervals$start[[i]], tolerance = 0)) ||
      !isTRUE(all.equal(as.numeric(local$end), expected_intervals$end[[i]], tolerance = 0)) ||
      !identical(as.logical(local$include_start), expected_intervals$include_start[[i]]) ||
      !identical(as.logical(local$include_end), expected_intervals$include_end[[i]])) {
    figure7_stop("Computed interval mismatch for ", expected_intervals$interval_id[[i]])
  }
}
if (any(as.character(intervals$derivation_analysis_id) !=
        as.character(interval_localization_test$analysis_id[[1L]])) ||
    any(as.character(intervals$derivation_support_type) !=
        "positive_pointwise_two_sided") ||
    any(as.numeric(intervals$derivation_pointwise_alpha) != 0.05) ||
    any(as.integer(intervals$derivation_n_permutations) != 4900L)) {
  figure7_stop("Computed interval definition is detached from localization inference")
}

gsea_path <- source_paths[["primary_gsea"]]
gsea <- read_csv(gsea_path)
assert_columns(
  gsea,
  c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leading_edge", "collection",
    "collection_label", "ranking_id", "pathway_label", "direction",
    "nPermSimple", "retry_round", "model_id"),
  "primary-adjacent GSEA"
)
if (!nrow(gsea) || any(gsea$ranking_id != "primary_adjacent_state") || any(gsea$model_id != model_id)) {
  figure7_stop("Primary-adjacent GSEA does not match the initial-ploidy model")
}
if (nrow(gsea) < 100L) {
  figure7_stop("Generated primary-adjacent GSEA table is implausibly small")
}
collection_ids <- c("H", "C2:CP:REACTOME", "C5:GO:BP")
if (!setequal(unique(as.character(gsea$collection)), collection_ids)) {
  figure7_stop("Generated primary-adjacent GSEA collections are incomplete")
}
finite_gsea_columns <- c("pval", "padj", "ES", "NES", "size")
if (any(!vapply(
      gsea[finite_gsea_columns],
      function(column) all(is.finite(as.numeric(column))),
      logical(1L)
    ))) {
  figure7_stop(
    "Generated primary-adjacent GSEA contains unresolved nonfinite results"
  )
}
if (anyDuplicated(gsea[c("collection", "pathway")])) {
  figure7_stop("Generated primary-adjacent GSEA contains duplicate pathways")
}
for (collection_id in collection_ids) {
  local <- gsea[gsea$collection == collection_id, , drop = FALSE]
  expected_padj <- stats::p.adjust(local$pval, method = "BH")
  if (!isTRUE(all.equal(
        as.numeric(local$padj),
        expected_padj,
        tolerance = 1e-12,
        check.attributes = FALSE
      ))) {
    figure7_stop(
      "Generated primary-adjacent GSEA adjusted P values are not ",
      "collection-wide BH values for ",
      collection_id
    )
  }
}

collection_labels <- c(
  "H" = "Hallmark",
  "C2:CP:REACTOME" = "Reactome",
  "C5:GO:BP" = "GO biological process"
)
selection_source <- gsea
selection_source$collection_id <- as.character(gsea$collection)
selection_source$pathway_id <- as.character(gsea$pathway)
selected <- figure7_select_generated_pathways(selection_source, config)
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
selected_counts <- table(factor(
  selected_export$collection_id,
  levels = collection_ids
))
max_per_collection <-
  as.integer(config$state_pathways$top_positive_per_collection) +
  as.integer(config$state_pathways$top_negative_per_collection)
if (!nrow(selected_export) || anyDuplicated(selected_keys) ||
    any(selected_counts < 1L) ||
    any(selected_counts > max_per_collection) ||
    any(!is.finite(selected_export$padj)) ||
    any(selected_export$padj >
      figure7_generated_pathway_fdr_threshold())) {
  figure7_stop(
    "Generated pathway selection must contain only unique ",
    "FDR-significant pathways within configured per-collection limits"
  )
}

activity_path <- source_paths[["activity"]]
if (!file.exists(activity_path)) figure7_stop("Missing state-pathway activity source: ", activity_path)
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
  figure7_stop("Activity rows do not match the selected initial-ploidy model")
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
expected_activity_rows <-
  nrow(selected_export) * as.integer(config$state_pathways$grid_size)
if (nrow(activity_export) != expected_activity_rows ||
    any(lengths(activity_groups) !=
      as.integer(config$state_pathways$grid_size)) ||
    !all(vapply(activity_groups, identical, logical(1L), activity_groups[[1L]]))) {
  figure7_stop(
    "Generated activity table must contain every selected pathway on ",
    "the identical configured grid"
  )
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
if (nrow(contrast) < 10000L || nrow(resolution) != nrow(contrast)) {
  figure7_stop("Generated complete gene tables are implausibly small or misaligned")
}
figure7_assert_human_feature_names(
  as.character(contrast$gene),
  analysis = "Generated panel-7F contrast",
  human_prefix = as.character(config$feature_species$human_prefix),
  mouse_prefix = as.character(config$feature_species$mouse_prefix)
)
retained <- resolution[as.logical(resolution$retained_for_gsea), , drop = FALSE]
contrast_index <- match(retained$gene, contrast$gene)
if (anyNA(contrast_index)) figure7_stop("Retained GSEA genes are missing from the contrast table")
ranked <- contrast[contrast_index, , drop = FALSE]
ranked$resolution_rank <- retained$resolution_rank
ranked <- ranked[order(-ranked$t_statistic, ranked$gene_symbol), , drop = FALSE]
if (nrow(ranked) < 10000L || anyDuplicated(ranked$gene_symbol)) {
  figure7_stop(
    "Generated human-only GSEA ranking is implausibly small or has ",
    "duplicate gene symbols"
  )
}
figure7_assert_human_feature_names(
  as.character(ranked$gene),
  analysis = "Generated panel-7F ranking",
  human_prefix = as.character(config$feature_species$human_prefix),
  mouse_prefix = as.character(config$feature_species$mouse_prefix)
)
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
  nPermSimple = gsea$nPermSimple,
  retry_round = gsea$retry_round,
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
    "dose", "dose_mg", "initial_ploidy", "retained_for_model", "exclusion_reason", "library_size"),
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
  c("model_id", "covariate_mode", "initial_ploidy_levels", "include_dose", "n_observations",
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
retained_design_columns <- strsplit(
  as.character(design$retained_design_columns[[1L]]),
  ";",
  fixed = TRUE
)[[1L]]
if (!identical(as.character(design$model_id[[1L]]), model_id) ||
    !identical(as.character(design$covariate_mode[[1L]]), "initial_ploidy") ||
    !identical(as.character(design$initial_ploidy_levels[[1L]]), "2N;4N") ||
    !"initial_ploidy_factor4N" %in% retained_design_columns ||
    any(grepl("ETP|endpoint|copy.number|cn_score", retained_design_columns,
              ignore.case = TRUE))) {
  figure7_stop(
    "Design audit must include injected initial ploidy and prohibit endpoint-CN-score covariates"
  )
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
write_reference(intervals, "state_pathway_interval_definition.tsv")

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
gene_set_value <- function(key) {
  hit <- which(as.character(gene_set_contract$key) == key)
  if (length(hit) != 1L) {
    figure7_stop("Gene-set contract must contain exactly one row for ", key)
  }
  as.character(gene_set_contract$value[[hit]])
}
expected_gene_set_contract <- c(
  provider = "msigdbr",
  msigdbr_package_version = as.character(config$gene_sets$package_version),
  database_release = as.character(config$gene_sets$database_release),
  species = as.character(config$gene_sets$species),
  collections = paste(
    as.character(unlist(config$state_pathways$collections)),
    collapse = ","
  ),
  membership_sha256 = figure7_sha256(gene_set_membership_path)
)
observed_gene_set_contract <- vapply(
  names(expected_gene_set_contract),
  gene_set_value,
  character(1L)
)
if (!identical(
  unname(observed_gene_set_contract),
  unname(expected_gene_set_contract)
)) {
  figure7_stop("Generated state-pathway gene-set contract is incompatible")
}

data_table_files <- c(
  "state_pathway_interval_definition.tsv",
  "panel_7F_pathway_activity_plot_data.tsv",
  "panel_7F_selected_pathway_gsea.tsv",
  "panel_7F_leading_edge_genes.tsv",
  "state_pathway_gene_ranking_complete.tsv",
  "state_pathway_gsea_complete.tsv",
  "state_pathway_sample_bin_coverage.tsv",
  "state_pathway_design_qc.tsv"
)
data_table_hashes <- vapply(
  file.path(staging_dir, data_table_files),
  figure7_sha256,
  character(1L)
)
names(data_table_hashes) <- paste0(
  sub("[.]tsv$", "", data_table_files),
  "_sha256"
)
gsea_retry_usage <- table(as.integer(gsea$nPermSimple))
gsea_retry_usage_value <- paste(
  names(gsea_retry_usage),
  as.integer(gsea_retry_usage),
  sep = ":",
  collapse = ","
)
provenance <- c(
  reference_kind = as.character(config$state_pathways$generated_reference_kind),
  canonical_publication_allowed = "false",
  source_results_id = "generated_pseudotime_state_pathways",
  source_code_revision = recorded_source_code_revision,
  exporter_script_sha256 =
    figure7_sha256(file.path(script_dir, "export_state_pathway_reference.R")),
  exporter_common_io_sha256 =
    figure7_sha256(file.path(script_dir, "src", "common_io.R")),
  seurat_rds_sha256 = input_value("seurat_rds", "sha256"),
  cellcycle_metadata_sha256 = input_value("cell_metadata", "sha256"),
  noncellcycle_metadata_sha256 = input_value("noncell_metadata", "sha256"),
  interval_definition_sha256 =
    observed_source_sha256[["interval_definition"]],
  interval_localization_grid_sha256 =
    observed_source_sha256[["interval_localization_grid"]],
  interval_localization_support_sha256 =
    observed_source_sha256[["interval_localization_support"]],
  interval_localization_test_sha256 =
    observed_source_sha256[["interval_localization_test"]],
  cell_expression_match_audit_sha256 =
    observed_source_sha256[["cell_expression_match_audit"]],
  localization_and_expression_n_cells = "2881",
  metadata_to_expression_match_rate = "1",
  interval_selection_rule =
    as.character(config$state_pathways$interval_selection_rule),
  interval_pointwise_p_min =
    format(min(pointwise_p[modeled_grid]), digits = 17),
  interval_pointwise_p_max =
    format(max(pointwise_p[modeled_grid]), digits = 17),
  assay = "RNA",
  counts_layer = "counts",
  nuisance_policy = "injected_initial_ploidy_only_no_endpoint_cn_score",
  initial_ploidy_levels = "2N,4N",
  endpoint_cn_score_covariate_prohibited = "true",
  spline_df = "5",
  pseudotime_bins = "20",
  minimum_cells_per_sample_bin = "5",
  grid_size = "501",
  seed = "1",
  gene_set_source = "pinned msigdbr runtime query",
  gene_set_release = gene_set_value("database_release"),
  gene_set_species = gene_set_value("species"),
  gene_set_collections = gene_set_value("collections"),
  gene_set_membership_sha256 = gene_set_value("membership_sha256"),
  gene_set_min_size = "15",
  gene_set_max_size = "500",
  expression_filter = "CPM > 1 in at least 14 retained sample-bin observations",
  normalization = "TMM",
  observation_model = "voom",
  mouse_block = "duplicateCorrelation",
  nuisance_terms = "dose_mg_factor,initial_ploidy_factor",
  treatment_by_pseudotime_interaction = "false",
  empirical_bayes = "robust",
  contrast = "mean(primary grid) - 0.5 * mean(left grid) - 0.5 * mean(right grid)",
  gsea_rank_statistic = "moderated_t",
  gsea_nperm_simple =
    as.character(config$state_pathways$gsea_nperm_simple),
  gsea_nperm_simple_max =
    as.character(config$state_pathways$gsea_nperm_simple_max),
  gsea_nperm_simple_multiplier =
    as.character(config$state_pathways$gsea_nperm_simple_multiplier),
  gsea_nperm_simple_usage = gsea_retry_usage_value,
  gsea_adaptive_retry_rule =
    "retry only unresolved pathways at geometric nPermSimple increments; merge by pathway; recompute collection-wide BH; fail closed at cap",
  pathway_activity = "row-standardize fitted expression per leading-edge gene, then average genes",
  feature_species_policy_id =
    as.character(config$feature_species$policy_id),
  feature_species_policy =
    as.character(config$state_pathways$generated_feature_species_policy),
  human_feature_prefix =
    as.character(config$feature_species$human_prefix),
  mouse_feature_prefix =
    as.character(config$feature_species$mouse_prefix),
  unknown_feature_policy =
    as.character(config$feature_species$unknown_feature_policy),
  n_input_features = as.character(species_counts[["input"]]),
  n_human_features_retained =
    as.character(species_counts[["human"]]),
  n_mouse_features_excluded =
    as.character(species_counts[["mouse"]]),
  n_ambiguous_features =
    as.character(species_counts[["ambiguous"]]),
  feature_species_audit_sha256 =
    observed_source_sha256[["feature_species_audit"]],
  feature_species_policy_code_sha256 =
    input_value("feature_species_policy_code", "sha256"),
  support_common_io_sha256 =
    input_value("common_io_code", "sha256"),
  density_localization_config_sha256 =
    input_value("density_localization_config", "sha256"),
  density_localization_code_sha256 =
    input_value("density_localization_code", "sha256"),
  gsea_ranking_rule = "one row per cleaned gene symbol; descending moderated t; ties retain cleaned-symbol order",
  pathway_selection_fdr_threshold =
    as.character(figure7_generated_pathway_fdr_threshold()),
  pathway_selection_rule =
    figure7_generated_pathway_selection_rule(),
  activity_table_sha256 =
    data_table_hashes[["panel_7F_pathway_activity_plot_data_sha256"]],
  data_table_hashes,
  figure7_config_sha256 = input_value("config", "sha256"),
  figure7_config_contract_sha256 =
    figure7_state_config_contract_sha256(config),
  generated_reference_id = reference_id,
  workflow_id = workflow_id,
  model_id = model_id,
  accumulated_interval = sprintf(
    "[%.3f,%.3f]",
    derived_intervals$primary_accumulated_state$start,
    derived_intervals$primary_accumulated_state$end
  ),
  left_neighbor_interval = sprintf(
    "[%.3f,%.3f)",
    derived_intervals$left_neighbor$start,
    derived_intervals$left_neighbor$end
  ),
  right_neighbor_interval = sprintf(
    "(%.3f,%.3f]",
    derived_intervals$right_neighbor$start,
    derived_intervals$right_neighbor$end
  ),
  recorded_cell_metadata_locator = input_value("cell_metadata", "locator"),
  recorded_noncell_metadata_locator = input_value("noncell_metadata", "locator"),
  recorded_seurat_rds_locator = input_value("seurat_rds", "locator"),
  msigdbr_package_version = package_value("msigdbr"),
  fgsea_package_version = package_value("fgsea"),
  source_activity_sha256 = observed_source_sha256[["activity"]],
  source_primary_gsea_sha256 = observed_source_sha256[["primary_gsea"]],
  source_leading_edge_sha256 = observed_source_sha256[["leading_edge"]],
  source_gene_contrast_sha256 = observed_source_sha256[["gene_contrast"]],
  source_gene_resolution_sha256 = observed_source_sha256[["gene_resolution"]],
  source_sample_bin_metadata_sha256 = observed_source_sha256[["sample_bin_metadata"]],
  source_design_audit_sha256 = observed_source_sha256[["design_audit"]],
  source_primary_coverage_sha256 = observed_source_sha256[["primary_coverage"]]
)
provenance_export <- data.frame(key = names(provenance), value = unname(provenance), stringsAsFactors = FALSE)
write_reference(provenance_export, "state_pathway_provenance.tsv")

reference_files <- c(
  "state_pathway_interval_definition.tsv",
  "panel_7F_pathway_activity_plot_data.tsv",
  "panel_7F_selected_pathway_gsea.tsv",
  "panel_7F_leading_edge_genes.tsv",
  "state_pathway_gene_ranking_complete.tsv",
  "state_pathway_gsea_complete.tsv",
  "state_pathway_sample_bin_coverage.tsv",
  "state_pathway_design_qc.tsv",
  "state_pathway_provenance.tsv"
)
expected_files <- reference_files
observed_files <- sort(list.files(staging_dir, all.files = FALSE, recursive = FALSE))
if (!identical(observed_files, sort(expected_files)) || any(file.info(file.path(staging_dir, expected_files))$size <= 0)) {
  figure7_stop(
    "Generated staging directory does not contain exactly nine reference ",
    "tables"
  )
}

dir.create(dirname(output_dir), recursive = TRUE, showWarnings = FALSE)
if (!file.rename(staging_dir, output_dir)) figure7_stop("Could not atomically materialize generated reference: ", output_dir)
completed <- TRUE
message(
  "Exported generated GRCh-only state-pathway panel-7F reference: ",
  output_dir
)
message("Activity rows: ", nrow(activity_export), "; selected pathways: ", nrow(selected_export))
