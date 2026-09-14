#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 280)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

format_value <- function(x) {
  if (length(x) == 0L) return("")
  paste(as.character(x), collapse = ";")
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05s_validate_latest_composition_association.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
p_dir <- file.path(results_root, "05p_latest_human_mouse_sample_matrix")
q_dir <- file.path(results_root, "05q_mouse_celltype_composition_2N_vs_4N")
r_dir <- file.path(results_root, "05r_human_cluster_mouse_celltype_association")
out_dir <- file.path(results_root, "05s_latest_composition_association_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

audit_rows <- list()
audit_index <- 0L
add_check <- function(check_id, passed, actual, expected, severity = "ERROR", note = "") {
  audit_index <<- audit_index + 1L
  audit_rows[[audit_index]] <<- data.frame(
    check_id = check_id,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    severity = severity,
    actual = format_value(actual),
    expected = format_value(expected),
    note = note,
    stringsAsFactors = FALSE
  )
}

required_files <- c(
  file.path(p_dir, "completion.txt"),
  file.path(p_dir, "latest_human_mouse_sample_matrix.rds"),
  file.path(p_dir, "input_manifest_sha256.csv"),
  file.path(q_dir, "completion.txt"),
  file.path(q_dir, "mouse_cell_type_2N_vs_4N_tests.csv"),
  file.path(q_dir, "mouse_fine_category_2N_vs_4N_tests.csv"),
  file.path(q_dir, "global_mouse_composition_2N_vs_4N_tests.csv"),
  file.path(r_dir, "completion.txt"),
  file.path(r_dir, "human_cluster_mouse_category_associations.csv"),
  file.path(r_dir, "human_cluster_mouse_category_ploidy_interactions.csv")
)
missing_required <- required_files[!file.exists(required_files)]
add_check("required_inputs_exist", length(missing_required) == 0L, basename(missing_required), "all required 05p-05r outputs")
if (length(missing_required) > 0L) {
  audit <- do.call(rbind, audit_rows)
  utils::write.csv(audit, file.path(out_dir, "validation_audit.csv"), row.names = FALSE)
  stop("Missing required 05p-05r outputs: ", paste(missing_required, collapse = ", "), call. = FALSE)
}

dat <- readRDS(file.path(p_dir, "latest_human_mouse_sample_matrix.rds"))
sample_metadata <- dat$sample_metadata
sample_summary <- dat$sample_summary
human_counts <- dat$human_cluster_counts
primary_counts <- dat$mouse_cell_type_counts
fine_counts <- dat$mouse_fine_category_counts
crosswalk <- dat$mouse_cluster_annotation

expected_samples <- names(unlist(cfg$sample_sources, use.names = TRUE))
excluded_samples <- as.character(unlist(cfg$immune_composition$excluded_samples, use.names = FALSE))
add_check("sample_set", setequal(sample_metadata$sample, expected_samples) && nrow(sample_metadata) == 16L, paste(sort(sample_metadata$sample), collapse = ","), "16 configured tissue samples")
add_check("cell_culture_exclusion", !any(sample_metadata$sample %in% excluded_samples), intersect(sample_metadata$sample, excluded_samples), "no excluded cell-culture samples")
ploidy_counts <- table(factor(as.character(sample_metadata$Ploidy), levels = c("2N", "4N")))
add_check("ploidy_replicates", identical(as.integer(ploidy_counts), c(8L, 8L)), as.integer(ploidy_counts), c(8L, 8L))
dose_ploidy <- table(factor(as.character(sample_metadata$Dose), levels = c("0", "30", "120")), factor(as.character(sample_metadata$Ploidy), levels = c("2N", "4N")))
add_check("dose_by_ploidy_replicates", identical(as.integer(dose_ploidy), c(4L, 2L, 2L, 4L, 2L, 2L)), as.integer(dose_ploidy), c(4L, 2L, 2L, 4L, 2L, 2L))
add_check("sample_is_statistical_unit", identical(dat$contract$biological_unit, "sample"), dat$contract$biological_unit, "sample")
add_check("cluster_level_annotation", identical(dat$contract$annotation_level, "existing 05d cluster only") && identical(dat$contract$cell_level_classifier_used, FALSE), paste(dat$contract$annotation_level, dat$contract$cell_level_classifier_used), "existing-cluster labels; no cell-level classifier")
add_check("no_joint_expression_analysis", identical(dat$contract$joint_expression_integration_or_clustering, FALSE), dat$contract$joint_expression_integration_or_clustering, FALSE)
add_check("mouse_labels_complete", !anyNA(primary_counts$mouse_cell_type_latest) && !any(primary_counts$mouse_cell_type_latest == ""), sum(is.na(primary_counts$mouse_cell_type_latest) | primary_counts$mouse_cell_type_latest == ""), 0)
add_check("fine_labels_complete", !anyNA(fine_counts$mouse_cell_type_state_latest) && !any(fine_counts$mouse_cell_type_state_latest == ""), sum(is.na(fine_counts$mouse_cell_type_state_latest) | fine_counts$mouse_cell_type_state_latest == ""), 0)
add_check("crosswalk_one_row_per_cluster", !anyDuplicated(as.character(crosswalk$seurat_clusters)) && nrow(crosswalk) == 20L, nrow(crosswalk), 20L)

fraction_sums <- list(
  human = stats::aggregate(fraction_within_human_tumor ~ sample, human_counts, sum)$fraction_within_human_tumor,
  mouse_primary = stats::aggregate(fraction_within_mouse ~ sample, primary_counts, sum)$fraction_within_mouse,
  mouse_fine = stats::aggregate(fraction_within_mouse ~ sample, fine_counts, sum)$fraction_within_mouse
)
for (name_i in names(fraction_sums)) {
  max_error <- max(abs(fraction_sums[[name_i]] - 1))
  add_check(paste0(name_i, "_fractions_sum_to_one"), is.finite(max_error) && max_error < 1e-10, signif(max_error, 5), "<1e-10")
}
mouse_total_from_summary <- sum(sample_summary$sample_total_mouse_cells)
mouse_total_from_counts <- sum(primary_counts$n_cells)
add_check("mouse_cell_total_consistent", identical(as.numeric(mouse_total_from_summary), as.numeric(mouse_total_from_counts)) && mouse_total_from_counts > 0, c(mouse_total_from_summary, mouse_total_from_counts), "equal positive totals")
human_total_from_summary <- sum(sample_summary$sample_total_human_tumor_cells)
human_total_from_counts <- sum(human_counts$n_cells)
add_check("human_cell_total_consistent", identical(as.numeric(human_total_from_summary), as.numeric(human_total_from_counts)) && human_total_from_counts > 0, c(human_total_from_summary, human_total_from_counts), "equal positive totals")

p_input_manifest <- utils::read.csv(file.path(p_dir, "input_manifest_sha256.csv"), check.names = FALSE)
add_check("latest_mouse_annotation_source", any(p_input_manifest$role == "latest_05h_cluster_level_mouse_annotation") && any(grepl("/05h_existing_immune_cluster_annotation/", p_input_manifest$path, fixed = TRUE)), paste(p_input_manifest$path, collapse = ";"), "05h existing-cluster annotation")
all_input_paths <- c(
  p_input_manifest$path,
  utils::read.csv(file.path(q_dir, "input_manifest_sha256.csv"), check.names = FALSE)$path,
  utils::read.csv(file.path(r_dir, "input_manifest_sha256.csv"), check.names = FALSE)$path
)
add_check("no_stale_05m_input", !any(grepl("/05m_", all_input_paths, fixed = TRUE)), all_input_paths[grepl("/05m_", all_input_paths, fixed = TRUE)], "no 05m input")

primary_tests <- utils::read.csv(file.path(q_dir, "mouse_cell_type_2N_vs_4N_tests.csv"), check.names = FALSE)
fine_tests <- utils::read.csv(file.path(q_dir, "mouse_fine_category_2N_vs_4N_tests.csv"), check.names = FALSE)
global_tests <- utils::read.csv(file.path(q_dir, "global_mouse_composition_2N_vs_4N_tests.csv"), check.names = FALSE)
for (test_name in c("primary_tests", "fine_tests")) {
  z <- get(test_name)
  add_check(paste0(test_name, "_sample_sizes"), all(z$n_samples_2N == 8L & z$n_samples_4N == 8L), paste(unique(paste(z$n_samples_2N, z$n_samples_4N, sep = "/")), collapse = ","), "8/8")
  add_check(paste0(test_name, "_exact_assignments"), all(z$n_exact_stratified_assignments == 2520L), unique(z$n_exact_stratified_assignments), 2520L)
  add_check(paste0(test_name, "_p_values_valid"), all(is.finite(z$p_value) & z$p_value >= 0 & z$p_value <= 1), range(z$p_value, na.rm = TRUE), "[0,1]")
  add_check(paste0(test_name, "_fdr_valid"), all(is.finite(z$p_value_bh) & z$p_value_bh >= 0 & z$p_value_bh <= 1), range(z$p_value_bh, na.rm = TRUE), "[0,1]")
}
add_check("global_tests_complete", nrow(global_tests) == 2L && all(global_tests$n_samples == 16L) && all(global_tests$n_exact_stratified_assignments == 2520L), paste(global_tests$n_samples, global_tests$n_exact_stratified_assignments, sep = "/"), "two resolutions; 16 samples; 2520 assignments")
add_check("global_p_values_valid", all(is.finite(global_tests$p_value) & global_tests$p_value >= 0 & global_tests$p_value <= 1), range(global_tests$p_value), "[0,1]")

associations <- utils::read.csv(file.path(r_dir, "human_cluster_mouse_category_associations.csv"), check.names = FALSE, na.strings = "NA")
interactions <- utils::read.csv(file.path(r_dir, "human_cluster_mouse_category_ploidy_interactions.csv"), check.names = FALSE, na.strings = "NA")
tested_assoc <- associations$test_status == "tested"
add_check("association_sample_sizes", all(associations$n_samples == 8L), unique(associations$n_samples), 8L)
add_check("association_exact_permutations", all(associations$n_exact_within_dose_permutations == 96L), unique(associations$n_exact_within_dose_permutations), 96L)
add_check("association_p_values_valid", any(tested_assoc) && all(is.finite(associations$p_value[tested_assoc]) & associations$p_value[tested_assoc] >= 0 & associations$p_value[tested_assoc] <= 1), range(associations$p_value[tested_assoc], na.rm = TRUE), "at least one tested pair; p in [0,1]")
add_check("association_p_values_exact_grid", all(abs(associations$p_value[tested_assoc] * 96 - round(associations$p_value[tested_assoc] * 96)) < 1e-8), max(abs(associations$p_value[tested_assoc] * 96 - round(associations$p_value[tested_assoc] * 96))), "<1e-8")
tested_assoc_fdr <- associations$p_value_bh[tested_assoc]
add_check("association_fdr_valid", all(is.finite(tested_assoc_fdr) & tested_assoc_fdr >= 0 & tested_assoc_fdr <= 1), range(tested_assoc_fdr, na.rm = TRUE), "[0,1]")
tested_interactions <- interactions$test_status == "tested_exploratory"
add_check("interaction_sample_sizes", all(interactions$n_samples == 16L), unique(interactions$n_samples), 16L)
add_check("interaction_p_values_valid", any(tested_interactions) && all(is.finite(interactions$interaction_p_value[tested_interactions]) & interactions$interaction_p_value[tested_interactions] >= 0 & interactions$interaction_p_value[tested_interactions] <= 1), range(interactions$interaction_p_value[tested_interactions], na.rm = TRUE), "at least one tested interaction; p in [0,1]")

plot_files <- c(
  file.path(q_dir, "mouse_cell_type_2N_vs_4N_grouped_bar.png"),
  file.path(q_dir, "mouse_cell_type_composition_2N_vs_4N_stacked.png"),
  file.path(q_dir, "mouse_fine_category_2N_vs_4N_sample_distribution.png"),
  file.path(r_dir, "heatmap_human_cluster_mouse_cell_type_composition_association.png"),
  file.path(r_dir, "heatmap_human_cluster_mouse_cell_type_recruitment_proxy_association.png"),
  file.path(r_dir, "heatmap_human_cluster_mouse_fine_category_association.png"),
  file.path(r_dir, "heatmap_primary_association_difference_4N_minus_2N.png"),
  file.path(r_dir, "top_human_cluster_mouse_cell_type_association_scatterplots.png")
)
plot_ok <- file.exists(plot_files) & file.info(plot_files)$size > 10000
add_check("plot_files_nonempty", all(plot_ok), basename(plot_files[!plot_ok]), "all eight PNG files >10 KB")

for (stage_dir in c(q_dir, r_dir)) {
  manifest_path <- file.path(stage_dir, "output_manifest_sha256.csv")
  manifest <- utils::read.csv(manifest_path, check.names = FALSE)
  paths <- file.path(stage_dir, manifest$file)
  exist_ok <- file.exists(paths)
  size_ok <- exist_ok & file.info(paths)$size == manifest$size_bytes
  sha_ok <- rep(FALSE, length(paths))
  sha_ok[exist_ok] <- vapply(paths[exist_ok], digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE) == manifest$sha256[exist_ok]
  stage_name <- basename(stage_dir)
  add_check(paste0(stage_name, "_manifest_files_exist"), all(exist_ok), basename(paths[!exist_ok]), "all manifested files")
  add_check(paste0(stage_name, "_manifest_sizes_match"), all(size_ok), basename(paths[!size_ok]), "recorded sizes")
  add_check(paste0(stage_name, "_manifest_sha256_match"), all(sha_ok), basename(paths[!sha_ok]), "recorded SHA256")
}

audit <- do.call(rbind, audit_rows)
utils::write.csv(audit, file.path(out_dir, "validation_audit.csv"), row.names = FALSE)
failed <- audit[audit$status == "FAIL" & audit$severity == "ERROR", , drop = FALSE]
if (nrow(failed) > 0L) {
  writeLines(c("VALIDATION FAILED", paste0(failed$check_id, ": ", failed$actual, " (expected ", failed$expected, ")")), file.path(out_dir, "validation_report.txt"))
  stop("05s validation failed: ", paste(failed$check_id, collapse = ", "), call. = FALSE)
}

readiness <- list(
  status = "ready_to_share_with_caveats",
  validated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  input_endpoint = "current 05h existing-cluster-level mouse annotation",
  statistical_unit = "biological sample",
  n_samples = 16L,
  n_samples_per_ploidy = 8L,
  dose_control = "exact restricted permutations within Dose strata",
  annotation_level = "existing cluster only",
  cell_level_classifier_used = FALSE,
  joint_human_mouse_expression_analysis = FALSE,
  caveats = c(
    "Eight samples per Ploidy limits power for cluster-by-cell-type associations.",
    "The 96 within-Dose permutations per Ploidy limit p-value resolution.",
    "Captured scRNA-seq fractions depend on dissociation and recovery and are not absolute histologic abundances.",
    "Associations may be consistent with differential recruitment but do not establish causal recruitment."
  )
)
yaml::write_yaml(readiness, file.path(out_dir, "analysis_readiness.yaml"))
writeLines(c(
  "05s validation passed.",
  paste0("Checks passed: ", nrow(audit), "; failed: 0."),
  paste0("Mouse cells represented in sample matrices: ", mouse_total_from_counts, "."),
  paste0("Human tumor cells represented in sample matrices: ", human_total_from_counts, "."),
  paste0("Primary mouse cell types: ", nrow(primary_tests), "; fine categories: ", nrow(fine_tests), "."),
  paste0("Human tumor clusters: ", length(unique(human_counts$human_cluster)), "."),
  paste0("Primary composition BH-FDR<0.05: ", sum(primary_tests$p_value_bh < 0.05), "."),
  paste0("Fine composition BH-FDR<0.05: ", sum(fine_tests$p_value_bh < 0.05), "."),
  paste0("Primary association BH-FDR<0.05: ", sum(associations$resolution == "primary_cell_type" & associations$p_value_bh < 0.05, na.rm = TRUE), "."),
  paste0("Primary exploratory Ploidy-interaction BH-FDR<0.05: ", sum(interactions$resolution == "primary_cell_type" & interactions$interaction_p_value_bh < 0.05, na.rm = TRUE), "."),
  "Status: ready_to_share_with_caveats."
), file.path(out_dir, "validation_report.txt"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05s.txt"))

input_files <- c(
  file.path(p_dir, "latest_human_mouse_sample_matrix.rds"),
  file.path(q_dir, "completion.txt"),
  file.path(r_dir, "completion.txt"),
  normalizePath(config_path, mustWork = TRUE), script_path
)
input_manifest <- data.frame(
  role = c("05p_current_sample_matrix", "05q_completion", "05r_completion", "analysis_config", "validation_script"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)
writeLines("05s latest composition and association validation completed successfully.", file.path(out_dir, "completion.txt"))
output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files), size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
message("[05s] Validation passed. Output: ", out_dir)
