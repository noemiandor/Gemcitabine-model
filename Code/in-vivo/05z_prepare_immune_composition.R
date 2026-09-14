#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

stop_if_missing_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
}

complete_counts <- function(counts, samples, subtypes) {
  grid <- expand.grid(sample = samples, immune_subtype = subtypes, stringsAsFactors = FALSE)
  out <- merge(grid, counts, by = c("sample", "immune_subtype"), all.x = TRUE, sort = FALSE)
  out$n_cells[is.na(out$n_cells)] <- 0L
  out
}

percent_label <- function(x, accuracy = 1) {
  paste0(format(round(100 * x / accuracy) * accuracy, trim = TRUE), "%")
}

exact_two_sample_permutation <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) < 2L || length(y) < 2L) {
    stop("Each Ploidy group must contain at least two biological samples.", call. = FALSE)
  }
  pooled <- c(x, y)
  combinations <- utils::combn(seq_along(pooled), length(x))
  permuted_differences <- apply(combinations, 2L, function(idx) {
    mean(pooled[idx]) - mean(pooled[-idx])
  })
  observed_difference <- mean(x) - mean(y)
  tolerance <- sqrt(.Machine$double.eps)
  p_value <- mean(abs(permuted_differences) >= abs(observed_difference) - tolerance)
  list(
    observed_difference = observed_difference,
    p_value = p_value,
    n_permutations = ncol(combinations)
  )
}

p_value_stars <- function(p) {
  ifelse(
    !is.finite(p), "NA",
    ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
  )
}

script_dir <- resolve_script_dir()
project_root <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_csv <- file.path(results_root, "05m_myeloid_states_final", "mouse_cell_annotations_with_myeloid_states.csv")
census_csv <- file.path(results_root, "05n_unified_human_mouse_census", "unified_human_mouse_captured_cell_census.csv")
out_dir <- file.path(results_root, "05z_immune_composition")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_csv)) stop("Missing 05m final annotation table: ", input_csv, call. = FALSE)
if (!file.exists(census_csv)) stop("Missing 05n unified census: ", census_csv, call. = FALSE)

# The sample-to-Ploidy/Dose mapping is an explicit transcription of the same
# sample_info.xlsx used by the human-cell workflow. Its hash is checked before
# analysis so Dose is never inferred from animal or harvest suffixes.
metadata_source <- file.path(project_root, as.character(cfg$immune_composition$sample_metadata_source))
if (!file.exists(metadata_source)) stop("Missing authoritative sample metadata workbook: ", metadata_source, call. = FALSE)
metadata_source_sha256 <- digest::digest(metadata_source, algo = "sha256", file = TRUE, serialize = FALSE)
expected_metadata_sha256 <- as.character(cfg$immune_composition$sample_metadata_source_sha256)
if (!identical(metadata_source_sha256, expected_metadata_sha256)) {
  stop("sample_info.xlsx SHA256 mismatch. Expected ", expected_metadata_sha256, ", observed ", metadata_source_sha256, call. = FALSE)
}

metadata_cfg <- cfg$immune_composition$sample_metadata
if (is.null(metadata_cfg) || is.null(names(metadata_cfg)) || any(!nzchar(names(metadata_cfg)))) {
  stop("immune_composition.sample_metadata must be a named mapping.", call. = FALSE)
}
sample_metadata <- do.call(rbind, lapply(names(metadata_cfg), function(sample_id) {
  entry <- metadata_cfg[[sample_id]]
  data.frame(
    sample = sample_id,
    Ploidy = as.character(entry$Ploidy),
    Dose = as.character(entry$Dose),
    stringsAsFactors = FALSE
  )
}))
rownames(sample_metadata) <- NULL
if (anyDuplicated(sample_metadata$sample)) stop("Duplicated sample identifiers in sample metadata mapping.", call. = FALSE)

expected_samples <- names(cfg$sample_sources)
if (!setequal(sample_metadata$sample, expected_samples)) {
  stop(
    "Sample metadata keys must exactly match sample_sources. Missing: ",
    paste(setdiff(expected_samples, sample_metadata$sample), collapse = ", "),
    "; unexpected: ", paste(setdiff(sample_metadata$sample, expected_samples), collapse = ", "),
    call. = FALSE
  )
}
ploidy_order <- as.character(unlist(cfg$immune_composition$ploidy_order, use.names = FALSE))
dose_order <- as.character(unlist(cfg$immune_composition$dose_order, use.names = FALSE))
if (!setequal(unique(sample_metadata$Ploidy), ploidy_order)) stop("Ploidy groups do not match configured ploidy_order.", call. = FALSE)
if (!setequal(unique(sample_metadata$Dose), dose_order)) stop("Dose groups do not match configured dose_order.", call. = FALSE)
sample_metadata$sample_order <- match(sample_metadata$sample, expected_samples)
sample_metadata$ploidy_order <- match(sample_metadata$Ploidy, ploidy_order)
sample_metadata$dose_order <- match(sample_metadata$Dose, dose_order)
sample_metadata <- sample_metadata[order(sample_metadata$sample_order), , drop = FALSE]

stratum_map <- unique(sample_metadata[c("Ploidy", "Dose")])
expected_strata <- expand.grid(Ploidy = ploidy_order, Dose = dose_order, stringsAsFactors = FALSE)
if (nrow(merge(expected_strata, stratum_map, by = c("Ploidy", "Dose"))) != nrow(expected_strata)) {
  stop("At least one configured Dose x Ploidy combination has no biological sample.", call. = FALSE)
}

annotations <- utils::read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
stop_if_missing_columns(
  annotations,
  c("cell", "sample", "mouse_final_cell_type_v2", "mouse_final_immune_type", "mouse_neutrophil_state_final", "mouse_final_cell_type_state"),
  "05m final annotation table"
)
if (anyDuplicated(annotations$cell)) stop("The 05m final annotation table contains duplicated cell identifiers.", call. = FALSE)
if (any(is.na(annotations$sample) | !nzchar(annotations$sample))) stop("One or more cells have a missing sample identifier.", call. = FALSE)

excluded_samples <- unlist(cfg$immune_composition$excluded_samples, use.names = FALSE)
excluded_found <- intersect(unique(annotations$sample), excluded_samples)
if (length(excluded_found) > 0L) {
  stop("Excluded cell-line samples are present in the final mouse-cell table: ", paste(excluded_found, collapse = ", "), call. = FALSE)
}
unexpected_samples <- setdiff(unique(annotations$sample), expected_samples)
missing_samples <- setdiff(expected_samples, unique(annotations$sample))
if (length(unexpected_samples) > 0L) stop("Unexpected samples in 05m final table: ", paste(unexpected_samples, collapse = ", "), call. = FALSE)
if (length(missing_samples) > 0L) stop("Expected samples missing from 05m final table: ", paste(missing_samples, collapse = ", "), call. = FALSE)

label_prefix <- as.character(cfg$immune_composition$final_label_prefix)
immune_pattern <- paste0("^", gsub("([][{}()+*^$|\\?.])", "\\\\\\1", label_prefix), "\\s*")
is_immune <- !is.na(annotations$mouse_final_immune_type) & nzchar(annotations$mouse_final_immune_type)
immune <- annotations[is_immune, , drop = FALSE]
if (nrow(immune) == 0L) stop("No cells have a non-missing 05m final immune type.", call. = FALSE)

immune$immune_subtype <- trimws(as.character(immune$mouse_final_immune_type))
immune$immune_subtype[is.na(immune$immune_subtype) | !nzchar(immune$immune_subtype)] <- "Ambiguous"
configured_subtypes <- unlist(cfg$immune_composition$subtype_order, use.names = FALSE)
unexpected_subtypes <- setdiff(unique(immune$immune_subtype), configured_subtypes)
if (length(unexpected_subtypes) > 0L) {
  stop("Final immune annotations not represented in immune_composition.subtype_order: ", paste(unexpected_subtypes, collapse = ", "), call. = FALSE)
}
observed_subtypes <- configured_subtypes[configured_subtypes %in% unique(immune$immune_subtype)]

metadata_idx <- match(immune$sample, sample_metadata$sample)
immune$Ploidy <- sample_metadata$Ploidy[metadata_idx]
immune$Dose <- sample_metadata$Dose[metadata_idx]
if (anyNA(immune$Ploidy) || anyNA(immune$Dose)) stop("One or more immune cells could not be mapped to Ploidy and Dose.", call. = FALSE)

raw_counts <- stats::aggregate(
  list(n_cells = rep.int(1L, nrow(immune))),
  by = list(sample = immune$sample, immune_subtype = immune$immune_subtype),
  FUN = sum
)
sample_counts <- complete_counts(raw_counts, expected_samples, observed_subtypes)
metadata_idx <- match(sample_counts$sample, sample_metadata$sample)
sample_counts$Ploidy <- sample_metadata$Ploidy[metadata_idx]
sample_counts$Dose <- sample_metadata$Dose[metadata_idx]
sample_counts$sample_order <- sample_metadata$sample_order[metadata_idx]
sample_counts$ploidy_order <- sample_metadata$ploidy_order[metadata_idx]
sample_counts$dose_order <- sample_metadata$dose_order[metadata_idx]
sample_totals <- stats::aggregate(n_cells ~ sample, sample_counts, sum)
names(sample_totals)[2] <- "sample_total_immune_cells"
sample_counts <- merge(sample_counts, sample_totals, by = "sample", all.x = TRUE, sort = FALSE)
if (any(sample_counts$sample_total_immune_cells <= 0L)) stop("At least one expected sample has no final immune cells.", call. = FALSE)
sample_counts$fraction <- sample_counts$n_cells / sample_counts$sample_total_immune_cells
sample_counts$immune_subtype_order <- match(sample_counts$immune_subtype, observed_subtypes)
sample_counts <- sample_counts[order(sample_counts$ploidy_order, sample_counts$dose_order, sample_counts$sample_order, sample_counts$immune_subtype_order), , drop = FALSE]

sample_fraction_sums <- stats::aggregate(fraction ~ sample, sample_counts, sum)
if (any(abs(sample_fraction_sums$fraction - 1) > 1e-12)) stop("Per-sample immune fractions do not sum to 1.", call. = FALSE)

stratum_pooled <- stats::aggregate(n_cells ~ Ploidy + Dose + immune_subtype, sample_counts, sum)
stratum_totals <- stats::aggregate(n_cells ~ Ploidy + Dose, stratum_pooled, sum)
names(stratum_totals)[3] <- "stratum_total_immune_cells"
stratum_pooled <- merge(stratum_pooled, stratum_totals, by = c("Ploidy", "Dose"), all.x = TRUE, sort = FALSE)
stratum_pooled$pooled_cell_fraction <- stratum_pooled$n_cells / stratum_pooled$stratum_total_immune_cells
stratum_pooled$ploidy_order <- match(stratum_pooled$Ploidy, ploidy_order)
stratum_pooled$dose_order <- match(stratum_pooled$Dose, dose_order)
stratum_pooled$immune_subtype_order <- match(stratum_pooled$immune_subtype, observed_subtypes)
stratum_pooled <- stratum_pooled[order(stratum_pooled$ploidy_order, stratum_pooled$dose_order, stratum_pooled$immune_subtype_order), , drop = FALSE]

stratum_mean <- stats::aggregate(fraction ~ Ploidy + Dose + immune_subtype, sample_counts, mean)
names(stratum_mean)[4] <- "mean_sample_fraction"
stratum_sd <- stats::aggregate(fraction ~ Ploidy + Dose + immune_subtype, sample_counts, stats::sd)
names(stratum_sd)[4] <- "sd_sample_fraction"
stratum_n <- stats::aggregate(sample ~ Ploidy + Dose + immune_subtype, sample_counts, function(x) length(unique(x)))
names(stratum_n)[4] <- "n_samples"
stratum_equal <- Reduce(
  function(x, y) merge(x, y, by = c("Ploidy", "Dose", "immune_subtype"), all = TRUE, sort = FALSE),
  list(stratum_mean, stratum_sd, stratum_n)
)
stratum_equal$sem_sample_fraction <- stratum_equal$sd_sample_fraction / sqrt(stratum_equal$n_samples)
stratum_equal$ploidy_order <- match(stratum_equal$Ploidy, ploidy_order)
stratum_equal$dose_order <- match(stratum_equal$Dose, dose_order)
stratum_equal$immune_subtype_order <- match(stratum_equal$immune_subtype, observed_subtypes)
stratum_equal <- stratum_equal[order(stratum_equal$ploidy_order, stratum_equal$dose_order, stratum_equal$immune_subtype_order), , drop = FALSE]
stratum_mean_sums <- stats::aggregate(mean_sample_fraction ~ Ploidy + Dose, stratum_equal, sum)
if (any(abs(stratum_mean_sums$mean_sample_fraction - 1) > 1e-12)) stop("Equal-weight Dose x Ploidy mean fractions do not sum to 1.", call. = FALSE)

if (!identical(ploidy_order, c("2N", "4N"))) {
  stop("The exact permutation comparison requires ploidy_order to be 2N, 4N.", call. = FALSE)
}
test_grid <- expand.grid(Dose = dose_order, immune_subtype = observed_subtypes, stringsAsFactors = FALSE)
test_rows <- lapply(seq_len(nrow(test_grid)), function(i) {
  dose_i <- test_grid$Dose[[i]]
  subtype_i <- test_grid$immune_subtype[[i]]
  dat_i <- sample_counts[
    as.character(sample_counts$Dose) == dose_i & sample_counts$immune_subtype == subtype_i,
    , drop = FALSE
  ]
  x <- dat_i$fraction[dat_i$Ploidy == "2N"]
  y <- dat_i$fraction[dat_i$Ploidy == "4N"]
  test_i <- exact_two_sample_permutation(x, y)
  data.frame(
    Dose = dose_i,
    immune_subtype = subtype_i,
    n_samples_2N = length(x),
    n_samples_4N = length(y),
    mean_fraction_2N = mean(x),
    mean_fraction_4N = mean(y),
    difference_2N_minus_4N = test_i$observed_difference,
    test = "two-sided exact permutation test of mean sample fractions",
    n_permutations = test_i$n_permutations,
    p_value = test_i$p_value,
    stringsAsFactors = FALSE
  )
})
ploidy_tests <- do.call(rbind, test_rows)
ploidy_tests$p_value_bh_within_dose <- NA_real_
for (dose_i in dose_order) {
  idx <- which(as.character(ploidy_tests$Dose) == dose_i)
  ploidy_tests$p_value_bh_within_dose[idx] <- stats::p.adjust(ploidy_tests$p_value[idx], method = "BH")
}
ploidy_tests$significance_fdr <- p_value_stars(ploidy_tests$p_value_bh_within_dose)
ploidy_tests$dose_order <- match(as.character(ploidy_tests$Dose), dose_order)
ploidy_tests$immune_subtype_order <- match(ploidy_tests$immune_subtype, observed_subtypes)
ploidy_tests <- ploidy_tests[order(ploidy_tests$dose_order, ploidy_tests$immune_subtype_order), , drop = FALSE]

palette_all <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
if (is.null(names(palette_all)) || any(!nzchar(names(palette_all)))) {
  stop("immune_composition.subtype_colors must be a named mapping.", call. = FALSE)
}
if (!all(configured_subtypes %in% names(palette_all))) {
  stop("The fixed palette does not cover all configured immune subtypes.", call. = FALSE)
}
palette_df <- data.frame(
  immune_subtype = observed_subtypes,
  color_hex = unname(palette_all[observed_subtypes]),
  immune_subtype_order = seq_along(observed_subtypes),
  stringsAsFactors = FALSE
)

sample_map <- merge(sample_metadata, sample_totals, by = "sample", all.x = TRUE, sort = FALSE)
sample_map <- sample_map[order(sample_map$sample_order), , drop = FALSE]
immune_metadata <- immune[c("cell", "sample", "mouse_final_cell_type_v2", "mouse_final_cell_type_state", "immune_subtype", "mouse_neutrophil_state_final", "Ploidy", "Dose")]
immune_metadata$sample_order <- match(immune_metadata$sample, expected_samples)

utils::write.csv(sample_counts, file.path(out_dir, "immune_composition_by_sample.csv"), row.names = FALSE)
utils::write.csv(sample_map, file.path(out_dir, "immune_sample_ploidy_dose_map.csv"), row.names = FALSE)
utils::write.csv(stratum_equal, file.path(out_dir, "immune_composition_by_dose_ploidy_equal_sample_weight.csv"), row.names = FALSE)
utils::write.csv(stratum_pooled, file.path(out_dir, "immune_composition_by_dose_ploidy_pooled_cells.csv"), row.names = FALSE)
utils::write.csv(ploidy_tests, file.path(out_dir, "immune_subtype_ploidy_tests_by_dose.csv"), row.names = FALSE)
utils::write.csv(palette_df, file.path(out_dir, "immune_subtype_palette.csv"), row.names = FALSE)
metadata_con <- gzfile(file.path(out_dir, "final_immune_cell_metadata.csv.gz"), open = "wt")
utils::write.csv(immune_metadata, metadata_con, row.names = FALSE, na = "NA")
close(metadata_con)

# Final Neutrophil-state composition. The denominator is all final v2
# Neutrophils within each biological sample. Colors are identical to the 05m
# final state UMAP.
neutrophils <- immune[immune$immune_subtype == "Neutrophil", , drop = FALSE]
if (nrow(neutrophils) == 0L) stop("05m final annotation contains no v2 Neutrophils.", call. = FALSE)
if (any(is.na(neutrophils$mouse_neutrophil_state_final) | !nzchar(neutrophils$mouse_neutrophil_state_final))) {
  stop("At least one final v2 Neutrophil lacks mouse_neutrophil_state_final.", call. = FALSE)
}
configured_state_order <- as.character(unlist(cfg$myeloid$state_order, use.names = FALSE))
unexpected_states <- setdiff(unique(neutrophils$mouse_neutrophil_state_final), configured_state_order)
if (length(unexpected_states) > 0L) stop("Unexpected final Neutrophil states: ", paste(unexpected_states, collapse = ", "), call. = FALSE)
observed_states <- configured_state_order[configured_state_order %in% unique(neutrophils$mouse_neutrophil_state_final)]
state_palette_all <- unlist(cfg$myeloid$state_colors, use.names = TRUE)
if (!all(configured_state_order %in% names(state_palette_all))) stop("myeloid.state_colors does not cover state_order.", call. = FALSE)

state_raw_counts <- stats::aggregate(
  list(n_cells = rep.int(1L, nrow(neutrophils))),
  by = list(sample = neutrophils$sample, neutrophil_state = neutrophils$mouse_neutrophil_state_final),
  FUN = sum
)
state_grid <- expand.grid(sample = expected_samples, neutrophil_state = observed_states, stringsAsFactors = FALSE)
state_sample_counts <- merge(state_grid, state_raw_counts, by = c("sample", "neutrophil_state"), all.x = TRUE, sort = FALSE)
state_sample_counts$n_cells[is.na(state_sample_counts$n_cells)] <- 0L
state_meta_idx <- match(state_sample_counts$sample, sample_metadata$sample)
state_sample_counts$Ploidy <- sample_metadata$Ploidy[state_meta_idx]
state_sample_counts$Dose <- sample_metadata$Dose[state_meta_idx]
state_sample_counts$sample_order <- sample_metadata$sample_order[state_meta_idx]
state_sample_counts$state_order <- match(state_sample_counts$neutrophil_state, observed_states)
state_totals <- stats::aggregate(n_cells ~ sample, state_sample_counts, sum)
names(state_totals)[2] <- "sample_total_neutrophils"
state_sample_counts <- merge(state_sample_counts, state_totals, by = "sample", all.x = TRUE, sort = FALSE)
if (any(state_sample_counts$sample_total_neutrophils <= 0L)) stop("At least one sample has no final v2 Neutrophils.", call. = FALSE)
state_sample_counts$fraction <- state_sample_counts$n_cells / state_sample_counts$sample_total_neutrophils
state_sample_counts <- state_sample_counts[order(state_sample_counts$sample_order, state_sample_counts$state_order), , drop = FALSE]
state_fraction_sums <- stats::aggregate(fraction ~ sample, state_sample_counts, sum)
if (any(abs(state_fraction_sums$fraction - 1) > 1e-12)) stop("Per-sample Neutrophil-state fractions do not sum to 1.", call. = FALSE)

state_equal_mean <- stats::aggregate(fraction ~ Ploidy + Dose + neutrophil_state, state_sample_counts, mean)
names(state_equal_mean)[4] <- "mean_sample_fraction"
state_equal_sd <- stats::aggregate(fraction ~ Ploidy + Dose + neutrophil_state, state_sample_counts, stats::sd)
names(state_equal_sd)[4] <- "sd_sample_fraction"
state_equal_n <- stats::aggregate(sample ~ Ploidy + Dose + neutrophil_state, state_sample_counts, function(x) length(unique(x)))
names(state_equal_n)[4] <- "n_samples"
state_equal <- Reduce(
  function(x, y) merge(x, y, by = c("Ploidy", "Dose", "neutrophil_state"), all = TRUE, sort = FALSE),
  list(state_equal_mean, state_equal_sd, state_equal_n)
)
state_equal$sem_sample_fraction <- state_equal$sd_sample_fraction / sqrt(state_equal$n_samples)
state_equal$state_order <- match(state_equal$neutrophil_state, observed_states)
state_equal$ploidy_order <- match(state_equal$Ploidy, ploidy_order)
state_equal$dose_order <- match(state_equal$Dose, dose_order)
state_equal <- state_equal[order(state_equal$ploidy_order, state_equal$dose_order, state_equal$state_order), , drop = FALSE]

state_test_grid <- expand.grid(Dose = dose_order, neutrophil_state = observed_states, stringsAsFactors = FALSE)
state_test_rows <- lapply(seq_len(nrow(state_test_grid)), function(i) {
  dose_i <- state_test_grid$Dose[[i]]
  state_i <- state_test_grid$neutrophil_state[[i]]
  dat_i <- state_sample_counts[as.character(state_sample_counts$Dose) == dose_i & state_sample_counts$neutrophil_state == state_i, , drop = FALSE]
  x <- dat_i$fraction[dat_i$Ploidy == "2N"]
  y <- dat_i$fraction[dat_i$Ploidy == "4N"]
  test_i <- exact_two_sample_permutation(x, y)
  data.frame(
    Dose = dose_i, neutrophil_state = state_i, n_samples_2N = length(x), n_samples_4N = length(y),
    mean_fraction_2N = mean(x), mean_fraction_4N = mean(y), difference_2N_minus_4N = test_i$observed_difference,
    test = "two-sided exact permutation test of mean sample fractions", n_permutations = test_i$n_permutations,
    p_value = test_i$p_value, stringsAsFactors = FALSE
  )
})
state_ploidy_tests <- do.call(rbind, state_test_rows)
state_ploidy_tests$p_value_bh_within_dose <- NA_real_
for (dose_i in dose_order) {
  idx <- which(as.character(state_ploidy_tests$Dose) == dose_i)
  state_ploidy_tests$p_value_bh_within_dose[idx] <- stats::p.adjust(state_ploidy_tests$p_value[idx], method = "BH")
}
state_ploidy_tests$significance_fdr <- p_value_stars(state_ploidy_tests$p_value_bh_within_dose)
state_ploidy_tests$state_order <- match(state_ploidy_tests$neutrophil_state, observed_states)
state_ploidy_tests$dose_order <- match(as.character(state_ploidy_tests$Dose), dose_order)
state_ploidy_tests <- state_ploidy_tests[order(state_ploidy_tests$dose_order, state_ploidy_tests$state_order), , drop = FALSE]
state_palette_df <- data.frame(
  neutrophil_state = observed_states, color_hex = unname(state_palette_all[observed_states]),
  state_order = seq_along(observed_states), stringsAsFactors = FALSE
)
utils::write.csv(state_sample_counts, file.path(out_dir, "neutrophil_state_composition_by_sample.csv"), row.names = FALSE)
utils::write.csv(state_equal, file.path(out_dir, "neutrophil_state_composition_by_dose_ploidy_equal_sample_weight.csv"), row.names = FALSE)
utils::write.csv(state_ploidy_tests, file.path(out_dir, "neutrophil_state_ploidy_tests_by_dose.csv"), row.names = FALSE)
utils::write.csv(state_palette_df, file.path(out_dir, "neutrophil_state_palette.csv"), row.names = FALSE)

# Unified human-mouse captured-cell denominators. No cross-species expression
# integration is performed here; 05n has already linked retained cells by
# sample plus raw Cell Ranger barcode.
census <- utils::read.csv(census_csv, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
stop_if_missing_columns(
  census,
  c("sample", "cell_key", "captured_species", "mouse_final_cell_type_v2", "mouse_final_cell_type_state"),
  "05n unified census"
)
if (anyDuplicated(census$cell_key)) stop("05n unified census contains duplicate sample-barcode keys.", call. = FALSE)
if (!setequal(unique(census$sample), expected_samples)) stop("05n unified census sample universe differs from the frozen 16-sample contract.", call. = FALSE)
if (!setequal(unique(census$captured_species), c("Human tumor", "Mouse host"))) stop("05n census must contain Human tumor and Mouse host categories.", call. = FALSE)
census_meta_idx <- match(census$sample, sample_metadata$sample)
census$Ploidy <- sample_metadata$Ploidy[census_meta_idx]
census$Dose <- sample_metadata$Dose[census_meta_idx]

species_grid <- expand.grid(sample = expected_samples, captured_species = c("Human tumor", "Mouse host"), stringsAsFactors = FALSE)
species_raw_counts <- stats::aggregate(list(n_cells = rep.int(1L, nrow(census))), by = list(sample = census$sample, captured_species = census$captured_species), FUN = sum)
species_counts <- merge(species_grid, species_raw_counts, by = c("sample", "captured_species"), all.x = TRUE, sort = FALSE)
species_counts$n_cells[is.na(species_counts$n_cells)] <- 0L
species_meta_idx <- match(species_counts$sample, sample_metadata$sample)
species_counts$Ploidy <- sample_metadata$Ploidy[species_meta_idx]
species_counts$Dose <- sample_metadata$Dose[species_meta_idx]
species_counts$sample_order <- sample_metadata$sample_order[species_meta_idx]
unified_totals <- stats::aggregate(n_cells ~ sample, species_counts, sum)
names(unified_totals)[2] <- "sample_total_human_plus_mouse_retained_singlets"
species_counts <- merge(species_counts, unified_totals, by = "sample", all.x = TRUE, sort = FALSE)
species_counts$fraction <- species_counts$n_cells / species_counts$sample_total_human_plus_mouse_retained_singlets

mouse_census <- census[census$captured_species == "Mouse host", , drop = FALSE]
if (any(is.na(mouse_census$mouse_final_cell_type_v2))) stop("A retained mouse census cell lacks mouse_final_cell_type_v2.", call. = FALSE)
mouse_type_levels <- sort(unique(mouse_census$mouse_final_cell_type_v2))
mouse_type_grid <- expand.grid(sample = expected_samples, mouse_final_cell_type_v2 = mouse_type_levels, stringsAsFactors = FALSE)
mouse_type_raw <- stats::aggregate(
  list(n_cells = rep.int(1L, nrow(mouse_census))),
  by = list(sample = mouse_census$sample, mouse_final_cell_type_v2 = mouse_census$mouse_final_cell_type_v2), FUN = sum
)
mouse_type_counts <- merge(mouse_type_grid, mouse_type_raw, by = c("sample", "mouse_final_cell_type_v2"), all.x = TRUE, sort = FALSE)
mouse_type_counts$n_cells[is.na(mouse_type_counts$n_cells)] <- 0L
mouse_type_meta_idx <- match(mouse_type_counts$sample, sample_metadata$sample)
mouse_type_counts$Ploidy <- sample_metadata$Ploidy[mouse_type_meta_idx]
mouse_type_counts$Dose <- sample_metadata$Dose[mouse_type_meta_idx]
mouse_type_counts$sample_order <- sample_metadata$sample_order[mouse_type_meta_idx]
mouse_totals <- stats::aggregate(n_cells ~ sample, mouse_type_counts, sum)
names(mouse_totals)[2] <- "sample_total_mouse_retained_singlets"
mouse_type_counts <- merge(mouse_type_counts, mouse_totals, by = "sample", all.x = TRUE, sort = FALSE)
mouse_type_counts <- merge(mouse_type_counts, unified_totals, by = "sample", all.x = TRUE, sort = FALSE)
mouse_type_counts$fraction_of_mouse_cells <- mouse_type_counts$n_cells / mouse_type_counts$sample_total_mouse_retained_singlets
mouse_type_counts$fraction_of_human_plus_mouse_retained_singlets <- mouse_type_counts$n_cells / mouse_type_counts$sample_total_human_plus_mouse_retained_singlets

unified_categories <- data.frame(
  sample = census$sample,
  captured_category = ifelse(census$captured_species == "Human tumor", "Human tumor", census$mouse_final_cell_type_v2),
  stringsAsFactors = FALSE
)
unified_category_levels <- c("Human tumor", mouse_type_levels)
unified_category_grid <- expand.grid(sample = expected_samples, captured_category = unified_category_levels, stringsAsFactors = FALSE)
unified_category_raw <- stats::aggregate(list(n_cells = rep.int(1L, nrow(unified_categories))), by = list(sample = unified_categories$sample, captured_category = unified_categories$captured_category), FUN = sum)
unified_category_counts <- merge(unified_category_grid, unified_category_raw, by = c("sample", "captured_category"), all.x = TRUE, sort = FALSE)
unified_category_counts$n_cells[is.na(unified_category_counts$n_cells)] <- 0L
unified_category_counts$sample_total_human_plus_mouse_retained_singlets <- unname(unified_totals$sample_total_human_plus_mouse_retained_singlets[match(unified_category_counts$sample, unified_totals$sample)])
unified_category_counts$fraction <- unified_category_counts$n_cells / unified_category_counts$sample_total_human_plus_mouse_retained_singlets
unified_category_counts$Ploidy <- sample_metadata$Ploidy[match(unified_category_counts$sample, sample_metadata$sample)]
unified_category_counts$Dose <- sample_metadata$Dose[match(unified_category_counts$sample, sample_metadata$sample)]
unified_category_counts$sample_order <- sample_metadata$sample_order[match(unified_category_counts$sample, sample_metadata$sample)]

utils::write.csv(species_counts, file.path(out_dir, "human_vs_mouse_captured_composition_by_sample.csv"), row.names = FALSE)
utils::write.csv(mouse_type_counts, file.path(out_dir, "mouse_cell_types_with_mouse_and_unified_denominators_by_sample.csv"), row.names = FALSE)
utils::write.csv(unified_category_counts, file.path(out_dir, "human_tumor_plus_mouse_cell_type_composition_by_sample.csv"), row.names = FALSE)

contract <- list(
  sources = list(
    final_mouse_cell_annotation = normalizePath(input_csv, mustWork = TRUE),
    unified_human_mouse_census = normalizePath(census_csv, mustWork = TRUE),
    sample_metadata = normalizePath(metadata_source, mustWork = TRUE),
    sample_metadata_sha256 = metadata_source_sha256
  ),
  immune_cell_rule = "05m mouse_final_immune_type is non-missing",
  sample_composition_denominator = "all cells with a final immune annotation within that sample",
  primary_dose_ploidy_summary = "arithmetic mean of per-sample immune-subtype fractions within each Dose x Ploidy stratum; every biological sample has equal weight",
  supplementary_dose_ploidy_summary = "pooled immune-cell counts within each Dose x Ploidy stratum",
  ploidy_levels = ploidy_order,
  dose_levels_mg_per_kg = dose_order,
  excluded_samples = excluded_samples,
  ambiguous_policy = "Immune: Ambiguous is retained as an explicit immune subtype",
  inferential_statistics = list(
    unit = "biological sample fraction; cells are not treated as independent replicates",
    comparison = "2N versus 4N separately for every observed immune subtype within each Dose",
    test = "two-sided exact permutation test of the difference in mean sample fractions",
    multiple_testing = "Benjamini-Hochberg correction across all observed immune-subtype tests within each Dose",
    significance_labels = "**** FDR<0.0001; *** FDR<0.001; ** FDR<0.01; * FDR<0.05; ns otherwise"
  ),
  immune_color_source = "immune_composition.subtype_colors in 05_mouse_cell_analysis_config.yaml; shared with the 05g/05m final UMAPs",
  neutrophil_state_color_source = "myeloid.state_colors in 05_mouse_cell_analysis_config.yaml; shared with the 05m final state UMAP",
  unified_denominator_interpretation = "scRNA-seq captured-cell fraction, not an absolute histologic tissue fraction"
)
yaml::write_yaml(contract, file.path(out_dir, "immune_composition_contract.yaml"))

input_files <- c(input_csv, census_csv, metadata_source)
input_info <- file.info(input_files)
input_manifest <- data.frame(
  role = c("05m final mouse-cell/myeloid-state annotations", "05n unified retained human-mouse census", "authoritative sample Ploidy/Dose metadata"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = input_info$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)

unique_sample_strata <- unique(sample_counts[c("sample", "Ploidy", "Dose")])
stratum_sample_counts <- stats::aggregate(sample ~ Ploidy + Dose, unique_sample_strata, length)
stratum_metrics <- paste0("samples_", stratum_sample_counts$Ploidy, "_dose_", stratum_sample_counts$Dose, "mg_per_kg")
audit <- data.frame(
  metric = c(
    "final_mouse_cells", "final_immune_cells", "immune_fraction_of_final_mouse_cells",
    "samples", "ploidy_groups", "dose_groups", "observed_immune_subtypes",
    "excluded_cell_line_samples_present", "sample_metadata_sha256_matches",
    "max_abs_sample_fraction_sum_error", "max_abs_dose_ploidy_mean_fraction_sum_error",
    "ploidy_tests", "ploidy_tests_fdr_lt_0.05",
    "final_neutrophils", "observed_neutrophil_states", "neutrophil_state_ploidy_tests",
    "neutrophil_state_tests_fdr_lt_0.05", "unified_human_mouse_denominator_cells",
    "unified_human_cells", "unified_mouse_cells",
    stratum_metrics
  ),
  value = c(
    nrow(annotations), nrow(immune), format(nrow(immune) / nrow(annotations), digits = 12),
    length(unique(immune$sample)), length(unique(immune$Ploidy)), length(unique(immune$Dose)),
    length(observed_subtypes), length(excluded_found), identical(metadata_source_sha256, expected_metadata_sha256),
    format(max(abs(sample_fraction_sums$fraction - 1)), scientific = TRUE),
    format(max(abs(stratum_mean_sums$mean_sample_fraction - 1)), scientific = TRUE),
    nrow(ploidy_tests), sum(ploidy_tests$p_value_bh_within_dose < 0.05),
    nrow(neutrophils), length(observed_states), nrow(state_ploidy_tests),
    sum(state_ploidy_tests$p_value_bh_within_dose < 0.05), nrow(census),
    sum(census$captured_species == "Human tumor"), sum(census$captured_species == "Mouse host"),
    stratum_sample_counts$sample
  ),
  stringsAsFactors = FALSE
)
utils::write.table(audit, file.path(out_dir, "immune_composition_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

# All descriptive immune-composition figures are generated in this 05z step.
# The sample-level table is the direct plotting source so the denominator and
# experimental mapping used in every panel remain auditable.
sample_levels <- expected_samples
subtype_levels <- observed_subtypes
palette_values <- stats::setNames(palette_df$color_hex, palette_df$immune_subtype)

sample_plot <- sample_counts
sample_plot$sample <- factor(sample_plot$sample, levels = sample_levels)
sample_plot$immune_subtype <- factor(sample_plot$immune_subtype, levels = subtype_levels)
sample_plot$Ploidy <- factor(sample_plot$Ploidy, levels = ploidy_order)
sample_plot$Dose <- factor(as.character(sample_plot$Dose), levels = dose_order)

equal_plot <- stratum_equal
equal_plot$Ploidy <- factor(equal_plot$Ploidy, levels = ploidy_order)
equal_plot$Dose <- factor(as.character(equal_plot$Dose), levels = dose_order)
equal_plot$immune_subtype <- factor(equal_plot$immune_subtype, levels = subtype_levels)

sample_total_map <- stats::setNames(sample_totals$sample_total_immune_cells, sample_totals$sample)
sample_labels <- stats::setNames(
  paste0(sample_levels, "\n(n=", format(sample_total_map[sample_levels], big.mark = ",", scientific = FALSE), ")"),
  sample_levels
)

n_by_stratum <- stats::aggregate(sample ~ Ploidy + Dose, unique_sample_strata, length)
names(n_by_stratum)[3] <- "n_samples"
n_by_stratum$Ploidy <- factor(n_by_stratum$Ploidy, levels = ploidy_order)
n_by_stratum$Dose <- factor(as.character(n_by_stratum$Dose), levels = dose_order)
n_by_stratum$label <- paste0("n=", n_by_stratum$n_samples)

base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom",
    legend.title = ggplot2::element_blank(),
    plot.title.position = "plot",
    plot.caption = ggplot2::element_text(hjust = 0)
  )

state_palette_values <- stats::setNames(state_palette_df$color_hex, state_palette_df$neutrophil_state)
state_sample_plot <- state_sample_counts
state_sample_plot$sample <- factor(state_sample_plot$sample, levels = sample_levels)
state_sample_plot$neutrophil_state <- factor(state_sample_plot$neutrophil_state, levels = observed_states)
state_sample_plot$Ploidy <- factor(state_sample_plot$Ploidy, levels = ploidy_order)
state_sample_plot$Dose <- factor(as.character(state_sample_plot$Dose), levels = dose_order)
state_equal_plot <- state_equal
state_equal_plot$neutrophil_state <- factor(state_equal_plot$neutrophil_state, levels = observed_states)
state_equal_plot$Ploidy <- factor(state_equal_plot$Ploidy, levels = ploidy_order)
state_equal_plot$Dose <- factor(as.character(state_equal_plot$Dose), levels = dose_order)

p_state_sample <- ggplot2::ggplot(state_sample_plot, ggplot2::aes(x = sample, y = fraction, fill = neutrophil_state)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"), Ploidy = function(x) paste0("Ploidy ", x))
  ) +
  ggplot2::scale_fill_manual(values = state_palette_values, breaks = observed_states, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_label, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Final Neutrophil-state composition by sample",
    subtitle = "Denominator = all final v2 Neutrophils within each sample",
    x = NULL, y = "Fraction of Neutrophils",
    caption = "State colors are identical to the 05m final UMAP. Cell-line samples are excluded."
  ) + base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1))

p_state_equal <- ggplot2::ggplot(state_equal_plot, ggplot2::aes(x = Dose, y = mean_sample_fraction, fill = neutrophil_state)) +
  ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  ggplot2::facet_grid(~Ploidy, labeller = ggplot2::labeller(Ploidy = function(x) paste0("Ploidy ", x))) +
  ggplot2::scale_fill_manual(values = state_palette_values, breaks = observed_states, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_label, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Final Neutrophil-state composition by Dose and Ploidy",
    subtitle = "Arithmetic mean of per-sample state fractions; every biological sample has equal weight",
    x = "Gemcitabine dose (mg/kg)", y = "Mean fraction of Neutrophils",
    caption = "Denominator = final v2 Neutrophils within each sample."
  ) + base_theme

state_test_plot <- state_ploidy_tests
state_test_plot$Dose <- factor(as.character(state_test_plot$Dose), levels = dose_order)
state_test_plot$neutrophil_state <- factor(state_test_plot$neutrophil_state, levels = observed_states)
state_test_plot$y_min <- mapply(function(dose_i, state_i) {
  min(state_sample_counts$fraction[as.character(state_sample_counts$Dose) == as.character(dose_i) & state_sample_counts$neutrophil_state == as.character(state_i)])
}, state_test_plot$Dose, state_test_plot$neutrophil_state)
state_test_plot$y_max <- mapply(function(dose_i, state_i) {
  max(state_sample_counts$fraction[as.character(state_sample_counts$Dose) == as.character(dose_i) & state_sample_counts$neutrophil_state == as.character(state_i)])
}, state_test_plot$Dose, state_test_plot$neutrophil_state)
state_test_plot$y_pad <- pmax((state_test_plot$y_max - state_test_plot$y_min) * 0.12, 0.0001)
state_test_plot$y_bracket <- state_test_plot$y_max + state_test_plot$y_pad
state_test_plot$y_label <- state_test_plot$y_max + state_test_plot$y_pad * 1.9
p_state_tests <- ggplot2::ggplot(
  state_sample_plot,
  ggplot2::aes(x = Ploidy, y = fraction, shape = Ploidy, color = neutrophil_state, fill = neutrophil_state)
) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1234), size = 2.4, stroke = 0.5, alpha = 0.9) +
  ggplot2::stat_summary(ggplot2::aes(group = Ploidy), fun = mean, geom = "point", shape = 95, size = 8) +
  ggplot2::geom_segment(data = state_test_plot, ggplot2::aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket), inherit.aes = FALSE, color = "#343A40", linewidth = 0.4) +
  ggplot2::geom_text(data = state_test_plot, ggplot2::aes(x = 1.5, y = y_label, label = significance_fdr), inherit.aes = FALSE, color = "#111827", fontface = "bold", size = 3.3) +
  ggplot2::facet_wrap(
    ggplot2::vars(Dose, neutrophil_state), ncol = length(observed_states), scales = "free_y",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"))
  ) +
  ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
  ggplot2::scale_color_manual(values = state_palette_values, breaks = observed_states, drop = FALSE) +
  ggplot2::scale_fill_manual(values = state_palette_values, breaks = observed_states, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = function(x) percent_label(x, accuracy = 0.001), expand = ggplot2::expansion(mult = c(0.04, 0.06))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Neutrophil-state fractions: 2N versus 4N within each Dose",
    subtitle = "Each point is one biological sample; stars use within-Dose BH-adjusted exact-permutation p-values",
    x = "Ploidy", y = "Fraction of Neutrophils",
    caption = "**** FDR<0.0001; *** <0.001; ** <0.01; * <0.05; ns otherwise."
  ) + base_theme + ggplot2::guides(shape = "none", color = "none", fill = "none")

species_plot <- species_counts
species_plot$sample <- factor(species_plot$sample, levels = sample_levels)
species_plot$Ploidy <- factor(species_plot$Ploidy, levels = ploidy_order)
species_plot$Dose <- factor(as.character(species_plot$Dose), levels = dose_order)
species_plot$captured_species <- factor(species_plot$captured_species, levels = c("Human tumor", "Mouse host"))
species_colors <- c(`Human tumor` = "#4C78A8", `Mouse host` = "#E67E22")
p_species <- ggplot2::ggplot(species_plot, ggplot2::aes(x = sample, y = fraction, fill = captured_species)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"), Ploidy = function(x) paste0("Ploidy ", x))
  ) +
  ggplot2::scale_fill_manual(values = species_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_label, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Human tumor versus mouse host captured-cell composition",
    subtitle = "Denominator = retained human tumor plus retained mouse QC singlets in the same Cell Ranger libraries",
    x = NULL, y = "Fraction of retained human + mouse singlets",
    caption = "Captured-cell fractions are affected by dissociation and recovery; they are not absolute histologic tissue fractions."
  ) + base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1))

nonimmune_colors <- unlist(cfg$annotation$final_nonimmune_colors, use.names = TRUE)
immune_lineage_colors <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
myeloid_lineage_colors <- unlist(cfg$myeloid$lineage_colors, use.names = TRUE)
mouse_type_colors <- stats::setNames(rep(NA_character_, length(mouse_type_levels)), mouse_type_levels)
for (label in mouse_type_levels) {
  if (startsWith(label, "Immune: ")) {
    subtype <- sub("^Immune: ", "", label)
    if (subtype %in% names(immune_lineage_colors)) mouse_type_colors[[label]] <- immune_lineage_colors[[subtype]]
    if (subtype %in% names(myeloid_lineage_colors)) mouse_type_colors[[label]] <- myeloid_lineage_colors[[subtype]]
  } else if (label %in% names(nonimmune_colors)) {
    mouse_type_colors[[label]] <- nonimmune_colors[[label]]
  }
}
if (anyNA(mouse_type_colors)) stop("Mouse composition palette lacks labels: ", paste(names(mouse_type_colors)[is.na(mouse_type_colors)], collapse = ", "), call. = FALSE)

mouse_type_plot <- mouse_type_counts
mouse_type_plot$sample <- factor(mouse_type_plot$sample, levels = sample_levels)
mouse_type_plot$mouse_final_cell_type_v2 <- factor(mouse_type_plot$mouse_final_cell_type_v2, levels = mouse_type_levels)
mouse_type_plot$Ploidy <- factor(mouse_type_plot$Ploidy, levels = ploidy_order)
mouse_type_plot$Dose <- factor(as.character(mouse_type_plot$Dose), levels = dose_order)
p_mouse_types <- ggplot2::ggplot(mouse_type_plot, ggplot2::aes(x = sample, y = fraction_of_mouse_cells, fill = mouse_final_cell_type_v2)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"), Ploidy = function(x) paste0("Ploidy ", x))
  ) +
  ggplot2::scale_fill_manual(values = mouse_type_colors, breaks = mouse_type_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_label, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Mouse host cell-type composition by sample",
    subtitle = "Denominator = all retained mouse QC singlets within each sample",
    x = NULL, y = "Fraction of mouse host cells",
    caption = "Cell-line samples are excluded."
  ) + base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1))

unified_plot <- unified_category_counts
unified_plot$sample <- factor(unified_plot$sample, levels = sample_levels)
unified_plot$captured_category <- factor(unified_plot$captured_category, levels = unified_category_levels)
unified_plot$Ploidy <- factor(unified_plot$Ploidy, levels = ploidy_order)
unified_plot$Dose <- factor(as.character(unified_plot$Dose), levels = dose_order)
unified_colors <- c(`Human tumor` = "#4C78A8", mouse_type_colors)
p_unified_categories <- ggplot2::ggplot(unified_plot, ggplot2::aes(x = sample, y = fraction, fill = captured_category)) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"), Ploidy = function(x) paste0("Ploidy ", x))
  ) +
  ggplot2::scale_fill_manual(values = unified_colors, breaks = unified_category_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), labels = percent_label, expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Human tumor plus mouse host cell-type composition",
    subtitle = "Denominator = retained human tumor plus retained mouse QC singlets",
    x = NULL, y = "Fraction of retained human + mouse singlets",
    caption = "The human category and mouse cell types share one sample-level captured-cell denominator; no cross-species expression clustering was performed."
  ) + base_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1))

# One 100% stacked bar per biological sample. The denominator is all cells with
# a final immune annotation within that sample.
p_sample_fraction <- ggplot2::ggplot(
  sample_plot,
  ggplot2::aes(x = sample, y = fraction, fill = immune_subtype)
) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(
      Dose = function(x) paste0(x, " mg/kg"),
      Ploidy = function(x) paste0("Ploidy ", x)
    )
  ) +
  ggplot2::scale_fill_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = sample_labels) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.2), labels = percent_label,
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Immune-cell subtype composition by sample",
    subtitle = "Samples are arranged by gemcitabine Dose and Ploidy; denominator = immune cells within each sample",
    x = NULL, y = "Fraction of immune cells",
    caption = "Immune: Ambiguous is retained as an explicit category. Cell-line samples are excluded."
  ) +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1)) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

p_sample_count <- ggplot2::ggplot(
  sample_plot,
  ggplot2::aes(x = sample, y = n_cells, fill = immune_subtype)
) +
  ggplot2::geom_col(width = 0.82, color = "white", linewidth = 0.18) +
  ggplot2::facet_wrap(
    ggplot2::vars(Ploidy, Dose), ncol = 3, scales = "free_x",
    labeller = ggplot2::labeller(
      Dose = function(x) paste0(x, " mg/kg"),
      Ploidy = function(x) paste0("Ploidy ", x)
    )
  ) +
  ggplot2::scale_fill_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_x_discrete(labels = stats::setNames(sample_levels, sample_levels)) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.04))) +
  ggplot2::labs(
    title = "Final immune-cell counts by sample",
    subtitle = "Absolute cell counts are shown as a sampling-depth companion to the composition plot",
    x = NULL, y = "Number of immune cells",
    caption = "Counts reflect retained cells and are descriptive; they are not normalized for sampling depth."
  ) +
  base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 52, hjust = 1, vjust = 1)) +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

# Each Dose x Ploidy bar is the arithmetic mean of biological-sample fractions,
# giving every sample equal weight regardless of its recovered cell count.
p_equal <- ggplot2::ggplot(
  equal_plot,
  ggplot2::aes(x = Dose, y = mean_sample_fraction, fill = immune_subtype)
) +
  ggplot2::geom_col(width = 0.72, color = "white", linewidth = 0.28) +
  ggplot2::geom_text(
    data = n_by_stratum,
    ggplot2::aes(x = Dose, y = 1.025, label = label),
    inherit.aes = FALSE, size = 3.3, color = "#343A40"
  ) +
  ggplot2::facet_grid(
    ~Ploidy,
    labeller = ggplot2::labeller(Ploidy = function(x) paste0("Ploidy ", x))
  ) +
  ggplot2::scale_fill_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 1, 0.2), labels = percent_label,
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1.06), clip = "off") +
  ggplot2::labs(
    title = "Immune-cell subtype composition by Dose and Ploidy",
    subtitle = "Mean of per-sample subtype fractions within each stratum; every biological sample has equal weight",
    x = "Gemcitabine dose (mg/kg)", y = "Mean fraction of immune cells",
    caption = "Final mouse immune cells only. This composition panel is descriptive; 2N-versus-4N tests are reported separately."
  ) +
  base_theme +
  ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2, byrow = TRUE))

# Companion view exposes all biological-sample observations. Shape and fill
# encode Ploidy so that the comparison does not rely on color alone.
p_distribution <- ggplot2::ggplot(
  sample_plot,
  ggplot2::aes(
    x = Dose, y = fraction, shape = Ploidy,
    color = immune_subtype, fill = immune_subtype, group = Ploidy
  )
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitterdodge(
      jitter.width = 0.08, jitter.height = 0, dodge.width = 0.52, seed = 1234
    ),
    size = 2.2, stroke = 0.45, alpha = 0.9
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = Ploidy), fun = mean, geom = "point",
    position = ggplot2::position_dodge(width = 0.52),
    shape = 95, size = 7
  ) +
  ggplot2::facet_wrap(~immune_subtype, scales = "free_y", ncol = 3) +
  ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
  ggplot2::scale_color_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_fill_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = function(x) percent_label(x, accuracy = 0.01),
    expand = ggplot2::expansion(mult = c(0.04, 0.1))
  ) +
  ggplot2::labs(
    title = "Between-sample immune-cell subtype fractions by Dose and Ploidy",
    subtitle = "Each point is one biological sample; horizontal ticks mark stratum means; facet y-axis ranges vary",
    x = "Gemcitabine dose (mg/kg)", y = "Fraction of immune cells",
    caption = "The denominator for each point is all final immune cells within that sample."
  ) +
  base_theme +
  ggplot2::guides(
    shape = ggplot2::guide_legend(title = NULL),
    color = "none",
    fill = "none"
  )

test_plot <- ploidy_tests
test_plot$Dose <- factor(as.character(test_plot$Dose), levels = dose_order)
test_plot$immune_subtype <- factor(test_plot$immune_subtype, levels = subtype_levels)
test_plot$y_min <- mapply(function(dose_i, subtype_i) {
  min(sample_counts$fraction[
    as.character(sample_counts$Dose) == as.character(dose_i) &
      sample_counts$immune_subtype == as.character(subtype_i)
  ])
}, test_plot$Dose, test_plot$immune_subtype)
test_plot$y_max <- mapply(function(dose_i, subtype_i) {
  max(sample_counts$fraction[
    as.character(sample_counts$Dose) == as.character(dose_i) &
      sample_counts$immune_subtype == as.character(subtype_i)
  ])
}, test_plot$Dose, test_plot$immune_subtype)
test_plot$y_pad <- pmax(
  (test_plot$y_max - test_plot$y_min) * 0.12,
  0.0001
)
test_plot$y_bracket <- test_plot$y_max + test_plot$y_pad
test_plot$y_label <- test_plot$y_max + test_plot$y_pad * 1.9

# Faceted biological-sample comparison with one exact 2N-versus-4N test per
# Dose and immune subtype. Stars use within-Dose BH-adjusted p-values.
p_ploidy_tests <- ggplot2::ggplot(
  sample_plot,
  ggplot2::aes(
    x = Ploidy, y = fraction, shape = Ploidy,
    color = immune_subtype, fill = immune_subtype
  )
) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 1234),
    size = 2.5, stroke = 0.5, alpha = 0.9
  ) +
  ggplot2::stat_summary(
    ggplot2::aes(group = Ploidy), fun = mean, geom = "point",
    shape = 95, size = 8
  ) +
  ggplot2::geom_segment(
    data = test_plot,
    ggplot2::aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
    inherit.aes = FALSE, color = "#343A40", linewidth = 0.4
  ) +
  ggplot2::geom_text(
    data = test_plot,
    ggplot2::aes(x = 1.5, y = y_label, label = significance_fdr),
    inherit.aes = FALSE, color = "#111827", fontface = "bold", size = 3.3
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(Dose, immune_subtype), ncol = length(subtype_levels), scales = "free_y",
    labeller = ggplot2::labeller(Dose = function(x) paste0(x, " mg/kg"))
  ) +
  ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
  ggplot2::scale_color_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_fill_manual(values = palette_values, breaks = subtype_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(
    labels = function(x) percent_label(x, accuracy = 0.001),
    expand = ggplot2::expansion(mult = c(0.04, 0.06))
  ) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Immune-cell subtype fractions: 2N versus 4N within each Dose",
    subtitle = "Each point is one biological sample; horizontal ticks show means; stars use within-Dose BH-adjusted p-values",
    x = "Ploidy", y = "Fraction of immune cells",
    caption = "Two-sided exact permutation tests of mean sample fractions. **** FDR<0.0001; *** <0.001; ** <0.01; * <0.05; ns otherwise."
  ) +
  base_theme +
  ggplot2::guides(shape = "none", color = "none", fill = "none")

sample_fraction_pdf <- file.path(out_dir, "immune_subtype_composition_by_sample.pdf")
sample_fraction_png <- file.path(out_dir, "immune_subtype_composition_by_sample.png")
sample_count_pdf <- file.path(out_dir, "immune_subtype_counts_by_sample.pdf")
sample_count_png <- file.path(out_dir, "immune_subtype_counts_by_sample.png")
equal_pdf <- file.path(out_dir, "immune_subtype_composition_by_dose_and_ploidy_equal_sample_weight.pdf")
equal_png <- file.path(out_dir, "immune_subtype_composition_by_dose_and_ploidy_equal_sample_weight.png")
distribution_pdf <- file.path(out_dir, "immune_subtype_fraction_distribution_by_dose_and_ploidy.pdf")
distribution_png <- file.path(out_dir, "immune_subtype_fraction_distribution_by_dose_and_ploidy.png")
ploidy_tests_pdf <- file.path(out_dir, "immune_subtype_ploidy_comparison_by_dose.pdf")
ploidy_tests_png <- file.path(out_dir, "immune_subtype_ploidy_comparison_by_dose.png")
state_sample_pdf <- file.path(out_dir, "neutrophil_state_composition_by_sample.pdf")
state_sample_png <- file.path(out_dir, "neutrophil_state_composition_by_sample.png")
state_equal_pdf <- file.path(out_dir, "neutrophil_state_composition_by_dose_and_ploidy_equal_sample_weight.pdf")
state_equal_png <- file.path(out_dir, "neutrophil_state_composition_by_dose_and_ploidy_equal_sample_weight.png")
state_tests_pdf <- file.path(out_dir, "neutrophil_state_ploidy_comparison_by_dose.pdf")
state_tests_png <- file.path(out_dir, "neutrophil_state_ploidy_comparison_by_dose.png")
species_pdf <- file.path(out_dir, "human_vs_mouse_captured_composition_by_sample.pdf")
species_png <- file.path(out_dir, "human_vs_mouse_captured_composition_by_sample.png")
mouse_types_pdf <- file.path(out_dir, "mouse_cell_type_composition_within_mouse_by_sample.pdf")
mouse_types_png <- file.path(out_dir, "mouse_cell_type_composition_within_mouse_by_sample.png")
unified_categories_pdf <- file.path(out_dir, "human_tumor_plus_mouse_cell_type_composition_by_sample.pdf")
unified_categories_png <- file.path(out_dir, "human_tumor_plus_mouse_cell_type_composition_by_sample.png")

ggplot2::ggsave(sample_fraction_pdf, p_sample_fraction, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(sample_fraction_png, p_sample_fraction, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(sample_count_pdf, p_sample_count, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(sample_count_png, p_sample_count, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(equal_pdf, p_equal, width = 11.2, height = 7.6, units = "in", bg = "white")
ggplot2::ggsave(equal_png, p_equal, width = 11.2, height = 7.6, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(distribution_pdf, p_distribution, width = 11.5, height = 8.5, units = "in", bg = "white")
ggplot2::ggsave(distribution_png, p_distribution, width = 11.5, height = 8.5, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(ploidy_tests_pdf, p_ploidy_tests, width = 18, height = 10.5, units = "in", bg = "white")
ggplot2::ggsave(ploidy_tests_png, p_ploidy_tests, width = 18, height = 10.5, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(state_sample_pdf, p_state_sample, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(state_sample_png, p_state_sample, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(state_equal_pdf, p_state_equal, width = 11.2, height = 7.6, units = "in", bg = "white")
ggplot2::ggsave(state_equal_png, p_state_equal, width = 11.2, height = 7.6, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(state_tests_pdf, p_state_tests, width = max(18, 2.8 * length(observed_states)), height = 10.5, units = "in", bg = "white")
ggplot2::ggsave(state_tests_png, p_state_tests, width = max(18, 2.8 * length(observed_states)), height = 10.5, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(species_pdf, p_species, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(species_png, p_species, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(mouse_types_pdf, p_mouse_types, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(mouse_types_png, p_mouse_types, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")
ggplot2::ggsave(unified_categories_pdf, p_unified_categories, width = 15, height = 10.2, units = "in", bg = "white")
ggplot2::ggsave(unified_categories_png, p_unified_categories, width = 15, height = 10.2, units = "in", dpi = 300, bg = "white")

output_files <- c(
  file.path(out_dir, "immune_composition_by_sample.csv"),
  file.path(out_dir, "immune_sample_ploidy_dose_map.csv"),
  file.path(out_dir, "immune_composition_by_dose_ploidy_equal_sample_weight.csv"),
  file.path(out_dir, "immune_composition_by_dose_ploidy_pooled_cells.csv"),
  file.path(out_dir, "immune_subtype_ploidy_tests_by_dose.csv"),
  file.path(out_dir, "immune_subtype_palette.csv"),
  file.path(out_dir, "final_immune_cell_metadata.csv.gz"),
  file.path(out_dir, "immune_composition_contract.yaml"),
  file.path(out_dir, "input_manifest_sha256.csv"),
  file.path(out_dir, "immune_composition_audit.tsv"),
  file.path(out_dir, "neutrophil_state_composition_by_sample.csv"),
  file.path(out_dir, "neutrophil_state_composition_by_dose_ploidy_equal_sample_weight.csv"),
  file.path(out_dir, "neutrophil_state_ploidy_tests_by_dose.csv"),
  file.path(out_dir, "neutrophil_state_palette.csv"),
  file.path(out_dir, "human_vs_mouse_captured_composition_by_sample.csv"),
  file.path(out_dir, "mouse_cell_types_with_mouse_and_unified_denominators_by_sample.csv"),
  file.path(out_dir, "human_tumor_plus_mouse_cell_type_composition_by_sample.csv"),
  sample_fraction_pdf,
  sample_fraction_png,
  sample_count_pdf,
  sample_count_png,
  equal_pdf,
  equal_png,
  distribution_pdf,
  distribution_png,
  ploidy_tests_pdf,
  ploidy_tests_png,
  state_sample_pdf,
  state_sample_png,
  state_equal_pdf,
  state_equal_png,
  state_tests_pdf,
  state_tests_png,
  species_pdf,
  species_png,
  mouse_types_pdf,
  mouse_types_png,
  unified_categories_pdf,
  unified_categories_png
)
output_info <- file.info(output_files)
output_manifest <- data.frame(
  path = normalizePath(output_files, mustWork = TRUE),
  size_bytes = output_info$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05z.txt"))
writeLines(
  c(
    "05z final composition data preparation and plotting completed.",
    paste0("Final immune cells: ", nrow(immune)),
    paste0("Samples: ", length(unique(immune$sample))),
    paste0("Ploidy levels: ", paste(ploidy_order, collapse = ", ")),
    paste0("Dose levels (mg/kg): ", paste(dose_order, collapse = ", ")),
    paste0("Observed immune subtypes: ", paste(observed_subtypes, collapse = ", ")),
    "All sample fractions and equal-weight Dose x Ploidy means sum to 1.",
    paste0("Exact 2N-versus-4N permutation tests: ", nrow(ploidy_tests), "; within-Dose FDR<0.05: ", sum(ploidy_tests$p_value_bh_within_dose < 0.05)),
    paste0("Final Neutrophil cells: ", nrow(neutrophils), "; observed states: ", paste(observed_states, collapse = ", ")),
    paste0("Neutrophil-state 2N-versus-4N tests: ", nrow(state_ploidy_tests), "; within-Dose FDR<0.05: ", sum(state_ploidy_tests$p_value_bh_within_dose < 0.05)),
    paste0("Unified retained human + mouse denominator cells: ", nrow(census)),
    "Immune-lineage and Neutrophil-state colors match the fixed 05m final UMAP palettes.",
    "Unified proportions are scRNA-seq captured-cell fractions, not absolute histologic tissue fractions."
  ),
  file.path(out_dir, "completion.txt")
)
message("[05z] Completed. Output: ", out_dir)
