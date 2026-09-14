#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

stop_if_missing_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

exact_two_sample_permutation <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  if (length(x) < 2L || length(y) < 2L || any(!is.finite(c(x, y)))) {
    stop("Each Ploidy group must contain at least two finite biological-sample fractions.", call. = FALSE)
  }
  pooled <- c(x, y)
  combinations <- utils::combn(seq_along(pooled), length(x))
  permuted_differences <- apply(combinations, 2L, function(idx) mean(pooled[idx]) - mean(pooled[-idx]))
  observed_difference <- mean(x) - mean(y)
  tolerance <- sqrt(.Machine$double.eps)
  list(
    observed_difference = observed_difference,
    p_value = mean(abs(permuted_differences) >= abs(observed_difference) - tolerance),
    n_permutations = ncol(combinations)
  )
}

p_value_stars <- function(p) {
  ifelse(!is.finite(p), "NA", ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))))
}

run_group_tests <- function(sample_fractions, analysis_family, group_order, dose_i) {
  rows <- lapply(group_order, function(group_i) {
    z <- sample_fractions[sample_fractions$final_immune_group == group_i, , drop = FALSE]
    x <- z$fraction[z$Ploidy == "2N"]
    y <- z$fraction[z$Ploidy == "4N"]
    test <- exact_two_sample_permutation(x, y)
    data.frame(
      analysis_family = analysis_family,
      Dose = dose_i,
      final_immune_group = group_i,
      denominator = unique(z$denominator),
      n_samples_2N = length(x), n_samples_4N = length(y),
      mean_fraction_2N = mean(x), mean_fraction_4N = mean(y),
      median_fraction_2N = stats::median(x), median_fraction_4N = stats::median(y),
      difference_2N_minus_4N = test$observed_difference,
      difference_percentage_points = 100 * test$observed_difference,
      test = "two-sided exact permutation test of the difference in mean biological-sample fractions",
      n_permutations = test$n_permutations,
      p_value = test$p_value,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_value_bh <- stats::p.adjust(out$p_value, method = "BH")
  out$significance_fdr <- p_value_stars(out$p_value_bh)
  out$final_immune_group_order <- match(out$final_immune_group, group_order)
  out
}

make_comparison_plot <- function(sample_fractions, tests, group_colors, dose_i, title, y_title, output_stem, out_dir) {
  group_order <- names(group_colors)
  sample_fractions$final_immune_group <- factor(sample_fractions$final_immune_group, levels = group_order)
  sample_fractions$Ploidy <- factor(sample_fractions$Ploidy, levels = c("2N", "4N"))
  tests$final_immune_group <- factor(tests$final_immune_group, levels = group_order)
  tests$y_min <- vapply(as.character(tests$final_immune_group), function(group_i) {
    min(sample_fractions$fraction[as.character(sample_fractions$final_immune_group) == group_i])
  }, numeric(1))
  tests$y_max <- vapply(as.character(tests$final_immune_group), function(group_i) {
    max(sample_fractions$fraction[as.character(sample_fractions$final_immune_group) == group_i])
  }, numeric(1))
  tests$y_pad <- pmax((tests$y_max - tests$y_min) * 0.12, 0.0001)
  tests$y_bracket <- tests$y_max + tests$y_pad
  tests$y_label <- tests$y_max + 1.9 * tests$y_pad
  n_2n <- length(unique(sample_fractions$sample[sample_fractions$Ploidy == "2N"]))
  n_4n <- length(unique(sample_fractions$sample[sample_fractions$Ploidy == "4N"]))
  n_allocations <- unique(tests$n_permutations)
  if (length(n_allocations) != 1L) stop("Permutation counts differ across final immune groups.", call. = FALSE)

  p <- ggplot2::ggplot(
    sample_fractions,
    ggplot2::aes(x = Ploidy, y = fraction, shape = Ploidy, color = final_immune_group, fill = final_immune_group)
  ) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 1234), size = 2.8, stroke = 0.6, alpha = 0.9) +
    ggplot2::stat_summary(ggplot2::aes(group = Ploidy), fun = mean, geom = "point", shape = 95, size = 8) +
    ggplot2::geom_segment(
      data = tests, ggplot2::aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
      inherit.aes = FALSE, color = "#343A40", linewidth = 0.45
    ) +
    ggplot2::geom_text(
      data = tests, ggplot2::aes(x = 1.5, y = y_label, label = significance_fdr),
      inherit.aes = FALSE, color = "#111827", fontface = "bold", size = 3.8
    ) +
    ggplot2::facet_wrap(ggplot2::vars(final_immune_group), scales = "free_y", ncol = 3) +
    ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
    ggplot2::scale_color_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::scale_y_continuous(
      labels = function(x) paste0(format(round(100 * x, 2), trim = TRUE), "%"),
      expand = ggplot2::expansion(mult = c(0.04, 0.10))
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::labs(
      title = title,
      subtitle = paste0("Dose = ", dose_i, " mg/kg; each point is one biological sample (2N n=", n_2n, ", 4N n=", n_4n, ")"),
      x = "Ploidy", y = y_title,
      caption = paste0("Exact two-sided permutation tests (", n_allocations, " allocations); labels use BH-adjusted p-values across all 9 final immune groups within this Dose.")
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(), strip.background = ggplot2::element_rect(fill = "#E5E7EB"),
      strip.text = ggplot2::element_text(face = "bold"), plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "none", plot.margin = ggplot2::margin(8, 18, 8, 8)
    )
  ggplot2::ggsave(file.path(out_dir, paste0(output_stem, ".pdf")), p, width = 12, height = 10)
  ggplot2::ggsave(file.path(out_dir, paste0(output_stem, ".png")), p, width = 12, height = 10, dpi = 300)
}

script_dir <- resolve_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
annotation_csv <- file.path(results_root, "05m_myeloid_states_final", "mouse_cell_annotations_with_myeloid_states.csv")
census_summary_csv <- file.path(results_root, "05n_unified_human_mouse_census", "unified_census_sample_summary.csv")
out_dir <- file.path(results_root, "05o_dose_stratified_2N_vs_4N_final_immune_groups")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(annotation_csv)) stop("Missing 05m annotation CSV: ", annotation_csv, call. = FALSE)
if (!file.exists(census_summary_csv)) stop("Missing 05n census summary: ", census_summary_csv, call. = FALSE)

sample_metadata <- do.call(rbind, lapply(names(cfg$immune_composition$sample_metadata), function(sample_id) {
  z <- cfg$immune_composition$sample_metadata[[sample_id]]
  data.frame(sample = sample_id, Ploidy = as.character(z$Ploidy), Dose = as.character(z$Dose), stringsAsFactors = FALSE)
}))
dose_order <- as.character(unlist(cfg$immune_composition$dose_order, use.names = FALSE))
if (!identical(dose_order, c("0", "30", "120"))) stop("05o requires Dose order 0, 30, 120 mg/kg.", call. = FALSE)
if (any(!sample_metadata$Dose %in% dose_order)) stop("Sample metadata contains an unconfigured Dose.", call. = FALSE)
sample_metadata$Ploidy <- factor(sample_metadata$Ploidy, levels = c("2N", "4N"))
sample_metadata$Dose <- factor(sample_metadata$Dose, levels = dose_order)
replicate_table <- table(sample_metadata$Dose, sample_metadata$Ploidy)
if (any(replicate_table < 2L)) stop("Every Dose x Ploidy group must contain at least two biological samples.", call. = FALSE)
if (any(replicate_table[, "2N"] != replicate_table[, "4N"])) stop("2N and 4N sample counts must match within every Dose.", call. = FALSE)

cells <- utils::read.csv(annotation_csv, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("NA", ""))
stop_if_missing_columns(
  cells,
  c("sample", "mouse_final_cell_type_state", "mouse_final_immune_type", "mouse_myeloid_recluster", "scDblFinder.class", "mouse_fraction"),
  "05m final annotation"
)
cells <- cells[cells$sample %in% sample_metadata$sample, , drop = FALSE]
if (any(cells$scDblFinder.class != "singlet", na.rm = TRUE)) stop("05m cells contain non-singlets.", call. = FALSE)
if (min(cells$mouse_fraction, na.rm = TRUE) < as.numeric(cfg$species$mouse_fraction_min)) stop("05m cells violate the mouse-fraction threshold.", call. = FALSE)

is_final_immune <- !is.na(cells$mouse_final_immune_type) & startsWith(cells$mouse_final_cell_type_state, "Immune: ")
immune_cells <- cells[is_final_immune, c("sample", "mouse_final_cell_type_state", "mouse_myeloid_recluster"), drop = FALSE]
immune_cells$final_immune_group <- sub("^Immune: ", "", immune_cells$mouse_final_cell_type_state)

state_order <- as.character(unlist(cfg$myeloid$state_order, use.names = FALSE))
configured_group_order <- c("Monocyte", "Macrophage", paste0("Neutrophil — ", state_order))
observed_group_order <- configured_group_order[configured_group_order %in% unique(immune_cells$final_immune_group)]
unconfigured_groups <- setdiff(unique(immune_cells$final_immune_group), configured_group_order)
if (length(unconfigured_groups) > 0L) stop("Final immune groups are absent from configured order: ", paste(unconfigured_groups, collapse = ", "), call. = FALSE)
if (length(observed_group_order) != 9L) stop("Expected exactly 9 final immune groups; observed: ", paste(observed_group_order, collapse = ", "), call. = FALSE)

# This verifies that the fine labels being counted remain cluster-level labels.
cluster_uniformity <- tapply(immune_cells$final_immune_group, immune_cells$mouse_myeloid_recluster, function(x) length(unique(x)))
if (any(cluster_uniformity != 1L)) stop("A reclustered myeloid cluster contains multiple final immune-group labels.", call. = FALSE)

count_grid <- expand.grid(sample = sample_metadata$sample, final_immune_group = observed_group_order, stringsAsFactors = FALSE)
observed_counts <- as.data.frame(table(sample = immune_cells$sample, final_immune_group = immune_cells$final_immune_group), stringsAsFactors = FALSE)
names(observed_counts)[3L] <- "n_cells"
group_counts <- merge(count_grid, observed_counts, by = c("sample", "final_immune_group"), all.x = TRUE, sort = FALSE)
group_counts$n_cells[is.na(group_counts$n_cells)] <- 0L

census <- utils::read.csv(census_summary_csv, stringsAsFactors = FALSE, check.names = FALSE)
stop_if_missing_columns(census, c("sample", "retained_human_tumor", "retained_mouse_qc_singlet", "main_denominator_cells"), "05n census summary")
sample_audit <- merge(sample_metadata, census, by = "sample", all.x = TRUE, sort = FALSE)
mouse_counts <- as.data.frame(table(sample = cells$sample), stringsAsFactors = FALSE)
names(mouse_counts)[2L] <- "mouse_cells_in_05m"
sample_audit <- merge(sample_audit, mouse_counts, by = "sample", all.x = TRUE, sort = FALSE)
if (anyNA(sample_audit$main_denominator_cells) || anyNA(sample_audit$retained_mouse_qc_singlet)) stop("05n census is incomplete for configured samples.", call. = FALSE)
if (any(sample_audit$mouse_cells_in_05m != sample_audit$retained_mouse_qc_singlet)) stop("05m and 05n mouse denominators disagree.", call. = FALSE)
immune_totals <- stats::aggregate(n_cells ~ sample, group_counts, sum)
names(immune_totals)[2L] <- "total_final_mouse_immune_cells"
sample_audit <- merge(sample_audit, immune_totals, by = "sample", all.x = TRUE, sort = FALSE)
sample_audit <- sample_audit[match(sample_metadata$sample, sample_audit$sample), , drop = FALSE]

group_counts <- merge(
  group_counts,
  sample_audit[c("sample", "Ploidy", "Dose", "retained_mouse_qc_singlet", "main_denominator_cells", "total_final_mouse_immune_cells")],
  by = "sample", all.x = TRUE, sort = FALSE
)
group_counts$Ploidy <- factor(as.character(group_counts$Ploidy), levels = c("2N", "4N"))
group_counts$fraction_unified <- group_counts$n_cells / group_counts$main_denominator_cells
group_counts$fraction_mouse <- group_counts$n_cells / group_counts$retained_mouse_qc_singlet
group_counts$fraction_within_immune <- group_counts$n_cells / group_counts$total_final_mouse_immune_cells
group_counts$final_immune_group_order <- match(group_counts$final_immune_group, observed_group_order)
group_counts <- group_counts[order(match(group_counts$sample, sample_metadata$sample), group_counts$final_immune_group_order), , drop = FALSE]

make_fraction_table <- function(fraction_column, denominator_label) {
  data.frame(
    sample = group_counts$sample, Ploidy = group_counts$Ploidy, Dose = group_counts$Dose,
    final_immune_group = group_counts$final_immune_group, numerator_cells = group_counts$n_cells,
    fraction = group_counts[[fraction_column]], denominator = denominator_label,
    final_immune_group_order = group_counts$final_immune_group_order,
    stringsAsFactors = FALSE
  )
}
primary_fractions <- make_fraction_table("fraction_unified", "retained human tumor + retained mouse QC singlets")
mouse_fractions <- make_fraction_table("fraction_mouse", "retained mouse QC singlets")
within_immune_fractions <- make_fraction_table("fraction_within_immune", "all final mouse immune cells")

run_tests_by_dose <- function(fractions, analysis_family) {
  do.call(rbind, lapply(dose_order, function(dose_i) {
    run_group_tests(
      fractions[as.character(fractions$Dose) == dose_i, , drop = FALSE],
      analysis_family, observed_group_order, dose_i
    )
  }))
}
primary_tests <- run_tests_by_dose(primary_fractions, "primary_unified_human_mouse_denominator")
mouse_tests <- run_tests_by_dose(mouse_fractions, "sensitivity_mouse_only_denominator")
within_immune_tests <- run_tests_by_dose(within_immune_fractions, "sensitivity_within_final_mouse_immune_composition")

utils::write.csv(sample_audit, file.path(out_dir, "all_doses_sample_denominator_audit.csv"), row.names = FALSE, na = "NA")
utils::write.csv(group_counts, file.path(out_dir, "all_doses_final_immune_group_counts_and_fractions_by_sample.csv"), row.names = FALSE, na = "NA")
utils::write.csv(primary_fractions, file.path(out_dir, "all_doses_final_immune_groups_primary_unified_denominator_sample_fractions.csv"), row.names = FALSE, na = "NA")
utils::write.csv(primary_tests, file.path(out_dir, "all_doses_final_immune_groups_primary_unified_denominator_2N_vs_4N_tests.csv"), row.names = FALSE, na = "NA")
utils::write.csv(mouse_fractions, file.path(out_dir, "all_doses_final_immune_groups_mouse_only_denominator_sample_fractions.csv"), row.names = FALSE, na = "NA")
utils::write.csv(mouse_tests, file.path(out_dir, "all_doses_final_immune_groups_mouse_only_denominator_2N_vs_4N_tests.csv"), row.names = FALSE, na = "NA")
utils::write.csv(within_immune_fractions, file.path(out_dir, "all_doses_final_immune_groups_within_immune_sample_fractions.csv"), row.names = FALSE, na = "NA")
utils::write.csv(within_immune_tests, file.path(out_dir, "all_doses_final_immune_groups_within_immune_2N_vs_4N_tests.csv"), row.names = FALSE, na = "NA")

subtype_colors <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
state_colors <- unlist(cfg$myeloid$state_colors, use.names = TRUE)
group_colors <- c(
  Monocyte = subtype_colors[["Monocyte"]],
  Macrophage = subtype_colors[["Macrophage"]],
  stats::setNames(unname(state_colors[state_order]), paste0("Neutrophil — ", state_order))
)[observed_group_order]
if (anyNA(group_colors)) stop("Final immune-group palette is incomplete.", call. = FALSE)
palette_df <- data.frame(final_immune_group = observed_group_order, color_hex = unname(group_colors), final_immune_group_order = seq_along(observed_group_order), stringsAsFactors = FALSE)
utils::write.csv(palette_df, file.path(out_dir, "final_immune_group_palette.csv"), row.names = FALSE)

for (dose_i in dose_order) {
  dose_stem <- paste0("dose", dose_i)
  dose_primary_fractions <- primary_fractions[as.character(primary_fractions$Dose) == dose_i, , drop = FALSE]
  dose_within_fractions <- within_immune_fractions[as.character(within_immune_fractions$Dose) == dose_i, , drop = FALSE]
  dose_primary_tests <- primary_tests[primary_tests$Dose == dose_i, , drop = FALSE]
  dose_mouse_tests <- mouse_tests[mouse_tests$Dose == dose_i, , drop = FALSE]
  dose_within_tests <- within_immune_tests[within_immune_tests$Dose == dose_i, , drop = FALSE]
  utils::write.csv(dose_primary_tests, file.path(out_dir, paste0(dose_stem, "_final_immune_groups_primary_unified_denominator_2N_vs_4N_tests.csv")), row.names = FALSE, na = "NA")
  utils::write.csv(dose_mouse_tests, file.path(out_dir, paste0(dose_stem, "_final_immune_groups_mouse_only_denominator_2N_vs_4N_tests.csv")), row.names = FALSE, na = "NA")
  utils::write.csv(dose_within_tests, file.path(out_dir, paste0(dose_stem, "_final_immune_groups_within_immune_2N_vs_4N_tests.csv")), row.names = FALSE, na = "NA")
  make_comparison_plot(
    dose_primary_fractions, dose_primary_tests, group_colors, dose_i,
    paste0("Dose=", dose_i, " final immune groups: 2N versus 4N"),
    "Fraction of retained human + mouse captured cells",
    paste0(dose_stem, "_final_immune_groups_primary_unified_denominator_2N_vs_4N"), out_dir
  )
  make_comparison_plot(
    dose_within_fractions, dose_within_tests, group_colors, dose_i,
    paste0("Dose=", dose_i, " final immune-group composition: 2N versus 4N"),
    "Fraction within all final mouse immune cells",
    paste0(dose_stem, "_final_immune_groups_within_immune_2N_vs_4N"), out_dir
  )
}

contract <- list(
  dose_mg_per_kg = as.numeric(dose_order),
  biological_replicates = lapply(dose_order, function(dose_i) {
    counts <- replicate_table[dose_i, c("2N", "4N")]
    list(`2N` = as.integer(counts[["2N"]]), `4N` = as.integer(counts[["4N"]]))
  }),
  final_immune_groups = observed_group_order,
  biological_annotation_level = "cluster; each 05j myeloid cluster has exactly one final lineage/state label",
  primary_denominator = "retained human tumor cells plus retained mouse QC singlets in the same Cell Ranger library",
  sensitivity_denominators = c("retained mouse QC singlets", "all final mouse immune cells"),
  unit_of_inference = "biological sample; cells are not treated as independent replicates",
  test = "two-sided exact permutation test of the difference in mean sample fractions within each Dose; 70 allocations at Dose=0 and 6 allocations at Dose=30/120",
  multiple_testing = "Benjamini-Hochberg across all 9 final immune groups separately within each Dose and analysis family",
  low_replication_warning = "Dose=30 and Dose=120 each have n=2 per Ploidy and only 6 exact allocations, so inferential resolution and power are very limited",
  interpretation = "captured-cell fractions, not absolute histologic tissue fractions"
)
names(contract$biological_replicates) <- paste0("Dose_", dose_order)
yaml::write_yaml(contract, file.path(out_dir, "05o_analysis_contract.yaml"))

summary_lines <- c(
  "05o Dose-stratified final immune-group comparisons completed for Dose=0, 30, and 120 mg/kg.",
  paste0("Final cluster-level immune groups tested (n=", length(observed_group_order), "): ", paste(observed_group_order, collapse = "; "), "."),
  unlist(lapply(dose_order, function(dose_i) {
    p <- primary_tests[primary_tests$Dose == dose_i, , drop = FALSE]
    w <- within_immune_tests[within_immune_tests$Dose == dose_i, , drop = FALSE]
    nominal_primary <- p$final_immune_group[p$p_value < 0.05]
    sig_primary <- p$final_immune_group[p$p_value_bh < 0.05]
    sig_within <- w$final_immune_group[w$p_value_bh < 0.05]
    n_2n <- replicate_table[dose_i, "2N"]
    n_4n <- replicate_table[dose_i, "4N"]
    c(
      paste0("Dose=", dose_i, " mg/kg: 2N n=", n_2n, "; 4N n=", n_4n, "; exact allocations=", unique(p$n_permutations), "."),
      paste0("Dose=", dose_i, " primary unified-denominator groups with nominal p < 0.05: ", if (length(nominal_primary) == 0L) "none" else paste(nominal_primary, collapse = ", "), "."),
      paste0("Dose=", dose_i, " primary unified-denominator groups with BH-FDR < 0.05: ", if (length(sig_primary) == 0L) "none" else paste(sig_primary, collapse = ", "), "."),
      paste0("Dose=", dose_i, " within-final-immune groups with BH-FDR < 0.05: ", if (length(sig_within) == 0L) "none" else paste(sig_within, collapse = ", "), ".")
    )
  })),
  "Dose=30 and Dose=120 each have only two biological samples per Ploidy and six exact allocations; inferential resolution and power are very limited.",
  "The primary denominator includes retained human tumor cells and retained mouse QC singlets; no human-mouse joint expression clustering was performed.",
  "Fractions are captured-cell fractions and are not absolute histologic tissue fractions."
)
writeLines(summary_lines, file.path(out_dir, "results_summary.txt"))

input_files <- c(annotation_csv, census_summary_csv, normalizePath(config_path, mustWork = TRUE))
input_info <- file.info(input_files)
input_manifest <- data.frame(
  role = c("05m final cluster-derived mouse annotation", "05n unified human-mouse sample denominators", "05 analysis config"),
  path = normalizePath(input_files, mustWork = TRUE), size_bytes = input_info$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05o.txt"))
writeLines(summary_lines, file.path(out_dir, "completion.txt"))

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_info <- file.info(output_files)
output_manifest <- data.frame(
  path = normalizePath(output_files, mustWork = TRUE), size_bytes = output_info$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
message("[05o] Completed. Output: ", out_dir)
