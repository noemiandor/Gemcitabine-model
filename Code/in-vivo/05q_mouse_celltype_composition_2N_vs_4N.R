#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 260)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

percent_label <- function(x, accuracy = 1) {
  paste0(format(round(100 * x / accuracy) * accuracy, trim = TRUE), "%")
}

p_value_stars <- function(p) {
  ifelse(
    !is.finite(p), "NA",
    ifelse(p < 0.0001, "****", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
  )
}

build_stratified_ploidy_assignments <- function(sample_metadata) {
  sample_metadata <- sample_metadata[order(sample_metadata$sample_order), , drop = FALSE]
  blocks <- split(seq_len(nrow(sample_metadata)), as.character(sample_metadata$Dose))
  choices <- lapply(blocks, function(idx) {
    n_2n <- sum(as.character(sample_metadata$Ploidy[idx]) == "2N")
    combn(idx, n_2n, simplify = FALSE)
  })
  choice_grid <- expand.grid(lapply(choices, seq_along), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  assignments <- matrix(FALSE, nrow = nrow(sample_metadata), ncol = nrow(choice_grid))
  for (j in seq_len(nrow(choice_grid))) {
    for (b in seq_along(choices)) assignments[choices[[b]][[choice_grid[j, b]]], j] <- TRUE
  }
  rownames(assignments) <- sample_metadata$sample
  assignments
}

run_category_tests <- function(counts, category_column, assignments) {
  categories <- unique(as.character(counts[[category_column]])[order(counts$category_order)])
  rows <- lapply(categories, function(category_i) {
    z <- counts[as.character(counts[[category_column]]) == category_i, , drop = FALSE]
    z <- z[match(rownames(assignments), z$sample), , drop = FALSE]
    if (anyNA(z$sample)) stop("A category table is missing one or more samples.", call. = FALSE)
    y <- as.numeric(z$fraction_within_mouse)
    observed_2n <- as.character(z$Ploidy) == "2N"
    observed_difference <- mean(y[observed_2n]) - mean(y[!observed_2n])
    n_2n <- colSums(assignments)
    sum_2n <- as.numeric(crossprod(y, assignments))
    permuted_difference <- sum_2n / n_2n - (sum(y) - sum_2n) / (length(y) - n_2n)
    p_value <- mean(abs(permuted_difference) >= abs(observed_difference) - sqrt(.Machine$double.eps))

    total_column <- grep("^sample_total_mouse_cells$", names(z), value = TRUE)
    sensitivity_beta <- sensitivity_p <- NA_real_
    if (length(total_column) == 1L && length(unique(z$n_cells)) > 1L && any(z$n_cells > 0L)) {
      fit <- stats::glm(
        cbind(n_cells, z[[total_column]] - n_cells) ~ Ploidy + Dose,
        family = stats::quasibinomial(), data = z
      )
      coefficient_table <- summary(fit)$coefficients
      if ("Ploidy4N" %in% rownames(coefficient_table)) {
        sensitivity_beta <- unname(coefficient_table["Ploidy4N", "Estimate"])
        sensitivity_p <- unname(coefficient_table["Ploidy4N", "Pr(>|t|)"])
      }
    }
    data.frame(
      category = category_i,
      n_samples_2N = sum(observed_2n),
      n_samples_4N = sum(!observed_2n),
      nonzero_samples_2N = sum(z$n_cells[observed_2n] > 0L),
      nonzero_samples_4N = sum(z$n_cells[!observed_2n] > 0L),
      mean_fraction_2N = mean(y[observed_2n]),
      mean_fraction_4N = mean(y[!observed_2n]),
      difference_2N_minus_4N = observed_difference,
      fold_change_2N_over_4N_with_0.5_cell_pseudocount =
        ((sum(z$n_cells[observed_2n]) + 0.5) / (sum(z$sample_total_mouse_cells[observed_2n]) + 0.5)) /
        ((sum(z$n_cells[!observed_2n]) + 0.5) / (sum(z$sample_total_mouse_cells[!observed_2n]) + 0.5)),
      test = "two-sided exact Dose-stratified permutation test of the difference in mean sample fractions",
      n_exact_stratified_assignments = ncol(assignments),
      p_value = p_value,
      quasibinomial_log_odds_4N_vs_2N = sensitivity_beta,
      quasibinomial_p_value = sensitivity_p,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$p_value_bh <- stats::p.adjust(out$p_value, method = "BH")
  out$significance_fdr <- p_value_stars(out$p_value_bh)
  out$quasibinomial_p_value_bh <- stats::p.adjust(out$quasibinomial_p_value, method = "BH")
  out$category_order <- match(out$category, categories)
  out[order(out$category_order), , drop = FALSE]
}

composition_summary <- function(counts, category_column) {
  categories <- unique(as.character(counts[[category_column]])[order(counts$category_order)])
  groups <- split(counts, interaction(counts$Ploidy, counts[[category_column]], drop = TRUE, lex.order = TRUE))
  rows <- lapply(groups, function(z) {
    values <- as.numeric(z$fraction_within_mouse)
    n <- length(values)
    sem <- if (n > 1L) stats::sd(values) / sqrt(n) else NA_real_
    margin <- if (n > 1L) stats::qt(0.975, df = n - 1L) * sem else NA_real_
    data.frame(
      Ploidy = as.character(z$Ploidy[[1L]]),
      category = as.character(z[[category_column]][[1L]]),
      n_samples = n,
      mean_sample_fraction = mean(values),
      sd_sample_fraction = stats::sd(values),
      sem_sample_fraction = sem,
      ci95_lower = max(0, mean(values) - margin),
      ci95_upper = min(1, mean(values) + margin),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$Ploidy <- factor(out$Ploidy, levels = c("2N", "4N"))
  out$category_order <- match(out$category, categories)
  out[order(out$category_order, out$Ploidy), , drop = FALSE]
}

multivariate_partial_f <- function(Y, dose, ploidy) {
  reduced <- stats::model.matrix(~ factor(dose))
  full <- stats::model.matrix(~ factor(dose) + factor(ploidy))
  sse <- function(X) {
    fitted <- X %*% qr.coef(qr(X), Y)
    sum((Y - fitted)^2)
  }
  sse_reduced <- sse(reduced)
  sse_full <- sse(full)
  df_effect <- qr(full)$rank - qr(reduced)$rank
  df_residual <- nrow(Y) - qr(full)$rank
  ((sse_reduced - sse_full) / df_effect) / (sse_full / df_residual)
}

run_global_composition_test <- function(counts, category_column, assignments) {
  samples <- rownames(assignments)
  categories <- unique(as.character(counts[[category_column]])[order(counts$category_order)])
  count_matrix <- matrix(0, nrow = length(samples), ncol = length(categories), dimnames = list(samples, categories))
  count_matrix[cbind(match(counts$sample, samples), match(as.character(counts[[category_column]]), categories))] <- counts$n_cells
  adjusted <- count_matrix + 0.5
  composition <- adjusted / rowSums(adjusted)
  log_composition <- log(composition)
  clr <- log_composition - rowMeans(log_composition)
  metadata <- unique(counts[c("sample", "Ploidy", "Dose", "sample_order")])
  metadata <- metadata[match(samples, metadata$sample), , drop = FALSE]
  observed_f <- multivariate_partial_f(clr, metadata$Dose, metadata$Ploidy)
  permuted_f <- vapply(seq_len(ncol(assignments)), function(j) {
    permuted_ploidy <- ifelse(assignments[, j], "2N", "4N")
    multivariate_partial_f(clr, metadata$Dose, permuted_ploidy)
  }, numeric(1))
  data.frame(
    resolution = category_column,
    n_samples = nrow(clr),
    n_categories = ncol(clr),
    pseudocount_cells = 0.5,
    transformation = "centered log-ratio",
    test = "exact Dose-stratified permutation test of multivariate partial F for Ploidy",
    observed_partial_f = observed_f,
    n_exact_stratified_assignments = ncol(assignments),
    p_value = mean(permuted_f >= observed_f - sqrt(.Machine$double.eps)),
    stringsAsFactors = FALSE
  )
}

write_plot_pair <- function(plot, stem, out_dir, width, height) {
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white", limitsize = FALSE)
  ggplot2::ggsave(file.path(out_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05q_mouse_celltype_composition_2N_vs_4N.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05p_latest_human_mouse_sample_matrix", "latest_human_mouse_sample_matrix.rds")
input_completion <- file.path(results_root, "05p_latest_human_mouse_sample_matrix", "completion.txt")
out_dir <- file.path(results_root, "05q_mouse_celltype_composition_2N_vs_4N")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds) || !file.exists(input_completion)) stop("05p current sample matrix is incomplete.", call. = FALSE)

dat <- readRDS(input_rds)
sample_metadata <- dat$sample_metadata
primary_counts <- dat$mouse_cell_type_counts
fine_counts <- dat$mouse_fine_category_counts
ploidy_counts <- table(factor(as.character(sample_metadata$Ploidy), levels = c("2N", "4N")))
if (nrow(sample_metadata) != 16L || !identical(as.integer(ploidy_counts), c(8L, 8L))) {
  stop("05q requires eight 2N and eight 4N biological samples.", call. = FALSE)
}
if (any(sample_metadata$sample %in% as.character(unlist(cfg$immune_composition$excluded_samples, use.names = FALSE)))) {
  stop("A cell-culture sample entered 05q.", call. = FALSE)
}

assignments <- build_stratified_ploidy_assignments(sample_metadata)
if (ncol(assignments) != 2520L) stop("Expected 2520 exact Dose-stratified Ploidy assignments.", call. = FALSE)
primary_tests <- run_category_tests(primary_counts, "mouse_cell_type_latest", assignments)
fine_tests <- run_category_tests(fine_counts, "mouse_cell_type_state_latest", assignments)
primary_summary <- composition_summary(primary_counts, "mouse_cell_type_latest")
fine_summary <- composition_summary(fine_counts, "mouse_cell_type_state_latest")
global_tests <- rbind(
  run_global_composition_test(primary_counts, "mouse_cell_type_latest", assignments),
  run_global_composition_test(fine_counts, "mouse_cell_type_state_latest", assignments)
)

primary_levels <- unique(primary_counts$mouse_cell_type_latest[order(primary_counts$category_order)])
fine_levels <- unique(fine_counts$mouse_cell_type_state_latest[order(fine_counts$category_order)])
primary_counts$Ploidy <- factor(as.character(primary_counts$Ploidy), levels = c("2N", "4N"))
primary_counts$mouse_cell_type_latest <- factor(primary_counts$mouse_cell_type_latest, levels = primary_levels)
fine_counts$Ploidy <- factor(as.character(fine_counts$Ploidy), levels = c("2N", "4N"))
fine_counts$mouse_cell_type_state_latest <- factor(fine_counts$mouse_cell_type_state_latest, levels = fine_levels)
primary_summary$category <- factor(primary_summary$category, levels = primary_levels)
fine_summary$category <- factor(fine_summary$category, levels = fine_levels)

ploidy_colors <- c(`2N` = "#3B6FB6", `4N` = "#E67E22")
broad_colors <- unlist(cfg$annotation$manual_broad_celltype_colors, use.names = TRUE)
immune_lineage_colors <- unlist(cfg$immune_fine_annotation$lineage_colors, use.names = TRUE)
primary_colors <- broad_colors
primary_colors["Neutrophil"] <- unname(immune_lineage_colors["Neutrophil"])
primary_colors["Macrophage"] <- unname(immune_lineage_colors["Macrophage"])
primary_colors["Classical monocyte"] <- unname(immune_lineage_colors["Classical monocyte"])
primary_colors["Ambiguous mixed myeloid"] <- unname(immune_lineage_colors["Ambiguous"])
for (label in primary_levels) {
  if (!label %in% names(primary_colors)) primary_colors[[label]] <- "#8D99AE"
}
primary_colors <- primary_colors[primary_levels]

plot_tests <- primary_tests
plot_tests$category <- factor(plot_tests$category, levels = primary_levels)
plot_tests$y <- vapply(primary_levels, function(label) {
  z <- primary_summary[as.character(primary_summary$category) == label, , drop = FALSE]
  max(z$ci95_upper, na.rm = TRUE) + 0.035
}, numeric(1))

base_theme <- ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"), legend.position = "top"
  )

p_grouped <- ggplot2::ggplot(primary_summary, ggplot2::aes(x = category, y = mean_sample_fraction, fill = Ploidy)) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.74), width = 0.68, color = "#343A40", linewidth = 0.25) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = ci95_lower, ymax = ci95_upper),
    position = ggplot2::position_dodge(width = 0.74), width = 0.16, linewidth = 0.45
  ) +
  ggplot2::geom_point(
    data = primary_counts,
    ggplot2::aes(x = mouse_cell_type_latest, y = fraction_within_mouse, shape = Ploidy),
    inherit.aes = FALSE,
    position = ggplot2::position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.74, seed = 20260911),
    color = "#111827", fill = "white", size = 1.9, stroke = 0.45, alpha = 0.85
  ) +
  ggplot2::geom_text(
    data = plot_tests,
    ggplot2::aes(x = category, y = y, label = significance_fdr),
    inherit.aes = FALSE, fontface = "bold", size = 3.5, color = "#111827"
  ) +
  ggplot2::scale_fill_manual(values = ploidy_colors, drop = FALSE) +
  ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
  ggplot2::scale_y_continuous(labels = function(x) percent_label(x, 0.01), expand = ggplot2::expansion(mult = c(0, 0.12))) +
  ggplot2::coord_cartesian(ylim = c(0, max(plot_tests$y, na.rm = TRUE) * 1.04), clip = "off") +
  ggplot2::labs(
    title = "Mouse host cell-type fractions: 2N versus 4N",
    subtitle = "Bars are equal-weight biological-sample means; points are samples; 95% t intervals; Dose-stratified exact permutation FDR",
    x = NULL, y = "Fraction of all retained mouse cells", fill = "Ploidy", shape = "Ploidy",
    caption = "Denominator is all retained mouse cells within each sample. Cell-culture samples are excluded."
  ) + base_theme +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 42, hjust = 1, vjust = 1))

stacked_summary <- primary_summary
stacked_summary$category <- factor(stacked_summary$category, levels = rev(primary_levels))
p_stacked <- ggplot2::ggplot(stacked_summary, ggplot2::aes(x = Ploidy, y = mean_sample_fraction, fill = category)) +
  ggplot2::geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  ggplot2::scale_fill_manual(values = primary_colors, breaks = primary_levels, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = percent_label, breaks = seq(0, 1, 0.2), expand = ggplot2::expansion(mult = c(0, 0.01))) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    title = "Mouse host cell-type composition by tumor-cell ploidy",
    subtitle = "Each segment is the arithmetic mean of per-sample fractions; every biological sample has equal weight",
    x = "Ploidy", y = "Mean fraction of all retained mouse cells", fill = "Mouse cell type",
    caption = "Composition is descriptive; the global and cell-type-specific tests control for Dose by restricted permutation."
  ) + base_theme

fine_plot_tests <- fine_tests
fine_plot_tests$category <- factor(fine_plot_tests$category, levels = fine_levels)
fine_plot_tests$mouse_cell_type_state_latest <- factor(as.character(fine_plot_tests$category), levels = fine_levels)
fine_plot_tests$y <- vapply(fine_levels, function(label) {
  z <- fine_counts[as.character(fine_counts$mouse_cell_type_state_latest) == label, , drop = FALSE]
  max(z$fraction_within_mouse, na.rm = TRUE) + max(0.003, diff(range(z$fraction_within_mouse)) * 0.12)
}, numeric(1))
p_fine <- ggplot2::ggplot(fine_counts, ggplot2::aes(x = Ploidy, y = fraction_within_mouse, shape = Ploidy, fill = Ploidy)) +
  ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = 20260911), size = 2.2, color = "#343A40", alpha = 0.9) +
  ggplot2::stat_summary(ggplot2::aes(group = Ploidy), fun = mean, geom = "point", shape = 95, size = 7, color = "#111827") +
  ggplot2::geom_text(
    data = fine_plot_tests, ggplot2::aes(x = 1.5, y = y, label = significance_fdr),
    inherit.aes = FALSE, fontface = "bold", size = 3.1, color = "#111827"
  ) +
  ggplot2::facet_wrap(~mouse_cell_type_state_latest, scales = "free_y", ncol = 4) +
  ggplot2::scale_shape_manual(values = c(`2N` = 21, `4N` = 22), drop = FALSE) +
  ggplot2::scale_fill_manual(values = ploidy_colors, drop = FALSE) +
  ggplot2::scale_y_continuous(labels = function(x) percent_label(x, 0.001), expand = ggplot2::expansion(mult = c(0.04, 0.16))) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::labs(
    title = "Fine mouse cell-type/state fractions: 2N versus 4N",
    subtitle = "Each point is one biological sample; horizontal ticks show means; labels use BH-adjusted Dose-stratified exact permutation p-values",
    x = "Ploidy", y = "Fraction of all retained mouse cells",
    caption = "Fine state labels remain existing-cluster-level annotations. Free facet scales are used only to make rare categories visible."
  ) + base_theme + ggplot2::guides(shape = "none", fill = "none")

utils::write.csv(primary_counts, file.path(out_dir, "mouse_cell_type_sample_fractions.csv"), row.names = FALSE)
utils::write.csv(primary_summary, file.path(out_dir, "mouse_cell_type_ploidy_equal_sample_summary.csv"), row.names = FALSE)
utils::write.csv(primary_tests, file.path(out_dir, "mouse_cell_type_2N_vs_4N_tests.csv"), row.names = FALSE, na = "NA")
utils::write.csv(fine_counts, file.path(out_dir, "mouse_fine_category_sample_fractions.csv"), row.names = FALSE)
utils::write.csv(fine_summary, file.path(out_dir, "mouse_fine_category_ploidy_equal_sample_summary.csv"), row.names = FALSE)
utils::write.csv(fine_tests, file.path(out_dir, "mouse_fine_category_2N_vs_4N_tests.csv"), row.names = FALSE, na = "NA")
utils::write.csv(global_tests, file.path(out_dir, "global_mouse_composition_2N_vs_4N_tests.csv"), row.names = FALSE)
utils::write.csv(data.frame(mouse_cell_type = names(primary_colors), color_hex = unname(primary_colors)), file.path(out_dir, "mouse_cell_type_palette.csv"), row.names = FALSE)
write_plot_pair(p_grouped, "mouse_cell_type_2N_vs_4N_grouped_bar", out_dir, 13.5, 7.8)
write_plot_pair(p_stacked, "mouse_cell_type_composition_2N_vs_4N_stacked", out_dir, 8.5, 7.2)
write_plot_pair(p_fine, "mouse_fine_category_2N_vs_4N_sample_distribution", out_dir, 15, 12)

contract <- list(
  question = "Do current cluster-level mouse host cell-type fractions differ between tumors initiated from 2N and 4N cells?",
  biological_unit = "sample",
  numerator = "retained current 05h mouse cells assigned to one cluster-level category",
  denominator = "all retained current 05h mouse cells within the same sample",
  primary_resolution = "mouse_cell_type_latest: broad non-immune types plus manually approved 05h immune lineages",
  supplementary_resolution = "mouse_cell_type_state_latest: mutually exclusive manually approved fine states where available",
  equal_sample_weight = TRUE,
  ploidy_sample_sizes = list(`2N` = 8L, `4N` = 8L),
  dose_adjustment = "Ploidy labels are permuted only within Dose blocks (0, 30, 120 mg/kg)",
  exact_assignments = ncol(assignments),
  global_test = "CLR multivariate partial-F test for Ploidy conditional on Dose",
  category_test = "two-sided exact restricted permutation test of the difference in mean sample fractions",
  multiple_testing = "Benjamini-Hochberg across all categories separately within each annotation resolution",
  sensitivity_model = "quasibinomial cell-type-versus-other count model with Ploidy and Dose",
  cell_level_classifier_used = FALSE,
  interpretation = "captured-cell composition, not absolute histologic abundance"
)
yaml::write_yaml(contract, file.path(out_dir, "analysis_contract.yaml"))

input_files <- c(input_rds, input_completion, normalizePath(config_path, mustWork = TRUE), script_path)
input_manifest <- data.frame(
  role = c("05p_current_sample_matrix", "05p_completion", "analysis_config", "analysis_script"),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05q.txt"))
writeLines(c(
  "05q sample-level 2N-versus-4N mouse composition analysis completed.",
  paste0("Biological samples: 16; 2N=8; 4N=8; exact Dose-stratified assignments=", ncol(assignments), "."),
  paste0("Primary mouse cell types tested: ", nrow(primary_tests), "; BH-FDR<0.05: ", sum(primary_tests$p_value_bh < 0.05), "."),
  paste0("Fine mutually exclusive categories tested: ", nrow(fine_tests), "; BH-FDR<0.05: ", sum(fine_tests$p_value_bh < 0.05), "."),
  paste0("Global primary-composition p-value: ", signif(global_tests$p_value[global_tests$resolution == "mouse_cell_type_latest"], 6), "."),
  "Cells were used only to calculate within-sample counts; statistical replication is at the biological-sample level.",
  "Dose was controlled by restricted permutation, and cell-culture samples were excluded."
), file.path(out_dir, "completion.txt"))

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files), size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
message("[05q] Completed. Output: ", out_dir)
