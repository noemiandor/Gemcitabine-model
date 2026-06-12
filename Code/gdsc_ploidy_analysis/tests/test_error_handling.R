#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("--file=", args, value = TRUE)
  if (length(hit) > 0) {
    return(dirname(normalizePath(sub("--file=", "", hit[1], fixed = TRUE))))
  }
  normalizePath(getwd())
}

tests_dir <- script_dir()
base_dir <- normalizePath(file.path(tests_dir, ".."))
source(file.path(base_dir, "src", "analysis_helpers.R"))

expect_error_contains <- function(expr, expected, label) {
  msg <- tryCatch(
    {
      force(expr)
      NA_character_
    },
    error = function(e) conditionMessage(e)
  )
  if (is.na(msg)) {
    stop(label, ": expected an error but none was thrown", call. = FALSE)
  }
  missing <- expected[!vapply(expected, grepl, logical(1), x = msg, fixed = TRUE)]
  if (length(missing) > 0) {
    stop(label, ": missing expected text in error message: ", paste(missing, collapse = ", "), "\nObserved: ", msg, call. = FALSE)
  }
}

annotations <- data.frame(
  drug = c("A", "B"),
  group = c("SIGNALING", "CYTOTOXIC"),
  stringsAsFactors = FALSE
)
values <- matrix(c(0.1, -0.2), dimnames = list(c("A", "B"), "TEST"))

expect_error_contains(
  run_enrichment_or_stop(
    values,
    annotations,
    cancer = "TEST",
    direction = "low_ploidy_sensitive",
    permute_n = 1,
    pvalue_cutoff = 0.05,
    enrichment_fn = function(...) stop("forced enrichment failure")
  ),
  c("cancer=TEST", "direction=low_ploidy_sensitive", "forced enrichment failure"),
  "enrichment context"
)

expect_error_contains(
  plot_barplot_or_stop(
    c(A = 0.1),
    colors = c("red"),
    cancer = "TEST",
    xlab = "test metric",
    plotting_fn = function(...) stop("forced plot failure")
  ),
  c("cancer=TEST", "plotted_drugs=1", "forced plot failure"),
  "plot context"
)

plot_page_table <- build_ploidy_sensitivity_plot_table(
  cancer_type = "TEST",
  plot_values = c(DrugA = -0.2, DrugB = 0.3),
  group_by_drug = c(DrugA = "CYTOTOXIC", DrugB = "SIGNALING"),
  color_by_group = c(CYTOTOXIC = "#111111", SIGNALING = "#222222"),
  metric = "Z_SCORE",
  metric_label = "GDSC Z-score",
  abs_r_threshold = 0.1,
  page_index = 1,
  page_tsv_file = "tables/ploidyVsDrugSensitivity_pages/ploidyVsDrugSensitivity_page_01_TEST.tsv"
)
if (!identical(plot_page_table$barplot_order, 1:2)) {
  stop("Ploidy sensitivity page table should preserve barplot order", call. = FALSE)
}
if (!identical(plot_page_table$plot_color, c("#111111", "#222222"))) {
  stop("Ploidy sensitivity page table should record group colors", call. = FALSE)
}

expect_error_contains(
  validate_required_groups(c("SIGNALING"), required_groups = c("SIGNALING", "CYTOTOXIC")),
  c("Missing required annotation groups", "CYTOTOXIC"),
  "required group validation"
)

dups <- data.frame(
  DRUG_NAME = c("DrugA", "DrugA", "DrugB"),
  CELL_LINE_NAME = c("Cell1", "Cell1", "Cell2"),
  Z_SCORE = c(10, 20, 30),
  RMSE = c(0.5, 0.1, 0.2),
  NLME_RESULT_ID = c(2, 1, 3),
  stringsAsFactors = FALSE
)
resolved <- resolve_duplicate_drug_cell_lines(dups, metric = "Z_SCORE", strategy = "lowest_rmse")
if (nrow(resolved) != 2) {
  stop("duplicate resolution should return two rows", call. = FALSE)
}
if (resolved$Z_SCORE[resolved$DRUG_NAME == "DrugA"] != 20) {
  stop("lowest_rmse duplicate resolution did not select the lower-RMSE row", call. = FALSE)
}

dups_no_rmse <- dups[, setdiff(colnames(dups), "RMSE"), drop = FALSE]
expect_error_contains(
  resolve_duplicate_drug_cell_lines(dups_no_rmse, metric = "Z_SCORE", strategy = "lowest_rmse"),
  c("missing required columns", "RMSE"),
  "lowest_rmse requires RMSE"
)

dup_qc_dir <- tempfile()
resolved_with_qc <- resolve_duplicate_drug_cell_lines(dups, metric = "Z_SCORE", strategy = "lowest_rmse", qc_dir = dup_qc_dir)
selection_audit <- read.table(
  file.path(dup_qc_dir, "duplicate_drug_cell_line_selection.tsv"),
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!any(selection_audit$selected)) {
  stop("Duplicate selection audit should mark selected rows", call. = FALSE)
}
if (!identical(nrow(resolved_with_qc), 2L)) {
  stop("Duplicate resolution with QC should still return resolved rows", call. = FALSE)
}

cor_stats <- compute_drug_ploidy_correlation(
  response = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  ploidy = c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29),
  min_n = 5
)
if (!identical(cor_stats$spearman_ci_method, "approximate_fisher_transform")) {
  stop("Spearman CI method should be labeled as approximate_fisher_transform", call. = FALSE)
}
if (!is.finite(cor_stats$pearson_r)) {
  stop("Correlation helper should return a finite Pearson r for valid input", call. = FALSE)
}

low_n_stats <- compute_drug_ploidy_correlation(
  response = c(1, 2, 3),
  ploidy = c(1, 2, 3),
  min_n = 5
)
if (!identical(low_n_stats$spearman_ci_method, "not_estimated_n_below_minimum")) {
  stop("Low-n Spearman CI method should be not_estimated_n_below_minimum", call. = FALSE)
}

gem_summary_path <- tempfile(fileext = ".tsv")
gem_summary <- write_gemcitabine_rank_summary(
  data.frame(
    cancer_type = c("A", "A"),
    drug = c("Gemcitabine", "OtherDrug"),
    metric = c("Z_SCORE", "Z_SCORE"),
    n = c(10, 10),
    stringsAsFactors = FALSE
  ),
  gem_summary_path
)
if (nrow(gem_summary) != 1 || gem_summary$drug != "Gemcitabine") {
  stop("Gemcitabine summary should neutrally subset gemcitabine aliases only", call. = FALSE)
}

tissue_dr <- data.frame(
  CELL_LINE_NAME = paste0("CELL", seq_len(12)),
  DRUG_NAME = rep("DrugA", 12),
  TCGA_DESC = rep(c("T1", "T2", "T3"), each = 4),
  Z_SCORE = seq_len(12),
  stringsAsFactors = FALSE
)
tissue_ploidy <- data.frame(
  `Cell iname` = paste0("CELL", seq_len(12)),
  ploidy = seq(2, 3.1, length.out = 12),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
tissue_models <- build_tissue_adjusted_ploidy_models(
  tissue_dr,
  tissue_ploidy,
  metrics = "Z_SCORE",
  min_n = 10,
  min_tissues = 3,
  min_rows_per_tissue = 2
)
if (!identical(tissue_models$model_status, "ok")) {
  stop("Tissue-adjusted model fixture should fit with model_status ok", call. = FALSE)
}

cat("Error handling tests passed.\n")
