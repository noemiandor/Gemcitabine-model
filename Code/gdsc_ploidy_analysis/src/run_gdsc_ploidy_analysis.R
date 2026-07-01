options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

module_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  hit <- grep(needle, args, value = TRUE)
  if (length(hit) > 0) {
    script_path <- dirname(normalizePath(sub(needle, "", hit[1])))
    if (basename(script_path) == "src") {
      return(dirname(script_path))
    }
    return(script_path)
  }
  normalizePath(getwd())
}

arg_value <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0("--", name, "=")
  hit <- grep(prefix, args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[1], fixed = TRUE)
}

base_dir <- module_dir()
src_dir <- file.path(base_dir, "src")
source(file.path(src_dir, "dependencies.R"))
setup_local_lib(base_dir)
require_packages(c("EnrichIntersect", "openxlsx", "data.table"))
suppressPackageStartupMessages({
  library(EnrichIntersect)
})
source(file.path(src_dir, "analysis_helpers.R"))
source(file.path(src_dir, "drug_class_workbook.R"))

find_python_for_plotting <- function() {
  python <- Sys.which("python3")
  if (!nzchar(python)) {
    python <- Sys.which("python")
  }
  if (!nzchar(python)) {
    stop("Could not find python3 or python to generate enrichment heatmaps.", call. = FALSE)
  }
  unname(python)
}

run_python_plot_script <- function(script, args, description) {
  if (!file.exists(script)) {
    stop("Missing plotting script for ", description, ": ", script, call. = FALSE)
  }
  python <- find_python_for_plotting()
  message("Generating ", description, " with ", basename(script))
  status <- system2(
    python,
    args = c(normalizePath(script), args),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L)) {
    stop("Failed to generate ", description, " with exit status ", status, ".", call. = FALSE)
  }
}

run_enrichment_heatmap_plots <- function(workbook, output_dir, source_dir, category_mode = "primary_secondary") {
  if (!file.exists(workbook)) {
    stop("Cannot generate enrichment heatmaps because workbook is missing: ", workbook, call. = FALSE)
  }
  suffix <- paste0("_", category_mode)
  run_python_plot_script(
    file.path(source_dir, "plot_ploidy_enrichment_panels.py"),
    c(
      normalizePath(workbook),
      "--out-prefix",
      file.path(output_dir, paste0("ploidy_enrichment_panels", suffix, "_ABC"))
    ),
    "manuscript-style ploidy-enrichment panels"
  )
  run_python_plot_script(
    file.path(source_dir, "plot_ploidy_enrichment_clustered_heatmaps.py"),
    c(
      normalizePath(workbook),
      "--out-prefix",
      file.path(output_dir, paste0("ploidy_enrichment_clustered", suffix))
    ),
    "clustered ploidy-enrichment heatmaps"
  )
}

data_dir <- file.path(base_dir, "data")
raw_data_dir <- file.path(data_dir, "raw")
repo_root <- normalizePath(file.path(base_dir, "..", ".."), mustWork = TRUE)
default_out_dir <- file.path(
  repo_root,
  "Results",
  "public_data",
  "gdsc_ploidy_analysis",
  "runs",
  paste0(format(Sys.time(), "%Y%m%dT%H%M%S"), "_gdsc")
)
out_dir <- normalizePath(arg_value("output-dir", default_out_dir), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gdsc_file <- file.path(raw_data_dir, "GDSC2_fitted_dose_response_24Jul22.txt")
ploidy_file <- file.path(raw_data_dir, "ploidyAcrossCellLines_V1.txt")
drug_class_workbook_file <- normalizePath(
  arg_value("drug-class-workbook", file.path(data_dir, "manual", "drug_class_final_used_with_primary_secondary_corrected.xlsx")),
  mustWork = FALSE
)
category_mode <- arg_value("category-mode", "primary_secondary")
if (!identical(category_mode, "primary_secondary")) {
  stop("Only --category-mode=primary_secondary is supported. Deprecated values curated, legacy, and proposal were removed.", call. = FALSE)
}

required_input_files <- c(gdsc_file, ploidy_file, drug_class_workbook_file)
if (!all(file.exists(required_input_files))) {
  stop("Missing required input files: ", paste(required_input_files[!file.exists(required_input_files)], collapse = ", "), call. = FALSE)
}
write_input_manifest(
  unique(required_input_files),
  file.path(out_dir, "metadata", "input_manifest.tsv")
)

dr <- read.table(gdsc_file, sep = "\t", header = TRUE)
dr$CELL_LINE_NAME_RAW <- dr$CELL_LINE_NAME
dr$CELL_LINE_NAME <- toupper(gsub("-", "", dr$CELL_LINE_NAME))
dr$CELL_LINE_KEY <- normalize_cell_line_name(dr$CELL_LINE_NAME_RAW)

appCL <- read.table(ploidy_file, sep = "\t", check.names = FALSE, header = TRUE)
appCL <- appCL[!is.na(appCL$ploidy), ]
appCL$CELL_LINE_KEY <- normalize_cell_line_name(appCL$`Cell iname`)

R <- list()
canonical_metric <- "Z_SCORE"
metric <- canonical_metric
metric_label <- "GDSC Z-score"
supplemental_metrics <- intersect(c("Z_SCORE", "LN_IC50", "AUC"), colnames(dr))
canonical_matching_mode <- "raw_name_matching"
supplemental_matching_mode <- "normalized_cell_line_key"
ploidy_collision_tolerance <- as.numeric(arg_value("ploidy-collision-tolerance", "0.05"))
correlation_min_n <- as.integer(arg_value("correlation-min-n", "10"))
analysis_mode <- arg_value("analysis-mode", "dev")
if (!analysis_mode %in% c("dev", "manuscript")) {
  stop("--analysis-mode must be either 'dev' or 'manuscript'.", call. = FALSE)
}
default_enrichment_permute_n <- 300L
enrichment_permute_n <- as.integer(arg_value("enrichment-permute-n", as.character(default_enrichment_permute_n)))
low_ploidy_pvalue_cutoff <- as.numeric(arg_value("low-ploidy-pvalue-cutoff", "0.05"))
high_ploidy_pvalue_cutoff <- as.numeric(arg_value("high-ploidy-pvalue-cutoff", "0.1"))
tissue_model_min_n <- as.integer(arg_value("tissue-model-min-n", "20"))
tissue_model_min_tissues <- as.integer(arg_value("tissue-model-min-tissues", "3"))
tissue_model_min_rows_per_tissue <- as.integer(arg_value("tissue-model-min-rows-per-tissue", "2"))
ploidy_sensitivity_plot_abs_r_threshold <- as.numeric(arg_value("ploidy-sensitivity-plot-abs-r-threshold", "0.185"))

metadata_dir <- file.path(out_dir, "metadata")
tables_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

gdsc_collisions <- audit_normalized_key_collisions(
  dr,
  raw_col = "CELL_LINE_NAME_RAW",
  key_col = "CELL_LINE_KEY",
  source_name = "GDSC"
)
ploidy_collisions <- audit_normalized_key_collisions(
  appCL,
  raw_col = "Cell iname",
  key_col = "CELL_LINE_KEY",
  source_name = "CellPassports_ploidy"
)
write_tsv(as.data.frame(gdsc_collisions), file.path(tables_dir, "cell_line_key_collisions_gdsc.tsv"))
write_tsv(as.data.frame(ploidy_collisions), file.path(tables_dir, "cell_line_key_collisions_ploidy.tsv"))
validate_ploidy_key_collisions(
  appCL,
  tolerance = ploidy_collision_tolerance,
  tables_dir = tables_dir
)
write_cell_line_matching_delta(
  dr,
  appCL,
  file.path(tables_dir, "cell_line_matching_delta_raw_vs_normalized.tsv")
)

appCL <- appCL[!duplicated(appCL$`Cell iname`), ]
rownames(appCL) <- appCL$`Cell iname`
duplicate_strategy <- "lowest_rmse"
duplicate_strategy_required_columns <- paste(
  c("DRUG_NAME", "CELL_LINE_NAME", metric, if (duplicate_strategy == "lowest_rmse") "RMSE" else character()),
  collapse = ";"
)
write_run_metadata(
  data.frame(
    key = c(
      "metric",
      "metric_label",
      "canonical_metric",
      "supplemental_metrics",
      "canonical_matching_mode",
      "supplemental_matching_mode",
      "ploidy_collision_tolerance",
      "correlation_min_n",
      "spearman_ci_method",
      "analysis_mode",
      "category_mode",
      "drug_class_workbook_file",
      "drug_class_workbook_md5",
      "drug_class_workbook_sha256",
      "enrichment_permute_n",
      "low_ploidy_pvalue_cutoff",
      "high_ploidy_pvalue_cutoff",
      "tissue_model_min_n",
      "tissue_model_min_tissues",
      "tissue_model_min_rows_per_tissue",
      "ploidy_sensitivity_plot_abs_r_threshold",
      "duplicate_strategy",
      "duplicate_strategy_required_columns"
    ),
    value = c(
      metric,
      metric_label,
      canonical_metric,
      paste(supplemental_metrics, collapse = ";"),
      canonical_matching_mode,
      supplemental_matching_mode,
      as.character(ploidy_collision_tolerance),
      as.character(correlation_min_n),
      "approximate_fisher_transform_for_estimable_spearman_ci",
      analysis_mode,
      category_mode,
      drug_class_workbook_file,
      file_checksum(drug_class_workbook_file),
      if (file.exists(drug_class_workbook_file)) unname(tools::sha256sum(drug_class_workbook_file)) else NA_character_,
      as.character(enrichment_permute_n),
      as.character(low_ploidy_pvalue_cutoff),
      as.character(high_ploidy_pvalue_cutoff),
      as.character(tissue_model_min_n),
      as.character(tissue_model_min_tissues),
      as.character(tissue_model_min_rows_per_tissue),
      as.character(ploidy_sensitivity_plot_abs_r_threshold),
      duplicate_strategy,
      duplicate_strategy_required_columns
    ),
    stringsAsFactors = FALSE
  ),
  metadata_dir
)
dr <- resolve_duplicate_drug_cell_lines(dr, metric = metric, strategy = duplicate_strategy, qc_dir = qc_dir)

canonical_correlations <- build_drug_ploidy_correlation_table(
  dr,
  appCL,
  metrics = canonical_metric,
  matching_mode = "raw",
  min_n = correlation_min_n
)
write_tsv(
  canonical_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_Z_SCORE.tsv")
)
write_correlations_xlsx(
  canonical_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_Z_SCORE.xlsx")
)

all_metric_correlations <- build_drug_ploidy_correlation_table(
  dr,
  appCL,
  metrics = supplemental_metrics,
  matching_mode = "raw",
  min_n = correlation_min_n
)
write_tsv(
  all_metric_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_all_metrics.tsv")
)
write_correlations_xlsx(
  all_metric_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_all_metrics.xlsx")
)

normalized_key_correlations <- build_drug_ploidy_correlation_table(
  dr,
  appCL,
  metrics = canonical_metric,
  matching_mode = "normalized",
  min_n = correlation_min_n
)
write_tsv(
  normalized_key_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.tsv")
)
write_correlations_xlsx(
  normalized_key_correlations,
  file.path(tables_dir, "drug_ploidy_correlations_by_cancer_Z_SCORE_normalized_key_SUPPLEMENTAL.xlsx")
)
write_correlation_delta(
  canonical_correlations,
  normalized_key_correlations,
  file.path(tables_dir, "correlation_delta_raw_vs_normalized_Z_SCORE.tsv")
)
write_gemcitabine_rank_summary(
  canonical_correlations,
  file.path(tables_dir, "gemcitabine_rank_summary_Z_SCORE.tsv")
)
write_gemcitabine_rank_summary(
  all_metric_correlations,
  file.path(tables_dir, "gemcitabine_rank_summary_all_metrics.tsv")
)

tissue_models_z <- build_tissue_adjusted_ploidy_models(
  dr,
  appCL,
  metrics = canonical_metric,
  min_n = tissue_model_min_n,
  min_tissues = tissue_model_min_tissues,
  min_rows_per_tissue = tissue_model_min_rows_per_tissue
)
write_tsv(
  tissue_models_z,
  file.path(tables_dir, "all_cancers_tissue_adjusted_ploidy_models_Z_SCORE.tsv")
)
tissue_models_all <- build_tissue_adjusted_ploidy_models(
  dr,
  appCL,
  metrics = supplemental_metrics,
  min_n = tissue_model_min_n,
  min_tissues = tissue_model_min_tissues,
  min_rows_per_tissue = tissue_model_min_rows_per_tissue
)
write_tsv(
  tissue_models_all,
  file.path(tables_dir, "all_cancers_tissue_adjusted_ploidy_models_all_metrics.tsv")
)

pdf(file.path(out_dir, "drugsVsPloidyCorr.pdf"), width = 15, height = 7)
par(mfrow = c(3, 7))
for (can in c("allcancers", unique(dr$TCGA_DESC))) {
  dr_sub <- dr
  if (can != "allcancers") {
    dr_sub <- dr[dr$TCGA_DESC == can, ]
  }
  ii <- intersect(dr_sub$CELL_LINE_NAME, rownames(appCL))
  if (length(ii) < 10) {
    next
  }

  r <- list()
  for (drug in unique(dr_sub$DRUG_NAME)) {
    dr_drug <- dr_sub[dr_sub$DRUG_NAME == drug, ]
    rownames(dr_drug) <- dr_drug$CELL_LINE_NAME
    if (sum(!is.na(dr_drug[ii, metric])) < 10) {
      next
    }
    r[[drug]] <- cor(dr_drug[ii, metric], appCL[ii, "ploidy"], use = "pairwise.complete.obs")
    if (abs(r[[drug]]) > 0.6) {
      plot(
        dr_drug[ii, metric],
        appCL[ii, "ploidy"],
        main = paste(can, drug),
        xlab = metric_label,
        ylab = "Ploidy"
      )
    }
  }
  R[[can]] <- sort(unlist(r))
}
dev.off()
save(R, file = file.path(out_dir, "drugsVsPloidyCorr.RData"))

correlation_drugs <- derive_correlation_eligible_drugs(R, dr)
write_tsv(correlation_drugs, file.path(tables_dir, "drug_class_correlation_eligible_drugs.tsv"))

primary_secondary_rows <- read_primary_secondary_drug_classes(drug_class_workbook_file)
validation_rows <- validate_primary_secondary_class_coverage(correlation_drugs, primary_secondary_rows)
write_tsv(validation_rows, file.path(tables_dir, "drug_class_workbook_validation.tsv"))
write_tsv(
  primary_secondary_used_table(primary_secondary_rows),
  file.path(tables_dir, "drug_class_primary_secondary_used.tsv")
)
primary_secondary_counts <- primary_secondary_class_counts(primary_secondary_rows)
write_tsv(
  primary_secondary_counts,
  file.path(tables_dir, "drug_class_primary_secondary_counts.tsv")
)
coxIn <- make_enrichment_class_table(primary_secondary_rows)
group_order <- attr(primary_secondary_rows, "class_order")
category_suffix <- "primary_secondary"
category_source_label <- "reviewed_primary_anticancer_class_workbook"

save(coxIn, file = file.path(out_dir, "coxIn.RData"))

for (can in names(R)) {
  names(R[[can]]) <- toupper(names(R[[can]]))
}
R_ <- sapply(R, function(x) x[names(x) %in% coxIn$drug])

x <- sapply(names(R_), function(can) matrix(R_[[can]], dimnames = list(names(R_[[can]]), can)))
lowpIsSens <- highpIsSens <- list()
selected_drugs_for_enrichment <- list()
enrichment_metadata <- list()
category_id_by_drug <- if ("category_id" %in% colnames(coxIn)) {
  setNames(coxIn$category_id, rownames(coxIn))
} else {
  setNames(coxIn$group, rownames(coxIn))
}
coxIn_for_enrichment <- coxIn[, c("drug", "group"), drop = FALSE]
# Stabilize the permutation-based enrichment step for reproducibility testing.
set.seed(1)
for (can in names(x)) {
  enrichment_drugs <- rownames(x[[can]])
  enrichment_values <- as.numeric(x[[can]][, 1])
  selected_drugs_for_enrichment[[paste(can, "low_ploidy_sensitive", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "low_ploidy_sensitive",
    metric = metric,
    drug = enrichment_drugs,
    correlation_value = enrichment_values,
    group = coxIn[enrichment_drugs, "group"],
    category_id = unname(category_id_by_drug[enrichment_drugs]),
    stringsAsFactors = FALSE
  )
  selected_drugs_for_enrichment[[paste(can, "high_ploidy_sensitive", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "high_ploidy_sensitive",
    metric = metric,
    drug = enrichment_drugs,
    correlation_value = -enrichment_values,
    group = coxIn[enrichment_drugs, "group"],
    category_id = unname(category_id_by_drug[enrichment_drugs]),
    stringsAsFactors = FALSE
  )
  enrichment_metadata[[paste(can, "low_ploidy_sensitive", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "low_ploidy_sensitive",
    metric = metric,
    analysis_mode = analysis_mode,
    permute_n = enrichment_permute_n,
    pvalue_cutoff = low_ploidy_pvalue_cutoff,
    input_drugs = length(enrichment_drugs),
    input_groups = length(unique(coxIn[enrichment_drugs, "group"])),
    category_source = category_source_label,
    category_mode = category_mode,
    stringsAsFactors = FALSE
  )
  enrichment_metadata[[paste(can, "high_ploidy_sensitive", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "high_ploidy_sensitive",
    metric = metric,
    analysis_mode = analysis_mode,
    permute_n = enrichment_permute_n,
    pvalue_cutoff = high_ploidy_pvalue_cutoff,
    input_drugs = length(enrichment_drugs),
    input_groups = length(unique(coxIn[enrichment_drugs, "group"])),
    category_source = category_source_label,
    category_mode = category_mode,
    stringsAsFactors = FALSE
  )

  lowpIsSens[[can]] <- run_enrichment_or_stop(
    x[[can]],
    coxIn_for_enrichment,
    cancer = can,
    direction = "low_ploidy_sensitive",
    permute_n = enrichment_permute_n,
    pvalue_cutoff = low_ploidy_pvalue_cutoff
  )
  highpIsSens[[can]] <- run_enrichment_or_stop(
    -x[[can]],
    coxIn_for_enrichment,
    cancer = can,
    direction = "high_ploidy_sensitive",
    permute_n = enrichment_permute_n,
    pvalue_cutoff = high_ploidy_pvalue_cutoff
  )
}

write_tsv(
  do.call(rbind, selected_drugs_for_enrichment),
  file.path(tables_dir, sprintf("class_enrichment_selected_drugs_%s_%s.tsv", category_suffix, metric))
)
write_tsv(
  do.call(rbind, enrichment_metadata),
  file.path(tables_dir, sprintf("class_enrichment_%s_metadata.tsv", category_suffix))
)

if (!is.null(group_order)) {
  groups <- group_order
} else {
  groups <- sort(unique(coxIn$group))
}
groups <- groups[groups %in% unique(coxIn$group)]
lowpIsSens <- sapply(lowpIsSens, function(x) as.data.frame(x)[groups, ])
highpIsSens <- sapply(highpIsSens, function(x) as.data.frame(x)[groups, ])
rownames(lowpIsSens) <- rownames(highpIsSens) <- groups

enrichment_long <- rbind(
  data.frame(
    direction = "low_ploidy_sensitive",
    metric = metric,
    cancer_type = rep(colnames(lowpIsSens), each = nrow(lowpIsSens)),
    group = rep(rownames(lowpIsSens), times = ncol(lowpIsSens)),
    pvalue = as.vector(lowpIsSens),
    stringsAsFactors = FALSE
  ),
  data.frame(
    direction = "high_ploidy_sensitive",
    metric = metric,
    cancer_type = rep(colnames(highpIsSens), each = nrow(highpIsSens)),
    group = rep(rownames(highpIsSens), times = ncol(highpIsSens)),
    pvalue = as.vector(highpIsSens),
    stringsAsFactors = FALSE
  )
)
write_tsv(enrichment_long, file.path(tables_dir, sprintf("class_enrichment_%s_%s.tsv", category_suffix, metric)))
if (category_mode == "primary_secondary" && exists("primary_secondary_counts")) {
  write_tsv(
    primary_secondary_review_summary(
      primary_secondary_counts,
      enrichment_long,
      low_ploidy_pvalue_cutoff = low_ploidy_pvalue_cutoff,
      high_ploidy_pvalue_cutoff = high_ploidy_pvalue_cutoff
    ),
    file.path(tables_dir, "drug_class_primary_secondary_review_summary.tsv")
  )
}

mode_workbook <- file.path(out_dir, sprintf("drugsVsPloidyCorr_%s_%s.xlsx", category_suffix, metric))
mode_wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(mode_wb, "lowpIsSens")
openxlsx::writeData(mode_wb, "lowpIsSens", t(lowpIsSens), rowNames = TRUE)
openxlsx::addWorksheet(mode_wb, "highpIsSens")
openxlsx::writeData(mode_wb, "highpIsSens", t(highpIsSens), rowNames = TRUE)
openxlsx::saveWorkbook(mode_wb, mode_workbook, overwrite = TRUE)
invisible(file.copy(mode_workbook, file.path(out_dir, "drugsVsPloidyCorr.xlsx"), overwrite = TRUE))
invisible(file.copy(
  mode_workbook,
  file.path(tables_dir, sprintf("drugsVsPloidyCorr_%s_%s.xlsx", category_suffix, metric)),
  overwrite = TRUE
))
write_tsv(
  data.frame(
    category_mode = category_mode,
    category_source = category_source_label,
    canonical_workbook = mode_workbook,
    compatibility_workbook = file.path(out_dir, "drugsVsPloidyCorr.xlsx"),
    canonical_workbook_md5 = file_checksum(mode_workbook),
    drug_class_workbook_file = drug_class_workbook_file,
    drug_class_workbook_md5 = file_checksum(drug_class_workbook_file),
    drug_class_workbook_sha256 = if (file.exists(drug_class_workbook_file)) unname(tools::sha256sum(drug_class_workbook_file)) else NA_character_,
    stringsAsFactors = FALSE
  ),
  file.path(metadata_dir, "category_mode_artifacts.tsv")
)

tmp <- sort(unique(coxIn$group))
col <- rainbow(length(tmp) * 1.3)[1:length(tmp)]
names(col) <- tmp
ploidy_sensitivity_pages_dir <- file.path(tables_dir, "ploidyVsDrugSensitivity_pages")
dir.create(ploidy_sensitivity_pages_dir, recursive = TRUE, showWarnings = FALSE)
ploidy_sensitivity_plot_tables <- list()
pdf(file.path(out_dir, "ploidyVsDrugSensitivity.pdf"), width = 3, height = 6)
for (page_index in seq_along(colnames(lowpIsSens))) {
  sheet <- colnames(lowpIsSens)[page_index]
  plot_vals <- R_[[sheet]][abs(R_[[sheet]]) >= ploidy_sensitivity_plot_abs_r_threshold]
  page_tsv_name <- sprintf(
    "ploidyVsDrugSensitivity_page_%02d_%s.tsv",
    page_index,
    safe_file_stem(sheet)
  )
  page_tsv_rel <- file.path("tables", "ploidyVsDrugSensitivity_pages", page_tsv_name)
  page_table <- build_ploidy_sensitivity_plot_table(
    cancer_type = sheet,
    plot_values = plot_vals,
    group_by_drug = setNames(coxIn$group, rownames(coxIn)),
    color_by_group = col,
    metric = metric,
    metric_label = metric_label,
    abs_r_threshold = ploidy_sensitivity_plot_abs_r_threshold,
    page_index = page_index,
    page_tsv_file = page_tsv_rel
  )
  write_tsv(page_table, file.path(ploidy_sensitivity_pages_dir, page_tsv_name))
  ploidy_sensitivity_plot_tables[[sheet]] <- page_table
  plot_barplot_or_stop(
    plot_vals,
    colors = col[coxIn[names(plot_vals), "group"]],
    cancer = sheet,
    xlab = paste("Pearson r between ploidy and drug sensitivity", paste0("(", metric_label, ")"))
  )
}
dev.off()
write_tsv(
  do.call(rbind, ploidy_sensitivity_plot_tables),
  file.path(tables_dir, "ploidyVsDrugSensitivity_plot_values_Z_SCORE.tsv")
)

run_enrichment_heatmap_plots(
  workbook = mode_workbook,
  output_dir = out_dir,
  source_dir = src_dir,
  category_mode = category_mode
)

write_session_metadata(file.path(metadata_dir, "session_info.txt"))
message("Analysis completed. Outputs written to: ", out_dir)
