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
require_packages(c("EnrichIntersect", "xlsx", "plyr", "openxlsx", "data.table"))
suppressPackageStartupMessages({
  library(EnrichIntersect)
  library(xlsx)
  library(plyr)
})
source(file.path(src_dir, "annotations.R"))
source(file.path(src_dir, "analysis_helpers.R"))
data_dir <- file.path(base_dir, "data")
raw_data_dir <- file.path(data_dir, "raw")
out_dir <- normalizePath(arg_value("output-dir", file.path(base_dir, "output")), mustWork = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

gdsc_file <- file.path(raw_data_dir, "GDSC2_fitted_dose_response_24Jul22.txt")
ploidy_file <- file.path(raw_data_dir, "ploidyAcrossCellLines_V1.txt")
cmap_file <- file.path(raw_data_dir, "small_molecule_20200407234909.csv")
custom_file <- file.path(data_dir, "manual", "custom_set_candidate.tsv")
annotation_file <- normalizePath(
  arg_value("annotation-file", file.path(data_dir, "derived", "pubchem_drug_annotations.tsv")),
  mustWork = FALSE
)

stopifnot(file.exists(gdsc_file), file.exists(ploidy_file), file.exists(cmap_file), file.exists(custom_file), file.exists(annotation_file))
write_input_manifest(
  c(gdsc_file, ploidy_file, cmap_file, custom_file, annotation_file),
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
canonical_matching_mode <- "legacy_raw_name_matching"
supplemental_matching_mode <- "normalized_cell_line_key"
ploidy_collision_tolerance <- as.numeric(arg_value("ploidy-collision-tolerance", "0.05"))
correlation_min_n <- as.integer(arg_value("correlation-min-n", "10"))
analysis_mode <- arg_value("analysis-mode", "dev")
if (!analysis_mode %in% c("dev", "manuscript")) {
  stop("--analysis-mode must be either 'dev' or 'manuscript'.", call. = FALSE)
}
default_enrichment_permute_n <- if (analysis_mode == "manuscript") 10000L else 300L
enrichment_permute_n <- as.integer(arg_value("enrichment-permute-n", as.character(default_enrichment_permute_n)))
low_ploidy_pvalue_cutoff <- as.numeric(arg_value("low-ploidy-pvalue-cutoff", "0.05"))
high_ploidy_pvalue_cutoff <- as.numeric(arg_value("high-ploidy-pvalue-cutoff", "0.1"))
tissue_model_min_n <- as.integer(arg_value("tissue-model-min-n", "20"))
tissue_model_min_tissues <- as.integer(arg_value("tissue-model-min-tissues", "3"))
tissue_model_min_rows_per_tissue <- as.integer(arg_value("tissue-model-min-rows-per-tissue", "2"))

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
      "enrichment_permute_n",
      "low_ploidy_pvalue_cutoff",
      "high_ploidy_pvalue_cutoff",
      "tissue_model_min_n",
      "tissue_model_min_tissues",
      "tissue_model_min_rows_per_tissue",
      "duplicate_strategy"
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
      as.character(enrichment_permute_n),
      as.character(low_ploidy_pvalue_cutoff),
      as.character(high_ploidy_pvalue_cutoff),
      as.character(tissue_model_min_n),
      as.character(tissue_model_min_tissues),
      as.character(tissue_model_min_rows_per_tissue),
      duplicate_strategy
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

coxIn <- data.frame(drug = unique(unlist(sapply(R, names))), stringsAsFactors = FALSE)
coxIn <- load_drug_annotations(coxIn$drug, annotation_file, custom_file)
coxIn$group <- coxIn$drugCategory_Pubchem
annotation_audit <- build_drug_annotation_audit(coxIn, custom_file)
write_tsv(annotation_audit, file.path(tables_dir, "drug_annotation_audit.tsv"))
write_curated_mapping_template(
  annotation_audit,
  file.path(tables_dir, "drug_class_final_curated_TEMPLATE.tsv")
)
warning("Enrichment still uses legacy_group_used_for_enrichment. Review drug_annotation_audit.tsv before switching to curated final_category.")

coxIn <- coxIn[, intersect(c("drug", "drugName", "group"), colnames(coxIn))]
coxIn <- coxIn[!duplicated(coxIn$drug), ]
rownames(coxIn) <- toupper(coxIn$drug)

save(coxIn, file = file.path(out_dir, "coxIn.RData"))

coxIn_other <- coxIn[is.na(coxIn$group), , drop = FALSE]
coxIn_other$group <- "NOTCLASSIFIED"
coxIn <- coxIn[!is.na(coxIn$group), , drop = FALSE]

coxIn$group <- normalize_legacy_drug_group(coxIn$group)
coxIn <- coxIn[nchar(coxIn$group) > 0, , drop = FALSE]

fr <- plyr::count(coxIn$group)
write_tsv(fr, file.path(tables_dir, "drug_category_counts_before_filter.tsv"))
coxIn <- coxIn[coxIn$group %in% fr$x[fr$freq > 1], , drop = FALSE]
fr_after <- plyr::count(coxIn$group)
write_tsv(fr_after, file.path(tables_dir, "drug_category_counts_after_filter.tsv"))
validate_required_groups(coxIn$group)
coxIn <- coxIn[, c("drug", "group"), drop = FALSE]
coxIn$drug <- toupper(coxIn$drug)
rownames(coxIn) <- coxIn$drug

for (can in names(R)) {
  names(R[[can]]) <- toupper(names(R[[can]]))
}
R_ <- sapply(R, function(x) x[names(x) %in% coxIn$drug])

x <- sapply(names(R_), function(can) matrix(R_[[can]], dimnames = list(names(R_[[can]]), can)))
lowpIsSens <- highpIsSens <- list()
selected_drugs_for_enrichment <- list()
enrichment_metadata <- list()
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
    stringsAsFactors = FALSE
  )
  selected_drugs_for_enrichment[[paste(can, "high_ploidy_sensitive", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "high_ploidy_sensitive",
    metric = metric,
    drug = enrichment_drugs,
    correlation_value = -enrichment_values,
    group = coxIn[enrichment_drugs, "group"],
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
    category_source = "legacy_group_used_for_enrichment",
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
    category_source = "legacy_group_used_for_enrichment",
    stringsAsFactors = FALSE
  )

  lowpIsSens[[can]] <- run_enrichment_or_stop(
    x[[can]],
    coxIn,
    cancer = can,
    direction = "low_ploidy_sensitive",
    permute_n = enrichment_permute_n,
    pvalue_cutoff = low_ploidy_pvalue_cutoff
  )
  highpIsSens[[can]] <- run_enrichment_or_stop(
    -x[[can]],
    coxIn,
    cancer = can,
    direction = "high_ploidy_sensitive",
    permute_n = enrichment_permute_n,
    pvalue_cutoff = high_ploidy_pvalue_cutoff
  )
}

write_tsv(
  do.call(rbind, selected_drugs_for_enrichment),
  file.path(tables_dir, "class_enrichment_selected_drugs_legacy_Z_SCORE.tsv")
)
write_tsv(
  do.call(rbind, enrichment_metadata),
  file.path(tables_dir, "class_enrichment_legacy_metadata.tsv")
)

groups <- unique(coxIn$group)
lowpIsSens <- sapply(lowpIsSens, function(x) as.data.frame(x)[groups, ])
highpIsSens <- sapply(highpIsSens, function(x) as.data.frame(x)[groups, ])
rownames(lowpIsSens) <- rownames(highpIsSens) <- groups
lowpIsSens <- lowpIsSens[, order(lowpIsSens["SIGNALING", ])]
lowpIsSens <- lowpIsSens[, order(lowpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["SIGNALING", ])]

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
write_tsv(enrichment_long, file.path(tables_dir, "class_enrichment_legacy_Z_SCORE.tsv"))

write.xlsx(t(lowpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr.xlsx"), sheetName = "lowpIsSens")
write.xlsx(t(highpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr.xlsx"), sheetName = "highpIsSens", append = TRUE)
invisible(file.copy(
  file.path(out_dir, "drugsVsPloidyCorr.xlsx"),
  file.path(tables_dir, "drugsVsPloidyCorr_legacy_Z_SCORE.xlsx"),
  overwrite = TRUE
))

tmp <- sort(unique(coxIn$group))
col <- rainbow(length(tmp) * 1.3)[1:length(tmp)]
names(col) <- tmp
pdf(file.path(out_dir, "ploidyVsDrugSensitivity.pdf"), width = 3, height = 6)
for (sheet in colnames(lowpIsSens)) {
  plot_vals <- R_[[sheet]][abs(R_[[sheet]]) >= 0.1]
  plot_barplot_or_stop(
    plot_vals,
    colors = col[coxIn[names(plot_vals), "group"]],
    cancer = sheet,
    xlab = paste("Pearson r between ploidy and drug sensitivity", paste0("(", metric_label, ")"))
  )
}
dev.off()

write_session_metadata(file.path(metadata_dir, "session_info.txt"))
message("Analysis completed. Outputs written to: ", out_dir)
