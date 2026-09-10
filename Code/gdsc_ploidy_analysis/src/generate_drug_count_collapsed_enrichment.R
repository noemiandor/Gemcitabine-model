#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(paste0("^", file_arg), args)])
if (length(script_path) == 0) {
  script_path <- normalizePath("Code/gdsc_ploidy_analysis/src/generate_drug_count_collapsed_enrichment.R")
}
src_dir <- dirname(normalizePath(script_path))
module_dir <- normalizePath(file.path(src_dir, ".."))
repo_root <- normalizePath(file.path(module_dir, "..", ".."))

source(file.path(src_dir, "dependencies.R"))
setup_local_lib(module_dir)
require_packages(c("openxlsx"))
source(file.path(src_dir, "analysis_helpers.R"))
source(file.path(src_dir, "drug_class_workbook.R"))

arg_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  hit <- grep(paste0("^", prefix), commandArgs(trailingOnly = TRUE), value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(prefix, "", hit[[1]])
}

input_run <- normalizePath(arg_value(
  "input-run",
  file.path(repo_root, "Results", "public_data", "gdsc_ploidy_analysis", "runs", "post_cleanup_check_gdsc")
))
drug_class_workbook <- normalizePath(arg_value(
  "drug-class-workbook",
  file.path(repo_root, "Data", "public_data", "gdsc", "manual", "drug_class_final_used_with_primary_secondary_corrected.xlsx")
))
output_dir <- arg_value(
  "output-dir",
  file.path(input_run, "drug_count_collapsed")
)
if (!grepl("^/", output_dir)) {
  output_dir <- file.path(repo_root, output_dir)
}
output_dir <- normalizePath(output_dir, mustWork = FALSE)

metric <- arg_value("metric", "Z_SCORE")
permute_n <- as.integer(arg_value("permute-n", "1000"))
significance_cutoff <- as.numeric(arg_value("significance-cutoff", "0.05"))
min_group_contexts <- as.integer(arg_value("min-group-contexts", "2"))
max_bidirectional_imbalance <- as.integer(arg_value("max-bidirectional-imbalance", "2"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
tables_dir <- file.path(output_dir, "tables")
figures_dir <- file.path(output_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(
  tables_dir,
  c(
    sprintf("drug_count_collapsed_direction_summary_p%g.tsv", significance_cutoff),
    sprintf("drug_count_collapsed_category_groups_p%g.tsv", significance_cutoff)
  )
))

rdata_file <- file.path(input_run, "drugsVsPloidyCorr.RData")
if (!file.exists(rdata_file)) {
  stop("Missing correlation RData: ", rdata_file, call. = FALSE)
}
load(rdata_file)
if (!exists("R")) {
  stop("Expected object `R` in ", rdata_file, call. = FALSE)
}

class_rows <- read_primary_secondary_drug_classes(drug_class_workbook)
drug_counts <- openxlsx::read.xlsx(drug_class_workbook, sheet = "Drug counts", detectDates = FALSE)
required_count_cols <- c("primary_anticancer_class", "count", "Combine.if.needed")
missing_count_cols <- setdiff(required_count_cols, colnames(drug_counts))
if (length(missing_count_cols) > 0) {
  stop("Drug counts sheet is missing required column(s): ", paste(missing_count_cols, collapse = ", "), call. = FALSE)
}
drug_counts$primary_anticancer_class <- trimws(as.character(drug_counts$primary_anticancer_class))
drug_counts$Combine.if.needed <- trimws(as.character(drug_counts$Combine.if.needed))
drug_counts$Combine.if.needed[is.na(drug_counts$Combine.if.needed)] <- ""
drug_counts$collapsed_class <- ifelse(
  nonempty(drug_counts$Combine.if.needed),
  drug_counts$Combine.if.needed,
  drug_counts$primary_anticancer_class
)

collapse_map <- drug_counts[, c("primary_anticancer_class", "count", "Combine.if.needed", "collapsed_class")]
collapse_map$collapse_rule <- ifelse(
  nonempty(collapse_map$Combine.if.needed),
  "Drug counts sheet Combine.if.needed",
  "unchanged primary_anticancer_class"
)
write_tsv(collapse_map, file.path(tables_dir, "drug_count_collapsed_class_map.tsv"))

included <- class_rows[class_rows$include_in_enrichment, , drop = FALSE]
missing_primary <- setdiff(sort(unique(included$primary_anticancer_class)), collapse_map$primary_anticancer_class)
if (length(missing_primary) > 0) {
  stop(
    "Included primary classes missing from Drug counts collapse map: ",
    paste(missing_primary, collapse = ", "),
    call. = FALSE
  )
}
included$collapsed_class <- collapse_map$collapsed_class[
  match(included$primary_anticancer_class, collapse_map$primary_anticancer_class)
]

coxIn <- data.frame(
  drug = normalize_excel_drug_key(included$drug_key),
  group = trimws(included$collapsed_class),
  primary_anticancer_class = trimws(included$primary_anticancer_class),
  stringsAsFactors = FALSE
)
rownames(coxIn) <- coxIn$drug

for (can in names(R)) {
  names(R[[can]]) <- normalize_drug_key(names(R[[can]]))
}
R_ <- lapply(R, function(x) x[names(x) %in% rownames(coxIn)])
empty_contexts <- names(R_)[vapply(R_, length, integer(1)) == 0]
if (length(empty_contexts) > 0) {
  stop("No collapsed-class-annotated drugs for context(s): ", paste(empty_contexts, collapse = ", "), call. = FALSE)
}

collapsed_order <- unique(collapse_map$collapsed_class)
collapsed_order <- collapsed_order[collapsed_order %in% coxIn$group]
collapsed_counts <- aggregate(drug ~ group, coxIn, length)
names(collapsed_counts) <- c("category_label", "n_drugs")
collapsed_counts$source_primary_classes <- vapply(collapsed_counts$category_label, function(label) {
  paste(sort(unique(coxIn$primary_anticancer_class[coxIn$group == label])), collapse = "; ")
}, character(1))
collapsed_counts$.order <- match(collapsed_counts$category_label, collapsed_order)
collapsed_counts <- collapsed_counts[order(collapsed_counts$.order, collapsed_counts$category_label), , drop = FALSE]
collapsed_counts$.order <- NULL
write_tsv(collapsed_counts, file.path(tables_dir, "drug_count_collapsed_class_counts.tsv"))

assignment_table <- data.frame(
  drug = included$drug,
  drug_key = normalize_excel_drug_key(included$drug_key),
  primary_anticancer_class = included$primary_anticancer_class,
  category_label = included$collapsed_class,
  secondary_mechanism_note = included$secondary_mechanism_note,
  stringsAsFactors = FALSE
)
assignment_table <- assignment_table[order(assignment_table$category_label, assignment_table$primary_anticancer_class, assignment_table$drug_key), , drop = FALSE]
write_tsv(assignment_table, file.path(tables_dir, "drug_count_collapsed_assignments.tsv"))

one_sided_ks_enrichment <- function(in_group) {
  in_group <- as.logical(in_group)
  n_in <- sum(in_group)
  n_out <- length(in_group) - n_in
  if (n_in == 0 || n_out == 0) {
    return(NA_real_)
  }
  max(cumsum(in_group) / n_in - cumsum(!in_group) / n_out, na.rm = TRUE)
}

collapsed_enrichment_pvalues <- function(values, annotations, groups, permute_n) {
  values <- values[is.finite(values)]
  values <- values[names(values) %in% annotations$drug]
  values <- sort(values, decreasing = TRUE)
  drug_order <- names(values)
  group_by_drug <- annotations$group[match(drug_order, annotations$drug)]
  n <- length(values)

  out <- rep(NA_real_, length(groups))
  names(out) <- groups
  for (group in groups) {
    in_group <- group_by_drug == group
    observed <- one_sided_ks_enrichment(in_group)
    if (!is.finite(observed)) {
      next
    }
    permuted <- replicate(permute_n, {
      one_sided_ks_enrichment(in_group[sample.int(n, n)])
    })
    out[[group]] <- sum(permuted >= observed, na.rm = TRUE) / sum(is.finite(permuted))
  }
  out
}

groups <- collapsed_counts$category_label
lowp <- highp <- list()
selected_drugs <- list()
set.seed(1)
for (can in names(R_)) {
  values <- R_[[can]]
  annotations <- coxIn[names(values), c("drug", "group", "primary_anticancer_class"), drop = FALSE]
  selected_drugs[[paste(can, "low", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "low_ploidy_sensitive",
    drug = names(values),
    correlation_value = as.numeric(values),
    category_label = annotations$group,
    primary_anticancer_class = annotations$primary_anticancer_class,
    stringsAsFactors = FALSE
  )
  selected_drugs[[paste(can, "high", sep = "::")]] <- data.frame(
    cancer_type = can,
    direction = "high_ploidy_sensitive",
    drug = names(values),
    correlation_value = -as.numeric(values),
    category_label = annotations$group,
    primary_anticancer_class = annotations$primary_anticancer_class,
    stringsAsFactors = FALSE
  )
  lowp[[can]] <- collapsed_enrichment_pvalues(values, annotations, groups, permute_n)
  highp[[can]] <- collapsed_enrichment_pvalues(-values, annotations, groups, permute_n)
}
write_tsv(
  do.call(rbind, selected_drugs),
  file.path(tables_dir, sprintf("drug_count_collapsed_selected_drugs_%s.tsv", metric))
)

lowp <- sapply(lowp, function(x) x[groups])
highp <- sapply(highp, function(x) x[groups])
rownames(lowp) <- rownames(highp) <- groups

enrichment_long <- rbind(
  data.frame(
    direction = "low_ploidy_sensitive",
    metric = metric,
    cancer_type = rep(colnames(lowp), each = nrow(lowp)),
    category_label = rep(rownames(lowp), times = ncol(lowp)),
    pvalue = as.vector(lowp),
    stringsAsFactors = FALSE
  ),
  data.frame(
    direction = "high_ploidy_sensitive",
    metric = metric,
    cancer_type = rep(colnames(highp), each = nrow(highp)),
    category_label = rep(rownames(highp), times = ncol(highp)),
    pvalue = as.vector(highp),
    stringsAsFactors = FALSE
  )
)
zero_pvalue_floor <- 1 / (permute_n + 1)
enrichment_long$pvalue_for_fdr <- ifelse(
  is.finite(enrichment_long$pvalue) & enrichment_long$pvalue <= 0,
  zero_pvalue_floor,
  enrichment_long$pvalue
)
enrichment_long$qvalue_bh <- p.adjust(enrichment_long$pvalue_for_fdr, method = "BH")
write_tsv(
  enrichment_long,
  file.path(tables_dir, sprintf("drug_count_collapsed_enrichment_%s.tsv", metric))
)

fdr_matrix <- function(direction) {
  subset <- enrichment_long[enrichment_long$direction == direction, , drop = FALSE]
  out <- matrix(
    NA_real_,
    nrow = length(groups),
    ncol = length(colnames(lowp)),
    dimnames = list(groups, colnames(lowp))
  )
  for (idx in seq_len(nrow(subset))) {
    out[subset$category_label[[idx]], subset$cancer_type[[idx]]] <- subset$qvalue_bh[[idx]]
  }
  out
}
lowq <- fdr_matrix("low_ploidy_sensitive")
highq <- fdr_matrix("high_ploidy_sensitive")

summary_rows <- lapply(groups, function(group) {
  low <- enrichment_long[enrichment_long$direction == "low_ploidy_sensitive" & enrichment_long$category_label == group, ]
  high <- enrichment_long[enrichment_long$direction == "high_ploidy_sensitive" & enrichment_long$category_label == group, ]
  merged <- merge(
    low[, c("cancer_type", "pvalue", "qvalue_bh")],
    high[, c("cancer_type", "pvalue", "qvalue_bh")],
    by = "cancer_type",
    suffixes = c("_low", "_high")
  )
  low_sig <- merged$qvalue_bh_low <= significance_cutoff
  high_sig <- merged$qvalue_bh_high <= significance_cutoff
  data.frame(
    category_label = group,
    n_drugs = collapsed_counts$n_drugs[match(group, collapsed_counts$category_label)],
    source_primary_classes = collapsed_counts$source_primary_classes[match(group, collapsed_counts$category_label)],
    n_low_contexts = sum(low_sig, na.rm = TRUE),
    n_high_contexts = sum(high_sig, na.rm = TRUE),
    n_low_only_contexts = sum(low_sig & !high_sig, na.rm = TRUE),
    n_high_only_contexts = sum(high_sig & !low_sig, na.rm = TRUE),
    n_both_contexts = sum(low_sig & high_sig, na.rm = TRUE),
    min_low_pvalue = min(low$pvalue, na.rm = TRUE),
    min_high_pvalue = min(high$pvalue, na.rm = TRUE),
    min_low_qvalue_bh = min(low$qvalue_bh, na.rm = TRUE),
    min_high_qvalue_bh = min(high$qvalue_bh, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
direction_summary <- do.call(rbind, summary_rows)
direction_summary$direction_bias_contexts <- direction_summary$n_high_contexts - direction_summary$n_low_contexts
direction_summary$total_significant_contexts <- direction_summary$n_low_contexts + direction_summary$n_high_contexts
direction_summary$category_group <- ifelse(
  direction_summary$total_significant_contexts < min_group_contexts,
  "weak_or_not_recurrent",
  ifelse(
    direction_summary$n_low_contexts > 0 &
      direction_summary$n_high_contexts > 0 &
      abs(direction_summary$direction_bias_contexts) <= max_bidirectional_imbalance,
    "bidirectional_context_dependent",
    ifelse(
      direction_summary$n_low_contexts > direction_summary$n_high_contexts,
      "low_ploidy_biased",
      ifelse(
        direction_summary$n_high_contexts > direction_summary$n_low_contexts,
        "high_ploidy_biased",
        "weak_or_not_recurrent"
      )
    )
  )
)
direction_summary$category_group_definition <- ifelse(
  direction_summary$category_group == "low_ploidy_biased",
  sprintf(">= %d total BH-FDR significant contexts and more low- than high-ploidy contexts after the bidirectional rule", min_group_contexts),
  ifelse(
    direction_summary$category_group == "high_ploidy_biased",
    sprintf(">= %d total BH-FDR significant contexts and more high- than low-ploidy contexts after the bidirectional rule", min_group_contexts),
    ifelse(
      direction_summary$category_group == "bidirectional_context_dependent",
      sprintf(">= %d total BH-FDR significant contexts, at least one context in each direction, and absolute high-minus-low context difference <= %d", min_group_contexts, max_bidirectional_imbalance),
      sprintf("< %d total BH-FDR significant low/high enrichment contexts", min_group_contexts)
    )
  )
)
direction_summary <- direction_summary[order(
  direction_summary$category_group,
  -abs(direction_summary$direction_bias_contexts),
  direction_summary$category_label
), , drop = FALSE]
write_tsv(
  direction_summary,
  file.path(tables_dir, sprintf("drug_count_collapsed_direction_summary_fdr%g.tsv", significance_cutoff))
)
write_tsv(
  direction_summary[, c(
    "category_label",
    "category_group",
    "category_group_definition",
    "source_primary_classes",
    "n_drugs",
    "n_low_contexts",
    "n_high_contexts",
    "total_significant_contexts",
    "n_low_only_contexts",
    "n_high_only_contexts",
    "n_both_contexts",
    "min_low_pvalue",
    "min_high_pvalue",
    "min_low_qvalue_bh",
    "min_high_qvalue_bh"
  )],
  file.path(tables_dir, sprintf("drug_count_collapsed_category_groups_fdr%g.tsv", significance_cutoff))
)

cancer_type_count_file <- file.path(input_run, "tables", sprintf("drug_ploidy_correlations_by_cancer_%s.tsv", metric))
cancer_type_counts <- data.frame(
  cancer_type = colnames(lowp),
  n_cell_lines = NA_integer_,
  count_method = "maximum matched cell lines across retained drug correlations",
  source_file = cancer_type_count_file,
  stringsAsFactors = FALSE
)
if (file.exists(cancer_type_count_file)) {
  correlation_counts <- read.delim(cancer_type_count_file, sep = "\t", header = TRUE, check.names = FALSE)
  if (all(c("cancer_type", "n") %in% colnames(correlation_counts))) {
    if ("metric" %in% colnames(correlation_counts)) {
      correlation_counts <- correlation_counts[correlation_counts$metric == metric, , drop = FALSE]
    }
    correlation_counts$n <- suppressWarnings(as.numeric(correlation_counts$n))
    correlation_counts <- correlation_counts[is.finite(correlation_counts$n), , drop = FALSE]
    if (nrow(correlation_counts) > 0) {
      max_counts <- aggregate(n ~ cancer_type, correlation_counts, max, na.rm = TRUE)
      names(max_counts)[names(max_counts) == "n"] <- "n_cell_lines"
      cancer_type_counts$n_cell_lines <- max_counts$n_cell_lines[
        match(cancer_type_counts$cancer_type, max_counts$cancer_type)
      ]
    }
  } else {
    warning("Cancer-type count source is missing cancer_type/n columns: ", cancer_type_count_file, call. = FALSE)
  }
} else {
  warning("Cancer-type count source missing: ", cancer_type_count_file, call. = FALSE)
}
write_tsv(
  cancer_type_counts,
  file.path(tables_dir, sprintf("drug_count_collapsed_cancer_type_cell_line_counts_%s.tsv", metric))
)

workbook <- file.path(output_dir, sprintf("drug_count_collapsed_drugsVsPloidyCorr_%s.xlsx", metric))
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "rawLowpIsSens")
openxlsx::writeData(wb, "rawLowpIsSens", t(lowp), rowNames = TRUE)
openxlsx::addWorksheet(wb, "rawHighpIsSens")
openxlsx::writeData(wb, "rawHighpIsSens", t(highp), rowNames = TRUE)
openxlsx::addWorksheet(wb, "lowpIsSens")
openxlsx::writeData(wb, "lowpIsSens", t(lowq), rowNames = TRUE)
openxlsx::addWorksheet(wb, "highpIsSens")
openxlsx::writeData(wb, "highpIsSens", t(highq), rowNames = TRUE)
openxlsx::addWorksheet(wb, "drugClassCounts")
openxlsx::writeData(wb, "drugClassCounts", collapsed_counts)
openxlsx::addWorksheet(wb, "cancerTypeCounts")
openxlsx::writeData(wb, "cancerTypeCounts", cancer_type_counts)
openxlsx::saveWorkbook(wb, workbook, overwrite = TRUE)

plot_script <- file.path(src_dir, "plot_ploidy_enrichment_clustered_heatmaps.py")
plot_prefix <- file.path(figures_dir, "drug_count_collapsed_enrichment")
cmd <- sprintf(
  "MPLCONFIGDIR=%s python3 %s %s --out-prefix %s",
  shQuote(file.path(tempdir(), "mplconfig_gdsc_collapsed_classes")),
  shQuote(plot_script),
  shQuote(workbook),
  shQuote(plot_prefix)
)
status <- system(cmd)
if (status != 0) {
  stop("Drug-count collapsed heatmap plotting failed with status ", status, call. = FALSE)
}

write_tsv(
  data.frame(
    input_run = input_run,
    drug_class_workbook = drug_class_workbook,
    rdata_file = rdata_file,
    metric = metric,
    permute_n = permute_n,
    significance_cutoff = significance_cutoff,
    significance_measure = "BH FDR q-value; exact zero permutation p-values are floored at 1/(permute_n + 1) before BH correction",
    zero_pvalue_floor = zero_pvalue_floor,
    min_group_contexts = min_group_contexts,
    max_bidirectional_imbalance = max_bidirectional_imbalance,
    collapse_rule = "Use Drug counts sheet Combine.if.needed where populated; otherwise primary_anticancer_class",
    workbook = workbook,
    plot_prefix = plot_prefix,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "drug_count_collapsed_run_config.tsv")
)

cat("Wrote drug-count collapsed enrichment outputs to:", output_dir, "\n")
