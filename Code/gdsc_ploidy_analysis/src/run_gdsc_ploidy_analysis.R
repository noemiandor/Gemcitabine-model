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
dr$CELL_LINE_NAME <- toupper(gsub("-", "", dr$CELL_LINE_NAME))

appCL <- read.table(ploidy_file, sep = "\t", check.names = FALSE, header = TRUE)
appCL <- appCL[!is.na(appCL$ploidy), ]
appCL <- appCL[!duplicated(appCL$`Cell iname`), ]
rownames(appCL) <- appCL$`Cell iname`

R <- list()
metric <- "Z_SCORE"
metric_label <- "GDSC Z-score"

metadata_dir <- file.path(out_dir, "metadata")
tables_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")
dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
duplicate_strategy <- "lowest_rmse"
write_run_metadata(
  data.frame(
    key = c("metric", "metric_label", "duplicate_strategy"),
    value = c(metric, metric_label, duplicate_strategy),
    stringsAsFactors = FALSE
  ),
  metadata_dir
)
dr <- resolve_duplicate_drug_cell_lines(dr, metric = metric, strategy = duplicate_strategy, qc_dir = qc_dir)

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

coxIn <- coxIn[, intersect(c("drug", "drugName", "group"), colnames(coxIn))]
coxIn <- coxIn[!duplicated(coxIn$drug), ]
rownames(coxIn) <- toupper(coxIn$drug)

save(coxIn, file = file.path(out_dir, "coxIn.RData"))

coxIn_other <- coxIn[is.na(coxIn$group), , drop = FALSE]
coxIn_other$group <- "NOTCLASSIFIED"
coxIn <- coxIn[!is.na(coxIn$group), , drop = FALSE]

tmp <- strsplit(coxIn$group, "; ", fixed = TRUE)
coxIn$group <- vapply(tmp, function(x) x[length(x)], character(1))
coxIn$group <- gsub("Cytotoxic medicines", "Cytotoxic", gsub(";", "", gsub(",", "", coxIn$group)))
coxIn$group[grep("Alkylating", coxIn$group)] <- "Alkylating"
coxIn$group[grep("Topoisomerase", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Tubulin", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Antimitotic", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Antineoplastic Agents", coxIn$group)] <- "Antineoplastic Agents"
coxIn$group <- toupper(coxIn$group)
coxIn$group <- gsub("(ANTI-)INFLAMMATORY", "IMMUNOSUPPRESSIVE AGENTS", coxIn$group, fixed = TRUE)
coxIn$group[coxIn$group %in% c("PARP INHIBITORS", "SIGNAL TRANSDUCTION INHIBITORS", "JAK INHIBITORS", "ENZYME INHIBITORS", "TARGETED THERAPIES")] <- "SIGNALING"
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
# Stabilize the permutation-based enrichment step for reproducibility testing.
set.seed(1)
for (can in names(x)) {
  lowpIsSens[[can]] <- run_enrichment_or_stop(
    x[[can]],
    coxIn,
    cancer = can,
    direction = "low_ploidy_sensitive",
    permute_n = 300,
    pvalue_cutoff = 0.05
  )
  highpIsSens[[can]] <- run_enrichment_or_stop(
    -x[[can]],
    coxIn,
    cancer = can,
    direction = "high_ploidy_sensitive",
    permute_n = 300,
    pvalue_cutoff = 0.1
  )
}

groups <- unique(coxIn$group)
lowpIsSens <- sapply(lowpIsSens, function(x) as.data.frame(x)[groups, ])
highpIsSens <- sapply(highpIsSens, function(x) as.data.frame(x)[groups, ])
rownames(lowpIsSens) <- rownames(highpIsSens) <- groups
lowpIsSens <- lowpIsSens[, order(lowpIsSens["SIGNALING", ])]
lowpIsSens <- lowpIsSens[, order(lowpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["SIGNALING", ])]

write.xlsx(t(lowpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr.xlsx"), sheetName = "lowpIsSens")
write.xlsx(t(highpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr.xlsx"), sheetName = "highpIsSens", append = TRUE)

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
