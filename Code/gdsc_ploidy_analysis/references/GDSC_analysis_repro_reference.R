options(java.parameters = "-Xmx7g")

library("EnrichIntersect")
library(xlsx)
library(matlab)
library(plyr)

dr <- read.table(
  "/Users/4470246/Projects/Resources/data/Databases/GDSC/GDSC2_fitted_dose_response_24Jul22.txt",
  sep = "\t",
  header = TRUE
)
dr$CELL_LINE_NAME <- toupper(gsub("-", "", dr$CELL_LINE_NAME))

appCL <- read.table(
  "/Users/4470246/Projects/Resources/data/ploidyAcrossCellLines/ploidyAcrossCellLines_V1.txt",
  sep = "\t",
  check.names = FALSE,
  header = TRUE
)
appCL <- appCL[!is.na(appCL$ploidy), ]
appCL <- appCL[!duplicated(appCL$`Cell iname`), ]
rownames(appCL) <- appCL$`Cell iname`

out_dir <- path.expand("~/Downloads")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

metric <- "Z_SCORE"
R <- list()

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
    dr_drug <- dr_drug[!duplicated(dr_drug$CELL_LINE_NAME), ]
    rownames(dr_drug) <- dr_drug$CELL_LINE_NAME
    if (sum(!is.na(dr_drug[ii, metric])) < 10) {
      next
    }
    r[[drug]] <- cor(dr_drug[ii, metric], appCL[ii, "ploidy"], use = "pairwise.complete.obs")
    if (abs(r[[drug]]) > 0.6) {
      plot(dr_drug[ii, metric], appCL[ii, "ploidy"], main = paste(can, drug))
    }
  }

  R[[can]] <- sort(unlist(r))
}
dev.off()

save(file = file.path(out_dir, "drugsVsPloidyCorr.RObj"), "R")

source("/Users/4470246/Projects/code/RCode/scripts/annotateFromDrugBank.R")
source("/Users/4470246/Projects/code/RCode/scripts/annotateFromPubchem.R")

coxIn <- data.frame(
  drug = unique(unlist(sapply(R, names))),
  stringsAsFactors = FALSE
)
coxIn$drugName <- coxIn$drug
coxIn <- annotateFromPubchem(coxIn)

if ("drugCategory_Pubchem" %in% colnames(coxIn)) {
  coxIn$group <- coxIn$drugCategory_Pubchem
} else if ("group" %in% colnames(coxIn)) {
  coxIn$group <- coxIn$group
} else {
  stop("Could not identify annotated drug category column.")
}

save(file = file.path(out_dir, "coxIn.RObj"), "coxIn")

if (!exists("custom.set")) {
  custom.set <- data.frame(drug = character(), group = character(), stringsAsFactors = FALSE)
}

load(file = file.path(out_dir, "drugsVsPloidyCorr.RObj"))
load(file = file.path(out_dir, "coxIn.RObj"))

coxIn <- coxIn[, intersect(c("drug", "group"), colnames(coxIn))]
coxIn <- coxIn[!duplicated(coxIn$drug), ]
rownames(coxIn) <- toupper(coxIn$drug)

if (nrow(custom.set) > 0) {
  rownames(custom.set) <- toupper(custom.set$drug)
  ii <- intersect(rownames(custom.set), rownames(coxIn[is.na(coxIn$group), , drop = FALSE]))
  coxIn[ii, "group"] <- custom.set[ii, "group"]
}

coxIn_other <- coxIn[is.na(coxIn$group), , drop = FALSE]
coxIn_other$group <- "NOTCLASSIFIED"
coxIn <- coxIn[!is.na(coxIn$group), , drop = FALSE]

tmp <- sapply(coxIn$group, function(x) strsplit(x, "; ")[[1]])
coxIn$group <- sapply(tmp, function(x) x[length(x)])
coxIn$group <- gsub("Cytotoxic medicines", "Cytotoxic", gsub(";", "", gsub(",", "", coxIn$group)))
coxIn$group[grep("Alkylating", coxIn$group)] <- "Alkylating"
coxIn$group[grep("Topoisomerase", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Tubulin", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Antimitotic", coxIn$group)] <- "Cytotoxic"
coxIn$group[grep("Antineoplastic Agents", coxIn$group)] <- "Antineoplastic Agents"
coxIn$group <- toupper(coxIn$group)
coxIn$group <- gsub("(ANTI-)INFLAMMATORY", "IMMUNOSUPPRESSIVE AGENTS", coxIn$group, fixed = TRUE)
coxIn$group[coxIn$group %in% c(
  "PARP INHIBITORS",
  "SIGNAL TRANSDUCTION INHIBITORS",
  "JAK INHIBITORS",
  "ENZYME INHIBITORS",
  "TARGETED THERAPIES"
)] <- "SIGNALING"
coxIn <- coxIn[nchar(coxIn$group) > 0, , drop = FALSE]

fr <- plyr::count(coxIn$group)
coxIn <- coxIn[coxIn$group %in% fr$x[fr$freq > 1], , drop = FALSE]
print(fr)

coxIn$drug <- toupper(coxIn$drug)
rownames(coxIn) <- coxIn$drug
for (can in names(R)) {
  names(R[[can]]) <- toupper(names(R[[can]]))
}
R_ <- sapply(R, function(x) x[names(x) %in% coxIn$drug])

x <- sapply(names(R_), function(can) matrix(R_[[can]], dimnames = list(names(R_[[can]]), can)))
lowpIsSens <- highpIsSens <- list()
for (can in names(x)) {
  lowpIsSens[[can]] <- try(enrichment(x[[can]], coxIn, permute.n = 200, normalize = FALSE, pvalue.cutoff = 0.05)$pvalue)
  highpIsSens[[can]] <- try(enrichment(-x[[can]], coxIn, permute.n = 200, normalize = FALSE, pvalue.cutoff = 0.1)$pvalue)
}

groups <- unique(coxIn$group)
lowpIsSens <- sapply(lowpIsSens, function(x) as.data.frame(x)[groups, ])
highpIsSens <- sapply(highpIsSens, function(x) as.data.frame(x)[groups, ])
rownames(lowpIsSens) <- rownames(highpIsSens) <- groups
lowpIsSens <- lowpIsSens[, order(lowpIsSens["SIGNALING", ])]
lowpIsSens <- lowpIsSens[, order(lowpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["CYTOTOXIC", ])]
highpIsSens <- highpIsSens[, order(highpIsSens["SIGNALING", ])]

write.xlsx(t(lowpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr_repro.xlsx"), sheetName = "lowpIsSens")
write.xlsx(t(highpIsSens), file = file.path(out_dir, "drugsVsPloidyCorr_repro.xlsx"), sheetName = "highpIsSens", append = TRUE)

tmp <- sort(unique(coxIn$group))
col <- rainbow(length(tmp) * 1.3)[1:length(tmp)]
names(col) <- tmp
pdf(file.path(out_dir, "ploidyVsDrugSensitivity_repro.pdf"), width = 3, height = 6)
for (sheet in colnames(lowpIsSens)) {
  try(barplot(
    R_[[sheet]],
    col = col[coxIn[names(R_[[sheet]]), "group"]],
    main = sheet,
    horiz = TRUE,
    las = 2,
    cex.lab = 0.7,
    cex.names = 0.35,
    xlab = "Pearson r between ploidy and drug sensitivity (IC50)"
  ))
}
dev.off()
