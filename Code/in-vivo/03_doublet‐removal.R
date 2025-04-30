# ───────────────────────────────────────────────────────────────────────────
# Pure R doublet‐removal pipeline (no LayerData dependencies)
# Methods: scDblFinder, computeDoubletDensity, QC outlier on nFeature_RNA
# ───────────────────────────────────────────────────────────────────────────

# Uncomment and run these lines if you haven’t installed the Bioconductor packages:
# if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
# BiocManager::install(c("scDblFinder", "scater"))

# 1. Load libraries ---------------------------------------------------------
library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)   # Method 1 & 2
library(scater)        # Method 3

# 2. Read data & convert to SCE --------------------------------------------
# Read your merged Seurat object
load('/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/data/SUM-159/C03_Integration/integrated_2025-04-17.RData')

integrated

######### all ########
# Convert to SingleCellExperiment for Bioconductor tools
sce <- as.SingleCellExperiment(integrated)

# 3. Estimate expected doublets --------------------------------------------
n_cells            <- ncol(sce)
doublet_rate       <- 0.075                    # adjust to your library
n_expected_doublet <- round(n_cells * doublet_rate)

# 4. Method 1: scDblFinder classification ----------------------------------
#   Simulates doublets and fits a mixture model
sce     <- scDblFinder(sce)
is_dbl1 <- colData(sce)$scDblFinder.class == "doublet"

# 5. Method 2: computeDoubletDensity ---------------------------------------
#   k-NN density scoring of doublet likelihood
density_scores <- computeDoubletDensity(sce)

#   Call top N cells as doublets
cutoff2 <- sort(density_scores, decreasing = TRUE)[n_expected_doublet]
is_dbl2 <- density_scores >= cutoff2

# 6. Method 3: QC outlier on nFeature_RNA -----------------------------------
# 6a. Compute per-cell QC metrics (adds 'detected' = number of genes)
sce <- addPerCellQC(sce)

# 6b. Flag cells with abnormally high gene counts as doublets
is_dbl3 <- isOutlier(
  sce$detected,   # pass the vector directly
  type  = "higher",
  log   = TRUE,
  nmads = 2
)

# 7. Intersection & filtering ----------------------------------------------
#   Identify barcodes called doublets by ALL three methods
doublet_barcodes <- colnames(sce)[ is_dbl1 & is_dbl2 & is_dbl3 ]

#   Subset the original Seurat object to remove those doublets
singlets <- subset(
  integrated,
  cells = setdiff(Cells(integrated), doublet_barcodes)
)

# 8. Save results -----------------------------------------------------------
saveRDS(
  singlets,
  file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds"
)

######### only cell line ########
cell_line<-c("2N-Cell-Culture","4N-Cell-Culture")

integrated_cell_line<-subset(integrated, subset= orig.ident %in% cell_line)

# Convert to SingleCellExperiment for Bioconductor tools
sce_cl <- as.SingleCellExperiment(integrated_cell_line)

# 3. Estimate expected doublets --------------------------------------------
n_cells            <- ncol(sce_cl)
doublet_rate       <- 0.075                    # adjust to your library
n_expected_doublet <- round(n_cells * doublet_rate)

# 4. Method 1: scDblFinder classification ----------------------------------
#   Simulates doublets and fits a mixture model
sce_cl     <- scDblFinder(sce_cl)
is_dbl1 <- colData(sce_cl)$scDblFinder.class == "doublet"

# 5. Method 2: computeDoubletDensity ---------------------------------------
#   k-NN density scoring of doublet likelihood
density_scores <- computeDoubletDensity(sce_cl)

#   Call top N cells as doublets
cutoff2 <- sort(density_scores, decreasing = TRUE)[n_expected_doublet]
is_dbl2 <- density_scores >= cutoff2

# 6. Method 3: QC outlier on nFeature_RNA -----------------------------------
# 6a. Compute per-cell QC metrics (adds 'detected' = number of genes)
sce_cl <- addPerCellQC(sce_cl)

# 6b. Flag cells with abnormally high gene counts as doublets
is_dbl3 <- isOutlier(
  sce_cl$detected,   # pass the vector directly
  type  = "higher",
  log   = TRUE,
  nmads = 2
)

# 7. Intersection & filtering ----------------------------------------------
#   Identify barcodes called doublets by ALL three methods
doublet_barcodes <- colnames(sce_cl)[ is_dbl1 & is_dbl2 & is_dbl3 ]

#   Subset the original Seurat object to remove those doublets
singlets_cl <- subset(
  integrated_cell_line,
  cells = setdiff(Cells(integrated_cell_line), doublet_barcodes)
)

# 8. Save results -----------------------------------------------------------
saveRDS(
  singlets_cl,
  file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds"
)




