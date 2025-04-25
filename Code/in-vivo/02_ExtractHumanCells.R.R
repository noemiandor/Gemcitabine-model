# 02_ExtractHumanCells.R

# Load required packages
library(Seurat)
library(DropletUtils)

# 1. Read in the combined Seurat object
combined <- readRDS(
  "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/combined.Rds"
)

# 2. Get all gene names and split into human vs mouse
all_genes   <- rownames(combined@assays$RNA)
human_genes <- grep("^GRCh38", all_genes, value = TRUE)
mouse_genes <- grep("^GRCm39", all_genes, value = TRUE)

# 3. Pull the raw (or normalized) count matrix
#    (counts is usually a sparse Matrix, so we use colSums)
counts      <- combined@assays$RNA@counts

# 4. Compute per-cell total for each species
human_counts <- Matrix::colSums(counts[human_genes, , drop = FALSE])
mouse_counts <- Matrix::colSums(counts[mouse_genes, , drop = FALSE])

# 5. Define “human” cells as those with more human UMIs than mouse UMIs
human_cells  <- names(which(human_counts > mouse_counts))

# 6. Subset the Seurat object
human <- subset(combined, cells = human_cells)

# 7. Save the human-only Seurat object
saveRDS(
  human,
  file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/01_HumanCells/human_cells.Rds"
)
