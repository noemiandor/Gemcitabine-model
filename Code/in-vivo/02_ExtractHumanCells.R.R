# 02_ExtractHumanCells.R

# Load required packages
library(Seurat)
library(DropletUtils)

# 1. Read in the combined Seurat object
combined <- readRDS(
  "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/combined.Rds"
)

# 2. Define where your CellRanger outputs live
parent_dir <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/data/SUM-159/A02_cellRanger"

# 3. List each sample folder (full paths)
sample_dirs <- list.dirs(
  path      = parent_dir,
  full.names = TRUE,
  recursive  = FALSE
)

# 4. For each sample, read its molecule_info.h5 and pull barcodes assigned to GRCh38
human_barcodes_list <- lapply(sample_dirs, function(sample_path) {
  # Path to the H5 file
  h5_file <- file.path(sample_path, "outs", "molecule_info.h5")
  
  # Read per-molecule info; returns a DataFrame with at least 'barcode' and 'genome' columns
  mol_info <- read10xMolInfo(h5_file)
  
  # Keep only barcodes where genome == "GRCh38"
  human <- mol_info$barcode[mol_info$genome == "GRCh38"]
  unique(human)
})

# Combine into one vector of barcodes
human_barcodes <- unique(unlist(human_barcodes_list))

# 5. Subset the merged Seurat object to human cells only
human_cells <- subset(combined, cells = human_barcodes)

# 6. Save the human-only Seurat object
saveRDS(
  human_cells,
  file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/human_cells.Rds"
)