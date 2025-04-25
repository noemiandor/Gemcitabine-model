# Load Seurat package
library(Seurat)

# 1. Define the parent directory containing all Cell Ranger outputs
parent_dir <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/data/SUM-159/A02_cellRanger"

# 2. List only the first‐level subdirectories (one per sample), returning full paths
sample_dirs <- list.dirs(
  path = parent_dir,
  full.names = TRUE,
  recursive = FALSE
)

# 3. Extract sample names from directory names, e.g. "2N-A1-0-Count-HM"
sample_names <- basename(sample_dirs)

# 4. Loop over each sample directory to:
#    a) point to the filtered feature-barcode matrix
#    b) read the 10X data
#    c) create a Seurat object
#    d) add sample name as metadata
seurat_list <- lapply(seq_along(sample_dirs), function(i) {
  # Path to the filtered matrix output by Cell Ranger
  data_dir <- file.path(sample_dirs[i], "outs", "filtered_feature_bc_matrix")
  
  # Read in counts matrix
  counts <- Read10X(data.dir = data_dir)
  
  # Create a Seurat object with basic filtering thresholds
  so <- CreateSeuratObject(
    counts = counts,
    project = sample_names[i],
    min.cells = 3,      # keep genes expressed in at least 3 cells
    min.features = 200  # keep cells with at least 200 detected features
  )
  
  # Store the sample identifier in metadata
  so$sample <- sample_names[i]
  
  return(so)
})

# Name each element of the list by its sample name
names(seurat_list) <- sample_names
saveRDS(seurat_list,    file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/seurat_list.Rds")

# 5. (Optional) Merge all individual Seurat objects into one combined object
#    Here we use Reduce() with merge(); cell IDs are prefixed by sample name
combined <- Reduce(
  f = function(x, y) merge(x, y, add.cell.ids = c(x$sample[1], y$sample[1])),
  x = seurat_list
)

saveRDS(combined,    file = "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/combined.Rds")


# Inspect the combined object
combined






