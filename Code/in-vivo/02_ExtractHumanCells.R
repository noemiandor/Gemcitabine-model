#@@@@@@@@@@ QC and integration
# Load Seurat
library(Seurat)
library(dplyr)
library(Matrix)

seurat_list<-readRDS('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/00_ReadData/seurat_list.Rds')

# 1. Extract human cells from each sample
seurat_list <- lapply(seurat_list, function(obj) {
  # Identify human vs mouse genes by prefix
  all_genes   <- rownames(obj@assays$RNA)
  human_idx   <- grep("^GRCh38", all_genes)
  mouse_idx   <- grep("^GRCm39", all_genes)
  counts      <- obj@assays$RNA@counts
  # Sum UMIs per barcode for human and mouse genes
  human_counts <- Matrix::colSums(counts[human_idx, , drop=FALSE])
  mouse_counts <- Matrix::colSums(counts[mouse_idx, , drop=FALSE])
  # Keep barcodes with more human UMIs
  human_cells  <- names(which(human_counts > mouse_counts))
  # Skip samples with no human cells
  if (length(human_cells) == 0) {
    warning("No human cells found in sample: ", obj@project.name)
    return(NULL)
  }
  # Subset object to human cells only
  subset(obj, cells = human_cells)
})
# Remove any NULL entries (samples with no human cells)
seurat_list <- Filter(Negate(is.null), seurat_list)

# 2. Perform QC on human-only objects
qc_seurat_list <- lapply(seurat_list, function(obj) {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern="^MT-")
  subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 10)
})

# 2.1) Merge them into one Seurat object, preserving the sample ID in orig.ident
qc_combined <- merge(
  x = qc_seurat_list[[1]],
  y = qc_seurat_list[-1],
  add.cell.ids = names(qc_seurat_list),
  project      = "QC_Combined"
)

# 2.2) Now run a violin plot on the combined object
#    We’ll show three key QC metrics side by side, one violin per sample
pdf(
  file   = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/01_HumanCells/QC_after_filtering_violin.pdf",
  width  = 15,    # inches
  height = 6      # inches
)

print(
  VlnPlot(
    object   = qc_combined,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident",
    pt.size  = 0
  ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title  = element_text(hjust = 0.5)
    ) +
    ggtitle("QC Metrics After Filtering\n(one violin per sample)")
)

dev.off()

# 3. Preprocess each QC-filtered object for integration
preprocessed_list <- lapply(qc_seurat_list, function(obj) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method="vst", nfeatures=2000)
  obj <- ScaleData(obj, features=VariableFeatures(obj), verbose=FALSE)
  RunPCA(obj, features=VariableFeatures(obj), npcs=30, verbose=FALSE)
})

# 4. Select integration features and run CCA integration
features <- SelectIntegrationFeatures(preprocessed_list, nfeatures=3000)
anchors  <- FindIntegrationAnchors(preprocessed_list, anchor.features=features,
                                   reduction="cca", dims=1:30)
integrated <- IntegrateData(anchorset=anchors, dims=1:30)
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose=FALSE)
integrated <- RunPCA(integrated, npcs=30, verbose=FALSE)
integrated <- RunUMAP(integrated, reduction="pca", dims=1:30)
integrated <- FindNeighbors(integrated, reduction="pca", dims=1:30)
integrated <- FindClusters(integrated, resolution=0.5)




# 5. Save integrated, human-only object
saveRDS(integrated, file="/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/01_HumanCells/integrated_human_RC.Rds")