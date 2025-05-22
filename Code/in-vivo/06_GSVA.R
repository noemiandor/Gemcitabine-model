library(Seurat)

singlets <-readRDS("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds")

# Load required packages for GSVA
library(msigdbr)
library(GSVA)
# Load future for parallelization
library(future)
# Increase maximum allowed size for future globals (e.g., for NormalizeData)
options(future.globals.maxSize = 8000 * 1024^2)  # ~8 GB
# Choose appropriate parallel processing plan: multicore if supported, otherwise multisession
if (future::supportsMulticore()) {
  plan("multicore", workers = 4)
} else {
  plan("multisession", workers = 4)
}

singlets <- NormalizeData(
  object = singlets,
  normalization.method = "LogNormalize",  
  scale.factor = 10000                    
)

# Prepare normalized expression matrix for GSVA
expr_mat <- as.matrix(GetAssayData(singlets, slot = "data"))

# Retrieve Hallmark gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
gs_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Run GSVA
gsva_res <- gsva(
  expr = expr_mat,
  gset.idx.list = gs_list,
  method = "gsva",
  kcdf = "Gaussian",
  verbose = TRUE
)

# Add GSVA scores back into Seurat metadata
gsva_df <- as.data.frame(t(gsva_res))
singlets <- AddMetaData(singlets, metadata = gsva_df)

# Save GSVA results and updated Seurat object
saveRDS(gsva_res, file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_results.Rds")
saveRDS(singlets, file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/singlets_with_gsva.Rds")

# Define features of interest: key Hallmark pathways and MKI67 gene
pathways_to_plot <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_APOPTOSIS"
)

# UMAP visualization of GSVA scores for selected pathways
FeaturePlot(
  object = singlets,
  features = pathways_to_plot,
  reduction = "umap",
  pt.size = 0.5
)

# UMAP visualization of MKI67 gene expression
FeaturePlot(
  object = singlets,
  features = "MKI67",
  reduction = "umap",
  pt.size = 0.5
)

# Violin plots for GSVA pathway activity across clusters
VlnPlot(
  object = singlets,
  features = pathways_to_plot,
  group.by = "seurat_clusters",
  pt.size = 0
)

# Violin plot for MKI67 expression across clusters
VlnPlot(
  object = singlets,
  features = "MKI67",
  group.by = "seurat_clusters",
  pt.size = 0
)

# Save UMAP FeaturePlots to PDF
pdf("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/GSVA_UMAP_FeaturePlots.pdf", width = 10, height = 8)
FeaturePlot(
  object = singlets,
  features = pathways_to_plot,
  reduction = "umap",
  pt.size = 0.5
)
FeaturePlot(
  object = singlets,
  features = "MKI67",
  reduction = "umap",
  pt.size = 0.5
)
dev.off()

# Save Violin plots to PDF
pdf("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/GSVA_ViolinPlots.pdf", width = 10, height = 8)
VlnPlot(
  object = singlets,
  features = pathways_to_plot,
  group.by = "seurat_clusters",
  pt.size = 0
)
VlnPlot(
  object = singlets,
  features = "MKI67",
  group.by = "seurat_clusters",
  pt.size = 0
)
dev.off()
