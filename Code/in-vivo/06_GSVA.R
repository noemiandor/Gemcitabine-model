library(Seurat)

singlets <-readRDS("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds")

# Load required packages for GSVA
library(msigdbr)
# library(GSVA)

singlets <- NormalizeData(
  object = singlets,
  normalization.method = "LogNormalize",  
  scale.factor = 10000                    
)

# Prepare normalized expression matrix for GSVA
#expr_mat <- as.matrix(GetAssayData(singlets, slot = "data"))

# Retrieve Hallmark gene sets from MSigDB
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H")
gs_list <- split(hallmark_sets$gene_symbol, hallmark_sets$gs_name)

# Convert gene list to annotation dataframe for scGSVA
annot_df <- data.frame(
  gene = unlist(gs_list),
  geneSet = rep(names(gs_list), lengths(gs_list)),
  weight = 1,
  stringsAsFactors = FALSE
)

hsko<-readRDS('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/hsko.Rds')



# Load scGSVA and BiocParallel for parallel execution
library(scGSVA)

# Run scGSVA on the Seurat object using custom Hallmark gene sets
gsva_res <- scgsva(
  obj = singlets,
  annot = hsko,
  method = "gsva",        # or "ssgsea" if preferred
  kcdf = "Gaussian",
  cores = 16,
  verbose = TRUE
)





res_pathway<-findPathway(gsva_res,group = "Subpopulation")
res_pathway_sig<-sigPathway(gsva_res,group = "Subpopulation")

res_pathway_sig<-res_pathway_sig[order(res_pathway_sig$p.adj),]
#res_pathway_sig<-res_pathway_sig[which(res_pathway_sig$adj.P.Val < 0.05),]
sig_pathway<-unique(res_pathway_sig$Path)


res_pathway_sig_annova<-sigPathway(gsva_res,group = "Subpopulation",test.use = "anova")
res_pathway_sig_annova<-res_pathway_sig_annova[order(res_pathway_sig_annova$p.adj),]

write.table(res_pathway_sig_annova,'/Volumes/Work/projects/BJ/active/20220330-GGT211227006/runtime_ver2/output/04_subpopulations/res_pathway_sig_annova.txt',row.names=F,col.names=T,quote=F,sep="\t")






# Add GSVA scores back into Seurat metadata
gsva_df <- as.data.frame(gsva_res@gsva)
colnames(gsva_df)<-"GSVA"
singlets <- AddMetaData(singlets, metadata = gsva_df)

# Save GSVA results and updated Seurat object
saveRDS(gsva_res, file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_results.Rds")
saveRDS(singlets, file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/singlets_with_gsva.Rds")

