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


# Convert gene list to annotation dataframe for scGSVA
annot_df <- data.frame(
  gene = unlist(gs_list),
  geneSet = rep(names(gs_list), lengths(gs_list)),
  weight = 1,
  stringsAsFactors = FALSE
)

hsko<-readRDS('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/hsko.Rds')

# Prepare Hallmark gene sets
m_df <- msigdbr(species = "Homo sapiens", category = "H")
# Prepare TERM2GENE as a two-column data.frame for GSEA
term2gene_df <- m_df[, c("gene_symbol","gs_id","gs_name")]
term2gene_df<-as.data.frame(term2gene_df)
colnames(term2gene_df)<-c("GeneID","PATH","Annot")



hallmark_sets<-hsko

hallmark_sets@annot<-term2gene_df
hallmark_sets@anntype<-"HALLMARK"

saveRDS(hallmark_sets,file = '/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/hallmark_sets.Rds')


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





gsva_res<-readRDS('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_res_brest_cancer.Rds')


Dose<-read.table('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/IDs_Dose.txt',header = T)


meta <- singlets@meta.data %>%
  rownames_to_column("cell") %>%                  # preserve barcode
  left_join(Dose, by = c("orig.ident" = "IDs")) %>%
  column_to_rownames("cell")

singlets@meta.data <- meta

# 2. Ensure metadata columns are factors with the desired order
singlets$sample_type       <- factor(singlets$sample_type,
                                     levels = c("2N-cellline","4N-cellline","2N-tumor","4N-tumor"))
singlets$Dose              <- factor(singlets$Dose,
                                     levels = c("0mg/kg","30mg/kg","120mg/kg"))
singlets$seurat_clusters   <- factor(singlets$seurat_clusters)

gsva_res@obj@meta.data<-singlets@meta.data

gsva_res_use<-gsva_res

orig_meta <- singlets@meta.data
empty_mat <- matrix(0, nrow = 0, ncol = nrow(orig_meta))
colnames(empty_mat) <- rownames(orig_meta) 
rownames(empty_mat) <- character(0) 

seurat_meta_only <- CreateSeuratObject(
  counts = empty_mat,
  meta.data = orig_meta,
  project = "MetaOnly"
)

gsva_res_use@obj<-seurat_meta_only

saveRDS(gsva_res_use,file = '/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_res_brest_cancer_use.Rds')

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





####################################### hallmark






gsva_res_HM<-readRDS('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_res_brest_cancer_hallmarks.Rds')


gsva_res_HM@obj@meta.data<-singlets@meta.data

gsva_res_HM_use<-gsva_res_HM

gsva_res_HM_use@obj<-seurat_meta_only

saveRDS(gsva_res_HM_use,file = '/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/gsva_res_brest_cancer_hallmarks_use.Rds')






