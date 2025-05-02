# ───────────────────────────────────────────────────────────────────────────
# Differential Expression Analyses & Plots for singlets
# 1) sample_type comparisons
# 2) dose comparisons (within 2N-tumor & 4N-tumor)
# 3) cluster × sample_type comparisons (2N-tumor vs 4N-tumor in each cluster)
# ───────────────────────────────────────────────────────────────────────────

library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)

setwd('/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/03_DiffExp_genes')
# 2. Read in Seurat object
singlets<-readRDS('/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds')
# 3. Read in dose info
Dose<-read.table('/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/data/SUM-159/IDs_Dose.txt',header = T)

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

# output dirs
out_dir_de   <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/03_DiffExp_genes/DE_Results"
out_dir_pdf  <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/03_DiffExp_genes/DE_Plots"
dir.create(out_dir_de,   recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_pdf,  recursive = TRUE, showWarnings = FALSE)
# ───────────────────────────────────────────────────────────────────────────
# 1) sample_type comparisons
# ───────────────────────────────────────────────────────────────────────────
Idents(singlets) <- singlets$sample_type

comparisons1 <- list(
  `2N-cellline_vs_4N-cellline` = c("2N-cellline","4N-cellline"),
  `2N-tumor_vs_4N-tumor`       = c("2N-tumor",    "4N-tumor")
)

for (nm in names(comparisons1)) {
  grp <- comparisons1[[nm]]
  markers <- FindMarkers(
    object       = singlets,
    ident.1      = grp[1],
    ident.2      = grp[2],
    min.pct      = 0.1,
    logfc.threshold = 0.25,
    test.use     = "wilcox"
  )
  # save table
  write.csv(
    markers,
    file = file.path(out_dir_de, paste0("DE_sampletype_", nm, ".csv"))
  )
  
  # Volcano
  p1<-EnhancedVolcano(
    markers,
    lab        = rownames(markers),
    x          = "avg_log2FC",
    y          = "p_val_adj",
    pCutoff    = 0.05,
    FCcutoff   = 0.5,
    title      = nm
  )
  pdf(file.path(out_dir_pdf, paste0("Volcano_sampletype_", nm, ".pdf")), 
      width=6, height=6)
  on.exit(dev.off(), add = TRUE)
  print(p1)
  dev.off()
  
  # Violin for top 5 up/down
  top5_up   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
  top5_down <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
  
  p2<-VlnPlot(
    singlets,
    features = c(top5_up, top5_down),
    group.by = "sample_type",
    pt.size  = 0.1
  ) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Vln_sampletype_", nm, ".pdf")), 
      width=6, height=6)
  on.exit(dev.off(), add = TRUE)
  print(p2)
  dev.off()
  
  # Heatmap for top 10 (by absolute logFC)
  top10 <- rownames(markers %>% 
                      mutate(absFC=abs(avg_log2FC)) %>% 
                      arrange(desc(absFC)) %>% 
                      head(10))
  p3<-DoHeatmap(
    singlets,
    features = top10,
    group.by = "sample_type",
    assay    = "RNA",
    slot     = "data"
  ) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Heatmap_sampletype_", nm, ".pdf")), 
      width=6, height=6)
  on.exit(dev.off(), add = TRUE)
  print(p3)
  dev.off()
}

# ───────────────────────────────────────────────────────────────────────────
# 2) dose comparisons within 2N-tumor and 4N-tumor
# ───────────────────────────────────────────────────────────────────────────
# subset to tumors
tumor <- subset(singlets, subset = sample_type %in% c("2N-tumor","4N-tumor"))

# for each tumor type separately:
for (stype in c("2N-tumor","4N-tumor")) {
  so <- subset(tumor, subset = sample_type==stype)
  Idents(so) <- so$Dose
  doses <- unique(so$Dose)
  combos <- combn(doses, 2, simplify = FALSE)
  for (cmb in combos) {
    nm <- paste0(gsub("-","",stype),"_", gsub("mg/kg","",cmb[1]),
                 "_vs_", gsub("mg/kg","",cmb[2]))
    markers <- FindMarkers(
      object          = so,
      ident.1         = cmb[1],
      ident.2         = cmb[2],
      min.pct         = 0.1,
      logfc.threshold = 0.25,
      test.use        = "wilcox"
    )
    write.csv(
      markers,
      file = file.path(out_dir_de, paste0("DE_dose_", nm, ".csv"))
    )
    # Volcano
    
    p1<-EnhancedVolcano(
      markers,
      lab      = rownames(markers),
      x        = "avg_log2FC",
      y        = "p_val_adj",
      pCutoff  = 0.05,
      FCcutoff = 0.5,
      title    = nm
    )
    pdf(file.path(out_dir_pdf, paste0("Volcano_dose_", nm, ".pdf")),6,6)
    on.exit(dev.off(), add = TRUE)
    print(p1)
    dev.off()
    # Violin top genes
    genes_up   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
    genes_down <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
    p2<-VlnPlot(so, features=c(genes_up, genes_down), group.by="Dose", pt.size=0.1) + NoLegend()
    pdf(file.path(out_dir_pdf, paste0("Vln_dose_", nm, ".pdf")),6,6)
    on.exit(dev.off(), add = TRUE)
    print(p2)
    dev.off()
    # Heatmap top10
    top10 <- rownames(markers %>% mutate(absFC=abs(avg_log2FC)) %>% arrange(desc(absFC)) %>% head(10))
    p3<-DoHeatmap(so, features=top10, group.by="Dose",assay    = "RNA",
              slot     = "data") + NoLegend()
    pdf(file.path(out_dir_pdf, paste0("Heatmap_dose_", nm, ".pdf")),6,6)
    on.exit(dev.off(), add = TRUE)
    print(p3)
    dev.off()
  }
}

# ───────────────────────────────────────────────────────────────────────────
# 3) cluster × sample_type comparisons (2N-tumor vs 4N-tumor in each cluster)
# ───────────────────────────────────────────────────────────────────────────
# ensure tumor still defined
tumor <- subset(singlets, subset = sample_type %in% c("2N-tumor","4N-tumor"))
clusters <- levels(tumor$seurat_clusters)

for (cl in clusters) {
  so <- subset(tumor, subset=seurat_clusters==cl)
  n2N <- sum(so$sample_type == "2N-tumor")
  n4N <- sum(so$sample_type == "4N-tumor")
  
  #  Skip if either group has <3 cells
  if (n2N < 3 || n4N < 3) {
    message("Skipping cluster ", cl, 
            ": only ", n2N, " 2N-tumor and ", 
            n4N, " 4N-tumor cells.")
    next
  }
  Idents(so) <- so$sample_type
  markers <- FindMarkers(
    object          = so,
    ident.1         = "2N-tumor",
    ident.2         = "4N-tumor",
    min.pct         = 0.1,
    logfc.threshold = 0.25,
    test.use        = "wilcox"
  )
  nm <- paste0("cluster", cl)
  write.csv(
    markers,
    file = file.path(out_dir_de, paste0("DE_cluster_", nm, ".csv"))
  )
  # Volcano
  p1<-EnhancedVolcano(
    markers,
    lab      = rownames(markers),
    x        = "avg_log2FC",
    y        = "p_val_adj",
    pCutoff  = 0.05,
    FCcutoff = 0.5,
    title    = nm
  )
  pdf(file.path(out_dir_pdf, paste0("Volcano_cluster_", nm, ".pdf")),6,6)
  on.exit(dev.off(), add = TRUE)
  print(p1)
  dev.off()
  # Violin top genes
  up5   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
  down5 <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
  p2<-VlnPlot(so, features=c(up5, down5), group.by="sample_type", pt.size=0.1) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Vln_cluster_", nm, ".pdf")),6,6)
  on.exit(dev.off(), add = TRUE)
  print(p2)
  dev.off()
  # Heatmap top10
  top10 <- rownames(markers %>% mutate(absFC=abs(avg_log2FC)) %>% arrange(desc(absFC)) %>% head(10))
  p3<-DoHeatmap(so, features = top10, group.by="sample_type",assay    = "RNA",
            slot     = "data") + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Heatmap_cluster_", nm, ".pdf")),6,6)
  on.exit(dev.off(), add = TRUE)
  print(p3)
  dev.off()
}















