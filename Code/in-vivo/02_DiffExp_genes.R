#!/usr/bin/env Rscript

script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_match <- grep(file_arg, cmd_args, value = TRUE)
if (length(file_match) > 0) {
  script_path <- normalizePath(sub(file_arg, "", file_match[1]), mustWork = FALSE)
}
if (is.null(script_path) || !nzchar(script_path)) {
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = FALSE) else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    script_path <- frame_files[length(frame_files)]
  }
}
script_dir <- if (!is.null(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

# ───────────────────────────────────────────────────────────────────────────
# Differential Expression Analyses & Plots for obj
# 1) sample_type comparisons
# 2) dose comparisons (within 2N-tumor & 4N-tumor)
# 3) cluster × sample_type comparisons (2N-tumor vs 4N-tumor in each cluster)
# ───────────────────────────────────────────────────────────────────────────

library(Seurat)
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(tibble)

output_root <- file.path(results_root, "02_DiffExp_genes")
out_dir_de <- file.path(output_root, "DE_Results")
out_dir_pdf <- file.path(output_root, "DE_Plots")
.ensure_dir(output_root)
.ensure_dir(out_dir_de)
.ensure_dir(out_dir_pdf)

# 2. Read in Seurat object
obj <- readRDS(file.path(results_root, "01_data", "integrated_sct_cca_seurat.rds"))

obj$sample_type <- dplyr::case_when(
  obj$sample == "2N-Cell-Culture" ~ "2N-cellline",
  obj$sample == "4N-Cell-Culture" ~ "4N-cellline",
  grepl("^2N", obj$sample) ~ "2N-tumor",
  grepl("^4N", obj$sample) ~ "4N-tumor",
  grepl("^A5", obj$sample) ~ "4N-tumor",
  grepl("^A6", obj$sample) ~ "4N-tumor",
  TRUE ~ NA_character_
)

na_samples <- sort(unique(obj$sample[is.na(obj$sample_type)]))
if (length(na_samples) > 0) {
  warning("sample_type is NA for samples: ", paste(na_samples, collapse = ", "))
}

# 2. Ensure metadata columns are factors with the desired order
obj$sample_type       <- factor(obj$sample_type,
                                     levels = c("2N-cellline","4N-cellline","2N-tumor","4N-tumor"))

obj$Dose <- dplyr::case_when(
  obj$Dose == "0" ~ "0mg/kg",
  obj$Dose == "30" ~ "30mg/kg",
  obj$Dose == "120" ~ "120mg/kg",
  TRUE ~ NA_character_
)

obj$Dose              <- factor(obj$Dose,
                                     levels = c("0mg/kg","30mg/kg","120mg/kg"))

obj$seurat_clusters   <- factor(obj$seurat_clusters)

# ───────────────────────────────────────────────────────────────────────────
# 1) sample_type comparisons
# ───────────────────────────────────────────────────────────────────────────
Idents(obj) <- obj$sample_type

comparisons1 <- list(
  `2N-cellline_vs_4N-cellline` = c("2N-cellline","4N-cellline"),
  `2N-tumor_vs_4N-tumor`       = c("2N-tumor",    "4N-tumor"),
  `2N-cellline_vs_4N-tumor`       = c("2N-cellline",    "4N-tumor"),
  `2N-tumor_vs_4N-cellline`       = c("2N-tumor",    "4N-cellline"),
  `2N-cellline_vs_2N-tumorr`       = c("2N-cellline",    "2N-tumor"),
  `4N-tumor_vs_4N-cellline`       = c("4N-tumor",    "4N-cellline")
)

DefaultAssay(obj) <- "RNA"

if (nrow(obj[["RNA"]]@data) == 0) {
  obj <- NormalizeData(obj, verbose = FALSE)
}

obj <- ScaleData(
  object = obj,
  assay = "RNA",
  features = rownames(obj),
  verbose = FALSE
)

strip_prefix_all_slots <- function(obj, prefix_regex = "^GRCh38[-_]") {
  assays <- Assays(obj)
  
  old_feats <- unique(unlist(lapply(assays, function(a) rownames(obj[[a]]))))
  new_feats <- make.unique(sub(prefix_regex, "", old_feats, perl = TRUE))
  map <- setNames(new_feats, old_feats)
  
  remap <- function(x) {
    y <- unname(map[as.character(x)])
    y[is.na(y)] <- x[is.na(y)]
    y
  }
  
  for (a in assays) {
    assay <- obj[[a]]
    
    if (nrow(assay@counts) > 0) rownames(assay@counts) <- remap(rownames(assay@counts))
    if (nrow(assay@data) > 0) rownames(assay@data) <- remap(rownames(assay@data))
    if (nrow(assay@scale.data) > 0) rownames(assay@scale.data) <- remap(rownames(assay@scale.data))
    if (nrow(assay@meta.features) > 0) rownames(assay@meta.features) <- remap(rownames(assay@meta.features))
    
    vf <- VariableFeatures(assay)
    if (length(vf) > 0) VariableFeatures(assay) <- remap(vf)
    
    if (inherits(assay, "SCTAssay") && length(assay@SCTModel.list) > 0) {
      for (m in names(assay@SCTModel.list)) {
        model <- assay@SCTModel.list[[m]]
        if ("feature.attributes" %in% slotNames(model) && nrow(model@feature.attributes) > 0) {
          rownames(model@feature.attributes) <- remap(rownames(model@feature.attributes))
        }
        assay@SCTModel.list[[m]] <- model
      }
    }
    
    obj[[a]] <- assay
  }
  
  for (dr in Reductions(obj)) {
    red <- obj[[dr]]
    if (nrow(red@feature.loadings) > 0) {
      rownames(red@feature.loadings) <- remap(rownames(red@feature.loadings))
    }
    if (nrow(red@feature.loadings.projected) > 0) {
      rownames(red@feature.loadings.projected) <- remap(rownames(red@feature.loadings.projected))
    }
    obj[[dr]] <- red
  }
  
  obj
}

obj <- strip_prefix_all_slots(obj, "^GRCh38[-_]")


for (nm in names(comparisons1)) {
  grp <- comparisons1[[nm]]
  # Subset to only the two sample types being compared
  so <- subset(obj, subset = sample_type %in% grp)
  markers <- FindMarkers(
    object       = so,
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
      width=12, height=12)
  on.exit(dev.off(), add = TRUE)
  print(p1)
  dev.off()
  
  # Violin for top 5 up/down
  top5_up   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
  top5_down <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
  
  p2<-VlnPlot(
    so,
    features = c(top5_up, top5_down),
    group.by = "sample_type",
    pt.size  = 0.1
  ) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Vln_sampletype_", nm, ".pdf")),
      width=12, height=12)
  on.exit(dev.off(), add = TRUE)
  print(p2)
  dev.off()
  
  # Heatmap for top 30 (by absolute logFC)
  top30 <- rownames(markers %>%
                      mutate(absFC=abs(avg_log2FC)) %>%
                      arrange(desc(absFC)) %>%
                      head(30))
  p3<-DoHeatmap(
    so,
    features = top30,
    group.by = "sample_type",
    assay    = "RNA",
    slot     = "scale.data",
    raster = FALSE
  ) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Heatmap_sampletype_", nm, ".pdf")),
      width=6, height=6)
  on.exit(dev.off(), add = TRUE)
  print(p3)
  dev.off()
  rm(so)
  gc()
}

# ───────────────────────────────────────────────────────────────────────────
# 2) dose comparisons within 2N-tumor and 4N-tumor
# ───────────────────────────────────────────────────────────────────────────
# subset to tumors
tumor <- subset(obj, subset = sample_type %in% c("2N-tumor","4N-tumor"))

# for each tumor type separately:
for (stype in c("2N-tumor","4N-tumor")) {
  so <- subset(tumor, subset = sample_type==stype)
  # Drop unused levels of sample_type so that only tumor levels remain
  so$sample_type <- droplevels(so$sample_type)
  Idents(so) <- so$Dose
  dose_order <- c("0mg/kg", "30mg/kg", "120mg/kg")
  doses <- dose_order[dose_order %in% unique(as.character(so$Dose))]
  if (length(doses) < 2) {
    message("Skipping ", stype, ": fewer than 2 dose levels present.")
    next
  }
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
    pdf(file.path(out_dir_pdf, paste0("Volcano_dose_", nm, ".pdf")),12,12)
    on.exit(dev.off(), add = TRUE)
    print(p1)
    dev.off()
    # Violin top genes
    genes_up   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
    genes_down <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
    p2<-VlnPlot(so, features=c(genes_up, genes_down), group.by="Dose", pt.size=0.1) + NoLegend()
    pdf(file.path(out_dir_pdf, paste0("Vln_dose_", nm, ".pdf")),12,12)
    on.exit(dev.off(), add = TRUE)
    print(p2)
    dev.off()
    # Heatmap top30
    top30 <- rownames(markers %>% mutate(absFC=abs(avg_log2FC)) %>% arrange(desc(absFC)) %>% head(30))
    p3<-DoHeatmap(so, features=top30, group.by="Dose",assay    = "RNA",
              slot     = "scale.data",raster = FALSE) + NoLegend()
    pdf(file.path(out_dir_pdf, paste0("Heatmap_dose_", nm, ".pdf")),6,6)
    on.exit(dev.off(), add = TRUE)
    print(p3)
    dev.off()
  }
}

# ───────────────────────────────────────────────────────────────────────────
# 2.5) tumor-type comparisons at the same Dose
# ───────────────────────────────────────────────────────────────────────────
for (d in levels(obj$Dose)) {
  so_d <- subset(tumor, subset = Dose == d)
  # Drop unused levels of sample_type so that only tumor levels remain
  so_d$sample_type <- droplevels(so_d$sample_type)
  Idents(so_d) <- so_d$sample_type
  markers_dt <- FindMarkers(
    object          = so_d,
    ident.1         = "2N-tumor",
    ident.2         = "4N-tumor",
    min.pct         = 0.1,
    logfc.threshold = 0.25,
    test.use        = "wilcox"
  )
  nm_dt <- paste0("tumorType_", gsub("mg/kg","",d))
  write.csv(
    markers_dt,
    file = file.path(out_dir_de, paste0("DE_", nm_dt, ".csv"))
  )
  # Volcano plot
  p_dt <- EnhancedVolcano(
    markers_dt,
    lab      = rownames(markers_dt),
    x        = "avg_log2FC",
    y        = "p_val_adj",
    pCutoff  = 0.05,
    FCcutoff = 0.5,
    title    = nm_dt
  )
  ggsave(
    filename = file.path(out_dir_pdf, paste0("Volcano_", nm_dt, ".pdf")),
    plot     = p_dt,
    width    = 12,
    height   = 12,
    units    = "in"
  )
  # Violin plot for top 5 up/down genes
  top5_up_dt   <- rownames(markers_dt %>% filter(avg_log2FC > 0) %>% head(5))
  top5_down_dt <- rownames(markers_dt %>% filter(avg_log2FC < 0) %>% head(5))
  p_vln_dt <- VlnPlot(
    so_d,
    features = c(top5_up_dt, top5_down_dt),
    group.by = "sample_type",
    pt.size  = 0.1
  ) + NoLegend()
  ggsave(
    filename = file.path(out_dir_pdf, paste0("Vln_", nm_dt, ".pdf")),
    plot     = p_vln_dt,
    width    = 12,
    height   = 12,
    units    = "in"
  )

  # Heatmap for top 30 genes by absolute logFC
  top30_dt <- rownames(markers_dt %>%
                        mutate(absFC = abs(avg_log2FC)) %>%
                        arrange(desc(absFC)) %>%
                        head(30))
  pdf(
    file   = file.path(out_dir_pdf, paste0("Heatmap_", nm_dt, ".pdf")),
    width  = 6,
    height = 6
  )
  p_hm_dt <- DoHeatmap(
    so_d,
    features = top30_dt,
    group.by = "sample_type",
    assay    = "RNA",
    slot     = "scale.data",
    raster = FALSE
  ) + NoLegend()
  print(p_hm_dt)
  dev.off()
}

# ───────────────────────────────────────────────────────────────────────────
# 3) cluster × sample_type comparisons (2N-tumor vs 4N-tumor in each cluster)
# ───────────────────────────────────────────────────────────────────────────
# ensure tumor still defined
tumor <- subset(obj, subset = sample_type %in% c("2N-tumor","4N-tumor"))
clusters <- levels(tumor$seurat_clusters)

for (cl in clusters) {
  so <- subset(tumor, subset=seurat_clusters==cl)
  # Drop unused levels of sample_type so that only tumor levels remain
  so$sample_type <- droplevels(so$sample_type)
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
  pdf(file.path(out_dir_pdf, paste0("Volcano_cluster_", nm, ".pdf")),12,12)
  on.exit(dev.off(), add = TRUE)
  print(p1)
  dev.off()
  # Violin top genes
  up5   <- rownames(markers %>% filter(avg_log2FC>0) %>% head(5))
  down5 <- rownames(markers %>% filter(avg_log2FC<0) %>% head(5))
  p2<-VlnPlot(so, features=c(up5, down5), group.by="sample_type", pt.size=0.1) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Vln_cluster_", nm, ".pdf")),12,12)
  on.exit(dev.off(), add = TRUE)
  print(p2)
  dev.off()
  # Heatmap top30
  top30 <- rownames(markers %>% mutate(absFC=abs(avg_log2FC)) %>% arrange(desc(absFC)) %>% head(30))
  p3<-DoHeatmap(so, features = top30, group.by="sample_type",assay    = "RNA",
            slot     = "scale.data",raster = FALSE) + NoLegend()
  pdf(file.path(out_dir_pdf, paste0("Heatmap_cluster_", nm, ".pdf")),6,6)
  on.exit(dev.off(), add = TRUE)
  print(p3)
  dev.off()
}


table(obj@meta.data$seurat_clusters,obj@meta.data$sample_type)

pdf(file.path(out_dir_pdf, "UMAP.pdf"),7,6)
on.exit(dev.off(), add = TRUE)
DimPlot(
  object    = obj,
  reduction = "umap",
  group.by  = "seurat_clusters",
  pt.size   = 0.5
) + ggtitle("UMAP colored by Cluster")
dev.off()






tumor <- subset(obj, subset = sample_type %in% c("2N-tumor", "4N-tumor"))
tumor$sample_type <- droplevels(factor(tumor$sample_type, levels = c("2N-tumor", "4N-tumor")))


out_dir <- out_dir_pdf

p_cluster <- DimPlot(
  tumor,
  reduction = "umap",
  group.by = "seurat_clusters",
  pt.size = 0.4
) + ggtitle("Tumor only UMAP - by cluster")

ggsave(
  filename = file.path(out_dir, "UMAP_tumor_by_cluster.pdf"),
  plot = p_cluster,
  width = 7, height = 6
)


p_type <- DimPlot(
  tumor,
  reduction = "umap",
  group.by = "sample_type",
  pt.size = 0.4
) + ggtitle("Tumor only UMAP - 2N-tumor vs 4N-tumor")

ggsave(
  filename = file.path(out_dir, "UMAP_tumor_by_sample_type.pdf"),
  plot = p_type,
  width = 7, height = 6
)



DefaultAssay(obj) <- "RNA"
if (nrow(obj[["RNA"]]@data) == 0) {
  obj <- NormalizeData(obj, verbose = FALSE)
}



dose_chr <- as.character(obj$Dose)
dose_chr[dose_chr == "0"] <- "0mg/kg"
dose_chr[dose_chr == "30"] <- "30mg/kg"
dose_chr[dose_chr == "120"] <- "120mg/kg"
obj$Dose <- dose_chr


obj_sub <- subset(
  obj,
  subset = sample_type %in% c("2N-cellline", "4N-cellline") |
    (sample_type %in% c("2N-tumor", "4N-tumor") & Dose == "0mg/kg")
)


obj_sub$ploidy <- ifelse(grepl("^2N", obj_sub$sample_type), "2N", "4N")
obj_sub$ploidy <- factor(obj_sub$ploidy, levels = c("2N", "4N"))
Idents(obj_sub) <- "ploidy"


de_2N_vs_4N <- FindMarkers(
  object = obj_sub,
  ident.1 = "2N",
  ident.2 = "4N",
  min.pct = 0.1,
  logfc.threshold = 0.25,
  test.use = "wilcox",
  assay = "RNA",
  slot = "data"
)


out_dir <- out_dir_pdf


p_ploidy <- DimPlot(
  obj_sub,
  reduction = "umap",
  group.by = "ploidy",
  pt.size = 0.4
) + ggtitle("CellLine + Dose0 tumor subset: 2N vs 4N")

ggsave(
  filename = file.path(out_dir, "UMAP_subset_CellLinePlusDose0Tumor_by_2N4N.pdf"),
  plot = p_ploidy,
  width = 7, height = 6
)

# UMAP 2: 按原sample_type着色（查看组成）
p_type <- DimPlot(
  obj_sub,
  reduction = "umap",
  group.by = "sample_type",
  pt.size = 0.4
) + ggtitle("CellLine + Dose0 tumor subset: sample_type")

ggsave(
  filename = file.path(out_dir, "UMAP_subset_CellLinePlusDose0Tumor_by_sample_type.pdf"),
  plot = p_type,
  width = 8, height = 6
)

# UMAP 3: 按cluster着色
p_cluster <- DimPlot(
  obj_sub,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.4
) + ggtitle("CellLine + Dose0 tumor subset: clusters")

ggsave(
  filename = file.path(out_dir, "UMAP_subset_CellLinePlusDose0Tumor_by_cluster.pdf"),
  plot = p_cluster,
  width = 8, height = 6
)

# 可选保存DE结果
# write.csv(
#   de_2N_vs_4N,
#   file.path(out_dir_de, "DE_2N_vs_4N_in_CellLinePlusDose0Tumor.csv")
# )



















saveRDS(obj, file.path(output_root, "seu_obj.Rds"))
