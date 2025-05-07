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
library(scater)
library(ggplot2)
library(VennDiagram)
library(tibble)
library(grid)
library(scds)          # For SCDS hybrid doublet detection
# 2. Read data & convert to SCE --------------------------------------------
# Read your merged Seurat object
load('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/C03_Integration/integrated_2025-04-17.RData')

integrated

DefaultAssay(integrated)

integrated <- ScaleData(integrated, verbose = FALSE)

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

# 6. Method 3: scds hybrid (CXDS + BCDS)
sce <- cxds(sce)
sce <- bcds(sce)
sce <- cxds_bcds_hybrid(sce, estNdbl = FALSE)  # just compute scores
scores <- colData(sce)$hybrid_score
# take top n_expected_doublet cells as doublets
cutoff <- sort(scores, decreasing=TRUE)[n_expected_doublet]
is_dbl3 <- scores >= cutoff

# 7. Intersection & filtering ----------------------------------------------
#   Identify barcodes called doublets by ALL three methods
doublet_barcodes <- colnames(sce)[ is_dbl1 & is_dbl2 & is_dbl3 ]

#   Subset the original Seurat object to remove those doublets
singlets <- subset(
  integrated,
  cells = setdiff(Cells(integrated), doublet_barcodes)
)

integrated$doublet_status <- ifelse(
  Cells(integrated) %in% doublet_barcodes,
  "doublet",
  "singlet"
)

# ---- visualize method-wise doublet distribution for full integration ----
library(tidyr)
# build a data.frame of flags per method
df_int <- integrated@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    scDblFinder = is_dbl1,
    density     = is_dbl2,
    hybrid      = as.logical(is_dbl3)
  ) %>%
  pivot_longer(
    cols = c("scDblFinder", "density", "hybrid"),
    names_to  = "method",
    values_to = "flag"
  )
# count by sample and method
dist_int <- df_int %>%
  filter(flag) %>%
  group_by(orig.ident, method) %>%
  summarise(n_doublets = n(), .groups="drop")
# barplot
gsea_plot_dir <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal"
p_int <- ggplot(dist_int, aes(x=orig.ident, y=n_doublets, fill=method)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Doublets per Sample by Method (integrated)",
       x="Sample", y="Number of doublets") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
pdf(file.path(gsea_plot_dir,"doublet_dist_integrated.pdf"), width=8, height=4)
print(p_int)
dev.off()
# build named list of doublet cells per method
cells_int <- df_int %>%
  filter(flag)                        # keep only rows where flag is TRUE
venn_list_int <- list(
  scDblFinder = unique(cells_int$cell[cells_int$method == "scDblFinder"]),
  density     = unique(cells_int$cell[cells_int$method == "density"]),
  hybrid      = unique(cells_int$cell[cells_int$method == "hybrid"])
)
# Temporarily change working directory to avoid log file creation in read-only FS
old_wd <- getwd()
tmp_wd <- tempdir()
setwd(tmp_wd)
# Render integrated Venn to PDF via grid
venn_plot_int <- venn.diagram(
  x        = venn_list_int,
  filename = NULL,
  fill     = c("red","green","blue"),
  alpha    = 0.5,
  cat.cex  = 0.8,
  cex      = 1.0,
  main     = "Venn: Methods (integrated)",
  logger   = FALSE      # disable logging
)
# Restore original working directory
setwd(old_wd)
pdf(
  file   = file.path(gsea_plot_dir, "venn_integrated.pdf"),
  width  = 5,
  height = 5
)
grid.newpage()
grid.draw(venn_plot_int)
dev.off()

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
is_dbl1_sc <- colData(sce_cl)$scDblFinder.class == "doublet"

# 5. Method 2: computeDoubletDensity ---------------------------------------
#   k-NN density scoring of doublet likelihood
density_scores <- computeDoubletDensity(sce_cl)

#   Call top N cells as doublets
cutoff2 <- sort(density_scores, decreasing = TRUE)[n_expected_doublet]
is_dbl2_sc <- density_scores >= cutoff2

# 6. Method 3: scds hybrid (CXDS + BCDS) for cell-line subset
sce_cl <- cxds(sce_cl)
sce_cl <- bcds(sce_cl)
sce_cl <- cxds_bcds_hybrid(sce_cl, estNdbl = FALSE)  # just compute scores
scores_cl <- colData(sce_cl)$hybrid_score
# take top n_expected_doublet cells as doublets
cutoff_cl <- sort(scores_cl, decreasing=TRUE)[n_expected_doublet]
is_dbl3_sc <- scores_cl >= cutoff_cl


# 7. Intersection & filtering ----------------------------------------------
#   Identify barcodes called doublets by ALL three methods
doublet_barcodes_sc <- colnames(sce_cl)[ is_dbl1_sc & is_dbl2_sc & is_dbl3_sc ]

#   Subset the original Seurat object to remove those doublets
singlets_cl <- subset(
  integrated_cell_line,
  cells = setdiff(Cells(integrated_cell_line), doublet_barcodes_sc)
)


integrated$doublet_status_cell_line_only <- ifelse(
  Cells(integrated) %in% doublet_barcodes_sc,
  "doublet",
  "singlet"
)

# ---- visualize method-wise doublet distribution for cell-line subset ----
df_cl <- integrated_cell_line@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(
    scDblFinder = is_dbl1_sc,
    density     = is_dbl2_sc,
    hybrid      = as.logical(is_dbl3_sc)
  ) %>%
  pivot_longer(
    cols = c("scDblFinder", "density", "hybrid"),
    names_to  = "method",
    values_to = "flag"
  )
dist_cl <- df_cl %>%
  filter(flag) %>%
  group_by(orig.ident, method) %>%
  summarise(n_doublets = n(), .groups="drop")
p_cl <- ggplot(dist_cl, aes(x=orig.ident, y=n_doublets, fill=method)) +
  geom_bar(stat="identity", position="dodge") +
  labs(title="Doublets per Sample by Method (cell-line only)",
       x="Sample", y="Number of doublets") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45,hjust=1))
pdf(file.path(gsea_plot_dir,"doublet_dist_cellline.pdf"), width=8, height=4)
print(p_cl)
dev.off()
# build named list of doublet cells per method for cell-line subset
cells_cl <- df_cl %>%
  filter(flag)
venn_list_cl <- list(
  scDblFinder = unique(cells_cl$cell[cells_cl$method == "scDblFinder"]),
  density     = unique(cells_cl$cell[cells_cl$method == "density"]),
  hybrid      = unique(cells_cl$cell[cells_cl$method == "hybrid"])
)
# Render cell-line subset Venn to PDF via grid
 # Temporarily change working directory to avoid log file creation in read-only FS
old_wd <- getwd()
tmp_wd <- tempdir()
setwd(tmp_wd)
venn_plot_cl <- venn.diagram(
  x        = venn_list_cl,
  filename = NULL,
  fill     = c("red","green","blue"),
  alpha    = 0.5,
  cat.cex  = 0.8,
  cex      = 1.0,
  main     = "Venn: Methods (cell-line only)",
  logger   = FALSE
)
# Restore original working directory
setwd(old_wd)
pdf(
  file   = file.path(gsea_plot_dir, "venn_cellline.pdf"),
  width  = 5,
  height = 5
)
grid.newpage()
grid.draw(venn_plot_cl)
dev.off()



###### Figures ######

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

# 1. Extract metadata
meta <- integrated@meta.data %>%
  as.data.frame() %>%
  tibble::rownames_to_column("cell")

# 2. Compute per-sample totals and doublet counts/percentages
sum_df <- meta %>%
  group_by(orig.ident) %>%
  summarise(
    total   = n(),
    doublet = sum(doublet_status == "doublet")
  ) %>%
  mutate(
    pct_label = sprintf("%d (%.1f%%)", doublet, 100 * doublet / total)
  )

# 3. Plot 1: stacked bar + annotation of doublet count (percentage)
p1 <- ggplot(meta, aes(x = orig.ident, fill = doublet_status)) +
  geom_bar() +
  geom_text(
    data = sum_df,
    aes(
      x     = orig.ident,
      y     = total,
      label = pct_label
    ),
    inherit.aes = FALSE,   
    vjust = -0.5,
    size  = 3.5
  ) +
  scale_fill_manual(
    values = c("singlet" = "steelblue", "doublet" = "firebrick")
  ) +
  labs(
    title = "Doublet vs Singlet by Sample",
    x     = "Sample (orig.ident)",
    y     = "Number of cells",
    fill  = "Status"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  expand_limits(y = max(sum_df$total) * 1.05)


pdf(
  file   = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/p1_doublet_comparison.pdf",
  width  = 15,
  height = 5
)
print(p1)
dev.off()

#  Prepare data for Plot 2: filter and pivot ------------------------------
meta2 <- meta %>%
  filter(orig.ident %in% c("2N-Cell-Culture", "4N-Cell-Culture")) %>%
  pivot_longer(
    cols      = c(doublet_status, doublet_status_cell_line_only),
    names_to  = "method",
    values_to = "status"
  )

# 3. Summarize counts & percentages ----------------------------------------
sum2 <- meta2 %>%
  group_by(orig.ident, method) %>%
  summarise(
    total   = n(),
    doublet = sum(status == "doublet"),
    .groups = "drop"
  ) %>%
  mutate(
    pct_label = sprintf("%d (%.1f%%)", doublet, 100 * doublet / total)
  )

# 4. Plot 2: stacked bar + annotation --------------------------------------
p2 <- ggplot(meta2, aes(x = method, fill = status)) +
  geom_bar(position = "stack") +
  facet_wrap(~ orig.ident, nrow = 1) +
  # add the count (percent) labels above each bar
  geom_text(
    data = sum2,
    aes(x = method, y = total, label = pct_label),
    inherit.aes = FALSE,
    vjust = -0.5,
    size  = 3.5
  ) +
  scale_fill_manual(
    values = c("singlet" = "steelblue", "doublet" = "firebrick")
  ) +
  labs(
    title = "Comparison of Doublet Calling Methods",
    x = "Method",
    y = "Number of cells",
    fill = "Status"
  ) +
  theme_classic() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "gray90", colour = NA)
  ) +
  # extend y-axis to make room for labels
  expand_limits(y = max(sum2$total) * 1.05)

pdf(
  file   = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/p2_doublet_comparison.pdf",
  width  = 5,
  height = 5
)
print(p2)
dev.off()


#  Save results -----------------------------------------------------------

singlets_cl<-subset(integrated, subset= doublet_status_cell_line_only == 'singlet')

singlets<-subset(integrated, subset= doublet_status == 'singlet')


saveRDS(
  singlets_cl,
  file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets_cl.Rds"
)

saveRDS(
  singlets,
  file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds"
)

saveRDS(
  integrated,
  file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/integrated.Rds"
)

save.image('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/02_doublet‐removal.RData')
