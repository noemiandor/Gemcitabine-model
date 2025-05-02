# ───────────────────────────────────────────────────────────────────────────
# Complete R script: 
#  • GSEA on MSigDB Hallmark, Cell Cycle & Stress gene sets 
#  • Save FGSEA result tables
#  • DSEA (Disease Ontology GSEA via DOSE::gseDO) 
#  • Save DSEA result tables
#  • Heatmap of NES (FGSEA) → PDF
#  • Dose‐response NES line plot → PDF
# ───────────────────────────────────────────────────────────────────────────
library(msigdbr)    # fetch MSigDB
library(fgsea)      # fast GSEA
library(DOSE)       # gseDO for Disease Ontology
library(pheatmap)   # NES heatmap
library(ggplot2)    # plotting
library(dplyr)      # data manipulation
library(tibble)     # rowname utilities
library(tidyr)      # data reshaping

# 1. Define directories -----------------------------------------------------
dose_dir       <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/03_DiffExp_genes/DE_results"
gsea_res_dir   <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/04_GSEA/GSEA_Results"
dsea_res_dir   <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/04_GSEA/DSEA_Results"
plot_dir       <- "/Users/4482173/Documents/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/04_GSEA/GSEA_Plots"
dir.create(gsea_res_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(dsea_res_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(plot_dir,       recursive=TRUE, showWarnings=FALSE)

# 2. Read DE result files & build ranked lists ------------------------------
dose_files <- list.files(dose_dir, pattern="^DE_dose_.*\\.csv$", full.names=TRUE)
gene_lists <- lapply(dose_files, function(f) {
  df <- read.csv(f, row.names=1, stringsAsFactors=FALSE)
  ranks <- df$avg_log2FC
  names(ranks) <- rownames(df)
  sort(ranks, decreasing=TRUE)
})
names(gene_lists) <- dose_files %>%
  basename() %>%
  gsub("^DE_dose_|\\.csv$", "", .)

# 3. Retrieve MSigDB gene sets ------------------------------------------------
# 3.1 Hallmark (category H)
hallmark_sets <- msigdbr(species="Homo sapiens", category="H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# 3.2 GO cell cycle & stress response (C5 GO)
# fetch only GO Biological Process gene sets
go_bp <- msigdbr(
  species    = "Homo sapiens", 
  category   = "C5", 
  subcategory= "BP"
)  

# then filter to the terms and split into a list
go_sets <- go_bp %>%
  filter(grepl("CELL_CYCLE|RESPONSE_TO_STRESS", gs_name)) %>%
  split(x = .$gene_symbol, f = .$gs_name)

# combine all
pathways <- c(hallmark_sets, go_sets)

# 4. Run fgsea for each comparison -------------------------------------------
gsea_results <- lapply(gene_lists, function(ranks) {
  fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = 15,
    maxSize  = 500,
    nperm    = 10000
  )
})

# Save raw GSEA tables to CSV --------------------------------------------
for (comp in names(gsea_results)) {
  # pull out the fgsea result and coerce to data.frame
  df <- as.data.frame(gsea_results[[comp]])
  
  # find any list‐columns
  is_list_col <- sapply(df, is.list)
  if (any(is_list_col)) {
    # collapse each list element to a single string
    df[is_list_col] <- lapply(
      df[is_list_col],
      function(col) sapply(col, function(x) paste(x, collapse = ";"))
    )
  }
  
  # write it out
  write.csv(
    df,
    file      = file.path(gsea_res_dir, paste0("GSEA_table_", comp, ".csv")),
    row.names = FALSE,
    quote     = TRUE
  )
}


# 5. Combine into long format ------------------------------------------------
gsea_long <- bind_rows(
  lapply(names(gsea_results), function(comp) {
    res <- gsea_results[[comp]]
    res %>% select(pathway, NES, padj) %>% mutate(comparison = comp)
  })
)

# 6. Build NES matrix & plot heatmap -----------------------------------------
nes_matrix <- gsea_long %>%
  select(pathway, comparison, NES) %>%
  pivot_wider(names_from = comparison, values_from = NES) %>%
  column_to_rownames("pathway") %>%
  as.matrix()
# identify rows with all finite values
good_rows <- apply(nes_matrix, 1, function(x) all(is.finite(x)))

# subset to only those
nes_matrix_clean <- nes_matrix[good_rows, ]

# now plot
p1<-pheatmap(
  nes_matrix_clean,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main         = "GSEA NES Heatmap (filtered)",
  fontsize_row = 8
)
# save heatmap as PDF
pdf(
  file   = file.path(gsea_plot_dir, "GSEA_NES_Heatmap.pdf"),
  width  = 8,
  height = 10
)
print(p1)
dev.off()

# 7. Prepare dose-response data & plot line graph ----------------------------
# select pathways to visualize
selected_pathways <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "GO_CELL_CYCLE",
  "GO_RESPONSE_TO_STRESS"
)

plot_df <- gsea_long %>%
  filter(pathway %in% selected_pathways) %>%
  mutate(
    comparison = factor(comparison, levels = names(gene_lists))
  )

dose_response_plot <- ggplot(plot_df, aes(x = comparison, y = NES, color = pathway, group = pathway)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Dose-Response of GSEA NES",
    x     = "Comparison",
    y     = "Normalized Enrichment Score (NES)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# save dose-response plot as PDF
ggsave(
  filename = file.path(gsea_plot_dir, "Dose_Response_GSEA_NES.pdf"),
  plot     = dose_response_plot,
  width    = 8,
  height   = 5,
  units    = "in"
)