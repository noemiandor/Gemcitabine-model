library(Seurat)
library(msigdbr)
library(scGSVA)
library(dplyr)
library(methods)  # ensure make.names is available
library(tibble)
library(limma)
output_dir<-'/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA'
setwd(output_dir)

# Define a helper function to run scGSVA, replace obj slot, and save results
run_scgsva_and_save <- function(singlets, gs_obj, gs_name, cores = 4, output_dir) {
  # Run scGSVA on the full Seurat object
  gsva_res <- scgsva(
    obj = singlets,
    annot = gs_obj,
    method = "ssgsea",       # or "gsva" / "UCell"
    kcdf = "Gaussian",
    cores = cores,
    verbose = TRUE
  )
  
  # Prepare a "meta-only" Seurat object with no expression data
  orig_meta <- singlets@meta.data
  empty_mat <- matrix(0, nrow = 0, ncol = nrow(orig_meta))
  colnames(empty_mat) <- rownames(orig_meta)
  rownames(empty_mat) <- character(0)
  seurat_meta_only <- CreateSeuratObject(
    counts = empty_mat,
    meta.data = orig_meta,
    project = "MetaOnly"
  )
  
  # Replace the obj slot in a copy to reduce file size
  gsva_res_use <- gsva_res
  gsva_res_use@obj <- seurat_meta_only
  
  # Save both the full GSVA result and the "use" version
  saveRDS(gsva_res, file = file.path(output_dir, paste0(gs_name, "_gsva_res.Rds")))
  saveRDS(gsva_res_use, file = file.path(output_dir, paste0(gs_name, "_gsva_res_use.Rds")))
  
  return(invisible(NULL))
}

# Load the Seurat object once
singlets <- readRDS("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/02_doublet‐removal/singlets.Rds")

# Ensure singlets is log-normalized before GSVA
singlets <- NormalizeData(
  object = singlets,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# Read or define the gene sets
# 1) hsko (already in the correct format)
hsko <- readRDS("/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/hsko.Rds")

# 2) hallmark_sets (create TERM2GENE + annotate)
m_df <- msigdbr(species = "Homo sapiens", category = "H")
term2gene_df <- data.frame(
  GeneID = m_df$gene_symbol,
  PATH = m_df$gs_id,
  Annot = m_df$gs_name,
  stringsAsFactors = FALSE
)
hallmark_sets <- hsko
hallmark_sets@annot <- term2gene_df
hallmark_sets@anntype <- "HALLMARK"
saveRDS(hallmark_sets, file = "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA/hallmark_sets.Rds")

# Define where to save outputs
output_dir <- '/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/ScRNA_Seq/05_GSVA'

# Run GSVA for both gene sets in parallel (two separate calls)
run_scgsva_and_save(singlets, hsko, "hsko", cores = 16, output_dir = output_dir)
run_scgsva_and_save(singlets, hallmark_sets, "hallmark", cores = 16, output_dir = output_dir)

Dose<-read.table('/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/IDs_Dose.txt',header = T)
meta <- singlets@meta.data %>%
  rownames_to_column("cell") %>%                  # preserve barcode
  left_join(Dose, by = c("orig.ident" = "IDs")) %>%
  column_to_rownames("cell")

singlets@meta.data <- meta


singlets$sample_type       <- factor(singlets$sample_type,
                                     levels = c("2N-cellline","4N-cellline","2N-tumor","4N-tumor"))
singlets$Dose              <- factor(singlets$Dose,
                                     levels = c("0mg/kg","30mg/kg","120mg/kg"))
singlets$seurat_clusters   <- factor(singlets$seurat_clusters)





# read back and add GSVA scores to singlets metadata
gsva_res_hsko <- readRDS(file.path(output_dir, "hsko_gsva_res_use.Rds"))
gsva_df_hsko <- as.data.frame(gsva_res_hsko@gsva)
colnames(gsva_df_hsko) <- paste0("GSVA_hsko_", colnames(gsva_df_hsko))
singlets <- AddMetaData(singlets, metadata = gsva_df_hsko)

gsva_res_hm <- readRDS(file.path(output_dir, "hallmarks_gsva_res_use.Rds"))
gsva_df_hm <- as.data.frame(gsva_res_hm@gsva)
colnames(gsva_df_hm) <- paste0("GSVA_", colnames(gsva_df_hm))
singlets <- AddMetaData(singlets, metadata = gsva_df_hm)

saveRDS(singlets, file = file.path(output_dir, "singlets_with_gsva.Rds"))


############################################################
meta_data<-singlets@meta.data

rm(singlets)
gc()#release RAM

saveRDS(meta_data, file = file.path(output_dir, "meta_data.Rds"))

# Create directory for mixed-effects model results
mixed_dir <- file.path(output_dir, "Mixed-Effects_Model")
dir.create(mixed_dir, showWarnings = FALSE)

meta_data<-readRDS(file.path(output_dir, "meta_data.Rds"))

# Function to perform all mixed-effects analyses and save results in mixed_dir
perform_mixed_effects_analysis <- function(meta_data, mixed_dir) {
  library(lme4); library(lmerTest); library(broom.mixed)
  library(tidyr); library(dplyr); library(purrr); library(ggpubr); library(ggplot2)
  
  # Identify GSVA score columns
  gsva_cols <- grep("^GSVA_", colnames(meta_data), value = TRUE)

  # Long-format helper
  to_long <- function(df) {
    df %>%
      rownames_to_column("cell") %>%
      pivot_longer(cols = all_of(gsva_cols), names_to = "pathway", values_to = "score") %>%
      mutate(sample = orig.ident)
  }

  # Helper to run one mixed-effects test and produce boxplots
  run_and_plot <- function(df, group_col, ref_level, prefix) {
    df <- df %>% filter(!is.na(.data[[group_col]]))
    df[[group_col]] <- relevel(factor(df[[group_col]]), ref_level)
    long_df <- to_long(df)

    formula_str <- paste0("score ~ ", group_col, " + (1|sample)")
    res <- long_df %>%
      group_by(pathway) %>%
      nest() %>%
      mutate(
        model = map(data, ~ lmer(as.formula(formula_str), data = .x,
                                 control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))),
        tidy = map(model, ~ if (inherits(.x, "lmerMod")) tidy(.x, effects="fixed", conf.int=TRUE) else NULL)
      ) %>%
      select(pathway, tidy) %>%
      unnest(tidy) %>%
      filter(term != "(Intercept)")

    # Save results
    write.csv(res, file = file.path(mixed_dir, paste0(prefix, "_mixed_effects_all.csv")), row.names = FALSE)
    sig <- res %>% filter(p.value < 0.05)
    write.csv(sig, file = file.path(mixed_dir, paste0(prefix, "_mixed_effects_sig.csv")), row.names = FALSE)

    # Boxplots
    pdf(file.path(mixed_dir, paste0(prefix, "_boxplots.pdf")), width = 6, height = 4)
    for (pw in unique(sig$pathway)) {
      df_plot <- df %>%
        rownames_to_column("cell") %>%
        select(cell, pathway_score = !!sym(pw), group = .data[[group_col]]) %>%
        filter(!is.na(group))
      plt <- ggplot(df_plot, aes(x = group, y = pathway_score)) +
        geom_boxplot() +
        stat_compare_means(method = "wilcox.test") +
        labs(title = paste(prefix, pw), x = group_col, y = "GSVA score")
      print(plt)
    }
    dev.off()
  return(res)
    }

  # Define scenarios
  scenarios <- list(
    list(df = meta_data %>% filter(sample_type %in% c("2N-cellline","4N-cellline")), group="sample_type", ref="2N-cellline", prefix="type_4N_vs_2N_cellline"),
    list(df = meta_data %>% filter(sample_type %in% c("2N-tumor","4N-tumor")), group="sample_type", ref="2N-tumor", prefix="type_4N_vs_2N_tumor"),
    list(df = { meta_data$sample_group <- factor(meta_data$sample_group, levels=c("2N","4N")); meta_data }, group="sample_group", ref="2N", prefix="sample_group_4N_vs_2N"),
    list(df = meta_data, group="Dose", ref="0mg/kg", prefix="dose_others_vs_0"),
    list(df = meta_data %>% filter(Dose %in% c("30mg/kg","120mg/kg")), group="Dose", ref="30mg/kg",prefix="dose_120_vs_30"),
    list(df = meta_data %>% filter(sample_type=="2N-tumor"), group="Dose", ref="0mg/kg",  prefix="2T_dose_others_vs_0"),
    list(df = meta_data %>% filter(sample_type=="2N-tumor", Dose %in% c("30mg/kg","120mg/kg")), group="Dose",  ref="30mg/kg",  prefix="2T_dose_120_vs_30"),
    list(df = meta_data %>% filter(sample_type=="4N-tumor"), group="Dose", ref="0mg/kg", prefix="4T_dose_others_vs_0"),
    list(df = meta_data %>% filter(sample_type=="4N-tumor", Dose %in% c("30mg/kg","120mg/kg")), group="Dose",  ref="30mg/kg", prefix="4T_dose_120_vs_30")
  )

  # Run all
  lapply(scenarios, function(s) run_and_plot(s$df, s$group, s$ref, s$prefix))
}


# Execute mixed-effects workflow
perform_mixed_effects_analysis(meta_data, mixed_dir)


# ----------------------------------------------------------------------
# True Pseudo-Bulk GSVA Differential Analysis
# ----------------------------------------------------------------------

# Create directory for pseudo-bulk results
pseudo_dir <- file.path(output_dir, "psudo-bulk")
dir.create(pseudo_dir, showWarnings = FALSE)

library(dplyr)
library(limma)
library(dplyr)
library(stringr)
library(methods)

meta_data_fixed <- meta_data %>%
  mutate(
    across(
      where(is.character),
      ~ str_replace_all(.x, "-", ".")
    ),
    across(
      where(is.factor),
      ~ factor(str_replace_all(as.character(.x), "-", "."))
    )
  )

# Identify all GSVA score columns (both hsko and hallmark)
gsva_cols <- grep("^GSVA", colnames(meta_data_fixed), value = TRUE)

# Helper to run pseudobulk for a given subset and grouping
run_pseudobulk <- function(df, group_col, levels_vec, prefix) {
  # Sanitize grouping values to syntactically valid names
  df[[group_col]] <- make.names(as.character(df[[group_col]]))
  san_levels <- make.names(levels_vec)
  df[[group_col]] <- factor(df[[group_col]], levels = san_levels)

  # 2) Aggregate (mean) per sample (orig.ident)
  pseudo <- df %>%
    select(orig.ident, all_of(gsva_cols), !!sym(group_col)) %>%
    group_by(orig.ident, !!sym(group_col)) %>%
    summarise(across(all_of(gsva_cols), mean), .groups = "drop")

  # 3) Sample annotation and design
  sample_anno <- pseudo %>%
    select(orig.ident, !!sym(group_col)) %>%
    distinct() %>%
    column_to_rownames("orig.ident")
  design <- model.matrix(~ 0 + ., data = sample_anno)
  # Ensure design column names match sanitized levels
  colnames(design) <- san_levels
  # colnames(design) <- levels_vec

  # 4) Expression matrix: pathways × samples
  expr_mat <- t(as.matrix(pseudo %>% column_to_rownames("orig.ident") %>% select(all_of(gsva_cols))))
  colnames(expr_mat) <- rownames(sample_anno)

  # 5) Define pairwise contrasts using sanitized level names
  lvl <- san_levels
  contrast_list <- combn(lvl, 2, function(x) paste0(x[1], "-", x[2]), simplify = FALSE)
  contrasts <- setNames(
    lapply(combn(lvl, 2, simplify = FALSE),
           function(x) paste0(x[1], "-", x[2])),
    contrast_list
  )
  cont_matrix <- makeContrasts(contrasts = unlist(contrasts), levels = design)

  # 6) Fit limma and eBayes
  fit <- lmFit(expr_mat, design)
  fit <- contrasts.fit(fit, cont_matrix)
  fit <- eBayes(fit)

  # 7) Assemble results
  res <- data.frame(Pathway = rownames(fit$coefficients), stringsAsFactors = FALSE)
  for (ctr in names(contrasts)) {
    res[[paste0("logFC_", ctr)]]   <- fit$coefficients[, ctr]
    res[[paste0("Pvalue_", ctr)]]  <- fit$p.value[, ctr]
  }
  res <- res %>%
    rename_with(~ .x %>% 
                  gsub("X", "", .) %>%     
                  gsub("\\.", "/", .)     
    )

  # 8) Save full and significant
  write.csv(res, file = file.path(pseudo_dir, paste0(prefix, "_all.csv")), row.names = FALSE)
  sig <- res %>% filter(if_any(starts_with("Pvalue_"), ~ . < 0.05))
  write.csv(sig, file = file.path(pseudo_dir, paste0(prefix, "_sig.csv")), row.names = FALSE)
}

# Run pseudobulk for each scenario, preserving the same comparisons:
# 2N-cellline vs 2N-tumor (all doses)
run_pseudobulk(
  meta_data_fixed %>% filter(sample_type %in% c("2N.cellline", "2N.tumor")),
  group_col = "sample_type",
  levels_vec = c("2N.cellline", "2N.tumor"),
  prefix = "2N.cellline_vs_2N.tumor"
)

# 4N-cellline vs 4N-tumor (all doses)
run_pseudobulk(
  meta_data_fixed %>% filter(sample_type %in% c("4N.cellline", "4N.tumor")),
  group_col = "sample_type",
  levels_vec = c("4N.cellline", "4N.tumor"),
  prefix = "4N.cellline_vs_4N.tumor"
)

# 2N-tumor dose comparisons
run_pseudobulk(
  meta_data_fixed %>% filter(sample_type == "2N.tumor"),
  group_col = "Dose",
  levels_vec = c("0mg/kg", "30mg/kg", "120mg/kg"),
  prefix = "2N_tumor_dose"
)


# 4N-tumor dose comparisons
run_pseudobulk(
  meta_data_fixed %>% filter(sample_type == "4N.tumor"),
  group_col = "Dose",
  levels_vec = c("0mg/kg", "30mg/kg", "120mg/kg"),
  prefix = "4N_tumor_dose"
)


# ----------------------------------------------------------------------
# Bootstrap Down-sampling Analysis (1000 cells × 1000 iterations)
# ----------------------------------------------------------------------
library(purrr)
library(dplyr)

# Create directory for bootstrap results
bootstrap_dir <- file.path(output_dir, "bootstrap_results")
dir.create(bootstrap_dir, showWarnings = FALSE)

# Parallel setup
library(parallel)
ncores <- max(1, parallel::detectCores() - 1)

# Define bootstrap scenarios matching pseudobulk prefixes
bs_scenarios <- list(
  list(prefix = "2N_cellline_vs_tumor", df = meta_data_fixed %>% filter(sample_type %in% c("2N.cellline","2N.tumor")), group_col = "sample_type", levels = c("2N.cellline","2N.tumor")),
  list(prefix = "4N_cellline_vs_tumor", df = meta_data_fixed %>% filter(sample_type %in% c("4N.cellline","4N.tumor")), group_col = "sample_type", levels = c("4N.cellline","4N.tumor")),
  list(prefix = "2N_tumor_dose",        df = meta_data_fixed %>% filter(sample_type == "2N.tumor"),                              group_col = "Dose",        levels = c("0mg/kg","30mg/kg","120mg/kg")),
  list(prefix = "4N_tumor_dose",        df = meta_data_fixed %>% filter(sample_type == "4N.tumor"),                              group_col = "Dose",        levels = c("0mg/kg","30mg/kg","120mg/kg"))
)

# Run bootstrap for each scenario
for (sc in bs_scenarios) {
  df <- sc$df
  # Precompute valid cells per group
  cells_by_group <- split(rownames(df), df[[sc$group_col]])
  # Initialize counts storage
  pathway_counts <- setNames(rep(0, length(gsva_cols)), gsva_cols)

  # Run bootstrap iterations in parallel, each returning a logical vector of significant flags
  sig_mat <- parallel::mclapply(
    seq_len(1000),
    function(i) {
      # Sample 1000 cells per group (with replacement if needed)
      samp_cells <- lapply(cells_by_group, function(cg) sample(cg, size = min(1000, length(cg)), replace = length(cg) < 1000))
      # Build sampled dataframe
      df_samp <- df[unlist(samp_cells), ]
      # Test each pathway
      sapply(gsva_cols, function(pw) {
        g1 <- df_samp[samp_cells[[1]], pw]
        g2 <- df_samp[samp_cells[[2]], pw]
        wilcox.test(g1, g2)$p.value < 0.05
      })
    },
    mc.cores = ncores
  )
  # Sum across iterations to count significant calls per pathway
  pathway_counts <- colSums(do.call(rbind, sig_mat))

  # Save counts: number of times p < 0.05 out of 1000
  bs_res <- data.frame(
    Pathway = names(pathway_counts),
    SignificantCount = as.integer(pathway_counts),
    stringsAsFactors = FALSE
  )
  write.csv(bs_res, file = file.path(bootstrap_dir, paste0(sc$prefix, "_bootstrap_counts.csv")), row.names = FALSE)
}


save.image(paste0(output_dir,'/05_GSVA.RData'))
