#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Seurat", "scDblFinder", "SingleCellExperiment", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05a_species_extraction", "tentative_mouse_seurat_list_pre_qc.rds")
out_dir <- file.path(results_root, "05b_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05a input: ", input_rds, call. = FALSE)

qc <- cfg$qc
min_features <- as.numeric(qc$min_features)
min_counts <- as.numeric(qc$min_counts)
max_percent_mt <- as.numeric(qc$max_percent_mt)
max_features_quantile <- as.numeric(qc$max_features_quantile)
max_counts_quantile <- as.numeric(qc$max_counts_quantile)
scdblfinder_min_cells <- as.integer(qc$scdblfinder_min_cells)
seed <- as.integer(qc$seed)
mouse_fraction_min <- as.numeric(cfg$species$mouse_fraction_min)

mouse_list <- readRDS(input_rds)
if (!is.list(mouse_list) || length(mouse_list) == 0L) stop("05a input is not a non-empty Seurat list.", call. = FALSE)

post_qc <- list()
decision_tables <- list()
summary_tables <- list()
qc_metric_tables <- list()

for (sample_id in names(mouse_list)) {
  obj <- mouse_list[[sample_id]]
  if (!inherits(obj, "Seurat")) stop("Non-Seurat object in 05a list: ", sample_id, call. = FALSE)
  if (any(obj$mouse_fraction < mouse_fraction_min, na.rm = TRUE)) {
    stop("05a threshold contract violated in sample: ", sample_id, call. = FALSE)
  }
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, assay = "RNA", pattern = "^mt-")
  md <- obj@meta.data
  feature_upper <- as.numeric(stats::quantile(md$nFeature_RNA, probs = max_features_quantile, na.rm = TRUE, names = FALSE))
  count_upper <- as.numeric(stats::quantile(md$nCount_RNA, probs = max_counts_quantile, na.rm = TRUE, names = FALSE))

  fail_min_features <- md$nFeature_RNA < min_features
  fail_min_counts <- md$nCount_RNA < min_counts
  fail_percent_mt <- md$percent.mt > max_percent_mt
  fail_max_features <- md$nFeature_RNA > feature_upper
  fail_max_counts <- md$nCount_RNA > count_upper
  pass_qc <- !(fail_min_features | fail_min_counts | fail_percent_mt | fail_max_features | fail_max_counts)

  decision <- data.frame(
    cell = rownames(md), sample = sample_id,
    nFeature_RNA = md$nFeature_RNA, nCount_RNA = md$nCount_RNA,
    percent.mt = md$percent.mt, mouse_fraction = md$mouse_fraction,
    feature_upper_99pct = feature_upper, count_upper_99pct = count_upper,
    fail_min_features = fail_min_features, fail_min_counts = fail_min_counts,
    fail_percent_mt = fail_percent_mt, fail_max_features = fail_max_features,
    fail_max_counts = fail_max_counts, pass_qc = pass_qc,
    stringsAsFactors = FALSE
  )
  qc_metric_tables[[sample_id]] <- decision
  passing_cells <- rownames(md)[pass_qc]
  n_after_qc <- length(passing_cells)

  if (n_after_qc == 0L) {
    obj <- NULL
    decision$scDblFinder.class <- NA_character_
    decision$scDblFinder.score <- NA_real_
    decision$keep_final <- FALSE
  } else {
    obj <- subset(obj, cells = passing_cells)
  }

  if (n_after_qc > 0L && n_after_qc < scdblfinder_min_cells) {
    obj$scDblFinder.class <- "singlet_low_cell_skip"
    obj$scDblFinder.score <- NA_real_
    decision$scDblFinder.class <- ifelse(decision$pass_qc, "singlet_low_cell_skip", NA_character_)
    decision$scDblFinder.score <- NA_real_
    decision$keep_final <- decision$pass_qc
  } else if (n_after_qc >= scdblfinder_min_cells) {
    sce <- SingleCellExperiment::SingleCellExperiment(
      assays = list(counts = Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
    )
    set.seed(seed)
    sce <- scDblFinder::scDblFinder(sce, verbose = FALSE)
    dbl <- as.data.frame(SingleCellExperiment::colData(sce))
    obj$scDblFinder.class <- as.character(dbl$scDblFinder.class)
    obj$scDblFinder.score <- as.numeric(dbl$scDblFinder.score)
    decision$scDblFinder.class <- NA_character_
    decision$scDblFinder.score <- NA_real_
    hit <- match(rownames(dbl), decision$cell)
    decision$scDblFinder.class[hit] <- as.character(dbl$scDblFinder.class)
    decision$scDblFinder.score[hit] <- as.numeric(dbl$scDblFinder.score)
    decision$keep_final <- decision$pass_qc & decision$scDblFinder.class == "singlet"
    obj <- subset(obj, cells = rownames(dbl)[dbl$scDblFinder.class == "singlet"])
  }

  decision_tables[[sample_id]] <- decision
  n_final <- if (is.null(obj)) 0L else ncol(obj)
  summary_tables[[sample_id]] <- data.frame(
    sample = sample_id,
    tentative_mouse_pre_qc = nrow(md),
    after_qc = n_after_qc,
    after_scdblfinder = n_final,
    removed_by_qc = nrow(md) - n_after_qc,
    removed_as_doublet = n_after_qc - n_final,
    feature_upper_99pct = feature_upper,
    count_upper_99pct = count_upper,
    scdblfinder_skipped_low_cells = n_after_qc > 0L && n_after_qc < scdblfinder_min_cells,
    stringsAsFactors = FALSE
  )
  if (n_final > 0L) post_qc[[sample_id]] <- obj
  rm(obj)
  invisible(gc())
}

if (length(post_qc) < 2L) stop("Fewer than two non-empty samples remain after mouse-cell QC; CCA cannot proceed.", call. = FALSE)
decision_df <- do.call(rbind, decision_tables)
summary_df <- do.call(rbind, summary_tables)
qc_metrics_df <- do.call(rbind, qc_metric_tables)
rownames(decision_df) <- rownames(summary_df) <- rownames(qc_metrics_df) <- NULL

if (any(vapply(post_qc, function(x) any(x$mouse_fraction < mouse_fraction_min, na.rm = TRUE), logical(1)))) {
  stop("A post-QC cell violates the mouse-fraction threshold.", call. = FALSE)
}

utils::write.csv(decision_df, file.path(out_dir, "cell_qc_decisions.csv"), row.names = FALSE, na = "NA")
utils::write.csv(summary_df, file.path(out_dir, "sample_qc_summary.csv"), row.names = FALSE, na = "NA")
utils::write.csv(qc_metrics_df, file.path(out_dir, "qc_metrics_pre_filter.csv"), row.names = FALSE, na = "NA")
saveRDS(post_qc, file.path(out_dir, "mouse_seurat_list_post_qc.rds"), compress = FALSE)

plot_df <- qc_metrics_df
plot_long <- rbind(
  data.frame(sample = plot_df$sample, metric = "nFeature_RNA", value = plot_df$nFeature_RNA),
  data.frame(sample = plot_df$sample, metric = "nCount_RNA", value = plot_df$nCount_RNA),
  data.frame(sample = plot_df$sample, metric = "percent.mt", value = plot_df$percent.mt)
)
p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = sample, y = value)) +
  ggplot2::geom_violin(scale = "width", fill = "grey85", color = "grey35") +
  ggplot2::geom_boxplot(width = 0.12, outlier.shape = NA) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
ggplot2::ggsave(file.path(out_dir, "mouse_qc_metrics_pre_filter.pdf"), p, width = 12, height = 12)
ggplot2::ggsave(file.path(out_dir, "mouse_qc_metrics_pre_filter.png"), p, width = 12, height = 12, dpi = 220)

writeLines(
  c(
    "05b mouse-cell QC completed.",
    paste0("Samples retained: ", length(post_qc)),
    paste0("Cells before QC: ", sum(summary_df$tentative_mouse_pre_qc)),
    paste0("Cells after QC and scDblFinder: ", sum(summary_df$after_scdblfinder))
  ),
  file.path(out_dir, "completion.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05b.txt"))
message("[05b] Completed. Output: ", out_dir)
