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
required_packages <- c("yaml", "Seurat", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05e_broad_annotation", "mouse_clustered_broad_annotated.rds")
marker_file <- file.path(results_root, "05d_clustering", "cluster_positive_markers.csv")
out_dir <- file.path(results_root, "05f_immune_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05e input: ", input_rds, call. = FALSE)
if (!file.exists(marker_file)) stop("Missing 05d cluster marker table: ", marker_file, call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05e input is not a Seurat object.", call. = FALSE)
if (!all(c("seurat_clusters", "mouse_broad_cell_type") %in% colnames(obj@meta.data))) {
  stop("05e object lacks cluster-level broad annotation.", call. = FALSE)
}
cluster_broad_n <- tapply(obj$mouse_broad_cell_type, as.character(obj$seurat_clusters), function(x) length(unique(x)))
if (any(cluster_broad_n != 1L)) stop("05e broad annotation is not uniform within every cluster.", call. = FALSE)
immune_clusters <- unique(as.character(obj$seurat_clusters)[obj$mouse_broad_cell_type == "Immune"])
if (length(immune_clusters) == 0L) stop("No clusters were annotated as Immune by 05e.", call. = FALSE)

markers <- utils::read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
if (!all(c("cluster", "gene") %in% colnames(markers)) || length(fc_column) != 1L) {
  stop("05d marker table must contain cluster, gene, and exactly one avg_log2FC/avg_logFC column.", call. = FALSE)
}
markers$cluster <- as.character(markers$cluster)
markers$gene_upper <- toupper(trimws(as.character(markers$gene)))
markers$rank_fc <- as.numeric(markers[[fc_column]])
markers <- markers[nzchar(markers$gene_upper) & is.finite(markers$rank_fc), , drop = FALSE]

top_marker_n <- as.integer(cfg$annotation$top_cluster_markers)
if (!is.finite(top_marker_n) || top_marker_n < 1L) stop("annotation.top_cluster_markers must be a positive integer.", call. = FALSE)

immune_signatures <- list(
  `T cell` = c(
    "Cd3d", "Cd3e", "Cd3g", "Trac", "Trbc1", "Trbc2", "Cd247", "Lck", "Lat", "Skap1",
    "Itk", "Bcl11b", "Tcf7", "Il7r", "Cd2", "Cd5", "Cd6", "Themis", "Mal", "Ltb"
  ),
  `NK cell` = c(
    "Nkg7", "Klrd1", "Klrk1", "Prf1", "Gzmb", "Gzma", "Ccl5", "Xcl1", "Xcl2", "Klrc1",
    "Klra4", "Klra7", "Klra8", "Klra9", "Klra10", "Klra3", "Eomes", "Tbx21", "Fasl", "S1pr5"
  ),
  `B cell` = c(
    "Cd79a", "Cd79b", "Ms4a1", "Cd37", "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cd22", "Cd19",
    "Bank1", "Pax5", "Ebf1", "Cd72", "Cd180", "Fcmr", "Cd83", "Blnk", "Spib", "Hvcn1"
  ),
  `Plasma cell` = c(
    "Jchain", "Mzb1", "Sdc1", "Xbp1", "Derl3", "Sec11c", "Ighm", "Igha", "Ighg1", "Ighg2b",
    "Igkc", "Slpi", "Creld2", "Pdia4", "Manf", "Txndc5", "Ero1l", "Hspa5", "Tent5c", "Pou2af1"
  ),
  Monocyte = c(
    "Ly6c2", "Ccr2", "Sell", "Plac8", "Chil3", "Fn1", "Vcan", "Ms4a6c", "Lyz2", "Ctss",
    "S100a10", "Fcgr1", "Lst1", "Tyrobp", "Fcer1g", "Ctsb", "Cst3", "Ifitm3", "Lgals3", "Ccr1"
  ),
  Macrophage = c(
    "Adgre1", "Csf1r", "Mertk", "Fcgr1", "Cd68", "Mafb", "Apoe", "C1qa", "C1qb", "C1qc",
    "Lgals3", "Ms4a7", "Trem2", "Marco", "Mrc1", "Cd163", "Folr2", "Ctsd", "Ctsb", "Lpl",
    "Acp5", "Cd86", "Aif1", "Lgmn"
  ),
  Neutrophil = c(
    "S100a8", "S100a9", "Ly6g", "Csf3r", "Retnlg", "Ngp", "Lcn2", "Mmp8", "Mmp9", "Camp",
    "Wfdc21", "Wfdc17", "Il1b", "Fpr1", "Fpr2", "Cxcr2", "Clec4d", "Clec4e", "Slpi", "Stfa2l1",
    "G0s2", "Hcar2", "Mreg", "Fnip2", "Lsmem1", "Dock10", "P2rx7", "Nceh1", "Il1r2", "Chil3",
    "Cxcr4", "Nfam1", "Adgre5", "Ccrl2", "Nlrp3", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl3", "Ccl4",
    "Cstdc4", "Ccl6", "Cd14", "Saa1", "Saa3", "Hp", "Chil1", "Lrg1", "Rsad2", "Ifit1",
    "Isg15", "Slfn4", "Gbp2", "Plac8", "Morrbid", "Retreg1", "Fmnl2", "Lyst", "Trem1", "Mcemp1",
    "Alox5ap", "Adam8", "Mrgpra2b", "Pglyrp1", "Ltf", "Prok2", "Cd177", "Serpinb1a"
  ),
  `Dendritic cell` = c(
    "Itgax", "Flt3", "Zbtb46", "Xcr1", "Clec9a", "Ccr7", "Fscn1", "Cd209a", "Clec10a", "Sirpa",
    "Cd74", "H2-Aa", "H2-Ab1", "H2-Eb1", "Cst3", "Batf3", "Irf8", "Relb", "Traf1", "Il12b"
  ),
  `Mast cell` = c(
    "Kit", "Fcer1a", "Cpa3", "Mcpt4", "Tpsb2", "Ms4a2", "Gata2", "Hdc", "Srgn", "Cma1",
    "Mcpt1", "Mcpt2", "Mcpt5", "Mcpt6", "Il1rl1", "Rgs13", "Socs2", "Cd200r3", "Gzmb", "Alox5ap"
  )
)

expressed_upper <- unique(toupper(rownames(obj[["RNA"]])))
immune_signatures <- lapply(immune_signatures, function(x) intersect(unique(toupper(x)), expressed_upper))
if (any(vapply(immune_signatures, length, integer(1)) == 0L)) stop("At least one immune signature has no genes in the RNA assay.", call. = FALSE)

make_top_markers <- function(cluster_id) {
  z <- markers[markers$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$rank_fc, z$gene_upper), , drop = FALSE]
  z <- z[!duplicated(z$gene_upper), , drop = FALSE]
  z <- utils::head(z, top_marker_n)
  if (nrow(z) == 0L) stop("Immune cluster has no positive markers: ", cluster_id, call. = FALSE)
  data.frame(
    cluster = cluster_id, marker_rank = seq_len(nrow(z)), gene = z$gene,
    gene_upper = z$gene_upper, avg_log2FC = z$rank_fc, stringsAsFactors = FALSE
  )
}
top_markers <- do.call(rbind, lapply(immune_clusters, make_top_markers))
rownames(top_markers) <- NULL

score_one <- function(cluster_id, signature_name) {
  marker_genes <- unique(top_markers$gene_upper[top_markers$cluster == cluster_id])
  sig_genes <- immune_signatures[[signature_name]]
  overlap_genes <- intersect(marker_genes, sig_genes)
  overlap <- length(overlap_genes)
  precision <- overlap / length(marker_genes)
  recall <- overlap / length(sig_genes)
  f1 <- if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  union_len <- length(unique(c(marker_genes, sig_genes)))
  data.frame(
    cluster = cluster_id, signature = signature_name, cluster_marker_n = length(marker_genes),
    signature_gene_n = length(sig_genes), overlap = overlap, precision = precision, recall = recall,
    f1 = f1, jaccard = if (union_len == 0L) 0 else overlap / union_len,
    overlap_genes = paste(overlap_genes, collapse = ";"), stringsAsFactors = FALSE
  )
}
score_grid <- expand.grid(cluster = immune_clusters, signature = names(immune_signatures), stringsAsFactors = FALSE)
scores <- do.call(rbind, Map(score_one, score_grid$cluster, score_grid$signature))
rownames(scores) <- NULL

select_cluster_label <- function(cluster_id) {
  z <- scores[scores$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$f1, -z$overlap, -z$jaccard, z$signature), , drop = FALSE]
  best <- z[1L, , drop = FALSE]
  second <- if (nrow(z) >= 2L) z[2L, , drop = FALSE] else best[0L, , drop = FALSE]
  label <- if (best$overlap[[1]] > 0L) best$signature[[1]] else "Ambiguous"
  confidence <- if (best$overlap[[1]] == 0L) {
    "ambiguous"
  } else if (best$overlap[[1]] >= 3L && (nrow(second) == 0L || best$f1[[1]] > second$f1[[1]])) {
    "high"
  } else if (best$overlap[[1]] >= 2L) {
    "medium"
  } else {
    "low"
  }
  data.frame(
    cluster = cluster_id, immune_cell_type = label, winning_signature = best$signature[[1]],
    overlap = best$overlap[[1]], precision = best$precision[[1]], recall = best$recall[[1]],
    f1 = best$f1[[1]], jaccard = best$jaccard[[1]], overlap_genes = best$overlap_genes[[1]],
    runner_up_signature = if (nrow(second) == 0L) NA_character_ else second$signature[[1]],
    runner_up_f1 = if (nrow(second) == 0L) NA_real_ else second$f1[[1]],
    f1_margin = if (nrow(second) == 0L) NA_real_ else best$f1[[1]] - second$f1[[1]],
    confidence = confidence, stringsAsFactors = FALSE
  )
}
cluster_annotation <- do.call(rbind, lapply(immune_clusters, select_cluster_label))
cluster_annotation$n_cells <- as.integer(table(as.character(obj$seurat_clusters))[cluster_annotation$cluster])
rownames(cluster_annotation) <- NULL

obj$mouse_immune_candidate <- as.character(obj$seurat_clusters) %in% immune_clusters
obj$mouse_immune_cell_type <- NA_character_
obj$mouse_immune_subtype <- NA_character_
obj$mouse_immune_confidence <- NA_character_
obj$mouse_immune_signature <- NA_character_
obj$mouse_immune_overlap <- NA_integer_
obj$mouse_immune_f1 <- NA_real_
obj$mouse_immune_overlap_genes <- NA_character_

immune_index <- which(obj$mouse_immune_candidate)
immune_cluster_vector <- as.character(obj$seurat_clusters)[immune_index]
lookup <- function(column) stats::setNames(cluster_annotation[[column]], cluster_annotation$cluster)[immune_cluster_vector]
obj$mouse_immune_cell_type[immune_index] <- unname(lookup("immune_cell_type"))
obj$mouse_immune_subtype[immune_index] <- unname(lookup("immune_cell_type"))
obj$mouse_immune_confidence[immune_index] <- unname(lookup("confidence"))
obj$mouse_immune_signature[immune_index] <- unname(lookup("winning_signature"))
obj$mouse_immune_overlap[immune_index] <- as.integer(unname(lookup("overlap")))
obj$mouse_immune_f1[immune_index] <- as.numeric(unname(lookup("f1")))
obj$mouse_immune_overlap_genes[immune_index] <- unname(lookup("overlap_genes"))

immune_labels_per_cluster <- tapply(obj$mouse_immune_cell_type[immune_index], immune_cluster_vector, function(x) length(unique(x)))
if (any(immune_labels_per_cluster != 1L)) stop("Immune annotation is not uniform within every immune cluster.", call. = FALSE)

signature_library <- do.call(rbind, lapply(names(immune_signatures), function(nm) {
  data.frame(signature = nm, gene_upper = immune_signatures[[nm]], stringsAsFactors = FALSE)
}))
utils::write.csv(signature_library, file.path(out_dir, "immune_signature_library.csv"), row.names = FALSE)
utils::write.csv(top_markers, file.path(out_dir, "immune_cluster_top_markers.csv"), row.names = FALSE)
utils::write.csv(scores, file.path(out_dir, "immune_cluster_signature_scores.csv"), row.names = FALSE)
utils::write.csv(cluster_annotation, file.path(out_dir, "immune_cluster_annotation.csv"), row.names = FALSE, na = "NA")

immune_df <- data.frame(
  cell = colnames(obj)[immune_index], sample = obj$sample[immune_index], cluster = immune_cluster_vector,
  broad_cell_type = obj$mouse_broad_cell_type[immune_index], immune_cell_type = obj$mouse_immune_cell_type[immune_index],
  immune_subtype = obj$mouse_immune_subtype[immune_index], winning_signature = obj$mouse_immune_signature[immune_index],
  overlap = obj$mouse_immune_overlap[immune_index], f1 = obj$mouse_immune_f1[immune_index],
  overlap_genes = obj$mouse_immune_overlap_genes[immune_index], confidence = obj$mouse_immune_confidence[immune_index],
  stringsAsFactors = FALSE
)
utils::write.csv(immune_df, file.path(out_dir, "immune_annotation_per_cell.csv"), row.names = FALSE, na = "NA")

obj$mouse_immune_plot_group <- ifelse(
  obj$mouse_immune_candidate, paste0("Immune: ", obj$mouse_immune_cell_type), paste0("Non-immune: ", obj$mouse_broad_cell_type)
)
p <- Seurat::DimPlot(obj, reduction = "umap", group.by = "mouse_immune_plot_group", label = TRUE, repel = TRUE, raster = FALSE) +
  ggplot2::labs(title = "Mouse immune cell types: cluster-DEG signature overlap")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_immune_cell_type.pdf"), p, width = 12, height = 9)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_immune_cell_type.png"), p, width = 12, height = 9, dpi = 300)

obj@misc$mouse05_immune_annotation_contract <- list(
  annotation_level = "cluster", immune_cluster_rule = "05e broad cluster annotation equals Immune",
  marker_source = normalizePath(marker_file, mustWork = TRUE), top_positive_cluster_markers = top_marker_n,
  signature_source = "curated canonical mouse immune-lineage signatures embedded in 05f_mouse_immune_annotation.R",
  score_metrics = c("overlap", "precision", "recall", "f1", "jaccard"),
  label_ranking = "descending F1, overlap, Jaccard; zero-overlap clusters are Ambiguous",
  method_reference = "PANcanKFLs/code/packages/Annotation/GBM_annotation.R cluster-level annotation via marker overlaps",
  immune_clusters = immune_clusters, immune_cells = length(immune_index)
)
saveRDS(obj, file.path(out_dir, "mouse_broad_and_immune_annotated.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05f.txt"))
writeLines(
  c(
    "05f immune cluster-level annotation completed.", paste0("Immune clusters: ", length(immune_clusters)),
    paste0("Immune cells: ", length(immune_index)), paste0("Top positive DE markers used per cluster: up to ", top_marker_n),
    "No cell-level classifier was used."
  ),
  file.path(out_dir, "completion.txt")
)
message("[05f] Completed. Output: ", out_dir)
