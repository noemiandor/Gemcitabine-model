#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 240)
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
input_rds <- file.path(results_root, "05k_myeloid_cluster_degs", "mouse_myeloid_reclustered_with_degs.rds")
marker_file <- file.path(results_root, "05k_myeloid_cluster_degs", "myeloid_cluster_positive_degs_full.csv")
out_dir <- file.path(results_root, "05l_myeloid_annotation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05k input: ", input_rds, call. = FALSE)
if (!file.exists(marker_file)) stop("Missing 05k DEG table: ", marker_file, call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05k input is not a Seurat object.", call. = FALSE)
if (!"seurat_clusters" %in% colnames(obj@meta.data)) stop("05k input lacks seurat_clusters.", call. = FALSE)
markers <- utils::read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC", "rank_fc"), colnames(markers))
if (!all(c("cluster", "gene") %in% colnames(markers)) || length(fc_column) < 1L) stop("05k DEG table lacks cluster/gene/fold-change columns.", call. = FALSE)
fc_column <- fc_column[[1L]]
markers$cluster <- as.character(markers$cluster)
markers$gene_upper <- toupper(trimws(as.character(markers$gene)))
markers$rank_fc <- as.numeric(markers[[fc_column]])
markers <- markers[nzchar(markers$gene_upper) & is.finite(markers$rank_fc), , drop = FALSE]
cluster_ids <- sort(unique(as.character(obj$seurat_clusters)))
if (length(setdiff(cluster_ids, unique(markers$cluster))) > 0L) stop("At least one reclustered myeloid cluster lacks DEGs.", call. = FALSE)

lineage_signatures <- list(
  Neutrophil = c("Csf3r", "Ly6g", "Retnlg", "Ngp", "Camp", "Ltf", "Mmp8", "Mmp9", "Cd177", "Cxcr2", "Fpr1", "Fpr2", "S100a8", "S100a9", "Lcn2", "Pglyrp1", "Wfdc21", "Slpi", "Trem1", "Clec4d", "Clec4e"),
  Monocyte = c("Ly6c2", "Ccr2", "Vcan", "Sell", "Plac8", "Fn1", "Ms4a6c", "Lyz2", "Ctss", "S100a10", "Fcgr1", "Lst1", "Fcer1g", "Ctsb", "Cst3", "Ifitm3", "Lgals3", "Ccr1"),
  Macrophage = c("Adgre1", "Csf1r", "Mertk", "Cd68", "Mafb", "Apoe", "C1qa", "C1qb", "C1qc", "Ms4a7", "Trem2", "Marco", "Mrc1", "Cd163", "Folr2", "Ctsd", "Lpl", "Acp5", "Aif1", "Lgmn")
)
state_signatures <- list(
  Hypoxic = c("Ero1a", "Bnip3", "Bnip3l", "P4ha1", "P4ha2", "Hmox1", "Ndrg1", "Slc2a1", "Vegfa", "Pdk1", "Adm", "Ldha", "Pgk1", "Eno1", "Aldoa", "Gpi", "Tpi1", "Anxa2"),
  `IFN-responsive` = c("Rsad2", "Ifit1", "Ifit2", "Ifit3", "Isg15", "Cxcl10", "Stat1", "Irf7", "Mx1", "Mx2", "Oas1a", "Oas2", "Ifitm3", "Gbp2", "Ifi47", "Usp18"),
  Inflammatory = c("Il1b", "Nlrp3", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl3", "Ccl4", "Tnf", "Ptgs2", "Nfkbia", "Icam1", "Clec4d", "Clec4e", "Trem1", "Cd14", "Fpr1", "Fpr2", "Hcar2"),
  `Granule-high` = c("Ngp", "Ltf", "Camp", "Ly6g", "Cd177", "Mmp8", "Mmp9", "Retnlg", "Lcn2", "Pglyrp1", "Wfdc21", "Wfdc17", "Slpi", "S100a8", "S100a9", "Serpinb1a", "Stfa2l1", "Cstdc4"),
  `Immature-like` = c("Saa1", "Saa3", "Hp", "Chil1", "Chil3", "Lrg1", "Plac8", "Slfn4", "Vcan", "G0s2", "Sell", "Csf3r", "Lyz2", "Morrbid", "Retreg1", "Lsmem1", "Dock10", "P2rx7"),
  `CXCR4-high/aged-like` = c("Cxcr4", "Itgam", "Icam1", "Cd63", "Lgals3", "Apoe", "Ctsb", "Ctsd", "Fcer1g", "Tyrobp", "Fcgr3", "Mmp9", "Socs3", "Junb", "Fos", "Klf6")
)
expressed_upper <- unique(toupper(rownames(obj[["RNA"]])))
lineage_signatures <- lapply(lineage_signatures, function(x) intersect(unique(toupper(x)), expressed_upper))
state_signatures <- lapply(state_signatures, function(x) intersect(unique(toupper(x)), expressed_upper))
if (any(vapply(c(lineage_signatures, state_signatures), length, integer(1)) == 0L)) stop("At least one lineage/state signature has no expressed genes.", call. = FALSE)

ranked <- do.call(rbind, lapply(cluster_ids, function(cluster_id) {
  z <- markers[markers$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$rank_fc, z$gene_upper), , drop = FALSE]
  z <- z[!duplicated(z$gene_upper), , drop = FALSE]
  z$marker_rank <- seq_len(nrow(z))
  z
}))
rownames(ranked) <- NULL

score_signatures <- function(signatures, top_sizes, score_class) {
  rows <- list()
  for (cluster_id in cluster_ids) {
    for (top_n in top_sizes) {
      marker_genes <- unique(ranked$gene_upper[ranked$cluster == cluster_id & ranked$marker_rank <= top_n])
      for (signature_name in names(signatures)) {
        signature <- signatures[[signature_name]]
        overlap_genes <- intersect(marker_genes, signature)
        overlap <- length(overlap_genes)
        precision <- overlap / length(marker_genes)
        recall <- overlap / length(signature)
        rows[[length(rows) + 1L]] <- data.frame(
          score_class = score_class, cluster = cluster_id, top_n = top_n, signature = signature_name,
          cluster_marker_n = length(marker_genes), signature_gene_n = length(signature), overlap = overlap,
          precision = precision, recall = recall,
          f1 = if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall),
          jaccard = overlap / length(unique(c(marker_genes, signature))),
          overlap_genes = paste(overlap_genes, collapse = ";"), stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

lineage_top_n <- as.integer(cfg$myeloid$lineage_top_markers)
state_top_sizes <- sort(unique(as.integer(unlist(cfg$myeloid$state_top_markers, use.names = FALSE))))
lineage_scores <- score_signatures(lineage_signatures, lineage_top_n, "lineage")
state_scores <- score_signatures(state_signatures, state_top_sizes, "state")

aggregate_scores <- function(scores) {
  keys <- unique(scores[c("cluster", "signature")])
  do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    z <- scores[scores$cluster == keys$cluster[[i]] & scores$signature == keys$signature[[i]], , drop = FALSE]
    genes <- unique(unlist(strsplit(z$overlap_genes[nzchar(z$overlap_genes)], ";", fixed = TRUE)))
    data.frame(
      cluster = keys$cluster[[i]], signature = keys$signature[[i]], max_overlap = max(z$overlap),
      mean_f1 = mean(z$f1), mean_jaccard = mean(z$jaccard),
      overlap_genes_union = paste(genes, collapse = ";"), stringsAsFactors = FALSE
    )
  }))
}
lineage_aggregate <- aggregate_scores(lineage_scores)
state_aggregate <- aggregate_scores(state_scores)
min_overlap <- as.integer(cfg$myeloid$min_signature_overlap)
margin_cutoff <- as.numeric(cfg$myeloid$mixed_f1_margin)

select_best <- function(aggregate_table, cluster_id, allow_mixed = TRUE) {
  z <- aggregate_table[aggregate_table$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$mean_f1, -z$max_overlap, -z$mean_jaccard, z$signature), , drop = FALSE]
  best <- z[1L, , drop = FALSE]
  second <- z[2L, , drop = FALSE]
  margin <- best$mean_f1[[1]] - second$mean_f1[[1]]
  label <- if (best$max_overlap[[1]] < min_overlap) {
    "Unresolved"
  } else if (allow_mixed && second$max_overlap[[1]] >= min_overlap && margin <= margin_cutoff) {
    "Mixed"
  } else {
    best$signature[[1]]
  }
  confidence <- if (label == "Unresolved") "unresolved" else if (label == "Mixed") "mixed" else if (best$max_overlap[[1]] >= 3L && margin > margin_cutoff) "high" else "medium"
  data.frame(
    cluster = cluster_id, label = label, winning_signature = best$signature[[1]],
    max_overlap = best$max_overlap[[1]], mean_f1 = best$mean_f1[[1]], mean_jaccard = best$mean_jaccard[[1]],
    overlap_genes = best$overlap_genes_union[[1]], runner_up_signature = second$signature[[1]],
    runner_up_mean_f1 = second$mean_f1[[1]], f1_margin = margin, confidence = confidence,
    stringsAsFactors = FALSE
  )
}

lineage_deg_annotation <- do.call(rbind, lapply(cluster_ids, function(x) select_best(lineage_aggregate, x, allow_mixed = TRUE)))
colnames(lineage_deg_annotation)[colnames(lineage_deg_annotation) == "label"] <- "mouse_myeloid_lineage_deg_candidate"
colnames(lineage_deg_annotation)[colnames(lineage_deg_annotation) == "confidence"] <- "deg_lineage_confidence"

# Within a compartment that was deliberately assembled from already annotated
# Neutrophil/Monocyte/Macrophage cells, state-driven DEGs can hide shared lineage
# genes (for example, an IFN-high Neutrophil cluster). Preserve lineage from the
# dominant source lineage at the *reclustered cluster level*. The DEG-overlap
# lineage call remains in the output as an explicit audit, while cluster DEGs are
# the authoritative evidence for the downstream Neutrophil state label.
source_lineage_min <- as.numeric(cfg$myeloid$lineage_source_dominance_min)
if (!is.finite(source_lineage_min) || source_lineage_min <= 0.5 || source_lineage_min > 1) {
  stop("myeloid.lineage_source_dominance_min must be in (0.5, 1].", call. = FALSE)
}
source_counts <- as.data.frame(table(
  cluster = as.character(obj$seurat_clusters),
  source_lineage = as.character(obj$mouse_myeloid_source_lineage)
), stringsAsFactors = FALSE)
colnames(source_counts)[3L] <- "n_cells"
source_counts <- source_counts[source_counts$n_cells > 0L, , drop = FALSE]
source_summary <- do.call(rbind, lapply(cluster_ids, function(cluster_id) {
  z <- source_counts[source_counts$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$n_cells, z$source_lineage), , drop = FALSE]
  total <- sum(z$n_cells)
  dominant_fraction <- z$n_cells[[1L]] / total
  final_lineage <- if (dominant_fraction >= source_lineage_min) z$source_lineage[[1L]] else "Mixed"
  data.frame(
    cluster = cluster_id,
    mouse_myeloid_lineage_v2 = final_lineage,
    confidence = if (dominant_fraction >= 0.90) "high" else if (dominant_fraction >= source_lineage_min) "medium" else "mixed",
    source_dominant_lineage = z$source_lineage[[1L]],
    source_dominant_fraction = dominant_fraction,
    source_lineage_counts = paste0(z$source_lineage, ":", z$n_cells, collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
lineage_annotation <- merge(source_summary, lineage_deg_annotation, by = "cluster", all.x = TRUE, sort = FALSE)
lineage_annotation <- lineage_annotation[match(cluster_ids, lineage_annotation$cluster), , drop = FALSE]
neutrophil_clusters <- lineage_annotation$cluster[lineage_annotation$mouse_myeloid_lineage_v2 == "Neutrophil"]
state_annotation <- if (length(neutrophil_clusters) > 0L) {
  do.call(rbind, lapply(neutrophil_clusters, function(x) select_best(state_aggregate, x, allow_mixed = TRUE)))
} else {
  data.frame()
}
if (nrow(state_annotation) > 0L) colnames(state_annotation)[colnames(state_annotation) == "label"] <- "mouse_neutrophil_state_reclustered"

cluster_annotation <- lineage_annotation
cluster_annotation$mouse_neutrophil_state_reclustered <- NA_character_
cluster_annotation$state_winning_signature <- NA_character_
cluster_annotation$state_max_overlap <- NA_integer_
cluster_annotation$state_mean_f1 <- NA_real_
cluster_annotation$state_overlap_genes <- NA_character_
cluster_annotation$state_runner_up_signature <- NA_character_
cluster_annotation$state_f1_margin <- NA_real_
cluster_annotation$state_confidence <- NA_character_
if (nrow(state_annotation) > 0L) {
  idx <- match(state_annotation$cluster, cluster_annotation$cluster)
  cluster_annotation$mouse_neutrophil_state_reclustered[idx] <- state_annotation$mouse_neutrophil_state_reclustered
  cluster_annotation$state_winning_signature[idx] <- state_annotation$winning_signature
  cluster_annotation$state_max_overlap[idx] <- state_annotation$max_overlap
  cluster_annotation$state_mean_f1[idx] <- state_annotation$mean_f1
  cluster_annotation$state_overlap_genes[idx] <- state_annotation$overlap_genes
  cluster_annotation$state_runner_up_signature[idx] <- state_annotation$runner_up_signature
  cluster_annotation$state_f1_margin[idx] <- state_annotation$f1_margin
  cluster_annotation$state_confidence[idx] <- state_annotation$confidence
}
cluster_annotation$n_cells <- as.integer(table(as.character(obj$seurat_clusters))[cluster_annotation$cluster])

lookup <- function(column) stats::setNames(cluster_annotation[[column]], cluster_annotation$cluster)[as.character(obj$seurat_clusters)]
obj$mouse_myeloid_recluster <- as.character(obj$seurat_clusters)
obj$mouse_myeloid_lineage_v2 <- unname(lookup("mouse_myeloid_lineage_v2"))
obj$mouse_myeloid_lineage_v2_confidence <- unname(lookup("confidence"))
obj$mouse_myeloid_lineage_v2_evidence <- paste0(
  "cluster=", obj$mouse_myeloid_recluster,
  ";source_dominant_lineage=", unname(lookup("source_dominant_lineage")),
  ";source_dominant_fraction=", signif(as.numeric(unname(lookup("source_dominant_fraction"))), 5),
  ";source_lineage_counts=", unname(lookup("source_lineage_counts")),
  ";deg_lineage_candidate=", unname(lookup("mouse_myeloid_lineage_deg_candidate")),
  ";deg_winning_signature=", unname(lookup("winning_signature")),
  ";max_overlap=", unname(lookup("max_overlap")),
  ";mean_f1=", signif(as.numeric(unname(lookup("mean_f1"))), 5),
  ";runner_up=", unname(lookup("runner_up_signature")),
  ";f1_margin=", signif(as.numeric(unname(lookup("f1_margin"))), 5),
  ";overlap_genes=", unname(lookup("overlap_genes"))
)
obj$mouse_neutrophil_state_reclustered <- unname(lookup("mouse_neutrophil_state_reclustered"))
obj$mouse_neutrophil_state_reclustered_confidence <- unname(lookup("state_confidence"))
obj$mouse_neutrophil_state_reclustered_evidence <- ifelse(
  obj$mouse_myeloid_lineage_v2 == "Neutrophil",
  paste0(
    "cluster=", obj$mouse_myeloid_recluster,
    ";winning_state=", unname(lookup("state_winning_signature")),
    ";max_overlap=", unname(lookup("state_max_overlap")),
    ";mean_f1=", signif(as.numeric(unname(lookup("state_mean_f1"))), 5),
    ";runner_up=", unname(lookup("state_runner_up_signature")),
    ";f1_margin=", signif(as.numeric(unname(lookup("state_f1_margin"))), 5),
    ";overlap_genes=", unname(lookup("state_overlap_genes"))
  ), NA_character_
)

lineage_uniformity <- tapply(obj$mouse_myeloid_lineage_v2, obj$mouse_myeloid_recluster, function(x) length(unique(x)))
if (any(lineage_uniformity != 1L)) stop("Reclustered myeloid lineage labels are not cluster-uniform.", call. = FALSE)
state_uniformity <- tapply(obj$mouse_neutrophil_state_reclustered[obj$mouse_myeloid_lineage_v2 == "Neutrophil"], obj$mouse_myeloid_recluster[obj$mouse_myeloid_lineage_v2 == "Neutrophil"], function(x) length(unique(x)))
if (length(state_uniformity) > 0L && any(state_uniformity != 1L)) stop("Reclustered Neutrophil state labels are not cluster-uniform.", call. = FALSE)

signature_library <- rbind(
  do.call(rbind, lapply(names(lineage_signatures), function(x) data.frame(score_class = "lineage", signature = x, gene_upper = lineage_signatures[[x]], stringsAsFactors = FALSE))),
  do.call(rbind, lapply(names(state_signatures), function(x) data.frame(score_class = "state", signature = x, gene_upper = state_signatures[[x]], stringsAsFactors = FALSE)))
)
utils::write.csv(signature_library, file.path(out_dir, "myeloid_lineage_and_neutrophil_state_signature_library.csv"), row.names = FALSE)
utils::write.csv(lineage_scores, file.path(out_dir, "myeloid_lineage_scores.csv"), row.names = FALSE)
utils::write.csv(lineage_aggregate, file.path(out_dir, "myeloid_lineage_scores_aggregated.csv"), row.names = FALSE)
utils::write.csv(source_counts, file.path(out_dir, "myeloid_cluster_source_lineage_counts.csv"), row.names = FALSE)
utils::write.csv(state_scores, file.path(out_dir, "neutrophil_state_scores.csv"), row.names = FALSE)
utils::write.csv(state_aggregate, file.path(out_dir, "neutrophil_state_scores_aggregated.csv"), row.names = FALSE)
utils::write.csv(cluster_annotation, file.path(out_dir, "myeloid_cluster_lineage_state_annotation.csv"), row.names = FALSE, na = "NA")

cell_annotation <- data.frame(
  cell = colnames(obj), sample = as.character(obj$sample), barcode_raw = as.character(obj$barcode_raw),
  parent_cluster_05d = as.character(obj$mouse_parent_cluster_05d), source_lineage = as.character(obj$mouse_myeloid_source_lineage),
  mouse_myeloid_recluster = obj$mouse_myeloid_recluster, mouse_myeloid_lineage_v2 = obj$mouse_myeloid_lineage_v2,
  mouse_myeloid_lineage_v2_confidence = obj$mouse_myeloid_lineage_v2_confidence,
  mouse_myeloid_lineage_v2_evidence = obj$mouse_myeloid_lineage_v2_evidence,
  mouse_neutrophil_state_existing = as.character(obj$mouse_neutrophil_state_existing),
  mouse_neutrophil_state_reclustered = obj$mouse_neutrophil_state_reclustered,
  mouse_neutrophil_state_reclustered_confidence = obj$mouse_neutrophil_state_reclustered_confidence,
  mouse_neutrophil_state_reclustered_evidence = obj$mouse_neutrophil_state_reclustered_evidence,
  stringsAsFactors = FALSE
)
utils::write.csv(cell_annotation, file.path(out_dir, "myeloid_cell_lineage_state_annotation.csv"), row.names = FALSE, na = "NA")

lineage_colors <- unlist(cfg$myeloid$lineage_colors, use.names = TRUE)
lineage_levels <- intersect(names(lineage_colors), unique(obj$mouse_myeloid_lineage_v2))
obj$mouse_myeloid_lineage_v2 <- factor(obj$mouse_myeloid_lineage_v2, levels = lineage_levels)
p_lineage <- Seurat::DimPlot(obj, reduction = "umap", group.by = "mouse_myeloid_lineage_v2", cols = unname(lineage_colors[lineage_levels]), label = TRUE, repel = TRUE, raster = FALSE) +
  ggplot2::labs(title = "Reclustered mouse myeloid lineages", subtitle = "Cluster-level dominant source lineage; DEG-overlap lineage candidate retained as an audit")

state_colors <- unlist(cfg$myeloid$state_colors, use.names = TRUE)
obj$mouse_myeloid_state_plot <- ifelse(as.character(obj$mouse_myeloid_lineage_v2) == "Neutrophil", obj$mouse_neutrophil_state_reclustered, as.character(obj$mouse_myeloid_lineage_v2))
state_plot_levels <- c(
  intersect(as.character(unlist(cfg$myeloid$state_order, use.names = FALSE)), unique(obj$mouse_myeloid_state_plot)),
  intersect(c("Monocyte", "Macrophage", "Mixed", "Unresolved"), unique(obj$mouse_myeloid_state_plot))
)
state_plot_levels <- unique(state_plot_levels)
obj$mouse_myeloid_state_plot <- factor(obj$mouse_myeloid_state_plot, levels = state_plot_levels)
state_plot_colors <- c(state_colors, lineage_colors)[state_plot_levels]
p_state <- Seurat::DimPlot(obj, reduction = "umap", group.by = "mouse_myeloid_state_plot", cols = unname(state_plot_colors), label = TRUE, repel = TRUE, raster = FALSE) +
  ggplot2::labs(title = "Reclustered mouse myeloid lineages and Neutrophil states")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_lineage_v2.pdf"), p_lineage, width = 11, height = 8)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_myeloid_lineage_v2.png"), p_lineage, width = 11, height = 8, dpi = 300)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_neutrophil_states_reclustered.pdf"), p_state, width = 12, height = 9)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_neutrophil_states_reclustered.png"), p_state, width = 12, height = 9, dpi = 300)

obj@misc$mouse05_myeloid_annotation_contract <- list(
  annotation_level = "reclustered_myeloid_cluster", lineage_first = TRUE,
  lineage_source = "dominant existing lineage within each reclustered cluster",
  lineage_source_dominance_min = source_lineage_min,
  lineage_deg_overlap_role = "audit only because state-driven within-myeloid DEGs can suppress shared lineage markers",
  state_source = "top positive cluster DEGs vs Neutrophil state signatures, evaluated only after Neutrophil lineage assignment",
  state_genes_used_for_lineage = FALSE, score_metrics = c("overlap", "precision", "recall", "f1", "jaccard"),
  min_signature_overlap = min_overlap, mixed_f1_margin = margin_cutoff
)
saveRDS(obj, file.path(out_dir, "mouse_myeloid_reclustered_annotated.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05l.txt"))
writeLines(c(
  "05l reclustered mouse myeloid annotation completed.", paste0("Cells: ", ncol(obj)),
  paste0("Clusters: ", length(cluster_ids)), paste0("Neutrophil clusters: ", length(neutrophil_clusters)),
  "Lineage was fixed per reclustered cluster from its dominant pre-existing lineage; DEG-overlap lineage calls were retained as an audit.",
  "Neutrophil state was assigned from reclustered-cluster DEGs only within lineage-confirmed Neutrophil clusters."
), file.path(out_dir, "completion.txt"))
message("[05l] Completed. Output: ", out_dir)
