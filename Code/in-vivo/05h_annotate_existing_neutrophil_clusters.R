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
required_packages <- c("yaml", "Seurat", "Matrix", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
input_rds <- file.path(results_root, "05g_final", "mouse_integrated_clustered_annotated_final.rds")
marker_file <- file.path(results_root, "05d_clustering", "cluster_positive_markers.csv")
out_dir <- file.path(results_root, "05h_existing_neutrophil_states")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(input_rds)) stop("Missing 05g input: ", input_rds, call. = FALSE)
if (!file.exists(marker_file)) stop("Missing 05d marker table: ", marker_file, call. = FALSE)

obj <- readRDS(input_rds)
if (!inherits(obj, "Seurat")) stop("05g input is not a Seurat object.", call. = FALSE)
required_meta <- c("sample", "seurat_clusters", "mouse_final_cell_type", "scDblFinder.class", "barcode_raw")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0L) stop("05g object lacks metadata: ", paste(missing_meta, collapse = ", "), call. = FALSE)
if (any(as.character(obj$scDblFinder.class) != "singlet")) stop("05h input contains non-singlets.", call. = FALSE)

neutrophil_label <- as.character(cfg$myeloid$existing_neutrophil_label)
neutrophil_cells <- colnames(obj)[as.character(obj$mouse_final_cell_type) == neutrophil_label]
if (length(neutrophil_cells) == 0L) stop("No existing Neutrophil cells were found in 05g.", call. = FALSE)
existing_clusters <- sort(unique(as.character(obj$seurat_clusters[neutrophil_cells])))
cluster_labels <- tapply(as.character(obj$mouse_final_cell_type), as.character(obj$seurat_clusters), function(x) unique(x))
if (any(vapply(cluster_labels[existing_clusters], length, integer(1)) != 1L)) {
  stop("05g final labels are not uniform within every existing Neutrophil cluster.", call. = FALSE)
}

markers <- utils::read.csv(marker_file, stringsAsFactors = FALSE, check.names = FALSE)
fc_column <- intersect(c("avg_log2FC", "avg_logFC"), colnames(markers))
if (!all(c("cluster", "gene") %in% colnames(markers)) || length(fc_column) != 1L) {
  stop("05d markers require cluster, gene, and exactly one avg_log2FC/avg_logFC column.", call. = FALSE)
}
markers$cluster <- as.character(markers$cluster)
markers$gene_upper <- toupper(trimws(as.character(markers$gene)))
markers$rank_fc <- as.numeric(markers[[fc_column]])
markers <- markers[nzchar(markers$gene_upper) & is.finite(markers$rank_fc), , drop = FALSE]
if (length(setdiff(existing_clusters, unique(markers$cluster))) > 0L) stop("At least one existing Neutrophil cluster lacks DEGs.", call. = FALSE)

state_signatures <- list(
  Hypoxic = c("Ero1a", "Bnip3", "Bnip3l", "P4ha1", "P4ha2", "Hmox1", "Ndrg1", "Slc2a1", "Vegfa", "Pdk1", "Adm", "Ldha", "Pgk1", "Eno1", "Aldoa", "Gpi", "Tpi1", "Anxa2"),
  `IFN-responsive` = c("Rsad2", "Ifit1", "Ifit2", "Ifit3", "Isg15", "Cxcl10", "Stat1", "Irf7", "Mx1", "Mx2", "Oas1a", "Oas2", "Ifitm3", "Gbp2", "Ifi47", "Usp18"),
  Inflammatory = c("Il1b", "Nlrp3", "Cxcl1", "Cxcl2", "Cxcl3", "Ccl3", "Ccl4", "Tnf", "Ptgs2", "Nfkbia", "Icam1", "Clec4d", "Clec4e", "Trem1", "Cd14", "Fpr1", "Fpr2", "Hcar2"),
  `Granule-high` = c("Ngp", "Ltf", "Camp", "Ly6g", "Cd177", "Mmp8", "Mmp9", "Retnlg", "Lcn2", "Pglyrp1", "Wfdc21", "Wfdc17", "Slpi", "S100a8", "S100a9", "Serpinb1a", "Stfa2l1", "Cstdc4"),
  `Immature-like` = c("Saa1", "Saa3", "Hp", "Chil1", "Chil3", "Lrg1", "Plac8", "Slfn4", "Vcan", "G0s2", "Sell", "Csf3r", "Lyz2", "Morrbid", "Retreg1", "Lsmem1", "Dock10", "P2rx7"),
  `CXCR4-high/aged-like` = c("Cxcr4", "Itgam", "Icam1", "Cd63", "Lgals3", "Apoe", "Ctsb", "Ctsd", "Fcer1g", "Tyrobp", "Fcgr3", "Mmp9", "Socs3", "Junb", "Fos", "Klf6")
)
expressed_upper <- unique(toupper(rownames(obj[["RNA"]])))
state_signatures <- lapply(state_signatures, function(x) intersect(unique(toupper(x)), expressed_upper))
if (any(vapply(state_signatures, length, integer(1)) == 0L)) stop("A Neutrophil-state signature has no genes in RNA.", call. = FALSE)

top_sizes <- sort(unique(as.integer(unlist(cfg$myeloid$existing_state_top_markers, use.names = FALSE))))
if (length(top_sizes) == 0L || any(!is.finite(top_sizes) | top_sizes < 1L)) stop("Invalid existing_state_top_markers.", call. = FALSE)
top_marker_rows <- list()
score_rows <- list()
row_i <- 1L
for (cluster_id in existing_clusters) {
  z <- markers[markers$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$rank_fc, z$gene_upper), , drop = FALSE]
  z <- z[!duplicated(z$gene_upper), , drop = FALSE]
  for (top_n in top_sizes) {
    local <- utils::head(z, top_n)
    marker_genes <- unique(local$gene_upper)
    top_marker_rows[[row_i]] <- data.frame(
      cluster = cluster_id, top_n = top_n, marker_rank = seq_len(nrow(local)),
      gene = local$gene, gene_upper = local$gene_upper, avg_log2FC = local$rank_fc,
      stringsAsFactors = FALSE
    )
    for (state in names(state_signatures)) {
      signature <- state_signatures[[state]]
      overlap_genes <- intersect(marker_genes, signature)
      overlap <- length(overlap_genes)
      precision <- overlap / length(marker_genes)
      recall <- overlap / length(signature)
      score_rows[[length(score_rows) + 1L]] <- data.frame(
        cluster = cluster_id, top_n = top_n, state = state,
        cluster_marker_n = length(marker_genes), signature_gene_n = length(signature), overlap = overlap,
        precision = precision, recall = recall,
        f1 = if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall),
        jaccard = overlap / length(unique(c(marker_genes, signature))),
        overlap_genes = paste(overlap_genes, collapse = ";"), stringsAsFactors = FALSE
      )
    }
    row_i <- row_i + 1L
  }
}
top_markers <- do.call(rbind, top_marker_rows)
scores <- do.call(rbind, score_rows)

aggregate_one <- function(cluster_id, state) {
  z <- scores[scores$cluster == cluster_id & scores$state == state, , drop = FALSE]
  genes <- unique(unlist(strsplit(z$overlap_genes[nzchar(z$overlap_genes)], ";", fixed = TRUE)))
  data.frame(
    cluster = cluster_id, state = state, max_overlap = max(z$overlap), mean_f1 = mean(z$f1),
    mean_jaccard = mean(z$jaccard), overlap_genes_union = paste(genes, collapse = ";"),
    stringsAsFactors = FALSE
  )
}
score_grid <- expand.grid(cluster = existing_clusters, state = names(state_signatures), stringsAsFactors = FALSE)
aggregate_scores <- do.call(rbind, Map(aggregate_one, score_grid$cluster, score_grid$state))
min_overlap <- as.integer(cfg$myeloid$min_signature_overlap)
margin_cutoff <- as.numeric(cfg$myeloid$mixed_f1_margin)
select_label <- function(cluster_id) {
  z <- aggregate_scores[aggregate_scores$cluster == cluster_id, , drop = FALSE]
  z <- z[order(-z$mean_f1, -z$max_overlap, -z$mean_jaccard, z$state), , drop = FALSE]
  best <- z[1L, , drop = FALSE]
  second <- z[2L, , drop = FALSE]
  label <- if (best$max_overlap[[1]] < min_overlap) "Unresolved" else best$state[[1]]
  margin <- best$mean_f1[[1]] - second$mean_f1[[1]]
  confidence <- if (label == "Unresolved") "unresolved" else if (best$max_overlap[[1]] >= 3L && margin > margin_cutoff) "high" else if (margin > 0) "medium" else "low"
  data.frame(
    cluster = cluster_id, mouse_neutrophil_state_existing = label,
    winning_state = best$state[[1]], max_overlap = best$max_overlap[[1]], mean_f1 = best$mean_f1[[1]],
    mean_jaccard = best$mean_jaccard[[1]], overlap_genes = best$overlap_genes_union[[1]],
    runner_up_state = second$state[[1]], runner_up_mean_f1 = second$mean_f1[[1]], f1_margin = margin,
    confidence = confidence, stringsAsFactors = FALSE
  )
}
cluster_annotation <- do.call(rbind, lapply(existing_clusters, select_label))
cluster_annotation$n_cells <- as.integer(table(as.character(obj$seurat_clusters[neutrophil_cells]))[cluster_annotation$cluster])

# Raw-expression lineage audit. This is validation only; state labels above are
# determined exclusively from cluster DEGs.
lineage_audit_signatures <- list(
  Neutrophil_core = c("S100a8", "S100a9", "Csf3r", "Retnlg", "Lcn2", "Cxcr2", "Fpr1", "Clec4d", "Clec4e", "Ly6g", "Ngp", "Camp"),
  Monocyte_core = c("Ly6c2", "Ccr2", "Vcan", "Sell", "Plac8", "Fn1", "Ms4a6c"),
  Macrophage_core = c("Adgre1", "Csf1r", "Mertk", "C1qa", "C1qb", "C1qc", "Mafb", "Apoe", "Mrc1")
)
Seurat::DefaultAssay(obj) <- "RNA"
expr <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
gene_lookup <- stats::setNames(rownames(expr), toupper(rownames(expr)))
audit_rows <- list()
for (cluster_id in existing_clusters) {
  cells <- neutrophil_cells[as.character(obj$seurat_clusters[neutrophil_cells]) == cluster_id]
  for (signature_name in names(lineage_audit_signatures)) {
    genes_upper <- intersect(toupper(lineage_audit_signatures[[signature_name]]), names(gene_lookup))
    genes <- unname(gene_lookup[genes_upper])
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      cluster = cluster_id, signature = signature_name, n_cells = length(cells), genes_detected = length(genes),
      mean_fraction_detected = mean(Matrix::rowMeans(expr[genes, cells, drop = FALSE] > 0)),
      mean_log_normalized_expression = mean(Matrix::rowMeans(expr[genes, cells, drop = FALSE])),
      stringsAsFactors = FALSE
    )
  }
}
lineage_audit <- do.call(rbind, audit_rows)

lookup <- function(column) stats::setNames(cluster_annotation[[column]], cluster_annotation$cluster)[as.character(obj$seurat_clusters)]
is_neutrophil <- colnames(obj) %in% neutrophil_cells
obj$mouse_neutrophil_state_existing <- NA_character_
obj$mouse_neutrophil_state_existing_confidence <- NA_character_
obj$mouse_neutrophil_state_existing_evidence <- NA_character_
obj$mouse_neutrophil_state_existing[is_neutrophil] <- unname(lookup("mouse_neutrophil_state_existing"))[is_neutrophil]
obj$mouse_neutrophil_state_existing_confidence[is_neutrophil] <- unname(lookup("confidence"))[is_neutrophil]
evidence_lookup <- stats::setNames(
  paste0("existing_cluster=", cluster_annotation$cluster, ";winning_state=", cluster_annotation$winning_state,
         ";max_overlap=", cluster_annotation$max_overlap, ";mean_f1=", signif(cluster_annotation$mean_f1, 5),
         ";runner_up=", cluster_annotation$runner_up_state, ";f1_margin=", signif(cluster_annotation$f1_margin, 5),
         ";overlap_genes=", cluster_annotation$overlap_genes),
  cluster_annotation$cluster
)
obj$mouse_neutrophil_state_existing_evidence[is_neutrophil] <- unname(evidence_lookup[as.character(obj$seurat_clusters[is_neutrophil])])
uniformity <- tapply(obj$mouse_neutrophil_state_existing[is_neutrophil], as.character(obj$seurat_clusters[is_neutrophil]), function(x) length(unique(x)))
if (any(uniformity != 1L)) stop("Existing-cluster state annotation is not cluster-uniform.", call. = FALSE)

state_library <- do.call(rbind, lapply(names(state_signatures), function(x) data.frame(state = x, gene_upper = state_signatures[[x]], stringsAsFactors = FALSE)))
utils::write.csv(state_library, file.path(out_dir, "neutrophil_state_signature_library.csv"), row.names = FALSE)
utils::write.csv(top_markers, file.path(out_dir, "existing_neutrophil_cluster_top_markers.csv"), row.names = FALSE)
utils::write.csv(scores, file.path(out_dir, "existing_neutrophil_state_scores_by_top_n.csv"), row.names = FALSE)
utils::write.csv(aggregate_scores, file.path(out_dir, "existing_neutrophil_state_scores_aggregated.csv"), row.names = FALSE)
utils::write.csv(cluster_annotation, file.path(out_dir, "existing_neutrophil_cluster_state_annotation.csv"), row.names = FALSE)
utils::write.csv(lineage_audit, file.path(out_dir, "existing_neutrophil_lineage_expression_audit.csv"), row.names = FALSE)

obj$mouse_neutrophil_state_existing_plot <- ifelse(is_neutrophil, obj$mouse_neutrophil_state_existing, "Other mouse cells")
state_colors <- unlist(cfg$myeloid$state_colors, use.names = TRUE)
plot_levels <- c(intersect(as.character(unlist(cfg$myeloid$state_order, use.names = FALSE)), unique(obj$mouse_neutrophil_state_existing[is_neutrophil])), "Other mouse cells")
obj$mouse_neutrophil_state_existing_plot <- factor(obj$mouse_neutrophil_state_existing_plot, levels = plot_levels)
plot_colors <- c(state_colors[setdiff(plot_levels, "Other mouse cells")], `Other mouse cells` = "#D9D9D9")
p <- Seurat::DimPlot(obj, reduction = "umap", group.by = "mouse_neutrophil_state_existing_plot", cols = unname(plot_colors), raster = FALSE) +
  ggplot2::labs(title = "Existing Neutrophil clusters: DEG-based state annotation")
ggplot2::ggsave(file.path(out_dir, "umap_existing_neutrophil_states.pdf"), p, width = 11, height = 8)
ggplot2::ggsave(file.path(out_dir, "umap_existing_neutrophil_states.png"), p, width = 11, height = 8, dpi = 300)

obj@misc$mouse05_existing_neutrophil_state_contract <- list(
  annotation_level = "existing_05d_cluster", state_label_source = "top positive cluster DEGs only",
  top_marker_sets = top_sizes, score_metrics = c("overlap", "precision", "recall", "f1", "jaccard"),
  raw_expression_role = "lineage validation only", existing_neutrophil_clusters = existing_clusters,
  existing_neutrophil_cells = length(neutrophil_cells)
)
saveRDS(obj, file.path(out_dir, "mouse_with_existing_neutrophil_states.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05h.txt"))
writeLines(c(
  "05h direct annotation of existing Neutrophil clusters completed.",
  paste0("Existing Neutrophil clusters: ", length(existing_clusters)),
  paste0("Existing Neutrophil cells: ", length(neutrophil_cells)),
  "State labels were derived from cluster DEGs; raw expression was used only for lineage audit."
), file.path(out_dir, "completion.txt"))
message("[05h] Completed. Output: ", out_dir)
