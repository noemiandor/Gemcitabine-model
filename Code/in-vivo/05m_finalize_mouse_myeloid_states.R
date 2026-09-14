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
required_packages <- c("yaml", "Seurat", "ggplot2", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
full_input <- file.path(results_root, "05h_existing_neutrophil_states", "mouse_with_existing_neutrophil_states.rds")
myeloid_input <- file.path(results_root, "05l_myeloid_annotation", "mouse_myeloid_reclustered_annotated.rds")
out_dir <- file.path(results_root, "05m_myeloid_states_final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(full_input)) stop("Missing 05h full-object input: ", full_input, call. = FALSE)
if (!file.exists(myeloid_input)) stop("Missing 05l myeloid input: ", myeloid_input, call. = FALSE)

full_obj <- readRDS(full_input)
myeloid <- readRDS(myeloid_input)
if (!inherits(full_obj, "Seurat") || !inherits(myeloid, "Seurat")) stop("05m inputs must both be Seurat objects.", call. = FALSE)
if (!all(colnames(myeloid) %in% colnames(full_obj))) stop("05l contains cells absent from the full 05h object.", call. = FALSE)
source_labels <- as.character(unlist(cfg$myeloid$source_final_labels, use.names = FALSE))
expected_myeloid <- colnames(full_obj)[as.character(full_obj$mouse_final_cell_type) %in% source_labels]
if (!setequal(colnames(myeloid), expected_myeloid)) stop("05l myeloid universe differs from configured 05h source labels.", call. = FALSE)

required_myeloid <- c(
  "mouse_myeloid_recluster", "mouse_myeloid_lineage_v2", "mouse_myeloid_lineage_v2_confidence",
  "mouse_myeloid_lineage_v2_evidence", "mouse_neutrophil_state_reclustered",
  "mouse_neutrophil_state_reclustered_confidence", "mouse_neutrophil_state_reclustered_evidence"
)
missing_myeloid <- setdiff(required_myeloid, colnames(myeloid@meta.data))
if (length(missing_myeloid) > 0L) stop("05l input lacks annotation fields: ", paste(missing_myeloid, collapse = ", "), call. = FALSE)

mapping_index <- match(colnames(full_obj), colnames(myeloid))
is_myeloid <- !is.na(mapping_index)
map_field <- function(field, mode = c("character", "numeric")) {
  mode <- match.arg(mode)
  out <- if (mode == "character") rep(NA_character_, ncol(full_obj)) else rep(NA_real_, ncol(full_obj))
  values <- myeloid@meta.data[[field]]
  if (is.factor(values)) values <- as.character(values)
  out[is_myeloid] <- values[mapping_index[is_myeloid]]
  out
}

full_obj$mouse_myeloid_recluster <- map_field("mouse_myeloid_recluster")
full_obj$mouse_myeloid_lineage_v2 <- map_field("mouse_myeloid_lineage_v2")
full_obj$mouse_myeloid_lineage_v2_confidence <- map_field("mouse_myeloid_lineage_v2_confidence")
full_obj$mouse_myeloid_lineage_v2_evidence <- map_field("mouse_myeloid_lineage_v2_evidence")
full_obj$mouse_neutrophil_state_final <- map_field("mouse_neutrophil_state_reclustered")
full_obj$mouse_neutrophil_state_confidence <- map_field("mouse_neutrophil_state_reclustered_confidence")
full_obj$mouse_neutrophil_state_evidence <- map_field("mouse_neutrophil_state_reclustered_evidence")
full_obj$mouse_neutrophil_state_final[full_obj$mouse_myeloid_lineage_v2 != "Neutrophil" | is.na(full_obj$mouse_myeloid_lineage_v2)] <- NA_character_
full_obj$mouse_neutrophil_state_confidence[is.na(full_obj$mouse_neutrophil_state_final)] <- NA_character_
full_obj$mouse_neutrophil_state_evidence[is.na(full_obj$mouse_neutrophil_state_final)] <- NA_character_

valid_lineages <- c("Neutrophil", "Monocyte", "Macrophage", "Mixed", "Unresolved")
bad_lineages <- setdiff(unique(stats::na.omit(full_obj$mouse_myeloid_lineage_v2)), valid_lineages)
if (length(bad_lineages) > 0L) stop("Unexpected 05l lineage labels: ", paste(bad_lineages, collapse = ", "), call. = FALSE)
configured_states <- as.character(unlist(cfg$myeloid$state_order, use.names = FALSE))
bad_states <- setdiff(unique(stats::na.omit(full_obj$mouse_neutrophil_state_final)), configured_states)
if (length(bad_states) > 0L) stop("Unexpected 05l Neutrophil states: ", paste(bad_states, collapse = ", "), call. = FALSE)

base_v2 <- as.character(full_obj$mouse_final_cell_type)
base_v2[is_myeloid] <- paste0("Immune: ", full_obj$mouse_myeloid_lineage_v2[is_myeloid])
base_v2[is_myeloid & (is.na(full_obj$mouse_myeloid_lineage_v2) | !nzchar(full_obj$mouse_myeloid_lineage_v2))] <- "Immune: Unresolved"
full_obj$mouse_final_cell_type_v2 <- base_v2

immune_prefix <- paste0("^", as.character(cfg$immune_composition$final_label_prefix), "\\s*")
immune_type <- rep(NA_character_, ncol(full_obj))
old_immune <- grepl(immune_prefix, as.character(full_obj$mouse_final_cell_type))
immune_type[old_immune] <- trimws(sub(immune_prefix, "", as.character(full_obj$mouse_final_cell_type[old_immune])))
immune_type[is_myeloid & full_obj$mouse_myeloid_lineage_v2 %in% c("Neutrophil", "Monocyte", "Macrophage")] <- full_obj$mouse_myeloid_lineage_v2[is_myeloid & full_obj$mouse_myeloid_lineage_v2 %in% c("Neutrophil", "Monocyte", "Macrophage")]
immune_type[is_myeloid & full_obj$mouse_myeloid_lineage_v2 %in% c("Mixed", "Unresolved")] <- "Other immune"
full_obj$mouse_final_immune_type <- immune_type

final_state_label <- base_v2
neutrophil_v2 <- is_myeloid & full_obj$mouse_myeloid_lineage_v2 == "Neutrophil"
if (any(neutrophil_v2 & is.na(full_obj$mouse_neutrophil_state_final))) stop("A v2 Neutrophil cell lacks a final state.", call. = FALSE)
final_state_label[neutrophil_v2] <- paste0("Immune: Neutrophil — ", full_obj$mouse_neutrophil_state_final[neutrophil_v2])
full_obj$mouse_final_cell_type_state <- final_state_label

comparison <- data.frame(
  cell = colnames(full_obj)[is_myeloid], sample = as.character(full_obj$sample[is_myeloid]),
  parent_cluster_05d = as.character(full_obj$seurat_clusters[is_myeloid]),
  source_final_cell_type = as.character(full_obj$mouse_final_cell_type[is_myeloid]),
  mouse_neutrophil_state_existing = as.character(full_obj$mouse_neutrophil_state_existing[is_myeloid]),
  mouse_myeloid_recluster = full_obj$mouse_myeloid_recluster[is_myeloid],
  mouse_myeloid_lineage_v2 = full_obj$mouse_myeloid_lineage_v2[is_myeloid],
  mouse_neutrophil_state_final = full_obj$mouse_neutrophil_state_final[is_myeloid],
  stringsAsFactors = FALSE
)
cluster_crosswalk <- as.data.frame(table(
  parent_cluster_05d = comparison$parent_cluster_05d,
  myeloid_recluster = comparison$mouse_myeloid_recluster,
  lineage_v2 = comparison$mouse_myeloid_lineage_v2
), stringsAsFactors = FALSE)
colnames(cluster_crosswalk)[4] <- "n_cells"
cluster_crosswalk <- cluster_crosswalk[cluster_crosswalk$n_cells > 0L, , drop = FALSE]
parent_totals <- tapply(cluster_crosswalk$n_cells, cluster_crosswalk$parent_cluster_05d, sum)
new_totals <- tapply(cluster_crosswalk$n_cells, cluster_crosswalk$myeloid_recluster, sum)
cluster_crosswalk$fraction_within_parent_cluster <- cluster_crosswalk$n_cells / unname(parent_totals[cluster_crosswalk$parent_cluster_05d])
cluster_crosswalk$fraction_within_myeloid_recluster <- cluster_crosswalk$n_cells / unname(new_totals[cluster_crosswalk$myeloid_recluster])

state_crosswalk <- as.data.frame(table(
  existing_state = comparison$mouse_neutrophil_state_existing,
  final_state = comparison$mouse_neutrophil_state_final,
  useNA = "ifany"
), stringsAsFactors = FALSE)
colnames(state_crosswalk)[3] <- "n_cells"
state_crosswalk <- state_crosswalk[state_crosswalk$n_cells > 0L, , drop = FALSE]
utils::write.csv(cluster_crosswalk, file.path(out_dir, "existing_cluster_to_myeloid_recluster_crosswalk.csv"), row.names = FALSE)
utils::write.csv(state_crosswalk, file.path(out_dir, "existing_vs_reclustered_neutrophil_state_crosswalk.csv"), row.names = FALSE, na = "NA")
utils::write.csv(comparison, file.path(out_dir, "myeloid_existing_vs_reclustered_annotation_per_cell.csv"), row.names = FALSE, na = "NA")

cell_df <- data.frame(
  cell = colnames(full_obj), sample = as.character(full_obj$sample), barcode_raw = as.character(full_obj$barcode_raw),
  mouse_fraction = as.numeric(full_obj$mouse_fraction), scDblFinder.class = as.character(full_obj$scDblFinder.class),
  parent_cluster_05d = as.character(full_obj$seurat_clusters),
  mouse_final_cell_type = as.character(full_obj$mouse_final_cell_type),
  mouse_neutrophil_state_existing = as.character(full_obj$mouse_neutrophil_state_existing),
  mouse_myeloid_recluster = full_obj$mouse_myeloid_recluster,
  mouse_myeloid_lineage_v2 = full_obj$mouse_myeloid_lineage_v2,
  mouse_myeloid_lineage_v2_confidence = full_obj$mouse_myeloid_lineage_v2_confidence,
  mouse_myeloid_lineage_v2_evidence = full_obj$mouse_myeloid_lineage_v2_evidence,
  mouse_neutrophil_state_final = full_obj$mouse_neutrophil_state_final,
  mouse_neutrophil_state_confidence = full_obj$mouse_neutrophil_state_confidence,
  mouse_neutrophil_state_evidence = full_obj$mouse_neutrophil_state_evidence,
  mouse_final_immune_type = full_obj$mouse_final_immune_type,
  mouse_final_cell_type_v2 = full_obj$mouse_final_cell_type_v2,
  mouse_final_cell_type_state = full_obj$mouse_final_cell_type_state,
  stringsAsFactors = FALSE
)
utils::write.csv(cell_df, file.path(out_dir, "mouse_cell_annotations_with_myeloid_states.csv"), row.names = FALSE, na = "NA")

state_colors <- unlist(cfg$myeloid$state_colors, use.names = TRUE)
immune_colors <- unlist(cfg$immune_composition$subtype_colors, use.names = TRUE)
nonimmune_colors <- unlist(cfg$annotation$final_nonimmune_colors, use.names = TRUE)
lineage_colors <- unlist(cfg$myeloid$lineage_colors, use.names = TRUE)
labels <- sort(unique(final_state_label))
palette <- stats::setNames(rep(NA_character_, length(labels)), labels)
for (label in labels) {
  if (startsWith(label, "Immune: Neutrophil — ")) {
    state <- sub("^Immune: Neutrophil — ", "", label)
    palette[[label]] <- state_colors[[state]]
  } else if (startsWith(label, "Immune: ")) {
    subtype <- sub("^Immune: ", "", label)
    if (subtype %in% names(immune_colors)) palette[[label]] <- immune_colors[[subtype]]
    if (subtype %in% names(lineage_colors)) palette[[label]] <- lineage_colors[[subtype]]
  } else if (label %in% names(nonimmune_colors)) {
    palette[[label]] <- nonimmune_colors[[label]]
  }
}
if (anyNA(palette)) stop("Final state UMAP palette lacks labels: ", paste(names(palette)[is.na(palette)], collapse = ", "), call. = FALSE)
full_obj$mouse_final_cell_type_state <- factor(final_state_label, levels = labels)
palette_df <- data.frame(mouse_final_cell_type_state = labels, color_hex = unname(palette), stringsAsFactors = FALSE)
utils::write.csv(palette_df, file.path(out_dir, "mouse_final_cell_type_state_palette.csv"), row.names = FALSE)
p <- Seurat::DimPlot(
  full_obj, reduction = "umap", group.by = "mouse_final_cell_type_state",
  cols = unname(palette), label = TRUE, repel = TRUE, raster = FALSE
) + ggplot2::labs(title = "Final mouse cell types and Neutrophil states")
ggplot2::ggsave(file.path(out_dir, "umap_mouse_final_cell_type_state.pdf"), p, width = 13, height = 9)
ggplot2::ggsave(file.path(out_dir, "umap_mouse_final_cell_type_state.png"), p, width = 13, height = 9, dpi = 300)

full_obj@misc$mouse05_myeloid_state_final_contract <- list(
  existing_state_field = "mouse_neutrophil_state_existing",
  authoritative_state_field = "mouse_neutrophil_state_final",
  authoritative_state_source = "05l independently reclustered myeloid compartment cluster DEGs",
  myeloid_lineage_field = "mouse_myeloid_lineage_v2", final_display_field = "mouse_final_cell_type_state",
  human_mouse_joint_clustering_used = FALSE, final_state_palette = as.list(palette)
)
final_rds <- file.path(out_dir, "mouse_integrated_annotated_myeloid_states_final.rds")
tmp_rds <- file.path(out_dir, paste0(".mouse_integrated_annotated_myeloid_states_final.rds.tmp.", Sys.getpid()))
on.exit(unlink(tmp_rds, force = TRUE), add = TRUE)
saveRDS(full_obj, tmp_rds, compress = FALSE)
if (!file.rename(tmp_rds, final_rds)) stop("Failed to atomically replace 05m final RDS.", call. = FALSE)

output_files <- c(
  final_rds, file.path(out_dir, "mouse_cell_annotations_with_myeloid_states.csv"),
  file.path(out_dir, "existing_cluster_to_myeloid_recluster_crosswalk.csv"),
  file.path(out_dir, "existing_vs_reclustered_neutrophil_state_crosswalk.csv"),
  file.path(out_dir, "mouse_final_cell_type_state_palette.csv"),
  file.path(out_dir, "umap_mouse_final_cell_type_state.pdf"),
  file.path(out_dir, "umap_mouse_final_cell_type_state.png")
)
fi <- file.info(output_files)
manifest <- data.frame(
  path = normalizePath(output_files, mustWork = TRUE), size_bytes = fi$size,
  sha256 = vapply(output_files, function(x) digest::digest(x, algo = "sha256", file = TRUE, serialize = FALSE), character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(manifest, file.path(out_dir, "05m_output_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05m.txt"))
writeLines(c(
  "05m final mouse myeloid/Neutrophil state mapping completed.",
  paste0("Full mouse singlets: ", ncol(full_obj)), paste0("Mapped myeloid cells: ", sum(is_myeloid)),
  paste0("Final v2 Neutrophil cells: ", sum(neutrophil_v2)),
  "mouse_neutrophil_state_final is authoritative; existing-cluster state is retained for comparison."
), file.path(out_dir, "completion.txt"))
message("[05m] Completed. Output: ", out_dir)
