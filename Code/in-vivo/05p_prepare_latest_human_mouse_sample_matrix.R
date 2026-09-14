#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 260)
Sys.setenv(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1", NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

cluster_order <- function(x) {
  x <- unique(as.character(x))
  numeric_x <- suppressWarnings(as.numeric(x))
  if (all(is.finite(numeric_x))) x[order(numeric_x)] else sort(x)
}

read_mex_triplet <- function(matrix_dir) {
  files <- file.path(matrix_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
  if (!all(file.exists(files))) stop("Missing Cell Ranger MEX triplet in: ", matrix_dir, call. = FALSE)
  features <- utils::read.delim(
    gzfile(files[[2L]]), header = FALSE, sep = "\t", quote = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  barcodes <- scan(gzfile(files[[3L]]), what = character(), quiet = TRUE)
  counts <- methods::as(Matrix::readMM(gzfile(files[[1L]])), "CsparseMatrix")
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
    stop("MEX dimensions disagree in: ", matrix_dir, call. = FALSE)
  }
  list(counts = counts, features = features, barcodes = barcodes, files = files)
}

complete_category_counts <- function(cells, category_column, category_levels, sample_metadata, total_name) {
  counts <- stats::aggregate(
    list(n_cells = rep.int(1L, nrow(cells))),
    by = list(sample = cells$sample, category = as.character(cells[[category_column]])),
    FUN = sum
  )
  grid <- expand.grid(
    sample = sample_metadata$sample,
    category = category_levels,
    stringsAsFactors = FALSE
  )
  out <- merge(grid, counts, by = c("sample", "category"), all.x = TRUE, sort = FALSE)
  out$n_cells[is.na(out$n_cells)] <- 0L
  totals <- stats::aggregate(n_cells ~ sample, out, sum)
  names(totals)[2L] <- total_name
  out <- merge(out, totals, by = "sample", all.x = TRUE, sort = FALSE)
  meta_idx <- match(out$sample, sample_metadata$sample)
  out$Ploidy <- sample_metadata$Ploidy[meta_idx]
  out$Dose <- sample_metadata$Dose[meta_idx]
  out$sample_order <- sample_metadata$sample_order[meta_idx]
  out$category_order <- match(out$category, category_levels)
  names(out)[names(out) == "category"] <- category_column
  out[order(out$sample_order, out$category_order), , drop = FALSE]
}

write_csv_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "NA")
}

script_dir <- resolve_script_dir()
script_path <- file.path(script_dir, "05p_prepare_latest_human_mouse_sample_matrix.R")
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1L]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Matrix", "SeuratObject", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
input_root <- normalizePath(Sys.getenv("MOUSE05_INPUT_ROOT", unset = cfg$input_root), mustWork = TRUE)
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
mouse_rds <- file.path(results_root, "05h_existing_immune_cluster_annotation", "mouse_existing_immune_clusters_annotated.rds")
mouse_completion <- file.path(results_root, "05h_existing_immune_cluster_annotation", "completion.txt")
human_rdata <- normalizePath(as.character(cfg$unified_census$human_integrated_rdata), mustWork = TRUE)
metadata_source <- normalizePath(file.path(dirname(script_dir), "..", as.character(cfg$immune_composition$sample_metadata_source)), mustWork = TRUE)
out_dir <- file.path(results_root, "05p_latest_human_mouse_sample_matrix")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_inputs <- c(mouse_rds, mouse_completion, human_rdata, config_path, metadata_source)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) stop("Missing input(s): ", paste(missing_inputs, collapse = ", "), call. = FALSE)

metadata_sha256 <- digest::digest(metadata_source, algo = "sha256", file = TRUE, serialize = FALSE)
if (!identical(metadata_sha256, as.character(cfg$immune_composition$sample_metadata_source_sha256))) {
  stop("Authoritative sample metadata SHA256 mismatch.", call. = FALSE)
}

sample_sources <- unlist(cfg$sample_sources, use.names = TRUE)
sample_ids <- names(sample_sources)
excluded_samples <- as.character(unlist(cfg$immune_composition$excluded_samples, use.names = FALSE))
if (length(sample_ids) != 16L || any(sample_ids %in% excluded_samples)) {
  stop("05p requires exactly 16 in-vivo tissue samples and no cell-culture samples.", call. = FALSE)
}
metadata_cfg <- cfg$immune_composition$sample_metadata
sample_metadata <- do.call(rbind, lapply(names(metadata_cfg), function(sample_id) {
  z <- metadata_cfg[[sample_id]]
  data.frame(sample = sample_id, Ploidy = as.character(z$Ploidy), Dose = as.character(z$Dose), stringsAsFactors = FALSE)
}))
rownames(sample_metadata) <- NULL
if (!setequal(sample_metadata$sample, sample_ids) || anyDuplicated(sample_metadata$sample)) {
  stop("Configured sample metadata does not exactly match the 16 tissue samples.", call. = FALSE)
}
sample_metadata$sample_order <- match(sample_metadata$sample, sample_ids)
sample_metadata <- sample_metadata[order(sample_metadata$sample_order), , drop = FALSE]
sample_metadata$Ploidy <- factor(sample_metadata$Ploidy, levels = c("2N", "4N"))
sample_metadata$Dose <- factor(sample_metadata$Dose, levels = as.character(unlist(cfg$immune_composition$dose_order, use.names = FALSE)))
replicate_table <- table(sample_metadata$Dose, sample_metadata$Ploidy)
if (!identical(as.integer(replicate_table), c(4L, 2L, 2L, 4L, 2L, 2L))) {
  stop("Unexpected Dose-by-Ploidy replicate structure.", call. = FALSE)
}

message("[05p] Reading latest 05h mouse object.")
mouse <- readRDS(mouse_rds)
if (!inherits(mouse, "Seurat")) stop("05h input is not a Seurat object.", call. = FALSE)
required_mouse_meta <- c(
  "sample", "barcode_raw", "seurat_clusters", "mouse_fraction", "scDblFinder.class",
  "mouse_broad_cell_type", "mouse_parent_immune_review_candidate",
  "mouse_parent_immune_lineage", "mouse_parent_immune_state",
  "mouse_parent_immune_review_status", "mouse_parent_immune_annotation_level"
)
missing_mouse_meta <- setdiff(required_mouse_meta, colnames(mouse@meta.data))
if (length(missing_mouse_meta) > 0L) stop("05h mouse object lacks: ", paste(missing_mouse_meta, collapse = ", "), call. = FALSE)
human_feature_n <- sum(startsWith(rownames(mouse[["RNA"]]), as.character(cfg$species$human_prefix)))
if (human_feature_n != 0L) stop("Human features remain in the latest mouse RNA assay.", call. = FALSE)

mouse_meta <- mouse@meta.data
mouse_meta$cell <- colnames(mouse)
mouse_meta$sample <- as.character(mouse_meta$sample)
mouse_meta$barcode_raw <- as.character(mouse_meta$barcode_raw)
mouse_meta$seurat_clusters <- as.character(mouse_meta$seurat_clusters)
review_mask <- as.logical(mouse_meta$mouse_parent_immune_review_candidate)
if (anyNA(review_mask)) stop("05h immune-review mask contains NA.", call. = FALSE)
if (any(as.character(mouse_meta$scDblFinder.class) != "singlet")) stop("Latest mouse object contains non-singlets.", call. = FALSE)
if (any(as.numeric(mouse_meta$mouse_fraction) < as.numeric(cfg$species$mouse_fraction_min))) {
  stop("Latest mouse object violates the >=90% mouse-transcript rule.", call. = FALSE)
}
if (!setequal(unique(mouse_meta$sample), sample_ids) || any(mouse_meta$sample %in% excluded_samples)) {
  stop("Latest mouse object does not contain exactly the 16 tissue samples.", call. = FALSE)
}

mouse_meta$mouse_cell_type_latest <- as.character(mouse_meta$mouse_broad_cell_type)
immune_lineage <- as.character(mouse_meta$mouse_parent_immune_lineage)
immune_lineage[immune_lineage == "Ambiguous"] <- "Ambiguous mixed myeloid"
immune_lineage[immune_lineage == "Unresolved"] <- "Unresolved immune/myeloid"
mouse_meta$mouse_cell_type_latest[review_mask] <- immune_lineage[review_mask]

immune_state <- as.character(mouse_meta$mouse_parent_immune_state)
immune_state[is.na(immune_state)] <- ""
mouse_meta$mouse_cell_type_state_latest <- mouse_meta$mouse_cell_type_latest
specific_state <- review_mask & nzchar(immune_state) & immune_state != "Unresolved"
unresolved_state <- review_mask & immune_state == "Unresolved"
mouse_meta$mouse_cell_type_state_latest[specific_state] <- immune_state[specific_state]
mouse_meta$mouse_cell_type_state_latest[unresolved_state] <- paste0(
  mouse_meta$mouse_cell_type_latest[unresolved_state], " / Unresolved state"
)
mouse_meta$mouse_annotation_source <- ifelse(review_mask, "05h_existing_cluster_manual_annotation", "05f_broad_cluster_manual_annotation")

if (any(is.na(mouse_meta$mouse_cell_type_latest) | !nzchar(mouse_meta$mouse_cell_type_latest))) {
  stop("At least one mouse cell lacks a current primary cell type.", call. = FALSE)
}
if (any(is.na(mouse_meta$mouse_cell_type_state_latest) | !nzchar(mouse_meta$mouse_cell_type_state_latest))) {
  stop("At least one mouse cell lacks a current mutually exclusive fine category.", call. = FALSE)
}
type_uniformity <- tapply(mouse_meta$mouse_cell_type_latest, mouse_meta$seurat_clusters, function(x) length(unique(x)))
fine_uniformity <- tapply(mouse_meta$mouse_cell_type_state_latest, mouse_meta$seurat_clusters, function(x) length(unique(x)))
if (any(type_uniformity != 1L) || any(fine_uniformity != 1L)) {
  stop("Latest mouse labels are not uniform within the existing 05d clusters.", call. = FALSE)
}
if (!identical(mouse@misc$mouse05_existing_immune_annotation_contract$cell_level_classifier_used, FALSE)) {
  stop("05h annotation contract does not explicitly prohibit cell-level classification.", call. = FALSE)
}

mouse_cluster_order <- cluster_order(mouse_meta$seurat_clusters)
type_by_cluster <- vapply(mouse_cluster_order, function(k) unique(mouse_meta$mouse_cell_type_latest[mouse_meta$seurat_clusters == k]), character(1))
fine_by_cluster <- vapply(mouse_cluster_order, function(k) unique(mouse_meta$mouse_cell_type_state_latest[mouse_meta$seurat_clusters == k]), character(1))
mouse_type_levels <- unique(type_by_cluster)
mouse_fine_levels <- unique(fine_by_cluster)
mouse_key <- paste(mouse_meta$sample, mouse_meta$barcode_raw, sep = "|||")
if (anyDuplicated(mouse_key)) stop("Mouse sample-plus-barcode keys are duplicated.", call. = FALSE)
mouse_meta$cell_key <- mouse_key
mouse_cell_metadata <- mouse_meta[c(
  "cell", "cell_key", "sample", "barcode_raw", "seurat_clusters", "mouse_fraction", "scDblFinder.class",
  "mouse_broad_cell_type", "mouse_parent_immune_review_candidate", "mouse_parent_immune_lineage",
  "mouse_parent_immune_state", "mouse_cell_type_latest", "mouse_cell_type_state_latest", "mouse_annotation_source"
)]
rm(mouse)
invisible(gc())

message("[05p] Reading established human integrated object.")
human_env <- new.env(parent = emptyenv())
loaded <- load(human_rdata, envir = human_env)
human_object_name <- as.character(cfg$unified_census$human_object_name)
if (!human_object_name %in% loaded) stop("Configured human object is absent from RData: ", human_object_name, call. = FALSE)
human <- get(human_object_name, envir = human_env)
if (!inherits(human, "Seurat")) stop("Configured human object is not Seurat.", call. = FALSE)
sample_column <- as.character(cfg$unified_census$human_sample_column)
if (!sample_column %in% colnames(human@meta.data) || !"seurat_clusters" %in% colnames(human@meta.data)) {
  stop("Human object lacks sample or seurat_clusters metadata.", call. = FALSE)
}
human_sample <- as.character(human@meta.data[[sample_column]])
human_keep <- human_sample %in% sample_ids
human_raw_barcode <- sub("_[0-9]+$", "", colnames(human))
human_meta <- data.frame(
  cell = colnames(human)[human_keep],
  sample = human_sample[human_keep],
  barcode_raw = human_raw_barcode[human_keep],
  human_cluster = as.character(human$seurat_clusters[human_keep]),
  stringsAsFactors = FALSE
)
human_meta$cell_key <- paste(human_meta$sample, human_meta$barcode_raw, sep = "|||")
if (anyDuplicated(human_meta$cell_key)) stop("Human sample-plus-barcode keys are duplicated.", call. = FALSE)
if (!setequal(unique(human_meta$sample), sample_ids) || any(human_meta$sample %in% excluded_samples)) {
  stop("Established human object does not cover exactly the 16 tissue samples after restriction.", call. = FALSE)
}
human_cluster_levels <- cluster_order(human_meta$human_cluster)
rm(human, human_env)
invisible(gc())

human_prefix <- as.character(cfg$species$human_prefix)
mouse_prefix <- as.character(cfg$species$mouse_prefix)
expected_human <- as.integer(cfg$species$expected_human_features)
expected_mouse <- as.integer(cfg$species$expected_mouse_features)
species_threshold <- as.numeric(cfg$species$mouse_fraction_min)
if (!identical(species_threshold, as.numeric(cfg$unified_census$human_fraction_min))) {
  stop("Human and mouse species thresholds must be identical.", call. = FALSE)
}

raw_tables <- vector("list", length(sample_ids))
mex_manifest <- vector("list", length(sample_ids))
names(raw_tables) <- names(mex_manifest) <- sample_ids
for (sample_id in sample_ids) {
  matrix_dir <- normalizePath(file.path(input_root, sample_sources[[sample_id]]), mustWork = TRUE)
  message("[05p] Reading joint-reference Cell Ranger matrix: ", sample_id)
  dat <- read_mex_triplet(matrix_dir)
  feature_id <- as.character(dat$features[[1L]])
  is_human <- startsWith(feature_id, human_prefix)
  is_mouse <- startsWith(feature_id, mouse_prefix)
  if (sum(is_human) != expected_human || sum(is_mouse) != expected_mouse || sum(is_human | is_mouse) != nrow(dat$counts)) {
    stop("Joint-reference feature contract failed for ", sample_id, call. = FALSE)
  }
  human_umi <- as.numeric(Matrix::colSums(dat$counts[is_human, , drop = FALSE]))
  mouse_umi <- as.numeric(Matrix::colSums(dat$counts[is_mouse, , drop = FALSE]))
  species_total_umi <- human_umi + mouse_umi
  human_fraction <- ifelse(species_total_umi > 0, human_umi / species_total_umi, NA_real_)
  mouse_fraction <- ifelse(species_total_umi > 0, mouse_umi / species_total_umi, NA_real_)
  raw_species_class <- ifelse(
    species_total_umi == 0, "invalid_zero_species_umi",
    ifelse(human_fraction >= species_threshold, "tentative_human",
      ifelse(mouse_fraction >= species_threshold, "tentative_mouse", "mixed_cross_species_candidate")
    )
  )
  raw_tables[[sample_id]] <- data.frame(
    sample = sample_id,
    barcode_raw = dat$barcodes,
    cell_key = paste(sample_id, dat$barcodes, sep = "|||"),
    human_umi = human_umi,
    mouse_umi = mouse_umi,
    species_total_umi = species_total_umi,
    human_fraction = human_fraction,
    mouse_fraction = mouse_fraction,
    raw_species_class = raw_species_class,
    stringsAsFactors = FALSE
  )
  file_info <- file.info(dat$files)
  mex_manifest[[sample_id]] <- data.frame(
    sample = sample_id,
    role = c("matrix", "features", "barcodes"),
    path = normalizePath(dat$files, mustWork = TRUE),
    size_bytes = file_info$size,
    modified_time = format(file_info$mtime, "%Y-%m-%d %H:%M:%S %z"),
    stringsAsFactors = FALSE
  )
  rm(dat)
  invisible(gc())
}
raw <- do.call(rbind, raw_tables)
rownames(raw) <- NULL
if (anyDuplicated(raw$cell_key)) stop("Raw Cell Ranger sample-plus-barcode keys are duplicated.", call. = FALSE)
if (anyNA(match(mouse_meta$cell_key, raw$cell_key))) stop("At least one current mouse cell is absent from raw Cell Ranger matrices.", call. = FALSE)
if (anyNA(match(human_meta$cell_key, raw$cell_key))) stop("At least one retained human-object cell is absent from raw Cell Ranger matrices.", call. = FALSE)

human_idx <- match(raw$cell_key, human_meta$cell_key)
mouse_idx <- match(raw$cell_key, mouse_meta$cell_key)
raw$human_pipeline_retained <- !is.na(human_idx)
raw$mouse_pipeline_retained <- !is.na(mouse_idx)
if (any(raw$human_pipeline_retained & raw$mouse_pipeline_retained)) {
  stop("A raw barcode maps to both current human and mouse retained-cell universes.", call. = FALSE)
}
raw$census_status <- ifelse(
  raw$raw_species_class == "tentative_human" & raw$human_pipeline_retained, "retained_human_tumor",
  ifelse(raw$raw_species_class == "tentative_mouse" & raw$mouse_pipeline_retained, "retained_mouse_qc_singlet",
    ifelse(raw$raw_species_class == "tentative_human", "human_not_retained_by_established_pipeline",
      ifelse(raw$raw_species_class == "tentative_mouse", "mouse_qc_or_doublet_removed",
        ifelse(raw$raw_species_class == "mixed_cross_species_candidate", "mixed_cross_species_candidate", "invalid_zero_species_umi")
      )
    )
  )
)
human_species_mismatch <- raw$human_pipeline_retained & raw$raw_species_class != "tentative_human"
mouse_species_mismatch <- raw$mouse_pipeline_retained & raw$raw_species_class != "tentative_mouse"
raw$census_status[human_species_mismatch] <- "human_pipeline_cell_fails_90pct_human_rule"
raw$census_status[mouse_species_mismatch] <- "mouse_pipeline_cell_fails_90pct_mouse_rule"
raw$main_denominator_member <- raw$census_status %in% c("retained_human_tumor", "retained_mouse_qc_singlet")

retained_raw <- raw[raw$main_denominator_member, , drop = FALSE]
retained_human_idx <- match(retained_raw$cell_key, human_meta$cell_key)
retained_mouse_idx <- match(retained_raw$cell_key, mouse_meta$cell_key)
captured_species <- ifelse(retained_raw$census_status == "retained_human_tumor", "Human tumor", "Mouse host")
meta_idx <- match(retained_raw$sample, sample_metadata$sample)
unified <- data.frame(
  sample = retained_raw$sample,
  Ploidy = as.character(sample_metadata$Ploidy[meta_idx]),
  Dose = as.character(sample_metadata$Dose[meta_idx]),
  barcode_raw = retained_raw$barcode_raw,
  cell_key = retained_raw$cell_key,
  captured_species = captured_species,
  pipeline_cell = ifelse(
    captured_species == "Human tumor",
    human_meta$cell[retained_human_idx],
    mouse_meta$cell[retained_mouse_idx]
  ),
  human_cluster = human_meta$human_cluster[retained_human_idx],
  mouse_cluster = mouse_meta$seurat_clusters[retained_mouse_idx],
  mouse_cell_type_latest = mouse_meta$mouse_cell_type_latest[retained_mouse_idx],
  mouse_cell_type_state_latest = mouse_meta$mouse_cell_type_state_latest[retained_mouse_idx],
  human_umi = retained_raw$human_umi,
  mouse_umi = retained_raw$mouse_umi,
  human_fraction = retained_raw$human_fraction,
  mouse_fraction = retained_raw$mouse_fraction,
  stringsAsFactors = FALSE
)
if (nrow(unified) != sum(raw$main_denominator_member) || anyDuplicated(unified$cell_key)) {
  stop("Unified current human-mouse census construction failed.", call. = FALSE)
}

human_cells <- unified[unified$captured_species == "Human tumor", , drop = FALSE]
mouse_cells <- unified[unified$captured_species == "Mouse host", , drop = FALSE]
if (any(is.na(human_cells$human_cluster)) || any(is.na(mouse_cells$mouse_cell_type_latest))) {
  stop("A retained species-specific cell lacks its current cluster/type label.", call. = FALSE)
}
if (nrow(mouse_cells) != nrow(mouse_meta)) stop("Not every current 05h mouse cell entered the unified census.", call. = FALSE)

human_counts <- complete_category_counts(human_cells, "human_cluster", human_cluster_levels, sample_metadata, "sample_total_human_tumor_cells")
human_counts$fraction_within_human_tumor <- human_counts$n_cells / human_counts$sample_total_human_tumor_cells
mouse_type_counts <- complete_category_counts(mouse_cells, "mouse_cell_type_latest", mouse_type_levels, sample_metadata, "sample_total_mouse_cells")
mouse_type_counts$fraction_within_mouse <- mouse_type_counts$n_cells / mouse_type_counts$sample_total_mouse_cells
mouse_fine_counts <- complete_category_counts(mouse_cells, "mouse_cell_type_state_latest", mouse_fine_levels, sample_metadata, "sample_total_mouse_cells")
mouse_fine_counts$fraction_within_mouse <- mouse_fine_counts$n_cells / mouse_fine_counts$sample_total_mouse_cells

human_totals <- unique(human_counts[c("sample", "sample_total_human_tumor_cells")])
unified_totals <- as.data.frame(table(sample = unified$sample), stringsAsFactors = FALSE)
names(unified_totals)[2L] <- "sample_total_unified_cells"
for (table_name in c("human_counts", "mouse_type_counts", "mouse_fine_counts")) {
  x <- get(table_name)
  if (!"sample_total_human_tumor_cells" %in% names(x)) {
    x <- merge(x, human_totals, by = "sample", all.x = TRUE, sort = FALSE)
  }
  x <- merge(x, unified_totals, by = "sample", all.x = TRUE, sort = FALSE)
  x$fraction_of_unified_cells <- x$n_cells / x$sample_total_unified_cells
  if (startsWith(table_name, "mouse")) x$mouse_type_to_human_tumor_ratio <- x$n_cells / x$sample_total_human_tumor_cells
  x <- x[order(x$sample_order, x$category_order), , drop = FALSE]
  assign(table_name, x)
}

species_counts <- as.data.frame(table(sample = unified$sample, captured_species = unified$captured_species), stringsAsFactors = FALSE)
names(species_counts)[3L] <- "n_cells"
species_counts <- merge(species_counts, unified_totals, by = "sample", all.x = TRUE, sort = FALSE)
species_counts$fraction_of_unified_cells <- species_counts$n_cells / species_counts$sample_total_unified_cells
species_counts <- merge(species_counts, sample_metadata[c("sample", "Ploidy", "Dose", "sample_order")], by = "sample", all.x = TRUE, sort = FALSE)

status_counts <- as.data.frame(table(sample = raw$sample, census_status = raw$census_status), stringsAsFactors = FALSE)
names(status_counts)[3L] <- "n_barcodes"
status_counts <- status_counts[status_counts$n_barcodes > 0L, , drop = FALSE]
sample_summary <- merge(
  unique(human_counts[c("sample", "sample_total_human_tumor_cells")]),
  unique(mouse_type_counts[c("sample", "sample_total_mouse_cells")]),
  by = "sample", all = TRUE, sort = FALSE
)
sample_summary <- merge(sample_summary, unified_totals, by = "sample", all = TRUE, sort = FALSE)
sample_summary <- merge(sample_summary, sample_metadata[c("sample", "Ploidy", "Dose", "sample_order")], by = "sample", all.x = TRUE, sort = FALSE)
sample_summary$mouse_to_human_tumor_ratio <- sample_summary$sample_total_mouse_cells / sample_summary$sample_total_human_tumor_cells
sample_summary <- sample_summary[order(sample_summary$sample_order), , drop = FALSE]

mouse_cluster_annotation <- unique(mouse_meta[c(
  "seurat_clusters", "mouse_broad_cell_type", "mouse_parent_immune_review_candidate",
  "mouse_parent_immune_lineage", "mouse_parent_immune_state",
  "mouse_cell_type_latest", "mouse_cell_type_state_latest", "mouse_annotation_source"
)])
mouse_cluster_annotation$cluster_order <- match(mouse_cluster_annotation$seurat_clusters, mouse_cluster_order)
mouse_cluster_annotation <- mouse_cluster_annotation[order(mouse_cluster_annotation$cluster_order), , drop = FALSE]

if (any(abs(stats::aggregate(fraction_within_human_tumor ~ sample, human_counts, sum)$fraction_within_human_tumor - 1) > 1e-12)) {
  stop("Within-human cluster fractions do not sum to one for every sample.", call. = FALSE)
}
if (any(abs(stats::aggregate(fraction_within_mouse ~ sample, mouse_type_counts, sum)$fraction_within_mouse - 1) > 1e-12)) {
  stop("Primary within-mouse cell-type fractions do not sum to one for every sample.", call. = FALSE)
}
if (any(abs(stats::aggregate(fraction_within_mouse ~ sample, mouse_fine_counts, sum)$fraction_within_mouse - 1) > 1e-12)) {
  stop("Fine within-mouse category fractions do not sum to one for every sample.", call. = FALSE)
}

write_csv_gz(unified, file.path(out_dir, "latest_unified_human_mouse_cell_census.csv.gz"))
write_csv_gz(mouse_cell_metadata, file.path(out_dir, "latest_mouse_cell_metadata.csv.gz"))
utils::write.csv(sample_metadata, file.path(out_dir, "sample_ploidy_dose_metadata.csv"), row.names = FALSE)
utils::write.csv(sample_summary, file.path(out_dir, "sample_human_mouse_totals.csv"), row.names = FALSE)
utils::write.csv(human_counts, file.path(out_dir, "sample_human_tumor_cluster_counts_fractions.csv"), row.names = FALSE)
utils::write.csv(mouse_type_counts, file.path(out_dir, "sample_mouse_cell_type_counts_fractions.csv"), row.names = FALSE)
utils::write.csv(mouse_fine_counts, file.path(out_dir, "sample_mouse_fine_category_counts_fractions.csv"), row.names = FALSE)
utils::write.csv(species_counts, file.path(out_dir, "sample_human_mouse_species_composition.csv"), row.names = FALSE)
utils::write.csv(status_counts, file.path(out_dir, "cellranger_barcode_status_counts.csv"), row.names = FALSE)
utils::write.csv(mouse_cluster_annotation, file.path(out_dir, "latest_mouse_cluster_annotation_crosswalk.csv"), row.names = FALSE, na = "")
utils::write.csv(do.call(rbind, mex_manifest), file.path(out_dir, "joint_cellranger_input_manifest.csv"), row.names = FALSE)

analysis_object <- list(
  sample_metadata = sample_metadata,
  sample_summary = sample_summary,
  human_cluster_counts = human_counts,
  mouse_cell_type_counts = mouse_type_counts,
  mouse_fine_category_counts = mouse_fine_counts,
  species_counts = species_counts,
  mouse_cluster_annotation = mouse_cluster_annotation,
  contract = list(
    biological_unit = "sample",
    n_tissue_samples = length(sample_ids),
    excluded_cell_culture_samples = excluded_samples,
    human_cluster_source = "seurat_clusters in integrated_2025-04-17.RData, restricted to tissue samples",
    mouse_primary_label = "05f broad type for non-immune clusters and 05h manually approved immune lineage for reviewed clusters",
    mouse_fine_label = "05h manually approved immune state where available; otherwise current lineage/broad type",
    annotation_level = "existing 05d cluster only",
    cell_level_classifier_used = FALSE,
    species_threshold = species_threshold,
    joint_expression_integration_or_clustering = FALSE,
    main_interpretation = "sample-level scRNA-seq captured-cell composition, not absolute histologic abundance"
  )
)
saveRDS(analysis_object, file.path(out_dir, "latest_human_mouse_sample_matrix.rds"), compress = FALSE)

input_files <- c(mouse_rds, mouse_completion, human_rdata, normalizePath(config_path, mustWork = TRUE), script_path, metadata_source)
input_manifest <- data.frame(
  role = c(
    "latest_05h_cluster_level_mouse_annotation", "05h_completion_contract",
    "established_human_integrated_object", "analysis_config", "analysis_script", "authoritative_sample_metadata"
  ),
  path = normalizePath(input_files, mustWork = TRUE),
  size_bytes = file.info(input_files)$size,
  sha256 = vapply(input_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(input_manifest, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05p.txt"))
writeLines(c(
  "05p latest human-mouse sample matrix completed.",
  paste0("Tissue samples: ", length(sample_ids), "; 2N=", sum(sample_metadata$Ploidy == "2N"), "; 4N=", sum(sample_metadata$Ploidy == "4N"), "."),
  paste0("Current mouse cells: ", nrow(mouse_meta), "; current primary mouse types: ", length(mouse_type_levels), "; current mutually exclusive fine categories: ", length(mouse_fine_levels), "."),
  paste0("Retained human tumor cells passing the shared species rule: ", nrow(human_cells), "; human tumor clusters: ", length(human_cluster_levels), "."),
  paste0("Human pipeline cells failing the >=90% human rule: ", sum(human_species_mismatch), "."),
  paste0("Mouse pipeline cells failing the >=90% mouse rule: ", sum(mouse_species_mismatch), "."),
  "All current mouse labels are existing-cluster-level labels; no cell-level classifier was used.",
  "Human and mouse expression matrices were not jointly integrated or clustered.",
  "Fractions are sample-level scRNA-seq captured-cell fractions, not absolute histologic fractions."
), file.path(out_dir, "completion.txt"))

output_files <- setdiff(list.files(out_dir, full.names = TRUE), file.path(out_dir, "output_manifest_sha256.csv"))
output_manifest <- data.frame(
  file = basename(output_files),
  size_bytes = file.info(output_files)$size,
  sha256 = vapply(output_files, digest::digest, character(1), algo = "sha256", file = TRUE, serialize = FALSE),
  stringsAsFactors = FALSE
)
utils::write.csv(output_manifest, file.path(out_dir, "output_manifest_sha256.csv"), row.names = FALSE)
message("[05p] Completed. Output: ", out_dir)
