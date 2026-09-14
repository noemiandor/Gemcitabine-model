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

read_mex_triplet <- function(matrix_dir) {
  files <- file.path(matrix_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
  if (!all(file.exists(files))) stop("Missing Cell Ranger MEX triplet in: ", matrix_dir, call. = FALSE)
  features <- utils::read.delim(gzfile(files[[2]]), header = FALSE, sep = "\t", quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  barcodes <- scan(gzfile(files[[3]]), what = character(), quiet = TRUE)
  counts <- methods::as(Matrix::readMM(gzfile(files[[1]])), "CsparseMatrix")
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) stop("MEX dimensions disagree in: ", matrix_dir, call. = FALSE)
  list(counts = counts, features = features, barcodes = barcodes, files = files)
}

script_dir <- resolve_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")
required_packages <- c("yaml", "Matrix", "SeuratObject")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
input_root <- normalizePath(Sys.getenv("MOUSE05_INPUT_ROOT", unset = cfg$input_root), mustWork = TRUE)
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
mouse_annotation_file <- file.path(results_root, "05m_myeloid_states_final", "mouse_cell_annotations_with_myeloid_states.csv")
human_rdata <- normalizePath(as.character(cfg$unified_census$human_integrated_rdata), mustWork = TRUE)
out_dir <- file.path(results_root, "05n_unified_human_mouse_census")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(mouse_annotation_file)) stop("Missing 05m annotation table: ", mouse_annotation_file, call. = FALSE)

sample_sources <- unlist(cfg$sample_sources, use.names = TRUE)
sample_ids <- names(sample_sources)
excluded <- as.character(unlist(cfg$immune_composition$excluded_samples, use.names = FALSE))
if (length(sample_sources) != 16L || any(excluded %in% sample_ids)) stop("05n requires the frozen 16 tissue-sample contract.", call. = FALSE)
human_prefix <- as.character(cfg$species$human_prefix)
mouse_prefix <- as.character(cfg$species$mouse_prefix)
expected_human <- as.integer(cfg$species$expected_human_features)
expected_mouse <- as.integer(cfg$species$expected_mouse_features)
human_fraction_min <- as.numeric(cfg$unified_census$human_fraction_min)
mouse_fraction_min <- as.numeric(cfg$species$mouse_fraction_min)
if (human_fraction_min != mouse_fraction_min) stop("05n requires identical human and mouse species thresholds.", call. = FALSE)

raw_tables <- list()
manifest <- list()
for (sample_id in sample_ids) {
  matrix_dir <- normalizePath(file.path(input_root, sample_sources[[sample_id]]), mustWork = TRUE)
  message("[05n] Reading joint-reference Cell Ranger matrix: ", sample_id)
  dat <- read_mex_triplet(matrix_dir)
  feature_id <- as.character(dat$features[[1]])
  is_human <- startsWith(feature_id, human_prefix)
  is_mouse <- startsWith(feature_id, mouse_prefix)
  if (sum(is_human) != expected_human || sum(is_mouse) != expected_mouse || sum(is_human | is_mouse) != nrow(dat$counts)) {
    stop("Joint-reference feature contract failed for ", sample_id, call. = FALSE)
  }
  human_umi <- as.numeric(Matrix::colSums(dat$counts[is_human, , drop = FALSE]))
  mouse_umi <- as.numeric(Matrix::colSums(dat$counts[is_mouse, , drop = FALSE]))
  total <- human_umi + mouse_umi
  human_fraction <- ifelse(total > 0, human_umi / total, NA_real_)
  mouse_fraction <- ifelse(total > 0, mouse_umi / total, NA_real_)
  species_class <- ifelse(
    total == 0, "invalid_zero_species_umi",
    ifelse(human_fraction >= human_fraction_min, "tentative_human",
      ifelse(mouse_fraction >= mouse_fraction_min, "tentative_mouse", "mixed_cross_species_candidate")
    )
  )
  raw_tables[[sample_id]] <- data.frame(
    sample = sample_id, barcode_raw = dat$barcodes, cell_key = paste(sample_id, dat$barcodes, sep = "|||"),
    human_umi = human_umi, mouse_umi = mouse_umi, species_total_umi = total,
    human_fraction = human_fraction, mouse_fraction = mouse_fraction, raw_species_class = species_class,
    stringsAsFactors = FALSE
  )
  fi <- file.info(dat$files)
  manifest[[sample_id]] <- data.frame(
    sample = sample_id, role = c("matrix", "features", "barcodes"), path = normalizePath(dat$files, mustWork = TRUE),
    size_bytes = fi$size, modified_time = format(fi$mtime, "%Y-%m-%d %H:%M:%S %z"), stringsAsFactors = FALSE
  )
  rm(dat)
  invisible(gc())
}
raw <- do.call(rbind, raw_tables)
rownames(raw) <- NULL
if (anyDuplicated(raw$cell_key)) stop("Raw sample plus barcode keys are not unique.", call. = FALSE)

# The previously established human object supplies the retained human-tumor
# universe. It is linked by sample and raw barcode; its expression matrix is not
# combined with the mouse expression matrix and is never used for joint PCA/CCA.
human_env <- new.env(parent = emptyenv())
loaded <- load(human_rdata, envir = human_env)
human_object_name <- as.character(cfg$unified_census$human_object_name)
if (!human_object_name %in% loaded) stop("Configured human object is absent from RData: ", human_object_name, call. = FALSE)
human <- get(human_object_name, envir = human_env)
if (!inherits(human, "Seurat")) stop("Configured human object is not Seurat.", call. = FALSE)
sample_column <- as.character(cfg$unified_census$human_sample_column)
if (!sample_column %in% colnames(human@meta.data)) stop("Human object lacks sample column: ", sample_column, call. = FALSE)
human_sample <- as.character(human@meta.data[[sample_column]])
if (any(excluded %in% sample_ids)) stop("Excluded cell-line samples leaked into the configured tissue sample list.", call. = FALSE)
human_keep <- human_sample %in% sample_ids
human_raw_barcode <- sub("_[0-9]+$", "", colnames(human))
human_key <- paste(human_sample, human_raw_barcode, sep = "|||")
human_meta <- data.frame(
  cell_key = human_key[human_keep], human_pipeline_cell = colnames(human)[human_keep],
  sample = human_sample[human_keep], barcode_raw = human_raw_barcode[human_keep],
  human_pipeline_retained = TRUE,
  human_cluster = if ("seurat_clusters" %in% colnames(human@meta.data)) as.character(human$seurat_clusters[human_keep]) else NA_character_,
  human_sample_type = if ("sample_type" %in% colnames(human@meta.data)) as.character(human$sample_type[human_keep]) else NA_character_,
  human_blueprint_main = if ("blueprint.main" %in% colnames(human@meta.data)) as.character(human@meta.data$blueprint.main[human_keep]) else NA_character_,
  human_hpca_main = if ("hpca.main" %in% colnames(human@meta.data)) as.character(human@meta.data$hpca.main[human_keep]) else NA_character_,
  stringsAsFactors = FALSE
)
if (anyDuplicated(human_meta$cell_key)) stop("Human sample plus raw-barcode keys are not unique.", call. = FALSE)
if (!setequal(unique(human_meta$sample), sample_ids)) stop("Established human object does not cover all 16 tissue samples.", call. = FALSE)

mouse_meta <- utils::read.csv(mouse_annotation_file, stringsAsFactors = FALSE, check.names = FALSE)
required_mouse <- c(
  "cell", "sample", "barcode_raw", "mouse_fraction", "scDblFinder.class", "mouse_final_cell_type_v2",
  "mouse_final_immune_type", "mouse_neutrophil_state_final", "mouse_final_cell_type_state"
)
missing_mouse <- setdiff(required_mouse, colnames(mouse_meta))
if (length(missing_mouse) > 0L) stop("05m annotation table lacks fields: ", paste(missing_mouse, collapse = ", "), call. = FALSE)
mouse_meta$cell_key <- paste(mouse_meta$sample, mouse_meta$barcode_raw, sep = "|||")
if (anyDuplicated(mouse_meta$cell_key)) stop("Mouse sample plus raw-barcode keys are not unique.", call. = FALSE)
if (any(mouse_meta$scDblFinder.class != "singlet") || any(mouse_meta$mouse_fraction < mouse_fraction_min)) stop("05m mouse table violates QC/singlet/species contracts.", call. = FALSE)

human_idx <- match(raw$cell_key, human_meta$cell_key)
mouse_idx <- match(raw$cell_key, mouse_meta$cell_key)
raw$human_pipeline_retained <- !is.na(human_idx)
raw$mouse_pipeline_retained <- !is.na(mouse_idx)
if (any(raw$human_pipeline_retained & raw$mouse_pipeline_retained)) stop("A raw barcode maps to both retained human and retained mouse cells.", call. = FALSE)

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
raw$main_denominator_member <- raw$census_status %in% c("retained_human_tumor", "retained_mouse_qc_singlet")

# A retained pipeline cell that fails the common 90-percent raw species rule is
# explicitly excluded from the main denominator and reported as a mismatch.
human_species_mismatch <- raw$human_pipeline_retained & raw$raw_species_class != "tentative_human"
mouse_species_mismatch <- raw$mouse_pipeline_retained & raw$raw_species_class != "tentative_mouse"
raw$census_status[human_species_mismatch] <- "human_pipeline_cell_fails_90pct_human_rule"
raw$census_status[mouse_species_mismatch] <- "mouse_pipeline_cell_fails_90pct_mouse_rule"
raw$main_denominator_member[human_species_mismatch | mouse_species_mismatch] <- FALSE

retained_raw <- raw[raw$main_denominator_member, , drop = FALSE]
retained_human_idx <- match(retained_raw$cell_key, human_meta$cell_key)
retained_mouse_idx <- match(retained_raw$cell_key, mouse_meta$cell_key)
retained_species <- ifelse(retained_raw$census_status == "retained_human_tumor", "Human tumor", "Mouse host")
unified <- data.frame(
  sample = retained_raw$sample, barcode_raw = retained_raw$barcode_raw, cell_key = retained_raw$cell_key,
  captured_species = retained_species, human_umi = retained_raw$human_umi, mouse_umi = retained_raw$mouse_umi,
  human_fraction = retained_raw$human_fraction, mouse_fraction = retained_raw$mouse_fraction,
  pipeline_cell = ifelse(
    retained_species == "Human tumor", human_meta$human_pipeline_cell[retained_human_idx], mouse_meta$cell[retained_mouse_idx]
  ),
  human_cluster = human_meta$human_cluster[retained_human_idx],
  human_sample_type = human_meta$human_sample_type[retained_human_idx],
  human_blueprint_main = human_meta$human_blueprint_main[retained_human_idx],
  human_hpca_main = human_meta$human_hpca_main[retained_human_idx],
  mouse_final_cell_type_v2 = mouse_meta$mouse_final_cell_type_v2[retained_mouse_idx],
  mouse_final_immune_type = mouse_meta$mouse_final_immune_type[retained_mouse_idx],
  mouse_neutrophil_state_final = mouse_meta$mouse_neutrophil_state_final[retained_mouse_idx],
  mouse_final_cell_type_state = mouse_meta$mouse_final_cell_type_state[retained_mouse_idx],
  stringsAsFactors = FALSE
)
if (nrow(unified) != sum(raw$main_denominator_member) || anyDuplicated(unified$cell_key)) stop("Unified main denominator construction failed.", call. = FALSE)

status_counts <- as.data.frame(table(sample = raw$sample, census_status = raw$census_status), stringsAsFactors = FALSE)
colnames(status_counts)[3] <- "n_cells"
status_counts <- status_counts[status_counts$n_cells > 0L, , drop = FALSE]
species_counts <- as.data.frame(table(sample = unified$sample, captured_species = unified$captured_species), stringsAsFactors = FALSE)
colnames(species_counts)[3] <- "n_cells"
species_totals <- tapply(species_counts$n_cells, species_counts$sample, sum)
species_counts$fraction_of_human_plus_mouse_retained_singlets <- species_counts$n_cells / unname(species_totals[species_counts$sample])

sample_summary <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
for (status in sort(unique(raw$census_status))) {
  tab <- stats::setNames(status_counts$n_cells[status_counts$census_status == status], status_counts$sample[status_counts$census_status == status])
  sample_summary[[status]] <- as.integer(ifelse(sample_summary$sample %in% names(tab), tab[sample_summary$sample], 0L))
}
sample_summary$main_denominator_cells <- as.integer(species_totals[sample_summary$sample])
sample_summary$main_denominator_definition <- "retained human tumor cells from established human object plus 05m mouse QC singlets; both pass >=90% raw species fraction"

utils::write.csv(raw, file.path(out_dir, "all_filtered_barcodes_species_and_pipeline_status.csv"), row.names = FALSE, na = "NA")
utils::write.csv(unified, file.path(out_dir, "unified_human_mouse_captured_cell_census.csv"), row.names = FALSE, na = "NA")
utils::write.csv(status_counts, file.path(out_dir, "census_status_counts_by_sample.csv"), row.names = FALSE)
utils::write.csv(species_counts, file.path(out_dir, "human_vs_mouse_retained_composition_by_sample.csv"), row.names = FALSE)
utils::write.csv(sample_summary, file.path(out_dir, "unified_census_sample_summary.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, manifest), file.path(out_dir, "joint_cellranger_input_manifest.csv"), row.names = FALSE)
utils::write.csv(human_meta, file.path(out_dir, "established_human_pipeline_cell_crosswalk.csv"), row.names = FALSE, na = "NA")

human_info <- file.info(human_rdata)
provenance <- data.frame(
  role = c("established human retained-cell universe", "05m final mouse retained-cell universe"),
  path = c(human_rdata, normalizePath(mouse_annotation_file, mustWork = TRUE)),
  size_bytes = c(human_info$size, file.info(mouse_annotation_file)$size),
  modified_time = format(c(human_info$mtime, file.info(mouse_annotation_file)$mtime), "%Y-%m-%d %H:%M:%S %z"),
  stringsAsFactors = FALSE
)
utils::write.csv(provenance, file.path(out_dir, "unified_census_pipeline_input_provenance.csv"), row.names = FALSE)

contract <- list(
  cell_ranger_reference = list(human_features = expected_human, mouse_features = expected_mouse),
  species_threshold = human_fraction_min, species_denominator = "human UMI plus mouse UMI within each raw filtered barcode",
  human_universe = "membership in integrated_2025-04-17.RData established human object, restricted to 16 tissue samples",
  mouse_universe = "05m cells retained after 05b QC/scDblFinder and subsequent annotation",
  main_denominator = "retained human tumor plus retained mouse host cells that also pass their >=90% species threshold",
  mixed_and_pipeline_filtered_cells = "excluded from main denominator but fully reported",
  joint_expression_clustering = FALSE,
  interpretation = "scRNA-seq captured-cell fraction; not an absolute histologic tissue fraction"
)
saveRDS(list(census = unified, sample_summary = sample_summary, contract = contract), file.path(out_dir, "unified_human_mouse_census.rds"), compress = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05n.txt"))
writeLines(c(
  "05n unified human-mouse captured-cell census completed.",
  paste0("Raw filtered barcodes classified: ", nrow(raw)),
  paste0("Retained human tumor cells in main denominator: ", sum(unified$captured_species == "Human tumor")),
  paste0("Retained mouse host cells in main denominator: ", sum(unified$captured_species == "Mouse host")),
  paste0("Human pipeline cells failing >=90% human rule: ", sum(human_species_mismatch)),
  paste0("Mouse pipeline cells failing >=90% mouse rule: ", sum(mouse_species_mismatch)),
  "Human and mouse expression matrices were not jointly integrated or clustered.",
  "Reported proportions are captured-cell fractions, not absolute tissue fractions."
), file.path(out_dir, "completion.txt"))
message("[05n] Completed. Output: ", out_dir)
