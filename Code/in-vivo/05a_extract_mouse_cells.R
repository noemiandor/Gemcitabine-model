#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, width = 220)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

resolve_script_dir <- function() {
  hit <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(hit) != 1L) stop("Cannot resolve script directory.", call. = FALSE)
  dirname(normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE))
}

script_dir <- resolve_script_dir()
args <- commandArgs(trailingOnly = TRUE)
config_path <- if (length(args) >= 1L) args[[1]] else file.path(script_dir, "05_mouse_cell_analysis_config.yaml")

required_packages <- c("yaml", "Matrix", "Seurat", "digest")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

cfg <- yaml::read_yaml(normalizePath(config_path, mustWork = TRUE))
input_root <- normalizePath(Sys.getenv("MOUSE05_INPUT_ROOT", unset = cfg$input_root), mustWork = TRUE)
results_root <- Sys.getenv("MOUSE05_RESULTS_ROOT", unset = cfg$results_root)
out_dir <- file.path(results_root, "05a_species_extraction")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

human_prefix <- as.character(cfg$species$human_prefix)
mouse_prefix <- as.character(cfg$species$mouse_prefix)
mouse_fraction_min <- as.numeric(cfg$species$mouse_fraction_min)
expected_human <- as.integer(cfg$species$expected_human_features)
expected_mouse <- as.integer(cfg$species$expected_mouse_features)
sample_sources <- unlist(cfg$sample_sources, use.names = TRUE)
excluded_cell_lines <- c("2N-Cell-Culture", "4N-Cell-Culture")

if (
  length(sample_sources) != 16L ||
  anyDuplicated(names(sample_sources)) ||
  any(excluded_cell_lines %in% names(sample_sources))
) {
  stop(
    "The frozen input contract must contain exactly 16 uniquely named xenograft-tissue samples and exclude both cell-line samples.",
    call. = FALSE
  )
}
if (!is.finite(mouse_fraction_min) || mouse_fraction_min < 0 || mouse_fraction_min > 1) {
  stop("species.mouse_fraction_min must be between 0 and 1.", call. = FALSE)
}

read_mex_triplet <- function(matrix_dir) {
  files <- file.path(matrix_dir, c("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"))
  if (!all(file.exists(files))) {
    stop("Missing Cell Ranger MEX triplet in: ", matrix_dir, call. = FALSE)
  }
  features <- utils::read.delim(
    gzfile(files[[2]]), header = FALSE, sep = "\t", quote = "",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (ncol(features) < 2L) stop("features.tsv.gz must contain at least two columns: ", files[[2]], call. = FALSE)
  barcodes <- scan(gzfile(files[[3]]), what = character(), quiet = TRUE)
  counts <- Matrix::readMM(gzfile(files[[1]]))
  counts <- methods::as(counts, "CsparseMatrix")
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
    stop("Matrix/features/barcodes dimensions do not agree in: ", matrix_dir, call. = FALSE)
  }
  list(counts = counts, features = features, barcodes = barcodes, files = files)
}

sha256 <- function(path) digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE)

mouse_objects <- list()
classification_tables <- list()
sample_summaries <- list()
input_manifest <- list()
feature_mapping_reference <- NULL

for (sample_id in names(sample_sources)) {
  matrix_dir <- normalizePath(file.path(input_root, sample_sources[[sample_id]]), mustWork = TRUE)
  message("[05a] Reading ", sample_id, " from ", matrix_dir)
  dat <- read_mex_triplet(matrix_dir)
  feature_id <- as.character(dat$features[[1]])
  feature_name <- as.character(dat$features[[2]])
  is_human <- startsWith(feature_id, human_prefix)
  is_mouse <- startsWith(feature_id, mouse_prefix)

  if (sum(is_human) != expected_human || sum(is_mouse) != expected_mouse) {
    stop(
      "Unexpected feature contract for ", sample_id,
      ": GRCh38=", sum(is_human), ", GRCm39=", sum(is_mouse),
      "; expected ", expected_human, " and ", expected_mouse, call. = FALSE
    )
  }
  if (any(is_human & is_mouse) || sum(is_human | is_mouse) != nrow(dat$counts)) {
    stop("Every feature must belong uniquely to GRCh38 or GRCm39 for sample: ", sample_id, call. = FALSE)
  }

  human_umi <- as.numeric(Matrix::colSums(dat$counts[is_human, , drop = FALSE]))
  mouse_umi <- as.numeric(Matrix::colSums(dat$counts[is_mouse, , drop = FALSE]))
  species_total_umi <- human_umi + mouse_umi
  mouse_fraction <- ifelse(species_total_umi > 0, mouse_umi / species_total_umi, NA_real_)
  tentative_mouse <- !is.na(mouse_fraction) & mouse_fraction >= mouse_fraction_min
  cell_name <- paste0(sample_id, "_", dat$barcodes)

  class_df <- data.frame(
    cell = cell_name,
    sample = sample_id,
    barcode_raw = dat$barcodes,
    human_umi = human_umi,
    mouse_umi = mouse_umi,
    species_total_umi = species_total_umi,
    mouse_fraction = mouse_fraction,
    mouse_percent = 100 * mouse_fraction,
    species_class = ifelse(
      species_total_umi == 0,
      "invalid_zero_species_umi",
      ifelse(tentative_mouse, "tentative_mouse", "not_mouse_below_90pct")
    ),
    stringsAsFactors = FALSE
  )
  classification_tables[[sample_id]] <- class_df

  sample_summaries[[sample_id]] <- data.frame(
    sample = sample_id,
    source_matrix_dir = matrix_dir,
    raw_filtered_barcodes = length(dat$barcodes),
    tentative_mouse_cells = sum(tentative_mouse),
    tentative_mouse_percent = 100 * mean(tentative_mouse),
    human_features = sum(is_human),
    mouse_features = sum(is_mouse),
    pooled_mouse_fraction = sum(mouse_umi) / sum(species_total_umi),
    stringsAsFactors = FALSE
  )

  fi <- file.info(dat$files)
  input_manifest[[sample_id]] <- data.frame(
    sample = sample_id,
    role = c("matrix", "features", "barcodes"),
    path = normalizePath(dat$files, mustWork = TRUE),
    size_bytes = fi$size,
    modified_time = format(fi$mtime, "%Y-%m-%d %H:%M:%S %z"),
    sha256 = vapply(dat$files, sha256, character(1)),
    stringsAsFactors = FALSE
  )

  mouse_feature_id <- sub(paste0("^", mouse_prefix), "", feature_id[is_mouse])
  mouse_feature_name <- sub(paste0("^", mouse_prefix), "", feature_name[is_mouse])
  mouse_feature_key <- make.unique(mouse_feature_name)
  feature_map <- data.frame(
    feature_key = mouse_feature_key,
    source_feature_id = feature_id[is_mouse],
    source_feature_name = feature_name[is_mouse],
    mouse_ensembl_id = mouse_feature_id,
    mouse_gene_symbol = mouse_feature_name,
    source_genome = "GRCm39",
    stringsAsFactors = FALSE
  )
  if (is.null(feature_mapping_reference)) {
    feature_mapping_reference <- feature_map
  } else if (!identical(feature_mapping_reference, feature_map)) {
    stop("GRCm39 feature mapping differs across samples; refusing to merge.", call. = FALSE)
  }

  if (any(tentative_mouse)) {
    mouse_counts <- dat$counts[is_mouse, tentative_mouse, drop = FALSE]
    rownames(mouse_counts) <- mouse_feature_key
    colnames(mouse_counts) <- cell_name[tentative_mouse]
    obj <- Seurat::CreateSeuratObject(
      counts = mouse_counts,
      project = sample_id,
      min.cells = 0,
      min.features = 0
    )
    selected <- class_df[tentative_mouse, , drop = FALSE]
    rownames(selected) <- selected$cell
    obj <- Seurat::AddMetaData(obj, selected[colnames(obj), setdiff(colnames(selected), "cell"), drop = FALSE])
    assay_obj <- obj[["RNA"]]
    if ("meta.features" %in% methods::slotNames(assay_obj)) {
      assay_obj@meta.features <- feature_map[match(rownames(assay_obj), feature_map$feature_key), , drop = FALSE]
      rownames(assay_obj@meta.features) <- rownames(assay_obj)
      obj[["RNA"]] <- assay_obj
    }
    obj@misc$mouse05_species_contract <- list(
      mouse_fraction_min = mouse_fraction_min,
      human_features_retained = 0L,
      mouse_features_retained = nrow(mouse_counts),
      source_matrix_dir = matrix_dir
    )
    mouse_objects[[sample_id]] <- obj
  }

  rm(dat)
  invisible(gc())
}

classification_df <- do.call(rbind, classification_tables)
sample_summary_df <- do.call(rbind, sample_summaries)
manifest_df <- do.call(rbind, input_manifest)
rownames(classification_df) <- NULL
rownames(sample_summary_df) <- NULL
rownames(manifest_df) <- NULL

if (length(mouse_objects) == 0L) stop("No tentative mouse cells were found at the configured threshold.", call. = FALSE)
if (any(classification_df$mouse_fraction[classification_df$species_class == "tentative_mouse"] < mouse_fraction_min)) {
  stop("A retained cell violates the mouse-fraction threshold.", call. = FALSE)
}
if (any(startsWith(feature_mapping_reference$source_feature_id, human_prefix))) {
  stop("Human features leaked into the mouse feature mapping.", call. = FALSE)
}

utils::write.csv(classification_df, file.path(out_dir, "cell_species_classification.csv"), row.names = FALSE, na = "NA")
utils::write.csv(sample_summary_df, file.path(out_dir, "sample_species_summary.csv"), row.names = FALSE, na = "NA")
utils::write.csv(manifest_df, file.path(out_dir, "input_manifest_sha256.csv"), row.names = FALSE, na = "NA")
utils::write.csv(feature_mapping_reference, file.path(out_dir, "GRCm39_feature_mapping.csv"), row.names = FALSE, na = "NA")
saveRDS(mouse_objects, file.path(out_dir, "tentative_mouse_seurat_list_pre_qc.rds"), compress = FALSE)

writeLines(
  c(
    "05a species extraction completed.",
    paste0("Samples in frozen contract: ", length(sample_sources)),
    paste0("Samples with tentative mouse cells: ", length(mouse_objects)),
    paste0("Tentative mouse cells: ", sum(classification_df$species_class == "tentative_mouse")),
    paste0("Mouse fraction threshold: ", mouse_fraction_min),
    "Human features retained: 0"
  ),
  file.path(out_dir, "completion.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_05a.txt"))
message("[05a] Completed. Output: ", out_dir)
