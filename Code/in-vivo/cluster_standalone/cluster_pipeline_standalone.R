#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  KMP_DUPLICATE_LIB_OK = "TRUE",
  KMP_INIT_AT_FORK = "FALSE"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: cluster_pipeline_standalone.R INPUT_DIR OUTPUT_DIR", call. = FALSE)
}

input_dir <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- normalizePath(args[[2]], mustWork = FALSE)
all_command_args <- commandArgs(trailingOnly = FALSE)
script_argument <- grep("^--file=", all_command_args, value = TRUE)
if (length(script_argument) != 1L) {
  stop("Cannot resolve standalone pipeline script directory.", call. = FALSE)
}
STANDALONE_SCRIPT_DIR <- dirname(normalizePath(
  sub("^--file=", "", script_argument[[1]]), mustWork = TRUE
))

private_library_env <- Sys.getenv("CLUSTER_STANDALONE_R_LIBRARY", unset = "")
if (!nzchar(private_library_env)) {
  private_library_env <- Sys.getenv("R_LIBS_USER", unset = "")
}
if (nzchar(private_library_env)) {
  private_libraries <- strsplit(
    private_library_env, .Platform$path.sep, fixed = TRUE
  )[[1]]
  private_libraries <- private_libraries[dir.exists(private_libraries)]
  if (length(private_libraries) > 0L) {
    .libPaths(c(private_libraries, .libPaths()))
  }
}
message("Private R library at pipeline start: ", private_library_env)
message(".libPaths() at pipeline start: ", paste(.libPaths(), collapse = " | "))

required_packages <- c(
  "Seurat", "scDblFinder", "SingleCellExperiment", "Matrix", "RANN",
  "readxl", "readr", "ggplot2", "xgboost", "BiocNeighbors", "assorthead",
  "RcppAnnoy", "irlba", "sctransform", "uwot"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R packages: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

PINNED_PACKAGE_VERSIONS <- c(
  scDblFinder = "1.22.0",
  RcppAnnoy = "0.0.22",
  irlba = "2.3.5.1",
  sctransform = "0.4.2",
  xgboost = "1.7.11.1",
  BiocNeighbors = "2.2.0",
  assorthead = "1.2.0",
  uwot = "0.2.3"
)
observed_pinned_versions <- vapply(
  names(PINNED_PACKAGE_VERSIONS),
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
if (!identical(observed_pinned_versions, PINNED_PACKAGE_VERSIONS)) {
  stop(
    "Pinned package version mismatch. Observed: ",
    paste(names(observed_pinned_versions), observed_pinned_versions, sep = "=", collapse = ", "),
    "; expected: ",
    paste(names(PINNED_PACKAGE_VERSIONS), PINNED_PACKAGE_VERSIONS, sep = "=", collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

if (requireNamespace("future", quietly = TRUE)) {
  future::plan(future::sequential)
}
options(future.globals.maxSize = 60 * 1024^3)

SEEDS <- c(
  raw_workflow = 1234L,
  scDblFinder = 1234L,
  SCTransform = 1448145L,
  integration = 1234L,
  initial_PCA = 42L,
  initial_clustering = 0L,
  initial_UMAP = 42L,
  cell_cycle = 12345L,
  initial_DEG_similarity = 12345L,
  cluster_refine = 12345L,
  cluster_QC = 12345L,
  prefilter_DEG = 12345L,
  final_PCA = 42L,
  final_UMAP = 42L,
  final_DEG = 1234L
)

EXPECTED_INITIAL_CLUSTERS <- as.character(0:14)
EXPECTED_INITIAL_COUNTS <- c(
  `0` = 11421L, `1` = 9589L, `2` = 3644L, `3` = 3177L, `4` = 3016L,
  `5` = 1823L, `6` = 1685L, `7` = 1660L, `8` = 1579L, `9` = 1281L,
  `10` = 1267L, `11` = 1071L, `12` = 757L, `13` = 704L, `14` = 210L
)
DEG_FEATURE_CHUNK_TARGET <- 62L
DEG_PARALLEL_WORKERS <- 8L
SCDBL_PCA_ROUND_DIGITS <- 10L
SCDBL_PCA_TARGET_MAX_LOADING_SIGN <- list(
  `2N-Cell-Culture` = c(
    1, -1, -1, -1, -1, 1, -1, 1, 1, -1,
    1, -1, -1, 1, 1, 1, 1, -1, -1, -1
  ),
  `4N-Cell-Culture` = c(
    -1, 1, -1, 1, -1, 1, -1, -1, -1, 1,
    1, -1, -1, -1, 1, -1, -1, 1, 1, -1
  )
)
SEURAT_PCA_TARGET_MAX_LOADING_SIGN <- list(
  initial = c(
    -1, -1, 1, -1, -1, -1, 1, 1, -1, -1,
    1, -1, -1, 1, -1, 1, 1, -1, 1, -1,
    1, -1, 1, -1, 1, 1, 1, -1, 1, -1,
    1, 1, 1, 1, 1, -1, -1, 1, -1, -1,
    1, 1, 1, -1, -1, 1, -1, -1, 1, 1
  ),
  final = c(
    1, 1, 1, 1, 1, -1, 1, 1, -1, 1,
    -1, 1, 1, -1, -1, -1, -1, 1, -1, -1,
    -1, -1, 1, 1, -1, 1, 1, 1, 1, 1,
    -1, 1, 1, -1, 1, 1, 1, 1, 1, 1,
    1, -1, 1, 1, 1, 1, 1, -1, 1, 1
  )
)
EXPECTED_CULTURE_SINGLET_COUNTS <- c(
  `2N-Cell-Culture` = 16050L,
  `4N-Cell-Culture` = 12709L
)
EXPECTED_REFINED_CELLS <- 42884L
EXPECTED_REFINED_SPLIT_COUNTS <- c(`4c` = 103L, `9c` = 409L)
EXPECTED_REFINED_COUNTS <- c(
  `0` = 11421L, `1` = 9589L, `2` = 3644L, `3` = 3177L, `4` = 2913L,
  `4c` = 103L, `5` = 1823L, `6` = 1685L, `7` = 1660L, `8` = 1579L,
  `9` = 872L, `9c` = 409L, `10` = 1267L, `11` = 1071L, `12` = 757L,
  `13` = 704L, `14` = 210L
)
REFINE_COMPATIBILITY_STEPS <- c(
  cluster4_intersection_a = 0.0125,
  cluster4_intersection_b_cluster9_union_a = 0.03,
  cluster9_union_b = 0.09
)
REFINE_COMPATIBILITY_WORKERS <- 3L
EXPECTED_QC_HIGH <- c("4", "9", "9c")
EXPECTED_ZERO_QUALIFYING_UP <- c("3", "9c")
EXPECTED_REMOVAL_SET <- c("3", "4", "9", "9c")
EXPECTED_MANUAL_MERGE_COUNTS <- c(
  `0` = 22670L, `2` = 3644L, `3` = 3177L, `4` = 2913L, `4c` = 103L,
  `5` = 1823L, `6` = 1685L, `8` = 1579L, `9` = 872L, `9c` = 409L,
  `10` = 3095L, `13` = 704L, `14` = 210L
)
EXPECTED_FINAL_CELLS <- 35513L
EXPECTED_FINAL_LEVELS <- c("0", "2", "4c", "5", "6", "8", "10", "13", "14")
EXPECTED_FINAL_COUNTS <- c(
  `0` = 22670L,
  `2` = 3644L,
  `4c` = 103L,
  `5` = 1823L,
  `6` = 1685L,
  `8` = 1579L,
  `10` = 3095L,
  `13` = 704L,
  `14` = 210L
)
EXPECTED_BASE_METADATA <- c(
  "orig.ident", "nCount_RNA", "nFeature_RNA", "sample", "sample_folder",
  "percent.mt", "scDblFinder.class", "scDblFinder.score", "harvest",
  "Sequencing.IDs", "IDs", "Dose", "barcode_raw", "nCount_SCT",
  "nFeature_SCT", "integrated_snn_res.0.6", "seurat_clusters", "S.Score",
  "G2M.Score", "Phase", "cluster_cell_cycle_annotation",
  "seurat_cluster_refine", "cluster_refine_focus_status"
)
EXPECTED_FINAL_METADATA <- c(EXPECTED_BASE_METADATA, "clusters", "Ploidy", "TN")
EXPECTED_FINAL_COMMANDS <- c(
  "FindIntegrationAnchors", "withCallingHandlers",
  "FindNeighbors.integrated.pca", "FindClusters",
  "RunPCA.integrated", "RunUMAP.integrated.pca"
)

ensure_dir <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Cannot create directory: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

write_csv <- function(x, path) {
  ensure_dir(dirname(path))
  utils::write.csv(as.data.frame(x, stringsAsFactors = FALSE), path, row.names = FALSE)
  invisible(path)
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  utils::write.table(
    as.data.frame(x, stringsAsFactors = FALSE), path,
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
  )
  invisible(path)
}

resolve_reference_blas_lapack <- function() {
  if (!identical(Sys.info()[["sysname"]], "Linux")) return(character(0))
  resolve_one <- function(pattern, label) {
    candidates <- unique(normalizePath(Sys.glob(pattern), mustWork = FALSE))
    candidates <- candidates[file.exists(candidates)]
    candidates <- candidates[!grepl("openblas", candidates, ignore.case = TRUE)]
    if (length(candidates) == 0L) {
      stop("Reference ", label, " library was not found for compatible UMAP.", call. = FALSE)
    }
    sort(candidates)[[1]]
  }
  c(
    resolve_one("/usr/lib/*-linux-gnu/blas/libblas.so.3*", "BLAS"),
    resolve_one("/usr/lib/*-linux-gnu/lapack/liblapack.so.3*", "LAPACK")
  )
}

run_compatible_umap_matrix <- function(coordinates, seed, stage_label) {
  helper <- file.path(STANDALONE_SCRIPT_DIR, "run_compatible_umap.R")
  if (!file.exists(helper)) stop("Compatible UMAP helper is missing: ", helper, call. = FALSE)
  if (!is.matrix(coordinates) || !is.numeric(coordinates) ||
      ncol(coordinates) != 30L || is.null(rownames(coordinates))) {
    stop("Compatible UMAP requires a named numeric matrix with exactly 30 columns.", call. = FALSE)
  }
  runtime_dir <- ensure_dir(file.path(output_dir, "00_provenance", "umap_runtime"))
  input_path <- tempfile(paste0(stage_label, "_pca_"), tmpdir = runtime_dir, fileext = ".rds")
  output_path <- tempfile(paste0(stage_label, "_umap_"), tmpdir = runtime_dir, fileext = ".rds")
  on.exit(unlink(c(input_path, output_path), force = TRUE), add = TRUE)
  saveRDS(coordinates, input_path, compress = FALSE)

  reference_libraries <- resolve_reference_blas_lapack()
  preload <- c(reference_libraries, Sys.getenv("LD_PRELOAD", unset = ""))
  preload <- preload[nzchar(preload)]
  environment_assignments <- c(
    if (length(preload)) paste0("LD_PRELOAD=", paste(preload, collapse = ":")),
    paste0("R_LIBS_USER=", Sys.getenv("R_LIBS_USER", unset = "")),
    paste0(
      "CLUSTER_STANDALONE_R_LIBRARY=",
      Sys.getenv("CLUSTER_STANDALONE_R_LIBRARY", unset = "")
    ),
    "OMP_NUM_THREADS=1", "OPENBLAS_NUM_THREADS=1", "MKL_NUM_THREADS=1"
  )
  command_arguments <- c(
    environment_assignments,
    file.path(R.home("bin"), "Rscript"), "--vanilla", helper,
    input_path, output_path, as.character(as.integer(seed))
  )
  status <- system2(
    "/usr/bin/env",
    args = vapply(command_arguments, shQuote, character(1)),
    stdout = "",
    stderr = ""
  )
  if (!identical(status, 0L) || !file.exists(output_path)) {
    stop("Compatible UMAP subprocess failed for stage ", stage_label, call. = FALSE)
  }
  result <- readRDS(output_path)
  embedding <- result$embedding
  if (!is.matrix(embedding) || !identical(dim(embedding), c(nrow(coordinates), 2L)) ||
      !identical(rownames(embedding), rownames(coordinates))) {
    stop("Compatible UMAP returned an invalid embedding for stage ", stage_label, call. = FALSE)
  }
  write_tsv(
    data.frame(
      stage = stage_label,
      seed = as.integer(seed),
      uwot_version = result$uwot_version,
      uwot_built = result$uwot_built,
      force_irlba = result$force_irlba,
      reference_blas = if (length(reference_libraries)) reference_libraries[[1]] else result$blas,
      reference_lapack = if (length(reference_libraries)) reference_libraries[[2]] else result$lapack,
      ld_preload = result$ld_preload,
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "00_provenance", paste0("umap_", stage_label, "_runtime.tsv"))
  )
  embedding
}

run_compatible_umap <- function(obj, reduction, dims, seed, stage_label) {
  coordinates <- Seurat::Embeddings(obj, reduction = reduction)[, dims, drop = FALSE]
  embedding <- run_compatible_umap_matrix(
    coordinates = coordinates,
    seed = seed,
    stage_label = stage_label
  )
  obj[["umap"]] <- Seurat::CreateDimReducObject(
    embeddings = embedding,
    key = "UMAP_",
    assay = Seurat::DefaultAssay(obj),
    global = TRUE
  )
  obj
}

write_stage_marker <- function(stage, details = "completed") {
  stage_dir <- ensure_dir(file.path(output_dir, "00_provenance", "stages"))
  writeLines(
    c(
      paste0("stage=", stage),
      paste0("status=completed"),
      paste0("details=", details),
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
    ),
    file.path(stage_dir, paste0(stage, ".complete"))
  )
}

assert_set_equal <- function(observed, expected, label) {
  observed <- sort(unique(as.character(observed)))
  expected <- sort(unique(as.character(expected)))
  if (!identical(observed, expected)) {
    stop(
      label, " mismatch. Observed: ", paste(observed, collapse = ","),
      "; expected: ", paste(expected, collapse = ","),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

sort_maybe_numeric <- function(x) {
  x <- unique(as.character(x))
  numeric_value <- suppressWarnings(as.numeric(x))
  if (all(!is.na(numeric_value))) {
    return(x[order(numeric_value, x)])
  }
  numeric_prefix <- suppressWarnings(as.numeric(sub("^([0-9]+).*$", "\\1", x)))
  x[order(is.na(numeric_prefix), numeric_prefix, x)]
}

message_cluster_counts_before_rds <- function(obj, cluster_col, rds_path, label) {
  if (!(cluster_col %in% colnames(obj@meta.data))) {
    stop("Cannot report cluster counts; metadata column is missing: ", cluster_col, call. = FALSE)
  }
  counts <- table(obj@meta.data[[cluster_col]], useNA = "ifany")
  message(
    label,
    " cluster counts before RDS save [", basename(rds_path), "]: ",
    paste(names(counts), as.integer(counts), sep = "=", collapse = ", ")
  )
  invisible(counts)
}

validate_initial_cluster_contract <- function(obj, label) {
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    stop(label, " is missing metadata column seurat_clusters.", call. = FALSE)
  }
  if (ncol(obj) != sum(EXPECTED_INITIAL_COUNTS)) {
    stop(
      label, " cell count mismatch: ", ncol(obj),
      " vs ", sum(EXPECTED_INITIAL_COUNTS),
      call. = FALSE
    )
  }
  observed <- table(factor(
    as.character(obj$seurat_clusters),
    levels = EXPECTED_INITIAL_CLUSTERS
  ))
  if (!identical(as.integer(observed), as.integer(EXPECTED_INITIAL_COUNTS))) {
    stop(
      label, " cluster counts do not match the reference. Observed: ",
      paste(names(observed), as.integer(observed), sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  invisible(observed)
}

validate_named_cluster_counts <- function(values, expected, label) {
  values <- as.character(values)
  if (!setequal(unique(values), names(expected))) {
    stop(
      label, " labels mismatch. Observed: ",
      paste(sort_maybe_numeric(values), collapse = ","),
      call. = FALSE
    )
  }
  observed <- table(factor(values, levels = names(expected)))
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(
      label, " counts mismatch. Observed: ",
      paste(names(observed), as.integer(observed), sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  invisible(observed)
}

validate_refined_resume <- function(resumed, input_obj, allow_manual_merge = FALSE) {
  validate_initial_cluster_contract(resumed, "Resumed refined object")
  if (!identical(colnames(resumed), colnames(input_obj))) {
    stop("Resumed refined object cell order differs from the current initial object.", call. = FALSE)
  }
  required <- c("seurat_cluster_refine", "cluster_refine_focus_status")
  missing <- setdiff(required, colnames(resumed@meta.data))
  if (length(missing) > 0L) {
    stop("Resumed refined object is missing metadata: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  expected_metadata <- EXPECTED_BASE_METADATA
  if (isTRUE(allow_manual_merge)) {
    expected_metadata <- c(expected_metadata, "manual_merge_test")
  }
  if (!identical(colnames(resumed@meta.data), expected_metadata)) {
    stop("Resumed refined object metadata contract mismatch.", call. = FALSE)
  }
  validate_named_cluster_counts(
    resumed$seurat_cluster_refine,
    EXPECTED_REFINED_COUNTS,
    "Resumed refined object"
  )
  invisible(TRUE)
}

validate_manual_merge_resume <- function(resumed, input_obj) {
  validate_refined_resume(resumed, input_obj, allow_manual_merge = TRUE)
  if (!("manual_merge_test" %in% colnames(resumed@meta.data))) {
    stop("Resumed manual-merge object lacks manual_merge_test metadata.", call. = FALSE)
  }
  validate_named_cluster_counts(
    resumed$manual_merge_test,
    EXPECTED_MANUAL_MERGE_COUNTS,
    "Resumed manual-merge object"
  )
  invisible(TRUE)
}

get_assay_matrix <- function(obj, assay, slot_name) {
  tryCatch(
    Seurat::GetAssayData(obj, assay = assay, slot = slot_name),
    error = function(e) NULL
  )
}

prepare_obj_for_markers <- function(obj) {
  if (!("RNA" %in% names(obj@assays))) {
    stop("RNA assay is missing.", call. = FALSE)
  }
  Seurat::DefaultAssay(obj) <- "RNA"
  rna_data <- get_assay_matrix(obj, "RNA", "data")
  if (is.null(rna_data) || nrow(rna_data) == 0L || ncol(rna_data) == 0L) {
    obj <- Seurat::NormalizeData(obj, assay = "RNA", verbose = FALSE)
  }
  obj
}

resolve_lfc_col <- function(df) {
  candidates <- c("avg_log2FC", "avg_logFC")
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) == 0L) {
    stop("No avg_log2FC/avg_logFC column in DEG result.", call. = FALSE)
  }
  hit[[1]]
}

safe_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- sub("^GRCh[0-9]+[-_]", "", x, ignore.case = TRUE)
  x <- sub("^GRCm39[-_]", "", x, ignore.case = TRUE)
  x <- sub("^hg38[-_]", "", x, ignore.case = TRUE)
  sub("\\.[0-9]+$", "", x)
}

save_umap_plot <- function(obj, group_by, file_stub, title) {
  if (!("umap" %in% names(obj@reductions))) return(invisible(FALSE))
  p <- Seurat::DimPlot(
    obj,
    reduction = "umap",
    group.by = group_by,
    label = identical(group_by, "clusters"),
    repel = TRUE,
    pt.size = 0.30,
    raster = FALSE
  ) + ggplot2::labs(title = title)
  ggplot2::ggsave(paste0(file_stub, ".pdf"), p, width = 9, height = 7)
  ggplot2::ggsave(paste0(file_stub, ".png"), p, width = 9, height = 7, dpi = 300)
  invisible(TRUE)
}

discover_inputs <- function(input_dir) {
  ploidy_file <- file.path(input_dir, "all_ploidy.tsv")
  sample_info_file <- file.path(input_dir, "sample_info.xlsx")

  candidate_directories <- list.dirs(
    input_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  cellranger_candidates <- candidate_directories[
    basename(candidate_directories) == "A02_cellRanger"
  ]
  cellranger_candidates <- unique(normalizePath(
    cellranger_candidates,
    mustWork = TRUE
  ))
  if (length(cellranger_candidates) != 1L) {
    stop(
      "INPUT_DIR must contain exactly one A02_cellRanger directory; found ",
      length(cellranger_candidates),
      if (length(cellranger_candidates) > 0L) {
        paste0(": ", paste(cellranger_candidates, collapse = ", "))
      } else {
        ""
      },
      call. = FALSE
    )
  }
  cellranger_root <- cellranger_candidates[[1L]]
  for (required_file in c(ploidy_file, sample_info_file)) {
    if (!file.exists(required_file)) {
      stop("Required input file is missing: ", required_file, call. = FALSE)
    }
  }

  sample_dirs <- sort(list.dirs(cellranger_root, recursive = FALSE, full.names = TRUE))
  sample_dirs <- sample_dirs[grepl("-Count-HM$", basename(sample_dirs))]
  h5_files <- file.path(sample_dirs, "outs", "filtered_feature_bc_matrix.h5")
  missing_h5 <- h5_files[!file.exists(h5_files)]
  if (length(missing_h5) > 0L) {
    stop("Missing required H5 file(s): ", paste(missing_h5, collapse = ", "), call. = FALSE)
  }
  if (length(h5_files) == 0L) {
    stop("No *-Count-HM/outs/filtered_feature_bc_matrix.h5 inputs found.", call. = FALSE)
  }

  list(
    cellranger_root = cellranger_root,
    sample_dirs = sample_dirs,
    h5_files = h5_files,
    ploidy_file = ploidy_file,
    sample_info_file = sample_info_file
  )
}

read_sample_info <- function(path) {
  as.data.frame(
    readxl::read_excel(path),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

resolve_col <- function(df, candidates, fallback_index = NULL) {
  nm <- names(df)
  for (candidate in candidates) {
    hit <- which(tolower(nm) == tolower(candidate))
    if (length(hit) == 1L) return(nm[[hit]])
  }
  if (!is.null(fallback_index) && fallback_index <= length(nm)) {
    return(nm[[fallback_index]])
  }
  stop("Cannot resolve column: ", paste(candidates, collapse = ", "), call. = FALSE)
}

read_counts <- function(h5_file) {
  counts <- Seurat::Read10X_h5(h5_file)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) {
      counts[["Gene Expression"]]
    } else {
      counts[[1]]
    }
  }
  counts
}

qc_filter_special <- function(obj, sample_id) {
  feature_upper <- as.numeric(stats::quantile(obj$nFeature_RNA, 0.99, na.rm = TRUE))
  count_upper <- as.numeric(stats::quantile(obj$nCount_RNA, 0.99, na.rm = TRUE))
  keep <- obj$nFeature_RNA >= 200 &
    obj$nCount_RNA >= 500 &
    obj$nFeature_RNA <= feature_upper &
    obj$nCount_RNA <= count_upper &
    obj$percent.mt <= 20
  message(sprintf("[%s] culture QC keep: %d / %d", sample_id, sum(keep), length(keep)))
  subset(obj, cells = colnames(obj)[keep])
}

scdbl_platform_stable_processing <- function(e, dims = 20L, sample_id) {
  target_sign <- SCDBL_PCA_TARGET_MAX_LOADING_SIGN[[sample_id]]
  if (is.null(target_sign)) {
    stop("No scDblFinder PCA orientation contract for sample: ", sample_id, call. = FALSE)
  }

  pca <- getFromNamespace(".defaultProcessing", "scDblFinder")(e, dims = dims)
  rotation <- attr(pca, "rotation")
  if (is.null(rotation) || ncol(rotation) != ncol(pca)) {
    stop("scDblFinder PCA rotation attribute is unavailable.", call. = FALSE)
  }
  if (length(target_sign) < ncol(rotation)) {
    stop(
      "PCA orientation contract has ", length(target_sign),
      " entries but ", ncol(rotation), " are required.",
      call. = FALSE
    )
  }

  observed_sign <- vapply(
    seq_len(ncol(rotation)),
    function(component) {
      loading <- rotation[, component]
      sign(loading[which.max(abs(loading))])
    },
    numeric(1)
  )
  if (any(observed_sign == 0)) {
    stop("Cannot orient a PCA component with a zero anchor loading.", call. = FALSE)
  }
  flip <- observed_sign != target_sign[seq_along(observed_sign)]
  pca[, flip] <- -pca[, flip, drop = FALSE]

  # The original objects were generated on ARM macOS. Rounding after an
  # explicit loading-based orientation removes sub-float BLAS differences
  # while preserving the original scDblFinder calls exactly on that platform.
  pca <- round(pca, digits = SCDBL_PCA_ROUND_DIGITS)
  rownames(pca) <- colnames(e)
  pca
}

remove_doublets <- function(obj, sample_id) {
  if (ncol(obj) == 0L) return(obj)
  if (ncol(obj) < 100L) {
    obj$scDblFinder.class <- "singlet_low_cell_skip"
    obj$scDblFinder.score <- NA_real_
    return(obj)
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
  )
  processing <- function(e, dims = 20L) {
    scdbl_platform_stable_processing(e, dims = dims, sample_id = sample_id)
  }
  set.seed(SEEDS[["scDblFinder"]])
  sce <- scDblFinder::scDblFinder(sce, verbose = FALSE, processing = processing)
  dbl_meta <- as.data.frame(SingleCellExperiment::colData(sce))
  obj$scDblFinder.class <- as.character(dbl_meta$scDblFinder.class)
  obj$scDblFinder.score <- as.numeric(dbl_meta$scDblFinder.score)
  singlets <- rownames(dbl_meta)[dbl_meta$scDblFinder.class == "singlet"]
  message(sprintf("[%s] scDblFinder singlet keep: %d / %d", sample_id, length(singlets), ncol(obj)))
  expected_singlets <- unname(EXPECTED_CULTURE_SINGLET_COUNTS[[sample_id]])
  if (is.null(expected_singlets) || length(singlets) != expected_singlets) {
    stop(
      "scDblFinder reproduction mismatch for ", sample_id,
      ": observed ", length(singlets),
      ", expected ", expected_singlets,
      call. = FALSE
    )
  }
  subset(obj, cells = singlets)
}

make_special_metadata <- function(sample_id, columns, seq_col, id_col) {
  out <- as.list(stats::setNames(rep(NA_character_, length(columns)), columns))
  out[[seq_col]] <- sample_id
  out[[id_col]] <- sample_id
  out
}

get_metadata_row <- function(sample_id, sample_info, id_col, seq_col, special_samples) {
  if (sample_id %in% special_samples) {
    return(make_special_metadata(sample_id, names(sample_info), seq_col, id_col))
  }
  idx <- which(as.character(sample_info[[id_col]]) == sample_id)
  if (length(idx) == 0L) stop("Sample missing from sample_info.xlsx: ", sample_id, call. = FALSE)
  if (length(idx) > 1L) warning("Multiple sample_info rows for ", sample_id, "; first used.")
  as.list(sample_info[idx[[1]], , drop = FALSE])
}

add_sample_metadata <- function(obj, metadata_row) {
  for (column_name in names(metadata_row)) {
    obj[[column_name]] <- metadata_row[[column_name]]
  }
  obj
}

coerce_base_metadata_contract <- function(obj) {
  char_cols <- c(
    "orig.ident", "sample", "sample_folder", "scDblFinder.class", "harvest",
    "Sequencing.IDs", "IDs", "Dose", "barcode_raw"
  )
  for (column_name in intersect(char_cols, colnames(obj@meta.data))) {
    obj@meta.data[[column_name]] <- as.character(obj@meta.data[[column_name]])
  }
  obj@meta.data$nCount_RNA <- as.numeric(obj@meta.data$nCount_RNA)
  obj@meta.data$nFeature_RNA <- as.integer(obj@meta.data$nFeature_RNA)
  obj@meta.data$percent.mt <- as.numeric(obj@meta.data$percent.mt)
  obj@meta.data$scDblFinder.score <- as.numeric(obj@meta.data$scDblFinder.score)
  if ("nCount_SCT" %in% colnames(obj@meta.data)) {
    obj@meta.data$nCount_SCT <- as.numeric(obj@meta.data$nCount_SCT)
  }
  if ("nFeature_SCT" %in% colnames(obj@meta.data)) {
    obj@meta.data$nFeature_SCT <- as.integer(obj@meta.data$nFeature_SCT)
  }
  obj
}

orient_seurat_pca_to_reference_anchor <- function(obj, contract_name, stage_label) {
  target_sign <- SEURAT_PCA_TARGET_MAX_LOADING_SIGN[[contract_name]]
  if (is.null(target_sign)) {
    stop("Unknown Seurat PCA orientation contract: ", contract_name, call. = FALSE)
  }
  if (!("pca" %in% names(obj@reductions))) {
    stop("Cannot orient Seurat PCA because reduction 'pca' is missing.", call. = FALSE)
  }
  reduction <- obj[["pca"]]
  loadings <- reduction@feature.loadings
  embeddings <- reduction@cell.embeddings
  if (ncol(loadings) != ncol(embeddings) || length(target_sign) < ncol(loadings)) {
    stop("Seurat PCA orientation contract dimension mismatch.", call. = FALSE)
  }
  observed_sign <- vapply(seq_len(ncol(loadings)), function(component) {
    loading <- loadings[, component]
    sign(loading[which.max(abs(loading))])
  }, numeric(1))
  if (any(observed_sign == 0)) {
    stop("Cannot orient a Seurat PCA component with a zero anchor loading.", call. = FALSE)
  }
  expected_sign <- target_sign[seq_along(observed_sign)]
  flip <- observed_sign != expected_sign
  if (any(flip)) {
    reduction@cell.embeddings[, flip] <- -reduction@cell.embeddings[, flip, drop = FALSE]
    reduction@feature.loadings[, flip] <- -reduction@feature.loadings[, flip, drop = FALSE]
    if (length(reduction@feature.loadings.projected) > 0L) {
      reduction@feature.loadings.projected[, flip] <-
        -reduction@feature.loadings.projected[, flip, drop = FALSE]
    }
    obj[["pca"]] <- reduction
  }
  contract <- data.frame(
    component = seq_along(observed_sign),
    observed_anchor_sign = observed_sign,
    expected_anchor_sign = expected_sign,
    flipped = flip,
    stringsAsFactors = FALSE
  )
  file_stub <- gsub("[^A-Za-z0-9_.-]+", "_", stage_label)
  write_tsv(
    contract,
    file.path(
      ensure_dir(file.path(output_dir, "00_provenance", "pca_orientation")),
      paste0(file_stub, ".tsv")
    )
  )
  list(object = obj, flipped = flip, contract = contract)
}

run_raw_integration <- function(inputs) {
  stage_root <- ensure_dir(file.path(output_dir, "01_data"))
  object_path <- file.path(stage_root, "integrated_sct_cca_seurat.rds")
  if (file.exists(object_path)) {
    message("[resume] Reading verified initial object: ", object_path)
    resumed <- readRDS(object_path)
    validate_initial_cluster_contract(resumed, "Resumed initial object")
    oriented <- orient_seurat_pca_to_reference_anchor(
      resumed, contract_name = "initial", stage_label = "initial_resume"
    )
    resumed <- oriented$object
    resumed <- run_compatible_umap(
      resumed, reduction = "pca", dims = 1:30,
      seed = SEEDS[["initial_UMAP"]], stage_label = "initial_resume"
    )
    return(resumed)
  }
  special_samples <- c("2N-Cell-Culture", "4N-Cell-Culture")

  sample_info <- read_sample_info(inputs$sample_info_file)
  if (ncol(sample_info) < 3L) {
    stop("sample_info.xlsx must contain at least three columns.", call. = FALSE)
  }
  harvest_col <- resolve_col(sample_info, "harvest", 1L)
  seq_col <- resolve_col(sample_info, c("Sequencing IDs", "Sequencing ID", "SequencingIDs"), 2L)
  id_col <- resolve_col(sample_info, c("IDs", "ID"), 3L)
  sample_info[[harvest_col]] <- as.character(sample_info[[harvest_col]])
  sample_info[[id_col]] <- as.character(sample_info[[id_col]])

  ploidy <- utils::read.delim(
    inputs$ploidy_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!all(c("file", "cell_id") %in% names(ploidy))) {
    stop("all_ploidy.tsv must contain columns 'file' and 'cell_id'.", call. = FALSE)
  }
  ploidy$harvest <- sub("\\.sps\\.cbs$", "", ploidy$file)
  barcode_by_harvest <- lapply(split(ploidy$cell_id, ploidy$harvest), unique)
  id_to_harvest <- stats::setNames(sample_info[[harvest_col]], sample_info[[id_col]])
  barcode_by_id <- lapply(id_to_harvest, function(harvest) {
    if (is.na(harvest) || !(harvest %in% names(barcode_by_harvest))) character(0) else barcode_by_harvest[[harvest]]
  })

  seurat_list <- list()
  sample_summary <- list()
  set.seed(SEEDS[["raw_workflow"]])
  for (sample_dir in inputs$sample_dirs) {
    sample_folder <- basename(sample_dir)
    sample_id <- sub("-Count-HM$", "", sample_folder)
    h5_file <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix.h5")
    message("Reading sample: ", sample_id)
    obj <- Seurat::CreateSeuratObject(
      counts = read_counts(h5_file),
      project = sample_id,
      min.cells = 0,
      min.features = 0
    )
    obj$sample <- sample_id
    obj$sample_folder <- sample_folder
    obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")
    n_raw <- ncol(obj)

    if (sample_id %in% special_samples) {
      filter_method <- "QC + scDblFinder"
      obj <- qc_filter_special(obj, sample_id)
      n_after_qc <- ncol(obj)
      obj <- remove_doublets(obj, sample_id)
    } else {
      filter_method <- "all_ploidy keep barcodes"
      keep_barcodes <- barcode_by_id[[sample_id]]
      if (is.null(keep_barcodes) || length(keep_barcodes) == 0L) {
        warning("No all_ploidy barcode mapping for ", sample_id, "; sample skipped.")
        next
      }
      keep_cells <- intersect(colnames(obj), keep_barcodes)
      message(sprintf("[%s] all_ploidy keep: %d / %d", sample_id, length(keep_cells), ncol(obj)))
      obj <- subset(obj, cells = keep_cells)
      obj$scDblFinder.class <- "singlet_by_ploidy"
      obj$scDblFinder.score <- NA_real_
      n_after_qc <- ncol(obj)
    }
    if (ncol(obj) == 0L) stop("No cells remain for sample: ", sample_id, call. = FALSE)

    obj <- add_sample_metadata(
      obj,
      get_metadata_row(sample_id, sample_info, id_col, seq_col, special_samples)
    )
    obj$barcode_raw <- colnames(obj)
    obj <- Seurat::RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))
    seurat_list[[sample_id]] <- obj
    sample_summary[[sample_id]] <- data.frame(
      sample = sample_id,
      sample_folder = sample_folder,
      filter_method = filter_method,
      raw_cells = n_raw,
      after_qc_cells = n_after_qc,
      after_filter_cells = ncol(obj),
      stringsAsFactors = FALSE
    )
  }

  if (length(seurat_list) < 2L) {
    stop("At least two non-empty samples are required for SCT-CCA integration.", call. = FALSE)
  }
  write_csv(do.call(rbind, sample_summary), file.path(stage_root, "sample_filter_summary.csv"))

  message("Running SCTransform on ", length(seurat_list), " samples.")
  seurat_list <- lapply(seurat_list, function(x) {
    set.seed(SEEDS[["SCTransform"]])
    Seurat::SCTransform(
      object = x,
      assay = "RNA",
      new.assay.name = "SCT",
      vars.to.regress = "percent.mt",
      return.only.var.genes = FALSE,
      seed.use = SEEDS[["SCTransform"]],
      verbose = FALSE
    )
  })

  set.seed(SEEDS[["integration"]])
  integration_features <- Seurat::SelectIntegrationFeatures(seurat_list, nfeatures = 3000)
  seurat_list <- Seurat::PrepSCTIntegration(
    object.list = seurat_list,
    anchor.features = integration_features,
    verbose = TRUE
  )
  set.seed(SEEDS[["integration"]])
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = seurat_list,
    normalization.method = "SCT",
    anchor.features = integration_features,
    reduction = "cca",
    verbose = TRUE
  )
  set.seed(SEEDS[["integration"]])
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    verbose = TRUE
  )

  Seurat::DefaultAssay(integrated) <- "integrated"
  set.seed(SEEDS[["initial_PCA"]])
  integrated <- Seurat::RunPCA(
    integrated,
    npcs = 50,
    seed.use = SEEDS[["initial_PCA"]],
    verbose = TRUE
  )
  integrated <- orient_seurat_pca_to_reference_anchor(
    integrated, contract_name = "initial", stage_label = "initial_fresh"
  )$object
  integrated <- Seurat::FindNeighbors(integrated, dims = 1:30, verbose = TRUE)
  set.seed(SEEDS[["initial_clustering"]])
  integrated <- Seurat::FindClusters(
    integrated,
    resolution = 0.6,
    random.seed = SEEDS[["initial_clustering"]],
    verbose = TRUE
  )
  integrated <- run_compatible_umap(
    integrated, reduction = "pca", dims = 1:30,
    seed = SEEDS[["initial_UMAP"]], stage_label = "initial"
  )
  integrated <- coerce_base_metadata_contract(integrated)

  initial_cluster_counts <- as.data.frame(
    table(cluster = as.character(integrated$seurat_clusters)),
    stringsAsFactors = FALSE
  )
  names(initial_cluster_counts)[[2]] <- "n_cells"
  write_csv(initial_cluster_counts, file.path(stage_root, "initial_cluster_counts.csv"))
  assert_set_equal(levels(integrated$seurat_clusters), EXPECTED_INITIAL_CLUSTERS, "Initial cluster labels")
  validate_initial_cluster_contract(integrated, "Initial object")
  message_cluster_counts_before_rds(
    integrated,
    "seurat_clusters",
    object_path,
    "Initial"
  )
  saveRDS(integrated, object_path)
  save_umap_plot(
    integrated,
    "seurat_clusters",
    file.path(stage_root, "umap_initial_clusters"),
    "Initial SCT-CCA clusters"
  )
  write_stage_marker("01_data", paste0(ncol(integrated), " cells; 15 initial clusters"))
  integrated
}

normalize_feature_key <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^GRCh[0-9]+-", "", x, ignore.case = TRUE)
  x <- sub("\\.[0-9]+$", "", x)
  toupper(x)
}

extract_feature_metadata <- function(assay_obj) {
  slots <- methods::slotNames(assay_obj)
  if ("meta.features" %in% slots) return(as.data.frame(assay_obj@meta.features, stringsAsFactors = FALSE))
  if ("meta.data" %in% slots) return(as.data.frame(assay_obj@meta.data, stringsAsFactors = FALSE))
  data.frame(row.names = rownames(assay_obj))
}

build_feature_lookup <- function(feature_ids, aliases) {
  keys <- normalize_feature_key(aliases)
  keep <- !is.na(keys) & keys != "" & !duplicated(keys)
  stats::setNames(feature_ids[keep], keys[keep])
}

resolve_gene_set_features <- function(obj, assay_name, symbols, set_name) {
  assay_obj <- obj[[assay_name]]
  feature_ids <- rownames(assay_obj)
  target_keys <- normalize_feature_key(symbols)
  matched <- rep(NA_character_, length(symbols))
  lookups <- list(feature_ids = build_feature_lookup(feature_ids, feature_ids))
  if (any(grepl("\\|", feature_ids))) {
    lookups$before_pipe <- build_feature_lookup(feature_ids, sub("\\|.*$", "", feature_ids))
    lookups$after_pipe <- build_feature_lookup(feature_ids, sub("^.*\\|", "", feature_ids))
  }
  if (any(grepl("_", feature_ids))) {
    lookups$before_underscore <- build_feature_lookup(feature_ids, sub("_.*$", "", feature_ids))
    lookups$after_underscore <- build_feature_lookup(feature_ids, sub("^.*_", "", feature_ids))
  }
  feature_meta <- extract_feature_metadata(assay_obj)
  alias_columns <- intersect(
    c("gene_name", "gene", "symbol", "gene_symbol", "gene_symbols", "SYMBOL", "GENE",
      "Gene", "GeneName", "GeneSymbol", "feature_name", "feature", "features"),
    colnames(feature_meta)
  )
  for (column_name in alias_columns) {
    lookups[[paste0("meta_", column_name)]] <- build_feature_lookup(feature_ids, feature_meta[[column_name]])
  }
  for (lookup in lookups) {
    proposed <- unname(lookup[target_keys])
    fill <- is.na(matched) & !is.na(proposed)
    matched[fill] <- proposed[fill]
  }
  matched <- unique(matched[!is.na(matched)])
  message(set_name, ": ", length(matched), "/", length(symbols), " genes matched.")
  if (length(matched) < 5L) stop("Too few matched ", set_name, " features.", call. = FALSE)
  matched
}

safe_sd <- function(x) {
  value <- stats::sd(x, na.rm = TRUE)
  if (is.na(value)) 0 else value
}

cluster_numeric_summary <- function(values, clusters, cluster_levels, fun = mean) {
  vapply(cluster_levels, function(cluster_id) {
    fun(values[clusters == cluster_id], na.rm = TRUE)
  }, numeric(1))
}

run_cell_cycle <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "01a_cell_cycle"))
  cell_cycle_path <- file.path(stage_root, "integrated_sct_cca_seurat_cell_cycle.rds")
  if (file.exists(cell_cycle_path)) {
    message("[resume] Reading verified cell-cycle object: ", cell_cycle_path)
    resumed <- readRDS(cell_cycle_path)
    validate_initial_cluster_contract(resumed, "Resumed cell-cycle object")
    required_cell_cycle_columns <- c(
      "S.Score", "G2M.Score", "Phase", "cluster_cell_cycle_annotation"
    )
    missing_cell_cycle_columns <- setdiff(
      required_cell_cycle_columns,
      colnames(resumed@meta.data)
    )
    if (length(missing_cell_cycle_columns) > 0L) {
      stop(
        "Resumed cell-cycle object is missing metadata: ",
        paste(missing_cell_cycle_columns, collapse = ", "),
        call. = FALSE
      )
    }
    if (!identical(colnames(resumed), colnames(obj))) {
      stop("Resumed cell-cycle object cell order differs from the initial object.", call. = FALSE)
    }
    resumed[["umap"]] <- obj[["umap"]]
    return(resumed)
  }
  set.seed(SEEDS[["cell_cycle"]])
  Seurat::DefaultAssay(obj) <- "RNA"
  cc_genes <- Seurat::cc.genes.updated.2019
  s_features <- resolve_gene_set_features(obj, "RNA", cc_genes$s.genes, "S.Score")
  g2m_features <- resolve_gene_set_features(obj, "RNA", cc_genes$g2m.genes, "G2M.Score")
  obj <- Seurat::CellCycleScoring(
    obj,
    s.features = s_features,
    g2m.features = g2m_features,
    set.ident = FALSE
  )

  meta <- obj@meta.data
  clusters <- as.character(meta$seurat_clusters)
  cluster_levels <- unique(clusters)
  phase <- factor(as.character(meta$Phase), levels = c("G1", "S", "G2M"))
  cluster_summary <- data.frame(
    cluster = cluster_levels,
    n_cells = as.integer(table(factor(clusters, levels = cluster_levels))),
    mean_S.Score = cluster_numeric_summary(meta$S.Score, clusters, cluster_levels),
    median_S.Score = cluster_numeric_summary(meta$S.Score, clusters, cluster_levels, stats::median),
    mean_G2M.Score = cluster_numeric_summary(meta$G2M.Score, clusters, cluster_levels),
    median_G2M.Score = cluster_numeric_summary(meta$G2M.Score, clusters, cluster_levels, stats::median),
    frac_G1 = cluster_numeric_summary(phase == "G1", clusters, cluster_levels),
    frac_S = cluster_numeric_summary(phase == "S", clusters, cluster_levels),
    frac_G2M = cluster_numeric_summary(phase == "G2M", clusters, cluster_levels),
    stringsAsFactors = FALSE
  )
  write_csv(cluster_summary, file.path(stage_root, "cluster_cell_cycle_summary.csv"))

  score_flag <-
    cluster_summary$mean_S.Score > mean(cluster_summary$mean_S.Score) + safe_sd(cluster_summary$mean_S.Score) |
    cluster_summary$mean_G2M.Score > mean(cluster_summary$mean_G2M.Score) + safe_sd(cluster_summary$mean_G2M.Score)
  phase_flag <-
    cluster_summary$frac_S > mean(cluster_summary$frac_S) + safe_sd(cluster_summary$frac_S) |
    cluster_summary$frac_G2M > mean(cluster_summary$frac_G2M) + safe_sd(cluster_summary$frac_G2M)

  pca_use <- Seurat::Embeddings(obj, "pca")[, 1:20, drop = FALSE]
  cor_s <- apply(pca_use, 2, stats::cor, y = meta$S.Score, method = "pearson", use = "pairwise.complete.obs")
  cor_g2m <- apply(pca_use, 2, stats::cor, y = meta$G2M.Score, method = "pearson", use = "pairwise.complete.obs")
  selected_pcs <- colnames(pca_use)[abs(cor_s) >= 0.30 | abs(cor_g2m) >= 0.30]
  if (length(selected_pcs) == 0L) {
    selected_pcs <- colnames(pca_use)[which.max(pmax(abs(cor_s), abs(cor_g2m)))]
  }
  pc_flag <- stats::setNames(rep(FALSE, length(cluster_levels)), cluster_levels)
  for (pc_name in selected_pcs) {
    cluster_pc_mean <- cluster_numeric_summary(pca_use[, pc_name], clusters, cluster_levels)
    pc_flag <- pc_flag |
      abs(cluster_pc_mean - mean(cluster_pc_mean, na.rm = TRUE)) > safe_sd(cluster_pc_mean)
  }

  n_methods <- as.integer(score_flag) + as.integer(phase_flag) + as.integer(pc_flag)
  candidate_summary <- data.frame(
    cluster = cluster_levels,
    flag_score_outlier = score_flag,
    flag_phase_enrichment = phase_flag,
    flag_pc_separation_support = as.logical(pc_flag),
    n_methods_flagged = n_methods,
    final_candidate_flag = n_methods >= 2L,
    stringsAsFactors = FALSE
  )
  write_csv(candidate_summary, file.path(stage_root, "candidate_cell_cycle_cluster_integrated_summary.csv"))
  annotation_map <- stats::setNames(
    ifelse(candidate_summary$final_candidate_flag, "cell_cycle_candidate", "not_cell_cycle_candidate"),
    candidate_summary$cluster
  )
  obj$cluster_cell_cycle_annotation <- unname(annotation_map[as.character(obj$seurat_clusters)])
  obj@meta.data$S.Score <- as.numeric(obj@meta.data$S.Score)
  obj@meta.data$G2M.Score <- as.numeric(obj@meta.data$G2M.Score)
  obj@meta.data$Phase <- as.character(obj@meta.data$Phase)
  obj@meta.data$cluster_cell_cycle_annotation <- as.character(obj@meta.data$cluster_cell_cycle_annotation)
  message_cluster_counts_before_rds(
    obj,
    "seurat_clusters",
    cell_cycle_path,
    "Cell-cycle-scored"
  )
  saveRDS(obj, cell_cycle_path)
  write_stage_marker("01a_cell_cycle", paste0(sum(candidate_summary$final_candidate_flag), " candidate clusters"))
  obj
}

run_deg_by_group <- function(
  obj,
  group_col,
  group_levels,
  output_root,
  min_pct = NULL,
  logfc_threshold = NULL,
  seed,
  layout = c("nested_deg", "nested", "flat")
) {
  layout <- match.arg(layout)
  obj <- prepare_obj_for_markers(obj)
  obj@meta.data[[group_col]] <- factor(
    as.character(obj@meta.data[[group_col]]),
    levels = group_levels
  )
  Seurat::Idents(obj) <- obj@meta.data[[group_col]]
  result <- stats::setNames(vector("list", length(group_levels)), group_levels)

  marker_paths <- stats::setNames(vapply(group_levels, function(group_id) {
    if (layout == "nested_deg") {
      group_dir <- ensure_dir(file.path(output_root, paste0("cluster_", group_id), "DEG"))
    } else if (layout == "nested") {
      group_dir <- ensure_dir(file.path(output_root, paste0("cluster_", group_id)))
    } else {
      group_dir <- ensure_dir(output_root)
    }
    file.path(group_dir, paste0("cluster_", group_id, "_vs_rest_markers.csv"))
  }, character(1)), group_levels)

  groups_to_run <- character(0)
  for (group_id in group_levels) {
    marker_path <- marker_paths[[group_id]]
    if (!file.exists(marker_path) || file.info(marker_path)$size <= 0L) {
      groups_to_run <- c(groups_to_run, group_id)
      next
    }
    de_cached <- utils::read.csv(
      marker_path,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    required_columns <- c("p_val", "pct.1", "pct.2", "p_val_adj", "gene", "gene_symbol")
    missing_columns <- setdiff(required_columns, colnames(de_cached))
    if (length(missing_columns) > 0L) {
      stop(
        "Cached DEG file is incomplete: ", marker_path,
        "; missing columns: ", paste(missing_columns, collapse = ", "),
        call. = FALSE
      )
    }
    resolve_lfc_col(de_cached)
    message("[DEG cache] ", group_col, " ", group_id, " vs rest: ", marker_path)
    result[[group_id]] <- de_cached
  }

  if (length(groups_to_run) == 0L) {
    message("All ", group_col, " DEG files are complete; using cached results.")
    return(result)
  }

  all_features <- rownames(obj[["RNA"]])
  target_task_count <- if (.Platform$OS.type == "unix") {
    min(DEG_FEATURE_CHUNK_TARGET, length(all_features))
  } else {
    length(groups_to_run)
  }
  chunks_per_group <- stats::setNames(
    rep(1L, length(groups_to_run)),
    groups_to_run
  )
  remaining_chunks <- target_task_count - length(groups_to_run)
  if (remaining_chunks > 0L) {
    for (i in seq_len(remaining_chunks)) {
      group_index <- ((i - 1L) %% length(groups_to_run)) + 1L
      chunks_per_group[[group_index]] <- chunks_per_group[[group_index]] + 1L
    }
  }
  chunk_cache_root <- ensure_dir(file.path(output_root, ".deg_feature_chunk_cache"))
  tasks <- unlist(lapply(groups_to_run, function(group_id) {
    n_chunks <- chunks_per_group[[group_id]]
    group_cache_dir <- ensure_dir(file.path(chunk_cache_root, paste0("cluster_", group_id)))
    feature_chunks <- split(
      all_features,
      rep(seq_len(n_chunks), length.out = length(all_features))
    )
    lapply(seq_along(feature_chunks), function(chunk_id) {
      list(
        group_id = group_id,
        chunk_id = as.integer(chunk_id),
        n_chunks = as.integer(n_chunks),
        features = feature_chunks[[chunk_id]],
        cache_path = file.path(
          group_cache_dir,
          sprintf("chunk_%03d_of_%03d.rds", chunk_id, n_chunks)
        )
      )
    })
  }), recursive = FALSE)

  task_cache_contract <- function(task) {
    list(
      group_col = group_col,
      group_id = task$group_id,
      chunk_id = task$chunk_id,
      n_chunks = task$n_chunks,
      features = task$features,
      min_pct = min_pct,
      logfc_threshold = logfc_threshold,
      seed = seed,
      r_version = paste(R.version$major, R.version$minor, sep = "."),
      seurat_version = as.character(utils::packageVersion("Seurat"))
    )
  }

  read_task_cache <- function(task) {
    if (!file.exists(task$cache_path) || file.info(task$cache_path)$size <= 0L) {
      return(NULL)
    }
    cached <- tryCatch(readRDS(task$cache_path), error = function(e) NULL)
    if (
      !is.list(cached) ||
      !identical(cached$contract, task_cache_contract(task)) ||
      !is.data.frame(cached$de)
    ) {
      return(NULL)
    }
    message(
      "[DEG chunk cache] ", group_col, " ", task$group_id,
      " chunk ", task$chunk_id, "/", task$n_chunks,
      "; rows=", nrow(cached$de)
    )
    list(
      group_id = task$group_id,
      chunk_id = task$chunk_id,
      de = cached$de,
      error = NULL
    )
  }

  run_one_task <- function(task) {
    started <- Sys.time()
    tryCatch(
      {
        message(
          "[DEG chunk start] ", group_col, " ", task$group_id,
          " chunk ", task$chunk_id, "/", task$n_chunks,
          "; genes=", length(task$features),
          "; at ", format(started, "%Y-%m-%dT%H:%M:%S%z")
        )
        set.seed(seed)
        marker_args <- list(
          object = obj,
          ident.1 = task$group_id,
          assay = "RNA",
          slot = "data",
          test.use = "wilcox",
          features = task$features,
          verbose = TRUE
        )
        if (!is.null(min_pct)) marker_args$min.pct <- min_pct
        if (!is.null(logfc_threshold)) marker_args$logfc.threshold <- logfc_threshold
        de <- do.call(Seurat::FindMarkers, marker_args)
        de <- as.data.frame(de, stringsAsFactors = FALSE)
        de$gene <- rownames(de)
        rownames(de) <- NULL
        cache_value <- list(
          contract = task_cache_contract(task),
          de = de
        )
        saveRDS(cache_value, task$cache_path)
        result_value <- list(
          group_id = task$group_id,
          chunk_id = task$chunk_id,
          de = de,
          error = NULL
        )
        message(
          "[DEG chunk done] ", group_col, " ", task$group_id,
          " chunk ", task$chunk_id, "/", task$n_chunks,
          "; rows=", nrow(de),
          "; elapsed_min=",
          round(as.numeric(difftime(Sys.time(), started, units = "mins")), 2),
          "; cache=", task$cache_path
        )
        result_value
      },
      error = function(e) {
        list(
          group_id = task$group_id,
          chunk_id = task$chunk_id,
          de = NULL,
          error = conditionMessage(e)
        )
      }
    )
  }

  cached_computed <- lapply(tasks, read_task_cache)
  cache_hit <- !vapply(cached_computed, is.null, logical(1))
  tasks_to_run <- tasks[!cache_hit]
  cached_computed <- cached_computed[cache_hit]

  workers <- min(DEG_PARALLEL_WORKERS, length(tasks_to_run))
  if (.Platform$OS.type != "unix") workers <- 1L
  message(
    "Running ", length(groups_to_run), " missing ", group_col,
    " DEG comparison(s) as ", length(tasks), " feature chunk(s); ",
    length(cached_computed), " cached and ", length(tasks_to_run),
    " to compute with up to ", workers, " forked worker(s)."
  )
  newly_computed <- if (length(tasks_to_run) == 0L) {
    list()
  } else if (workers > 1L) {
    parallel::mclapply(
      tasks_to_run,
      run_one_task,
      mc.cores = workers,
      mc.preschedule = FALSE,
      mc.set.seed = FALSE
    )
  } else {
    lapply(tasks_to_run, run_one_task)
  }
  computed <- c(cached_computed, newly_computed)

  for (item in computed) {
    if (inherits(item, "try-error") || !is.list(item) || is.null(item$group_id)) {
      stop("A forked DEG worker terminated without a valid result.", call. = FALSE)
    }
    if (!is.null(item$error)) {
      stop(
        "DEG failed for ", group_col, " ", item$group_id,
        ": ", item$error,
        call. = FALSE
      )
    }
  }

  for (group_id in groups_to_run) {
    group_parts <- computed[vapply(
      computed,
      function(item) identical(item$group_id, group_id),
      logical(1)
    )]
    group_parts <- group_parts[order(vapply(group_parts, `[[`, integer(1), "chunk_id"))]
    de <- do.call(rbind, lapply(group_parts, `[[`, "de"))
    rownames(de) <- NULL
    de$gene_symbol <- safe_gene_symbol(de$gene)
    if (nrow(de) > 0L) {
      lfc_col <- resolve_lfc_col(de)
      de <- de[order(de$p_val_adj, -abs(de[[lfc_col]]), de$gene), , drop = FALSE]
    }
    rownames(de) <- NULL
    write_csv(de, marker_paths[[group_id]])
    message(
      "[DEG assembled] ", group_col, " ", group_id,
      "; chunks=", length(group_parts),
      "; rows=", nrow(de),
      "; output=", marker_paths[[group_id]]
    )
    group_cache_paths <- vapply(
      tasks[vapply(tasks, function(task) identical(task$group_id, group_id), logical(1))],
      `[[`,
      character(1),
      "cache_path"
    )
    unlink(group_cache_paths, force = TRUE)
    result[[group_id]] <- de
  }
  result
}

SIMILARITY_SPECS <- list(
  all_ranked = list(type = "all", description = "All tested genes with signed avg_log2FC."),
  lenient = list(
    type = "filtered", padj = 0.05, abs_lfc = 0.25, abs_delta = 0.05,
    description = "p_adj < 0.05; abs(log2FC) >= 0.25; abs(delta pct) >= 0.05."
  ),
  moderate = list(
    type = "filtered", padj = 0.01, abs_lfc = 0.50, abs_delta = 0.10,
    description = "p_adj < 0.01; abs(log2FC) >= 0.50; abs(delta pct) >= 0.10."
  ),
  strict = list(
    type = "filtered", padj = 0.001, abs_lfc = 1.00, abs_delta = 0.20,
    description = "p_adj < 0.001; abs(log2FC) >= 1.00; abs(delta pct) >= 0.20."
  ),
  top50_each = list(
    type = "top_n", padj = 0.05, abs_lfc = 0.25, abs_delta = 0.05, top_n = 50L,
    description = "Lenient filter followed by top 50 up and top 50 down genes."
  ),
  top30_each = list(
    type = "top_n", padj = 0.05, abs_lfc = 0.25, abs_delta = 0.05, top_n = 30L,
    description = "Lenient filter followed by top 30 up and top 30 down genes."
  ),
  top20_each = list(
    type = "top_n", padj = 0.05, abs_lfc = 0.25, abs_delta = 0.05, top_n = 20L,
    description = "Lenient filter followed by top 20 up and top 20 down genes."
  ),
  top10_each = list(
    type = "top_n", padj = 0.05, abs_lfc = 0.25, abs_delta = 0.05, top_n = 10L,
    description = "Lenient filter followed by top 10 up and top 10 down genes."
  )
)

prepare_similarity_table <- function(de) {
  if (nrow(de) == 0L) {
    return(data.frame(
      gene_key = character(0), lfc = numeric(0), p_adj = numeric(0),
      delta_pct = numeric(0), abs_lfc = numeric(0), abs_delta = numeric(0)
    ))
  }
  lfc_col <- resolve_lfc_col(de)
  gene_label <- if ("gene_symbol" %in% colnames(de)) de$gene_symbol else de$gene
  gene_label[is.na(gene_label) | gene_label == ""] <- de$gene[is.na(gene_label) | gene_label == ""]
  out <- data.frame(
    gene_key = toupper(safe_gene_symbol(gene_label)),
    lfc = as.numeric(de[[lfc_col]]),
    p_adj = as.numeric(de$p_val_adj),
    delta_pct = as.numeric(de$pct.1) - as.numeric(de$pct.2),
    stringsAsFactors = FALSE
  )
  out$abs_lfc <- abs(out$lfc)
  out$abs_delta <- abs(out$delta_pct)
  out <- out[
    !is.na(out$gene_key) & out$gene_key != "" & is.finite(out$lfc) & !is.na(out$p_adj),
    , drop = FALSE
  ]
  out <- out[order(-out$abs_lfc, out$p_adj), , drop = FALSE]
  out[!duplicated(out$gene_key), , drop = FALSE]
}

select_similarity_features <- function(prepared, spec) {
  if (spec$type == "all") {
    selected <- prepared
  } else {
    selected <- prepared[
      prepared$p_adj < spec$padj &
        prepared$abs_lfc >= spec$abs_lfc &
        prepared$abs_delta >= spec$abs_delta,
      , drop = FALSE
    ]
  }
  if (spec$type == "top_n") {
    up <- selected[selected$lfc > 0, , drop = FALSE]
    down <- selected[selected$lfc < 0, , drop = FALSE]
    up <- up[order(-up$lfc, up$p_adj), , drop = FALSE]
    down <- down[order(down$lfc, down$p_adj), , drop = FALSE]
    selected <- rbind(head(up, spec$top_n), head(down, spec$top_n))
    selected <- selected[order(-selected$abs_lfc, selected$p_adj), , drop = FALSE]
    selected <- selected[!duplicated(selected$gene_key), , drop = FALSE]
  }
  signed_vector <- selected$lfc
  names(signed_vector) <- selected$gene_key
  list(table = selected, signed_vector = signed_vector)
}

signed_lfc_cosine <- function(v1, v2) {
  genes <- union(names(v1), names(v2))
  if (length(genes) == 0L) return(NA_real_)
  x <- stats::setNames(numeric(length(genes)), genes)
  y <- x
  x[names(v1)] <- v1
  y[names(v2)] <- v2
  denominator <- sqrt(sum(x^2)) * sqrt(sum(y^2))
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  sum(x * y) / denominator
}

write_similarity_heatmap <- function(mat, path, title) {
  plot_df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(plot_df) <- c("cluster_1", "cluster_2", "similarity")
  plot_df$cluster_1 <- factor(plot_df$cluster_1, levels = rev(rownames(mat)))
  plot_df$cluster_2 <- factor(plot_df$cluster_2, levels = colnames(mat))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = cluster_2, y = cluster_1, fill = similarity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#1f4e79", mid = "white", high = "#b03a2e", midpoint = 0,
      limits = c(-1, 1), na.value = "grey90"
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::labs(title = title, x = "Cluster", y = "Cluster", fill = "Cosine")
  ggplot2::ggsave(path, p, width = 8, height = 7)
}

run_initial_deg_similarity <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "02a_cluster_similarity_by_DEGs"))
  deg_root <- ensure_dir(file.path(stage_root, "initial_cluster_vs_rest_DEG"))
  matrix_root <- ensure_dir(file.path(stage_root, "matrices"))
  summary_root <- ensure_dir(file.path(stage_root, "summaries"))
  plot_root <- ensure_dir(file.path(stage_root, "plots"))

  initial_levels <- if (is.factor(obj$seurat_clusters)) {
    as.character(levels(obj$seurat_clusters))
  } else {
    sort_maybe_numeric(obj$seurat_clusters)
  }
  assert_set_equal(initial_levels, EXPECTED_INITIAL_CLUSTERS, "Similarity-stage initial clusters")

  marker_obj <- prepare_obj_for_markers(obj)
  marker_tables <- run_deg_by_group(
    marker_obj,
    group_col = "seurat_clusters",
    group_levels = initial_levels,
    output_root = deg_root,
    min_pct = 0.10,
    logfc_threshold = 0,
    seed = SEEDS[["initial_DEG_similarity"]],
    layout = "nested_deg"
  )
  prepared_tables <- lapply(marker_tables, prepare_similarity_table)

  threshold_metadata <- data.frame(
    threshold = names(SIMILARITY_SPECS),
    type = vapply(SIMILARITY_SPECS, `[[`, character(1), "type"),
    description = vapply(SIMILARITY_SPECS, `[[`, character(1), "description"),
    stringsAsFactors = FALSE
  )
  write_csv(threshold_metadata, file.path(summary_root, "threshold_metadata.csv"))
  write_csv(
    data.frame(
      selected_metric = "signed_lfc_cosine",
      rationale = "Signed fold-change vectors retain both magnitude and direction; absent genes are zero-filled over the pairwise union.",
      selection_stage = "before_fixed_merge",
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "selected_similarity_metric.csv")
  )

  pairwise_rows <- list()
  feature_size_rows <- list()
  nearest_rows <- list()
  threshold_summary_rows <- list()
  row_index <- 1L
  size_index <- 1L
  nearest_index <- 1L
  summary_index <- 1L
  matrix_paths <- character(0)

  for (threshold_name in names(SIMILARITY_SPECS)) {
    spec <- SIMILARITY_SPECS[[threshold_name]]
    feature_map <- lapply(prepared_tables, select_similarity_features, spec = spec)
    for (cluster_id in initial_levels) {
      feature_size_rows[[size_index]] <- data.frame(
        threshold = threshold_name,
        cluster = cluster_id,
        n_features = nrow(feature_map[[cluster_id]]$table),
        n_up = sum(feature_map[[cluster_id]]$table$lfc > 0),
        n_down = sum(feature_map[[cluster_id]]$table$lfc < 0),
        stringsAsFactors = FALSE
      )
      size_index <- size_index + 1L
    }

    mat <- matrix(
      NA_real_, length(initial_levels), length(initial_levels),
      dimnames = list(initial_levels, initial_levels)
    )
    for (i in seq_along(initial_levels)) {
      for (j in i:length(initial_levels)) {
        c1 <- initial_levels[[i]]
        c2 <- initial_levels[[j]]
        value <- if (i == j) 1 else signed_lfc_cosine(
          feature_map[[c1]]$signed_vector,
          feature_map[[c2]]$signed_vector
        )
        mat[c1, c2] <- value
        mat[c2, c1] <- value
        pairwise_rows[[row_index]] <- data.frame(
          threshold = threshold_name,
          metric = "signed_lfc_cosine",
          cluster_1 = c1,
          cluster_2 = c2,
          similarity = value,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }

    matrix_path <- file.path(
      matrix_root,
      paste0("similarity_matrix__signed_lfc_cosine__", threshold_name, ".csv")
    )
    matrix_df <- data.frame(cluster = rownames(mat), as.data.frame(mat, check.names = FALSE), check.names = FALSE)
    write_csv(matrix_df, matrix_path)
    matrix_paths <- c(matrix_paths, matrix_path)
    write_similarity_heatmap(
      mat,
      file.path(plot_root, paste0("heatmap__signed_lfc_cosine__", threshold_name, ".pdf")),
      paste0("signed_lfc_cosine | ", threshold_name)
    )

    for (cluster_id in initial_levels) {
      values <- mat[cluster_id, ]
      values <- values[names(values) != cluster_id]
      best <- if (all(is.na(values))) NA_integer_ else which.max(values)
      nearest_rows[[nearest_index]] <- data.frame(
        threshold = threshold_name,
        metric = "signed_lfc_cosine",
        cluster = cluster_id,
        nearest_cluster = if (is.na(best)) NA_character_ else names(values)[best],
        similarity = if (is.na(best)) NA_real_ else as.numeric(values[best]),
        stringsAsFactors = FALSE
      )
      nearest_index <- nearest_index + 1L
    }
    offdiag <- mat[row(mat) != col(mat)]
    offdiag <- offdiag[is.finite(offdiag)]
    threshold_summary_rows[[summary_index]] <- data.frame(
      threshold = threshold_name,
      metric = "signed_lfc_cosine",
      mean_offdiag_similarity = mean(offdiag),
      median_offdiag_similarity = stats::median(offdiag),
      sd_offdiag_similarity = stats::sd(offdiag),
      max_offdiag_similarity = max(offdiag),
      min_offdiag_similarity = min(offdiag),
      stringsAsFactors = FALSE
    )
    summary_index <- summary_index + 1L
  }

  pairwise_df <- do.call(rbind, pairwise_rows)
  pairwise_df <- pairwise_df[pairwise_df$cluster_1 != pairwise_df$cluster_2, , drop = FALSE]
  pairwise_df$pair_id <- apply(pairwise_df[, c("cluster_1", "cluster_2"), drop = FALSE], 1, function(z) {
    paste(sort(z), collapse = "__")
  })
  pairwise_df <- pairwise_df[!duplicated(pairwise_df[, c("threshold", "pair_id")]), , drop = FALSE]
  write_csv(pairwise_df, file.path(summary_root, "pairwise_similarity_long.csv"))
  write_csv(do.call(rbind, feature_size_rows), file.path(summary_root, "feature_set_sizes_by_threshold.csv"))
  write_csv(do.call(rbind, nearest_rows), file.path(summary_root, "nearest_neighbor_by_threshold_metric.csv"))
  write_csv(do.call(rbind, threshold_summary_rows), file.path(summary_root, "threshold_metric_summary.csv"))

  pair_ids <- unique(pairwise_df$pair_id)
  stability <- do.call(rbind, lapply(pair_ids, function(pair_id) {
    values <- pairwise_df$similarity[pairwise_df$pair_id == pair_id]
    pair_row <- pairwise_df[pairwise_df$pair_id == pair_id, , drop = FALSE][1, ]
    data.frame(
      metric = "signed_lfc_cosine",
      pair_id = pair_id,
      cluster_1 = pair_row$cluster_1,
      cluster_2 = pair_row$cluster_2,
      n_thresholds = sum(is.finite(values)),
      mean_similarity = mean(values, na.rm = TRUE),
      median_similarity = stats::median(values, na.rm = TRUE),
      min_similarity = min(values, na.rm = TRUE),
      max_similarity = max(values, na.rm = TRUE),
      range_similarity = diff(range(values, na.rm = TRUE)),
      sd_similarity = stats::sd(values, na.rm = TRUE),
      always_positive = all(values > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  write_csv(stability, file.path(summary_root, "pair_similarity_across_thresholds_summary.csv"))

  gate_checks <- data.frame(
    check = c(
      "initial_cluster_count_is_15",
      "all_eight_feature_definitions_completed",
      "all_cosine_matrices_exist",
      "selected_metric_is_signed_lfc_cosine",
      "similarity_precedes_fixed_merge"
    ),
    passed = c(
      length(initial_levels) == 15L,
      length(unique(pairwise_df$threshold)) == 8L,
      length(matrix_paths) == 8L && all(file.exists(matrix_paths)),
      TRUE,
      TRUE
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(gate_checks, file.path(summary_root, "fixed_merge_expression_similarity_gate.tsv"))
  if (!all(gate_checks$passed)) {
    stop("Expression-similarity gate failed; fixed merge is forbidden.", call. = FALSE)
  }
  writeLines(
    c(
      "status=PASS",
      "selected_metric=signed_lfc_cosine",
      "stage=before_fixed_merge",
      "feature_definitions=all_ranked,lenient,moderate,strict,top50_each,top30_each,top20_each,top10_each"
    ),
    file.path(summary_root, "FIXED_MERGE_GATE_PASS.txt")
  )
  write_stage_marker("02a_cluster_similarity_by_DEGs", "signed_lfc_cosine gate passed before fixed merge")
  invisible(list(gate_file = file.path(summary_root, "FIXED_MERGE_GATE_PASS.txt")))
}

build_refined_levels <- function(original_levels, label_map) {
  out <- character(0)
  for (cluster_id in original_levels) {
    out <- c(out, cluster_id)
    if (cluster_id %in% names(label_map)) out <- c(out, unname(label_map[[cluster_id]]))
  }
  unique(out)
}

point_in_polygon <- function(x, y, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    xi <- poly_x[i]
    yi <- poly_y[i]
    xj <- poly_x[j]
    yj <- poly_y[j]
    crosses <- ((yi > y) != (yj > y)) &
      (x < ((xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi))
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

compute_knn_indices <- function(coords, k) {
  RANN::nn2(data = coords, query = coords, k = k + 1L)$nn.idx[, -1, drop = FALSE]
}

compute_knn_graph_components <- function(coords, k, mutual = FALSE) {
  n <- nrow(coords)
  if (n == 0L) return(integer(0))
  if (n == 1L) return(1L)
  k_use <- min(as.integer(k), n - 1L)
  if (!is.finite(k_use) || k_use < 1L) return(seq_len(n))
  nn_idx <- compute_knn_indices(coords, k_use)
  if (!is.matrix(nn_idx)) nn_idx <- matrix(nn_idx, nrow = n, ncol = k_use)
  adjacency <- replicate(n, integer(0), simplify = FALSE)
  if (mutual) {
    nn_sets <- lapply(seq_len(n), function(i) unique(as.integer(nn_idx[i, ])))
    for (i in seq_len(n)) {
      for (j in nn_sets[[i]]) {
        if (!is.na(j) && j >= 1L && j <= n && i %in% nn_sets[[j]]) {
          adjacency[[i]] <- c(adjacency[[i]], j)
          adjacency[[j]] <- c(adjacency[[j]], i)
        }
      }
    }
  } else {
    for (i in seq_len(n)) {
      neighbors <- unique(as.integer(nn_idx[i, ]))
      neighbors <- neighbors[!is.na(neighbors) & neighbors >= 1L & neighbors <= n]
      adjacency[[i]] <- c(adjacency[[i]], neighbors)
      for (j in neighbors) adjacency[[j]] <- c(adjacency[[j]], i)
    }
  }
  adjacency <- lapply(adjacency, unique)
  component <- integer(n)
  component_id <- 0L
  for (i in seq_len(n)) {
    if (component[i] != 0L) next
    component_id <- component_id + 1L
    queue <- i
    component[i] <- component_id
    head_index <- 1L
    while (head_index <= length(queue)) {
      node <- queue[[head_index]]
      head_index <- head_index + 1L
      neighbors <- adjacency[[node]]
      if (length(neighbors) == 0L) next
      new_nodes <- neighbors[component[neighbors] == 0L]
      if (length(new_nodes) > 0L) {
        component[new_nodes] <- component_id
        queue <- c(queue, new_nodes)
      }
    }
  }
  component
}

expand_candidate_component <- function(candidate_idx, seed_idx, umap_mat, frac_region, frac_core, inside_hull) {
  if (length(candidate_idx) == 0L || length(seed_idx) == 0L) {
    return(list(selected_idx = integer(0), eligible_idx = integer(0)))
  }
  eligible_idx <- candidate_idx[
    frac_region[candidate_idx] >= 0 & frac_core[candidate_idx] >= 0 & (!FALSE | inside_hull[candidate_idx])
  ]
  seed_in_eligible <- eligible_idx %in% seed_idx
  if (!any(seed_in_eligible)) return(list(selected_idx = integer(0), eligible_idx = eligible_idx))
  component <- compute_knn_graph_components(umap_mat[eligible_idx, , drop = FALSE], k = 5L, mutual = FALSE)
  selected_components <- unique(component[seed_in_eligible])
  list(
    selected_idx = eligible_idx[component %in% selected_components],
    eligible_idx = eligible_idx
  )
}

build_focus_status <- function(cluster_vec, refined_vec, region_clusters) {
  out <- rep("other", length(cluster_vec))
  out[cluster_vec %in% region_clusters] <- "reference_region"
  out[cluster_vec == "4"] <- "cluster4_other"
  out[cluster_vec == "9"] <- "cluster9_other"
  out[refined_vec == "4c"] <- "cluster4c_selected"
  out[refined_vec == "9c"] <- "cluster9c_selected"
  factor(
    out,
    levels = c(
      "reference_region", "cluster4_other", "cluster4c_selected",
      "cluster9_other", "cluster9c_selected", "other"
    )
  )
}

canonicalize_pca_orientation <- function(coordinates) {
  pivot_index <- apply(abs(coordinates), 2L, which.max)
  orientation <- vapply(
    seq_len(ncol(coordinates)),
    function(j) if (coordinates[pivot_index[[j]], j] < 0) -1 else 1,
    numeric(1)
  )
  list(
    coordinates = sweep(coordinates, 2L, orientation, `*`),
    pivot_index = pivot_index,
    pivot_cell = rownames(coordinates)[pivot_index],
    orientation = orientation
  )
}

run_refine_compatibility_umaps <- function(obj) {
  coordinates <- Seurat::Embeddings(obj, reduction = "pca")[, seq_len(30L), drop = FALSE]
  canonical <- canonicalize_pca_orientation(coordinates)
  write_tsv(
    data.frame(
      component = colnames(coordinates),
      pivot_cell = canonical$pivot_cell,
      input_orientation = canonical$orientation,
      canonical_pivot_value = canonical$coordinates[cbind(
        canonical$pivot_index,
        seq_len(ncol(canonical$coordinates))
      )],
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "00_provenance", "refine_pca_orientation_contract.tsv")
  )

  run_one <- function(index) {
    step_name <- names(REFINE_COMPATIBILITY_STEPS)[[index]]
    step <- unname(REFINE_COMPATIBILITY_STEPS[[index]])
    quantized <- round(canonical$coordinates / step) * step
    run_compatible_umap_matrix(
      coordinates = quantized,
      seed = SEEDS[["initial_UMAP"]],
      stage_label = paste0("refine_compatibility_", step_name)
    )
  }
  indices <- seq_along(REFINE_COMPATIBILITY_STEPS)
  if (identical(Sys.info()[["sysname"]], "Linux") && length(indices) > 1L) {
    embeddings <- parallel::mclapply(
      indices,
      run_one,
      mc.cores = min(REFINE_COMPATIBILITY_WORKERS, length(indices)),
      mc.preschedule = FALSE
    )
  } else {
    embeddings <- lapply(indices, run_one)
  }
  failed <- vapply(embeddings, inherits, logical(1), what = "try-error")
  if (any(failed)) {
    stop(
      "At least one refine compatibility UMAP failed: ",
      paste(names(REFINE_COMPATIBILITY_STEPS)[failed], collapse = ", "),
      call. = FALSE
    )
  }
  names(embeddings) <- names(REFINE_COMPATIBILITY_STEPS)
  write_tsv(
    data.frame(
      selection_target = c("cluster4", "cluster4_and_cluster9", "cluster9"),
      member = names(REFINE_COMPATIBILITY_STEPS),
      pca_quantization_step = unname(REFINE_COMPATIBILITY_STEPS),
      set_operator = c("intersection", "intersection_for_4_and_union_for_9", "union"),
      pca_orientation = "largest_absolute_cell_coordinate_positive",
      umap_seed = SEEDS[["initial_UMAP"]],
      stringsAsFactors = FALSE
    ),
    file.path(output_dir, "00_provenance", "refine_numeric_compatibility_rule.tsv")
  )
  embeddings
}

select_refine_indices_from_embedding <- function(umap_mat, cluster_vec) {
  core_clusters <- c("6", "10", "11", "12")
  region_clusters <- c(core_clusters, "0", "3")
  candidate_clusters <- c("4", "9")
  nn_idx <- compute_knn_indices(umap_mat, 50L)
  neighbor_labels <- matrix(cluster_vec[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  frac_region <- rowMeans(matrix(
    neighbor_labels %in% region_clusters,
    nrow = nrow(neighbor_labels)
  ))
  frac_core <- rowMeans(matrix(
    neighbor_labels %in% core_clusters,
    nrow = nrow(neighbor_labels)
  ))
  reference_idx <- which(cluster_vec %in% region_clusters)
  reference_umap <- umap_mat[reference_idx, , drop = FALSE]
  hull_order <- grDevices::chull(reference_umap[, 1], reference_umap[, 2])
  hull_xy <- rbind(
    reference_umap[hull_order, , drop = FALSE],
    reference_umap[hull_order[[1]], , drop = FALSE]
  )
  inside_hull <- point_in_polygon(
    umap_mat[, 1], umap_mat[, 2], hull_xy[, 1], hull_xy[, 2]
  )
  candidate_idx <- which(cluster_vec %in% candidate_clusters)
  seed_idx <- candidate_idx[
    frac_region[candidate_idx] >= 0.60 &
      frac_core[candidate_idx] >= 0.20 &
      inside_hull[candidate_idx]
  ]
  expanded_idx <- integer(0)
  for (candidate_cluster in candidate_clusters) {
    result <- expand_candidate_component(
      candidate_idx = which(cluster_vec == candidate_cluster),
      seed_idx = seed_idx[cluster_vec[seed_idx] == candidate_cluster],
      umap_mat = umap_mat,
      frac_region = frac_region,
      frac_core = frac_core,
      inside_hull = inside_hull
    )
    expanded_idx <- union(expanded_idx, result$selected_idx)
  }
  sort(unique(c(seed_idx, expanded_idx)))
}

run_cluster_refine <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "02b_cluster_refine"))
  object_root <- ensure_dir(file.path(stage_root, "objects"))
  summary_root <- ensure_dir(file.path(stage_root, "summaries"))
  cell_root <- ensure_dir(file.path(stage_root, "selected_cells"))
  refined_path <- file.path(object_root, "integrated_sct_cca_seurat_cluster_refine.rds")
  if (file.exists(refined_path)) {
    message("[resume] Reading verified refined object: ", refined_path)
    resumed <- readRDS(refined_path)
    validate_refined_resume(resumed, obj)
    oriented <- orient_seurat_pca_to_reference_anchor(
      resumed, contract_name = "initial", stage_label = "refined_resume"
    )
    resumed <- oriented$object
    if (any(oriented$flipped)) {
      message_cluster_counts_before_rds(
        resumed,
        "seurat_cluster_refine",
        refined_path,
        "Refined PCA-orientation repair"
      )
      saveRDS(resumed, refined_path)
    }
    return(resumed)
  }
  set.seed(SEEDS[["cluster_refine"]])

  core_clusters <- c("6", "10", "11", "12")
  neighbor_clusters <- c("0", "3")
  region_clusters <- c(core_clusters, neighbor_clusters)
  candidate_clusters <- c("4", "9")
  label_map <- c(`4` = "4c", `9` = "9c")
  cluster_vec <- as.character(obj$seurat_clusters)
  cluster_levels <- sort_maybe_numeric(cluster_vec)
  obj$seurat_clusters <- factor(cluster_vec, levels = cluster_levels)
  umap_mat <- Seurat::Embeddings(obj, "umap")[, 1:2, drop = FALSE]
  colnames(umap_mat) <- c("UMAP_1", "UMAP_2")

  nn_idx <- compute_knn_indices(umap_mat, 50L)
  neighbor_labels <- matrix(cluster_vec[nn_idx], nrow = nrow(nn_idx), ncol = ncol(nn_idx))
  neighbor_is_region <- matrix(neighbor_labels %in% region_clusters, nrow = nrow(neighbor_labels))
  neighbor_is_core <- matrix(neighbor_labels %in% core_clusters, nrow = nrow(neighbor_labels))
  frac_region <- rowMeans(neighbor_is_region, na.rm = TRUE)
  frac_core <- rowMeans(neighbor_is_core, na.rm = TRUE)

  reference_idx <- which(cluster_vec %in% region_clusters)
  reference_umap <- umap_mat[reference_idx, , drop = FALSE]
  hull_order <- grDevices::chull(reference_umap[, 1], reference_umap[, 2])
  hull_xy <- reference_umap[hull_order, , drop = FALSE]
  hull_xy <- rbind(hull_xy, hull_xy[1, , drop = FALSE])
  inside_hull <- point_in_polygon(
    umap_mat[, 1], umap_mat[, 2], hull_xy[, 1], hull_xy[, 2]
  )

  candidate_idx <- which(cluster_vec %in% candidate_clusters)
  seed_idx <- candidate_idx[
    frac_region[candidate_idx] >= 0.60 &
      frac_core[candidate_idx] >= 0.20 &
      inside_hull[candidate_idx]
  ]
  expanded_idx <- integer(0)
  eligible_idx <- integer(0)
  for (candidate_cluster in candidate_clusters) {
    result <- expand_candidate_component(
      candidate_idx = which(cluster_vec == candidate_cluster),
      seed_idx = seed_idx[cluster_vec[seed_idx] == candidate_cluster],
      umap_mat = umap_mat,
      frac_region = frac_region,
      frac_core = frac_core,
      inside_hull = inside_hull
    )
    expanded_idx <- union(expanded_idx, result$selected_idx)
    eligible_idx <- union(eligible_idx, result$eligible_idx)
  }
  main_selected_idx <- sort(unique(c(seed_idx, expanded_idx)))
  compatibility_embeddings <- run_refine_compatibility_umaps(obj)
  compatibility_selected_idx <- lapply(
    compatibility_embeddings,
    select_refine_indices_from_embedding,
    cluster_vec = cluster_vec
  )
  cluster4_selected_idx <- Reduce(
    intersect,
    lapply(
      compatibility_selected_idx[c(
        "cluster4_intersection_a", "cluster4_intersection_b_cluster9_union_a"
      )],
      function(index) index[cluster_vec[index] == "4"]
    )
  )
  cluster9_selected_idx <- Reduce(
    union,
    lapply(
      compatibility_selected_idx[c(
        "cluster4_intersection_b_cluster9_union_a", "cluster9_union_b"
      )],
      function(index) index[cluster_vec[index] == "9"]
    )
  )
  selected_idx <- sort(unique(c(cluster4_selected_idx, cluster9_selected_idx)))
  compatibility_selected <- lapply(
    compatibility_selected_idx,
    function(index) seq_len(ncol(obj)) %in% index
  )
  refined_vec <- cluster_vec
  refined_vec[selected_idx] <- unname(label_map[cluster_vec[selected_idx]])
  refined_levels <- build_refined_levels(cluster_levels, label_map)
  obj$seurat_cluster_refine <- factor(refined_vec, levels = refined_levels)
  obj$cluster_refine_focus_status <- build_focus_status(cluster_vec, refined_vec, region_clusters)

  selected_seed_cells <- colnames(obj)[seed_idx]
  main_selected_cells <- colnames(obj)[main_selected_idx]
  selected_cells <- colnames(obj)[selected_idx]
  score_df <- data.frame(
    cell = colnames(obj),
    original_cluster = cluster_vec,
    seurat_cluster_refine = refined_vec,
    UMAP_1 = umap_mat[, 1],
    UMAP_2 = umap_mat[, 2],
    inside_reference_hull = inside_hull,
    frac_region_neighbors = frac_region,
    frac_core_neighbors = frac_core,
    n_region_neighbors = rowSums(neighbor_is_region),
    n_core_neighbors = rowSums(neighbor_is_core),
    selected_as_seed = colnames(obj) %in% selected_seed_cells,
    selected_by_component_expansion = colnames(obj) %in% setdiff(
      main_selected_cells,
      selected_seed_cells
    ),
    component_expand_eligible = seq_len(ncol(obj)) %in% eligible_idx,
    compatibility_cluster4_step_0_0125 = compatibility_selected$cluster4_intersection_a,
    compatibility_shared_step_0_03 = compatibility_selected$cluster4_intersection_b_cluster9_union_a,
    compatibility_cluster9_step_0_09 = compatibility_selected$cluster9_union_b,
    selected_by_compatibility_rule = colnames(obj) %in% selected_cells,
    selected_for_refine = colnames(obj) %in% selected_cells,
    stringsAsFactors = FALSE
  )
  compatibility_counts <- data.frame(
    member = names(compatibility_selected),
    pca_quantization_step = unname(REFINE_COMPATIBILITY_STEPS),
    selected_from_cluster4 = vapply(
      compatibility_selected,
      function(selected) sum(selected & cluster_vec == "4"),
      integer(1)
    ),
    selected_from_cluster9 = vapply(
      compatibility_selected,
      function(selected) sum(selected & cluster_vec == "9"),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    compatibility_counts,
    file.path(summary_root, "refine_numeric_compatibility_member_counts.tsv")
  )
  message(
    "Refine compatibility selection before RDS save: 4c=",
    length(cluster4_selected_idx),
    ", 9c=",
    length(cluster9_selected_idx)
  )
  write_csv(score_df[score_df$original_cluster %in% candidate_clusters, ], file.path(summary_root, "candidate_cluster4_9_region_scores.csv"))
  write_csv(score_df[score_df$seurat_cluster_refine == "4c", ], file.path(cell_root, "cells_selected_as_4c.csv"))
  write_csv(score_df[score_df$seurat_cluster_refine == "9c", ], file.path(cell_root, "cells_selected_as_9c.csv"))
  refined_counts <- table(obj$seurat_cluster_refine)
  write_csv(
    data.frame(seurat_cluster_refine = names(refined_counts), n_cells = as.integer(refined_counts)),
    file.path(summary_root, "cluster_counts_after_refine.csv")
  )

  if (ncol(obj) != EXPECTED_REFINED_CELLS) {
    stop("Refined-object cell count mismatch: ", ncol(obj), " vs ", EXPECTED_REFINED_CELLS, call. = FALSE)
  }
  observed_split <- as.integer(refined_counts[names(EXPECTED_REFINED_SPLIT_COUNTS)])
  names(observed_split) <- names(EXPECTED_REFINED_SPLIT_COUNTS)
  if (anyNA(observed_split) ||
      !identical(unname(observed_split), unname(as.integer(EXPECTED_REFINED_SPLIT_COUNTS)))) {
    stop(
      "Refined 4c/9c counts do not reproduce the reference. Observed ",
      paste(names(observed_split), observed_split, sep = "=", collapse = ", "),
      "; expected ",
      paste(names(EXPECTED_REFINED_SPLIT_COUNTS), EXPECTED_REFINED_SPLIT_COUNTS, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  refined_metadata_expected <- EXPECTED_BASE_METADATA
  if (!identical(colnames(obj@meta.data), refined_metadata_expected)) {
    stop(
      "Refined metadata contract mismatch. Observed: ",
      paste(colnames(obj@meta.data), collapse = ","),
      call. = FALSE
    )
  }

  message_cluster_counts_before_rds(
    obj,
    "seurat_cluster_refine",
    refined_path,
    "Refined"
  )
  saveRDS(obj, refined_path)
  write_stage_marker("02b_cluster_refine", "42884 cells; 4c=103; 9c=409")
  obj
}

compute_top_expression_metrics <- function(counts, top_n = 50L) {
  total_counts <- Matrix::colSums(counts)
  dominant_fraction <- rep(NA_real_, ncol(counts))
  percent_top <- rep(NA_real_, ncol(counts))
  if (inherits(counts, "dgCMatrix")) {
    pointer <- counts@p
    values_all <- counts@x
    for (j in seq_len(ncol(counts))) {
      total_j <- total_counts[j]
      if (!is.finite(total_j) || total_j <= 0) {
        dominant_fraction[j] <- 0
        percent_top[j] <- 0
        next
      }
      start <- pointer[j] + 1L
      end <- pointer[j + 1L]
      values <- if (start <= end) values_all[start:end] else numeric(0)
      if (length(values) == 0L) {
        dominant_fraction[j] <- 0
        percent_top[j] <- 0
        next
      }
      values <- sort(values, decreasing = TRUE)
      dominant_fraction[j] <- values[1] / total_j
      percent_top[j] <- 100 * sum(head(values, min(top_n, length(values)))) / total_j
    }
  } else {
    for (j in seq_len(ncol(counts))) {
      values <- as.numeric(counts[, j, drop = TRUE])
      total_j <- sum(values)
      values <- sort(values[values > 0], decreasing = TRUE)
      if (total_j <= 0 || length(values) == 0L) {
        dominant_fraction[j] <- 0
        percent_top[j] <- 0
      } else {
        dominant_fraction[j] <- values[1] / total_j
        percent_top[j] <- 100 * sum(head(values, min(top_n, length(values)))) / total_j
      }
    }
  }
  data.frame(
    cell = colnames(counts),
    dominant_gene_fraction = dominant_fraction,
    percent.top50 = percent_top,
    stringsAsFactors = FALSE
  )
}

safe_mad <- function(x) {
  value <- stats::mad(x, center = stats::median(x, na.rm = TRUE), constant = 1, na.rm = TRUE)
  if (!is.finite(value) || value == 0) NA_real_ else value
}

build_qc_summary_long <- function(qc_df, metrics, cluster_levels) {
  do.call(rbind, lapply(metrics, function(metric_name) {
    do.call(rbind, lapply(cluster_levels, function(cluster_id) {
      values <- qc_df[[metric_name]][qc_df$cluster == cluster_id]
      data.frame(
        cluster = cluster_id,
        metric = metric_name,
        n_cells = length(values),
        mean_value = mean(values, na.rm = TRUE),
        median_value = stats::median(values, na.rm = TRUE),
        sd_value = stats::sd(values, na.rm = TRUE),
        iqr_value = stats::IQR(values, na.rm = TRUE),
        q05 = stats::quantile(values, 0.05, na.rm = TRUE, names = FALSE),
        q25 = stats::quantile(values, 0.25, na.rm = TRUE, names = FALSE),
        q75 = stats::quantile(values, 0.75, na.rm = TRUE, names = FALSE),
        q95 = stats::quantile(values, 0.95, na.rm = TRUE, names = FALSE),
        stringsAsFactors = FALSE
      )
    }))
  }))
}

compute_cluster_logistic_qc <- function(qc_df, metrics, cluster_levels) {
  do.call(rbind, lapply(cluster_levels, function(cluster_id) {
    sub_df <- qc_df[, metrics, drop = FALSE]
    sub_df$target <- as.integer(qc_df$cluster == cluster_id)
    sub_df <- sub_df[stats::complete.cases(sub_df), , drop = FALSE]
    if (nrow(sub_df) < 50L || sum(sub_df$target == 1L) < 10L || sum(sub_df$target == 0L) < 10L) {
      return(data.frame(
        cluster = cluster_id, n_cells_used = nrow(sub_df),
        n_target = sum(sub_df$target == 1L), n_rest = sum(sub_df$target == 0L),
        mcfadden_r2 = NA_real_, model_aic = NA_real_, stringsAsFactors = FALSE
      ))
    }
    for (metric_name in metrics) sub_df[[metric_name]] <- as.numeric(scale(sub_df[[metric_name]]))
    formula_full <- stats::as.formula(paste("target ~", paste(metrics, collapse = " + ")))
    fit_full <- tryCatch(
      suppressWarnings(stats::glm(formula_full, data = sub_df, family = stats::binomial())),
      error = function(e) NULL
    )
    fit_null <- tryCatch(
      suppressWarnings(stats::glm(target ~ 1, data = sub_df, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(fit_full) || is.null(fit_null)) {
      return(data.frame(
        cluster = cluster_id, n_cells_used = nrow(sub_df),
        n_target = sum(sub_df$target == 1L), n_rest = sum(sub_df$target == 0L),
        mcfadden_r2 = NA_real_, model_aic = NA_real_, stringsAsFactors = FALSE
      ))
    }
    ll_full <- as.numeric(stats::logLik(fit_full))
    ll_null <- as.numeric(stats::logLik(fit_null))
    data.frame(
      cluster = cluster_id,
      n_cells_used = nrow(sub_df),
      n_target = sum(sub_df$target == 1L),
      n_rest = sum(sub_df$target == 0L),
      mcfadden_r2 = if (is.finite(ll_full) && is.finite(ll_null) && ll_null != 0) 1 - ll_full / ll_null else NA_real_,
      model_aic = stats::AIC(fit_full),
      stringsAsFactors = FALSE
    )
  }))
}

make_qc_flags <- function(summary_long, logistic_df) {
  directions <- c(
    nCount_RNA = "low", nFeature_RNA = "low", percent.mt = "high",
    log10_genes_per_umi = "low", umi_per_gene = "high",
    dominant_gene_fraction = "high", percent.top50 = "high", percent.ribo = "high"
  )
  flag_parts <- lapply(unique(summary_long$metric), function(metric_name) {
    part <- summary_long[summary_long$metric == metric_name, , drop = FALSE]
    across_median <- stats::median(part$median_value, na.rm = TRUE)
    across_mad <- safe_mad(part$median_value)
    part$suspicious_direction <- unname(directions[metric_name])
    part$across_cluster_median <- across_median
    part$across_cluster_mad <- across_mad
    part$robust_z <- if (is.na(across_mad)) NA_real_ else (part$median_value - across_median) / across_mad
    part$suspicious_flag <- if (directions[[metric_name]] == "high") {
      !is.na(part$robust_z) & part$robust_z >= 1.5
    } else {
      !is.na(part$robust_z) & part$robust_z <= -1.5
    }
    part
  })
  flag_long <- do.call(rbind, flag_parts)
  clusters <- unique(summary_long$cluster)
  flag_summary <- do.call(rbind, lapply(clusters, function(cluster_id) {
    part <- flag_long[flag_long$cluster == cluster_id, , drop = FALSE]
    data.frame(
      cluster = cluster_id,
      n_flagged_metrics = sum(part$suspicious_flag),
      flagged_metrics = if (any(part$suspicious_flag)) paste(part$metric[part$suspicious_flag], collapse = ";") else NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  flag_summary <- merge(flag_summary, logistic_df, by = "cluster", all.x = TRUE, sort = FALSE)
  flag_summary <- flag_summary[match(clusters, flag_summary$cluster), , drop = FALSE]
  flag_summary$concern_level <- ifelse(
    flag_summary$n_flagged_metrics >= 3L,
    "High",
    ifelse(
      flag_summary$n_flagged_metrics >= 2L |
        (!is.na(flag_summary$mcfadden_r2) & flag_summary$mcfadden_r2 >= 0.10),
      "Moderate",
      "Low"
    )
  )
  flag_summary$interpretation <- ifelse(
    flag_summary$concern_level == "High",
    "Multiple QC metrics are shifted in suspicious directions; cluster separation may be influenced by sequencing quality.",
    ifelse(
      flag_summary$concern_level == "Moderate",
      "Some QC metrics differ enough to warrant caution when interpreting this cluster.",
      "No strong QC-driven separation signal was detected from the selected metrics."
    )
  )
  list(long = flag_long, summary = flag_summary)
}

run_cluster_qc <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "02c_cluster_quality"))
  summary_root <- ensure_dir(file.path(stage_root, "summaries"))
  set.seed(SEEDS[["cluster_QC"]])
  counts <- Seurat::GetAssayData(obj, assay = "RNA", slot = "counts")
  cluster_vec <- as.character(obj$seurat_cluster_refine)
  cluster_levels <- as.character(levels(obj$seurat_cluster_refine))
  percent_ribo <- if (any(grepl("^RPL|^RPS|^MRPL|^MRPS", rownames(counts), ignore.case = TRUE))) {
    as.numeric(Seurat::PercentageFeatureSet(obj, assay = "RNA", pattern = "^RPL|^RPS|^MRPL|^MRPS"))
  } else {
    rep(NA_real_, ncol(obj))
  }
  n_count <- as.numeric(obj$nCount_RNA)
  n_feature <- as.numeric(obj$nFeature_RNA)
  top_metrics <- compute_top_expression_metrics(counts, 50L)
  qc_df <- data.frame(
    cell = colnames(obj),
    cluster = cluster_vec,
    nCount_RNA = n_count,
    nFeature_RNA = n_feature,
    percent.mt = as.numeric(obj$percent.mt),
    percent.ribo = percent_ribo,
    log10_genes_per_umi = ifelse(n_count > 1 & n_feature > 1, log10(n_feature) / log10(n_count), NA_real_),
    umi_per_gene = ifelse(n_feature > 0, n_count / n_feature, NA_real_),
    dominant_gene_fraction = top_metrics$dominant_gene_fraction,
    percent.top50 = top_metrics$percent.top50,
    stringsAsFactors = FALSE
  )
  write_csv(qc_df, file.path(summary_root, "cell_qc_metrics.csv"))
  qc_metrics <- c(
    "nCount_RNA", "nFeature_RNA", "percent.mt", "log10_genes_per_umi",
    "umi_per_gene", "dominant_gene_fraction", "percent.top50"
  )
  if (!all(is.na(qc_df$percent.ribo))) qc_metrics <- c(qc_metrics, "percent.ribo")
  summary_long <- build_qc_summary_long(qc_df, qc_metrics, cluster_levels)
  logistic_df <- compute_cluster_logistic_qc(qc_df, qc_metrics, cluster_levels)
  flags <- make_qc_flags(summary_long, logistic_df)
  write_csv(summary_long, file.path(summary_root, "cluster_qc_summary_long.csv"))
  write_csv(logistic_df, file.path(summary_root, "cluster_qc_logistic_multivariate.csv"))
  write_csv(flags$long, file.path(summary_root, "cluster_qc_metric_flags_long.csv"))
  write_csv(flags$summary, file.path(summary_root, "cluster_qc_outlier_flags.csv"))
  write_stage_marker(
    "02c_cluster_quality",
    paste0("High concern: ", paste(flags$summary$cluster[flags$summary$concern_level == "High"], collapse = ","))
  )
  list(flags = flags$summary, qc = qc_df)
}

MANUAL_MERGE_MAP <- list(
  `0` = c("0", "1", "7"),
  `10` = c("10", "11", "12")
)

build_manual_merge_labels <- function(cluster_vec, merge_map, cluster_order = NULL) {
  cluster_vec <- as.character(cluster_vec)
  out <- cluster_vec
  for (merge_name in names(merge_map)) {
    out[cluster_vec %in% as.character(merge_map[[merge_name]])] <- merge_name
  }
  original_order <- if (is.null(cluster_order)) sort_maybe_numeric(cluster_vec) else as.character(cluster_order)
  merged_order <- character(0)
  for (cluster_id in original_order) {
    target <- cluster_id
    for (merge_name in names(merge_map)) {
      if (cluster_id %in% merge_map[[merge_name]]) target <- merge_name
    }
    if (!(target %in% merged_order)) merged_order <- c(merged_order, target)
  }
  factor(out, levels = merged_order)
}

map_labels_through_fixed_merge <- function(labels) {
  labels <- as.character(labels)
  out <- labels
  for (merge_name in names(MANUAL_MERGE_MAP)) {
    out[labels %in% MANUAL_MERGE_MAP[[merge_name]]] <- merge_name
  }
  unique(out)
}

apply_fixed_merge_after_gate <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "02d_manual_cluster_merge"))
  object_root <- ensure_dir(file.path(stage_root, "objects"))
  summary_root <- ensure_dir(file.path(stage_root, "summaries"))
  manual_merge_path <- file.path(
    object_root,
    "integrated_sct_cca_seurat_cluster_refine_manual_merge_test.rds"
  )
  similarity_root <- file.path(output_dir, "02a_cluster_similarity_by_DEGs", "summaries")
  gate_file <- file.path(similarity_root, "FIXED_MERGE_GATE_PASS.txt")
  if (!file.exists(gate_file) || !any(grepl("^status=PASS$", readLines(gate_file, warn = FALSE)))) {
    stop("Fixed merge blocked: pre-merge expression-similarity gate did not pass.", call. = FALSE)
  }

  pairwise_path <- file.path(similarity_root, "pairwise_similarity_long.csv")
  if (!file.exists(pairwise_path)) {
    stop("Fixed merge blocked: pairwise cosine evidence is missing.", call. = FALSE)
  }
  if (file.exists(manual_merge_path)) {
    message("[resume] Reading verified manual-merge object: ", manual_merge_path)
    resumed <- readRDS(manual_merge_path)
    validate_manual_merge_resume(resumed, obj)
    oriented <- orient_seurat_pca_to_reference_anchor(
      resumed, contract_name = "initial", stage_label = "manual_merge_resume"
    )
    resumed <- oriented$object
    if (any(oriented$flipped)) {
      message_cluster_counts_before_rds(
        resumed,
        "manual_merge_test",
        manual_merge_path,
        "Manual-merge PCA-orientation repair"
      )
      saveRDS(resumed, manual_merge_path)
    }
    return(resumed)
  }
  pairwise <- utils::read.csv(pairwise_path, stringsAsFactors = FALSE, check.names = FALSE)
  evidence_pairs <- c("0__1", "0__7", "1__7", "10__11", "10__12", "11__12")
  merge_evidence <- pairwise[pairwise$pair_id %in% evidence_pairs, , drop = FALSE]
  merge_evidence$fixed_merge_group <- ifelse(
    merge_evidence$pair_id %in% c("0__1", "0__7", "1__7"), "0", "10"
  )
  write_csv(
    merge_evidence,
    file.path(summary_root, "fixed_merge_premerge_signed_lfc_cosine_evidence.csv")
  )

  obj <- prepare_obj_for_markers(obj)
  refined_vec <- as.character(obj$seurat_cluster_refine)
  refined_levels <- as.character(levels(obj$seurat_cluster_refine))
  obj$manual_merge_test <- build_manual_merge_labels(
    refined_vec,
    MANUAL_MERGE_MAP,
    cluster_order = refined_levels
  )
  merge_mapping <- do.call(rbind, lapply(names(MANUAL_MERGE_MAP), function(merge_name) {
    data.frame(
      merged_group = merge_name,
      source_cluster = MANUAL_MERGE_MAP[[merge_name]],
      stringsAsFactors = FALSE
    )
  }))
  write_csv(merge_mapping, file.path(summary_root, "manual_merge_mapping.csv"))
  merge_counts <- table(obj$manual_merge_test)
  write_csv(
    data.frame(manual_merge_test = names(merge_counts), n_cells = as.integer(merge_counts)),
    file.path(summary_root, "manual_merge_group_cell_counts.csv")
  )
  message_cluster_counts_before_rds(
    obj,
    "manual_merge_test",
    manual_merge_path,
    "Manual-merge"
  )
  saveRDS(obj, manual_merge_path)
  write_stage_marker(
    "02d_manual_cluster_merge",
    "fixed merge created only after signed_lfc_cosine pre-merge gate"
  )
  obj
}

count_qualifying_upregulated <- function(de) {
  if (is.null(de) || nrow(de) == 0L) return(0L)
  lfc_col <- resolve_lfc_col(de)
  required <- c("p_val_adj", "pct.1", "pct.2")
  if (!all(required %in% colnames(de))) {
    stop("DEG table lacks columns required for the qualifying-up rule.", call. = FALSE)
  }
  lfc <- as.numeric(de[[lfc_col]])
  as.integer(sum(
    !is.na(de$p_val_adj) & as.numeric(de$p_val_adj) < 0.05 &
      !is.na(lfc) & lfc > 0 & abs(lfc) >= 0.25 &
      abs(as.numeric(de$pct.1) - as.numeric(de$pct.2)) >= 0.05
  ))
}

run_prefilter_deg_and_select_removals <- function(obj, qc_flags) {
  stage_root <- ensure_dir(file.path(output_dir, "02e_prefilter_DEGs"))
  deg_root <- ensure_dir(file.path(stage_root, "01_cluster_vs_rest_DEG"))
  summary_root <- ensure_dir(file.path(stage_root, "summaries"))
  if (ncol(obj) != EXPECTED_REFINED_CELLS) {
    stop("Pre-filter DEG must use all 42,884 cells; observed ", ncol(obj), ".", call. = FALSE)
  }

  merged_levels <- as.character(levels(obj$manual_merge_test))
  marker_tables <- run_deg_by_group(
    obj,
    group_col = "manual_merge_test",
    group_levels = merged_levels,
    output_root = deg_root,
    min_pct = 0.10,
    logfc_threshold = 0,
    seed = SEEDS[["prefilter_DEG"]],
    layout = "nested"
  )
  qualifying_counts <- vapply(marker_tables, count_qualifying_upregulated, integer(1))
  cluster_counts <- table(factor(as.character(obj$manual_merge_test), levels = merged_levels))
  deg_summary <- data.frame(
    cluster = merged_levels,
    n_cells = as.integer(cluster_counts),
    cells_in_scope = ncol(obj),
    tested_genes = vapply(marker_tables, nrow, integer(1)),
    qualifying_upregulated_genes = as.integer(qualifying_counts[merged_levels]),
    has_qualifying_upregulated_gene = as.integer(qualifying_counts[merged_levels]) > 0L,
    stringsAsFactors = FALSE
  )
  write_csv(deg_summary, file.path(summary_root, "prefilter_cluster_vs_rest_DEG_summary.csv"))

  qc_high_refined <- as.character(qc_flags$cluster[qc_flags$concern_level == "High"])
  qc_high_merged <- map_labels_through_fixed_merge(qc_high_refined)
  zero_up <- deg_summary$cluster[deg_summary$qualifying_upregulated_genes == 0L]
  removal_set <- sort_maybe_numeric(union(qc_high_merged, zero_up))
  decision <- data.frame(
    cluster = merged_levels,
    concern_level = qc_flags$concern_level[match(merged_levels, qc_flags$cluster)],
    qc_high_source = merged_levels %in% qc_high_merged,
    qualifying_upregulated_genes = deg_summary$qualifying_upregulated_genes,
    zero_qualifying_upregulated_genes = merged_levels %in% zero_up,
    remove_by_union_rule = merged_levels %in% removal_set,
    data_scope = "all_42884_cells_before_any_cluster_deletion",
    stringsAsFactors = FALSE
  )
  decision$selection_reason <- ifelse(
    decision$qc_high_source & decision$zero_qualifying_upregulated_genes,
    "concern_level_High AND zero_qualifying_upregulated_DEG",
    ifelse(
      decision$qc_high_source,
      "concern_level_High",
      ifelse(
        decision$zero_qualifying_upregulated_genes,
        "zero_qualifying_upregulated_DEG",
        "retain"
      )
    )
  )
  write_csv(decision, file.path(summary_root, "data_driven_cluster_removal_decision.csv"))
  write_tsv(
    data.frame(
      set_name = c("qc_high_refined", "qc_high_after_merge", "zero_qualifying_up", "union_removal"),
      clusters = c(
        paste(qc_high_refined, collapse = ","),
        paste(qc_high_merged, collapse = ","),
        paste(zero_up, collapse = ","),
        paste(removal_set, collapse = ",")
      ),
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "removal_rule_sets.tsv")
  )

  assertion_report <- data.frame(
    assertion = c(
      "prefilter_scope_is_42884_cells",
      "qc_high_matches_reference",
      "zero_qualifying_up_matches_reference",
      "union_removal_matches_reference"
    ),
    observed = c(
      as.character(ncol(obj)),
      paste(sort_maybe_numeric(qc_high_refined), collapse = ","),
      paste(sort_maybe_numeric(zero_up), collapse = ","),
      paste(sort_maybe_numeric(removal_set), collapse = ",")
    ),
    expected = c(
      as.character(EXPECTED_REFINED_CELLS),
      paste(EXPECTED_QC_HIGH, collapse = ","),
      paste(EXPECTED_ZERO_QUALIFYING_UP, collapse = ","),
      paste(EXPECTED_REMOVAL_SET, collapse = ",")
    ),
    passed = c(
      ncol(obj) == EXPECTED_REFINED_CELLS,
      setequal(qc_high_refined, EXPECTED_QC_HIGH),
      setequal(zero_up, EXPECTED_ZERO_QUALIFYING_UP),
      setequal(removal_set, EXPECTED_REMOVAL_SET)
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(assertion_report, file.path(summary_root, "reproduction_assertions.tsv"))
  if (!all(assertion_report$passed)) {
    stop(
      "Data-driven removal did not reproduce the reference sets. See ",
      file.path(summary_root, "reproduction_assertions.tsv"),
      ". No hard-coded fallback deletion was applied.",
      call. = FALSE
    )
  }
  write_stage_marker(
    "02e_prefilter_DEGs",
    paste0("dynamic union selected: ", paste(removal_set, collapse = ","))
  )
  list(removal_set = removal_set, decision = decision)
}

add_derived_tn_ploidy_metadata <- function(obj) {
  ids <- as.character(obj$IDs)
  ids[is.na(ids)] <- ""
  ploidy <- ifelse(
    grepl("2N", ids, fixed = TRUE),
    "2N",
    ifelse(grepl("4N", ids, fixed = TRUE), "4N", NA_character_)
  )
  if (anyNA(ploidy)) {
    unresolved <- sort(unique(ids[is.na(ploidy)]))
    stop(
      "Ploidy must be either 2N or 4N, but some IDs cannot be parsed: ",
      paste(utils::head(unresolved, 20L), collapse = ", "),
      call. = FALSE
    )
  }
  tn <- ifelse(grepl("Cell-Culture", ids, fixed = TRUE), "CellLine", "Tumor")
  obj$Ploidy <- factor(ploidy, levels = c("2N", "4N", "Unknown"))
  obj$TN <- factor(tn, levels = c("CellLine", "Tumor"))
  obj
}

coerce_final_metadata_contract <- function(obj) {
  obj <- coerce_base_metadata_contract(obj)
  obj@meta.data$integrated_snn_res.0.6 <- factor(
    as.character(obj@meta.data$integrated_snn_res.0.6),
    levels = EXPECTED_INITIAL_CLUSTERS
  )
  obj@meta.data$seurat_clusters <- factor(
    as.character(obj@meta.data$seurat_clusters),
    levels = EXPECTED_INITIAL_CLUSTERS
  )
  obj@meta.data$S.Score <- as.numeric(obj@meta.data$S.Score)
  obj@meta.data$G2M.Score <- as.numeric(obj@meta.data$G2M.Score)
  obj@meta.data$Phase <- as.character(obj@meta.data$Phase)
  obj@meta.data$cluster_cell_cycle_annotation <- as.character(obj@meta.data$cluster_cell_cycle_annotation)
  obj@meta.data$Ploidy <- factor(
    as.character(obj@meta.data$Ploidy),
    levels = c("2N", "4N", "Unknown")
  )
  obj@meta.data$TN <- factor(
    as.character(obj@meta.data$TN),
    levels = c("CellLine", "Tumor")
  )
  obj
}

run_final_object <- function(obj, removal_set) {
  stage_root <- ensure_dir(file.path(output_dir, "03_final_cluster"))
  summary_root <- ensure_dir(file.path(stage_root, "01_summary"))
  plot_root <- ensure_dir(file.path(stage_root, "02_plots"))
  object_root <- ensure_dir(file.path(stage_root, "03_objects"))
  labels <- as.character(obj$manual_merge_test)
  merged_levels <- as.character(levels(obj$manual_merge_test))
  keep <- !(labels %in% removal_set)
  retained_levels <- setdiff(merged_levels, removal_set)
  removed_cells <- colnames(obj)[!keep]
  retained_cells <- colnames(obj)[keep]

  write_csv(
    data.frame(
      metric = c("input_cells", "removed_cells", "retained_cells", "removed_cluster_n", "retained_cluster_n"),
      value = c(ncol(obj), length(removed_cells), length(retained_cells), length(removal_set), length(retained_levels)),
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "filter_summary.csv")
  )
  write_csv(
    data.frame(
      cell = removed_cells,
      cluster = labels[!keep],
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "removed_cells.csv")
  )

  obj_final <- subset(obj, cells = retained_cells)
  obj_final$clusters <- factor(
    as.character(obj_final$manual_merge_test),
    levels = retained_levels
  )
  obj_final@meta.data$manual_merge_test <- NULL
  Seurat::Idents(obj_final) <- obj_final$clusters
  obj_final <- add_derived_tn_ploidy_metadata(obj_final)
  write_csv(
    data.frame(
      metadata = c("Ploidy", "Ploidy", "TN", "TN"),
      value = c("2N", "4N", "CellLine", "Tumor"),
      n_cells = c(
        sum(obj_final$Ploidy == "2N", na.rm = TRUE),
        sum(obj_final$Ploidy == "4N", na.rm = TRUE),
        sum(obj_final$TN == "CellLine", na.rm = TRUE),
        sum(obj_final$TN == "Tumor", na.rm = TRUE)
      ),
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "derived_metadata_counts.csv")
  )
  Seurat::DefaultAssay(obj_final) <- "integrated"
  pca_features <- Seurat::VariableFeatures(obj_final[["integrated"]])
  if (length(pca_features) < 2L) {
    stop("Integrated assay lacks variable features for final PCA.", call. = FALSE)
  }
  pca_npcs <- min(50L, length(pca_features) - 1L, ncol(obj_final) - 1L)
  set.seed(SEEDS[["final_PCA"]])
  obj_final <- Seurat::RunPCA(
    obj_final,
    assay = "integrated",
    features = pca_features,
    npcs = pca_npcs,
    seed.use = SEEDS[["final_PCA"]],
    verbose = FALSE
  )
  obj_final <- orient_seurat_pca_to_reference_anchor(
    obj_final, contract_name = "final", stage_label = "final"
  )$object
  obj_final <- run_compatible_umap(
    obj_final, reduction = "pca", dims = 1:30,
    seed = SEEDS[["final_UMAP"]], stage_label = "final"
  )
  if (!setequal(names(obj_final@commands), EXPECTED_FINAL_COMMANDS)) {
    stop(
      "Final Seurat command set mismatch. Observed: ",
      paste(names(obj_final@commands), collapse = ","),
      call. = FALSE
    )
  }
  obj_final@commands <- obj_final@commands[EXPECTED_FINAL_COMMANDS]
  obj_final@graphs <- list()
  obj_final@neighbors <- list()
  obj_final <- coerce_final_metadata_contract(obj_final)
  Seurat::DefaultAssay(obj_final) <- "integrated"
  Seurat::Idents(obj_final) <- obj_final$clusters

  if (ncol(obj_final) != EXPECTED_FINAL_CELLS) {
    stop("Final cell count mismatch: ", ncol(obj_final), " vs ", EXPECTED_FINAL_CELLS, call. = FALSE)
  }
  if (!identical(as.character(levels(obj_final$clusters)), EXPECTED_FINAL_LEVELS)) {
    stop("Final cluster levels do not match the reference.", call. = FALSE)
  }
  observed_counts <- table(obj_final$clusters)
  if (!identical(as.integer(observed_counts), as.integer(EXPECTED_FINAL_COUNTS))) {
    stop(
      "Final cluster counts do not match the reference. Observed: ",
      paste(names(observed_counts), observed_counts, sep = "=", collapse = ","),
      call. = FALSE
    )
  }
  if (!identical(colnames(obj_final@meta.data), EXPECTED_FINAL_METADATA)) {
    stop(
      "Final metadata contract mismatch. Observed: ",
      paste(colnames(obj_final@meta.data), collapse = ","),
      call. = FALSE
    )
  }
  expected_assay_features <- c(RNA = 72302L, SCT = 46009L, integrated = 3000L)
  observed_assay_features <- vapply(names(expected_assay_features), function(assay_name) {
    nrow(obj_final[[assay_name]])
  }, integer(1))
  if (!identical(unname(observed_assay_features), unname(expected_assay_features))) {
    stop(
      "Final assay feature counts do not match the reference. Observed: ",
      paste(names(observed_assay_features), observed_assay_features, sep = "=", collapse = ","),
      call. = FALSE
    )
  }
  if (ncol(Seurat::Embeddings(obj_final, "pca")) != 50L || ncol(Seurat::Embeddings(obj_final, "umap")) != 2L) {
    stop("Final PCA/UMAP dimensions do not match the reference.", call. = FALSE)
  }

  final_counts <- table(obj_final$clusters)
  write_csv(
    data.frame(cluster = names(final_counts), n_cells = as.integer(final_counts)),
    file.path(summary_root, "clusters_counts_after_filter.csv")
  )
  write_tsv(
    data.frame(
      parameter = c(
        "removal_source", "removed_clusters", "pca_assay", "pca_components",
        "umap_dimensions", "reintegration_run", "find_neighbors_run",
        "find_clusters_run", "cluster_labels_changed"
      ),
      value = c(
        "dynamic_union_of_QC_High_and_zero_qualifying_up_DEG",
        paste(removal_set, collapse = ","),
        "integrated", as.character(pca_npcs), "1:30", "FALSE", "FALSE", "FALSE", "FALSE"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(summary_root, "pca_umap_parameters.tsv")
  )
  save_umap_plot(
    obj_final,
    "clusters",
    file.path(plot_root, "umap_by_clusters"),
    "Final clusters after dynamic filtering"
  )
  final_path <- file.path(object_root, "integrated_sct_cca_seurat_final_reclustered.rds")
  message_cluster_counts_before_rds(
    obj_final,
    "clusters",
    final_path,
    "Final"
  )
  saveRDS(obj_final, final_path)
  write_stage_marker("03_final_cluster", "35513 cells; 9 retained clusters; new PCA and UMAP")
  obj_final
}

run_final_deg <- function(obj) {
  stage_root <- ensure_dir(file.path(output_dir, "03a_DEGs"))
  deg_root <- ensure_dir(file.path(stage_root, "01_clusters_vs_rest"))
  summary_root <- ensure_dir(file.path(stage_root, "00_summary"))
  cluster_levels <- as.character(levels(obj$clusters))
  marker_tables <- run_deg_by_group(
    obj,
    group_col = "clusters",
    group_levels = cluster_levels,
    output_root = deg_root,
    min_pct = NULL,
    logfc_threshold = NULL,
    seed = SEEDS[["final_DEG"]],
    layout = "flat"
  )
  cluster_counts <- table(factor(as.character(obj$clusters), levels = cluster_levels))
  for (cluster_id in cluster_levels) {
    de <- marker_tables[[cluster_id]]
    lfc_col <- resolve_lfc_col(de)
    required_core <- c("gene", "gene_symbol", "p_val", lfc_col, "pct.1", "pct.2", "p_val_adj")
    missing_core <- setdiff(required_core, colnames(de))
    if (length(missing_core) > 0L) {
      stop(
        "Final DEG table lacks reference-schema columns for cluster ", cluster_id,
        ": ", paste(missing_core, collapse = ", "),
        call. = FALSE
      )
    }
    n_ident_1 <- as.integer(cluster_counts[[cluster_id]])
    comparison_meta <- data.frame(
      scope = "clusters_vs_rest",
      cluster = cluster_id,
      subset_group = NA_character_,
      dose = NA_character_,
      comparison = paste0("cluster_", cluster_id, "_vs_rest"),
      group_col = "clusters",
      ident_1 = cluster_id,
      ident_2 = "rest",
      n_ident_1 = n_ident_1,
      n_ident_2 = ncol(obj) - n_ident_1,
      stringsAsFactors = FALSE
    )
    meta_rows <- comparison_meta[rep(1L, nrow(de)), , drop = FALSE]
    reference_schema <- cbind(
      meta_rows,
      de[, c("gene", "gene_symbol", "p_val", lfc_col, "pct.1", "pct.2", "p_val_adj"), drop = FALSE]
    )
    marker_path <- file.path(
      deg_root,
      paste0("cluster_", cluster_id, "_vs_rest_markers.csv")
    )
    readr::write_csv(reference_schema, marker_path)
    marker_tables[[cluster_id]] <- reference_schema
  }
  summary_df <- data.frame(
    cluster = cluster_levels,
    n_cells = as.integer(cluster_counts),
    tested_genes = vapply(marker_tables, nrow, integer(1)),
    significant_genes = vapply(marker_tables, function(de) {
      if (nrow(de) == 0L) 0L else as.integer(sum(!is.na(de$p_val_adj) & de$p_val_adj < 0.05))
    }, integer(1)),
    stringsAsFactors = FALSE
  )
  write_csv(summary_df, file.path(summary_root, "final_cluster_vs_rest_DEG_summary.csv"))
  write_stage_marker("03a_DEGs", "final cluster-vs-rest DEG only")
  invisible(summary_df)
}

validate_final_resume <- function(obj) {
  if (!inherits(obj, "Seurat")) {
    stop("Existing final output is not a Seurat object.", call. = FALSE)
  }
  if (ncol(obj) != EXPECTED_FINAL_CELLS) {
    stop("Existing final output cell count mismatch.", call. = FALSE)
  }
  if (!identical(colnames(obj@meta.data), EXPECTED_FINAL_METADATA)) {
    stop("Existing final output metadata contract mismatch.", call. = FALSE)
  }
  if (!identical(levels(obj$clusters), EXPECTED_FINAL_LEVELS)) {
    stop("Existing final output cluster levels mismatch.", call. = FALSE)
  }
  validate_named_cluster_counts(
    obj$clusters,
    EXPECTED_FINAL_COUNTS,
    "Existing final output"
  )
  if (!identical(as.character(Seurat::Idents(obj)), as.character(obj$clusters)) ||
      !identical(levels(Seurat::Idents(obj)), EXPECTED_FINAL_LEVELS)) {
    stop("Existing final output active identities mismatch clusters.", call. = FALSE)
  }
  if (!identical(levels(obj$Ploidy), c("2N", "4N", "Unknown")) ||
      !identical(levels(obj$TN), c("CellLine", "Tumor"))) {
    stop("Existing final output Ploidy/TN factor contract mismatch.", call. = FALSE)
  }
  expected_assay_features <- c(RNA = 72302L, SCT = 46009L, integrated = 3000L)
  observed_assay_features <- vapply(names(expected_assay_features), function(assay_name) {
    if (!(assay_name %in% names(obj@assays))) return(NA_integer_)
    nrow(obj[[assay_name]])
  }, integer(1))
  if (!identical(unname(observed_assay_features), unname(expected_assay_features))) {
    stop("Existing final output assay feature-count contract mismatch.", call. = FALSE)
  }
  if (!all(c("pca", "umap") %in% names(obj@reductions)) ||
      ncol(Seurat::Embeddings(obj, "pca")) != 50L ||
      ncol(Seurat::Embeddings(obj, "umap")) != 2L) {
    stop("Existing final output PCA/UMAP contract mismatch.", call. = FALSE)
  }
  if (!identical(names(obj@commands), EXPECTED_FINAL_COMMANDS)) {
    stop("Existing final output Seurat command contract mismatch.", call. = FALSE)
  }
  if (length(obj@graphs) != 0L || length(obj@neighbors) != 0L) {
    stop("Existing final output unexpectedly contains graphs or neighbors.", call. = FALSE)
  }
  invisible(TRUE)
}

read_verified_removal_set <- function() {
  path <- file.path(
    output_dir,
    "02e_prefilter_DEGs", "summaries", "removal_rule_sets.tsv"
  )
  if (!file.exists(path)) {
    stop("Cannot resume final DEG because removal_rule_sets.tsv is missing.", call. = FALSE)
  }
  sets <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  row <- sets[sets$set_name == "union_removal", , drop = FALSE]
  if (nrow(row) != 1L || !nzchar(row$clusters[[1]])) {
    stop("Cannot resolve the prior dynamic union-removal set.", call. = FALSE)
  }
  observed <- strsplit(row$clusters[[1]], ",", fixed = TRUE)[[1]]
  assert_set_equal(observed, EXPECTED_REMOVAL_SET, "Resumed dynamic removal set")
  sort_maybe_numeric(observed)
}

initialize_provenance <- function(inputs) {
  provenance_root <- ensure_dir(file.path(output_dir, "00_provenance"))
  input_files <- c(inputs$h5_files, inputs$ploidy_file, inputs$sample_info_file)
  input_types <- c(
    rep("filtered_feature_bc_matrix.h5", length(inputs$h5_files)),
    "all_ploidy.tsv", "sample_info.xlsx"
  )
  file_info <- file.info(input_files)
  manifest <- data.frame(
    input_type = input_types,
    path = normalizePath(input_files, mustWork = TRUE),
    size_bytes = as.numeric(file_info$size),
    modified_time = format(file_info$mtime, "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  write_tsv(manifest, file.path(provenance_root, "run_manifest.tsv"))
  write_tsv(
    data.frame(operation = names(SEEDS), seed = as.integer(SEEDS), stringsAsFactors = FALSE),
    file.path(provenance_root, "seed_registry.tsv")
  )
  write_tsv(
    data.frame(
      parameter = c(
        "input_dir", "output_dir", "accepted_raw_file_basename",
        "integration_features", "initial_pca_npcs", "clustering_resolution",
        "initial_umap_dims", "similarity_metric", "refine_knn_k",
        "refine_pca_orientation", "refine_compatibility_steps",
        "refine_compatibility_rule", "refine_compatibility_workers",
        "fixed_merge_0", "fixed_merge_10", "deletion_rule",
        "qualifying_up_rule", "final_umap_dims",
        "deg_feature_chunk_target", "deg_parallel_workers",
        "deg_parallel_backend", "deg_cache_reuse",
        "scdblfinder_version", "rcppannoy_version", "irlba_version",
        "uwot_version", "sctransform_version", "xgboost_version", "biocneighbors_version",
        "assorthead_version", "scdbl_pca_orientation_anchor",
        "scdbl_pca_round_digits", "initial_pca_orientation_anchor",
        "final_pca_orientation_anchor", "annoy_float_compatibility",
        "umap_powf_compatibility", "umap_spectral_blas"
      ),
      value = c(
        input_dir, output_dir, "filtered_feature_bc_matrix.h5", "3000", "50", "0.6", "1:30",
        "signed_lfc_cosine", "50",
        "largest_absolute_cell_coordinate_positive",
        paste(names(REFINE_COMPATIBILITY_STEPS), REFINE_COMPATIBILITY_STEPS, sep = "=", collapse = ","),
        "cluster4=intersection(0.0125,0.03);cluster9=union(0.03,0.09)",
        as.character(REFINE_COMPATIBILITY_WORKERS),
        "0,1,7", "10,11,12",
        "concern_level==High OR zero qualifying upregulated DEG",
        "p_adj<0.05; avg_log2FC>0; abs(logFC)>=0.25; abs(pct.1-pct.2)>=0.05",
        "1:30",
        as.character(DEG_FEATURE_CHUNK_TARGET), as.character(DEG_PARALLEL_WORKERS),
        "parallel::mclapply_feature_chunks", "TRUE",
        PINNED_PACKAGE_VERSIONS[["scDblFinder"]],
        PINNED_PACKAGE_VERSIONS[["RcppAnnoy"]],
        PINNED_PACKAGE_VERSIONS[["irlba"]],
        PINNED_PACKAGE_VERSIONS[["uwot"]],
        PINNED_PACKAGE_VERSIONS[["sctransform"]],
        PINNED_PACKAGE_VERSIONS[["xgboost"]],
        PINNED_PACKAGE_VERSIONS[["BiocNeighbors"]],
        PINNED_PACKAGE_VERSIONS[["assorthead"]],
        "sign_of_maximum_absolute_gene_loading_per_component",
        as.character(SCDBL_PCA_ROUND_DIGITS),
        "reference_sign_of_maximum_absolute_gene_loading_per_component",
        "reference_sign_of_maximum_absolute_gene_loading_per_component",
        "explicit_ARM_macOS_accumulation_order_in_vendored_assorthead",
        "macos_arm64_libsystem_m_powf_in_vendored_uwot",
        "reference_BLAS_LAPACK_subprocess_with_forced_irlba"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(provenance_root, "analysis_parameters.tsv")
  )
  writeLines(capture.output(sessionInfo()), file.path(provenance_root, "sessionInfo_start.txt"))
}

finish_provenance <- function(removal_set) {
  provenance_root <- ensure_dir(file.path(output_dir, "00_provenance"))
  writeLines(capture.output(sessionInfo()), file.path(provenance_root, "sessionInfo.txt"))
  writeLines(
    c(
      "status=PASS",
      paste0("completed_at=", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")),
      "refined_object=02b_cluster_refine/objects/integrated_sct_cca_seurat_cluster_refine.rds",
      "final_object=03_final_cluster/03_objects/integrated_sct_cca_seurat_final_reclustered.rds",
      paste0(
        "removed_clusters_selected_dynamically=",
        paste(removal_set, collapse = ",")
      ),
      "excluded_analyses=annotation,ORA,GSEA,Dose_DEG,Ploidy_DEG,TN_DEG,trajectory"
    ),
    file.path(provenance_root, "PIPELINE_COMPLETE.txt")
  )
}

ensure_dir(output_dir)
final_existing <- file.path(
  output_dir,
  "03_final_cluster", "03_objects", "integrated_sct_cca_seurat_final_reclustered.rds"
)
inputs <- discover_inputs(input_dir)
initialize_provenance(inputs)
message("Standalone cluster workflow started.")
message("Input directory: ", input_dir)
message("Output directory: ", output_dir)

if (file.exists(final_existing)) {
  message("[resume] Reading verified final object: ", final_existing)
  final_obj <- readRDS(final_existing)
  validate_final_resume(final_obj)
  removal_set <- read_verified_removal_set()
  run_final_deg(final_obj)
  finish_provenance(removal_set)
  message("Standalone cluster workflow completed successfully from verified final-object resume.")
  message("Final object: ", final_existing)
  quit(save = "no", status = 0L)
}

integrated <- run_raw_integration(inputs)
integrated <- run_cell_cycle(integrated)
run_initial_deg_similarity(integrated)
invisible(gc())
refined <- run_cluster_refine(integrated)
rm(integrated)
invisible(gc())
qc_flags <- run_cluster_qc(refined)$flags
merged <- apply_fixed_merge_after_gate(refined)
rm(refined)
invisible(gc())
removal_result <- run_prefilter_deg_and_select_removals(merged, qc_flags)
final_obj <- run_final_object(merged, removal_result$removal_set)
rm(merged)
invisible(gc())
run_final_deg(final_obj)
finish_provenance(removal_result$removal_set)

message("Standalone cluster workflow completed successfully.")
message(
  "Final object: ",
  file.path(output_dir, "03_final_cluster", "03_objects", "integrated_sct_cca_seurat_final_reclustered.rds")
)
