#!/usr/bin/env Rscript

script_path <- NULL
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
file_match <- grep(file_arg, cmd_args, value = TRUE)
if (length(file_match) > 0) {
  script_path <- normalizePath(sub(file_arg, "", file_match[1]), mustWork = FALSE)
}
if (is.null(script_path) || !nzchar(script_path)) {
  frame_files <- vapply(
    sys.frames(),
    function(x) {
      if (!is.null(x$ofile)) normalizePath(x$ofile, mustWork = FALSE) else NA_character_
    },
    character(1)
  )
  frame_files <- frame_files[!is.na(frame_files)]
  if (length(frame_files) > 0) {
    script_path <- frame_files[length(frame_files)]
  }
}
script_dir <- if (!is.null(script_path) && nzchar(script_path)) dirname(script_path) else getwd()
source(file.path(script_dir, "Utils.R"))
config <- load_in_vivo_config(file.path(script_dir, "in_vivo_config.yaml"))
results_root <- get_results_root(config)

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(1234)

# -----------------------------
# User-configurable parameters
# -----------------------------
input_root <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/A02_cellRanger"
ploidy_file <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/all_ploidy.tsv"
sample_info_file <- "/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/data/SUM-159/sample_info.xlsx"
output_dir <- file.path(results_root, "01_data")

special_samples <- c("2N-Cell-Culture", "4N-Cell-Culture")

# QC for special samples only
qc_min_features <- 200
qc_min_counts <- 500
qc_max_percent_mt <- 20
qc_max_features_quantile <- 0.99
qc_max_counts_quantile <- 0.99

# Integration settings
nfeatures_integration <- 3000
npcs <- 30
cluster_resolution <- 0.6
future_globals_maxsize_gb <- 60

assert_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required. Please install it first.")
  }
}

for (pkg in c("Seurat", "scDblFinder", "SingleCellExperiment")) {
  assert_pkg(pkg)
}

read_sample_info <- function(xlsx_path) {
  if (requireNamespace("readxl", quietly = TRUE)) {
    df <- readxl::read_excel(xlsx_path)
    return(as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Please install 'readxl' (recommended) or 'xml2' to read sample_info.xlsx")
  }

  tmp_dir <- tempfile("sample_info_xlsx_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  utils::unzip(
    zipfile = xlsx_path,
    files = c("xl/sharedStrings.xml", "xl/worksheets/sheet1.xml"),
    exdir = tmp_dir
  )

  shared_file <- file.path(tmp_dir, "xl/sharedStrings.xml")
  sheet_file <- file.path(tmp_dir, "xl/worksheets/sheet1.xml")
  if (!file.exists(sheet_file)) {
    stop("Cannot find sheet1.xml in sample_info.xlsx: ", xlsx_path)
  }

  ns <- c(x = "http://schemas.openxmlformats.org/spreadsheetml/2006/main")
  shared_strings <- character(0)
  if (file.exists(shared_file)) {
    shared_doc <- xml2::read_xml(shared_file)
    si_nodes <- xml2::xml_find_all(shared_doc, ".//x:si", ns = ns)
    shared_strings <- vapply(
      si_nodes,
      function(node) {
        txt <- xml2::xml_find_all(node, ".//x:t", ns = ns)
        paste(xml2::xml_text(txt), collapse = "")
      },
      FUN.VALUE = character(1)
    )
  }

  excel_col_to_index <- function(cell_ref) {
    col <- gsub("[0-9]+", "", cell_ref)
    chars <- utf8ToInt(col) - utf8ToInt("A") + 1L
    sum(chars * (26L ^ rev(seq_along(chars) - 1L)))
  }

  read_cell <- function(cell_node) {
    cell_type <- xml2::xml_attr(cell_node, "t")
    v_node <- xml2::xml_find_first(cell_node, "x:v", ns = ns)
    if (inherits(v_node, "xml_missing")) {
      return(NA_character_)
    }

    value <- xml2::xml_text(v_node)
    if (!is.na(cell_type) && cell_type == "s") {
      idx <- as.integer(value) + 1L
      if (is.na(idx) || idx < 1L || idx > length(shared_strings)) {
        return(NA_character_)
      }
      return(shared_strings[idx])
    }
    value
  }

  sheet_doc <- xml2::read_xml(sheet_file)
  row_nodes <- xml2::xml_find_all(sheet_doc, ".//x:sheetData/x:row", ns = ns)
  if (length(row_nodes) < 2L) {
    stop("sample_info.xlsx has no data rows: ", xlsx_path)
  }

  parsed_rows <- lapply(row_nodes, function(row_node) {
    cell_nodes <- xml2::xml_find_all(row_node, "x:c", ns = ns)
    if (length(cell_nodes) == 0L) {
      return(character(0))
    }
    refs <- xml2::xml_attr(cell_nodes, "r")
    idx <- vapply(refs, excel_col_to_index, FUN.VALUE = integer(1))
    vals <- vapply(cell_nodes, read_cell, FUN.VALUE = character(1))
    out <- rep(NA_character_, max(idx))
    out[idx] <- vals
    out
  })

  max_cols <- max(vapply(parsed_rows, length, FUN.VALUE = integer(1)))
  mat <- do.call(
    rbind,
    lapply(parsed_rows, function(x) {
      length(x) <- max_cols
      x
    })
  )

  header <- mat[1, ]
  keep <- !is.na(header) & header != ""
  data_mat <- mat[-1, keep, drop = FALSE]
  header <- header[keep]

  out <- as.data.frame(data_mat, stringsAsFactors = FALSE, check.names = FALSE)
  names(out) <- header
  out
}

resolve_col <- function(df, candidates, fallback_index = NULL) {
  nm <- names(df)
  for (cand in candidates) {
    idx <- which(tolower(nm) == tolower(cand))
    if (length(idx) == 1L) {
      return(nm[idx])
    }
  }
  if (!is.null(fallback_index) && fallback_index <= length(nm)) {
    return(nm[fallback_index])
  }
  stop("Cannot resolve column: ", paste(candidates, collapse = ", "))
}

read_counts <- function(h5_file) {
  counts <- Seurat::Read10X_h5(h5_file)
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }
  counts
}

qc_filter_special <- function(obj, sample_id) {
  feature_upper <- as.numeric(stats::quantile(obj$nFeature_RNA, probs = qc_max_features_quantile, na.rm = TRUE))
  count_upper <- as.numeric(stats::quantile(obj$nCount_RNA, probs = qc_max_counts_quantile, na.rm = TRUE))

  keep <- obj$nFeature_RNA >= qc_min_features &
    obj$nCount_RNA >= qc_min_counts &
    obj$nFeature_RNA <= feature_upper &
    obj$nCount_RNA <= count_upper &
    obj$percent.mt <= qc_max_percent_mt

  message(sprintf("[%s] QC keep: %d / %d", sample_id, sum(keep), length(keep)))
  subset(obj, cells = colnames(obj)[keep])
}

remove_doublets <- function(obj, sample_id) {
  if (ncol(obj) == 0L) {
    return(obj)
  }

  if (ncol(obj) < 100L) {
    warning(sprintf("[%s] <100 cells after QC, skip scDblFinder and keep all cells.", sample_id))
    obj$scDblFinder.class <- "singlet_low_cell_skip"
    obj$scDblFinder.score <- NA_real_
    return(obj)
  }

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = Seurat::GetAssayData(obj, assay = "RNA", slot = "counts"))
  )
  set.seed(1234)
  sce <- scDblFinder::scDblFinder(sce, verbose = FALSE)
  dbl_meta <- as.data.frame(SingleCellExperiment::colData(sce))

  obj$scDblFinder.class <- as.character(dbl_meta$scDblFinder.class)
  obj$scDblFinder.score <- as.numeric(dbl_meta$scDblFinder.score)

  singlets <- rownames(dbl_meta)[dbl_meta$scDblFinder.class == "singlet"]
  message(sprintf("[%s] Singlet keep: %d / %d", sample_id, length(singlets), ncol(obj)))
  subset(obj, cells = singlets)
}

make_special_metadata <- function(sample_id, cols, seq_col, id_col) {
  out <- as.list(stats::setNames(rep(NA_character_, length(cols)), cols))
  out[[seq_col]] <- sample_id
  out[[id_col]] <- sample_id
  out
}

get_metadata_row <- function(sample_id, sample_info, id_col, seq_col, special_samples) {
  if (sample_id %in% special_samples) {
    return(make_special_metadata(sample_id, names(sample_info), seq_col, id_col))
  }

  idx <- which(as.character(sample_info[[id_col]]) == sample_id)
  if (length(idx) == 0L) {
    stop("Sample not found in sample_info IDs: ", sample_id)
  }
  if (length(idx) > 1L) {
    warning("Multiple sample_info rows for sample ", sample_id, ". Use first row.")
  }

  as.list(sample_info[idx[1], , drop = FALSE])
}

add_sample_metadata <- function(obj, metadata_row) {
  for (col_name in names(metadata_row)) {
    obj[[col_name]] <- metadata_row[[col_name]]
  }
  obj
}

if (!dir.exists(input_root)) {
  stop("Input folder does not exist: ", input_root)
}
if (!file.exists(ploidy_file)) {
  stop("Cannot find all_ploidy.tsv: ", ploidy_file)
}
if (!file.exists(sample_info_file)) {
  stop("Cannot find sample_info.xlsx: ", sample_info_file)
}
.ensure_dir(output_dir)

sample_info <- read_sample_info(sample_info_file)
if (ncol(sample_info) < 3L) {
  stop("sample_info.xlsx must contain at least 3 columns.")
}

harvest_col <- resolve_col(sample_info, c("harvest"), fallback_index = 1L)
seq_col <- resolve_col(sample_info, c("Sequencing IDs", "Sequencing ID", "SequencingIDs"), fallback_index = 2L)
id_col <- resolve_col(sample_info, c("IDs", "ID"), fallback_index = 3L)

sample_info[[harvest_col]] <- as.character(sample_info[[harvest_col]])
sample_info[[id_col]] <- as.character(sample_info[[id_col]])

ploidy <- utils::read.delim(
  ploidy_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!all(c("file", "cell_id") %in% names(ploidy))) {
  stop("all_ploidy.tsv must contain columns: file, cell_id")
}

ploidy$harvest <- sub("\\.sps\\.cbs$", "", ploidy$file)
barcode_by_harvest <- split(ploidy$cell_id, ploidy$harvest)
barcode_by_harvest <- lapply(barcode_by_harvest, unique)

id_to_harvest <- stats::setNames(sample_info[[harvest_col]], sample_info[[id_col]])
barcode_by_id <- lapply(id_to_harvest, function(harvest) {
  if (is.na(harvest) || !harvest %in% names(barcode_by_harvest)) {
    return(character(0))
  }
  barcode_by_harvest[[harvest]]
})

sample_dirs <- sort(list.dirs(input_root, recursive = FALSE, full.names = TRUE))
sample_dirs <- sample_dirs[grepl("-Count-HM$", basename(sample_dirs))]
if (length(sample_dirs) == 0L) {
  stop("No sample directories ending with '-Count-HM' in: ", input_root)
}

seurat_list <- list()
sample_summary <- list()

for (sample_dir in sample_dirs) {
  sample_folder <- basename(sample_dir)
  sample_id <- sub("-Count-HM$", "", sample_folder)
  h5_file <- file.path(sample_dir, "outs", "filtered_feature_bc_matrix.h5")

  if (!file.exists(h5_file)) {
    warning("Skip missing h5: ", sample_folder)
    next
  }

  message("\n========== ", sample_id, " ==========")
  counts <- read_counts(h5_file)
  obj <- Seurat::CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 0,
    min.features = 0
  )
  obj$sample <- sample_id
  obj$sample_folder <- sample_folder
  obj[["percent.mt"]] <- Seurat::PercentageFeatureSet(obj, pattern = "^MT-")

  n_raw <- ncol(obj)
  n_after_qc <- NA_integer_
  n_after_filter <- NA_integer_
  filter_method <- NA_character_

  if (sample_id %in% special_samples) {
    filter_method <- "QC + scDblFinder"
    obj <- qc_filter_special(obj, sample_id)
    n_after_qc <- ncol(obj)
    obj <- remove_doublets(obj, sample_id)
    n_after_filter <- ncol(obj)
  } else {
    filter_method <- "all_ploidy keep barcodes"
    keep_barcodes <- barcode_by_id[[sample_id]]
    if (is.null(keep_barcodes) || length(keep_barcodes) == 0L) {
      warning("No all_ploidy barcode mapping for sample: ", sample_id, ". Skip sample.")
      next
    }
    keep_cells <- intersect(colnames(obj), keep_barcodes)
    message(sprintf("[%s] all_ploidy keep: %d / %d", sample_id, length(keep_cells), ncol(obj)))
    obj <- subset(obj, cells = keep_cells)
    obj$scDblFinder.class <- "singlet_by_ploidy"
    obj$scDblFinder.score <- NA_real_
    n_after_qc <- ncol(obj)
    n_after_filter <- ncol(obj)
  }

  if (ncol(obj) == 0L) {
    warning("No cells remain after filtering: ", sample_id)
    next
  }

  meta_row <- get_metadata_row(sample_id, sample_info, id_col, seq_col, special_samples)
  obj <- add_sample_metadata(obj, meta_row)

  obj$barcode_raw <- colnames(obj)
  obj <- Seurat::RenameCells(obj, new.names = paste0(sample_id, "_", colnames(obj)))

  seurat_list[[sample_id]] <- obj
  sample_summary[[sample_id]] <- data.frame(
    sample = sample_id,
    sample_folder = sample_folder,
    filter_method = filter_method,
    raw_cells = n_raw,
    after_qc_cells = n_after_qc,
    after_filter_cells = n_after_filter,
    stringsAsFactors = FALSE
  )
}

if (length(seurat_list) < 2L) {
  stop("At least 2 non-empty samples are required for SCT-CCA integration.")
}

summary_df <- do.call(rbind, sample_summary)
utils::write.csv(summary_df, file.path(output_dir, "sample_filter_summary.csv"), row.names = FALSE)

message("\nRunning SCTransform ...")
seurat_list <- lapply(
  seurat_list,
  function(x) {
    Seurat::SCTransform(
      object = x,
      assay = "RNA",
      new.assay.name = "SCT",
      vars.to.regress = "percent.mt",
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
  }
)

message("Selecting features and preparing integration ...")
integration_features <- Seurat::SelectIntegrationFeatures(
  object.list = seurat_list,
  nfeatures = nfeatures_integration
)
configure_future_for_seurat(max_size_gb = future_globals_maxsize_gb)
invisible(gc())
seurat_list <- Seurat::PrepSCTIntegration(
  object.list = seurat_list,
  anchor.features = integration_features,
  verbose = FALSE
)

message("Finding SCT-CCA anchors and integrating ...")
anchors <- Seurat::FindIntegrationAnchors(
  object.list = seurat_list,
  normalization.method = "SCT",
  anchor.features = integration_features,
  reduction = "cca",
  verbose = FALSE
)
integrated <- Seurat::IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT",
  verbose = FALSE
)

Seurat::DefaultAssay(integrated) <- "integrated"

message("Running PCA / neighbors / clusters / UMAP ...")
integrated <- Seurat::RunPCA(integrated, npcs = max(50, npcs), verbose = FALSE)
dims_use <- 1:min(npcs, ncol(Seurat::Embeddings(integrated, reduction = "pca")))
integrated <- Seurat::FindNeighbors(integrated, dims = dims_use, verbose = FALSE)
integrated <- Seurat::FindClusters(integrated, resolution = cluster_resolution, verbose = FALSE)
integrated <- Seurat::RunUMAP(integrated, dims = dims_use, verbose = FALSE)

saveRDS(seurat_list, file.path(output_dir, "seurat_list_postfilter.rds"))
saveRDS(integrated, file.path(output_dir, "integrated_sct_cca_seurat.rds"))
writeLines(capture.output(sessionInfo()), con = file.path(output_dir, "sessionInfo_01_data.txt"))

pdf(file.path(output_dir, "umap_overview.pdf"), width = 9, height = 7)
print(Seurat::DimPlot(integrated, reduction = "umap", group.by = "sample"))
print(Seurat::DimPlot(integrated, reduction = "umap", group.by = "seurat_clusters", label = TRUE))
dev.off()

message("Writing stack figures ...")
write_stackfig_outputs(
  obj = integrated,
  output_dir = output_dir
)

message("\nDone. Output directory: ", output_dir)
