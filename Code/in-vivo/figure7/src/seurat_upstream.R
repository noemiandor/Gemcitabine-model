# Minimal Seurat reconstruction helpers used only by Figure 7 and
# Supplementary Figures 4-7.
#
# These functions preserve the scientific mutations performed by Tao's
# published 01_data.R, 01a_cell_cycle.R, 02b_cluster_refine.R,
# 02d_manual_cluster_merge.R, and 03_final_cluster.R scripts.  Diagnostic
# plots, marker scans, and unrelated exports from those scripts are omitted.

figure7_upstream_require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    figure7_stop("R package '", package, "' is required for Seurat reconstruction")
  }
}

figure7_upstream_require_columns <- function(data, columns, label) {
  missing <- setdiff(columns, colnames(data))
  if (length(missing)) {
    figure7_stop(
      label,
      " is missing column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  invisible(data)
}

figure7_upstream_sort_maybe_numeric <- function(values) {
  values <- as.character(values)
  observed <- unique(values)
  numeric_values <- suppressWarnings(as.numeric(observed))
  if (!anyNA(numeric_values)) {
    observed[order(numeric_values)]
  } else {
    sort(observed)
  }
}

figure7_upstream_set_single_thread <- function() {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1",
    KMP_DUPLICATE_LIB_OK = "TRUE",
    KMP_INIT_AT_FORK = "FALSE"
  )
  invisible(TRUE)
}

figure7_upstream_configure_future <- function(
  jobs = 1L,
  max_size_gb = 60
) {
  jobs <- suppressWarnings(as.integer(jobs))
  if (length(jobs) != 1L || is.na(jobs) || jobs < 1L) {
    figure7_stop("Seurat reconstruction jobs must be a positive integer")
  }
  max_size_bytes <- as.numeric(max_size_gb) * 1024^3
  current_max <- getOption("future.globals.maxSize")
  if (is.null(current_max) || is.na(current_max) ||
      current_max < max_size_bytes) {
    options(future.globals.maxSize = max_size_bytes)
  }
  figure7_upstream_require_namespace("future")
  future::plan(future::sequential)
  invisible(list(
    strategy = "sequential",
    workers = 1L,
    requested_workers = jobs
  ))
}

figure7_upstream_read_sample_info <- function(path) {
  if (!file.exists(path)) {
    figure7_stop("Missing Seurat reconstruction sample workbook: ", path)
  }
  if (requireNamespace("readxl", quietly = TRUE)) {
    return(as.data.frame(
      readxl::read_excel(path),
      stringsAsFactors = FALSE,
      check.names = FALSE
    ))
  }
  if (!requireNamespace("xml2", quietly = TRUE)) {
    figure7_stop(
      "R package 'readxl' or 'xml2' is required to read sample_info.xlsx"
    )
  }

  temporary <- tempfile("figure7_sample_info_")
  dir.create(temporary, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  utils::unzip(
    path,
    files = c("xl/sharedStrings.xml", "xl/worksheets/sheet1.xml"),
    exdir = temporary
  )
  shared_path <- file.path(temporary, "xl", "sharedStrings.xml")
  sheet_path <- file.path(temporary, "xl", "worksheets", "sheet1.xml")
  if (!file.exists(sheet_path)) {
    figure7_stop("sample_info.xlsx does not contain sheet1.xml")
  }
  namespace <- c(
    x = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
  )
  shared_strings <- character()
  if (file.exists(shared_path)) {
    document <- xml2::read_xml(shared_path)
    nodes <- xml2::xml_find_all(document, ".//x:si", ns = namespace)
    shared_strings <- vapply(nodes, function(node) {
      paste(
        xml2::xml_text(
          xml2::xml_find_all(node, ".//x:t", ns = namespace)
        ),
        collapse = ""
      )
    }, character(1L))
  }
  column_index <- function(reference) {
    letters <- gsub("[0-9]+", "", reference)
    values <- utf8ToInt(letters) - utf8ToInt("A") + 1L
    sum(values * (26L ^ rev(seq_along(values) - 1L)))
  }
  cell_value <- function(node) {
    type <- xml2::xml_attr(node, "t")
    value_node <- xml2::xml_find_first(node, "x:v", ns = namespace)
    if (inherits(value_node, "xml_missing")) return(NA_character_)
    value <- xml2::xml_text(value_node)
    if (!is.na(type) && type == "s") {
      index <- suppressWarnings(as.integer(value)) + 1L
      if (is.na(index) || index < 1L || index > length(shared_strings)) {
        return(NA_character_)
      }
      return(shared_strings[[index]])
    }
    value
  }
  document <- xml2::read_xml(sheet_path)
  row_nodes <- xml2::xml_find_all(
    document,
    ".//x:sheetData/x:row",
    ns = namespace
  )
  if (length(row_nodes) < 2L) {
    figure7_stop("sample_info.xlsx has no data rows")
  }
  parsed <- lapply(row_nodes, function(row_node) {
    cells <- xml2::xml_find_all(row_node, "x:c", ns = namespace)
    if (!length(cells)) return(character())
    indexes <- vapply(
      xml2::xml_attr(cells, "r"),
      column_index,
      integer(1L)
    )
    values <- vapply(cells, cell_value, character(1L))
    output <- rep(NA_character_, max(indexes))
    output[indexes] <- values
    output
  })
  width <- max(vapply(parsed, length, integer(1L)))
  matrix <- do.call(rbind, lapply(parsed, function(row) {
    length(row) <- width
    row
  }))
  header <- matrix[1L, ]
  keep <- !is.na(header) & nzchar(header)
  output <- as.data.frame(
    matrix[-1L, keep, drop = FALSE],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(output) <- header[keep]
  output
}

figure7_upstream_resolve_column <- function(
  data,
  candidates,
  fallback_index = NULL
) {
  for (candidate in candidates) {
    match <- which(tolower(names(data)) == tolower(candidate))
    if (length(match) == 1L) return(names(data)[match])
  }
  if (!is.null(fallback_index) &&
      fallback_index >= 1L &&
      fallback_index <= ncol(data)) {
    return(names(data)[fallback_index])
  }
  figure7_stop(
    "Cannot resolve upstream metadata column: ",
    paste(candidates, collapse = ", ")
  )
}

figure7_upstream_expected_samples <- function(sample_info) {
  id_column <- figure7_upstream_resolve_column(
    sample_info,
    c("IDs", "ID"),
    3L
  )
  ids <- as.character(sample_info[[id_column]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    figure7_stop("sample_info.xlsx contains invalid or duplicate sample IDs")
  }
  sort(c(ids, "2N-Cell-Culture", "4N-Cell-Culture"))
}

figure7_upstream_h5_inventory <- function(
  cellranger_root,
  sample_info
) {
  if (!dir.exists(cellranger_root)) {
    figure7_stop("Missing Cell Ranger root: ", cellranger_root)
  }
  expected_samples <- figure7_upstream_expected_samples(sample_info)
  sample_directories <- sort(
    list.dirs(cellranger_root, recursive = FALSE, full.names = TRUE)
  )
  sample_directories <- sample_directories[
    grepl("-Count-HM$", basename(sample_directories))
  ]
  observed_samples <- sub("-Count-HM$", "", basename(sample_directories))
  missing <- setdiff(expected_samples, observed_samples)
  unexpected <- setdiff(observed_samples, expected_samples)
  if (length(missing) || length(unexpected) ||
      anyDuplicated(observed_samples)) {
    figure7_stop(
      "Cell Ranger sample inventory differs from the reviewed 18-sample ",
      "contract; missing=[",
      paste(missing, collapse = ","),
      "]; unexpected=[",
      paste(unexpected, collapse = ","),
      "]"
    )
  }
  paths <- file.path(
    sample_directories,
    "outs",
    "filtered_feature_bc_matrix.h5"
  )
  absent <- paths[!file.exists(paths)]
  if (length(absent)) {
    figure7_stop(
      "Missing Cell Ranger filtered matrix file(s): ",
      paste(absent, collapse = ", ")
    )
  }
  data.frame(
    sample = observed_samples,
    path = normalizePath(paths, mustWork = TRUE),
    sha256 = vapply(paths, figure7_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
}

figure7_upstream_read_counts <- function(path) {
  figure7_upstream_require_namespace("Seurat")
  counts <- Seurat::Read10X_h5(path)
  if (is.list(counts)) {
    counts <- if ("Gene Expression" %in% names(counts)) {
      counts[["Gene Expression"]]
    } else {
      counts[[1L]]
    }
  }
  counts
}

figure7_upstream_qc_special <- function(
  object,
  min_features = 200,
  min_counts = 500,
  max_percent_mt = 20,
  upper_quantile = 0.99
) {
  feature_upper <- as.numeric(stats::quantile(
    object$nFeature_RNA,
    probs = upper_quantile,
    na.rm = TRUE
  ))
  count_upper <- as.numeric(stats::quantile(
    object$nCount_RNA,
    probs = upper_quantile,
    na.rm = TRUE
  ))
  keep <- object$nFeature_RNA >= min_features &
    object$nCount_RNA >= min_counts &
    object$nFeature_RNA <= feature_upper &
    object$nCount_RNA <= count_upper &
    object$percent.mt <= max_percent_mt
  subset(object, cells = colnames(object)[keep])
}

figure7_upstream_remove_doublets <- function(object) {
  if (!ncol(object)) return(object)
  if (ncol(object) < 100L) {
    object$scDblFinder.class <- "singlet_low_cell_skip"
    object$scDblFinder.score <- NA_real_
    return(object)
  }
  figure7_upstream_require_namespace("SingleCellExperiment")
  figure7_upstream_require_namespace("scDblFinder")
  single_cell <- SingleCellExperiment::SingleCellExperiment(
    assays = list(
      counts = Seurat::GetAssayData(
        object,
        assay = "RNA",
        slot = "counts"
      )
    )
  )
  set.seed(1234)
  single_cell <- scDblFinder::scDblFinder(single_cell, verbose = FALSE)
  metadata <- as.data.frame(SingleCellExperiment::colData(single_cell))
  object$scDblFinder.class <- as.character(metadata$scDblFinder.class)
  object$scDblFinder.score <- as.numeric(metadata$scDblFinder.score)
  singlets <- rownames(metadata)[metadata$scDblFinder.class == "singlet"]
  subset(object, cells = singlets)
}

figure7_upstream_build_integrated <- function(
  cellranger_root,
  all_ploidy_path,
  sample_info_path,
  jobs = 1L
) {
  for (package in c("Seurat", "SingleCellExperiment", "scDblFinder")) {
    figure7_upstream_require_namespace(package)
  }
  sample_info <- figure7_upstream_read_sample_info(sample_info_path)
  if (ncol(sample_info) < 3L) {
    figure7_stop("sample_info.xlsx must contain at least three columns")
  }
  harvest_column <- figure7_upstream_resolve_column(
    sample_info,
    "harvest",
    1L
  )
  sequencing_column <- figure7_upstream_resolve_column(
    sample_info,
    c("Sequencing IDs", "Sequencing ID", "SequencingIDs"),
    2L
  )
  id_column <- figure7_upstream_resolve_column(
    sample_info,
    c("IDs", "ID"),
    3L
  )
  sample_info[[harvest_column]] <- as.character(
    sample_info[[harvest_column]]
  )
  sample_info[[id_column]] <- as.character(sample_info[[id_column]])

  ploidy <- figure7_read_tsv(
    all_ploidy_path,
    c("file", "cell_id"),
    "all_ploidy.tsv"
  )
  ploidy$harvest <- sub("\\.sps\\.cbs$", "", ploidy$file)
  barcode_by_harvest <- lapply(
    split(ploidy$cell_id, ploidy$harvest),
    unique
  )
  harvest_by_id <- stats::setNames(
    sample_info[[harvest_column]],
    sample_info[[id_column]]
  )
  barcode_by_id <- lapply(harvest_by_id, function(harvest) {
    if (is.na(harvest) || !harvest %in% names(barcode_by_harvest)) {
      character()
    } else {
      barcode_by_harvest[[harvest]]
    }
  })

  inventory <- figure7_upstream_h5_inventory(
    cellranger_root,
    sample_info
  )
  special_samples <- c("2N-Cell-Culture", "4N-Cell-Culture")
  objects <- list()
  summaries <- list()
  set.seed(1234)
  for (index in seq_len(nrow(inventory))) {
    sample_id <- inventory$sample[[index]]
    counts <- figure7_upstream_read_counts(inventory$path[[index]])
    object <- Seurat::CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = 0,
      min.features = 0
    )
    object$sample <- sample_id
    object$sample_folder <- paste0(sample_id, "-Count-HM")
    object[["percent.mt"]] <- Seurat::PercentageFeatureSet(
      object,
      pattern = "^MT-"
    )
    raw_cells <- ncol(object)
    if (sample_id %in% special_samples) {
      method <- "QC + scDblFinder"
      object <- figure7_upstream_qc_special(object)
      after_qc <- ncol(object)
      object <- figure7_upstream_remove_doublets(object)
    } else {
      method <- "all_ploidy keep barcodes"
      keep_barcodes <- barcode_by_id[[sample_id]]
      if (is.null(keep_barcodes) || !length(keep_barcodes)) {
        figure7_stop(
          "No endpoint-ploidy barcode mapping for sample: ",
          sample_id
        )
      }
      keep_cells <- intersect(colnames(object), keep_barcodes)
      if (!length(keep_cells)) {
        figure7_stop(
          "No Cell Ranger barcodes survived endpoint-ploidy filtering for ",
          sample_id
        )
      }
      object <- subset(object, cells = keep_cells)
      object$scDblFinder.class <- "singlet_by_ploidy"
      object$scDblFinder.score <- NA_real_
      after_qc <- ncol(object)
    }
    if (!ncol(object)) {
      figure7_stop("No cells remain after filtering sample ", sample_id)
    }
    if (sample_id %in% special_samples) {
      metadata <- as.list(stats::setNames(
        rep(NA_character_, ncol(sample_info)),
        names(sample_info)
      ))
      metadata[[sequencing_column]] <- sample_id
      metadata[[id_column]] <- sample_id
    } else {
      match <- which(sample_info[[id_column]] == sample_id)
      if (length(match) != 1L) {
        figure7_stop(
          "sample_info.xlsx must contain exactly one row for ",
          sample_id
        )
      }
      metadata <- as.list(sample_info[match, , drop = FALSE])
    }
    for (column in names(metadata)) {
      object[[column]] <- metadata[[column]]
    }
    object$barcode_raw <- colnames(object)
    object <- Seurat::RenameCells(
      object,
      new.names = paste0(sample_id, "_", colnames(object))
    )
    objects[[sample_id]] <- object
    summaries[[sample_id]] <- data.frame(
      sample = sample_id,
      sample_folder = paste0(sample_id, "-Count-HM"),
      filter_method = method,
      raw_cells = raw_cells,
      after_qc_cells = after_qc,
      after_filter_cells = ncol(object),
      stringsAsFactors = FALSE
    )
  }
  if (length(objects) != 18L) {
    figure7_stop(
      "Expected 18 non-empty samples after filtering; observed ",
      length(objects)
    )
  }

  objects <- lapply(objects, function(object) {
    Seurat::SCTransform(
      object,
      assay = "RNA",
      new.assay.name = "SCT",
      vars.to.regress = "percent.mt",
      return.only.var.genes = FALSE,
      verbose = FALSE
    )
  })
  features <- Seurat::SelectIntegrationFeatures(
    object.list = objects,
    nfeatures = 3000
  )
  figure7_upstream_configure_future(jobs, max_size_gb = 60)
  invisible(gc())
  objects <- Seurat::PrepSCTIntegration(
    object.list = objects,
    anchor.features = features,
    verbose = FALSE
  )
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = objects,
    normalization.method = "SCT",
    anchor.features = features,
    reduction = "cca",
    verbose = FALSE
  )
  integrated <- Seurat::IntegrateData(
    anchorset = anchors,
    normalization.method = "SCT",
    verbose = FALSE
  )
  Seurat::DefaultAssay(integrated) <- "integrated"
  integrated <- Seurat::RunPCA(
    integrated,
    npcs = 50L,
    verbose = FALSE
  )
  dimensions <- seq_len(min(
    30L,
    ncol(Seurat::Embeddings(integrated, reduction = "pca"))
  ))
  integrated <- Seurat::FindNeighbors(
    integrated,
    dims = dimensions,
    verbose = FALSE
  )
  integrated <- Seurat::FindClusters(
    integrated,
    resolution = 0.6,
    verbose = FALSE
  )
  integrated <- Seurat::RunUMAP(
    integrated,
    dims = dimensions,
    verbose = FALSE
  )
  list(
    object = integrated,
    sample_summary = do.call(rbind, summaries),
    h5_inventory = inventory
  )
}

figure7_upstream_normalize_feature_key <- function(values) {
  values <- trimws(as.character(values))
  values <- sub("^GRCh[0-9]+-", "", values, ignore.case = TRUE)
  values <- sub("\\.[0-9]+$", "", values)
  toupper(values)
}

figure7_upstream_feature_metadata <- function(assay) {
  slots <- methods::slotNames(assay)
  if ("meta.features" %in% slots) {
    return(as.data.frame(assay@meta.features, stringsAsFactors = FALSE))
  }
  if ("meta.data" %in% slots) {
    return(as.data.frame(assay@meta.data, stringsAsFactors = FALSE))
  }
  data.frame(row.names = rownames(assay))
}

figure7_upstream_feature_lookup <- function(feature_ids, aliases) {
  keys <- figure7_upstream_normalize_feature_key(aliases)
  keep <- !is.na(keys) & nzchar(keys) & !duplicated(keys)
  stats::setNames(feature_ids[keep], keys[keep])
}

figure7_upstream_resolve_gene_set <- function(
  object,
  assay_name,
  symbols,
  label,
  minimum = 5L
) {
  assay <- object[[assay_name]]
  feature_ids <- rownames(assay)
  if (!length(feature_ids)) {
    figure7_stop("No features are present in assay ", assay_name)
  }
  keys <- figure7_upstream_normalize_feature_key(symbols)
  matches <- rep(NA_character_, length(symbols))
  lookups <- list(
    feature_ids = figure7_upstream_feature_lookup(
      feature_ids,
      feature_ids
    )
  )
  if (any(grepl("\\|", feature_ids))) {
    lookups$before_pipe <- figure7_upstream_feature_lookup(
      feature_ids,
      sub("\\|.*$", "", feature_ids)
    )
    lookups$after_pipe <- figure7_upstream_feature_lookup(
      feature_ids,
      sub("^.*\\|", "", feature_ids)
    )
  }
  if (any(grepl("_", feature_ids))) {
    lookups$before_underscore <- figure7_upstream_feature_lookup(
      feature_ids,
      sub("_.*$", "", feature_ids)
    )
    lookups$after_underscore <- figure7_upstream_feature_lookup(
      feature_ids,
      sub("^.*_", "", feature_ids)
    )
  }
  feature_metadata <- figure7_upstream_feature_metadata(assay)
  metadata_columns <- intersect(
    c(
      "gene_name", "gene", "symbol", "gene_symbol", "gene_symbols",
      "SYMBOL", "GENE", "Gene", "GeneName", "GeneSymbol",
      "feature_name", "feature", "features"
    ),
    colnames(feature_metadata)
  )
  for (column in metadata_columns) {
    lookups[[paste0("meta_", column)]] <-
      figure7_upstream_feature_lookup(
        feature_ids,
        feature_metadata[[column]]
      )
  }
  for (lookup in lookups) {
    proposed <- unname(lookup[keys])
    fill <- is.na(matches) & !is.na(proposed)
    matches[fill] <- proposed[fill]
  }
  features <- unique(matches[!is.na(matches)])
  if (length(features) < minimum) {
    figure7_stop(
      label,
      " has only ",
      length(features),
      " matched features in assay ",
      assay_name
    )
  }
  features
}

figure7_upstream_safe_sd <- function(values) {
  output <- stats::sd(values, na.rm = TRUE)
  if (is.na(output)) 0 else output
}

figure7_upstream_cell_cycle_candidate_table <- function(
  metadata,
  pca_embeddings
) {
  figure7_upstream_require_columns(
    metadata,
    c("seurat_clusters", "S.Score", "G2M.Score", "Phase"),
    "Cell-cycle metadata"
  )
  if (nrow(pca_embeddings) != nrow(metadata)) {
    figure7_stop("PCA embeddings do not align with cell-cycle metadata")
  }
  if (ncol(pca_embeddings) < 20L) {
    figure7_stop("At least 20 PCs are required for cell-cycle annotation")
  }
  clusters <- as.character(metadata$seurat_clusters)
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    figure7_stop("Cell-cycle metadata contains missing cluster labels")
  }
  cluster_levels <- unique(clusters)
  summaries <- do.call(rbind, lapply(cluster_levels, function(cluster) {
    selected <- clusters == cluster
    data.frame(
      cluster = cluster,
      n_cells = sum(selected),
      mean_S.Score = mean(metadata$S.Score[selected], na.rm = TRUE),
      mean_G2M.Score = mean(metadata$G2M.Score[selected], na.rm = TRUE),
      frac_S = mean(as.character(metadata$Phase[selected]) == "S",
                    na.rm = TRUE),
      frac_G2M = mean(as.character(metadata$Phase[selected]) == "G2M",
                      na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  score_threshold_s <- mean(
    summaries$mean_S.Score,
    na.rm = TRUE
  ) + figure7_upstream_safe_sd(summaries$mean_S.Score)
  score_threshold_g2m <- mean(
    summaries$mean_G2M.Score,
    na.rm = TRUE
  ) + figure7_upstream_safe_sd(summaries$mean_G2M.Score)
  score_flag <- summaries$mean_S.Score > score_threshold_s |
    summaries$mean_G2M.Score > score_threshold_g2m

  phase_threshold_s <- mean(
    summaries$frac_S,
    na.rm = TRUE
  ) + figure7_upstream_safe_sd(summaries$frac_S)
  phase_threshold_g2m <- mean(
    summaries$frac_G2M,
    na.rm = TRUE
  ) + figure7_upstream_safe_sd(summaries$frac_G2M)
  phase_flag <- summaries$frac_S > phase_threshold_s |
    summaries$frac_G2M > phase_threshold_g2m

  pca_use <- pca_embeddings[, seq_len(20L), drop = FALSE]
  if (is.null(colnames(pca_use)) || any(!nzchar(colnames(pca_use)))) {
    colnames(pca_use) <- paste0("PC", seq_len(ncol(pca_use)))
  }
  cor_s <- apply(pca_use, 2L, function(values) {
    stats::cor(
      values,
      metadata$S.Score,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  })
  cor_g2m <- apply(pca_use, 2L, function(values) {
    stats::cor(
      values,
      metadata$G2M.Score,
      method = "pearson",
      use = "pairwise.complete.obs"
    )
  })
  selected_pc <- which(abs(cor_s) >= 0.30 | abs(cor_g2m) >= 0.30)
  if (!length(selected_pc)) {
    selection_score <- pmax(abs(cor_s), abs(cor_g2m))
    selected_pc <- which.max(selection_score)
  }
  pc_flag <- rep(FALSE, length(cluster_levels))
  for (pc in selected_pc) {
    cluster_means <- vapply(cluster_levels, function(cluster) {
      mean(pca_use[clusters == cluster, pc], na.rm = TRUE)
    }, numeric(1L))
    across_mean <- mean(cluster_means, na.rm = TRUE)
    across_sd <- figure7_upstream_safe_sd(cluster_means)
    pc_flag <- pc_flag |
      abs(cluster_means - across_mean) > across_sd
  }
  method_count <- as.integer(score_flag) +
    as.integer(phase_flag) +
    as.integer(pc_flag)
  data.frame(
    cluster = cluster_levels,
    flag_score_outlier = score_flag,
    flag_phase_enrichment = phase_flag,
    flag_pc_separation_support = pc_flag,
    n_methods_flagged = method_count,
    final_candidate_flag = method_count >= 2L,
    stringsAsFactors = FALSE
  )
}

figure7_upstream_annotate_cell_cycle <- function(object) {
  set.seed(12345)
  figure7_upstream_require_namespace("Seurat")
  figure7_upstream_require_columns(
    object@meta.data,
    c("seurat_clusters", "sample", "Dose", "orig.ident"),
    "Integrated Seurat metadata"
  )
  if (!"RNA" %in% names(object@assays)) {
    figure7_stop("Integrated Seurat object is missing RNA assay")
  }
  if (!all(c("pca", "umap") %in% names(object@reductions))) {
    figure7_stop("Integrated Seurat object requires PCA and UMAP reductions")
  }
  Seurat::DefaultAssay(object) <- "RNA"
  genes <- Seurat::cc.genes.updated.2019
  s_features <- figure7_upstream_resolve_gene_set(
    object,
    "RNA",
    genes$s.genes,
    "S.Score"
  )
  g2m_features <- figure7_upstream_resolve_gene_set(
    object,
    "RNA",
    genes$g2m.genes,
    "G2M.Score"
  )
  object <- Seurat::CellCycleScoring(
    object,
    s.features = s_features,
    g2m.features = g2m_features,
    set.ident = FALSE
  )
  candidates <- figure7_upstream_cell_cycle_candidate_table(
    object@meta.data,
    Seurat::Embeddings(object, reduction = "pca")
  )
  annotation <- ifelse(
    candidates$final_candidate_flag,
    "cell_cycle_candidate",
    "not_cell_cycle_candidate"
  )
  names(annotation) <- candidates$cluster
  object$cluster_cell_cycle_annotation <- unname(
    annotation[as.character(object$seurat_clusters)]
  )
  if (anyNA(object$cluster_cell_cycle_annotation)) {
    figure7_stop("Cell-cycle annotation could not be assigned to every cell")
  }
  list(object = object, candidate_table = candidates)
}

figure7_upstream_build_refined_levels <- function(
  original_levels,
  label_map
) {
  output <- character()
  for (cluster in original_levels) {
    output <- c(output, cluster)
    if (cluster %in% names(label_map)) {
      output <- c(output, unname(label_map[[cluster]]))
    }
  }
  unique(output)
}

figure7_upstream_point_in_polygon <- function(
  x,
  y,
  polygon_x,
  polygon_y
) {
  n <- length(polygon_x)
  inside <- rep(FALSE, length(x))
  j <- n
  for (i in seq_len(n)) {
    xi <- polygon_x[[i]]
    yi <- polygon_y[[i]]
    xj <- polygon_x[[j]]
    yj <- polygon_y[[j]]
    crosses <- ((yi > y) != (yj > y)) &
      (x < ((xj - xi) * (y - yi) / (yj - yi + 1e-12) + xi))
    inside <- xor(inside, crosses)
    j <- i
  }
  inside
}

figure7_upstream_knn_indices <- function(coordinates, k) {
  figure7_upstream_require_namespace("RANN")
  n <- nrow(coordinates)
  if (n < 2L) {
    figure7_stop("At least two cells are required for UMAP kNN")
  }
  k <- min(as.integer(k), n - 1L)
  nearest <- RANN::nn2(
    data = coordinates,
    query = coordinates,
    k = k + 1L
  )
  nearest$nn.idx[, -1L, drop = FALSE]
}

figure7_upstream_knn_components <- function(
  coordinates,
  k,
  mutual = FALSE
) {
  n <- nrow(coordinates)
  if (!n) return(integer())
  if (n == 1L) return(1L)
  k <- min(as.integer(k), n - 1L)
  if (!is.finite(k) || k < 1L) return(seq_len(n))
  nearest <- figure7_upstream_knn_indices(coordinates, k)
  if (!is.matrix(nearest)) {
    nearest <- matrix(nearest, nrow = n, ncol = k)
  }
  adjacency <- replicate(n, integer(), simplify = FALSE)
  if (isTRUE(mutual)) {
    sets <- lapply(seq_len(n), function(index) {
      unique(as.integer(nearest[index, ]))
    })
    for (i in seq_len(n)) {
      for (j in sets[[i]]) {
        if (is.na(j) || j < 1L || j > n) next
        if (i %in% sets[[j]]) {
          adjacency[[i]] <- c(adjacency[[i]], j)
          adjacency[[j]] <- c(adjacency[[j]], i)
        }
      }
    }
  } else {
    for (i in seq_len(n)) {
      neighbors <- unique(as.integer(nearest[i, ]))
      neighbors <- neighbors[
        !is.na(neighbors) & neighbors >= 1L & neighbors <= n
      ]
      if (!length(neighbors)) next
      adjacency[[i]] <- c(adjacency[[i]], neighbors)
      for (j in neighbors) {
        adjacency[[j]] <- c(adjacency[[j]], i)
      }
    }
  }
  adjacency <- lapply(adjacency, unique)
  component <- integer(n)
  component_id <- 0L
  for (i in seq_len(n)) {
    if (component[[i]] != 0L) next
    component_id <- component_id + 1L
    queue <- i
    component[[i]] <- component_id
    head <- 1L
    while (head <= length(queue)) {
      node <- queue[[head]]
      head <- head + 1L
      neighbors <- adjacency[[node]]
      if (!length(neighbors)) next
      new_nodes <- neighbors[component[neighbors] == 0L]
      if (!length(new_nodes)) next
      component[new_nodes] <- component_id
      queue <- c(queue, new_nodes)
    }
  }
  component
}

figure7_upstream_expand_candidate_component <- function(
  candidate_index,
  seed_index,
  coordinates,
  fraction_region,
  fraction_core,
  inside_hull,
  minimum_region_fraction = 0,
  minimum_core_fraction = 0,
  require_inside_hull = FALSE,
  graph_k = 5L,
  mutual_knn = FALSE
) {
  if (!length(candidate_index) || !length(seed_index)) {
    return(list(selected_index = integer(), eligible_index = integer()))
  }
  eligible <- candidate_index[
    fraction_region[candidate_index] >= minimum_region_fraction &
      fraction_core[candidate_index] >= minimum_core_fraction &
      (!require_inside_hull | inside_hull[candidate_index])
  ]
  if (!length(eligible)) {
    return(list(selected_index = integer(), eligible_index = integer()))
  }
  seed_in_eligible <- eligible %in% seed_index
  if (!any(seed_in_eligible)) {
    return(list(selected_index = integer(), eligible_index = eligible))
  }
  components <- figure7_upstream_knn_components(
    coordinates[eligible, , drop = FALSE],
    graph_k,
    mutual_knn
  )
  selected_components <- unique(components[seed_in_eligible])
  list(
    selected_index = eligible[components %in% selected_components],
    eligible_index = eligible
  )
}

figure7_upstream_refine_cluster_values <- function(
  cluster_values,
  coordinates,
  knn_k = 50L,
  core_clusters = c("6", "10", "11", "12"),
  neighbor_clusters = c("0", "3"),
  candidate_label_map = c("4" = "4c", "9" = "9c")
) {
  cluster_values <- as.character(cluster_values)
  if (length(cluster_values) != nrow(coordinates)) {
    figure7_stop("Cluster labels do not align with UMAP coordinates")
  }
  required_clusters <- unique(c(
    core_clusters,
    neighbor_clusters,
    names(candidate_label_map)
  ))
  missing <- setdiff(required_clusters, unique(cluster_values))
  if (length(missing)) {
    figure7_stop(
      "Cluster refinement is missing configured cluster(s): ",
      paste(missing, collapse = ", ")
    )
  }
  region_clusters <- c(core_clusters, neighbor_clusters)
  nearest <- figure7_upstream_knn_indices(coordinates, knn_k)
  neighbor_labels <- matrix(
    cluster_values[nearest],
    nrow = nrow(nearest),
    ncol = ncol(nearest)
  )
  region_matrix <- matrix(
    neighbor_labels %in% region_clusters,
    nrow = nrow(neighbor_labels),
    ncol = ncol(neighbor_labels)
  )
  core_matrix <- matrix(
    neighbor_labels %in% core_clusters,
    nrow = nrow(neighbor_labels),
    ncol = ncol(neighbor_labels)
  )
  fraction_region <- rowMeans(region_matrix, na.rm = TRUE)
  fraction_core <- rowMeans(core_matrix, na.rm = TRUE)
  reference <- coordinates[
    cluster_values %in% region_clusters,
    seq_len(2L),
    drop = FALSE
  ]
  hull_order <- grDevices::chull(reference[, 1L], reference[, 2L])
  hull <- reference[hull_order, , drop = FALSE]
  hull <- rbind(hull, hull[1L, , drop = FALSE])
  inside_hull <- figure7_upstream_point_in_polygon(
    coordinates[, 1L],
    coordinates[, 2L],
    hull[, 1L],
    hull[, 2L]
  )
  candidate_index <- which(cluster_values %in% names(candidate_label_map))
  seed_index <- candidate_index[
    fraction_region[candidate_index] >= 0.60 &
      fraction_core[candidate_index] >= 0.20 &
      inside_hull[candidate_index]
  ]
  expanded_index <- integer()
  eligible_index <- integer()
  for (cluster in names(candidate_label_map)) {
    cluster_index <- which(cluster_values == cluster)
    cluster_seed <- seed_index[cluster_values[seed_index] == cluster]
    expansion <- figure7_upstream_expand_candidate_component(
      candidate_index = cluster_index,
      seed_index = cluster_seed,
      coordinates = coordinates,
      fraction_region = fraction_region,
      fraction_core = fraction_core,
      inside_hull = inside_hull,
      minimum_region_fraction = 0,
      minimum_core_fraction = 0,
      require_inside_hull = FALSE,
      graph_k = 5L,
      mutual_knn = FALSE
    )
    expanded_index <- union(
      expanded_index,
      expansion$selected_index
    )
    eligible_index <- union(
      eligible_index,
      expansion$eligible_index
    )
  }
  selected_index <- sort(unique(c(seed_index, expanded_index)))
  refined <- cluster_values
  refined[selected_index] <- unname(
    candidate_label_map[cluster_values[selected_index]]
  )
  levels <- figure7_upstream_build_refined_levels(
    figure7_upstream_sort_maybe_numeric(cluster_values),
    candidate_label_map
  )
  list(
    refined = factor(refined, levels = levels),
    seed_index = seed_index,
    selected_index = selected_index,
    eligible_index = eligible_index,
    fraction_region = fraction_region,
    fraction_core = fraction_core,
    inside_hull = inside_hull
  )
}

figure7_upstream_refine_clusters <- function(object) {
  figure7_upstream_set_single_thread()
  set.seed(12345)
  figure7_upstream_require_namespace("Seurat")
  figure7_upstream_require_columns(
    object@meta.data,
    "seurat_clusters",
    "Integrated Seurat metadata"
  )
  if (!"umap" %in% names(object@reductions)) {
    figure7_stop("Integrated Seurat object is missing UMAP reduction")
  }
  coordinates <- Seurat::Embeddings(object, reduction = "umap")
  if (ncol(coordinates) < 2L) {
    figure7_stop("UMAP reduction requires at least two dimensions")
  }
  refinement <- figure7_upstream_refine_cluster_values(
    object$seurat_clusters,
    coordinates[, seq_len(2L), drop = FALSE]
  )
  object$seurat_cluster_refine <- refinement$refined
  list(object = object, audit = refinement)
}

figure7_upstream_manual_merge_values <- function(
  cluster_values,
  cluster_order = NULL,
  merge_map = list(
    `0` = c("0", "1", "7"),
    `10` = c("10", "11", "12")
  )
) {
  cluster_values <- as.character(cluster_values)
  output <- cluster_values
  members <- unlist(merge_map, use.names = FALSE)
  if (anyDuplicated(members)) {
    figure7_stop("Manual cluster merge map contains overlapping members")
  }
  missing <- setdiff(members, unique(cluster_values))
  if (length(missing)) {
    figure7_stop(
      "Manual cluster merge is missing configured cluster(s): ",
      paste(missing, collapse = ", ")
    )
  }
  for (label in names(merge_map)) {
    output[cluster_values %in% merge_map[[label]]] <- label
  }
  if (is.null(cluster_order)) {
    cluster_order <- figure7_upstream_sort_maybe_numeric(cluster_values)
  }
  levels <- character()
  for (cluster in as.character(cluster_order)) {
    target <- if (cluster %in% members) {
      names(merge_map)[vapply(
        merge_map,
        function(group) cluster %in% group,
        logical(1L)
      )][1L]
    } else {
      cluster
    }
    if (!target %in% levels) levels <- c(levels, target)
  }
  factor(output, levels = levels)
}

figure7_upstream_join_layers <- function(object, assay = "RNA") {
  if (!"JoinLayers" %in% getNamespaceExports("Seurat") ||
      !assay %in% names(object@assays)) {
    return(object)
  }
  tryCatch(
    Seurat::JoinLayers(object, assay = assay),
    error = function(error) object
  )
}

figure7_upstream_assay_data <- function(
  object,
  assay = "RNA",
  slot = "data"
) {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, slot = slot),
    error = function(first_error) {
      tryCatch(
        Seurat::GetAssayData(object, assay = assay, layer = slot),
        error = function(second_error) NULL
      )
    }
  )
}

figure7_upstream_merge_clusters <- function(object, jobs = 1L) {
  figure7_upstream_set_single_thread()
  set.seed(12345)
  figure7_upstream_configure_future(jobs, max_size_gb = 30)
  figure7_upstream_require_namespace("Seurat")
  figure7_upstream_require_columns(
    object@meta.data,
    "seurat_cluster_refine",
    "Refined Seurat metadata"
  )
  if (!"RNA" %in% names(object@assays)) {
    figure7_stop("Refined Seurat object is missing RNA assay")
  }
  Seurat::DefaultAssay(object) <- "RNA"
  object <- figure7_upstream_join_layers(object, "RNA")
  data <- figure7_upstream_assay_data(object, "RNA", "data")
  if (is.null(data) || !nrow(data) || !ncol(data)) {
    object <- Seurat::NormalizeData(
      object,
      assay = "RNA",
      verbose = FALSE
    )
  }
  cluster_values <- as.character(object$seurat_cluster_refine)
  cluster_order <- if (is.factor(object$seurat_cluster_refine)) {
    levels(object$seurat_cluster_refine)
  } else {
    figure7_upstream_sort_maybe_numeric(cluster_values)
  }
  object$seurat_cluster_refine <- factor(
    cluster_values,
    levels = cluster_order
  )
  object$manual_merge_test <- figure7_upstream_manual_merge_values(
    cluster_values,
    cluster_order
  )
  object
}

figure7_upstream_derive_tn_ploidy <- function(metadata) {
  figure7_upstream_require_columns(
    metadata,
    "IDs",
    "Final Seurat metadata"
  )
  ids <- as.character(metadata$IDs)
  ids[is.na(ids)] <- ""
  ploidy <- ifelse(
    grepl("2N", ids, fixed = TRUE),
    "2N",
    ifelse(grepl("4N", ids, fixed = TRUE), "4N", NA_character_)
  )
  if (anyNA(ploidy)) {
    unresolved <- sort(unique(ids[is.na(ploidy)]))
    figure7_stop(
      "Cannot derive 2N/4N ploidy from IDs: ",
      paste(utils::head(unresolved, 20L), collapse = ", ")
    )
  }
  context <- ifelse(
    grepl("Cell-Culture", ids, fixed = TRUE),
    "CellLine",
    "Tumor"
  )
  metadata$Ploidy <- factor(ploidy, levels = c("2N", "4N"))
  metadata$TN <- factor(context, levels = c("CellLine", "Tumor"))
  metadata
}

figure7_upstream_choose_pca_assay <- function(
  object,
  preferred = c("integrated", "SCT", "RNA")
) {
  available <- preferred[preferred %in% names(object@assays)]
  if (!length(available)) {
    figure7_stop(
      "No preferred PCA assay is present: ",
      paste(preferred, collapse = ", ")
    )
  }
  available[[1L]]
}

figure7_upstream_prepare_pca_assay <- function(
  object,
  assay,
  fallback_features = 3000L
) {
  Seurat::DefaultAssay(object) <- assay
  if (assay == "RNA") {
    object <- figure7_upstream_join_layers(object, assay)
    data <- figure7_upstream_assay_data(object, assay, "data")
    if (is.null(data) || !nrow(data) || !ncol(data)) {
      object <- Seurat::NormalizeData(
        object,
        assay = assay,
        verbose = FALSE
      )
    }
    if (length(Seurat::VariableFeatures(object)) < 2L) {
      object <- Seurat::FindVariableFeatures(
        object,
        assay = assay,
        nfeatures = fallback_features,
        verbose = FALSE
      )
    }
    object <- Seurat::ScaleData(
      object,
      assay = assay,
      features = Seurat::VariableFeatures(object),
      verbose = FALSE
    )
    return(object)
  }
  if (length(Seurat::VariableFeatures(object)) < 2L) {
    scaled <- figure7_upstream_assay_data(object, assay, "scale.data")
    data <- figure7_upstream_assay_data(object, assay, "data")
    candidates <- if (!is.null(scaled) && nrow(scaled) > 1L) {
      rownames(scaled)
    } else if (!is.null(data) && nrow(data) > 1L) {
      rownames(data)
    } else {
      character()
    }
    if (length(candidates) < 2L) {
      figure7_stop("PCA assay has fewer than two usable features: ", assay)
    }
    Seurat::VariableFeatures(object) <- utils::head(
      candidates,
      min(length(candidates), fallback_features)
    )
  }
  scaled <- figure7_upstream_assay_data(object, assay, "scale.data")
  if (is.null(scaled) || !nrow(scaled) || !ncol(scaled)) {
    object <- Seurat::ScaleData(
      object,
      assay = assay,
      features = Seurat::VariableFeatures(object),
      verbose = FALSE
    )
  }
  object
}

figure7_upstream_finalize_clusters <- function(
  object,
  jobs = 1L,
  rerun_reductions = TRUE
) {
  figure7_upstream_set_single_thread()
  set.seed(1234)
  figure7_upstream_require_namespace("Seurat")
  figure7_upstream_require_columns(
    object@meta.data,
    c("manual_merge_test", "sample", "IDs"),
    "Merged Seurat metadata"
  )
  cluster_values <- as.character(object$manual_merge_test)
  if (anyNA(cluster_values) || any(!nzchar(cluster_values))) {
    figure7_stop("Merged Seurat object contains missing cluster labels")
  }
  remove <- c("9", "4", "3", "9c")
  missing <- setdiff(remove, unique(cluster_values))
  if (length(missing)) {
    figure7_stop(
      "Final filtering is missing configured cluster(s): ",
      paste(missing, collapse = ", ")
    )
  }
  source_levels <- if (is.factor(object$manual_merge_test)) {
    levels(object$manual_merge_test)
  } else {
    figure7_upstream_sort_maybe_numeric(cluster_values)
  }
  retained_levels <- setdiff(source_levels, remove)
  keep <- !cluster_values %in% remove
  cells <- rownames(object@meta.data)[keep]
  object <- subset(object, cells = cells)
  if ("clusters" %in% colnames(object@meta.data)) {
    object$clusters_orig <- object$clusters
  }
  object$clusters <- factor(
    as.character(object$manual_merge_test),
    levels = retained_levels
  )
  object$manual_merge_test <- NULL
  Seurat::Idents(object) <- object$clusters
  object@meta.data <- figure7_upstream_derive_tn_ploidy(
    object@meta.data
  )
  if (!isTRUE(rerun_reductions)) return(object)

  figure7_upstream_configure_future(jobs, max_size_gb = 60)
  assay <- figure7_upstream_choose_pca_assay(object)
  object <- figure7_upstream_prepare_pca_assay(
    object,
    assay,
    fallback_features = 3000L
  )
  Seurat::DefaultAssay(object) <- assay
  features <- Seurat::VariableFeatures(object)
  if (length(features) < 2L) {
    figure7_stop("Fewer than two variable features remain for final PCA")
  }
  number_pcs <- min(
    50L,
    length(features) - 1L,
    ncol(object) - 1L
  )
  if (number_pcs < 2L) {
    figure7_stop("Too few cells/features remain for final PCA")
  }
  object <- Seurat::RunPCA(
    object,
    assay = assay,
    features = features,
    npcs = number_pcs,
    verbose = FALSE
  )
  dimensions <- seq_len(min(
    30L,
    ncol(Seurat::Embeddings(object, reduction = "pca"))
  ))
  Seurat::RunUMAP(
    object,
    reduction = "pca",
    dims = dimensions,
    reduction.name = "umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )
}

figure7_upstream_expected_counts <- function(stage) {
  counts <- list(
    integrated = c(
      `0` = 11421L, `1` = 9589L, `2` = 3644L, `3` = 3177L,
      `4` = 3016L, `5` = 1823L, `6` = 1685L, `7` = 1660L,
      `8` = 1579L, `9` = 1281L, `10` = 1267L, `11` = 1071L,
      `12` = 757L, `13` = 704L, `14` = 210L
    ),
    refined = c(
      `0` = 11421L, `1` = 9589L, `2` = 3644L, `3` = 3177L,
      `4` = 2913L, `4c` = 103L, `5` = 1823L, `6` = 1685L,
      `7` = 1660L, `8` = 1579L, `9` = 872L, `9c` = 409L,
      `10` = 1267L, `11` = 1071L, `12` = 757L, `13` = 704L,
      `14` = 210L
    ),
    merged = c(
      `0` = 22670L, `2` = 3644L, `3` = 3177L, `4` = 2913L,
      `4c` = 103L, `5` = 1823L, `6` = 1685L, `8` = 1579L,
      `9` = 872L, `9c` = 409L, `10` = 3095L, `13` = 704L,
      `14` = 210L
    ),
    final = c(
      `0` = 22670L, `2` = 3644L, `4c` = 103L, `5` = 1823L,
      `6` = 1685L, `8` = 1579L, `10` = 3095L, `13` = 704L,
      `14` = 210L
    )
  )
  if (!stage %in% names(counts)) {
    figure7_stop("Unknown upstream acceptance stage: ", stage)
  }
  counts[[stage]]
}

figure7_upstream_assert_counts <- function(values, stage) {
  observed <- table(as.character(values))
  expected <- figure7_upstream_expected_counts(stage)
  observed <- observed[names(expected)]
  if (anyNA(observed) ||
      !identical(as.integer(observed), as.integer(expected)) ||
      !setequal(names(table(as.character(values))), names(expected))) {
    observed_text <- paste(
      paste(names(table(as.character(values))),
            as.integer(table(as.character(values))),
            sep = "="),
      collapse = ","
    )
    expected_text <- paste(
      paste(names(expected), expected, sep = "="),
      collapse = ","
    )
    figure7_stop(
      "Seurat ",
      stage,
      " cluster counts differ from the reviewed contract; expected ",
      expected_text,
      "; observed ",
      observed_text
    )
  }
  invisible(TRUE)
}

figure7_upstream_expected_tumor_samples <- function() {
  counts <- c(
    "2N-A1-0" = 413L,
    "2N-A1-LR" = 369L,
    "2N-A1-R" = 196L,
    "2N-A1-RR" = 1495L,
    "2N-A2-0" = 317L,
    "2N-A2-L" = 505L,
    "2N-A4-R" = 305L,
    "2N-A4-RL" = 1280L,
    "4N-A5-0" = 888L,
    "4N-A5-RR" = 393L,
    "A5-4N-L" = 385L,
    "A5-4N-R" = 358L,
    "A6-4N-O" = 189L,
    "A6-4N-RR" = 660L,
    "4N-A8-RL" = 1832L,
    "4N-A8-RR" = 247L
  )
  data.frame(
    sample = names(counts),
    n_cells = unname(counts),
    Ploidy = c(rep("2N", 8L), rep("4N", 8L)),
    Dose = c(
      rep(0, 4L), rep(30, 2L), rep(120, 2L),
      rep(0, 4L), rep(30, 2L), rep(120, 2L)
    ),
    stringsAsFactors = FALSE
  )
}

figure7_upstream_validate_final_metadata <- function(metadata) {
  if (!is.data.frame(metadata)) {
    figure7_stop("Final Seurat metadata must be a data frame")
  }
  figure7_upstream_require_columns(
    metadata,
    c("clusters", "TN", "Ploidy", "Dose", "sample"),
    "Final Seurat metadata"
  )
  if (nrow(metadata) != 35513L) {
    figure7_stop(
      "Final Seurat metadata cell count differs from reviewed value 35513: ",
      nrow(metadata)
    )
  }

  clusters <- as.character(metadata$clusters)
  discarded <- intersect(unique(clusters), c("3", "4", "9", "9c"))
  if (length(discarded)) {
    figure7_stop(
      "Final Seurat metadata retains discarded QC cluster(s): ",
      paste(discarded, collapse = ", ")
    )
  }
  figure7_upstream_assert_counts(clusters, "final")

  assert_exact_counts <- function(values, expected, label) {
    values <- as.character(values)
    if (anyNA(values) || any(!nzchar(values))) {
      figure7_stop("Final Seurat ", label, " contains missing values")
    }
    observed <- table(values)
    if (!setequal(names(observed), names(expected)) ||
        !identical(
          as.integer(observed[names(expected)]),
          as.integer(expected)
        )) {
      figure7_stop(
        "Final Seurat ",
        label,
        " counts differ from reviewed values"
      )
    }
    invisible(TRUE)
  }

  context <- as.character(metadata$TN)
  ploidy <- as.character(metadata$Ploidy)
  assert_exact_counts(
    context,
    c(CellLine = 25681L, Tumor = 9832L),
    "TN"
  )
  assert_exact_counts(
    ploidy,
    c(`2N` = 19716L, `4N` = 15797L),
    "Ploidy"
  )

  context_ploidy <- table(
    factor(context, levels = c("CellLine", "Tumor")),
    factor(ploidy, levels = c("2N", "4N"))
  )
  expected_context_ploidy <- matrix(
    c(14836L, 4880L, 10845L, 4952L),
    nrow = 2L,
    dimnames = list(c("CellLine", "Tumor"), c("2N", "4N"))
  )
  if (!identical(
    unname(as.integer(context_ploidy)),
    unname(as.integer(expected_context_ploidy))
  )) {
    figure7_stop(
      "Final Seurat TN-by-Ploidy counts differ from reviewed values"
    )
  }

  dose_raw <- metadata$Dose
  cell_line <- context == "CellLine"
  tumor <- context == "Tumor"
  if (any(!is.na(dose_raw[cell_line]))) {
    figure7_stop("Final Seurat CellLine Dose values must be missing")
  }
  tumor_dose <- suppressWarnings(as.numeric(as.character(dose_raw[tumor])))
  if (any(!is.finite(tumor_dose)) ||
      any(!tumor_dose %in% c(0, 30, 120))) {
    figure7_stop(
      "Final Seurat Tumor Dose values must be numeric 0, 30, or 120"
    )
  }

  treated <- tumor_dose %in% c(30, 120)
  if (sum(treated) != 5335L) {
    figure7_stop(
      "Final Seurat treated Tumor count differs from reviewed value 5335: ",
      sum(treated)
    )
  }
  treated_ploidy <- ploidy[tumor][treated]
  assert_exact_counts(
    treated_ploidy,
    c(`2N` = 2407L, `4N` = 2928L),
    "treated Tumor Ploidy"
  )

  expected_samples <- figure7_upstream_expected_tumor_samples()
  tumor_sample <- as.character(metadata$sample[tumor])
  observed_sample_counts <- table(factor(
    tumor_sample,
    levels = expected_samples$sample
  ))
  if (anyNA(tumor_sample) || any(!nzchar(tumor_sample)) ||
      !setequal(unique(tumor_sample), expected_samples$sample) ||
      !identical(
        as.integer(observed_sample_counts),
        as.integer(expected_samples$n_cells)
      )) {
    figure7_stop(
      "Final Seurat Tumor sample counts differ from the reviewed ",
      "16-sample universe"
    )
  }
  sample_match <- match(tumor_sample, expected_samples$sample)
  if (any(ploidy[tumor] != expected_samples$Ploidy[sample_match]) ||
      any(tumor_dose != expected_samples$Dose[sample_match])) {
    figure7_stop(
      "Final Seurat Tumor sample-to-Ploidy/Dose mapping differs from ",
      "reviewed values"
    )
  }
  invisible(TRUE)
}

figure7_upstream_validate_final <- function(
  object,
  strict_counts = TRUE
) {
  if (!inherits(object, "Seurat")) {
    figure7_stop("Final reconstruction output is not a Seurat object")
  }
  required_metadata <- c(
    "clusters", "cluster_cell_cycle_annotation", "S.Score", "G2M.Score",
    "Phase", "Ploidy", "TN", "sample", "Dose", "IDs", "harvest",
    "barcode_raw", "nCount_RNA", "nFeature_RNA", "percent.mt",
    "scDblFinder.class", "scDblFinder.score"
  )
  figure7_upstream_require_columns(
    object@meta.data,
    required_metadata,
    "Final reconstructed Seurat metadata"
  )
  if (!all(c("pca", "umap") %in% names(object@reductions))) {
    figure7_stop("Final reconstructed Seurat object lacks PCA/UMAP")
  }
  if (!"RNA" %in% names(object@assays)) {
    figure7_stop("Final reconstructed Seurat object lacks RNA assay")
  }
  cells <- rownames(object@meta.data)
  if (anyNA(cells) || any(!nzchar(cells)) || anyDuplicated(cells)) {
    figure7_stop("Final reconstructed Seurat object has invalid cell IDs")
  }
  annotations <- unique(as.character(
    object$cluster_cell_cycle_annotation
  ))
  expected_annotations <- c(
    "cell_cycle_candidate",
    "not_cell_cycle_candidate"
  )
  if (!setequal(annotations, expected_annotations)) {
    figure7_stop("Final Seurat cluster cell-cycle annotation is invalid")
  }
  if (isTRUE(strict_counts)) {
    figure7_upstream_validate_final_metadata(object@meta.data)
  }
  invisible(TRUE)
}
