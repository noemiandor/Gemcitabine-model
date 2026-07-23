#!/usr/bin/env Rscript

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (grepl("^--[^=]+=", token)) {
      key <- sub("^--([^=]+)=.*$", "\\1", token)
      out[[key]] <- sub("^--[^=]+=", "", token)
      i <- i + 1L
    } else if (grepl("^--", token)) {
      key <- sub("^--", "", token)
      if (i < length(args) && !grepl("^--", args[[i + 1L]])) {
        out[[key]] <- args[[i + 1L]]
        i <- i + 2L
      } else {
        out[[key]] <- "TRUE"
        i <- i + 1L
      }
    } else {
      i <- i + 1L
    }
  }
  out
}

arg_value <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (!is.null(value) && length(value) == 1L && nzchar(value)) value else default
}

arg_flag <- function(args, name, default = FALSE) {
  value <- arg_value(args, name, NULL)
  if (is.null(value)) return(default)
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript Code/in-vivo/figure7/generate_scvelo_cell_metrics.R \\",
      "    --seurat_rds \"/path/to/integrated_sct_cca_seurat_final_reclustered.rds\" \\",
      "    --loom_root \"/path/to/velocyto_loom\" \\",
      "    --output Data/in-vivo/scvelo_cell_metrics.csv \\",
      "    --seurat_metadata_output Data/in-vivo/seurat_metadata.csv \\",
      "    --seurat_metadata_provenance_output Data/in-vivo/seurat_metadata_provenance.tsv",
      "",
      "Alternatively, use the legacy --input_root layout:",
      "  seurat_obj_annotated/integrated_sct_cca_seurat_final_reclustered.rds",
      "  velocyto_loom/<sample_folder>/<sample_folder>.loom",
      "",
      "Common options:",
      "  --python /path/to/python",
      "  --work_dir /path/to/work_dir",
      "  --keep_work",
      "  --root_clusters 6",
      "  --n_jobs 16",
      sep = "\n"
    ),
    "\n"
  )
}

resolve_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_match <- grep("--file=", cmd_args, value = TRUE)
  if (length(file_match) > 0) {
    return(normalizePath(sub("--file=", "", file_match[1]), mustWork = TRUE))
  }
  NA_character_
}

is_absolute_path <- function(path) {
  grepl("^/", path)
}

create_scvelo_multiprocessing_tmpdir <- function() {
  tmp_root <- "/tmp"
  if (!dir.exists(tmp_root) || file.access(tmp_root, mode = 2L) != 0L) {
    stop("A writable /tmp directory is required for scVelo multiprocessing", call. = FALSE)
  }
  path <- tempfile("f7mp_", tmpdir = tmp_root)
  if (!dir.create(path, recursive = FALSE, showWarnings = FALSE)) {
    stop("Cannot create short scVelo multiprocessing temp directory: ", path, call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

resolve_output_path <- function(path, repo_root) {
  if (is_absolute_path(path)) {
    normalizePath(path, mustWork = FALSE)
  } else {
    normalizePath(file.path(repo_root, path), mustWork = FALSE)
  }
}

require_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing required R package: ", pkg, call. = FALSE)
  }
}

file_sha256_local <- function(path) {
  require_package("digest")
  unname(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
}

write_tsv_checked <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    data, path, sep = "\t", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = "NA"
  )
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to write TSV output: ", path, call. = FALSE)
  }
  invisible(path)
}

git_revision_local <- function(repo_root) {
  result <- suppressWarnings(system2(
    "git", c("-C", shQuote(repo_root), "rev-parse", "HEAD"),
    stdout = TRUE, stderr = FALSE
  ))
  if (length(result) == 1L && grepl("^[0-9a-f]{40}$", result[[1L]])) result[[1L]] else "not_available"
}

split_values <- function(value) {
  value <- trimws(as.character(value))
  if (!nzchar(value)) return(character(0))
  out <- unlist(strsplit(value, "[,;]", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  out[!is.na(out) & nzchar(out)]
}

sanitize_path_component_local <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == ""] <- "NA"
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x[x == ""] <- "NA"
  x
}

standardize_in_vivo_dose_local <- function(dose_values) {
  dose_chr <- trimws(as.character(dose_values))
  dose_chr[dose_chr %in% c("0", "0mg", "0 mg/kg", "0mg/kg")] <- "0mg/kg"
  dose_chr[dose_chr %in% c("30", "30mg", "30 mg/kg", "30mg/kg")] <- "30mg/kg"
  dose_chr[dose_chr %in% c("120", "120mg", "120 mg/kg", "120mg/kg")] <- "120mg/kg"
  dose_chr[dose_chr == ""] <- NA_character_
  dose_chr
}

resolve_col <- function(df, candidates) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  lower_names <- tolower(colnames(df))
  for (candidate in candidates) {
    idx <- match(tolower(candidate), lower_names)
    if (!is.na(idx)) return(colnames(df)[idx])
  }
  NA_character_
}

infer_ploidy_from_text <- function(values) {
  values <- trimws(as.character(values))
  values_upper <- toupper(values)
  out <- rep(NA_character_, length(values))
  out[!is.na(values_upper) & grepl("2N", values_upper, fixed = TRUE)] <- "2N"
  out[is.na(out) & !is.na(values_upper) & grepl("4N", values_upper, fixed = TRUE)] <- "4N"
  out
}

infer_ploidy_from_sample_label <- function(values) {
  values <- trimws(as.character(values))
  out <- infer_ploidy_from_text(values)
  out[is.na(out) & !is.na(values) & grepl("^A5", values, ignore.case = TRUE)] <- "4N"
  out[is.na(out) & !is.na(values) & grepl("^A6", values, ignore.case = TRUE)] <- "4N"
  out
}

infer_trajectory_context_from_ids <- function(values) {
  values <- trimws(as.character(values))
  values_lower <- tolower(values)
  ifelse(!is.na(values_lower) & grepl("cell-culture", values_lower, fixed = TRUE), "CellLine", "Tumor")
}

infer_sample_type <- function(values) {
  values <- trimws(as.character(values))
  out <- rep(NA_character_, length(values))
  out[grepl("^2N", values)] <- "2N-tumor"
  out[grepl("^4N", values)] <- "4N-tumor"
  out[grepl("^A5", values)] <- "4N-tumor"
  out[grepl("^A6", values)] <- "4N-tumor"
  out[values == "2N-Cell-Culture"] <- "2N-cellline"
  out[values == "4N-Cell-Culture"] <- "4N-cellline"
  out
}

make_branch_label <- function(root_clusters, end_clusters) {
  root_values <- split_values(root_clusters)
  end_values <- split_values(end_clusters)
  if (length(root_values) == 0) {
    stop("--root_clusters must contain at least one cluster label.", call. = FALSE)
  }
  root_label <- paste(sanitize_path_component_local(root_values), collapse = "_")
  end_label <- if (length(end_values) == 0) "NULL" else paste(sanitize_path_component_local(end_values), collapse = "_")
  paste0("ROOT_", root_label, "_END_", end_label)
}

build_all_cells_metadata <- function(
  seurat_rds,
  sample_folder_col = "sample_folder",
  dose_col = "Dose",
  id_col = "ID",
  cluster_col = "cluster_final",
  expected_doses = c("0mg/kg", "30mg/kg", "120mg/kg"),
  n_pcs = 30L,
  pca_prefix = "PCA_",
  umap_reduction = "umap",
  pca_reduction = "pca",
  cluster_annotation_col = "cluster_cell_cycle_annotation",
  base_cluster_col = "integrated_snn_res.0.6",
  source_qc_fields = character(),
  root_clusters = "6",
  end_clusters = "",
  seurat_metadata_output = NULL
) {
  message("Reading Seurat object: ", seurat_rds)
  obj <- readRDS(seurat_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS is not a Seurat object: ", seurat_rds, call. = FALSE)
  }
  require_package("Seurat")

  seurat_metadata_raw <- obj@meta.data
  meta <- seurat_metadata_raw
  if (!(dose_col %in% colnames(meta))) {
    stop("Seurat metadata is missing dose column: ", dose_col, call. = FALSE)
  }
  if (!(umap_reduction %in% names(obj@reductions))) {
    stop("Seurat object is missing reduction: ", umap_reduction, call. = FALSE)
  }
  if (!(pca_reduction %in% names(obj@reductions))) {
    stop("Seurat object is missing reduction: ", pca_reduction, call. = FALSE)
  }

  resolved_sample_folder_col <- resolve_col(
    meta,
    c(sample_folder_col, "sample_folder", "Sequencing.IDs", "sample", "orig.ident", "IDs")
  )
  if (is.na(resolved_sample_folder_col)) {
    stop("Cannot resolve a sample folder metadata column.", call. = FALSE)
  }
  sample_folder_col <- resolved_sample_folder_col

  if (!(cluster_col %in% colnames(meta))) {
    cluster_col <- resolve_col(meta, c("clusters", "seurat_clusters"))
  }
  if (is.na(cluster_col)) {
    stop("Cannot find cluster_final, clusters, or seurat_clusters metadata column.", call. = FALSE)
  }
  for (required_field in c(cluster_annotation_col, base_cluster_col)) {
    if (!(required_field %in% colnames(meta))) {
      stop("Seurat metadata is missing reviewed SI Figure 4 field: ", required_field, call. = FALSE)
    }
  }
  missing_qc_fields <- setdiff(source_qc_fields, colnames(meta))
  if (length(missing_qc_fields)) {
    stop("Seurat metadata is missing reviewed QC field(s): ", paste(missing_qc_fields, collapse = ", "), call. = FALSE)
  }

  dose_values <- standardize_in_vivo_dose_local(meta[[dose_col]])
  observed_doses <- unique(dose_values[!is.na(dose_values)])
  missing_doses <- setdiff(expected_doses, observed_doses)
  if (length(missing_doses) > 0) {
    stop("Missing expected Dose group(s): ", paste(missing_doses, collapse = ", "), call. = FALSE)
  }
  dose_levels <- c(expected_doses, sort(setdiff(observed_doses, expected_doses)))
  dose_levels <- dose_levels[dose_levels %in% observed_doses]
  meta[[dose_col]] <- factor(dose_values, levels = dose_levels)

  cells <- rownames(meta)
  umap <- as.data.frame(Seurat::Embeddings(obj, umap_reduction), check.names = FALSE)
  if (ncol(umap) < 2) {
    stop("UMAP embedding has fewer than 2 dimensions.", call. = FALSE)
  }
  umap <- umap[cells, seq_len(2), drop = FALSE]
  colnames(umap) <- c("UMAP_1", "UMAP_2")

  pca <- as.data.frame(Seurat::Embeddings(obj, pca_reduction), check.names = FALSE)
  if (ncol(pca) < 1) {
    stop("PCA embedding has no usable dimensions.", call. = FALSE)
  }
  pca <- pca[cells, seq_len(min(ncol(pca), n_pcs)), drop = FALSE]
  colnames(pca) <- paste0(pca_prefix, seq_len(ncol(pca)))

  id_candidates <- unique(c("IDs", id_col, "ID"))
  meta_cols <- unique(c(
    dose_col,
    id_candidates,
    cluster_col,
    "clusters",
    "sample",
    sample_folder_col,
    "Sequencing.IDs",
    "orig.ident",
    "sample_type",
    cluster_annotation_col,
    base_cluster_col,
    "cluster_final_annotation_primary",
    "cluster_final_annotation_multi",
    "seurat_clusters",
    "nCount_RNA",
    "nFeature_RNA",
    "percent.mt"
  ))
  meta_cols <- meta_cols[meta_cols %in% colnames(meta)]

  meta_all <- data.frame(cell = cells, meta[, meta_cols, drop = FALSE], check.names = FALSE)
  meta_all <- cbind(meta_all, umap, pca)
  meta_all$Dose <- factor(as.character(meta[[dose_col]]), levels = dose_levels)
  meta_all$dose_safe <- sanitize_path_component_local(as.character(meta[[dose_col]]))

  if (!("sample_type" %in% colnames(meta_all)) && "sample" %in% colnames(meta_all)) {
    meta_all$sample_type <- infer_sample_type(meta_all$sample)
  }

  resolved_id_col <- resolve_col(meta_all, id_candidates)
  if (is.na(resolved_id_col)) {
    stop("Cannot resolve an ID metadata column. Tried: ", paste(id_candidates, collapse = ", "), call. = FALSE)
  }

  n_cells <- nrow(meta_all)
  ploidy_from_id <- infer_ploidy_from_text(meta_all[[resolved_id_col]])
  ploidy_from_sample_folder <- if (sample_folder_col %in% colnames(meta_all)) {
    infer_ploidy_from_sample_label(meta_all[[sample_folder_col]])
  } else {
    rep(NA_character_, n_cells)
  }
  ploidy_from_sample <- if ("sample" %in% colnames(meta_all)) {
    infer_ploidy_from_sample_label(meta_all$sample)
  } else {
    rep(NA_character_, n_cells)
  }

  meta_all$ploidy <- ploidy_from_id
  meta_all$ploidy_source <- ifelse(!is.na(ploidy_from_id), resolved_id_col, NA_character_)
  missing_ploidy <- is.na(meta_all$ploidy)
  meta_all$ploidy[missing_ploidy] <- ploidy_from_sample_folder[missing_ploidy]
  meta_all$ploidy_source[missing_ploidy & !is.na(ploidy_from_sample_folder)] <- sample_folder_col
  missing_ploidy <- is.na(meta_all$ploidy)
  meta_all$ploidy[missing_ploidy] <- ploidy_from_sample[missing_ploidy]
  meta_all$ploidy_source[missing_ploidy & !is.na(ploidy_from_sample)] <- "sample"

  if (any(is.na(meta_all$ploidy))) {
    bad_cells <- unique(meta_all$cell[is.na(meta_all$ploidy)])
    stop("Could not infer ploidy. Example unresolved cells: ", paste(head(bad_cells, 20), collapse = ", "), call. = FALSE)
  }

  ploidy_levels <- c("2N", "4N")
  ploidy_levels <- ploidy_levels[ploidy_levels %in% unique(meta_all$ploidy)]
  meta_all$ploidy <- factor(meta_all$ploidy, levels = ploidy_levels)
  meta_all$trajectory_context <- factor(
    infer_trajectory_context_from_ids(meta_all[[resolved_id_col]]),
    levels = c("CellLine", "Tumor")
  )
  meta_all$Ploidy <- factor(as.character(meta_all$ploidy), levels = c("2N", "4N"))
  meta_all$TN <- factor(as.character(meta_all$trajectory_context), levels = c("CellLine", "Tumor"))
  meta_all$Dose_DEG <- factor(as.character(meta_all$Dose), levels = expected_doses)
  if (!("clusters" %in% colnames(meta_all))) {
    meta_all$clusters <- as.character(meta_all[[cluster_col]])
  }
  meta_all$trajectory_group <- factor(as.character(meta_all$TN), levels = c("CellLine", "Tumor"))
  meta_all$trajectory_shape_group <- factor(as.character(meta_all$TN), levels = c("CellLine", "Tumor"))
  meta_all$trajectory_analysis_group <- "All_cells"
  meta_all$trajectory_tn_scope <- "All"
  meta_all$trajectory_ploidy_scope <- "ploidy_all"
  meta_all$trajectory_branch <- make_branch_label(root_clusters, end_clusters)

  meta_all$cell <- enc2utf8(as.character(meta_all$cell))
  if (anyNA(meta_all$cell) || any(!nzchar(meta_all$cell))) {
    stop("Metadata contains missing or empty cell IDs.", call. = FALSE)
  }
  if (anyDuplicated(meta_all$cell)) {
    stop("Metadata contains duplicated cell IDs.", call. = FALSE)
  }
  meta_all <- meta_all[
    order(meta_all$cell, method = "radix"),
    ,
    drop = FALSE
  ]
  rownames(meta_all) <- NULL

  if (!is.null(seurat_metadata_output) && nzchar(seurat_metadata_output)) {
    metadata_index <- match(meta_all$cell, rownames(meta))
    if (anyNA(metadata_index)) {
      stop("Could not align Seurat metadata rows to exported cell IDs.", call. = FALSE)
    }
    seurat_metadata <- data.frame(
      cell = meta_all$cell,
      seurat_metadata_raw[metadata_index, , drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    exported_umap <- umap[meta_all$cell, c("UMAP_1", "UMAP_2"), drop = FALSE]
    seurat_metadata$UMAP_1 <- as.numeric(exported_umap$UMAP_1)
    seurat_metadata$UMAP_2 <- as.numeric(exported_umap$UMAP_2)
    metadata_output_dir <- dirname(seurat_metadata_output)
    if (!dir.exists(metadata_output_dir)) {
      dir.create(metadata_output_dir, recursive = TRUE, showWarnings = FALSE)
    }
    readr::write_csv(seurat_metadata, seurat_metadata_output, na = "NA")
    if (!file.exists(seurat_metadata_output) || file.info(seurat_metadata_output)$size <= 0) {
      stop("Failed to write Seurat metadata output: ", seurat_metadata_output, call. = FALSE)
    }
    message("Wrote Seurat metadata: ", seurat_metadata_output)
  }

  attr(meta_all, "sample_folder_col") <- sample_folder_col
  attr(meta_all, "cluster_col") <- cluster_col
  attr(meta_all, "source_metadata_columns") <- colnames(seurat_metadata_raw)
  meta_all
}

resolve_loom_files <- function(loom_root, metadata_df, sample_folder_col) {
  if (!dir.exists(loom_root)) {
    stop("Missing loom root: ", loom_root, call. = FALSE)
  }
  if (!(sample_folder_col %in% colnames(metadata_df))) {
    stop("Metadata does not contain sample folder column: ", sample_folder_col, call. = FALSE)
  }
  sample_folders <- sort(unique(as.character(metadata_df[[sample_folder_col]])))
  sample_folders <- sample_folders[!is.na(sample_folders) & nzchar(sample_folders)]
  nested_files <- file.path(loom_root, sample_folders, paste0(sample_folders, ".loom"))
  flat_files <- file.path(loom_root, paste0(sample_folders, ".loom"))
  loom_files <- ifelse(file.exists(nested_files), nested_files, flat_files)
  missing <- loom_files[!file.exists(loom_files)]
  if (length(missing) > 0) {
    missing_samples <- sample_folders[!file.exists(loom_files)]
    expected <- vapply(missing_samples, function(sample) {
      paste0(
        file.path(loom_root, paste0(sample, ".loom")), " or ",
        file.path(loom_root, sample, paste0(sample, ".loom"))
      )
    }, character(1L))
    stop("Missing expected loom file(s):\n", paste(expected, collapse = "\n"), call. = FALSE)
  }
  normalizePath(loom_files, mustWork = TRUE)
}

embedded_scvelo_worker <- function() {
  r"---(
#!/usr/bin/env python3

import argparse
import inspect
import re
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Standalone scVelo metrics worker embedded by generate_scvelo_cell_metrics.R."
    )
    parser.add_argument("--input", required=True, help="Comma-separated .loom/.h5ad velocity inputs.")
    parser.add_argument("--metadata", required=True, help="Seurat-exported cell metadata CSV.")
    parser.add_argument("--output", required=True, help="Output scvelo_cell_metrics.csv path.")
    parser.add_argument("--cluster-col", default="clusters")
    parser.add_argument("--mode", default="stochastic", choices=["deterministic", "stochastic", "dynamical"])
    parser.add_argument("--min-shared-counts", type=int, default=20)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--min-matched-cells", type=int, default=50)
    parser.add_argument("--no-use-r-umap", action="store_true")
    parser.add_argument("--use-metadata-pca", action="store_true")
    parser.add_argument("--pca-prefix", default="PCA_")
    parser.add_argument("--root-clusters", default="")
    parser.add_argument("--end-clusters", default="")
    parser.add_argument("--root-key-name", default="root_cells_manual")
    parser.add_argument("--end-key-name", default="end_points_manual")
    return parser.parse_args()


def import_required():
    try:
        import numpy as np
        import pandas as pd
        import anndata as ad
        import scanpy as sc
        import scvelo as scv
    except ImportError as exc:
        raise SystemExit(
            "Missing Python package for scVelo metrics extraction: "
            f"{exc}. Use --python to point to an environment with scvelo, scanpy, "
            "anndata, pandas, numpy."
        ) from exc
    return np, pd, ad, sc, scv


def clean_cell_id(value):
    value = str(value).strip()
    value = re.sub(r":", "_", value)
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"x$", "", value)
    return value


def extract_10x_barcode(value):
    value = str(value).strip()
    if not value:
        return ""
    parts = re.split(r"[:_]", value)
    for part in reversed(parts):
        part = part.strip()
        match = re.search(r"([ACGTN]{8,})(?:[-_]\d+)?x?$", part, flags=re.IGNORECASE)
        if match:
            return match.group(1).upper()
    match = re.search(r"([ACGTN]{8,})(?:[-_]\d+)?x?", value, flags=re.IGNORECASE)
    return match.group(1).upper() if match else ""


def normalize_sample_value(value):
    value = str(value).strip()
    if not value or value in {"nan", "None", "NA"}:
        return ""
    return value


def infer_sample_from_cell(value, known_samples):
    value = str(value).strip()
    if not value:
        return ""
    for sample in known_samples:
        sample = str(sample)
        if value == sample or value.startswith(f"{sample}:") or value.startswith(f"{sample}_"):
            return sample
    return ""


def make_match_keys(value, sample="", known_samples=None):
    known_samples = known_samples or []
    value = str(value).strip()
    sample = normalize_sample_value(sample)
    barcode = extract_10x_barcode(value)
    inferred_sample = sample or infer_sample_from_cell(value, known_samples)
    keys = []
    for candidate in [value, clean_cell_id(value)]:
        if candidate:
            keys.append(("cell", f"cell:{candidate}"))
            keys.append(("cell_upper", f"cell_upper:{candidate.upper()}"))
    if barcode:
        keys.append(("barcode", f"barcode:{barcode}"))
        if inferred_sample:
            keys.append(("sample_barcode", f"sample_barcode:{inferred_sample}|{barcode}"))
    return keys


def build_unique_key_to_obs(obs_names, known_samples):
    key_to_obs = {}
    duplicated = set()
    for obs in obs_names:
        for _, key in make_match_keys(obs, known_samples=known_samples):
            if key in key_to_obs and key_to_obs[key] != obs:
                duplicated.add(key)
            elif key not in duplicated:
                key_to_obs[key] = obs
    for key in duplicated:
        key_to_obs.pop(key, None)
    return key_to_obs


def parse_input_paths(input_arg):
    paths = [Path(x.strip()) for x in str(input_arg).split(",") if x.strip()]
    if not paths:
        raise SystemExit("--input did not contain usable file paths.")
    return paths


def read_velocity_input(path, scv, sc):
    suffix = path.suffix.lower()
    if suffix == ".loom":
        if hasattr(scv, "read"):
            return scv.read(str(path), cache=False)
        if hasattr(scv, "read_loom"):
            return scv.read_loom(str(path), cache=False)
        return sc.read_loom(str(path))
    if suffix == ".h5ad":
        return sc.read_h5ad(str(path))
    raise SystemExit(f"Unsupported velocity input format: {path}")


def read_velocity_inputs(input_arg, scv, sc, ad):
    paths = parse_input_paths(input_arg)
    adatas = []
    keys = []
    for idx, path in enumerate(paths):
        adata_i = read_velocity_input(path, scv, sc)
        adata_i.var_names_make_unique()
        adata_i.obs_names_make_unique()
        adata_i.obs["velocity_source_file"] = str(path)
        adata_i.obs["velocity_source_index"] = idx
        key = path.stem
        keys.append(key)
        adatas.append(adata_i)
    if len(adatas) == 1:
        return adatas[0], [str(paths[0])]
    merged = ad.concat(
        adatas,
        axis=0,
        join="outer",
        merge="same",
        label="velocity_batch",
        keys=keys,
        index_unique=None,
    )
    if not merged.obs_names.is_unique:
        merged.obs["velocity_original_cell"] = merged.obs_names.astype(str)
        merged.obs_names = [
            f"{batch}:{cell}"
            for batch, cell in zip(
                merged.obs["velocity_batch"].astype(str),
                merged.obs["velocity_original_cell"].astype(str),
            )
        ]
    return merged, [str(path) for path in paths]


def choose_metadata_sample_column(meta):
    for col in ["sample_folder", "Sequencing.IDs", "sample", "orig.ident"]:
        if col in meta.columns:
            return col
    return ""


def load_metadata(path, pd):
    meta = pd.read_csv(path)
    if "cell" not in meta.columns:
        raise SystemExit(f"Metadata file is missing required column 'cell': {path}")
    meta["cell"] = meta["cell"].astype(str)
    return meta


def match_cells(adata, meta, pd):
    sample_col = choose_metadata_sample_column(meta)
    known_samples = []
    if sample_col:
        known_samples = sorted(meta[sample_col].dropna().astype(str).unique().tolist(), key=len, reverse=True)
    key_to_obs = build_unique_key_to_obs([str(x) for x in adata.obs_names], known_samples)
    rows = []
    for _, row in meta.iterrows():
        cell = str(row["cell"])
        sample = str(row[sample_col]) if sample_col else ""
        velocity_hint = str(row["velocity_cell"]) if "velocity_cell" in meta.columns and pd.notna(row["velocity_cell"]) else cell
        matched_obs = None
        match_method = "none"
        for method, key in make_match_keys(velocity_hint, sample=sample, known_samples=known_samples):
            if key in key_to_obs:
                matched_obs = key_to_obs[key]
                match_method = method
                break
        if matched_obs is None and velocity_hint != cell:
            for method, key in make_match_keys(cell, sample=sample, known_samples=known_samples):
                if key in key_to_obs:
                    matched_obs = key_to_obs[key]
                    match_method = method
                    break
        if matched_obs is None:
            continue
        out = row.to_dict()
        out["velocity_obs"] = matched_obs
        out["velocity_cell"] = matched_obs
        out["match_method"] = match_method
        rows.append(out)
    return pd.DataFrame(rows)


def ensure_velocity_layers(adata):
    missing = [layer for layer in ["spliced", "unspliced"] if layer not in adata.layers]
    if missing:
        raise SystemExit(f"Velocity input is missing required layer(s): {', '.join(missing)}")


def set_metadata_and_umap(adata, matched_meta, args, np):
    adata = adata[matched_meta["velocity_obs"].astype(str).tolist()].copy()
    adata.obs_names = matched_meta["cell"].astype(str).values
    skip_cols = {"velocity_obs", "match_method"}
    for col in matched_meta.columns:
        if col in skip_cols:
            continue
        adata.obs[col] = matched_meta[col].values
    if not args.no_use_r_umap and {"UMAP_1", "UMAP_2"}.issubset(adata.obs.columns):
        coords = adata.obs[["UMAP_1", "UMAP_2"]].astype(float).to_numpy()
        if np.isfinite(coords).all():
            adata.obsm["X_umap"] = coords.astype("float32", copy=False)
    return adata


def scvelo_filter_and_normalize(adata, args, scv):
    signature = inspect.signature(scv.pp.filter_and_normalize)
    kwargs = {"min_shared_counts": args.min_shared_counts}
    if "n_top_genes" in signature.parameters:
        kwargs["n_top_genes"] = args.n_top_genes
    try:
        scv.pp.filter_and_normalize(adata, **kwargs)
    except TypeError:
        scv.pp.filter_and_normalize(adata, min_shared_counts=args.min_shared_counts)


def get_pca_columns(columns, prefix):
    pca_cols = []
    for col in columns:
        if str(col).startswith(prefix):
            suffix = str(col)[len(prefix):]
            if suffix.isdigit():
                pca_cols.append((int(suffix), col))
    pca_cols = [col for _, col in sorted(pca_cols)]
    return pca_cols


def apply_metadata_pca_neighbors(adata, args, sc, np):
    if not args.use_metadata_pca:
        return False
    pca_cols = get_pca_columns(adata.obs.columns, args.pca_prefix)
    if len(pca_cols) == 0:
        raise SystemExit(f"No metadata PCA columns found with prefix {args.pca_prefix!r}")
    pca_cols = pca_cols[: args.n_pcs]
    pca = adata.obs[pca_cols].astype(float).to_numpy()
    if not np.isfinite(pca).all():
        raise SystemExit("Metadata PCA contains non-finite values.")
    adata.obsm["X_pca"] = pca.astype("float32", copy=False)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=min(args.n_pcs, pca.shape[1]), use_rep="X_pca")
    return True


def run_scvelo_moments(adata, args, scv):
    if args.use_metadata_pca:
        try:
            scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
        except TypeError:
            scv.pp.moments(adata)
        return
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)


def parse_comma_arg(value):
    return [x.strip() for x in str(value).split(",") if x.strip()]


def set_manual_terminal_states(adata, args):
    root_clusters = parse_comma_arg(args.root_clusters)
    end_clusters = parse_comma_arg(args.end_clusters)
    clusters = adata.obs[args.cluster_col].astype(str) if args.cluster_col in adata.obs.columns else None
    if clusters is None and (root_clusters or end_clusters):
        raise SystemExit(f"Manual root/end clusters requested, but cluster column {args.cluster_col!r} is absent.")
    if root_clusters:
        root_mask = clusters.isin(root_clusters)
        if int(root_mask.sum()) == 0:
            observed = sorted(clusters.dropna().unique().tolist())
            raise SystemExit(f"No root cells found in {root_clusters}. Observed cluster labels include: {', '.join(observed[:30])}")
        adata.obs[args.root_key_name] = root_mask.astype(float).values
    if end_clusters:
        end_mask = clusters.isin(end_clusters)
        if int(end_mask.sum()) == 0:
            observed = sorted(clusters.dropna().unique().tolist())
            raise SystemExit(f"No endpoint cells found in {end_clusters}. Observed cluster labels include: {', '.join(observed[:30])}")
        adata.obs[args.end_key_name] = end_mask.astype(float).values


def run_velocity_pseudotime(adata, args, scv):
    kwargs = {}
    signature = inspect.signature(scv.tl.velocity_pseudotime)
    if args.root_clusters and "root_key" in signature.parameters:
        kwargs["root_key"] = args.root_key_name
    if args.end_clusters and "end_key" in signature.parameters:
        kwargs["end_key"] = args.end_key_name
    if args.root_clusters and "root_key" not in signature.parameters:
        adata.obs["root_cells"] = adata.obs[args.root_key_name].values
    if args.end_clusters and "end_key" not in signature.parameters:
        adata.obs["end_points"] = adata.obs[args.end_key_name].values
    scv.tl.velocity_pseudotime(adata, **kwargs)


def write_metrics(adata, args, output, pd):
    metric_cols = [
        "velocity_cell",
        "Dose",
        "ploidy",
        "Ploidy",
        "TN",
        "Dose_DEG",
        "trajectory_context",
        "trajectory_group",
        "trajectory_analysis_group",
        "trajectory_shape_group",
        "trajectory_tn_scope",
        "trajectory_ploidy_scope",
        "trajectory_branch",
        args.cluster_col,
        "clusters",
        "sample",
        "sample_type",
        "cluster_final_annotation_primary",
        "velocity_source_file",
        "velocity_batch",
        "velocity_length",
        "velocity_confidence",
        "velocity_confidence_transition",
        args.root_key_name,
        args.end_key_name,
        "velocity_pseudotime",
        "latent_time",
    ]
    metric_cols = [x for x in dict.fromkeys(metric_cols) if x in adata.obs.columns]
    metrics = adata.obs[metric_cols].copy()
    metrics.insert(0, "cell", adata.obs_names.astype(str))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(output, index=False)


def main():
    args = parse_args()
    np, pd, ad, sc, scv = import_required()
    scv.settings.verbosity = 3
    adata, input_paths = read_velocity_inputs(args.input, scv, sc, ad)
    adata.var_names_make_unique()
    ensure_velocity_layers(adata)
    meta = load_metadata(args.metadata, pd)
    matched_meta = match_cells(adata, meta, pd)
    if matched_meta.shape[0] < args.min_matched_cells:
        raise SystemExit(
            f"Only {matched_meta.shape[0]} cells matched between Seurat metadata and velocity input; "
            f"minimum required is {args.min_matched_cells}."
        )
    adata = set_metadata_and_umap(adata, matched_meta, args, np)
    scvelo_filter_and_normalize(adata, args, scv)
    apply_metadata_pca_neighbors(adata, args, sc, np)
    run_scvelo_moments(adata, args, scv)
    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode=args.mode)
    scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
    scv.tl.velocity_confidence(adata)
    set_manual_terminal_states(adata, args)
    try:
        run_velocity_pseudotime(adata, args, scv)
    except Exception as exc:
        if args.root_clusters or args.end_clusters:
            raise
        print(f"WARNING: velocity_pseudotime failed: {exc}", file=sys.stderr)
    try:
        scv.tl.latent_time(adata)
    except Exception as exc:
        print(f"WARNING: latent_time failed: {exc}", file=sys.stderr)
    write_metrics(adata, args, args.output, pd)


if __name__ == "__main__":
    main()
)---"
}

write_embedded_scvelo_worker <- function(path) {
  writeLines(embedded_scvelo_worker(), con = path)
  Sys.chmod(path, mode = "0755")
  invisible(path)
}

main <- function() {
  args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
  if (arg_flag(args, "help", FALSE) || arg_flag(args, "h", FALSE)) {
    usage()
    quit(save = "no", status = 0)
  }

  script_path <- resolve_script_path()
  script_dir <- dirname(script_path)
  in_vivo_dir <- normalizePath(dirname(script_dir), mustWork = TRUE)
  repo_root <- normalizePath(dirname(dirname(in_vivo_dir)), mustWork = TRUE)
  config_path <- resolve_output_path(
    arg_value(args, "config", file.path(script_dir, "figure7_config.yaml")),
    repo_root
  )
  if (!file.exists(config_path)) stop("Missing Figure 7 config: ", config_path, call. = FALSE)
  require_package("yaml")
  config <- yaml::read_yaml(config_path)
  si_config <- config$si_figure4
  if (is.null(si_config)) stop("Figure 7 config is missing si_figure4", call. = FALSE)

  input_root_arg <- arg_value(args, "input_root", NULL)
  seurat_rds_arg <- arg_value(args, "seurat_rds", NULL)
  loom_root_arg <- arg_value(args, "loom_root", NULL)
  if (is.null(input_root_arg) && (is.null(seurat_rds_arg) || is.null(loom_root_arg))) {
    usage()
    stop("Pass both --seurat_rds and --loom_root, or use legacy --input_root.", call. = FALSE)
  }
  input_root <- if (is.null(input_root_arg)) NULL else normalizePath(input_root_arg, mustWork = TRUE)

  output_path <- resolve_output_path(arg_value(args, "output", "Data/in-vivo/scvelo_cell_metrics.csv"), repo_root)
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  seurat_metadata_output <- resolve_output_path(
    arg_value(args, "seurat_metadata_output", "Data/in-vivo/seurat_metadata.csv"),
    repo_root
  )
  seurat_metadata_provenance_output <- resolve_output_path(
    arg_value(args, "seurat_metadata_provenance_output", "Data/in-vivo/seurat_metadata_provenance.tsv"),
    repo_root
  )

  work_dir_arg <- arg_value(args, "work_dir", NULL)
  work_dir <- if (is.null(work_dir_arg)) {
    tempfile("scvelo_cell_metrics_")
  } else if (is_absolute_path(work_dir_arg)) {
    work_dir_arg
  } else {
    file.path(repo_root, work_dir_arg)
  }
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  keep_work <- arg_flag(args, "keep_work", FALSE) || !is.null(work_dir_arg)

  seurat_rds <- if (!is.null(seurat_rds_arg)) {
    resolve_output_path(seurat_rds_arg, repo_root)
  } else {
    file.path(
      input_root,
      arg_value(args, "seurat_rds_rel", file.path("seurat_obj_annotated", "integrated_sct_cca_seurat_final_reclustered.rds"))
    )
  }
  if (!file.exists(seurat_rds)) {
    stop("Missing Seurat RDS: ", seurat_rds, call. = FALSE)
  }

  require_package("Seurat")
  require_package("readr")

  n_pcs <- as.integer(arg_value(args, "n_pcs", "30"))
  n_neighbors <- as.integer(arg_value(args, "n_neighbors", "30"))
  n_jobs <- as.integer(arg_value(args, "n_jobs", "16"))
  min_shared_counts <- as.integer(arg_value(args, "min_shared_counts", "20"))
  n_top_genes <- as.integer(arg_value(args, "n_top_genes", "2000"))
  min_matched_cells <- as.integer(arg_value(args, "min_matched_cells", "50"))
  root_clusters <- arg_value(args, "root_clusters", "6")
  end_clusters <- arg_value(args, "end_clusters", "")
  mode <- arg_value(args, "mode", "stochastic")
  pca_prefix <- arg_value(args, "pca_prefix", "PCA_")
  expected_doses <- split_values(arg_value(args, "expected_doses", "0mg/kg,30mg/kg,120mg/kg"))

  metadata_df <- build_all_cells_metadata(
    seurat_rds = seurat_rds,
    sample_folder_col = arg_value(args, "sample_folder_col", "sample_folder"),
    dose_col = arg_value(args, "dose_col", "Dose"),
    id_col = arg_value(args, "id_col", "ID"),
    cluster_col = arg_value(args, "cluster_col", as.character(si_config$cluster_id_field)),
    expected_doses = expected_doses,
    n_pcs = n_pcs,
    pca_prefix = pca_prefix,
    umap_reduction = arg_value(args, "umap_reduction", as.character(si_config$umap_reduction)),
    pca_reduction = arg_value(args, "pca_reduction", as.character(si_config$pca_reduction)),
    cluster_annotation_col = as.character(si_config$cluster_annotation_field),
    base_cluster_col = as.character(si_config$base_cluster_field),
    source_qc_fields = as.character(unlist(si_config$source_qc_fields)),
    root_clusters = root_clusters,
    end_clusters = end_clusters,
    seurat_metadata_output = seurat_metadata_output
  )
  sample_folder_col <- attr(metadata_df, "sample_folder_col")
  cluster_col <- attr(metadata_df, "cluster_col")
  loom_root <- if (!is.null(loom_root_arg)) {
    resolve_output_path(loom_root_arg, repo_root)
  } else {
    file.path(input_root, "velocyto_loom")
  }
  loom_files <- resolve_loom_files(loom_root, metadata_df, sample_folder_col)

  metadata_file <- file.path(work_dir, "cells_metadata_umap.csv")
  readr::write_csv(as.data.frame(metadata_df, stringsAsFactors = FALSE), metadata_file)

  python_bin <- arg_value(args, "python", Sys.which("python"))
  if (!nzchar(python_bin)) {
    stop("Cannot find python. Pass --python /path/to/python.", call. = FALSE)
  }
  python_script <- file.path(work_dir, "standalone_scvelo_metrics_worker.py")
  write_embedded_scvelo_worker(python_script)

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )

  log_file <- file.path(work_dir, "scvelo.log")
  multiprocessing_tmpdir <- create_scvelo_multiprocessing_tmpdir()
  on.exit(unlink(multiprocessing_tmpdir, recursive = TRUE, force = TRUE), add = TRUE)

  scvelo_args <- c(
    python_script,
    "--input", paste(loom_files, collapse = ","),
    "--metadata", metadata_file,
    "--output", output_path,
    "--cluster-col", cluster_col,
    "--mode", mode,
    "--min-shared-counts", as.character(min_shared_counts),
    "--n-top-genes", as.character(n_top_genes),
    "--n-pcs", as.character(n_pcs),
    "--n-neighbors", as.character(n_neighbors),
    "--n-jobs", as.character(n_jobs),
    "--min-matched-cells", as.character(min_matched_cells),
    "--use-metadata-pca",
    "--pca-prefix", pca_prefix,
    "--root-clusters", root_clusters
  )
  if (nzchar(end_clusters)) {
    scvelo_args <- c(scvelo_args, "--end-clusters", end_clusters)
  }

  message(
    "Running scVelo metrics extraction with ", length(loom_files),
    " loom file(s); multiprocessing tmpdir: ", multiprocessing_tmpdir
  )
  status <- system2(
    python_bin,
    args = scvelo_args,
    stdout = log_file,
    stderr = log_file,
    env = c(
      paste0("TMPDIR=", multiprocessing_tmpdir),
      paste0("TMP=", multiprocessing_tmpdir),
      paste0("TEMP=", multiprocessing_tmpdir)
    )
  )
  if (!identical(status, 0L)) {
    log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character(0)
    if (length(log_lines) > 0) {
      message("scVelo log tail:")
      message(paste(tail(log_lines, 80), collapse = "\n"))
    }
    stop("scVelo metrics extraction failed. Log: ", log_file, call. = FALSE)
  }

  if (!file.exists(output_path)) {
    stop("scVelo completed but did not write output: ", output_path, call. = FALSE)
  }
  out_header <- readLines(output_path, n = 1L, warn = FALSE)
  if (!grepl("^cell,velocity_cell,", out_header)) {
    stop("Unexpected scVelo metrics header in output: ", output_path, call. = FALSE)
  }

  cellcycle_mapping <- unlist(si_config$cellcycle_mapping, use.names = TRUE)
  inclusion <- si_config$inclusion
  provenance <- data.frame(
    key = c(
      "schema_version", "artifact", "source_seurat_rds", "source_seurat_rds_sha256",
      "source_object_class", "source_object_cells", "source_metadata_columns",
      "seurat_metadata_output", "seurat_metadata_sha256", "scvelo_metrics_output",
      "scvelo_metrics_sha256", "umap_reduction", "umap_dimensions", "pca_reduction",
      "cluster_id_field", "base_cluster_field", "clustering_resolution",
      "cluster_annotation_field", "sample_field", "dose_field", "ploidy_field",
      "context_field", "cellcycle_mapping", "inclusion_context", "inclusion_ploidy_levels",
      "inclusion_dose_levels", "inclusion_required_nonmissing", "source_qc_fields",
      "source_qc_policy", "figure7_config", "figure7_config_sha256",
      "export_script", "export_script_sha256", "source_code_revision"
    ),
    value = c(
      "1", "seurat_metadata_export", normalizePath(seurat_rds, mustWork = TRUE),
      file_sha256_local(seurat_rds), "Seurat", as.character(nrow(metadata_df)),
      paste(attr(metadata_df, "source_metadata_columns"), collapse = ","),
      normalizePath(seurat_metadata_output, mustWork = TRUE), file_sha256_local(seurat_metadata_output),
      normalizePath(output_path, mustWork = TRUE), file_sha256_local(output_path),
      arg_value(args, "umap_reduction", as.character(si_config$umap_reduction)), "UMAP_1,UMAP_2",
      arg_value(args, "pca_reduction", as.character(si_config$pca_reduction)),
      attr(metadata_df, "cluster_col"), as.character(si_config$base_cluster_field),
      as.character(si_config$clustering_resolution), as.character(si_config$cluster_annotation_field),
      as.character(si_config$sample_field), as.character(si_config$dose_field),
      as.character(si_config$ploidy_field), as.character(si_config$context_field),
      paste(paste(names(cellcycle_mapping), cellcycle_mapping, sep = "->"), collapse = ";"),
      as.character(inclusion$context), paste(unlist(inclusion$ploidy_levels), collapse = ","),
      paste(unlist(inclusion$dose_levels), collapse = ","),
      paste(unlist(inclusion$required_nonmissing), collapse = ","),
      paste(unlist(si_config$source_qc_fields), collapse = ","),
      as.character(si_config$source_qc_policy), normalizePath(config_path, mustWork = TRUE),
      file_sha256_local(config_path), normalizePath(script_path, mustWork = TRUE),
      file_sha256_local(script_path), git_revision_local(repo_root)
    ),
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(provenance$key) || anyNA(provenance$value) || any(!nzchar(provenance$value))) {
    stop("Seurat metadata provenance contains missing or duplicated fields", call. = FALSE)
  }
  write_tsv_checked(provenance, seurat_metadata_provenance_output)

  message("Wrote scVelo cell metrics: ", output_path)
  message("Seurat metadata output: ", seurat_metadata_output)
  message("Seurat metadata provenance: ", seurat_metadata_provenance_output)
  if (!keep_work) {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  } else {
    message("Kept work directory: ", normalizePath(work_dir, mustWork = FALSE))
  }
}

main()
