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
      "    --output Results/in-vivo/figure7/intermediates/scvelo_cell_metrics.csv",
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
      "  --velocity_embedding_output /path/to/scvelo_velocity_umap.tsv",
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

split_values <- function(value) {
  value <- trimws(as.character(value))
  if (!nzchar(value)) return(character(0))
  out <- unlist(strsplit(value, "[,;]", perl = TRUE), use.names = FALSE)
  out <- trimws(out)
  out[!is.na(out) & nzchar(out)]
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

build_all_cells_metadata <- function(
  seurat_rds,
  sample_folder_col = "sample_folder",
  dose_col = "Dose",
  ploidy_col = "Ploidy",
  context_col = "TN",
  id_col = "ID",
  cluster_col = "cluster_final",
  expected_doses = c("0mg/kg", "30mg/kg", "120mg/kg"),
  n_pcs = 30L,
  pca_prefix = "PCA_",
  umap_reduction = "umap",
  pca_reduction = "pca"
) {
  message("Reading Seurat object: ", seurat_rds)
  obj <- readRDS(seurat_rds)
  if (!inherits(obj, "Seurat")) {
    stop("Input RDS is not a Seurat object: ", seurat_rds, call. = FALSE)
  }
  require_package("Seurat")

  meta <- obj@meta.data
  if (!(dose_col %in% colnames(meta))) {
    stop("Seurat metadata is missing dose column: ", dose_col, call. = FALSE)
  }
  for (required_col in c(ploidy_col, context_col)) {
    if (!(required_col %in% colnames(meta))) {
      stop(
        "Seurat metadata is missing reviewed field: ",
        required_col,
        call. = FALSE
      )
    }
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
    ploidy_col,
    context_col,
    id_candidates,
    cluster_col,
    "clusters",
    "sample",
    sample_folder_col,
    "Sequencing.IDs",
    "orig.ident",
    "seurat_clusters"
  ))
  meta_cols <- meta_cols[meta_cols %in% colnames(meta)]

  meta_all <- data.frame(cell = cells, meta[, meta_cols, drop = FALSE], check.names = FALSE)
  meta_all <- cbind(meta_all, umap, pca)
  meta_all$Dose <- factor(as.character(meta[[dose_col]]), levels = dose_levels)

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

  inferred_ploidy <- ploidy_from_id
  missing_ploidy <- is.na(inferred_ploidy)
  inferred_ploidy[missing_ploidy] <-
    ploidy_from_sample_folder[missing_ploidy]
  missing_ploidy <- is.na(inferred_ploidy)
  inferred_ploidy[missing_ploidy] <- ploidy_from_sample[missing_ploidy]

  reviewed_ploidy <- toupper(trimws(as.character(meta_all[[ploidy_col]])))
  reviewed_ploidy[!reviewed_ploidy %in% c("2N", "4N")] <- NA_character_
  reviewed_context <- tolower(trimws(as.character(
    meta_all[[context_col]]
  )))
  reviewed_context <- ifelse(
    reviewed_context %in% c("cellline", "cell line"),
    "CellLine",
    ifelse(reviewed_context == "tumor", "Tumor", NA_character_)
  )
  if (anyNA(reviewed_ploidy) || anyNA(reviewed_context)) {
    stop(
      "Reviewed Seurat Ploidy/TN fields contain missing or invalid values",
      call. = FALSE
    )
  }
  ploidy_conflict <- !is.na(inferred_ploidy) &
    inferred_ploidy != reviewed_ploidy
  inferred_context <- infer_trajectory_context_from_ids(
    meta_all[[resolved_id_col]]
  )
  context_conflict <- !is.na(inferred_context) &
    inferred_context != reviewed_context
  if (any(ploidy_conflict) || any(context_conflict)) {
    stop(
      "Reviewed Seurat Ploidy/TN fields contradict ID/sample inference",
      call. = FALSE
    )
  }
  meta_all$Ploidy <- factor(reviewed_ploidy, levels = c("2N", "4N"))
  meta_all$TN <- factor(
    reviewed_context,
    levels = c("CellLine", "Tumor")
  )
  if (!("clusters" %in% colnames(meta_all))) {
    meta_all$clusters <- as.character(meta_all[[cluster_col]])
  }

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

  attr(meta_all, "sample_folder_col") <- sample_folder_col
  attr(meta_all, "cluster_col") <- cluster_col
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
    parser.add_argument(
        "--velocity-embedding-output",
        default="",
        help=(
            "Optional tab-delimited per-cell UMAP velocity-vector output. "
            "The embedding is calculated from the same stochastic velocity graph "
            "and reviewed Seurat UMAP used for velocity pseudotime."
        ),
    )
    parser.add_argument("--cluster-col", default="clusters")
    parser.add_argument("--mode", default="stochastic", choices=["deterministic", "stochastic", "dynamical"])
    parser.add_argument("--min-shared-counts", type=int, default=20)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--n-jobs", type=int, default=4)
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
        "Dose",
        "Ploidy",
        "TN",
        args.cluster_col,
        "clusters",
        "sample",
        "velocity_pseudotime",
    ]
    metric_cols = [x for x in dict.fromkeys(metric_cols) if x in adata.obs.columns]
    metrics = adata.obs[metric_cols].copy()
    metrics.insert(0, "cell", adata.obs_names.astype(str))
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    metrics.to_csv(output, index=False)


def write_velocity_embedding(adata, args, output, pd, np, scv):
    if not output:
        return
    if "X_umap" not in adata.obsm:
        raise SystemExit(
            "Cannot export the velocity embedding because the reviewed Seurat "
            "UMAP is absent from adata.obsm['X_umap']."
        )
    scv.tl.velocity_embedding(adata, basis="umap")
    if "velocity_umap" not in adata.obsm:
        raise SystemExit(
            "scVelo did not create adata.obsm['velocity_umap']."
        )
    coords = np.asarray(adata.obsm["X_umap"], dtype=float)
    vectors = np.asarray(adata.obsm["velocity_umap"], dtype=float)
    if (
        coords.ndim != 2
        or vectors.ndim != 2
        or coords.shape != vectors.shape
        or coords.shape[1] != 2
    ):
        raise SystemExit(
            "Expected matched two-dimensional X_umap and velocity_umap arrays."
        )
    if not np.isfinite(coords).all() or not np.isfinite(vectors).all():
        raise SystemExit("UMAP coordinates and velocity vectors must be finite.")
    root_values = (
        adata.obs[args.root_key_name].astype(float).to_numpy()
        if args.root_key_name in adata.obs.columns
        else np.zeros(adata.n_obs, dtype=float)
    )
    cluster_values = (
        adata.obs[args.cluster_col].astype(str).to_numpy()
        if args.cluster_col in adata.obs.columns
        else np.repeat("", adata.n_obs)
    )
    context_values = (
        adata.obs["TN"].astype(str).to_numpy()
        if "TN" in adata.obs.columns
        else np.repeat("", adata.n_obs)
    )
    sample_values = (
        adata.obs["sample"].astype(str).to_numpy()
        if "sample" in adata.obs.columns
        else np.repeat("", adata.n_obs)
    )
    table = pd.DataFrame(
        {
            "cell": adata.obs_names.astype(str),
            "UMAP_1": coords[:, 0],
            "UMAP_2": coords[:, 1],
            "velocity_UMAP_1": vectors[:, 0],
            "velocity_UMAP_2": vectors[:, 1],
            "velocity_pseudotime": adata.obs["velocity_pseudotime"].astype(float).to_numpy(),
            "cluster": cluster_values,
            "context": context_values,
            "sample": sample_values,
            "is_root": root_values > 0,
        }
    )
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, sep="\t", index=False)


def main():
    args = parse_args()
    np, pd, ad, sc, scv = import_required()
    scv.settings.verbosity = 3
    adata, input_paths = read_velocity_inputs(args.input, scv, sc, ad)
    adata.var_names_make_unique()
    ensure_velocity_layers(adata)
    meta = load_metadata(args.metadata, pd)
    matched_meta = match_cells(adata, meta, pd)
    if matched_meta.shape[0] != meta.shape[0]:
        raise SystemExit(
            f"Matched {matched_meta.shape[0]} of {meta.shape[0]} Seurat cells "
            "to velocity observations; the Figure 7 contract requires one "
            "velocity observation for every deposited Seurat cell."
        )
    if matched_meta["velocity_obs"].astype(str).duplicated().any():
        raise SystemExit(
            "Cell matching reused a velocity observation for more than one "
            "Seurat cell."
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
    set_manual_terminal_states(adata, args)
    try:
        run_velocity_pseudotime(adata, args, scv)
    except Exception as exc:
        if args.root_clusters or args.end_clusters:
            raise
        print(f"WARNING: velocity_pseudotime failed: {exc}", file=sys.stderr)
    write_metrics(adata, args, args.output, pd)
    write_velocity_embedding(
        adata,
        args,
        args.velocity_embedding_output,
        pd,
        np,
        scv,
    )


if __name__ == "__main__":
    main()
)---"
}

write_embedded_scvelo_worker <- function(path) {
  writeLines(embedded_scvelo_worker(), con = path)
  Sys.chmod(path, mode = "0755")
  invisible(path)
}

validate_scvelo_metrics_output <- function(
  output_path,
  metadata_df,
  cluster_col
) {
  metrics <- utils::read.csv(
    output_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  expected_columns <- c(
    "cell", "Dose", "Ploidy", "TN",
    cluster_col, "sample", "velocity_pseudotime"
  )
  if (!identical(names(metrics), expected_columns)) {
    stop(
      "Unexpected scVelo metrics columns: ",
      paste(names(metrics), collapse = ","),
      call. = FALSE
    )
  }
  expected_cells <- as.character(metadata_df$cell)
  observed_cells <- as.character(metrics$cell)
  if (nrow(metrics) != nrow(metadata_df) ||
      anyNA(observed_cells) ||
      any(!nzchar(observed_cells)) ||
      anyDuplicated(observed_cells) ||
      !setequal(observed_cells, expected_cells)) {
    stop(
      "scVelo output must contain exactly one row for every deposited ",
      "Seurat cell",
      call. = FALSE
    )
  }
  pseudotime <- suppressWarnings(as.numeric(metrics$velocity_pseudotime))
  if (any(!is.finite(pseudotime)) ||
      any(pseudotime < 0 | pseudotime > 1)) {
    stop(
      "scVelo output velocity_pseudotime must be finite and within [0,1] ",
      "for every cell",
      call. = FALSE
    )
  }
  order_index <- match(expected_cells, observed_cells)
  metrics <- metrics[order_index, , drop = FALSE]
  for (field in c("Dose", "Ploidy", "TN", cluster_col, "sample")) {
    if (!(field %in% names(metadata_df)) ||
        !identical(
          as.character(metrics[[field]]),
          as.character(metadata_df[[field]])
        )) {
      stop(
        "scVelo output field differs from deposited Seurat metadata: ",
        field,
        call. = FALSE
      )
    }
  }
  invisible(metrics)
}

validate_scvelo_velocity_embedding_output <- function(
  output_path,
  metadata_df,
  cluster_col,
  root_clusters = "6"
) {
  embedding <- utils::read.delim(
    output_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA", "NaN")
  )
  expected_columns <- c(
    "cell", "UMAP_1", "UMAP_2", "velocity_UMAP_1",
    "velocity_UMAP_2", "velocity_pseudotime", "cluster",
    "context", "sample", "is_root"
  )
  if (!identical(names(embedding), expected_columns)) {
    stop(
      "Unexpected scVelo velocity-embedding columns: ",
      paste(names(embedding), collapse = ","),
      call. = FALSE
    )
  }
  expected_cells <- as.character(metadata_df$cell)
  observed_cells <- as.character(embedding$cell)
  if (nrow(embedding) != nrow(metadata_df) ||
      anyNA(observed_cells) || any(!nzchar(observed_cells)) ||
      anyDuplicated(observed_cells) ||
      !setequal(observed_cells, expected_cells)) {
    stop(
      "scVelo velocity embedding must contain exactly one row for every ",
      "deposited Seurat cell",
      call. = FALSE
    )
  }
  numeric_columns <- c(
    "UMAP_1", "UMAP_2", "velocity_UMAP_1",
    "velocity_UMAP_2", "velocity_pseudotime"
  )
  numeric_values <- lapply(
    embedding[numeric_columns],
    function(value) suppressWarnings(as.numeric(value))
  )
  if (any(!vapply(numeric_values, function(value) {
    all(is.finite(value))
  }, logical(1L)))) {
    stop(
      "scVelo velocity embedding coordinates, vectors, and pseudotime ",
      "must be finite",
      call. = FALSE
    )
  }
  pseudotime <- numeric_values$velocity_pseudotime
  if (any(pseudotime < 0 | pseudotime > 1)) {
    stop("scVelo velocity-embedding pseudotime must be within [0,1]", call. = FALSE)
  }
  order_index <- match(expected_cells, observed_cells)
  embedding <- embedding[order_index, , drop = FALSE]
  metadata_cluster <- as.character(metadata_df[[cluster_col]])
  if (!identical(as.character(embedding$cluster), metadata_cluster) ||
      !identical(as.character(embedding$context), as.character(metadata_df$TN)) ||
      !identical(as.character(embedding$sample), as.character(metadata_df$sample))) {
    stop(
      "scVelo velocity embedding metadata differs from deposited Seurat metadata",
      call. = FALSE
    )
  }
  root_values <- tolower(as.character(embedding$is_root)) %in%
    c("true", "t", "1")
  expected_roots <- metadata_cluster %in% split_values(root_clusters)
  if (!identical(root_values, expected_roots) || !any(root_values)) {
    stop(
      "scVelo velocity embedding root flags do not exactly identify the ",
      "configured root cluster(s)",
      call. = FALSE
    )
  }
  embedding$is_root <- root_values
  invisible(embedding)
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
  si_config <- config$si_figures
  if (is.null(si_config)) stop("Figure 7 config is missing si_figures", call. = FALSE)

  input_root_arg <- arg_value(args, "input_root", NULL)
  seurat_rds_arg <- arg_value(args, "seurat_rds", NULL)
  loom_root_arg <- arg_value(args, "loom_root", NULL)
  if (is.null(input_root_arg) && (is.null(seurat_rds_arg) || is.null(loom_root_arg))) {
    usage()
    stop("Pass both --seurat_rds and --loom_root, or use legacy --input_root.", call. = FALSE)
  }
  input_root <- if (is.null(input_root_arg)) NULL else normalizePath(input_root_arg, mustWork = TRUE)

  output_path <- resolve_output_path(arg_value(args, "output", "Data/in-vivo/scvelo_cell_metrics.csv"), repo_root)
  velocity_embedding_arg <- arg_value(args, "velocity_embedding_output", "")
  velocity_embedding_path <- if (nzchar(velocity_embedding_arg)) {
    resolve_output_path(velocity_embedding_arg, repo_root)
  } else {
    ""
  }
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (nzchar(velocity_embedding_path) && !dir.exists(dirname(velocity_embedding_path))) {
    dir.create(dirname(velocity_embedding_path), recursive = TRUE, showWarnings = FALSE)
  }
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
  root_clusters <- arg_value(args, "root_clusters", "6")
  end_clusters <- arg_value(args, "end_clusters", "")
  mode <- arg_value(args, "mode", "stochastic")
  pca_prefix <- arg_value(args, "pca_prefix", "PCA_")
  expected_doses <- split_values(arg_value(args, "expected_doses", "0mg/kg,30mg/kg,120mg/kg"))

  metadata_df <- build_all_cells_metadata(
    seurat_rds = seurat_rds,
    sample_folder_col = arg_value(args, "sample_folder_col", "sample_folder"),
    dose_col = arg_value(args, "dose_col", "Dose"),
    ploidy_col = arg_value(
      args,
      "ploidy_col",
      as.character(si_config$ploidy_field)
    ),
    context_col = arg_value(
      args,
      "context_col",
      as.character(si_config$context_field)
    ),
    id_col = arg_value(args, "id_col", "ID"),
    cluster_col = arg_value(args, "cluster_col", as.character(si_config$cluster_id_field)),
    expected_doses = expected_doses,
    n_pcs = n_pcs,
    pca_prefix = pca_prefix,
    umap_reduction = arg_value(args, "umap_reduction", as.character(si_config$umap_reduction)),
    pca_reduction = arg_value(args, "pca_reduction", as.character(si_config$pca_reduction))
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
    "--use-metadata-pca",
    "--pca-prefix", pca_prefix,
    "--root-clusters", root_clusters
  )
  if (nzchar(velocity_embedding_path)) {
    scvelo_args <- c(
      scvelo_args,
      "--velocity-embedding-output", velocity_embedding_path
    )
  }
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
  validate_scvelo_metrics_output(
    output_path,
    metadata_df,
    attr(metadata_df, "cluster_col")
  )
  if (nzchar(velocity_embedding_path)) {
    if (!file.exists(velocity_embedding_path)) {
      stop(
        "scVelo completed but did not write velocity embedding: ",
        velocity_embedding_path,
        call. = FALSE
      )
    }
    validate_scvelo_velocity_embedding_output(
      velocity_embedding_path,
      metadata_df,
      attr(metadata_df, "cluster_col"),
      root_clusters
    )
  }

  message("Wrote scVelo cell metrics: ", output_path)
  if (nzchar(velocity_embedding_path)) {
    message("Wrote scVelo UMAP velocity embedding: ", velocity_embedding_path)
  }
  if (!keep_work) {
    unlink(work_dir, recursive = TRUE, force = TRUE)
  } else {
    message("Kept work directory: ", normalizePath(work_dir, mustWork = FALSE))
  }
}

if (identical(environment(), globalenv())) main()
