#!/usr/bin/env python3

import argparse
import inspect
import json
import os
import re
import sys
import traceback
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run scVelo for one trajectory subset using Seurat-exported metadata."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input .loom or .h5ad file with spliced/unspliced layers. Multiple files may be comma-separated.",
    )
    parser.add_argument("--metadata", required=True, help="CSV exported by 04_trajectory.R for one trajectory analysis.")
    parser.add_argument("--output-dir", required=True, help="Output directory for this trajectory analysis.")
    parser.add_argument("--dose-label", required=True, help="Human-readable trajectory analysis label.")
    parser.add_argument("--basis", default="umap", help="Embedding basis to use for velocity plots.")
    parser.add_argument("--cluster-col", default="cluster_final", help="Categorical column for cluster coloring.")
    parser.add_argument("--color-cols", default="cluster_final,Dose,sample_type", help="Comma-separated metadata columns to plot.")
    parser.add_argument("--shape-col", default="", help="Optional metadata column used to vary point shape in metric UMAPs.")
    parser.add_argument("--shape-order", default="", help="Comma-separated order for --shape-col categories.")
    parser.add_argument("--shape-markers", default="", help="Comma-separated matplotlib markers for --shape-col categories.")
    parser.add_argument("--cell-map", default="", help="Optional CSV with seurat_cell and velocity_cell columns.")
    parser.add_argument("--mode", default="stochastic", choices=["deterministic", "stochastic", "dynamical"])
    parser.add_argument("--min-shared-counts", type=int, default=20)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--min-matched-cells", type=int, default=50)
    parser.add_argument("--stream-density", type=float, default=1.2)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--no-use-r-umap", action="store_true", help="Do not overwrite X_umap with Seurat-exported UMAP.")
    parser.add_argument(
        "--use-metadata-pca",
        action="store_true",
        help="Use PCA coordinates exported in metadata to build Scanpy neighbors before scVelo moments.",
    )
    parser.add_argument("--pca-prefix", default="PCA_", help="Prefix for metadata PCA columns, e.g. PCA_1.")
    parser.add_argument("--root-clusters", default="", help="Comma-separated cluster labels to force as velocity pseudotime roots.")
    parser.add_argument("--end-clusters", default="", help="Comma-separated cluster labels to force as velocity pseudotime endpoints.")
    parser.add_argument("--root-key-name", default="root_cells_manual", help="adata.obs column for manual root cells.")
    parser.add_argument("--end-key-name", default="end_points_manual", help="adata.obs column for manual endpoint cells.")
    return parser.parse_args()


def import_required():
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import anndata as ad
        import scanpy as sc
        import scvelo as scv
    except ImportError as exc:
        raise SystemExit(
            "Missing Python package for velocity analysis: "
            f"{exc}. Install scvelo, scanpy, anndata, pandas, numpy, matplotlib."
        ) from exc

    return plt, np, pd, sc, scv, ad


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


def choose_metadata_sample_column(meta):
    for col in ["sample_folder", "Sequencing.IDs", "sample", "orig.ident"]:
        if col in meta.columns:
            return col
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
            keys.append(f"cell:{candidate}")
            keys.append(f"cell_upper:{candidate.upper()}")
    if barcode:
        keys.append(f"barcode:{barcode}")
        if inferred_sample:
            keys.append(f"sample_barcode:{inferred_sample}|{barcode}")
    return keys


def build_unique_key_to_obs(obs_names, known_samples):
    key_to_obs = {}
    duplicated = set()
    for obs in obs_names:
        for key in make_match_keys(obs, known_samples=known_samples):
            if key in key_to_obs and key_to_obs[key] != obs:
                duplicated.add(key)
            elif key not in duplicated:
                key_to_obs[key] = obs
    for key in duplicated:
        key_to_obs.pop(key, None)
    return key_to_obs


def read_velocity_input(path, scv, sc):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".loom":
        if hasattr(scv, "read"):
            return scv.read(str(path), cache=False)
        if hasattr(scv, "read_loom"):
            return scv.read_loom(str(path), cache=False)
        return sc.read_loom(str(path))
    if suffix == ".h5ad":
        return sc.read_h5ad(str(path))
    raise SystemExit(f"Unsupported velocity input format: {path}. Expected .loom or .h5ad.")


def parse_input_paths(input_arg):
    paths = [x.strip() for x in str(input_arg).split(",") if x.strip()]
    if not paths:
        raise SystemExit("--input did not contain any usable file paths.")
    return [Path(x) for x in paths]


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


def load_metadata(metadata_path, cell_map_path, pd):
    meta = pd.read_csv(metadata_path)
    if "cell" not in meta.columns:
        raise SystemExit(f"Metadata file is missing required column 'cell': {metadata_path}")
    meta["cell"] = meta["cell"].astype(str)

    if cell_map_path:
        cell_map = pd.read_csv(cell_map_path)
        required = {"seurat_cell", "velocity_cell"}
        missing = required.difference(cell_map.columns)
        if missing:
            raise SystemExit(
                f"Cell map is missing required columns: {', '.join(sorted(missing))}"
            )
        cell_map = cell_map[["seurat_cell", "velocity_cell"]].copy()
        cell_map["seurat_cell"] = cell_map["seurat_cell"].astype(str)
        cell_map["velocity_cell"] = cell_map["velocity_cell"].astype(str)
        meta = meta.merge(cell_map, left_on="cell", right_on="seurat_cell", how="left")

    return meta


def match_cells(adata, meta, pd):
    obs_names = pd.Index(adata.obs_names.astype(str))
    meta = meta.copy()
    meta["velocity_cell"] = meta.get("velocity_cell", meta["cell"])
    meta["velocity_cell"] = meta["velocity_cell"].fillna(meta["cell"]).astype(str)

    exact_keep = meta["velocity_cell"].isin(obs_names)
    if exact_keep.sum() > 0:
        matched = meta.loc[exact_keep].copy()
        matched["matched_velocity_cell"] = matched["velocity_cell"]
        matched["match_method"] = "exact"
        return matched

    known_samples = []
    sample_col = choose_metadata_sample_column(meta)
    if sample_col:
        known_samples = [
            normalize_sample_value(x)
            for x in meta[sample_col].dropna().astype(str).unique().tolist()
        ]
        known_samples = sorted([x for x in known_samples if x], key=len, reverse=True)

    key_to_obs = build_unique_key_to_obs(obs_names.astype(str).tolist(), known_samples)
    matched_rows = []
    matched_velocity_cells = []
    match_keys = []
    for _, row in meta.iterrows():
        sample = row[sample_col] if sample_col else ""
        candidates = []
        candidates.extend(make_match_keys(row["velocity_cell"], sample=sample, known_samples=known_samples))
        candidates.extend(make_match_keys(row["cell"], sample=sample, known_samples=known_samples))
        seen = set()
        for key in candidates:
            if key in seen:
                continue
            seen.add(key)
            if key in key_to_obs:
                matched_rows.append(row)
                matched_velocity_cells.append(key_to_obs[key])
                match_keys.append(key)
                break

    if matched_rows:
        matched = pd.DataFrame(matched_rows).copy()
        matched["matched_velocity_cell"] = matched_velocity_cells
        matched["match_method"] = "normalized_barcode"
        matched["match_key"] = match_keys
        return matched

    matched = meta.iloc[0:0].copy()
    matched["matched_velocity_cell"] = []
    matched["match_method"] = []
    return matched


def set_metadata_and_umap(adata, matched_meta, args, np):
    velocity_cells = matched_meta["matched_velocity_cell"].astype(str).tolist()
    adata_sub = adata[velocity_cells].copy()
    adata_sub.obs_names = matched_meta["cell"].astype(str).tolist()

    for col in matched_meta.columns:
        if col in {"matched_velocity_cell", "match_method"}:
            continue
        adata_sub.obs[col] = matched_meta[col].values

    adata_sub.obs["velocity_cell"] = velocity_cells
    adata_sub.obs["match_method"] = matched_meta["match_method"].values

    if not args.no_use_r_umap and {"UMAP_1", "UMAP_2"}.issubset(matched_meta.columns):
        coords = matched_meta[["UMAP_1", "UMAP_2"]].astype(float).to_numpy()
        if np.isfinite(coords).all():
            adata_sub.obsm[f"X_{args.basis}"] = coords

    return adata_sub


def ensure_velocity_layers(adata):
    required_layers = {"spliced", "unspliced"}
    missing = required_layers.difference(set(adata.layers.keys()))
    if missing:
        raise SystemExit(
            "Velocity input is missing required layer(s): "
            f"{', '.join(sorted(missing))}. A true RNA velocity run needs spliced/unspliced counts."
        )


def save_plot_both(plt, output_stub, dpi):
    output_stub = Path(output_stub)
    output_stub.parent.mkdir(parents=True, exist_ok=True)
    try:
        plt.savefig(str(output_stub.with_suffix(".pdf")), bbox_inches="tight")
    except Exception:
        err_file = output_stub.with_suffix(".pdf_error.txt")
        err_file.write_text(traceback.format_exc())
        print(
            f"WARNING: Failed to save PDF plot {output_stub.with_suffix('.pdf')}; "
            f"see {err_file}",
            file=sys.stderr,
        )
    try:
        plt.savefig(str(output_stub.with_suffix(".png")), bbox_inches="tight", dpi=dpi)
    except Exception:
        err_file = output_stub.with_suffix(".png_error.txt")
        err_file.write_text(traceback.format_exc())
        print(
            f"WARNING: Failed to save PNG plot {output_stub.with_suffix('.png')}; "
            f"see {err_file}",
            file=sys.stderr,
        )
    finally:
        plt.close("all")


def parse_comma_arg(value):
    return [x.strip() for x in str(value).split(",") if x.strip()]


def metadata_pca_columns(columns, prefix):
    hits = []
    pattern = re.compile(rf"^{re.escape(prefix)}(\d+)$")
    for col in columns:
        match = pattern.match(str(col))
        if match:
            hits.append((int(match.group(1)), col))
    hits.sort(key=lambda x: x[0])
    return [col for _, col in hits]


def apply_metadata_pca_neighbors(adata, args, sc, np):
    if not args.use_metadata_pca:
        return {"use_metadata_pca": False, "metadata_pca_dims": 0}

    pca_cols = metadata_pca_columns(adata.obs.columns, args.pca_prefix)
    if not pca_cols:
        raise SystemExit(
            f"--use-metadata-pca was set, but no metadata PCA columns matched prefix {args.pca_prefix!r}."
        )
    pca_cols = pca_cols[: args.n_pcs]
    pca = adata.obs[pca_cols].astype(float).to_numpy()
    if not np.isfinite(pca).all():
        bad_cells = adata.obs_names[~np.isfinite(pca).all(axis=1)][:20].tolist()
        raise SystemExit(
            "Metadata PCA contains non-finite values for matched cells. "
            f"Example cells: {', '.join(map(str, bad_cells))}"
        )

    adata.obsm["X_pca"] = pca.astype("float32", copy=False)
    n_pcs = min(args.n_pcs, pca.shape[1])
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=n_pcs, use_rep="X_pca")
    return {
        "use_metadata_pca": True,
        "metadata_pca_dims": int(pca.shape[1]),
        "metadata_pca_columns": pca_cols,
        "neighbors_use_rep": "X_pca",
        "neighbors_n_pcs": int(n_pcs),
        "neighbors_n_neighbors": int(args.n_neighbors),
    }


def run_scvelo_moments(adata, args, scv):
    if args.use_metadata_pca:
        try:
            scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
        except TypeError:
            scv.pp.moments(adata)
        return
    scv.pp.moments(adata, n_pcs=args.n_pcs, n_neighbors=args.n_neighbors)


def set_manual_terminal_states(adata, args):
    root_clusters = parse_comma_arg(args.root_clusters)
    end_clusters = parse_comma_arg(args.end_clusters)
    terminal_report = {
        "root_clusters": root_clusters,
        "end_clusters": end_clusters,
        "root_key": "",
        "end_key": "",
        "n_root_cells_manual": 0,
        "n_end_cells_manual": 0,
    }
    if not root_clusters and not end_clusters:
        return terminal_report

    if args.cluster_col not in adata.obs.columns:
        raise SystemExit(
            f"Manual root/end clusters were requested, but cluster column {args.cluster_col!r} is absent."
        )
    clusters = adata.obs[args.cluster_col].astype(str)

    if root_clusters:
        root_mask = clusters.isin(root_clusters)
        n_root = int(root_mask.sum())
        if n_root == 0:
            observed = sorted(clusters.dropna().unique().tolist())
            raise SystemExit(
                f"No root cells found for {args.cluster_col} in {root_clusters}. "
                f"Observed cluster labels include: {', '.join(observed[:30])}"
            )
        adata.obs[args.root_key_name] = root_mask.astype(float).values
        terminal_report["root_key"] = args.root_key_name
        terminal_report["n_root_cells_manual"] = n_root

    if end_clusters:
        end_mask = clusters.isin(end_clusters)
        n_end = int(end_mask.sum())
        if n_end == 0:
            observed = sorted(clusters.dropna().unique().tolist())
            raise SystemExit(
                f"No endpoint cells found for {args.cluster_col} in {end_clusters}. "
                f"Observed cluster labels include: {', '.join(observed[:30])}"
            )
        adata.obs[args.end_key_name] = end_mask.astype(float).values
        terminal_report["end_key"] = args.end_key_name
        terminal_report["n_end_cells_manual"] = n_end

    return terminal_report


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


def get_shape_marker_map(adata, args):
    if not args.shape_col or args.shape_col not in adata.obs.columns:
        return [], {}

    shape_values = adata.obs[args.shape_col].astype(str)
    observed = [x for x in shape_values.dropna().unique().tolist() if x and x != "nan"]
    requested_order = parse_comma_arg(args.shape_order)
    order = [x for x in requested_order if x in observed]
    order.extend([x for x in observed if x not in order])

    default_markers = ["o", "^", "*", "s", "D", "P", "X", "v", "<", ">"]
    requested_markers = parse_comma_arg(args.shape_markers)
    markers = requested_markers + default_markers
    marker_map = {group: markers[idx % len(markers)] for idx, group in enumerate(order)}
    return order, marker_map


def plot_metric_shape_umap(adata, args, plt, np):
    if not args.shape_col or args.shape_col not in adata.obs.columns:
        return
    basis_key = f"X_{args.basis}"
    if basis_key not in adata.obsm:
        return

    order, marker_map = get_shape_marker_map(adata, args)
    if not order:
        return

    coords = adata.obsm[basis_key]
    shape_values = adata.obs[args.shape_col].astype(str)
    plots_dir = Path(args.output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    for metric in ["velocity_pseudotime", "latent_time", "velocity_confidence"]:
        if metric not in adata.obs.columns:
            continue
        values = adata.obs[metric].astype(float).to_numpy()
        finite = np.isfinite(values)
        if not finite.any():
            continue

        fig, ax = plt.subplots(figsize=(6.5, 5.5))
        norm = plt.Normalize(vmin=float(np.nanmin(values[finite])), vmax=float(np.nanmax(values[finite])))
        scatter_for_colorbar = None
        for group in order:
            idx = (shape_values == group).to_numpy() & finite
            if not idx.any():
                continue
            scatter_for_colorbar = ax.scatter(
                coords[idx, 0],
                coords[idx, 1],
                c=values[idx],
                cmap="viridis",
                norm=norm,
                marker=marker_map[group],
                s=9 if marker_map[group] != "*" else 18,
                alpha=0.85,
                linewidths=0,
                label=group,
            )
        if scatter_for_colorbar is None:
            plt.close(fig)
            continue
        ax.set_title(f"{args.dose_label} | {metric} by {args.shape_col}")
        ax.set_aspect("equal", adjustable="box")
        ax.set_axis_off()
        ax.legend(title=args.shape_col, loc="best", frameon=False, markerscale=1.6)
        fig.colorbar(scatter_for_colorbar, ax=ax, fraction=0.046, pad=0.04, label=metric)
        save_plot_both(plt, plots_dir / f"umap_{metric}_by_shape", args.dpi)


def prepare_obs_columns_for_plotting(adata, columns):
    metric_cols = {
        "velocity_length",
        "velocity_confidence",
        "velocity_confidence_transition",
        "velocity_pseudotime",
        "latent_time",
    }
    for col in columns:
        if not col or col not in adata.obs.columns or col in metric_cols:
            continue
        series = adata.obs[col]
        if getattr(series.dtype, "kind", "") in {"b", "i", "u", "f", "c"}:
            continue
        values = series.astype(str)
        values = values.where(series.notna(), "NA")
        adata.obs[col] = values.astype("category")


def plot_scvelo_outputs(adata, args, plt, scv, np):
    color_cols = [x.strip() for x in args.color_cols.split(",") if x.strip()]
    color_cols = [x for x in color_cols if x in adata.obs.columns]
    if args.cluster_col in adata.obs.columns and args.cluster_col not in color_cols:
        color_cols.insert(0, args.cluster_col)
    prepare_obs_columns_for_plotting(
        adata,
        color_cols + [args.cluster_col, args.shape_col, "trajectory_shape_group"],
    )

    plots_dir = Path(args.output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    stream_color = args.cluster_col if args.cluster_col in adata.obs.columns else None
    scv.pl.velocity_embedding_stream(
        adata,
        basis=args.basis,
        color=stream_color,
        density=args.stream_density,
        title=f"{args.dose_label} velocity stream",
        legend_loc="right margin",
        frameon=False,
        show=False,
    )
    save_plot_both(plt, plots_dir / "velocity_stream", args.dpi)

    scv.pl.velocity_embedding_grid(
        adata,
        basis=args.basis,
        color=stream_color,
        title=f"{args.dose_label} velocity grid",
        legend_loc="right margin",
        frameon=False,
        show=False,
    )
    save_plot_both(plt, plots_dir / "velocity_grid", args.dpi)

    for col in color_cols:
        scv.pl.scatter(
            adata,
            basis=args.basis,
            color=col,
            title=f"{args.dose_label} | {col}",
            frameon=False,
            show=False,
        )
        safe_col = re.sub(r"[^A-Za-z0-9._-]+", "_", col).strip("_")
        save_plot_both(plt, plots_dir / f"umap_{safe_col}", args.dpi)

    for metric in ["velocity_pseudotime", "latent_time", "velocity_confidence"]:
        if metric in adata.obs.columns:
            scv.pl.scatter(
                adata,
                basis=args.basis,
                color=metric,
                color_map="viridis",
                title=f"{args.dose_label} | {metric}",
                frameon=False,
                show=False,
            )
            save_plot_both(plt, plots_dir / f"umap_{metric}", args.dpi)

    plot_metric_shape_umap(adata, args, plt, np)


def scvelo_filter_and_normalize(adata, args, scv):
    kwargs = {"min_shared_counts": args.min_shared_counts}
    signature = inspect.signature(scv.pp.filter_and_normalize)
    if "n_top_genes" in signature.parameters:
        kwargs["n_top_genes"] = args.n_top_genes
    else:
        print(
            "WARNING: scv.pp.filter_and_normalize does not expose n_top_genes; "
            "running without n_top_genes.",
            file=sys.stderr,
        )
    try:
        scv.pp.filter_and_normalize(adata, **kwargs)
    except TypeError as exc:
        if "n_top_genes" not in str(exc):
            raise
        kwargs.pop("n_top_genes", None)
        print(
            "WARNING: scv.pp.filter_and_normalize rejected n_top_genes; "
            "retrying without n_top_genes.",
            file=sys.stderr,
        )
        scv.pp.filter_and_normalize(adata, **kwargs)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt, np, pd, sc, scv, ad = import_required()
    scv.settings.verbosity = 3
    scv.settings.set_figure_params(
        "scvelo",
        dpi=args.dpi,
        dpi_save=args.dpi,
        fontsize=10,
        frameon=False,
        facecolor="white",
    )

    adata, input_paths = read_velocity_inputs(args.input, scv, sc, ad)
    adata.var_names_make_unique()
    ensure_velocity_layers(adata)

    meta = load_metadata(args.metadata, args.cell_map, pd)
    matched_meta = match_cells(adata, meta, pd)
    report = {
        "dose_label": args.dose_label,
        "input": input_paths,
        "n_input_files": len(input_paths),
        "metadata_cells": int(meta.shape[0]),
        "velocity_cells": int(adata.n_obs),
        "matched_cells": int(matched_meta.shape[0]),
        "match_method": matched_meta["match_method"].iloc[0] if matched_meta.shape[0] else "none",
        "shape_col": args.shape_col,
        "shape_order": parse_comma_arg(args.shape_order),
        "shape_markers": parse_comma_arg(args.shape_markers),
        "use_metadata_pca": bool(args.use_metadata_pca),
        "pca_prefix": args.pca_prefix,
        "root_clusters": parse_comma_arg(args.root_clusters),
        "end_clusters": parse_comma_arg(args.end_clusters),
    }
    (output_dir / "cell_matching_report.json").write_text(json.dumps(report, indent=2))
    matched_meta.to_csv(output_dir / "matched_cells.csv", index=False)

    if matched_meta.shape[0] < args.min_matched_cells:
        sample_col = choose_metadata_sample_column(meta)
        debug = {
            "metadata_cell_examples": meta["cell"].astype(str).head(20).tolist(),
            "metadata_velocity_cell_examples": meta.get("velocity_cell", meta["cell"]).astype(str).head(20).tolist(),
            "metadata_sample_column": sample_col,
            "metadata_sample_examples": meta[sample_col].astype(str).head(20).tolist() if sample_col else [],
            "velocity_cell_examples": [str(x) for x in adata.obs_names[:20]],
            "metadata_barcode_examples": [extract_10x_barcode(x) for x in meta["cell"].astype(str).head(20).tolist()],
            "velocity_barcode_examples": [extract_10x_barcode(x) for x in adata.obs_names[:20]],
        }
        (output_dir / "cell_matching_debug.json").write_text(json.dumps(debug, indent=2))
        raise SystemExit(
            f"Only {matched_meta.shape[0]} cells matched between Seurat metadata and velocity input; "
            f"minimum required is {args.min_matched_cells}. See cell_matching_debug.json. "
            "Provide --cell-map if cell IDs differ."
        )

    adata = set_metadata_and_umap(adata, matched_meta, args, np)

    scvelo_filter_and_normalize(adata, args, scv)
    pca_neighbor_report = apply_metadata_pca_neighbors(adata, args, sc, np)
    run_scvelo_moments(adata, args, scv)

    if args.mode == "dynamical":
        scv.tl.recover_dynamics(adata, n_jobs=args.n_jobs)
        scv.tl.velocity(adata, mode="dynamical")
    else:
        scv.tl.velocity(adata, mode=args.mode)

    scv.tl.velocity_graph(adata, n_jobs=args.n_jobs)
    scv.tl.velocity_confidence(adata)
    terminal_report = set_manual_terminal_states(adata, args)
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
        "TrajectoryComparisonID",
        "TrajectoryScope",
        "TrajectoryComparison",
        "TrajectoryIdent1",
        "TrajectoryIdent2",
        "TrajectoryGroupCol",
        "TrajectoryCluster",
        "TrajectorySubsetGroup",
        "TrajectoryDose",
        "TrajectoryComparisonGroup",
        "TrajectoryOriginalGroup",
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
    metrics.to_csv(output_dir / "scvelo_cell_metrics.csv", index=False)
    adata.write_h5ad(output_dir / "scvelo_result.h5ad")

    try:
        plot_scvelo_outputs(adata, args, plt, scv, np)
    except Exception:
        plot_error_file = output_dir / "plotting_error.txt"
        plot_error_file.write_text(traceback.format_exc())
        print(
            f"WARNING: Plotting failed after scVelo metrics were written; see {plot_error_file}",
            file=sys.stderr,
        )

    report["post_filter_cells"] = int(adata.n_obs)
    report["post_filter_genes"] = int(adata.n_vars)
    report["mode"] = args.mode
    report.update(pca_neighbor_report)
    report.update(terminal_report)
    (output_dir / "run_summary.json").write_text(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
