#!/usr/bin/env python3

import argparse
import json
import re
import sys
import traceback
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Scanpy PAGA for one trajectory subset using Seurat-exported metadata."
    )
    parser.add_argument("--metadata", required=True, help="CSV exported by 04_trajectory.R for one trajectory analysis.")
    parser.add_argument("--output-dir", required=True, help="Output directory for this PAGA analysis.")
    parser.add_argument("--analysis-label", required=True, help="Human-readable trajectory analysis label.")
    parser.add_argument("--basis", default="umap", help="Embedding basis name used for output plots.")
    parser.add_argument("--cluster-col", default="cluster_final", help="Categorical column used as PAGA groups.")
    parser.add_argument("--color-cols", default="cluster_final,Dose,TN,Ploidy", help="Comma-separated obs columns to plot on UMAP.")
    parser.add_argument("--pca-prefix", default="PCA_", help="Prefix for metadata PCA columns, e.g. PCA_1.")
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-cells", type=int, default=20)
    parser.add_argument("--paga-threshold", type=float, default=0.03)
    parser.add_argument("--root-clusters", default="", help="Optional comma-separated cluster labels used to root DPT.")
    parser.add_argument("--dpi", type=int, default=300)
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
    except ImportError as exc:
        raise SystemExit(
            "Missing Python package for PAGA analysis: "
            f"{exc}. Install scanpy, anndata, pandas, numpy, scipy, scikit-learn, "
            "matplotlib, python-igraph, and leidenalg."
        ) from exc

    return plt, np, pd, ad, sc


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


def safe_filename(value):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_") or "value"


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2))


def save_plot_both(plt, output_stub, dpi):
    output_stub = Path(output_stub)
    output_stub.parent.mkdir(parents=True, exist_ok=True)
    try:
        plt.savefig(str(output_stub.with_suffix(".pdf")), bbox_inches="tight")
    except Exception:
        err_file = output_stub.with_suffix(".pdf_error.txt")
        err_file.write_text(traceback.format_exc())
        print(
            f"WARNING: Failed to save PDF plot {output_stub.with_suffix('.pdf')}; see {err_file}",
            file=sys.stderr,
        )
    try:
        plt.savefig(str(output_stub.with_suffix(".png")), bbox_inches="tight", dpi=dpi)
    except Exception:
        err_file = output_stub.with_suffix(".png_error.txt")
        err_file.write_text(traceback.format_exc())
        print(
            f"WARNING: Failed to save PNG plot {output_stub.with_suffix('.png')}; see {err_file}",
            file=sys.stderr,
        )
    finally:
        plt.close("all")


def load_metadata(metadata_path, args, pd, np):
    meta = pd.read_csv(metadata_path)
    required = {"cell", "UMAP_1", "UMAP_2", args.cluster_col}
    missing = required.difference(set(meta.columns))
    if missing:
        raise SystemExit(
            f"Metadata file is missing required column(s): {', '.join(sorted(missing))}. "
            f"File: {metadata_path}"
        )

    pca_cols = metadata_pca_columns(meta.columns, args.pca_prefix)
    if not pca_cols:
        raise SystemExit(
            f"No PCA columns matched prefix {args.pca_prefix!r} in metadata: {metadata_path}"
        )
    pca_cols = pca_cols[: args.n_pcs]

    meta = meta.copy()
    meta["cell"] = meta["cell"].astype(str)
    keep = meta["cell"].notna()
    keep &= meta[args.cluster_col].notna()
    keep &= meta[["UMAP_1", "UMAP_2"]].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    keep &= meta[pca_cols].apply(pd.to_numeric, errors="coerce").notna().all(axis=1)
    meta = meta.loc[keep].copy()
    if meta.shape[0] == 0:
        raise SystemExit(f"No usable cells remained after filtering metadata: {metadata_path}")

    meta = meta.drop_duplicates(subset=["cell"], keep="first")
    pca = meta[pca_cols].astype(float).to_numpy()
    umap = meta[["UMAP_1", "UMAP_2"]].astype(float).to_numpy()
    if not np.isfinite(pca).all() or not np.isfinite(umap).all():
        raise SystemExit("Metadata PCA/UMAP contains non-finite values after filtering.")

    return meta, pca, umap, pca_cols


def prepare_obs_columns(adata, columns):
    for col in columns:
        if not col or col not in adata.obs.columns:
            continue
        values = adata.obs[col]
        if getattr(values.dtype, "kind", "") in {"b", "i", "u", "f", "c"}:
            continue
        values = values.astype(str).where(values.notna(), "NA")
        adata.obs[col] = values.astype("category")


def build_adata(meta, pca, umap, args, ad, np):
    obs = meta.copy()
    obs.index = obs["cell"].astype(str)
    obs.index.name = None
    adata = ad.AnnData(X=np.zeros((obs.shape[0], 1), dtype="float32"), obs=obs)
    adata.obsm["X_pca"] = pca.astype("float32", copy=False)
    adata.obsm[f"X_{args.basis}"] = umap.astype("float32", copy=False)
    adata.obs[args.cluster_col] = adata.obs[args.cluster_col].astype(str).astype("category")
    return adata


def write_paga_tables(adata, args, pd, output_dir):
    groups = adata.obs[args.cluster_col].cat.categories.astype(str).tolist()
    conn = adata.uns["paga"]["connectivities"]
    conn_dense = conn.toarray() if hasattr(conn, "toarray") else conn

    conn_df = pd.DataFrame(conn_dense, index=groups, columns=groups)
    conn_df.index.name = args.cluster_col
    conn_df.to_csv(output_dir / "paga_connectivities.csv")

    edge_rows = []
    for i, group_i in enumerate(groups):
        for j in range(i + 1, len(groups)):
            weight = float(conn_dense[i, j])
            if weight > 0:
                edge_rows.append(
                    {
                        "group_1": group_i,
                        "group_2": groups[j],
                        "connectivity": weight,
                        "above_threshold": bool(weight >= args.paga_threshold),
                    }
                )
    edge_df = pd.DataFrame(edge_rows, columns=["group_1", "group_2", "connectivity", "above_threshold"])
    edge_df.to_csv(output_dir / "paga_edges.csv", index=False)
    return edge_df, conn_df


def plot_scanpy_umaps(adata, args, plt, sc):
    plots_dir = Path(args.output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    color_cols = [x for x in parse_comma_arg(args.color_cols) if x in adata.obs.columns]
    if args.cluster_col in adata.obs.columns and args.cluster_col not in color_cols:
        color_cols.insert(0, args.cluster_col)
    prepare_obs_columns(adata, color_cols + ["dpt_pseudotime"])

    for col in color_cols:
        sc.pl.embedding(
            adata,
            basis=args.basis,
            color=col,
            title=f"{args.analysis_label} | {col}",
            frameon=False,
            show=False,
        )
        save_plot_both(plt, plots_dir / f"umap_{safe_filename(col)}", args.dpi)

    if "dpt_pseudotime" in adata.obs.columns:
        sc.pl.embedding(
            adata,
            basis=args.basis,
            color="dpt_pseudotime",
            color_map="viridis",
            title=f"{args.analysis_label} | DPT pseudotime",
            frameon=False,
            show=False,
        )
        save_plot_both(plt, plots_dir / "umap_dpt_pseudotime", args.dpi)


def plot_paga_graph(adata, args, plt, sc):
    plots_dir = Path(args.output_dir) / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    sc.pl.paga(
        adata,
        threshold=args.paga_threshold,
        color=args.cluster_col,
        title=f"{args.analysis_label} PAGA",
        frameon=False,
        show=False,
    )
    save_plot_both(plt, plots_dir / "paga_graph", args.dpi)


def plot_paga_umap_overlay(adata, args, plt, np, output_dir):
    plots_dir = output_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    coords = adata.obsm[f"X_{args.basis}"]
    clusters = adata.obs[args.cluster_col].astype(str)
    groups = adata.obs[args.cluster_col].cat.categories.astype(str).tolist()
    conn = adata.uns["paga"]["connectivities"]
    conn_dense = conn.toarray() if hasattr(conn, "toarray") else conn

    cmap = plt.get_cmap("tab20")
    group_colors = {group: cmap(i % 20) for i, group in enumerate(groups)}
    centroids = {}
    for group in groups:
        mask = (clusters == group).to_numpy()
        if mask.any():
            centroids[group] = np.median(coords[mask, :], axis=0)

    fig, ax = plt.subplots(figsize=(6.5, 5.6))
    for group in groups:
        mask = (clusters == group).to_numpy()
        if not mask.any():
            continue
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=4,
            alpha=0.55,
            linewidths=0,
            color=group_colors[group],
            label=group,
        )

    max_weight = float(conn_dense.max()) if conn_dense.size else 0.0
    if max_weight > 0:
        for i, group_i in enumerate(groups):
            for j in range(i + 1, len(groups)):
                weight = float(conn_dense[i, j])
                if weight < args.paga_threshold:
                    continue
                group_j = groups[j]
                if group_i not in centroids or group_j not in centroids:
                    continue
                p0 = centroids[group_i]
                p1 = centroids[group_j]
                ax.plot(
                    [p0[0], p1[0]],
                    [p0[1], p1[1]],
                    color="#222222",
                    alpha=0.75,
                    linewidth=0.4 + 4.0 * weight / max_weight,
                    zorder=3,
                )

    for group, xy in centroids.items():
        ax.text(
            xy[0],
            xy[1],
            str(group),
            ha="center",
            va="center",
            fontsize=8,
            fontweight="bold",
            color="black",
            bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.75},
            zorder=4,
        )

    ax.set_title(f"{args.analysis_label} PAGA edges on UMAP")
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    if len(groups) <= 25:
        ax.legend(
            title=args.cluster_col,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            markerscale=2.0,
            fontsize=7,
        )
    save_plot_both(plt, plots_dir / "paga_umap_overlay", args.dpi)


def run_optional_dpt(adata, args, sc, np):
    root_clusters = parse_comma_arg(args.root_clusters)
    if not root_clusters:
        return {"dpt_status": "skipped_no_root", "root_clusters": [], "n_root_cells": 0}

    clusters = adata.obs[args.cluster_col].astype(str)
    root_mask = clusters.isin(root_clusters).to_numpy()
    n_root = int(root_mask.sum())
    if n_root == 0:
        observed = sorted(clusters.dropna().unique().tolist())
        return {
            "dpt_status": "skipped_missing_root",
            "root_clusters": root_clusters,
            "n_root_cells": 0,
            "observed_clusters": observed,
        }

    try:
        adata.uns["iroot"] = int(np.flatnonzero(root_mask)[0])
        sc.tl.diffmap(adata)
        sc.tl.dpt(adata)
    except Exception as exc:
        return {
            "dpt_status": "failed",
            "root_clusters": root_clusters,
            "n_root_cells": n_root,
            "message": str(exc),
        }

    return {
        "dpt_status": "completed",
        "root_clusters": root_clusters,
        "n_root_cells": n_root,
        "iroot": int(adata.uns["iroot"]),
    }


def write_cell_metrics(adata, args, output_dir):
    cols = [
        "cell",
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
        "dpt_pseudotime",
    ]
    cols = [x for x in dict.fromkeys(cols) if x in adata.obs.columns]
    metrics = adata.obs[cols].copy()
    if "cell" not in metrics.columns:
        metrics.insert(0, "cell", adata.obs_names.astype(str))
    metrics.to_csv(output_dir / "paga_cell_metrics.csv", index=False)


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plt, np, pd, ad, sc = import_required()
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(
        dpi=args.dpi,
        dpi_save=args.dpi,
        fontsize=10,
        frameon=False,
        facecolor="white",
    )

    meta, pca, umap, pca_cols = load_metadata(args.metadata, args, pd, np)
    adata = build_adata(meta, pca, umap, args, ad, np)
    groups = adata.obs[args.cluster_col].cat.categories.astype(str).tolist()

    effective_neighbors = min(args.n_neighbors, max(1, adata.n_obs - 1))
    effective_pcs = min(args.n_pcs, pca.shape[1])
    report = {
        "analysis_label": args.analysis_label,
        "metadata": str(args.metadata),
        "n_cells": int(adata.n_obs),
        "n_groups": int(len(groups)),
        "groups": groups,
        "cluster_col": args.cluster_col,
        "pca_prefix": args.pca_prefix,
        "pca_columns": pca_cols,
        "requested_n_pcs": int(args.n_pcs),
        "effective_n_pcs": int(effective_pcs),
        "requested_n_neighbors": int(args.n_neighbors),
        "effective_n_neighbors": int(effective_neighbors),
        "paga_threshold": float(args.paga_threshold),
    }

    if adata.n_obs < args.min_cells:
        report["status"] = "skipped_too_few_cells"
        write_json(output_dir / "paga_run_summary.json", report)
        write_cell_metrics(adata, args, output_dir)
        adata.write_h5ad(output_dir / "paga_result.h5ad")
        return
    if len(groups) < 2:
        report["status"] = "skipped_single_group"
        write_json(output_dir / "paga_run_summary.json", report)
        write_cell_metrics(adata, args, output_dir)
        adata.write_h5ad(output_dir / "paga_result.h5ad")
        return

    sc.pp.neighbors(adata, n_neighbors=effective_neighbors, n_pcs=effective_pcs, use_rep="X_pca")
    sc.tl.paga(adata, groups=args.cluster_col)
    edge_df, _ = write_paga_tables(adata, args, pd, output_dir)
    dpt_report = run_optional_dpt(adata, args, sc, np)
    write_cell_metrics(adata, args, output_dir)
    adata.write_h5ad(output_dir / "paga_result.h5ad")

    try:
        plot_scanpy_umaps(adata, args, plt, sc)
        plot_paga_graph(adata, args, plt, sc)
        plot_paga_umap_overlay(adata, args, plt, np, output_dir)
    except Exception:
        plot_error_file = output_dir / "plotting_error.txt"
        plot_error_file.write_text(traceback.format_exc())
        print(
            f"WARNING: PAGA plotting failed after tables were written; see {plot_error_file}",
            file=sys.stderr,
        )

    report["status"] = "completed"
    report["n_edges"] = int(edge_df.shape[0])
    report["n_edges_above_threshold"] = int(edge_df["above_threshold"].sum()) if edge_df.shape[0] else 0
    report.update(dpt_report)
    write_json(output_dir / "paga_run_summary.json", report)


if __name__ == "__main__":
    main()
