#!/usr/bin/env python3
"""Clustered heatmaps for GDSC ploidy-enrichment panels.

This exploratory figure keeps the existing Figure 1 workflow untouched while
making clustered versions of the low- and high-ploidy enrichment heatmaps.
"""

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "gemcitabine_model_mpl"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "gemcitabine_model_xdg"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.path.join(os.environ["XDG_CACHE_HOME"], "fontconfig"), exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import leaves_list, linkage
from scipy.spatial.distance import pdist

from plot_ploidy_enrichment_panels import (
    HEATMAP_CMAP,
    infer_zero_replacement,
    ordered_by_low_ploidy_significance,
    p_to_neglog10,
    pretty_label,
    read_enrichment_workbook,
    star_color_for_value,
)

def finite_for_clustering(values):
    """Replace missing values with no-enrichment score for clustering."""
    values = np.asarray(values, dtype=float)
    return np.where(np.isfinite(values), values, 0.0)


def clustered_order(matrix, axis=0, method="ward", metric="euclidean"):
    """Return hierarchical-clustering leaf order for rows or columns."""
    values = finite_for_clustering(matrix)
    if axis == 1:
        values = values.T
    if values.shape[0] <= 1:
        return np.arange(values.shape[0])
    distances = pdist(values, metric=metric)
    if len(distances) == 0 or np.allclose(distances, 0):
        return np.arange(values.shape[0])
    return leaves_list(linkage(distances, method=method))


def display_label(label):
    """Readable drug-class label for compact clustered heatmaps."""
    label = pretty_label(label).replace(" AND ", " & ")
    compact_labels = {
        "Serine/threonine\nkinase inhibitors": "Ser/Thr kinase\ninhibitors",
        "Receptor tyrosine\nkinase inhibitors": "Receptor tyrosine\nkinase inhibitors",
        "Non-receptor tyrosine\nkinase inhibitors": "Non-receptor\ntyrosine kinase\ninhibitors",
        "Chaperone/protein-\nhomeostasis inhibitors": "Chaperone/protein\nhomeostasis\ninhibitors",
    }
    if label in compact_labels:
        return compact_labels[label]
    label = label.replace("Epigenetics & Transcription", "Epigenetics &\ntranscription")
    label = label.replace("Hormones & Antihormones", "Hormones &\nantihormones")
    label = label.replace("Tyrosine Kinase Inhibitors", "Tyrosine\nkinase\ninhibitors")
    return label


def ordered_labels(labels):
    return [display_label(label) for label in labels]


def read_optional_sheet(xlsx_path, sheet_name):
    """Read an optional workbook sheet, returning None if absent."""
    try:
        return pd.read_excel(xlsx_path, sheet_name=sheet_name)
    except ValueError:
        return None


def ordered_numeric_values(table, key_col, value_col, labels):
    """Return numeric values from a metadata table in heatmap label order."""
    if table is None or key_col not in table.columns or value_col not in table.columns:
        return None
    keys = table[key_col].astype(str)
    values = pd.to_numeric(table[value_col], errors="coerce")
    mapping = dict(zip(keys, values))
    ordered = np.array([mapping.get(str(label), np.nan) for label in labels], dtype=float)
    if np.all(~np.isfinite(ordered)):
        return None
    return ordered


def set_seaborn_x_labels(ax, labels, fontsize=7):
    """Place x labels at seaborn heatmap cell centers."""
    ax.set_xticks(np.arange(len(labels)) + 0.5)
    ax.set_xticklabels(labels, fontsize=fontsize)
    for label in ax.get_xticklabels():
        label.set_rotation(55)
        label.set_ha("right")
        label.set_va("top")
        label.set_rotation_mode("anchor")
    ax.tick_params(axis="x", which="major", pad=2)


def add_significance_stars(ax, qvalues, row_order, col_order, scores=None, vmin=0, vmax=None):
    """Mark cells with BH-FDR q <= 0.05."""
    ordered = qvalues.iloc[row_order, col_order].to_numpy(dtype=float)
    ordered_scores = None
    if scores is not None:
        ordered_scores = scores.iloc[row_order, col_order].to_numpy(dtype=float)
    if vmax is None and ordered_scores is not None:
        vmax = np.nanmax(ordered_scores)
    for row_idx in range(ordered.shape[0]):
        for col_idx in range(ordered.shape[1]):
            qvalue = ordered[row_idx, col_idx]
            if np.isfinite(qvalue) and qvalue <= 0.05:
                score = ordered_scores[row_idx, col_idx] if ordered_scores is not None else np.nan
                ax.text(
                    col_idx + 0.5,
                    row_idx + 0.5,
                    "*",
                    ha="center",
                    va="center",
                    color=star_color_for_value(score, vmin, vmax),
                    fontsize=9,
                    fontweight="bold",
                )


def save_shared_order_heatmap(
    low,
    high,
    low_scores,
    high_scores,
    out_prefix,
    drug_counts=None,
    cancer_type_counts=None,
):
    """Save two heatmaps using one shared row order and a low-to-high class order."""
    low_scores = pd.DataFrame(low_scores, index=low.index, columns=low.columns)
    high_scores = pd.DataFrame(high_scores, index=high.index, columns=high.columns)

    ordered_low, ordered_high = ordered_by_low_ploidy_significance(low, high)
    row_order = np.array([low.index.get_loc(label) for label in ordered_low.index])
    col_order = np.array([low.columns.get_loc(label) for label in ordered_low.columns])
    vmax = np.nanmax([np.nanmax(low_scores), np.nanmax(high_scores), -np.log10(0.05)])
    ordered_row_labels = list(ordered_low.index)
    ordered_col_labels = list(ordered_low.columns)
    ordered_drug_counts = ordered_numeric_values(
        drug_counts,
        key_col="category_label",
        value_col="n_drugs",
        labels=ordered_col_labels,
    )
    ordered_cell_line_counts = ordered_numeric_values(
        cancer_type_counts,
        key_col="cancer_type",
        value_col="n_cell_lines",
        labels=ordered_row_labels,
    )
    show_marginals = ordered_drug_counts is not None and ordered_cell_line_counts is not None

    if show_marginals:
        fig = plt.figure(figsize=(17.5, 10.0))
        grid = fig.add_gridspec(
            2,
            4,
            height_ratios=[0.22, 1.0],
            width_ratios=[0.26, 1.0, 1.0, 0.045],
            hspace=0.05,
            wspace=0.08,
        )
        top_axes = [fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[0, 2])]
        cell_count_ax = fig.add_subplot(grid[1, 0])
        axes = [fig.add_subplot(grid[1, 1]), fig.add_subplot(grid[1, 2]), fig.add_subplot(grid[1, 3])]

        n_cols = len(ordered_col_labels)
        y_max = max(1.0, np.nanmax(ordered_drug_counts) * 1.12)
        for idx, ax in enumerate(top_axes):
            ax.bar(
                np.arange(n_cols) + 0.5,
                ordered_drug_counts,
                width=0.80,
                color="#4C4C4C",
                edgecolor="none",
            )
            ax.set_xlim(0, n_cols)
            ax.set_ylim(0, y_max)
            ax.set_xticks([])
            ax.tick_params(axis="y", labelsize=7, length=2)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if idx == 0:
                ax.set_ylabel("Drugs\nper class", fontsize=8)
            else:
                ax.set_yticklabels([])
                ax.spines["left"].set_visible(False)

        n_rows = len(ordered_row_labels)
        finite_counts = np.where(np.isfinite(ordered_cell_line_counts), ordered_cell_line_counts, 0)
        cell_count_ax.barh(
            np.arange(n_rows) + 0.5,
            finite_counts,
            height=0.78,
            color="#6F6F6F",
            edgecolor="none",
        )
        cell_count_ax.set_ylim(n_rows, 0)
        cell_count_ax.set_yticks(np.arange(n_rows) + 0.5)
        cell_count_ax.set_yticklabels(
            [
                f"{label} (n={int(count)})" if np.isfinite(count) else f"{label} (n=NA)"
                for label, count in zip(ordered_row_labels, ordered_cell_line_counts)
            ],
            fontsize=8,
        )
        cell_count_ax.set_xlabel("Cell lines", fontsize=8)
        cell_count_ax.tick_params(axis="x", labelsize=7, length=2)
        cell_count_ax.tick_params(axis="y", length=0, pad=2)
        cell_count_ax.spines["top"].set_visible(False)
        cell_count_ax.spines["right"].set_visible(False)
    else:
        fig, axes = plt.subplots(
            1,
            3,
            figsize=(16, 8.5),
            gridspec_kw={"width_ratios": [1, 1, 0.04], "wspace": 0.08},
        )

    cmap = sns.color_palette(HEATMAP_CMAP, as_cmap=True)
    panels = [
        (
            axes[0],
            low_scores,
            low,
            "Low-ploidy-selective enrichment\n(ordered by chemotherapy-agent FDR significance)",
        ),
        (
            axes[1],
            high_scores,
            high,
            "High-ploidy-selective enrichment\n(same row and drug-class order)",
        ),
    ]

    heatmap = None
    for panel_idx, (ax, scores, pvalues, title) in enumerate(panels):
        ordered_scores = scores.iloc[row_order, col_order]
        heatmap = sns.heatmap(
            ordered_scores,
            ax=ax,
            cmap=cmap,
            vmin=0,
            vmax=vmax,
            cbar=False,
            linewidths=0.25,
            linecolor="#E8E8E8",
            xticklabels=False,
            yticklabels=False if show_marginals else ordered_scores.index,
        )
        set_seaborn_x_labels(ax, ordered_labels(ordered_scores.columns), fontsize=7)
        add_significance_stars(ax, pvalues, row_order, col_order, scores=scores, vmin=0, vmax=vmax)
        if show_marginals:
            top_axes[panel_idx].set_title(title, fontsize=12, pad=8)
        else:
            ax.set_title(title, fontsize=12, pad=10)
        ax.set_xlabel("Drug category")
        ax.tick_params(axis="y", labelsize=8)

    axes[1].set_ylabel("")
    axes[1].set_yticklabels([])
    cbar = fig.colorbar(heatmap.collections[0], cax=axes[2])
    cbar.set_label("-log10(BH-FDR q-value)")
    fig.suptitle(
        "GDSC ploidy-enrichment heatmaps ordered by low-ploidy chemotherapy FDR signal",
        fontsize=14,
    )
    if show_marginals:
        fig.subplots_adjust(left=0.14, right=0.93, top=0.90, bottom=0.28)
    else:
        fig.subplots_adjust(left=0.09, right=0.92, top=0.86, bottom=0.30)

    for suffix in (".png", ".pdf"):
        fig.savefig(f"{out_prefix}_shared_order{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_clustermap(scores, qvalues, title, out_prefix, vmax):
    """Save a single clustered heatmap with row and column dendrograms."""
    score_df = pd.DataFrame(scores, index=qvalues.index, columns=qvalues.columns).fillna(0.0)
    score_df.columns = ordered_labels(score_df.columns)

    grid = sns.clustermap(
        score_df,
        method="ward",
        metric="euclidean",
        cmap=HEATMAP_CMAP,
        vmin=0,
        vmax=vmax,
        linewidths=0.25,
        linecolor="#E8E8E8",
        figsize=(10.5, 9.5),
        cbar_kws={"label": "-log10(BH-FDR q-value)"},
        dendrogram_ratio=(0.16, 0.14),
        cbar_pos=(0.91, 0.34, 0.022, 0.30),
        xticklabels=False,
    )
    grid.fig.subplots_adjust(right=0.82, top=0.92, bottom=0.16)
    grid.ax_cbar.set_position((0.91, 0.34, 0.022, 0.30))
    grid.fig.suptitle(title, fontsize=13, y=1.02)
    grid.ax_heatmap.set_xlabel("Drug class")
    grid.ax_heatmap.set_ylabel("")
    grid.ax_heatmap.tick_params(axis="y", labelsize=8)

    row_order = grid.dendrogram_row.reordered_ind
    col_order = grid.dendrogram_col.reordered_ind
    set_seaborn_x_labels(
        grid.ax_heatmap,
        [score_df.columns[i] for i in col_order],
        fontsize=7,
    )
    add_significance_stars(grid.ax_heatmap, qvalues, row_order, col_order, scores=score_df, vmin=0, vmax=vmax)

    for suffix in (".png", ".pdf"):
        grid.fig.savefig(f"{out_prefix}{suffix}", dpi=300, bbox_inches="tight")
    plt.close(grid.fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx", help="Input drugsVsPloidyCorr.xlsx workbook")
    parser.add_argument(
        "--out-prefix",
        default="ploidy_enrichment_clustered",
        help="Output prefix for clustered heatmap PNG/PDF files",
    )
    args = parser.parse_args()

    xlsx_path = Path(args.xlsx)
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    low, high = read_enrichment_workbook(xlsx_path)
    drug_counts = read_optional_sheet(xlsx_path, "drugClassCounts")
    cancer_type_counts = read_optional_sheet(xlsx_path, "cancerTypeCounts")
    zero_replacement = infer_zero_replacement(low, high)
    low_scores = p_to_neglog10(low, zero_replacement)
    high_scores = p_to_neglog10(high, zero_replacement)
    vmax = np.nanmax([np.nanmax(low_scores), np.nanmax(high_scores), -np.log10(0.05)])

    save_shared_order_heatmap(
        low,
        high,
        low_scores,
        high_scores,
        out_prefix,
        drug_counts=drug_counts,
        cancer_type_counts=cancer_type_counts,
    )
    save_clustermap(
        low_scores,
        low,
        "Low-ploidy-selective enrichment by BH FDR, clustered rows and columns",
        f"{out_prefix}_lowpIsSens_clustermap",
        vmax,
    )
    save_clustermap(
        high_scores,
        high,
        "High-ploidy-selective enrichment by BH FDR, clustered rows and columns",
        f"{out_prefix}_highpIsSens_clustermap",
        vmax,
    )

    print("Wrote clustered heatmaps with prefix:", out_prefix)


if __name__ == "__main__":
    main()
