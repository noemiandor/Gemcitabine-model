
#!/usr/bin/env python3
"""
Create manuscript-style Panels A-C from drugsVsPloidyCorr.xlsx.

Input workbook requirements:
- Sheet "lowpIsSens": enrichment p-values for drug classes enriched among low-ploidy-selective drugs.
- Sheet "highpIsSens": enrichment p-values for drug classes enriched among high-ploidy-selective drugs.
- Rows are cancer types.
- Columns are drug classes.
- Cell values are enrichment p-values.

Outputs:
- Combined figure with:
  A. Workflow schematic
  B. Heatmap of low-ploidy-selective enrichment
  C. Heatmap of high-ploidy-selective enrichment
- PNG and PDF versions.
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D


def read_enrichment_workbook(xlsx_path, low_sheet="lowpIsSens", high_sheet="highpIsSens"):
    """Read the two expected enrichment-p-value sheets."""
    low = pd.read_excel(xlsx_path, sheet_name=low_sheet, index_col=0)
    high = pd.read_excel(xlsx_path, sheet_name=high_sheet, index_col=0)

    # Ensure all values are numeric p-values.
    low = low.apply(pd.to_numeric, errors="coerce")
    high = high.apply(pd.to_numeric, errors="coerce")

    # Align row and column order so both heatmaps are directly comparable.
    row_order = list(dict.fromkeys(list(low.index) + list(high.index)))
    col_order = list(dict.fromkeys(list(low.columns) + list(high.columns)))
    low = low.reindex(index=row_order, columns=col_order)
    high = high.reindex(index=row_order, columns=col_order)

    return low, high


def infer_zero_replacement(*dfs):
    """
    Replace p=0 values for visualization.

    Permutation-derived p-values can be exactly zero when no permutation is
    as extreme as the observed statistic. For plotting -log10(p), replace
    zeros with the smallest positive p-value observed in the workbook.
    """
    vals = []
    for df in dfs:
        arr = df.to_numpy(dtype=float).ravel()
        vals.extend(arr[np.isfinite(arr) & (arr > 0)])
    if not vals:
        return 1e-6
    return float(np.min(vals))


def p_to_neglog10(df, zero_replacement):
    """Convert p-values to -log10(p), replacing zeros and preserving NaNs."""
    arr = df.to_numpy(dtype=float)
    arr = np.where(arr == 0, zero_replacement, arr)
    arr = np.where(arr > 1, np.nan, arr)
    arr = np.where(arr < 0, np.nan, arr)
    with np.errstate(divide="ignore", invalid="ignore"):
        out = -np.log10(arr)
    out[~np.isfinite(out)] = np.nan
    return out


def pretty_label(x):
    """Make drug-class labels more readable while preserving common acronyms."""
    x = str(x).replace(".", " ").replace("_", " ")
    x = " ".join(x.split())
    acronyms = {"MEK", "DNA", "RNA", "GDSC"}
    words = []
    for w in x.split():
        if w.upper() in acronyms:
            words.append(w.upper())
        elif len(w) <= 3 and w.isupper():
            words.append(w)
        else:
            words.append(w.capitalize())

    label = " ".join(words)

    # Manual wrapping for compact heatmap labels.
    label = label.replace("Tyrosine Kinase Inhibitors", "Tyrosine\nkinase\ninhibitors")
    label = label.replace("Antineoplastic Agents", "Antineoplastic\nagents")
    label = label.replace("Immunosuppressive Agents", "Immunosuppressive\nagents")
    label = label.replace("Epigenetics And Transcription", "Epigenetics &\ntranscription")
    label = label.replace("Hormones And Antihormones", "Hormones &\nantihormones")
    label = label.replace("MEK Inhibitors", "MEK\ninhibitors")
    label = label.replace("Metabolicinhibitor", "Metabolic\ninhibitor")
    return label


def set_imshow_x_labels(ax, labels, fontsize=8):
    """Place x labels at imshow cell centers and anchor rotated text to ticks."""
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, fontsize=fontsize)
    for label in ax.get_xticklabels():
        label.set_rotation(35)
        label.set_ha("right")
        label.set_va("top")
        label.set_rotation_mode("anchor")
    ax.tick_params(axis="x", which="major", pad=3)


def add_workflow_panel(ax):
    """Draw Panel A: workflow schematic."""
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    boxes = [
        (0.03, 0.35, 0.24, 0.36, "GDSC drug\nresponse"),
        (0.285, 0.35, 0.20, 0.36, "Cell-line\nploidy"),
        (0.535, 0.35, 0.20, 0.36, "Drug-level\ncorrelations"),
        (0.785, 0.35, 0.19, 0.36, "Drug-class\nenrichment"),
    ]

    for x, y, w, h, label in boxes:
        rect = patches.FancyBboxPatch(
            (x, y), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.03",
            linewidth=1.2,
            facecolor="white",
            edgecolor="black"
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=11)

    arrow_y = 0.53
    arrowprops = dict(arrowstyle="->", linewidth=1.2, shrinkA=0, shrinkB=0)
    ax.annotate("", xy=(0.285, arrow_y), xytext=(0.27, arrow_y), arrowprops=arrowprops)
    ax.annotate("", xy=(0.535, arrow_y), xytext=(0.485, arrow_y), arrowprops=arrowprops)
    ax.annotate("", xy=(0.785, arrow_y), xytext=(0.735, arrow_y), arrowprops=arrowprops)

    ax.text(0.41, 0.20, "For each cancer type and drug: correlate response with ploidy",
            ha="center", va="center", fontsize=9)
    ax.text(0.76, 0.20, "Then test whether drug classes are enriched by direction",
            ha="center", va="center", fontsize=9)

    ax.text(0.0, 0.98, "A", fontsize=16, fontweight="bold", ha="left", va="top")


def add_heatmap_panel(ax, data, pvalues, title, panel_letter, vmin, vmax):
    """Draw one heatmap panel."""
    masked = np.ma.masked_invalid(data)
    cmap = plt.get_cmap().copy()
    cmap.set_bad("#F2F2F2")

    im = ax.imshow(masked, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)

    ax.set_title(title, fontsize=12, pad=12)
    set_imshow_x_labels(ax, [pretty_label(c) for c in pvalues.columns], fontsize=8)
    ax.set_yticks(np.arange(pvalues.shape[0]))
    ax.set_yticklabels([str(i) for i in pvalues.index], fontsize=8)

    # Light gridlines around cells.
    ax.set_xticks(np.arange(-0.5, pvalues.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, pvalues.shape[0], 1), minor=True)
    ax.grid(which="minor", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    # Mark nominally significant cells to make interpretation immediate.
    arr = pvalues.to_numpy(dtype=float)
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            p = arr[i, j]
            if np.isfinite(p) and p <= 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=9, fontweight="bold")

    ax.text(-0.12, 1.04, panel_letter, transform=ax.transAxes,
            fontsize=16, fontweight="bold", ha="left", va="bottom")

    return im


def make_figure(xlsx_path, out_png, out_pdf=None):
    low, high = read_enrichment_workbook(xlsx_path)
    zero_replacement = infer_zero_replacement(low, high)

    low_z = p_to_neglog10(low, zero_replacement)
    high_z = p_to_neglog10(high, zero_replacement)

    vmax = np.nanmax([np.nanmax(low_z), np.nanmax(high_z), -np.log10(0.05)])
    vmin = 0

    fig = plt.figure(figsize=(18, 14), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=2,
        height_ratios=[1.15, 6.0],
        width_ratios=[1, 1],
        hspace=0.34,
        wspace=0.25
    )

    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    add_workflow_panel(ax_a)

    im_b = add_heatmap_panel(
        ax_b, low_z, low,
        "Low-ploidy-selective enrichment",
        "B", vmin, vmax
    )
    im_c = add_heatmap_panel(
        ax_c, high_z, high,
        "High-ploidy-selective enrichment",
        "C", vmin, vmax
    )

    # Leave generous margins for rotated x labels and an external colorbar.
    fig.subplots_adjust(left=0.07, right=0.885, top=0.955, bottom=0.185)

    cbar_ax = fig.add_axes([0.91, 0.31, 0.014, 0.40])
    cbar = fig.colorbar(im_c, cax=cbar_ax)
    cbar.set_label("-log10(enrichment p-value)", fontsize=10)

    # Add concise explanatory footnote.
    fig.text(
        0.5, 0.035,
        f"Stars mark nominal enrichment p ≤ 0.05. Exact p=0 entries were plotted at the permutation floor "
        f"(smallest positive p = {zero_replacement:g}) to avoid infinite -log10 values.",
        ha="center", va="bottom", fontsize=9
    )

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    if out_pdf is not None:
        fig.savefig(out_pdf, bbox_inches="tight")

    return {
        "low_shape": low.shape,
        "high_shape": high.shape,
        "zero_replacement": zero_replacement,
        "vmax": vmax,
        "out_png": str(out_png),
        "out_pdf": str(out_pdf) if out_pdf is not None else None
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("xlsx", help="Input enrichment p-value workbook")
    parser.add_argument("--out-prefix", default="ploidy_enrichment_panels_ABC",
                        help="Output prefix for PNG/PDF")
    args = parser.parse_args()

    xlsx_path = Path(args.xlsx)
    out_prefix = Path(args.out_prefix)
    out_png = out_prefix.with_suffix(".png")
    out_pdf = out_prefix.with_suffix(".pdf")

    summary = make_figure(xlsx_path, out_png, out_pdf)
    print(summary)


if __name__ == "__main__":
    main()
