
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
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "gemcitabine_model_mpl"))
os.environ.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "gemcitabine_model_xdg"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.path.join(os.environ["XDG_CACHE_HOME"], "fontconfig"), exist_ok=True)

import matplotlib

matplotlib.use("Agg")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.lines import Line2D

HEATMAP_CMAP = "viridis"
CHEMOTHERAPY_CATEGORY = "Chemotherapy agents"
CHEMOTHERAPY_PRIMARY_CLASSES = {
    "Alkylating agents",
    "Antimetabolites",
    "Antimitotic agents",
    "Topoisomerase inhibitors",
    "Tumor antibiotics",
}
ANTIMITOTIC_AGENTS_LABEL = "Antimitotic agents"
TOP_ROW_LABELS = ("allcancers",)
ROW_AFTER_LABELS = {"unclassified": "stad"}
ANTIMITOTIC_PRIORITY_ROWS = ("BRCA", "ALL")
ANTIMETABOLITES_LABEL = "antimetabolites"
OTHER_LABEL = "other"
NON_RECEPTOR_TK_LABEL = "non-receptor tyrosine kinase inhibitors"


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


def chemotherapy_columns(columns):
    """Find the aggregate chemotherapy column or its primary-class members."""
    columns = list(columns)
    if CHEMOTHERAPY_CATEGORY in columns:
        return [CHEMOTHERAPY_CATEGORY]
    return [col for col in columns if col in CHEMOTHERAPY_PRIMARY_CLASSES]


def norm_label(label):
    """Normalize display labels for explicit ordering rules."""
    return str(label).casefold()


def place_special_rows(row_order):
    """Pin allcancers first and UNCLASSIFIED immediately below STAD."""
    rows = list(row_order)
    top_norms = {norm_label(label) for label in TOP_ROW_LABELS}
    after_norms = set(ROW_AFTER_LABELS)
    top_rows = [row for row in rows if norm_label(row) in top_norms]
    moved_after_rows = [row for row in rows if norm_label(row) in after_norms]
    remaining = [
        row for row in rows
        if norm_label(row) not in top_norms and norm_label(row) not in after_norms
    ]

    ordered = list(top_rows)
    inserted_after_rows = set()
    for row in remaining:
        ordered.append(row)
        anchor_norm = norm_label(row)
        for moved_norm, after_norm in ROW_AFTER_LABELS.items():
            if after_norm == anchor_norm:
                ordered.extend(
                    moved_row for moved_row in moved_after_rows
                    if norm_label(moved_row) == moved_norm
                )
                inserted_after_rows.add(moved_norm)

    ordered.extend(
        moved_row for moved_row in moved_after_rows
        if norm_label(moved_row) not in inserted_after_rows
    )
    return pd.Index(ordered)


def place_antimitotic_rows(row_order):
    """Pin allcancers first, then BRCA above ALL, before antimitotic-ranked rows."""
    rows = list(row_order)
    top_norms = {norm_label(label) for label in TOP_ROW_LABELS}
    priority_norms = {norm_label(label) for label in ANTIMITOTIC_PRIORITY_ROWS}
    top_rows = [row for row in rows if norm_label(row) in top_norms]
    priority_rows = [
        row for label in ANTIMITOTIC_PRIORITY_ROWS
        for row in rows
        if norm_label(row) == norm_label(label)
    ]
    remaining = [
        row for row in rows
        if norm_label(row) not in top_norms and norm_label(row) not in priority_norms
    ]
    return pd.Index(top_rows + priority_rows + remaining)


def is_epigenetic_drug_category(column):
    """Flag epigenetic categories requested next to the right edge."""
    return "epigenetic" in norm_label(column)


def is_kinase_drug_category(column):
    """Flag kinase categories requested as the rightmost columns."""
    return "kinase" in norm_label(column)


def is_non_receptor_tyrosine_kinase_category(column):
    """Flag the kinase subclass requested before epigenetic inhibitors."""
    return norm_label(column) == NON_RECEPTOR_TK_LABEL


def move_columns_before(columns, moving_labels, anchor_label):
    """Move named columns immediately before an anchor when both are present."""
    moving_norms = {norm_label(label) for label in moving_labels}
    anchor_norm = norm_label(anchor_label)
    moving_cols = [col for col in columns if norm_label(col) in moving_norms]
    if not moving_cols:
        return list(columns)

    remaining = [col for col in columns if norm_label(col) not in moving_norms]
    ordered = []
    inserted = False
    for col in remaining:
        if norm_label(col) == anchor_norm and not inserted:
            ordered.extend(moving_cols)
            inserted = True
        ordered.append(col)
    if not inserted:
        ordered.extend(moving_cols)
    return ordered


def place_special_columns(col_order):
    """Apply requested drug-category placements after significance sorting."""
    cols = list(col_order)
    left_cols = [
        col for col in cols
        if (
            not is_epigenetic_drug_category(col)
            and not is_kinase_drug_category(col)
            and not is_non_receptor_tyrosine_kinase_category(col)
        )
    ]
    left_cols = move_columns_before(left_cols, [ANTIMETABOLITES_LABEL], OTHER_LABEL)
    non_receptor_tk_cols = [
        col for col in cols
        if is_non_receptor_tyrosine_kinase_category(col)
    ]
    epigenetic_cols = [col for col in cols if is_epigenetic_drug_category(col)]
    kinase_cols = [
        col for col in cols
        if is_kinase_drug_category(col) and not is_non_receptor_tyrosine_kinase_category(col)
    ]
    return pd.Index(left_cols + non_receptor_tk_cols + epigenetic_cols + kinase_cols)


def ordered_by_low_ploidy_significance(low, high):
    """
    Order both heatmaps from the low-ploidy panel.

    Rows are sorted by antimitotic-agent enrichment significance when that
    primary class is present; otherwise by chemotherapy-agent enrichment
    significance in the low-ploidy-selective heatmap. Columns are sorted by their strongest
    low-ploidy-selective enrichment, with epigenetic categories placed near
    the right edge and kinase categories after them. The same row/column order
    is then reused for the high-ploidy-selective heatmap.
    """
    if ANTIMITOTIC_AGENTS_LABEL in low.columns:
        row_pvalue = low[ANTIMITOTIC_AGENTS_LABEL]
        row_signal = p_to_neglog10(
            low[[ANTIMITOTIC_AGENTS_LABEL]],
            infer_zero_replacement(low[[ANTIMITOTIC_AGENTS_LABEL]]),
        ).sum(axis=1)
        row_placer = place_antimitotic_rows
    else:
        chemo_cols = chemotherapy_columns(low.columns)
        if chemo_cols:
            row_pvalue = low[chemo_cols].min(axis=1, skipna=True)
        else:
            row_pvalue = low.min(axis=1, skipna=True)
        row_signal = p_to_neglog10(low, infer_zero_replacement(low)).sum(axis=1)
        row_placer = place_special_rows
    row_order = (
        pd.DataFrame(
            {
                "pvalue": row_pvalue,
                "signal": row_signal,
                "label": low.index.astype(str),
            },
            index=low.index,
        )
        .assign(pvalue=lambda df: df["pvalue"].fillna(np.inf))
        .sort_values(["pvalue", "signal", "label"], ascending=[True, False, True], kind="mergesort")
        .index
    )
    row_order = row_placer(row_order)

    col_pvalue = low.min(axis=0, skipna=True)
    col_signal = p_to_neglog10(low, infer_zero_replacement(low)).sum(axis=0)
    col_order = (
        pd.DataFrame(
            {
                "pvalue": col_pvalue,
                "signal": col_signal,
                "label": low.columns.astype(str),
            },
            index=low.columns,
        )
        .assign(pvalue=lambda df: df["pvalue"].fillna(np.inf))
        .sort_values(["pvalue", "signal", "label"], ascending=[True, False, True], kind="mergesort")
        .index
    )
    col_order = place_special_columns(col_order)

    return low.loc[row_order, col_order], high.loc[row_order, col_order]


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


def star_color_for_value(value, vmin, vmax, cmap_name=HEATMAP_CMAP):
    """Choose black or white significance stars based on heatmap-cell luminance."""
    if not np.isfinite(value) or not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        return "black"
    scaled = np.clip((value - vmin) / (vmax - vmin), 0.0, 1.0)
    r, g, b, _ = plt.get_cmap(cmap_name)(scaled)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "black" if luminance >= 0.55 else "white"


PRIMARY_CLASS_LABELS = {
    "Serine/threonine kinase inhibitors": "Serine/threonine\nkinase inhibitors",
    "Epigenetic inhibitors": "Epigenetic\ninhibitors",
    "Other": "Other",
    "Receptor tyrosine kinase inhibitors": "Receptor tyrosine\nkinase inhibitors",
    "DNA damage repair inhibitors": "DNA damage\nrepair inhibitors",
    "Metabolic/redox agents": "Metabolic/redox\nagents",
    "Proapoptotic agents": "Proapoptotic\nagents",
    "Non-receptor tyrosine kinase inhibitors": "Non-receptor tyrosine\nkinase inhibitors",
    "WNT-pathway modulators": "WNT-pathway\nmodulators",
    "Antimetabolites": "Antimetabolites",
    "Topoisomerase inhibitors": "Topoisomerase\ninhibitors",
    "Chaperone/protein-homeostasis inhibitors": "Chaperone/protein-\nhomeostasis inhibitors",
    "Alkylating agents": "Alkylating\nagents",
    "Antimitotic agents": "Antimitotic\nagents",
    "p53/MDM2 pathway": "p53/MDM2\npathway",
    "Hormone therapy": "Hormone\ntherapy",
    "Tumor antibiotics": "Tumor\nantibiotics",
    "Proteasome inhibitors": "Proteasome\ninhibitors",
    "Hedgehog pathway inhibitors": "Hedgehog pathway\ninhibitors",
}


def pretty_label(x):
    """Make drug-class labels readable without changing biological labels."""
    label = " ".join(str(x).replace(".", " ").replace("_", " ").split())
    if label in PRIMARY_CLASS_LABELS:
        return PRIMARY_CLASS_LABELS[label]

    replacements = {
        "Wnt": "WNT",
        "Dna": "DNA",
        "Rna": "RNA",
        "Gdsc": "GDSC",
        "P53/mdm2": "p53/MDM2",
        "Mdm2": "MDM2",
    }
    for old, new in replacements.items():
        label = label.replace(old, new)
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


def add_heatmap_panel(ax, data, pvalues, title, panel_letter, vmin, vmax, show_ylabels=True):
    """Draw one heatmap panel."""
    masked = np.ma.masked_invalid(data)
    cmap = plt.get_cmap(HEATMAP_CMAP).copy()
    cmap.set_bad("#F2F2F2")

    im = ax.imshow(masked, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)

    ax.set_title(title, fontsize=12, pad=12)
    set_imshow_x_labels(ax, [pretty_label(c) for c in pvalues.columns], fontsize=8)
    ax.set_yticks(np.arange(pvalues.shape[0]))
    if show_ylabels:
        ax.set_yticklabels([str(i) for i in pvalues.index], fontsize=8)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", labelleft=False)

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
                ax.text(
                    j,
                    i,
                    "*",
                    ha="center",
                    va="center",
                    color=star_color_for_value(data[i, j], vmin, vmax),
                    fontsize=9,
                    fontweight="bold",
                )

    ax.text(-0.12, 1.04, panel_letter, transform=ax.transAxes,
            fontsize=16, fontweight="bold", ha="left", va="bottom")

    return im


def make_figure(xlsx_path, out_png, out_pdf=None):
    low, high = read_enrichment_workbook(xlsx_path)
    low, high = ordered_by_low_ploidy_significance(low, high)
    zero_replacement = infer_zero_replacement(low, high)

    low_z = p_to_neglog10(low, zero_replacement)
    high_z = p_to_neglog10(high, zero_replacement)

    vmax = np.nanmax([np.nanmax(low_z), np.nanmax(high_z), -np.log10(0.05)])
    vmin = 0

    fig = plt.figure(figsize=(18, 14), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=2,
        ncols=3,
        height_ratios=[1.15, 6.0],
        width_ratios=[1, 1, 0.04],
        hspace=0.34,
        wspace=0.08
    )

    ax_a = fig.add_subplot(gs[0, :2])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_cbar = fig.add_subplot(gs[1, 2])

    add_workflow_panel(ax_a)

    im_b = add_heatmap_panel(
        ax_b, low_z, low,
        "Low-ploidy-selective enrichment\n(ordered by chemotherapy-agent significance)",
        "B", vmin, vmax
    )
    im_c = add_heatmap_panel(
        ax_c, high_z, high,
        "High-ploidy-selective enrichment\n(same row and drug-class order)",
        "C", vmin, vmax, show_ylabels=False
    )

    # Leave generous margins for rotated x labels and align the colorbar with the heatmap row.
    fig.subplots_adjust(left=0.07, right=0.925, top=0.955, bottom=0.185)

    cbar = fig.colorbar(im_c, cax=ax_cbar)
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
