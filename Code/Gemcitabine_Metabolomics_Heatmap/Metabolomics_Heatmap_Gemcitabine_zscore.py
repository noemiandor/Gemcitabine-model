#!/usr/bin/env python3
"""
Generate the ordered row-z-score heatmap:
"Strong gemcitabine-response metabolites ordered by 2N/4N response class"

This script generates the heatmap of metabolomics results for gemcitabine response.

Input file requirements:
- A metabolite identifier column named "row identity (all IDs)" or, if absent, the first column.
- Replicate abundance columns named:
  2N_C1, 2N_C2, 2N_C3, 2N_C4
  2N_G1, 2N_G2, 2N_G3, 2N_G4
  4N_C1, 4N_C2, 4N_C3, 4N_C4
  4N_G1, 4N_G2, 4N_G3, 4N_G4

Methods:
- Zeros/missing values are imputed as half of the row-wise minimum positive value.
- Intensities are log2-transformed.
- Welch t-tests identify 2N and 4N gemcitabine responses.
- Two-way ANOVA estimates the ploidy × treatment interaction.
- Strong-response metabolites are selected with p < 0.05 and >=2-fold change,
  plus interaction-only hits where applicable.
- Heatmap values are row z-scores across all 16 individual samples, clipped at ±3.
"""

import argparse
from datetime import datetime
import os
import re
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ttest_ind
from statsmodels.formula.api import ols
import statsmodels.api as sm
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LinearSegmentedColormap


def load_and_preprocess(input_file: str):
    df_raw = pd.read_excel(input_file, sheet_name=0)

    id_candidates = ["row identity (all IDs)", "Metabolite ID", "Metabolite", "metabolite"]
    id_col = next((c for c in id_candidates if c in df_raw.columns), df_raw.columns[0])

    replicate_pattern = re.compile(r"^(2N|4N)_[CG][1-4]$")
    rep_cols = [c for c in df_raw.columns if replicate_pattern.match(str(c))]

    expected_cols = [f"{p}_{t}{i}" for p in ["2N", "4N"] for t in ["C", "G"] for i in range(1, 5)]
    ordered_cols = [c for c in expected_cols if c in rep_cols]

    required = set(expected_cols)
    missing = [c for c in expected_cols if c not in ordered_cols]
    if missing:
        raise ValueError(f"Missing expected replicate columns: {missing}")

    df = df_raw[[id_col] + ordered_cols].copy()
    df[id_col] = df[id_col].astype(str).str.strip()

    # Preserve duplicate metabolite names as separate measured features.
    counts = df[id_col].value_counts()
    df["Feature_ID"] = [
        f"{name} | row_{i}" if counts[name] > 1 else name
        for i, name in zip(df_raw.index, df[id_col])
    ]

    X_raw = df[ordered_cols].apply(pd.to_numeric, errors="coerce")

    # Drop features with no positive signal.
    keep = (X_raw.fillna(0) > 0).sum(axis=1) > 0
    df = df.loc[keep].reset_index(drop=True)
    X_raw = X_raw.loc[keep].reset_index(drop=True)

    # Impute zeros/missing as half of the row-wise minimum positive intensity.
    X_imp = X_raw.copy()
    row_min_positive = X_imp.where(X_imp > 0).min(axis=1)
    global_min_positive = X_imp.where(X_imp > 0).min().min()
    row_min_positive = row_min_positive.fillna(global_min_positive)

    for i in range(X_imp.shape[0]):
        fill_value = row_min_positive.iloc[i] / 2.0
        X_imp.iloc[i, :] = X_imp.iloc[i, :].replace(0, np.nan).fillna(fill_value)

    X_log2 = np.log2(X_imp)
    X_log2.index = df["Feature_ID"]

    meta = []
    for col in ordered_cols:
        ploidy, rest = col.split("_")
        treatment = "Gemcitabine" if rest.startswith("G") else "Control"
        group = f"{ploidy}_{'G' if treatment == 'Gemcitabine' else 'C'}"
        meta.append({
            "sample": col,
            "ploidy": ploidy,
            "treatment": treatment,
            "group": group,
            "replicate": rest[1:]
        })
    meta = pd.DataFrame(meta)

    group_cols = {
        "2N_C": [c for c in ordered_cols if c.startswith("2N_C")],
        "2N_G": [c for c in ordered_cols if c.startswith("2N_G")],
        "4N_C": [c for c in ordered_cols if c.startswith("4N_C")],
        "4N_G": [c for c in ordered_cols if c.startswith("4N_G")],
    }

    return df, X_log2, meta, group_cols, ordered_cols, id_col


def compute_statistics(df, X_log2, meta, group_cols, id_col):
    alpha = 0.05
    fc_cut = 1.0  # log2(2-fold)

    log2fc_2n = X_log2[group_cols["2N_G"]].mean(axis=1) - X_log2[group_cols["2N_C"]].mean(axis=1)
    log2fc_4n = X_log2[group_cols["4N_G"]].mean(axis=1) - X_log2[group_cols["4N_C"]].mean(axis=1)

    p_2n = ttest_ind(
        X_log2[group_cols["2N_G"]].values,
        X_log2[group_cols["2N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit"
    ).pvalue

    p_4n = ttest_ind(
        X_log2[group_cols["4N_G"]].values,
        X_log2[group_cols["4N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit"
    ).pvalue

    p_interaction = []
    for idx in range(X_log2.shape[0]):
        dat = meta.copy()
        dat["value"] = X_log2.iloc[idx].values
        try:
            model = ols("value ~ C(ploidy) + C(treatment) + C(ploidy):C(treatment)", data=dat).fit()
            aov = sm.stats.anova_lm(model, typ=2)
            p_interaction.append(aov.loc["C(ploidy):C(treatment)", "PR(>F)"])
        except Exception:
            p_interaction.append(np.nan)

    p_interaction = np.array(p_interaction, dtype=float)
    diff_response = log2fc_2n.values - log2fc_4n.values

    results = pd.DataFrame({
        "Feature_ID": df["Feature_ID"],
        "Metabolite": df[id_col],
        "log2FC_2N_G_vs_2N_C": log2fc_2n.values,
        "p_2N_G_vs_2N_C": p_2n,
        "log2FC_4N_G_vs_4N_C": log2fc_4n.values,
        "p_4N_G_vs_4N_C": p_4n,
        "Differential_response_2N_minus_4N": diff_response,
        "p_interaction": p_interaction,
    })

    order_category = []
    response_category = []

    for _, r in results.iterrows():
        s2 = (r["p_2N_G_vs_2N_C"] < alpha) and (abs(r["log2FC_2N_G_vs_2N_C"]) >= fc_cut)
        s4 = (r["p_4N_G_vs_4N_C"] < alpha) and (abs(r["log2FC_4N_G_vs_4N_C"]) >= fc_cut)
        inter = (r["p_interaction"] < alpha) and (abs(r["Differential_response_2N_minus_4N"]) >= fc_cut)

        if s2 and s4:
            if np.sign(r["log2FC_2N_G_vs_2N_C"]) == np.sign(r["log2FC_4N_G_vs_4N_C"]):
                order_category.append("Common up" if r["log2FC_2N_G_vs_2N_C"] > 0 else "Common down")
                response_category.append("Common gem response")
            else:
                order_category.append("Opposite")
                response_category.append("Opposite gem response")
        elif s2:
            order_category.append("2N up" if r["log2FC_2N_G_vs_2N_C"] > 0 else "2N down")
            response_category.append("2N-specific response")
        elif s4:
            order_category.append("4N up" if r["log2FC_4N_G_vs_4N_C"] > 0 else "4N down")
            response_category.append("4N-specific response")
        elif inter:
            order_category.append("Interaction only")
            response_category.append("Differential response only")
        else:
            order_category.append("No strong response")
            response_category.append("No strong response")

    category_order = [
        "2N up", "2N down", "4N up", "4N down",
        "Common up", "Common down", "Opposite", "Interaction only"
    ]
    order_map = {c: i for i, c in enumerate(category_order)}

    results["Order_category"] = order_category
    results["Response_category"] = response_category
    results["cat_order"] = results["Order_category"].map(order_map).fillna(99).astype(int)

    return results


def make_ordered_zscore_heatmap(X_log2, ordered_cols, results, output_png, output_pdf=None, matrix_csv=None, metabolites_csv=None):
    strong = results[results["Order_category"] != "No strong response"].copy()

    # Match the earlier heatmap ordering: response class first, then interaction p-value, then Feature_ID.
    strong = strong.sort_values(["cat_order", "p_interaction", "Feature_ID"], ascending=[True, True, True])

    z_scores = pd.DataFrame(
        StandardScaler(with_mean=True, with_std=True).fit_transform(X_log2.T).T,
        index=X_log2.index,
        columns=ordered_cols
    ).clip(-3, 3)

    heatmap_matrix = z_scores.loc[strong["Feature_ID"], ordered_cols].copy()

    labels = []
    seen = {}
    for metabolite in strong["Metabolite"].astype(str):
        seen[metabolite] = seen.get(metabolite, 0) + 1
        labels.append(metabolite if seen[metabolite] == 1 else f"{metabolite} ({seen[metabolite]})")
    heatmap_matrix.index = labels

    if matrix_csv:
        heatmap_matrix.to_csv(matrix_csv)
    if metabolites_csv:
        strong.to_csv(metabolites_csv, index=False)

    custom_bwr = LinearSegmentedColormap.from_list("custom_bwr", ["#001a4d", "#ffffff", "#800000"])

    plt.figure(figsize=(12, max(9, 0.22 * len(heatmap_matrix))))
    ax = sns.heatmap(
        heatmap_matrix,
        cmap=custom_bwr,
        vmin=-3,
        vmax=3,
        center=0,
        linewidths=0.25,
        annot=False,
        cbar_kws={"label": "Row z-score (clipped at ±3)"}
    )
    ax.set_title("Strong gemcitabine-response metabolites ordered by 2N/4N response class", fontsize=14)
    ax.set_xlabel("Samples", fontsize=11)
    ax.set_ylabel("Metabolites", fontsize=11)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=7)

    plt.tight_layout()
    plt.savefig(output_png, dpi=300, bbox_inches="tight")
    if output_pdf:
        plt.savefig(output_pdf, bbox_inches="tight")
    plt.close()

    return strong, heatmap_matrix


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input .xlsm/.xlsx metabolomics file")
    parser.add_argument("--outdir", help="Output directory")
    parser.add_argument("--output-dir", help="Canonical output directory; writes plots to figures/ and tables to tables/.")
    args = parser.parse_args()

    canonical_output = args.output_dir is not None or args.outdir is None
    if canonical_output:
        if args.output_dir:
            args.outdir = args.output_dir
        else:
            repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            args.outdir = os.path.join(
                repo_root,
                "Results",
                "in-vitro",
                "metabolomics",
                "runs",
                f"{datetime.now().strftime('%Y%m%dT%H%M%S')}_metabolomics_zscore",
            )

    os.makedirs(args.outdir, exist_ok=True)
    figdir = os.path.join(args.outdir, "figures") if canonical_output else args.outdir
    tabledir = os.path.join(args.outdir, "tables") if canonical_output else args.outdir
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabledir, exist_ok=True)

    df, X_log2, meta, group_cols, ordered_cols, id_col = load_and_preprocess(args.input)
    results = compute_statistics(df, X_log2, meta, group_cols, id_col)

    output_png = os.path.join(figdir, "ordered_response_heatmap_reproduced.png")
    output_pdf = os.path.join(figdir, "ordered_response_heatmap_reproduced.pdf")
    matrix_csv = os.path.join(tabledir, "ordered_response_heatmap_zscore_matrix.csv")
    metabolites_csv = os.path.join(tabledir, "ordered_response_heatmap_metabolites.csv")
    results_csv = os.path.join(tabledir, "ordered_response_heatmap_statistics.csv")

    strong, heatmap_matrix = make_ordered_zscore_heatmap(
        X_log2,
        ordered_cols,
        results,
        output_png,
        output_pdf=output_pdf,
        matrix_csv=matrix_csv,
        metabolites_csv=metabolites_csv,
    )

    results.to_csv(results_csv, index=False)

    print(f"Loaded {len(df)} metabolite features.")
    print(f"Selected {len(strong)} strong-response metabolites for the heatmap.")
    print("Saved:")
    print(" -", output_png)
    print(" -", output_pdf)
    print(" -", matrix_csv)
    print(" -", metabolites_csv)
    print(" -", results_csv)


if __name__ == "__main__":
    main()
