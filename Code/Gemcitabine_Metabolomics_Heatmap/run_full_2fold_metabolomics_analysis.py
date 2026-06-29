#!/usr/bin/env python3
"""
Full metabolomics analysis for 2N/4N control and gemcitabine-treated samples.

This script performs the 2-fold / p<0.05 version of the workflow:
- Reads replicate columns: 2N_C1-C4, 2N_G1-G4, 4N_C1-C4, 4N_G1-G4
- Imputes zero/missing values as half of the row-wise minimum positive intensity
- Log2 transforms metabolite intensities
- Runs pairwise Welch t-tests
- Runs two-way ANOVA: log2 abundance ~ ploidy + treatment + ploidy:treatment
- Runs one-way ANOVA and Tukey HSD for top hits
- Classifies response categories using p<0.05 and |log2FC|>=1
- Creates PCA figures with and without sample labels
- Creates volcano plots with black insignificant points and red significant points,
  with and without metabolite labels
- Creates z-score and replicate-level fold-change heatmaps
- Creates pathway-class enrichment tables and figure

Run:
    python run_full_2fold_metabolomics_analysis.py \
        --input /path/to/Metabolomics_2N_4N_Full.xlsm \
        --out /path/to/output_folder
"""

import argparse
import os
import re
import zipfile
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind, f_oneway, hypergeom
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.patches import Patch


ALPHA = 0.05
FC_CUT = np.log2(2.0)
HEATMAP_CLIP = 3.0


def safe_neglog10(p):
    return -np.log10(np.clip(np.asarray(p, dtype=float), 1e-300, None))


def truncate_label(label, max_len=80):
    label = str(label)
    return label if len(label) <= max_len else label[: max_len - 3] + "..."


def unique_labels(labels, max_len=80):
    out, seen = [], {}
    for lab in labels:
        lab = truncate_label(lab, max_len=max_len)
        seen[lab] = seen.get(lab, 0) + 1
        out.append(lab if seen[lab] == 1 else f"{lab} ({seen[lab]})")
    return out


def savefig(path, dpi=300):
    plt.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close()


def load_and_preprocess(input_file):
    df_raw = pd.read_excel(input_file, sheet_name=0)
    id_candidates = ["row identity (all IDs)", "Metabolite ID", "Metabolite", "metabolite"]
    id_col = next((c for c in id_candidates if c in df_raw.columns), df_raw.columns[0])

    replicate_pattern = re.compile(r"^(2N|4N)_[CG][1-4]$")
    rep_cols = [c for c in df_raw.columns if replicate_pattern.match(str(c))]
    expected_cols = [f"{p}_{t}{i}" for p in ["2N", "4N"] for t in ["C", "G"] for i in range(1, 5)]
    ordered_cols = [c for c in expected_cols if c in rep_cols]

    df = df_raw[[id_col] + ordered_cols].copy()
    df[id_col] = df[id_col].astype(str).str.strip()

    counts = df[id_col].value_counts()
    df["Feature_ID"] = [
        f"{name} | row_{i}" if counts[name] > 1 else name
        for i, name in zip(df_raw.index, df[id_col])
    ]

    X_raw = df[ordered_cols].apply(pd.to_numeric, errors="coerce")
    keep = (X_raw.fillna(0) > 0).sum(axis=1) > 0
    df = df.loc[keep].reset_index(drop=True)
    X_raw = X_raw.loc[keep].reset_index(drop=True)

    X_imp = X_raw.copy()
    row_min_positive = X_imp.where(X_imp > 0).min(axis=1)
    global_min_positive = X_imp.where(X_imp > 0).min().min()
    row_min_positive = row_min_positive.fillna(global_min_positive)
    for i in range(X_imp.shape[0]):
        fill_value = row_min_positive.iloc[i] / 2.0
        X_imp.iloc[i, :] = X_imp.iloc[i, :].replace(0, np.nan).fillna(fill_value)

    X_log2 = np.log2(X_imp)
    X_log2.index = df["Feature_ID"]

    meta_rows = []
    for col in ordered_cols:
        ploidy, rest = col.split("_")
        treatment = "Gemcitabine" if rest.startswith("G") else "Control"
        group = f"{ploidy}_{'G' if treatment == 'Gemcitabine' else 'C'}"
        meta_rows.append({"sample": col, "ploidy": ploidy, "treatment": treatment, "group": group, "replicate": rest[1:]})
    meta = pd.DataFrame(meta_rows)

    group_cols = {
        "2N_C": [c for c in ordered_cols if c.startswith("2N_C")],
        "2N_G": [c for c in ordered_cols if c.startswith("2N_G")],
        "4N_C": [c for c in ordered_cols if c.startswith("4N_C")],
        "4N_G": [c for c in ordered_cols if c.startswith("4N_G")],
    }
    return df, X_log2, meta, group_cols, ordered_cols, id_col


def compute_statistics(df, X_log2, meta, group_cols, id_col):
    def pairwise_stats(group_a, group_b, label):
        A = X_log2[group_cols[group_a]]
        B = X_log2[group_cols[group_b]]
        log2fc = B.mean(axis=1) - A.mean(axis=1)
        pvals = ttest_ind(B.values, A.values, axis=1, equal_var=False, nan_policy="omit").pvalue
        fdr = multipletests(np.nan_to_num(pvals, nan=1.0), method="fdr_bh")[1]
        return pd.DataFrame({
            "Feature_ID": df["Feature_ID"],
            "Metabolite": df[id_col],
            "comparison": label,
            "group_a": group_a,
            "group_b": group_b,
            "log2FC": log2fc.values,
            "fold_change": (2 ** log2fc).values,
            "p_value": pvals,
            "FDR_BH": fdr,
        })

    pairwise_results = pd.concat([
        pairwise_stats("2N_C", "2N_G", "Gem_response_2N_G_vs_2N_C"),
        pairwise_stats("4N_C", "4N_G", "Gem_response_4N_G_vs_4N_C"),
        pairwise_stats("2N_C", "4N_C", "Baseline_4N_C_vs_2N_C"),
        pairwise_stats("2N_G", "4N_G", "Treated_4N_G_vs_2N_G"),
    ], ignore_index=True)

    pw_wide = pairwise_results.pivot_table(
        index=["Feature_ID", "Metabolite"],
        columns="comparison",
        values=["log2FC", "fold_change", "p_value", "FDR_BH"],
        aggfunc="first",
    )
    pw_wide.columns = [f"{metric}_{comparison}" for metric, comparison in pw_wide.columns]
    pw_wide = pw_wide.reset_index()

    anova_rows = []
    for idx in range(X_log2.shape[0]):
        dat = meta.copy()
        dat["value"] = X_log2.iloc[idx].values
        try:
            model = ols("value ~ C(ploidy) + C(treatment) + C(ploidy):C(treatment)", data=dat).fit()
            aov = sm.stats.anova_lm(model, typ=2)
            anova_rows.append({
                "Feature_ID": df.loc[idx, "Feature_ID"],
                "Metabolite": df.loc[idx, id_col],
                "p_ploidy": aov.loc["C(ploidy)", "PR(>F)"],
                "p_treatment": aov.loc["C(treatment)", "PR(>F)"],
                "p_interaction": aov.loc["C(ploidy):C(treatment)", "PR(>F)"],
                "F_ploidy": aov.loc["C(ploidy)", "F"],
                "F_treatment": aov.loc["C(treatment)", "F"],
                "F_interaction": aov.loc["C(ploidy):C(treatment)", "F"],
            })
        except Exception:
            anova_rows.append({
                "Feature_ID": df.loc[idx, "Feature_ID"], "Metabolite": df.loc[idx, id_col],
                "p_ploidy": np.nan, "p_treatment": np.nan, "p_interaction": np.nan,
                "F_ploidy": np.nan, "F_treatment": np.nan, "F_interaction": np.nan,
            })

    anova_df = pd.DataFrame(anova_rows)
    for col in ["p_ploidy", "p_treatment", "p_interaction"]:
        anova_df[col.replace("p_", "FDR_")] = multipletests(np.nan_to_num(anova_df[col], nan=1.0), method="fdr_bh")[1]

    oneway_rows = []
    for idx in range(X_log2.shape[0]):
        vals = [X_log2.loc[df.loc[idx, "Feature_ID"], group_cols[g]].values for g in ["2N_C", "2N_G", "4N_C", "4N_G"]]
        try:
            stat, p = f_oneway(*vals)
        except Exception:
            stat, p = np.nan, np.nan
        oneway_rows.append({"Feature_ID": df.loc[idx, "Feature_ID"], "Metabolite": df.loc[idx, id_col], "oneway_F": stat, "oneway_p": p})
    oneway_df = pd.DataFrame(oneway_rows)
    oneway_df["oneway_FDR_BH"] = multipletests(np.nan_to_num(oneway_df["oneway_p"], nan=1.0), method="fdr_bh")[1]

    log2fc_2n = X_log2[group_cols["2N_G"]].mean(axis=1) - X_log2[group_cols["2N_C"]].mean(axis=1)
    log2fc_4n = X_log2[group_cols["4N_G"]].mean(axis=1) - X_log2[group_cols["4N_C"]].mean(axis=1)
    diff_response = log2fc_2n - log2fc_4n

    diff_df = pd.DataFrame({
        "Feature_ID": df["Feature_ID"],
        "Metabolite": df[id_col],
        "log2FC_2N_G_vs_C": log2fc_2n.values,
        "FC_2N_G_vs_C": (2 ** log2fc_2n).values,
        "log2FC_4N_G_vs_C": log2fc_4n.values,
        "FC_4N_G_vs_C": (2 ** log2fc_4n).values,
        "Differential_response_delta_log2FC_2N_minus_4N": diff_response.values,
        "p_interaction": anova_df["p_interaction"].values,
        "FDR_interaction_BH": anova_df["FDR_interaction"].values,
    })

    results = (
        df[["Feature_ID", id_col]].rename(columns={id_col: "Metabolite"})
        .merge(pw_wide, on=["Feature_ID", "Metabolite"], how="left")
        .merge(anova_df, on=["Feature_ID", "Metabolite"], how="left")
        .merge(oneway_df, on=["Feature_ID", "Metabolite"], how="left")
        .merge(diff_df, on=["Feature_ID", "Metabolite"], how="left")
    )

    if "p_interaction" not in results.columns:
        results["p_interaction"] = results["p_interaction_x"] if "p_interaction_x" in results.columns else results["p_interaction_y"]

    return pairwise_results, anova_df, oneway_df, diff_df, results


def categorize_results(results):
    category_order = ["2N up", "2N down", "4N up", "4N down", "Common up", "Common down", "Opposite", "Interaction only"]
    p_2n = results["p_value_Gem_response_2N_G_vs_2N_C"].values
    p_4n = results["p_value_Gem_response_4N_G_vs_4N_C"].values
    a_2n = results["log2FC_Gem_response_2N_G_vs_2N_C"].values
    a_4n = results["log2FC_Gem_response_4N_G_vs_4N_C"].values
    p_int = results["p_interaction"].values
    delta = results["Differential_response_delta_log2FC_2N_minus_4N"].values

    order_cats, response_cats = [], []
    for i in range(len(results)):
        s2 = (p_2n[i] < ALPHA) and (abs(a_2n[i]) >= FC_CUT)
        s4 = (p_4n[i] < ALPHA) and (abs(a_4n[i]) >= FC_CUT)
        inter = (p_int[i] < ALPHA) and (abs(delta[i]) >= FC_CUT)
        if s2 and s4:
            if np.sign(a_2n[i]) == np.sign(a_4n[i]):
                order_cats.append("Common up" if a_2n[i] > 0 else "Common down")
                response_cats.append("Common gem response")
            else:
                order_cats.append("Opposite")
                response_cats.append("Opposite gem response")
        elif s2:
            order_cats.append("2N up" if a_2n[i] > 0 else "2N down")
            response_cats.append("2N-specific response")
        elif s4:
            order_cats.append("4N up" if a_4n[i] > 0 else "4N down")
            response_cats.append("4N-specific response")
        elif inter:
            order_cats.append("Interaction only")
            response_cats.append("Differential response only")
        else:
            order_cats.append("No strong response")
            response_cats.append("No strong response")

    results["Order_category_2fold_p05"] = order_cats
    results["Response_category_2fold_p05"] = response_cats
    results["cat_order"] = results["Order_category_2fold_p05"].map({c: i for i, c in enumerate(category_order)}).fillna(99).astype(int)
    return results


def make_figures(df, X_log2, meta, group_cols, ordered_cols, results, figdir):
    os.makedirs(figdir, exist_ok=True)

    # PCA
    scaled = StandardScaler().fit_transform(X_log2.T)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(scaled)
    pca_df = meta.copy()
    pca_df["PC1"] = pcs[:, 0]
    pca_df["PC2"] = pcs[:, 1]

    def make_pca(label_points, filename):
        fig, ax = plt.subplots(figsize=(7, 6))
        for group in ["2N_C", "2N_G", "4N_C", "4N_G"]:
            sub = pca_df[pca_df["group"] == group]
            marker = "x" if group.endswith("_G") else "o"
            ax.scatter(sub["PC1"], sub["PC2"], s=95, marker=marker, label=group)
            if label_points:
                for _, r in sub.iterrows():
                    ax.text(r["PC1"], r["PC2"], r["sample"], fontsize=8, ha="left", va="bottom")
        ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}% variance)")
        ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}% variance)")
        ax.set_title("PCA of log2-imputed metabolite intensities" + (" with sample labels" if label_points else ""))
        ax.legend(frameon=False)
        plt.tight_layout()
        savefig(os.path.join(figdir, filename))

    make_pca(False, "01_PCA_unlabeled.png")
    make_pca(True, "02_PCA_labeled.png")

    # Volcano
    def get_v(kind):
        if kind == "2N":
            v = results[["Feature_ID", "Metabolite", "log2FC_Gem_response_2N_G_vs_2N_C", "p_value_Gem_response_2N_G_vs_2N_C"]].copy()
            v = v.rename(columns={"log2FC_Gem_response_2N_G_vs_2N_C": "x", "p_value_Gem_response_2N_G_vs_2N_C": "p_value"})
            return v, "log2FC 2N-G / 2N-C", "2N: Gemcitabine vs control"
        if kind == "4N":
            v = results[["Feature_ID", "Metabolite", "log2FC_Gem_response_4N_G_vs_4N_C", "p_value_Gem_response_4N_G_vs_4N_C"]].copy()
            v = v.rename(columns={"log2FC_Gem_response_4N_G_vs_4N_C": "x", "p_value_Gem_response_4N_G_vs_4N_C": "p_value"})
            return v, "log2FC 4N-G / 4N-C", "4N: Gemcitabine vs control"
        v = results[["Feature_ID", "Metabolite", "Differential_response_delta_log2FC_2N_minus_4N", "p_interaction"]].copy()
        v = v.rename(columns={"Differential_response_delta_log2FC_2N_minus_4N": "x", "p_interaction": "p_value"})
        return v, "Δlog2FC = 2N response − 4N response", "Differential response: ploidy × gemcitabine"

    def add_volcano(ax, kind, labels=False):
        v, xlabel, title = get_v(kind)
        v["neglog10p"] = safe_neglog10(v["p_value"])
        v["significant"] = (v["p_value"] < ALPHA) & (v["x"].abs() >= FC_CUT)
        sig = v["significant"]
        ax.scatter(v.loc[~sig, "x"], v.loc[~sig, "neglog10p"], s=18, alpha=0.65, c="black", label="Not significant")
        ax.scatter(v.loc[sig, "x"], v.loc[sig, "neglog10p"], s=28, alpha=0.90, c="red", label="p<0.05 & ≥2-fold")
        ax.axhline(-np.log10(ALPHA), linestyle="--", linewidth=1, color="gray")
        ax.axvline(FC_CUT, linestyle="--", linewidth=1, color="gray")
        ax.axvline(-FC_CUT, linestyle="--", linewidth=1, color="gray")
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("-log10(p-value)")
        ax.legend(fontsize=8, frameon=False)
        if labels:
            top = v.loc[sig].copy()
            top["score"] = top["neglog10p"] + 0.5 * top["x"].abs()
            top = top.sort_values("score", ascending=False).head(12)
            for _, r in top.iterrows():
                ax.text(r["x"], r["neglog10p"], truncate_label(r["Metabolite"], 34), fontsize=7, ha="left", va="bottom")

    for labels, fname in [(False, "03_volcano_black_red_unlabeled.png"), (True, "04_volcano_black_red_labeled_top_hits.png")]:
        fig, axes = plt.subplots(1, 3, figsize=(20, 6))
        for ax, kind in zip(axes, ["2N", "4N", "interaction"]):
            add_volcano(ax, kind, labels=labels)
        plt.tight_layout()
        savefig(os.path.join(figdir, fname))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--out", default="metabolomics_2fold_full_package")
    parser.add_argument("--output-dir", help="Canonical output directory alias for --out.")
    args = parser.parse_args()
    if args.output_dir:
        args.out = args.output_dir

    figdir = os.path.join(args.out, "figures")
    tabledir = os.path.join(args.out, "tables")
    os.makedirs(figdir, exist_ok=True)
    os.makedirs(tabledir, exist_ok=True)

    df, X_log2, meta, group_cols, ordered_cols, id_col = load_and_preprocess(args.input)
    pairwise_results, anova_df, oneway_df, diff_df, results = compute_statistics(df, X_log2, meta, group_cols, id_col)
    results = categorize_results(results)

    pairwise_results.to_csv(os.path.join(tabledir, "pairwise_comparisons_2fold.csv"), index=False)
    anova_df.to_csv(os.path.join(tabledir, "two_way_anova.csv"), index=False)
    oneway_df.to_csv(os.path.join(tabledir, "one_way_anova.csv"), index=False)
    diff_df.to_csv(os.path.join(tabledir, "differential_response_2fold.csv"), index=False)
    results.to_csv(os.path.join(tabledir, "metabolomics_full_statistics_2fold.csv"), index=False)
    results["Response_category_2fold_p05"].value_counts().reset_index().to_csv(os.path.join(tabledir, "response_category_counts_2fold.csv"), index=False)

    make_figures(df, X_log2, meta, group_cols, ordered_cols, results, figdir)

    print("Analysis complete.")
    print("Output:", args.out)


if __name__ == "__main__":
    main()
