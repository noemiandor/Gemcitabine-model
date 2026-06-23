#!/usr/bin/env python3
"""

"Curated pathway-class enrichment across 2-fold response sets"

This script generates the curated pathway-class enrichment heatmap.
It performs the 2-fold response-set definition, curated pathway-class annotation,
hypergeometric over-representation testing, BH FDR correction, and plots a heatmap
of -log10(FDR).

Input format:
- Metabolite ID/name column: "row identity (all IDs)" or first column
- Replicate columns:
  2N_C1, 2N_C2, 2N_C3, 2N_C4
  2N_G1, 2N_G2, 2N_G3, 2N_G4
  4N_C1, 4N_C2, 4N_C3, 4N_C4
  4N_G1, 4N_G2, 4N_G3, 4N_G4

Run:
python generate_curated_pathway_enrichment_heatmap_2fold.py \
    --input Metabolomics_2N_4N_Full.xlsm \
    --outdir pathway_enrichment_heatmap_output
"""

import argparse
import os
import re
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import ttest_ind, hypergeom
from statsmodels.stats.multitest import multipletests
from statsmodels.formula.api import ols
import statsmodels.api as sm


ALPHA = 0.05
LOG2FC_CUTOFF = 1.0  # log2(2), i.e. 2-fold


def load_and_preprocess(input_file):
    """Load metabolomics data, detect replicates, impute zeros/missing, log2-transform."""
    df_raw = pd.read_excel(input_file, sheet_name=0)

    id_candidates = ["row identity (all IDs)", "Metabolite ID", "Metabolite", "metabolite"]
    id_col = next((c for c in id_candidates if c in df_raw.columns), df_raw.columns[0])

    replicate_pattern = re.compile(r"^(2N|4N)_[CG][1-4]$")
    rep_cols = [c for c in df_raw.columns if replicate_pattern.match(str(c))]

    expected_cols = [
        f"{p}_{t}{i}"
        for p in ["2N", "4N"]
        for t in ["C", "G"]
        for i in range(1, 5)
    ]
    ordered_cols = [c for c in expected_cols if c in rep_cols]

    if len(ordered_cols) != 16:
        raise ValueError(f"Expected 16 replicate columns, found {len(ordered_cols)}: {ordered_cols}")

    df = df_raw[[id_col] + ordered_cols].copy()
    df[id_col] = df[id_col].astype(str).str.strip()

    # Keep duplicate metabolite labels as separate measured features.
    counts = df[id_col].value_counts()
    df["Feature_ID"] = [
        f"{name} | row_{i}" if counts[name] > 1 else name
        for i, name in zip(df_raw.index, df[id_col])
    ]

    X_raw = df[ordered_cols].apply(pd.to_numeric, errors="coerce")

    # Drop rows with no positive signal.
    keep = (X_raw.fillna(0) > 0).sum(axis=1) > 0
    df = df.loc[keep].reset_index(drop=True)
    X_raw = X_raw.loc[keep].reset_index(drop=True)

    # Impute zero/missing values by half of the row-wise minimum positive value.
    X_imp = X_raw.copy()
    row_min_positive = X_imp.where(X_imp > 0).min(axis=1)
    global_min_positive = X_imp.where(X_imp > 0).min().min()
    row_min_positive = row_min_positive.fillna(global_min_positive)

    for i in range(X_imp.shape[0]):
        fill_value = row_min_positive.iloc[i] / 2.0
        X_imp.iloc[i, :] = X_imp.iloc[i, :].replace(0, np.nan).fillna(fill_value)

    X_log2 = np.log2(X_imp)
    X_log2.index = df["Feature_ID"]

    # Sample metadata for interaction ANOVA.
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
            "replicate": rest[1:],
        })
    meta = pd.DataFrame(meta)

    group_cols = {
        "2N_C": [c for c in ordered_cols if c.startswith("2N_C")],
        "2N_G": [c for c in ordered_cols if c.startswith("2N_G")],
        "4N_C": [c for c in ordered_cols if c.startswith("4N_C")],
        "4N_G": [c for c in ordered_cols if c.startswith("4N_G")],
    }

    return df, X_log2, meta, group_cols, ordered_cols, id_col


def compute_response_sets(df, X_log2, meta, group_cols, id_col):
    """Compute 2N/4N gemcitabine responses and the interaction p-value."""
    log2fc_2n = (
        X_log2[group_cols["2N_G"]].mean(axis=1)
        - X_log2[group_cols["2N_C"]].mean(axis=1)
    )
    log2fc_4n = (
        X_log2[group_cols["4N_G"]].mean(axis=1)
        - X_log2[group_cols["4N_C"]].mean(axis=1)
    )

    p_2n = ttest_ind(
        X_log2[group_cols["2N_G"]].values,
        X_log2[group_cols["2N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    ).pvalue

    p_4n = ttest_ind(
        X_log2[group_cols["4N_G"]].values,
        X_log2[group_cols["4N_C"]].values,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    ).pvalue

    p_interaction = []
    for idx in range(X_log2.shape[0]):
        dat = meta.copy()
        dat["value"] = X_log2.iloc[idx].values
        try:
            model = ols(
                "value ~ C(ploidy) + C(treatment) + C(ploidy):C(treatment)",
                data=dat,
            ).fit()
            aov = sm.stats.anova_lm(model, typ=2)
            p_interaction.append(aov.loc["C(ploidy):C(treatment)", "PR(>F)"])
        except Exception:
            p_interaction.append(np.nan)

    p_interaction = np.asarray(p_interaction, dtype=float)
    diff_response = log2fc_2n.values - log2fc_4n.values

    results = pd.DataFrame({
        "Feature_ID": df["Feature_ID"],
        "Metabolite": df[id_col],
        "log2FC_2N_G_vs_C": log2fc_2n.values,
        "p_2N_G_vs_C": p_2n,
        "log2FC_4N_G_vs_C": log2fc_4n.values,
        "p_4N_G_vs_C": p_4n,
        "Differential_response_2N_minus_4N": diff_response,
        "p_interaction": p_interaction,
    })

    sets = {
        "2N_up_after_gem_2fold": (log2fc_2n.values >= LOG2FC_CUTOFF) & (p_2n < ALPHA),
        "2N_down_after_gem_2fold": (log2fc_2n.values <= -LOG2FC_CUTOFF) & (p_2n < ALPHA),
        "4N_up_after_gem_2fold": (log2fc_4n.values >= LOG2FC_CUTOFF) & (p_4n < ALPHA),
        "4N_down_after_gem_2fold": (log2fc_4n.values <= -LOG2FC_CUTOFF) & (p_4n < ALPHA),
        "Differential_response_2fold": (np.abs(diff_response) >= LOG2FC_CUTOFF) & (p_interaction < ALPHA),
        "Interaction_p05": (p_interaction < ALPHA),
    }

    return results, sets


def curated_pathway_annotation(metabolite_name):
    """
    Curated keyword-based pathway/metabolite-class annotation.
    This is the same class-level strategy used for the displayed heatmap.
    """
    pathways = {
        "Kennedy/choline-phospholipid": [
            "phosphocholine", "citicoline", "cdp-choline", "diphosphocholine",
            "phosphorylethanolamine", "cdp-ethanolamine", "choline"
        ],
        "Purine metabolism": [
            "adenosine", "adenine", "amp", "adp", "atp", "guanine", "guanosine",
            "gmp", "gdp", "gtp", "inosine", "imp", "hypoxanthine", "xanthosine",
            "xanthine", "damp", "dgmp"
        ],
        "Pyrimidine metabolism": [
            "uridine", "uracil", "ump", "udp", "utp", "cytidine", "cytosine",
            "cmp", "cdp", "ctp", "thymidine", "tmp", "dtmp", "orotate",
            "dihydroorotate", "dfdcmp"
        ],
        "Folate/one-carbon/methylation": [
            "folic", "methyltetrahydrofolic", "methionine", "s-adenosyl",
            "homocysteine", "sah", "sam", "betaine"
        ],
        "PPP/glycolysis/carbohydrate": [
            "glucose", "fructose", "phosphogluconic", "sedoheptulose", "ribose",
            "lactate", "pyruvic", "phosphoenolpyruvic", "galactitol", "mannitol",
            "sorbitol", "trehalose", "maltose", "melibiose"
        ],
        "Carnitine/FAO": [
            "carnitine", "acylcarnitine", "isovaleryl", "lauroyl", "decanoyl",
            "hexanoyl", "myristoyl", "stearoyl", "propionylcarnitine"
        ],
        "Fatty acids/lipids": [
            "palmitate", "palmitic", "stearidonic", "arachidonic", "myristic",
            "oleoyl", "glycerol", "dihome", "eicosenoic", "docosadienoate",
            "caprylic", "linoleic"
        ],
        "Amino acid/nitrogen": [
            "glutamine", "glutamate", "aspartate", "asparagine", "tryptophan",
            "phenylalanine", "serine", "alanine", "citrulline", "guanidino",
            "taurine", "hypotaurine", "aminobutanoate", "beta-alanine"
        ],
        "Cofactors/vitamins/redox": [
            "nadp", "nad", "pantothenic", "riboflavin", "nicotinamide",
            "coenzyme", "pterin"
        ],
    }

    n = str(metabolite_name).lower()
    hits = []
    for pathway, keywords in pathways.items():
        if any(k.lower() in n for k in keywords):
            hits.append(pathway)

    return "; ".join(hits) if hits else "Other"


def hypergeometric_enrichment(results, sets):
    """Run curated pathway-class over-representation enrichment."""
    results = results.copy()
    results["Curated_pathway_class"] = results["Metabolite"].apply(curated_pathway_annotation)

    pathway_classes = [
        "Kennedy/choline-phospholipid",
        "Purine metabolism",
        "Pyrimidine metabolism",
        "Folate/one-carbon/methylation",
        "PPP/glycolysis/carbohydrate",
        "Carnitine/FAO",
        "Fatty acids/lipids",
        "Amino acid/nitrogen",
        "Cofactors/vitamins/redox",
        "Other",
    ]

    rows = []
    M = len(results)

    for set_name, mask in sets.items():
        selected = results.loc[mask].copy()
        N = len(selected)
        if N == 0:
            continue

        for pathway in pathway_classes:
            if pathway == "Other":
                K = (results["Curated_pathway_class"] == "Other").sum()
                x = (selected["Curated_pathway_class"] == "Other").sum()
            else:
                K = results["Curated_pathway_class"].str.contains(re.escape(pathway), na=False).sum()
                x = selected["Curated_pathway_class"].str.contains(re.escape(pathway), na=False).sum()

            p_value = hypergeom.sf(x - 1, M, K, N) if x > 0 and K > 0 else 1.0

            rows.append({
                "Set": set_name,
                "Pathway_class": pathway,
                "Selected_count": int(x),
                "Background_count": int(K),
                "Set_size": int(N),
                "p_value": p_value,
            })

    enrichment = pd.DataFrame(rows)
    enrichment["FDR_BH"] = np.nan

    # Match the prior workflow: BH correction within each response set.
    for set_name in enrichment["Set"].unique():
        idx = enrichment["Set"] == set_name
        enrichment.loc[idx, "FDR_BH"] = multipletests(
            enrichment.loc[idx, "p_value"].astype(float),
            method="fdr_bh"
        )[1]

    return enrichment, results


def make_heatmap(enrichment, output_png, output_pdf=None, output_matrix_csv=None):
    """Plot the pathway enrichment heatmap exactly in the prior style."""
    enr = enrichment[
        (enrichment["Pathway_class"] != "Other")
        & (enrichment["Selected_count"] > 0)
    ].copy()

    enr["score"] = -np.log10(np.clip(enr["FDR_BH"].astype(float), 1e-300, None))

    # This filtering reproduces the displayed figure: show only pathways with at
    # least nominal enrichment or FDR trend in at least one response set.
    enr_show = enr[(enr["p_value"] < 0.25) | (enr["FDR_BH"] < 0.25)].copy()
    if len(enr_show) == 0:
        enr_show = enr.sort_values("p_value").head(25)

    pivot = enr_show.pivot_table(
        index="Pathway_class",
        columns="Set",
        values="score",
        aggfunc="max"
    ).fillna(0)

    set_order = [
        "2N_up_after_gem_2fold",
        "2N_down_after_gem_2fold",
        "4N_up_after_gem_2fold",
        "4N_down_after_gem_2fold",
        "Differential_response_2fold",
        "Interaction_p05",
    ]
    pivot = pivot[[s for s in set_order if s in pivot.columns]]

    # Alphabetical row ordering matches the attached heatmap.
    pivot = pivot.sort_index()

    if output_matrix_csv:
        pivot.to_csv(output_matrix_csv)

    fig, ax = plt.subplots(figsize=(12, max(5, 0.45 * len(pivot))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="viridis", interpolation="nearest")

    ax.set_title("Curated pathway-class enrichment across 2-fold response sets", fontsize=14)
    ax.set_xticks(np.arange(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=90)
    ax.set_yticks(np.arange(pivot.shape[0]))
    ax.set_yticklabels(pivot.index)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("-log10(FDR)")

    ax.set_xlabel("Response set")
    ax.set_ylabel("Pathway / metabolite class")

    plt.tight_layout()
    fig.savefig(output_png, dpi=300, bbox_inches="tight")
    if output_pdf:
        fig.savefig(output_pdf, bbox_inches="tight")
    plt.close(fig)

    return pivot


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input metabolomics .xlsm/.xlsx file")
    parser.add_argument("--outdir", default="pathway_enrichment_heatmap_output", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df, X_log2, meta, group_cols, ordered_cols, id_col = load_and_preprocess(args.input)
    results, sets = compute_response_sets(df, X_log2, meta, group_cols, id_col)
    enrichment, annotated_results = hypergeometric_enrichment(results, sets)

    annotated_results_csv = os.path.join(args.outdir, "annotated_metabolites_with_response_stats.csv")
    enrichment_csv = os.path.join(args.outdir, "curated_pathway_enrichment_2fold.csv")
    matrix_csv = os.path.join(args.outdir, "curated_pathway_enrichment_heatmap_matrix.csv")
    output_png = os.path.join(args.outdir, "curated_pathway_enrichment_heatmap_2fold_reproduced.png")
    output_pdf = os.path.join(args.outdir, "curated_pathway_enrichment_heatmap_2fold_reproduced.pdf")

    annotated_results.to_csv(annotated_results_csv, index=False)
    enrichment.to_csv(enrichment_csv, index=False)

    pivot = make_heatmap(
        enrichment,
        output_png=output_png,
        output_pdf=output_pdf,
        output_matrix_csv=matrix_csv
    )

    print(f"Loaded {len(df)} metabolite features.")
    print(f"Heatmap matrix shape: {pivot.shape[0]} pathway classes x {pivot.shape[1]} response sets.")
    print("Saved:")
    print(" -", output_png)
    print(" -", output_pdf)
    print(" -", enrichment_csv)
    print(" -", matrix_csv)
    print(" -", annotated_results_csv)


if __name__ == "__main__":
    main()
