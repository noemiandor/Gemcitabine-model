#!/usr/bin/env python3
"""Recompute primary-class enrichment with joint Monte Carlo permutations."""

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd


def enrichment_scores(membership):
    """One-sided unweighted KS-like maximum running-sum scores."""
    membership = np.asarray(membership, dtype=np.int16)
    n = membership.shape[-2]
    sizes = membership.sum(axis=-2)
    cumulative = membership.cumsum(axis=-2)
    positions = np.arange(1, n + 1, dtype=float)
    shape = (1,) * (membership.ndim - 2) + (n, 1)
    positions = positions.reshape(shape)
    with np.errstate(divide="ignore", invalid="ignore"):
        running = cumulative / sizes[..., None, :] - (
            positions - cumulative
        ) / (n - sizes[..., None, :])
    scores = np.nanmax(running, axis=-2)
    invalid = (sizes == 0) | (sizes == n)
    return np.where(invalid, np.nan, scores)


def bh_adjust(pvalues):
    """Benjamini-Hochberg adjusted p-values, preserving NaNs."""
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full_like(pvalues, np.nan)
    finite_idx = np.flatnonzero(np.isfinite(pvalues))
    finite = pvalues[finite_idx]
    order = np.argsort(finite, kind="mergesort")
    ranked = finite[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[finite_idx] = restored
    return out


def simes_pvalue(pvalues):
    finite = np.sort(np.asarray(pvalues, dtype=float)[np.isfinite(pvalues)])
    if not len(finite):
        return np.nan
    return min(1.0, np.min(len(finite) * finite / np.arange(1, len(finite) + 1)))


def analyze_context(table, groups, permutations, rng, batch_size):
    table = table.sort_values("correlation_value", ascending=False, kind="mergesort")
    labels = table["group"].astype(str).to_numpy()
    membership = np.column_stack([labels == group for group in groups])
    observed_low = enrichment_scores(membership)
    observed_high = enrichment_scores(membership[::-1])
    exceed_low = np.zeros(len(groups), dtype=np.int64)
    exceed_high = np.zeros(len(groups), dtype=np.int64)

    completed = 0
    while completed < permutations:
        current = min(batch_size, permutations - completed)
        permutation_order = np.argsort(
            rng.random((current, len(table))), axis=1, kind="quicksort"
        )
        permuted = membership[permutation_order]
        permuted_low = enrichment_scores(permuted)
        permuted_high = enrichment_scores(permuted[:, ::-1, :])
        exceed_low += np.sum(permuted_low >= observed_low, axis=0)
        exceed_high += np.sum(permuted_high >= observed_high, axis=0)
        completed += current

    denominator = permutations + 1
    return {
        "low_pvalue": (exceed_low + 1) / denominator,
        "high_pvalue": (exceed_high + 1) / denominator,
        "low_score": observed_low,
        "high_score": observed_high,
    }


def matrix_from_long(table, value, direction, cancers, groups):
    subset = table[table["direction"] == direction]
    return subset.pivot(index="cancer_type", columns="group", values=value).reindex(
        index=cancers, columns=groups
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-drugs", required=True)
    parser.add_argument("--template-workbook", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--permutations", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument(
        "--analysis-kind",
        choices=("primary", "collapsed"),
        default="primary",
    )
    args = parser.parse_args()

    started = time.time()
    out_dir = Path(args.out_dir)
    tables_dir = out_dir / "tables"
    metadata_dir = out_dir / "metadata"
    tables_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    selected = pd.read_csv(args.selected_drugs, sep="\t")
    selected = selected[selected["direction"] == "low_ploidy_sensitive"].copy()
    if args.analysis_kind == "collapsed":
        group_column = "category_label"
        statistics_name = "drug_count_collapsed_enrichment_Z_SCORE.tsv"
        significant_name = "drug_count_collapsed_significant_cells_fdr0.05_Z_SCORE.tsv"
        effect_name = "drug_count_collapsed_effect_scores_Z_SCORE.tsv"
        family_name = "drug_count_collapsed_family_simes_bh_sensitivity_Z_SCORE.tsv"
        workbook_name = "drug_count_collapsed_drugsVsPloidyCorr_Z_SCORE.xlsx"
    else:
        group_column = "group"
        statistics_name = "class_enrichment_primary_secondary_Z_SCORE.tsv"
        significant_name = "class_enrichment_primary_secondary_significant_cells_fdr0.05_Z_SCORE.tsv"
        effect_name = "class_enrichment_primary_secondary_effect_scores_Z_SCORE.tsv"
        family_name = "class_enrichment_primary_secondary_family_simes_bh_sensitivity_Z_SCORE.tsv"
        workbook_name = "drugsVsPloidyCorr_primary_secondary_Z_SCORE.xlsx"
    if group_column not in selected.columns:
        raise ValueError(f"Selected-drug table is missing required column: {group_column}")
    selected["group"] = selected[group_column].astype(str)
    template_low = pd.read_excel(args.template_workbook, sheet_name="lowpIsSens", index_col=0)
    class_counts = pd.read_excel(args.template_workbook, sheet_name="drugClassCounts")
    count_key = next(
        key for key in ("primary_anticancer_class", "category_label", "group")
        if key in class_counts.columns
    )
    groups = class_counts[count_key].astype(str).tolist()
    cancers = template_low.index.astype(str).tolist()

    seed_sequences = np.random.SeedSequence(args.seed).spawn(len(cancers))
    results = {}
    for cancer, seed_sequence in zip(cancers, seed_sequences):
        context = selected[selected["cancer_type"].astype(str) == cancer]
        results[cancer] = analyze_context(
            context,
            groups,
            args.permutations,
            np.random.default_rng(seed_sequence),
            args.batch_size,
        )

    rows = []
    for direction, prefix in (
        ("low_ploidy_sensitive", "low"),
        ("high_ploidy_sensitive", "high"),
    ):
        for cancer in cancers:
            result = results[cancer]
            for group_idx, group in enumerate(groups):
                rows.append(
                    {
                        "direction": direction,
                        "metric": "Z_SCORE",
                        "cancer_type": cancer,
                        "group": group,
                        "pvalue": result[f"{prefix}_pvalue"][group_idx],
                        "enrichment_score": result[f"{prefix}_score"][group_idx],
                    }
                )
    enrichment = pd.DataFrame(rows)
    enrichment["pvalue_for_fdr"] = enrichment["pvalue"]
    enrichment["qvalue_bh"] = bh_adjust(enrichment["pvalue"].to_numpy())
    enrichment.to_csv(
        tables_dir / statistics_name,
        sep="\t",
        index=False,
    )
    enrichment.loc[enrichment["qvalue_bh"] <= 0.05].to_csv(
        tables_dir / significant_name,
        sep="\t",
        index=False,
    )
    enrichment[["direction", "cancer_type", "group", "enrichment_score"]].to_csv(
        tables_dir / effect_name,
        sep="\t",
        index=False,
    )

    family_rows = []
    for (direction, cancer), family in enrichment.groupby(["direction", "cancer_type"], sort=False):
        family_rows.append(
            {"direction": direction, "cancer_type": cancer, "simes_pvalue": simes_pvalue(family["pvalue"])}
        )
    family_table = pd.DataFrame(family_rows)
    family_table["family_qvalue_bh"] = bh_adjust(family_table["simes_pvalue"].to_numpy())
    family_table.to_csv(
        tables_dir / family_name,
        sep="\t",
        index=False,
    )

    matrices = {
        "rawLowpIsSens": matrix_from_long(enrichment, "pvalue", "low_ploidy_sensitive", cancers, groups),
        "rawHighpIsSens": matrix_from_long(enrichment, "pvalue", "high_ploidy_sensitive", cancers, groups),
        "lowpIsSens": matrix_from_long(enrichment, "qvalue_bh", "low_ploidy_sensitive", cancers, groups),
        "highpIsSens": matrix_from_long(enrichment, "qvalue_bh", "high_ploidy_sensitive", cancers, groups),
        "lowEnrichmentScore": matrix_from_long(
            enrichment, "enrichment_score", "low_ploidy_sensitive", cancers, groups
        ),
        "highEnrichmentScore": matrix_from_long(
            enrichment, "enrichment_score", "high_ploidy_sensitive", cancers, groups
        ),
    }
    minimum_p = 1 / (args.permutations + 1)
    scope = "both directions, all cancer types, and all drug categories"
    metadata = pd.DataFrame(
        {
            "key": [
                "permutation_count", "minimum_permutation_pvalue", "zero_pvalue_floor",
                "multiple_testing_method", "adjustment_scope", "significance_cutoff",
                "heatmap_value", "permutation_pvalue_formula", "permutation_strategy",
                "effect_size", "source_selected_drugs", "seed",
            ],
            "value": [
                args.permutations, minimum_p, minimum_p, "Benjamini-Hochberg", scope, 0.05,
                "BH-FDR q-value", "(b + 1) / (B + 1)",
                "synchronized drug-class label permutations within each cancer type and across low/high directions",
                "one-sided unweighted KS-like maximum running-sum enrichment score",
                str(Path(args.selected_drugs)), args.seed,
            ],
        }
    )
    workbook = out_dir / workbook_name
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        for sheet, matrix in matrices.items():
            matrix.to_excel(writer, sheet_name=sheet)
        metadata.to_excel(writer, sheet_name="analysisMetadata", index=False)
        class_counts.to_excel(writer, sheet_name="drugClassCounts", index=False)
        pd.read_excel(args.template_workbook, sheet_name="cancerTypeCounts").to_excel(
            writer, sheet_name="cancerTypeCounts", index=False
        )
    matrices["lowpIsSens"].to_csv(tables_dir / "lowpIsSens_bh_qvalues.tsv", sep="\t")
    matrices["highpIsSens"].to_csv(tables_dir / "highpIsSens_bh_qvalues.tsv", sep="\t")
    (out_dir / "drugsVsPloidyCorr.xlsx").write_bytes(workbook.read_bytes())
    (tables_dir / workbook.name).write_bytes(workbook.read_bytes())

    summary = {
        "elapsed_seconds": time.time() - started,
        "tests": int(len(enrichment)),
        "permutations": args.permutations,
        "minimum_pvalue": float(enrichment["pvalue"].min()),
        "minimum_qvalue": float(enrichment["qvalue_bh"].min()),
        "significant_cells_q_le_0.05": int((enrichment["qvalue_bh"] <= 0.05).sum()),
        "significant_families_q_le_0.05": int((family_table["family_qvalue_bh"] <= 0.05).sum()),
    }
    (metadata_dir / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    metadata.to_csv(metadata_dir / "run_config.tsv", sep="\t", index=False)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
