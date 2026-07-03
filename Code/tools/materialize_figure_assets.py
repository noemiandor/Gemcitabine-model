#!/usr/bin/env python3
"""Copy selected source panels into manuscript-facing figures/ folders."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from figure_output_contract import FIGURE_MANIFEST_COLUMNS, repo_root_from, write_tsv


PANEL_SPECS = [
    {
        "module": "gdsc",
        "source": "ploidy_enrichment_panels_primary_secondary_ABC.png",
        "figure": "Figure1",
        "panel": "1A-1C",
        "asset": "panel_1ABC_gdsc_ploidy_enrichment.png",
        "caption_role": "GDSC workflow and primary-class enrichment overview source panel",
    },
    {
        "module": "gdsc",
        "source": "drug_count_collapsed/figures/drug_count_collapsed_enrichment_shared_order.png",
        "figure": "Figure1",
        "panel": "1B-1C",
        "asset": "panel_1BC_gdsc_collapsed_drug_class_enrichment_heatmap.png",
        "caption_role": "Workbook-collapsed drug-category enrichment heatmap ordered by low-ploidy chemotherapy-agent significance",
    },
    {
        "module": "ccle",
        "source": "ccle_drug_ploidy_correlations_ic50.pdf",
        "figure": "Figure2",
        "panel": "2B",
        "asset": "panel_2B_ccle_drug_ploidy_correlations_ic50.pdf",
        "caption_role": "Local CCLE IC50 ploidy-support plot",
    },
    {
        "module": "drug_response",
        "source": "figures/gemcitabine_normalized_dose_response_fits.png",
        "figure": "Figure3",
        "panel": "3D-3F",
        "asset": "panel_3DEF_gemcitabine_normalized_dose_response_fits.png",
        "caption_role": "Gemcitabine normalized dose-response fits",
    },
    {
        "module": "drug_response",
        "source": "figures/gemcitabine_absolute_auc_vs_ploidy_mean.png",
        "figure": "Figure3",
        "panel": "3H",
        "asset": "panel_3H_gemcitabine_absolute_auc_vs_ploidy_mean.png",
        "caption_role": "Absolute AUC versus mean ploidy source plot",
    },
    {
        "module": "lci_overlays",
        "source": "figures/D5_2_timecourse_panel.png",
        "figure": "Figure3",
        "panel": "3I",
        "asset": "panel_3I_lci_timecourse_overlay.png",
        "caption_role": "Live-cell imaging overlay time-course panel",
        "optional": True,
    },
    {
        "module": "pkpd",
        "source": "dfdctp_signal_curve_combined_ploidy.png",
        "figure": "Figure4",
        "panel": "4B",
        "asset": "panel_4B_dfdctp_signal_driver_combined_ploidy.png",
        "caption_role": "Baseline-subtracted PK-derived dFdCTP signal driver",
    },
    {
        "module": "pkpd",
        "source": "dose_50_nm_2n_vs_4n.png",
        "figure": "Figure5",
        "panel": "5B",
        "asset": "panel_5B_dose_50_nm_2n_vs_4n.png",
        "caption_role": "50 nM live/dead fitted model source panel comparing 2N and 4N",
    },
    {
        "module": "pkpd",
        "source": "dfdctp_signal_curve_2n.png",
        "figure": "Figure5",
        "panel": "5C",
        "asset": "panel_5C_dfdctp_signal_curve_2n.png",
        "caption_role": "2N baseline-subtracted PK-derived dFdCTP signal driver",
    },
    {
        "module": "pkpd",
        "source": "dfdctp_signal_curve_4n.png",
        "figure": "Figure5",
        "panel": "5C",
        "asset": "panel_5C_dfdctp_signal_curve_4n.png",
        "caption_role": "4N baseline-subtracted PK-derived dFdCTP signal driver",
    },
    {
        "module": "pkpd",
        "source": "effective_dfdctp_signal_curve_combined_ploidy.png",
        "figure": "Figure5",
        "panel": "5C",
        "asset": "panel_5C_effective_dfdctp_signal_combined_ploidy.png",
        "caption_role": "Fitted beta/Hill-corrected effective dFdCTP signal",
    },
    {
        "module": "pkpd",
        "source": "ploidy_parameter_log2_fold_change.png",
        "figure": "Figure5",
        "panel": "5D",
        "asset": "panel_5D_ploidy_parameter_log2_fold_change.png",
        "caption_role": "4N versus 2N fitted parameter log2 fold-change",
    },
    {
        "module": "pkpd",
        "source": "dose_response_ploidy_comparison.png",
        "figure": "Figure5",
        "panel": "5D",
        "asset": "panel_5D_dose_response_ploidy_comparison.png",
        "caption_role": "Dose-wise fitted 2N/4N sensitivity comparison",
    },
    {
        "module": "metabolomics",
        "source": "figures/01_PCA_unlabeled.png",
        "figure": "Figure6",
        "panel": "6A",
        "asset": "panel_6A_metabolomics_pca_unlabeled.png",
        "caption_role": "Metabolomics PCA source plot",
    },
    {
        "module": "metabolomics_pathway",
        "source": "figures/corrected_curated_pathway_enrichment_heatmap_2fold.png",
        "figure": "Figure6",
        "panel": "6C",
        "asset": "panel_6C_corrected_curated_pathway_enrichment_heatmap_2fold.png",
        "caption_role": "Curated pathway enrichment heatmap",
    },
    {
        "module": "metabolomics_zscore",
        "source": "figures/ordered_response_heatmap_reproduced.png",
        "figure": "Figure6",
        "panel": "6D",
        "asset": "panel_6D_ordered_response_heatmap_reproduced.png",
        "caption_role": "Ordered metabolite response z-score heatmap",
    },
    {
        "module": "metabolomics",
        "source": "figures/03_volcano_black_red_unlabeled.png",
        "figure": "Figure6",
        "panel": "6E",
        "asset": "panel_6E_volcano_black_red_unlabeled.png",
        "caption_role": "Metabolomics volcano source plot",
    },
    {
        "module": "pkpd",
        "source": "cohort_joint_fit_2n.png",
        "figure": "Supplementary",
        "panel": "SuppFig1A",
        "asset": "panel_SuppFig1A_cohort_joint_fit_2n.png",
        "caption_role": "Supplementary 2N cohort fit source panel",
    },
    {
        "module": "pkpd",
        "source": "cohort_joint_fit_4n.png",
        "figure": "Supplementary",
        "panel": "SuppFig1B",
        "asset": "panel_SuppFig1B_cohort_joint_fit_4n.png",
        "caption_role": "Supplementary 4N cohort fit source panel",
    },
]

EXTERNAL_ROWS = [
    {
        "figure": "Figure2",
        "panel": "2C",
        "source_kind": "external",
        "caption_role": "NCI-60 / standard-therapy evidence cited externally",
        "not_regenerated_reason": "No local code or input data for this cited external panel.",
        "citation_or_uri": "Choudhary et al.",
    },
    {
        "figure": "Figure4",
        "panel": "4A",
        "source_kind": "manual_composite",
        "caption_role": "Gemcitabine metabolism schematic",
        "not_regenerated_reason": "Schematic/manual panel; no code-generated asset exists locally.",
    },
    {
        "figure": "Figure4",
        "panel": "4C",
        "source_kind": "external",
        "caption_role": "Representative imaging panel",
        "not_regenerated_reason": "Representative image selection/composite is external to the current reproducible code workflow.",
    },
    {
        "figure": "Figure5",
        "panel": "5E",
        "source_kind": "external",
        "caption_role": "Checkpoint immunoblot panel",
        "not_regenerated_reason": "No immunoblot quantification or plotting code exists locally.",
    },
]


def parse_module_run(values: list[str], repo_root: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"--module-run must be NAME=PATH, got: {value}")
        name, raw_path = value.split("=", 1)
        path = Path(raw_path)
        if not path.is_absolute():
            path = repo_root / path
        out[name] = path.resolve()
    return out


def rel(path: Path, repo_root: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def generated_row(spec: dict[str, object], source: Path, asset: Path, repo_root: Path, run_id: str) -> dict[str, str]:
    return {
        "figure": str(spec["figure"]),
        "panel": str(spec["panel"]),
        "asset_path": rel(asset, repo_root),
        "source_file": rel(source, repo_root),
        "source_kind": "generated_panel",
        "generated_by": "Manager.sh",
        "command": "materialize_figure_assets.py",
        "input_data": rel(source, repo_root),
        "result_run_dir": rel(source.parent if source.parent.name != "figures" else source.parent.parent, repo_root),
        "run_id": run_id,
        "caption_role": str(spec["caption_role"]),
        "asset_status": "generated",
        "not_regenerated_reason": "",
        "local_provenance_path": "",
        "citation_or_uri": "",
        "notes": "",
    }


def external_row(row: dict[str, str]) -> dict[str, str]:
    out = {column: "" for column in FIGURE_MANIFEST_COLUMNS}
    out.update(row)
    out["asset_status"] = out.get("source_kind", "external")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-root", type=Path, default=Path("figures"))
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--module-run", action="append", default=[], help="Module run path as NAME=PATH. May be repeated.")
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root else repo_root_from(Path.cwd())
    figure_root = args.figure_root
    if not figure_root.is_absolute():
        figure_root = repo_root / figure_root
    figure_root = figure_root.resolve()
    module_runs = parse_module_run(args.module_run, repo_root=repo_root)

    rows_by_figure: dict[str, list[dict[str, str]]] = {}
    for spec in PANEL_SPECS:
        module = str(spec["module"])
        if module not in module_runs:
            if spec.get("optional"):
                continue
            continue
        source = module_runs[module] / str(spec["source"])
        if not source.exists():
            if spec.get("optional"):
                continue
            raise FileNotFoundError(f"Missing source for {spec['figure']} {spec['panel']}: {source}")
        out_dir = figure_root / str(spec["figure"])
        out_dir.mkdir(parents=True, exist_ok=True)
        asset = out_dir / str(spec["asset"])
        if asset.exists() and not args.overwrite:
            raise FileExistsError(f"Asset exists; use --overwrite to replace: {asset}")
        shutil.copy2(source, asset)
        rows_by_figure.setdefault(str(spec["figure"]), []).append(
            generated_row(spec, source=source, asset=asset, repo_root=repo_root, run_id=args.run_id)
        )

    touched_figures = set(rows_by_figure)
    for row in EXTERNAL_ROWS:
        if row["figure"] in touched_figures:
            rows_by_figure.setdefault(row["figure"], []).append(external_row(row))

    for figure, rows in sorted(rows_by_figure.items()):
        manifest_path = figure_root / figure / "manifest.tsv"
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(manifest_path, rows, FIGURE_MANIFEST_COLUMNS)
        print(f"Wrote {manifest_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
