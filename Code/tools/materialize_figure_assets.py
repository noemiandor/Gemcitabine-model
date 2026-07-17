#!/usr/bin/env python3
"""Copy selected source panels into manuscript-facing figures/ folders."""

from __future__ import annotations

import argparse
import shutil
from collections import Counter
from pathlib import Path

from figure_output_contract import (
    FIGURE_MANIFEST_COLUMNS,
    first_local_path,
    read_tsv,
    repo_root_from,
    sha256_file,
    validate_expected_panel_set,
    validate_module_manifest,
    validate_run_id,
    write_tsv,
)


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
        "caption_role": "Workbook-collapsed drug-category enrichment heatmap with drug-class and cancer-type cell-line count marginals",
    },
    {
        "module": "ccle",
        "source": "ccle_drug_ploidy_correlations_z_score.pdf",
        "figure": "Figure2",
        "panel": "2B",
        "asset": "panel_2B_ccle_drug_ploidy_correlations_z_score.pdf",
        "caption_role": "Local CCLE BreastCancerDrugSensitivity Z-score ploidy-support plot",
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
        "source": "dfdctp_signal_curve_stacked_2n_4n.png",
        "figure": "Figure5",
        "panel": "5B",
        "asset": "panel_5B_dfdctp_signal_curve_stacked_2n_4n.png",
        "caption_role": "Stacked baseline-subtracted PK-derived dFdCTP signal driver curves for 2N and 4N",
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
        "source": "dose_25_nm_2n_vs_4n.png",
        "figure": "Figure5",
        "panel": "5E",
        "asset": "panel_5E_dose_25_nm_2n_vs_4n.png",
        "caption_role": "25 nM live/dead fitted model source panel comparing 2N and 4N",
    },
    {
        "module": "pkpd",
        "source": "effective_dfdctp_signal_curve_combined_ploidy.png",
        "figure": "Figure5",
        "panel": "5F",
        "asset": "panel_5F_effective_dfdctp_signal_combined_ploidy.png",
        "caption_role": "Fitted beta/Hill-corrected effective dFdCTP signal",
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
    {
        "module": "drug_response",
        "source": "figures/gemcitabine_absolute_ec50_vs_ploidy_mean.png",
        "figure": "Supplementary",
        "panel": "SuppFig2A",
        "asset": "panel_SuppFig2A_gemcitabine_absolute_ec50_vs_ploidy_mean.png",
        "caption_role": "Absolute gemcitabine EC50 versus mean ploidy",
    },
    {
        "module": "drug_response",
        "source": "figures/gemcitabine_absolute_ic50_vs_ploidy.png",
        "figure": "Supplementary",
        "panel": "SuppFig2B",
        "asset": "panel_SuppFig2B_gemcitabine_absolute_ic50_vs_ploidy.png",
        "caption_role": "Absolute gemcitabine IC50 versus ploidy",
    },
    {
        "module": "drug_response",
        "source": "figures/gemcitabine_delta_auc_vs_delta_ploidy.png",
        "figure": "Supplementary",
        "panel": "SuppFig2C",
        "asset": "panel_SuppFig2C_gemcitabine_delta_auc_vs_delta_ploidy.png",
        "caption_role": "Paired gemcitabine delta AUC versus delta ploidy",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7A_day17_tgi_calculation.pdf",
        "figure": "Figure7",
        "panel": "7A",
        "asset": "panel_7A_day17_tgi_calculation.pdf",
        "caption_role": "Tumor-growth trajectories explaining the Day-17 TGI calculation",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7B_cellcycle_selected_ecdf_comparisons.pdf",
        "figure": "Figure7",
        "panel": "7B",
        "asset": "panel_7B_cellcycle_selected_ecdf_comparisons.pdf",
        "caption_role": "Selected CellCycle mean-ECDF comparisons",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7C_day17_tgi_by_initial_ploidy.pdf",
        "figure": "Figure7",
        "panel": "7C",
        "asset": "panel_7C_day17_tgi_by_initial_ploidy.pdf",
        "caption_role": "Day-17 TGI in initial 2N versus 4N treated tumors",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf",
        "figure": "Figure7",
        "panel": "7D",
        "asset": "panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf",
        "caption_role": "Within-dose-centered TGI and CellCycle ECDF-shift association",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7E_day17_tgi_vs_mean_etp.pdf",
        "figure": "Figure7",
        "panel": "7E",
        "asset": "panel_7E_day17_tgi_vs_mean_etp.pdf",
        "caption_role": "Day-17 TGI versus sample mean endpoint ploidy",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7F_pseudotime_state_pathway_activity.pdf",
        "figure": "Figure7",
        "panel": "7F",
        "asset": "panel_7F_pseudotime_state_pathway_activity.pdf",
        "caption_role": "Pathway activity across the accumulated CellCycle pseudotime state",
        "optional": True,
    },
]

STRICT_FIGURE_MODULES = {"in_vivo_figure7"}
FIGURE_SUFFIXES = {".pdf", ".png", ".jpg", ".jpeg", ".svg", ".tif", ".tiff"}

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
        "caption_role": "dFdU time-course measurements",
        "not_regenerated_reason": "dFdU measurement/composite panel is external to the current reproducible code workflow.",
    },
    {
        "figure": "Figure4",
        "panel": "4D",
        "source_kind": "external",
        "caption_role": "Representative imaging panel",
        "not_regenerated_reason": "Representative image selection/composite is external to the current reproducible code workflow.",
    },
    {
        "figure": "Figure5",
        "panel": "5A",
        "source_kind": "manual_composite",
        "caption_role": "Live/dead imaging and segmentation example panel",
        "not_regenerated_reason": "Panel is manually selected/assembled from live-cell imaging outputs and is not regenerated by the manuscript asset materializer.",
    },
    {
        "figure": "Figure5",
        "panel": "5C",
        "source_kind": "manual_composite",
        "caption_role": "Delay-aware live/dead model schematic",
        "not_regenerated_reason": "Schematic/manual panel; executable model code exists locally, but this graphic is not code-generated.",
    },
]


def parse_module_run(values: list[str], repo_root: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for value in values:
        if "=" not in value:
            raise ValueError(f"--module-run must be NAME=PATH, got: {value}")
        name, raw_path = value.split("=", 1)
        if not name or not raw_path:
            raise ValueError(f"--module-run must be NAME=PATH, got: {value}")
        if name in out:
            raise ValueError(f"Duplicate --module-run for module: {name}")
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


def generated_row(
    spec: dict[str, object],
    source: Path,
    asset: Path,
    result_run_root: Path,
    repo_root: Path,
    source_run_id: str,
    operation_id: str,
) -> dict[str, str]:
    return {
        "figure": str(spec["figure"]),
        "panel": str(spec["panel"]),
        "asset_path": rel(asset, repo_root),
        "source_file": rel(source, repo_root),
        "source_kind": "generated_panel",
        "generated_by": "Manager.sh",
        "command": "materialize_figure_assets.py",
        "input_data": rel(source, repo_root),
        "result_run_dir": rel(result_run_root, repo_root),
        "run_id": source_run_id,
        "caption_role": str(spec["caption_role"]),
        "asset_status": "generated",
        "not_regenerated_reason": "",
        "local_provenance_path": "",
        "citation_or_uri": "",
        "notes": f"materialization_operation_id={operation_id}",
    }


def external_row(row: dict[str, str]) -> dict[str, str]:
    out = {column: "" for column in FIGURE_MANIFEST_COLUMNS}
    out.update(row)
    out["asset_status"] = out.get("source_kind", "external")
    return out


def validate_strict_source_run(
    module: str,
    run_root: Path,
    selected_specs: list[dict[str, object]],
    repo_root: Path,
    source_run_id: str,
) -> None:
    if module not in STRICT_FIGURE_MODULES:
        return
    expected_name = f"{source_run_id}_figure7"
    if run_root.name != expected_name:
        raise ValueError(f"{module} source run must be named {expected_name}, got {run_root.name}")

    input_manifest = run_root / "metadata" / "input_manifest.tsv"
    if not input_manifest.is_file():
        raise FileNotFoundError(f"Missing source input manifest: {input_manifest}")
    input_errors = validate_module_manifest(input_manifest, repo_root=repo_root)
    if input_errors:
        raise ValueError("Invalid source input manifest:\n" + "\n".join(input_errors))
    _, input_rows = read_tsv(input_manifest)
    if any(
        row.get("module") != module or row.get("command_id") != source_run_id
        for row in input_rows
    ):
        raise ValueError(
            f"Source input-manifest provenance must use module={module} "
            f"and command_id={source_run_id}"
        )

    output_manifest = run_root / "metadata" / "output_manifest.tsv"
    if not output_manifest.is_file():
        raise FileNotFoundError(f"Missing source output manifest: {output_manifest}")
    errors = validate_module_manifest(output_manifest, repo_root=repo_root, output_root=run_root)
    if errors:
        raise ValueError("Invalid source output manifest:\n" + "\n".join(errors))

    expected_sources = {
        (run_root / str(spec["source"])).resolve()
        for spec in selected_specs
        if str(spec["module"]) == module
        and (not spec.get("optional") or (run_root / str(spec["source"])).is_file())
    }
    has_panel_f = any(path.name == "panel_7F_pseudotime_state_pathway_activity.pdf" for path in expected_sources)
    expected_panel_set = "a-f" if has_panel_f else "a-e"
    run_config = run_root / "metadata" / "run_config.tsv"
    if not run_config.is_file():
        raise FileNotFoundError(f"Missing source Figure 7 run config: {run_config}")
    _, config_rows = read_tsv(run_config)
    panel_set_rows = [row for row in config_rows if row.get("key") == "panel_set"]
    if len(panel_set_rows) != 1 or panel_set_rows[0].get("value") != expected_panel_set:
        raise ValueError(
            f"Source Figure 7 run must explicitly record panel_set={expected_panel_set}"
        )

    panel_contract = run_root / "metadata" / "panel_contract.tsv"
    if not panel_contract.is_file():
        raise FileNotFoundError(f"Missing source Figure 7 panel contract: {panel_contract}")
    _, contract_rows = read_tsv(panel_contract)
    expected_contract = [
        (str(spec["panel"]), Path(str(spec["source"])).name)
        for spec in selected_specs
        if str(spec["module"]) == module
        and (not spec.get("optional") or (run_root / str(spec["source"])).is_file())
    ]
    observed_contract = [(row.get("panel_id", ""), row.get("filename", "")) for row in contract_rows]
    if observed_contract != expected_contract:
        raise ValueError(
            f"Source Figure 7 panel contract mismatch: expected={expected_contract}; "
            f"observed={observed_contract}"
        )
    observed_figures = {
        path.resolve()
        for path in run_root.rglob("*")
        if path.is_file() and path.suffix.lower() in FIGURE_SUFFIXES
    }
    if observed_figures != expected_sources:
        missing = sorted(str(path) for path in expected_sources - observed_figures)
        unexpected = sorted(str(path) for path in observed_figures - expected_sources)
        raise ValueError(
            f"{module} violates its exact panel inventory: missing={missing}; unexpected={unexpected}"
        )

    _, rows = read_tsv(output_manifest)
    rows_by_path: dict[Path, list[dict[str, str]]] = {}
    for row in rows:
        path = first_local_path(row, repo_root, ("repo_relative_path", "absolute_path", "path"))
        if path is not None:
            rows_by_path.setdefault(path.resolve(), []).append(row)
    for source in sorted(expected_sources):
        matches = rows_by_path.get(source, [])
        if len(matches) != 1:
            raise ValueError(f"Expected one output-manifest row for {source}; found {len(matches)}")
        row = matches[0]
        if row.get("role") != "output_figure" or row.get("source_kind") != "generated_panel":
            raise ValueError(f"Output-manifest row is not a generated figure: {source}")
        if row.get("module") != module or row.get("command_id") != source_run_id:
            raise ValueError(f"Output-manifest provenance mismatch for {source}")
        if row.get("sha256", "").strip() != sha256_file(source):
            raise ValueError(f"Output-manifest checksum mismatch for {source}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-root", type=Path, default=Path("figures"))
    run_id_group = parser.add_mutually_exclusive_group()
    run_id_group.add_argument("--source-run-id")
    run_id_group.add_argument(
        "--run-id",
        dest="legacy_run_id",
        help="Backward-compatible alias for --source-run-id.",
    )
    parser.add_argument("--operation-id", default="")
    parser.add_argument("--module-run", action="append", default=[], help="Module run path as NAME=PATH. May be repeated.")
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    source_run_id = args.source_run_id or args.legacy_run_id
    if not source_run_id:
        parser.error("one of --source-run-id or --run-id is required")

    repo_root = args.repo_root.resolve() if args.repo_root else repo_root_from(Path.cwd())
    figure_root = args.figure_root
    if not figure_root.is_absolute():
        figure_root = repo_root / figure_root
    figure_root = figure_root.resolve()
    module_runs = parse_module_run(args.module_run, repo_root=repo_root)
    operation_id = args.operation_id or source_run_id
    validate_run_id(source_run_id)
    validate_run_id(operation_id)
    selected_specs = [spec for spec in PANEL_SPECS if str(spec["module"]) in module_runs]
    duplicate_specs = [
        key
        for key, count in Counter(
            (str(spec["figure"]), str(spec["panel"])) for spec in selected_specs
        ).items()
        if count > 1
    ]
    if duplicate_specs:
        raise ValueError(f"Duplicate selected panel contract(s): {duplicate_specs}")
    for module, run_root in module_runs.items():
        if not run_root.is_dir():
            raise FileNotFoundError(f"Missing module run directory for {module}: {run_root}")
        validate_strict_source_run(
            module, run_root, selected_specs, repo_root, source_run_id
        )

    rows_by_figure: dict[str, list[dict[str, str]]] = {}
    expected_by_figure: dict[str, set[str]] = {}
    for spec in selected_specs:
        module = str(spec["module"])
        source = module_runs[module] / str(spec["source"])
        if not source.exists():
            if spec.get("optional"):
                continue
            raise FileNotFoundError(f"Missing source for {spec['figure']} {spec['panel']}: {source}")
        expected_by_figure.setdefault(str(spec["figure"]), set()).add(str(spec["panel"]))
        out_dir = figure_root / str(spec["figure"])
        out_dir.mkdir(parents=True, exist_ok=True)
        asset = out_dir / str(spec["asset"])
        if asset.exists() and not args.overwrite:
            raise FileExistsError(f"Asset exists; use --overwrite to replace: {asset}")
        shutil.copy2(source, asset)
        if sha256_file(asset) != sha256_file(source):
            raise OSError(f"Copied asset checksum does not match source: {asset}")
        rows_by_figure.setdefault(str(spec["figure"]), []).append(
            generated_row(
                spec,
                source=source,
                asset=asset,
                result_run_root=module_runs[module],
                repo_root=repo_root,
                source_run_id=source_run_id,
                operation_id=operation_id,
            )
        )

    touched_figures = set(rows_by_figure)
    for row in EXTERNAL_ROWS:
        if row["figure"] in touched_figures:
            rows_by_figure.setdefault(row["figure"], []).append(external_row(row))

    for figure, rows in sorted(rows_by_figure.items()):
        rows = sorted(rows, key=lambda row: row.get("panel", ""))
        manifest_path = figure_root / figure / "manifest.tsv"
        panel_errors = validate_expected_panel_set(
            rows, expected_by_figure.get(figure, set()), manifest_path
        )
        if panel_errors:
            raise ValueError("Invalid materialized panel set:\n" + "\n".join(panel_errors))
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(manifest_path, rows, FIGURE_MANIFEST_COLUMNS)
        print(f"Wrote {manifest_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
