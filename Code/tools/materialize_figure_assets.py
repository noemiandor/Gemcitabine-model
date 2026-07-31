#!/usr/bin/env python3
"""Copy selected source panels into manuscript-facing figures/ folders."""

from __future__ import annotations

import argparse
import os
import shutil
from collections import Counter
from pathlib import Path

from figure_output_contract import (
    FIGURE_MANIFEST_COLUMNS,
    GENERATED_SOURCE_KINDS,
    first_local_path,
    module_manifest_local_path,
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
        "module": "si_figures",
        "source": "figures/panel_SuppFig4_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4",
        "asset": "panel_SuppFig4_composite.pdf",
        "caption_role": "Supplementary Figure 4 final composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig4_png",
        "asset": "panel_SuppFig4_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 4 composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig5_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig5",
        "asset": "panel_SuppFig5_composite.pdf",
        "caption_role": "Supplementary Figure 5 final composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig5_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig5_png",
        "asset": "panel_SuppFig5_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 5 composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig6_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig6",
        "asset": "panel_SuppFig6_composite.pdf",
        "caption_role": "Supplementary Figure 6 final composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig6_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig6_png",
        "asset": "panel_SuppFig6_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 6 composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig7_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig7",
        "asset": "panel_SuppFig7_composite.pdf",
        "caption_role": "Supplementary Figure 7 final composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig7_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig7_png",
        "asset": "panel_SuppFig7_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 7 composite",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7A_day17_tgi_calculation.pdf",
        "figure": "Figure7",
        "panel": "7A",
        "asset": "panel_7A_day17_tgi_calculation.pdf",
        "caption_role": "Tumor-growth trajectories explaining the Day-17 TGI calculation",
        "variant": "pdf",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7B_cellcycle_selected_ecdf_comparisons.pdf",
        "figure": "Figure7",
        "panel": "7B",
        "asset": "panel_7B_cellcycle_selected_ecdf_comparisons.pdf",
        "caption_role": "Selected CellCycle mean-ECDF comparisons",
        "variant": "pdf",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7C_day17_tgi_by_initial_ploidy.pdf",
        "figure": "Figure7",
        "panel": "7C",
        "asset": "panel_7C_day17_tgi_by_initial_ploidy.pdf",
        "caption_role": "Day-17 TGI in initial 2N versus 4N treated tumors",
        "variant": "pdf",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf",
        "figure": "Figure7",
        "panel": "7D",
        "asset": "panel_7D_day17_tgi_vs_centered_ecdf_shift.pdf",
        "caption_role": "Within-dose-centered TGI and CellCycle ECDF-shift association",
        "variant": "pdf",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7E_day17_tgi_vs_mean_etp.pdf",
        "figure": "Figure7",
        "panel": "7E",
        "asset": "panel_7E_day17_tgi_vs_mean_etp.pdf",
        "caption_role": "Day-17 TGI versus sample mean endpoint ploidy",
        "variant": "pdf",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7F_pseudotime_state_pathway_activity.pdf",
        "figure": "Figure7",
        "panel": "7F",
        "asset": "panel_7F_pseudotime_state_pathway_activity.pdf",
        "caption_role": "Pathway activity across the accumulated CellCycle pseudotime state",
        "variant": "pdf",
        "optional": True,
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7A_day17_tgi_calculation.png",
        "figure": "Figure7",
        "panel": "7A_png",
        "asset": "panel_7A_day17_tgi_calculation.png",
        "caption_role": "PNG derivative of the Day-17 TGI calculation panel",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7B_cellcycle_selected_ecdf_comparisons.png",
        "figure": "Figure7",
        "panel": "7B_png",
        "asset": "panel_7B_cellcycle_selected_ecdf_comparisons.png",
        "caption_role": "PNG derivative of the selected CellCycle mean-ECDF panel",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7C_day17_tgi_by_initial_ploidy.png",
        "figure": "Figure7",
        "panel": "7C_png",
        "asset": "panel_7C_day17_tgi_by_initial_ploidy.png",
        "caption_role": "PNG derivative of the initial-ploidy Day-17 TGI panel",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7D_day17_tgi_vs_centered_ecdf_shift.png",
        "figure": "Figure7",
        "panel": "7D_png",
        "asset": "panel_7D_day17_tgi_vs_centered_ecdf_shift.png",
        "caption_role": "PNG derivative of the centered ECDF-shift association panel",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7E_day17_tgi_vs_mean_etp.png",
        "figure": "Figure7",
        "panel": "7E_png",
        "asset": "panel_7E_day17_tgi_vs_mean_etp.png",
        "caption_role": "PNG derivative of the endpoint-ploidy association panel",
        "variant": "png",
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/panel_7F_pseudotime_state_pathway_activity.png",
        "figure": "Figure7",
        "panel": "7F_png",
        "asset": "panel_7F_pseudotime_state_pathway_activity.png",
        "caption_role": "PNG derivative of the pseudotime state-pathway activity panel",
        "variant": "png",
        "optional": True,
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/Figure7_reviewed_GRCh.png",
        "figure": "Figure7",
        "panel": "7A-7K_composite",
        "asset": "Figure7_reviewed_GRCh.png",
        "caption_role": (
            "Main Figure 7 A-K composite assembled in first-citation order"
        ),
        "variant": "png",
        "optional": True,
        "contract": True,
    },
]

STRICT_FIGURE_MODULES = {"in_vivo_figure7", "si_figures"}
FIGURE_SUFFIXES = {".pdf", ".png", ".jpg", ".jpeg", ".svg", ".tif", ".tiff"}
FIGURE7_REVIEWED_REFERENCE_ID = (
    "state_pathway_grch_human_only_etp2_24_day17_v2"
)
FIGURE7_REVIEWED_REFERENCE_KIND = "reviewed_human_only_frozen"
FIGURE7_REVIEWED_REFERENCE_ROOT = (
    Path("Data/in-vivo/figure7/saved_state_pathway")
    / FIGURE7_REVIEWED_REFERENCE_ID
)
FIGURE7_REVIEWED_FILES = {
    "panel_7F_pathway_activity_plot_data.tsv": "c2a6d455fd11c1d95bbc31441419953d4cc0d88eb8acf611fc3488f04af8d0f8",
    "panel_7F_selected_pathway_gsea.tsv": "ca86833bb6159cfde1e1218d897183809eba00a0cf096fc5c4951b5065767460",
    "panel_7F_leading_edge_genes.tsv": "6e86e697b64565cf0281e8684d29001eed84c83e1a3b9e45f7ffd37b78fdf9dd",
    "state_pathway_gene_ranking_complete.tsv": "ac9ced0c9bb691b490b0197c9962a72a6d272781e6883631c76d7c9f2397d34c",
    "state_pathway_gsea_complete.tsv": "4d1283e89b28e5100586a591621b5522b1698859d5d64b69054c4c470356a8d6",
    "state_pathway_sample_bin_coverage.tsv": "a0a408a91968f172ea7a30d9432c63571f74cea5fe27f4c5ad1c8f6044172a47",
    "state_pathway_design_qc.tsv": "66a1cf1a5362539f50777ae482b37b5b30e4d92bd673eee9aa4c842e91160151",
    "state_pathway_provenance.tsv": "406c4fa97b96dc001b1744e01658a533fc571897724a485080e8bbbc787293a4",
}
SI7_REVIEWED_FEATURE_POLICY = (
    "Human tumor/cell-line analysis: retain exact GRCh38-prefixed features "
    "and exclude exact GRCm39-prefixed features before RNA normalization, "
    "differential expression, symbol cleanup, deduplication, ORA, and GSEA; "
    "reject unrecognized feature prefixes."
)
SI7_REVIEWED_GENE_SET_DATABASE = (
    "MSigDB 2026.1.Hs Hallmark (Homo sapiens symbols)"
)
SI7_REVIEWED_FROZEN_MATRIX_NOTE = (
    "Reviewed SI7 matrices are the exact outputs from raw-refit run "
    "grch_human_only_v2_20260729_raw_refit_retry3_si_figures: exact "
    "GRCh38-prefixed RNA counts were retained before fresh RNA normalization, "
    "differential expression, symbol cleanup, deduplication, ORA, and GSEA."
)
SI_REVIEWED_MANIFEST_SHA256 = (
    "b624c3f3ff945c51f09b9e6e512a97df57eb4e514b3fba28a65e97a38207f135"
)
FIGURE7_PROCESSED_INPUTS = {
    Path(
        "Data/in-vivo/figure7/processed/"
        "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ): "bd0c6fcf3f6691114445a7c8a9fb8de7e55ef6580a9737284cd38f5bba216186",
    Path(
        "Data/in-vivo/figure7/processed/"
        "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
    ): "1e44794fb2b2c402b49cf9d17abcc54cd7690dc5d1d611e38a6aa3b996093d3c",
}

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


def panel_specs_for_figure7_variant(
    tgi_day: int,
    figure_name: str,
) -> list[dict[str, object]]:
    if tgi_day < 0:
        raise ValueError("--figure7-tgi-day must be a non-negative integer")
    if not (
        figure_name == "Figure7"
        or (
            figure_name.startswith(("Figure7_", "Figure7-", "Figure7."))
            and all(char.isalnum() or char in "._-" for char in figure_name)
        )
    ):
        raise ValueError(
            "--figure7-figure-name must be Figure7 or a Figure7-prefixed folder name"
        )
    specs: list[dict[str, object]] = []
    for original in PANEL_SPECS:
        spec = dict(original)
        if str(spec["module"]) == "in_vivo_figure7":
            spec["figure"] = figure_name
            for key in ("source", "asset", "caption_role"):
                spec[key] = str(spec[key]).replace("day17", f"day{tgi_day}").replace(
                    "Day-17", f"Day-{tgi_day}"
                )
        specs.append(spec)
    return specs


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


def manifest_relative(path: Path, manifest_parent: Path) -> str:
    relative = Path(os.path.relpath(path.resolve(), manifest_parent.resolve()))
    return f"manifest:{relative.as_posix()}"


def portable_locator(path: Path, repo_root: Path, manifest_parent: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return manifest_relative(path, manifest_parent)


def generated_row(
    spec: dict[str, object],
    source: Path,
    asset: Path,
    result_run_root: Path,
    repo_root: Path,
    source_run_id: str,
    operation_id: str,
) -> dict[str, str]:
    manifest_parent = asset.parent
    published_asset = portable_locator(asset, repo_root, manifest_parent)
    source_locator = portable_locator(source, repo_root, manifest_parent)
    source_run_locator = portable_locator(
        result_run_root,
        repo_root,
        manifest_parent,
    )
    return {
        "figure": str(spec["figure"]),
        "panel": str(spec["panel"]),
        "asset_path": published_asset,
        "source_file": published_asset,
        "source_kind": "generated_panel",
        "generated_by": "Manager.sh",
        "command": "materialize_figure_assets.py",
        "input_data": published_asset,
        "result_run_dir": portable_locator(
            manifest_parent,
            repo_root,
            manifest_parent,
        ),
        "run_id": source_run_id,
        "caption_role": str(spec["caption_role"]),
        "asset_status": "generated",
        "not_regenerated_reason": "",
        "local_provenance_path": "",
        "citation_or_uri": "",
        "notes": ";".join(
            (
                f"materialization_operation_id={operation_id}",
                f"validated_source_file={source_locator}",
                f"validated_source_run={source_run_locator}",
                f"validated_source_sha256={sha256_file(source)}",
            )
        ),
    }


def external_row(row: dict[str, str]) -> dict[str, str]:
    out = {column: "" for column in FIGURE_MANIFEST_COLUMNS}
    out.update(row)
    out["asset_status"] = out.get("source_kind", "external")
    return out


def read_unique_key_values(path: Path, label: str) -> dict[str, str]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    _, rows = read_tsv(path)
    keys = [row.get("key", "") for row in rows]
    values = [row.get("value", "") for row in rows]
    if (
        not rows
        or any(not key or not value for key, value in zip(keys, values))
        or len(set(keys)) != len(keys)
    ):
        raise ValueError(f"{label} is malformed")
    return dict(zip(keys, values))


def validate_si_publication_contract(run_root: Path, repo_root: Path) -> None:
    provenance = read_unique_key_values(
        run_root / "metadata" / "si_figures_provenance.tsv",
        "source SI Figures provenance",
    )
    run_config = read_unique_key_values(
        run_root / "metadata" / "run_config.tsv",
        "source SI Figures run config",
    )
    canonical_manifest = repo_root / "Data" / "in-vivo" / "SIfigures" / "manifest.tsv"
    if not canonical_manifest.is_file():
        raise FileNotFoundError(
            f"Missing reviewed SI Figures table manifest: {canonical_manifest}"
        )
    if sha256_file(canonical_manifest) != SI_REVIEWED_MANIFEST_SHA256:
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the repository "
            "does not contain the exact reviewed 11-table manifest"
        )
    manifest_headers, manifest_rows = read_tsv(canonical_manifest)
    if manifest_headers != [
        "filename",
        "bytes",
        "sha256",
        "source_revision",
        "notes",
    ]:
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the reviewed "
            "cache manifest schema is invalid"
        )
    reviewed_names = [row.get("filename", "") for row in manifest_rows]
    if (
        len(reviewed_names) != 11
        or len(set(reviewed_names)) != 11
        or any(
            not name
            or Path(name).name != name
            or name == "manifest.tsv"
            for name in reviewed_names
        )
    ):
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the reviewed "
            "cache manifest must name exactly 11 safe tables"
        )
    run_tables = run_root / "tables"
    run_manifest = run_tables / "manifest.tsv"
    observed_run_names = sorted(
        path.name for path in run_tables.iterdir() if path.is_file()
    ) if run_tables.is_dir() else []
    expected_run_names = sorted(["manifest.tsv", *reviewed_names])
    if observed_run_names != expected_run_names:
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the source "
            "run does not contain the exact reviewed 11-table cache inventory"
        )
    if sha256_file(run_manifest) != SI_REVIEWED_MANIFEST_SHA256:
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the source "
            "run table manifest is not the reviewed manifest"
        )
    canonical_cache = canonical_manifest.parent
    for row in manifest_rows:
        filename = row["filename"]
        canonical_table = canonical_cache / filename
        run_table = run_tables / filename
        if (
            not canonical_table.is_file()
            or not run_table.is_file()
            or row.get("bytes") != str(canonical_table.stat().st_size)
            or row.get("sha256") != sha256_file(canonical_table)
            or sha256_file(run_table) != row.get("sha256")
            or run_table.stat().st_size != canonical_table.stat().st_size
        ):
            raise ValueError(
                "Canonical SI Figures materialization is prohibited: source "
                f"run table {filename} is not the exact reviewed cache table"
            )
    si7_expected = {
        "si7_canonical_publication_allowed": "true",
        "si7_feature_species_policy": SI7_REVIEWED_FEATURE_POLICY,
        "si7_gene_set_database": SI7_REVIEWED_GENE_SET_DATABASE,
    }
    composition_expected = {
        "composition_normalization": (
            "within-sample cluster proportions averaged with equal sample "
            "weights within group"
        ),
        "composition_enrichment_test": (
            "exact independent-sample label permutation; one group versus "
            "exchangeable remaining samples"
        ),
        "composition_permutation_strata": (
            "SI4E=initial_ploidy;SI4G=context;SI5F=descriptive;"
            "SI5G=initial_ploidy;SI5H=initial_ploidy;SI5I=dose"
        ),
        "composition_multiple_testing": (
            "Benjamini-Hochberg across all group-by-cluster contrasts "
            "within each panel"
        ),
        "composition_fdr_threshold": "0.05",
    }
    renderer = (
        repo_root
        / "Code"
        / "in-vivo"
        / "SI_figures"
        / "generate_supplementary_figures.R"
    )
    if not renderer.is_file():
        raise FileNotFoundError(f"Missing SI Figures renderer: {renderer}")
    composition_helper = (
        repo_root
        / "Code"
        / "in-vivo"
        / "SI_figures"
        / "normalized_composition.R"
    )
    if not composition_helper.is_file():
        raise FileNotFoundError(
            f"Missing SI Figures normalized-composition helper: {composition_helper}"
        )
    shared_context_helper = (
        repo_root
        / "Code"
        / "in-vivo"
        / "SI_figures"
        / "shared_context_panels.R"
    )
    if not shared_context_helper.is_file():
        raise FileNotFoundError(
            f"Missing SI Figures shared-context helper: {shared_context_helper}"
        )
    if (
        any(
            provenance.get(key) != value
            for key, value in si7_expected.items()
        )
        or any(
            run_config.get(key) != value
            for key, value in {**si7_expected, **composition_expected}.items()
        )
        or provenance.get("si7_frozen_matrix_note")
        != SI7_REVIEWED_FROZEN_MATRIX_NOTE
        or run_config.get("module") != "si_figures"
        or run_config.get("figures") != "4,5,6,7"
        or run_config.get("figure_file_count") != "8"
        or run_config.get("cache_file_count") != "11"
        or run_config.get("table_mode") != "frozen_plot_tables_only"
        or provenance.get("artifact") != "supplementary_figures_4_7"
        or provenance.get("entrypoint") != (
            "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        )
        or provenance.get("entrypoint_sha256") != sha256_file(renderer)
        or provenance.get("normalized_composition_helper") != (
            "Code/in-vivo/SI_figures/normalized_composition.R"
        )
        or provenance.get("normalized_composition_helper_sha256")
        != sha256_file(composition_helper)
        or provenance.get("shared_context_panels_helper") != (
            "Code/in-vivo/SI_figures/shared_context_panels.R"
        )
        or provenance.get("shared_context_panels_helper_sha256")
        != sha256_file(shared_context_helper)
        or provenance.get("table_cache_manifest_sha256")
        != sha256_file(canonical_manifest)
    ):
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the run "
            "does not carry the reviewed frozen-table publication contract"
        )


def validate_strict_source_run(
    module: str,
    run_root: Path,
    selected_specs: list[dict[str, object]],
    repo_root: Path,
    source_run_id: str,
    figure7_tgi_day: int,
) -> None:
    if module not in STRICT_FIGURE_MODULES:
        return
    if module == "si_figures":
        validate_si_publication_contract(run_root, repo_root)
    expected_name = (
        f"{source_run_id}_figure7"
        if module == "in_vivo_figure7"
        else f"{source_run_id}_si_figures"
    )
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
    if module == "si_figures":
        cache_root = repo_root / "Data" / "in-vivo" / "SIfigures"
        cache_manifest = cache_root / "manifest.tsv"
        cache_headers, cache_rows = read_tsv(cache_manifest)
        cache_names = [row.get("filename", "") for row in cache_rows]
        if (
            cache_headers
            != ["filename", "bytes", "sha256", "source_revision", "notes"]
            or sha256_file(cache_manifest) != SI_REVIEWED_MANIFEST_SHA256
            or
            len(cache_names) != 11
            or len(set(cache_names)) != 11
            or any(
                not name
                or Path(name).name != name
                or name == "manifest.tsv"
                for name in cache_names
            )
        ):
            raise ValueError(
                "Reviewed SI Figures cache manifest must list exactly 11 "
                "safe plot-facing tables"
            )
        for row in cache_rows:
            table = cache_root / row["filename"]
            if not table.is_file():
                raise FileNotFoundError(
                    f"Missing reviewed SI Figures cache table: {table}"
                )
            if (
                row.get("bytes") != str(table.stat().st_size)
                or row.get("sha256") != sha256_file(table)
            ):
                raise ValueError(
                    "Reviewed SI Figures cache manifest does not match table "
                    f"bytes: {row['filename']}"
                )
        required_inputs = {
            (
                repo_root
                / "Code"
                / "in-vivo"
                / "figure7"
                / "figure7_config.yaml"
            ).resolve(),
            (
                repo_root
                / "Code"
                / "in-vivo"
                / "SI_figures"
                / "run_supplementary_figures.R"
            ).resolve(),
            (
                repo_root
                / "Code"
                / "in-vivo"
                / "SI_figures"
                / "generate_supplementary_figures.R"
            ).resolve(),
            (
                repo_root
                / "Code"
                / "in-vivo"
                / "SI_figures"
                / "shared_context_panels.R"
            ).resolve(),
            (
                repo_root
                / "Code"
                / "in-vivo"
                / "SI_figures"
                / "normalized_composition.R"
            ).resolve(),
            (
                repo_root
                / "Code"
                / "tools"
                / "validate_si_figures_table_cache.py"
            ).resolve(),
            cache_manifest.resolve(),
            *{(cache_root / name).resolve() for name in cache_names},
        }
        observed_inputs = {
            path.resolve()
            for row in input_rows
            if (
                path := first_local_path(
                    row,
                    repo_root,
                    ("repo_relative_path", "absolute_path", "path"),
                )
            )
            is not None
        }
        missing_inputs = sorted(
            str(path) for path in required_inputs - observed_inputs
        )
        if missing_inputs:
            raise ValueError(
                "SI Figures input manifest does not bind the reviewed "
                f"11-table cache: missing={missing_inputs}"
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
    if module == "in_vivo_figure7":
        panel_f_paths = {
            (run_root / str(spec["source"])).resolve()
            for spec in selected_specs
            if str(spec["module"]) == module
            and str(spec["panel"]).startswith("7F")
        }
        present_panel_f_paths = {path for path in panel_f_paths if path.is_file()}
        if present_panel_f_paths and present_panel_f_paths != panel_f_paths:
            raise ValueError(
                "Source Figure 7 run must contain both PDF and PNG panel-F "
                "assets or neither"
            )
        has_panel_f = present_panel_f_paths == panel_f_paths
        composite_specs = [
            spec
            for spec in selected_specs
            if str(spec["module"]) == module
            and str(spec["panel"]) == "7A-7K_composite"
        ]
        if len(composite_specs) != 1:
            raise ValueError("Figure 7 must define one A-K composite contract")
        composite_path = (
            run_root / str(composite_specs[0]["source"])
        ).resolve()
        if composite_path.is_file() != has_panel_f:
            requirement = "present" if has_panel_f else "absent"
            raise ValueError(
                "Source Figure 7 A-K composite must be "
                f"{requirement} exactly when panel F is present"
            )
        expected_panel_set = "a-f" if has_panel_f else "a-e"
        run_config = read_unique_key_values(
            run_root / "metadata" / "run_config.tsv",
            "source Figure 7 run config",
        )
        if run_config.get("panel_set") != expected_panel_set:
            raise ValueError(
                "Source Figure 7 run must explicitly record "
                f"panel_set={expected_panel_set}"
            )
        if run_config.get("tgi_day") != str(figure7_tgi_day):
            raise ValueError(
                "Source Figure 7 run must explicitly record "
                f"tgi_day={figure7_tgi_day}"
            )
        observed_input_paths = {
            path.resolve()
            for row in input_rows
            if (
                path := module_manifest_local_path(
                    row,
                    repo_root,
                )
            ) is not None
        }
        for relative_path, expected_hash in FIGURE7_PROCESSED_INPUTS.items():
            processed_path = (repo_root / relative_path).resolve()
            if processed_path not in observed_input_paths:
                raise ValueError(
                    "Figure 7 source input manifest does not bind the exact "
                    f"processed input: {relative_path}"
                )
            if (
                not processed_path.is_file()
                or sha256_file(processed_path) != expected_hash
            ):
                raise ValueError(
                    "Figure 7 processed input differs from its reviewed hash: "
                    f"{relative_path}"
                )
        if has_panel_f:
            figure7_config = (
                repo_root
                / "Code/in-vivo/figure7/figure7_config.yaml"
            )
            if (
                not figure7_config.is_file()
                or run_config.get("config_sha256")
                != sha256_file(figure7_config)
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "the run does not bind the active reviewed config"
                )
            reviewed_identity = {
                "state_pathway_reference_id": FIGURE7_REVIEWED_REFERENCE_ID,
                "state_pathway_reference_kind": FIGURE7_REVIEWED_REFERENCE_KIND,
                "canonical_publication_allowed": "true",
            }
            if any(
                run_config.get(key) != value
                for key, value in reviewed_identity.items()
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "run metadata does not carry the exact reviewed panel-7F "
                    "publication identity"
                )
            for filename, expected_hash in FIGURE7_REVIEWED_FILES.items():
                source_file = (
                    run_root / "metadata" / filename
                    if filename == "state_pathway_provenance.tsv"
                    else run_root / "tables" / filename
                )
                if not source_file.is_file():
                    raise FileNotFoundError(
                        f"Missing reviewed panel-7F contract file: {source_file}"
                    )
                if sha256_file(source_file) != expected_hash:
                    raise ValueError(
                        "Canonical Figure 7 materialization is prohibited: "
                        f"{filename} is not the reviewed "
                        f"{FIGURE7_REVIEWED_REFERENCE_ID} artifact"
                    )
            required_input_paths = {
                (
                    repo_root
                    / "Code/in-vivo/figure7/run_figure7.R"
                ).resolve(),
                (
                    repo_root
                    / "Code/in-vivo/figure7/figure7_config.yaml"
                ).resolve(),
                *{
                    (
                        repo_root
                        / FIGURE7_REVIEWED_REFERENCE_ROOT
                        / filename
                    ).resolve()
                    for filename in FIGURE7_REVIEWED_FILES
                },
            }
            missing_input_paths = sorted(
                str(path)
                for path in required_input_paths - observed_input_paths
            )
            if missing_input_paths:
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "the source input manifest does not bind the reviewed "
                    "renderer, config, and eight-file panel-7F reference: "
                    f"missing={missing_input_paths}"
                )
            provenance = read_unique_key_values(
                run_root / "metadata" / "state_pathway_provenance.tsv",
                "source Figure 7 state-pathway provenance",
            )
            if (
                provenance.get("canonical_reference_id")
                != FIGURE7_REVIEWED_REFERENCE_ID
                or provenance.get("reference_kind")
                != FIGURE7_REVIEWED_REFERENCE_KIND
                or provenance.get("canonical_publication_allowed")
                != "true"
                or provenance.get("reviewed_source_provenance_sha256")
                != (
                    "c325359b1d0fe67c3ad1f13523e993cb80e5168d164d06e1"
                    "b2fa1b1b584be888"
                )
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "reviewed provenance lineage is invalid"
                )
            expected_composite = {
                "main_composite_panel_set": "a-k",
                "main_composite_filename": "Figure7_reviewed_GRCh.png",
                "main_composite_panel_order": (
                    "A=7A;B=7C;C=SI4A;D=SI4B;E=SI4C;F=SI4E;"
                    "G=SI7B;H=7B;I=7F;J=7D;K=7E"
                ),
                "reviewed_si_cache_manifest": (
                    "Data/in-vivo/SIfigures/manifest.tsv"
                ),
                "reviewed_si_cache_manifest_sha256": (
                    SI_REVIEWED_MANIFEST_SHA256
                ),
                "si_context_cache_policy": "reviewed",
                "si_context_cache_kind": "reviewed_human_only_frozen",
                "si_context_cache_manifest": (
                    "Data/in-vivo/SIfigures/manifest.tsv"
                ),
                "si_context_cache_manifest_sha256": (
                    SI_REVIEWED_MANIFEST_SHA256
                ),
                "si_context_cache_canonical_publication_allowed": "true",
            }
            if any(
                run_config.get(key) != value
                for key, value in expected_composite.items()
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "the A-K composite/cache contract is invalid"
                )

            cache_root = repo_root / "Data/in-vivo/SIfigures"
            cache_manifest = cache_root / "manifest.tsv"
            cache_headers, cache_rows = read_tsv(cache_manifest)
            cache_names = [row.get("filename", "") for row in cache_rows]
            if (
                cache_headers
                != [
                    "filename",
                    "bytes",
                    "sha256",
                    "source_revision",
                    "notes",
                ]
                or sha256_file(cache_manifest)
                != SI_REVIEWED_MANIFEST_SHA256
                or len(cache_names) != 11
                or len(set(cache_names)) != 11
            ):
                raise ValueError(
                    "Canonical Figure 7 requires the exact reviewed SI cache"
                )
            for row in cache_rows:
                table = cache_root / row["filename"]
                if (
                    not table.is_file()
                    or row.get("bytes") != str(table.stat().st_size)
                    or row.get("sha256") != sha256_file(table)
                ):
                    raise ValueError(
                        "Canonical Figure 7 reviewed SI cache mismatch: "
                        f"{row['filename']}"
                    )
            observed_inputs = {
                path.resolve()
                for row in input_rows
                if (
                    path := first_local_path(
                        row,
                        repo_root,
                        ("repo_relative_path", "absolute_path", "path"),
                    )
                )
                is not None
            }
            required_context_inputs = {
                (repo_root / "Code/in-vivo/figure7/run_figure7.R").resolve(),
                (
                    repo_root
                    / "Code/in-vivo/figure7/src/context_panels.R"
                ).resolve(),
                (
                    repo_root
                    / "Code/in-vivo/SI_figures/shared_context_panels.R"
                ).resolve(),
                (
                    repo_root
                    / "Code/in-vivo/SI_figures/normalized_composition.R"
                ).resolve(),
                (
                    repo_root
                    / "Code/tools/validate_si_figures_table_cache.py"
                ).resolve(),
                cache_manifest.resolve(),
                *{(cache_root / name).resolve() for name in cache_names},
            }
            missing_context_inputs = sorted(
                str(path)
                for path in required_context_inputs - observed_inputs
            )
            if missing_context_inputs:
                raise ValueError(
                    "Canonical Figure 7 input manifest does not bind the "
                    "reviewed SI context inputs: "
                    f"missing={missing_context_inputs}"
                )

    panel_contract = run_root / "metadata" / "panel_contract.tsv"
    if not panel_contract.is_file():
        raise FileNotFoundError(f"Missing source panel contract: {panel_contract}")
    _, contract_rows = read_tsv(panel_contract)
    expected_contract = [
        (str(spec["panel"]), Path(str(spec["source"])).name)
        for spec in selected_specs
        if str(spec["module"]) == module
        and (
            module == "si_figures"
            or spec.get("contract", spec.get("variant", "pdf") == "pdf")
        )
        and (not spec.get("optional") or (run_root / str(spec["source"])).is_file())
    ]
    observed_contract = [(row.get("panel_id", ""), row.get("filename", "")) for row in contract_rows]
    if observed_contract != expected_contract:
        raise ValueError(
            f"Source panel contract mismatch: expected={expected_contract}; "
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
        path = module_manifest_local_path(
            row,
            repo_root,
            output_root=run_root,
        )
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
    parser.add_argument("--figure7-tgi-day", type=int, default=17)
    parser.add_argument("--figure7-figure-name", default="Figure7")
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
    panel_specs = panel_specs_for_figure7_variant(
        args.figure7_tgi_day, args.figure7_figure_name
    )
    selected_specs = [spec for spec in panel_specs if str(spec["module"]) in module_runs]
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
            module,
            run_root,
            selected_specs,
            repo_root,
            source_run_id,
            args.figure7_tgi_day,
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
    for figure in touched_figures:
        manifest_path = figure_root / figure / "manifest.tsv"
        if not manifest_path.is_file():
            continue
        headers, existing_rows = read_tsv(manifest_path)
        if headers != list(FIGURE_MANIFEST_COLUMNS):
            raise ValueError(
                f"Existing figure manifest has an invalid schema: {manifest_path}"
            )
        replaced_panels = set(expected_by_figure.get(figure, set()))
        replaced_panels.update(
            str(row["panel"])
            for row in EXTERNAL_ROWS
            if row["figure"] == figure
        )
        preserved_rows = [
            row
            for row in existing_rows
            if row.get("panel", "") not in replaced_panels
        ]
        rows_by_figure[figure].extend(preserved_rows)
        expected_by_figure.setdefault(figure, set()).update(
            row.get("panel", "")
            for row in preserved_rows
            if row.get("source_kind", "") in GENERATED_SOURCE_KINDS
            and row.get("panel", "")
        )
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
