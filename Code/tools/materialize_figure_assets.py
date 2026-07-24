#!/usr/bin/env python3
"""Copy selected source panels into manuscript-facing figures/ folders."""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import statistics
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
        "module": "si_figures",
        "source": "figures/panel_SuppFig4A_umap_cluster.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4A",
        "asset": "panel_SuppFig4A_umap_cluster.pdf",
        "caption_role": "In-vivo tumor-cell UMAP by cluster",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4A_umap_cluster.png",
        "figure": "Supplementary",
        "panel": "SuppFig4A_png",
        "asset": "panel_SuppFig4A_umap_cluster.png",
        "caption_role": "PNG derivative of the in-vivo tumor-cell UMAP by cluster",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4B_umap_initial_ploidy.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4B",
        "asset": "panel_SuppFig4B_umap_initial_ploidy.pdf",
        "caption_role": "In-vivo tumor-cell UMAP by initial ploidy",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4B_umap_initial_ploidy.png",
        "figure": "Supplementary",
        "panel": "SuppFig4B_png",
        "asset": "panel_SuppFig4B_umap_initial_ploidy.png",
        "caption_role": "PNG derivative of the in-vivo tumor-cell UMAP by initial ploidy",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4C_umap_treatment_dose.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4C",
        "asset": "panel_SuppFig4C_umap_treatment_dose.pdf",
        "caption_role": "In-vivo tumor-cell UMAP by treatment dose",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4C_umap_treatment_dose.png",
        "figure": "Supplementary",
        "panel": "SuppFig4C_png",
        "asset": "panel_SuppFig4C_umap_treatment_dose.png",
        "caption_role": "PNG derivative of the in-vivo tumor-cell UMAP by treatment dose",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4D_umap_mouse_facets.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4D",
        "asset": "panel_SuppFig4D_umap_mouse_facets.pdf",
        "caption_role": "In-vivo tumor-cell UMAP faceted by mouse",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4D_umap_mouse_facets.png",
        "figure": "Supplementary",
        "panel": "SuppFig4D_png",
        "asset": "panel_SuppFig4D_umap_mouse_facets.png",
        "caption_role": "PNG derivative of the in-vivo tumor-cell UMAP faceted by mouse",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4E_cluster_composition_by_mouse.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4E",
        "asset": "panel_SuppFig4E_cluster_composition_by_mouse.pdf",
        "caption_role": "Tumor-cluster composition by mouse",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4E_cluster_composition_by_mouse.png",
        "figure": "Supplementary",
        "panel": "SuppFig4E_png",
        "asset": "panel_SuppFig4E_cluster_composition_by_mouse.png",
        "caption_role": "PNG derivative of tumor-cluster composition by mouse",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4F_cluster_composition_by_ploidy_dose.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4F",
        "asset": "panel_SuppFig4F_cluster_composition_by_ploidy_dose.pdf",
        "caption_role": "Mouse-weighted cluster composition by initial ploidy and dose",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4F_cluster_composition_by_ploidy_dose.png",
        "figure": "Supplementary",
        "panel": "SuppFig4F_png",
        "asset": "panel_SuppFig4F_cluster_composition_by_ploidy_dose.png",
        "caption_role": "PNG derivative of mouse-weighted cluster composition by initial ploidy and dose",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig4_composite",
        "asset": "panel_SuppFig4_composite.pdf",
        "caption_role": "Supplementary Figures 4-7 composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig4_composite_png",
        "asset": "panel_SuppFig4_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figures 4-7 composite",
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
]

PANEL_SPECS = [spec for spec in PANEL_SPECS if spec["module"] != "si_figures"]
SI_FIGURE_PANEL_FILES = [
    ("SuppFig4A", "panel_SuppFig4A_umap_cluster", "All-cell UMAP by cluster"),
    ("SuppFig4B", "panel_SuppFig4B_umap_initial_ploidy", "All-cell UMAP by initial ploidy"),
    ("SuppFig4C", "panel_SuppFig4C_umap_context", "All-cell UMAP by Tumor/CellLine context"),
    ("SuppFig4D", "panel_SuppFig4D_umap_s_phase_score", "All-cell UMAP by S phase score"),
    ("SuppFig4E", "panel_SuppFig4E_cluster_context_proportion", "Tumor and CellLine proportions by cluster"),
    ("SuppFig4F", "panel_SuppFig4F_cluster_context_count", "Tumor and CellLine counts by cluster"),
    ("SuppFig4G", "panel_SuppFig4G_cluster_initial_ploidy_proportion", "Initial-ploidy proportions by cluster"),
    ("SuppFig4H", "panel_SuppFig4H_cluster_initial_ploidy_count", "Initial-ploidy counts by cluster"),
    ("SuppFig4_composite", "panel_SuppFig4_composite", "Supplementary Figure 4 composite"),
    ("SuppFig5A", "panel_SuppFig5A_umap_cluster", "Tumor UMAP by cluster"),
    ("SuppFig5B", "panel_SuppFig5B_umap_initial_ploidy", "Tumor UMAP by initial ploidy"),
    ("SuppFig5C", "panel_SuppFig5C_umap_treatment_dose", "Tumor UMAP by Gemcitabine dose"),
    ("SuppFig5D", "panel_SuppFig5D_umap_s_phase_score", "Tumor UMAP by S phase score"),
    ("SuppFig5E", "panel_SuppFig5E_umap_mouse_facets", "Mouse-faceted UMAP colored by initial ploidy"),
    ("SuppFig5F", "panel_SuppFig5F_cluster_composition_by_mouse", "Per-mouse cluster composition"),
    ("SuppFig5G", "panel_SuppFig5G_mouse_weighted_cluster_composition", "Mouse-weighted cluster composition"),
    ("SuppFig5H", "panel_SuppFig5H_cluster_composition_by_dose", "Dose composition by cluster"),
    ("SuppFig5I", "panel_SuppFig5I_cluster_composition_by_initial_ploidy", "Initial-ploidy composition by cluster"),
    ("SuppFig5_composite", "panel_SuppFig5_composite", "Supplementary Figure 5 composite"),
    ("SuppFig6A", "panel_SuppFig6A_umap_endpoint_ploidy_all_tumors", "All-tumor UMAP by endpoint ploidy"),
    ("SuppFig6B", "panel_SuppFig6B_umap_endpoint_ploidy_initial_2N", "Initial-2N tumor UMAP by endpoint ploidy"),
    ("SuppFig6C", "panel_SuppFig6C_umap_endpoint_ploidy_initial_4N", "Initial-4N tumor UMAP by endpoint ploidy"),
    ("SuppFig6D", "panel_SuppFig6D_umap_endpoint_ploidy_mouse_facets", "Mouse-faceted UMAP by endpoint ploidy"),
    ("SuppFig6_composite", "panel_SuppFig6_composite", "Supplementary Figure 6 composite"),
    ("SuppFig7A", "panel_SuppFig7A_cluster_Hallmark_ORA_annotation_score_heatmap_top20", "Cluster Hallmark ORA annotation-score heatmap"),
    ("SuppFig7B", "panel_SuppFig7B_cluster_Hallmark_GSEA_NES_heatmap_top20", "Cluster Hallmark GSEA NES heatmap"),
    ("SuppFig7_composite", "panel_SuppFig7_composite", "Supplementary Figure 7 composite"),
]
SI_FIGURE_CLUSTERS = ("0", "2", "4c", "5", "6", "8", "10", "13", "14")
SI_FIGURE_TABLE_CACHE_FILES = {
    "si_figures_cell_metadata.csv",
    "si_figures_cluster_key.tsv",
    "si_figure4_cluster_context_composition.csv",
    "si_figure4_cluster_initial_ploidy_composition.csv",
    "si_figure5_cluster_composition_by_mouse.csv",
    "si_figure5_cluster_dose_composition.csv",
    "si_figure5_cluster_initial_ploidy_composition.csv",
    "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
    "si_figure6_endpoint_ploidy_join_audit.csv",
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    "si_figure7_cluster_Hallmark_GSEA_all.csv",
    "si_figure7_cluster_Hallmark_ORA_all.csv",
    "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
    "si_figure7_hallmark_ORA_universe.csv",
}
SI_FIGURE_TABLE_CACHE_FILES.update(
    {
        filename
        for cluster in SI_FIGURE_CLUSTERS
        for filename in (
            f"si_figure7_cluster_{cluster}_top100_up_ORA_input.csv",
            f"si_figure7_cluster_{cluster}_vs_rest_DEG.csv",
        )
    }
)
for panel_id, filename_stub, caption_role in SI_FIGURE_PANEL_FILES:
    for extension in ("pdf", "png"):
        panel = panel_id if extension == "pdf" else f"{panel_id}_png"
        filename = f"{filename_stub}.{extension}"
        PANEL_SPECS.append(
            {
                "module": "si_figures",
                "source": f"figures/{filename}",
                "figure": "Supplementary",
                "panel": panel,
                "asset": filename,
                "caption_role": (
                    caption_role
                    if extension == "pdf"
                    else f"PNG derivative of {caption_role.lower()}"
                ),
                "variant": extension,
            }
        )

STRICT_FIGURE_MODULES = {"in_vivo_figure7", "si_figures"}
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


def read_delimited_rows(path: Path, delimiter: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing or empty SI Figures 4-7 contract file: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        rows = list(reader)
    if not fields or len(fields) != len(set(fields)) or not rows:
        raise ValueError(f"Invalid SI Figures 4-7 contract table: {path}")
    return fields, rows


def validate_si_figures_data_contract(
    run_root: Path,
    repo_root: Path,
    source_run_id: str,
    input_rows: list[dict[str, str]],
    output_rows_by_path: dict[Path, list[dict[str, str]]],
) -> None:
    required_artifacts = (
        run_root / "tables/si_figures_cell_metadata.csv",
        run_root / "tables/si_figures_cluster_key.tsv",
        run_root / "tables/si_figure4_cluster_context_composition.csv",
        run_root / "tables/si_figure4_cluster_initial_ploidy_composition.csv",
        run_root / "tables/si_figure5_cluster_composition_by_mouse.csv",
        run_root / "tables/si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
        run_root / "tables/si_figure5_cluster_dose_composition.csv",
        run_root / "tables/si_figure5_cluster_initial_ploidy_composition.csv",
        run_root / "tables/si_figure6_endpoint_ploidy_join_audit.csv",
        run_root / "tables/si_figure7_cluster_Hallmark_ORA_all.csv",
        run_root / "tables/si_figure7_cluster_Hallmark_GSEA_all.csv",
        run_root / "tables/si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
        run_root / "tables/si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
        run_root / "metadata/si_figures_provenance.tsv",
        run_root / "metadata/input_qc.tsv",
    )
    for path in required_artifacts:
        if not path.is_file() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"Missing SI Figures contract artifact: {path}")
        matches = output_rows_by_path.get(path.resolve(), [])
        if len(matches) != 1:
            raise ValueError(
                f"Expected one output-manifest row for SI Figures artifact {path}; "
                f"found {len(matches)}"
            )
        row = matches[0]
        if (
            row.get("module") != "si_figures"
            or row.get("command_id") != source_run_id
            or row.get("sha256", "").strip() != sha256_file(path)
        ):
            raise ValueError(f"SI Figures output provenance mismatch for {path}")

    canonical_path = run_root / "tables/si_figures_cell_metadata.csv"
    canonical_fields, canonical = read_delimited_rows(canonical_path, ",")
    required_canonical_fields = {
        "cell_id",
        "UMAP_1",
        "UMAP_2",
        "sample_id",
        "cluster_id",
        "initial_ploidy",
        "s_phase_score",
        "endpoint_ploidy",
        "endpoint_file",
        "endpoint_cell_id",
        "context",
        "included_in_si_figures",
    }
    if not required_canonical_fields <= set(canonical_fields):
        raise ValueError(
            "SI Figures canonical cell schema is incomplete: "
            f"{sorted(required_canonical_fields - set(canonical_fields))}"
        )
    cell_ids = [row["cell_id"].strip() for row in canonical]
    if any(not cell_id for cell_id in cell_ids) or len(cell_ids) != len(set(cell_ids)):
        raise ValueError("SI Figures canonical cell IDs must be nonempty and unique")
    tumor = [row for row in canonical if row["context"] == "Tumor"]
    cellline = [row for row in canonical if row["context"] == "CellLine"]
    if len(tumor) + len(cellline) != len(canonical):
        raise ValueError("SI Figures canonical context must be Tumor or CellLine")
    if any(not row["endpoint_ploidy"].strip() for row in tumor):
        raise ValueError("Every tumor cell must have endpoint ploidy")
    if any(row["endpoint_ploidy"].strip() for row in cellline):
        raise ValueError("CellLine cells must not have endpoint ploidy")

    endpoint_path = run_root / "tables/si_figure6_endpoint_ploidy_join_audit.csv"
    _, endpoint_rows = read_delimited_rows(endpoint_path, ",")
    if len(endpoint_rows) != len(canonical):
        raise ValueError("Endpoint-ploidy join audit row count differs from canonical cells")
    endpoint_tumor = [row for row in endpoint_rows if row["context"] == "Tumor"]
    if len(endpoint_tumor) != len(tumor) or any(
        row["matched"].upper() != "TRUE" for row in endpoint_tumor
    ):
        raise ValueError("Endpoint-ploidy join audit does not match every tumor cell")

    run_config_path = run_root / "metadata/run_config.tsv"
    _, run_config_rows = read_tsv(run_config_path)
    run_config = {row["key"]: row["value"] for row in run_config_rows}
    table_mode = run_config.get("table_mode")
    if table_mode not in {"full_reanalysis", "canonical_cache"}:
        raise ValueError(f"Unsupported SI Figures table mode: {table_mode}")

    input_names = {
        Path(row.get("absolute_path") or row.get("path", "")).name for row in input_rows
    }
    required_inputs = {"figure7_config.yaml", "generate_supplementary_figures.R"}
    if table_mode == "full_reanalysis":
        required_inputs.update(
            {
                "seurat_metadata.csv",
                "scvelo_cell_metrics.csv",
                "all_ploidy.tsv",
                "integrated_sct_cca_seurat_final_reclustered.rds",
            }
        )
    else:
        required_inputs.update(SI_FIGURE_TABLE_CACHE_FILES)
    if not required_inputs <= input_names:
        raise ValueError(
            "SI Figures input manifest is incomplete: "
            f"{sorted(required_inputs - input_names)}"
        )

    expected_config = {
        "module": "si_figures",
        "figures": "4,5,6,7",
        "logical_panels": str(len(SI_FIGURE_PANEL_FILES)),
        "figure_file_count": str(2 * len(SI_FIGURE_PANEL_FILES)),
        "si_figure7_from_raw_rds": (
            "true" if table_mode == "full_reanalysis" else "false"
        ),
        "raw_seurat_loaded": (
            "true" if table_mode == "full_reanalysis" else "false"
        ),
    }
    if table_mode == "canonical_cache":
        expected_config["deg_analysis_executed"] = "false"
    for key, expected in expected_config.items():
        if run_config.get(key) != expected:
            raise ValueError(
                f"SI Figures run config mismatch for {key}: "
                f"expected={expected}; observed={run_config.get(key)}"
            )

    _, qc_rows = read_tsv(run_root / "metadata/input_qc.tsv")
    qc = {row["key"]: row["value"] for row in qc_rows}
    if (
        int(qc.get("seurat_cells", "-1")) != len(canonical)
        or int(qc.get("tumor_cells", "-1")) != len(tumor)
        or int(qc.get("cellline_cells", "-1")) != len(cellline)
        or int(qc.get("tumor_cells_missing_endpoint_ploidy", "-1")) != 0
        or int(qc.get("cellline_cells_with_endpoint_ploidy", "-1")) != 0
        or int(qc.get("panel_files", "-1")) != 2 * len(SI_FIGURE_PANEL_FILES)
    ):
        raise ValueError("SI Figures QC counts do not reconcile to canonical cells/panels")

    _, provenance_rows = read_tsv(run_root / "metadata/si_figures_provenance.tsv")
    provenance = {row["key"]: row["value"] for row in provenance_rows}
    sha_keys = ["script_sha256", "figure7_config_sha256"]
    if table_mode == "full_reanalysis":
        sha_keys.extend(
            (
                "seurat_metadata_sha256",
                "scvelo_metrics_sha256",
                "all_ploidy_sha256",
                "seurat_rds_sha256",
            )
        )
    else:
        expected_not_read = {
            "seurat_metadata_sha256": "not_read_in_plot_only_mode",
            "scvelo_metrics_sha256": "not_read_in_plot_only_mode",
            "all_ploidy_sha256": "not_read_in_plot_only_mode",
            "seurat_rds_sha256": "not_run",
        }
        for key, expected in expected_not_read.items():
            if provenance.get(key) != expected:
                raise ValueError(
                    f"SI Figures plot-only provenance mismatch for {key}: "
                    f"expected={expected}; observed={provenance.get(key)}"
                )
    for key in sha_keys:
        value = provenance.get(key, "")
        if len(value) != 64 or any(char not in "0123456789abcdef" for char in value):
            raise ValueError(f"SI Figures provenance has an invalid SHA-256: {key}")
    return

    canonical_path = run_root / "tables/si_figures_cell_metadata.csv"
    cluster_key_path = run_root / "tables/si_figures_cluster_key.tsv"
    composition_mouse_path = run_root / "tables/cluster_composition_by_mouse.csv"
    composition_group_path = run_root / "tables/cluster_composition_by_ploidy_dose.csv"
    provenance_path = run_root / "metadata/si_figures_provenance.tsv"
    formal_paths = (
        canonical_path,
        cluster_key_path,
        composition_mouse_path,
        composition_group_path,
        provenance_path,
    )

    for path in formal_paths:
        matches = output_rows_by_path.get(path.resolve(), [])
        if len(matches) != 1:
            raise ValueError(f"Expected one output-manifest row for {path}; found {len(matches)}")
        row = matches[0]
        if row.get("role") != "output_table" or row.get("source_kind") != "generated_table":
            raise ValueError(f"SI Figures 4-7 contract artifact is not a generated table: {path}")
        if row.get("module") != "si_figures" or row.get("command_id") != source_run_id:
            raise ValueError(f"SI Figures 4-7 contract provenance mismatch for {path}")
        if row.get("sha256", "").strip() != sha256_file(path):
            raise ValueError(f"SI Figures 4-7 contract checksum mismatch for {path}")

    input_names = {
        Path(row.get("absolute_path") or row.get("path", "")).name for row in input_rows
    }
    required_inputs = {
        "figure7_config.yaml",
        "seurat_metadata.csv",
        "scvelo_cell_metrics.csv",
        "seurat_metadata_provenance.tsv",
    }
    if not required_inputs <= input_names:
        raise ValueError(
            "SI Figures 4-7 input manifest is missing its reviewed input bundle/config: "
            f"{sorted(required_inputs - input_names)}"
        )

    canonical_fields, canonical = read_delimited_rows(canonical_path, ",")
    expected_canonical_fields = [
        "cell_id",
        "UMAP_1",
        "UMAP_2",
        "sample_id",
        "cluster_id",
        "cluster_annotation",
        "cluster_order",
        "cluster_color",
        "initial_ploidy",
        "treatment",
        "dose",
        "dose_mg_per_kg",
        "cellcycle_classification",
        "context",
        "included_in_si_figures",
        "exclusion_reason",
    ]
    if canonical_fields != expected_canonical_fields:
        raise ValueError(
            "SI Figures 4-7 canonical cell schema mismatch: "
            f"expected={expected_canonical_fields}; observed={canonical_fields}"
        )
    cell_ids = [row["cell_id"].strip() for row in canonical]
    if any(not value for value in cell_ids) or len(cell_ids) != len(set(cell_ids)):
        raise ValueError("SI Figures 4-7 canonical cell IDs must be nonempty and unique")

    cluster_fields, cluster_rows = read_delimited_rows(cluster_key_path, "\t")
    expected_cluster_fields = [
        "cluster_id",
        "cluster_annotation",
        "cellcycle_classification",
        "cluster_order",
        "color",
        "n_all_cells",
        "n_included_tumor_cells",
    ]
    if cluster_fields != expected_cluster_fields:
        raise ValueError("SI Figures 4-7 cluster-key schema mismatch")
    cluster_by_id = {row["cluster_id"]: row for row in cluster_rows}
    if len(cluster_by_id) != len(cluster_rows):
        raise ValueError("SI Figures 4-7 cluster IDs are duplicated")
    orders = sorted(int(row["cluster_order"]) for row in cluster_rows)
    colors = [row["color"] for row in cluster_rows]
    if orders != list(range(1, len(cluster_rows) + 1)):
        raise ValueError("SI Figures 4-7 cluster order must be unique and continuous from one")
    if len(colors) != len(set(colors)) or any(
        len(color) != 7 or not color.startswith("#") for color in colors
    ):
        raise ValueError("SI Figures 4-7 cluster colors must be unique six-digit hex values")
    if any(
        not row["cluster_annotation"].strip()
        or row["cellcycle_classification"] not in {"CellCycle", "NonCellCycle"}
        for row in cluster_rows
    ):
        raise ValueError("SI Figures 4-7 cluster annotation/classification is incomplete")

    all_counts: Counter[str] = Counter()
    included_counts: Counter[tuple[str, str]] = Counter()
    included_sample_totals: Counter[str] = Counter()
    included_sample_design: dict[str, tuple[str, str]] = {}
    for row in canonical:
        cluster_id = row["cluster_id"]
        if cluster_id not in cluster_by_id:
            raise ValueError(f"Canonical cell uses cluster absent from cluster key: {cluster_id}")
        key = cluster_by_id[cluster_id]
        if (
            row["cluster_annotation"] != key["cluster_annotation"]
            or row["cluster_order"] != key["cluster_order"]
            or row["cluster_color"] != key["color"]
            or row["cellcycle_classification"] != key["cellcycle_classification"]
        ):
            raise ValueError(f"Canonical cell disagrees with cluster key: {row['cell_id']}")
        try:
            coordinates = (float(row["UMAP_1"]), float(row["UMAP_2"]))
        except ValueError as exc:
            raise ValueError("Canonical UMAP coordinate is nonnumeric") from exc
        if not all(math.isfinite(value) for value in coordinates):
            raise ValueError("Canonical UMAP coordinate is nonfinite")
        included = row["included_in_si_figures"].upper()
        if included not in {"TRUE", "FALSE"}:
            raise ValueError("Canonical inclusion flag must be TRUE or FALSE")
        if included == "TRUE" and row["exclusion_reason"].strip():
            raise ValueError("Included canonical cell has an exclusion reason")
        if included == "FALSE" and not row["exclusion_reason"].strip():
            raise ValueError("Excluded canonical cell lacks an exclusion reason")
        all_counts[cluster_id] += 1
        if included == "TRUE":
            if row["context"] != "Tumor":
                raise ValueError("An included canonical cell is not a tumor cell")
            if (
                row["initial_ploidy"] not in {"2N", "4N"}
                or row["dose"] not in {"0mg/kg", "30mg/kg", "120mg/kg"}
                or row["treatment"]
                != ("Control" if row["dose"] == "0mg/kg" else "Gemcitabine")
                or not row["sample_id"].strip()
            ):
                raise ValueError("An included canonical cell has an invalid design field")
            expected_dose_numeric = float(row["dose"].removesuffix("mg/kg"))
            if not math.isclose(
                float(row["dose_mg_per_kg"]),
                expected_dose_numeric,
                rel_tol=0,
                abs_tol=1e-12,
            ):
                raise ValueError("Canonical dose label and numeric dose disagree")
            design = (row["initial_ploidy"], row["dose"])
            previous_design = included_sample_design.setdefault(row["sample_id"], design)
            if previous_design != design:
                raise ValueError("An included sample has inconsistent ploidy or dose")
            included_counts[(row["sample_id"], cluster_id)] += 1
            included_sample_totals[row["sample_id"]] += 1

    for cluster_id, key in cluster_by_id.items():
        if int(key["n_all_cells"]) != all_counts[cluster_id]:
            raise ValueError(f"Cluster-key all-cell count mismatch for {cluster_id}")
        observed_included = sum(
            count for (sample, cluster), count in included_counts.items() if cluster == cluster_id
        )
        if int(key["n_included_tumor_cells"]) != observed_included:
            raise ValueError(f"Cluster-key included-cell count mismatch for {cluster_id}")

    mouse_fields, mouse_rows = read_delimited_rows(composition_mouse_path, ",")
    expected_mouse_fields = [
        "mouse",
        "initial_ploidy",
        "dose",
        "cluster_final",
        "cluster_annotation",
        "cluster_order",
        "n_cells",
        "total_cells",
        "proportion",
        "denominator_definition",
    ]
    if mouse_fields != expected_mouse_fields:
        raise ValueError("SI Figures 4-7 mouse-composition schema mismatch")
    proportion_sums: Counter[str] = Counter()
    mouse_rows_by_key: dict[tuple[str, str], dict[str, str]] = {}
    for row in mouse_rows:
        key = (row["mouse"], row["cluster_final"])
        if key in mouse_rows_by_key:
            raise ValueError(f"Composition row is duplicated for {key}")
        mouse_rows_by_key[key] = row
        if row["cluster_final"] not in cluster_by_id:
            raise ValueError(f"Composition uses an unknown cluster: {key}")
        cluster_key = cluster_by_id[row["cluster_final"]]
        if (
            row["cluster_annotation"] != cluster_key["cluster_annotation"]
            or row["cluster_order"] != cluster_key["cluster_order"]
            or (row["initial_ploidy"], row["dose"])
            != included_sample_design.get(row["mouse"])
            or not row["denominator_definition"].strip()
        ):
            raise ValueError(f"Composition metadata mismatch for {key}")
        if int(row["n_cells"]) != included_counts[key]:
            raise ValueError(f"Composition count mismatch for {key}")
        if int(row["total_cells"]) != included_sample_totals[row["mouse"]]:
            raise ValueError(f"Composition denominator mismatch for {row['mouse']}")
        expected = included_counts[key] / included_sample_totals[row["mouse"]]
        if not math.isclose(float(row["proportion"]), expected, rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Composition proportion mismatch for {key}")
        proportion_sums[row["mouse"]] += float(row["proportion"])
    expected_mouse_keys = {
        (sample, cluster_id)
        for sample in included_sample_totals
        for cluster_id in cluster_by_id
    }
    if set(mouse_rows_by_key) != expected_mouse_keys:
        raise ValueError("SI Figures 4-7 mouse composition is not a complete sample-by-cluster grid")
    if any(not math.isclose(value, 1.0, rel_tol=0, abs_tol=1e-12) for value in proportion_sums.values()):
        raise ValueError("SI Figures 4-7 mouse composition does not sum to one")

    group_fields, group_rows = read_delimited_rows(composition_group_path, ",")
    required_group_fields = {
        "initial_ploidy",
        "dose",
        "cluster_final",
        "cluster_annotation",
        "cluster_order",
        "n_mice",
        "sum_n_cells",
        "sum_total_cells",
        "mean_proportion",
        "sd_proportion",
        "min_proportion",
        "max_proportion",
    }
    if set(group_fields) != required_group_fields or not group_rows:
        raise ValueError("SI Figures 4-7 group-composition schema mismatch")
    observed_group_rows: dict[tuple[str, str, str], dict[str, str]] = {}
    expected_group_pieces: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for row in mouse_rows:
        key = (row["initial_ploidy"], row["dose"], row["cluster_final"])
        expected_group_pieces.setdefault(key, []).append(row)
    for row in group_rows:
        key = (row["initial_ploidy"], row["dose"], row["cluster_final"])
        if key in observed_group_rows:
            raise ValueError(f"Group composition row is duplicated for {key}")
        observed_group_rows[key] = row
    if set(observed_group_rows) != set(expected_group_pieces):
        raise ValueError("SI Figures 4-7 group composition does not cover the mouse composition")
    for key, pieces in expected_group_pieces.items():
        row = observed_group_rows[key]
        proportions = [float(piece["proportion"]) for piece in pieces]
        expected_sd = statistics.stdev(proportions) if len(proportions) > 1 else math.nan
        observed_sd_text = row["sd_proportion"].strip()
        observed_sd = (
            math.nan
            if observed_sd_text.upper() in {"", "NA", "NAN"}
            else float(observed_sd_text)
        )
        cluster_key = cluster_by_id[key[2]]
        if (
            row["cluster_annotation"] != cluster_key["cluster_annotation"]
            or row["cluster_order"] != cluster_key["cluster_order"]
            or int(row["n_mice"]) != len(pieces)
            or int(row["sum_n_cells"]) != sum(int(piece["n_cells"]) for piece in pieces)
            or int(row["sum_total_cells"]) != sum(int(piece["total_cells"]) for piece in pieces)
            or not math.isclose(
                float(row["mean_proportion"]),
                statistics.mean(proportions),
                rel_tol=0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(row["min_proportion"]), min(proportions), rel_tol=0, abs_tol=1e-12
            )
            or not math.isclose(
                float(row["max_proportion"]), max(proportions), rel_tol=0, abs_tol=1e-12
            )
            or (
                len(proportions) == 1
                and not math.isnan(observed_sd)
            )
            or (
                len(proportions) > 1
                and not math.isclose(observed_sd, expected_sd, rel_tol=0, abs_tol=1e-12)
            )
        ):
            raise ValueError(f"Group composition summary mismatch for {key}")

    _, provenance_rows = read_delimited_rows(provenance_path, "\t")
    provenance = {row["key"]: row["value"] for row in provenance_rows}
    required_provenance = {
        "source_seurat_rds_sha256",
        "umap_reduction",
        "cluster_id_field",
        "base_cluster_field",
        "clustering_resolution",
        "cluster_annotation_field",
        "source_qc_fields",
        "source_qc_policy",
        "cell_inclusion_context",
        "cell_inclusion_required_nonmissing",
        "canonical_cell_rows",
        "included_cell_rows",
        "excluded_cell_rows",
        "canonical_cell_table_sha256",
        "cluster_key_sha256",
        "composition_by_mouse_sha256",
        "composition_by_ploidy_dose_sha256",
        "upstream_analysis_documentation_doi",
        "upstream_analysis_documentation_url",
        "upstream_analysis_documentation_statement",
    }
    missing = sorted(required_provenance - provenance.keys())
    if missing:
        raise ValueError(f"SI Figures 4-7 provenance is missing keys: {missing}")
    if any(not provenance[key].strip() for key in required_provenance):
        raise ValueError("SI Figures 4-7 provenance contains an empty required value")
    source_hash = provenance["source_seurat_rds_sha256"]
    if len(source_hash) != 64 or any(char not in "0123456789abcdef" for char in source_hash):
        raise ValueError("SI Figures 4-7 source Seurat RDS checksum is invalid")
    if (
        provenance["umap_reduction"] != "umap"
        or provenance["cluster_id_field"] != "clusters"
        or provenance["base_cluster_field"] != "integrated_snn_res.0.6"
        or provenance["clustering_resolution"] != "0.6"
        or provenance["cluster_annotation_field"] != "cluster_cell_cycle_annotation"
    ):
        raise ValueError("SI Figures 4-7 provenance does not match the reviewed analysis contract")
    if (
        provenance["upstream_analysis_documentation_doi"]
        != "10.5281/zenodo.21463392"
        or provenance["upstream_analysis_documentation_url"]
        != "https://zenodo.org/records/21463392"
    ):
        raise ValueError("SI Figures 4-7 upstream-analysis documentation reference is invalid")
    if int(provenance["canonical_cell_rows"]) != len(canonical):
        raise ValueError("SI Figures 4-7 provenance canonical row count mismatch")
    if int(provenance["included_cell_rows"]) != sum(included_sample_totals.values()):
        raise ValueError("SI Figures 4-7 provenance included row count mismatch")
    if int(provenance["excluded_cell_rows"]) != len(canonical) - sum(
        included_sample_totals.values()
    ):
        raise ValueError("SI Figures 4-7 provenance excluded row count mismatch")
    expected_hashes = {
        "canonical_cell_table_sha256": canonical_path,
        "cluster_key_sha256": cluster_key_path,
        "composition_by_mouse_sha256": composition_mouse_path,
        "composition_by_ploidy_dose_sha256": composition_group_path,
    }
    for key, path in expected_hashes.items():
        if provenance[key] != sha256_file(path):
            raise ValueError(f"SI Figures 4-7 provenance checksum mismatch: {key}")


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
    figure7_tgi_day: int,
) -> None:
    if module not in STRICT_FIGURE_MODULES:
        return
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
    run_config = run_root / "metadata" / "run_config.tsv"
    if not run_config.is_file():
        raise FileNotFoundError(f"Missing source run config: {run_config}")
    _, config_rows = read_tsv(run_config)
    panel_contract = run_root / "metadata" / "panel_contract.tsv"
    if not panel_contract.is_file():
        raise FileNotFoundError(f"Missing source panel contract: {panel_contract}")
    _, contract_rows = read_tsv(panel_contract)
    if module == "in_vivo_figure7":
        panel_f_paths = {
            (run_root / str(spec["source"])).resolve()
            for spec in selected_specs
            if str(spec["module"]) == module and str(spec["panel"]).startswith("7F")
        }
        present_panel_f_paths = {path for path in panel_f_paths if path.is_file()}
        if present_panel_f_paths and present_panel_f_paths != panel_f_paths:
            raise ValueError("Source Figure 7 run must contain both PDF and PNG panel-F assets or neither")
        has_panel_f = present_panel_f_paths == panel_f_paths
        expected_panel_set = "a-f" if has_panel_f else "a-e"
        panel_set_rows = [row for row in config_rows if row.get("key") == "panel_set"]
        if len(panel_set_rows) != 1 or panel_set_rows[0].get("value") != expected_panel_set:
            raise ValueError(
                f"Source Figure 7 run must explicitly record panel_set={expected_panel_set}"
            )
        tgi_day_rows = [row for row in config_rows if row.get("key") == "tgi_day"]
        if len(tgi_day_rows) != 1 or tgi_day_rows[0].get("value") != str(figure7_tgi_day):
            raise ValueError(
                f"Source Figure 7 run must explicitly record tgi_day={figure7_tgi_day}"
            )
        expected_contract = [
            (str(spec["panel"]), Path(str(spec["source"])).name)
            for spec in selected_specs
            if str(spec["module"]) == module
            and spec.get("variant", "pdf") == "pdf"
            and (not spec.get("optional") or (run_root / str(spec["source"])).is_file())
        ]
    else:
        module_rows = [row for row in config_rows if row.get("key") == "module"]
        if len(module_rows) != 1 or module_rows[0].get("value") != "si_figures":
            raise ValueError("Source SI Figures 4-7 run must explicitly record module=si_figures")
        expected_contract = [
            (str(spec["panel"]), Path(str(spec["source"])).name)
            for spec in selected_specs
            if str(spec["module"]) == module
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
    if module == "si_figures":
        validate_si_figures_data_contract(
            run_root,
            repo_root,
            source_run_id,
            input_rows,
            rows_by_path,
        )


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
    parser.add_argument(
        "--touched-manifest-list",
        type=Path,
        help="Optional file receiving the absolute manifest paths updated by this invocation.",
    )
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

    if "si_figures" in module_runs and args.overwrite:
        supplementary_dir = figure_root / "Supplementary"
        selected_si_assets = {
            str(spec["asset"])
            for spec in selected_specs
            if str(spec["module"]) == "si_figures"
        }
        if supplementary_dir.is_dir():
            for pattern in ("panel_SuppFig4*", "panel_SuppFig5*", "panel_SuppFig6*", "panel_SuppFig7*"):
                for existing in supplementary_dir.glob(pattern):
                    if (
                        existing.is_file()
                        and existing.suffix.lower() in FIGURE_SUFFIXES
                        and existing.name not in selected_si_assets
                    ):
                        existing.unlink()

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
        staged_asset = out_dir / f".{asset.name}.tmp.{operation_id}"
        if staged_asset.exists():
            staged_asset.unlink()
        shutil.copy2(source, staged_asset)
        if sha256_file(staged_asset) != sha256_file(source):
            staged_asset.unlink(missing_ok=True)
            raise OSError(f"Staged asset checksum does not match source: {asset}")
        staged_asset.replace(asset)
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
        if manifest_path.is_file():
            existing_fields, existing_rows = read_tsv(manifest_path)
            if tuple(existing_fields) != tuple(FIGURE_MANIFEST_COLUMNS):
                raise ValueError(f"Existing manifest schema mismatch: {manifest_path}")
            replaced_panels = {row.get("panel", "") for row in rows}
            preserved_rows = [
                row for row in existing_rows if row.get("panel", "") not in replaced_panels
            ]
            rows = sorted(preserved_rows + rows, key=lambda row: row.get("panel", ""))
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        staged_manifest = manifest_path.with_name(
            f".{manifest_path.name}.tmp.{operation_id}"
        )
        write_tsv(staged_manifest, rows, FIGURE_MANIFEST_COLUMNS)
        staged_manifest.replace(manifest_path)
        print(f"Wrote {manifest_path}")

    if args.touched_manifest_list:
        touched_manifest_list = args.touched_manifest_list
        if not touched_manifest_list.is_absolute():
            touched_manifest_list = repo_root / touched_manifest_list
        touched_manifest_list.parent.mkdir(parents=True, exist_ok=True)
        touched_manifest_list.write_text(
            "".join(
                f"{figure_root / figure / 'manifest.tsv'}\n"
                for figure in sorted(touched_figures)
            ),
            encoding="utf-8",
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
