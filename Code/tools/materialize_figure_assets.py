#!/usr/bin/env python3
"""Copy selected source panels into manuscript-facing figures/ folders."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
import shutil
import statistics
import struct
from collections import Counter
from pathlib import Path

from figure_output_contract import (
    FIGURE_MANIFEST_COLUMNS,
    GENERATED_SOURCE_KINDS,
    first_local_path,
    module_manifest_local_path,
    path_within,
    read_tsv,
    repo_root_from,
    resolve_repo_path,
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
        "caption_role": "Supplementary Figure 4 final D/F/I-only composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig4_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig4_png",
        "asset": "panel_SuppFig4_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 4 D/F/I-only composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig5_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig5",
        "asset": "panel_SuppFig5_composite.pdf",
        "caption_role": "Supplementary Figure 5 final E/G-only composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig5_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig5_png",
        "asset": "panel_SuppFig5_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 5 E/G-only composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig6_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig6",
        "asset": "panel_SuppFig6_composite.pdf",
        "caption_role": "Supplementary Figure 6 final A-E composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig6_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig6_png",
        "asset": "panel_SuppFig6_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 6 A-E composite",
        "variant": "png",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig7_composite.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig7",
        "asset": "panel_SuppFig7_composite.pdf",
        "caption_role": "Supplementary Figure 7 final A-only composite",
        "variant": "pdf",
    },
    {
        "module": "si_figures",
        "source": "figures/panel_SuppFig7_composite.png",
        "figure": "Supplementary",
        "panel": "SuppFig7_png",
        "asset": "panel_SuppFig7_composite.png",
        "caption_role": "PNG derivative of the Supplementary Figure 7 A-only composite",
        "variant": "png",
    },
    {
        "module": "in_vivo_endpoint_flow",
        "source": "figures/panel_SuppFig9_endpoint_flow_cytometry.pdf",
        "figure": "Supplementary",
        "panel": "SuppFig9",
        "asset": "panel_SuppFig9_endpoint_flow_cytometry.pdf",
        "caption_role": (
            "Supplementary Figure 9 panels A-B: endpoint-flow gating provenance "
            "and DNA-content distributions"
        ),
        "variant": "pdf",
    },
    {
        "module": "in_vivo_endpoint_flow",
        "source": "figures/panel_SuppFig9_endpoint_flow_cytometry.png",
        "figure": "Supplementary",
        "panel": "SuppFig9_png",
        "asset": "panel_SuppFig9_endpoint_flow_cytometry.png",
        "caption_role": "PNG derivative of Supplementary Figure 9 panels A-B",
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
        "caption_role": (
            "Equal-mouse CellCycle mean-ECDF comparisons and exact "
            "treated-cell-excess localization"
        ),
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
        "caption_role": (
            "Day-17 TGI versus mean endpoint tumor-cell ploidy using the "
            "unadjusted mouse-level Pearson association"
        ),
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
        "caption_role": (
            "PNG derivative of the CellCycle mean-ECDF and "
            "treated-cell-excess localization panel"
        ),
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
        "caption_role": (
            "PNG derivative of the raw endpoint-ploidy association panel"
        ),
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
        "panel": "7A-7L_composite",
        "asset": "Figure7_reviewed_GRCh.png",
        "caption_role": (
            "Main Figure 7 A-L composite assembled in manuscript order"
        ),
        "variant": "png",
        "optional": True,
        "contract": True,
    },
    {
        "module": "in_vivo_figure7",
        "source": "figures/Figure7_reviewed_GRCh.pdf",
        "figure": "Figure7",
        "panel": "7A-7L_composite_pdf",
        "asset": "Figure7_reviewed_GRCh.pdf",
        "caption_role": (
            "Vector PDF of the publication-scale main Figure 7 A-L composite"
        ),
        "variant": "pdf",
        "optional": True,
        "contract": False,
    },
]

STRICT_FIGURE_MODULES = {"in_vivo_figure7", "si_figures"}
FIGURE_SUFFIXES = {".pdf", ".png", ".jpg", ".jpeg", ".svg", ".tif", ".tiff"}
FIGURE7_REVIEWED_REFERENCE_ID = (
    "state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4"
)
FIGURE7_REVIEWED_REFERENCE_KIND = (
    "reviewed_human_only_initial_ploidy_computed_pointwise_interval"
)


def read_png_geometry(path: Path) -> tuple[int, int, float | None, float | None]:
    """Return PNG width, height, and optional x/y DPI from its pHYs chunk."""
    with path.open("rb") as handle:
        if handle.read(8) != b"\x89PNG\r\n\x1a\n":
            raise ValueError(f"Invalid PNG signature: {path}")
        width = height = None
        dpi_x = dpi_y = None
        while True:
            length_bytes = handle.read(4)
            if not length_bytes:
                break
            if len(length_bytes) != 4:
                raise ValueError(f"Truncated PNG chunk length: {path}")
            length = struct.unpack(">I", length_bytes)[0]
            chunk_type = handle.read(4)
            payload = handle.read(length)
            checksum = handle.read(4)
            if len(chunk_type) != 4 or len(payload) != length or len(checksum) != 4:
                raise ValueError(f"Truncated PNG chunk: {path}")
            if chunk_type == b"IHDR":
                if length != 13:
                    raise ValueError(f"Invalid PNG IHDR length: {path}")
                width, height = struct.unpack(">II", payload[:8])
            elif chunk_type == b"pHYs" and length == 9:
                pixels_x, pixels_y, unit = struct.unpack(">IIB", payload)
                if unit == 1:
                    dpi_x = pixels_x * 0.0254
                    dpi_y = pixels_y * 0.0254
            elif chunk_type == b"IEND":
                break
    if width is None or height is None:
        raise ValueError(f"PNG has no IHDR geometry: {path}")
    return width, height, dpi_x, dpi_y
FIGURE7_REVIEWED_REFERENCE_ROOT = (
    Path("Data/in-vivo/figure7/saved_state_pathway")
    / FIGURE7_REVIEWED_REFERENCE_ID
)
FIGURE7_REVIEWED_FILES = {
    "state_pathway_interval_definition.tsv": "8f51b4c24ebdbb3da6c391be03acaacded60c01acf29c6dde77387bd2a6c8b99",
    "panel_7F_pathway_activity_plot_data.tsv": "7d16e7e0ef4a3258c413a43307437282d5d0cdd7f0ec84baabd15c2bbbeefae3",
    "panel_7F_selected_pathway_gsea.tsv": "1350e2ed472d0be0a9c809eb85ebb7a4ca05a8d8fb704567ddd65b2e87ebf5bf",
    "panel_7F_leading_edge_genes.tsv": "888c75f988c29f39271db03aae1089df894f847b2ec656320abe21e999569a2e",
    "state_pathway_gene_ranking_complete.tsv": "150d8f3ce96eade2926b6dd27fae859f36c5bcf3dbaa83f856cdcd73aa4c67b0",
    "state_pathway_gsea_complete.tsv": "9659891bc273677832355dbf85a17805e62bbcb10884f6a698ee637d389f540d",
    "state_pathway_sample_bin_coverage.tsv": "0ef5470eca70272bacdbbd0d416b21a9ad4f9e448ac172cad4ba06cb7f2719aa",
    "state_pathway_design_qc.tsv": "b430aafcb4f00c6e1ac42bd380ed540eb84c585c5556e3d724680ddd99a8300d",
    "state_pathway_provenance.tsv": "2a5f15402cb54ca67b368ffba4adcf931b642fbb6cec1ba083780fc28f1bce9b",
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
SI_REVIEWED_CBS_MANIFEST_SHA256 = (
    "756c644df06c95f95ccd7a6a1d7bfbcc972b7873ebc1188aac7da5b72f1876f9"
)
SI_INJECTED_REFERENCE_MANIFEST_SHA256 = (
    "be07657b14493523f982498dba74599a3b8185b6f6d6e7e4ab56a13d838d0492"
)
SI_ENDPOINT_PLOIDY_SHA256 = (
    "6db48ee5f196b37b58aa71d0472dd3deb06aaacb4b637070af1b27d9425db2b3"
)
FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256 = (
    "80f4e6b78e7b6d8b73030da4889ecb5c09ee97c9f83fb771aec4d3908511b569"
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
FIGURE7_PANEL_L_SOURCE_INPUTS = (
    Path("Code/in-vivo/figure7/src/tgi_data.R"),
    Path("Code/in-vivo/figure7/src/tgi_statistics.R"),
    Path("Code/in-vivo/figure7/src/tgi_panels.R"),
)
FIGURE7_PANEL_L_RESULTS = {
    17: {
        "estimate": -0.6984010192530142,
        "asymptotic_p": 0.0540069781511847,
        "permutation_p_two_sided": 0.059126984126984125,
    },
    24: {
        "estimate": -0.387348978978662,
        "asymptotic_p": 0.343097626516948,
        "permutation_p_two_sided": 0.34692460317460316,
    },
    31: {
        "estimate": -0.2969484557892882,
        "asymptotic_p": 0.475086351653848,
        "permutation_p_two_sided": 0.4857142857142857,
    },
}
FIGURE7_PANEL_L_SAMPLE_DESIGN = {
    ("2N-A2-0", "2N", "30mg/kg", "30"),
    ("2N-A2-L", "2N", "30mg/kg", "30"),
    ("2N-A4-R", "2N", "120mg/kg", "120"),
    ("2N-A4-RL", "2N", "120mg/kg", "120"),
    ("A6-4N-O", "4N", "30mg/kg", "30"),
    ("A6-4N-RR", "4N", "30mg/kg", "30"),
    ("4N-A8-RL", "4N", "120mg/kg", "120"),
    ("4N-A8-RR", "4N", "120mg/kg", "120"),
}
FIGURE7_PANEL_L_SAMPLE_FILES = {
    "2N-A2-0": "SUM159-2N-30-0_harvest.sps.cbs",
    "2N-A2-L": "SUM159-2N-30-L_harvest.sps.cbs",
    "2N-A4-R": "SUM159-2N-120-R_harvest.sps.cbs",
    "2N-A4-RL": "SUM159-2N-120-RL_harvest.sps.cbs",
    "A6-4N-O": "SUM159-4N-30-0_harvest.sps.cbs",
    "A6-4N-RR": "SUM159-4N-30-RR_harvest.sps.cbs",
    "4N-A8-RL": "SUM159-4N-120-RL_harvest.sps.cbs",
    "4N-A8-RR": "SUM159-4N-120-RR_harvest.sps.cbs",
}
FIGURE7_CURATED_SAMPLE_COUNTS = {
    "2N-A1-0": 413,
    "2N-A1-LR": 369,
    "2N-A1-R": 196,
    "2N-A1-RR": 1495,
    "2N-A2-0": 317,
    "2N-A2-L": 505,
    "2N-A4-R": 305,
    "2N-A4-RL": 1280,
    "4N-A5-0": 888,
    "4N-A5-RR": 393,
    "A5-4N-L": 385,
    "A5-4N-R": 358,
    "A6-4N-O": 189,
    "A6-4N-RR": 660,
    "4N-A8-RL": 1832,
    "4N-A8-RR": 247,
}
FIGURE7J_COPY_NUMBER_ANNOTATION_COLUMNS = [
    "heatmap_row_id",
    "file",
    "cell_id",
    "sample_id",
    "initial_ploidy",
    "dose_mg_per_kg",
    "endpoint_ploidy",
    "frac_covered",
    "display_order",
]

# Panel identifiers retired by a current compositor contract must not survive
# targeted materialization as stale aliases in the same figure manifest.
SUPERSEDED_PANEL_IDS = {
    "Figure7": {"7A-7K_composite", "7A-7K_composite_pdf"},
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


def finite_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} is not numeric: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{label} is not finite: {value!r}")
    return parsed


def validate_figure7_copy_number_outputs(run_root: Path, repo_root: Path) -> None:
    """Bind main Figure 7J to the exact final-QC tumor universe."""
    annotation_path = (
        run_root
        / "metadata"
        / "figure7J_copy_number_cell_annotations.tsv"
    ).resolve()
    matrix_path = (
        run_root
        / "metadata"
        / "figure7J_copy_number_heatmap_matrix.rds"
    ).resolve()
    for label, path in (
        ("Figure 7J copy-number cell annotations", annotation_path),
        ("Figure 7J copy-number heatmap matrix", matrix_path),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"Missing {label}: {path}")

    annotation_headers, annotation_rows = read_tsv(annotation_path)
    if (
        annotation_headers != FIGURE7J_COPY_NUMBER_ANNOTATION_COLUMNS
        or len(annotation_rows) != 9832
    ):
        raise ValueError(
            "Figure 7J copy-number cell annotations must have the exact reviewed "
            "schema and 9,832 final-QC tumor rows"
        )

    audit_path = (
        repo_root
        / "Data"
        / "in-vivo"
        / "SIfigures"
        / "si_figure6_endpoint_ploidy_join_audit.csv"
    )
    if not audit_path.is_file():
        raise FileNotFoundError(
            f"Missing canonical SI6 endpoint-ploidy audit: {audit_path}"
        )
    canonical_by_key: dict[tuple[str, str], dict[str, str]] = {}
    with audit_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != [
            "cell",
            "context",
            "initial_ploidy",
            "endpoint_file",
            "endpoint_cell_id",
            "endpoint_ploidy",
            "matched",
        ]:
            raise ValueError("Canonical SI6 endpoint-ploidy audit schema is invalid")
        for row in reader:
            matched = row["matched"].strip().lower() in {
                "1",
                "true",
                "t",
                "yes",
            }
            if row["context"].strip() != "Tumor" or not matched:
                continue
            key = (
                row["endpoint_file"].strip(),
                row["endpoint_cell_id"].strip(),
            )
            if not all(key) or key in canonical_by_key:
                raise ValueError(
                    "Canonical SI6 endpoint-ploidy audit has missing or "
                    "duplicated matched CBS cell keys"
                )
            canonical_by_key[key] = row
    if len(canonical_by_key) != 9832:
        raise ValueError(
            "Canonical SI6 endpoint-ploidy audit must identify exactly "
            "9,832 matched final-QC tumor cells"
        )

    observed_by_key: dict[tuple[str, str], dict[str, str]] = {}
    heatmap_row_ids: set[str] = set()
    display_orders: set[int] = set()
    sample_counts: Counter[str] = Counter()
    origin_counts: Counter[str] = Counter()
    treated_count = 0
    for row in annotation_rows:
        file_name = row["file"].strip()
        cell_id = row["cell_id"].strip()
        sample_id = row["sample_id"].strip()
        origin = row["initial_ploidy"].strip()
        dose = row["dose_mg_per_kg"].strip()
        heatmap_row_id = row["heatmap_row_id"].strip()
        key = (file_name, cell_id)
        canonical = canonical_by_key.get(key)
        expected_heatmap_row_id = f"{file_name}::{cell_id}"
        full_cell_id = f"{sample_id}_{cell_id}"
        file_parts = file_name.split("-")
        expected_dose = file_parts[2] if len(file_parts) >= 4 else ""
        endpoint_ploidy = finite_float(
            row["endpoint_ploidy"],
            "SI6 annotated endpoint ploidy",
        )
        frac_covered = finite_float(
            row["frac_covered"],
            "SI6 annotated covered fraction",
        )
        try:
            display_order = int(row["display_order"])
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "Figure 7J copy-number display_order must contain integers"
            ) from exc
        if (
            canonical is None
            or key in observed_by_key
            or not file_name.endswith(".sps.cbs")
            or Path(file_name).name != file_name
            or not sample_id
            or sample_id not in FIGURE7_CURATED_SAMPLE_COUNTS
            or heatmap_row_id != expected_heatmap_row_id
            or heatmap_row_id in heatmap_row_ids
            or full_cell_id != canonical["cell"].strip()
            or origin != canonical["initial_ploidy"].strip()
            or len(file_parts) < 4
            or file_parts[0] != "SUM159"
            or file_parts[1] != origin
            or dose != expected_dose
            or dose not in {"0", "30", "120"}
            or not math.isclose(
                endpoint_ploidy,
                finite_float(
                    canonical["endpoint_ploidy"],
                    "canonical SI6 endpoint ploidy",
                ),
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not 0 < frac_covered <= 1
            or display_order < 1
            or display_order in display_orders
        ):
            raise ValueError(
            "Figure 7J copy-number annotations do not preserve the exact "
                "canonical cell/file/barcode/ploidy identity and display contract"
            )
        observed_by_key[key] = row
        heatmap_row_ids.add(heatmap_row_id)
        display_orders.add(display_order)
        sample_counts[sample_id] += 1
        origin_counts[origin] += 1
        treated_count += int(dose != "0")

    if (
        set(observed_by_key) != set(canonical_by_key)
        or display_orders != set(range(1, 9833))
        or sample_counts != Counter(FIGURE7_CURATED_SAMPLE_COUNTS)
        or origin_counts != Counter({"2N": 4880, "4N": 4952})
        or treated_count != 5335
    ):
        raise ValueError(
            "Figure 7J copy-number annotations must contain the exact 9,832-cell "
            "final-QC universe, including 5,335 treated cells and the "
            "reviewed per-sample/origin counts"
        )

    output_manifest = run_root / "metadata" / "output_manifest.tsv"
    if not output_manifest.is_file():
        raise FileNotFoundError(f"Missing source output manifest: {output_manifest}")
    errors = validate_module_manifest(
        output_manifest,
        repo_root=repo_root,
        output_root=run_root,
    )
    if errors:
        raise ValueError("Invalid source output manifest:\n" + "\n".join(errors))
    _, output_rows = read_tsv(output_manifest)
    expected_outputs = {
        annotation_path: ("output_table", "generated_table"),
        matrix_path: ("output_model", "generated_model"),
    }
    expected_command_id = run_root.name.removesuffix("_si_figures")
    for path, (role, source_kind) in expected_outputs.items():
        matches = [
            row
            for row in output_rows
            if (
                (local_path := module_manifest_local_path(
                    row,
                    repo_root,
                    output_root=run_root,
                ))
                is not None
                and local_path.resolve() == path
            )
        ]
        if len(matches) != 1:
            raise ValueError(
                "Figure 7J output manifest must bind exactly one row for "
                f"{path.name}; found {len(matches)}"
            )
        manifest_row = matches[0]
        if (
            manifest_row.get("role") != role
            or manifest_row.get("source_kind") != source_kind
            or manifest_row.get("module") != "si_figures"
            or manifest_row.get("command_id") != expected_command_id
            or manifest_row.get("sha256", "").strip() != sha256_file(path)
        ):
            raise ValueError(
                "Figure 7J output-manifest provenance is invalid for "
                f"{path.name}"
            )


def validate_si4i_density_localization_outputs(
    run_root: Path,
    repo_root: Path,
) -> None:
    """Validate the relocated 2,881-cell density-localization panel."""
    input_manifest = run_root / "metadata" / "input_manifest.tsv"
    _, input_rows = read_tsv(input_manifest)
    observed_inputs = {
        path.resolve()
        for row in input_rows
        if (path := module_manifest_local_path(row, repo_root)) is not None
    }
    required_inputs = {
        (
            repo_root
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        ).resolve(),
        (repo_root / "Code/in-vivo/figure7/src/common_io.R").resolve(),
        (repo_root / "Code/in-vivo/figure7/src/tgi_statistics.R").resolve(),
        (repo_root / "Code/in-vivo/figure7/src/tgi_panels.R").resolve(),
        (
            repo_root / "Code/in-vivo/figure7/density_localization_config.yaml"
        ).resolve(),
        (
            repo_root
            / "Data/in-vivo/figure7/processed/"
            "CellCycleCells_pseudotime_distribution_per_sample_"
            "cell_level_with_ploidy_dose_tgi.csv"
        ).resolve(),
    }
    missing_inputs = sorted(str(path) for path in required_inputs - observed_inputs)
    if missing_inputs:
        raise ValueError(
            "Supplementary Figure 4I input manifest does not bind its exact "
            f"renderer, method, and 2,881-cell input: missing={missing_inputs}"
        )

    table_paths = {
        "grid": run_root / "metadata/si_figure4I_density_localization_grid.tsv",
        "intervals": (
            run_root / "metadata/si_figure4I_density_localization_intervals.tsv"
        ),
        "test": run_root / "metadata/si_figure4I_density_localization_test.tsv",
    }
    output_manifest = run_root / "metadata" / "output_manifest.tsv"
    _, output_rows = read_tsv(output_manifest)
    expected_command_id = run_root.name.removesuffix("_si_figures")
    for label, path in table_paths.items():
        path = path.resolve()
        if not path.is_file():
            raise FileNotFoundError(
                f"Missing Supplementary Figure 4I {label} table: {path}"
            )
        matches = [
            row
            for row in output_rows
            if (
                (local_path := module_manifest_local_path(
                    row, repo_root, output_root=run_root
                ))
                is not None
                and local_path.resolve() == path
            )
        ]
        if len(matches) != 1:
            raise ValueError(
                "Supplementary Figure 4I output manifest must bind exactly "
                f"one {label} table"
            )
        row = matches[0]
        if (
            row.get("role") != "output_table"
            or row.get("source_kind") != "generated_table"
            or row.get("module") != "si_figures"
            or row.get("command_id") != expected_command_id
            or row.get("sha256", "").strip() != sha256_file(path)
        ):
            raise ValueError(
                "Supplementary Figure 4I output-manifest provenance is "
                f"invalid for {path.name}"
            )

    grid_headers, grid_rows = read_tsv(table_paths["grid"])
    interval_headers, interval_rows = read_tsv(table_paths["intervals"])
    test_headers, test_rows = read_tsv(table_paths["test"])
    if (
        len(grid_rows) != 501
        or len(interval_rows) != 2
        or len(test_rows) != 1
        or "pseudotime" not in grid_headers
        or not {"support_type", "start", "end"}.issubset(interval_headers)
        or not {
            "analysis_id", "cell_universe", "n_cells", "n_samples",
            "pointwise_start", "pointwise_end", "simultaneous_start",
            "simultaneous_end", "raw_excess_peak_pseudotime",
            "max_abs_t_pseudotime", "global_max_t_p_two_sided",
        }.issubset(test_headers)
    ):
        raise ValueError(
            "Supplementary Figure 4I tables have an invalid schema or row count"
        )
    test_row = test_rows[0]
    if (
        test_row.get("analysis_id")
        != "equal_mouse_kde_exact_origin_stratified_max_t_v1"
        or test_row.get("cell_universe")
        != "reviewed_qc_retained_cellcycle_2881"
    ):
        raise ValueError(
            "Supplementary Figure 4I does not use the reviewed analysis identity"
        )
    expected_numeric = {
        "n_cells": 2881.0,
        "n_samples": 16.0,
        "pointwise_start": 0.296,
        "pointwise_end": 0.486,
        "simultaneous_start": 0.414,
        "simultaneous_end": 0.426,
        "raw_excess_peak_pseudotime": 0.452,
        "max_abs_t_pseudotime": 0.420,
        "global_max_t_p_two_sided": 232.0 / 4900.0,
    }
    if any(
        not math.isclose(
            finite_float(test_row.get(key, ""), f"SI4I {key}"),
            value,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for key, value in expected_numeric.items()
    ):
        raise ValueError(
            "Supplementary Figure 4I does not reproduce the reviewed intervals, "
            "peaks, and global exact result"
        )


def validate_figure7_panel_l_contract(
    run_root: Path,
    repo_root: Path,
    input_rows: list[dict[str, str]],
    output_rows: list[dict[str, str]],
    source_run_id: str,
    tgi_day: int,
) -> None:
    expected_result = FIGURE7_PANEL_L_RESULTS.get(tgi_day)
    if expected_result is None:
        raise ValueError(
            "Canonical Figure 7 panel L has reviewed results only for "
            "TGI days 17, 24, and 31"
        )

    run_config = read_unique_key_values(
        run_root / "metadata" / "run_config.tsv",
        "source Figure 7 run config",
    )
    expected_endpoint_metadata = {
        "endpoint_ploidy_sha256": FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256,
        "endpoint_ploidy_n_cells": "14125",
        "endpoint_ploidy_n_files": "16",
        "endpoint_ploidy_score_universe_n_cells": "9832",
        "endpoint_ploidy_treated_score_n_cells": "5335",
        "endpoint_ploidy_score_policy": (
            "arithmetic_mean_of_finite_postprocessed_cell_ploidy_in_exact_"
            "qc_passed_cellcycle_noncellcycle_union_per_sample"
        ),
        "endpoint_ploidy_mapping_policy": (
            "exact_processed_sample_barcode_to_canonical_cbs_file_cell_and_value;"
            "score_universe=reviewed_final_seurat_tumor_cells"
        ),
    }
    if any(
        run_config.get(key) != value
        for key, value in expected_endpoint_metadata.items()
    ):
        raise ValueError(
            "Figure 7 panel L run metadata does not bind the reviewed "
            "14,125-cell endpoint-CBS inventory and exact 9,832-cell "
            "curated score-universe contract"
        )
    endpoint_locator = run_config.get("endpoint_ploidy_input", "")
    endpoint_path = resolve_repo_path(endpoint_locator, repo_root)
    if (
        endpoint_path is None
        or not path_within(endpoint_path, repo_root)
        or not endpoint_path.is_file()
        or sha256_file(endpoint_path)
        != FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256
    ):
        raise ValueError(
            "Figure 7 panel L does not bind the exact portable reviewed "
            "endpoint-ploidy table"
        )
    endpoint_input_rows = [
        row
        for row in input_rows
        if (
            (path := module_manifest_local_path(row, repo_root)) is not None
            and path.resolve() == endpoint_path.resolve()
        )
    ]
    if (
        len(endpoint_input_rows) != 1
        or endpoint_input_rows[0].get("sha256")
        != FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256
    ):
        raise ValueError(
            "Figure 7 source input manifest must bind exactly one reviewed "
            "14,125-cell endpoint-ploidy table"
        )

    endpoint_headers, endpoint_rows = read_tsv(endpoint_path)
    if endpoint_headers != [
        "file",
        "cell_id",
        "ploidy",
        "total_chromosomes",
        "frac_covered",
        "format",
    ] or len(endpoint_rows) != 14125:
        raise ValueError(
            "Figure 7 endpoint-ploidy source does not have the reviewed "
            "14,125-cell schema"
        )
    endpoint_keys: set[tuple[str, str]] = set()
    endpoint_scores: dict[str, list[float]] = {}
    endpoint_score_by_key: dict[tuple[str, str], float] = {}
    for row in endpoint_rows:
        file_name = row.get("file", "")
        cell_id = row.get("cell_id", "")
        key = (file_name, cell_id)
        score = finite_float(row.get("ploidy", ""), "endpoint cell ploidy")
        coverage = finite_float(
            row.get("frac_covered", ""),
            "endpoint cell covered fraction",
        )
        total_chromosomes = finite_float(
            row.get("total_chromosomes", ""),
            "endpoint cell total chromosome count",
        )
        if (
            not file_name.endswith(".sps.cbs")
            or Path(file_name).name != file_name
            or not cell_id
            or key in endpoint_keys
            or score <= 0
            or total_chromosomes <= 0
            or not 0 < coverage <= 1
            or row.get("format") != "wide"
        ):
            raise ValueError(
                "Figure 7 endpoint-ploidy source contains an invalid or "
                "duplicated CBS cell"
            )
        endpoint_keys.add(key)
        endpoint_score_by_key[key] = score
        endpoint_scores.setdefault(file_name, []).append(score)
    if len(endpoint_scores) != 16:
        raise ValueError(
            "Figure 7 endpoint-ploidy source must contain exactly 16 CBS files"
        )

    curated_scores_by_sample: dict[str, list[float]] = {}
    curated_file_by_sample: dict[str, str] = {}
    curated_keys: set[tuple[str, str]] = set()
    compartment_counts: Counter[str] = Counter()
    for relative_path, expected_hash in FIGURE7_PROCESSED_INPUTS.items():
        role = "cellcycle" if relative_path.name.startswith("CellCycle") else "noncellcycle"
        locator = run_config.get(f"{role}_input", "")
        processed_path = resolve_repo_path(locator, repo_root)
        if (
            processed_path is None
            or not path_within(processed_path, repo_root)
            or processed_path.resolve() != (repo_root / relative_path).resolve()
            or not processed_path.is_file()
            or run_config.get(f"{role}_sha256") != expected_hash
            or sha256_file(processed_path) != expected_hash
        ):
            raise ValueError(
                "Figure 7 panel L does not bind the exact portable reviewed "
                f"{role} score-universe table"
            )
        processed_input_rows = [
            row
            for row in input_rows
            if (
                (path := module_manifest_local_path(row, repo_root)) is not None
                and path.resolve() == processed_path.resolve()
            )
        ]
        if (
            len(processed_input_rows) != 1
            or processed_input_rows[0].get("sha256") != expected_hash
        ):
            raise ValueError(
                "Figure 7 source input manifest must bind exactly one reviewed "
                f"{role} score-universe table"
            )
        with processed_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            required = {
                "cell_id",
                "sample_id",
                "growth_curve_harvest",
                "cell_ploidy",
            }
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError(
                    f"Figure 7 {role} score-universe table has an invalid schema"
                )
            for row in reader:
                sample_id = row["sample_id"].strip()
                full_cell_id = row["cell_id"].strip()
                prefix = f"{sample_id}_"
                harvest = row["growth_curve_harvest"].strip()
                file_name = f"{harvest}.sps.cbs"
                if not sample_id or not full_cell_id.startswith(prefix) or not harvest:
                    raise ValueError(
                        f"Figure 7 {role} score-universe table violates the "
                        "sample-prefixed cell-key contract"
                    )
                key = (file_name, full_cell_id[len(prefix) :])
                processed_score = finite_float(
                    row["cell_ploidy"],
                    f"Figure 7 {role} cell ploidy",
                )
                canonical_score = endpoint_score_by_key.get(key)
                if (
                    key in curated_keys
                    or canonical_score is None
                    or not math.isclose(
                        processed_score,
                        canonical_score,
                        rel_tol=0.0,
                        abs_tol=1e-12,
                    )
                ):
                    raise ValueError(
                        "Figure 7 curated score-universe cells do not map "
                        "uniquely and exactly to the canonical CBS inventory"
                    )
                previous_file = curated_file_by_sample.setdefault(sample_id, file_name)
                if previous_file != file_name:
                    raise ValueError(
                        "Figure 7 curated score universe maps one sample to "
                        "multiple CBS files"
                    )
                curated_keys.add(key)
                curated_scores_by_sample.setdefault(sample_id, []).append(canonical_score)
                compartment_counts[role] += 1
    curated_counts = {
        sample_id: len(scores)
        for sample_id, scores in curated_scores_by_sample.items()
    }
    if (
        len(curated_keys) != 9832
        or compartment_counts != Counter({"cellcycle": 2881, "noncellcycle": 6951})
        or curated_counts != FIGURE7_CURATED_SAMPLE_COUNTS
        or sum(
            curated_counts[sample_id]
            for sample_id in FIGURE7_PANEL_L_SAMPLE_FILES
        )
        != 5335
    ):
        raise ValueError(
            "Figure 7 panel L score universe must be the exact 9,832-cell "
            "QC-passed CellCycle + NonCellCycle union (5,335 treated cells)"
        )

    table_paths = {
        "plot": (run_root / "tables" / "panel_7E_plot_data.tsv").resolve(),
        "test": (run_root / "tables" / "panel_7E_test.tsv").resolve(),
    }
    manifest_rows_by_path: dict[Path, list[dict[str, str]]] = {}
    for row in output_rows:
        path = module_manifest_local_path(
            row,
            repo_root,
            output_root=run_root,
        )
        if path is not None:
            manifest_rows_by_path.setdefault(path.resolve(), []).append(row)
    for label, table_path in table_paths.items():
        if not table_path.is_file():
            raise FileNotFoundError(
                f"Missing Figure 7 panel L {label} table: {table_path}"
            )
        matches = manifest_rows_by_path.get(table_path, [])
        if len(matches) != 1:
            raise ValueError(
                "Figure 7 panel L output manifest must bind exactly one "
                f"{label} table row; found {len(matches)}"
            )
        row = matches[0]
        if (
            row.get("role") != "output_table"
            or row.get("source_kind") != "generated_table"
            or row.get("module") != "in_vivo_figure7"
            or row.get("command_id") != source_run_id
            or row.get("sha256", "").strip() != sha256_file(table_path)
        ):
            raise ValueError(
                "Figure 7 panel L output-manifest provenance is invalid for "
                f"{table_path.name}"
            )

    plot_headers, plot_rows = read_tsv(table_paths["plot"])
    test_headers, test_rows = read_tsv(table_paths["test"])
    tgi_measure = f"TGI_percent_Day_{tgi_day}"
    required_plot_headers = {
        "sample_id",
        "initial_ploidy",
        "dose",
        "dose_mg",
        "endpoint_ploidy_file",
        "sample_mean_endpoint_ploidy",
        "n_endpoint_ploidy_cells",
        "endpoint_ploidy_source_total_cells",
        "endpoint_ploidy_score_universe_total_cells",
        "endpoint_ploidy_source_file_count",
        "endpoint_ploidy_source_sha256",
        "endpoint_ploidy_score_policy",
        "endpoint_ploidy_mapping_policy",
        tgi_measure,
        "tgi_outcome",
        "tgi_day",
        "tgi_measure",
        "matched_control_summary",
        "matched_control_group",
    }
    required_test_headers = {
        "n",
        "estimate",
        "asymptotic_p",
        "permutation_p_two_sided",
        "n_permutations",
        "permutation_mode",
        "permutation_strata",
        "association_type",
        "score_variable",
        "score_source_sha256",
        "score_source_n_cells",
        "score_inventory_n_cells",
        "score_source_n_files",
        "treated_score_n_cells",
        "score_aggregation_policy",
        "sample_mapping_policy",
        "score_standardization",
        "adjustment_terms",
        "outcome_variable",
        "plot_x",
        "plot_y",
        "tgi_outcome",
        "tgi_day",
        "tgi_measure",
        "matched_control_summary",
        "matched_control_group",
    }
    adjusted_plot_headers = {
        "terminal_postprocessed_cn_score",
        "terminal_cn_score_within_origin_z",
        "permutation_stratum",
        "terminal_cn_score_nuisance_residual",
        "tgi_origin_dose_residual",
    }
    adjusted_test_headers = {
        "partial_correlation",
        "effect_per_within_origin_sd",
    }
    if (
        len(plot_rows) != 8
        or len(test_rows) != 1
        or not required_plot_headers.issubset(plot_headers)
        or not required_test_headers.issubset(test_headers)
        or adjusted_plot_headers.intersection(plot_headers)
        or adjusted_test_headers.intersection(test_headers)
        or "etp_group" in plot_headers
        or "etp_group" in test_headers
    ):
        raise ValueError(
            "Figure 7 panel L tables do not have the reviewed raw-analysis "
            "schema and eight-tumor scope"
        )

    observed_design = {
        (
            row["sample_id"],
            row["initial_ploidy"],
            row["dose"],
            row["dose_mg"],
        )
        for row in plot_rows
    }
    if observed_design != FIGURE7_PANEL_L_SAMPLE_DESIGN:
        raise ValueError(
            "Figure 7 panel L does not contain the exact reviewed treated-"
            "tumor origin-by-dose design"
        )

    metadata = {
        "tgi_outcome": "day",
        "tgi_day": str(tgi_day),
        "tgi_measure": tgi_measure,
        "matched_control_summary": "mean",
        "matched_control_group": "initial_ploidy",
    }
    if any(
        row.get(key) != value
        for row in [*plot_rows, *test_rows]
        for key, value in metadata.items()
    ):
        raise ValueError(
            "Figure 7 panel L tables do not bind the selected TGI day and "
            "initial-ploidy-matched control contract"
        )

    treated_endpoint_cell_count = 0
    raw_scores: list[float] = []
    outcome: list[float] = []
    for row in plot_rows:
        source_score = finite_float(
            row["sample_mean_endpoint_ploidy"],
            "panel L source terminal CN score",
        )
        outcome_value = finite_float(
            row[tgi_measure],
            f"panel L {tgi_measure}",
        )
        expected_file = FIGURE7_PANEL_L_SAMPLE_FILES.get(row["sample_id"])
        file_scores = curated_scores_by_sample.get(row["sample_id"], [])
        endpoint_cell_count_value = finite_float(
            row["n_endpoint_ploidy_cells"],
            "panel L endpoint CBS cell count",
        )
        endpoint_cell_count = int(endpoint_cell_count_value)
        if (
            row["endpoint_ploidy_file"] != expected_file
            or curated_file_by_sample.get(row["sample_id"]) != expected_file
            or endpoint_cell_count_value != endpoint_cell_count
            or endpoint_cell_count != len(file_scores)
            or not file_scores
            or not math.isclose(
                source_score,
                sum(file_scores) / len(file_scores),
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or row["endpoint_ploidy_source_total_cells"] != "14125"
            or row["endpoint_ploidy_score_universe_total_cells"] != "9832"
            or row["endpoint_ploidy_source_file_count"] != "16"
            or row["endpoint_ploidy_source_sha256"]
            != FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256
            or row["endpoint_ploidy_score_policy"]
            != expected_endpoint_metadata["endpoint_ploidy_score_policy"]
            or row["endpoint_ploidy_mapping_policy"]
            != expected_endpoint_metadata["endpoint_ploidy_mapping_policy"]
        ):
            raise ValueError(
                "Figure 7 panel L score is not the per-mouse mean over the "
                "exact QC-passed curated cells from its mapped CBS file"
            )
        treated_endpoint_cell_count += endpoint_cell_count
        raw_scores.append(source_score)
        outcome.append(outcome_value)
    if treated_endpoint_cell_count != 5335:
        raise ValueError(
            "Figure 7 panel L must use the exact 5,335 QC-passed curated "
            "CBS cells from the eight treated tumors"
        )

    def pearson(left: list[float], right: list[float]) -> float:
        left_mean = statistics.mean(left)
        right_mean = statistics.mean(right)
        left_centered = [value - left_mean for value in left]
        right_centered = [value - right_mean for value in right]
        left_sum_squares = sum(value * value for value in left_centered)
        right_sum_squares = sum(value * value for value in right_centered)
        denominator = math.sqrt(left_sum_squares * right_sum_squares)
        if denominator <= 0 or not math.isfinite(denominator):
            raise ValueError(
                "Figure 7 panel L raw Pearson association is not estimable"
            )
        return sum(
            left_value * right_value
            for left_value, right_value in zip(left_centered, right_centered)
        ) / denominator

    recomputed_r = pearson(raw_scores, outcome)
    n_permutations = math.factorial(len(outcome))
    extreme_permutations = 0
    for permutation in itertools.permutations(range(len(outcome))):
        permuted_outcome = [outcome[index] for index in permutation]
        permuted_r = pearson(raw_scores, permuted_outcome)
        if abs(permuted_r) >= abs(recomputed_r) - 1e-15:
            extreme_permutations += 1
    recomputed_permutation_p = extreme_permutations / n_permutations

    test = test_rows[0]
    expected_method = {
        "n": "8",
        "n_permutations": "40320",
        "permutation_mode": "exact_TGI_label_enumeration",
        "permutation_strata": "none",
        "association_type": "unadjusted_mouse_level_pearson",
        "score_variable": "sample_mean_qc_passed_curated_cbs_cell_ploidy",
        "score_source_sha256": FIGURE7_PANEL_L_ENDPOINT_PLOIDY_SHA256,
        "score_source_n_cells": "9832",
        "score_inventory_n_cells": "14125",
        "score_source_n_files": "16",
        "treated_score_n_cells": "5335",
        "score_aggregation_policy": expected_endpoint_metadata[
            "endpoint_ploidy_score_policy"
        ],
        "sample_mapping_policy": expected_endpoint_metadata[
            "endpoint_ploidy_mapping_policy"
        ],
        "score_standardization": "none",
        "adjustment_terms": "none",
        "outcome_variable": tgi_measure,
        "plot_x": "sample_mean_endpoint_ploidy",
        "plot_y": tgi_measure,
    }
    if any(test.get(key) != value for key, value in expected_method.items()):
        raise ValueError(
            "Figure 7 panel L does not use the reviewed raw mouse-level "
            "Pearson and exact unrestricted permutation method"
        )
    numeric_results = {
        "estimate": recomputed_r,
        "asymptotic_p": expected_result["asymptotic_p"],
        "permutation_p_two_sided": recomputed_permutation_p,
    }
    if any(
        not math.isclose(
            finite_float(test.get(key, ""), f"panel L {key}"),
            value,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for key, value in numeric_results.items()
    ):
        raise ValueError(
            "Figure 7 panel L does not reproduce the exact reviewed "
            f"Day-{tgi_day} raw Pearson correlation and exact P"
        )
    recomputed_reviewed = {
        "estimate": recomputed_r,
        "permutation_p_two_sided": recomputed_permutation_p,
    }
    if any(
        not math.isclose(
            recomputed_reviewed[key],
            expected_result[key],
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for key in recomputed_reviewed
    ):
        raise ValueError(
            "Figure 7 panel L curated-CBS inputs no longer reproduce the "
            f"reviewed Day-{tgi_day} raw Pearson contract"
        )


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
        "displayed_panel_sets": (
            "SuppFig4=D,F,I;SuppFig5=E,G;SuppFig6=A,B,C,D,E;SuppFig7=A"
        ),
    }
    displayed_contract = run_root / "metadata" / "displayed_panel_contract.tsv"
    if not displayed_contract.is_file():
        raise FileNotFoundError(
            f"Missing SI displayed-panel contract: {displayed_contract}"
        )
    displayed_headers, displayed_rows = read_tsv(displayed_contract)
    expected_displayed_rows = [
        {
            "figure_id": "SuppFig4",
            "panel_ids": "D,F,I",
            "selection_policy": "panels_cited_in_manuscript_results",
        },
        {
            "figure_id": "SuppFig5",
            "panel_ids": "E,G",
            "selection_policy": "panels_cited_in_manuscript_results",
        },
        {
            "figure_id": "SuppFig6",
            "panel_ids": "A,B,C,D,E",
            "selection_policy": "panels_cited_in_manuscript_results",
        },
        {
            "figure_id": "SuppFig7",
            "panel_ids": "A",
            "selection_policy": "panels_cited_in_manuscript_results",
        },
    ]
    if (
        displayed_headers
        != ["figure_id", "panel_ids", "selection_policy"]
        or displayed_rows != expected_displayed_rows
    ):
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the source "
            "run does not contain the exact manuscript-displayed panel sets"
        )
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
    copy_number_helper = (
        repo_root
        / "Code"
        / "in-vivo"
        / "SI_figures"
        / "copy_number_heatmap.R"
    )
    if not copy_number_helper.is_file():
        raise FileNotFoundError(
            f"Missing SI Figures copy-number helper: {copy_number_helper}"
        )
    endpoint_ploidy_locator = provenance.get("endpoint_ploidy_source", "")
    all_ploidy = resolve_repo_path(endpoint_ploidy_locator, repo_root)
    if (
        all_ploidy is None
        or not path_within(all_ploidy, repo_root)
        or not all_ploidy.is_file()
    ):
        raise FileNotFoundError(
            "Missing portable SI Figures postprocessed copy-number table: "
            f"{endpoint_ploidy_locator}"
        )
    if sha256_file(all_ploidy) != SI_ENDPOINT_PLOIDY_SHA256:
        raise ValueError(
            "Canonical SI Figures require the exact reviewed 14,125-cell "
            "endpoint-ploidy table, whether tracked or reconstructed into "
            "the source run"
        )
    cbs_root = repo_root / "Data" / "in-vivo" / "scRNAseq_Numbat"
    cbs_manifest = cbs_root / "cbs_manifest.tsv"
    if not cbs_manifest.is_file():
        raise FileNotFoundError(
            f"Missing reviewed CBS manifest: {cbs_manifest}"
        )
    if sha256_file(cbs_manifest) != SI_REVIEWED_CBS_MANIFEST_SHA256:
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the "
            "repository does not contain the exact reviewed CBS manifest"
        )
    cbs_headers, cbs_rows = read_tsv(cbs_manifest)
    cbs_names = [row.get("filename", "") for row in cbs_rows]
    if (
        cbs_headers != ["filename", "bytes", "sha256", "notes"]
        or len(cbs_names) != 16
        or len(set(cbs_names)) != 16
        or any(not name or Path(name).name != name for name in cbs_names)
    ):
        raise ValueError(
            "Reviewed CBS manifest must contain exactly 16 unique matrices"
        )
    cbs_matrices = [cbs_root / name for name in cbs_names]
    observed_cbs_names = sorted(
        path.name for path in cbs_root.glob("*.sps.cbs")
    )
    if observed_cbs_names != sorted(cbs_names):
        raise ValueError(
            "Canonical SI Figures require the exact reviewed CBS inventory"
        )
    for path, row in zip(cbs_matrices, cbs_rows):
        if (
            not path.is_file()
            or row.get("bytes") != str(path.stat().st_size)
            or row.get("sha256") != sha256_file(path)
        ):
            raise ValueError(
                "Canonical SI Figures require the reviewed CBS matrix: "
                f"{path.name}"
            )
    cbs_hashes = ";".join(
        f"{path.name}={sha256_file(path)}" for path in cbs_matrices
    )
    reference_root = cbs_root / "injected_reference"
    reference_manifest = reference_root / "reference_manifest.tsv"
    if (
        not reference_manifest.is_file()
        or sha256_file(reference_manifest)
        != SI_INJECTED_REFERENCE_MANIFEST_SHA256
    ):
        raise ValueError(
            "Canonical SI Figures require the reviewed injected-cell "
            "reference manifest"
        )
    reference_headers, reference_rows = read_tsv(reference_manifest)
    reference_names = [row.get("filename", "") for row in reference_rows]
    if (
        reference_headers
        != [
            "filename",
            "injected_origin",
            "reference_label",
            "n_cells",
            "bytes",
            "sha256",
            "source_repository",
            "source_commit",
            "designation_basis",
            "chr999_unit",
            "chr999_interpretation_basis",
            "ploidy_policy",
        ]
        or len(reference_names) != 2
        or {row.get("injected_origin") for row in reference_rows}
        != {"2N", "4N"}
        or len(set(reference_names)) != 2
        or any(not name or Path(name).name != name for name in reference_names)
    ):
        raise ValueError(
            "Injected-cell reference manifest must pin one 2N and one 4N matrix"
        )
    reference_matrices = [reference_root / name for name in reference_names]
    if sorted(path.name for path in reference_root.glob("*.sps.cbs")) != sorted(
        reference_names
    ):
        raise ValueError(
            "Canonical SI Figures require the exact injected-reference inventory"
        )
    for path, row in zip(reference_matrices, reference_rows):
        if (
            not path.is_file()
            or row.get("bytes") != str(path.stat().st_size)
            or row.get("sha256") != sha256_file(path)
        ):
            raise ValueError(
                "Canonical SI Figures require the injected-cell reference: "
                f"{path.name}"
            )
    reference_hashes = ";".join(
        f"{path.name}={sha256_file(path)}" for path in reference_matrices
    )
    copy_number_expected = {
        "figure7j_copy_number_source": (
            "postprocessed NUMBAT-derived cell-by-segment CBS matrices"
        ),
        "figure7j_column_statistic": (
            "per-cell length-weighted mean across available CBS segments "
            "within each chromosome and file-specific schema"
        ),
        "figure7j_cbs_matrix_count": "16",
        "figure7j_source_cell_count": "14125",
        "figure7j_qc_passed_tumor_cell_count": "9832",
        "figure7j_treated_tumor_cell_count": "5335",
        "figure7j_qc_selection_policy": (
            "exact cells matched by the frozen endpoint-ploidy audit from the "
            "final QC-curated Seurat object; clusters 3, 4, 9, and 9c excluded"
        ),
        "figure7j_chromosome_count": "22",
        "figure7j_row_order": (
            "injected origin, dose, mouse, post-processed copy-number score, "
            "cell ID; no row clustering"
        ),
        "figure7j_column_order": (
            "chromosomes 1-22 in genomic order; no column clustering"
        ),
        "si6_postprocessed_copy_number_score_unit": (
            "sequenced mouse/CBS file"
        ),
        "si6_injected_reference_ploidy_policy": (
            "autosomal length-weighted estimate plus chr999 "
            "haploid-genome-equivalent unassigned DNA"
        ),
        "si6_injected_reference_cell_counts": "2N=20;4N=16",
        "si6_injected_reference_source_repository": "miningcloneid",
        "si6_injected_reference_source_commit": (
            "c505cd9159fa2a8c0974c7379f6aacd09fe19abc"
        ),
        "si6_injected_reference_chr999_unit": (
            "haploid-genome-equivalent unassigned DNA"
        ),
        "si6_injected_reference_chr999_interpretation_basis": (
            "project-confirmed 2026-08-01"
        ),
        "si6_endpoint_summary_analysis_type": "descriptive_only",
        "si6_2n_reference_mean_ploidy": "2.151242953153243",
        "si6_2n_endpoint_mouse_balanced_mean_ploidy": "2.135513297009228",
        "si6_2n_relative_change_percent": "-0.7311892002230036",
        "si6_4n_reference_mean_ploidy": "3.94651299970957",
        "si6_4n_endpoint_mouse_balanced_mean_ploidy": "2.318552949440266",
        "si6_4n_relative_change_percent": "-41.250593888557",
        "si6_reference_4n_minus_2n_mean_ploidy": "1.795270046556327",
        "si6_endpoint_4n_minus_2n_mouse_balanced_mean_ploidy": (
            "0.1830396524310385"
        ),
        "si6_separation_contraction_percent": "89.80433875214797",
    }
    si4i_expected = {
        "si4i_cell_universe": "reviewed_qc_retained_cellcycle_2881",
        "si4i_n_cells": "2881",
        "si4i_n_mice": "16",
        "si4i_pointwise_positive_interval": "0.296-0.486",
        "si4i_simultaneous_positive_interval": "0.414-0.426",
        "si4i_global_max_abs_t_p_two_sided": "0.0473469387755102",
    }
    si4i_cellcycle = (
        repo_root
        / "Data/in-vivo/figure7/processed/"
        "CellCycleCells_pseudotime_distribution_per_sample_"
        "cell_level_with_ploidy_dose_tgi.csv"
    )
    si4i_config = repo_root / "Code/in-vivo/figure7/density_localization_config.yaml"
    si4i_provenance_expected = {
        "si4i_cellcycle_pseudotime": (
            "Data/in-vivo/figure7/processed/"
            "CellCycleCells_pseudotime_distribution_per_sample_"
            "cell_level_with_ploidy_dose_tgi.csv"
        ),
        "si4i_cellcycle_pseudotime_sha256": sha256_file(si4i_cellcycle),
        "si4i_density_localization_config": (
            "Code/in-vivo/figure7/density_localization_config.yaml"
        ),
        "si4i_density_localization_config_sha256": sha256_file(si4i_config),
        "si4i_density_localization_method": (
            "equal-mouse Gaussian-kernel treated-minus-vehicle density "
            "contrast; exact injected-origin-stratified pointwise and "
            "studentized max-|T| permutation inference on the reviewed "
            "2,881-cell subset"
        ),
    }
    qc_selection_note = (
        "the complete 14,125-cell CBS source is checksum/value validated, "
        "then restricted by the frozen endpoint-ploidy audit to the exact "
        "9,832 QC-passed tumor cells (5,335 treated); clusters 3, 4, 9, "
        "and 9c remain excluded"
    )
    ploidy_reduction_note = (
        "project-designated lineage-matched 2N-A7M/4N-A5M karyotype "
        "reference distributions, including the chr999 unassigned-extra-DNA "
        "haploid-genome-equivalent term added to autosomal ploidy, compared "
        "descriptively with one postprocessed endpoint mean per mouse; no "
        "formal P value because each reference is one culture-level biological "
        "unit and origin-specific endpoint runs use different "
        "schemas/calibration"
    )
    if (
        any(
            provenance.get(key) != value
            for key, value in si7_expected.items()
        )
        or any(
            run_config.get(key) != value
            for key, value in {
                **si7_expected,
                **composition_expected,
                **copy_number_expected,
                **si4i_expected,
            }.items()
        )
        or any(
            provenance.get(key) != value
            for key, value in si4i_provenance_expected.items()
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
        or provenance.get("copy_number_heatmap_helper") != (
            "Code/in-vivo/SI_figures/copy_number_heatmap.R"
        )
        or provenance.get("copy_number_heatmap_helper_sha256")
        != sha256_file(copy_number_helper)
        or provenance.get("endpoint_ploidy_derivation_helper") != (
            "Data/in-vivo/weighted_ploidy.py"
        )
        or provenance.get("endpoint_ploidy_derivation_helper_sha256")
        != sha256_file(repo_root / "Data/in-vivo/weighted_ploidy.py")
        or provenance.get("endpoint_ploidy_source")
        != endpoint_ploidy_locator
        or provenance.get("endpoint_ploidy_source_sha256")
        != sha256_file(all_ploidy)
        or provenance.get("numbat_cbs_manifest") != (
            "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"
        )
        or provenance.get("numbat_cbs_manifest_sha256")
        != SI_REVIEWED_CBS_MANIFEST_SHA256
        or provenance.get("numbat_cbs_matrix_count") != "16"
        or provenance.get("numbat_cbs_matrix_hashes") != cbs_hashes
        or provenance.get("figure7j_qc_selection") != qc_selection_note
        or provenance.get("injected_cell_reference_manifest") != (
            "Data/in-vivo/scRNAseq_Numbat/injected_reference/"
            "reference_manifest.tsv"
        )
        or provenance.get("injected_cell_reference_manifest_sha256")
        != SI_INJECTED_REFERENCE_MANIFEST_SHA256
        or provenance.get("injected_cell_reference_matrix_hashes")
        != reference_hashes
        or provenance.get("si6e_ploidy_reduction_comparison")
        != ploidy_reduction_note
        or provenance.get("si6e_summary_analysis_type") != (
            "descriptive_only; no endpoint cross-origin or "
            "reference-to-endpoint test"
        )
        or provenance.get("table_cache_manifest_sha256")
        != sha256_file(canonical_manifest)
    ):
        raise ValueError(
            "Canonical SI Figures materialization is prohibited: the run "
            "does not carry the reviewed frozen-table publication contract"
        )
    validate_figure7_copy_number_outputs(run_root, repo_root)
    validate_si4i_density_localization_outputs(run_root, repo_root)


def validate_figure7_density_localization_contract(
    run_root: Path,
    repo_root: Path,
    input_rows: list[dict[str, str]],
    output_rows: list[dict[str, str]],
    source_run_id: str,
) -> None:
    """Validate the exact reviewed 2,881-cell localization displayed as SI4I."""
    figure7_config = repo_root / "Code/in-vivo/figure7/figure7_config.yaml"
    localization_config = (
        repo_root / "Code/in-vivo/figure7/density_localization_config.yaml"
    )
    run_config = read_unique_key_values(
        run_root / "metadata" / "run_config.tsv",
        "source Figure 7 run config",
    )
    if (
        not figure7_config.is_file()
        or not localization_config.is_file()
        or run_config.get("config_sha256") != sha256_file(figure7_config)
        or run_config.get("density_localization_config")
        != "Code/in-vivo/figure7/density_localization_config.yaml"
        or run_config.get("density_localization_config_sha256")
        != sha256_file(localization_config)
    ):
        raise ValueError(
            "Figure 7 density localization does not bind the active reviewed "
            "configuration"
        )

    required_inputs = {
        (
            repo_root / "Code/in-vivo/figure7/run_figure7.R"
        ).resolve(),
        figure7_config.resolve(),
        localization_config.resolve(),
        (
            repo_root / "Code/in-vivo/figure7/src/common_io.R"
        ).resolve(),
        (
            repo_root / "Code/in-vivo/figure7/src/tgi_statistics.R"
        ).resolve(),
        (
            repo_root / "Code/in-vivo/figure7/src/tgi_panels.R"
        ).resolve(),
    }
    observed_inputs = {
        path.resolve()
        for row in input_rows
        if (path := module_manifest_local_path(row, repo_root)) is not None
    }
    missing_inputs = sorted(str(path) for path in required_inputs - observed_inputs)
    if missing_inputs:
        raise ValueError(
            "Figure 7 density localization input manifest does not bind its "
            f"renderer, configs, statistics, and plotting code: missing={missing_inputs}"
        )

    table_paths = {
        "grid": (
            run_root / "tables/panel_7B_density_localization_grid.tsv"
        ).resolve(),
        "intervals": (
            run_root / "tables/panel_7B_density_localization_intervals.tsv"
        ).resolve(),
        "test": (
            run_root / "tables/panel_7B_density_localization_test.tsv"
        ).resolve(),
    }
    manifest_rows_by_path: dict[Path, list[dict[str, str]]] = {}
    for row in output_rows:
        path = module_manifest_local_path(
            row,
            repo_root,
            output_root=run_root,
        )
        if path is not None:
            manifest_rows_by_path.setdefault(path.resolve(), []).append(row)
    for label, table_path in table_paths.items():
        if not table_path.is_file():
            raise FileNotFoundError(
                f"Missing Figure 7 density-localization {label} table: {table_path}"
            )
        matches = manifest_rows_by_path.get(table_path, [])
        if len(matches) != 1:
            raise ValueError(
                "Figure 7 density-localization output manifest must bind "
                f"exactly one {label} table row; found {len(matches)}"
            )
        row = matches[0]
        if (
            row.get("role") != "output_table"
            or row.get("source_kind") != "generated_table"
            or row.get("module") != "in_vivo_figure7"
            or row.get("command_id") != source_run_id
            or row.get("sha256", "").strip() != sha256_file(table_path)
        ):
            raise ValueError(
                "Figure 7 density-localization output-manifest provenance is "
                f"invalid for {table_path.name}"
            )

    grid_headers, grid_rows = read_tsv(table_paths["grid"])
    interval_headers, interval_rows = read_tsv(table_paths["intervals"])
    test_headers, test_rows = read_tsv(table_paths["test"])
    required_grid_headers = {
        "pseudotime",
        "vehicle_mean_density",
        "treated_mean_density",
        "treated_minus_vehicle_density",
        "permutation_sd",
        "observed_studentized",
        "pointwise_p_two_sided",
        "max_t_adjusted_p_two_sided",
        "simultaneous_critical",
        "simultaneous_lower_envelope",
        "simultaneous_upper_envelope",
        "pointwise_positive_supported",
        "simultaneous_positive_supported",
        "modeled_state_interval",
    }
    required_interval_headers = {"support_type", "start", "end", "width", "alpha"}
    required_test_headers = {
        "analysis_id",
        "cell_universe",
        "n_cells",
        "n_samples",
        "n_vehicle_samples",
        "n_treated_samples",
        "sample_weighting",
        "density_estimator",
        "bandwidth_method",
        "bandwidth",
        "grid_start",
        "grid_end",
        "grid_points",
        "contrast",
        "permutation_strata",
        "n_permutations",
        "pointwise_test",
        "pointwise_alpha",
        "pointwise_start",
        "pointwise_end",
        "simultaneous_test",
        "simultaneous_alpha",
        "simultaneous_critical",
        "simultaneous_start",
        "simultaneous_end",
        "raw_excess_peak_pseudotime",
        "raw_excess_peak_density_difference",
        "max_abs_t_pseudotime",
        "max_abs_t_observed_statistic",
        "global_max_abs_t",
        "global_max_t_p_two_sided",
    }
    if (
        not required_grid_headers.issubset(grid_headers)
        or not required_interval_headers.issubset(interval_headers)
        or not required_test_headers.issubset(test_headers)
        or len(grid_rows) != 501
        or len(interval_rows) != 2
        or len(test_rows) != 1
    ):
        raise ValueError(
            "Figure 7 density-localization tables have an invalid schema or row count"
        )

    test_row = test_rows[0]
    expected_strings = {
        "analysis_id": "equal_mouse_kde_exact_origin_stratified_max_t_v1",
        "cell_universe": "reviewed_qc_retained_cellcycle_2881",
        "sample_weighting": "equal_mouse",
        "density_estimator": "stats::density Gaussian kernel",
        "bandwidth_method": "pooled label-invariant stats::bw.nrd0",
        "contrast": "treated_minus_vehicle",
        "permutation_strata": "initial_ploidy",
        "pointwise_test": "two-sided exact permutation at each grid point",
        "simultaneous_test": (
            "studentized max-absolute-T exact permutation null envelope"
        ),
    }
    if any(test_row.get(key) != value for key, value in expected_strings.items()):
        raise ValueError(
            "Figure 7 density localization does not use the reviewed "
            "equal-mouse exact origin-stratified method"
        )

    expected_numeric = {
        "n_cells": 2881.0,
        "n_samples": 16.0,
        "n_vehicle_samples": 8.0,
        "n_treated_samples": 8.0,
        "bandwidth": 0.0506470455660707,
        "grid_start": 0.0,
        "grid_end": 1.0,
        "grid_points": 501.0,
        "n_permutations": 4900.0,
        "pointwise_alpha": 0.05,
        "pointwise_start": 0.296,
        "pointwise_end": 0.486,
        "simultaneous_alpha": 0.05,
        "simultaneous_critical": 2.6737373166169203,
        "simultaneous_start": 0.414,
        "simultaneous_end": 0.426,
        "raw_excess_peak_pseudotime": 0.452,
        "raw_excess_peak_density_difference": 0.283948126904316,
        "max_abs_t_pseudotime": 0.420,
        "max_abs_t_observed_statistic": 2.68138128387289,
        "global_max_abs_t": 2.68138128387289,
        "global_max_t_p_two_sided": 232.0 / 4900.0,
    }
    if any(
        not math.isclose(
            finite_float(test_row.get(key, ""), f"density-localization {key}"),
            value,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for key, value in expected_numeric.items()
    ):
        raise ValueError(
            "Figure 7 density localization does not reproduce the reviewed "
            "2,881-cell intervals, peaks, or global exact result"
        )

    expected_intervals = (
        ("positive_pointwise_two_sided", 0.296, 0.486, 0.190, 0.05),
        ("positive_simultaneous_max_abs_t", 0.414, 0.426, 0.012, 0.05),
    )
    for row, expected in zip(interval_rows, expected_intervals, strict=True):
        support_type, start, end, width, alpha = expected
        if row.get("support_type") != support_type or any(
            not math.isclose(
                finite_float(row.get(key, ""), f"density-localization {key}"),
                value,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            for key, value in (
                ("start", start),
                ("end", end),
                ("width", width),
                ("alpha", alpha),
            )
        ):
            raise ValueError(
                "Figure 7 density-localization interval table is not the "
                "reviewed pointwise and simultaneous result"
            )

    raw_values: list[float] = []
    studentized_values: list[float] = []
    for index, row in enumerate(grid_rows):
        x = finite_float(row.get("pseudotime", ""), "density grid pseudotime")
        vehicle = finite_float(
            row.get("vehicle_mean_density", ""), "vehicle mean density"
        )
        treated = finite_float(
            row.get("treated_mean_density", ""), "treated mean density"
        )
        difference = finite_float(
            row.get("treated_minus_vehicle_density", ""), "density difference"
        )
        permutation_sd = finite_float(
            row.get("permutation_sd", ""), "density permutation SD"
        )
        studentized = finite_float(
            row.get("observed_studentized", ""), "studentized density difference"
        )
        pointwise_p = finite_float(
            row.get("pointwise_p_two_sided", ""), "pointwise density P"
        )
        adjusted_p = finite_float(
            row.get("max_t_adjusted_p_two_sided", ""), "adjusted density P"
        )
        critical = finite_float(
            row.get("simultaneous_critical", ""), "density critical value"
        )
        lower = finite_float(
            row.get("simultaneous_lower_envelope", ""), "density lower envelope"
        )
        upper = finite_float(
            row.get("simultaneous_upper_envelope", ""), "density upper envelope"
        )
        booleans = {
            key: row.get(key, "")
            for key in (
                "pointwise_positive_supported",
                "simultaneous_positive_supported",
                "modeled_state_interval",
            )
        }
        if any(value not in {"TRUE", "FALSE"} for value in booleans.values()):
            raise ValueError(
                "Figure 7 density-localization support masks must be strict booleans"
            )
        pointwise_supported = booleans["pointwise_positive_supported"] == "TRUE"
        simultaneous_supported = (
            booleans["simultaneous_positive_supported"] == "TRUE"
        )
        modeled = booleans["modeled_state_interval"] == "TRUE"
        expected_x = index / 500.0
        expected_pointwise = 0.296 <= expected_x <= 0.486
        expected_simultaneous = 0.414 <= expected_x <= 0.426
        expected_modeled = expected_pointwise
        if (
            not math.isclose(x, expected_x, rel_tol=0.0, abs_tol=1e-12)
            or permutation_sd <= 0
            or not math.isclose(treated - vehicle, difference, rel_tol=0.0, abs_tol=1e-12)
            or not math.isclose(
                difference / permutation_sd,
                studentized,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                critical,
                expected_numeric["simultaneous_critical"],
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                lower,
                -critical * permutation_sd,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                upper,
                critical * permutation_sd,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not 0 <= pointwise_p <= 1
            or not 0 <= adjusted_p <= 1
            or pointwise_supported != (difference > 0 and pointwise_p <= 0.05)
            or simultaneous_supported != (studentized > critical)
            or pointwise_supported != expected_pointwise
            or simultaneous_supported != expected_simultaneous
            or modeled != expected_modeled
        ):
            raise ValueError(
                "Figure 7 density-localization grid, null envelope, or support "
                "masks disagree with the reviewed contract"
            )
        raw_values.append(difference)
        studentized_values.append(studentized)

    raw_peak_index = max(range(len(raw_values)), key=raw_values.__getitem__)
    max_t_index = max(
        range(len(studentized_values)),
        key=lambda index: abs(studentized_values[index]),
    )
    if (
        raw_peak_index != 226
        or max_t_index != 210
        or not math.isclose(
            raw_values[raw_peak_index],
            expected_numeric["raw_excess_peak_density_difference"],
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        or not math.isclose(
            abs(studentized_values[max_t_index]),
            expected_numeric["global_max_abs_t"],
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    ):
        raise ValueError(
            "Figure 7 density-localization grid does not distinguish the raw "
            "0.452 peak from the 0.420 maximum standardized evidence"
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
        si_provenance = read_unique_key_values(
            run_root / "metadata" / "si_figures_provenance.tsv",
            "source SI Figures provenance",
        )
        endpoint_ploidy = resolve_repo_path(
            si_provenance.get("endpoint_ploidy_source", ""),
            repo_root,
        )
        if (
            endpoint_ploidy is None
            or not path_within(endpoint_ploidy, repo_root)
            or not endpoint_ploidy.is_file()
            or sha256_file(endpoint_ploidy) != SI_ENDPOINT_PLOIDY_SHA256
        ):
            raise ValueError(
                "SI Figures source provenance does not bind the exact portable "
                "endpoint-ploidy table"
            )
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
                / "in-vivo"
                / "SI_figures"
                / "copy_number_heatmap.R"
            ).resolve(),
            (repo_root / "Data" / "in-vivo" / "weighted_ploidy.py").resolve(),
            endpoint_ploidy.resolve(),
            (
                repo_root
                / "Data"
                / "in-vivo"
                / "scRNAseq_Numbat"
                / "cbs_manifest.tsv"
            ).resolve(),
            *{
                path.resolve()
                for path in (
                    repo_root / "Data" / "in-vivo" / "scRNAseq_Numbat"
                ).glob("*.sps.cbs")
            },
            (
                repo_root
                / "Data"
                / "in-vivo"
                / "scRNAseq_Numbat"
                / "injected_reference"
                / "reference_manifest.tsv"
            ).resolve(),
            *{
                path.resolve()
                for path in (
                    repo_root
                    / "Data"
                    / "in-vivo"
                    / "scRNAseq_Numbat"
                    / "injected_reference"
                ).glob("*.sps.cbs")
            },
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
    _, output_rows = read_tsv(output_manifest)

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
            and str(spec["panel"]) == "7A-7L_composite"
        ]
        if len(composite_specs) != 1:
            raise ValueError("Figure 7 must define one A-L composite contract")
        composite_path = (
            run_root / str(composite_specs[0]["source"])
        ).resolve()
        if composite_path.is_file() != has_panel_f:
            requirement = "present" if has_panel_f else "absent"
            raise ValueError(
                "Source Figure 7 A-L composite must be "
                f"{requirement} exactly when panel F is present"
            )
        composite_pdf_specs = [
            spec
            for spec in selected_specs
            if str(spec["module"]) == module
            and str(spec["panel"]) == "7A-7L_composite_pdf"
        ]
        if len(composite_pdf_specs) != 1:
            raise ValueError("Figure 7 must define one vector A-L composite")
        composite_pdf_path = (
            run_root / str(composite_pdf_specs[0]["source"])
        ).resolve()
        if composite_pdf_path.is_file() != has_panel_f:
            requirement = "present" if has_panel_f else "absent"
            raise ValueError(
                "Source Figure 7 vector A-L composite must be "
                f"{requirement} exactly when panel F is present"
            )
        if has_panel_f:
            width_px, height_px, dpi_x, dpi_y = read_png_geometry(composite_path)
            if (
                (width_px, height_px) != (2130, 3193)
                or dpi_x is None
                or dpi_y is None
                or abs(dpi_x - 300.0) > 0.5
                or abs(dpi_y - 300.0) > 0.5
            ):
                raise ValueError(
                    "Source Figure 7 A-L PNG must be exactly 2130x3193 pixels "
                    "with 300-DPI metadata (7.1x10.645 inches)"
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
        figure7_config = repo_root / "Code/in-vivo/figure7/figure7_config.yaml"
        density_localization_config = (
            repo_root / "Code/in-vivo/figure7/density_localization_config.yaml"
        )
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
        required_panel_l_sources = {
            (repo_root / relative_path).resolve()
            for relative_path in FIGURE7_PANEL_L_SOURCE_INPUTS
        }
        missing_panel_l_sources = sorted(
            str(path)
            for path in required_panel_l_sources - observed_input_paths
        )
        if missing_panel_l_sources:
            raise ValueError(
                "Figure 7 source input manifest does not bind the panel L "
                "statistics and plotting helpers: "
                f"missing={missing_panel_l_sources}"
            )
        validate_figure7_panel_l_contract(
            run_root,
            repo_root,
            input_rows,
            output_rows,
            source_run_id,
            figure7_tgi_day,
        )
        validate_figure7_density_localization_contract(
            run_root,
            repo_root,
            input_rows,
            output_rows,
            source_run_id,
        )
        if has_panel_f:
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
            copy_cbs_root = repo_root / "Data/in-vivo/scRNAseq_Numbat"
            copy_cbs_manifest = copy_cbs_root / "cbs_manifest.tsv"
            if (
                not copy_cbs_manifest.is_file()
                or sha256_file(copy_cbs_manifest)
                != SI_REVIEWED_CBS_MANIFEST_SHA256
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "panel J does not bind the exact reviewed CBS manifest"
                )
            copy_cbs_headers, copy_cbs_rows = read_tsv(copy_cbs_manifest)
            copy_cbs_names = [row.get("filename", "") for row in copy_cbs_rows]
            if (
                copy_cbs_headers != ["filename", "bytes", "sha256", "notes"]
                or len(copy_cbs_names) != 16
                or len(set(copy_cbs_names)) != 16
                or any(
                    not name or Path(name).name != name
                    for name in copy_cbs_names
                )
            ):
                raise ValueError(
                    "Canonical Figure 7 panel J requires exactly 16 safe "
                    "reviewed CBS matrix names"
                )
            copy_cbs_paths = {
                (copy_cbs_root / name).resolve() for name in copy_cbs_names
            }
            for row in copy_cbs_rows:
                path = (copy_cbs_root / row["filename"]).resolve()
                if (
                    not path.is_file()
                    or row.get("bytes") != str(path.stat().st_size)
                    or row.get("sha256") != sha256_file(path)
                ):
                    raise ValueError(
                        "Canonical Figure 7 panel J requires reviewed CBS "
                        f"matrix {path.name}"
                    )
            copy_ploidy = (repo_root / "Data/in-vivo/all_ploidy.tsv").resolve()
            if (
                not copy_ploidy.is_file()
                or sha256_file(copy_ploidy) != SI_ENDPOINT_PLOIDY_SHA256
            ):
                raise ValueError(
                    "Canonical Figure 7 panel J requires the exact reviewed "
                    "14,125-cell endpoint-ploidy table"
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
                (
                    repo_root
                    / "Code/in-vivo/figure7/density_localization_config.yaml"
                ).resolve(),
                (
                    repo_root
                    / "Code/in-vivo/figure7/src/copy_number_panel.R"
                ).resolve(),
                (
                    repo_root
                    / "Code/in-vivo/SI_figures/copy_number_heatmap.R"
                ).resolve(),
                (repo_root / "Data/in-vivo/all_ploidy.tsv").resolve(),
                (
                    repo_root
                    / "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"
                ).resolve(),
                (
                    repo_root
                    / "Data/in-vivo/SIfigures/si_figure6_endpoint_ploidy_join_audit.csv"
                ).resolve(),
                (
                    repo_root
                    / "Data/in-vivo/SIfigures/si_figures_cell_metadata.csv"
                ).resolve(),
                *copy_cbs_paths,
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
                    "renderer, config, and nine-file panel-7F reference: "
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
                    "c5192325e7c6ebd8531c0dfc2c8343c7032817e891f5f27"
                    "a69efd2c42fdb384c"
                )
            ):
                raise ValueError(
                    "Canonical Figure 7 materialization is prohibited: "
                    "reviewed provenance lineage is invalid"
                )
            expected_composite = {
                "main_composite_panel_set": "a-l",
                "main_composite_filename": "Figure7_reviewed_GRCh.png",
                "main_composite_panel_order": (
                    "A=7A;B=7C;C=SI4A;D=SI4B;E=SI4C;F=SI4E;"
                    "G=SI7B;H=7B;I=7F;J=7J;K=7D;L=7E"
                ),
                "main_composite_width_in": "7.1",
                "main_composite_height_in": "10.645",
                "main_composite_png_dpi": "300",
                "main_composite_layout_rows": "A/B;C/D/E;F/G;H;I/J;K/L",
                "copy_number_panel_qc_cells": "9832",
                "copy_number_panel_treated_cells": "5335",
                "copy_number_panel_mice": "16",
                "copy_number_panel_chromosomes": "22",
                "copy_number_panel_annotation_bars": (
                    "injected_origin;gemcitabine_dose;mouse"
                ),
                "copy_number_panel_row_ordering_policy": (
                    "hierarchical_clustering_separately_within_each_mouse"
                ),
                "copy_number_panel_row_distance_method": "euclidean",
                "copy_number_panel_row_linkage_method": "ward.D2",
                "copy_number_panel_row_tie_break_method": (
                    "canonical_heatmap_row_id_input_order"
                ),
                "density_localization_config": (
                    "Code/in-vivo/figure7/density_localization_config.yaml"
                ),
                "density_localization_config_sha256": sha256_file(
                    density_localization_config
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
                    "the A-L composite/cache contract is invalid"
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

    rows = output_rows
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
        replaced_asset_paths = {
            row.get("asset_path", "")
            for row in rows_by_figure.get(figure, [])
            if row.get("asset_path", "")
        }
        superseded_panels = SUPERSEDED_PANEL_IDS.get(figure, set())
        preserved_rows = [
            row
            for row in existing_rows
            if row.get("panel", "") not in replaced_panels
            and row.get("panel", "") not in superseded_panels
            and row.get("asset_path", "") not in replaced_asset_paths
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
