from __future__ import annotations

import csv
import os
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"

import sys

sys.path.insert(0, str(TOOLS_DIR))
from figure_output_contract import MODULE_MANIFEST_COLUMNS, sha256_file, write_tsv  # noqa: E402
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    panel_specs_for_figure7_variant,
)


class ManagerFigure7CliTest(unittest.TestCase):
    def _run(self, *args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPO_ROOT / "Manager.sh"), *args],
            cwd=REPO_ROOT,
            env=env,
            text=True,
            capture_output=True,
        )

    @staticmethod
    def _write_si_figure4_input_pair(root: Path) -> tuple[Path, Path]:
        root.mkdir(parents=True, exist_ok=True)
        seurat = root / "seurat_metadata.csv"
        scvelo = root / "scvelo_cell_metrics.csv"
        seurat.write_text(
            "cell,UMAP_1,UMAP_2,Dose,clusters,sample\n"
            "c1,0,0,0mg/kg,6,s1\n"
        )
        scvelo.write_text(
            "cell,velocity_cell,velocity_pseudotime,TN,clusters,sample,Ploidy,Dose\n"
            "c1,c1,0.1,Tumor,6,s1,2N,0mg/kg\n"
        )
        return seurat, scvelo

    @staticmethod
    def _write_fake_full_workflow_rscript(path: Path) -> None:
        path.write_text(
            """#!/usr/bin/env bash
set -euo pipefail
entrypoint="$1"
shift
[[ "$entrypoint" == *run_figure7.R ]] || exit 99
output_dir=""
intermediate_dir=""
seurat_metadata=""
for arg in "$@"; do
  case "$arg" in
    --output-dir=*) output_dir="${arg#*=}" ;;
    --intermediate-dir=*) intermediate_dir="${arg#*=}" ;;
    --seurat-metadata-output=*) seurat_metadata="${arg#*=}" ;;
  esac
done
mkdir -p "$output_dir/figures" "$output_dir/metadata" "$output_dir/tables" "$intermediate_dir"
scvelo="$intermediate_dir/scvelo_cell_metrics.csv"
cellcycle="$intermediate_dir/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
noncellcycle="$intermediate_dir/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv"
[[ -n "$seurat_metadata" ]] || seurat_metadata="$intermediate_dir/seurat_metadata.csv"
printf 'cell,velocity_cell,velocity_pseudotime,TN,clusters,sample,Ploidy,Dose\nc1,c1,0.1,Tumor,6,s1,2N,0mg/kg\n' > "$scvelo"
printf 'cell,UMAP_1,UMAP_2,Dose,clusters,sample\nc1,0,0,0mg/kg,6,s1\n' > "$seurat_metadata"
for table in "$cellcycle" "$noncellcycle"; do
  printf 'cell_id,sample_id,initial_ploidy,gemcitabine_dose,gemcitabine_dose_mg_per_kg,pseudotime,cell_ploidy\nc1,s1,2N,0mg/kg,0,0.1,2.0\n' > "$table"
done
scvelo_sha="$(shasum -a 256 "$scvelo" | awk '{print $1}')"
cellcycle_sha="$(shasum -a 256 "$cellcycle" | awk '{print $1}')"
noncellcycle_sha="$(shasum -a 256 "$noncellcycle" | awk '{print $1}')"
seurat_metadata_sha="$(shasum -a 256 "$seurat_metadata" | awk '{print $1}')"
if [[ "${BAD_FIGURE7_HASHES:-false}" == true ]]; then
  scvelo_sha="$(printf '0%.0s' {1..64})"
fi
for name in \
  panel_7A_day17_tgi_calculation \
  panel_7B_cellcycle_selected_ecdf_comparisons \
  panel_7C_day17_tgi_by_initial_ploidy \
  panel_7D_day17_tgi_vs_centered_ecdf_shift \
  panel_7E_day17_tgi_vs_mean_etp \
  panel_7F_pseudotime_state_pathway_activity; do
  printf 'fake pdf\n' > "$output_dir/figures/$name.pdf"
  printf 'fake png\n' > "$output_dir/figures/$name.png"
done
{
  printf 'key\tvalue\n'
  printf 'module\tin_vivo_figure7\n'
  printf 'mode\tfull-workflow\n'
  printf 'panel_set\ta-f\n'
  printf 'tgi_day\t17\n'
  printf 'workflow_executed_stages\tscvelo_metrics,celllevel_inputs\n'
  printf 'workflow_scvelo_metrics\t%s\n' "$scvelo"
  printf 'workflow_scvelo_sha256\t%s\n' "$scvelo_sha"
  printf 'workflow_cellcycle_input\t%s\n' "$cellcycle"
  printf 'workflow_cellcycle_sha256\t%s\n' "$cellcycle_sha"
  printf 'workflow_noncellcycle_input\t%s\n' "$noncellcycle"
  printf 'workflow_noncellcycle_sha256\t%s\n' "$noncellcycle_sha"
  printf 'workflow_seurat_metadata\t%s\n' "$seurat_metadata"
  printf 'workflow_seurat_metadata_sha256\t%s\n' "$seurat_metadata_sha"
} > "$output_dir/metadata/run_config.tsv"
{
  printf 'panel_id\tfilename\n'
  printf '7A\tpanel_7A_day17_tgi_calculation.pdf\n'
  printf '7B\tpanel_7B_cellcycle_selected_ecdf_comparisons.pdf\n'
  printf '7C\tpanel_7C_day17_tgi_by_initial_ploidy.pdf\n'
  printf '7D\tpanel_7D_day17_tgi_vs_centered_ecdf_shift.pdf\n'
  printf '7E\tpanel_7E_day17_tgi_vs_mean_etp.pdf\n'
  printf '7F\tpanel_7F_pseudotime_state_pathway_activity.pdf\n'
} > "$output_dir/metadata/panel_contract.tsv"
"""
        )
        path.chmod(0o755)

    def test_panels_only_requires_source_run_id(self) -> None:
        result = self._run(
            "--mode", "panels-only", "--modules", "in_vivo_figure7", "--run-id", "operation"
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--source-run-id is required", result.stderr)

    def test_default_manuscript_modules_include_figure7(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        default_line = next(
            line for line in manager_text.splitlines() if line.startswith('modules="')
        )
        self.assertIn("in_vivo_figure7", default_line)
        self.assertIn("si_figure4", default_line)
        self.assertLess(default_line.index("si_figure4"), default_line.index("in_vivo_figure7"))

    def test_si_figure4_dry_run_prefers_published_input_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data_root = tmp_path / "Data/in-vivo"
            seurat, scvelo = self._write_si_figure4_input_pair(data_root)
            env = os.environ.copy()
            env["FIGURE7_DATA_ROOT"] = str(data_root)
            result = self._run(
                "--dry-run", "--modules", "si_figure4", "--run-id", "si_published",
                "--output-root", str(tmp_path / "Results"), env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn(f"--seurat-metadata {seurat.resolve()}", result.stdout)
            self.assertIn(f"--scvelo-metrics {scvelo.resolve()}", result.stdout)
            self.assertNotIn("[si_figure4_input_prep]", result.stdout)

    def test_si_figure4_dry_run_falls_back_to_intermediate_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data_root = tmp_path / "Data/in-vivo"
            intermediate = tmp_path / "intermediates"
            seurat, scvelo = self._write_si_figure4_input_pair(intermediate)
            env = os.environ.copy()
            env["FIGURE7_DATA_ROOT"] = str(data_root)
            result = self._run(
                "--dry-run", "--modules", "si_figure4", "--run-id", "si_intermediate",
                "--output-root", str(tmp_path / "Results"),
                "--figure7-intermediate-dir", str(intermediate), env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("[si_figure4_input_publish]", result.stdout)
            self.assertIn(f"--seurat-metadata {seurat.resolve()}", result.stdout)
            self.assertIn(f"--scvelo-metrics {scvelo.resolve()}", result.stdout)
            self.assertFalse(data_root.exists())

    def test_si_figure4_dry_run_plans_scvelo_pair_generation_when_both_sources_absent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            data_root = tmp_path / "Data/in-vivo"
            intermediate = tmp_path / "intermediates"
            env = os.environ.copy()
            env["FIGURE7_DATA_ROOT"] = str(data_root)
            result = self._run(
                "--dry-run", "--modules", "si_figure4", "--run-id", "si_generate",
                "--output-root", str(tmp_path / "Results"),
                "--figure7-intermediate-dir", str(intermediate), env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("[si_figure4_input_prep]", result.stdout)
            self.assertIn("--mode=prepare-scvelo-inputs", result.stdout)
            self.assertIn(str(intermediate.resolve() / "seurat_metadata.csv"), result.stdout)
            self.assertIn(str(intermediate.resolve() / "scvelo_cell_metrics.csv"), result.stdout)
            self.assertFalse(data_root.exists())
            self.assertFalse(intermediate.exists())

    def test_default_mode_runs_figure7_full_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            inputs = [tmp_path / name for name in ("all_ploidy.tsv", "sample_info.xlsx", "growth.xlsx")]
            for path in inputs:
                path.write_text("fixture\n")
            result = self._run(
                "--dry-run", "--modules", "in_vivo_figure7", "--run-id", "default_full",
                "--figure7-cell-ploidy-input", str(inputs[0]),
                "--figure7-sample-info-input", str(inputs[1]),
                "--figure7-growth-curve-input", str(inputs[2]),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--mode=full-workflow", result.stdout)
            self.assertIn("Manager full-refit", result.stdout)

    def test_source_run_id_is_rejected_outside_panels_only(self) -> None:
        result = self._run(
            "--mode", "standard", "--modules", "in_vivo_figure7",
            "--run-id", "operation", "--source-run-id", "source",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("valid only with --mode panels-only", result.stderr)

    def test_full_analysis_requires_pinned_gene_set_artifact(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check", "--figure7-full-analysis",
            "--figure7-seurat-rds", "/tmp/object.rds",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--figure7-gene-set-artifact", result.stderr)

    def test_refresh_inputs_requires_full_refit_or_check_only(self) -> None:
        result = self._run(
            "--mode", "standard", "--modules", "in_vivo_figure7",
            "--run-id", "bad_refresh_mode", "--figure7-refresh-inputs",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires --mode full-refit or check-only", result.stderr)

    def test_refresh_check_forwards_full_workflow_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            cell_ploidy = tmp_path / "all_ploidy.tsv"
            sample_info = tmp_path / "sample_info.xlsx"
            growth_curve = tmp_path / "growth.xlsx"
            for path in (cell_ploidy, sample_info, growth_curve):
                path.write_text("fixture\n")
            intermediate = tmp_path / "intermediates"
            result = self._run(
                "--mode", "check-only", "--modules", "in_vivo_figure7",
                "--run-id", "refresh_check", "--figure7-refresh-inputs",
                "--figure7-intermediate-dir", str(intermediate),
                "--figure7-python", "/opt/scvelo/bin/python",
                "--figure7-cell-ploidy-input", str(cell_ploidy),
                "--figure7-sample-info-input", str(sample_info),
                "--figure7-growth-curve-input", str(growth_curve),
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--mode=full-workflow", result.stdout)
            self.assertIn(f"--intermediate-dir={intermediate.resolve()}", result.stdout)
            self.assertIn("--python=/opt/scvelo/bin/python", result.stdout)

    def test_refresh_publishes_four_validated_csvs_after_success(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            intermediate = tmp_path / "intermediates"
            data_root = tmp_path / "Data/in-vivo"
            cell_ploidy = tmp_path / "all_ploidy.tsv"
            sample_info = tmp_path / "sample_info.xlsx"
            growth_curve = tmp_path / "growth.xlsx"
            for path in (cell_ploidy, sample_info, growth_curve):
                path.write_text("fixture\n")
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            self._write_fake_full_workflow_rscript(fake_bin / "Rscript")
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["FIGURE7_DATA_ROOT"] = str(data_root)

            result = self._run(
                "--mode", "full-refit", "--modules", "in_vivo_figure7",
                "--run-id", "refresh_publish", "--figure7-refresh-inputs",
                "--output-root", str(output_root), "--figure-root", str(figure_root),
                "--figure7-intermediate-dir", str(intermediate),
                "--figure7-cell-ploidy-input", str(cell_ploidy),
                "--figure7-sample-info-input", str(sample_info),
                "--figure7-growth-curve-input", str(growth_curve),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            expected_pairs = (
                (
                    intermediate / "seurat_metadata.csv",
                    data_root / "seurat_metadata.csv",
                ),
                (
                    intermediate / "scvelo_cell_metrics.csv",
                    data_root / "scvelo_cell_metrics.csv",
                ),
                (
                    intermediate / "CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
                    data_root / "figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
                ),
                (
                    intermediate / "NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
                    data_root / "figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv",
                ),
            )
            for source, target in expected_pairs:
                self.assertEqual(source.read_bytes(), target.read_bytes())

            metadata_path = (
                output_root
                / "manager/runs/refresh_publish/metadata/figure7_input_materialization.tsv"
            )
            with metadata_path.open(newline="") as handle:
                metadata = {
                    row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(metadata["status"], "ok")
            self.assertEqual(metadata["materialized_file_count"], "4")

    def test_refresh_hash_failure_does_not_publish_any_csv(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            intermediate = tmp_path / "intermediates"
            data_root = tmp_path / "Data/in-vivo"
            inputs = [tmp_path / name for name in ("all_ploidy.tsv", "sample_info.xlsx", "growth.xlsx")]
            for path in inputs:
                path.write_text("fixture\n")
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            self._write_fake_full_workflow_rscript(fake_bin / "Rscript")
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["FIGURE7_DATA_ROOT"] = str(data_root)
            env["BAD_FIGURE7_HASHES"] = "true"

            result = self._run(
                "--mode", "full-refit", "--modules", "in_vivo_figure7",
                "--run-id", "refresh_bad_hash", "--figure7-refresh-inputs",
                "--output-root", str(output_root), "--figure-root", str(figure_root),
                "--figure7-intermediate-dir", str(intermediate),
                "--figure7-cell-ploidy-input", str(inputs[0]),
                "--figure7-sample-info-input", str(inputs[1]),
                "--figure7-growth-curve-input", str(inputs[2]),
                env=env,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn("changed after validation", result.stderr)
            self.assertFalse(data_root.exists())

    def test_ae_only_check_omits_panel_f_inputs_and_sets_panel_contract(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check_ae", "--figure7-panels-ae-only",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--panel-set=a-e", result.stdout)
        self.assertNotIn("--saved-state-pathway-dir", result.stdout)

    def test_tgi24_and_supplement_destination_are_forwarded(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check_tgi24",
            "--figure7-tgi-day", "24",
            "--figure7-figure-name", "Figure7_Supplement",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--tgi-day=24", result.stdout)

        specs = [
            spec
            for spec in panel_specs_for_figure7_variant(24, "Figure7_Supplement")
            if spec["module"] == "in_vivo_figure7"
        ]
        self.assertTrue(all(spec["figure"] == "Figure7_Supplement" for spec in specs))
        self.assertTrue(
            any(spec["asset"] == "panel_7A_day24_tgi_calculation.pdf" for spec in specs)
        )
        self.assertFalse(any("day17" in str(spec["asset"]) for spec in specs))

    def test_invalid_figure7_tgi_day_is_rejected(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "bad_tgi_day", "--figure7-tgi-day", "TGI_24",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("non-negative integer", result.stderr)

    def test_ae_only_rejects_full_panel_f_analysis(self) -> None:
        result = self._run(
            "--mode", "full-refit", "--modules", "in_vivo_figure7",
            "--run-id", "bad_combo", "--figure7-panels-ae-only",
            "--figure7-full-analysis", "--figure7-seurat-rds", "/tmp/object.rds",
            "--figure7-gene-set-artifact", "/tmp/gene_sets.rds",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("mutually exclusive", result.stderr)

    def test_state_pathway_export_rejects_ae_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            result = self._run(
                "--mode", "check-only", "--modules", "in_vivo_figure7",
                "--run-id", "bad_export_combo", "--figure7-panels-ae-only",
                "--figure7-state-pathway-results-root", tmp,
            )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot be used with --figure7-panels-ae-only", result.stderr)

    def test_state_pathway_results_root_runs_exporter_and_is_recorded(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_results_root = tmp_path / "state pathway results"
            source_results_root.mkdir()
            canonical_reference_root = tmp_path / "canonical Data" / "taoli_state_pathway_etp2_24_day17_v1"
            run_id = "integrated_state_export"

            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            fake_rscript = fake_bin / "Rscript"
            fake_rscript.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
entrypoint="$1"
shift
output_dir=""
saved_state_dir=""
source_results_root=""
for arg in "$@"; do
  case "$arg" in
    --output-dir=*) output_dir="${arg#*=}" ;;
    --saved-state-pathway-dir=*) saved_state_dir="${arg#*=}" ;;
    --results-root=*|--state-pathway-results-root=*) source_results_root="${arg#*=}" ;;
  esac
done
if [[ "$entrypoint" == *export_state_pathway_reference.R ]]; then
  mkdir -p "$output_dir"
  for name in \
    panel_7F_pathway_activity_plot_data.tsv \
    panel_7F_selected_pathway_gsea.tsv \
    panel_7F_leading_edge_genes.tsv \
    state_pathway_gene_ranking_complete.tsv \
    state_pathway_gsea_complete.tsv \
    state_pathway_sample_bin_coverage.tsv \
    state_pathway_design_qc.tsv; do
    printf 'value\\nfixture\\n' > "$output_dir/$name"
  done
  printf 'key\\tvalue\\nsource_results_id\\tpseudotime_state_pathways\\n' \
    > "$output_dir/state_pathway_provenance.tsv"
  exit 0
fi
if [[ "$entrypoint" == *run_figure7.R ]]; then
  mkdir -p "$output_dir/figures" "$output_dir/metadata" "$output_dir/tables"
  for name in \
    panel_7A_day17_tgi_calculation \
    panel_7B_cellcycle_selected_ecdf_comparisons \
    panel_7C_day17_tgi_by_initial_ploidy \
    panel_7D_day17_tgi_vs_centered_ecdf_shift \
    panel_7E_day17_tgi_vs_mean_etp \
    panel_7F_pseudotime_state_pathway_activity; do
    printf 'fake pdf\\n' > "$output_dir/figures/$name.pdf"
    printf 'fake png\\n' > "$output_dir/figures/$name.png"
  done
  printf 'key\\tvalue\\npanel_set\\ta-f\\ntgi_day\\t17\\nstate_pathway_source_results_root\\t%s\\n' \
    "$source_results_root" > "$output_dir/metadata/run_config.tsv"
  printf 'panel_id\\tfilename\\n' > "$output_dir/metadata/panel_contract.tsv"
  printf '7A\\tpanel_7A_day17_tgi_calculation.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7B\\tpanel_7B_cellcycle_selected_ecdf_comparisons.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7C\\tpanel_7C_day17_tgi_by_initial_ploidy.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7D\\tpanel_7D_day17_tgi_vs_centered_ecdf_shift.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7E\\tpanel_7E_day17_tgi_vs_mean_etp.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7F\\tpanel_7F_pseudotime_state_pathway_activity.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  cp "$saved_state_dir/state_pathway_provenance.tsv" "$output_dir/metadata/state_pathway_provenance.tsv"
  exit 0
fi
exit 99
"""
            )
            fake_rscript.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["FIGURE7_CANONICAL_REFERENCE_ROOT"] = str(canonical_reference_root)

            result = self._run(
                "--mode", "standard", "--modules", "in_vivo_figure7",
                "--run-id", run_id,
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                "--figure7-state-pathway-results-root", str(source_results_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            manager_root = output_root / "manager/runs" / run_id
            export_metadata_path = manager_root / "metadata/figure7_state_pathway_export.tsv"
            with export_metadata_path.open(newline="") as handle:
                export_metadata = {
                    row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(export_metadata["status"], "ok")
            self.assertEqual(export_metadata["source_results_root"], str(source_results_root.resolve()))
            reference_root = Path(export_metadata["exported_reference_dir"])
            self.assertEqual(reference_root.name, "taoli_state_pathway_etp2_24_day17_v1")
            self.assertEqual(len(list(reference_root.glob("*.tsv"))), 8)
            self.assertEqual(
                export_metadata["canonical_data_reference_dir"], str(canonical_reference_root)
            )
            self.assertEqual(len(list(canonical_reference_root.glob("*.tsv"))), 8)
            canonical_provenance = (
                canonical_reference_root / "state_pathway_provenance.tsv"
            ).read_text()
            self.assertNotIn(str(source_results_root.resolve()), canonical_provenance)
            for exported_path in reference_root.glob("*.tsv"):
                self.assertEqual(
                    exported_path.read_bytes(),
                    (canonical_reference_root / exported_path.name).read_bytes(),
                )

            with (manager_root / "metadata/figure7_state_pathway_materialization.tsv").open(
                newline=""
            ) as handle:
                materialization_metadata = {
                    row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(materialization_metadata["status"], "ok")
            self.assertEqual(
                materialization_metadata["target_data_reference_dir"],
                str(canonical_reference_root),
            )
            self.assertEqual(materialization_metadata["materialized_file_count"], "8")

            run_root = output_root / "in-vivo/figure7/runs" / f"{run_id}_figure7"
            with (run_root / "metadata/run_config.tsv").open(newline="") as handle:
                run_config = {row["key"]: row["value"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(
                run_config["state_pathway_source_results_root"], str(source_results_root.resolve())
            )
            module_runs = (manager_root / "metadata/module_runs.tsv").read_text()
            self.assertIn(f"state_pathway_source_results_root={source_results_root.resolve()}", module_runs)
            self.assertIn(f"state_pathway_canonical_data_dir={canonical_reference_root}", module_runs)
            self.assertTrue((figure_root / "Figure7/manifest.tsv").is_file())

    def test_failed_figure7_run_does_not_materialize_canonical_reference(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_results_root = tmp_path / "state_pathway_results"
            canonical_reference_root = tmp_path / "canonical" / "taoli_state_pathway_etp2_24_day17_v1"
            source_results_root.mkdir()

            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            fake_rscript = fake_bin / "Rscript"
            fake_rscript.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
entrypoint="$1"
shift
output_dir=""
for arg in "$@"; do
  case "$arg" in
    --output-dir=*) output_dir="${arg#*=}" ;;
  esac
done
if [[ "$entrypoint" == *export_state_pathway_reference.R ]]; then
  mkdir -p "$output_dir"
  for name in \
    panel_7F_pathway_activity_plot_data.tsv \
    panel_7F_selected_pathway_gsea.tsv \
    panel_7F_leading_edge_genes.tsv \
    state_pathway_gene_ranking_complete.tsv \
    state_pathway_gsea_complete.tsv \
    state_pathway_sample_bin_coverage.tsv \
    state_pathway_design_qc.tsv \
    state_pathway_provenance.tsv; do
    printf 'value\\nfixture\\n' > "$output_dir/$name"
  done
  exit 0
fi
exit 7
"""
            )
            fake_rscript.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["FIGURE7_CANONICAL_REFERENCE_ROOT"] = str(canonical_reference_root)

            result = self._run(
                "--mode", "standard", "--modules", "in_vivo_figure7",
                "--run-id", "failed_before_materialization",
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                "--figure7-state-pathway-results-root", str(source_results_root),
                env=env,
            )

            self.assertEqual(result.returncode, 7)
            self.assertFalse(canonical_reference_root.exists())
            self.assertFalse(
                (
                    output_root
                    / "manager/runs/failed_before_materialization/metadata"
                    / "figure7_state_pathway_materialization.tsv"
                ).exists()
            )

    def test_panels_only_uses_source_run_without_invoking_r(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_id = "source_analysis"
            operation_id = "materialization_operation"
            run_root = output_root / "in-vivo/figure7/runs" / f"{source_id}_figure7"
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()

            rows = []
            for spec in PANEL_SPECS:
                if spec["module"] != "in_vivo_figure7":
                    continue
                source = run_root / str(spec["source"])
                source.write_bytes(f"fake PDF {spec['panel']}\n".encode())
                rows.append(
                    {
                        "path": str(source),
                        "repo_relative_path": "",
                        "absolute_path": str(source),
                        "role": "output_figure",
                        "source_kind": "generated_panel",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(source),
                        "checksum_unavailable_reason": "",
                        "byte_size": source.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "Figure7",
                        "panel": spec["panel"],
                        "notes": "fixture",
                    }
                )
            write_tsv(run_root / "metadata/output_manifest.tsv", rows, MODULE_MANIFEST_COLUMNS)
            write_tsv(
                run_root / "metadata/run_config.tsv",
                [
                    {"key": "panel_set", "value": "a-f"},
                    {"key": "tgi_day", "value": "17"},
                ],
                ["key", "value"],
            )
            write_tsv(
                run_root / "metadata/panel_contract.tsv",
                [
                    {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                    for spec in PANEL_SPECS
                    if spec["module"] == "in_vivo_figure7"
                    and spec.get("variant", "pdf") == "pdf"
                ],
                ["panel_id", "filename"],
            )
            input_path = tmp_path / "figure7_input.tsv"
            input_path.write_text("value\n1\n")
            write_tsv(
                run_root / "metadata/input_manifest.tsv",
                [
                    {
                        "path": str(input_path),
                        "repo_relative_path": "",
                        "absolute_path": str(input_path),
                        "role": "input_data",
                        "source_kind": "input_file",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(input_path),
                        "checksum_unavailable_reason": "",
                        "byte_size": input_path.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "fixture",
                    }
                ],
                MODULE_MANIFEST_COLUMNS,
            )

            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            sentinel = tmp_path / "rscript_was_called"
            fake_rscript = fake_bin / "Rscript"
            fake_rscript.write_text(f"#!/usr/bin/env bash\ntouch '{sentinel}'\nexit 99\n")
            fake_rscript.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"

            result = self._run(
                "--mode", "panels-only",
                "--modules", "in_vivo_figure7",
                "--source-run-id", source_id,
                "--run-id", operation_id,
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(sentinel.exists(), "panels-only unexpectedly invoked Rscript")

            with (figure_root / "Figure7/manifest.tsv").open(newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["run_id"] for row in manifest_rows}, {source_id})
            self.assertTrue(
                all(f"materialization_operation_id={operation_id}" in row["notes"] for row in manifest_rows)
            )
            module_runs = (
                output_root / "manager/runs" / operation_id / "metadata/module_runs.tsv"
            ).read_text()
            self.assertIn(f"source_run_id={source_id}", module_runs)
            self.assertIn(f"materialization_operation_id={operation_id}", module_runs)


if __name__ == "__main__":
    unittest.main()
