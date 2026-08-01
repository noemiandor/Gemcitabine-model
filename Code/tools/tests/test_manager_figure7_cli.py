from __future__ import annotations

import csv
import os
import statistics
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

    def _figure7_input_lineage_paths(
        self,
        run_dir: Path,
        *,
        manager_mode: str = "full-refit",
    ) -> list[str]:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("figure7_runtime_source_paths()")
        end = manager_text.index("\nrequired_input_paths_for_module()", start)
        function_block = manager_text[start:end]
        script = f"""
set -euo pipefail
{function_block}
metadata_value() {{
  local path="$1" key="$2"
  awk -F '\\t' -v expected="${{key}}" '
    $1 == expected {{ count += 1; value = $2 }}
    END {{ if (count != 1 || value == "") exit 1; print value }}
  ' "${{path}}"
}}
mode="$3"
figure7_full_analysis=false
figure7_panels_ae_only=true
figure7_state_pathway_results_root=""
figure7_reference_root=unused
figure7_seurat_rds=""
figure7_gene_set_artifact=""
figure7_raw_seurat_rds=""
figure7_raw_data_dir="$2"
figure7_loom_root=""
figure7_cell_ploidy_input=cell_ploidy.tsv
figure7_endpoint_cbs_score_input=Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv
figure7_endpoint_cbs_score_will_be_derived=false
si_figures_cbs_dir=Data/in-vivo/scRNAseq_Numbat
figure7_sample_info_input=sample_info.xlsx
figure7_growth_curve_input=growth_curve.xlsx
input_paths_for_module in_vivo_figure7 "$1"
"""
        result = subprocess.run(
            [
                "bash",
                "-c",
                script,
                "fixture",
                str(run_dir),
                str(run_dir / "raw"),
                manager_mode,
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return [line for line in result.stdout.splitlines() if line]

    @staticmethod
    def _write_run_config(run_dir: Path, rows: dict[str, str]) -> None:
        metadata = run_dir / "metadata"
        metadata.mkdir(parents=True)
        with (metadata / "run_config.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["key", "value"])
            writer.writerows(rows.items())

    @staticmethod
    def _copy_reviewed_reference(run_dir: Path) -> None:
        reference_root = (
            REPO_ROOT
            / "Data/in-vivo/figure7/saved_state_pathway"
            / "taoli_04i_etp2_24_day17_v1"
        )
        (run_dir / "tables").mkdir(exist_ok=True)
        for source in reference_root.glob("*.tsv"):
            destination = (
                run_dir / "metadata" / source.name
                if source.name == "state_pathway_provenance.tsv"
                else run_dir / "tables" / source.name
            )
            destination.write_bytes(source.read_bytes())

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

    def test_ae_only_check_omits_panel_f_inputs_and_sets_panel_contract(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check_ae", "--figure7-panels-ae-only",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--panel-set=a-e", result.stdout)
        self.assertIn(
            "--endpoint-cbs-score-input="
            "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv",
            result.stdout,
        )
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

    def test_input_lineage_omits_unconsumed_stale_scvelo_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            stale_scvelo = run_dir / "stale_scvelo.csv"
            stale_manifest = run_dir / "scvelo_stage_manifest.tsv"
            stale_scvelo.write_text("stale\n")
            stale_manifest.write_text("key\tvalue\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(stale_scvelo),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertNotIn(str(stale_scvelo), paths)
            self.assertNotIn(str(stale_manifest), paths)

    def test_standard_input_lineage_uses_generic_processed_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            processed = run_dir / "processed"
            processed.mkdir()
            cellcycle = processed / "cellcycle.csv"
            noncellcycle = processed / "noncellcycle.csv"
            legacy_cellcycle = run_dir / "legacy-cellcycle.csv"
            legacy_noncellcycle = run_dir / "legacy-noncellcycle.csv"
            self._write_run_config(
                run_dir,
                {
                    "mode": "standard",
                    "cellcycle_input": str(cellcycle),
                    "cellcycle_sha256": "a" * 64,
                    "noncellcycle_input": str(noncellcycle),
                    "noncellcycle_sha256": "b" * 64,
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(legacy_cellcycle),
                    "workflow_noncellcycle_input": str(legacy_noncellcycle),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )

            paths = self._figure7_input_lineage_paths(
                run_dir,
                manager_mode="standard",
            )

            self.assertEqual(paths.count(str(cellcycle)), 1)
            self.assertEqual(paths.count(str(noncellcycle)), 1)
            self.assertEqual(
                paths.count(
                    "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
                ),
                1,
            )
            self.assertNotIn("cell_ploidy.tsv", paths)
            self.assertNotIn(str(legacy_cellcycle), paths)
            self.assertNotIn(str(legacy_noncellcycle), paths)

    def test_input_lineage_uses_exact_explicit_workflow_seurat_rds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            explicit_rds = run_dir / "external inputs" / "selected_object.rds"
            explicit_rds.parent.mkdir()
            explicit_rds.write_text("fixture rds\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "state_pathway_support",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(run_dir / "unused_scvelo.csv"),
                    "workflow_seurat_rds": str(explicit_rds),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "a" * 64,
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertIn(str(explicit_rds), paths)
            self.assertNotIn(
                str(run_dir / "raw/integrated_sct_cca_seurat_final_reclustered.rds"),
                paths,
            )

    def test_input_lineage_records_present_reconstruction_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            upstream = run_dir / "seurat_upstream"
            final_stage = upstream / "03_final_cluster/stage_manifest.tsv"
            reconstruction = upstream / "reconstruction_manifest.tsv"
            final_stage.parent.mkdir(parents=True)
            final_stage.write_text("key\tvalue\nstage\tfinal\n")
            reconstruction.write_text("key\tvalue\nschema_version\t1\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(run_dir / "unused_scvelo.csv"),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "workflow_seurat_reconstruction_manifest": str(reconstruction),
                    "workflow_seurat_reconstruction_manifest_sha256": "a" * 64,
                    "workflow_seurat_final_stage_manifest": str(final_stage),
                    "workflow_seurat_final_stage_manifest_sha256": "b" * 64,
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertIn(str(reconstruction), paths)
            self.assertIn(str(final_stage), paths)
            self.assertIn(
                "Code/in-vivo/figure7/"
                "generate_final_seurat_from_cellranger.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/figure7/src/seurat_upstream.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/figure7/src/seurat_upstream_selection.R",
                paths,
            )

            reconstruction.unlink()
            final_stage.unlink()
            archived_paths = self._figure7_input_lineage_paths(run_dir)
            self.assertNotIn(str(reconstruction), archived_paths)
            self.assertNotIn(str(final_stage), archived_paths)
            self.assertIn(
                "Code/in-vivo/figure7/"
                "generate_final_seurat_from_cellranger.R",
                archived_paths,
            )

    def test_state_pathway_results_root_is_retained_for_historical_audit(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_results_root = tmp_path / "04i results"
            source_results_root.mkdir()
            canonical_reference_root = (
                tmp_path
                / "canonical Data"
                / "taoli_04i_etp2_24_day17_v1"
            )
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
if [[ "$entrypoint" == *export_04i_state_pathway_reference.R ]]; then
  mkdir -p "$output_dir"
  cp "$FIGURE7_TEST_REFERENCE_ROOT"/*.tsv "$output_dir/"
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
  printf 'key\\tvalue\\nmode\\tstandard\\npanel_set\\ta-f\\ntgi_day\\t17\\ncellcycle_input\\tData/in-vivo/figure7/processed/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv\\ncellcycle_sha256\\tfixture\\nnoncellcycle_input\\tData/in-vivo/figure7/processed/NonCellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv\\nnoncellcycle_sha256\\tfixture\\nstate_pathway_reference_id\\ttaoli_04i_etp2_24_day17_v1\\nstate_pathway_reference_kind\\treviewed_frozen\\ncanonical_publication_allowed\\ttrue\\nstate_pathway_source_results_root\\t%s\\n' \
    "$source_results_root" > "$output_dir/metadata/run_config.tsv"
  printf 'panel_id\\tfilename\\n' > "$output_dir/metadata/panel_contract.tsv"
  printf '7A\\tpanel_7A_day17_tgi_calculation.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7B\\tpanel_7B_cellcycle_selected_ecdf_comparisons.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7C\\tpanel_7C_day17_tgi_by_initial_ploidy.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7D\\tpanel_7D_day17_tgi_vs_centered_ecdf_shift.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7E\\tpanel_7E_day17_tgi_vs_mean_etp.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  printf '7F\\tpanel_7F_pseudotime_state_pathway_activity.pdf\\n' >> "$output_dir/metadata/panel_contract.tsv"
  for source in "$saved_state_dir"/*.tsv; do
    name="${source##*/}"
    if [[ "$name" == "state_pathway_provenance.tsv" ]]; then
      cp "$source" "$output_dir/metadata/$name"
    else
      cp "$source" "$output_dir/tables/$name"
    fi
  done
  exit 0
fi
exit 99
"""
            )
            fake_rscript.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"
            env["FIGURE7_CANONICAL_REFERENCE_ROOT"] = str(
                canonical_reference_root
            )
            env["FIGURE7_TEST_REFERENCE_ROOT"] = str(
                REPO_ROOT
                / "Data/in-vivo/figure7/saved_state_pathway"
                / "taoli_04i_etp2_24_day17_v1"
            )

            result = self._run(
                "--mode", "standard", "--modules", "in_vivo_figure7",
                "--run-id", run_id,
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                "--figure7-state-pathway-results-root",
                str(source_results_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)

            manager_root = output_root / "manager/runs" / run_id
            export_metadata_path = (
                manager_root
                / "metadata/figure7_state_pathway_export.tsv"
            )
            with export_metadata_path.open(newline="") as handle:
                export_metadata = {
                    row["key"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(export_metadata["status"], "ok")
            self.assertEqual(
                export_metadata["source_results_root"],
                str(source_results_root.resolve()),
            )
            reference_root = Path(
                export_metadata["exported_audit_reference_dir"]
            )
            self.assertEqual(
                reference_root.name,
                "taoli_04i_etp2_24_day17_v1",
            )
            self.assertEqual(len(list(reference_root.glob("*.tsv"))), 8)
            self.assertEqual(
                export_metadata["historical_reference_id"],
                "taoli_04i_etp2_24_day17_v1",
            )
            self.assertEqual(
                export_metadata["canonical_data_materialization"],
                "prohibited",
            )
            self.assertFalse(canonical_reference_root.exists())
            audit_provenance = (
                reference_root / "state_pathway_provenance.tsv"
            ).read_text()
            self.assertNotIn(
                str(source_results_root.resolve()),
                audit_provenance,
            )

            materialization_path = (
                manager_root
                / "metadata/figure7_state_pathway_materialization.tsv"
            )
            self.assertFalse(materialization_path.exists())

            run_root = (
                output_root
                / "in-vivo/figure7/runs"
                / f"{run_id}_figure7"
            )
            with (
                run_root / "metadata/run_config.tsv"
            ).open(newline="") as handle:
                run_config = {
                    row["key"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(
                run_config["state_pathway_source_results_root"],
                str(source_results_root.resolve()),
            )
            module_runs = (
                manager_root / "metadata/module_runs.tsv"
            ).read_text()
            self.assertIn(
                "state_pathway_source_results_root="
                f"{source_results_root.resolve()}",
                module_runs,
            )
            self.assertIn(
                "historical_mixed_v1_audit_only=true",
                module_runs,
            )
            self.assertIn(
                "canonical_data_materialization=prohibited",
                module_runs,
            )
            self.assertIn("ok_noncanonical", module_runs)
            self.assertFalse(
                (figure_root / "Figure7/manifest.tsv").exists()
            )

    def test_failed_figure7_run_does_not_materialize_canonical_reference(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_results_root = tmp_path / "04i_results"
            canonical_reference_root = (
                tmp_path
                / "canonical"
                / "taoli_04i_etp2_24_day17_v1"
            )
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
if [[ "$entrypoint" == *export_04i_state_pathway_reference.R ]]; then
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
            env["FIGURE7_CANONICAL_REFERENCE_ROOT"] = str(
                canonical_reference_root
            )

            result = self._run(
                "--mode", "standard", "--modules", "in_vivo_figure7",
                "--run-id", "failed_before_materialization",
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                "--figure7-state-pathway-results-root",
                str(source_results_root),
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
            (run_root / "tables").mkdir()
            self._copy_reviewed_reference(run_root)

            rows = []
            for spec in PANEL_SPECS:
                if (
                    spec["module"] != "in_vivo_figure7"
                    or str(spec["panel"]).startswith("7F")
                    or str(spec["panel"]) == "7A-7K_composite"
                ):
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
            design = [
                ("2N-A2-0", "2N", "30mg/kg", "30", "SUM159-2N-30-0_harvest.sps.cbs"),
                ("2N-A2-L", "2N", "30mg/kg", "30", "SUM159-2N-30-L_harvest.sps.cbs"),
                ("2N-A4-R", "2N", "120mg/kg", "120", "SUM159-2N-120-R_harvest.sps.cbs"),
                ("2N-A4-RL", "2N", "120mg/kg", "120", "SUM159-2N-120-RL_harvest.sps.cbs"),
                ("A6-4N-O", "4N", "30mg/kg", "30", "SUM159-4N-30-0_harvest.sps.cbs"),
                ("A6-4N-RR", "4N", "30mg/kg", "30", "SUM159-4N-30-RR_harvest.sps.cbs"),
                ("4N-A8-RL", "4N", "120mg/kg", "120", "SUM159-4N-120-RL_harvest.sps.cbs"),
                ("4N-A8-RR", "4N", "120mg/kg", "120", "SUM159-4N-120-RR_harvest.sps.cbs"),
            ]
            endpoint_path = Path(
                "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
            )
            cellcycle_path = Path(
                "Data/in-vivo/figure7/processed/"
                "CellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            )
            noncellcycle_path = Path(
                "Data/in-vivo/figure7/processed/"
                "NonCellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            )
            endpoint_score_by_key: dict[tuple[str, str], float] = {}
            endpoint_files: set[str] = set()
            with (REPO_ROOT / endpoint_path).open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    key = (row["file"], row["cell_id"])
                    self.assertNotIn(key, endpoint_score_by_key)
                    endpoint_score_by_key[key] = float(row["ploidy"])
                    endpoint_files.add(row["file"])

            curated_scores_by_sample: dict[str, list[float]] = {}
            curated_file_by_sample: dict[str, str] = {}
            curated_keys: set[tuple[str, str]] = set()
            compartment_counts: dict[str, int] = {}
            tgi_by_sample: dict[str, float] = {}
            for compartment, relative_path in (
                ("CellCycle", cellcycle_path),
                ("NonCellCycle", noncellcycle_path),
            ):
                compartment_count = 0
                with (REPO_ROOT / relative_path).open(newline="") as handle:
                    for row in csv.DictReader(handle):
                        sample_id = row["sample_id"]
                        prefix = f"{sample_id}_"
                        full_cell_id = row["cell_id"]
                        self.assertTrue(full_cell_id.startswith(prefix))
                        endpoint_file = f"{row['growth_curve_harvest']}.sps.cbs"
                        key = (endpoint_file, full_cell_id[len(prefix) :])
                        self.assertNotIn(key, curated_keys)
                        self.assertIn(key, endpoint_score_by_key)
                        canonical_score = endpoint_score_by_key[key]
                        self.assertLess(
                            abs(float(row["cell_ploidy"]) - canonical_score),
                            1e-12,
                        )
                        previous_file = curated_file_by_sample.setdefault(
                            sample_id,
                            endpoint_file,
                        )
                        self.assertEqual(previous_file, endpoint_file)
                        curated_keys.add(key)
                        curated_scores_by_sample.setdefault(sample_id, []).append(
                            canonical_score
                        )
                        compartment_count += 1
                        if compartment == "CellCycle":
                            tgi_by_sample.setdefault(
                                sample_id,
                                float(row["TGI_percent_Day_17"]),
                            )
                compartment_counts[compartment] = compartment_count

            expected_curated_counts = {
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
            observed_curated_counts = {
                sample_id: len(values)
                for sample_id, values in curated_scores_by_sample.items()
            }
            treated_sample_ids = {row[0] for row in design}
            self.assertEqual(len(endpoint_score_by_key), 14125)
            self.assertEqual(len(endpoint_files), 16)
            self.assertEqual(compartment_counts, {
                "CellCycle": 2881,
                "NonCellCycle": 6951,
            })
            self.assertEqual(observed_curated_counts, expected_curated_counts)
            self.assertEqual(len(curated_keys), 9832)
            self.assertEqual(
                sum(observed_curated_counts[sample_id] for sample_id in treated_sample_ids),
                5335,
            )
            scores = [
                statistics.mean(curated_scores_by_sample[sample_id])
                for sample_id, *_ in design
            ]
            z_scores = []
            for origin in ("2N", "4N"):
                indices = [
                    index
                    for index, row in enumerate(design)
                    if row[1] == origin
                ]
                local = [scores[index] for index in indices]
                center = statistics.mean(local)
                spread = statistics.stdev(local)
                z_scores.extend((value - center) / spread for value in local)
            panel_k_plot = run_root / "tables/panel_7E_plot_data.tsv"
            plot_rows = []
            for index, (
                sample_id,
                origin,
                dose,
                dose_mg,
                endpoint_file,
            ) in enumerate(design):
                score_z = z_scores[index]
                score = scores[index]
                plot_rows.append(
                    {
                        "sample_id": sample_id,
                        "initial_ploidy": origin,
                        "dose": dose,
                        "dose_mg": dose_mg,
                        "endpoint_ploidy_file": endpoint_file,
                        "sample_mean_endpoint_ploidy": score,
                        "n_endpoint_ploidy_cells": len(
                            curated_scores_by_sample[sample_id]
                        ),
                        "endpoint_ploidy_source_total_cells": "14125",
                        "endpoint_ploidy_score_universe_total_cells": "9832",
                        "endpoint_ploidy_source_file_count": "16",
                        "endpoint_ploidy_source_sha256": sha256_file(
                            REPO_ROOT / endpoint_path
                        ),
                        "endpoint_ploidy_score_policy": (
                            "arithmetic_mean_of_finite_postprocessed_cell_"
                            "ploidy_in_exact_qc_passed_cellcycle_"
                            "noncellcycle_union_per_sample"
                        ),
                        "endpoint_ploidy_mapping_policy": (
                            "exact_processed_sample_barcode_to_canonical_cbs_"
                            "file_cell_and_value;score_universe=reviewed_final_"
                            "seurat_tumor_cells"
                        ),
                        "terminal_postprocessed_cn_score": score,
                        "terminal_cn_score_within_origin_z": score_z,
                        "permutation_stratum": f"{origin}|{dose_mg}",
                        "terminal_cn_score_nuisance_residual": "pending",
                        "tgi_origin_dose_residual": "pending",
                        "TGI_percent_Day_17": tgi_by_sample[sample_id],
                        "tgi_outcome": "day",
                        "tgi_day": "17",
                        "tgi_measure": "TGI_percent_Day_17",
                        "matched_control_summary": "mean",
                        "matched_control_group": "initial_ploidy",
                    }
                )

            def additive_residual(column: str) -> list[float]:
                values = [float(row[column]) for row in plot_rows]
                grand = statistics.mean(values)
                origin_means = {
                    origin: statistics.mean(
                        value
                        for value, row in zip(values, plot_rows)
                        if row["initial_ploidy"] == origin
                    )
                    for origin in ("2N", "4N")
                }
                dose_means = {
                    dose_mg: statistics.mean(
                        value
                        for value, row in zip(values, plot_rows)
                        if row["dose_mg"] == dose_mg
                    )
                    for dose_mg in ("30", "120")
                }
                return [
                    value
                    - origin_means[row["initial_ploidy"]]
                    - dose_means[row["dose_mg"]]
                    + grand
                    for value, row in zip(values, plot_rows)
                ]

            x_residual = additive_residual(
                "terminal_cn_score_within_origin_z"
            )
            y_residual = additive_residual("TGI_percent_Day_17")
            for row, x_value, y_value in zip(
                plot_rows,
                x_residual,
                y_residual,
            ):
                row["terminal_cn_score_nuisance_residual"] = x_value
                row["tgi_origin_dose_residual"] = y_value
            write_tsv(panel_k_plot, plot_rows, list(plot_rows[0]))
            panel_k_test = run_root / "tables/panel_7E_test.tsv"
            test_row = {
                "n": "8",
                "estimate": "-0.3431355546551816",
                "partial_correlation": "-0.3431355546551816",
                "effect_per_within_origin_sd": "-8.792298734079608",
                "permutation_p_two_sided": "0.75",
                "n_permutations": "16",
                "permutation_mode": (
                    "exact_TGI_label_enumeration_within_initial_ploidy_x_dose"
                ),
                "permutation_strata": "initial_ploidy:dose_mg",
                "score_variable": (
                    "sample_mean_qc_passed_curated_cbs_cell_ploidy"
                ),
                "score_source_sha256": sha256_file(REPO_ROOT / endpoint_path),
                "score_source_n_cells": "9832",
                "score_inventory_n_cells": "14125",
                "score_source_n_files": "16",
                "treated_score_n_cells": "5335",
                "score_aggregation_policy": (
                    "arithmetic_mean_of_finite_postprocessed_cell_"
                    "ploidy_in_exact_qc_passed_cellcycle_"
                    "noncellcycle_union_per_sample"
                ),
                "sample_mapping_policy": (
                    "exact_processed_sample_barcode_to_canonical_cbs_"
                    "file_cell_and_value;score_universe=reviewed_final_"
                    "seurat_tumor_cells"
                ),
                "score_standardization": "z_score_within_initial_ploidy",
                "adjustment_terms": "initial_ploidy+dose_mg",
                "outcome_variable": "TGI_percent_Day_17",
                "plot_x": "terminal_cn_score_nuisance_residual",
                "plot_y": "tgi_origin_dose_residual",
                "tgi_outcome": "day",
                "tgi_day": "17",
                "tgi_measure": "TGI_percent_Day_17",
                "matched_control_summary": "mean",
                "matched_control_group": "initial_ploidy",
            }
            write_tsv(panel_k_test, [test_row], list(test_row))
            for source in (panel_k_plot, panel_k_test):
                rows.append(
                    {
                        "path": str(source),
                        "repo_relative_path": "",
                        "absolute_path": str(source),
                        "role": "output_table",
                        "source_kind": "generated_table",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(source),
                        "checksum_unavailable_reason": "",
                        "byte_size": source.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "panel K fixture",
                    }
                )
            write_tsv(run_root / "metadata/output_manifest.tsv", rows, MODULE_MANIFEST_COLUMNS)
            write_tsv(
                run_root / "metadata/run_config.tsv",
                [
                    {"key": "mode", "value": "standard"},
                    {"key": "panel_set", "value": "a-e"},
                    {"key": "tgi_day", "value": "17"},
                    {
                        "key": "cellcycle_input",
                        "value": str(cellcycle_path),
                    },
                    {
                        "key": "cellcycle_sha256",
                        "value": sha256_file(REPO_ROOT / cellcycle_path),
                    },
                    {
                        "key": "noncellcycle_input",
                        "value": str(noncellcycle_path),
                    },
                    {
                        "key": "noncellcycle_sha256",
                        "value": sha256_file(REPO_ROOT / noncellcycle_path),
                    },
                    {
                        "key": "endpoint_ploidy_input",
                        "value": str(endpoint_path),
                    },
                    {
                        "key": "endpoint_ploidy_sha256",
                        "value": sha256_file(REPO_ROOT / endpoint_path),
                    },
                    {
                        "key": "endpoint_ploidy_n_cells",
                        "value": "14125",
                    },
                    {
                        "key": "endpoint_ploidy_n_files",
                        "value": "16",
                    },
                    {
                        "key": "endpoint_ploidy_score_universe_n_cells",
                        "value": "9832",
                    },
                    {
                        "key": "endpoint_ploidy_treated_score_n_cells",
                        "value": "5335",
                    },
                    {
                        "key": "endpoint_ploidy_score_policy",
                        "value": (
                            "arithmetic_mean_of_finite_postprocessed_cell_"
                            "ploidy_in_exact_qc_passed_cellcycle_"
                            "noncellcycle_union_per_sample"
                        ),
                    },
                    {
                        "key": "endpoint_ploidy_mapping_policy",
                        "value": (
                            "exact_processed_sample_barcode_to_canonical_cbs_"
                            "file_cell_and_value;score_universe=reviewed_final_"
                            "seurat_tumor_cells"
                        ),
                    },
                ],
                ["key", "value"],
            )
            write_tsv(
                run_root / "metadata/panel_contract.tsv",
                [
                    {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                    for spec in PANEL_SPECS
                    if spec["module"] == "in_vivo_figure7"
                    and not str(spec["panel"]).startswith("7F")
                    and str(spec["panel"]) != "7A-7K_composite"
                    and spec.get("variant", "pdf") == "pdf"
                ],
                ["panel_id", "filename"],
            )
            processed_input_paths = [
                cellcycle_path,
                noncellcycle_path,
                endpoint_path,
                Path("Code/in-vivo/figure7/src/tgi_data.R"),
                Path("Code/in-vivo/figure7/src/tgi_statistics.R"),
                Path("Code/in-vivo/figure7/src/tgi_panels.R"),
            ]
            write_tsv(
                run_root / "metadata/input_manifest.tsv",
                [
                    {
                        "path": str(relative_path),
                        "repo_relative_path": str(relative_path),
                        "absolute_path": str(REPO_ROOT / relative_path),
                        "role": "input_data",
                        "source_kind": "input_file",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(REPO_ROOT / relative_path),
                        "checksum_unavailable_reason": "",
                        "byte_size": (
                            REPO_ROOT / relative_path
                        ).stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "fixture",
                    }
                    for relative_path in processed_input_paths
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
