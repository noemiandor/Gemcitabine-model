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
from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class ManagerFigure7CliTest(unittest.TestCase):
    def _run(self, *args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPO_ROOT / "Manager.sh"), *args],
            cwd=REPO_ROOT,
            env=env,
            text=True,
            capture_output=True,
        )

    def test_panels_only_requires_source_run_id(self) -> None:
        result = self._run(
            "--mode", "panels-only", "--modules", "in_vivo_figure7", "--run-id", "operation"
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--source-run-id is required", result.stderr)

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
        self.assertNotIn("--saved-state-pathway-dir", result.stdout)

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
            source_results_root = tmp_path / "04i results"
            source_results_root.mkdir()
            canonical_reference_root = tmp_path / "canonical Data" / "taoli_04i_etp2_24_day17_v1"
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
  printf 'key\\tvalue\\nexport_source_results_root\\t%s\\n' "$source_results_root" \
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
  printf 'key\\tvalue\\npanel_set\\ta-f\\nstate_pathway_source_results_root\\t%s\\n' \
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
            self.assertEqual(reference_root.name, "taoli_04i_etp2_24_day17_v1")
            self.assertEqual(len(list(reference_root.glob("*.tsv"))), 8)
            self.assertEqual(
                export_metadata["canonical_data_reference_dir"], str(canonical_reference_root)
            )
            self.assertEqual(len(list(canonical_reference_root.glob("*.tsv"))), 8)
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
            source_results_root = tmp_path / "04i_results"
            canonical_reference_root = tmp_path / "canonical" / "taoli_04i_etp2_24_day17_v1"
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
                [{"key": "panel_set", "value": "a-f"}],
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
