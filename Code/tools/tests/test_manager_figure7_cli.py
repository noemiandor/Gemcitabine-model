from __future__ import annotations

import csv
import hashlib
import os
import statistics
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"
TESTS_DIR = TOOLS_DIR / "tests"

import sys

sys.path.insert(0, str(TOOLS_DIR))
sys.path.insert(0, str(TESTS_DIR))
from figure_output_contract import MODULE_MANIFEST_COLUMNS, sha256_file, write_tsv  # noqa: E402
from figure7_density_fixture import write_density_localization_fixture  # noqa: E402
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    panel_specs_for_figure7_variant,
    validate_generated_candidate_source_run,
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
            / "state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4"
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

    def test_removed_include_in_vivo_flag_has_migration_error(self) -> None:
        result = self._run("--include-in-vivo")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--include-in-vivo was removed", result.stderr)
        self.assertIn("--modules in_vivo_figure7", result.stderr)

    def test_help_does_not_advertise_removed_in_vivo_flag(self) -> None:
        result = self._run("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("--include-in-vivo", result.stdout)

    def test_removed_in_vivo_module_has_migration_error(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo",
            "--run-id", "removed_in_vivo",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("Module in_vivo was removed", result.stderr)
        self.assertIn("TGI-AUC and CellCycle-subset ploidy analyses", result.stderr)
        self.assertIn("--modules in_vivo_figure7", result.stderr)

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

    def test_scrna_source_rejects_unknown_value(self) -> None:
        result = self._run("--figure7-scrna-source", "automatic")
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be rds or h5", result.stderr)

    def test_generated_candidate_requires_h5_full_refit_and_isolated_root(self) -> None:
        result = self._run("--publish-generated-candidate")
        self.assertEqual(result.returncode, 2)
        self.assertIn("requires --mode full-refit", result.stderr)

        result = self._run(
            "--mode", "full-refit",
            "--figure7-scrna-source", "h5",
            "--publish-generated-candidate",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("explicit isolated --figure-root", result.stderr)

        result = self._run(
            "--mode", "full-refit",
            "--figure7-scrna-source", "h5",
            "--publish-generated-candidate",
            "--figure-root", str(REPO_ROOT / "figures"),
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot write to the canonical figures root", result.stderr)

    def test_generated_candidate_specs_use_noncanonical_composite_name(self) -> None:
        specs = panel_specs_for_figure7_variant(
            17,
            "Figure7",
            "generated-candidate",
        )
        composite_assets = {
            str(spec["asset"])
            for spec in specs
            if str(spec["panel"]).startswith("7A-7L_composite")
        }
        self.assertEqual(
            composite_assets,
            {
                "Figure7_generated_GRCh_candidate.png",
                "Figure7_generated_GRCh_candidate.pdf",
            },
        )
        self.assertNotIn("Figure7_reviewed_GRCh.png", composite_assets)

    def test_rds_scrna_source_plans_full_sif_preflight_and_downstream(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            full_sif = tmp_path / "full.sif"
            full_sif.write_text("fixture full sif")
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            apptainer = fake_bin / "apptainer"
            apptainer.write_text("#!/usr/bin/env bash\nexit 99\n")
            apptainer.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
            env["CLUSTER_STANDALONE_POST_FILTER_SIF"] = str(full_sif)
            result = self._run(
                "--mode", "check-only",
                "--modules", "in_vivo_figure7",
                "--figure7-panels-ae-only",
                "--figure7-scrna-source", "rds",
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("[figure7_scrna_raw_preflight]", result.stdout)
            self.assertIn("--roles=loom\\,seurat_rds\\,support", result.stdout)
            self.assertNotIn("cellranger_h5", result.stdout)
            self.assertNotIn("[figure7_scrna_semantic_audit]", result.stdout)
            self.assertIn("integrated_sct_cca_seurat_final_reclustered.rds", result.stdout)
            self.assertIn("apptainer exec --cleanenv", result.stdout)

    def test_h5_scrna_source_plans_semantic_audit_before_downstream(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            full_sif = tmp_path / "full.sif"
            cluster_sif = tmp_path / "cluster.sif"
            full_sif.write_text("fixture full sif")
            cluster_sif.write_text("fixture cluster sif")
            cellranger_root = tmp_path / "A02_cellRanger"
            for index in range(18):
                sample = f"sample{index:02d}-Count-HM"
                outs = cellranger_root / sample / "outs"
                outs.mkdir(parents=True)
                (outs / f"{sample}_filtered_feature_bc_matrix.h5").write_text(
                    "fixture h5"
                )
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            apptainer = fake_bin / "apptainer"
            apptainer.write_text("#!/usr/bin/env bash\nexit 99\n")
            apptainer.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
            env["CLUSTER_STANDALONE_POST_FILTER_SIF"] = str(full_sif)
            env["CLUSTER_STANDALONE_PRE_FILTER_SIF"] = str(cluster_sif)
            result = self._run(
                "--mode",
                "check-only",
                "--modules",
                "in_vivo_figure7",
                "--figure7-panels-ae-only",
                "--figure7-scrna-source",
                "h5",
                "--figure7-cellranger-root",
                str(cellranger_root),
                "--figure7-raw-data-dir",
                str(tmp_path / "raw"),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("[figure7_scrna_h5]", result.stdout)
            self.assertIn("[figure7_scrna_semantic_audit]", result.stdout)
            self.assertIn(
                "Code/in-vivo/scRNA_Seq_analysis/audit_generated_seurat_rds.R",
                result.stdout,
            )
            self.assertLess(
                result.stdout.index("[figure7_scrna_semantic_audit]"),
                result.stdout.index("[in_vivo_figure7]"),
            )
            downstream_line = next(
                line
                for line in result.stdout.splitlines()
                if line.startswith("[in_vivo_figure7]")
            )
            self.assertIn(
                "scRNA_Seq_analysis/03_final_cluster/03_objects/"
                "integrated_sct_cca_seurat_final_reclustered.rds",
                downstream_line,
            )

    def test_h5_scrna_source_fails_preflight_when_canonical_h5_are_absent(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            full_sif = tmp_path / "full.sif"
            cluster_sif = tmp_path / "cluster.sif"
            full_sif.write_text("fixture full sif")
            cluster_sif.write_text("fixture cluster sif")
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            apptainer = fake_bin / "apptainer"
            apptainer.write_text("#!/usr/bin/env bash\nexit 99\n")
            apptainer.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
            env["CLUSTER_STANDALONE_POST_FILTER_SIF"] = str(full_sif)
            env["CLUSTER_STANDALONE_PRE_FILTER_SIF"] = str(cluster_sif)
            result = self._run(
                "--mode", "check-only",
                "--modules", "in_vivo_figure7",
                "--figure7-panels-ae-only",
                "--figure7-scrna-source", "h5",
                "--figure7-raw-data-dir", str(tmp_path / "raw"),
                env=env,
            )
            self.assertEqual(result.returncode, 1)
            self.assertIn(
                "Expected 18 downloaded and checksum-validated Cell Ranger H5 files",
                result.stderr,
            )
            self.assertIn("cellranger_h5", result.stdout)

    def test_generated_rds_binding_accepts_a_passed_audit_and_detects_tampering(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("portable_md5()")
        end = manager_text.index("\nprepare_figure7_scrna_source()", start)
        function_block = manager_text[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            generated = tmp_path / "generated.rds"
            generated.write_bytes(b"generated serialized fixture")
            marker = tmp_path / "AUDIT_COMPLETE.txt"
            generated_md5 = hashlib.md5(generated.read_bytes()).hexdigest()
            generated_sha256 = sha256_file(generated)
            marker.write_text(
                "\n".join(
                    [
                        "schema_version=semantic_rds_audit_v1",
                        "status=PASS",
                        f"generated_rds={generated.resolve()}",
                        f"generated_rds_size_bytes={generated.stat().st_size}",
                        f"generated_rds_md5={generated_md5}",
                        f"generated_rds_sha256={generated_sha256}",
                        f"downstream_rds={generated.resolve()}",
                    ]
                )
                + "\n"
            )
            script = f"""
set -euo pipefail
require_file() {{ [[ -f \"$1\" ]] || return 1; }}
{function_block}
figure7_scrna_source=h5
bind_figure7_audited_generated_rds \"$1\" \"$2\"
verify_figure7_audited_generated_rds
printf '%s\\n' \"$figure7_audited_generated_rds_sha256\"
"""
            passed = subprocess.run(
                ["bash", "-c", script, "fixture", str(generated), str(marker)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(passed.returncode, 0, passed.stderr)
            self.assertIn(generated_sha256, passed.stdout)

            generated.write_bytes(b"tampered serialized fixture")
            failed = subprocess.run(
                ["bash", "-c", script, "fixture", str(generated), str(marker)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("no longer matches", failed.stderr)

    def test_standalone_uses_explicit_inputs_and_sample_prefixed_h5(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        shell_text = (
            REPO_ROOT / "Code/in-vivo/scRNA_Seq_analysis/run_cluster_standalone.sh"
        ).read_text()
        pipeline_text = (
            REPO_ROOT / "Code/in-vivo/scRNA_Seq_analysis/cluster_pipeline_standalone.R"
        ).read_text()
        self.assertIn(
            "CELLRANGER_ROOT ALL_PLOIDY_TSV SAMPLE_INFO_XLSX OUTPUT_DIR",
            shell_text,
        )
        self.assertIn(
            'paste0(sample_folder, "_filtered_feature_bc_matrix.h5")',
            pipeline_text,
        )
        self.assertNotIn(
            'file.path(sample_dir, "outs", "filtered_feature_bc_matrix.h5")',
            pipeline_text,
        )
        self.assertIn(
            '${figure7_seurat_upstream_dir}/scRNA_Seq_analysis/00_provenance',
            manager_text,
        )
        self.assertIn(
            'scRNA_Seq_analysis/00_validation/rds_semantic_audit',
            manager_text,
        )
        self.assertIn(
            "Code/in-vivo/scRNA_Seq_analysis/audit_generated_seurat_rds.R",
            manager_text,
        )
        self.assertIn('figure7_active_scrna_rds="${figure7_audited_generated_rds:-${generated_rds}}"', manager_text)
        self.assertNotIn("write_figure7_rds_md5_gate", manager_text)

    def test_generated_candidate_validator_requires_bound_semantic_audit(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            source_id = "semantic_candidate"
            run_root = tmp_path / f"{source_id}_si_figures"
            metadata = run_root / "metadata"
            provenance = run_root / "upstream/00_provenance"
            audit_dir = run_root / "upstream/00_validation/rds_semantic_audit"
            metadata.mkdir(parents=True)
            provenance.mkdir(parents=True)
            audit_dir.mkdir(parents=True)

            generated = run_root / "generated.rds"
            reference = run_root / "zenodo_reference.rds"
            generated.write_bytes(b"generated RDS fixture")
            reference.write_bytes(b"different Zenodo RDS fixture")
            for name in (
                "run_manifest.tsv",
                "final_artifact_runtime.tsv",
                "PIPELINE_COMPLETE.txt",
            ):
                (provenance / name).write_text(f"fixture={name}\n")

            summary = audit_dir / "audit_summary.tsv"
            write_tsv(
                summary,
                [
                    {
                        "group": "object",
                        "check": "fixture",
                        "status": "PASS",
                        "observed": "equivalent",
                        "reference": "equivalent",
                        "threshold": "fixture",
                        "details": "fixture",
                    }
                ],
                [
                    "group",
                    "check",
                    "status",
                    "observed",
                    "reference",
                    "threshold",
                    "details",
                ],
            )
            report_names = {
                "audit_summary.tsv",
                "metadata_comparison.tsv",
                "cluster_comparison.tsv",
                "graph_comparison.tsv",
                "assay_numeric_comparison.tsv",
                "pca_comparison.tsv",
                "umap_displacement_summary.tsv",
                "umap_largest_displacements.tsv",
                "command_comparison.tsv",
                "rds_file_identity.tsv",
            }
            for name in report_names - {
                "audit_summary.tsv",
                "rds_file_identity.tsv",
            }:
                (audit_dir / name).write_text("check\tstatus\nfixture\tPASS\n")
            generated_md5 = hashlib.md5(generated.read_bytes()).hexdigest()
            reference_md5 = hashlib.md5(reference.read_bytes()).hexdigest()
            write_tsv(
                audit_dir / "rds_file_identity.tsv",
                [
                    {
                        "artifact": "generated",
                        "path": str(generated.resolve()),
                        "size_bytes": generated.stat().st_size,
                        "md5": generated_md5,
                        "sha256": sha256_file(generated),
                        "role": "candidate_downstream_input",
                    },
                    {
                        "artifact": "zenodo_reference",
                        "path": str(reference.resolve()),
                        "size_bytes": reference.stat().st_size,
                        "md5": reference_md5,
                        "sha256": sha256_file(reference),
                        "role": "semantic_reference",
                    },
                ],
                ["artifact", "path", "size_bytes", "md5", "sha256", "role"],
            )
            marker_rows = [
                "schema_version=semantic_rds_audit_v1",
                "status=PASS",
                "checks_total=1",
                "checks_passed=1",
                "checks_failed=0",
                f"generated_rds={generated.resolve()}",
                f"generated_rds_size_bytes={generated.stat().st_size}",
                f"generated_rds_md5={generated_md5}",
                f"generated_rds_sha256={sha256_file(generated)}",
                f"zenodo_reference_rds={reference.resolve()}",
                f"zenodo_reference_rds_size_bytes={reference.stat().st_size}",
                f"zenodo_reference_rds_md5={reference_md5}",
                f"zenodo_reference_rds_sha256={sha256_file(reference)}",
                f"downstream_rds={generated.resolve()}",
                f"audit_summary={summary.resolve()}",
                f"audit_summary_sha256={sha256_file(summary)}",
            ]
            marker_rows.extend(
                f"report_sha256.{name}={sha256_file(audit_dir / name)}"
                for name in sorted(report_names)
            )
            (audit_dir / "AUDIT_COMPLETE.txt").write_text(
                "\n".join(marker_rows) + "\n"
            )

            write_tsv(
                metadata / "run_config.tsv",
                [{"key": "si7_canonical_publication_allowed", "value": "false"}],
                ["key", "value"],
            )
            write_tsv(metadata / "panel_contract.tsv", [], ["panel_id", "filename"])
            output_artifact = run_root / "output.txt"
            output_artifact.write_text("fixture output\n")

            def manifest_row(path: Path, role: str) -> dict[str, object]:
                return {
                    "path": str(path.resolve()),
                    "repo_relative_path": "",
                    "absolute_path": str(path.resolve()),
                    "role": role,
                    "source_kind": "input_file" if role == "input_data" else "generated_table",
                    "module": "si_figures",
                    "generated_by": "test",
                    "command_id": source_id,
                    "sha256": sha256_file(path),
                    "checksum_unavailable_reason": "",
                    "byte_size": path.stat().st_size,
                    "mtime_utc": "2026-08-05T00:00:00+00:00",
                    "figure": "",
                    "panel": "",
                    "notes": "fixture",
                }

            input_paths = [generated, reference, *provenance.iterdir(), *audit_dir.iterdir()]
            write_tsv(
                metadata / "input_manifest.tsv",
                [manifest_row(path, "input_data") for path in input_paths],
                MODULE_MANIFEST_COLUMNS,
            )
            write_tsv(
                metadata / "output_manifest.tsv",
                [manifest_row(output_artifact, "output_table")],
                MODULE_MANIFEST_COLUMNS,
            )

            validate_generated_candidate_source_run(
                "si_figures",
                run_root,
                [],
                REPO_ROOT,
                source_id,
                17,
            )
            (audit_dir / "graph_comparison.tsv").write_text(
                "check\tstatus\nfixture\tFAIL\n"
            )
            with self.assertRaisesRegex(
                ValueError,
                "sha256 does not match|report was modified",
            ):
                validate_generated_candidate_source_run(
                    "si_figures",
                    run_root,
                    [],
                    REPO_ROOT,
                    source_id,
                    17,
                )

    def test_ae_only_rejects_full_panel_f_analysis(self) -> None:
        result = self._run(
            "--mode", "full-refit", "--modules", "in_vivo_figure7",
            "--run-id", "bad_combo", "--figure7-panels-ae-only",
            "--figure7-full-analysis", "--figure7-seurat-rds", "/tmp/object.rds",
            "--figure7-gene-set-artifact", "/tmp/gene_sets.rds",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("mutually exclusive", result.stderr)

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
                    or str(spec["panel"]).startswith("7A-7L_composite")
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
            panel_l_plot = run_root / "tables/panel_7E_plot_data.tsv"
            plot_rows = []
            for index, (
                sample_id,
                origin,
                dose,
                dose_mg,
                endpoint_file,
            ) in enumerate(design):
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
                        "TGI_percent_Day_17": tgi_by_sample[sample_id],
                        "tgi_outcome": "day",
                        "tgi_day": "17",
                        "tgi_measure": "TGI_percent_Day_17",
                        "matched_control_summary": "mean",
                        "matched_control_group": "initial_ploidy",
                    }
                )

            write_tsv(panel_l_plot, plot_rows, list(plot_rows[0]))
            panel_l_test = run_root / "tables/panel_7E_test.tsv"
            test_row = {
                "n": "8",
                "estimate": "-0.698401019253014",
                "asymptotic_p": "0.0540069781511847",
                "permutation_p_two_sided": "0.0591269841269841",
                "n_permutations": "40320",
                "permutation_mode": "exact_TGI_label_enumeration",
                "permutation_strata": "none",
                "association_type": "unadjusted_mouse_level_pearson",
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
                "score_standardization": "none",
                "adjustment_terms": "none",
                "outcome_variable": "TGI_percent_Day_17",
                "plot_x": "sample_mean_endpoint_ploidy",
                "plot_y": "TGI_percent_Day_17",
                "tgi_outcome": "day",
                "tgi_day": "17",
                "tgi_measure": "TGI_percent_Day_17",
                "matched_control_summary": "mean",
                "matched_control_group": "initial_ploidy",
            }
            write_tsv(panel_l_test, [test_row], list(test_row))
            for source in (panel_l_plot, panel_l_test):
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
                        "notes": "panel L fixture",
                    }
                )
            for source in write_density_localization_fixture(run_root):
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
                        "notes": "density-localization fixture",
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
                        "key": "config_sha256",
                        "value": sha256_file(
                            REPO_ROOT
                            / "Code/in-vivo/figure7/figure7_config.yaml"
                        ),
                    },
                    {
                        "key": "density_localization_config",
                        "value": (
                            "Code/in-vivo/figure7/"
                            "density_localization_config.yaml"
                        ),
                    },
                    {
                        "key": "density_localization_config_sha256",
                        "value": sha256_file(
                            REPO_ROOT
                            / "Code/in-vivo/figure7/"
                            "density_localization_config.yaml"
                        ),
                    },
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
                    and not str(spec["panel"]).startswith("7A-7L_composite")
                    and spec.get("variant", "pdf") == "pdf"
                ],
                ["panel_id", "filename"],
            )
            processed_input_paths = [
                cellcycle_path,
                noncellcycle_path,
                endpoint_path,
                Path("Code/in-vivo/figure7/run_figure7.R"),
                Path("Code/in-vivo/figure7/figure7_config.yaml"),
                Path("Code/in-vivo/figure7/density_localization_config.yaml"),
                Path("Code/in-vivo/figure7/src/common_io.R"),
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
