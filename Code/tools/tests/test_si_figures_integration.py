from __future__ import annotations

import csv
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"
sys.path.insert(0, str(TOOLS_DIR))

from figure_output_contract import (  # noqa: E402
    FIGURE_MANIFEST_COLUMNS,
    MODULE_MANIFEST_COLUMNS,
    sha256_file,
    write_tsv,
)
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    FIGURE7J_COPY_NUMBER_ANNOTATION_COLUMNS,
    SI7_REVIEWED_FEATURE_POLICY,
    SI7_REVIEWED_FROZEN_MATRIX_NOTE,
    SI7_REVIEWED_GENE_SET_DATABASE,
)


class SiFiguresManagerTest(unittest.TestCase):
    @staticmethod
    def _dry_run_commands(stdout: str) -> dict[str, list[str]]:
        commands: dict[str, list[str]] = {}
        for line in stdout.splitlines():
            if not line.startswith("[") or "] " not in line:
                continue
            label, command = line.split("] ", 1)
            commands[label.removeprefix("[")] = shlex.split(command)
        return commands

    @staticmethod
    def _figure7_required_input_paths(
        output_root: Path,
        run_id: str,
    ) -> list[str]:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("module_run_dir()")
        end = manager_text.index("\ncheck_module_inputs()", start)
        function_block = manager_text[start:end]
        script = f"""
set -euo pipefail
{function_block}
mode=full-refit
dry_run=false
figure7_full_analysis=false
figure7_panels_ae_only=false
figure7_state_pathway_results_root=""
figure7_reference_root=unused
figure7_cell_ploidy_input=Data/in-vivo/all_ploidy.tsv
figure7_sample_info_input=Data/in-vivo/sample_info.xlsx
figure7_growth_curve_input=Data/in-vivo/dt_Gem_VT_20241223_v4.xlsx
output_root="$1"
run_id="$2"
required_input_paths_for_module in_vivo_figure7
"""
        result = subprocess.run(
            [
                "bash",
                "-c",
                script,
                "fixture",
                str(output_root),
                run_id,
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )
        if result.returncode != 0:
            raise AssertionError(result.stderr)
        return [line for line in result.stdout.splitlines() if line]

    def test_default_modules_include_frozen_si_figures(self) -> None:
        default_line = next(
            line
            for line in (REPO_ROOT / "Manager.sh").read_text().splitlines()
            if line.startswith('modules="')
        )
        self.assertIn("si_figures", default_line)

    def test_reviewed_si_policy_matches_materialization_guard(self) -> None:
        config_text = (
            REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml"
        ).read_text()
        generator_text = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        ).read_text()
        self.assertIn(
            f'si7_feature_species_policy: "{SI7_REVIEWED_FEATURE_POLICY}"',
            config_text,
        )
        self.assertIn(
            "grch_human_only_v2_20260729_raw_refit_retry3_si_figures",
            generator_text,
        )
        self.assertNotIn(
            "Frozen matrices were recalculated from Tao's cluster DEG cache",
            generator_text,
        )

    def test_check_only_uses_frozen_cache_without_raw_workflow(self) -> None:
        result = subprocess.run(
            [
                "bash",
                str(REPO_ROOT / "Manager.sh"),
                "--mode",
                "check-only",
                "--modules",
                "si_figures",
                "--run-id",
                "check_si_figures",
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("run_supplementary_figures.R", result.stdout)
        self.assertIn("--mode=plot-only", result.stdout)
        self.assertIn("Data/in-vivo/SIfigures", result.stdout)
        self.assertNotIn("seurat-rds", result.stdout.lower())
        self.assertNotIn("force-reanalysis", result.stdout.lower())

    def test_missing_default_endpoint_table_uses_run_scoped_cbs_fallback(
        self,
    ) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("quote_args()")
        end = manager_text.index("\nregistry_module_name()", start)
        function_block = manager_text[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            derivation_script = tmp_path / "Data/in-vivo/weighted_ploidy.py"
            cbs_manifest = (
                tmp_path / "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"
            )
            derivation_script.parent.mkdir(parents=True)
            cbs_manifest.parent.mkdir(parents=True)
            derivation_script.write_text("# fixture\n")
            shutil.copy2(
                REPO_ROOT / "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
                cbs_manifest,
            )
            script = f"""
set -euo pipefail
{function_block}
module_list=(si_figures)
figure7_cell_ploidy_input=Data/in-vivo/all_ploidy.tsv
figure7_cell_ploidy_input_explicit=false
figure7_python=/fixture/python
si_figures_cbs_dir=Data/in-vivo/scRNAseq_Numbat
manager_run_dir=Results/manager/runs/fallback
mode=check-only
dry_run=false
prepare_endpoint_ploidy_if_missing
printf 'DERIVED=%s\\n' "${{figure7_cell_ploidy_input}}"
"""
            result = subprocess.run(
                ["bash", "-c", script],
                cwd=tmp_path,
                text=True,
                capture_output=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        expected = (
            "Results/manager/runs/fallback/artifacts/endpoint_ploidy/"
            "all_ploidy.tsv"
        )
        self.assertIn("[endpoint_ploidy_fallback]", result.stdout)
        self.assertIn("--manifest Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv", result.stdout)
        self.assertIn("--omit-total-chromosomes", result.stdout)
        self.assertIn(
            "--expected-sha256 "
            "6db48ee5f196b37b58aa71d0472dd3deb06aaacb4b637070af1b27d9425db2b3",
            result.stdout,
        )
        self.assertIn(f"DERIVED={expected}", result.stdout)

    def test_endpoint_fallback_rejects_tampered_cbs_manifest(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("quote_args()")
        end = manager_text.index("\nregistry_module_name()", start)
        function_block = manager_text[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            derivation_script = tmp_path / "Data/in-vivo/weighted_ploidy.py"
            cbs_manifest = (
                tmp_path / "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"
            )
            derivation_script.parent.mkdir(parents=True)
            cbs_manifest.parent.mkdir(parents=True)
            derivation_script.write_text("# fixture\n")
            cbs_manifest.write_text(
                "filename\tbytes\tsha256\tnotes\n"
                "tampered.sps.cbs\t1\tdeadbeef\ttampered\n"
            )
            script = f"""
set -euo pipefail
{function_block}
module_list=(si_figures)
figure7_cell_ploidy_input=Data/in-vivo/all_ploidy.tsv
figure7_cell_ploidy_input_explicit=false
figure7_python=/fixture/python
si_figures_cbs_dir=Data/in-vivo/scRNAseq_Numbat
manager_run_dir=Results/manager/runs/fallback
mode=check-only
dry_run=false
prepare_endpoint_ploidy_if_missing
"""
            result = subprocess.run(
                ["bash", "-c", script],
                cwd=tmp_path,
                text=True,
                capture_output=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("requires the reviewed CBS manifest", result.stderr)

    def test_generator_has_no_raw_analysis_entrypoint(self) -> None:
        generator = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        ).read_text()
        for forbidden in (
            "Seurat::",
            "msigdbr",
            "fgsea",
            "force-reanalysis",
            "seurat-rds",
            "deg-cache",
        ):
            self.assertNotIn(forbidden, generator)

    def test_generated_renderer_requires_and_labels_raw_lineage(self) -> None:
        renderer = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        ).read_text()
        self.assertIn(
            "allow_generated_human_only_si7 && "
            "is.na(upstream_manifest_path)",
            renderer,
        )
        for generated_label in (
            "si_figures_generated_cache_manifest",
            "si_figures_generated_table",
            "generated_table_source_revision",
            "si7_generated_matrix_note",
            "Mode: generated human-only plot-facing tables (noncanonical)",
        ):
            self.assertIn(generated_label, renderer)

    def test_full_refit_wires_only_figure7_and_si4_7_raw_stages(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "Results"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "full-refit",
                    "--modules",
                    "in_vivo_figure7,si_figures",
                    "--run-id",
                    "raw_figure7_contract",
                    "--output-root",
                    str(output_root),
                    "--figure-root",
                    str(Path(tmp) / "figures"),
                    "--no-update-latest",
                    "--dry-run",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--mode=full-workflow", result.stdout)
        self.assertNotIn("--mode=prepare-scvelo-inputs", result.stdout)
        self.assertIn("run_supplementary_figures.R", result.stdout)
        self.assertIn("build_raw_supplementary_tables.R", (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/run_supplementary_figures.R"
        ).read_text())
        self.assertIn("--download-missing-raw=true", result.stdout)
        self.assertIn(
            "--cell-ploidy-input=Data/in-vivo/all_ploidy.tsv",
            result.stdout,
        )
        self.assertNotIn("[figure1]", result.stdout)
        self.assertNotIn("[figure2]", result.stdout)
        commands = self._dry_run_commands(result.stdout)
        expected_upstream = (
            output_root
            / "in-vivo/figure7/intermediates/seurat_upstream"
        )
        for module in ("in_vivo_figure7", "si_figures"):
            self.assertIn(
                f"--seurat-upstream-dir={expected_upstream}",
                commands[module],
            )
        self.assertIn(
            "--sample-info=Data/in-vivo/sample_info.xlsx",
            commands["si_figures"],
        )
        generated_si_run = (
            output_root
            / "in-vivo/SI_figures/runs"
            / "raw_figure7_contract_si_figures"
        )
        figure7_command = commands["in_vivo_figure7"]
        self.assertIn(
            f"--si-table-cache-dir={generated_si_run / 'tables'}",
            figure7_command,
        )
        self.assertIn(
            "--si-cache-policy=generated-human-only",
            figure7_command,
        )
        self.assertIn(
            "--si-cache-upstream-input-manifest="
            f"{generated_si_run / 'metadata/analysis_input_manifest.tsv'}",
            figure7_command,
        )
        self.assertIn(
            "--si-cache-source-run-config="
            f"{generated_si_run / 'metadata/run_config.tsv'}",
            figure7_command,
        )
        self.assertIn(
            "--si-cache-source-provenance="
            f"{generated_si_run / 'metadata/si_figures_provenance.tsv'}",
            figure7_command,
        )
        self.assertFalse(
            any("Data/in-vivo/SIfigures" in arg for arg in figure7_command)
        )
        self.assertLess(
            result.stdout.index("[si_figures]"),
            result.stdout.index("[in_vivo_figure7]"),
        )

    def test_figure7_only_full_refit_adds_run_scoped_si_prerequisite(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            run_id = "figure7_with_automatic_si"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "full-refit",
                    "--modules",
                    "in_vivo_figure7",
                    "--run-id",
                    run_id,
                    "--output-root",
                    str(output_root),
                    "--figure-root",
                    str(tmp_path / "figures"),
                    "--no-update-latest",
                    "--dry-run",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        commands = self._dry_run_commands(result.stdout)
        self.assertEqual(
            [
                module
                for module in commands
                if module in {"si_figures", "in_vivo_figure7"}
            ],
            ["si_figures", "in_vivo_figure7"],
        )
        generated_si_run = (
            output_root
            / "in-vivo/SI_figures/runs"
            / f"{run_id}_si_figures"
        )
        figure7_command = commands["in_vivo_figure7"]
        self.assertIn(
            f"--si-table-cache-dir={generated_si_run / 'tables'}",
            figure7_command,
        )
        self.assertIn(
            "--si-cache-policy=generated-human-only",
            figure7_command,
        )

    def test_full_refit_figure7_required_inputs_exclude_reviewed_si_cache(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "Results"
            run_id = "generated_si_required_inputs"
            paths = self._figure7_required_input_paths(output_root, run_id)

        generated_si_run = (
            output_root
            / "in-vivo/SI_figures/runs"
            / f"{run_id}_si_figures"
        )
        self.assertIn(str(generated_si_run / "tables/manifest.tsv"), paths)
        self.assertIn(
            str(generated_si_run / "metadata/analysis_input_manifest.tsv"),
            paths,
        )
        self.assertFalse(
            any("Data/in-vivo/SIfigures" in path for path in paths)
        )

    def test_generated_figure7_input_paths_recover_external_run_scoped_si_lineage(
        self,
    ) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("module_run_dir()")
        end = manager_text.index(
            "\nrequired_input_paths_for_module()",
            start,
        )
        function_block = manager_text[start:end]

        with tempfile.TemporaryDirectory() as tmp:
            output_root = Path(tmp) / "external-results"
            run_id = "external_generated_si_handoff"
            figure7_run = (
                output_root
                / "in-vivo/figure7/runs"
                / f"{run_id}_figure7"
            )
            si_run = (
                output_root
                / "in-vivo/SI_figures/runs"
                / f"{run_id}_si_figures"
            )
            figure7_metadata = figure7_run / "metadata"
            cell_table_root = (
                output_root
                / "in-vivo/figure7/intermediates/celllevel_inputs"
            )
            cellcycle = cell_table_root / "cellcycle.csv"
            noncellcycle = cell_table_root / "noncellcycle.csv"
            si_metadata = si_run / "metadata"
            si_tables = si_run / "tables"
            figure7_metadata.mkdir(parents=True)
            cell_table_root.mkdir(parents=True)
            si_metadata.mkdir(parents=True)
            si_tables.mkdir(parents=True)
            cellcycle.write_text("cell\n")
            noncellcycle.write_text("cell\n")

            table_names = ["context.csv", "gsea.tsv"]
            (si_tables / "manifest.tsv").write_text(
                "filename\tbytes\tsha256\tsource_revision\tnotes\n"
                + "".join(
                    f"{name}\t1\t{'a' * 64}\traw-fixture\tgenerated\n"
                    for name in table_names
                )
            )
            for name in table_names:
                (si_tables / name).write_text("x")
            lineage_names = [
                "analysis_input_manifest.tsv",
                "run_config.tsv",
                "si_figures_provenance.tsv",
            ]
            for name in lineage_names:
                (si_metadata / name).write_text("fixture\n")

            write_tsv(
                figure7_metadata / "run_config.tsv",
                [
                    {"key": "mode", "value": "full-workflow"},
                    {
                        "key": "si_context_cache_policy",
                        "value": "generated-human-only",
                    },
                    {
                        "key": "si_context_cache_manifest",
                        "value": "external:manifest.tsv",
                    },
                    {
                        "key": "si_context_upstream_input_manifest",
                        "value": "external:analysis_input_manifest.tsv",
                    },
                    {
                        "key": "si_context_source_run_config",
                        "value": "external:run_config.tsv",
                    },
                    {
                        "key": "si_context_source_provenance",
                        "value": "external:si_figures_provenance.tsv",
                    },
                    {
                        "key": "cellcycle_input",
                        "value": "external:cellcycle.csv",
                    },
                    {
                        "key": "noncellcycle_input",
                        "value": "external:noncellcycle.csv",
                    },
                    {
                        "key": "workflow_cellcycle_input",
                        "value": str(cellcycle),
                    },
                    {
                        "key": "workflow_noncellcycle_input",
                        "value": str(noncellcycle),
                    },
                    {"key": "workflow_executed_stages", "value": "none"},
                    {
                        "key": "workflow_state_pathway_reference",
                        "value": "not_applicable",
                    },
                    {
                        "key": "workflow_state_pathway_results",
                        "value": "not_applicable",
                    },
                    {"key": "raw_data_dir", "value": "not_applicable"},
                    {
                        "key": "seurat_rds_sha256",
                        "value": "not_available",
                    },
                    {"key": "loom_file_count", "value": "0"},
                ],
                ["key", "value"],
            )

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
mode=full-refit
figure7_full_analysis=false
figure7_panels_ae_only=false
figure7_state_pathway_results_root=""
figure7_reference_root=unused
figure7_raw_data_dir=unused
figure7_loom_root=""
figure7_cell_ploidy_input=unused
figure7_sample_info_input=unused
figure7_growth_curve_input=unused
figure7_seurat_rds=unused
figure7_gene_set_artifact=unused
output_root="$1"
run_id="$2"
input_paths_for_module in_vivo_figure7 "$3"
"""
            result = subprocess.run(
                [
                    "bash",
                    "-c",
                    script,
                    "fixture",
                    str(output_root),
                    run_id,
                    str(figure7_run),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            paths = set(result.stdout.splitlines())
            expected_si_paths = {
                str(si_tables / "manifest.tsv"),
                *(str(si_tables / name) for name in table_names),
                *(str(si_metadata / name) for name in lineage_names),
                str(cellcycle),
                str(noncellcycle),
            }
            self.assertTrue(expected_si_paths.issubset(paths))
            self.assertFalse(any(path.startswith("external:") for path in paths))

    def test_full_refit_requires_raw_lineage_despite_reviewed_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            missing_ploidy = tmp_path / "absent-all-ploidy.tsv"
            missing_sample_info = tmp_path / "absent-sample-info.xlsx"
            missing_seurat = tmp_path / "absent-deposited-seurat.rds"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "full-refit",
                    "--modules",
                    "si_figures",
                    "--run-id",
                    "cache_first_si_precheck",
                    "--output-root",
                    str(tmp_path / "Results"),
                    "--figure-root",
                    str(tmp_path / "figures"),
                    "--figure7-cell-ploidy-input",
                    str(missing_ploidy),
                    "--figure7-sample-info-input",
                    str(missing_sample_info),
                    "--figure7-raw-seurat-rds",
                    str(missing_seurat),
                    "--no-update-latest",
                    "--dry-run",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            f"Missing explicitly selected endpoint-ploidy input: "
            f"{missing_ploidy}",
            result.stderr,
        )

    def test_explicit_shared_upstream_and_cellranger_are_forwarded_to_both(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            shared_upstream = tmp_path / "shared upstream"
            cellranger_root = tmp_path / "cell ranger matrices"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "full-refit",
                    "--modules",
                    "in_vivo_figure7,si_figures",
                    "--run-id",
                    "explicit_shared_upstream",
                    "--output-root",
                    str(tmp_path / "Results"),
                    "--figure-root",
                    str(tmp_path / "figures"),
                    "--figure7-seurat-upstream-dir",
                    str(shared_upstream),
                    "--figure7-cellranger-root",
                    str(cellranger_root),
                    "--dry-run",
                    "--no-update-latest",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        commands = self._dry_run_commands(result.stdout)
        for module in ("in_vivo_figure7", "si_figures"):
            self.assertIn(
                f"--seurat-upstream-dir={shared_upstream}",
                commands[module],
            )
            self.assertIn(
                f"--cellranger-root={cellranger_root}",
                commands[module],
            )

    def test_si_only_full_refit_uses_shared_figure7_upstream_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            shared_upstream = tmp_path / "shared"
            cellranger_root = tmp_path / "cellranger"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "full-refit",
                    "--modules",
                    "si_figures",
                    "--run-id",
                    "si_only_shared_upstream",
                    "--output-root",
                    str(tmp_path / "Results"),
                    "--figure7-seurat-upstream-dir",
                    str(shared_upstream),
                    "--figure7-cellranger-root",
                    str(cellranger_root),
                    "--dry-run",
                    "--no-update-latest",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        commands = self._dry_run_commands(result.stdout)
        self.assertNotIn("in_vivo_figure7", commands)
        self.assertIn(
            f"--seurat-upstream-dir={shared_upstream}",
            commands["si_figures"],
        )
        self.assertIn(
            f"--cellranger-root={cellranger_root}",
            commands["si_figures"],
        )
        self.assertIn(
            "--sample-info=Data/in-vivo/sample_info.xlsx",
            commands["si_figures"],
        )

    def test_manager_retains_and_binds_raw_analysis_lineage(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("input_paths_for_module()")
        end = manager_text.index(
            "\nrequired_input_paths_for_module()",
            start,
        )
        function_block = manager_text[start:end]
        deposited_hash = "a" * 64
        contract_hash = "b" * 64
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "raw_si_run"
            metadata = run_root / "metadata"
            metadata.mkdir(parents=True)
            authored = metadata / "input_manifest.tsv"
            with authored.open("w", newline="") as handle:
                writer = csv.writer(
                    handle,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writerow(
                    ["role", "repo_relative_path", "sha256", "bytes"]
                )
                writer.writerow(
                    [
                        "raw_build_deposited_seurat_rds",
                        "external:deposited.rds",
                        deposited_hash,
                        "123",
                    ]
                )
                writer.writerow(
                    [
                        "raw_build_si_raw_config_contract",
                        "contract:si_raw_config_contract",
                        contract_hash,
                        "0",
                    ]
                )
                writer.writerow(
                    [
                        "raw_build_audit_reconstruction_manifest",
                        "Results/in-vivo/figure7/intermediates/"
                        "seurat_upstream/reconstruction_manifest.tsv",
                        "c" * 64,
                        "123",
                    ]
                )
            (metadata / "run_config.tsv").write_text(
                "key\tvalue\n"
                "si7_canonical_publication_allowed\tfalse\n"
            )
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
mode=full-refit
figure7_raw_seurat_rds=""
retain_module_analysis_manifest si_figures "$1"
input_paths_for_module si_figures "$1"
"""
            result = subprocess.run(
                ["bash", "-c", script, "fixture", str(run_root)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            retained = metadata / "analysis_input_manifest.tsv"
            self.assertFalse(authored.exists())
            self.assertTrue(retained.is_file())
            self.assertIn(deposited_hash, retained.read_text())
            paths = [
                line for line in result.stdout.splitlines() if line
            ]
            self.assertIn(
                "Code/in-vivo/SI_figures/run_supplementary_figures.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/SI_figures/generate_supplementary_figures.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/SI_figures/normalized_composition.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/SI_figures/shared_context_panels.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/SI_figures/copy_number_heatmap.R",
                paths,
            )
            self.assertIn("Data/in-vivo/all_ploidy.tsv", paths)
            self.assertIn(
                "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
                paths,
            )
            self.assertIn(
                "Data/in-vivo/scRNAseq_Numbat/injected_reference/"
                "reference_manifest.tsv",
                paths,
            )
            self.assertIn(
                "Code/tools/validate_si_figures_table_cache.py",
                paths,
            )
            self.assertIn(str(retained), paths)
            self.assertNotIn(
                "Results/in-vivo/figure7/intermediates/"
                "seurat_upstream/reconstruction_manifest.tsv",
                paths,
            )
            expected_cbs = {
                str(path.relative_to(REPO_ROOT))
                for path in (
                    REPO_ROOT / "Data/in-vivo/scRNAseq_Numbat"
                ).glob("*.sps.cbs")
            }
            self.assertEqual(len(expected_cbs), 16)
            expected_injected_references = {
                str(path.relative_to(REPO_ROOT))
                for path in (
                    REPO_ROOT
                    / "Data/in-vivo/scRNAseq_Numbat/injected_reference"
                ).glob("*.sps.cbs")
            }
            self.assertEqual(len(expected_injected_references), 2)
            self.assertEqual(
                set(paths),
                {
                    "Code/in-vivo/SI_figures/"
                    "run_supplementary_figures.R",
                    "Code/in-vivo/SI_figures/"
                    "generate_supplementary_figures.R",
                    "Code/in-vivo/SI_figures/"
                    "normalized_composition.R",
                    "Code/in-vivo/SI_figures/"
                    "shared_context_panels.R",
                    "Code/in-vivo/SI_figures/"
                    "copy_number_heatmap.R",
                    "Code/in-vivo/figure7/src/common_io.R",
                    "Code/in-vivo/figure7/src/tgi_statistics.R",
                    "Code/in-vivo/figure7/src/tgi_panels.R",
                    "Code/in-vivo/figure7/density_localization_config.yaml",
                    "Data/in-vivo/figure7/processed/"
                    "CellCycleCells_pseudotime_distribution_per_sample_"
                    "cell_level_with_ploidy_dose_tgi.csv",
                    "Data/in-vivo/weighted_ploidy.py",
                    "Code/tools/validate_si_figures_table_cache.py",
                    "Data/in-vivo/all_ploidy.tsv",
                    "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
                    "Data/in-vivo/scRNAseq_Numbat/injected_reference/"
                    "reference_manifest.tsv",
                    str(retained),
                    *expected_cbs,
                    *expected_injected_references,
                },
            )

            manager_manifest = metadata / "input_manifest.tsv"
            command = [
                sys.executable,
                str(REPO_ROOT / "Code/tools/write_file_manifest.py"),
                "--manifest-type",
                "input",
                "--module",
                "si_figures",
                "--output",
                str(manager_manifest),
                "--generated-by",
                "Manager.sh",
                "--command-id",
                "raw_si_run",
                "--portable",
            ]
            for path in paths:
                command.extend(["--path", path])
            written = subprocess.run(
                command,
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(written.returncode, 0, written.stderr)
            with manager_manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            retained_rows = [
                row for row in rows
                if row["sha256"] == sha256_file(retained)
            ]
            self.assertEqual(len(retained_rows), 1)
            self.assertEqual(
                retained_rows[0]["path"],
                "external:analysis_input_manifest.tsv",
            )
            self.assertIn(deposited_hash, retained.read_text())

            # Rebuild that Manager manifest inside a checkout-shaped tree,
            # deliberately omit the raw-build ancestor named by the retained
            # analysis manifest, relocate the tree, and validate there.
            moved_repo = Path(tmp) / "moved-checkout"
            moved_paths: list[Path] = []
            for relative in (
                Path(
                    "Code/in-vivo/SI_figures/"
                    "run_supplementary_figures.R"
                ),
                Path(
                    "Code/in-vivo/SI_figures/"
                    "generate_supplementary_figures.R"
                ),
                Path(
                    "Code/in-vivo/SI_figures/"
                    "normalized_composition.R"
                ),
                Path(
                    "Code/in-vivo/SI_figures/"
                    "shared_context_panels.R"
                ),
                Path("Code/tools/validate_si_figures_table_cache.py"),
            ):
                destination = moved_repo / relative
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(REPO_ROOT / relative, destination)
                moved_paths.append(destination)
            moved_retained = (
                moved_repo
                / "Results/in-vivo/SI_figures/runs/raw_si_run/metadata"
                / "analysis_input_manifest.tsv"
            )
            moved_retained.parent.mkdir(parents=True)
            shutil.copy2(retained, moved_retained)
            moved_paths.append(moved_retained)
            moved_manager_manifest = (
                moved_retained.parent / "input_manifest.tsv"
            )
            moved_command = [
                sys.executable,
                str(REPO_ROOT / "Code/tools/write_file_manifest.py"),
                "--manifest-type",
                "input",
                "--module",
                "si_figures",
                "--output",
                str(moved_manager_manifest),
                "--generated-by",
                "Manager.sh",
                "--command-id",
                "raw_si_run",
                "--portable",
                "--repo-root",
                str(moved_repo),
            ]
            for path in moved_paths:
                moved_command.extend(["--path", str(path)])
            moved_written = subprocess.run(
                moved_command,
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(
                moved_written.returncode,
                0,
                moved_written.stderr,
            )
            relocated_repo = Path(tmp) / "relocated-checkout"
            moved_repo.rename(relocated_repo)
            relocated_manifest = (
                relocated_repo
                / moved_manager_manifest.relative_to(moved_repo)
            )
            validated = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "Code/tools/validate_manifest.py"),
                    str(relocated_manifest),
                    "--repo-root",
                    str(relocated_repo),
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(
                validated.returncode,
                0,
                validated.stderr,
            )
            self.assertFalse(
                (
                    relocated_repo
                    / "Results/in-vivo/figure7/intermediates/"
                    "seurat_upstream/reconstruction_manifest.tsv"
                ).exists()
            )

    def test_manager_rejects_missing_retained_si_analysis_manifest(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("input_paths_for_module()")
        end = manager_text.index(
            "\nrequired_input_paths_for_module()",
            start,
        )
        function_block = manager_text[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "missing_analysis_manifest"
            (run_root / "metadata").mkdir(parents=True)
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
mode=full-refit
figure7_raw_seurat_rds=""
input_paths_for_module si_figures "$1"
"""
            result = subprocess.run(
                ["bash", "-c", script, "fixture", str(run_root)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Missing retained SI analysis input manifest",
            result.stderr,
        )

    def test_manager_stops_when_si_input_path_enumeration_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            fake_rscript = fake_bin / "Rscript"
            fake_rscript.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
output_dir=""
for argument in "$@"; do
  case "${argument}" in
    --output-dir=*) output_dir="${argument#*=}" ;;
  esac
done
if [[ -z "${output_dir}" ]]; then
  echo "Fake Rscript did not receive --output-dir" >&2
  exit 2
fi
mkdir -p "${output_dir}/metadata"
printf 'role\\trepo_relative_path\\tsha256\\tbytes\\n' \\
  > "${output_dir}/metadata/input_manifest.tsv"
printf 'figure7_config\\tCode/in-vivo/figure7/figure7_config.yaml\\t%s\\t1\\n' \\
  "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" \\
  >> "${output_dir}/metadata/input_manifest.tsv"
printf 'key\\tvalue\\nsi7_canonical_publication_allowed\\ttrue\\n' \\
  > "${output_dir}/metadata/run_config.tsv"
"""
            )
            fake_rscript.chmod(0o755)
            output_root = tmp_path / "Results"
            run_id = "bad_si_manifest_enumeration"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}{os.pathsep}{env['PATH']}"
            result = subprocess.run(
                [
                    "bash",
                    str(REPO_ROOT / "Manager.sh"),
                    "--mode",
                    "standard",
                    "--modules",
                    "si_figures",
                    "--run-id",
                    run_id,
                    "--output-root",
                    str(output_root),
                    "--figure-root",
                    str(tmp_path / "figures"),
                    "--no-update-latest",
                ],
                cwd=REPO_ROOT,
                env=env,
                text=True,
                capture_output=True,
            )
            metadata = (
                output_root
                / "in-vivo/SI_figures/runs"
                / f"{run_id}_si_figures/metadata"
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Canonical SI analysis manifest has an invalid "
                "reviewed-cache binding",
                result.stderr,
            )
            self.assertIn(
                "Failed to resolve input manifest paths for module: "
                "si_figures",
                result.stderr,
            )
            self.assertTrue(
                (metadata / "analysis_input_manifest.tsv").is_file()
            )
            self.assertFalse((metadata / "input_manifest.tsv").exists())
            self.assertFalse((metadata / "output_manifest.tsv").exists())

    def test_manager_expands_only_reviewed_cache_rows_for_canonical_si(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("input_paths_for_module()")
        end = manager_text.index(
            "\nrequired_input_paths_for_module()",
            start,
        )
        function_block = manager_text[start:end]
        with tempfile.TemporaryDirectory() as tmp:
            run_root = Path(tmp) / "canonical_si"
            metadata = run_root / "metadata"
            metadata.mkdir(parents=True)
            (metadata / "run_config.tsv").write_text(
                "key\tvalue\n"
                "si7_canonical_publication_allowed\ttrue\n"
            )
            with (
                REPO_ROOT / "Data/in-vivo/SIfigures/manifest.tsv"
            ).open(newline="") as handle:
                frozen_locators = [
                    "Data/in-vivo/SIfigures/" + row["filename"]
                    for row in csv.DictReader(handle, delimiter="\t")
                ]
            self.assertEqual(len(frozen_locators), 11)
            rows = [
                (
                    "figure7_config",
                    "Code/in-vivo/figure7/figure7_config.yaml",
                ),
                (
                    "normalized_composition_helper",
                    "Code/in-vivo/SI_figures/normalized_composition.R",
                ),
                (
                    "shared_context_panels_helper",
                    "Code/in-vivo/SI_figures/shared_context_panels.R",
                ),
                (
                    "copy_number_heatmap_helper",
                    "Code/in-vivo/SI_figures/copy_number_heatmap.R",
                ),
                (
                    "figure7_common_io_helper",
                    "Code/in-vivo/figure7/src/common_io.R",
                ),
                (
                    "figure7_tgi_statistics_helper",
                    "Code/in-vivo/figure7/src/tgi_statistics.R",
                ),
                (
                    "figure7_tgi_panels_helper",
                    "Code/in-vivo/figure7/src/tgi_panels.R",
                ),
                (
                    "density_localization_config",
                    "Code/in-vivo/figure7/density_localization_config.yaml",
                ),
                (
                    "si4i_cellcycle_pseudotime",
                    "Data/in-vivo/figure7/processed/"
                    "CellCycleCells_pseudotime_distribution_per_sample_"
                    "cell_level_with_ploidy_dose_tgi.csv",
                ),
                (
                    "endpoint_ploidy_source",
                    "Data/in-vivo/all_ploidy.tsv",
                ),
                (
                    "numbat_cbs_manifest",
                    "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
                ),
                (
                    "si_figures_cache_manifest",
                    "Data/in-vivo/SIfigures/manifest.tsv",
                ),
            ] + [
                (
                    "numbat_cbs_matrix_"
                    + "".join(
                        character if character.isalnum() else "_"
                        for character in path.name
                    ),
                    str(path.relative_to(REPO_ROOT)),
                )
                for path in sorted(
                    (
                        REPO_ROOT / "Data/in-vivo/scRNAseq_Numbat"
                    ).glob("*.sps.cbs")
                )
            ] + [
                (
                    "si_figures_frozen_table",
                    locator,
                )
                for locator in frozen_locators
            ] + [
                (
                    "raw_build_audit_reconstruction_manifest",
                    "Results/in-vivo/figure7/intermediates/"
                    "seurat_upstream/reconstruction_manifest.tsv",
                ),
            ]
            with (
                metadata / "analysis_input_manifest.tsv"
            ).open("w", newline="") as handle:
                writer = csv.writer(
                    handle,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writerow(
                    ["role", "repo_relative_path", "sha256", "bytes"]
                )
                for index, (role, locator) in enumerate(rows, start=1):
                    writer.writerow(
                        [role, locator, f"{index:x}" * 64, "1"]
                    )
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
mode=full-refit
figure7_raw_seurat_rds=""
input_paths_for_module si_figures "$1"
"""
            result = subprocess.run(
                ["bash", "-c", script, "fixture", str(run_root)],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
        self.assertEqual(result.returncode, 0, result.stderr)
        paths = set(result.stdout.splitlines())
        for _, locator in rows[:-1]:
            self.assertIn(locator, paths)
        self.assertNotIn(rows[-1][1], paths)


class SiFiguresMaterializationTest(unittest.TestCase):
    def _canonical_fixture(
        self,
        tmp: str,
    ) -> tuple[Path, Path, list[dict[str, object]]]:
        specs = [
            spec for spec in PANEL_SPECS
            if spec["module"] == "si_figures"
        ]
        repo = Path(tmp) / "repo"
        repo.mkdir()
        run_root = (
            repo
            / "Results/in-vivo/SI_figures/runs/source_si_figures"
        )
        (run_root / "figures").mkdir(parents=True)
        (run_root / "metadata").mkdir()
        for spec in specs:
            source = run_root / str(spec["source"])
            source.write_bytes(f"fixture {spec['panel']}\n".encode())
        self._write_publication_contract(repo, run_root, allowed=True)
        return repo, run_root, specs

    @staticmethod
    def _materializer_command(repo: Path, run_root: Path) -> list[str]:
        return [
            sys.executable,
            str(TOOLS_DIR / "materialize_figure_assets.py"),
            "--repo-root",
            str(repo),
            "--figure-root",
            str(repo / "figures"),
            "--source-run-id",
            "source",
            "--operation-id",
            "publish",
            "--module-run",
            f"si_figures={run_root}",
            "--overwrite",
        ]

    def test_contract_contains_only_four_composite_pairs(self) -> None:
        specs = [spec for spec in PANEL_SPECS if spec["module"] == "si_figures"]
        self.assertEqual(len(specs), 8)
        self.assertEqual(
            [spec["panel"] for spec in specs],
            [
                "SuppFig4",
                "SuppFig4_png",
                "SuppFig5",
                "SuppFig5_png",
                "SuppFig6",
                "SuppFig6_png",
                "SuppFig7",
                "SuppFig7_png",
            ],
        )
        self.assertTrue(
            all("composite" in str(spec["source"]) for spec in specs)
        )

    def test_materializer_copies_only_composites(self) -> None:
        specs = [spec for spec in PANEL_SPECS if spec["module"] == "si_figures"]
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            repo.mkdir()
            run_root = (
                repo
                / "Results/in-vivo/SI_figures/runs/source_si_figures"
            )
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()
            for spec in specs:
                source = run_root / str(spec["source"])
                source.write_bytes(f"fixture {spec['panel']}\n".encode())
            self._write_publication_contract(repo, run_root, allowed=True)

            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--source-run-id",
                    "source",
                    "--operation-id",
                    "publish",
                    "--module-run",
                    f"si_figures={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = repo / "figures/Supplementary/manifest.tsv"
            with manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 8)
            self.assertEqual(
                [row["panel"] for row in rows],
                [
                    "SuppFig4",
                    "SuppFig4_png",
                    "SuppFig5",
                    "SuppFig5_png",
                    "SuppFig6",
                    "SuppFig6_png",
                    "SuppFig7",
                    "SuppFig7_png",
                ],
            )
            for row in rows:
                self.assertEqual(row["source_file"], row["asset_path"])
                self.assertEqual(row["input_data"], row["asset_path"])
                self.assertEqual(
                    row["result_run_dir"],
                    "figures/Supplementary",
                )
                self.assertIn("validated_source_sha256=", row["notes"])
            published = {
                path.name
                for path in (repo / "figures/Supplementary").iterdir()
                if path.suffix in {".pdf", ".png"}
            }
            self.assertEqual(
                published,
                {Path(str(spec["asset"])).name for spec in specs},
            )

    def test_materializer_preserves_unselected_supplementary_rows(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, specs = self._canonical_fixture(tmp)
            supplementary = repo / "figures/Supplementary"
            supplementary.mkdir(parents=True)
            legacy_asset = supplementary / "panel_SuppFig1A_legacy.png"
            legacy_asset.write_bytes(b"legacy supplementary panel\n")
            legacy_row = {
                "figure": "Supplementary",
                "panel": "SuppFig1A",
                "asset_path": "manifest:panel_SuppFig1A_legacy.png",
                "source_file": "manifest:panel_SuppFig1A_legacy.png",
                "source_kind": "generated_panel",
                "generated_by": "legacy-test",
                "command": "legacy-test",
                "input_data": "manifest:panel_SuppFig1A_legacy.png",
                "result_run_dir": "manifest:.",
                "run_id": "legacy_run",
                "caption_role": "Unselected legacy panel",
                "asset_status": "generated",
                "not_regenerated_reason": "",
                "local_provenance_path": "",
                "citation_or_uri": "",
                "notes": "must survive targeted SI4-7 materialization",
            }
            write_tsv(
                supplementary / "manifest.tsv",
                [legacy_row],
                FIGURE_MANIFEST_COLUMNS,
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (supplementary / "manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), len(specs) + 1)
            self.assertEqual(
                [row for row in rows if row["panel"] == "SuppFig1A"],
                [legacy_row],
            )

    def test_materializer_requires_composition_helper_input_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            input_manifest = run_root / "metadata/input_manifest.tsv"
            with input_manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows = [
                row for row in rows
                if not row["path"].endswith("normalized_composition.R")
            ]
            write_tsv(
                input_manifest,
                rows,
                MODULE_MANIFEST_COLUMNS,
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("normalized_composition.R", result.stderr)

    def test_materializer_rejects_stale_composition_helper_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            helper = (
                repo
                / "Code/in-vivo/SI_figures/normalized_composition.R"
            )
            helper.write_text(helper.read_text() + "\n# tampered after render\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "does not carry the reviewed frozen-table publication contract",
                result.stderr,
            )

    def test_materializer_rejects_wrong_displayed_panel_set(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            contract = run_root / "metadata/displayed_panel_contract.tsv"
            with contract.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows[0]["panel_ids"] = "D,F"
            write_tsv(
                contract,
                rows,
                ["figure_id", "panel_ids", "selection_policy"],
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "exact manuscript-displayed panel sets",
                result.stderr,
            )

    def test_materializer_requires_shared_context_helper_input_row(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            input_manifest = run_root / "metadata/input_manifest.tsv"
            with input_manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows = [
                row for row in rows
                if not row["path"].endswith("shared_context_panels.R")
            ]
            write_tsv(
                input_manifest,
                rows,
                MODULE_MANIFEST_COLUMNS,
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("shared_context_panels.R", result.stderr)

    def test_materializer_rejects_stale_shared_context_helper_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            helper = (
                repo
                / "Code/in-vivo/SI_figures/shared_context_panels.R"
            )
            helper.write_text(helper.read_text() + "\n# tampered after render\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "does not carry the reviewed frozen-table publication contract",
                result.stderr,
            )

    def test_materializer_rejects_missing_reviewed_cbs_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            matrix = next(
                (repo / "Data/in-vivo/scRNAseq_Numbat").glob("*.sps.cbs")
            )
            matrix.unlink()

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exact reviewed CBS inventory", result.stderr)

    def test_materializer_rejects_tampered_reviewed_cbs_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            matrix = next(
                (repo / "Data/in-vivo/scRNAseq_Numbat").glob("*.sps.cbs")
            )
            with matrix.open("a") as handle:
                handle.write("\n# tampered after review\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("reviewed CBS matrix", result.stderr)

    def test_materializer_rejects_unreviewed_cbs_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            manifest = (
                repo / "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv"
            )
            with manifest.open("a") as handle:
                handle.write("\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("exact reviewed CBS manifest", result.stderr)

    def test_materializer_rejects_tampered_endpoint_ploidy_table(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            endpoint = repo / "Data/in-vivo/all_ploidy.tsv"
            with endpoint.open("a") as handle:
                handle.write("\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("endpoint-ploidy", result.stderr)

    def test_materializer_rejects_tampered_si6_cell_identity(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            annotation_path = (
                run_root
                / "metadata"
                / "figure7J_copy_number_cell_annotations.tsv"
            )
            with annotation_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            rows[0]["sample_id"] = "2N-A1-R"
            write_tsv(
                annotation_path,
                rows,
                FIGURE7J_COPY_NUMBER_ANNOTATION_COLUMNS,
            )
            output_manifest = run_root / "metadata/output_manifest.tsv"
            with output_manifest.open(newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            annotation_locator = str(annotation_path.relative_to(repo))
            matches = [
                row
                for row in output_rows
                if row["repo_relative_path"] == annotation_locator
            ]
            self.assertEqual(len(matches), 1)
            matches[0]["sha256"] = sha256_file(annotation_path)
            matches[0]["byte_size"] = str(annotation_path.stat().st_size)
            write_tsv(
                output_manifest,
                output_rows,
                MODULE_MANIFEST_COLUMNS,
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "exact canonical cell/file/barcode/ploidy identity",
                result.stderr,
            )

    def test_materializer_rejects_wrong_si6_matrix_manifest_role(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            output_manifest = run_root / "metadata/output_manifest.tsv"
            with output_manifest.open(newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            matches = [
                row
                for row in output_rows
                if row["repo_relative_path"].endswith(
                    "figure7J_copy_number_heatmap_matrix.rds"
                )
            ]
            self.assertEqual(len(matches), 1)
            matches[0]["role"] = "output_table"
            matches[0]["source_kind"] = "generated_table"
            write_tsv(
                output_manifest,
                output_rows,
                MODULE_MANIFEST_COLUMNS,
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Figure 7J output-manifest provenance is invalid",
                result.stderr,
            )

    def test_materializer_rejects_tampered_injected_reference_matrix(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            matrix = next(
                (
                    repo
                    / "Data/in-vivo/scRNAseq_Numbat/injected_reference"
                ).glob("*.sps.cbs")
            )
            with matrix.open("a") as handle:
                handle.write("\n# tampered after review\n")

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("injected-cell reference", result.stderr)

    def test_materializer_rejects_tampered_reference_policy_manifest(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo, run_root, _ = self._canonical_fixture(tmp)
            manifest = (
                repo
                / "Data/in-vivo/scRNAseq_Numbat/injected_reference/"
                "reference_manifest.tsv"
            )
            text = manifest.read_text()
            manifest.write_text(
                text.replace(
                    "chr999 haploid-genome-equivalent unassigned DNA",
                    "tampered policy",
                    1,
                )
            )

            result = subprocess.run(
                self._materializer_command(repo, run_root),
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("reference manifest", result.stderr)

    @staticmethod
    def _manifest_row(
        path: Path,
        repo: Path,
        source_id: str,
        *,
        role: str,
        source_kind: str,
        figure: str = "",
        panel: str = "",
    ) -> dict[str, object]:
        relative = str(path.relative_to(repo))
        return {
            "path": relative,
            "repo_relative_path": relative,
            "absolute_path": "",
            "role": role,
            "source_kind": source_kind,
            "module": "si_figures",
            "generated_by": "test",
            "command_id": source_id,
            "sha256": sha256_file(path),
            "checksum_unavailable_reason": "",
            "byte_size": path.stat().st_size,
            "mtime_utc": "2026-07-28T00:00:00+00:00",
            "figure": figure,
            "panel": panel,
            "notes": "strict publication fixture",
        }

    @classmethod
    def _write_publication_contract(
        cls,
        repo: Path,
        run_root: Path,
        allowed: bool,
    ) -> None:
        value = "true" if allowed else "false"
        cache_source = REPO_ROOT / "Data/in-vivo/SIfigures"
        cache_target = repo / "Data/in-vivo/SIfigures"
        shutil.copytree(cache_source, cache_target)
        shutil.copytree(cache_target, run_root / "tables")
        renderer_source = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        )
        renderer_target = (
            repo
            / "Code/in-vivo/SI_figures/generate_supplementary_figures.R"
        )
        renderer_target.parent.mkdir(parents=True)
        shutil.copy2(renderer_source, renderer_target)
        helper_source = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/normalized_composition.R"
        )
        helper_target = (
            repo
            / "Code/in-vivo/SI_figures/normalized_composition.R"
        )
        shutil.copy2(helper_source, helper_target)
        shared_context_source = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/shared_context_panels.R"
        )
        shared_context_target = (
            repo
            / "Code/in-vivo/SI_figures/shared_context_panels.R"
        )
        shutil.copy2(shared_context_source, shared_context_target)
        copy_number_source = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/copy_number_heatmap.R"
        )
        copy_number_target = (
            repo
            / "Code/in-vivo/SI_figures/copy_number_heatmap.R"
        )
        shutil.copy2(copy_number_source, copy_number_target)
        all_ploidy_source = REPO_ROOT / "Data/in-vivo/all_ploidy.tsv"
        all_ploidy_target = repo / "Data/in-vivo/all_ploidy.tsv"
        all_ploidy_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(all_ploidy_source, all_ploidy_target)
        ploidy_derivation_source = (
            REPO_ROOT / "Data/in-vivo/weighted_ploidy.py"
        )
        ploidy_derivation_target = (
            repo / "Data/in-vivo/weighted_ploidy.py"
        )
        shutil.copy2(ploidy_derivation_source, ploidy_derivation_target)
        cbs_source_root = REPO_ROOT / "Data/in-vivo/scRNAseq_Numbat"
        cbs_target_root = repo / "Data/in-vivo/scRNAseq_Numbat"
        cbs_target_root.mkdir(parents=True, exist_ok=True)
        cbs_manifest_target = cbs_target_root / "cbs_manifest.tsv"
        shutil.copy2(
            cbs_source_root / "cbs_manifest.tsv",
            cbs_manifest_target,
        )
        cbs_targets = []
        with cbs_manifest_target.open(newline="") as handle:
            reviewed_cbs_names = [
                row["filename"]
                for row in csv.DictReader(handle, delimiter="\t")
            ]
        for filename in reviewed_cbs_names:
            source = cbs_source_root / filename
            target = cbs_target_root / source.name
            shutil.copy2(source, target)
            cbs_targets.append(target)
        if len(cbs_targets) != 16:
            raise AssertionError("fixture requires exactly 16 CBS matrices")
        reference_source_root = cbs_source_root / "injected_reference"
        reference_target_root = cbs_target_root / "injected_reference"
        reference_target_root.mkdir(parents=True, exist_ok=True)
        reference_manifest_target = (
            reference_target_root / "reference_manifest.tsv"
        )
        shutil.copy2(
            reference_source_root / "reference_manifest.tsv",
            reference_manifest_target,
        )
        with reference_manifest_target.open(newline="") as handle:
            reviewed_reference_names = [
                row["filename"]
                for row in csv.DictReader(handle, delimiter="\t")
            ]
        reference_targets = []
        for filename in reviewed_reference_names:
            source = reference_source_root / filename
            target = reference_target_root / source.name
            shutil.copy2(source, target)
            reference_targets.append(target)
        if len(reference_targets) != 2:
            raise AssertionError(
                "fixture requires exactly two injected-cell references"
            )
        runner_source = (
            REPO_ROOT
            / "Code/in-vivo/SI_figures/run_supplementary_figures.R"
        )
        runner_target = (
            repo
            / "Code/in-vivo/SI_figures/run_supplementary_figures.R"
        )
        shutil.copy2(runner_source, runner_target)
        validator_source = (
            REPO_ROOT / "Code/tools/validate_si_figures_table_cache.py"
        )
        validator_target = (
            repo / "Code/tools/validate_si_figures_table_cache.py"
        )
        validator_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(validator_source, validator_target)
        config_source = (
            REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml"
        )
        config_target = (
            repo / "Code/in-vivo/figure7/figure7_config.yaml"
        )
        config_target.parent.mkdir(parents=True)
        shutil.copy2(config_source, config_target)
        figure7_shared_targets = []
        for relative in (
            "Code/in-vivo/figure7/src/common_io.R",
            "Code/in-vivo/figure7/src/tgi_statistics.R",
            "Code/in-vivo/figure7/src/tgi_panels.R",
            "Code/in-vivo/figure7/density_localization_config.yaml",
        ):
            source = REPO_ROOT / relative
            target = repo / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
            figure7_shared_targets.append(target)
        cellcycle_relative = (
            "Data/in-vivo/figure7/processed/"
            "CellCycleCells_pseudotime_distribution_per_sample_"
            "cell_level_with_ploidy_dose_tgi.csv"
        )
        cellcycle_source = REPO_ROOT / cellcycle_relative
        cellcycle_target = repo / cellcycle_relative
        cellcycle_target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(cellcycle_source, cellcycle_target)

        cache_manifest = cache_target / "manifest.tsv"
        provenance_rows = [
            {"key": "artifact", "value": "supplementary_figures_4_7"},
            {
                "key": "entrypoint",
                "value": (
                    "Code/in-vivo/SI_figures/"
                    "generate_supplementary_figures.R"
                ),
            },
            {
                "key": "entrypoint_sha256",
                "value": sha256_file(renderer_target),
            },
            {
                "key": "normalized_composition_helper",
                "value": (
                    "Code/in-vivo/SI_figures/"
                    "normalized_composition.R"
                ),
            },
            {
                "key": "normalized_composition_helper_sha256",
                "value": sha256_file(helper_target),
            },
            {
                "key": "shared_context_panels_helper",
                "value": (
                    "Code/in-vivo/SI_figures/shared_context_panels.R"
                ),
            },
            {
                "key": "shared_context_panels_helper_sha256",
                "value": sha256_file(shared_context_target),
            },
            {
                "key": "si4i_cellcycle_pseudotime",
                "value": cellcycle_relative,
            },
            {
                "key": "si4i_cellcycle_pseudotime_sha256",
                "value": sha256_file(cellcycle_target),
            },
            {
                "key": "si4i_density_localization_config",
                "value": (
                    "Code/in-vivo/figure7/density_localization_config.yaml"
                ),
            },
            {
                "key": "si4i_density_localization_config_sha256",
                "value": sha256_file(figure7_shared_targets[3]),
            },
            {
                "key": "si4i_density_localization_method",
                "value": (
                    "equal-mouse Gaussian-kernel treated-minus-vehicle "
                    "density contrast; exact injected-origin-stratified "
                    "pointwise and studentized max-|T| permutation "
                    "inference on the reviewed 2,881-cell subset"
                ),
            },
            {
                "key": "copy_number_heatmap_helper",
                "value": (
                    "Code/in-vivo/SI_figures/copy_number_heatmap.R"
                ),
            },
            {
                "key": "copy_number_heatmap_helper_sha256",
                "value": sha256_file(copy_number_target),
            },
            {
                "key": "endpoint_ploidy_derivation_helper",
                "value": "Data/in-vivo/weighted_ploidy.py",
            },
            {
                "key": "endpoint_ploidy_derivation_helper_sha256",
                "value": sha256_file(ploidy_derivation_target),
            },
            {
                "key": "endpoint_ploidy_source",
                "value": "Data/in-vivo/all_ploidy.tsv",
            },
            {
                "key": "endpoint_ploidy_source_sha256",
                "value": sha256_file(all_ploidy_target),
            },
            {
                "key": "numbat_cbs_manifest",
                "value": "Data/in-vivo/scRNAseq_Numbat/cbs_manifest.tsv",
            },
            {
                "key": "numbat_cbs_manifest_sha256",
                "value": sha256_file(cbs_manifest_target),
            },
            {"key": "numbat_cbs_matrix_count", "value": "16"},
            {
                "key": "numbat_cbs_matrix_hashes",
                "value": ";".join(
                    f"{path.name}={sha256_file(path)}"
                    for path in cbs_targets
                ),
            },
            {
                "key": "figure7j_qc_selection",
                "value": (
                    "the complete 14,125-cell CBS source is checksum/value "
                    "validated, then restricted by the frozen endpoint-ploidy "
                    "audit to the exact 9,832 QC-passed tumor cells (5,335 "
                    "treated); clusters 3, 4, 9, and 9c remain excluded"
                ),
            },
            {
                "key": "injected_cell_reference_manifest",
                "value": (
                    "Data/in-vivo/scRNAseq_Numbat/injected_reference/"
                    "reference_manifest.tsv"
                ),
            },
            {
                "key": "injected_cell_reference_manifest_sha256",
                "value": sha256_file(reference_manifest_target),
            },
            {
                "key": "injected_cell_reference_matrix_hashes",
                "value": ";".join(
                    f"{path.name}={sha256_file(path)}"
                    for path in reference_targets
                ),
            },
            {
                "key": "figure7j_harmonization",
                "value": (
                    "each separately postprocessed file-specific CBS schema "
                    "reduced independently to chr1-22 finite-segment "
                    "length-weighted means; exported and finite coverage "
                    "audited separately; no coordinate alignment; rows and "
                    "columns unclustered"
                ),
            },
            {
                "key": "si6e_reference_source",
                "value": (
                    "miningcloneid@"
                    "c505cd9159fa2a8c0974c7379f6aacd09fe19abc"
                ),
            },
            {
                "key": "si6e_reference_chr999_interpretation",
                "value": (
                    "haploid-genome-equivalent unassigned DNA; "
                    "project-confirmed 2026-08-01"
                ),
            },
            {
                "key": "si6e_reference_designation_basis",
                "value": (
                    "2N=project-designated 2N-lineage proxy; not the "
                    "same-passage A6M inoculum;4N=project-designated "
                    "4N-lineage proxy; not the same-passage A4M inoculum"
                ),
            },
            {
                "key": "si6e_reference_ploidy_policy",
                "value": (
                    "autosomal length-weighted estimate plus chr999 "
                    "haploid-genome-equivalent unassigned DNA"
                ),
            },
            {
                "key": "si6e_summary_analysis_type",
                "value": (
                    "descriptive_only; no endpoint cross-origin or "
                    "reference-to-endpoint test"
                ),
            },
            {
                "key": "si6e_ploidy_reduction_comparison",
                "value": (
                    "project-designated lineage-matched 2N-A7M/4N-A5M "
                    "karyotype reference distributions, including the "
                    "chr999 unassigned-extra-DNA haploid-genome-equivalent "
                    "term added to autosomal ploidy, compared descriptively "
                    "with one postprocessed endpoint mean per mouse; no "
                    "formal P value because each reference is one culture-level "
                    "biological unit and origin-specific endpoint runs use "
                    "different schemas/calibration"
                ),
            },
            {
                "key": "table_cache_manifest_sha256",
                "value": sha256_file(cache_manifest),
            },
            {
                "key": "si7_feature_species_policy",
                "value": SI7_REVIEWED_FEATURE_POLICY,
            },
            {
                "key": "si7_gene_set_database",
                "value": SI7_REVIEWED_GENE_SET_DATABASE,
            },
            {
                "key": "si7_canonical_publication_allowed",
                "value": value,
            },
            {
                "key": "si7_frozen_matrix_note",
                "value": SI7_REVIEWED_FROZEN_MATRIX_NOTE,
            },
        ]
        run_config_rows = [
            {"key": "module", "value": "si_figures"},
            {"key": "figures", "value": "4,5,6,7"},
            {
                "key": "displayed_panel_sets",
                "value": (
                    "SuppFig4=D,F,I;SuppFig5=E,G;"
                    "SuppFig6=A,B,C,D,E;SuppFig7=A"
                ),
            },
            {"key": "figure_file_count", "value": "8"},
            {"key": "table_mode", "value": "frozen_plot_tables_only"},
            {"key": "cache_file_count", "value": "11"},
            {
                "key": "composition_normalization",
                "value": (
                    "within-sample cluster proportions averaged with equal "
                    "sample weights within group"
                ),
            },
            {
                "key": "composition_enrichment_test",
                "value": (
                    "exact independent-sample label permutation; one group "
                    "versus exchangeable remaining samples"
                ),
            },
            {
                "key": "composition_permutation_strata",
                "value": (
                    "SI4E=initial_ploidy;SI4G=context;SI5F=descriptive;"
                    "SI5G=initial_ploidy;SI5H=initial_ploidy;SI5I=dose"
                ),
            },
            {
                "key": "composition_multiple_testing",
                "value": (
                    "Benjamini-Hochberg across all group-by-cluster "
                    "contrasts within each panel"
                ),
            },
            {"key": "composition_fdr_threshold", "value": "0.05"},
            {
                "key": "si4i_cell_universe",
                "value": "reviewed_qc_retained_cellcycle_2881",
            },
            {"key": "si4i_n_cells", "value": "2881"},
            {"key": "si4i_n_mice", "value": "16"},
            {
                "key": "si4i_pointwise_positive_interval",
                "value": "0.296-0.486",
            },
            {
                "key": "si4i_simultaneous_positive_interval",
                "value": "0.414-0.426",
            },
            {
                "key": "si4i_global_max_abs_t_p_two_sided",
                "value": "0.0473469387755102",
            },
            {
                "key": "figure7j_copy_number_source",
                "value": (
                    "postprocessed NUMBAT-derived cell-by-segment CBS matrices"
                ),
            },
            {
                "key": "figure7j_column_statistic",
                "value": (
                    "per-cell length-weighted mean across available CBS "
                    "segments within each chromosome and file-specific schema"
                ),
            },
            {"key": "figure7j_cbs_matrix_count", "value": "16"},
            {"key": "figure7j_source_cell_count", "value": "14125"},
            {
                "key": "figure7j_qc_passed_tumor_cell_count",
                "value": "9832",
            },
            {
                "key": "figure7j_treated_tumor_cell_count",
                "value": "5335",
            },
            {
                "key": "figure7j_qc_selection_policy",
                "value": (
                    "exact cells matched by the frozen endpoint-ploidy audit "
                    "from the final QC-curated Seurat object; clusters 3, 4, "
                    "9, and 9c excluded"
                ),
            },
            {"key": "figure7j_chromosome_count", "value": "22"},
            {
                "key": "figure7j_row_order",
                "value": (
                    "injected origin, dose, mouse, post-processed copy-number "
                    "score, cell ID; no row clustering"
                ),
            },
            {
                "key": "figure7j_column_order",
                "value": (
                    "chromosomes 1-22 in genomic order; no column clustering"
                ),
            },
            {
                "key": "si6_postprocessed_copy_number_score_unit",
                "value": "sequenced mouse/CBS file",
            },
            {
                "key": "si6_injected_reference_ploidy_policy",
                "value": (
                    "autosomal length-weighted estimate plus chr999 "
                    "haploid-genome-equivalent unassigned DNA"
                ),
            },
            {
                "key": "si6_injected_reference_cell_counts",
                "value": "2N=20;4N=16",
            },
            {
                "key": "si6_injected_reference_source_repository",
                "value": "miningcloneid",
            },
            {
                "key": "si6_injected_reference_source_commit",
                "value": "c505cd9159fa2a8c0974c7379f6aacd09fe19abc",
            },
            {
                "key": "si6_injected_reference_chr999_unit",
                "value": "haploid-genome-equivalent unassigned DNA",
            },
            {
                "key": "si6_injected_reference_chr999_interpretation_basis",
                "value": "project-confirmed 2026-08-01",
            },
            {
                "key": "si6_endpoint_summary_analysis_type",
                "value": "descriptive_only",
            },
            {
                "key": "si6_2n_reference_mean_ploidy",
                "value": "2.151242953153243",
            },
            {
                "key": "si6_2n_endpoint_mouse_balanced_mean_ploidy",
                "value": "2.135513297009228",
            },
            {
                "key": "si6_2n_relative_change_percent",
                "value": "-0.7311892002230036",
            },
            {
                "key": "si6_4n_reference_mean_ploidy",
                "value": "3.94651299970957",
            },
            {
                "key": "si6_4n_endpoint_mouse_balanced_mean_ploidy",
                "value": "2.318552949440266",
            },
            {
                "key": "si6_4n_relative_change_percent",
                "value": "-41.250593888557",
            },
            {
                "key": "si6_reference_4n_minus_2n_mean_ploidy",
                "value": "1.795270046556327",
            },
            {
                "key": "si6_endpoint_4n_minus_2n_mouse_balanced_mean_ploidy",
                "value": "0.1830396524310385",
            },
            {
                "key": "si6_separation_contraction_percent",
                "value": "89.80433875214797",
            },
            {
                "key": "si7_feature_species_policy",
                "value": SI7_REVIEWED_FEATURE_POLICY,
            },
            {
                "key": "si7_gene_set_database",
                "value": SI7_REVIEWED_GENE_SET_DATABASE,
            },
            {
                "key": "si7_canonical_publication_allowed",
                "value": value,
            },
        ]
        write_tsv(
            run_root / "metadata/si_figures_provenance.tsv",
            provenance_rows,
            ["key", "value"],
        )
        write_tsv(
            run_root / "metadata/run_config.tsv",
            run_config_rows,
            ["key", "value"],
        )
        write_tsv(
            run_root / "metadata/displayed_panel_contract.tsv",
            [
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
            ],
            ["figure_id", "panel_ids", "selection_policy"],
        )

        specs = [
            spec for spec in PANEL_SPECS
            if spec["module"] == "si_figures"
        ]
        write_tsv(
            run_root / "metadata/panel_contract.tsv",
            [
                {
                    "panel_id": spec["panel"],
                    "filename": Path(str(spec["source"])).name,
                }
                for spec in specs
            ],
            ["panel_id", "filename"],
        )

        annotation_path = (
            run_root
            / "metadata"
            / "figure7J_copy_number_cell_annotations.tsv"
        )
        matrix_path = (
            run_root
            / "metadata"
            / "figure7J_copy_number_heatmap_matrix.rds"
        )
        audit_path = (
            cache_target / "si_figure6_endpoint_ploidy_join_audit.csv"
        )
        annotation_rows = []
        with audit_path.open(newline="") as handle:
            audit_rows = csv.DictReader(handle)
            for row in audit_rows:
                if (
                    row["context"] != "Tumor"
                    or row["matched"].strip().lower() != "true"
                ):
                    continue
                cell_id = row["endpoint_cell_id"]
                suffix = f"_{cell_id}"
                if not row["cell"].endswith(suffix):
                    raise AssertionError(
                        "fixture audit cell does not have its CBS barcode suffix"
                    )
                sample_id = row["cell"][: -len(suffix)]
                file_name = row["endpoint_file"]
                annotation_rows.append(
                    {
                        "heatmap_row_id": f"{file_name}::{cell_id}",
                        "file": file_name,
                        "cell_id": cell_id,
                        "sample_id": sample_id,
                        "initial_ploidy": row["initial_ploidy"],
                        "dose_mg_per_kg": file_name.split("-")[2],
                        "endpoint_ploidy": row["endpoint_ploidy"],
                        "frac_covered": "1",
                        "display_order": str(len(annotation_rows) + 1),
                    }
                )
        if len(annotation_rows) != 9832:
            raise AssertionError(
                "fixture audit must contain exactly 9,832 matched tumor cells"
            )
        write_tsv(
            annotation_path,
            annotation_rows,
            FIGURE7J_COPY_NUMBER_ANNOTATION_COLUMNS,
        )
        matrix_path.write_bytes(b"strict SI6 matrix RDS fixture\n")

        si4i_grid_path = (
            run_root
            / "metadata/si_figure4I_density_localization_grid.tsv"
        )
        si4i_interval_path = (
            run_root
            / "metadata/si_figure4I_density_localization_intervals.tsv"
        )
        si4i_test_path = (
            run_root
            / "metadata/si_figure4I_density_localization_test.tsv"
        )
        write_tsv(
            si4i_grid_path,
            [
                {
                    "pseudotime": f"{index / 500:.3f}",
                    "treated_minus_vehicle_density": "0",
                }
                for index in range(501)
            ],
            ["pseudotime", "treated_minus_vehicle_density"],
        )
        write_tsv(
            si4i_interval_path,
            [
                {
                    "support_type": "positive_pointwise_two_sided",
                    "start": "0.296",
                    "end": "0.486",
                },
                {
                    "support_type": "positive_simultaneous_max_abs_t",
                    "start": "0.414",
                    "end": "0.426",
                },
            ],
            ["support_type", "start", "end"],
        )
        write_tsv(
            si4i_test_path,
            [
                {
                    "analysis_id": (
                        "equal_mouse_kde_exact_origin_stratified_max_t_v1"
                    ),
                    "cell_universe": "reviewed_qc_retained_cellcycle_2881",
                    "n_cells": "2881",
                    "n_samples": "16",
                    "pointwise_start": "0.296",
                    "pointwise_end": "0.486",
                    "simultaneous_start": "0.414",
                    "simultaneous_end": "0.426",
                    "raw_excess_peak_pseudotime": "0.452",
                    "max_abs_t_pseudotime": "0.420",
                    "global_max_t_p_two_sided": "0.0473469387755102",
                }
            ],
            [
                "analysis_id",
                "cell_universe",
                "n_cells",
                "n_samples",
                "pointwise_start",
                "pointwise_end",
                "simultaneous_start",
                "simultaneous_end",
                "raw_excess_peak_pseudotime",
                "max_abs_t_pseudotime",
                "global_max_t_p_two_sided",
            ],
        )

        source_id = run_root.name.removesuffix("_si_figures")
        with cache_manifest.open(newline="") as handle:
            cache_rows = list(csv.DictReader(handle, delimiter="\t"))
        input_paths = [
            config_target,
            runner_target,
            renderer_target,
            helper_target,
            shared_context_target,
            copy_number_target,
            *figure7_shared_targets,
            cellcycle_target,
            ploidy_derivation_target,
            all_ploidy_target,
            cbs_manifest_target,
            *cbs_targets,
            reference_manifest_target,
            *reference_targets,
            validator_target,
            cache_manifest,
            *[
                cache_target / row["filename"]
                for row in cache_rows
            ],
        ]
        write_tsv(
            run_root / "metadata/input_manifest.tsv",
            [
                cls._manifest_row(
                    path,
                    repo,
                    source_id,
                    role="input_data",
                    source_kind="input_file",
                )
                for path in input_paths
            ],
            MODULE_MANIFEST_COLUMNS,
        )
        write_tsv(
            run_root / "metadata/output_manifest.tsv",
            [
                cls._manifest_row(
                    run_root / str(spec["source"]),
                    repo,
                    source_id,
                    role="output_figure",
                    source_kind="generated_panel",
                    figure="Supplementary",
                    panel=str(spec["panel"]),
                )
                for spec in specs
            ]
            + [
                cls._manifest_row(
                    annotation_path,
                    repo,
                    source_id,
                    role="output_table",
                    source_kind="generated_table",
                ),
                cls._manifest_row(
                    matrix_path,
                    repo,
                    source_id,
                    role="output_model",
                    source_kind="generated_model",
                ),
                *[
                    cls._manifest_row(
                        path,
                        repo,
                        source_id,
                        role="output_table",
                        source_kind="generated_table",
                    )
                    for path in (
                        si4i_grid_path,
                        si4i_interval_path,
                        si4i_test_path,
                    )
                ],
            ],
            MODULE_MANIFEST_COLUMNS,
        )

    def test_materializer_rejects_noncanonical_raw_si_run(self) -> None:
        specs = [
            spec for spec in PANEL_SPECS
            if spec["module"] == "si_figures"
        ]
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            run_root = (
                repo
                / "Results/in-vivo/SI_figures/runs/source_si_figures"
            )
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()
            for spec in specs:
                source = run_root / str(spec["source"])
                source.write_bytes(b"noncanonical raw fixture\n")
            self._write_publication_contract(
                repo,
                run_root,
                allowed=False,
            )
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--source-run-id",
                    "source",
                    "--module-run",
                    f"si_figures={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "Canonical SI Figures materialization is prohibited",
                result.stderr,
            )
            self.assertFalse((repo / "figures").exists())
            (
                run_root
                / "metadata"
                / "si_figures_provenance.tsv"
            ).unlink()
            missing = subprocess.run(
                result.args,
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(missing.returncode, 0)
            self.assertIn(
                "Missing source SI Figures provenance",
                missing.stderr,
            )
            self.assertFalse((repo / "figures").exists())

    def test_materializer_rejects_tampered_run_scoped_table_cache(
        self,
    ) -> None:
        specs = [
            spec for spec in PANEL_SPECS
            if spec["module"] == "si_figures"
        ]
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            run_root = (
                repo
                / "Results/in-vivo/SI_figures/runs/source_si_figures"
            )
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()
            for spec in specs:
                source = run_root / str(spec["source"])
                source.write_bytes(b"canonical render fixture\n")
            self._write_publication_contract(repo, run_root, allowed=True)
            (
                run_root
                / "tables"
                / "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv"
            ).write_text("tampered\n")

            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--source-run-id",
                    "source",
                    "--module-run",
                    f"si_figures={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "is not the exact reviewed cache table",
                result.stderr,
            )
            self.assertFalse((repo / "figures").exists())

    def test_materializer_rejects_false_reviewed_si7_lineage(self) -> None:
        specs = [
            spec for spec in PANEL_SPECS
            if spec["module"] == "si_figures"
        ]
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            run_root = (
                repo
                / "Results/in-vivo/SI_figures/runs/source_si_figures"
            )
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()
            for spec in specs:
                source = run_root / str(spec["source"])
                source.write_bytes(b"canonical render fixture\n")
            self._write_publication_contract(repo, run_root, allowed=True)
            provenance_path = (
                run_root / "metadata/si_figures_provenance.tsv"
            )
            with provenance_path.open(newline="") as handle:
                provenance_rows = list(
                    csv.DictReader(handle, delimiter="\t")
                )
            for row in provenance_rows:
                if row["key"] == "si7_frozen_matrix_note":
                    row["value"] = (
                        "Frozen matrices came from the historical DEG cache."
                    )
            write_tsv(
                provenance_path,
                provenance_rows,
                ["key", "value"],
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--source-run-id",
                    "source",
                    "--module-run",
                    f"si_figures={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn(
                "does not carry the reviewed frozen-table publication contract",
                result.stderr,
            )
            self.assertFalse((repo / "figures").exists())


if __name__ == "__main__":
    unittest.main()
