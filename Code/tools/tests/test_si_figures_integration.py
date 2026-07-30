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
    MODULE_MANIFEST_COLUMNS,
    sha256_file,
    write_tsv,
)
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    SI7_REVIEWED_FEATURE_POLICY,
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
        self.assertIn(
            f'si7_feature_species_policy: "{SI7_REVIEWED_FEATURE_POLICY}"',
            config_text,
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
        self.assertLess(
            result.stdout.index("[in_vivo_figure7]"),
            result.stdout.index("[si_figures]"),
        )

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
            f"Missing required file: {missing_ploidy}",
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
                "Code/tools/validate_si_figures_table_cache.py",
                paths,
            )
            self.assertIn(str(retained), paths)
            self.assertNotIn(
                "Results/in-vivo/figure7/intermediates/"
                "seurat_upstream/reconstruction_manifest.tsv",
                paths,
            )
            self.assertEqual(
                set(paths),
                {
                    "Code/in-vivo/SI_figures/"
                    "run_supplementary_figures.R",
                    "Code/in-vivo/SI_figures/"
                    "generate_supplementary_figures.R",
                    "Code/tools/validate_si_figures_table_cache.py",
                    str(retained),
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
                    "si_figures_cache_manifest",
                    "Data/in-vivo/SIfigures/manifest.tsv",
                ),
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
            published = {
                path.name
                for path in (repo / "figures/Supplementary").iterdir()
                if path.suffix in {".pdf", ".png"}
            }
            self.assertEqual(
                published,
                {Path(str(spec["asset"])).name for spec in specs},
            )

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
        config_source = (
            REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml"
        )
        config_target = (
            repo / "Code/in-vivo/figure7/figure7_config.yaml"
        )
        config_target.parent.mkdir(parents=True)
        shutil.copy2(config_source, config_target)

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
        ]
        run_config_rows = [
            {"key": "module", "value": "si_figures"},
            {"key": "figures", "value": "4,5,6,7"},
            {"key": "figure_file_count", "value": "8"},
            {"key": "table_mode", "value": "frozen_plot_tables_only"},
            {"key": "cache_file_count", "value": "11"},
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

        source_id = run_root.name.removesuffix("_si_figures")
        with cache_manifest.open(newline="") as handle:
            cache_rows = list(csv.DictReader(handle, delimiter="\t"))
        input_paths = [
            config_target,
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


if __name__ == "__main__":
    unittest.main()
