from __future__ import annotations

import csv
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = TOOLS_DIR.parents[1]
sys.path.insert(0, str(TOOLS_DIR))

from figure_output_contract import (  # noqa: E402
    MODULE_MANIFEST_COLUMNS,
    sha256_file,
    validate_expected_panel_set,
    validate_figure_manifest,
    validate_module_manifest,
    write_tsv,
)
from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class Figure7MaterializationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name) / "repo"
        self.repo.mkdir()
        self.source_id = "analysis_run"
        self.operation_id = "materialize_run"
        self.run_root = (
            self.repo
            / "Results/in-vivo/figure7/runs"
            / f"{self.source_id}_figure7"
        )
        (self.run_root / "figures").mkdir(parents=True)
        (self.run_root / "metadata").mkdir()
        self.input_path = self.repo / "Data/in-vivo/figure7/fixture.tsv"
        self.input_path.parent.mkdir(parents=True)
        self.input_path.write_text("value\n1\n")
        self.figure7_specs = [
            spec for spec in PANEL_SPECS if spec["module"] == "in_vivo_figure7"
        ]
        for spec in self.figure7_specs:
            source = self.run_root / str(spec["source"])
            source.write_bytes(f"fake PDF for {spec['panel']}\n".encode())
        self._write_run_metadata(include_f=True)
        self._write_input_manifest()
        self._write_output_manifest()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _write_run_metadata(self, include_f: bool) -> None:
        selected = [
            spec for spec in self.figure7_specs
            if spec.get("variant", "pdf") == "pdf"
            and (include_f or not str(spec["panel"]).startswith("7F"))
        ]
        write_tsv(
            self.run_root / "metadata/run_config.tsv",
            [{"key": "panel_set", "value": "a-f" if include_f else "a-e"}],
            ["key", "value"],
        )
        write_tsv(
            self.run_root / "metadata/panel_contract.tsv",
            [
                {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                for spec in selected
            ],
            ["panel_id", "filename"],
        )

    def _write_input_manifest(self) -> None:
        relative = str(self.input_path.relative_to(self.repo))
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [
                {
                    "path": relative,
                    "repo_relative_path": relative,
                    "absolute_path": str(self.input_path),
                    "role": "input_data",
                    "source_kind": "input_file",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(self.input_path),
                    "checksum_unavailable_reason": "",
                    "byte_size": self.input_path.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "",
                    "panel": "",
                    "notes": "test fixture",
                }
            ],
            MODULE_MANIFEST_COLUMNS,
        )

    def _write_output_manifest(self) -> None:
        rows = []
        for spec in self.figure7_specs:
            source = self.run_root / str(spec["source"])
            relative = str(source.relative_to(self.repo))
            rows.append(
                {
                    "path": relative,
                    "repo_relative_path": relative,
                    "absolute_path": str(source),
                    "role": "output_figure",
                    "source_kind": "generated_panel",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(source),
                    "checksum_unavailable_reason": "",
                    "byte_size": source.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "Figure7",
                    "panel": spec["panel"],
                    "notes": "test",
                }
            )
        write_tsv(
            self.run_root / "metadata/output_manifest.tsv",
            rows,
            MODULE_MANIFEST_COLUMNS,
        )

    def _run_materializer(self) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            [
                sys.executable,
                str(TOOLS_DIR / "materialize_figure_assets.py"),
                "--repo-root",
                str(self.repo),
                "--figure-root",
                str(self.repo / "figures"),
                "--source-run-id",
                self.source_id,
                "--operation-id",
                self.operation_id,
                "--module-run",
                f"in_vivo_figure7={self.run_root}",
                "--overwrite",
            ],
            text=True,
            capture_output=True,
        )

    def test_materializes_exact_six_with_source_provenance(self) -> None:
        result = self._run_materializer()
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = self.repo / "figures/Figure7/manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected_panels = [value for letter in "ABCDEF" for value in (f"7{letter}", f"7{letter}_png")]
        self.assertEqual([row["panel"] for row in rows], expected_panels)
        self.assertEqual({row["run_id"] for row in rows}, {self.source_id})
        self.assertEqual(
            {row["result_run_dir"] for row in rows},
            {str(self.run_root.relative_to(self.repo))},
        )
        self.assertTrue(
            all(f"materialization_operation_id={self.operation_id}" in row["notes"] for row in rows)
        )

    def test_materializes_explicit_ae_run_without_optional_panel_f(self) -> None:
        panel_f = [spec for spec in self.figure7_specs if str(spec["panel"]).startswith("7F")]
        for spec in panel_f:
            (self.run_root / str(spec["source"])).unlink()
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = [
                row for row in csv.DictReader(handle, delimiter="\t")
                if not row["panel"].startswith("7F")
            ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        self._write_run_metadata(include_f=False)

        result = self._run_materializer()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.repo / "figures/Figure7/manifest.tsv").open(newline="") as handle:
            materialized = list(csv.DictReader(handle, delimiter="\t"))
        expected_panels = [value for letter in "ABCDE" for value in (f"7{letter}", f"7{letter}_png")]
        self.assertEqual([row["panel"] for row in materialized], expected_panels)
        self.assertFalse((self.repo / "figures/Figure7/panel_7F_pseudotime_state_pathway_activity.pdf").exists())

    def test_rejects_unrecorded_ae_omission(self) -> None:
        panel_f = [spec for spec in self.figure7_specs if str(spec["panel"]).startswith("7F")]
        for spec in panel_f:
            (self.run_root / str(spec["source"])).unlink()
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = [
                row for row in csv.DictReader(handle, delimiter="\t")
                if not row["panel"].startswith("7F")
            ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("explicitly record panel_set=a-e", result.stderr)

    def test_rejects_missing_or_unexpected_figure(self) -> None:
        (self.run_root / str(self.figure7_specs[0]["source"])).unlink()
        missing = self._run_materializer()
        self.assertNotEqual(missing.returncode, 0)

        (self.run_root / str(self.figure7_specs[0]["source"])).write_bytes(b"restored")
        self._write_output_manifest()
        (self.run_root / "tables").mkdir()
        (self.run_root / "tables/unexpected.png").write_bytes(b"unexpected")
        unexpected = self._run_materializer()
        self.assertNotEqual(unexpected.returncode, 0)
        self.assertIn("exact panel inventory", unexpected.stderr)

    def test_rejects_checksum_mismatch(self) -> None:
        source = self.run_root / str(self.figure7_specs[0]["source"])
        source.write_bytes(b"changed after manifest")
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sha256 does not match", result.stderr)

    def test_rejects_duplicate_source_manifest_row(self) -> None:
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows.append(dict(rows[0]))
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Expected one output-manifest row", result.stderr)

    def test_rejects_missing_or_invalid_input_manifest(self) -> None:
        input_manifest = self.run_root / "metadata/input_manifest.tsv"
        input_manifest.unlink()
        missing = self._run_materializer()
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("Missing source input manifest", missing.stderr)

        self._write_input_manifest()
        self.input_path.write_text("changed after manifest\n")
        invalid = self._run_materializer()
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("Invalid source input manifest", invalid.stderr)

    def test_rejects_input_manifest_from_another_run(self) -> None:
        input_manifest = self.run_root / "metadata/input_manifest.tsv"
        with input_manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[0]["command_id"] = "different_source_run"
        write_tsv(input_manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("input-manifest provenance", result.stderr)


class ManifestContractTest(unittest.TestCase):
    def test_expected_set_and_duplicate_checks(self) -> None:
        rows = [
            {"panel": "7A", "source_kind": "generated_panel"},
            {"panel": "7A", "source_kind": "external"},
        ]
        errors = validate_expected_panel_set(rows, {"7A", "7B"}, Path("manifest.tsv"))
        self.assertTrue(any("duplicate panel ID" in error for error in errors))
        self.assertTrue(any("missing expected generated panel" in error for error in errors))

    def test_repo_relative_path_precedes_stale_absolute_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            data = repo / "Results/run/value.tsv"
            data.parent.mkdir(parents=True)
            data.write_text("value\n1\n")
            manifest = repo / "manifest.tsv"
            write_tsv(
                manifest,
                [
                    {
                        "path": "Results/run/value.tsv",
                        "repo_relative_path": "Results/run/value.tsv",
                        "absolute_path": "/old/checkout/Results/run/value.tsv",
                        "role": "output_table",
                        "source_kind": "generated_table",
                        "module": "test",
                        "generated_by": "test",
                        "command_id": "run",
                        "sha256": sha256_file(data),
                        "checksum_unavailable_reason": "",
                        "byte_size": data.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "test",
                    }
                ],
                MODULE_MANIFEST_COLUMNS,
            )
            self.assertEqual(validate_module_manifest(manifest, repo, repo / "Results/run"), [])

    def test_external_absolute_path_follows_missing_repo_relative_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "checkout"
            repo.mkdir()
            external = Path(tmp) / "external.tsv"
            external.write_text("value\n1\n")
            manifest = repo / "manifest.tsv"
            write_tsv(
                manifest,
                [
                    {
                        "path": "Data/missing.tsv",
                        "repo_relative_path": "Data/missing.tsv",
                        "absolute_path": str(external),
                        "role": "input_data",
                        "source_kind": "input_file",
                        "module": "test",
                        "generated_by": "test",
                        "command_id": "run",
                        "sha256": sha256_file(external),
                        "checksum_unavailable_reason": "",
                        "byte_size": external.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "external fallback",
                    }
                ],
                MODULE_MANIFEST_COLUMNS,
            )
            self.assertEqual(validate_module_manifest(manifest, repo), [])

    def test_existing_figure_materialization_remains_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            run_root = repo / "Results/gdsc/runs/source_gdsc"
            (run_root / "drug_count_collapsed/figures").mkdir(parents=True)
            (run_root / "ploidy_enrichment_panels_primary_secondary_ABC.png").write_bytes(b"one")
            (
                run_root
                / "drug_count_collapsed/figures/drug_count_collapsed_enrichment_shared_order.png"
            ).write_bytes(b"two")
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--run-id",
                    "source",
                    "--module-run",
                    f"gdsc={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((repo / "figures/Figure1/manifest.tsv").is_file())

    def test_optional_lci_panel_remains_optional(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            drug_run = repo / "Results/drug_response"
            lci_run = repo / "Results/lci_overlays"
            lci_run.mkdir(parents=True)
            for spec in PANEL_SPECS:
                if spec["module"] != "drug_response":
                    continue
                source = drug_run / str(spec["source"])
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(b"source")
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root", str(repo),
                    "--figure-root", str(repo / "figures"),
                    "--run-id", "source",
                    "--module-run", f"drug_response={drug_run}",
                    "--module-run", f"lci_overlays={lci_run}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (repo / "figures/Figure3/manifest.tsv").open(newline="") as handle:
                panels = {row["panel"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertNotIn("3I", panels)

    def test_existing_manifests_validate_in_a_moved_checkout(self) -> None:
        manifest_paths = sorted((REPO_ROOT / "figures").glob("Figure[1-6]/manifest.tsv"))
        manifest_paths.append(REPO_ROOT / "figures/Supplementary/manifest.tsv")
        with tempfile.TemporaryDirectory() as tmp:
            moved_repo = Path(tmp) / "moved-checkout"
            for manifest in manifest_paths:
                with manifest.open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                for row in rows:
                    for key in ("asset_path", "source_file", "result_run_dir"):
                        value = row.get(key, "").strip()
                        if not value or Path(value).is_absolute():
                            continue
                        source = REPO_ROOT / value
                        destination = moved_repo / value
                        if source.is_file():
                            destination.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source, destination)
                        elif source.is_dir():
                            destination.mkdir(parents=True, exist_ok=True)
                moved_manifest = moved_repo / manifest.relative_to(REPO_ROOT)
                moved_manifest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(manifest, moved_manifest)
                self.assertEqual(
                    validate_figure_manifest(moved_manifest, moved_repo),
                    [],
                    f"failed after checkout move: {manifest}",
                )


if __name__ == "__main__":
    unittest.main()
