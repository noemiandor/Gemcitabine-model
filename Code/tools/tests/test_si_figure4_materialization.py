from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

from figure_output_contract import MODULE_MANIFEST_COLUMNS, sha256_file, write_tsv  # noqa: E402
from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class SiFigure4MaterializationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name) / "repo"
        self.repo.mkdir()
        self.source_id = "si_analysis"
        self.operation_id = "si_materialize"
        self.run_root = (
            self.repo
            / "Results/in-vivo/SI_figure4/runs"
            / f"{self.source_id}_si_figure4"
        )
        (self.run_root / "figures").mkdir(parents=True)
        (self.run_root / "metadata").mkdir()
        self.specs = [spec for spec in PANEL_SPECS if spec["module"] == "si_figure4"]
        for spec in self.specs:
            source = self.run_root / str(spec["source"])
            source.write_bytes(f"fixture {spec['panel']}\n".encode())
        self.input_path = self.repo / "Data/in-vivo/seurat_metadata.csv"
        self.input_path.parent.mkdir(parents=True)
        self.input_path.write_text("cell,UMAP_1,UMAP_2\nc1,0,0\n")
        self._write_metadata()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _write_metadata(self) -> None:
        write_tsv(
            self.run_root / "metadata/run_config.tsv",
            [{"key": "module", "value": "si_figure4"}],
            ["key", "value"],
        )
        write_tsv(
            self.run_root / "metadata/panel_contract.tsv",
            [
                {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                for spec in self.specs
            ],
            ["panel_id", "filename"],
        )
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [self._manifest_row(self.input_path, "input_data", "input_file", "")],
            MODULE_MANIFEST_COLUMNS,
        )
        write_tsv(
            self.run_root / "metadata/output_manifest.tsv",
            [
                self._manifest_row(
                    self.run_root / str(spec["source"]),
                    "output_figure",
                    "generated_panel",
                    str(spec["panel"]),
                )
                for spec in self.specs
            ],
            MODULE_MANIFEST_COLUMNS,
        )

    def _manifest_row(
        self, path: Path, role: str, source_kind: str, panel: str
    ) -> dict[str, str | int]:
        relative = str(path.relative_to(self.repo))
        return {
            "path": relative,
            "repo_relative_path": relative,
            "absolute_path": str(path),
            "role": role,
            "source_kind": source_kind,
            "module": "si_figure4",
            "generated_by": "Manager.sh",
            "command_id": self.source_id,
            "sha256": sha256_file(path),
            "checksum_unavailable_reason": "",
            "byte_size": path.stat().st_size,
            "mtime_utc": "2026-07-22T00:00:00+00:00",
            "figure": "Supplementary" if panel else "",
            "panel": panel,
            "notes": "fixture",
        }

    def _run(self) -> subprocess.CompletedProcess[str]:
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
                f"si_figure4={self.run_root}",
                "--overwrite",
            ],
            text=True,
            capture_output=True,
        )

    def test_materializes_exact_si_figure4_contract(self) -> None:
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.repo / "figures/Supplementary/manifest.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual([row["panel"] for row in rows], [str(spec["panel"]) for spec in self.specs])
        self.assertEqual({row["run_id"] for row in rows}, {self.source_id})
        for spec in self.specs:
            published = self.repo / "figures/Supplementary" / str(spec["asset"])
            self.assertEqual(published.read_bytes(), (self.run_root / str(spec["source"])).read_bytes())

    def test_rejects_an_unexpected_figure_file(self) -> None:
        (self.run_root / "figures/unexpected.png").write_bytes(b"unexpected")
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exact panel inventory", result.stderr)


if __name__ == "__main__":
    unittest.main()
