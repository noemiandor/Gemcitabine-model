from __future__ import annotations

import csv
import shutil
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"

import sys

sys.path.insert(0, str(TOOLS_DIR))
from validate_si_figures_table_cache import (  # noqa: E402
    EXPECTED_FILES,
    validate_cache,
)


class SiFiguresCacheTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.canonical = REPO_ROOT / "Data/in-vivo/SIfigures"

    def copy_cache(self, destination: Path) -> Path:
        cache = destination / "SIfigures"
        shutil.copytree(self.canonical, cache)
        return cache

    def test_canonical_cache_is_exact_and_valid(self) -> None:
        self.assertEqual(validate_cache(self.canonical), [])
        observed = {
            path.name
            for path in self.canonical.iterdir()
            if path.is_file() and path.name != "manifest.tsv"
        }
        self.assertEqual(observed, set(EXPECTED_FILES))
        self.assertEqual(len(observed), 11)

    def test_rejects_missing_and_unexpected_tables(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = self.copy_cache(Path(tmp))
            (cache / EXPECTED_FILES[0]).unlink()
            (cache / "raw_workflow_intermediate.csv").write_text("value\n1\n")
            errors = validate_cache(cache)
        self.assertTrue(any("missing expected" in error for error in errors))
        self.assertTrue(any("unexpected cache" in error for error in errors))

    def test_rejects_nonfinite_si7_matrix(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = self.copy_cache(Path(tmp))
            matrix = (
                cache
                / "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv"
            )
            with matrix.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                headers = list(rows[0])
            rows[0]["0"] = "not-a-number"
            with matrix.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=headers,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
            errors = validate_cache(cache)
        self.assertTrue(any("matrix values must be finite" in error for error in errors))
        self.assertTrue(any("SHA-256 mismatch" in error for error in errors))

    def test_rejects_nonportable_manifest_filename(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            cache = self.copy_cache(Path(tmp))
            manifest = cache / "manifest.tsv"
            with manifest.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
                headers = list(rows[0])
            rows[0]["filename"] = "/hpc/tao/" + rows[0]["filename"]
            with manifest.open("w", newline="") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=headers,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                writer.writerows(rows)
            errors = validate_cache(cache)
        self.assertTrue(any("portable basenames" in error for error in errors))


if __name__ == "__main__":
    unittest.main()
