from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"
sys.path.insert(0, str(TOOLS_DIR))

from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class SiFiguresManagerTest(unittest.TestCase):
    def test_default_modules_include_frozen_si_figures(self) -> None:
        default_line = next(
            line
            for line in (REPO_ROOT / "Manager.sh").read_text().splitlines()
            if line.startswith('modules="')
        )
        self.assertIn("si_figures", default_line)

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
        self.assertIn("generate_supplementary_figures.R", result.stdout)
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
            for spec in specs:
                source = run_root / str(spec["source"])
                source.write_bytes(f"fixture {spec['panel']}\n".encode())

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


if __name__ == "__main__":
    unittest.main()
