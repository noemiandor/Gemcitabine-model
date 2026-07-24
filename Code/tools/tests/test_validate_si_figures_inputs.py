from __future__ import annotations

import argparse
import tempfile
import unittest
from pathlib import Path
import sys


TOOLS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml"
sys.path.insert(0, str(TOOLS_DIR))
from validate_si_figures_inputs import validate  # noqa: E402


class ValidateSiFiguresInputsTest(unittest.TestCase):
    def _bundle(self, root: Path, metric_cell: str = "c1") -> argparse.Namespace:
        seurat = root / "seurat_metadata.csv"
        scvelo = root / "scvelo_cell_metrics.csv"
        all_ploidy = root / "all_ploidy.tsv"
        seurat.write_text(
            "cell,UMAP_1,UMAP_2,S.Score,Dose,clusters,sample,Ploidy,TN,"
            "harvest,barcode_raw,cluster_cell_cycle_annotation,"
            "integrated_snn_res.0.6\n"
            "c1,0.5,-0.5,0.1,0mg/kg,6,s1,2N,Tumor,h1,b1,"
            "cell_cycle_candidate,6\n"
            "c2,1,1,-0.2,,6,s2,4N,CellLine,h2,b2,"
            "cell_cycle_candidate,6\n"
        )
        scvelo.write_text(
            "cell,TN,clusters,sample,Ploidy,Dose\n"
            f"{metric_cell},Tumor,6,s1,2N,0mg/kg\n"
        )
        all_ploidy.write_text(
            "file\tcell_id\tploidy\n"
            "h1.sps.cbs\tb1\t2.1\n"
        )
        return argparse.Namespace(
            seurat_metadata=seurat,
            scvelo_metrics=scvelo,
            all_ploidy=all_ploidy,
            seurat_rds=None,
            seurat_metadata_provenance=None,
            config=CONFIG_PATH,
        )

    def test_accepts_complete_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            validate(self._bundle(Path(tmp)))

    def test_rejects_scvelo_cell_absent_from_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError, "absent from Seurat metadata"):
                validate(self._bundle(Path(tmp), metric_cell="missing"))

    def test_rejects_missing_tumor_endpoint_ploidy(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = self._bundle(Path(tmp))
            args.all_ploidy.write_text("file\tcell_id\tploidy\nother\tcell\t2.0\n")
            with self.assertRaisesRegex(ValueError, "missing tumor cells"):
                validate(args)

    def test_rejects_cellline_endpoint_match(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            args = self._bundle(Path(tmp))
            args.all_ploidy.write_text(
                "file\tcell_id\tploidy\n"
                "h1.sps.cbs\tb1\t2.1\n"
                "h2.sps.cbs\tb2\t2.2\n"
            )
            with self.assertRaisesRegex(ValueError, "CellLine"):
                validate(args)


if __name__ == "__main__":
    unittest.main()
