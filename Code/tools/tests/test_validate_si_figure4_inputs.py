from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sys

TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))
from validate_si_figure4_inputs import validate


class ValidateSiFigure4InputsTest(unittest.TestCase):
    @staticmethod
    def _write_pair(root: Path, scvelo_cell: str = "c1") -> tuple[Path, Path]:
        seurat = root / "seurat_metadata.csv"
        scvelo = root / "scvelo_cell_metrics.csv"
        seurat.write_text(
            "cell,UMAP_1,UMAP_2,Dose,cluster_final,sample\n"
            "c1,0.5,-0.5,0mg/kg,6,s1\n"
        )
        scvelo.write_text(
            "cell,TN,cluster_final,sample,Ploidy,Dose\n"
            f"{scvelo_cell},Tumor,6,s1,2N,0mg/kg\n"
        )
        return seurat, scvelo

    def test_accepts_a_complete_cell_aligned_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            seurat, scvelo = self._write_pair(Path(tmp))
            validate(seurat, scvelo)

    def test_rejects_scvelo_cells_missing_from_seurat_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            seurat, scvelo = self._write_pair(Path(tmp), scvelo_cell="missing")
            with self.assertRaisesRegex(ValueError, "missing from Seurat metadata"):
                validate(seurat, scvelo)

    def test_rejects_nonfinite_umap_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seurat, scvelo = self._write_pair(root)
            seurat.write_text(
                "cell,UMAP_1,UMAP_2,Dose,cluster_final,sample\n"
                "c1,NaN,-0.5,0mg/kg,6,s1\n"
            )
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                validate(seurat, scvelo)


if __name__ == "__main__":
    unittest.main()
