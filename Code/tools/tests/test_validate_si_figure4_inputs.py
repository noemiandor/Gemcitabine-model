from __future__ import annotations

import hashlib
import tempfile
import unittest
from pathlib import Path
import sys

TOOLS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
CONFIG_PATH = REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml"
sys.path.insert(0, str(TOOLS_DIR))
from validate_si_figure4_inputs import validate


class ValidateSiFigure4InputsTest(unittest.TestCase):
    @staticmethod
    def _sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @classmethod
    def _write_bundle(
        cls,
        root: Path,
        scvelo_cell: str = "c1",
    ) -> tuple[Path, Path, Path]:
        seurat = root / "seurat_metadata.csv"
        scvelo = root / "scvelo_cell_metrics.csv"
        provenance = root / "seurat_metadata_provenance.tsv"
        seurat.write_text(
            "cell,UMAP_1,UMAP_2,Dose,clusters,sample,"
            "cluster_cell_cycle_annotation,integrated_snn_res.0.6,"
            "nCount_RNA,nFeature_RNA,percent.mt,scDblFinder.class,scDblFinder.score\n"
            "c1,0.5,-0.5,0mg/kg,6,s1,cell_cycle_candidate,6,"
            "1000,500,2.5,singlet,0.01\n"
        )
        scvelo.write_text(
            "cell,TN,clusters,sample,Ploidy,Dose\n"
            f"{scvelo_cell},Tumor,6,s1,2N,0mg/kg\n"
        )
        rows = {
            "source_seurat_rds": "/data/source.rds",
            "source_seurat_rds_sha256": "a" * 64,
            "source_object_cells": "1",
            "seurat_metadata_sha256": cls._sha256(seurat),
            "scvelo_metrics_sha256": cls._sha256(scvelo),
            "umap_reduction": "umap",
            "cluster_id_field": "clusters",
            "base_cluster_field": "integrated_snn_res.0.6",
            "clustering_resolution": "0.6",
            "cluster_annotation_field": "cluster_cell_cycle_annotation",
            "sample_field": "sample",
            "dose_field": "Dose",
            "ploidy_field": "Ploidy",
            "context_field": "TN",
            "cellcycle_mapping": (
                "cell_cycle_candidate->CellCycle;"
                "not_cell_cycle_candidate->NonCellCycle"
            ),
            "inclusion_context": "Tumor",
            "inclusion_ploidy_levels": "2N,4N",
            "inclusion_dose_levels": "0mg/kg,30mg/kg,120mg/kg",
            "inclusion_required_nonmissing": (
                "cell_id,UMAP_1,UMAP_2,sample_id,cluster_id,"
                "cluster_annotation,initial_ploidy,dose,"
                "cellcycle_classification"
            ),
            "source_qc_fields": (
                "nCount_RNA,nFeature_RNA,percent.mt,"
                "scDblFinder.class,scDblFinder.score"
            ),
            "source_qc_policy": (
                "Use the reviewed Seurat object as provided; SI Figure 4 "
                "applies no additional expression, mitochondrial, or doublet threshold."
            ),
            "figure7_config_sha256": cls._sha256(CONFIG_PATH),
            "export_script_sha256": "b" * 64,
            "source_code_revision": "c" * 40,
        }
        provenance.write_text(
            "key\tvalue\n"
            + "".join(f"{key}\t{value}\n" for key, value in rows.items())
        )
        return seurat, scvelo, provenance

    def test_accepts_a_complete_cell_aligned_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            seurat, scvelo, provenance = self._write_bundle(Path(tmp))
            validate(seurat, scvelo, provenance, CONFIG_PATH)

    def test_rejects_scvelo_cells_missing_from_seurat_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            seurat, scvelo, provenance = self._write_bundle(
                Path(tmp), scvelo_cell="missing"
            )
            with self.assertRaisesRegex(ValueError, "missing from Seurat metadata"):
                validate(seurat, scvelo, provenance, CONFIG_PATH)

    def test_rejects_nonfinite_umap_coordinates(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            seurat, scvelo, provenance = self._write_bundle(root)
            seurat.write_text(
                "cell,UMAP_1,UMAP_2,Dose,clusters,sample,"
                "cluster_cell_cycle_annotation,integrated_snn_res.0.6,"
                "nCount_RNA,nFeature_RNA,percent.mt,scDblFinder.class,scDblFinder.score\n"
                "c1,NaN,-0.5,0mg/kg,6,s1,cell_cycle_candidate,6,"
                "1000,500,2.5,singlet,0.01\n"
            )
            with self.assertRaisesRegex(ValueError, "nonfinite"):
                validate(seurat, scvelo, provenance, CONFIG_PATH)

    def test_rejects_tampered_seurat_metadata(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            seurat, scvelo, provenance = self._write_bundle(Path(tmp))
            seurat.write_text(seurat.read_text().replace("0.5", "0.6"))
            with self.assertRaisesRegex(ValueError, "checksum mismatch"):
                validate(seurat, scvelo, provenance, CONFIG_PATH)


if __name__ == "__main__":
    unittest.main()
