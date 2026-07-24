from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = TOOLS_DIR / "validate_si_figures_table_cache.py"
SPEC = importlib.util.spec_from_file_location("validate_si_figures_table_cache", MODULE_PATH)
assert SPEC and SPEC.loader
CACHE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CACHE)


def write_table(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


def create_valid_cache(root: Path) -> None:
    write_table(
        root / "si_figures_cell_metadata.csv",
        [
            "cell_id", "UMAP_1", "UMAP_2", "sample_id", "cluster_id",
            "cluster_annotation", "initial_ploidy", "dose", "s_phase_score",
            "endpoint_ploidy", "context", "included_in_si_figures",
        ],
        [
            {
                "cell_id": "tumor1", "UMAP_1": 0, "UMAP_2": 1, "sample_id": "m1",
                "cluster_id": "0", "cluster_annotation": "Cycling", "initial_ploidy": "2N",
                "dose": "0mg/kg", "s_phase_score": 0.1, "endpoint_ploidy": 2.1,
                "context": "Tumor", "included_in_si_figures": "TRUE",
            },
            {
                "cell_id": "line1", "UMAP_1": 2, "UMAP_2": 3, "sample_id": "line",
                "cluster_id": "2", "cluster_annotation": "Stress", "initial_ploidy": "4N",
                "dose": "0mg/kg", "s_phase_score": -0.1, "endpoint_ploidy": "",
                "context": "CellLine", "included_in_si_figures": "FALSE",
            },
        ],
    )
    write_table(
        root / "si_figures_cluster_key.tsv",
        ["cluster_id", "cluster_annotation", "cellcycle_classification", "cluster_order", "color"],
        [
            {
                "cluster_id": cluster, "cluster_annotation": f"Cluster {cluster}",
                "cellcycle_classification": "CellCycle", "cluster_order": index,
                "color": "#000000",
            }
            for index, cluster in enumerate(CACHE.CLUSTERS, 1)
        ],
    )
    composition_rows = [
        {
            "group_value": "0", "fill_value": "Tumor", "n_cells": 1,
            "total_cells": 1, "proportion": 1,
        }
    ]
    for name in (
        "si_figure4_cluster_context_composition.csv",
        "si_figure4_cluster_initial_ploidy_composition.csv",
        "si_figure5_cluster_dose_composition.csv",
        "si_figure5_cluster_initial_ploidy_composition.csv",
    ):
        write_table(root / name, list(composition_rows[0]), composition_rows)
    write_table(
        root / "si_figure5_cluster_composition_by_mouse.csv",
        ["mouse", "initial_ploidy", "dose", "cluster_final", "n_cells", "total_cells", "proportion"],
        [{
            "mouse": "m1", "initial_ploidy": "2N", "dose": "0mg/kg",
            "cluster_final": "0", "n_cells": 1, "total_cells": 1, "proportion": 1,
        }],
    )
    write_table(
        root / "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
        ["initial_ploidy", "dose", "cluster_final", "n_mice", "mean_proportion"],
        [{
            "initial_ploidy": "2N", "dose": "0mg/kg", "cluster_final": "0",
            "n_mice": 1, "mean_proportion": 1,
        }],
    )
    write_table(
        root / "si_figure6_endpoint_ploidy_join_audit.csv",
        ["cell", "context", "initial_ploidy", "endpoint_file", "endpoint_cell_id", "endpoint_ploidy", "matched"],
        [
            {
                "cell": "tumor1", "context": "Tumor", "initial_ploidy": "2N",
                "endpoint_file": "m1.sps.cbs", "endpoint_cell_id": "tumor1",
                "endpoint_ploidy": 2.1, "matched": "TRUE",
            },
            {
                "cell": "line1", "context": "CellLine", "initial_ploidy": "4N",
                "endpoint_file": "", "endpoint_cell_id": "",
                "endpoint_ploidy": "", "matched": "FALSE",
            },
        ],
    )
    for cluster in CACHE.CLUSTERS:
        write_table(
            root / f"si_figure7_cluster_{cluster}_vs_rest_DEG.csv",
            ["gene", "gene_symbol", "cluster"],
            [{"gene": f"gene_{cluster}", "gene_symbol": f"GENE_{cluster}", "cluster": cluster}],
        )
        write_table(
            root / f"si_figure7_cluster_{cluster}_top100_up_ORA_input.csv",
            ["gene_symbol"],
            [{"gene_symbol": f"GENE_{cluster}"}],
        )
    matrix_rows = [
        {"pathway": f"PATHWAY_{index}", **{cluster: index / 10 for cluster in CACHE.CLUSTERS}}
        for index in range(1, 21)
    ]
    for name in (
        "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
        "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    ):
        write_table(root / name, ["pathway", *CACHE.CLUSTERS], matrix_rows)
    for name in (
        "si_figure7_cluster_Hallmark_ORA_all.csv",
        "si_figure7_cluster_Hallmark_GSEA_all.csv",
        "si_figure7_hallmark_ORA_universe.csv",
    ):
        write_table(root / name, ["key", "value"], [{"key": "fixture", "value": 1}])


class ValidateSiFiguresTableCacheTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        create_valid_cache(self.root)

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def test_accepts_exact_complete_cache(self) -> None:
        self.assertEqual(CACHE.validate_cache(self.root), [])

    def test_rejects_missing_and_unexpected_files(self) -> None:
        (self.root / CACHE.EXPECTED_FILES[0]).unlink()
        (self.root / "unexpected.csv").write_text("x\n1\n")
        errors = CACHE.validate_cache(self.root)
        self.assertTrue(any("missing expected" in error for error in errors))
        self.assertTrue(any("unexpected cache" in error for error in errors))

    def test_rejects_nonfinite_umap(self) -> None:
        path = self.root / "si_figures_cell_metadata.csv"
        text = path.read_text().replace("tumor1,0,1", "tumor1,not-a-number,1")
        path.write_text(text)
        self.assertTrue(any("UMAP_1 must be finite" in error for error in CACHE.validate_cache(self.root)))

    def test_writes_checksum_manifest(self) -> None:
        manifest = self.root / "manifest.tsv"
        CACHE.write_manifest(self.root, manifest, "source_run")
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 32)
        self.assertEqual({row["filename"] for row in rows}, set(CACHE.EXPECTED_FILES))
        self.assertEqual({row["source_run"] for row in rows}, {"source_run"})


if __name__ == "__main__":
    unittest.main()
