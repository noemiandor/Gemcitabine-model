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
        (self.run_root / "tables").mkdir()
        self.specs = [spec for spec in PANEL_SPECS if spec["module"] == "si_figure4"]
        for spec in self.specs:
            source = self.run_root / str(spec["source"])
            source.write_bytes(f"fixture {spec['panel']}\n".encode())
        input_root = self.repo / "Data/in-vivo"
        input_root.mkdir(parents=True)
        self.input_paths = [
            self.repo / "Code/in-vivo/figure7/figure7_config.yaml",
            input_root / "seurat_metadata.csv",
            input_root / "scvelo_cell_metrics.csv",
            input_root / "seurat_metadata_provenance.tsv",
        ]
        self.input_paths[0].parent.mkdir(parents=True)
        for path in self.input_paths:
            path.write_text(f"fixture for {path.name}\n")
        self._write_data_contract()
        self._write_metadata()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _write_data_contract(self) -> None:
        self.canonical_path = self.run_root / "tables/si_figure4_cell_metadata.csv"
        self.cluster_key_path = self.run_root / "tables/si_figure4_cluster_key.tsv"
        self.mouse_path = self.run_root / "tables/cluster_composition_by_mouse.csv"
        self.group_path = self.run_root / "tables/cluster_composition_by_ploidy_dose.csv"
        self.provenance_path = self.run_root / "metadata/si_figure4_provenance.tsv"
        self.canonical_path.write_text(
            "cell_id,UMAP_1,UMAP_2,sample_id,cluster_id,cluster_annotation,"
            "cluster_order,cluster_color,initial_ploidy,treatment,dose,"
            "dose_mg_per_kg,cellcycle_classification,context,"
            "included_in_si_figure4,exclusion_reason\n"
            "c1,0,0,s1,6,cell-cycle,1,#112233,2N,Control,0mg/kg,0,"
            "CellCycle,Tumor,TRUE,\n"
            "c2,1,1,s2,6,cell-cycle,1,#112233,2N,Control,0mg/kg,0,"
            "CellCycle,CellLine,FALSE,context_not_Tumor\n"
        )
        self.cluster_key_path.write_text(
            "cluster_id\tcluster_annotation\tcellcycle_classification\t"
            "cluster_order\tcolor\tn_all_cells\tn_included_tumor_cells\n"
            "6\tcell-cycle\tCellCycle\t1\t#112233\t2\t1\n"
        )
        self.mouse_path.write_text(
            "mouse,initial_ploidy,dose,cluster_final,cluster_annotation,"
            "cluster_order,n_cells,total_cells,proportion,denominator_definition\n"
            "s1,2N,0mg/kg,6,cell-cycle,1,1,1,1,"
            "all included SI Figure 4 tumor cells from the mouse\n"
        )
        self.group_path.write_text(
            "initial_ploidy,dose,cluster_final,cluster_annotation,cluster_order,"
            "n_mice,sum_n_cells,sum_total_cells,mean_proportion,sd_proportion,"
            "min_proportion,max_proportion\n"
            "2N,0mg/kg,6,cell-cycle,1,1,1,1,1,NA,1,1\n"
        )
        provenance = {
            "source_seurat_rds_sha256": "a" * 64,
            "umap_reduction": "umap",
            "cluster_id_field": "clusters",
            "base_cluster_field": "integrated_snn_res.0.6",
            "clustering_resolution": "0.6",
            "cluster_annotation_field": "cluster_cell_cycle_annotation",
            "source_qc_fields": (
                "nCount_RNA,nFeature_RNA,percent.mt,"
                "scDblFinder.class,scDblFinder.score"
            ),
            "source_qc_policy": "reviewed Seurat object as provided",
            "cell_inclusion_context": "Tumor",
            "cell_inclusion_required_nonmissing": (
                "cell,UMAP_1,UMAP_2,sample,cluster,annotation"
            ),
            "canonical_cell_rows": "2",
            "included_cell_rows": "1",
            "excluded_cell_rows": "1",
            "canonical_cell_table_sha256": sha256_file(self.canonical_path),
            "cluster_key_sha256": sha256_file(self.cluster_key_path),
            "composition_by_mouse_sha256": sha256_file(self.mouse_path),
            "composition_by_ploidy_dose_sha256": sha256_file(self.group_path),
            "upstream_analysis_documentation_doi": "10.5281/zenodo.21463392",
            "upstream_analysis_documentation_url": (
                "https://zenodo.org/records/21463392"
            ),
            "upstream_analysis_documentation_statement": (
                "The processed inputs have documented upstream analysis."
            ),
        }
        self.provenance_path.write_text(
            "key\tvalue\n"
            + "".join(f"{key}\t{value}\n" for key, value in provenance.items())
        )

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
            [
                self._manifest_row(path, "input_data", "input_file", "")
                for path in self.input_paths
            ],
            MODULE_MANIFEST_COLUMNS,
        )
        formal_outputs = [
            self.canonical_path,
            self.cluster_key_path,
            self.mouse_path,
            self.group_path,
            self.provenance_path,
        ]
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
            ]
            + [
                self._manifest_row(
                    path,
                    "output_table",
                    "generated_table",
                    "",
                )
                for path in formal_outputs
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

    def test_rejects_tampered_canonical_data(self) -> None:
        self.canonical_path.write_text(
            self.canonical_path.read_text().replace(
                "c1,0,0,s1", "c1,not-a-number,0,s1"
            )
        )
        self._write_metadata()
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Canonical UMAP coordinate is nonnumeric", result.stderr)

    def test_rejects_tampered_composition_denominator(self) -> None:
        self.mouse_path.write_text(
            self.mouse_path.read_text().replace(
                "cell-cycle,1,1,1,1,all included",
                "cell-cycle,1,1,2,0.5,all included",
            )
        )
        self._write_metadata()
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Composition denominator mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
