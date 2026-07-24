from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(TOOLS_DIR))

from figure_output_contract import (  # noqa: E402
    FIGURE_MANIFEST_COLUMNS,
    MODULE_MANIFEST_COLUMNS,
    sha256_file,
    write_tsv,
)
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    SI_FIGURE_PANEL_FILES,
    SI_FIGURE_TABLE_CACHE_FILES,
)


class SiFiguresMaterializationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name) / "repo"
        self.repo.mkdir()
        self.source_id = "si_analysis"
        self.operation_id = "si_materialize"
        self.run_root = (
            self.repo
            / "Results/in-vivo/SI_figures/runs"
            / f"{self.source_id}_si_figures"
        )
        for subdir in ("figures", "metadata", "tables"):
            (self.run_root / subdir).mkdir(parents=True, exist_ok=True)
        self.specs = [spec for spec in PANEL_SPECS if spec["module"] == "si_figures"]
        for spec in self.specs:
            source = self.run_root / str(spec["source"])
            source.write_bytes(f"fixture {spec['panel']}\n".encode())
        self.input_paths = self._write_inputs()
        self.output_tables = self._write_contract()
        self._write_manifests()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _write_inputs(self) -> list[Path]:
        paths = [
            self.repo / "Code/in-vivo/SI_figures/generate_supplementary_figures.R",
            self.repo / "Code/in-vivo/figure7/figure7_config.yaml",
            self.repo / "Data/in-vivo/seurat_metadata.csv",
            self.repo / "Data/in-vivo/scvelo_cell_metrics.csv",
            self.repo / "Data/in-vivo/all_ploidy.tsv",
            self.repo
            / "Data/in-vivo/figure7/raw/zenodo_21463392"
            / "integrated_sct_cca_seurat_final_reclustered.rds",
        ]
        for path in paths:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(f"fixture for {path.name}\n")
        return paths

    def _write_contract(self) -> list[Path]:
        canonical = self.run_root / "tables/si_figures_cell_metadata.csv"
        canonical.write_text(
            "cell_id,UMAP_1,UMAP_2,sample_id,cluster_id,initial_ploidy,"
            "s_phase_score,endpoint_ploidy,endpoint_file,endpoint_cell_id,"
            "context,included_in_si_figures\n"
            "c1,0,0,s1,6,2N,0.1,2.1,h1.sps.cbs,b1,Tumor,TRUE\n"
            "c2,1,1,s2,6,4N,-0.2,,,,CellLine,FALSE\n"
        )
        endpoint = self.run_root / "tables/si_figure6_endpoint_ploidy_join_audit.csv"
        endpoint.write_text(
            "cell,context,initial_ploidy,endpoint_file,endpoint_cell_id,"
            "endpoint_ploidy,matched\n"
            "c1,Tumor,2N,h1.sps.cbs,b1,2.1,TRUE\n"
            "c2,CellLine,4N,h2.sps.cbs,b2,,FALSE\n"
        )
        table_names = [
            "si_figures_cluster_key.tsv",
            "si_figure4_cluster_context_composition.csv",
            "si_figure4_cluster_initial_ploidy_composition.csv",
            "si_figure5_cluster_composition_by_mouse.csv",
            "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
            "si_figure5_cluster_dose_composition.csv",
            "si_figure5_cluster_initial_ploidy_composition.csv",
            "si_figure7_cluster_Hallmark_ORA_all.csv",
            "si_figure7_cluster_Hallmark_GSEA_all.csv",
            "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
            "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
        ]
        tables = [canonical, endpoint]
        for name in table_names:
            path = self.run_root / "tables" / name
            path.write_text("key\tvalue\nfixture\t1\n")
            tables.append(path)

        run_config = self.run_root / "metadata/run_config.tsv"
        write_tsv(
            run_config,
            [
                {"key": "module", "value": "si_figures"},
                {"key": "figures", "value": "4,5,6,7"},
                {"key": "logical_panels", "value": str(len(SI_FIGURE_PANEL_FILES))},
                {"key": "figure_file_count", "value": str(2 * len(SI_FIGURE_PANEL_FILES))},
                {"key": "si_figure7_from_raw_rds", "value": "true"},
                {"key": "table_mode", "value": "full_reanalysis"},
                {"key": "raw_seurat_loaded", "value": "true"},
                {"key": "deg_analysis_executed", "value": "true"},
            ],
            ["key", "value"],
        )
        qc = self.run_root / "metadata/input_qc.tsv"
        write_tsv(
            qc,
            [
                {"key": "seurat_cells", "value": "2"},
                {"key": "tumor_cells", "value": "1"},
                {"key": "cellline_cells", "value": "1"},
                {"key": "tumor_cells_missing_endpoint_ploidy", "value": "0"},
                {"key": "cellline_cells_with_endpoint_ploidy", "value": "0"},
                {"key": "panel_files", "value": str(2 * len(SI_FIGURE_PANEL_FILES))},
            ],
            ["key", "value"],
        )
        provenance = self.run_root / "metadata/si_figures_provenance.tsv"
        write_tsv(
            provenance,
            [
                {"key": key, "value": "a" * 64}
                for key in (
                    "script_sha256",
                    "seurat_metadata_sha256",
                    "scvelo_metrics_sha256",
                    "all_ploidy_sha256",
                    "seurat_rds_sha256",
                    "figure7_config_sha256",
                )
            ],
            ["key", "value"],
        )
        return tables + [qc, provenance]

    def _manifest_row(
        self, path: Path, role: str, source_kind: str, panel: str = ""
    ) -> dict[str, str | int]:
        relative = str(path.relative_to(self.repo))
        return {
            "path": relative,
            "repo_relative_path": relative,
            "absolute_path": str(path),
            "role": role,
            "source_kind": source_kind,
            "module": "si_figures",
            "generated_by": "Manager.sh",
            "command_id": self.source_id,
            "sha256": sha256_file(path),
            "checksum_unavailable_reason": "",
            "byte_size": path.stat().st_size,
            "mtime_utc": "2026-07-24T00:00:00+00:00",
            "figure": "Supplementary" if panel else "",
            "panel": panel,
            "notes": "fixture",
        }

    def _write_manifests(self) -> None:
        write_tsv(
            self.run_root / "metadata/panel_contract.tsv",
            [
                {
                    "panel_id": spec["panel"],
                    "filename": Path(str(spec["source"])).name,
                    "caption_role": spec["caption_role"],
                }
                for spec in self.specs
            ],
            ["panel_id", "filename", "caption_role"],
        )
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [self._manifest_row(path, "input_data", "input_file") for path in self.input_paths],
            MODULE_MANIFEST_COLUMNS,
        )
        output_rows = [
            self._manifest_row(
                self.run_root / str(spec["source"]),
                "output_figure",
                "generated_panel",
                str(spec["panel"]),
            )
            for spec in self.specs
        ]
        output_rows.extend(
            self._manifest_row(path, "output_table", "generated_table")
            for path in self.output_tables
        )
        write_tsv(
            self.run_root / "metadata/output_manifest.tsv",
            output_rows,
            MODULE_MANIFEST_COLUMNS,
        )

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
                f"si_figures={self.run_root}",
                "--touched-manifest-list",
                str(self.repo / "manager/touched_figure_manifests.txt"),
                "--overwrite",
            ],
            text=True,
            capture_output=True,
        )

    def test_materializes_exact_si_figures_contract(self) -> None:
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        manifest = self.repo / "figures/Supplementary/manifest.tsv"
        touched = self.repo / "manager/touched_figure_manifests.txt"
        self.assertEqual(touched.read_text().splitlines(), [str(manifest.resolve())])
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual({row["panel"] for row in rows}, {str(spec["panel"]) for spec in self.specs})
        for spec in self.specs:
            published = manifest.parent / str(spec["asset"])
            self.assertEqual(
                published.read_bytes(),
                (self.run_root / str(spec["source"])).read_bytes(),
            )

    def test_preserves_existing_supplementary_manifest_rows(self) -> None:
        out_dir = self.repo / "figures/Supplementary"
        out_dir.mkdir(parents=True)
        existing_asset = out_dir / "panel_SuppFig1A_existing.png"
        existing_asset.write_bytes(b"existing")
        existing_row = {column: "" for column in FIGURE_MANIFEST_COLUMNS}
        existing_row.update(
            {
                "figure": "Supplementary",
                "panel": "SuppFig1A",
                "asset_path": str(existing_asset.relative_to(self.repo)),
                "asset_status": "generated",
            }
        )
        write_tsv(
            out_dir / "manifest.tsv",
            [existing_row],
            FIGURE_MANIFEST_COLUMNS,
        )
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (out_dir / "manifest.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertIn("SuppFig1A", {row["panel"] for row in rows})
        self.assertTrue(existing_asset.is_file())

    def test_rejects_unexpected_figure_file(self) -> None:
        (self.run_root / "figures/unexpected.png").write_bytes(b"unexpected")
        result = self._run()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exact panel inventory", result.stderr)

    def test_accepts_plot_only_run_from_complete_table_cache(self) -> None:
        cache_root = self.repo / "Data/in-vivo/SIfigures"
        cache_root.mkdir(parents=True)
        cache_paths = []
        for name in sorted(SI_FIGURE_TABLE_CACHE_FILES):
            path = cache_root / name
            path.write_text(f"fixture for {name}\n")
            cache_paths.append(path)
        script_and_config = self.input_paths[:2]
        write_tsv(
            self.run_root / "metadata/run_config.tsv",
            [
                {"key": "module", "value": "si_figures"},
                {"key": "figures", "value": "4,5,6,7"},
                {"key": "logical_panels", "value": str(len(SI_FIGURE_PANEL_FILES))},
                {"key": "figure_file_count", "value": str(2 * len(SI_FIGURE_PANEL_FILES))},
                {"key": "si_figure7_from_raw_rds", "value": "false"},
                {"key": "table_mode", "value": "canonical_cache"},
                {"key": "raw_seurat_loaded", "value": "false"},
                {"key": "deg_analysis_executed", "value": "false"},
            ],
            ["key", "value"],
        )
        write_tsv(
            self.run_root / "metadata/si_figures_provenance.tsv",
            [
                {"key": "script_sha256", "value": "a" * 64},
                {"key": "figure7_config_sha256", "value": "a" * 64},
                {
                    "key": "seurat_metadata_sha256",
                    "value": "not_read_in_plot_only_mode",
                },
                {
                    "key": "scvelo_metrics_sha256",
                    "value": "not_read_in_plot_only_mode",
                },
                {
                    "key": "all_ploidy_sha256",
                    "value": "not_read_in_plot_only_mode",
                },
                {"key": "seurat_rds_sha256", "value": "not_run"},
            ],
            ["key", "value"],
        )
        self._write_manifests()
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [
                self._manifest_row(path, "input_data", "input_file")
                for path in [*script_and_config, *cache_paths]
            ],
            MODULE_MANIFEST_COLUMNS,
        )
        result = self._run()
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
