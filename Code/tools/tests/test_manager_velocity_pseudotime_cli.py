from __future__ import annotations

import csv
import hashlib
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]


class ManagerVelocityPseudotimeCliTest(unittest.TestCase):
    @staticmethod
    def sha256(path: Path) -> str:
        digest = hashlib.sha256()
        digest.update(path.read_bytes())
        return digest.hexdigest()

    @staticmethod
    def write_key_values(path: Path, rows: dict[str, str]) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["key", "value"],
                delimiter="\t",
                lineterminator="\n",
            )
            writer.writeheader()
            for key, value in rows.items():
                writer.writerow({"key": key, "value": value})

    def run_manager(self, *args: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPO_ROOT / "Manager.sh"), *args],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )

    def test_module_is_opt_in_until_the_frozen_table_is_reviewed(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text(encoding="utf-8")
        default_line = next(
            line for line in manager_text.splitlines() if line.startswith('modules="')
        )
        self.assertNotIn("in_vivo_velocity_pseudotime", default_line)

    def test_full_refit_check_only_wires_the_raw_fallback_without_downloading(self) -> None:
        result = self.run_manager(
            "--mode",
            "full-refit",
            "--modules",
            "in_vivo_velocity_pseudotime",
            "--run-id",
            "velocity_panel_contract",
            "--dry-run",
            "--figure7-no-download-missing-raw",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(
            "Code/in-vivo/velocity_pseudotime/run_velocity_pseudotime.R",
            result.stdout,
        )
        self.assertIn("--mode=full-workflow", result.stdout)
        self.assertIn("--download-missing-raw=false", result.stdout)
        self.assertIn("--jobs=1", result.stdout)
        self.assertIn(
            "Results/in-vivo/velocity_pseudotime/intermediates",
            result.stdout,
        )

    def test_materializer_registers_both_publication_derivatives(self) -> None:
        text = (
            REPO_ROOT / "Code/tools/materialize_figure_assets.py"
        ).read_text(encoding="utf-8")
        self.assertIn('"module": "in_vivo_velocity_pseudotime"', text)
        self.assertIn("panel_SuppFig10_velocity_pseudotime.pdf", text)
        self.assertIn("panel_SuppFig10_velocity_pseudotime.png", text)

    def run_materializer_with_config(
        self,
        rows: dict[str, str],
    ) -> subprocess.CompletedProcess[str]:
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        root = Path(temporary.name)
        run_id = "velocity_guard"
        run_root = root / f"{run_id}_velocity_pseudotime"
        metadata = run_root / "metadata"
        metadata.mkdir(parents=True)
        with (metadata / "run_config.tsv").open("w", encoding="utf-8") as handle:
            handle.write("key\tvalue\n")
            for key, value in rows.items():
                handle.write(f"{key}\t{value}\n")
        return subprocess.run(
            [
                "python3",
                str(REPO_ROOT / "Code/tools/materialize_figure_assets.py"),
                "--repo-root",
                str(REPO_ROOT),
                "--figure-root",
                str(root / "figures"),
                "--source-run-id",
                run_id,
                "--module-run",
                f"in_vivo_velocity_pseudotime={run_root}",
                "--overwrite",
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )

    def test_direct_materialization_rejects_raw_generated_output(self) -> None:
        result = self.run_materializer_with_config(
            {
                "module": "in_vivo_velocity_pseudotime",
                "source_kind": "raw_regenerated_unreviewed",
                "source_table": (
                    "Results/in-vivo/velocity_pseudotime/intermediates/"
                    "cellcycle_velocity_pseudotime_umap.tsv"
                ),
                "source_provenance": (
                    "Results/in-vivo/velocity_pseudotime/intermediates/"
                    "cellcycle_velocity_pseudotime_umap.provenance.tsv"
                ),
                "root_cluster": "6",
                "n_cells": "2881",
                "scvelo_mode": "stochastic",
                "canonical_publication_allowed": "false",
            }
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("exact reviewed repo-relative frozen bundle", result.stderr)

    def test_direct_materialization_rejects_absolute_or_overridden_bundle(self) -> None:
        canonical = {
            "module": "in_vivo_velocity_pseudotime",
            "source_kind": "reviewed_frozen_table",
            "root_cluster": "6",
            "n_cells": "2881",
            "scvelo_mode": "stochastic",
            "canonical_publication_allowed": "true",
        }
        cases = [
            ("/tmp/velocity.tsv", "/tmp/velocity.provenance.tsv"),
            (
                "Data/in-vivo/figure7/processed/alternate_velocity.tsv",
                "Data/in-vivo/figure7/processed/alternate_velocity.provenance.tsv",
            ),
        ]
        for table, provenance in cases:
            with self.subTest(table=table):
                result = self.run_materializer_with_config(
                    {
                        **canonical,
                        "source_table": table,
                        "source_provenance": provenance,
                    }
                )
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(
                    "exact reviewed repo-relative frozen bundle",
                    result.stderr,
                )

    def test_direct_materialization_accepts_only_a_complete_canonical_contract(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo_root = Path(temporary) / "repo"
            run_id = "velocity_canonical"
            run_root = (
                repo_root
                / "Results/in-vivo/velocity_pseudotime/runs"
                / f"{run_id}_velocity_pseudotime"
            )
            table = (
                repo_root
                / "Data/in-vivo/figure7/processed/"
                "cellcycle_velocity_pseudotime_umap.tsv"
            )
            provenance = table.with_name(
                "cellcycle_velocity_pseudotime_umap.provenance.tsv"
            )
            table.parent.mkdir(parents=True)
            columns = [
                "cell_id",
                "sample_id",
                "cluster",
                "UMAP_1",
                "UMAP_2",
                "velocity_UMAP_1",
                "velocity_UMAP_2",
                "velocity_pseudotime",
                "is_root",
            ]
            clusters = ["4c"] * 89 + ["6"] * 686 + ["10"] * 2106
            with table.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=columns,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for index, cluster in enumerate(clusters):
                    writer.writerow(
                        {
                            "cell_id": f"cell-{index:05d}",
                            "sample_id": f"mouse-{index % 16:02d}",
                            "cluster": cluster,
                            "UMAP_1": str(index / 100),
                            "UMAP_2": str((index % 97) / 10),
                            "velocity_UMAP_1": "0.1",
                            "velocity_UMAP_2": "-0.05",
                            "velocity_pseudotime": str(index / 2880),
                            "is_root": "TRUE" if cluster == "6" else "FALSE",
                        }
                    )
            frozen_rows = {
                "schema_version": "1",
                "artifact": "figure7_cellcycle_velocity_umap",
                "table_sha256": self.sha256(table),
                "root_cluster": "6",
                "cellcycle_clusters": "4c,6,10",
                "n_cells": "2881",
                "cluster_4c_cells": "89",
                "cluster_6_cells": "686",
                "cluster_10_cells": "2106",
                "scvelo_mode": "stochastic",
                "velocity_basis": "reviewed_seurat_umap",
                "velocity_source": (
                    "scvelo.tl.velocity_embedding_from_velocity_graph"
                ),
                "embedding_source_sha256": "a" * 64,
                "embedding_lineage_sha256": "b" * 64,
                "embedding_inventory_sha256": "c" * 64,
                "cellcycle_source_sha256": "d" * 64,
                "si_metadata_source_sha256": "e" * 64,
                "builder_sha256": "f" * 64,
                "panel_logic_sha256": "1" * 64,
                "scvelo_generator_sha256": "2" * 64,
            }
            self.write_key_values(provenance, frozen_rows)

            figures = run_root / "figures"
            tables = run_root / "tables"
            metadata = run_root / "metadata"
            figures.mkdir(parents=True)
            tables.mkdir(parents=True)
            pdf = figures / "panel_SuppFig10_velocity_pseudotime.pdf"
            png = figures / "panel_SuppFig10_velocity_pseudotime.png"
            pdf.write_bytes(b"%PDF-1.4\nvelocity contract\n")
            png.write_bytes(b"velocity PNG contract\n")
            plot_table = tables / "panel_velocity_pseudotime_plot_data.tsv"
            shutil.copy2(table, plot_table)
            grid = tables / "panel_velocity_pseudotime_vector_grid.tsv"
            grid_columns = [
                "UMAP_1",
                "UMAP_2",
                "velocity_UMAP_1",
                "velocity_UMAP_2",
                "n_cells",
                "vector_length",
                "UMAP_1_to",
                "UMAP_2_to",
            ]
            with grid.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=grid_columns,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()
                for index in range(12):
                    writer.writerow(
                        {
                            "UMAP_1": str(index),
                            "UMAP_2": str(index / 2),
                            "velocity_UMAP_1": "0.1",
                            "velocity_UMAP_2": "0.2",
                            "n_cells": "4",
                            "vector_length": "0.2236068",
                            "UMAP_1_to": str(index + 0.1),
                            "UMAP_2_to": str(index / 2 + 0.2),
                        }
                    )
            run_provenance = {
                **frozen_rows,
                "rendered_pdf_sha256": self.sha256(pdf),
                "rendered_png_sha256": self.sha256(png),
                "plot_table_sha256": self.sha256(plot_table),
                "vector_grid_sha256": self.sha256(grid),
                "vector_grid_rows": "12",
                "root_display_policy": "cluster_6_label_plus_90_percent_core_hull",
            }
            self.write_key_values(
                metadata / "velocity_pseudotime_provenance.tsv",
                run_provenance,
            )
            self.write_key_values(
                metadata / "run_config.tsv",
                {
                    "module": "in_vivo_velocity_pseudotime",
                    "mode": "plot-only",
                    "source_kind": "reviewed_frozen_table",
                    "source_table": table.relative_to(repo_root).as_posix(),
                    "source_table_sha256": self.sha256(table),
                    "source_provenance": provenance.relative_to(repo_root).as_posix(),
                    "source_provenance_sha256": self.sha256(provenance),
                    "root_cluster": "6",
                    "n_cells": "2881",
                    "scvelo_mode": "stochastic",
                    "canonical_publication_allowed": "true",
                },
            )
            manifest_script = REPO_ROOT / "Code/tools/write_file_manifest.py"
            subprocess.run(
                [
                    "python3",
                    str(manifest_script),
                    "--manifest-type",
                    "input",
                    "--module",
                    "in_vivo_velocity_pseudotime",
                    "--output",
                    str(metadata / "input_manifest.tsv"),
                    "--generated-by",
                    "Manager.sh",
                    "--command-id",
                    run_id,
                    "--repo-root",
                    str(repo_root),
                    "--portable",
                    "--path",
                    str(table),
                    "--path",
                    str(provenance),
                ],
                check=True,
                cwd=REPO_ROOT,
            )
            subprocess.run(
                [
                    "python3",
                    str(manifest_script),
                    "--manifest-type",
                    "output",
                    "--module",
                    "in_vivo_velocity_pseudotime",
                    "--output",
                    str(metadata / "output_manifest.tsv"),
                    "--generated-by",
                    "Manager.sh",
                    "--command-id",
                    run_id,
                    "--repo-root",
                    str(repo_root),
                    "--locator-root",
                    str(run_root),
                    "--scan-dir",
                    str(run_root),
                ],
                check=True,
                cwd=REPO_ROOT,
            )
            result = subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "Code/tools/materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo_root),
                    "--figure-root",
                    str(repo_root / "figures"),
                    "--source-run-id",
                    run_id,
                    "--module-run",
                    f"in_vivo_velocity_pseudotime={run_root}",
                    "--overwrite",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue(
                (
                    repo_root
                    / "figures/Supplementary/"
                    "panel_SuppFig10_velocity_pseudotime.pdf"
                ).is_file()
            )
            plot_lines = plot_table.read_text(encoding="utf-8").splitlines()
            altered = plot_lines[1].split("\t")
            altered[3] = "9.5"
            plot_lines[1] = "\t".join(altered)
            plot_table.write_text("\n".join(plot_lines) + "\n", encoding="utf-8")
            run_provenance["plot_table_sha256"] = self.sha256(plot_table)
            self.write_key_values(
                metadata / "velocity_pseudotime_provenance.tsv",
                run_provenance,
            )
            subprocess.run(
                [
                    "python3",
                    str(manifest_script),
                    "--manifest-type",
                    "output",
                    "--module",
                    "in_vivo_velocity_pseudotime",
                    "--output",
                    str(metadata / "output_manifest.tsv"),
                    "--generated-by",
                    "Manager.sh",
                    "--command-id",
                    run_id,
                    "--repo-root",
                    str(repo_root),
                    "--locator-root",
                    str(run_root),
                    "--scan-dir",
                    str(run_root),
                ],
                check=True,
                cwd=REPO_ROOT,
            )
            altered_result = subprocess.run(
                [
                    "python3",
                    str(REPO_ROOT / "Code/tools/materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo_root),
                    "--figure-root",
                    str(repo_root / "altered_figures"),
                    "--source-run-id",
                    run_id,
                    "--module-run",
                    f"in_vivo_velocity_pseudotime={run_root}",
                    "--overwrite",
                ],
                cwd=REPO_ROOT,
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(altered_result.returncode, 0)
            self.assertIn("not byte-identical", altered_result.stderr)


if __name__ == "__main__":
    unittest.main()
