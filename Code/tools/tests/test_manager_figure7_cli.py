from __future__ import annotations

import csv
import os
import statistics
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
TOOLS_DIR = REPO_ROOT / "Code/tools"
TESTS_DIR = TOOLS_DIR / "tests"

import sys

sys.path.insert(0, str(TOOLS_DIR))
sys.path.insert(0, str(TESTS_DIR))
from figure_output_contract import MODULE_MANIFEST_COLUMNS, sha256_file, write_tsv  # noqa: E402
from figure7_density_fixture import write_density_localization_fixture  # noqa: E402
from materialize_figure_assets import (  # noqa: E402
    PANEL_SPECS,
    panel_specs_for_figure7_variant,
)


class ManagerFigure7CliTest(unittest.TestCase):
    def _run(self, *args: str, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", str(REPO_ROOT / "Manager.sh"), *args],
            cwd=REPO_ROOT,
            env=env,
            text=True,
            capture_output=True,
        )

    def _figure7_input_lineage_paths(
        self,
        run_dir: Path,
        *,
        manager_mode: str = "full-refit",
    ) -> list[str]:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        start = manager_text.index("figure7_runtime_source_paths()")
        end = manager_text.index("\nrequired_input_paths_for_module()", start)
        function_block = manager_text[start:end]
        script = f"""
set -euo pipefail
{function_block}
metadata_value() {{
  local path="$1" key="$2"
  awk -F '\\t' -v expected="${{key}}" '
    $1 == expected {{ count += 1; value = $2 }}
    END {{ if (count != 1 || value == "") exit 1; print value }}
  ' "${{path}}"
}}
mode="$3"
figure7_full_analysis=false
figure7_panels_ae_only=true
figure7_state_pathway_results_root=""
figure7_reference_root=unused
figure7_seurat_rds=""
figure7_gene_set_artifact=""
figure7_raw_seurat_rds=""
figure7_raw_data_dir="$2"
figure7_loom_root=""
figure7_cell_ploidy_input=cell_ploidy.tsv
figure7_endpoint_cbs_score_input=Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv
figure7_endpoint_cbs_score_will_be_derived=false
si_figures_cbs_dir=Data/in-vivo/scRNAseq_Numbat
figure7_sample_info_input=sample_info.xlsx
figure7_growth_curve_input=growth_curve.xlsx
input_paths_for_module in_vivo_figure7 "$1"
"""
        result = subprocess.run(
            [
                "bash",
                "-c",
                script,
                "fixture",
                str(run_dir),
                str(run_dir / "raw"),
                manager_mode,
            ],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        return [line for line in result.stdout.splitlines() if line]

    @staticmethod
    def _write_run_config(run_dir: Path, rows: dict[str, str]) -> None:
        metadata = run_dir / "metadata"
        metadata.mkdir(parents=True)
        with (metadata / "run_config.tsv").open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["key", "value"])
            writer.writerows(rows.items())

    @staticmethod
    def _copy_reviewed_reference(run_dir: Path) -> None:
        reference_root = (
            REPO_ROOT
            / "Data/in-vivo/figure7/saved_state_pathway"
            / "state_pathway_grch_human_only_initial_ploidy_day17_pointwise_v4"
        )
        (run_dir / "tables").mkdir(exist_ok=True)
        for source in reference_root.glob("*.tsv"):
            destination = (
                run_dir / "metadata" / source.name
                if source.name == "state_pathway_provenance.tsv"
                else run_dir / "tables" / source.name
            )
            destination.write_bytes(source.read_bytes())

    def test_panels_only_requires_source_run_id(self) -> None:
        result = self._run(
            "--mode", "panels-only", "--modules", "in_vivo_figure7", "--run-id", "operation"
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--source-run-id is required", result.stderr)

    def test_default_manuscript_modules_include_figure7(self) -> None:
        manager_text = (REPO_ROOT / "Manager.sh").read_text()
        default_line = next(
            line for line in manager_text.splitlines() if line.startswith('modules="')
        )
        self.assertIn("in_vivo_figure7", default_line)

    def test_removed_include_in_vivo_flag_has_migration_error(self) -> None:
        result = self._run("--include-in-vivo")
        self.assertEqual(result.returncode, 2)
        self.assertIn("--include-in-vivo was removed", result.stderr)
        self.assertIn("--modules in_vivo_figure7", result.stderr)

    def test_help_does_not_advertise_removed_in_vivo_flag(self) -> None:
        result = self._run("--help")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("--include-in-vivo", result.stdout)

    def test_removed_in_vivo_module_has_migration_error(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo",
            "--run-id", "removed_in_vivo",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("Module in_vivo was removed", result.stderr)
        self.assertIn("TGI-AUC and CellCycle-subset ploidy analyses", result.stderr)
        self.assertIn("--modules in_vivo_figure7", result.stderr)

    def test_source_run_id_is_rejected_outside_panels_only(self) -> None:
        result = self._run(
            "--mode", "standard", "--modules", "in_vivo_figure7",
            "--run-id", "operation", "--source-run-id", "source",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("valid only with --mode panels-only", result.stderr)

    def test_full_analysis_requires_pinned_gene_set_artifact(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check", "--figure7-full-analysis",
            "--figure7-seurat-rds", "/tmp/object.rds",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--figure7-gene-set-artifact", result.stderr)

    def test_ae_only_check_omits_panel_f_inputs_and_sets_panel_contract(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check_ae", "--figure7-panels-ae-only",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--panel-set=a-e", result.stdout)
        self.assertIn(
            "--endpoint-cbs-score-input="
            "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv",
            result.stdout,
        )
        self.assertNotIn("--saved-state-pathway-dir", result.stdout)

    def test_tgi24_and_supplement_destination_are_forwarded(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "check_tgi24",
            "--figure7-tgi-day", "24",
            "--figure7-figure-name", "Figure7_Supplement",
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--tgi-day=24", result.stdout)

        specs = [
            spec
            for spec in panel_specs_for_figure7_variant(24, "Figure7_Supplement")
            if spec["module"] == "in_vivo_figure7"
        ]
        self.assertTrue(all(spec["figure"] == "Figure7_Supplement" for spec in specs))
        self.assertTrue(
            any(spec["asset"] == "panel_7A_day24_tgi_calculation.pdf" for spec in specs)
        )
        self.assertFalse(any("day17" in str(spec["asset"]) for spec in specs))

    def test_invalid_figure7_tgi_day_is_rejected(self) -> None:
        result = self._run(
            "--mode", "check-only", "--modules", "in_vivo_figure7",
            "--run-id", "bad_tgi_day", "--figure7-tgi-day", "TGI_24",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("non-negative integer", result.stderr)

    def test_ae_only_rejects_full_panel_f_analysis(self) -> None:
        result = self._run(
            "--mode", "full-refit", "--modules", "in_vivo_figure7",
            "--run-id", "bad_combo", "--figure7-panels-ae-only",
            "--figure7-full-analysis", "--figure7-seurat-rds", "/tmp/object.rds",
            "--figure7-gene-set-artifact", "/tmp/gene_sets.rds",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("mutually exclusive", result.stderr)

    def test_input_lineage_omits_unconsumed_stale_scvelo_cache(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            stale_scvelo = run_dir / "stale_scvelo.csv"
            stale_manifest = run_dir / "scvelo_stage_manifest.tsv"
            stale_scvelo.write_text("stale\n")
            stale_manifest.write_text("key\tvalue\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(stale_scvelo),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertNotIn(str(stale_scvelo), paths)
            self.assertNotIn(str(stale_manifest), paths)

    def test_standard_input_lineage_uses_generic_processed_pair(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            processed = run_dir / "processed"
            processed.mkdir()
            cellcycle = processed / "cellcycle.csv"
            noncellcycle = processed / "noncellcycle.csv"
            legacy_cellcycle = run_dir / "legacy-cellcycle.csv"
            legacy_noncellcycle = run_dir / "legacy-noncellcycle.csv"
            self._write_run_config(
                run_dir,
                {
                    "mode": "standard",
                    "cellcycle_input": str(cellcycle),
                    "cellcycle_sha256": "a" * 64,
                    "noncellcycle_input": str(noncellcycle),
                    "noncellcycle_sha256": "b" * 64,
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(legacy_cellcycle),
                    "workflow_noncellcycle_input": str(legacy_noncellcycle),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )

            paths = self._figure7_input_lineage_paths(
                run_dir,
                manager_mode="standard",
            )

            self.assertEqual(paths.count(str(cellcycle)), 1)
            self.assertEqual(paths.count(str(noncellcycle)), 1)
            self.assertEqual(
                paths.count(
                    "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
                ),
                1,
            )
            self.assertNotIn("cell_ploidy.tsv", paths)
            self.assertNotIn(str(legacy_cellcycle), paths)
            self.assertNotIn(str(legacy_noncellcycle), paths)

    def test_input_lineage_uses_exact_explicit_workflow_seurat_rds(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            explicit_rds = run_dir / "external inputs" / "selected_object.rds"
            explicit_rds.parent.mkdir()
            explicit_rds.write_text("fixture rds\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "state_pathway_support",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(run_dir / "unused_scvelo.csv"),
                    "workflow_seurat_rds": str(explicit_rds),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "a" * 64,
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertIn(str(explicit_rds), paths)
            self.assertNotIn(
                str(run_dir / "raw/integrated_sct_cca_seurat_final_reclustered.rds"),
                paths,
            )

    def test_input_lineage_records_present_reconstruction_manifests(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            upstream = run_dir / "seurat_upstream"
            final_stage = upstream / "03_final_cluster/stage_manifest.tsv"
            reconstruction = upstream / "reconstruction_manifest.tsv"
            final_stage.parent.mkdir(parents=True)
            final_stage.write_text("key\tvalue\nstage\tfinal\n")
            reconstruction.write_text("key\tvalue\nschema_version\t1\n")
            self._write_run_config(
                run_dir,
                {
                    "workflow_executed_stages": "none",
                    "workflow_cellcycle_input": str(run_dir / "cellcycle.csv"),
                    "workflow_noncellcycle_input": str(run_dir / "noncellcycle.csv"),
                    "workflow_scvelo_metrics": str(run_dir / "unused_scvelo.csv"),
                    "workflow_state_pathway_results": "not_applicable",
                    "workflow_state_pathway_reference": "not_applicable",
                    "workflow_seurat_reconstruction_manifest": str(reconstruction),
                    "workflow_seurat_reconstruction_manifest_sha256": "a" * 64,
                    "workflow_seurat_final_stage_manifest": str(final_stage),
                    "workflow_seurat_final_stage_manifest_sha256": "b" * 64,
                    "raw_data_dir": str(run_dir / "raw"),
                    "seurat_rds_sha256": "not_available",
                    "loom_file_count": "0",
                },
            )
            paths = self._figure7_input_lineage_paths(run_dir)
            self.assertIn(str(reconstruction), paths)
            self.assertIn(str(final_stage), paths)
            self.assertIn(
                "Code/in-vivo/figure7/"
                "generate_final_seurat_from_cellranger.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/figure7/src/seurat_upstream.R",
                paths,
            )
            self.assertIn(
                "Code/in-vivo/figure7/src/seurat_upstream_selection.R",
                paths,
            )

            reconstruction.unlink()
            final_stage.unlink()
            archived_paths = self._figure7_input_lineage_paths(run_dir)
            self.assertNotIn(str(reconstruction), archived_paths)
            self.assertNotIn(str(final_stage), archived_paths)
            self.assertIn(
                "Code/in-vivo/figure7/"
                "generate_final_seurat_from_cellranger.R",
                archived_paths,
            )

    def test_panels_only_uses_source_run_without_invoking_r(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            output_root = tmp_path / "Results"
            figure_root = tmp_path / "figures"
            source_id = "source_analysis"
            operation_id = "materialization_operation"
            run_root = output_root / "in-vivo/figure7/runs" / f"{source_id}_figure7"
            (run_root / "figures").mkdir(parents=True)
            (run_root / "metadata").mkdir()
            (run_root / "tables").mkdir()
            self._copy_reviewed_reference(run_root)

            rows = []
            for spec in PANEL_SPECS:
                if (
                    spec["module"] != "in_vivo_figure7"
                    or str(spec["panel"]).startswith("7F")
                    or str(spec["panel"]).startswith("7A-7L_composite")
                ):
                    continue
                source = run_root / str(spec["source"])
                source.write_bytes(f"fake PDF {spec['panel']}\n".encode())
                rows.append(
                    {
                        "path": str(source),
                        "repo_relative_path": "",
                        "absolute_path": str(source),
                        "role": "output_figure",
                        "source_kind": "generated_panel",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(source),
                        "checksum_unavailable_reason": "",
                        "byte_size": source.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "Figure7",
                        "panel": spec["panel"],
                        "notes": "fixture",
                    }
                )
            design = [
                ("2N-A2-0", "2N", "30mg/kg", "30", "SUM159-2N-30-0_harvest.sps.cbs"),
                ("2N-A2-L", "2N", "30mg/kg", "30", "SUM159-2N-30-L_harvest.sps.cbs"),
                ("2N-A4-R", "2N", "120mg/kg", "120", "SUM159-2N-120-R_harvest.sps.cbs"),
                ("2N-A4-RL", "2N", "120mg/kg", "120", "SUM159-2N-120-RL_harvest.sps.cbs"),
                ("A6-4N-O", "4N", "30mg/kg", "30", "SUM159-4N-30-0_harvest.sps.cbs"),
                ("A6-4N-RR", "4N", "30mg/kg", "30", "SUM159-4N-30-RR_harvest.sps.cbs"),
                ("4N-A8-RL", "4N", "120mg/kg", "120", "SUM159-4N-120-RL_harvest.sps.cbs"),
                ("4N-A8-RR", "4N", "120mg/kg", "120", "SUM159-4N-120-RR_harvest.sps.cbs"),
            ]
            endpoint_path = Path(
                "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
            )
            cellcycle_path = Path(
                "Data/in-vivo/figure7/processed/"
                "CellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            )
            noncellcycle_path = Path(
                "Data/in-vivo/figure7/processed/"
                "NonCellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            )
            endpoint_score_by_key: dict[tuple[str, str], float] = {}
            endpoint_files: set[str] = set()
            with (REPO_ROOT / endpoint_path).open(newline="") as handle:
                for row in csv.DictReader(handle, delimiter="\t"):
                    key = (row["file"], row["cell_id"])
                    self.assertNotIn(key, endpoint_score_by_key)
                    endpoint_score_by_key[key] = float(row["ploidy"])
                    endpoint_files.add(row["file"])

            curated_scores_by_sample: dict[str, list[float]] = {}
            curated_file_by_sample: dict[str, str] = {}
            curated_keys: set[tuple[str, str]] = set()
            compartment_counts: dict[str, int] = {}
            tgi_by_sample: dict[str, float] = {}
            for compartment, relative_path in (
                ("CellCycle", cellcycle_path),
                ("NonCellCycle", noncellcycle_path),
            ):
                compartment_count = 0
                with (REPO_ROOT / relative_path).open(newline="") as handle:
                    for row in csv.DictReader(handle):
                        sample_id = row["sample_id"]
                        prefix = f"{sample_id}_"
                        full_cell_id = row["cell_id"]
                        self.assertTrue(full_cell_id.startswith(prefix))
                        endpoint_file = f"{row['growth_curve_harvest']}.sps.cbs"
                        key = (endpoint_file, full_cell_id[len(prefix) :])
                        self.assertNotIn(key, curated_keys)
                        self.assertIn(key, endpoint_score_by_key)
                        canonical_score = endpoint_score_by_key[key]
                        self.assertLess(
                            abs(float(row["cell_ploidy"]) - canonical_score),
                            1e-12,
                        )
                        previous_file = curated_file_by_sample.setdefault(
                            sample_id,
                            endpoint_file,
                        )
                        self.assertEqual(previous_file, endpoint_file)
                        curated_keys.add(key)
                        curated_scores_by_sample.setdefault(sample_id, []).append(
                            canonical_score
                        )
                        compartment_count += 1
                        if compartment == "CellCycle":
                            tgi_by_sample.setdefault(
                                sample_id,
                                float(row["TGI_percent_Day_17"]),
                            )
                compartment_counts[compartment] = compartment_count

            expected_curated_counts = {
                "2N-A1-0": 413,
                "2N-A1-LR": 369,
                "2N-A1-R": 196,
                "2N-A1-RR": 1495,
                "2N-A2-0": 317,
                "2N-A2-L": 505,
                "2N-A4-R": 305,
                "2N-A4-RL": 1280,
                "4N-A5-0": 888,
                "4N-A5-RR": 393,
                "A5-4N-L": 385,
                "A5-4N-R": 358,
                "A6-4N-O": 189,
                "A6-4N-RR": 660,
                "4N-A8-RL": 1832,
                "4N-A8-RR": 247,
            }
            observed_curated_counts = {
                sample_id: len(values)
                for sample_id, values in curated_scores_by_sample.items()
            }
            treated_sample_ids = {row[0] for row in design}
            self.assertEqual(len(endpoint_score_by_key), 14125)
            self.assertEqual(len(endpoint_files), 16)
            self.assertEqual(compartment_counts, {
                "CellCycle": 2881,
                "NonCellCycle": 6951,
            })
            self.assertEqual(observed_curated_counts, expected_curated_counts)
            self.assertEqual(len(curated_keys), 9832)
            self.assertEqual(
                sum(observed_curated_counts[sample_id] for sample_id in treated_sample_ids),
                5335,
            )
            scores = [
                statistics.mean(curated_scores_by_sample[sample_id])
                for sample_id, *_ in design
            ]
            panel_l_plot = run_root / "tables/panel_7E_plot_data.tsv"
            plot_rows = []
            for index, (
                sample_id,
                origin,
                dose,
                dose_mg,
                endpoint_file,
            ) in enumerate(design):
                score = scores[index]
                plot_rows.append(
                    {
                        "sample_id": sample_id,
                        "initial_ploidy": origin,
                        "dose": dose,
                        "dose_mg": dose_mg,
                        "endpoint_ploidy_file": endpoint_file,
                        "sample_mean_endpoint_ploidy": score,
                        "n_endpoint_ploidy_cells": len(
                            curated_scores_by_sample[sample_id]
                        ),
                        "endpoint_ploidy_source_total_cells": "14125",
                        "endpoint_ploidy_score_universe_total_cells": "9832",
                        "endpoint_ploidy_source_file_count": "16",
                        "endpoint_ploidy_source_sha256": sha256_file(
                            REPO_ROOT / endpoint_path
                        ),
                        "endpoint_ploidy_score_policy": (
                            "arithmetic_mean_of_finite_postprocessed_cell_"
                            "ploidy_in_exact_qc_passed_cellcycle_"
                            "noncellcycle_union_per_sample"
                        ),
                        "endpoint_ploidy_mapping_policy": (
                            "exact_processed_sample_barcode_to_canonical_cbs_"
                            "file_cell_and_value;score_universe=reviewed_final_"
                            "seurat_tumor_cells"
                        ),
                        "TGI_percent_Day_17": tgi_by_sample[sample_id],
                        "tgi_outcome": "day",
                        "tgi_day": "17",
                        "tgi_measure": "TGI_percent_Day_17",
                        "matched_control_summary": "mean",
                        "matched_control_group": "initial_ploidy",
                    }
                )

            write_tsv(panel_l_plot, plot_rows, list(plot_rows[0]))
            panel_l_test = run_root / "tables/panel_7E_test.tsv"
            test_row = {
                "n": "8",
                "estimate": "-0.698401019253014",
                "asymptotic_p": "0.0540069781511847",
                "permutation_p_two_sided": "0.0591269841269841",
                "n_permutations": "40320",
                "permutation_mode": "exact_TGI_label_enumeration",
                "permutation_strata": "none",
                "association_type": "unadjusted_mouse_level_pearson",
                "score_variable": (
                    "sample_mean_qc_passed_curated_cbs_cell_ploidy"
                ),
                "score_source_sha256": sha256_file(REPO_ROOT / endpoint_path),
                "score_source_n_cells": "9832",
                "score_inventory_n_cells": "14125",
                "score_source_n_files": "16",
                "treated_score_n_cells": "5335",
                "score_aggregation_policy": (
                    "arithmetic_mean_of_finite_postprocessed_cell_"
                    "ploidy_in_exact_qc_passed_cellcycle_"
                    "noncellcycle_union_per_sample"
                ),
                "sample_mapping_policy": (
                    "exact_processed_sample_barcode_to_canonical_cbs_"
                    "file_cell_and_value;score_universe=reviewed_final_"
                    "seurat_tumor_cells"
                ),
                "score_standardization": "none",
                "adjustment_terms": "none",
                "outcome_variable": "TGI_percent_Day_17",
                "plot_x": "sample_mean_endpoint_ploidy",
                "plot_y": "TGI_percent_Day_17",
                "tgi_outcome": "day",
                "tgi_day": "17",
                "tgi_measure": "TGI_percent_Day_17",
                "matched_control_summary": "mean",
                "matched_control_group": "initial_ploidy",
            }
            write_tsv(panel_l_test, [test_row], list(test_row))
            for source in (panel_l_plot, panel_l_test):
                rows.append(
                    {
                        "path": str(source),
                        "repo_relative_path": "",
                        "absolute_path": str(source),
                        "role": "output_table",
                        "source_kind": "generated_table",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(source),
                        "checksum_unavailable_reason": "",
                        "byte_size": source.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "panel L fixture",
                    }
                )
            for source in write_density_localization_fixture(run_root):
                rows.append(
                    {
                        "path": str(source),
                        "repo_relative_path": "",
                        "absolute_path": str(source),
                        "role": "output_table",
                        "source_kind": "generated_table",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(source),
                        "checksum_unavailable_reason": "",
                        "byte_size": source.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "density-localization fixture",
                    }
                )
            write_tsv(run_root / "metadata/output_manifest.tsv", rows, MODULE_MANIFEST_COLUMNS)
            write_tsv(
                run_root / "metadata/run_config.tsv",
                [
                    {"key": "mode", "value": "standard"},
                    {"key": "panel_set", "value": "a-e"},
                    {"key": "tgi_day", "value": "17"},
                    {
                        "key": "config_sha256",
                        "value": sha256_file(
                            REPO_ROOT
                            / "Code/in-vivo/figure7/figure7_config.yaml"
                        ),
                    },
                    {
                        "key": "density_localization_config",
                        "value": (
                            "Code/in-vivo/figure7/"
                            "density_localization_config.yaml"
                        ),
                    },
                    {
                        "key": "density_localization_config_sha256",
                        "value": sha256_file(
                            REPO_ROOT
                            / "Code/in-vivo/figure7/"
                            "density_localization_config.yaml"
                        ),
                    },
                    {
                        "key": "cellcycle_input",
                        "value": str(cellcycle_path),
                    },
                    {
                        "key": "cellcycle_sha256",
                        "value": sha256_file(REPO_ROOT / cellcycle_path),
                    },
                    {
                        "key": "noncellcycle_input",
                        "value": str(noncellcycle_path),
                    },
                    {
                        "key": "noncellcycle_sha256",
                        "value": sha256_file(REPO_ROOT / noncellcycle_path),
                    },
                    {
                        "key": "endpoint_ploidy_input",
                        "value": str(endpoint_path),
                    },
                    {
                        "key": "endpoint_ploidy_sha256",
                        "value": sha256_file(REPO_ROOT / endpoint_path),
                    },
                    {
                        "key": "endpoint_ploidy_n_cells",
                        "value": "14125",
                    },
                    {
                        "key": "endpoint_ploidy_n_files",
                        "value": "16",
                    },
                    {
                        "key": "endpoint_ploidy_score_universe_n_cells",
                        "value": "9832",
                    },
                    {
                        "key": "endpoint_ploidy_treated_score_n_cells",
                        "value": "5335",
                    },
                    {
                        "key": "endpoint_ploidy_score_policy",
                        "value": (
                            "arithmetic_mean_of_finite_postprocessed_cell_"
                            "ploidy_in_exact_qc_passed_cellcycle_"
                            "noncellcycle_union_per_sample"
                        ),
                    },
                    {
                        "key": "endpoint_ploidy_mapping_policy",
                        "value": (
                            "exact_processed_sample_barcode_to_canonical_cbs_"
                            "file_cell_and_value;score_universe=reviewed_final_"
                            "seurat_tumor_cells"
                        ),
                    },
                ],
                ["key", "value"],
            )
            write_tsv(
                run_root / "metadata/panel_contract.tsv",
                [
                    {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                    for spec in PANEL_SPECS
                    if spec["module"] == "in_vivo_figure7"
                    and not str(spec["panel"]).startswith("7F")
                    and not str(spec["panel"]).startswith("7A-7L_composite")
                    and spec.get("variant", "pdf") == "pdf"
                ],
                ["panel_id", "filename"],
            )
            processed_input_paths = [
                cellcycle_path,
                noncellcycle_path,
                endpoint_path,
                Path("Code/in-vivo/figure7/run_figure7.R"),
                Path("Code/in-vivo/figure7/figure7_config.yaml"),
                Path("Code/in-vivo/figure7/density_localization_config.yaml"),
                Path("Code/in-vivo/figure7/src/common_io.R"),
                Path("Code/in-vivo/figure7/src/tgi_data.R"),
                Path("Code/in-vivo/figure7/src/tgi_statistics.R"),
                Path("Code/in-vivo/figure7/src/tgi_panels.R"),
            ]
            write_tsv(
                run_root / "metadata/input_manifest.tsv",
                [
                    {
                        "path": str(relative_path),
                        "repo_relative_path": str(relative_path),
                        "absolute_path": str(REPO_ROOT / relative_path),
                        "role": "input_data",
                        "source_kind": "input_file",
                        "module": "in_vivo_figure7",
                        "generated_by": "test",
                        "command_id": source_id,
                        "sha256": sha256_file(REPO_ROOT / relative_path),
                        "checksum_unavailable_reason": "",
                        "byte_size": (
                            REPO_ROOT / relative_path
                        ).stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "fixture",
                    }
                    for relative_path in processed_input_paths
                ],
                MODULE_MANIFEST_COLUMNS,
            )

            fake_bin = tmp_path / "bin"
            fake_bin.mkdir()
            sentinel = tmp_path / "rscript_was_called"
            fake_rscript = fake_bin / "Rscript"
            fake_rscript.write_text(f"#!/usr/bin/env bash\ntouch '{sentinel}'\nexit 99\n")
            fake_rscript.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env['PATH']}"

            result = self._run(
                "--mode", "panels-only",
                "--modules", "in_vivo_figure7",
                "--source-run-id", source_id,
                "--run-id", operation_id,
                "--output-root", str(output_root),
                "--figure-root", str(figure_root),
                env=env,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertFalse(sentinel.exists(), "panels-only unexpectedly invoked Rscript")

            with (figure_root / "Figure7/manifest.tsv").open(newline="") as handle:
                manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["run_id"] for row in manifest_rows}, {source_id})
            self.assertTrue(
                all(f"materialization_operation_id={operation_id}" in row["notes"] for row in manifest_rows)
            )
            module_runs = (
                output_root / "manager/runs" / operation_id / "metadata/module_runs.tsv"
            ).read_text()
            self.assertIn(f"source_run_id={source_id}", module_runs)
            self.assertIn(f"materialization_operation_id={operation_id}", module_runs)


if __name__ == "__main__":
    unittest.main()
