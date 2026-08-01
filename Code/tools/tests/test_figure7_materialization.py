from __future__ import annotations

import csv
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


TOOLS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = TOOLS_DIR.parents[1]
sys.path.insert(0, str(TOOLS_DIR))

from figure_output_contract import (  # noqa: E402
    MODULE_MANIFEST_COLUMNS,
    sha256_file,
    validate_expected_panel_set,
    validate_figure_manifest,
    validate_module_manifest,
    write_tsv,
)
from materialize_figure_assets import PANEL_SPECS  # noqa: E402


class Figure7MaterializationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.tmp = tempfile.TemporaryDirectory()
        self.repo = Path(self.tmp.name) / "repo"
        self.repo.mkdir()
        self.source_id = "analysis_run"
        self.operation_id = "materialize_run"
        self.run_root = (
            self.repo
            / "Results/in-vivo/figure7/runs"
            / f"{self.source_id}_figure7"
        )
        (self.run_root / "figures").mkdir(parents=True)
        (self.run_root / "metadata").mkdir()
        (self.run_root / "tables").mkdir()
        self.input_path = self.repo / "Data/in-vivo/figure7/fixture.tsv"
        self.input_path.parent.mkdir(parents=True)
        self.input_path.write_text("value\n1\n")
        self.reviewed_reference_source = (
            REPO_ROOT
            / "Data/in-vivo/figure7/saved_state_pathway"
            / "state_pathway_grch_human_only_initial_ploidy_day17_v3"
        )
        self.reviewed_reference_root = (
            self.repo
            / "Data/in-vivo/figure7/saved_state_pathway"
            / "state_pathway_grch_human_only_initial_ploidy_day17_v3"
        )
        shutil.copytree(
            self.reviewed_reference_source,
            self.reviewed_reference_root,
        )
        self.figure7_config = (
            self.repo / "Code/in-vivo/figure7/figure7_config.yaml"
        )
        self.figure7_config.parent.mkdir(parents=True)
        shutil.copy2(
            REPO_ROOT / "Code/in-vivo/figure7/figure7_config.yaml",
            self.figure7_config,
        )
        self.figure7_renderer = (
            self.repo / "Code/in-vivo/figure7/run_figure7.R"
        )
        shutil.copy2(
            REPO_ROOT / "Code/in-vivo/figure7/run_figure7.R",
            self.figure7_renderer,
        )
        self.processed_inputs = []
        for relative in (
            Path(
                "Data/in-vivo/figure7/processed/"
                "CellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            ),
            Path(
                "Data/in-vivo/figure7/processed/"
                "NonCellCycleCells_pseudotime_distribution_per_sample_"
                "cell_level_with_ploidy_dose_tgi.csv"
            ),
        ):
            destination = self.repo / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPO_ROOT / relative, destination)
            self.processed_inputs.append(destination)
        self.endpoint_ploidy = (
            self.repo / "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv"
        )
        self.endpoint_ploidy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(
            REPO_ROOT / "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv",
            self.endpoint_ploidy,
        )
        self.si_cache_root = self.repo / "Data/in-vivo/SIfigures"
        shutil.copytree(
            REPO_ROOT / "Data/in-vivo/SIfigures",
            self.si_cache_root,
        )
        self.context_inputs = []
        for relative in (
            "Code/in-vivo/figure7/src/context_panels.R",
            "Code/in-vivo/figure7/src/tgi_data.R",
            "Code/in-vivo/figure7/src/tgi_statistics.R",
            "Code/in-vivo/figure7/src/tgi_panels.R",
            "Code/in-vivo/SI_figures/shared_context_panels.R",
            "Code/in-vivo/SI_figures/normalized_composition.R",
            "Code/tools/validate_si_figures_table_cache.py",
        ):
            source = REPO_ROOT / relative
            destination = self.repo / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            self.context_inputs.append(destination)
        self.figure7_specs = [
            spec for spec in PANEL_SPECS if spec["module"] == "in_vivo_figure7"
        ]
        for spec in self.figure7_specs:
            source = self.run_root / str(spec["source"])
            source.write_bytes(f"fake PDF for {spec['panel']}\n".encode())
        self._write_run_metadata(include_f=True)
        self._write_state_provenance(canonical_publication_allowed="true")
        self._write_input_manifest()
        self._write_panel_k_tables()
        self._write_output_manifest()

    def tearDown(self) -> None:
        self.tmp.cleanup()

    def _write_run_metadata(self, include_f: bool) -> None:
        selected = [
            spec for spec in self.figure7_specs
            if spec.get("contract", spec.get("variant", "pdf") == "pdf")
            and (
                include_f
                or (
                    not str(spec["panel"]).startswith("7F")
                    and str(spec["panel"]) != "7A-7K_composite"
                )
            )
        ]
        run_config = [
            {"key": "panel_set", "value": "a-f" if include_f else "a-e"},
            {"key": "tgi_day", "value": "17"},
            {
                "key": "cellcycle_input",
                "value": str(self.processed_inputs[0].relative_to(self.repo)),
            },
            {
                "key": "cellcycle_sha256",
                "value": sha256_file(self.processed_inputs[0]),
            },
            {
                "key": "noncellcycle_input",
                "value": str(self.processed_inputs[1].relative_to(self.repo)),
            },
            {
                "key": "noncellcycle_sha256",
                "value": sha256_file(self.processed_inputs[1]),
            },
            {
                "key": "endpoint_ploidy_input",
                "value": "Data/in-vivo/scRNAseq_Numbat/all_ploidy.csv",
            },
            {
                "key": "endpoint_ploidy_sha256",
                "value": sha256_file(self.endpoint_ploidy),
            },
            {"key": "endpoint_ploidy_n_cells", "value": "14125"},
            {"key": "endpoint_ploidy_n_files", "value": "16"},
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
                    "arithmetic_mean_of_finite_postprocessed_cell_ploidy_"
                    "in_exact_qc_passed_cellcycle_noncellcycle_union_per_sample"
                ),
            },
            {
                "key": "endpoint_ploidy_mapping_policy",
                "value": (
                    "exact_processed_sample_barcode_to_canonical_cbs_file_"
                    "cell_and_value;score_universe=reviewed_final_seurat_tumor_cells"
                ),
            },
        ]
        if include_f:
            run_config.extend(
                [
                    {
                        "key": "state_pathway_reference_id",
                        "value": (
                            "state_pathway_grch_human_only_"
                            "initial_ploidy_day17_v3"
                        ),
                    },
                    {
                        "key": "state_pathway_reference_kind",
                        "value": (
                            "reviewed_human_only_initial_ploidy_frozen"
                        ),
                    },
                    {
                        "key": "canonical_publication_allowed",
                        "value": "true",
                    },
                    {
                        "key": "config_sha256",
                        "value": sha256_file(self.figure7_config),
                    },
                    {"key": "main_composite_panel_set", "value": "a-k"},
                    {
                        "key": "main_composite_filename",
                        "value": "Figure7_reviewed_GRCh.png",
                    },
                    {
                        "key": "main_composite_panel_order",
                        "value": (
                            "A=7A;B=7C;C=SI4A;D=SI4B;E=SI4C;F=SI4E;"
                            "G=SI7B;H=7B;I=7F;J=7D;K=7E"
                        ),
                    },
                    {
                        "key": "reviewed_si_cache_manifest",
                        "value": "Data/in-vivo/SIfigures/manifest.tsv",
                    },
                    {
                        "key": "reviewed_si_cache_manifest_sha256",
                        "value": sha256_file(
                            self.si_cache_root / "manifest.tsv"
                        ),
                    },
                    {
                        "key": "si_context_cache_policy",
                        "value": "reviewed",
                    },
                    {
                        "key": "si_context_cache_kind",
                        "value": "reviewed_human_only_frozen",
                    },
                    {
                        "key": "si_context_cache_manifest",
                        "value": "Data/in-vivo/SIfigures/manifest.tsv",
                    },
                    {
                        "key": "si_context_cache_manifest_sha256",
                        "value": sha256_file(
                            self.si_cache_root / "manifest.tsv"
                        ),
                    },
                    {
                        "key": (
                            "si_context_cache_canonical_"
                            "publication_allowed"
                        ),
                        "value": "true",
                    },
                ]
            )
        write_tsv(
            self.run_root / "metadata/run_config.tsv",
            run_config,
            ["key", "value"],
        )
        write_tsv(
            self.run_root / "metadata/panel_contract.tsv",
            [
                {"panel_id": spec["panel"], "filename": Path(str(spec["source"])).name}
                for spec in selected
            ],
            ["panel_id", "filename"],
        )

    def _write_input_manifest(self) -> None:
        input_paths = [
            self.input_path,
            self.figure7_config,
            self.figure7_renderer,
            self.endpoint_ploidy,
            *self.processed_inputs,
            *self.context_inputs,
            self.si_cache_root / "manifest.tsv",
            *sorted(
                path
                for path in self.si_cache_root.iterdir()
                if path.name != "manifest.tsv"
            ),
            *sorted(self.reviewed_reference_root.glob("*.tsv")),
        ]
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [
                {
                    "path": str(path.relative_to(self.repo)),
                    "repo_relative_path": str(path.relative_to(self.repo)),
                    "absolute_path": str(path),
                    "role": "input_data",
                    "source_kind": "input_file",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(path),
                    "checksum_unavailable_reason": "",
                    "byte_size": path.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "",
                    "panel": "",
                    "notes": "test fixture",
                }
                for path in input_paths
            ],
            MODULE_MANIFEST_COLUMNS,
        )

    def _write_state_provenance(
        self, canonical_publication_allowed: str
    ) -> None:
        if canonical_publication_allowed == "true":
            for source in self.reviewed_reference_root.glob("*.tsv"):
                destination = (
                    self.run_root / "metadata" / source.name
                    if source.name == "state_pathway_provenance.tsv"
                    else self.run_root / "tables" / source.name
                )
                shutil.copy2(source, destination)
            return
        write_tsv(
            self.run_root / "metadata/state_pathway_provenance.tsv",
            [
                {
                    "key": "reference_kind",
                    "value": "test_fixture",
                },
                {
                    "key": "canonical_publication_allowed",
                    "value": canonical_publication_allowed,
                },
            ],
            ["key", "value"],
        )

    def _write_panel_k_tables(self) -> None:
        plot_columns = [
            "sample_id",
            "initial_ploidy",
            "dose",
            "dose_mg",
            "endpoint_ploidy_file",
            "sample_mean_endpoint_ploidy",
            "n_endpoint_ploidy_cells",
            "endpoint_ploidy_source_total_cells",
            "endpoint_ploidy_score_universe_total_cells",
            "endpoint_ploidy_source_file_count",
            "endpoint_ploidy_source_sha256",
            "endpoint_ploidy_score_policy",
            "endpoint_ploidy_mapping_policy",
            "terminal_postprocessed_cn_score",
            "terminal_cn_score_within_origin_z",
            "permutation_stratum",
            "terminal_cn_score_nuisance_residual",
            "tgi_origin_dose_residual",
            "TGI_percent_Day_17",
            "tgi_outcome",
            "tgi_day",
            "tgi_measure",
            "matched_control_summary",
            "matched_control_group",
        ]
        design = [
            (
                "2N-A2-0", "2N", "30mg/kg", "30",
                "SUM159-2N-30-0_harvest.sps.cbs",
            ),
            (
                "2N-A2-L", "2N", "30mg/kg", "30",
                "SUM159-2N-30-L_harvest.sps.cbs",
            ),
            (
                "2N-A4-R", "2N", "120mg/kg", "120",
                "SUM159-2N-120-R_harvest.sps.cbs",
            ),
            (
                "2N-A4-RL", "2N", "120mg/kg", "120",
                "SUM159-2N-120-RL_harvest.sps.cbs",
            ),
            (
                "A6-4N-O", "4N", "30mg/kg", "30",
                "SUM159-4N-30-0_harvest.sps.cbs",
            ),
            (
                "A6-4N-RR", "4N", "30mg/kg", "30",
                "SUM159-4N-30-RR_harvest.sps.cbs",
            ),
            (
                "4N-A8-RL", "4N", "120mg/kg", "120",
                "SUM159-4N-120-RL_harvest.sps.cbs",
            ),
            (
                "4N-A8-RR", "4N", "120mg/kg", "120",
                "SUM159-4N-120-RR_harvest.sps.cbs",
            ),
        ]
        endpoint_scores: dict[str, list[float]] = {}
        endpoint_score_by_key: dict[tuple[str, str], float] = {}
        with self.endpoint_ploidy.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                score = float(row["ploidy"])
                endpoint_scores.setdefault(row["file"], []).append(score)
                endpoint_score_by_key[(row["file"], row["cell_id"])] = score
        curated_scores_by_sample: dict[str, list[float]] = {}
        for processed_path in self.processed_inputs:
            with processed_path.open(newline="") as handle:
                for row in csv.DictReader(handle):
                    sample_id = row["sample_id"]
                    prefix = f"{sample_id}_"
                    key = (
                        f"{row['growth_curve_harvest']}.sps.cbs",
                        row["cell_id"][len(prefix) :],
                    )
                    self.assertTrue(row["cell_id"].startswith(prefix))
                    self.assertAlmostEqual(
                        float(row["cell_ploidy"]),
                        endpoint_score_by_key[key],
                        places=12,
                    )
                    curated_scores_by_sample.setdefault(sample_id, []).append(
                        endpoint_score_by_key[key]
                    )
        self.assertEqual(sum(map(len, curated_scores_by_sample.values())), 9832)
        terminal_scores = [
            statistics.mean(curated_scores_by_sample[sample_id])
            for sample_id, *_ in design
        ]
        z_scores: list[float] = []
        for origin in ("2N", "4N"):
            indices = [
                index
                for index, row in enumerate(design)
                if row[1] == origin
            ]
            local_scores = [terminal_scores[index] for index in indices]
            local_mean = statistics.mean(local_scores)
            local_sd = statistics.stdev(local_scores)
            local_z = [
                (score - local_mean) / local_sd for score in local_scores
            ]
            z_scores.extend(local_z)
        tgi_by_sample: dict[str, float] = {}
        with self.processed_inputs[0].open(newline="") as handle:
            for row in csv.DictReader(handle):
                sample_id = row["sample_id"]
                if sample_id not in tgi_by_sample:
                    tgi_by_sample[sample_id] = float(
                        row["TGI_percent_Day_17"]
                    )
        plot_rows = []
        for index, (
            sample_id,
            origin,
            dose,
            dose_mg,
            endpoint_file,
        ) in enumerate(design):
            z_score = z_scores[index]
            terminal_score = terminal_scores[index]
            plot_rows.append(
                {
                    "sample_id": sample_id,
                    "initial_ploidy": origin,
                    "dose": dose,
                    "dose_mg": dose_mg,
                    "endpoint_ploidy_file": endpoint_file,
                    "sample_mean_endpoint_ploidy": terminal_score,
                    "n_endpoint_ploidy_cells": len(
                        curated_scores_by_sample[sample_id]
                    ),
                    "endpoint_ploidy_source_total_cells": "14125",
                    "endpoint_ploidy_score_universe_total_cells": "9832",
                    "endpoint_ploidy_source_file_count": "16",
                    "endpoint_ploidy_source_sha256": sha256_file(
                        self.endpoint_ploidy
                    ),
                    "endpoint_ploidy_score_policy": (
                        "arithmetic_mean_of_finite_postprocessed_cell_ploidy_"
                        "in_exact_qc_passed_cellcycle_noncellcycle_union_per_sample"
                    ),
                    "endpoint_ploidy_mapping_policy": (
                        "exact_processed_sample_barcode_to_canonical_cbs_file_"
                        "cell_and_value;score_universe=reviewed_final_seurat_tumor_cells"
                    ),
                    "terminal_postprocessed_cn_score": terminal_score,
                    "terminal_cn_score_within_origin_z": z_score,
                    "permutation_stratum": f"{origin}|{dose_mg}",
                    "terminal_cn_score_nuisance_residual": "pending",
                    "tgi_origin_dose_residual": "pending",
                    "TGI_percent_Day_17": tgi_by_sample[sample_id],
                    "tgi_outcome": "day",
                    "tgi_day": "17",
                    "tgi_measure": "TGI_percent_Day_17",
                    "matched_control_summary": "mean",
                    "matched_control_group": "initial_ploidy",
                }
            )

        def additive_residual(column: str) -> list[float]:
            values = [float(row[column]) for row in plot_rows]
            grand_mean = statistics.mean(values)
            origin_means = {
                origin: statistics.mean(
                    value
                    for value, row in zip(values, plot_rows)
                    if row["initial_ploidy"] == origin
                )
                for origin in ("2N", "4N")
            }
            dose_means = {
                dose_mg: statistics.mean(
                    value
                    for value, row in zip(values, plot_rows)
                    if row["dose_mg"] == dose_mg
                )
                for dose_mg in ("30", "120")
            }
            return [
                value
                - origin_means[row["initial_ploidy"]]
                - dose_means[row["dose_mg"]]
                + grand_mean
                for value, row in zip(values, plot_rows)
            ]

        x_residual = additive_residual("terminal_cn_score_within_origin_z")
        y_residual = additive_residual("TGI_percent_Day_17")
        for row, x_value, y_value in zip(
            plot_rows,
            x_residual,
            y_residual,
        ):
            row["terminal_cn_score_nuisance_residual"] = x_value
            row["tgi_origin_dose_residual"] = y_value
        write_tsv(
            self.run_root / "tables/panel_7E_plot_data.tsv",
            plot_rows,
            plot_columns,
        )
        test_row = {
            "n": "8",
            "estimate": "-0.3431355546551816",
            "partial_correlation": "-0.3431355546551816",
            "effect_per_within_origin_sd": "-8.792298734079608",
            "permutation_p_two_sided": "0.75",
            "n_permutations": "16",
            "permutation_mode": (
                "exact_TGI_label_enumeration_within_initial_ploidy_x_dose"
            ),
            "permutation_strata": "initial_ploidy:dose_mg",
            "score_variable": "sample_mean_qc_passed_curated_cbs_cell_ploidy",
            "score_source_sha256": sha256_file(self.endpoint_ploidy),
            "score_source_n_cells": "9832",
            "score_inventory_n_cells": "14125",
            "score_source_n_files": "16",
            "treated_score_n_cells": "5335",
            "score_aggregation_policy": (
                "arithmetic_mean_of_finite_postprocessed_cell_ploidy_"
                "in_exact_qc_passed_cellcycle_noncellcycle_union_per_sample"
            ),
            "sample_mapping_policy": (
                "exact_processed_sample_barcode_to_canonical_cbs_file_"
                "cell_and_value;score_universe=reviewed_final_seurat_tumor_cells"
            ),
            "score_standardization": "z_score_within_initial_ploidy",
            "adjustment_terms": "initial_ploidy+dose_mg",
            "outcome_variable": "TGI_percent_Day_17",
            "plot_x": "terminal_cn_score_nuisance_residual",
            "plot_y": "tgi_origin_dose_residual",
            "tgi_outcome": "day",
            "tgi_day": "17",
            "tgi_measure": "TGI_percent_Day_17",
            "matched_control_summary": "mean",
            "matched_control_group": "initial_ploidy",
        }
        write_tsv(
            self.run_root / "tables/panel_7E_test.tsv",
            [test_row],
            list(test_row),
        )

    def _write_output_manifest(self) -> None:
        rows = []
        for spec in self.figure7_specs:
            source = self.run_root / str(spec["source"])
            relative = str(source.relative_to(self.repo))
            rows.append(
                {
                    "path": relative,
                    "repo_relative_path": relative,
                    "absolute_path": str(source),
                    "role": "output_figure",
                    "source_kind": "generated_panel",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(source),
                    "checksum_unavailable_reason": "",
                    "byte_size": source.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "Figure7",
                    "panel": spec["panel"],
                    "notes": "test",
                }
            )
        for source in (
            self.run_root / "tables/panel_7E_plot_data.tsv",
            self.run_root / "tables/panel_7E_test.tsv",
        ):
            relative = str(source.relative_to(self.repo))
            rows.append(
                {
                    "path": relative,
                    "repo_relative_path": relative,
                    "absolute_path": str(source),
                    "role": "output_table",
                    "source_kind": "generated_table",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(source),
                    "checksum_unavailable_reason": "",
                    "byte_size": source.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "",
                    "panel": "",
                    "notes": "panel K test fixture",
                }
            )
        write_tsv(
            self.run_root / "metadata/output_manifest.tsv",
            rows,
            MODULE_MANIFEST_COLUMNS,
        )

    def _run_materializer(self) -> subprocess.CompletedProcess[str]:
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
                f"in_vivo_figure7={self.run_root}",
                "--overwrite",
            ],
            text=True,
            capture_output=True,
        )

    def test_materializes_exact_reviewed_human_only_initial_ploidy_v3(
        self,
    ) -> None:
        result = self._run_materializer()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(
            (
                self.repo
                / "figures/Figure7/"
                "panel_7F_pseudotime_state_pathway_activity.pdf"
            ).is_file()
        )

    def test_materialized_panel_7e_records_confound_safe_panel_7k_role(
        self,
    ) -> None:
        result = self._run_materializer()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (
            self.repo / "figures/Figure7/manifest.tsv"
        ).open(newline="") as handle:
            rows = {
                row["panel"]: row
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertEqual(
            rows["7E"]["caption_role"],
            (
                "Day-17 TGI versus within-origin standardized terminal "
                "postprocessed CN score, adjusted for origin and dose"
            ),
        )

        run_config_path = self.run_root / "metadata/run_config.tsv"
        with run_config_path.open(newline="") as handle:
            run_config = {
                row["key"]: row["value"]
                for row in csv.DictReader(handle, delimiter="\t")
            }
        self.assertTrue(
            run_config["main_composite_panel_order"].endswith("K=7E")
        )

    def test_rejects_missing_panel_k_source_code_binding(self) -> None:
        manifest = self.run_root / "metadata/input_manifest.tsv"
        with manifest.open(newline="") as handle:
            original_rows = list(csv.DictReader(handle, delimiter="\t"))
        for relative in (
            "Code/in-vivo/figure7/src/tgi_data.R",
            "Code/in-vivo/figure7/src/tgi_statistics.R",
            "Code/in-vivo/figure7/src/tgi_panels.R",
        ):
            with self.subTest(relative=relative):
                rows = [
                    row for row in original_rows
                    if row["repo_relative_path"] != relative
                ]
                write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
                result = self._run_materializer()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(
                    "does not bind the panel K statistics and plotting helpers",
                    result.stderr,
                )
                self.assertIn(relative, result.stderr)

    def test_rejects_missing_panel_k_output_manifest_binding(self) -> None:
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows = [
            row for row in rows
            if not row["repo_relative_path"].endswith("panel_7E_test.tsv")
        ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "output manifest must bind exactly one test table row",
            result.stderr,
        )

    def test_rejects_pooled_confounded_panel_k_method(self) -> None:
        table = self.run_root / "tables/panel_7E_test.tsv"
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            columns = list(rows[0])
        rows[0].update(
            {
                "n_permutations": "40320",
                "permutation_mode": "exact_TGI_label_enumeration",
                "permutation_strata": "",
                "score_standardization": "none",
                "adjustment_terms": "none",
                "plot_x": "sample_mean_endpoint_ploidy",
                "plot_y": "TGI_percent_Day_17",
            }
        )
        write_tsv(table, rows, columns)
        self._write_output_manifest()

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "does not use the reviewed within-origin standardization",
            result.stderr,
        )

    def test_rejects_nonreviewed_panel_k_numeric_result(self) -> None:
        table = self.run_root / "tables/panel_7E_test.tsv"
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            columns = list(rows[0])
        rows[0]["effect_per_within_origin_sd"] = "-8.7"
        write_tsv(table, rows, columns)
        self._write_output_manifest()

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "does not reproduce the exact reviewed Day-17 adjusted beta",
            result.stderr,
        )

    def test_rejects_panel_k_scores_not_standardized_within_origin(
        self,
    ) -> None:
        table = self.run_root / "tables/panel_7E_plot_data.tsv"
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            columns = list(rows[0])
        rows[0]["terminal_cn_score_within_origin_z"] = "0"
        write_tsv(table, rows, columns)
        self._write_output_manifest()

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "not the per-mouse mean over the exact QC-passed curated cells",
            result.stderr,
        )

    def test_rejects_missing_exact_processed_input_binding(self) -> None:
        manifest = self.run_root / "metadata/input_manifest.tsv"
        with manifest.open(newline="") as handle:
            original_rows = list(csv.DictReader(handle, delimiter="\t"))

        for processed_input in self.processed_inputs:
            with self.subTest(processed_input=processed_input.name):
                omitted = str(processed_input.relative_to(self.repo))
                rows = [
                    row for row in original_rows
                    if row["repo_relative_path"] != omitted
                ]
                write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)

                result = self._run_materializer()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn(
                    "does not bind the exact processed input",
                    result.stderr,
                )
                self.assertIn(omitted, result.stderr)
                self.assertFalse((self.repo / "figures").exists())

    def test_rejects_missing_complete_endpoint_ploidy_inventory_binding(self) -> None:
        manifest = self.run_root / "metadata/input_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        endpoint_relative = str(self.endpoint_ploidy.relative_to(self.repo))
        rows = [
            row
            for row in rows
            if row["repo_relative_path"] != endpoint_relative
        ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "input manifest must bind exactly one reviewed 14,125-cell",
            result.stderr,
        )

    def test_rejects_plot_subset_or_wrong_mouse_panel_k_score(self) -> None:
        table = self.run_root / "tables/panel_7E_plot_data.tsv"
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            columns = list(rows[0])
        rows[0]["sample_mean_endpoint_ploidy"] = str(
            float(rows[0]["sample_mean_endpoint_ploidy"]) + 0.01
        )
        rows[0]["terminal_postprocessed_cn_score"] = rows[0][
            "sample_mean_endpoint_ploidy"
        ]
        write_tsv(table, rows, columns)
        self._write_output_manifest()

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "not the per-mouse mean over the exact QC-passed curated cells",
            result.stderr,
        )

    def test_rejects_tampered_panel_k_adjusted_residual(self) -> None:
        table = self.run_root / "tables/panel_7E_plot_data.tsv"
        with table.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
            columns = list(rows[0])
        rows[0]["terminal_cn_score_nuisance_residual"] = str(
            float(rows[0]["terminal_cn_score_nuisance_residual"]) + 0.01
        )
        write_tsv(table, rows, columns)
        self._write_output_manifest()

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "plotted residuals do not reproduce the origin-and-dose nuisance adjustment",
            result.stderr,
        )

    def test_rejects_generated_si_context_identity_as_canonical(self) -> None:
        run_config_path = self.run_root / "metadata/run_config.tsv"
        with run_config_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        generated_values = {
            "si_context_cache_policy": "generated-human-only",
            "si_context_cache_kind": "generated_human_only_run_scoped",
            "si_context_cache_manifest": (
                "Results/in-vivo/SI_figures/intermediates/"
                "generated_human_only_fixture/tables/manifest.tsv"
            ),
            "si_context_cache_manifest_sha256": "a" * 64,
            "si_context_cache_canonical_publication_allowed": "false",
        }
        for row in rows:
            if row["key"] in generated_values:
                row["value"] = generated_values[row["key"]]
        write_tsv(run_config_path, rows, ["key", "value"])

        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "A-K composite/cache contract is invalid",
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())

    def test_rejects_historical_mixed_v1_spoofed_as_reviewed_v3(
        self,
    ) -> None:
        historical = (
            REPO_ROOT
            / "Data/in-vivo/figure7/saved_state_pathway"
            / "taoli_04i_etp2_24_day17_v1"
        )
        reviewed_filenames = {
            path.name for path in self.reviewed_reference_root.glob("*.tsv")
            if path.name != "state_pathway_provenance.tsv"
        }
        for target in self.run_root.joinpath("tables").glob("*.tsv"):
            if target.name in reviewed_filenames:
                target.unlink()
        (self.run_root / "metadata/state_pathway_provenance.tsv").unlink()
        for source in historical.glob("*.tsv"):
            destination = (
                self.run_root / "metadata" / source.name
                if source.name == "state_pathway_provenance.tsv"
                else self.run_root / "tables" / source.name
            )
            shutil.copy2(source, destination)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Canonical Figure 7 materialization is prohibited",
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())

    def test_rejects_explicit_historical_mixed_v1_for_af(self) -> None:
        run_config_path = self.run_root / "metadata/run_config.tsv"
        with run_config_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        replacements = {
            "state_pathway_reference_id": "taoli_04i_etp2_24_day17_v1",
            "state_pathway_reference_kind": "historical_mixed_frozen",
            "canonical_publication_allowed": "false",
        }
        for row in rows:
            if row["key"] in replacements:
                row["value"] = replacements[row["key"]]
        write_tsv(run_config_path, rows, ["key", "value"])
        self._write_state_provenance(
            canonical_publication_allowed="false"
        )
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Canonical Figure 7 materialization is prohibited",
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())

    def test_materializes_explicit_ae_run_without_optional_panel_f(self) -> None:
        full_only = [
            spec
            for spec in self.figure7_specs
            if str(spec["panel"]).startswith("7F")
            or str(spec["panel"]) == "7A-7K_composite"
        ]
        for spec in full_only:
            (self.run_root / str(spec["source"])).unlink()
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = [
                row for row in csv.DictReader(handle, delimiter="\t")
                if not row["panel"].startswith("7F")
                and row["panel"] != "7A-7K_composite"
            ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        self._write_run_metadata(include_f=False)

        result = self._run_materializer()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.repo / "figures/Figure7/manifest.tsv").open(newline="") as handle:
            materialized = list(csv.DictReader(handle, delimiter="\t"))
        expected_panels = [value for letter in "ABCDE" for value in (f"7{letter}", f"7{letter}_png")]
        self.assertEqual([row["panel"] for row in materialized], expected_panels)
        self.assertFalse((self.repo / "figures/Figure7/panel_7F_pseudotime_state_pathway_activity.pdf").exists())

    def test_rejects_unrecorded_ae_omission(self) -> None:
        full_only = [
            spec
            for spec in self.figure7_specs
            if str(spec["panel"]).startswith("7F")
            or str(spec["panel"]) == "7A-7K_composite"
        ]
        for spec in full_only:
            (self.run_root / str(spec["source"])).unlink()
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = [
                row for row in csv.DictReader(handle, delimiter="\t")
                if not row["panel"].startswith("7F")
                and row["panel"] != "7A-7K_composite"
            ]
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("explicitly record panel_set=a-e", result.stderr)

    def test_rejects_missing_or_unexpected_figure(self) -> None:
        (self.run_root / str(self.figure7_specs[0]["source"])).unlink()
        missing = self._run_materializer()
        self.assertNotEqual(missing.returncode, 0)

        (self.run_root / str(self.figure7_specs[0]["source"])).write_bytes(b"restored")
        self._write_output_manifest()
        (self.run_root / "tables").mkdir(exist_ok=True)
        (self.run_root / "tables/unexpected.png").write_bytes(b"unexpected")
        unexpected = self._run_materializer()
        self.assertNotEqual(unexpected.returncode, 0)
        self.assertIn("exact panel inventory", unexpected.stderr)

    def test_rejects_checksum_mismatch(self) -> None:
        source = self.run_root / str(self.figure7_specs[0]["source"])
        source.write_bytes(b"changed after manifest")
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("sha256 does not match", result.stderr)

    def test_rejects_noncanonical_generated_panel_f(self) -> None:
        self._write_state_provenance(
            canonical_publication_allowed="false"
        )
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Canonical Figure 7 materialization is prohibited",
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())

    def test_rejects_spoofed_legacy_canonical_identity(self) -> None:
        write_tsv(
            self.run_root / "metadata/state_pathway_provenance.tsv",
            [
                {
                    "key": "canonical_reference_id",
                    "value": "taoli_04i_etp2_24_day17_v1",
                },
                {
                    "key": "code_revision_04i",
                    "value": "spoof",
                },
            ],
            ["key", "value"],
        )
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "Canonical Figure 7 materialization is prohibited",
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())

    def test_rejects_duplicate_source_manifest_row(self) -> None:
        manifest = self.run_root / "metadata/output_manifest.tsv"
        with manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows.append(dict(rows[0]))
        write_tsv(manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("Expected one output-manifest row", result.stderr)

    def test_rejects_missing_or_invalid_input_manifest(self) -> None:
        input_manifest = self.run_root / "metadata/input_manifest.tsv"
        input_manifest.unlink()
        missing = self._run_materializer()
        self.assertNotEqual(missing.returncode, 0)
        self.assertIn("Missing source input manifest", missing.stderr)

        self._write_input_manifest()
        self.input_path.write_text("changed after manifest\n")
        invalid = self._run_materializer()
        self.assertNotEqual(invalid.returncode, 0)
        self.assertIn("Invalid source input manifest", invalid.stderr)

    def test_rejects_input_manifest_from_another_run(self) -> None:
        input_manifest = self.run_root / "metadata/input_manifest.tsv"
        with input_manifest.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        rows[0]["command_id"] = "different_source_run"
        write_tsv(input_manifest, rows, MODULE_MANIFEST_COLUMNS)
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn("input-manifest provenance", result.stderr)

    def test_rejects_canonical_identity_without_reviewed_input_binding(
        self,
    ) -> None:
        retained_paths = [
            self.input_path,
            self.endpoint_ploidy,
            *self.processed_inputs,
            self.repo / "Code/in-vivo/figure7/src/tgi_data.R",
            self.repo / "Code/in-vivo/figure7/src/tgi_statistics.R",
            self.repo / "Code/in-vivo/figure7/src/tgi_panels.R",
        ]
        write_tsv(
            self.run_root / "metadata/input_manifest.tsv",
            [
                {
                    "path": str(path.relative_to(self.repo)),
                    "repo_relative_path": str(path.relative_to(self.repo)),
                    "absolute_path": str(path),
                    "role": "input_data",
                    "source_kind": "input_file",
                    "module": "in_vivo_figure7",
                    "generated_by": "Manager.sh",
                    "command_id": self.source_id,
                    "sha256": sha256_file(path),
                    "checksum_unavailable_reason": "",
                    "byte_size": path.stat().st_size,
                    "mtime_utc": "2026-07-16T00:00:00+00:00",
                    "figure": "",
                    "panel": "",
                    "notes": "insufficient binding",
                }
                for path in retained_paths
            ],
            MODULE_MANIFEST_COLUMNS,
        )
        result = self._run_materializer()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            (
                "does not bind the reviewed renderer, config, and "
                "eight-file"
            ),
            result.stderr,
        )
        self.assertFalse((self.repo / "figures").exists())


class ManifestContractTest(unittest.TestCase):
    def test_expected_set_and_duplicate_checks(self) -> None:
        rows = [
            {"panel": "7A", "source_kind": "generated_panel"},
            {"panel": "7A", "source_kind": "external"},
        ]
        errors = validate_expected_panel_set(rows, {"7A", "7B"}, Path("manifest.tsv"))
        self.assertTrue(any("duplicate panel ID" in error for error in errors))
        self.assertTrue(any("missing expected generated panel" in error for error in errors))

    def test_repo_relative_path_precedes_stale_absolute_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp)
            data = repo / "Results/run/value.tsv"
            data.parent.mkdir(parents=True)
            data.write_text("value\n1\n")
            manifest = repo / "manifest.tsv"
            write_tsv(
                manifest,
                [
                    {
                        "path": "Results/run/value.tsv",
                        "repo_relative_path": "Results/run/value.tsv",
                        "absolute_path": "/old/checkout/Results/run/value.tsv",
                        "role": "output_table",
                        "source_kind": "generated_table",
                        "module": "test",
                        "generated_by": "test",
                        "command_id": "run",
                        "sha256": sha256_file(data),
                        "checksum_unavailable_reason": "",
                        "byte_size": data.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "test",
                    }
                ],
                MODULE_MANIFEST_COLUMNS,
            )
            self.assertEqual(validate_module_manifest(manifest, repo, repo / "Results/run"), [])

    def test_external_absolute_path_follows_missing_repo_relative_path(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "checkout"
            repo.mkdir()
            external = Path(tmp) / "external.tsv"
            external.write_text("value\n1\n")
            manifest = repo / "manifest.tsv"
            write_tsv(
                manifest,
                [
                    {
                        "path": "Data/missing.tsv",
                        "repo_relative_path": "Data/missing.tsv",
                        "absolute_path": str(external),
                        "role": "input_data",
                        "source_kind": "input_file",
                        "module": "test",
                        "generated_by": "test",
                        "command_id": "run",
                        "sha256": sha256_file(external),
                        "checksum_unavailable_reason": "",
                        "byte_size": external.stat().st_size,
                        "mtime_utc": "2026-07-16T00:00:00+00:00",
                        "figure": "",
                        "panel": "",
                        "notes": "external fallback",
                    }
                ],
                MODULE_MANIFEST_COLUMNS,
            )
            self.assertEqual(validate_module_manifest(manifest, repo), [])

    def test_existing_figure_materialization_remains_supported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            run_root = repo / "Results/gdsc/runs/source_gdsc"
            (run_root / "drug_count_collapsed/figures").mkdir(parents=True)
            (run_root / "ploidy_enrichment_panels_primary_secondary_ABC.png").write_bytes(b"one")
            (
                run_root
                / "drug_count_collapsed/figures/drug_count_collapsed_enrichment_shared_order.png"
            ).write_bytes(b"two")
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root",
                    str(repo),
                    "--figure-root",
                    str(repo / "figures"),
                    "--run-id",
                    "source",
                    "--module-run",
                    f"gdsc={run_root}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertTrue((repo / "figures/Figure1/manifest.tsv").is_file())

    def test_optional_lci_panel_remains_optional(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            repo = Path(tmp) / "repo"
            drug_run = repo / "Results/drug_response"
            lci_run = repo / "Results/lci_overlays"
            lci_run.mkdir(parents=True)
            for spec in PANEL_SPECS:
                if spec["module"] != "drug_response":
                    continue
                source = drug_run / str(spec["source"])
                source.parent.mkdir(parents=True, exist_ok=True)
                source.write_bytes(b"source")
            result = subprocess.run(
                [
                    sys.executable,
                    str(TOOLS_DIR / "materialize_figure_assets.py"),
                    "--repo-root", str(repo),
                    "--figure-root", str(repo / "figures"),
                    "--run-id", "source",
                    "--module-run", f"drug_response={drug_run}",
                    "--module-run", f"lci_overlays={lci_run}",
                    "--overwrite",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with (repo / "figures/Figure3/manifest.tsv").open(newline="") as handle:
                panels = {row["panel"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertNotIn("3I", panels)

    def test_existing_manifests_validate_in_a_moved_checkout(self) -> None:
        manifest_paths = sorted((REPO_ROOT / "figures").glob("Figure[1-6]/manifest.tsv"))
        manifest_paths.append(REPO_ROOT / "figures/Supplementary/manifest.tsv")
        manifest_paths.append(
            REPO_ROOT / "figures/Figure7_Supplement/si8_manifest.tsv"
        )
        with tempfile.TemporaryDirectory() as tmp:
            moved_repo = Path(tmp) / "moved-checkout"
            for manifest in manifest_paths:
                with manifest.open(newline="") as handle:
                    rows = list(csv.DictReader(handle, delimiter="\t"))
                for row in rows:
                    for key in (
                        "asset_path",
                        "source_file",
                        "result_run_dir",
                        "local_provenance_path",
                    ):
                        value = row.get(key, "").strip()
                        if not value or Path(value).is_absolute():
                            continue
                        source = REPO_ROOT / value
                        destination = moved_repo / value
                        if source.is_file():
                            destination.parent.mkdir(parents=True, exist_ok=True)
                            shutil.copy2(source, destination)
                        elif source.is_dir():
                            destination.mkdir(parents=True, exist_ok=True)
                moved_manifest = moved_repo / manifest.relative_to(REPO_ROOT)
                moved_manifest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(manifest, moved_manifest)
                self.assertEqual(
                    validate_figure_manifest(moved_manifest, moved_repo),
                    [],
                    f"failed after checkout move: {manifest}",
                )


if __name__ == "__main__":
    unittest.main()
