from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from figure_output_contract import sha256_file
from materialize_figure7_data import (
    AUDIT_FILES,
    MANIFEST_COLUMNS,
    REFERENCE_FILES,
    materialize_generated_candidate_data,
)


class Figure7DataMaterializationTest(unittest.TestCase):
    def test_publishes_fixed_inputs_and_run_scoped_audited_bundle(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            repo = Path(temporary) / "repo"
            repo.mkdir()
            data_root = repo / "Data/in-vivo/figure7"
            source_root = repo / "Results/source"
            reference_id = "state_pathway_fixture_candidate"
            reference = source_root / "saved_state_pathway" / reference_id
            reference.mkdir(parents=True)
            for name in REFERENCE_FILES:
                (reference / name).write_text(f"reference:{name}\n")

            processed = source_root / "processed"
            processed.mkdir()
            cellcycle = processed / "CellCycleCells.csv"
            noncellcycle = processed / "NonCellCycleCells.csv"
            cellcycle.write_text("cell\nA\n")
            noncellcycle.write_text("cell\nB\n")
            generated_rds = source_root / "generated.rds"
            generated_rds.write_bytes(b"generated-rds-fixture")

            audit = source_root / "audit"
            audit.mkdir()
            for name in AUDIT_FILES:
                content = "status=PASS\n" if name == "AUDIT_COMPLETE.txt" else f"audit:{name}\n"
                (audit / name).write_text(content)

            run_id = "fixture_h5_candidate"
            run = repo / "Results/runs" / f"{run_id}_figure7"
            (run / "metadata").mkdir(parents=True)
            (run / "tables").mkdir()
            (run / "tables/panel.tsv").write_text("value\n1\n")
            config = {
                "canonical_publication_allowed": "false",
                "state_pathway_reference_id": reference_id,
                "workflow_seurat_rds": str(generated_rds),
                "workflow_state_pathway_reference": str(reference),
                "workflow_cellcycle_input": str(cellcycle),
                "workflow_cellcycle_sha256": sha256_file(cellcycle),
                "workflow_noncellcycle_input": str(noncellcycle),
                "workflow_noncellcycle_sha256": sha256_file(noncellcycle),
                "seurat_rds_sha256": sha256_file(generated_rds),
            }
            with (run / "metadata/run_config.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["key", "value"])
                writer.writerows(config.items())
            with (run / "metadata/input_manifest.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(["path"])
                writer.writerow([str(audit / "AUDIT_COMPLETE.txt")])

            bundle = materialize_generated_candidate_data(
                repo_root=repo,
                data_root=data_root,
                source_run_id=run_id,
                figure7_run=run,
                si_run=None,
                tgi_day=17,
                overwrite=False,
                validate_sources=False,
            )

            self.assertEqual(
                (data_root / "processed" / cellcycle.name).read_bytes(),
                cellcycle.read_bytes(),
            )
            for name in REFERENCE_FILES:
                self.assertEqual(
                    (data_root / "saved_state_pathway" / reference_id / name).read_bytes(),
                    (reference / name).read_bytes(),
                )
            self.assertEqual(
                (bundle / "seurat" / generated_rds.name).read_bytes(),
                generated_rds.read_bytes(),
            )
            self.assertTrue((bundle / "semantic_rds_audit/AUDIT_COMPLETE.txt").is_file())
            self.assertTrue((bundle / "figure7/tables/panel.tsv").is_file())
            self.assertIn("status=PASS", (bundle / "BUNDLE_COMPLETE.txt").read_text())
            with (bundle / "publication_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertGreater(len(rows), len(REFERENCE_FILES) + len(AUDIT_FILES))
            self.assertEqual(list(rows[0]), list(MANIFEST_COLUMNS))
            self.assertTrue(
                all(".tmp-" not in row["published_path"] for row in rows)
            )


if __name__ == "__main__":
    unittest.main()
