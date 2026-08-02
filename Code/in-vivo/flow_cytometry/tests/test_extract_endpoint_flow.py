#!/usr/bin/env python3

from __future__ import annotations

import csv
from collections import Counter
import importlib.util
import io
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "extract_endpoint_flow.py"
REPO_ROOT = Path(__file__).resolve().parents[4]
SPEC = importlib.util.spec_from_file_location("extract_endpoint_flow", SCRIPT)
assert SPEC and SPEC.loader
extractor = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = extractor
SPEC.loader.exec_module(extractor)


def write_fcs(path: Path, total: int) -> None:
    delimiter = "|"
    values = {
        "$TOT": str(total),
        "$PAR": "2",
        "$BYTEORD": "4,3,2,1",
        "$DATATYPE": "F",
        "$FIL": path.name,
        "$CYT": "FACSCantoII",
        "$P1N": "450/50 Violet B-A",
        "$P2N": "FSC-A",
        "$DATE": "28-JAN-2025",
        "$BTIM": "12:00:00",
        "$ETIM": "12:01:00",
    }
    body = delimiter + delimiter.join(item for pair in values.items() for item in pair) + delimiter
    raw = body.encode("latin-1")
    text_start = 256
    text_end = text_start + len(raw) - 1
    header = (
        "FCS3.0"
        + " " * 4
        + f"{text_start:>8}"
        + f"{text_end:>8}"
        + f"{text_end + 1:>8}"
        + f"{text_end + 8:>8}"
        + f"{0:>8}"
        + f"{0:>8}"
    ).encode("ascii")
    assert len(header) == 58
    path.write_bytes(header + b" " * (text_start - len(header)) + raw + b"12345678")


def population(name: str, count: int) -> str:
    return f'<Population name="{name}" count="{count}" />'


def write_workspace(path: Path, samples: list[dict[str, object]]) -> None:
    chunks = ["<?xml version=\"1.0\"?><Workspace><SampleList>"]
    for index, sample in enumerate(samples, start=1):
        filename = str(sample["filename"])
        peak = str(sample["peak"])
        total = int(sample["total"])
        human = int(sample["human"])
        chunks.append(
            f'<Sample><DataSet uri="file:/old/machine/{filename}" sampleID="{index}" />'
            f'<SampleNode name="{filename}" count="{total}" sampleID="{index}">'
            f'<Subpopulations><Population name="Human Cell Enrichment" count="{total}">'
            f'<Subpopulations><Population name="HumanCells" count="{human}">'
            "<Subpopulations>"
            + population(peak, int(sample["peak_count"]))
            + population("2N", int(sample["two_n"]))
            + population("4N", int(sample["four_n"]))
            + population(r"mCh+\2N", int(sample["mch_two_n"]))
            + population(r"mCh+\4N", int(sample["mch_four_n"]))
            + "</Subpopulations></Population></Subpopulations></Population>"
            "</Subpopulations></SampleNode></Sample>"
        )
    chunks.append("</SampleList></Workspace>")
    path.write_text("".join(chunks), encoding="utf-8")


def write_crosswalk(path: Path, samples: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=extractor.REQUIRED_CROSSWALK_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for index, sample in enumerate(samples, start=1):
            writer.writerow(
                {
                    "mouse_id": f"mouse_{index}",
                    "injected_origin": "2N" if index == 1 else "4N",
                    "dose_mg_kg": "0" if index == 1 else "30",
                    "fcs_file": sample["filename"],
                    "wsp_sample_name": sample["filename"],
                    "peak_gate_name": sample["peak"],
                    "reviewed_peak_annotation_n": sample[
                        "reviewed_peak_annotation_n"
                    ],
                    "acquisition_date": "2025-01-28",
                }
            )


def write_expected_manifest(
    path: Path,
    crosswalk: Path,
    sensitivity: Path,
    workspace: Path,
    fcs_dir: Path,
) -> None:
    inputs = [
        ("crosswalk", crosswalk),
        ("selection_sensitivity", sensitivity),
        ("workspace", workspace),
        *(("fcs", item) for item in sorted(fcs_dir.glob("*.fcs"))),
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=extractor.EXPECTED_MANIFEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for role, input_path in inputs:
            writer.writerow(
                {
                    "input_role": role,
                    "relative_path": input_path.relative_to(path.parent).as_posix(),
                    "size_bytes": input_path.stat().st_size,
                    "sha256": extractor.sha256_file(input_path),
                }
            )


def write_sensitivity(path: Path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=extractor.SENSITIVITY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow(
            {
                "mouse_id": "mouse_2",
                "selected_wsp_sample_name": "sample_two.fcs",
                "selected_peak_gate_name": "2.20N",
                "selected_peak_annotation_n": "2.20",
                "paired_wsp_sample_name": "sample_two_M.fcs",
                "paired_peak_gate_name": "2.08N",
                "paired_peak_annotation_n": "2.08",
                "relationship_status": "unresolved",
            }
        )


SAMPLES = [
    {
        "filename": "sample_one.fcs",
        "peak": "2N 50",
        "reviewed_peak_annotation_n": "2.04",
        "total": 1000,
        "human": 800,
        "peak_count": 600,
        "two_n": 650,
        "four_n": 20,
        "mch_two_n": 300,
        "mch_four_n": 3,
    },
    {
        "filename": "sample_two.fcs",
        "peak": "2.20N",
        "reviewed_peak_annotation_n": "2.20",
        "total": 2000,
        "human": 1600,
        "peak_count": 1200,
        "two_n": 1300,
        "four_n": 40,
        "mch_two_n": 600,
        "mch_four_n": 12,
    },
]

PAIRED_SAMPLES = [
    {
        "filename": "sample_two_M.fcs",
        "peak": "2.08N",
        "reviewed_peak_annotation_n": "2.08",
        "total": 1800,
        "human": 1500,
        "peak_count": 1100,
        "two_n": 1250,
        "four_n": 30,
        "mch_two_n": 550,
        "mch_four_n": 10,
    }
]


class EndpointFlowExtractorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.fcs_dir = self.root / "fcs"
        self.fcs_dir.mkdir()
        for sample in SAMPLES:
            write_fcs(self.fcs_dir / str(sample["filename"]), int(sample["total"]))
        self.workspace = self.root / "reviewed.wsp"
        self.crosswalk = self.root / "crosswalk.tsv"
        self.sensitivity = self.root / "paired_acquisition_sensitivity.tsv"
        self.expected_manifest = self.root / "content_manifest.tsv"
        write_workspace(self.workspace, SAMPLES + PAIRED_SAMPLES)
        write_crosswalk(self.crosswalk, SAMPLES)
        write_sensitivity(self.sensitivity)
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def command(self, output: Path) -> list[str]:
        return [
            sys.executable,
            str(SCRIPT),
            "--crosswalk",
            str(self.crosswalk),
            "--paired-sensitivity",
            str(self.sensitivity),
            "--workspace",
            str(self.workspace),
            "--expected-input-manifest",
            str(self.expected_manifest),
            "--fcs-dir",
            str(self.fcs_dir),
            "--output-dir",
            str(output),
            "--expected-samples",
            "2",
        ]

    def test_fcs_metadata_reader(self) -> None:
        metadata = extractor.read_fcs_metadata(self.fcs_dir / "sample_one.fcs")
        self.assertEqual(metadata.version, "FCS3.0")
        self.assertEqual(metadata.keywords["$TOT"], "1000")
        self.assertEqual(metadata.keywords["$FIL"], "sample_one.fcs")

    def test_end_to_end_outputs_are_deterministic(self) -> None:
        first = self.root / "first"
        second = self.root / "second"
        result = subprocess.run(self.command(first), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        result = subprocess.run(self.command(second), text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stderr)
        for filename in (
            "endpoint_flow_per_sample.tsv",
            "endpoint_flow_group_summary.tsv",
            "input_hashes.tsv",
            "paired_workspace_sensitivity.tsv",
        ):
            self.assertEqual((first / filename).read_bytes(), (second / filename).read_bytes())

        with (first / "endpoint_flow_per_sample.tsv").open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["human_cells_pct_total"], "80.000000")
        self.assertEqual(rows[0]["human_cells_pct_enrichment"], "80.000000")
        self.assertEqual(rows[0]["human_cells_qc_pass"], "FALSE")
        self.assertEqual(
            rows[0]["mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts"],
            "0.990099",
        )
        self.assertEqual(rows[1]["reviewed_peak_annotation_n"], "2.200000")
        self.assertEqual(rows[0]["fcs_cytometer"], "FACSCantoII")
        self.assertEqual(rows[0]["dna_content_parameter_index"], "1")
        self.assertEqual(rows[0]["dna_content_parameter_name"], "450/50 Violet B-A")

        with (first / "endpoint_flow_group_summary.tsv").open(newline="") as handle:
            groups = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(groups[0]["group_level"], "all")
        self.assertEqual(groups[0]["reviewed_peak_gate_labels"], "2.20N;2N 50")
        self.assertEqual(groups[0]["n_samples_human_cells_qc_pass"], "1")
        self.assertEqual(len(groups), 5)

        with (first / "input_hashes.tsv").open(newline="") as handle:
            hashes = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(hashes), 5)
        self.assertTrue(all(len(row["sha256"]) == 64 for row in hashes))
        self.assertFalse(any(str(self.root) in row["logical_path"] for row in hashes))

    def test_expected_sample_count_is_enforced(self) -> None:
        command = self.command(self.root / "out")
        command[-1] = "16"
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("expected exactly 16", result.stderr)

    def test_duplicate_file_mapping_is_rejected(self) -> None:
        text = self.crosswalk.read_text(encoding="utf-8")
        text = text.replace("sample_two.fcs", "sample_one.fcs")
        self.crosswalk.write_text(text, encoding="utf-8")
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be unique", result.stderr)

    def test_parent_directory_traversal_is_rejected(self) -> None:
        text = self.crosswalk.read_text(encoding="utf-8")
        text = text.replace("sample_one.fcs", "../sample_one.fcs")
        self.crosswalk.write_text(text, encoding="utf-8")
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("cannot traverse outside", result.stderr)

    def test_missing_reviewed_gate_is_rejected(self) -> None:
        text = self.workspace.read_text(encoding="utf-8")
        text = text.replace('name="mCh+\\4N"', 'name="removed"', 1)
        self.workspace.write_text(text, encoding="utf-8")
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("missing direct HumanCells gate", result.stderr)

    def test_child_gate_count_above_human_cells_is_rejected(self) -> None:
        text = self.workspace.read_text(encoding="utf-8")
        text = text.replace('name="2N 50" count="600"', 'name="2N 50" count="900"', 1)
        self.workspace.write_text(text, encoding="utf-8")
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("count exceeds HumanCells count", result.stderr)

    def test_crosswalk_acquisition_date_must_match_fcs(self) -> None:
        text = self.crosswalk.read_text(encoding="utf-8")
        text = text.replace("2025-01-28", "2025-01-29")
        self.crosswalk.write_text(text, encoding="utf-8")
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("does not match FCS $DATE", result.stderr)

    def test_unselected_workspace_node_with_duplicate_gates_is_ignored(self) -> None:
        text = self.workspace.read_text(encoding="utf-8")
        unrelated = (
            '<Sample><DataSet uri="file:/old/control.fcs" sampleID="99" />'
            '<SampleNode name="control.fcs" count="10" sampleID="99">'
            '<Subpopulations><Population name="HumanCells" count="8">'
            '<Subpopulations><Population name="2N" count="4" />'
            '<Population name="2N" count="4" /></Subpopulations>'
            '</Population></Subpopulations></SampleNode></Sample>'
        )
        text = text.replace("</SampleList>", unrelated + "</SampleList>")
        self.workspace.write_text(text, encoding="utf-8")
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_fcs_workspace_total_mismatch_is_rejected(self) -> None:
        write_fcs(self.fcs_dir / "sample_one.fcs", 999)
        write_expected_manifest(
            self.expected_manifest,
            self.crosswalk,
            self.sensitivity,
            self.workspace,
            self.fcs_dir,
        )
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("does not equal workspace SampleNode count", result.stderr)

    def test_input_hash_mismatch_is_rejected(self) -> None:
        with (self.fcs_dir / "sample_one.fcs").open("ab") as handle:
            handle.write(b"tampered")
        result = subprocess.run(
            self.command(self.root / "out"), text=True, capture_output=True
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("does not match expected size/SHA-256", result.stderr)

    def test_unexpected_dna_channel_is_rejected(self) -> None:
        command = self.command(self.root / "out")
        command.extend(["--expected-dna-channel", "wrong channel"])
        result = subprocess.run(command, text=True, capture_output=True)
        self.assertEqual(result.returncode, 2)
        self.assertIn("expected exactly one FCS parameter", result.stderr)

    def test_text_parser_accepts_escaped_delimiter(self) -> None:
        parsed = extractor._parse_fcs_text(b"|KEY|a||b|NEXT|value|")
        self.assertEqual(parsed, {"KEY": "a|b", "NEXT": "value"})


class RealEndpointCohortContractTest(unittest.TestCase):
    DATASET = (
        REPO_ROOT
        / "Data"
        / "in-vivo"
        / "flow_cytometry"
        / "endpoint_tumors_20250128"
    )
    EXPECTED_ROWS = {
        "2N-A1-0": ("2N", "0", "Sample_A1-0-GENOMICS_035.fcs", "2N", 2.00),
        "2N-A1-LR": ("2N", "0", "Sample_A1-RL-GENOMICS_038.fcs", "2N", 2.00),
        "2N-A1-R": ("2N", "0", "Sample_A1-R_GENOMICS_036.fcs", "2N", 2.00),
        "2N-A1-RR": ("2N", "0", "Sample_A1-RR-GENOMICS_037.fcs", "2N", 2.00),
        "2N-A2-0": ("2N", "30", "Sample_A2-0_001.fcs", "2N 50", 2.00),
        "2N-A2-L": ("2N", "30", "Sample_A2-L_002.fcs", "2N 50", 2.00),
        "2N-A4-R": ("2N", "120", "Sample_A4-R_005.fcs", "2N 50", 2.00),
        "2N-A4-RL": ("2N", "120", "Sample_A4-RL_006.fcs", "2N 50", 2.00),
        "4N-A5-0": ("4N", "0", "Sample_A5-0_007.fcs", "1.88N", 1.88),
        "A5-4N-L": ("4N", "0", "Sample_A5-L_010.fcs", "2.08N", 2.08),
        "A5-4N-R": ("4N", "0", "Sample_A5-R_009.fcs", "2.2N", 2.20),
        "4N-A5-RR": ("4N", "0", "Sample_A5-RR_013.fcs", "1.88N", 1.88),
        "A6-4N-O": ("4N", "30", "Sample_A6-0_015.fcs", "2.12N", 2.12),
        "A6-4N-RR": ("4N", "30", "Sample_A6-RR_020.fcs", "2.12N", 2.12),
        "4N-A8-RL": ("4N", "120", "Sample_A8-RL_027.fcs", "2.04N", 2.04),
        "4N-A8-RR": ("4N", "120", "Sample_A8-RR_029.fcs", "2.08N", 2.08),
    }
    # total, HumanCells, peak, 2N, 4N, mCh+\2N, mCh+\4N
    EXPECTED_4N_COUNTS = {
        "4N-A5-0": (140295, 17969, 17519, 17684, 80, 10593, 142),
        "A5-4N-L": (177054, 14410, 14054, 14226, 11, 6612, 22),
        "A5-4N-R": (96925, 19997, 19606, 19430, 12, 10429, 84),
        "4N-A5-RR": (134017, 18658, 18402, 18534, 22, 10843, 62),
        "A6-4N-O": (218863, 18253, 18043, 18070, 19, 10439, 27),
        "A6-4N-RR": (176134, 14251, 13654, 14069, 24, 7949, 50),
        "4N-A8-RL": (153997, 18366, 17872, 18199, 71, 10447, 87),
        "4N-A8-RR": (192185, 18140, 17890, 18010, 22, 12108, 42),
    }

    def test_reviewed_cohort_identity_design_and_counts(self) -> None:
        crosswalk_path = self.DATASET / "crosswalk.tsv"
        sensitivity_path = self.DATASET / "paired_acquisition_sensitivity.tsv"
        workspace_path = self.DATASET / "workspace" / "20250129_TumorSamples.wsp"
        fcs_dir = self.DATASET / "fcs"
        expected_manifest_path = self.DATASET / "content_manifest.tsv"
        self.assertTrue(crosswalk_path.is_file(), "reviewed crosswalk is required")
        self.assertTrue(
            sensitivity_path.is_file(), "paired sensitivity mapping is required"
        )
        self.assertTrue(workspace_path.is_file(), "reviewed FlowJo workspace is required")
        self.assertTrue(
            expected_manifest_path.is_file(), "pinned input manifest is required"
        )

        rows = extractor.read_crosswalk(crosswalk_path, expected_samples=16)
        self.assertEqual({row.mouse_id for row in rows}, set(self.EXPECTED_ROWS))
        self.assertEqual(
            Counter((row.injected_origin, row.dose_mg_kg) for row in rows),
            Counter(
                {
                    ("2N", "0"): 4,
                    ("2N", "30"): 2,
                    ("2N", "120"): 2,
                    ("4N", "0"): 4,
                    ("4N", "30"): 2,
                    ("4N", "120"): 2,
                }
            ),
        )
        for row in rows:
            expected = self.EXPECTED_ROWS[row.mouse_id]
            observed = (
                row.injected_origin,
                row.dose_mg_kg,
                row.fcs_file,
                row.peak_gate_name,
                float(row.reviewed_peak_annotation_n),
            )
            self.assertEqual(observed, expected, row.mouse_id)
            self.assertEqual(row.wsp_sample_name, row.fcs_file, row.mouse_id)
            self.assertTrue((fcs_dir / row.fcs_file).is_file(), row.mouse_id)

        sensitivity = extractor.read_sensitivity_table(sensitivity_path)
        workspace = extractor.read_workspace(
            workspace_path,
            {row.wsp_sample_name for row in rows}
            | {row.paired_wsp_sample_name for row in sensitivity},
        )
        verified_hashes = extractor.verify_expected_input_manifest(
            expected_manifest_path,
            crosswalk_path,
            sensitivity_path,
            workspace_path,
            rows,
            fcs_dir,
        )
        self.assertEqual(len(verified_hashes), 19)
        extracted = extractor.extract_samples(
            rows,
            workspace,
            fcs_dir,
            expected_cytometer="FACSCantoII",
            expected_dna_channel="450/50 Violet B-A",
            min_human_cells=1000,
        )
        extracted_by_mouse = {str(row["mouse_id"]): row for row in extracted}
        self.assertEqual(len(extracted_by_mouse), 16)
        for mouse_id, expected in self.EXPECTED_4N_COUNTS.items():
            row = extracted_by_mouse[mouse_id]
            observed = (
                int(row["wsp_total_count"]),
                int(row["human_cells_count"]),
                int(row["peak_count"]),
                int(row["two_n_count"]),
                int(row["four_n_count"]),
                int(row["mch_pos_two_n_count"]),
                int(row["mch_pos_four_n_count"]),
            )
            self.assertEqual(observed, expected, mouse_id)

        four_n_rows = [
            extracted_by_mouse[mouse_id] for mouse_id in self.EXPECTED_4N_COUNTS
        ]
        peak_values = [
            float(row["reviewed_peak_annotation_n"]) for row in four_n_rows
        ]
        mch_four_n = [
            float(
                row[
                    "mch_pos_four_n_pct_of_summed_two_n_four_n_gate_counts"
                ]
            )
            for row in four_n_rows
        ]
        self.assertEqual(
            sorted(peak_values), [1.88, 1.88, 2.04, 2.08, 2.08, 2.12, 2.12, 2.2]
        )
        self.assertAlmostEqual(min(mch_four_n), 0.257978, places=6)
        self.assertAlmostEqual(max(mch_four_n), 1.322776, places=6)
        self.assertEqual(extracted_by_mouse["2N-A1-0"]["human_cells_qc_pass"], "FALSE")
        self.assertTrue(
            all(
                row["human_cells_qc_pass"] == "TRUE"
                for mouse_id, row in extracted_by_mouse.items()
                if mouse_id != "2N-A1-0"
            )
        )
        paired = extractor.extract_paired_sensitivity(sensitivity, rows, workspace)
        self.assertEqual(len(paired), 4)
        self.assertEqual(
            {
                str(row["mouse_id"]): (
                    int(row["paired_human_cells_count"]),
                    int(row["paired_peak_count"]),
                    str(row["paired_peak_gate_name"]),
                )
                for row in paired
            },
            {
                "4N-A5-0": (20516, 19258, "1.92N"),
                "4N-A5-RR": (19969, 19377, "1.88N"),
                "4N-A8-RL": (18263, 15163, "2.2N"),
                "4N-A8-RR": (18189, 14388, "2.2N"),
            },
        )


if __name__ == "__main__":
    unittest.main()
