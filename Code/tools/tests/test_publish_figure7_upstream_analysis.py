from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path
import sys


TOOLS_DIR = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[3]
SOURCE_MANIFEST = (
    REPO_ROOT
    / "Code/in-vivo/figure7/zenodo_upstream_analysis_documents.tsv"
)
sys.path.insert(0, str(TOOLS_DIR))

from publish_figure7_upstream_analysis import publish  # noqa: E402


class PublishFigure7UpstreamAnalysisTest(unittest.TestCase):
    def test_publishes_uniform_offline_index(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "upstream_analysis"
            publish("si_figures", SOURCE_MANIFEST, output)
            self.assertEqual(
                {path.name for path in output.iterdir()},
                {"README.md", "zenodo_document_manifest.tsv", "provenance.tsv"},
            )
            with (output / "zenodo_document_manifest.tsv").open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 11)
            self.assertEqual({row["doi"] for row in rows}, {"10.5281/zenodo.21463392"})
            self.assertEqual(
                {row["availability"] for row in rows},
                {"authoritative_copy_on_zenodo"},
            )
            with (output / "provenance.tsv").open(newline="") as handle:
                provenance = {
                    row["key"]: row["value"]
                    for row in csv.DictReader(handle, delimiter="\t")
                }
            self.assertEqual(provenance["module"], "si_figures")
            self.assertEqual(provenance["document_count"], "11")
            self.assertIn("upstream", provenance["documentation_statement"])

    def test_rejects_a_tampered_source_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = root / "documents.tsv"
            manifest.write_text(
                SOURCE_MANIFEST.read_text().replace(
                    "f90aedf2586b79b1cb98f8ac28a9880d",
                    "not-an-md5",
                )
            )
            with self.assertRaisesRegex(ValueError, "Invalid upstream-document"):
                publish("in_vivo_figure7", manifest, root / "out")


if __name__ == "__main__":
    unittest.main()
