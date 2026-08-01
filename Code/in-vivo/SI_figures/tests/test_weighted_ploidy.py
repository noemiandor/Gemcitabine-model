#!/usr/bin/env python3
"""Focused byte-contract tests for the reviewed SI6 ploidy derivation."""

from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest


REPO_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_PATH = REPO_ROOT / "Data" / "in-vivo" / "weighted_ploidy.py"
CBS_ROOT = REPO_ROOT / "Data" / "in-vivo" / "scRNAseq_Numbat"
MANIFEST_PATH = CBS_ROOT / "cbs_manifest.tsv"
CANONICAL_FULL = CBS_ROOT / "all_ploidy.csv"
CANONICAL_REDUCED = REPO_ROOT / "Data" / "in-vivo" / "all_ploidy.tsv"
FULL_SHA256 = "80f4e6b78e7b6d8b73030da4889ecb5c09ee97c9f83fb771aec4d3908511b569"
REDUCED_SHA256 = "6db48ee5f196b37b58aa71d0472dd3deb06aaacb4b637070af1b27d9425db2b3"

spec = importlib.util.spec_from_file_location("weighted_ploidy", SCRIPT_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"Could not load {SCRIPT_PATH}")
weighted_ploidy = importlib.util.module_from_spec(spec)
try:
    spec.loader.exec_module(weighted_ploidy)
    OPTIONAL_IMPORT_ERROR = None
except ModuleNotFoundError as error:
    if error.name not in {"numpy", "pandas"}:
        raise
    OPTIONAL_IMPORT_ERROR = error


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@unittest.skipIf(
    OPTIONAL_IMPORT_ERROR is not None,
    "weighted-ploidy byte tests require numpy and pandas",
)
class ReviewedPloidyDerivationTest(unittest.TestCase):
    def test_manifest_order_matches_canonical_first_occurrence(self) -> None:
        paths = weighted_ploidy.files_from_reviewed_manifest(str(MANIFEST_PATH))
        observed_order = [Path(path).name for path in paths]
        canonical = weighted_ploidy.pd.read_csv(
            CANONICAL_FULL, sep="\t", usecols=["file"]
        )
        canonical_order = list(dict.fromkeys(canonical["file"].tolist()))
        self.assertEqual(observed_order, canonical_order)
        self.assertEqual(len(observed_order), 16)

    def test_manifest_regenerates_full_table_byte_for_byte(self) -> None:
        with tempfile.TemporaryDirectory(prefix="si6-ploidy-full-") as temp_dir:
            output = Path(temp_dir) / "all_ploidy.csv"
            status = weighted_ploidy.main(
                [
                    "--manifest",
                    str(MANIFEST_PATH),
                    "--out",
                    str(output),
                    "--sep",
                    "tsv",
                    "--expected-sha256",
                    FULL_SHA256,
                ]
            )
            self.assertEqual(status, 0)
            self.assertEqual(sha256(output), FULL_SHA256)
            self.assertEqual(output.read_bytes(), CANONICAL_FULL.read_bytes())

    def test_manifest_regenerates_reduced_table_byte_for_byte(self) -> None:
        with tempfile.TemporaryDirectory(prefix="si6-ploidy-reduced-") as temp_dir:
            output = Path(temp_dir) / "all_ploidy.tsv"
            status = weighted_ploidy.main(
                [
                    "--manifest",
                    str(MANIFEST_PATH),
                    "--out",
                    str(output),
                    "--sep",
                    "tsv",
                    "--omit-total-chromosomes",
                    "--expected-sha256",
                    REDUCED_SHA256,
                ]
            )
            self.assertEqual(status, 0)
            self.assertEqual(sha256(output), REDUCED_SHA256)
            self.assertEqual(output.read_bytes(), CANONICAL_REDUCED.read_bytes())

    def test_manifest_checksum_tampering_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="si6-ploidy-manifest-") as temp_dir:
            fixture = Path(temp_dir)
            manifest = weighted_ploidy.pd.read_csv(
                MANIFEST_PATH, sep="\t", dtype=str, keep_default_na=False
            )
            for filename in manifest["filename"]:
                shutil.copy2(CBS_ROOT / filename, fixture / filename)
            manifest.loc[0, "sha256"] = "0" * 64
            fixture_manifest = fixture / "cbs_manifest.tsv"
            manifest.to_csv(fixture_manifest, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "checksum changed"):
                weighted_ploidy.files_from_reviewed_manifest(
                    str(fixture_manifest)
                )


if __name__ == "__main__":
    unittest.main()
