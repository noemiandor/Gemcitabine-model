#!/usr/bin/env python3
"""Publish a deterministic index of Figure 7 upstream-analysis documentation."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


DOI = "10.5281/zenodo.21463392"
RECORD_URL = "https://zenodo.org/records/21463392"
RECORD_TITLE = (
    "Resource limitation rewires chromosome instability and ploidy evolution "
    "across in vitro and in vivo cancer models"
)
PUBLICATION_DATE = "2026-07-21"
EXPECTED_COLUMNS = [
    "document_role",
    "filename",
    "size_bytes",
    "md5",
    "url",
    "description",
]
EXPECTED_DOCUMENTS = 11
ALLOWED_MODULES = {"si_figures", "in_vivo_figure7"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_source_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise FileNotFoundError(f"Missing upstream-document manifest: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(
                "Upstream-document manifest schema mismatch: "
                f"expected={EXPECTED_COLUMNS}; observed={reader.fieldnames}"
            )
        rows = list(reader)
    if len(rows) != EXPECTED_DOCUMENTS:
        raise ValueError(
            f"Expected {EXPECTED_DOCUMENTS} upstream documents; found {len(rows)}"
        )
    roles = [row["document_role"] for row in rows]
    filenames = [row["filename"] for row in rows]
    if len(roles) != len(set(roles)) or len(filenames) != len(set(filenames)):
        raise ValueError("Upstream-document roles and filenames must be unique")
    for row in rows:
        if (
            not row["description"].strip()
            or not row["filename"].strip()
            or int(row["size_bytes"]) <= 0
            or len(row["md5"]) != 32
            or any(char not in "0123456789abcdef" for char in row["md5"])
            or row["url"]
            != (
                "https://zenodo.org/api/records/21463392/files/"
                f"{row['filename']}/content"
            )
        ):
            raise ValueError(
                f"Invalid upstream-document manifest row: {row['document_role']}"
            )
    return rows


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def publish(module: str, source_manifest: Path, output_dir: Path) -> None:
    if module not in ALLOWED_MODULES:
        raise ValueError(f"Unsupported module for Figure 7 upstream provenance: {module}")
    rows = read_source_manifest(source_manifest)
    output_dir.mkdir(parents=True, exist_ok=True)
    existing = list(output_dir.iterdir())
    if existing:
        raise FileExistsError(
            f"Upstream-analysis publication directory is not empty: {output_dir}"
        )

    published_columns = [
        "doi",
        "record_url",
        "document_role",
        "filename",
        "size_bytes",
        "md5",
        "download_url",
        "availability",
        "description",
    ]
    published_rows = [
        {
            "doi": DOI,
            "record_url": RECORD_URL,
            "document_role": row["document_role"],
            "filename": row["filename"],
            "size_bytes": row["size_bytes"],
            "md5": row["md5"],
            "download_url": row["url"],
            "availability": "authoritative_copy_on_zenodo",
            "description": row["description"],
        }
        for row in rows
    ]
    document_manifest = output_dir / "zenodo_document_manifest.tsv"
    write_tsv(document_manifest, published_rows, published_columns)

    statement = (
        "The Seurat RDS and loom files consumed by this workflow are processed "
        "artifacts with upstream data-generation and analysis steps. The reviewed "
        "descriptions, parameters, software/runtime records, and checksums for "
        f"those preceding steps are archived at Zenodo DOI {DOI}."
    )
    provenance_rows = [
        {"key": "schema_version", "value": "1"},
        {"key": "artifact", "value": "figure7_upstream_analysis_documentation"},
        {"key": "module", "value": module},
        {"key": "zenodo_doi", "value": DOI},
        {"key": "zenodo_record_url", "value": RECORD_URL},
        {"key": "zenodo_record_title", "value": RECORD_TITLE},
        {"key": "zenodo_publication_date", "value": PUBLICATION_DATE},
        {"key": "documentation_scope", "value": "upstream_seurat_and_loom_generation"},
        {"key": "documentation_statement", "value": statement},
        {"key": "document_count", "value": str(len(rows))},
        {
            "key": "source_manifest",
            "value": "Code/in-vivo/figure7/zenodo_upstream_analysis_documents.tsv",
        },
        {"key": "source_manifest_sha256", "value": sha256_file(source_manifest)},
        {
            "key": "published_document_manifest_sha256",
            "value": sha256_file(document_manifest),
        },
        {
            "key": "publication_mode",
            "value": "offline_index_of_authoritative_zenodo_documents",
        },
    ]
    write_tsv(output_dir / "provenance.tsv", provenance_rows, ["key", "value"])

    lines = [
        "# Figure 7 upstream analysis documentation",
        "",
        statement,
        "",
        "This Manager run begins from the processed Seurat RDS and loom artifacts. "
        "It does not rerun Cell Ranger, Seurat integration/refinement, or velocyto "
        "loom generation. The indexed Zenodo documents describe those preceding "
        "steps and remain the authoritative copies.",
        "",
        f"- Zenodo DOI: [{DOI}](https://doi.org/{DOI})",
        f"- Record: [{RECORD_TITLE}]({RECORD_URL})",
        f"- Indexed documents: {len(rows)}",
        f"- Consuming module: `{module}`",
        "",
        "## Published index",
        "",
        "| Role | File | Description |",
        "|---|---|---|",
    ]
    lines.extend(
        f"| {row['document_role']} | [{row['filename']}]({row['url']}) | "
        f"{row['description']} |"
        for row in rows
    )
    lines.extend(
        [
            "",
            "The accompanying `zenodo_document_manifest.tsv` records the size and "
            "Zenodo MD5 for every indexed document. `provenance.tsv` records the "
            "pinned source-manifest SHA-256 and publication contract.",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", required=True, choices=sorted(ALLOWED_MODULES))
    parser.add_argument("--source-manifest", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    publish(
        args.module,
        args.source_manifest.resolve(),
        args.output_dir.resolve(),
    )
    print(f"Published Figure 7 upstream-analysis index: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
