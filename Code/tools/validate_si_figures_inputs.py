#!/usr/bin/env python3
"""Validate the input bundle used by Supplementary Figures 4-7."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(
    path: Path, label: str, delimiter: str = ","
) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"Missing or empty {label}: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"{label} has missing or duplicated column names")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{label} has no data rows")
    return fields, rows


def require_fields(fields: list[str], required: tuple[str, ...], label: str) -> None:
    missing = [name for name in required if fields.count(name) != 1]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def unique_values(rows: list[dict[str, str]], field: str, label: str) -> set[str]:
    values = [row.get(field, "").strip() for row in rows]
    if any(not value for value in values):
        raise ValueError(f"{label} contains an empty {field}")
    if len(values) != len(set(values)):
        raise ValueError(f"{label} contains duplicated {field} values")
    return set(values)


def finite_number(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} contains a nonnumeric value") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} contains a nonfinite value")
    return result


def read_provenance(path: Path) -> dict[str, str]:
    fields, rows = read_rows(path, "Seurat metadata provenance", "\t")
    if fields != ["key", "value"]:
        raise ValueError("Seurat metadata provenance must have key and value columns")
    keys = [row["key"].strip() for row in rows]
    if any(not key for key in keys) or len(keys) != len(set(keys)):
        raise ValueError("Seurat metadata provenance contains missing or duplicated keys")
    return {row["key"].strip(): row["value"].strip() for row in rows}


def validate(args: argparse.Namespace) -> None:
    seurat_fields, seurat_rows = read_rows(args.seurat_metadata, "Seurat metadata")
    metrics_fields, metrics_rows = read_rows(args.scvelo_metrics, "scVelo metrics")
    ploidy_fields, ploidy_rows = read_rows(args.all_ploidy, "endpoint ploidy", "\t")

    require_fields(
        seurat_fields,
        (
            "cell",
            "UMAP_1",
            "UMAP_2",
            "S.Score",
            "Dose",
            "clusters",
            "sample",
            "Ploidy",
            "TN",
            "harvest",
            "barcode_raw",
            "cluster_cell_cycle_annotation",
            "integrated_snn_res.0.6",
        ),
        "Seurat metadata",
    )
    require_fields(metrics_fields, ("cell", "clusters", "Ploidy", "TN"), "scVelo metrics")
    require_fields(ploidy_fields, ("file", "cell_id", "ploidy"), "endpoint ploidy")

    seurat_cells = unique_values(seurat_rows, "cell", "Seurat metadata")
    metrics_cells = unique_values(metrics_rows, "cell", "scVelo metrics")
    missing_metrics = metrics_cells - seurat_cells
    if missing_metrics:
        raise ValueError(
            "scVelo cells are absent from Seurat metadata: "
            + ", ".join(sorted(missing_metrics)[:10])
        )

    metrics_by_cell = {row["cell"].strip(): row for row in metrics_rows}
    for row in seurat_rows:
        finite_number(row["UMAP_1"], "Seurat UMAP_1")
        finite_number(row["UMAP_2"], "Seurat UMAP_2")
        finite_number(row["S.Score"], "Seurat S.Score")
        metric = metrics_by_cell.get(row["cell"].strip())
        if metric is not None:
            if row["clusters"].strip() != metric["clusters"].strip():
                raise ValueError("Seurat/scVelo cluster mismatch")
            if row["Ploidy"].strip() != metric["Ploidy"].strip():
                raise ValueError("Seurat/scVelo initial-ploidy mismatch")
            if row["TN"].strip() != metric["TN"].strip():
                raise ValueError("Seurat/scVelo context mismatch")

    endpoint: dict[tuple[str, str], float] = {}
    for row in ploidy_rows:
        key = (row["file"].strip(), row["cell_id"].strip())
        if not all(key):
            raise ValueError("Endpoint ploidy contains an empty file/cell_id key")
        if key in endpoint:
            raise ValueError(f"Endpoint ploidy contains a duplicated key: {key}")
        endpoint[key] = finite_number(row["ploidy"], "endpoint ploidy")

    tumor_rows = [row for row in seurat_rows if row["TN"].strip() == "Tumor"]
    cellline_rows = [row for row in seurat_rows if row["TN"].strip() == "CellLine"]
    if len(tumor_rows) + len(cellline_rows) != len(seurat_rows):
        raise ValueError("Seurat TN must contain only Tumor and CellLine")

    def endpoint_key(row: dict[str, str]) -> tuple[str, str]:
        return (f"{row['harvest'].strip()}.sps.cbs", row["barcode_raw"].strip())

    missing_tumor = [row["cell"] for row in tumor_rows if endpoint_key(row) not in endpoint]
    matched_cellline = [row["cell"] for row in cellline_rows if endpoint_key(row) in endpoint]
    if missing_tumor:
        raise ValueError(
            "Endpoint ploidy is missing tumor cells: " + ", ".join(missing_tumor[:10])
        )
    if matched_cellline:
        raise ValueError(
            "Endpoint ploidy unexpectedly matches CellLine cells: "
            + ", ".join(matched_cellline[:10])
        )

    if not args.config.is_file() or args.config.stat().st_size <= 0:
        raise ValueError(f"Missing Figure 7 config: {args.config}")
    if args.seurat_rds is not None:
        if not args.seurat_rds.is_file() or args.seurat_rds.stat().st_size <= 0:
            raise ValueError(f"Missing raw Seurat RDS: {args.seurat_rds}")

    if args.seurat_metadata_provenance is not None:
        provenance = read_provenance(args.seurat_metadata_provenance)
        required = {
            "seurat_metadata_sha256",
            "scvelo_metrics_sha256",
            "figure7_config_sha256",
            "source_seurat_rds_sha256",
        }
        missing = required - provenance.keys()
        if missing:
            raise ValueError(
                "Seurat metadata provenance is missing keys: " + ", ".join(sorted(missing))
            )
        expected = {
            "seurat_metadata_sha256": sha256_file(args.seurat_metadata),
            "scvelo_metrics_sha256": sha256_file(args.scvelo_metrics),
            "figure7_config_sha256": sha256_file(args.config),
        }
        for key, value in expected.items():
            if provenance[key] != value:
                raise ValueError(f"Seurat metadata provenance checksum mismatch: {key}")
        if args.seurat_rds is not None:
            if provenance["source_seurat_rds_sha256"] != sha256_file(args.seurat_rds):
                raise ValueError("Raw Seurat RDS checksum differs from provenance")

    endpoint_values = [endpoint[endpoint_key(row)] for row in tumor_rows]
    print(
        "Validated SI Figures input bundle: "
        f"seurat_cells={len(seurat_rows)}; scvelo_cells={len(metrics_rows)}; "
        f"tumor_cells={len(tumor_rows)}; cellline_cells={len(cellline_rows)}; "
        f"endpoint_matches={len(endpoint_values)}; "
        f"endpoint_range={min(endpoint_values):.15g},{max(endpoint_values):.15g}"
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seurat-metadata", required=True, type=Path)
    parser.add_argument("--scvelo-metrics", required=True, type=Path)
    parser.add_argument("--all-ploidy", required=True, type=Path)
    parser.add_argument("--seurat-rds", type=Path)
    parser.add_argument("--seurat-metadata-provenance", type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    validate(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
