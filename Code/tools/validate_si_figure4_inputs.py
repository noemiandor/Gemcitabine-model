#!/usr/bin/env python3
"""Validate the paired Seurat metadata and scVelo inputs used by SI Figure 4."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def choose_column(fieldnames: list[str], candidates: tuple[str, ...], label: str) -> str:
    for candidate in candidates:
        if fieldnames.count(candidate) == 1:
            return candidate
    lower = [name.lower() for name in fieldnames]
    for candidate in candidates:
        matches = [fieldnames[i] for i, name in enumerate(lower) if name == candidate.lower()]
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(f"{label} has ambiguous matches for {candidate}")
    raise ValueError(f"{label} is missing one of: {', '.join(candidates)}")


def read_rows(path: Path, label: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"Missing or empty {label}: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if not fieldnames or len(fieldnames) != len(set(fieldnames)):
            raise ValueError(f"{label} has missing or duplicated column names: {path}")
        rows = list(reader)
    if not rows:
        raise ValueError(f"{label} has no data rows: {path}")
    return fieldnames, rows


def require_exact(fieldnames: list[str], required: tuple[str, ...], label: str) -> None:
    missing = [name for name in required if fieldnames.count(name) != 1]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def unique_cells(rows: list[dict[str, str]], label: str) -> set[str]:
    cells = [row.get("cell", "").strip() for row in rows]
    if any(not cell for cell in cells):
        raise ValueError(f"{label} contains an empty cell identifier")
    if len(cells) != len(set(cells)):
        raise ValueError(f"{label} contains duplicated cell identifiers")
    return set(cells)


def validate(seurat_path: Path, metrics_path: Path) -> None:
    seurat_fields, seurat_rows = read_rows(seurat_path, "Seurat metadata")
    metrics_fields, metrics_rows = read_rows(metrics_path, "scVelo metrics")
    require_exact(seurat_fields, ("cell", "UMAP_1", "UMAP_2", "Dose"), "Seurat metadata")
    require_exact(metrics_fields, ("cell",), "scVelo metrics")
    cluster = choose_column(seurat_fields, ("cluster_final", "clusters"), "Seurat metadata cluster")
    choose_column(seurat_fields, ("sample", "IDs", "ID", "orig.ident", "Sequencing.IDs"), "Seurat metadata sample")
    choose_column(metrics_fields, ("Ploidy", "ploidy"), "scVelo ploidy")
    choose_column(metrics_fields, ("TN", "trajectory_context", "trajectory_group", "sample_type"), "scVelo context")
    choose_column(metrics_fields, ("Dose", "Dose_DEG"), "scVelo dose")
    choose_column(metrics_fields, (cluster,), "scVelo cluster")

    seurat_cells = unique_cells(seurat_rows, "Seurat metadata")
    metrics_cells = unique_cells(metrics_rows, "scVelo metrics")
    missing = metrics_cells - seurat_cells
    if missing:
        examples = ", ".join(sorted(missing)[:10])
        raise ValueError(f"scVelo cells are missing from Seurat metadata: {examples}")
    for row in seurat_rows:
        for column in ("UMAP_1", "UMAP_2"):
            try:
                value = float(row[column])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"Seurat {column} contains a nonnumeric value") from exc
            if not math.isfinite(value):
                raise ValueError(f"Seurat {column} contains a nonfinite value")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seurat-metadata", required=True, type=Path)
    parser.add_argument("--scvelo-metrics", required=True, type=Path)
    args = parser.parse_args()
    validate(args.seurat_metadata.resolve(), args.scvelo_metrics.resolve())
    print("SI Figure 4 input pair validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
