#!/usr/bin/env python3
"""Validate the paired Seurat metadata and scVelo inputs used by SI Figure 4."""

from __future__ import annotations

import argparse
import csv
import hashlib
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_provenance(path: Path) -> dict[str, str]:
    if not path.is_file() or path.stat().st_size <= 0:
        raise ValueError(f"Missing or empty Seurat metadata provenance: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or set(rows[0]) != {"key", "value"}:
        raise ValueError("Seurat metadata provenance must have key and value columns")
    keys = [row["key"].strip() for row in rows]
    if any(not key for key in keys) or len(keys) != len(set(keys)):
        raise ValueError("Seurat metadata provenance contains missing or duplicated keys")
    return {row["key"].strip(): row["value"].strip() for row in rows}


def validate(
    seurat_path: Path,
    metrics_path: Path,
    provenance_path: Path,
    config_path: Path,
) -> None:
    seurat_fields, seurat_rows = read_rows(seurat_path, "Seurat metadata")
    metrics_fields, metrics_rows = read_rows(metrics_path, "scVelo metrics")
    require_exact(seurat_fields, ("cell", "UMAP_1", "UMAP_2", "Dose"), "Seurat metadata")
    require_exact(metrics_fields, ("cell",), "scVelo metrics")
    require_exact(
        seurat_fields,
        (
            "clusters",
            "sample",
            "cluster_cell_cycle_annotation",
            "integrated_snn_res.0.6",
            "nCount_RNA",
            "nFeature_RNA",
            "percent.mt",
            "scDblFinder.class",
            "scDblFinder.score",
        ),
        "Seurat metadata",
    )
    cluster = "clusters"
    choose_column(metrics_fields, ("Ploidy", "ploidy"), "scVelo ploidy")
    choose_column(metrics_fields, ("TN", "trajectory_context", "trajectory_group", "sample_type"), "scVelo context")
    choose_column(metrics_fields, ("Dose", "Dose_DEG"), "scVelo dose")
    choose_column(metrics_fields, (cluster,), "scVelo cluster")
    observed_classifications = {
        row["cluster_cell_cycle_annotation"].strip() for row in seurat_rows
    }
    expected_classifications = {
        "cell_cycle_candidate",
        "not_cell_cycle_candidate",
    }
    if not observed_classifications or not observed_classifications <= expected_classifications:
        raise ValueError(
            "Seurat cluster_cell_cycle_annotation violates the reviewed classification contract"
        )

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
    provenance = read_provenance(provenance_path)
    required_provenance = {
        "source_seurat_rds",
        "source_seurat_rds_sha256",
        "source_object_cells",
        "seurat_metadata_sha256",
        "scvelo_metrics_sha256",
        "umap_reduction",
        "cluster_id_field",
        "base_cluster_field",
        "clustering_resolution",
        "cluster_annotation_field",
        "sample_field",
        "dose_field",
        "ploidy_field",
        "context_field",
        "cellcycle_mapping",
        "inclusion_context",
        "inclusion_ploidy_levels",
        "inclusion_dose_levels",
        "inclusion_required_nonmissing",
        "source_qc_fields",
        "source_qc_policy",
        "figure7_config_sha256",
        "export_script_sha256",
        "source_code_revision",
    }
    missing = sorted(required_provenance - provenance.keys())
    if missing:
        raise ValueError(
            "Seurat metadata provenance is missing required keys: " + ", ".join(missing)
        )
    if provenance["seurat_metadata_sha256"] != sha256_file(seurat_path):
        raise ValueError("Seurat metadata provenance checksum mismatch")
    if provenance["scvelo_metrics_sha256"] != sha256_file(metrics_path):
        raise ValueError("scVelo metrics provenance checksum mismatch")
    if provenance["figure7_config_sha256"] != sha256_file(config_path):
        raise ValueError("Figure 7 config provenance checksum mismatch")
    expected_contract = {
        "umap_reduction": "umap",
        "cluster_id_field": "clusters",
        "base_cluster_field": "integrated_snn_res.0.6",
        "clustering_resolution": "0.6",
        "cluster_annotation_field": "cluster_cell_cycle_annotation",
        "sample_field": "sample",
        "dose_field": "Dose",
        "ploidy_field": "Ploidy",
        "context_field": "TN",
        "cellcycle_mapping": (
            "cell_cycle_candidate->CellCycle;"
            "not_cell_cycle_candidate->NonCellCycle"
        ),
        "inclusion_context": "Tumor",
        "inclusion_ploidy_levels": "2N,4N",
        "inclusion_dose_levels": "0mg/kg,30mg/kg,120mg/kg",
        "inclusion_required_nonmissing": (
            "cell_id,UMAP_1,UMAP_2,sample_id,cluster_id,cluster_annotation,"
            "initial_ploidy,dose,cellcycle_classification"
        ),
        "source_qc_fields": (
            "nCount_RNA,nFeature_RNA,percent.mt,"
            "scDblFinder.class,scDblFinder.score"
        ),
        "source_qc_policy": (
            "Use the reviewed Seurat object as provided; SI Figure 4 applies "
            "no additional expression, mitochondrial, or doublet threshold."
        ),
    }
    if any(provenance[key] != value for key, value in expected_contract.items()):
        raise ValueError("Seurat metadata provenance does not match the reviewed SI Figure 4 contract")
    source_hash = provenance["source_seurat_rds_sha256"]
    if len(source_hash) != 64 or any(char not in "0123456789abcdef" for char in source_hash):
        raise ValueError("Source Seurat RDS provenance checksum is invalid")
    if int(provenance["source_object_cells"]) != len(seurat_rows):
        raise ValueError("Source-object cell count does not match Seurat metadata")
    for key, length in (("export_script_sha256", 64), ("source_code_revision", 40)):
        value = provenance[key]
        if len(value) != length or any(char not in "0123456789abcdef" for char in value):
            raise ValueError(f"Seurat metadata provenance has an invalid {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seurat-metadata", required=True, type=Path)
    parser.add_argument("--scvelo-metrics", required=True, type=Path)
    parser.add_argument("--seurat-metadata-provenance", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    validate(
        args.seurat_metadata.resolve(),
        args.scvelo_metrics.resolve(),
        args.seurat_metadata_provenance.resolve(),
        args.config.resolve(),
    )
    print("SI Figure 4 input bundle validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
