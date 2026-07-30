#!/usr/bin/env python3
"""Validate the minimal frozen-table cache used to plot Supplementary Figures 4-7."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


CLUSTERS = ("0", "2", "4c", "5", "6", "8", "10", "13", "14")
EXPECTED_FILES = (
    "si_figure4_cluster_context_composition.csv",
    "si_figure4_cluster_initial_ploidy_composition.csv",
    "si_figure5_cluster_composition_by_mouse.csv",
    "si_figure5_cluster_dose_composition.csv",
    "si_figure5_cluster_initial_ploidy_composition.csv",
    "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
    "si_figure6_endpoint_ploidy_join_audit.csv",
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
    "si_figures_cell_metadata.csv",
    "si_figures_cluster_key.tsv",
)
MANIFEST_COLUMNS = ("filename", "bytes", "sha256", "source_revision", "notes")
REVIEWED_MANIFEST_SHA256 = (
    "5713379814b457d470753ec92a4e9155ecf776fe8f22eed8c8d66eb56881167d"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return list(reader.fieldnames or ()), list(reader)


def require_columns(
    filename: str,
    headers: list[str],
    required: tuple[str, ...],
    errors: list[str],
) -> None:
    missing = [column for column in required if column not in headers]
    if missing:
        errors.append(f"{filename}: missing column(s): {', '.join(missing)}")


def finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def truthy(value: str) -> bool:
    return value.strip().lower() in {"true", "t", "1"}


def validate_proportions(
    filename: str,
    rows: list[dict[str, str]],
    value_column: str,
    errors: list[str],
) -> None:
    if any(
        not finite(row.get(value_column, ""))
        or not 0 <= float(row[value_column]) <= 1
        for row in rows
    ):
        errors.append(f"{filename}: {value_column} must be finite and within [0,1]")


def validate_cache(
    cache_dir: Path,
    si7_policy: str = "corrected-human-only",
) -> list[str]:
    errors: list[str] = []
    if not cache_dir.is_dir():
        return [f"cache directory does not exist: {cache_dir}"]

    observed = {
        path.name
        for path in cache_dir.iterdir()
        if path.is_file() and path.name != "manifest.tsv"
    }
    missing = sorted(set(EXPECTED_FILES) - observed)
    unexpected = sorted(observed - set(EXPECTED_FILES))
    if missing:
        errors.append("missing expected file(s): " + ", ".join(missing))
    if unexpected:
        errors.append("unexpected cache file(s): " + ", ".join(unexpected))
    if errors:
        return errors

    tables: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
    for filename in EXPECTED_FILES:
        path = cache_dir / filename
        if path.stat().st_size <= 0:
            errors.append(f"{filename}: file is empty")
            continue
        headers, rows = read_rows(path)
        if not headers or not rows:
            errors.append(f"{filename}: header or data rows are missing")
        tables[filename] = (headers, rows)
    if errors:
        return errors

    cell_name = "si_figures_cell_metadata.csv"
    cell_headers, cells = tables[cell_name]
    require_columns(
        cell_name,
        cell_headers,
        (
            "cell_id",
            "UMAP_1",
            "UMAP_2",
            "sample_id",
            "cluster_id",
            "cluster_annotation",
            "cluster_order",
            "cluster_color",
            "initial_ploidy",
            "dose",
            "s_phase_score",
            "endpoint_ploidy",
            "endpoint_file",
            "endpoint_cell_id",
            "context",
            "included_in_si_figures",
        ),
        errors,
    )
    cell_ids = [row.get("cell_id", "").strip() for row in cells]
    if any(not cell_id for cell_id in cell_ids) or len(set(cell_ids)) != len(cell_ids):
        errors.append(f"{cell_name}: cell_id must be nonempty and unique")
    for column in ("UMAP_1", "UMAP_2", "s_phase_score"):
        if any(not finite(row.get(column, "")) for row in cells):
            errors.append(f"{cell_name}: {column} must be finite")
    if {row.get("context", "") for row in cells} != {"Tumor", "CellLine"}:
        errors.append(f"{cell_name}: context must contain exactly Tumor and CellLine")
    observed_clusters = {row.get("cluster_id", "") for row in cells}
    if observed_clusters != set(CLUSTERS):
        errors.append(
            f"{cell_name}: expected clusters {','.join(CLUSTERS)}; "
            f"observed {','.join(sorted(observed_clusters))}"
        )

    key_name = "si_figures_cluster_key.tsv"
    key_headers, cluster_key = tables[key_name]
    require_columns(
        key_name,
        key_headers,
        (
            "cluster_id",
            "cluster_annotation",
            "cellcycle_classification",
            "cluster_order",
            "color",
            "n_all_cells",
            "n_included_tumor_cells",
        ),
        errors,
    )
    if tuple(row.get("cluster_id", "") for row in cluster_key) != CLUSTERS:
        errors.append(f"{key_name}: cluster order must be {','.join(CLUSTERS)}")
    for row in cluster_key:
        cluster = row.get("cluster_id", "")
        expected_all = sum(cell.get("cluster_id") == cluster for cell in cells)
        expected_included = sum(
            cell.get("cluster_id") == cluster
            and cell.get("context") == "Tumor"
            and truthy(cell.get("included_in_si_figures", ""))
            for cell in cells
        )
        if row.get("n_all_cells") != str(expected_all):
            errors.append(f"{key_name}: n_all_cells mismatch for cluster {cluster}")
        if row.get("n_included_tumor_cells") != str(expected_included):
            errors.append(
                f"{key_name}: n_included_tumor_cells mismatch for cluster {cluster}"
            )

    composition_contracts = {
        "si_figure4_cluster_context_composition.csv": (
            ("group_value", "fill_value", "n_cells", "total_cells", "proportion"),
            "proportion",
        ),
        "si_figure4_cluster_initial_ploidy_composition.csv": (
            ("group_value", "fill_value", "n_cells", "total_cells", "proportion"),
            "proportion",
        ),
        "si_figure5_cluster_composition_by_mouse.csv": (
            (
                "mouse",
                "initial_ploidy",
                "dose",
                "cluster_final",
                "n_cells",
                "total_cells",
                "proportion",
            ),
            "proportion",
        ),
        "si_figure5_cluster_dose_composition.csv": (
            ("group_value", "fill_value", "n_cells", "total_cells", "proportion"),
            "proportion",
        ),
        "si_figure5_cluster_initial_ploidy_composition.csv": (
            ("group_value", "fill_value", "n_cells", "total_cells", "proportion"),
            "proportion",
        ),
        "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv": (
            ("initial_ploidy", "dose", "cluster_final", "n_mice", "mean_proportion"),
            "mean_proportion",
        ),
    }
    for filename, (required, proportion_column) in composition_contracts.items():
        headers, rows = tables[filename]
        require_columns(filename, headers, required, errors)
        validate_proportions(filename, rows, proportion_column, errors)

    endpoint_name = "si_figure6_endpoint_ploidy_join_audit.csv"
    endpoint_headers, endpoint_rows = tables[endpoint_name]
    require_columns(
        endpoint_name,
        endpoint_headers,
        (
            "cell",
            "context",
            "initial_ploidy",
            "endpoint_file",
            "endpoint_cell_id",
            "endpoint_ploidy",
            "matched",
        ),
        errors,
    )
    endpoint_by_cell = {row.get("cell", ""): row for row in endpoint_rows}
    if len(endpoint_by_cell) != len(endpoint_rows) or set(endpoint_by_cell) != set(
        cell_ids
    ):
        errors.append(f"{endpoint_name}: cell inventory must match canonical metadata")
    else:
        for cell in cells:
            endpoint = endpoint_by_cell[cell["cell_id"]]
            if endpoint.get("context") != cell.get("context"):
                errors.append(f"{endpoint_name}: context differs from canonical metadata")
                break
            is_tumor = cell.get("context") == "Tumor"
            if is_tumor != truthy(endpoint.get("matched", "")):
                errors.append(f"{endpoint_name}: matched flag violates the context contract")
                break
            if is_tumor and not finite(endpoint.get("endpoint_ploidy", "")):
                errors.append(f"{endpoint_name}: every Tumor row needs endpoint ploidy")
                break
            if not is_tumor and endpoint.get("endpoint_ploidy", "").strip():
                errors.append(f"{endpoint_name}: CellLine rows cannot have endpoint ploidy")
                break

    for filename in (
        "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
        "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    ):
        headers, matrix_rows = tables[filename]
        if headers != ["pathway", *CLUSTERS]:
            errors.append(f"{filename}: expected pathway plus {','.join(CLUSTERS)}")
        if len(matrix_rows) != 20:
            errors.append(f"{filename}: expected exactly 20 pathway rows")
        pathways = [row.get("pathway", "").strip() for row in matrix_rows]
        if any(not pathway for pathway in pathways) or len(set(pathways)) != len(pathways):
            errors.append(f"{filename}: pathways must be nonempty and unique")
        if any(
            not finite(row.get(cluster, ""))
            for row in matrix_rows
            for cluster in CLUSTERS
        ):
            errors.append(f"{filename}: matrix values must be finite")

    manifest_path = cache_dir / "manifest.tsv"
    if not manifest_path.is_file():
        errors.append("manifest.tsv: missing")
    else:
        if (
            si7_policy == "corrected-human-only"
            and sha256(manifest_path) != REVIEWED_MANIFEST_SHA256
        ):
            errors.append(
                "manifest.tsv: canonical cache is not the exact reviewed manifest"
            )
        manifest_headers, manifest_rows = read_rows(manifest_path)
        if tuple(manifest_headers) != MANIFEST_COLUMNS:
            errors.append(
                "manifest.tsv: expected columns " + ",".join(MANIFEST_COLUMNS)
            )
        if any(
            Path(row.get("filename", "")).name != row.get("filename", "")
            for row in manifest_rows
        ):
            errors.append("manifest.tsv: filenames must be portable basenames")
        by_name = {row.get("filename", ""): row for row in manifest_rows}
        if len(by_name) != len(manifest_rows) or set(by_name) != set(EXPECTED_FILES):
            errors.append("manifest.tsv: inventory must match the exact 11-table cache")
        else:
            for filename in EXPECTED_FILES:
                path = cache_dir / filename
                row = by_name[filename]
                if row.get("bytes") != str(path.stat().st_size):
                    errors.append(f"manifest.tsv: byte count mismatch for {filename}")
                if row.get("sha256") != sha256(path):
                    errors.append(f"manifest.tsv: SHA-256 mismatch for {filename}")
                if not row.get("source_revision", "").strip():
                    errors.append(f"manifest.tsv: source_revision missing for {filename}")
        matrix_rows = [
            row
            for row in manifest_rows
            if row.get("filename", "").startswith("si_figure7_cluster_Hallmark")
        ]
        if len(matrix_rows) == 2:
            if si7_policy == "corrected-human-only" and any(
                "human-only GRCh" not in row.get("notes", "")
                for row in matrix_rows
            ):
                errors.append(
                    "manifest.tsv: SI7 matrices must record human-only GRCh policy"
                )
            if si7_policy == "generated-human-only" and any(
                "generated human-only GRCh policy"
                not in row.get("notes", "")
                or "not approved for canonical publication"
                not in row.get("notes", "")
                or not row.get("source_revision", "").startswith(
                    "raw-generated-human-only@"
                )
                for row in matrix_rows
            ):
                errors.append(
                    "manifest.tsv: generated SI7 matrices must record the "
                    "human-only policy, generated lineage, and canonical-"
                    "publication prohibition"
                )
    return errors


def write_manifest(
    cache_dir: Path,
    output: Path,
    source_revision: str,
) -> None:
    if not source_revision:
        raise ValueError("--source-revision is required with --write-manifest")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=MANIFEST_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for filename in EXPECTED_FILES:
            path = cache_dir / filename
            is_si7 = filename.startswith("si_figure7_cluster_Hallmark")
            writer.writerow(
                {
                    "filename": filename,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "source_revision": source_revision,
                    "notes": (
                        "human-only GRCh feature policy; mouse-aligned GRCm39 "
                        "features excluded; MSigDB 2026.1.Hs"
                        if is_si7
                        else "frozen plot-facing table"
                    ),
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument(
        "--si7-policy",
        choices=(
            "corrected-human-only",
            "generated-human-only",
        ),
        default="corrected-human-only",
        help=(
            "Expected SI7 feature policy. corrected-human-only is the exact "
            "reviewed Data cache; generated-human-only is unreviewed and "
            "noncanonical."
        ),
    )
    parser.add_argument("--write-manifest", type=Path)
    parser.add_argument("--source-revision", default="")
    args = parser.parse_args()

    cache_dir = args.cache_dir.resolve()
    if args.write_manifest is not None:
        write_manifest(
            cache_dir,
            args.write_manifest.resolve(),
            args.source_revision,
        )
    errors = validate_cache(cache_dir, args.si7_policy)
    if errors:
        for error in errors:
            print(error)
        return 1
    print(
        "Validated SI Figures table cache: "
        f"files={len(EXPECTED_FILES)}; clusters={len(CLUSTERS)}; cache={cache_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
