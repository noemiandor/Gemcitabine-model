#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import hashlib
import math
from pathlib import Path


CLUSTERS = ("0", "2", "4c", "5", "6", "8", "10", "13", "14")

COMMON_FILES = (
    "si_figures_cell_metadata.csv",
    "si_figures_cluster_key.tsv",
    "si_figure4_cluster_context_composition.csv",
    "si_figure4_cluster_initial_ploidy_composition.csv",
    "si_figure5_cluster_composition_by_mouse.csv",
    "si_figure5_cluster_dose_composition.csv",
    "si_figure5_cluster_initial_ploidy_composition.csv",
    "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv",
    "si_figure6_endpoint_ploidy_join_audit.csv",
    "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    "si_figure7_cluster_Hallmark_GSEA_all.csv",
    "si_figure7_cluster_Hallmark_ORA_all.csv",
    "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
    "si_figure7_hallmark_ORA_universe.csv",
)

CLUSTER_FILES = tuple(
    filename
    for cluster in CLUSTERS
    for filename in (
        f"si_figure7_cluster_{cluster}_top100_up_ORA_input.csv",
        f"si_figure7_cluster_{cluster}_vs_rest_DEG.csv",
    )
)

EXPECTED_FILES = tuple(sorted((*COMMON_FILES, *CLUSTER_FILES)))


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
        headers = list(reader.fieldnames or ())
        rows = list(reader)
    return headers, rows


def require_columns(
    path: Path,
    headers: list[str],
    required: tuple[str, ...],
    errors: list[str],
) -> None:
    missing = [column for column in required if column not in headers]
    if missing:
        errors.append(f"{path.name}: missing column(s): {', '.join(missing)}")


def parse_finite(value: str) -> bool:
    try:
        return math.isfinite(float(value))
    except (TypeError, ValueError):
        return False


def validate_cache(cache_dir: Path) -> list[str]:
    errors: list[str] = []
    if not cache_dir.is_dir():
        return [f"cache directory does not exist: {cache_dir}"]

    observed = {
        path.name
        for path in cache_dir.iterdir()
        if path.is_file() and path.name != "manifest.tsv"
    }
    expected = set(EXPECTED_FILES)
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
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

    cell_headers, cells = tables["si_figures_cell_metadata.csv"]
    require_columns(
        cache_dir / "si_figures_cell_metadata.csv",
        cell_headers,
        (
            "cell_id",
            "UMAP_1",
            "UMAP_2",
            "sample_id",
            "cluster_id",
            "cluster_annotation",
            "initial_ploidy",
            "dose",
            "s_phase_score",
            "endpoint_ploidy",
            "context",
            "included_in_si_figures",
        ),
        errors,
    )
    cell_ids = [row.get("cell_id", "").strip() for row in cells]
    if any(not value for value in cell_ids) or len(set(cell_ids)) != len(cell_ids):
        errors.append("si_figures_cell_metadata.csv: cell_id must be nonempty and unique")
    for column in ("UMAP_1", "UMAP_2", "s_phase_score"):
        if any(not parse_finite(row.get(column, "")) for row in cells):
            errors.append(f"si_figures_cell_metadata.csv: {column} must be finite")
    contexts = {row.get("context", "") for row in cells}
    if contexts != {"Tumor", "CellLine"}:
        errors.append(
            "si_figures_cell_metadata.csv: context must contain exactly Tumor and CellLine"
        )
    observed_clusters = {row.get("cluster_id", "") for row in cells}
    if not observed_clusters <= set(CLUSTERS):
        errors.append(
            "si_figures_cell_metadata.csv: unexpected cluster(s): "
            + ",".join(sorted(observed_clusters - set(CLUSTERS)))
        )

    key_headers, cluster_key = tables["si_figures_cluster_key.tsv"]
    require_columns(
        cache_dir / "si_figures_cluster_key.tsv",
        key_headers,
        (
            "cluster_id",
            "cluster_annotation",
            "cellcycle_classification",
            "cluster_order",
            "color",
        ),
        errors,
    )
    key_clusters = tuple(row.get("cluster_id", "") for row in cluster_key)
    if key_clusters != CLUSTERS:
        errors.append(
            "si_figures_cluster_key.tsv: expected ordered clusters "
            + ",".join(CLUSTERS)
            + "; observed "
            + ",".join(key_clusters)
        )

    composition_contracts = {
        "si_figure4_cluster_context_composition.csv": (
            "group_value",
            "fill_value",
            "n_cells",
            "total_cells",
            "proportion",
        ),
        "si_figure4_cluster_initial_ploidy_composition.csv": (
            "group_value",
            "fill_value",
            "n_cells",
            "total_cells",
            "proportion",
        ),
        "si_figure5_cluster_composition_by_mouse.csv": (
            "mouse",
            "initial_ploidy",
            "dose",
            "cluster_final",
            "n_cells",
            "total_cells",
            "proportion",
        ),
        "si_figure5_cluster_dose_composition.csv": (
            "group_value",
            "fill_value",
            "n_cells",
            "total_cells",
            "proportion",
        ),
        "si_figure5_cluster_initial_ploidy_composition.csv": (
            "group_value",
            "fill_value",
            "n_cells",
            "total_cells",
            "proportion",
        ),
        "si_figure5_mouse_weighted_composition_by_initial_ploidy_dose.csv": (
            "initial_ploidy",
            "dose",
            "cluster_final",
            "n_mice",
            "mean_proportion",
        ),
    }
    for filename, required in composition_contracts.items():
        headers, rows = tables[filename]
        require_columns(cache_dir / filename, headers, required, errors)
        proportion_field = (
            "mean_proportion"
            if "mean_proportion" in required
            else "proportion"
        )
        if any(
            not parse_finite(row.get(proportion_field, ""))
            or not 0 <= float(row[proportion_field]) <= 1
            for row in rows
        ):
            errors.append(f"{filename}: {proportion_field} must be within [0,1]")

    endpoint_headers, endpoint_rows = tables[
        "si_figure6_endpoint_ploidy_join_audit.csv"
    ]
    require_columns(
        cache_dir / "si_figure6_endpoint_ploidy_join_audit.csv",
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
    tumor_rows = [row for row in endpoint_rows if row.get("context") == "Tumor"]
    endpoint_by_cell = {row.get("cell", ""): row for row in endpoint_rows}
    if len(endpoint_by_cell) != len(endpoint_rows):
        errors.append(
            "si_figure6_endpoint_ploidy_join_audit.csv: cell values must be unique"
        )
    if set(endpoint_by_cell) != set(cell_ids):
        errors.append(
            "si_figure6_endpoint_ploidy_join_audit.csv: cell inventory must match canonical metadata"
        )
    else:
        for cell in cells:
            audit = endpoint_by_cell[cell["cell_id"]]
            if audit.get("context") != cell.get("context"):
                errors.append(
                    "si_figure6_endpoint_ploidy_join_audit.csv: context differs from canonical metadata"
                )
                break
    if not tumor_rows or any(
        row.get("matched", "").lower() != "true"
        or not parse_finite(row.get("endpoint_ploidy", ""))
        for row in tumor_rows
    ):
        errors.append(
            "si_figure6_endpoint_ploidy_join_audit.csv: every Tumor row must have a finite endpoint ploidy"
        )
    cellline_rows = [
        row for row in endpoint_rows if row.get("context") == "CellLine"
    ]
    if any(
        row.get("matched", "").lower() == "true"
        or row.get("endpoint_ploidy", "").strip()
        for row in cellline_rows
    ):
        errors.append(
            "si_figure6_endpoint_ploidy_join_audit.csv: CellLine rows must not have endpoint ploidy"
        )

    for cluster in CLUSTERS:
        deg_name = f"si_figure7_cluster_{cluster}_vs_rest_DEG.csv"
        deg_headers, deg_rows = tables[deg_name]
        require_columns(
            cache_dir / deg_name,
            deg_headers,
            ("gene", "gene_symbol", "cluster"),
            errors,
        )
        if any(row.get("cluster", "") != cluster for row in deg_rows):
            errors.append(f"{deg_name}: cluster column does not match {cluster}")
        genes = [row.get("gene", "") for row in deg_rows]
        if any(not gene for gene in genes) or len(set(genes)) != len(genes):
            errors.append(f"{deg_name}: gene values must be nonempty and unique")

        ora_name = f"si_figure7_cluster_{cluster}_top100_up_ORA_input.csv"
        ora_headers, ora_rows = tables[ora_name]
        require_columns(cache_dir / ora_name, ora_headers, ("gene_symbol",), errors)
        if len(ora_rows) > 100:
            errors.append(f"{ora_name}: expected no more than 100 marker rows")

    for filename in (
        "si_figure7_cluster_Hallmark_ORA_annotation_score_heatmap_top20_matrix.tsv",
        "si_figure7_cluster_Hallmark_GSEA_NES_heatmap_top20_matrix.tsv",
    ):
        headers, matrix_rows = tables[filename]
        if headers != ["pathway", *CLUSTERS]:
            errors.append(
                f"{filename}: expected pathway plus ordered clusters {','.join(CLUSTERS)}"
            )
        if len(matrix_rows) != 20:
            errors.append(f"{filename}: expected exactly 20 pathway rows")
        for row in matrix_rows:
            if not row.get("pathway", ""):
                errors.append(f"{filename}: pathway values must be nonempty")
                break
            if any(not parse_finite(row.get(cluster, "")) for cluster in CLUSTERS):
                errors.append(f"{filename}: matrix values must be finite")
                break
        pathways = [row.get("pathway", "") for row in matrix_rows]
        if len(pathways) != len(set(pathways)):
            errors.append(f"{filename}: pathway values must be unique")

    manifest_path = cache_dir / "manifest.tsv"
    if manifest_path.is_file():
        manifest_headers, manifest_rows = read_rows(manifest_path)
        require_columns(
            manifest_path,
            manifest_headers,
            ("filename", "bytes", "sha256", "source_run"),
            errors,
        )
        manifest_by_name = {
            row.get("filename", ""): row for row in manifest_rows
        }
        if (
            len(manifest_by_name) != len(manifest_rows)
            or set(manifest_by_name) != set(EXPECTED_FILES)
        ):
            errors.append(
                "manifest.tsv: filename inventory must match the exact 32-table cache"
            )
        else:
            for filename in EXPECTED_FILES:
                path = cache_dir / filename
                row = manifest_by_name[filename]
                if row.get("bytes") != str(path.stat().st_size):
                    errors.append(f"manifest.tsv: byte count mismatch for {filename}")
                if row.get("sha256") != sha256(path):
                    errors.append(f"manifest.tsv: SHA-256 mismatch for {filename}")

    return errors


def write_manifest(
    cache_dir: Path,
    output: Path,
    source_run: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("filename", "bytes", "sha256", "source_run"),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for filename in EXPECTED_FILES:
            path = cache_dir / filename
            writer.writerow(
                {
                    "filename": filename,
                    "bytes": path.stat().st_size,
                    "sha256": sha256(path),
                    "source_run": source_run,
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate the complete canonical SI Figures table cache."
    )
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--write-manifest", type=Path)
    parser.add_argument("--source-run", default="")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    cache_dir = args.cache_dir.resolve()
    errors = validate_cache(cache_dir)
    if errors:
        for error in errors:
            print(error)
        return 1
    if args.write_manifest is not None:
        write_manifest(
            cache_dir,
            args.write_manifest.resolve(),
            args.source_run,
        )
    print(
        "Validated SI Figures table cache: "
        f"files={len(EXPECTED_FILES)}; clusters={len(CLUSTERS)}; cache={cache_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
