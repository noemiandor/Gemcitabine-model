#!/usr/bin/env python3
"""Publish the reusable data contract for an audited Figure 7 H5 candidate."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from figure_output_contract import read_tsv, repo_root_from, sha256_file, validate_run_id
from materialize_figure_assets import (
    module_manifest_local_path,
    panel_specs_for_figure7_variant,
    validate_generated_candidate_source_run,
)


REFERENCE_FILES = (
    "state_pathway_interval_definition.tsv",
    "panel_7F_pathway_activity_plot_data.tsv",
    "panel_7F_selected_pathway_gsea.tsv",
    "panel_7F_leading_edge_genes.tsv",
    "state_pathway_gene_ranking_complete.tsv",
    "state_pathway_gsea_complete.tsv",
    "state_pathway_sample_bin_coverage.tsv",
    "state_pathway_design_qc.tsv",
    "state_pathway_provenance.tsv",
)
AUDIT_FILES = (
    "AUDIT_COMPLETE.txt",
    "audit_summary.tsv",
    "metadata_comparison.tsv",
    "cluster_comparison.tsv",
    "graph_comparison.tsv",
    "assay_numeric_comparison.tsv",
    "pca_comparison.tsv",
    "umap_displacement_summary.tsv",
    "umap_largest_displacements.tsv",
    "command_comparison.tsv",
    "rds_file_identity.tsv",
)
MANIFEST_COLUMNS = (
    "role",
    "source_path",
    "published_path",
    "sha256",
    "byte_size",
    "materialization",
)


def read_key_values(path: Path) -> dict[str, str]:
    headers, rows = read_tsv(path)
    if headers != ["key", "value"]:
        raise ValueError(f"Invalid key/value metadata schema: {path}")
    result: dict[str, str] = {}
    for row in rows:
        key = row.get("key", "")
        value = row.get("value", "")
        if not key or not value or key in result:
            raise ValueError(f"Invalid or duplicate key/value metadata row: {path}")
        result[key] = value
    return result


def resolve_path(value: str, repo_root: Path) -> Path:
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def portable_path(path: Path, repo_root: Path) -> str:
    try:
        return path.resolve().relative_to(repo_root.resolve()).as_posix()
    except ValueError:
        return str(path.resolve())


def require_file(path: Path, label: str) -> Path:
    if not path.is_file():
        raise FileNotFoundError(f"Missing {label}: {path}")
    return path


def exact_audit_dir(figure7_run: Path, repo_root: Path) -> Path:
    _, rows = read_tsv(figure7_run / "metadata" / "input_manifest.tsv")
    markers = []
    for row in rows:
        path = module_manifest_local_path(row, repo_root)
        if path is not None and path.name == "AUDIT_COMPLETE.txt":
            markers.append(path.resolve())
    if len(markers) != 1:
        raise ValueError(
            "Figure 7 input manifest must bind exactly one semantic audit marker; "
            f"found {len(markers)}"
        )
    marker = require_file(markers[0], "semantic RDS audit marker")
    marker_values = {}
    for line in marker.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            marker_values[key] = value
    if marker_values.get("status") != "PASS":
        raise ValueError(f"Semantic RDS audit did not pass: {marker}")
    return marker.parent


def iter_regular_files(root: Path) -> list[Path]:
    if not root.is_dir():
        raise FileNotFoundError(f"Missing publication source directory: {root}")
    return sorted(path for path in root.rglob("*") if path.is_file())


def copy_exact(source: Path, destination: Path, *, hardlink: bool = False) -> str:
    source = require_file(source, "publication source")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.resolve() == destination.resolve():
        return "existing-exact"
    temporary = destination.with_name(f".{destination.name}.tmp-{os.getpid()}")
    if temporary.exists():
        temporary.unlink()
    method = "copy"
    if hardlink:
        try:
            os.link(source, temporary)
            method = "hardlink"
        except OSError:
            shutil.copy2(source, temporary)
    else:
        shutil.copy2(source, temporary)
    if sha256_file(temporary) != sha256_file(source):
        temporary.unlink(missing_ok=True)
        raise OSError(f"Published file checksum mismatch: {destination}")
    os.replace(temporary, destination)
    return method


def materialize_generated_candidate_data(
    *,
    repo_root: Path,
    data_root: Path,
    source_run_id: str,
    figure7_run: Path,
    si_run: Path | None,
    tgi_day: int,
    overwrite: bool,
    validate_sources: bool = True,
) -> Path:
    validate_run_id(source_run_id)
    repo_root = repo_root.resolve()
    data_root = data_root.resolve()
    figure7_run = figure7_run.resolve()
    si_run = si_run.resolve() if si_run is not None else None
    specs = panel_specs_for_figure7_variant(
        tgi_day, "Figure7", "generated-candidate"
    )
    if validate_sources:
        validate_generated_candidate_source_run(
            "in_vivo_figure7",
            figure7_run,
            specs,
            repo_root,
            source_run_id,
            tgi_day,
        )
        if si_run is not None:
            validate_generated_candidate_source_run(
                "si_figures",
                si_run,
                specs,
                repo_root,
                source_run_id,
                tgi_day,
            )

    config = read_key_values(figure7_run / "metadata" / "run_config.tsv")
    if config.get("canonical_publication_allowed") != "false":
        raise ValueError("Data materializer accepts only generated candidates")
    reference_id = config.get("state_pathway_reference_id", "")
    if not re.fullmatch(r"[A-Za-z0-9._-]+_candidate", reference_id):
        raise ValueError(f"Invalid generated state-pathway reference id: {reference_id}")
    rds_source = resolve_path(config["workflow_seurat_rds"], repo_root)
    reference_source = resolve_path(
        config["workflow_state_pathway_reference"], repo_root
    )
    processed_sources = (
        resolve_path(config["workflow_cellcycle_input"], repo_root),
        resolve_path(config["workflow_noncellcycle_input"], repo_root),
    )
    for source in processed_sources:
        require_file(source, "processed Figure 7 input")
        expected = config[
            "workflow_cellcycle_sha256"
            if source == processed_sources[0]
            else "workflow_noncellcycle_sha256"
        ]
        if sha256_file(source) != expected:
            raise ValueError(f"Processed Figure 7 input hash mismatch: {source}")
    if sha256_file(require_file(rds_source, "generated Seurat RDS")) != config["seurat_rds_sha256"]:
        raise ValueError("Generated Seurat RDS hash disagrees with Figure 7 run metadata")

    bundle = data_root / "generated_candidates" / source_run_id
    if bundle.exists():
        if not overwrite:
            raise FileExistsError(f"Candidate data bundle exists: {bundle}")
        stale = bundle.with_name(
            f"{bundle.name}.stale.{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}"
        )
        os.replace(bundle, stale)
    bundle.parent.mkdir(parents=True, exist_ok=True)
    build_root = Path(
        tempfile.mkdtemp(prefix=f".{source_run_id}.tmp-", dir=bundle.parent)
    )
    rows: list[dict[str, str]] = []

    def publish(source: Path, destination: Path, role: str, *, hardlink: bool = False) -> None:
        method = copy_exact(source, destination, hardlink=hardlink)
        final_destination = destination
        if build_root == destination or build_root in destination.parents:
            final_destination = bundle / destination.relative_to(build_root)
        rows.append(
            {
                "role": role,
                "source_path": portable_path(source, repo_root),
                "published_path": portable_path(final_destination, repo_root),
                "sha256": sha256_file(destination),
                "byte_size": str(destination.stat().st_size),
                "materialization": method,
            }
        )

    try:
        for source in processed_sources:
            publish(source, data_root / "processed" / source.name, "processed_cell_table")
        for name in REFERENCE_FILES:
            publish(
                reference_source / name,
                data_root / "saved_state_pathway" / reference_id / name,
                "state_pathway_reference",
            )
        publish(
            rds_source,
            build_root / "seurat" / rds_source.name,
            "audited_generated_rds",
            hardlink=True,
        )
        audit_dir = exact_audit_dir(figure7_run, repo_root)
        for name in AUDIT_FILES:
            publish(
                audit_dir / name,
                build_root / "semantic_rds_audit" / name,
                "semantic_rds_audit",
            )
        for label, run_root in (("figure7", figure7_run), ("si_figures", si_run)):
            if run_root is None:
                continue
            for section in ("tables", "metadata"):
                source_root = run_root / section
                for source in iter_regular_files(source_root):
                    publish(
                        source,
                        build_root / label / section / source.relative_to(source_root),
                        f"{label}_{section}",
                    )
        manifest = build_root / "publication_manifest.tsv"
        with manifest.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=MANIFEST_COLUMNS, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda row: row["published_path"]))
        complete = build_root / "BUNDLE_COMPLETE.txt"
        complete.write_text(
            "schema_version=figure7_generated_candidate_data_v1\n"
            f"status=PASS\nsource_run_id={source_run_id}\n"
            f"files={len(rows)}\n"
            f"manifest_sha256={sha256_file(manifest)}\n"
        )
        os.replace(build_root, bundle)
    except Exception:
        shutil.rmtree(build_root, ignore_errors=True)
        raise
    print(f"Published Figure 7 reusable data bundle: {bundle}")
    print(f"Published files: {len(rows)}")
    return bundle


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--data-root", type=Path, default=Path("Data/in-vivo/figure7"))
    parser.add_argument("--source-run-id", required=True)
    parser.add_argument("--figure7-run-dir", type=Path, required=True)
    parser.add_argument("--si-run-dir", type=Path)
    parser.add_argument("--figure7-tgi-day", type=int, default=17)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    repo_root = args.repo_root.resolve() if args.repo_root else repo_root_from(Path.cwd())
    data_root = args.data_root
    if not data_root.is_absolute():
        data_root = repo_root / data_root
    materialize_generated_candidate_data(
        repo_root=repo_root,
        data_root=data_root,
        source_run_id=args.source_run_id,
        figure7_run=args.figure7_run_dir,
        si_run=args.si_run_dir,
        tgi_day=args.figure7_tgi_day,
        overwrite=args.overwrite,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
