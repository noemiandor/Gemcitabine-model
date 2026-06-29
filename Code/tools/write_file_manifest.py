#!/usr/bin/env python3
"""Write a standard input or output manifest for local files."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from figure_output_contract import (
    MODULE_MANIFEST_COLUMNS,
    path_within,
    repo_root_from,
    sha256_file,
    write_tsv,
)


FIGURE_SUFFIXES = {".png", ".pdf", ".jpg", ".jpeg", ".tif", ".tiff", ".svg"}
TABLE_SUFFIXES = {".tsv", ".csv", ".xlsx", ".xlsm", ".txt"}
MODEL_SUFFIXES = {".rdata", ".rds", ".pkl", ".pickle", ".parquet"}


def role_for(path: Path, manifest_type: str) -> str:
    if manifest_type == "input":
        return "input_data"
    suffix = path.suffix.lower()
    if suffix in FIGURE_SUFFIXES:
        return "output_figure"
    if suffix in TABLE_SUFFIXES:
        return "output_table"
    if suffix in MODEL_SUFFIXES:
        return "output_model"
    if "log" in path.parts:
        return "output_log"
    if "metadata" in path.parts:
        return "output_metadata"
    return "output_file"


def source_kind_for(path: Path, manifest_type: str) -> str:
    if manifest_type == "input":
        return "input_file"
    suffix = path.suffix.lower()
    if suffix in FIGURE_SUFFIXES:
        return "generated_panel"
    if suffix in TABLE_SUFFIXES:
        return "generated_table"
    if suffix in MODEL_SUFFIXES:
        return "generated_model"
    if "log" in path.parts:
        return "generated_log"
    if "metadata" in path.parts:
        return "generated_metadata"
    return "generated_file"


def iter_files(scan_dir: Path) -> Iterable[Path]:
    for path in sorted(scan_dir.rglob("*")):
        if not path.is_file():
            continue
        if path.name == "output_manifest.tsv":
            continue
        yield path


def manifest_row(
    path: Path,
    repo_root: Path,
    module: str,
    manifest_type: str,
    generated_by: str,
    command_id: str,
) -> dict[str, object]:
    resolved = path.resolve()
    repo_relative = str(resolved.relative_to(repo_root)) if path_within(resolved, repo_root) else ""
    checksum = sha256_file(resolved) if resolved.is_file() else ""
    return {
        "path": repo_relative or str(resolved),
        "repo_relative_path": repo_relative,
        "absolute_path": str(resolved),
        "role": role_for(resolved, manifest_type),
        "source_kind": source_kind_for(resolved, manifest_type),
        "module": module,
        "generated_by": generated_by,
        "command_id": command_id,
        "sha256": checksum,
        "checksum_unavailable_reason": "" if checksum else "not_a_regular_file",
        "byte_size": resolved.stat().st_size if resolved.exists() else "",
        "mtime_utc": datetime.fromtimestamp(resolved.stat().st_mtime, tz=timezone.utc).isoformat()
        if resolved.exists()
        else "",
        "figure": "",
        "panel": "",
        "notes": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-type", choices=("input", "output"), required=True)
    parser.add_argument("--module", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--generated-by", default="")
    parser.add_argument("--command-id", default="")
    parser.add_argument("--repo-root", type=Path)
    parser.add_argument("--scan-dir", type=Path, help="Directory to scan for output files.")
    parser.add_argument("--path", action="append", default=[], help="Input or output file path. May be repeated.")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root else repo_root_from(Path.cwd())
    paths: list[Path] = []
    if args.scan_dir is not None:
        scan_dir = args.scan_dir
        if not scan_dir.is_absolute():
            scan_dir = repo_root / scan_dir
        paths.extend(iter_files(scan_dir.resolve()))
    for raw_path in args.path:
        path = Path(raw_path)
        if not path.is_absolute():
            path = repo_root / path
        paths.append(path.resolve())

    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError("Manifest path(s) do not exist: " + ", ".join(missing))

    rows = [
        manifest_row(
            path,
            repo_root=repo_root,
            module=args.module,
            manifest_type=args.manifest_type,
            generated_by=args.generated_by,
            command_id=args.command_id,
        )
        for path in sorted(set(paths))
    ]
    if not rows:
        raise ValueError("No files were available to write to the manifest")
    write_tsv(args.output, rows, MODULE_MANIFEST_COLUMNS)
    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
