#!/usr/bin/env python3
"""Shared helpers for manuscript figure output manifests."""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import re
import tempfile
from typing import Iterable, Mapping, Sequence


RUN_ID_RE = re.compile(r"^[A-Za-z0-9._-]+$")

MODULE_MANIFEST_COLUMNS = (
    "path",
    "repo_relative_path",
    "absolute_path",
    "role",
    "source_kind",
    "module",
    "generated_by",
    "command_id",
    "sha256",
    "checksum_unavailable_reason",
    "byte_size",
    "mtime_utc",
    "figure",
    "panel",
    "notes",
)

FIGURE_MANIFEST_COLUMNS = (
    "figure",
    "panel",
    "asset_path",
    "source_file",
    "source_kind",
    "generated_by",
    "command",
    "input_data",
    "result_run_dir",
    "run_id",
    "caption_role",
    "asset_status",
    "not_regenerated_reason",
    "local_provenance_path",
    "citation_or_uri",
    "notes",
)

GENERATED_SOURCE_KINDS = {"generated_panel", "generated_table"}
NONLOCAL_SOURCE_KINDS = {"external", "placeholder"}
MANUAL_SOURCE_KINDS = {"manual_composite", "external", "placeholder"}


def validate_run_id(run_id: str) -> None:
    if not run_id:
        raise ValueError("run_id is required")
    if not RUN_ID_RE.match(run_id):
        raise ValueError(
            "run_id may contain only letters, numbers, dots, underscores, and hyphens"
        )


def repo_root_from(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).resolve()
    for candidate in (current, *current.parents):
        if (candidate / ".git").exists():
            return candidate
    raise FileNotFoundError(f"Could not find repository root from {current}")


def resolve_repo_path(value: str, repo_root: Path) -> Path | None:
    value = (value or "").strip()
    if not value:
        return None
    path = Path(value)
    if not path.is_absolute():
        path = repo_root / path
    return path.resolve()


def path_within(path: Path, root: Path) -> bool:
    path = path.resolve()
    root = root.resolve()
    return path == root or root in path.parents


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            return [], []
        rows = [{key: (value if value is not None else "") for key, value in row.items()} for row in reader]
        return list(reader.fieldnames), rows


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            out_row = {key: row.get(key, "") for key in columns}
            if columns and not str(out_row.get(columns[-1], "")).strip():
                out_row[columns[-1]] = "."
            writer.writerow(out_row)


def ensure_columns(headers: Sequence[str], required: Sequence[str], manifest_path: Path) -> list[str]:
    missing = [column for column in required if column not in headers]
    if missing:
        return [f"{manifest_path}: missing required column(s): {', '.join(missing)}"]
    return []


def first_local_path(row: Mapping[str, str], repo_root: Path, keys: Iterable[str]) -> Path | None:
    first_resolved: Path | None = None
    for key in keys:
        path = resolve_repo_path(row.get(key, ""), repo_root)
        if path is None:
            continue
        if first_resolved is None:
            first_resolved = path
        if path.exists():
            return path
    return first_resolved


def module_manifest_local_path(
    row: Mapping[str, str],
    repo_root: Path,
    output_root: Path | None = None,
) -> Path | None:
    """Resolve a module-manifest row, including output-root-relative locators."""
    for key in ("repo_relative_path", "absolute_path"):
        value = (row.get(key, "") or "").strip()
        if not value:
            continue
        path = Path(value)
        candidate = path if path.is_absolute() else repo_root / path
        if candidate.exists():
            return candidate.resolve()

    value = (row.get("path", "") or "").strip()
    if not value or value.startswith("external:"):
        return None
    path = Path(value)
    if path.is_absolute():
        return path.resolve()
    if output_root is not None and row.get("role", "").startswith("output"):
        return (output_root / path).resolve()
    return (repo_root / path).resolve()


def resolve_figure_locator(
    value: str,
    repo_root: Path,
    manifest_parent: Path,
) -> Path | None:
    value = (value or "").strip()
    if not value or value.startswith("external:"):
        return None
    prefix = "manifest:"
    if value.startswith(prefix):
        relative = value[len(prefix) :]
        if not relative or Path(relative).is_absolute():
            return None
        return (manifest_parent / relative).resolve()
    return resolve_repo_path(value, repo_root)


def validate_expected_panel_set(
    rows: Sequence[Mapping[str, str]],
    expected_generated_panels: Iterable[str] | None,
    manifest_path: Path,
) -> list[str]:
    """Validate unique panel IDs and an optional generated-panel allow-list."""
    errors: list[str] = []
    counts: dict[str, int] = {}
    for row in rows:
        panel = row.get("panel", "").strip()
        if panel:
            counts[panel] = counts.get(panel, 0) + 1
    for panel, count in sorted(counts.items()):
        if count > 1:
            errors.append(f"{manifest_path}: duplicate panel ID {panel!r} appears {count} times")
    if expected_generated_panels is None:
        return errors
    expected = set(expected_generated_panels)
    observed = {
        row.get("panel", "").strip()
        for row in rows
        if row.get("source_kind", "").strip() in GENERATED_SOURCE_KINDS
        and row.get("panel", "").strip()
    }
    missing = sorted(expected - observed)
    unexpected = sorted(observed - expected)
    if missing:
        errors.append(f"{manifest_path}: missing expected generated panel(s): {', '.join(missing)}")
    if unexpected:
        errors.append(f"{manifest_path}: unexpected generated panel(s): {', '.join(unexpected)}")
    return errors


def validate_module_manifest(path: Path, repo_root: Path, output_root: Path | None = None) -> list[str]:
    errors: list[str] = []
    headers, rows = read_tsv(path)
    errors.extend(ensure_columns(headers, MODULE_MANIFEST_COLUMNS, path))
    if errors:
        return errors
    if not rows:
        return [f"{path}: manifest has no data rows"]

    output_root = output_root.resolve() if output_root is not None else None
    for idx, row in enumerate(rows, start=2):
        for required in ("role", "source_kind", "module"):
            if not row.get(required, "").strip():
                errors.append(f"{path}:{idx}: {required} is required")

        # Prefer the portable checkout-relative location. The recorded absolute
        # path is only a fallback for genuinely external inputs.
        local_path = module_manifest_local_path(
            row,
            repo_root,
            output_root=output_root,
        )
        source_kind = row.get("source_kind", "").strip()
        if local_path is not None and source_kind not in NONLOCAL_SOURCE_KINDS:
            if not local_path.exists():
                errors.append(f"{path}:{idx}: local path does not exist: {local_path}")
            elif output_root is not None and row.get("role", "").startswith("output") and not path_within(local_path, output_root):
                errors.append(f"{path}:{idx}: output path is outside output root: {local_path}")

            declared_size = row.get("byte_size", "").strip()
            if declared_size and local_path.exists():
                try:
                    expected_size = int(declared_size)
                    if local_path.stat().st_size != expected_size:
                        errors.append(f"{path}:{idx}: byte_size does not match file size for {local_path}")
                except ValueError:
                    errors.append(f"{path}:{idx}: byte_size is not an integer")

            declared_sha = row.get("sha256", "").strip()
            if declared_sha and local_path.exists() and local_path.is_file():
                actual_sha = sha256_file(local_path)
                if actual_sha != declared_sha:
                    errors.append(f"{path}:{idx}: sha256 does not match file contents for {local_path}")
            elif not declared_sha and local_path is not None:
                reason = row.get("checksum_unavailable_reason", "").strip()
                if not reason:
                    errors.append(f"{path}:{idx}: sha256 or checksum_unavailable_reason is required")

    return errors


def validate_figure_manifest(path: Path, repo_root: Path) -> list[str]:
    errors: list[str] = []
    headers, rows = read_tsv(path)
    errors.extend(ensure_columns(headers, FIGURE_MANIFEST_COLUMNS, path))
    if errors:
        return errors
    if not rows:
        return [f"{path}: manifest has no data rows"]

    errors.extend(validate_expected_panel_set(rows, None, path))

    for idx, row in enumerate(rows, start=2):
        source_kind = row.get("source_kind", "").strip()
        for required in ("figure", "panel", "source_kind", "asset_status"):
            if not row.get(required, "").strip():
                errors.append(f"{path}:{idx}: {required} is required")

        asset_path = resolve_figure_locator(
            row.get("asset_path", ""),
            repo_root,
            path.parent,
        )
        if source_kind in GENERATED_SOURCE_KINDS:
            for required in ("asset_path", "source_file", "generated_by", "command", "input_data", "result_run_dir", "run_id"):
                if not row.get(required, "").strip():
                    errors.append(f"{path}:{idx}: {required} is required for {source_kind}")
            if asset_path is None or not asset_path.exists():
                errors.append(f"{path}:{idx}: generated asset is missing: {row.get('asset_path', '')}")
            elif not path_within(asset_path, path.parent):
                errors.append(
                    f"{path}:{idx}: generated asset escapes its figure directory: {asset_path}"
                )
            source_path = resolve_figure_locator(
                row.get("source_file", ""),
                repo_root,
                path.parent,
            )
            if source_path is None or not source_path.exists():
                errors.append(f"{path}:{idx}: generated source_file is missing: {row.get('source_file', '')}")
            result_dir = resolve_figure_locator(
                row.get("result_run_dir", ""),
                repo_root,
                path.parent,
            )
            if result_dir is None or not result_dir.exists():
                errors.append(f"{path}:{idx}: result_run_dir is missing: {row.get('result_run_dir', '')}")
            elif source_path is not None and not path_within(source_path, result_dir):
                errors.append(
                    f"{path}:{idx}: generated source_file is outside result_run_dir: {source_path}"
                )
            input_path = resolve_figure_locator(
                row.get("input_data", ""),
                repo_root,
                path.parent,
            )
            if (
                input_path is None
                or source_path is None
                or input_path != source_path
            ):
                errors.append(
                    f"{path}:{idx}: input_data must resolve to source_file"
                )
        elif source_kind in MANUAL_SOURCE_KINDS:
            if not row.get("not_regenerated_reason", "").strip():
                errors.append(f"{path}:{idx}: not_regenerated_reason is required for {source_kind}")
            if asset_path is not None and row.get("asset_status", "").strip() != "external" and not asset_path.exists():
                errors.append(f"{path}:{idx}: local manual asset is missing: {asset_path}")
        else:
            errors.append(f"{path}:{idx}: unknown source_kind: {source_kind}")

    return errors


def atomic_write_latest(latest_path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    latest_path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False, dir=latest_path.parent, newline="") as handle:
        tmp_path = Path(handle.name)
        writer = csv.DictWriter(handle, fieldnames=("key", "value"), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({"key": row.get("key", ""), "value": row.get("value", "")})
    os.replace(tmp_path, latest_path)
