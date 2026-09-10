#!/usr/bin/env python3
"""Check the curated dataset index. This is not an unused-file detector."""

import csv
from pathlib import Path
import subprocess
import sys


def main():
    root = Path(__file__).resolve().parents[2]
    index = root / "Data/figure_datasets.tsv"
    tracked = set(subprocess.check_output(
        ["git", "ls-files", "-z"], cwd=root
    ).decode().split("\0"))
    errors = []
    seen = set()
    with index.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        name, relative = row["dataset_id"], row["path"]
        path = root / relative
        if name in seen:
            errors.append(f"Duplicate dataset ID: {name}")
        seen.add(name)
        if Path(relative).is_absolute() or ".." in Path(relative).parts:
            errors.append(f"Nonportable dataset path: {relative}")
            continue
        if row["availability"] == "optional_local":
            print(f"OPTIONAL {name}: {'present' if path.exists() else 'absent'}")
            continue
        if row["availability"] != "git":
            errors.append(f"Unknown availability for {name}: {row['availability']}")
            continue
        if not path.exists():
            errors.append(f"Missing required dataset: {relative}")
            continue
        registered = [p for p in tracked if p == relative or p.startswith(relative + "/")]
        if not registered:
            errors.append(f"Dataset is not tracked/staged in Git: {relative}")
        for p in registered:
            if not (root / p).is_file():
                errors.append(f"Missing registered file: {p}")
        # A directory entry must not conceal local-only scientific tables.
        # Explicitly ignored caches and generated files remain outside this check.
        untracked = subprocess.check_output(
            ["git", "ls-files", "--others", "--exclude-standard", "-z", "--", relative],
            cwd=root,
        ).decode().split("\0")
        for p in untracked:
            if p and Path(p).suffix.lower() in {
                ".csv", ".tsv", ".txt", ".parquet", ".rds", ".xlsx", ".xlsm", ".fcs"
            }:
                errors.append(f"Untracked data inside registered dataset: {p}")
        print(f"GIT {name}: {len(registered)} registered file(s)")
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Dataset index OK: {len(rows)} datasets; required inputs are present in Git/index.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
