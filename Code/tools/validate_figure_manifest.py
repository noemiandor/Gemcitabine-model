#!/usr/bin/env python3
"""Validate a manuscript-facing figure manifest TSV."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

from figure_output_contract import repo_root_from, validate_figure_manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=None,
        help="Repository root. Defaults to nearest parent containing .git.",
    )
    args = parser.parse_args()

    repo_root = args.repo_root.resolve() if args.repo_root else repo_root_from(args.manifest)
    errors = validate_figure_manifest(args.manifest, repo_root=repo_root)
    if errors:
        for error in errors:
            print(error, file=sys.stderr)
        return 1
    print(f"Figure manifest OK: {args.manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
