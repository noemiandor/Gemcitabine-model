#!/usr/bin/env python3
"""Smoke-test saved PKPD/live-dead fit summaries and regenerated plot outputs."""

from __future__ import annotations

from pathlib import Path
import struct
import subprocess
import sys

import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_ROOT = REPO_ROOT / "Data" / "in-vitro" / "pkpd_live_dead_model"
FIT_OUTPUTS_ROOT = DATA_ROOT / "invitro_fitting_outputs"
PLOTTER = MODULE_ROOT / "plot_invitro_fit_outputs.py"
MANIFEST = FIT_OUTPUTS_ROOT / "saved_summary_output_manifest.tsv"

DEFAULT_FOLDER = FIT_OUTPUTS_ROOT / "alsoGoodFit_20260514T093906"
BEST_FIT_FOLDER = FIT_OUTPUTS_ROOT / "bestFitSoFar_20260513T164159"


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a valid PNG")
    width, height = struct.unpack(">II", header[16:24])
    return int(width), int(height)


def required_columns(raw_value: object) -> list[str]:
    if raw_value is None or pd.isna(raw_value):
        return []
    return [col.strip() for col in str(raw_value).split(",") if col.strip()]


def check_best_row_selectable(summary_path: Path) -> None:
    summary_df = pd.read_csv(summary_path, sep="\t")
    if summary_df.empty:
        raise AssertionError(f"{summary_path} has no rows")
    if "posterior_objective" not in summary_df.columns:
        raise AssertionError(f"{summary_path} lacks posterior_objective")
    objective = pd.to_numeric(summary_df["posterior_objective"], errors="coerce")
    if not objective.notna().any():
        raise AssertionError(f"{summary_path} has no numeric posterior_objective values")
    _ = summary_df.loc[objective.idxmin()]


def check_manifest_outputs(folder: Path, comparison_group: str, manifest_df: pd.DataFrame) -> None:
    active_groups = {"common", comparison_group}
    for _, row in manifest_df[manifest_df["requirement_group"].isin(active_groups)].iterrows():
        output_path = folder / str(row["output_name"])
        if not output_path.exists():
            raise AssertionError(f"Missing required output: {output_path}")
        if output_path.stat().st_size <= 0:
            raise AssertionError(f"Required output is empty: {output_path}")

        output_type = str(row["output_type"])
        if output_type == "tsv":
            table_columns = set(pd.read_csv(output_path, sep="\t", nrows=0).columns)
            missing = [col for col in required_columns(row.get("required_columns")) if col not in table_columns]
            if missing:
                raise AssertionError(f"{output_path} is missing required columns: {missing}")
        elif output_type == "png":
            width, height = png_dimensions(output_path)
            expected_width = row.get("width_px")
            expected_height = row.get("height_px")
            if pd.notna(expected_width) and int(expected_width) != width:
                raise AssertionError(f"{output_path} width {width} != expected {int(expected_width)}")
            if pd.notna(expected_height) and int(expected_height) != height:
                raise AssertionError(f"{output_path} height {height} != expected {int(expected_height)}")
        else:
            raise AssertionError(f"Unsupported manifest output_type={output_type!r}")


def run_plotter(args: list[str]) -> None:
    subprocess.run([sys.executable, str(PLOTTER), *args], cwd=REPO_ROOT, check=True)


def main() -> None:
    if not MANIFEST.exists():
        raise FileNotFoundError(f"Missing saved-summary output manifest: {MANIFEST}")
    manifest_df = pd.read_csv(MANIFEST, sep="\t")

    runs = [
        ([], DEFAULT_FOLDER, "dose_50"),
        ([str(BEST_FIT_FOLDER)], BEST_FIT_FOLDER, "dose_50"),
        ([str(BEST_FIT_FOLDER), "--comparison-dose", "25 nM"], BEST_FIT_FOLDER, "dose_25"),
    ]

    for args, output_folder, comparison_group in runs:
        run_plotter(args)
        check_best_row_selectable(output_folder / "joint_fit_summary.tsv")
        check_manifest_outputs(output_folder, comparison_group, manifest_df)

    print("Saved-summary smoke test passed.")


if __name__ == "__main__":
    main()
