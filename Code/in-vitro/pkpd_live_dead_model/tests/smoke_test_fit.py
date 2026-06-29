#!/usr/bin/env python3
"""Run a bounded PKPD/live-dead smoke fit and validate plotter compatibility."""

from __future__ import annotations

from pathlib import Path
import struct
import subprocess
import sys
import tempfile

import pandas as pd


MODULE_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = Path(__file__).resolve().parents[4]
DATA_ROOT = REPO_ROOT / "Data" / "in-vitro" / "pkpd_live_dead_model"
FIT_OUTPUTS_ROOT = DATA_ROOT / "invitro_fitting_outputs"
RUNNER = MODULE_ROOT / "run_invitro_fit.py"
PLOTTER = MODULE_ROOT / "plot_invitro_fit_outputs.py"
FITTING_SOURCE = MODULE_ROOT / "src" / "invitro_fitting.py"
MANIFEST = FIT_OUTPUTS_ROOT / "saved_summary_output_manifest.tsv"
IMPORTED_SUMMARY = FIT_OUTPUTS_ROOT / "alsoGoodFit_20260514T093906" / "joint_fit_summary.tsv"


def png_dimensions(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.read(24)
    if len(header) < 24 or header[:8] != b"\x89PNG\r\n\x1a\n":
        raise AssertionError(f"{path} is not a valid PNG")
    return struct.unpack(">II", header[16:24])


def required_columns(raw_value: object) -> list[str]:
    if raw_value is None or pd.isna(raw_value):
        return []
    return [value.strip() for value in str(raw_value).split(",") if value.strip()]


def assert_manifest_outputs(output_dir: Path) -> None:
    manifest_df = pd.read_csv(MANIFEST, sep="\t")
    active_groups = {"common", "dose_50", "dose_25"}
    for _, row in manifest_df[manifest_df["requirement_group"].isin(active_groups)].iterrows():
        output_path = output_dir / str(row["output_name"])
        if not output_path.exists():
            raise AssertionError(f"Missing required output: {output_path}")
        if output_path.stat().st_size <= 0:
            raise AssertionError(f"Required output is empty: {output_path}")

        if row["output_type"] == "tsv":
            columns = set(pd.read_csv(output_path, sep="\t", nrows=0).columns)
            missing = [col for col in required_columns(row["required_columns"]) if col not in columns]
            if missing:
                raise AssertionError(f"{output_path} is missing required columns: {missing}")
        elif row["output_type"] == "png":
            width, height = png_dimensions(output_path)
            expected_width = row.get("width_px")
            expected_height = row.get("height_px")
            if pd.notna(expected_width) and int(expected_width) != width:
                raise AssertionError(f"{output_path} width {width} != expected {int(expected_width)}")
            if pd.notna(expected_height) and int(expected_height) != height:
                raise AssertionError(f"{output_path} height {height} != expected {int(expected_height)}")


def assert_summary_compatible(output_dir: Path) -> None:
    imported_columns = set(pd.read_csv(IMPORTED_SUMMARY, sep="\t", nrows=0).columns)
    summary = pd.read_csv(output_dir / "joint_fit_summary.tsv", sep="\t")
    attempts = pd.read_csv(output_dir / "optimizer_attempts.tsv", sep="\t")
    if summary.empty:
        raise AssertionError("Smoke-fit joint_fit_summary.tsv has no rows")
    if attempts.empty:
        raise AssertionError("Smoke-fit optimizer_attempts.tsv has no rows")

    missing_imported_columns = sorted(imported_columns.difference(summary.columns))
    if missing_imported_columns:
        raise AssertionError(f"Smoke summary is missing imported-summary columns: {missing_imported_columns}")

    objective = pd.to_numeric(summary["posterior_objective"], errors="coerce")
    if not objective.notna().any():
        raise AssertionError("Smoke summary lacks numeric posterior_objective values")
    best_row = summary.loc[objective.idxmin()]
    required_best_row_fields = [
        "posterior_objective",
        "objective",
        "observation_channels",
        "model_variant",
        "n_tr",
        "theta_alive",
        "theta_dead",
        "2N_r",
        "4N_r",
        "2N_K",
        "4N_K",
        "2N_k_tr",
        "4N_k_tr",
        "2N_k_kill",
        "4N_k_kill",
        "2N_k_clear",
        "4N_k_clear",
    ]
    missing = [field for field in required_best_row_fields if field not in best_row.index or pd.isna(best_row[field])]
    if missing:
        raise AssertionError(f"Best smoke-fit row is missing required plotting fields: {missing}")

    required_attempt_columns = {"start_idx", "final_raw_objective", "optimizer_success", "success"}
    missing_attempt_columns = sorted(required_attempt_columns.difference(attempts.columns))
    if missing_attempt_columns:
        raise AssertionError(f"optimizer_attempts.tsv missing columns: {missing_attempt_columns}")


def assert_visual_contract_source_guards() -> None:
    plotter_source = PLOTTER.read_text()
    fitting_source = FITTING_SOURCE.read_text()
    required_plotter_snippets = [
        'ALIVE_OBS_COLOR = "#4D4D4D"',
        'ALIVE_MODEL_COLOR = "#000000"',
        'DEAD_OBS_COLOR = "#F28E8E"',
        'DEAD_MODEL_COLOR = "#D62728"',
        "cols = 2",
        "rows = max(5",
        "MaxNLocator(integer=True)",
        'ax.set_xscale("log")',
    ]
    for snippet in required_plotter_snippets:
        if snippet not in plotter_source:
            raise AssertionError(f"Plotter source no longer contains expected visual-contract snippet: {snippet}")
    if 'ax.set_xscale("log")' not in fitting_source:
        raise AssertionError("Single-ploidy dFdCTP signal plot no longer sets a log x-axis")


def run_command(args: list[str]) -> None:
    subprocess.run([sys.executable, *args], cwd=REPO_ROOT, check=True)


def main() -> None:
    assert_visual_contract_source_guards()
    with tempfile.TemporaryDirectory(prefix="gemcitabine_pkpd_smoke_fit_") as tmpdir:
        output_dir = Path(tmpdir) / "fit_output"
        run_command([
            str(RUNNER),
            "--output-dir",
            str(output_dir),
            "--smoke-test",
            "--n-starts",
            "1",
            "--max-parallel",
            "1",
        ])
        assert_summary_compatible(output_dir)
        run_command([str(PLOTTER), str(output_dir)])
        run_command([str(PLOTTER), str(output_dir), "--comparison-dose", "25 nM"])
        assert_manifest_outputs(output_dir)
    print("Smoke-fit regression test passed.")


if __name__ == "__main__":
    main()
