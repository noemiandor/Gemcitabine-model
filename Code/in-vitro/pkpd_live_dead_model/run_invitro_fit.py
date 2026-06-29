#!/usr/bin/env python3
"""Run or validate the imported in vitro PKPD/live-dead fitting workflow."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
import tempfile
from typing import Any, Dict, Iterable, Optional, Sequence

_CACHE_ROOT = Path(tempfile.gettempdir()) / "gemcitabine_model_pkpd_live_dead_model"
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg_cache"))

MODULE_ROOT = Path(__file__).resolve().parent
SRC_DIR = MODULE_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import pandas as pd

import invitro_fitting


REPO_ROOT = invitro_fitting.PROJECT_ROOT
DATA_ROOT = REPO_ROOT / "Data" / "in-vitro" / "pkpd_live_dead_model"
DEFAULT_SMOKE_MAX_NFEV = 25
DEFAULT_SMOKE_N_TR_VALUES = (2,)
DEFAULT_SMOKE_DOSE_LABELS = ("0 nM", "25 nM", "50 nM")
DEFAULT_SMOKE_FIT_T_MAX = 2.0


def path_within_repo(path: Path) -> bool:
    resolved = path.resolve()
    repo = REPO_ROOT.resolve()
    return resolved == repo or repo in resolved.parents


def repo_relative(path: Path) -> str:
    resolved = path.resolve()
    if not path_within_repo(resolved):
        raise ValueError(f"Path resolves outside this repository: {resolved}")
    return str(resolved.relative_to(REPO_ROOT.resolve()))


def parse_int_values(values: Optional[Sequence[str]]) -> Optional[tuple[int, ...]]:
    if values is None:
        return None
    parsed: list[int] = []
    for raw in values:
        for token in str(raw).split(","):
            token = token.strip()
            if token:
                parsed.append(int(token))
    if not parsed:
        raise ValueError("--n-tr-values must include at least one integer")
    return tuple(parsed)


def resolve_n_jobs(args: argparse.Namespace) -> Optional[int]:
    if args.n_jobs is not None and args.max_parallel is not None and args.n_jobs != args.max_parallel:
        raise ValueError("--n-jobs and --max-parallel both map to JointFitConfig.n_jobs and must agree")
    n_jobs = args.n_jobs if args.n_jobs is not None else args.max_parallel
    if args.smoke_test:
        if n_jobs not in (None, 1):
            print("Smoke test forces serial execution; using n_jobs=1.", file=sys.stderr)
        return 1
    return n_jobs


def build_paths(output_dir: Optional[Path]) -> invitro_fitting.ExperimentPaths:
    base = invitro_fitting.default_experiment_paths(project_root=REPO_ROOT)
    resolved_output_dir = base.output_dir
    if output_dir is not None:
        resolved_output_dir = output_dir
        if not resolved_output_dir.is_absolute():
            resolved_output_dir = REPO_ROOT / resolved_output_dir
        resolved_output_dir = resolved_output_dir.resolve()
    return invitro_fitting.ExperimentPaths(
        project_root=base.project_root,
        exp_base=base.exp_base,
        counts_raw=base.counts_raw,
        counts_agg=base.counts_agg,
        platemap=base.platemap,
        pkpd_constants=base.pkpd_constants,
        output_dir=resolved_output_dir,
    )


def guard_output_dir(output_dir: Path, overwrite: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()) and not overwrite:
        raise FileExistsError(f"Output directory already exists and is not empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)


def build_fit_config(args: argparse.Namespace) -> invitro_fitting.JointFitConfig:
    overrides: Dict[str, Any] = {}
    if args.max_days is not None:
        overrides["max_days"] = float(args.max_days)
    if args.fit_t_max is not None:
        overrides["fit_t_max"] = float(args.fit_t_max)
    if args.smoke_test and args.max_days is None:
        overrides["max_days"] = DEFAULT_SMOKE_FIT_T_MAX
    if args.smoke_test and args.fit_t_max is None:
        overrides["fit_t_max"] = DEFAULT_SMOKE_FIT_T_MAX

    n_tr_values = parse_int_values(args.n_tr_values)
    if args.smoke_test and n_tr_values is None:
        n_tr_values = DEFAULT_SMOKE_N_TR_VALUES
    if n_tr_values is not None:
        overrides["n_tr_values"] = n_tr_values

    max_nfev = args.max_nfev
    if args.smoke_test and max_nfev is None:
        max_nfev = DEFAULT_SMOKE_MAX_NFEV
    if max_nfev is not None:
        overrides["max_nfev"] = int(max_nfev)

    n_jobs = resolve_n_jobs(args)
    if n_jobs is not None:
        overrides["n_jobs"] = int(n_jobs)

    n_starts = args.n_starts
    if args.smoke_test and n_starts is None:
        n_starts = 1
    if n_starts is not None:
        n_starts = int(n_starts)
        if n_starts <= 0:
            raise ValueError("--n-starts must be positive")
        overrides["optimizer_start_indices"] = tuple(range(n_starts))

    if args.model_preset is not None:
        overrides["model_preset"] = args.model_preset
    if args.smoke_test:
        overrides["fit_dose_labels"] = DEFAULT_SMOKE_DOSE_LABELS
        overrides["max_replicates_per_ploidy_dose"] = 1
        overrides["accept_incomplete_optimizer"] = True

    return invitro_fitting.resolve_joint_fit_config(**overrides)


def check_expected_inputs(paths: invitro_fitting.ExperimentPaths, fit_config: invitro_fitting.JointFitConfig) -> Path:
    required_path_attrs = ("counts_agg", "platemap", "pkpd_constants")
    invitro_fitting.validate_paths(paths, required=required_path_attrs)
    for attr_name in ("project_root", "exp_base", *required_path_attrs):
        path_value = Path(getattr(paths, attr_name))
        if not path_within_repo(path_value):
            raise ValueError(f"{attr_name} resolves outside this repository: {path_value.resolve()}")

    counts_df = invitro_fitting.load_and_clean_counts(max_days=fit_config.max_days, paths=paths)
    platemap_df = invitro_fitting.load_platemap(paths=paths)
    pk_sheets = invitro_fitting.import_and_clean_pkpd(paths=paths)
    surfaces = invitro_fitting.build_preferred_dfdctp_signal_surfaces(
        pk_sheets,
        ploidy_keys=("2N", "4N"),
        fallback_half_life_days=invitro_fitting.PKConfig().fallback_half_life_days,
    )
    modeling_df = invitro_fitting.assemble_modeling_dataset(paths=paths, fit_config=fit_config)

    ploidies = sorted(str(value) for value in modeling_df["ploidy"].dropna().unique())
    dose_labels = sorted(
        (str(value) for value in modeling_df["gem"].dropna().unique()),
        key=lambda value: float(value.split()[0]),
    )
    required_ploidies = {"2N", "4N"}
    required_doses = {"0 nM", "25 nM", "50 nM"}
    if not required_ploidies.issubset(ploidies):
        raise ValueError(f"Missing expected ploidies: {sorted(required_ploidies.difference(ploidies))}")
    if not required_doses.issubset(dose_labels):
        raise ValueError(f"Missing expected dose labels: {sorted(required_doses.difference(dose_labels))}")

    rows: list[dict[str, Any]] = [
        {"resource": "counts_agg_path", "path": repo_relative(paths.counts_agg), "row_count": len(counts_df), "details": "cleaned aggregated count rows"},
        {"resource": "platemap_path", "path": repo_relative(paths.platemap), "row_count": len(platemap_df), "details": "standardized platemap wells"},
        {"resource": "pkpd_workbook_path", "path": repo_relative(paths.pkpd_constants), "row_count": sum(len(sheet) for sheet in pk_sheets.values()), "details": f"{len(pk_sheets)} workbook sheets"},
        {"resource": "modeling_dataset", "path": "", "row_count": len(modeling_df), "details": f"ploidies={','.join(ploidies)}; doses={','.join(dose_labels)}; max_days={fit_config.max_days:g}"},
    ]

    for ploidy, surface in surfaces.items():
        rows.append(
            {
                "resource": f"dfdctp_surface_{ploidy}",
                "path": "",
                "row_count": len(surface.calibration_profiles_by_dose),
                "details": f"calibrated_doses_uM={','.join(f'{dose:g}' for dose in surface.calibration_doses_uM)}",
            }
        )

    for ploidy in sorted(required_ploidies):
        for dose_label in dose_labels:
            aligned = invitro_fitting.get_aligned_live_dead_data(
                modeling_df,
                gem_dose=dose_label,
                ploidy=ploidy,
                t_max=fit_config.fit_t_max,
                count_transitional_as_alive=fit_config.count_transitional_as_alive,
            )
            replicate_count = len(aligned["replicate_columns"])
            if replicate_count <= 0:
                raise ValueError(f"No aligned replicates for {ploidy} {dose_label}")
            rows.append(
                {
                    "resource": f"aligned_live_dead_{ploidy}_{invitro_fitting.slugify_label(dose_label)}",
                    "path": "",
                    "row_count": len(aligned["t"]),
                    "details": (
                        f"replicates={replicate_count}; "
                        f"time_min={float(aligned['t'].min()):g}; "
                        f"time_max={float(aligned['t'].max()):g}; "
                        f"dropped_timepoints={aligned['dropped_timepoints']}; "
                        f"dropped_replicates={aligned['dropped_replicates']}"
                    ),
                }
            )

    manifest_path = DATA_ROOT / "input_manifest.tsv"
    pd.DataFrame(rows).to_csv(manifest_path, sep="\t", index=False)
    return manifest_path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, help="Output directory for fitting results. Defaults to a timestamped folder under Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs.")
    parser.add_argument("--max-days", type=float, help="Maximum live/dead observation time, mapped to JointFitConfig.max_days.")
    parser.add_argument("--fit-t-max", type=float, help="Maximum aligned fitting time, mapped to JointFitConfig.fit_t_max.")
    parser.add_argument("--n-starts", type=int, help="Use the first N optimizer starts from the existing start grid.")
    parser.add_argument("--max-parallel", type=int, help="Compatibility alias for JointFitConfig.n_jobs.")
    parser.add_argument("--n-jobs", type=int, help="Number of n_tr worker processes; use 1 for serial execution.")
    parser.add_argument("--n-tr-values", action="append", help="n_tr values to evaluate, as a comma-separated list or repeated option.")
    parser.add_argument("--max-nfev", type=int, help="Maximum optimizer iterations, mapped to JointFitConfig.max_nfev.")
    parser.add_argument("--model-preset", choices=sorted(invitro_fitting.ALLOWED_MODEL_PRESETS), help="Model preset resolved by src/invitro_fitting.py.")
    parser.add_argument("--smoke-test", action="store_true", help="Run a cheap serial fit: one n_tr value, one optimizer start by default, and reduced max_nfev.")
    parser.add_argument("--check-inputs", action="store_true", help="Validate input paths and data-loading/PKPD alignment, then write input_manifest.tsv without fitting.")
    parser.add_argument("--overwrite", action="store_true", help="Allow writing into a nonempty --output-dir.")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    paths = build_paths(args.output_dir)
    fit_config = build_fit_config(args)

    if args.check_inputs:
        manifest_path = check_expected_inputs(paths, fit_config)
        print(f"Input checks passed. Wrote {manifest_path}")
        return 0

    guard_output_dir(paths.output_dir, overwrite=args.overwrite)
    result = invitro_fitting.main(paths=paths, fit_config=fit_config)

    required_outputs = ["joint_fit_summary.tsv", "optimizer_attempts.tsv"]
    missing = [name for name in required_outputs if not (paths.output_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Fitting run did not produce required outputs: {missing}")
    print(f"Fit output directory: {paths.output_dir}")
    return 0 if result.get("success") else 2


if __name__ == "__main__":
    raise SystemExit(main())
