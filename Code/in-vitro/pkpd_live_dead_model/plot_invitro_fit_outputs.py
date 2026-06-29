#!/usr/bin/env python3
"""
Plot manuscript-oriented summaries from an invitro_fitting output folder.

Usage:
    python Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py
    python Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
        Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159
    python Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
        --fit-output Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906 \
        --output-dir Results/in-vitro/pkpd_live_dead_model/runs/saved_fit
"""

from __future__ import annotations

import argparse
from datetime import datetime
import os
from pathlib import Path
import sys
import tempfile
from typing import Any, Dict, Iterable, Optional, Tuple

_CACHE_ROOT = Path(tempfile.gettempdir()) / "gemcitabine_model_pkpd_live_dead_model"
os.environ.setdefault("MPLCONFIGDIR", str(_CACHE_ROOT / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_CACHE_ROOT / "xdg_cache"))

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
import numpy as np
import pandas as pd

MODULE_ROOT = Path(__file__).resolve().parent
SRC_DIR = MODULE_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import invitro_fitting

REPO_ROOT = invitro_fitting.PROJECT_ROOT
DATA_ROOT = REPO_ROOT / "Data" / "in-vitro" / "pkpd_live_dead_model"
FIT_OUTPUTS_ROOT = DATA_ROOT / "invitro_fitting_outputs"

DEFAULT_FIT_OUTPUT_FOLDER = (
    FIT_OUTPUTS_ROOT
    / "alsoGoodFit_20260514T093906"
)

# The archived alsoGoodFit dose-response summaries used the 4-day endpoint even
# though the current fitting default is 5 days.
DEFAULT_FIT_T_MAX = 4.0


DEFAULT_PARAMETERS = [
    "r",
    "K",
    "k_tr",
    "k_kill",
    "k_clear",
    "k_cyto",
    "beta_dose",
    "mu_base_death",
    "mu_confluence_death",
]


def default_output_folder() -> Path:
    return (
        REPO_ROOT
        / "Results"
        / "in-vitro"
        / "pkpd_live_dead_model"
        / "runs"
        / f"{datetime.now().strftime('%Y%m%dT%H%M%S')}_pkpd_saved_fit"
    )


DISPLAY_NAMES = {
    "r": "r",
    "K": "K",
    "k_tr": "k_tr",
    "k_kill": "k_kill",
    "k_clear": "k_clear",
    "k_cyto": "k_cyto",
    "beta_dose": "beta_dose",
    "mu_base_death": "mu_base_death",
    "mu_confluence_death": "mu_confluence_death",
}


PLOIDY_STYLES = {
    "2N": {"linestyle": "-", "marker": "o"},
    "4N": {"linestyle": "--", "marker": "s"},
}


ALIVE_OBS_COLOR = "#4D4D4D"
ALIVE_MODEL_COLOR = "#000000"
DEAD_OBS_COLOR = "#F28E8E"
DEAD_MODEL_COLOR = "#D62728"


def choose_best_row(summary_df: pd.DataFrame) -> pd.Series:
    if "posterior_objective" in summary_df.columns:
        values = pd.to_numeric(summary_df["posterior_objective"], errors="coerce")
        if values.notna().any():
            return summary_df.loc[values.idxmin()]
    if len(summary_df) == 0:
        raise ValueError("joint_fit_summary.tsv has no rows")
    return summary_df.iloc[0]


def row_bool(row: pd.Series, key: str, default: bool = False) -> bool:
    if key not in row.index or pd.isna(row[key]):
        return default
    value = row[key]
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def row_float(row: pd.Series, key: str, default: float = np.nan) -> float:
    if key not in row.index:
        return default
    value = pd.to_numeric(pd.Series([row[key]]), errors="coerce").iloc[0]
    return float(value) if np.isfinite(value) else default


def build_comparison_table(best_row: pd.Series, parameters: Iterable[str]) -> pd.DataFrame:
    rows = []
    for parameter in parameters:
        col_2n = f"2N_{parameter}"
        col_4n = f"4N_{parameter}"
        if col_2n not in best_row.index or col_4n not in best_row.index:
            continue
        value_2n = pd.to_numeric(pd.Series([best_row[col_2n]]), errors="coerce").iloc[0]
        value_4n = pd.to_numeric(pd.Series([best_row[col_4n]]), errors="coerce").iloc[0]
        if not np.isfinite(value_2n) or not np.isfinite(value_4n):
            continue
        if value_2n <= 0 or value_4n <= 0:
            continue
        ratio = float(value_4n / value_2n)
        rows.append(
            {
                "parameter": parameter,
                "display_name": DISPLAY_NAMES.get(parameter, parameter),
                "value_2N": float(value_2n),
                "value_4N": float(value_4n),
                "ratio_4N_over_2N": ratio,
                "log2_ratio_4N_over_2N": float(np.log2(ratio)),
            }
        )
    if not rows:
        raise ValueError("No positive paired 2N/4N parameter columns were found")
    return pd.DataFrame(rows)


def plot_fold_change(comparison_df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(comparison_df))
    values = comparison_df["log2_ratio_4N_over_2N"].to_numpy(dtype=float)
    colors = np.where(values >= 0.0, "#4C78A8", "#F58518")

    ax.bar(x, values, color=colors, edgecolor="black", linewidth=0.6)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels(comparison_df["display_name"], rotation=35, ha="right")
    ax.set_ylabel("log2(4N / 2N)")
    ax.set_title("Fitted Parameter Fold-Change: 4N vs 2N")
    ax.grid(axis="y", linestyle="--", alpha=0.35)

    for idx, value in enumerate(values):
        label_y = value + (0.06 if value >= 0 else -0.06)
        va = "bottom" if value >= 0 else "top"
        ax.text(idx, label_y, f"{2 ** value:.2g}x", ha="center", va=va, fontsize=8)

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_paired_values(comparison_df: pd.DataFrame, output_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5.5))
    x_2n = np.zeros(len(comparison_df))
    x_4n = np.ones(len(comparison_df))

    for idx, row in comparison_df.iterrows():
        y_values = [row["value_2N"], row["value_4N"]]
        ax.plot([0, 1], y_values, color="#888888", linewidth=1.2, alpha=0.8)
        ax.scatter([0], [row["value_2N"]], color="#4C78A8", s=45, zorder=3)
        ax.scatter([1], [row["value_4N"]], color="#F58518", s=45, zorder=3)
        geometric_midpoint = np.sqrt(row["value_2N"] * row["value_4N"])
        ax.text(1.05, geometric_midpoint, row["display_name"], va="center", fontsize=9)

    ax.set_yscale("log")
    ax.set_xlim(-0.25, 1.85)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["2N", "4N"])
    ax.set_ylabel("Fitted parameter value (log scale)")
    ax.set_title("Paired Fitted Parameters by Ploidy")
    ax.grid(axis="y", which="both", linestyle="--", alpha=0.35)

    # Keep explicit arrays in the plot code so future readers can see that all
    # parameters share the same two ploidy positions.
    _ = x_2n, x_4n

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_dfdctp_surfaces() -> Dict[str, invitro_fitting.DfdctpSignalSurface]:
    paths = invitro_fitting.default_experiment_paths()
    pk_config = invitro_fitting.PKConfig()
    pk_sheets = invitro_fitting.import_and_clean_pkpd(paths=paths, pk_config=pk_config)
    return invitro_fitting.build_preferred_dfdctp_signal_surfaces(
        pk_sheets,
        ploidy_keys=("2N", "4N"),
        fallback_half_life_days=pk_config.fallback_half_life_days,
    )


def config_from_best_row(
    best_row: pd.Series,
    fit_t_max: Optional[float] = None,
) -> invitro_fitting.JointFitConfig:
    model_preset = best_row.get("model_preset", None)
    if pd.isna(model_preset):
        model_preset = None
    kwargs: Dict[str, Any] = {
        "objective": str(best_row["objective"]),
        "observation_channels": str(best_row["observation_channels"]),
        "model_variant": str(best_row["model_variant"]),
        "model_preset": model_preset,
        "fit_beta_dose": row_bool(best_row, "fit_beta_dose", True),
        "fixed_beta_dose": row_float(best_row, "fixed_beta_dose", 1.0),
        "use_hill_dose_gate": row_bool(best_row, "use_hill_dose_gate", False),
        "fit_hill_dose_gate": row_bool(best_row, "fit_hill_dose_gate", True),
        "use_confluence_death": row_bool(best_row, "use_confluence_death", False),
        "fit_mu_confluence_death": row_bool(best_row, "fit_mu_confluence_death", True),
        "fixed_mu_confluence_death": row_float(best_row, "fixed_mu_confluence_death", 0.0),
        "confluence_death_exponent": row_float(best_row, "confluence_death_exponent", 4.0),
    }
    if fit_t_max is not None:
        kwargs["fit_t_max"] = float(fit_t_max)
        kwargs["max_days"] = float(fit_t_max)
    return invitro_fitting.resolve_joint_fit_config(
        invitro_fitting.JointFitConfig(**kwargs)
    )


def ploidy_parameter_dict(
    best_row: pd.Series,
    ploidy: str,
    fit_config: invitro_fitting.JointFitConfig,
) -> Dict[str, float]:
    params: Dict[str, float] = {}
    for name in invitro_fitting.get_parameter_names_for_config(fit_config):
        col = f"{ploidy}_{name}"
        if col in best_row.index:
            value = pd.to_numeric(pd.Series([best_row[col]]), errors="coerce").iloc[0]
            if np.isfinite(value):
                params[name] = float(value)
    return params


def interpolate_series_at_time(t: np.ndarray, values: np.ndarray, t_eval: float) -> float:
    if len(t) == 0:
        return np.nan
    return float(np.interp(float(t_eval), np.asarray(t, dtype=float), np.asarray(values, dtype=float)))


def mean_observed_series(matrix: np.ndarray) -> np.ndarray:
    arr = np.asarray(matrix, dtype=float)
    if arr.ndim == 1:
        return arr
    return np.nanmean(arr, axis=1)


def simulate_endpoint_metrics(
    fit_config: invitro_fitting.JointFitConfig,
    dfdctp_signal_curve: invitro_fitting.DfdctpSignalSurface,
    ploidy_params: Dict[str, float],
    dose_uM: float,
    N0: float,
    D0: float,
    endpoint_days: float,
    n_tr: int,
    dose_gate_ec50_uM: float,
    dose_gate_hill: float,
) -> Tuple[float, float]:
    treatment_params = [
        ploidy_params[name]
        for name in invitro_fitting.get_treatment_parameter_names_for_config(fit_config)
    ]
    alive, dead = invitro_fitting.simulate_joint_ext(
        np.array([0.0, float(endpoint_days)], dtype=float),
        treatment_params,
        float(N0),
        float(D0),
        float(ploidy_params["r"]),
        float(ploidy_params["K"]),
        float(dose_uM),
        dfdctp_signal_curve,
        int(n_tr),
        model_variant=fit_config.model_variant,
        fit_config=fit_config,
        beta_dose=float(ploidy_params.get("beta_dose", fit_config.fixed_beta_dose)),
        use_confluence_death=bool(fit_config.use_confluence_death),
        mu_confluence_death=float(
            ploidy_params.get(
                "mu_confluence_death",
                fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
            )
        ),
        confluence_death_exponent=float(fit_config.confluence_death_exponent),
        use_hill_dose_gate=bool(fit_config.use_hill_dose_gate),
        dose_gate_ec50_uM=float(dose_gate_ec50_uM) if fit_config.use_hill_dose_gate else None,
        dose_gate_hill=float(dose_gate_hill) if fit_config.use_hill_dose_gate else None,
    )
    return float(alive[-1]), float(dead[-1])


def build_dose_response_comparison_table(
    best_row: pd.Series,
    fit_config: invitro_fitting.JointFitConfig,
    curves_by_ploidy: Dict[str, invitro_fitting.DfdctpSignalSurface],
) -> pd.DataFrame:
    paths = invitro_fitting.default_experiment_paths()
    df = invitro_fitting.assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    gem_doses = sorted(
        [dose for dose in df["gem"].dropna().unique()],
        key=lambda x: float(str(x).split()[0]),
    )
    gate_ec50 = pd.to_numeric(pd.Series([best_row.get("dose_gate_ec50_uM", np.nan)]), errors="coerce").iloc[0]
    gate_hill = pd.to_numeric(pd.Series([best_row.get("dose_gate_hill", np.nan)]), errors="coerce").iloc[0]
    n_tr = int(best_row["n_tr"])
    ploidy_params = {ploidy: ploidy_parameter_dict(best_row, ploidy, fit_config) for ploidy in ("2N", "4N")}
    rows = []

    for dose_label in gem_doses:
        dose_uM = float(str(dose_label).split()[0]) / 1000.0
        if dose_uM <= 0:
            continue
        aligned = {}
        for ploidy in ("2N", "4N"):
            try:
                aligned[ploidy] = invitro_fitting.get_aligned_live_dead_data(
                    df,
                    gem_dose=dose_label,
                    ploidy=ploidy,
                    t_max=fit_config.fit_t_max,
                    count_transitional_as_alive=fit_config.count_transitional_as_alive,
                )
            except ValueError:
                aligned[ploidy] = None
        if aligned["2N"] is None or aligned["4N"] is None:
            continue
        endpoint_days = min(float(np.max(aligned["2N"]["t"])), float(np.max(aligned["4N"]["t"])))

        control_aligned = {}
        for ploidy in ("2N", "4N"):
            control_aligned[ploidy] = invitro_fitting.get_aligned_live_dead_data(
                df,
                gem_dose="0 nM",
                ploidy=ploidy,
                t_max=fit_config.fit_t_max,
                count_transitional_as_alive=fit_config.count_transitional_as_alive,
            )

        row = {
            "dose_label": dose_label,
            "dose_uM": dose_uM,
            "endpoint_days": endpoint_days,
        }
        pred_alive_vs_control = {}
        obs_alive_vs_control = {}
        for ploidy in ("2N", "4N"):
            dose_data = aligned[ploidy]
            control_data = control_aligned[ploidy]
            N0 = float(np.nanmean(dose_data["y_alive"][0, :]))
            D0 = float(np.nanmean(dose_data["y_dead"][0, :]))
            control_N0 = float(np.nanmean(control_data["y_alive"][0, :]))
            control_D0 = float(np.nanmean(control_data["y_dead"][0, :]))

            pred_alive_end, pred_dead_end = simulate_endpoint_metrics(
                fit_config=fit_config,
                dfdctp_signal_curve=curves_by_ploidy[ploidy],
                ploidy_params=ploidy_params[ploidy],
                dose_uM=dose_uM,
                N0=N0,
                D0=D0,
                endpoint_days=endpoint_days,
                n_tr=n_tr,
                dose_gate_ec50_uM=float(gate_ec50) if np.isfinite(gate_ec50) else float(fit_config.fixed_dose_gate_ec50_uM),
                dose_gate_hill=float(gate_hill) if np.isfinite(gate_hill) else float(fit_config.fixed_dose_gate_hill),
            )
            pred_alive_ctrl_end, pred_dead_ctrl_end = simulate_endpoint_metrics(
                fit_config=fit_config,
                dfdctp_signal_curve=curves_by_ploidy[ploidy],
                ploidy_params=ploidy_params[ploidy],
                dose_uM=0.0,
                N0=control_N0,
                D0=control_D0,
                endpoint_days=endpoint_days,
                n_tr=n_tr,
                dose_gate_ec50_uM=float(gate_ec50) if np.isfinite(gate_ec50) else float(fit_config.fixed_dose_gate_ec50_uM),
                dose_gate_hill=float(gate_hill) if np.isfinite(gate_hill) else float(fit_config.fixed_dose_gate_hill),
            )

            observed_alive_end = interpolate_series_at_time(
                dose_data["t"],
                mean_observed_series(dose_data["y_alive"]),
                endpoint_days,
            )
            observed_alive_ctrl_end = interpolate_series_at_time(
                control_data["t"],
                mean_observed_series(control_data["y_alive"]),
                endpoint_days,
            )

            pred_alive_vs_control[ploidy] = pred_alive_end / pred_alive_ctrl_end if pred_alive_ctrl_end > 0 else np.nan
            obs_alive_vs_control[ploidy] = (
                observed_alive_end / observed_alive_ctrl_end
                if np.isfinite(observed_alive_ctrl_end) and observed_alive_ctrl_end > 0 else np.nan
            )

            row.update(
                {
                    f"pred_alive_end_{ploidy}": pred_alive_end,
                    f"pred_dead_end_{ploidy}": pred_dead_end,
                    f"pred_alive_fraction_initial_{ploidy}": pred_alive_end / N0 if N0 > 0 else np.nan,
                    f"pred_dead_fraction_initial_{ploidy}": pred_dead_end / N0 if N0 > 0 else np.nan,
                    f"pred_alive_fraction_control_{ploidy}": pred_alive_vs_control[ploidy],
                    f"obs_alive_fraction_control_{ploidy}": obs_alive_vs_control[ploidy],
                    f"sensitivity_index_{ploidy}": 1.0 - pred_alive_vs_control[ploidy] if np.isfinite(pred_alive_vs_control[ploidy]) else np.nan,
                }
            )

        ratio = (
            pred_alive_vs_control["4N"] / pred_alive_vs_control["2N"]
            if np.isfinite(pred_alive_vs_control["2N"]) and pred_alive_vs_control["2N"] > 0
            and np.isfinite(pred_alive_vs_control["4N"]) and pred_alive_vs_control["4N"] > 0
            else np.nan
        )
        row["ratio_4N_over_2N_pred_alive_fraction_control"] = ratio
        row["log2_ratio_4N_over_2N_pred_alive_fraction_control"] = float(np.log2(ratio)) if np.isfinite(ratio) else np.nan
        row["delta_sensitivity_index_4N_minus_2N"] = (
            row["sensitivity_index_4N"] - row["sensitivity_index_2N"]
            if np.isfinite(row["sensitivity_index_2N"]) and np.isfinite(row["sensitivity_index_4N"])
            else np.nan
        )
        rows.append(row)

    if not rows:
        raise ValueError("No matched nonzero doses were available for dose-wise ploidy comparison")
    return pd.DataFrame(rows).sort_values("dose_uM").reset_index(drop=True)


def plot_dose_response_comparison(comparison_df: pd.DataFrame, output_path: Path) -> None:
    dose_nM = comparison_df["dose_uM"].to_numpy(dtype=float) * 1000.0
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 7.2), sharex=True, gridspec_kw={"height_ratios": [2.2, 1.2]})

    ax = axes[0]
    for ploidy, color in (("2N", "#4C78A8"), ("4N", "#F58518")):
        pred = comparison_df[f"pred_alive_fraction_control_{ploidy}"].to_numpy(dtype=float)
        obs = comparison_df[f"obs_alive_fraction_control_{ploidy}"].to_numpy(dtype=float)
        ax.plot(dose_nM, pred, color=color, marker="o", linewidth=2.0, label=f"{ploidy} model")
        ax.scatter(dose_nM, obs, color=color, marker="x", s=42, linewidths=1.2, label=f"{ploidy} observed")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.6)
    ax.set_xscale("log")
    ax.set_ylabel("Alive at endpoint / matched 0 nM")
    ax.set_title("Dose-wise ploidy sensitivity comparison from fitted cohort model")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    ax.legend(fontsize=8)

    ax = axes[1]
    log2_ratio = comparison_df["log2_ratio_4N_over_2N_pred_alive_fraction_control"].to_numpy(dtype=float)
    delta_sensitivity = comparison_df["delta_sensitivity_index_4N_minus_2N"].to_numpy(dtype=float)
    ax.plot(dose_nM, log2_ratio, color="#7A5195", marker="s", linewidth=2.0, label=r"$\log_2$ (4N / 2N)")
    ax.axhline(0.0, color="black", linewidth=1.0, linestyle="--", alpha=0.6)
    ax2 = ax.twinx()
    ax2.plot(dose_nM, delta_sensitivity, color="#2F9E44", marker="^", linewidth=1.6, label=r"$\Delta$ sensitivity index")
    ax.set_xscale("log")
    ax.set_xlabel("Administered gemcitabine dose (nM)")
    ax.set_ylabel(r"$\log_2$ alive/control ratio")
    ax2.set_ylabel(r"$\Delta(1-\mathrm{alive/control})$ 4N$-$2N")
    ax.grid(True, which="both", linestyle="--", alpha=0.3)
    lines_1, labels_1 = ax.get_legend_handles_labels()
    lines_2, labels_2 = ax2.get_legend_handles_labels()
    ax.legend(lines_1 + lines_2, labels_1 + labels_2, fontsize=8, loc="best")

    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_combined_dfdctp_signal_curves(
    curves_by_ploidy: Dict[str, invitro_fitting.DfdctpSignalSurface],
    output_path: Path,
) -> None:
    available_curves = {
        ploidy: curve
        for ploidy, curve in curves_by_ploidy.items()
        if len(curve.calibration_profiles_by_dose) > 0
    }
    if not available_curves:
        raise ValueError("No calibrated dFdCTP signal curves are available")

    max_time = max(
        float(profile.time_days.max())
        for curve in available_curves.values()
        for profile in curve.calibration_profiles_by_dose.values()
        if len(profile.time_days) > 0
    )
    positive_times = [
        float(time_value)
        for curve in available_curves.values()
        for profile in curve.calibration_profiles_by_dose.values()
        for time_value in np.asarray(profile.time_days, dtype=float)
        if np.isfinite(time_value) and time_value > 0
    ]
    min_positive_time = min(positive_times) if positive_times else 1.0 / 24.0
    log_time_floor = min_positive_time / 2.0
    t_grid = np.linspace(0.0, max(5.0, max_time), 400)
    t_plot_grid = np.maximum(t_grid, log_time_floor)
    all_doses = sorted(
        {
            float(dose)
            for curve in available_curves.values()
            for dose in curve.calibration_doses_uM
        }
    )
    dose_colors = {
        dose: color
        for dose, color in zip(all_doses, plt.cm.viridis(np.linspace(0.18, 0.82, len(all_doses))))
    }

    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for ploidy in ("2N", "4N"):
        if ploidy not in available_curves:
            continue
        curve = available_curves[ploidy]
        style = PLOIDY_STYLES.get(ploidy, {"linestyle": "-", "marker": "o"})
        for dose in curve.calibration_doses_uM:
            dose = float(dose)
            profile = curve.calibration_profiles_by_dose[dose]
            color = dose_colors[dose]
            ax.plot(
                t_plot_grid,
                curve(t_grid, dose),
                color=color,
                linestyle=style["linestyle"],
                linewidth=2.0,
                label=f"{ploidy}, {dose:g} uM",
            )
            time_plot = np.maximum(np.asarray(profile.time_days, dtype=float), log_time_floor)
            ax.scatter(
                time_plot,
                profile.induced_signal_uM_values,
                color=color,
                marker=style["marker"],
                s=34,
                alpha=0.82,
                edgecolor="black",
                linewidth=0.35,
            )

    ax.set_xscale("log")
    ax.set_xlim(log_time_floor, max(5.0, max_time))
    ax.set_xlabel("Time (days, log scale; day 0 shown at plotting floor)", fontsize=10)
    ax.set_ylabel("Baseline-subtracted intracellular dFdCTP (uM)", fontsize=10)
    ax.set_title("Intracellular dFdCTP signal drivers by ploidy", fontsize=11)
    ax.tick_params(axis="both", labelsize=9)
    ax.grid(True, alpha=0.25)
    ax.legend(
        loc="upper right",
        fontsize=7.5,
        frameon=True,
        title="Modeled profile",
        title_fontsize=7.5,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def get_live_dead_dose_items(
    fit_config: invitro_fitting.JointFitConfig,
) -> list[Tuple[str, float]]:
    paths = invitro_fitting.default_experiment_paths()
    df = invitro_fitting.assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    dose_labels = sorted(
        [dose for dose in df["gem"].dropna().unique()],
        key=lambda x: float(str(x).split()[0]),
    )
    dose_items = []
    for dose_label in dose_labels:
        dose_uM = float(str(dose_label).split()[0]) / 1000.0
        if dose_uM > 0:
            dose_items.append((str(dose_label), dose_uM))
    return dose_items


def format_nm_dose_label(dose_uM: float) -> str:
    return f"{float(dose_uM) * 1000.0:g} nM"


def plot_effective_dfdctp_signal_curves(
    best_row: pd.Series,
    fit_config: invitro_fitting.JointFitConfig,
    curves_by_ploidy: Dict[str, invitro_fitting.DfdctpSignalSurface],
    output_path: Path,
) -> None:
    available_curves = {
        ploidy: curve
        for ploidy, curve in curves_by_ploidy.items()
        if len(curve.calibration_profiles_by_dose) > 0
    }
    if not available_curves:
        raise ValueError("No calibrated dFdCTP signal curves are available")

    max_time = max(
        float(profile.time_days.max())
        for curve in available_curves.values()
        for profile in curve.calibration_profiles_by_dose.values()
        if len(profile.time_days) > 0
    )
    positive_times = [
        float(time_value)
        for curve in available_curves.values()
        for profile in curve.calibration_profiles_by_dose.values()
        for time_value in np.asarray(profile.time_days, dtype=float)
        if np.isfinite(time_value) and time_value > 0
    ]
    min_positive_time = min(positive_times) if positive_times else 1.0 / 24.0
    max_plot_time = max(5.0, max_time)
    t_grid = np.geomspace(min_positive_time, max_plot_time, 400)

    plotted_doses_uM = {
        float(dose_uM)
        for _dose_label, dose_uM in get_live_dead_dose_items(fit_config)
        if float(dose_uM) > 0
    }
    plotted_doses_uM.update(
        float(dose)
        for curve in available_curves.values()
        for dose in curve.calibration_doses_uM
        if float(dose) > 0
    )
    dose_items = [
        (format_nm_dose_label(dose_uM), dose_uM)
        for dose_uM in sorted(plotted_doses_uM)
    ]
    dose_colors = {
        dose_label: color
        for (dose_label, _dose_uM), color in zip(
            dose_items,
            plt.cm.viridis(np.linspace(0.12, 0.88, len(dose_items))),
        )
    }

    gate_ec50 = row_float(best_row, "dose_gate_ec50_uM", fit_config.fixed_dose_gate_ec50_uM)
    gate_hill = row_float(best_row, "dose_gate_hill", fit_config.fixed_dose_gate_hill)
    use_hill = bool(
        fit_config.use_hill_dose_gate
        and np.isfinite(gate_ec50)
        and gate_ec50 > 0
        and np.isfinite(gate_hill)
        and gate_hill > 0
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5), sharex=True, sharey=True)
    for ax, ploidy in zip(axes, ("2N", "4N")):
        if ploidy not in available_curves:
            ax.set_axis_off()
            continue
        curve = available_curves[ploidy]
        params = ploidy_parameter_dict(best_row, ploidy, fit_config)
        beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
        reference_dose = invitro_fitting.get_dose_scaling_reference_uM(curve)
        calibrated_doses = np.asarray(curve.calibration_doses_uM, dtype=float)
        for dose_label, dose_uM in dose_items:
            pk_signal = np.asarray(curve(t_grid, dose_uM), dtype=float)
            effective_signal = np.asarray(
                [
                    invitro_fitting.apply_effective_dose_correction(
                        signal_uM=float(signal_value),
                        dose_uM=float(dose_uM),
                        reference_dose_uM=reference_dose,
                        beta_dose=beta_dose,
                        use_hill_dose_gate=use_hill,
                        dose_gate_ec50_uM=gate_ec50 if use_hill else None,
                        dose_gate_hill=gate_hill if use_hill else None,
                    )
                    for signal_value in pk_signal
                ],
                dtype=float,
            )
            is_calibrated_pk_dose = bool(np.any(np.isclose(calibrated_doses, dose_uM, rtol=1e-6, atol=1e-9)))
            ax.plot(
                t_grid,
                effective_signal,
                color=dose_colors[dose_label],
                linestyle="-" if is_calibrated_pk_dose else "--",
                linewidth=1.8,
                alpha=0.95,
                label=dose_label,
            )

        ax.set_xscale("log")
        ax.set_xlim(min_positive_time, max_plot_time)
        ax.set_title(f"{ploidy} effective signal", fontsize=10)
        ax.grid(True, alpha=0.25)
        ax.tick_params(axis="both", labelsize=9)
        ax.text(
            0.03,
            0.96,
            f"beta={beta_dose:.3g}\nref={reference_dose * 1000.0:.3g} nM",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.82, "edgecolor": "0.75"},
        )

    axes[0].set_ylabel("Effective dFdCTP signal (uM-equivalent)", fontsize=10)
    for ax in axes:
        ax.set_xlabel("Time (days, log scale; t=0 omitted)", fontsize=10)
    correction_label = "beta + Hill corrected" if use_hill else "beta corrected"
    fig.suptitle(f"Fitted effective intracellular dFdCTP signal ({correction_label})", fontsize=11)
    handles, labels = axes[-1].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            labels,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            fontsize=8,
            title="Dose",
            title_fontsize=8,
        )
    if use_hill:
        fig.text(
            0.5,
            0.005,
            f"Shared Hill gate: EC50={gate_ec50 * 1000.0:.3g} nM, h={gate_hill:.3g}",
            ha="center",
            fontsize=8,
        )
    style_handles = [
        Line2D([0], [0], color="0.25", linestyle="-", linewidth=1.8, label="PK-calibrated dose"),
        Line2D([0], [0], color="0.25", linestyle="--", linewidth=1.8, label="Interpolated/scaled dose"),
    ]
    fig.legend(
        handles=style_handles,
        loc="lower left",
        bbox_to_anchor=(0.01, 0.005),
        fontsize=8,
        frameon=False,
    )
    fig.text(
        0.5,
        0.025,
        "Effective signal is zero at t=0; t=0 is omitted from the log-scale curves.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0.0, 0.08, 0.88, 0.95))
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def build_cohort_dose_data_list(
    df: pd.DataFrame,
    gem_doses: Iterable[str],
    ploidy: str,
    fit_config: invitro_fitting.JointFitConfig,
) -> list[Dict[str, Any]]:
    dose_data_list = []
    for gem_dose in gem_doses:
        try:
            aligned = invitro_fitting.get_aligned_live_dead_data(
                df,
                gem_dose=gem_dose,
                ploidy=ploidy,
                t_max=fit_config.fit_t_max,
                count_transitional_as_alive=fit_config.count_transitional_as_alive,
            )
        except ValueError:
            continue
        dose_data_list.append(
            {
                "dose_label": gem_dose,
                "dose_muM": float(str(gem_dose).split()[0]) / 1000.0,
                "t": aligned["t"],
                "y_alive": aligned["y_alive"],
                "y_dead": aligned["y_dead"],
                "N0": float(np.nanmean(aligned["y_alive"][0, :])),
                "D0": float(np.nanmean(aligned["y_dead"][0, :])),
            }
        )
    return dose_data_list


def treatment_params_for_plot(
    params: Dict[str, float],
    fit_config: invitro_fitting.JointFitConfig,
) -> list[float]:
    beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
    mu_confluence_death = float(
        params.get(
            "mu_confluence_death",
            fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
        )
    )
    treatment_value_by_name = {
        "k_tr": params["k_tr"],
        "k_kill": params["k_kill"],
        "k_clear": params["k_clear"],
        "k_cyto": params.get("k_cyto", 1e-8),
        "beta_dose": beta_dose,
        "mu_base_death": params["mu_base_death"],
        "mu_confluence_death": mu_confluence_death,
    }
    return [
        treatment_value_by_name[name]
        for name in invitro_fitting.get_treatment_parameter_names_for_config(fit_config)
    ]


def plot_cohort_fit_subplots(
    dose_data_list: list[Dict[str, Any]],
    ploidy_label: str,
    params: Dict[str, float],
    n_tr: int,
    dfdctp_signal_curve: invitro_fitting.DfdctpSignalSurface,
    fit_config: invitro_fitting.JointFitConfig,
    dose_gate_ec50_uM: float,
    dose_gate_hill: float,
    fit_summary: Dict[str, Any],
    output_folder: Path,
) -> None:
    n_doses = len(dose_data_list)
    cols = 2
    rows = max(5, int(np.ceil(n_doses / cols)))
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5.6, rows * 3.25 + 1.0), sharex=True)
    axes = np.asarray(axes).flatten()

    beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
    mu_base_death = float(params["mu_base_death"])
    mu_confluence_death = float(
        params.get(
            "mu_confluence_death",
            fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
        )
    )
    treatment_params = treatment_params_for_plot(params, fit_config)

    for i, data in enumerate(dose_data_list):
        ax = axes[i]
        t_data = data["t"]
        y_alive_raw = data["y_alive"]
        y_dead_raw = data["y_dead"]
        alive_sim, dead_sim = invitro_fitting.simulate_joint_ext(
            t_data,
            treatment_params,
            data["N0"],
            data["D0"],
            float(params["r"]),
            float(params["K"]),
            data["dose_muM"],
            dfdctp_signal_curve,
            n_tr,
            model_variant=fit_config.model_variant,
            fit_config=fit_config,
            beta_dose=beta_dose,
            use_confluence_death=fit_config.use_confluence_death,
            mu_confluence_death=mu_confluence_death,
            confluence_death_exponent=fit_config.confluence_death_exponent,
            use_hill_dose_gate=fit_config.use_hill_dose_gate,
            dose_gate_ec50_uM=dose_gate_ec50_uM,
            dose_gate_hill=dose_gate_hill,
        )

        if y_alive_raw.ndim > 1:
            for j in range(y_alive_raw.shape[1]):
                ax.scatter(t_data, y_alive_raw[:, j], color=ALIVE_OBS_COLOR, alpha=0.35, s=20)
        else:
            ax.scatter(t_data, y_alive_raw, color=ALIVE_OBS_COLOR, alpha=0.6, s=20)
        ax.plot(
            t_data,
            alive_sim,
            color=ALIVE_MODEL_COLOR,
            linewidth=2,
            label="Alive Model" if i == 0 else "",
        )

        if y_dead_raw.ndim > 1:
            for j in range(y_dead_raw.shape[1]):
                ax.scatter(t_data, y_dead_raw[:, j], color=DEAD_OBS_COLOR, alpha=0.35, s=20)
        else:
            ax.scatter(t_data, y_dead_raw, color=DEAD_OBS_COLOR, alpha=0.6, s=20)
        ax.plot(
            t_data,
            dead_sim,
            color=DEAD_MODEL_COLOR,
            linewidth=2,
            linestyle="--",
            label="Dead Model" if i == 0 else "",
        )

        ax.set_title(f"Dose: {data['dose_label']}", fontsize=12, fontweight="bold")
        ax.grid(True, linestyle="--", alpha=0.6)
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        if i % cols == 0:
            ax.set_ylabel("Cell Count", fontsize=11)
        if i >= (rows - 1) * cols:
            ax.set_xlabel("Time (Days)", fontsize=11)
        if i == 0:
            ax.legend(loc="upper left", fontsize=10)

    for j in range(n_doses, len(axes)):
        fig.delaxes(axes[j])

    param_lines = [
        f"Optima: n_tr = {n_tr} | Model: {fit_summary['model_variant']}",
        f"$k_{{tr}}$ = {params['k_tr']:.3f} | $k_{{kill}}$ = {params['k_kill']:.3f} | "
        f"$k_{{clear}}$ = {params['k_clear']:.3f}"
        + (f" | $k_{{cyto}}$ = {params['k_cyto']:.3f}" if fit_config.model_variant != "delayed_death_only" and "k_cyto" in params else "")
        + f" | $\\beta_{{dose}}$ = {beta_dose:.3f}",
        f"Hill EC50 = {dose_gate_ec50_uM * 1000.0:.3f} nM | Hill h = {dose_gate_hill:.3f} | "
        f"$\\mu_{{base death}}$ = {mu_base_death:.4f} day$^{{-1}}$ | "
        f"$\\mu_{{conf death}}$ = {mu_confluence_death:.4f} day$^{{-1}}$ | "
        f"Confluence exponent = {fit_config.confluence_death_exponent:.2f}",
        f"Objective: {fit_summary['objective']} | Channels: {fit_summary['observation_channels']} | "
        f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f} | "
        f"$\\theta_{{dead}}$ = {fit_summary['theta_dead']:.3f} | Data NLL = {fit_summary['nll']:.3f}",
    ]
    fig.suptitle(
        f"{ploidy_label} Cohort - Joint Live/Dead Kinetic Fit\n" + "\n".join(param_lines),
        fontsize=10,
        fontweight="bold",
        y=0.995,
    )
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.90])
    output_path = output_folder / f"cohort_joint_fit_{invitro_fitting.slugify_label(ploidy_label)}.png"
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_single_dose_ploidy_comparison(
    best_row: pd.Series,
    fit_config: invitro_fitting.JointFitConfig,
    curves_by_ploidy: Dict[str, invitro_fitting.DfdctpSignalSurface],
    output_folder: Path,
    dose_label: str = "50 nM",
) -> Path:
    paths = invitro_fitting.default_experiment_paths()
    df = invitro_fitting.assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    n_tr = int(best_row["n_tr"])
    dose_gate_ec50_uM = row_float(best_row, "dose_gate_ec50_uM", fit_config.fixed_dose_gate_ec50_uM)
    dose_gate_hill = row_float(best_row, "dose_gate_hill", fit_config.fixed_dose_gate_hill)

    plot_rows = []
    for ploidy in ("2N", "4N"):
        aligned = invitro_fitting.get_aligned_live_dead_data(
            df,
            gem_dose=dose_label,
            ploidy=ploidy,
            t_max=fit_config.fit_t_max,
            count_transitional_as_alive=fit_config.count_transitional_as_alive,
        )
        params = ploidy_parameter_dict(best_row, ploidy, fit_config)
        beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
        mu_confluence_death = float(
            params.get(
                "mu_confluence_death",
                fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
            )
        )
        dose_uM = float(str(dose_label).split()[0]) / 1000.0
        t_data = aligned["t"]
        alive_sim, dead_sim = invitro_fitting.simulate_joint_ext(
            t_data,
            treatment_params_for_plot(params, fit_config),
            float(np.nanmean(aligned["y_alive"][0, :])),
            float(np.nanmean(aligned["y_dead"][0, :])),
            float(params["r"]),
            float(params["K"]),
            dose_uM,
            curves_by_ploidy[ploidy],
            n_tr,
            model_variant=fit_config.model_variant,
            fit_config=fit_config,
            beta_dose=beta_dose,
            use_confluence_death=fit_config.use_confluence_death,
            mu_confluence_death=mu_confluence_death,
            confluence_death_exponent=fit_config.confluence_death_exponent,
            use_hill_dose_gate=fit_config.use_hill_dose_gate,
            dose_gate_ec50_uM=dose_gate_ec50_uM,
            dose_gate_hill=dose_gate_hill,
        )
        plot_rows.append(
            {
                "ploidy": ploidy,
                "t": t_data,
                "alive_obs": aligned["y_alive"],
                "dead_obs": aligned["y_dead"],
                "alive_sim": alive_sim,
                "dead_sim": dead_sim,
            }
        )

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), sharex=True, sharey=True)
    for ax, row in zip(axes, plot_rows):
        t_data = row["t"]
        for j in range(row["alive_obs"].shape[1]):
            ax.scatter(t_data, row["alive_obs"][:, j], color=ALIVE_OBS_COLOR, alpha=0.35, s=22)
        ax.plot(t_data, row["alive_sim"], color=ALIVE_MODEL_COLOR, linewidth=2.2, label="Alive Model")
        for j in range(row["dead_obs"].shape[1]):
            ax.scatter(t_data, row["dead_obs"][:, j], color=DEAD_OBS_COLOR, alpha=0.35, s=22)
        ax.plot(t_data, row["dead_sim"], color=DEAD_MODEL_COLOR, linewidth=2.2, linestyle="--", label="Dead Model")
        ax.set_title(f"{row['ploidy']} at {dose_label}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Time (Days)")
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.grid(True, linestyle="--", alpha=0.6)
    axes[0].set_ylabel("Cell Count")
    axes[0].legend(loc="upper left", fontsize=9)
    fig.suptitle(f"{dose_label} live/dead fit comparison: 2N vs 4N", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.93])

    output_path = output_folder / f"dose_{invitro_fitting.slugify_label(dose_label)}_2n_vs_4n.png"
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


def plot_cohort_fits_from_summary(
    best_row: pd.Series,
    fit_config: invitro_fitting.JointFitConfig,
    curves_by_ploidy: Dict[str, invitro_fitting.DfdctpSignalSurface],
    output_folder: Path,
) -> None:
    paths = invitro_fitting.default_experiment_paths()
    df = invitro_fitting.assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    gem_doses = sorted(
        [dose for dose in df["gem"].dropna().unique()],
        key=lambda x: float(str(x).split()[0]),
    )
    n_tr = int(best_row["n_tr"])
    dose_gate_ec50_uM = row_float(best_row, "dose_gate_ec50_uM", fit_config.fixed_dose_gate_ec50_uM)
    dose_gate_hill = row_float(best_row, "dose_gate_hill", fit_config.fixed_dose_gate_hill)

    for ploidy in ("2N", "4N"):
        params = ploidy_parameter_dict(best_row, ploidy, fit_config)
        dose_data_list = build_cohort_dose_data_list(df, gem_doses, ploidy, fit_config)
        if not dose_data_list:
            continue

        beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
        mu_base_death = float(params["mu_base_death"])
        mu_confluence_death = float(
            params.get(
                "mu_confluence_death",
                fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
            )
        )
        fit_summary = {
            "objective": best_row["objective"],
            "observation_channels": best_row["observation_channels"],
            "model_variant": best_row["model_variant"],
            "theta_alive": row_float(best_row, "theta_alive"),
            "theta_dead": row_float(best_row, "theta_dead"),
            "beta_dose": beta_dose,
            "mu_base_death": mu_base_death,
            "mu_confluence_death": mu_confluence_death,
            "use_confluence_death": bool(fit_config.use_confluence_death),
            "confluence_death_exponent": float(fit_config.confluence_death_exponent),
            "use_hill_dose_gate": bool(fit_config.use_hill_dose_gate),
            "dose_gate_ec50_uM": dose_gate_ec50_uM,
            "dose_gate_hill": dose_gate_hill,
            "nll": row_float(best_row, "data_nll"),
            "aic": None,
            "bic": None,
            "rmse": np.nan,
        }
        plot_cohort_fit_subplots(
            dose_data_list=dose_data_list,
            ploidy_label=ploidy,
            params=params,
            n_tr=n_tr,
            dfdctp_signal_curve=curves_by_ploidy[ploidy],
            fit_config=fit_config,
            dose_gate_ec50_uM=dose_gate_ec50_uM,
            dose_gate_hill=dose_gate_hill,
            fit_summary=fit_summary,
            output_folder=output_folder,
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot manuscript-oriented summaries from an invitro_fitting output folder."
    )
    parser.add_argument(
        "output_folder",
        type=Path,
        nargs="?",
        default=DEFAULT_FIT_OUTPUT_FOLDER,
        help=(
            "Legacy positional folder containing joint_fit_summary.tsv. "
            "Generated plots now default to a timestamped Results/ run unless --output-dir is supplied."
        ),
    )
    parser.add_argument(
        "--fit-output",
        type=Path,
        help="Folder containing joint_fit_summary.tsv. Overrides the positional folder.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Directory for generated plots/tables. If omitted, outputs are written "
            "to a timestamped Results/in-vitro/pkpd_live_dead_model/runs directory."
        ),
    )
    parser.add_argument(
        "--fit-t-max",
        type=float,
        default=DEFAULT_FIT_T_MAX,
        help=(
            "Maximum day used when reconstructing observed trajectories from raw data. "
            f"Defaults to {DEFAULT_FIT_T_MAX:g}, matching alsoGoodFit_20260514T093906."
        ),
    )
    parser.add_argument(
        "--parameters",
        nargs="+",
        default=DEFAULT_PARAMETERS,
        help="Parameter names to plot, without the 2N_/4N_ prefixes.",
    )
    parser.add_argument(
        "--skip-dfdctp",
        action="store_true",
        help="Do not regenerate dFdCTP signal-curve plots or PK tail diagnostics.",
    )
    parser.add_argument(
        "--skip-cohort",
        action="store_true",
        help="Do not regenerate cohort_joint_fit_2n.png and cohort_joint_fit_4n.png.",
    )
    parser.add_argument(
        "--comparison-dose",
        default="50 nM",
        help="Dose label for the single-dose 2N-vs-4N comparison plot.",
    )
    parser.add_argument(
        "--skip-dose-comparison",
        action="store_true",
        help="Do not regenerate the single-dose 2N-vs-4N comparison plot.",
    )
    return parser.parse_args()


def resolve_output_folder(output_folder: Path) -> Path:
    """Resolve target-repo and legacy miningcloneid-style fit-output paths."""
    output_folder = Path(output_folder)
    if output_folder.exists():
        return output_folder
    parts = output_folder.parts
    if "invitro_fitting_outputs" in parts:
        idx = parts.index("invitro_fitting_outputs")
        suffix = Path(*parts[idx + 1:]) if idx + 1 < len(parts) else Path()
        candidate = FIT_OUTPUTS_ROOT / suffix
        if candidate.exists():
            return candidate
    return output_folder


def main() -> None:
    args = parse_args()
    fit_output_folder = resolve_output_folder(args.fit_output or args.output_folder)
    summary_path = fit_output_folder / "joint_fit_summary.tsv"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing summary file: {summary_path}")

    output_folder = args.output_dir or default_output_folder()
    if not output_folder.is_absolute():
        output_folder = REPO_ROOT / output_folder
    output_folder = output_folder.resolve()
    output_folder.mkdir(parents=True, exist_ok=True)

    summary_df = pd.read_csv(summary_path, sep="\t")
    best_row = choose_best_row(summary_df)
    comparison_df = build_comparison_table(best_row, args.parameters)
    fit_config = config_from_best_row(best_row, fit_t_max=args.fit_t_max)

    table_path = output_folder / "ploidy_parameter_comparison.tsv"
    fold_change_path = output_folder / "ploidy_parameter_log2_fold_change.png"
    paired_path = output_folder / "ploidy_parameter_paired_values.png"
    dfdctp_path = output_folder / "dfdctp_signal_curve_combined_ploidy.png"
    effective_dfdctp_path = output_folder / "effective_dfdctp_signal_curve_combined_ploidy.png"
    dose_response_table_path = output_folder / "dose_response_ploidy_comparison.tsv"
    dose_response_plot_path = output_folder / "dose_response_ploidy_comparison.png"
    cohort_2n_path = output_folder / "cohort_joint_fit_2n.png"
    cohort_4n_path = output_folder / "cohort_joint_fit_4n.png"
    single_dose_comparison_path = (
        output_folder / f"dose_{invitro_fitting.slugify_label(args.comparison_dose)}_2n_vs_4n.png"
    )

    comparison_df.to_csv(table_path, sep="\t", index=False)
    plot_fold_change(comparison_df, fold_change_path)
    plot_paired_values(comparison_df, paired_path)
    curves_by_ploidy = build_dfdctp_surfaces()
    if not args.skip_dfdctp:
        invitro_fitting.save_pk_tail_diagnostics(curves_by_ploidy, output_folder)
        dfdctp_y_limits = invitro_fitting.get_dfdctp_signal_curve_y_limits(curves_by_ploidy)
        for ploidy in ("2N", "4N"):
            invitro_fitting.plot_dfdctp_signal_curve(
                ploidy,
                curves_by_ploidy[ploidy],
                output_dir=output_folder,
                y_limits=dfdctp_y_limits,
            )
            invitro_fitting.plot_dfdctp_amplitude_scaling(ploidy, curves_by_ploidy[ploidy], output_dir=output_folder)
        plot_combined_dfdctp_signal_curves(curves_by_ploidy, dfdctp_path)
        plot_effective_dfdctp_signal_curves(best_row, fit_config, curves_by_ploidy, effective_dfdctp_path)
    dose_response_df = build_dose_response_comparison_table(best_row, fit_config, curves_by_ploidy)
    dose_response_df.to_csv(dose_response_table_path, sep="\t", index=False)
    plot_dose_response_comparison(dose_response_df, dose_response_plot_path)
    if not args.skip_cohort:
        plot_cohort_fits_from_summary(best_row, fit_config, curves_by_ploidy, output_folder)
    if not args.skip_dose_comparison:
        plot_single_dose_ploidy_comparison(
            best_row,
            fit_config,
            curves_by_ploidy,
            output_folder,
            dose_label=args.comparison_dose,
        )

    print(f"Read {summary_path}")
    print(f"Wrote {table_path}")
    print(f"Wrote {fold_change_path}")
    print(f"Wrote {paired_path}")
    if not args.skip_dfdctp:
        print(f"Wrote {output_folder / 'pk_tail_fit_diagnostics.tsv'}")
        print(f"Wrote {output_folder / 'dfdctp_signal_curve_2n.png'}")
        print(f"Wrote {output_folder / 'dfdctp_signal_curve_4n.png'}")
        print(f"Wrote {output_folder / 'dfdctp_amplitude_scaling_2n.png'}")
        print(f"Wrote {output_folder / 'dfdctp_amplitude_scaling_4n.png'}")
        print(f"Wrote {dfdctp_path}")
        print(f"Wrote {effective_dfdctp_path}")
    print(f"Wrote {dose_response_table_path}")
    print(f"Wrote {dose_response_plot_path}")
    if not args.skip_cohort:
        print(f"Wrote {cohort_2n_path}")
        print(f"Wrote {cohort_4n_path}")
    if not args.skip_dose_comparison:
        print(f"Wrote {single_dose_comparison_path}")


if __name__ == "__main__":
    main()
