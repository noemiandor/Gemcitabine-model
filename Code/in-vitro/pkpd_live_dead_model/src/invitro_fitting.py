
import math
import os
import re
import sys
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass, field, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import odeint, solve_ivp
from scipy.optimize import least_squares, minimize
from scipy.special import gammaln
from scipy.stats import linregress

SCRIPT_PATH = Path(__file__).resolve()
MODULE_ROOT = SCRIPT_PATH.parents[1]
PROJECT_ROOT = MODULE_ROOT.parents[2]

# Current modeling convention:
# The preferred live/dead model is driven directly by intracellular dFdCTP.
# Parent Gemcitabine PK is not used to generate the live/dead drug signal.

# PK dFdCTP measurements are reported in ng/mL in the workbook. Dividing by the
# molecular weight converts ng/mL to nmol/mL, which is numerically equal to uM.
DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL = 503.1

# Current project guidance is that the main 2N/4N PK sheets correspond to
# 1.0 uM Gemcitabine reference dosing, while the explicitly named
# low-initial-Gemcitabine sheets correspond to 0.1 uM.
PKPD_REFERENCE_DOSE_BY_SHEET_UM = {
    "2N": 1.0,
    "4N": 1.0,
    "2N_lowInitialGemcitabine": 0.1,
    "4N_lowInitialGemcitabine": 0.1,
}

ALLOWED_MODEL_VARIANTS = {
    "delayed_death_only",
    "immediate_cytostasis_delayed_death",
    "delayed_cytostasis_delayed_death",
}

ALLOWED_MODEL_PRESETS = {
    "minimal",
    "beta_only",
    "hill_only",
    "beta_hill",
    "beta_baseline_confluence",
    "beta_hill_baseline_confluence",
}


def validate_model_variant(model_variant: str) -> str:
    if model_variant not in ALLOWED_MODEL_VARIANTS:
        raise ValueError(
            "model_variant must be one of "
            f"{sorted(ALLOWED_MODEL_VARIANTS)}, got {model_variant!r}"
        )
    return model_variant


BASE_PARAMETER_NAMES = ("r", "K", "k_tr", "k_kill", "k_clear")
DOSE_SCALING_PARAMETER_NAME = "beta_dose"
BASELINE_DEATH_PARAMETER_NAME = "mu_base_death"
CONFLUENCE_DEATH_PARAMETER_NAME = "mu_confluence_death"
OPTIONAL_DYNAMIC_PARAMETER_NAMES = (
    "k_cyto",
    DOSE_SCALING_PARAMETER_NAME,
    BASELINE_DEATH_PARAMETER_NAME,
    CONFLUENCE_DEATH_PARAMETER_NAME,
)
ALL_DYNAMIC_PARAMETER_NAMES = BASE_PARAMETER_NAMES + OPTIONAL_DYNAMIC_PARAMETER_NAMES


def get_parameter_names_for_config(fit_config: "JointFitConfig") -> Tuple[str, ...]:
    model_variant = validate_model_variant(fit_config.model_variant)
    names = list(BASE_PARAMETER_NAMES)
    if model_variant != "delayed_death_only":
        names.append("k_cyto")
    if fit_config.fit_beta_dose:
        names.append(DOSE_SCALING_PARAMETER_NAME)
    names.append(BASELINE_DEATH_PARAMETER_NAME)
    if fit_config.use_confluence_death and fit_config.fit_mu_confluence_death:
        names.append(CONFLUENCE_DEATH_PARAMETER_NAME)
    return tuple(names)


def get_parameter_names(model_variant: str) -> Tuple[str, ...]:
    return get_parameter_names_for_config(resolve_joint_fit_config(JointFitConfig(model_variant=model_variant)))


def get_treatment_parameter_names_for_config(fit_config: "JointFitConfig") -> Tuple[str, ...]:
    parameter_names = get_parameter_names_for_config(fit_config)
    return tuple(name for name in parameter_names if name not in {"r", "K"})


def get_treatment_parameter_names(model_variant: str) -> Tuple[str, ...]:
    return get_treatment_parameter_names_for_config(resolve_joint_fit_config(JointFitConfig(model_variant=model_variant)))


def unpack_treatment_params_for_config(params: Sequence[float], fit_config: "JointFitConfig") -> Dict[str, float]:
    treatment_names = get_treatment_parameter_names_for_config(fit_config)
    params_arr = np.asarray(params, dtype=float)
    if len(params_arr) != len(treatment_names):
        raise ValueError(
            f"Expected {len(treatment_names)} treatment parameters for {fit_config.model_variant}, got {len(params_arr)}"
        )
    out = {name: float(value) for name, value in zip(treatment_names, params_arr)}
    if DOSE_SCALING_PARAMETER_NAME not in out:
        out[DOSE_SCALING_PARAMETER_NAME] = float(fit_config.fixed_beta_dose)
    if CONFLUENCE_DEATH_PARAMETER_NAME not in out:
        out[CONFLUENCE_DEATH_PARAMETER_NAME] = (
            float(fit_config.fixed_mu_confluence_death)
            if fit_config.use_confluence_death else 0.0
        )
    return out


def unpack_treatment_params(params: Sequence[float], model_variant: str) -> Dict[str, float]:
    return unpack_treatment_params_for_config(params, resolve_joint_fit_config(JointFitConfig(model_variant=model_variant)))


@dataclass(frozen=True)
class ExperimentPaths:
    project_root: Path
    exp_base: Path
    counts_raw: Path
    counts_agg: Path
    platemap: Path
    pkpd_constants: Path
    output_dir: Path


@dataclass(frozen=True)
class PKConfig:
    dfdctp_molecular_weight_ng_per_nmol: float = 503.1
    censored_strategy: str = "nan"
    fallback_half_life_days: float = 1.0
    reference_dose_by_sheet_uM: Dict[str, float] = field(default_factory=lambda: dict(PKPD_REFERENCE_DOSE_BY_SHEET_UM))


@dataclass(frozen=True)
class JointFitConfig:
    max_days: float = 5.0
    fit_t_max: float = 5.0
    objective: str = "negative_binomial"
    observation_channels: str = "alive_dead"
    model_variant: str = "immediate_cytostasis_delayed_death"
    model_preset: Optional[str] = None
    count_transitional_as_alive: bool = False
    use_hill_dose_gate: bool = False
    fit_beta_dose: bool = True
    fixed_beta_dose: float = 1.0
    fit_hill_dose_gate: bool = True
    fixed_dose_gate_ec50_uM: float = 0.0125
    fixed_dose_gate_hill: float = 2.0
    use_confluence_death: bool = False
    fit_mu_confluence_death: bool = True
    fixed_mu_confluence_death: float = 0.0
    confluence_death_exponent: float = 4.0
    # For multi-process sweeps, set BLAS threads externally, e.g.
    # OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python invitro_fitting.py
    n_jobs: int = 7
    fit_means_only: bool = False
    high_dose_weight: float = 1.0
    n_tr_values: Tuple[int, ...] = tuple(range(2, 8))
    max_nfev: int = 3000
    optimizer_start_indices: Optional[Tuple[int, ...]] = None
    fit_dose_labels: Optional[Tuple[str, ...]] = None
    max_replicates_per_ploidy_dose: Optional[int] = None
    accept_incomplete_optimizer: bool = False
    large_objective_penalty: float = 1e30
    optimizer_method: str = "L-BFGS-B"
    solver_method: str = "LSODA"
    solver_rtol: float = 1e-6
    solver_atol: float = 1e-8
    prior_sd_log_r: float = 0.75
    prior_sd_log_K: float = 0.75
    prior_sd_log_k_tr: float = 1.00
    prior_sd_log_k_kill: float = 1.00
    prior_sd_log_k_clear: float = 1.00
    prior_sd_log_k_cyto: float = 1.00
    prior_sd_log_beta_dose: float = 0.50
    prior_sd_log_mu_base_death: float = 1.00
    prior_sd_log_mu_confluence_death: float = 1.00
    prior_sd_log_dose_gate_ec50_uM: float = 0.75
    prior_sd_log_dose_gate_hill: float = 0.75
    lower_bounds: Dict[str, float] = field(default_factory=lambda: {
        "r": 1e-8,
        "K": 1e-8,
        "k_tr": 1e-4,
        "k_kill": 1e-4,
        "k_clear": 1e-8,
        "k_cyto": 1e-8,
        "beta_dose": 0.25,
        "mu_base_death": 1e-8,
        "mu_confluence_death": 1e-8,
        "dose_gate_ec50_uM": 0.001,
        "dose_gate_hill": 0.5,
        "theta_alive": 1e-3,
        "theta_dead": 1e-3,
    })
    upper_bounds: Dict[str, float] = field(default_factory=lambda: {
        "r": 100.0,
        "K": 1e9,
        "k_tr": 50.0,
        "k_kill": 500.0,
        "k_clear": 5.0,
        "k_cyto": 1e6,
        "beta_dose": 4.0,
        "mu_base_death": 5.0,
        "mu_confluence_death": 5.0,
        "dose_gate_ec50_uM": 0.1,
        "dose_gate_hill": 8.0,
        "theta_alive": 1e6,
        "theta_dead": 1e6,
    })

    def __post_init__(self):
        validate_model_variant(self.model_variant)
        if self.model_preset is not None and self.model_preset not in ALLOWED_MODEL_PRESETS:
            raise ValueError(
                f"model_preset must be one of {sorted(ALLOWED_MODEL_PRESETS)}, got {self.model_preset!r}"
            )
        if not np.isfinite(self.fixed_beta_dose) or self.fixed_beta_dose <= 0:
            raise ValueError("fixed_beta_dose must be positive and finite")
        if not np.isfinite(self.fixed_mu_confluence_death) or self.fixed_mu_confluence_death < 0:
            raise ValueError("fixed_mu_confluence_death must be finite and >= 0")
        if not np.isfinite(self.confluence_death_exponent) or self.confluence_death_exponent <= 0:
            raise ValueError("confluence_death_exponent must be finite and > 0")
        if not np.isfinite(self.prior_sd_log_mu_confluence_death) or self.prior_sd_log_mu_confluence_death <= 0:
            raise ValueError("prior_sd_log_mu_confluence_death must be finite and > 0")
        if int(self.n_jobs) != self.n_jobs:
            raise ValueError("n_jobs must be an integer")
        if int(self.n_jobs) == 0:
            raise ValueError("n_jobs must be nonzero; use 1 for serial execution or a negative value for all cores")
        if int(self.max_nfev) <= 0:
            raise ValueError("max_nfev must be positive")
        if not self.n_tr_values:
            raise ValueError("n_tr_values must contain at least one value")
        if any(int(value) != value or int(value) < 1 for value in self.n_tr_values):
            raise ValueError("n_tr_values must be positive integers")
        if self.optimizer_start_indices is not None:
            if not self.optimizer_start_indices:
                raise ValueError("optimizer_start_indices must be None or contain at least one index")
            if any(int(value) != value or int(value) < 0 for value in self.optimizer_start_indices):
                raise ValueError("optimizer_start_indices must contain nonnegative integer indices")
        if self.fit_dose_labels is not None and not self.fit_dose_labels:
            raise ValueError("fit_dose_labels must be None or contain at least one dose label")
        if self.max_replicates_per_ploidy_dose is not None:
            if int(self.max_replicates_per_ploidy_dose) != self.max_replicates_per_ploidy_dose:
                raise ValueError("max_replicates_per_ploidy_dose must be an integer")
            if int(self.max_replicates_per_ploidy_dose) <= 0:
                raise ValueError("max_replicates_per_ploidy_dose must be positive")
        if self.use_hill_dose_gate:
            if not np.isfinite(self.fixed_dose_gate_ec50_uM) or self.fixed_dose_gate_ec50_uM <= 0:
                raise ValueError("fixed_dose_gate_ec50_uM must be positive and finite when Hill gate is enabled")
            if not np.isfinite(self.fixed_dose_gate_hill) or self.fixed_dose_gate_hill <= 0:
                raise ValueError("fixed_dose_gate_hill must be positive and finite when Hill gate is enabled")


def joint_fit_config_from_preset(preset: str, **overrides: Any) -> JointFitConfig:
    if preset not in ALLOWED_MODEL_PRESETS:
        raise ValueError(
            f"preset must be one of {sorted(ALLOWED_MODEL_PRESETS)}, got {preset!r}"
        )
    base: Dict[str, Any] = {"model_preset": preset}
    if preset == "minimal":
        base.update(
            use_hill_dose_gate=False,
            fit_beta_dose=False,
            fixed_beta_dose=1.0,
            use_confluence_death=False,
        )
    elif preset == "beta_only":
        base.update(
            use_hill_dose_gate=False,
            fit_beta_dose=True,
            use_confluence_death=False,
        )
    elif preset == "hill_only":
        base.update(
            use_hill_dose_gate=True,
            fit_hill_dose_gate=True,
            fit_beta_dose=False,
            fixed_beta_dose=1.0,
            use_confluence_death=False,
        )
    elif preset == "beta_hill":
        base.update(
            use_hill_dose_gate=True,
            fit_hill_dose_gate=True,
            fit_beta_dose=True,
            use_confluence_death=False,
        )
    elif preset == "beta_baseline_confluence":
        base.update(
            use_hill_dose_gate=False,
            fit_beta_dose=True,
            use_confluence_death=True,
            fit_mu_confluence_death=True,
        )
    elif preset == "beta_hill_baseline_confluence":
        base.update(
            use_hill_dose_gate=True,
            fit_hill_dose_gate=True,
            fit_beta_dose=True,
            use_confluence_death=True,
            fit_mu_confluence_death=True,
        )
    base.update(overrides)
    return JointFitConfig(**base)


def resolve_joint_fit_config(fit_config: Optional[JointFitConfig] = None, **overrides: Any) -> JointFitConfig:
    if fit_config is None:
        if overrides.get("model_preset") is not None:
            preset = overrides.pop("model_preset")
            return joint_fit_config_from_preset(preset, **overrides)
        return JointFitConfig(**overrides) if overrides else JointFitConfig()
    if fit_config.model_preset is not None:
        default_values = asdict(JointFitConfig())
        requested_values = asdict(fit_config)
        preset_overrides = {
            key: value
            for key, value in requested_values.items()
            if key != "model_preset" and value != default_values[key]
        }
        resolved = joint_fit_config_from_preset(
            fit_config.model_preset,
            **preset_overrides,
        )
    else:
        resolved = fit_config
    return replace(resolved, **overrides) if overrides else resolved


def resolve_beta_dose_for_params(
    params: Dict[str, float],
    fit_config: JointFitConfig,
) -> float:
    if fit_config.fit_beta_dose:
        if DOSE_SCALING_PARAMETER_NAME not in params:
            raise KeyError("beta_dose is fitted but missing from parameter dictionary")
        beta = float(params[DOSE_SCALING_PARAMETER_NAME])
    else:
        beta = float(fit_config.fixed_beta_dose)
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError(f"Effective beta_dose must be positive and finite, got {beta}")
    return beta


def resolve_mu_confluence_death_for_params(
    params: Dict[str, float],
    fit_config: JointFitConfig,
) -> float:
    if not fit_config.use_confluence_death:
        return 0.0
    if fit_config.fit_mu_confluence_death:
        if CONFLUENCE_DEATH_PARAMETER_NAME not in params:
            raise KeyError("mu_confluence_death is fitted but missing from parameter dictionary")
        mu_conf = float(params[CONFLUENCE_DEATH_PARAMETER_NAME])
    else:
        mu_conf = float(fit_config.fixed_mu_confluence_death)
    if not np.isfinite(mu_conf) or mu_conf < 0:
        raise ValueError(f"Effective mu_confluence_death must be finite and >= 0, got {mu_conf}")
    return mu_conf


def resolve_hill_dose_gate_params_for_unpacked(
    unpacked: Dict[str, Any],
    fit_config: JointFitConfig,
) -> Dict[str, Any]:
    if not fit_config.use_hill_dose_gate:
        return {
            "use_hill_dose_gate": False,
            "dose_gate_ec50_uM": np.nan,
            "dose_gate_hill": np.nan,
            "fit_hill_dose_gate": bool(fit_config.fit_hill_dose_gate),
        }
    if fit_config.fit_hill_dose_gate:
        dose_gate_params = unpacked.get("dose_gate_params", {})
        ec50 = float(dose_gate_params["dose_gate_ec50_uM"])
        hill = float(dose_gate_params["dose_gate_hill"])
    else:
        ec50 = float(fit_config.fixed_dose_gate_ec50_uM)
        hill = float(fit_config.fixed_dose_gate_hill)
    if not np.isfinite(ec50) or ec50 <= 0:
        raise ValueError(f"dose_gate_ec50_uM must be positive and finite, got {ec50}")
    if not np.isfinite(hill) or hill <= 0:
        raise ValueError(f"dose_gate_hill must be positive and finite, got {hill}")
    return {
        "use_hill_dose_gate": True,
        "dose_gate_ec50_uM": ec50,
        "dose_gate_hill": hill,
        "fit_hill_dose_gate": bool(fit_config.fit_hill_dose_gate),
    }


def resolve_simulation_effect_parameters(
    treatment_params: Dict[str, float],
    fit_config: JointFitConfig,
    beta_dose: Optional[float] = None,
    use_confluence_death: Optional[bool] = None,
    mu_confluence_death: Optional[float] = None,
    confluence_death_exponent: Optional[float] = None,
    use_hill_dose_gate: Optional[bool] = None,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
) -> Dict[str, Any]:
    beta = resolve_beta_dose_for_params(treatment_params, fit_config)
    if beta_dose is not None:
        beta = float(beta_dose)
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError(f"Effective beta_dose must be positive and finite, got {beta}")
    resolved_use_confluence = (
        fit_config.use_confluence_death if use_confluence_death is None else bool(use_confluence_death)
    )
    if mu_confluence_death is None:
        mu_conf = resolve_mu_confluence_death_for_params(treatment_params, fit_config)
    else:
        mu_conf = float(mu_confluence_death)
    if not resolved_use_confluence:
        mu_conf = 0.0
    if not np.isfinite(mu_conf) or mu_conf < 0:
        raise ValueError(f"Effective mu_confluence_death must be finite and >= 0, got {mu_conf}")
    exponent = (
        fit_config.confluence_death_exponent
        if confluence_death_exponent is None else float(confluence_death_exponent)
    )
    if not np.isfinite(exponent) or exponent <= 0:
        raise ValueError(f"confluence_death_exponent must be finite and > 0, got {exponent}")
    resolved_use_hill = fit_config.use_hill_dose_gate if use_hill_dose_gate is None else bool(use_hill_dose_gate)
    if resolved_use_hill:
        gate_ec50 = (
            float(fit_config.fixed_dose_gate_ec50_uM)
            if dose_gate_ec50_uM is None else float(dose_gate_ec50_uM)
        )
        gate_hill = (
            float(fit_config.fixed_dose_gate_hill)
            if dose_gate_hill is None else float(dose_gate_hill)
        )
        if not np.isfinite(gate_ec50) or gate_ec50 <= 0:
            raise ValueError(f"dose_gate_ec50_uM must be positive and finite, got {gate_ec50}")
        if not np.isfinite(gate_hill) or gate_hill <= 0:
            raise ValueError(f"dose_gate_hill must be positive and finite, got {gate_hill}")
    else:
        gate_ec50 = None
        gate_hill = None
    return {
        "beta_dose": beta,
        "use_confluence_death": resolved_use_confluence,
        "mu_confluence_death": mu_conf,
        "confluence_death_exponent": exponent,
        "use_hill_dose_gate": resolved_use_hill,
        "dose_gate_ec50_uM": gate_ec50,
        "dose_gate_hill": gate_hill,
    }


def confluence_death_hazard(
    alive_count: float,
    K: float,
    mu_confluence_death: float,
    confluence_death_exponent: float,
) -> float:
    if mu_confluence_death <= 0 or alive_count <= 0:
        return 0.0
    if not np.isfinite(K) or K <= 0:
        raise ValueError(f"K must be finite and > 0, got {K}")
    if not np.isfinite(confluence_death_exponent) or confluence_death_exponent <= 0:
        raise ValueError(
            "confluence_death_exponent must be finite and > 0, "
            f"got {confluence_death_exponent}"
        )
    density_fraction = np.clip(float(alive_count) / float(K), 0.0, 2.0)
    hazard = float(mu_confluence_death) * density_fraction ** float(confluence_death_exponent)
    if not np.isfinite(hazard) or hazard < 0:
        raise ValueError(
            "Confluence death hazard became invalid for "
            f"A={alive_count}, K={K}, mu={mu_confluence_death}, exponent={confluence_death_exponent}"
        )
    return hazard


def total_drug_independent_death_hazard(
    alive_count: float,
    K: float,
    mu_base_death: float,
    use_confluence_death: bool,
    mu_confluence_death: float = 0.0,
    confluence_death_exponent: float = 4.0,
) -> float:
    if not np.isfinite(mu_base_death) or mu_base_death < 0:
        raise ValueError(f"mu_base_death must be finite and >= 0, got {mu_base_death}")
    hazard = float(mu_base_death)
    if use_confluence_death:
        hazard += confluence_death_hazard(
            alive_count=alive_count,
            K=K,
            mu_confluence_death=mu_confluence_death,
            confluence_death_exponent=confluence_death_exponent,
        )
    if not np.isfinite(hazard) or hazard < 0:
        raise ValueError(
            "Drug-independent death hazard became invalid for "
            f"A={alive_count}, K={K}, mu_base={mu_base_death}, "
            f"use_confluence_death={use_confluence_death}, mu_confluence={mu_confluence_death}"
        )
    return hazard


@dataclass
class PKTailFitDiagnostics:
    sheet_name: str
    ploidy: str
    dose_uM: float
    analyte: str
    n_timepoints_total: int
    n_timepoints_positive_finite: int
    peak_idx: int
    peak_time_days: float
    peak_signal_uM: float
    tail_n_points: int
    tail_time_days: Tuple[float, ...]
    tail_signal_uM: Tuple[float, ...]
    tail_fit_used: bool
    fallback_used: bool
    tail_decay_per_day: float
    tail_half_life_days: float
    tail_fit_r2: Optional[float]
    warning_message: Optional[str]


@dataclass
class SimulationResult:
    alive: np.ndarray
    dead: np.ndarray
    success: bool
    message: str = ""


@dataclass
class ReplicateTrajectory:
    ploidy: str
    dose_label: str
    dose_uM: float
    replicate_id: Tuple[str, int]
    t: np.ndarray
    alive: np.ndarray
    dead: np.ndarray
    N0: float
    D0: float


def trim_finite_live_dead_observations(
    t: Sequence[float],
    alive: Sequence[float],
    dead: Sequence[float],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """
    Keep only timepoints with finite paired alive/dead observations.

    The joint count-likelihood fit requires paired finite observations for both
    channels at a given timepoint. Missing trailing observations are dropped
    here rather than allowed to propagate into the likelihood as NaNs.
    """
    t_arr = np.asarray(t, dtype=float)
    alive_arr = np.asarray(alive, dtype=float)
    dead_arr = np.asarray(dead, dtype=float)
    if t_arr.shape != alive_arr.shape or t_arr.shape != dead_arr.shape:
        raise ValueError("t, alive, and dead must have matching shapes")
    valid_mask = np.isfinite(t_arr) & np.isfinite(alive_arr) & np.isfinite(dead_arr)
    if not np.any(valid_mask):
        raise ValueError("No finite paired alive/dead observations remain")
    first_valid_idx = int(np.flatnonzero(valid_mask)[0])
    if first_valid_idx != 0:
        raise ValueError("Replicate is missing a finite baseline alive/dead observation at the first timepoint")
    dropped = int(np.count_nonzero(~valid_mask))
    return t_arr[valid_mask], alive_arr[valid_mask], dead_arr[valid_mask], dropped


def default_experiment_paths(project_root: Optional[Path] = None) -> ExperimentPaths:
    root = (project_root or PROJECT_ROOT).resolve()
    exp_base = root / "Data" / "in-vitro" / "pkpd_live_dead_model"
    count_base = root / "Data" / "GemDelayKillTerm" / "processed"
    output_dir = (
        root
        / "Results"
        / "in-vitro"
        / "pkpd_live_dead_model"
        / "runs"
        / f"{datetime.now().strftime('%Y%m%dT%H%M%S')}_pkpd_fit"
    )
    return ExperimentPaths(
        project_root=root,
        exp_base=exp_base,
        counts_raw=count_base / "counts_by_well_time.parquet",
        counts_agg=count_base / "counts_by_well_time_wellAggregated.parquet",
        platemap=exp_base / "raw" / "Gemcitabine_PlateMap_20240111.xlsx",
        pkpd_constants=exp_base / "raw" / "drugKinetics" / "GemcitabineExposure_PKPD.xlsx",
        output_dir=output_dir,
    )


def validate_paths(paths: ExperimentPaths, required: Sequence[str]) -> None:
    missing = []
    for attr_name in required:
        if not hasattr(paths, attr_name):
            raise ValueError(f"Unknown ExperimentPaths attribute '{attr_name}'")
        path_value = getattr(paths, attr_name)
        if not Path(path_value).exists():
            missing.append(f"{attr_name}={path_value}")
    if missing:
        raise FileNotFoundError("Missing required experiment paths: " + ", ".join(missing))


_DEFAULT_PATHS = default_experiment_paths()
PATHS = {
    "counts_raw": _DEFAULT_PATHS.counts_raw,
    "counts_agg": _DEFAULT_PATHS.counts_agg,
    "platemap": _DEFAULT_PATHS.platemap,
    "pkpd_constants": _DEFAULT_PATHS.pkpd_constants,
}


def slugify_label(value: str) -> str:
    """Converts labels to filesystem-safe filename fragments."""
    return re.sub(r'[^A-Za-z0-9]+', '_', value.strip()).strip('_').lower()


def normalize_column_name(name: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "_", str(name).strip().lower()).strip("_")


def parse_well_id(well: Any) -> Tuple[str, int]:
    match = re.fullmatch(r"\s*([A-Ha-h])\s*0*([1-9]|1[0-2])\s*", str(well))
    if match is None:
        raise ValueError(f"Unsupported well identifier '{well}'")
    return match.group(1).upper(), int(match.group(2))


def parse_ploidy_label(value: Any) -> Optional[str]:
    if pd.isna(value):
        return None
    match = re.search(r"(?i)\b([24]N)\b", str(value))
    if match is None:
        return None
    return match.group(1).upper()


def get_pk_reference_dose_uM(sheet_name: str, pk_config: Optional[PKConfig] = None) -> float:
    """Returns the documented administered Gemcitabine dose for a PK workbook sheet."""
    config = pk_config or PKConfig()
    if sheet_name not in config.reference_dose_by_sheet_uM:
        raise KeyError(
            f"No documented PK reference dose configured for sheet '{sheet_name}'. "
            "Update PKConfig.reference_dose_by_sheet_uM to match the workbook metadata."
        )
    return config.reference_dose_by_sheet_uM[sheet_name]


@dataclass
class DfdctpDoseProfile:
    """
    One calibrated dose-specific intracellular dFdCTP profile in uM.

    Each PK sheet keeps its own timing, peak value, and terminal tail. The
    higher-level surface interpolates between these dose-specific profiles in
    log-dose space.
    """
    sheet_name: str
    dose_uM: float
    analyte_column: str
    time_days: np.ndarray
    replicate_count_values: np.ndarray
    raw_signal_uM_values: np.ndarray
    raw_signal_uM_sd_values: np.ndarray
    induced_signal_uM_values: np.ndarray
    induced_signal_uM_sd_values: np.ndarray
    peak_signal_uM: float
    peak_time_days: float
    tail_half_life_days: float
    tail_decay_per_day: float
    tail_fit_r2: Optional[float] = None
    warning_message: Optional[str] = None
    tail_fit_diagnostics: Optional[PKTailFitDiagnostics] = None

    def evaluate(self, t_days):
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)
        if len(self.time_days) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        else:
            values = np.interp(
                t_eval,
                self.time_days,
                self.induced_signal_uM_values,
                left=float(self.induced_signal_uM_values[0]),
                right=np.nan,
            )
            tail_mask = t_eval > self.time_days[-1]
            if np.any(tail_mask):
                last_value = float(self.induced_signal_uM_values[-1])
                values[tail_mask] = last_value * np.exp(
                    -self.tail_decay_per_day * (t_eval[tail_mask] - float(self.time_days[-1]))
                )
            values = np.clip(values, 0.0, np.inf)
        if scalar_input:
            return float(values[0])
        return values


@dataclass
class DfdctpSignalSurface:
    """
    Dose-dependent intracellular dFdCTP surface C_dFdCTP(t, dose; ploidy).

    For each calibrated dose, the corresponding PK sheet keeps its own temporal
    profile. Evaluation proceeds by:
    1. evaluating each calibrated dose profile at time t, including its own tail
    2. interpolating the resulting dFdCTP concentrations across administered dose
       in log-dose space

    Outside the calibrated dose range:
    - below the minimum calibrated dose, the minimum-dose profile timing/shape is
      preserved and its amplitude is scaled proportionally with dose
    - above the maximum calibrated dose, the current default is explicit linear
      extrapolation in log-dose space, with warnings reported in summaries and
      plots
    """
    source_ploidy: str
    analyte_column: str
    calibration_profiles_by_dose: Dict[float, DfdctpDoseProfile]
    calibration_sheet_names: Tuple[str, ...]
    extrapolation_policy: str = (
        "within-range log-dose interpolation; below-range proportional scaling "
        "from minimum calibrated profile; above-range log-dose linear extrapolation"
    )
    warning_messages: Tuple[str, ...] = ()

    @property
    def calibration_doses_uM(self) -> np.ndarray:
        return np.array(sorted(self.calibration_profiles_by_dose.keys()), dtype=float)

    @property
    def min_calibration_dose_uM(self) -> float:
        doses = self.calibration_doses_uM
        return float(np.min(doses)) if len(doses) > 0 else np.nan

    @property
    def max_calibration_dose_uM(self) -> float:
        doses = self.calibration_doses_uM
        return float(np.max(doses)) if len(doses) > 0 else np.nan

    def get_matching_calibrated_dose(
        self,
        dose_uM: float,
        rtol: float = 1e-9,
        atol: float = 1e-12,
    ) -> Optional[float]:
        for calibrated_dose in self.calibration_doses_uM:
            if np.isclose(float(dose_uM), float(calibrated_dose), rtol=rtol, atol=atol):
                return float(calibrated_dose)
        return None

    def profile_for_dose(self, dose_uM: float) -> Optional[DfdctpDoseProfile]:
        matched_dose = self.get_matching_calibrated_dose(dose_uM)
        if matched_dose is None:
            return None
        return self.calibration_profiles_by_dose[matched_dose]

    def calibration_status_for_doses(self, dose_uM_values: Sequence[float]) -> str:
        dose_arr = np.asarray(dose_uM_values, dtype=float)
        if len(dose_arr) == 0 or len(self.calibration_profiles_by_dose) == 0:
            return "no live/dead doses provided"
        min_dose = self.min_calibration_dose_uM
        max_dose = self.max_calibration_dose_uM
        below = np.any(dose_arr < min_dose)
        above = np.any(dose_arr > max_dose)
        if below and above:
            return (
                "live/dead doses extend below and above the PK calibration range; "
                "below-range doses use proportional scaling of the minimum-dose profile"
            )
        if below:
            return (
                "live/dead doses extend below the PK calibration range; "
                "below-range doses use proportional scaling of the minimum-dose profile"
            )
        if above:
            return "live/dead doses extend above the PK calibration range"
        return "live/dead doses lie within the PK calibration range"

    def evaluate_at_calibrated_doses(self, t_days) -> Tuple[np.ndarray, np.ndarray]:
        doses = self.calibration_doses_uM
        values = np.vstack(
            [
                np.atleast_1d(
                    np.asarray(self.calibration_profiles_by_dose[float(dose)].evaluate(t_days), dtype=float)
                )
                for dose in doses
            ]
        )
        return doses, values

    def evaluate_below_range_dose(self, t_days, dose_muM):
        if dose_muM < 0:
            raise ValueError(f"dose_muM must be >= 0, got {dose_muM}")
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)
        if dose_muM == 0 or len(self.calibration_profiles_by_dose) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        else:
            min_dose = self.min_calibration_dose_uM
            min_profile = self.calibration_profiles_by_dose[min_dose]
            values = min_profile.evaluate(t_eval) * (float(dose_muM) / float(min_dose))
        if scalar_input:
            return float(np.asarray(values, dtype=float).reshape(-1)[0])
        return np.asarray(values, dtype=float)

    def __call__(self, t_days, dose_muM):
        if dose_muM < 0:
            raise ValueError(f"dose_muM must be >= 0, got {dose_muM}")
        doses = self.calibration_doses_uM
        t_arr = np.asarray(t_days, dtype=float)
        scalar_input = np.ndim(t_arr) == 0
        t_eval = np.atleast_1d(t_arr)

        if len(doses) == 0:
            values = np.zeros_like(t_eval, dtype=float)
        elif dose_muM == 0:
            values = np.zeros_like(t_eval, dtype=float)
        elif len(doses) == 1:
            only_dose = float(doses[-1])
            if dose_muM < only_dose:
                values = self.evaluate_below_range_dose(t_eval, dose_muM)
            else:
                values = self.calibration_profiles_by_dose[only_dose].evaluate(t_eval)
        else:
            matched_dose = self.get_matching_calibrated_dose(dose_muM)
            if matched_dose is not None:
                values = self.calibration_profiles_by_dose[matched_dose].evaluate(t_eval)
            elif dose_muM < self.min_calibration_dose_uM:
                values = self.evaluate_below_range_dose(t_eval, dose_muM)
            else:
                dose_grid, profile_values = self.evaluate_at_calibrated_doses(t_eval)
                log_doses = np.log(dose_grid)
                target_log_dose = np.log(max(float(dose_muM), 1e-12))
                if target_log_dose > log_doses[-1]:
                    idx0, idx1 = len(log_doses) - 2, len(log_doses) - 1
                else:
                    idx1 = int(np.searchsorted(log_doses, target_log_dose, side="right"))
                    idx0 = idx1 - 1
                x0, x1 = float(log_doses[idx0]), float(log_doses[idx1])
                y0 = profile_values[idx0, :]
                y1 = profile_values[idx1, :]
                if x1 == x0:
                    values = y0.copy()
                else:
                    frac = (target_log_dose - x0) / (x1 - x0)
                    values = y0 + frac * (y1 - y0)
                values = np.clip(values, 0.0, np.inf)

        if scalar_input:
            return float(np.asarray(values, dtype=float).reshape(-1)[0])
        return np.asarray(values, dtype=float)


def get_dose_scaling_reference_uM(dfdctp_signal_curve: DfdctpSignalSurface) -> float:
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    reference_dose = float(dfdctp_signal_curve.min_calibration_dose_uM)
    if not np.isfinite(reference_dose) or reference_dose <= 0:
        return 1.0
    return reference_dose


def apply_dose_power_correction(
    signal_uM: float,
    dose_uM: float,
    reference_dose_uM: float,
    beta_dose: float,
) -> float:
    signal = float(signal_uM)
    dose = float(dose_uM)
    ref = float(reference_dose_uM)
    beta = float(beta_dose)
    if dose <= 0 or signal <= 0:
        return 0.0
    if not np.isfinite(ref) or ref <= 0:
        raise ValueError(f"reference_dose_uM must be positive and finite, got {reference_dose_uM}")
    if not np.isfinite(beta) or beta <= 0:
        raise ValueError(f"beta_dose must be positive and finite, got {beta_dose}")
    corrected = signal * (dose / ref) ** (beta - 1.0)
    if not np.isfinite(corrected) or corrected <= 0:
        return 0.0
    return float(max(corrected, 0.0))


def hill_dose_gate(
    dose_uM: float,
    ec50_uM: float,
    hill_coefficient: float,
) -> float:
    dose = float(dose_uM)
    ec50 = float(ec50_uM)
    hill = float(hill_coefficient)
    if dose <= 0:
        return 0.0
    if not np.isfinite(ec50) or ec50 <= 0:
        raise ValueError(f"ec50_uM must be positive and finite, got {ec50_uM}")
    if not np.isfinite(hill) or hill <= 0:
        raise ValueError(f"hill_coefficient must be positive and finite, got {hill_coefficient}")
    log_ratio = hill * (np.log(dose) - np.log(ec50))
    if log_ratio >= 0:
        return float(1.0 / (1.0 + np.exp(-log_ratio)))
    exp_lr = np.exp(log_ratio)
    return float(exp_lr / (1.0 + exp_lr))


def normalized_hill_dose_gate(
    dose_uM: float,
    reference_dose_uM: float,
    ec50_uM: float,
    hill_coefficient: float,
) -> float:
    dose = float(dose_uM)
    ref = float(reference_dose_uM)
    if dose <= 0:
        return 0.0
    if not np.isfinite(ref) or ref <= 0:
        raise ValueError(f"reference_dose_uM must be positive and finite, got {reference_dose_uM}")
    gate = hill_dose_gate(dose, ec50_uM, hill_coefficient)
    ref_gate = hill_dose_gate(ref, ec50_uM, hill_coefficient)
    if not np.isfinite(ref_gate) or ref_gate <= 0:
        raise ValueError(f"reference Hill gate must be positive and finite, got {ref_gate}")
    out = gate / ref_gate
    if not np.isfinite(out) or out < 0:
        raise ValueError(f"normalized Hill gate must be finite and nonnegative, got {out}")
    return float(out)


def apply_effective_dose_correction(
    signal_uM: float,
    dose_uM: float,
    reference_dose_uM: float,
    beta_dose: float,
    use_hill_dose_gate: bool = False,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
) -> float:
    corrected = apply_dose_power_correction(
        signal_uM=signal_uM,
        dose_uM=dose_uM,
        reference_dose_uM=reference_dose_uM,
        beta_dose=beta_dose,
    )
    if corrected <= 0 or dose_uM <= 0:
        return 0.0
    if not use_hill_dose_gate:
        return corrected
    if dose_gate_ec50_uM is None or dose_gate_hill is None:
        raise ValueError("dose_gate_ec50_uM and dose_gate_hill are required when Hill gate is enabled")
    gate = normalized_hill_dose_gate(
        dose_uM=dose_uM,
        reference_dose_uM=reference_dose_uM,
        ec50_uM=dose_gate_ec50_uM,
        hill_coefficient=dose_gate_hill,
    )
    out = corrected * gate
    if not np.isfinite(out) or out <= 0:
        return 0.0
    return float(out)


def assert_preferred_live_dead_driver(driver: Any) -> DfdctpSignalSurface:
    """
    Guards the preferred live/dead modeling path.

    The preferred path must be driven directly by a `DfdctpSignalSurface`, not
    by parent Gemcitabine PK or any legacy extracellular-exposure object.
    """
    if not isinstance(driver, DfdctpSignalSurface):
        raise TypeError(
            "Preferred live/dead model driver must be a DfdctpSignalSurface; "
            f"got {type(driver).__name__}"
        )
    return driver


#################### eta, K_decay fitting utilities ####################

def parse_pk_concentration_value(
    value,
    censored_strategy: str = "nan",
) -> Tuple[float, bool, float]:
    """
    Parses PK concentration values while preserving censoring semantics.

    Returns `(numeric_value, was_censored, censor_upper_bound)`. Missing spreadsheet values
    (`np.nan`, blank cells) return `(np.nan, False)`. Censored strings such as
    `BDL`, `BQL`, `BLQ`, `N/F`, and `<0.5` return `was_censored=True`.

    For strings like `<0.5`, `censor_upper_bound` preserves the numeric upper
    bound in the original assay units. Token-only censored values retain
    `censor_upper_bound=np.nan` unless assay metadata provide a numeric LLOQ.
    """
    if censored_strategy not in {"nan", "zero", "half_lod"}:
        raise ValueError(
            "censored_strategy must be one of {'nan', 'zero', 'half_lod'}, "
            f"got {censored_strategy}"
        )

    if pd.isna(value):
        return np.nan, False, np.nan

    if isinstance(value, (int, float, np.integer, np.floating)):
        return float(value), False, np.nan

    value_str = str(value).strip()
    if value_str == "":
        return np.nan, False, np.nan

    try:
        return float(value_str), False, np.nan
    except ValueError:
        pass

    normalized = value_str.upper()
    censored_tokens = {
        "BDL",
        "BQL",
        "BLQ",
        "N/F",
        "ND",
        "N/D",
        "<LLOQ",
        "<LOQ",
        "BELOW DETECTION",
        "BELOW QUANTIFICATION",
    }
    missing_tokens = {"NA", "N/A", "NONE", "NULL"}

    if normalized in missing_tokens:
        return np.nan, False, np.nan

    threshold_match = re.search(r'<\s*([0-9]+(?:\.[0-9]+)?)', value_str)
    if threshold_match:
        threshold = float(threshold_match.group(1))
        if censored_strategy == "half_lod":
            return threshold / 2.0, True, threshold
        if censored_strategy == "zero":
            return 0.0, True, threshold
        return np.nan, True, threshold

    if normalized in censored_tokens:
        if censored_strategy == "zero":
            return 0.0, True, np.nan
        return np.nan, True, np.nan

    return np.nan, False, np.nan

def import_and_clean_pkpd(
    censored_strategy: Optional[str] = None,
    paths: Optional[ExperimentPaths] = None,
    pk_config: Optional[PKConfig] = None,
):
    """
    Imports PK sheets and parses analyte columns with explicit censored-value handling.

    `censored_strategy="nan"` is the conservative default: BDL/BQL/BLQ-style
    entries become `np.nan` and are ignored by downstream means/fits. Set
    `"zero"` to reproduce the prior zero-imputation behavior explicitly, or
    `"half_lod"` to impute half the reported threshold for `<x` values. The
    original censoring metadata and any numeric censor upper bounds are retained
    in companion columns for downstream PK fitting.
    """
    paths = paths or default_experiment_paths()
    pk_config = pk_config or PKConfig()
    validate_paths(paths, required=("pkpd_constants",))
    strategy = censored_strategy or pk_config.censored_strategy

    raw_sheets = pd.read_excel(paths.pkpd_constants, sheet_name=None)
    clean_data_dict = {}
    
    for sheet_name, df in raw_sheets.items():
        df.columns = [re.sub(r'\s+', ' ', str(c).strip()) for c in df.columns]
        
        if 'Timepoint' not in df.columns and 'Sample' in df.columns:
            time_regex, rep_regex = r'-(\d+\.?\d*)h', r'h-(\d+)'
            df['Timepoint'] = df['Sample'].apply(lambda x: float(re.search(time_regex, str(x)).group(1)) if re.search(time_regex, str(x)) else 0.0)
            df['replicate'] = df['Sample'].apply(lambda x: int(re.search(rep_regex, str(x)).group(1)) if re.search(rep_regex, str(x)) else 1)
        
        analyte_cols = [c for c in df.columns if '(ng/mL)' in c]
        for col in analyte_cols:
            parsed_values = df[col].apply(
                lambda value: parse_pk_concentration_value(
                    value,
                    censored_strategy=strategy,
                )
            )
            df[col] = parsed_values.apply(lambda item: item[0]).astype(float)
            df[f"{col}__was_censored"] = parsed_values.apply(lambda item: item[1]).astype(bool)
            df[f"{col}__censor_upper_bound"] = parsed_values.apply(lambda item: item[2]).astype(float)
            
        clean_data_dict[sheet_name] = df
    return clean_data_dict

def identify_pk_analyte_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Finds key PK analyte columns by normalized name."""
    analyte_map: Dict[str, str] = {}
    for col in df.columns:
        normalized = str(col).strip().lower()
        if '(ng/ml)' not in normalized:
            continue
        if 'gemcitabine' in normalized and 'gemcitabine' not in analyte_map:
            analyte_map['gemcitabine'] = col
        elif 'dfdu' in normalized and 'dfdu' not in analyte_map:
            analyte_map['dfdu'] = col
        elif 'dfdctp' in normalized and 'dfdctp' not in analyte_map:
            analyte_map['dfdctp'] = col
        elif 'dfdcdp' in normalized and 'dfdcdp' not in analyte_map:
            analyte_map['dfdcdp'] = col
        elif 'dfdcmp' in normalized and 'dfdcmp' not in analyte_map:
            analyte_map['dfdcmp'] = col
    return analyte_map

def fit_exponential_pk_decay_model(time_days, concentrations) -> Optional[Dict[str, float]]:
    """
    Fits C(t) = C0 * exp(-k * t) on positive, non-missing concentrations.

    The fit uses the decay phase beginning at the observed peak concentration.
    """
    time_arr = np.asarray(time_days, dtype=float)
    conc_arr = np.asarray(concentrations, dtype=float)
    valid_mask = np.isfinite(time_arr) & np.isfinite(conc_arr) & (conc_arr > 0)
    if np.count_nonzero(valid_mask) < 2:
        return None

    time_valid = time_arr[valid_mask]
    conc_valid = conc_arr[valid_mask]
    peak_idx = int(np.argmax(conc_valid))
    tail_time = time_valid[peak_idx:]
    tail_conc = conc_valid[peak_idx:]
    if len(tail_time) < 2:
        return None

    slope, intercept, r_val, _, _ = linregress(tail_time, np.log(tail_conc))
    k_ext_decay_per_day = -slope
    if k_ext_decay_per_day <= 0:
        return None

    c0 = float(np.exp(intercept))
    return {
        "C0": c0,
        "k_ext_decay_per_day": float(k_ext_decay_per_day),
        "half_life_days": float(np.log(2.0) / k_ext_decay_per_day),
        "r2": float(r_val ** 2),
        "tail_time_days": tuple(float(x) for x in tail_time),
        "tail_conc": tuple(float(x) for x in tail_conc),
        "tail_n_points": int(len(tail_time)),
        "peak_idx": int(peak_idx),
    }

def print_pk_workbook_summary(pk_sheets: Dict[str, pd.DataFrame]) -> None:
    print("PK workbook:")
    for sheet_name, df in pk_sheets.items():
        analyte_columns = identify_pk_analyte_columns(df)
        analyte_names = []
        for key in ["gemcitabine", "dfdu", "dfdctp", "dfdcdp", "dfdcmp"]:
            if key in analyte_columns:
                analyte_names.append(analyte_columns[key].replace(" (ng/mL)", ""))
        analyte_summary = ", ".join(analyte_names) if analyte_names else "no recognized analytes"
        print(f"  {sheet_name} sheet: found {analyte_summary}")

def print_dfdctp_signal_curve_summary(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    live_dead_dose_uM_values: Optional[Sequence[float]] = None,
) -> None:
    sheet_dose_summary = ", ".join(
        f"{sheet}={get_pk_reference_dose_uM(sheet):.4f} uM" for sheet in curve.calibration_sheet_names
    )
    print(f"{ploidy_label} intracellular dFdCTP signal:")
    print(f"  analyte driver: {curve.analyte_column}")
    print(f"  PK sheets used: {', '.join(curve.calibration_sheet_names)}")
    print(f"  PK sheet reference doses: {sheet_dose_summary}")
    print(
        "  calibration dose range: "
        f"{curve.min_calibration_dose_uM:.4f} to {curve.max_calibration_dose_uM:.4f} uM"
    )
    print(f"  extrapolation policy: {curve.extrapolation_policy}")
    for dose in curve.calibration_doses_uM:
        profile = curve.calibration_profiles_by_dose[float(dose)]
        print(
            f"  dose {dose:.4f} uM -> "
            f"peak {profile.peak_signal_uM:.6f} uM at {profile.peak_time_days:.4f} d; "
            f"tail t1/2 {profile.tail_half_life_days:.4f} d"
            + (
                f"; tail R^2 {profile.tail_fit_r2:.4f}"
                if profile.tail_fit_r2 is not None else
                "; tail R^2 n/a"
            )
        )
        if profile.warning_message:
            print(f"    warning: {profile.warning_message}")
    for warning_message in curve.warning_messages:
        print(f"  warning: {warning_message}")
    if live_dead_dose_uM_values:
        dose_status = curve.calibration_status_for_doses(live_dead_dose_uM_values)
        print(f"  live/dead dose coverage: {dose_status}")
        dose_arr = np.asarray(live_dead_dose_uM_values, dtype=float)
        if np.any(dose_arr < curve.min_calibration_dose_uM):
            print(
                "  below-range policy: doses below "
                f"{curve.min_calibration_dose_uM:.4f} uM use proportional scaling "
                f"of the {curve.min_calibration_dose_uM:.4f} uM dFdCTP profile "
                f"({np.min(dose_arr):.6f} to {curve.min_calibration_dose_uM:.4f} uM)"
            )
        if np.any(dose_arr > curve.max_calibration_dose_uM):
            print(
                "  above-range policy: doses above the PK calibration range use "
                "log-dose linear extrapolation "
                f"({curve.max_calibration_dose_uM:.4f} to {np.max(dose_arr):.6f} uM)"
            )


def collect_pk_tail_diagnostics(
    curves_by_ploidy: Dict[str, DfdctpSignalSurface]
) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for curve in curves_by_ploidy.values():
        for dose, profile in curve.calibration_profiles_by_dose.items():
            diagnostics = profile.tail_fit_diagnostics
            if diagnostics is None:
                continue
            rows.append(
                {
                    "sheet_name": diagnostics.sheet_name,
                    "ploidy": diagnostics.ploidy,
                    "dose_uM": diagnostics.dose_uM,
                    "analyte": diagnostics.analyte,
                    "n_timepoints_total": diagnostics.n_timepoints_total,
                    "n_timepoints_positive_finite": diagnostics.n_timepoints_positive_finite,
                    "peak_idx": diagnostics.peak_idx,
                    "peak_time_days": diagnostics.peak_time_days,
                    "peak_signal_uM": diagnostics.peak_signal_uM,
                    "tail_n_points": diagnostics.tail_n_points,
                    "tail_time_days": diagnostics.tail_time_days,
                    "tail_signal_uM": diagnostics.tail_signal_uM,
                    "tail_fit_used": diagnostics.tail_fit_used,
                    "fallback_used": diagnostics.fallback_used,
                    "tail_decay_per_day": diagnostics.tail_decay_per_day,
                    "tail_half_life_days": diagnostics.tail_half_life_days,
                    "tail_fit_r2": diagnostics.tail_fit_r2,
                    "warning_message": diagnostics.warning_message,
                }
            )
    return pd.DataFrame(rows)


def save_pk_tail_diagnostics(
    curves_by_ploidy: Dict[str, DfdctpSignalSurface],
    output_dir: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    diagnostics_df = collect_pk_tail_diagnostics(curves_by_ploidy)
    output_path = output_dir / "pk_tail_fit_diagnostics.tsv"
    diagnostics_df.to_csv(output_path, sep="\t", index=False)
    return output_path


def get_dfdctp_signal_curve_y_limits(
    curves_by_ploidy: Dict[str, DfdctpSignalSurface],
) -> Tuple[float, float]:
    y_max_values = []
    for curve in curves_by_ploidy.values():
        if len(curve.calibration_profiles_by_dose) == 0:
            continue
        max_time = max(
            float(profile.time_days.max())
            for profile in curve.calibration_profiles_by_dose.values()
            if len(profile.time_days) > 0
        )
        t_grid = np.linspace(0.0, max(5.0, max_time), 400)
        for dose in curve.calibration_doses_uM:
            modeled_signal = np.asarray(curve(t_grid, float(dose)), dtype=float)
            if np.any(np.isfinite(modeled_signal)):
                y_max_values.append(float(np.nanmax(modeled_signal)))
            profile = curve.calibration_profiles_by_dose[float(dose)]
            measured_upper = np.asarray(profile.induced_signal_uM_values, dtype=float) + np.where(
                np.isfinite(profile.induced_signal_uM_sd_values),
                np.asarray(profile.induced_signal_uM_sd_values, dtype=float),
                0.0,
            )
            if np.any(np.isfinite(measured_upper)):
                y_max_values.append(float(np.nanmax(measured_upper)))

    y_max = max(y_max_values) if y_max_values else 1.0
    if not np.isfinite(y_max) or y_max <= 0:
        y_max = 1.0
    return 0.0, y_max * 1.08


def plot_dfdctp_signal_curve(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
    y_limits: Optional[Tuple[float, float]] = None,
):
    """
    Plots measured and modeled baseline-subtracted intracellular dFdCTP signal.
    """
    if len(curve.calibration_profiles_by_dose) == 0:
        return
    max_time = max(float(profile.time_days.max()) for profile in curve.calibration_profiles_by_dose.values() if len(profile.time_days) > 0)
    positive_times = [
        float(time_value)
        for profile in curve.calibration_profiles_by_dose.values()
        for time_value in np.asarray(profile.time_days, dtype=float)
        if np.isfinite(time_value) and time_value > 0
    ]
    min_positive_time = min(positive_times) if positive_times else 1.0 / 24.0
    log_time_floor = min_positive_time / 2.0
    t_grid = np.linspace(0.0, max(5.0, max_time), 400)
    t_plot_grid = np.maximum(t_grid, log_time_floor)
    color_values = plt.cm.tab10(np.linspace(0.0, 1.0, max(len(curve.calibration_doses_uM), 1)))
    dose_colors = {
        float(dose): color
        for dose, color in zip(curve.calibration_doses_uM, color_values)
    }

    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for dose in curve.calibration_doses_uM:
        dose = float(dose)
        dose_label = f"{dose * 1000.0:g} nM"
        profile = curve.calibration_profiles_by_dose[float(dose)]
        modeled_signal = curve(t_grid, float(dose))
        color = dose_colors[dose]
        ax.plot(
            t_plot_grid,
            modeled_signal,
            color=color,
            linewidth=2.0,
            label=f"{dose_label} driver curve",
        )
        if len(profile.time_days) > 0:
            time_plot = np.maximum(np.asarray(profile.time_days, dtype=float), log_time_floor)
            yerr = np.where(
                np.isfinite(profile.induced_signal_uM_sd_values),
                profile.induced_signal_uM_sd_values,
                0.0,
            )
            ax.errorbar(
                time_plot,
                profile.induced_signal_uM_values,
                yerr=yerr,
                fmt="o",
                markersize=5.5,
                capsize=3.0,
                elinewidth=1.0,
                alpha=0.95,
                color=color,
                markeredgecolor="black",
                markeredgewidth=0.45,
                label=f"{dose_label} measured mean +/- SD ({profile.sheet_name})",
            )
    ax.set_xscale("log")
    ax.set_xlim(log_time_floor, max(5.0, max_time))
    if y_limits is not None:
        ax.set_ylim(*y_limits)
    ax.set_xlabel("Time (Days, log scale; day 0 shown at plotting floor)")
    ax.set_ylabel("Baseline-subtracted intracellular dFdCTP (uM)")
    ax.set_title(f"{ploidy_label} dFdCTP Signal Driver")
    ax.text(
        0.02,
        0.98,
        "\n".join([
            f"calibration range: {curve.min_calibration_dose_uM:.3f}-{curve.max_calibration_dose_uM:.3f} uM",
            f"policy: {curve.extrapolation_policy}",
        ]),
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right", fontsize=8)
    plt.tight_layout()
    if output_dir is not None:
        output_path = output_dir / f"dfdctp_signal_curve_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches="tight")
    elif close_fig:
        plt.close(fig)

def plot_dfdctp_amplitude_scaling(
    ploidy_label: str,
    curve: DfdctpSignalSurface,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
):
    """
    Plots calibrated PK-sheet peak dFdCTP signal by administered dose.
    """
    calibration_doses = curve.calibration_doses_uM
    calibration_peaks = np.array(
        [curve.calibration_profiles_by_dose[float(dose)].peak_signal_uM for dose in calibration_doses],
        dtype=float,
    )
    if len(calibration_doses) == 0:
        return

    fig, ax = plt.subplots(figsize=(6.5, 4.4))
    ax.scatter(
        calibration_doses,
        calibration_peaks,
        color="crimson",
        edgecolor="black",
        linewidth=0.4,
        s=60,
        zorder=3,
        label="PK-sheet peak dFdCTP",
    )
    ax.set_xlabel("Administered Gemcitabine Dose (uM)")
    ax.set_ylabel("Peak dFdCTP Signal (uM)")
    ax.set_title(f"{ploidy_label} dFdCTP Calibration Peaks")
    ax.text(
        0.02,
        0.98,
        curve.extrapolation_policy,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85, "edgecolor": "0.7"},
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=8)
    plt.tight_layout()
    if output_dir is not None:
        output_path = output_dir / f"dfdctp_amplitude_scaling_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches="tight")
    elif close_fig:
        plt.close(fig)

def extract_dfdctp_signal_summary(
    df: pd.DataFrame,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Returns per-timepoint dFdCTP replicate counts, means, and SDs in uM."""
    if "Timepoint" not in df.columns or analyte not in df.columns:
        return (
            np.array([], dtype=float),
            np.array([], dtype=int),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )

    data = df[["Timepoint", analyte]].copy()
    data[analyte] = pd.to_numeric(data[analyte], errors="coerce")
    data = data.dropna(subset=["Timepoint", analyte])
    if data.empty:
        return (
            np.array([], dtype=float),
            np.array([], dtype=int),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )

    grouped = data.groupby("Timepoint", sort=True)[analyte]
    summary = grouped.agg(["count", "mean", "std"]).dropna(subset=["mean"])
    if summary.empty:
        return (
            np.array([], dtype=float),
            np.array([], dtype=int),
            np.array([], dtype=float),
            np.array([], dtype=float),
        )

    time_days = summary.index.to_numpy(dtype=float) / 24.0
    mean_uM_values = dfdctp_ng_per_ml_to_uM(
        summary["mean"].to_numpy(dtype=float),
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    sd_uM_values = dfdctp_ng_per_ml_to_uM(
        summary["std"].to_numpy(dtype=float),
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    return (
        time_days,
        summary["count"].to_numpy(dtype=int),
        np.asarray(mean_uM_values, dtype=float),
        np.asarray(sd_uM_values, dtype=float),
    )

def baseline_subtract_treatment_induced_signal(signal_uM) -> np.ndarray:
    """
    Baseline-subtracts a dFdCTP signal profile so the curve represents the
    treatment-induced active metabolite above the 0h level.
    """
    signal_arr = np.asarray(signal_uM, dtype=float)
    if signal_arr.ndim != 1:
        raise ValueError("signal_uM must be a 1D array")
    if len(signal_arr) == 0:
        return signal_arr.copy()
    baseline_value = float(signal_arr[0]) if np.isfinite(signal_arr[0]) else 0.0
    return np.clip(signal_arr - baseline_value, 0.0, np.inf)

def build_dfdctp_profile_from_sheet(
    df: pd.DataFrame,
    sheet_name: str,
    reference_dose_uM: float,
    fallback_half_life_days: float = 1.0,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> Optional[DfdctpDoseProfile]:
    """
    Builds a mean baseline-subtracted dFdCTP profile from one PK sheet.
    """
    time_days, replicate_counts, raw_signal_uM, raw_signal_uM_sd = extract_dfdctp_signal_summary(
        df,
        analyte=analyte,
        molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
    )
    if len(time_days) == 0:
        return None

    induced_signal_uM = baseline_subtract_treatment_induced_signal(raw_signal_uM)
    induced_signal_uM_sd = np.asarray(raw_signal_uM_sd, dtype=float)
    if not np.any(np.isfinite(induced_signal_uM)):
        return None

    peak_idx = int(np.nanargmax(induced_signal_uM))
    peak_signal_uM = float(induced_signal_uM[peak_idx])
    time_of_peak_days = float(time_days[peak_idx])
    tail_fit = fit_exponential_pk_decay_model(time_days, induced_signal_uM)
    if tail_fit is not None:
        tail_decay_per_day = float(tail_fit["k_ext_decay_per_day"])
        tail_half_life_days = float(tail_fit["half_life_days"])
        fit_r2 = float(tail_fit["r2"])
        warning_message = None
        tail_fit_diagnostics = PKTailFitDiagnostics(
            sheet_name=sheet_name,
            ploidy=sheet_name.split("_")[0],
            dose_uM=float(reference_dose_uM),
            analyte=analyte,
            n_timepoints_total=int(len(time_days)),
            n_timepoints_positive_finite=int(np.count_nonzero(np.isfinite(induced_signal_uM) & (induced_signal_uM > 0))),
            peak_idx=int(peak_idx),
            peak_time_days=time_of_peak_days,
            peak_signal_uM=peak_signal_uM,
            tail_n_points=int(tail_fit["tail_n_points"]),
            tail_time_days=tuple(float(x) for x in tail_fit["tail_time_days"]),
            tail_signal_uM=tuple(float(x) for x in tail_fit["tail_conc"]),
            tail_fit_used=True,
            fallback_used=False,
            tail_decay_per_day=tail_decay_per_day,
            tail_half_life_days=tail_half_life_days,
            tail_fit_r2=fit_r2,
            warning_message=None,
        )
    else:
        tail_half_life_days = float(fallback_half_life_days)
        tail_decay_per_day = float(decay_rate_from_half_life_days(tail_half_life_days))
        fit_r2 = None
        warning_message = (
            f"Using fallback dFdCTP tail half-life of {tail_half_life_days:.3f} days "
            f"for {sheet_name} because the terminal fit was not identifiable."
        )
        positive_mask = np.isfinite(induced_signal_uM) & (induced_signal_uM > 0)
        tail_fit_diagnostics = PKTailFitDiagnostics(
            sheet_name=sheet_name,
            ploidy=sheet_name.split("_")[0],
            dose_uM=float(reference_dose_uM),
            analyte=analyte,
            n_timepoints_total=int(len(time_days)),
            n_timepoints_positive_finite=int(np.count_nonzero(positive_mask)),
            peak_idx=int(peak_idx),
            peak_time_days=time_of_peak_days,
            peak_signal_uM=peak_signal_uM,
            tail_n_points=int(np.count_nonzero(positive_mask[peak_idx:])),
            tail_time_days=tuple(float(x) for x in time_days[peak_idx:][positive_mask[peak_idx:]]),
            tail_signal_uM=tuple(float(x) for x in induced_signal_uM[peak_idx:][positive_mask[peak_idx:]]),
            tail_fit_used=False,
            fallback_used=True,
            tail_decay_per_day=tail_decay_per_day,
            tail_half_life_days=tail_half_life_days,
            tail_fit_r2=None,
            warning_message=warning_message,
        )

    return DfdctpDoseProfile(
        sheet_name=sheet_name,
        dose_uM=float(reference_dose_uM),
        analyte_column=analyte,
        time_days=np.asarray(time_days, dtype=float),
        replicate_count_values=np.asarray(replicate_counts, dtype=int),
        raw_signal_uM_values=np.asarray(raw_signal_uM, dtype=float),
        raw_signal_uM_sd_values=np.asarray(raw_signal_uM_sd, dtype=float),
        induced_signal_uM_values=np.asarray(induced_signal_uM, dtype=float),
        induced_signal_uM_sd_values=np.asarray(induced_signal_uM_sd, dtype=float),
        peak_signal_uM=peak_signal_uM,
        peak_time_days=time_of_peak_days,
        tail_half_life_days=tail_half_life_days,
        tail_decay_per_day=tail_decay_per_day,
        tail_fit_r2=fit_r2,
        warning_message=warning_message,
        tail_fit_diagnostics=tail_fit_diagnostics,
    )

def build_dfdctp_signal_curve_for_ploidy(
    pk_sheets: Dict[str, pd.DataFrame],
    ploidy_label: str,
    fallback_half_life_days: float = 1.0,
    analyte: str = "dFdCTP (ng/mL)",
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> DfdctpSignalSurface:
    """
    Builds a ploidy-specific dose-dependent dFdCTP surface from all usable PK
    sheets for that ploidy.
    """
    sheet_candidates = [ploidy_label, f"{ploidy_label}_lowInitialGemcitabine"]
    profiles_by_dose: Dict[float, DfdctpDoseProfile] = {}
    sheet_names_used: List[str] = []
    warning_messages: List[str] = []

    for sheet_name in sheet_candidates:
        sheet_df = pk_sheets.get(sheet_name)
        if sheet_df is None or analyte not in sheet_df.columns:
            continue
        reference_dose = get_pk_reference_dose_uM(sheet_name)
        profile = build_dfdctp_profile_from_sheet(
            sheet_df,
            sheet_name=sheet_name,
            reference_dose_uM=reference_dose,
            fallback_half_life_days=fallback_half_life_days,
            analyte=analyte,
            molecular_weight_ng_per_nmol=molecular_weight_ng_per_nmol,
        )
        if profile is None:
            continue
        profiles_by_dose[reference_dose] = profile
        sheet_names_used.append(sheet_name)
        if profile.warning_message:
            warning_messages.append(profile.warning_message)

    if len(profiles_by_dose) == 0:
        raise ValueError(f"No usable {analyte} PK profile found for ploidy {ploidy_label}")

    return DfdctpSignalSurface(
        source_ploidy=ploidy_label,
        analyte_column=analyte,
        calibration_profiles_by_dose=dict(sorted(profiles_by_dose.items())),
        calibration_sheet_names=tuple(sheet_names_used),
        extrapolation_policy="log-dose linear extrapolation",
        warning_messages=tuple(warning_messages),
    )


def build_preferred_dfdctp_signal_surfaces(
    pk_sheets: Dict[str, pd.DataFrame],
    ploidy_keys: Sequence[str] = ("2N", "4N"),
    fallback_half_life_days: float = 1.0,
) -> Dict[str, DfdctpSignalSurface]:
    """
    Builds the preferred live/dead model drivers from dFdCTP PK only.

    Parent Gemcitabine PK is intentionally not used here. If `Gemcitabine
    (ng/mL)` is absent but dFdCTP is present, the preferred driver still builds.
    """
    return {
        ploidy_key: build_dfdctp_signal_curve_for_ploidy(
            pk_sheets,
            ploidy_label=ploidy_key,
            fallback_half_life_days=fallback_half_life_days,
        )
        for ploidy_key in ploidy_keys
    }

def single_clone_dfdctp_ode_joint(
    y,
    t,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve: Callable[[Union[float, np.ndarray], float], Union[float, np.ndarray]],
    k_tr,
    k_kill,
    k_clear,
    n_tr,
    model_variant="immediate_cytostasis_delayed_death",
    k_cyto=None,
    beta_dose=1.0,
    mu_base_death=0.0,
    use_confluence_death: bool = False,
    mu_confluence_death: float = 0.0,
    confluence_death_exponent: float = 4.0,
    use_hill_dose_gate: bool = False,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
):
    """
    ODE system driven directly by intracellular dFdCTP signal.

    State order:
    - A: Alive cells
    - Z_1 ... Z_n: transit compartments smoothing/delaying dFdCTP signal
    - D_obs: observable dead cells
    """
    model_variant = validate_model_variant(model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    if not np.isfinite(mu_base_death) or mu_base_death < 0:
        raise ValueError(f"mu_base_death must be finite and >= 0, got {mu_base_death}")
    A = y[0]
    D_obs = y[-1]
    c_raw = float(np.asarray(dfdctp_signal_curve(t, dose_muM), dtype=float).reshape(-1)[0])
    c_now = apply_effective_dose_correction(
        signal_uM=c_raw,
        dose_uM=dose_muM,
        reference_dose_uM=get_dose_scaling_reference_uM(dfdctp_signal_curve),
        beta_dose=beta_dose,
        use_hill_dose_gate=use_hill_dose_gate,
        dose_gate_ec50_uM=dose_gate_ec50_uM,
        dose_gate_hill=dose_gate_hill,
    )

    dZ = np.zeros(n_tr)
    if n_tr > 0:
        delayed_signal = y[1:1+n_tr]
        dZ[0] = k_tr * (c_now - delayed_signal[0])
        for i in range(1, n_tr):
            dZ[i] = k_tr * (delayed_signal[i - 1] - delayed_signal[i])
        c_delayed = delayed_signal[-1]
    else:
        c_delayed = c_now

    kappa_drug = k_kill * c_delayed
    baseline_hazard = total_drug_independent_death_hazard(
        alive_count=A,
        K=K,
        mu_base_death=mu_base_death,
        use_confluence_death=use_confluence_death,
        mu_confluence_death=mu_confluence_death,
        confluence_death_exponent=confluence_death_exponent,
    )
    total_death_hazard = baseline_hazard + kappa_drug
    if model_variant == "delayed_death_only":
        growth_multiplier = 1.0
    else:
        if k_cyto is None or not np.isfinite(k_cyto) or k_cyto < 0:
            raise ValueError(
                f"k_cyto must be finite and >= 0 for model_variant={model_variant}, got {k_cyto}"
            )
        cyto_signal = c_now if model_variant == "immediate_cytostasis_delayed_death" else c_delayed
        growth_multiplier = cytostasis_multiplier(cyto_signal, k_cyto)

    dA = r * growth_multiplier * A * (1 - A / K) - total_death_hazard * A
    dD_obs = total_death_hazard * A - k_clear * D_obs
    return [dA] + dZ.tolist() + [dD_obs]


def simulate_joint_dfdctp(
    t,
    params,
    N0,
    D0,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve,
    n_tr,
    model_variant=JointFitConfig().model_variant,
    fit_config: Optional[JointFitConfig] = None,
    beta_dose: Optional[float] = None,
    use_confluence_death: Optional[bool] = None,
    mu_confluence_death: Optional[float] = None,
    confluence_death_exponent: Optional[float] = None,
    use_hill_dose_gate: Optional[bool] = None,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
):
    """
    Simulates the live/dead ODE with dFdCTP supplied directly as the drug-driver
    curve instead of integrating a separate intracellular drug state.
    """
    model_variant = validate_model_variant(model_variant)
    fit_config = resolve_joint_fit_config(fit_config, model_variant=model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    treatment_params = unpack_treatment_params_for_config(params, fit_config)
    effect_params = resolve_simulation_effect_parameters(
        treatment_params=treatment_params,
        fit_config=fit_config,
        beta_dose=beta_dose,
        use_confluence_death=use_confluence_death,
        mu_confluence_death=mu_confluence_death,
        confluence_death_exponent=confluence_death_exponent,
        use_hill_dose_gate=use_hill_dose_gate,
        dose_gate_ec50_uM=dose_gate_ec50_uM,
        dose_gate_hill=dose_gate_hill,
    )
    y0 = [N0] + [0.0] * n_tr + [D0]
    sol = odeint(
        single_clone_dfdctp_ode_joint,
        y0,
        t,
        args=(
            r,
            K,
            dose_muM,
            dfdctp_signal_curve,
            treatment_params["k_tr"],
            treatment_params["k_kill"],
            treatment_params["k_clear"],
            n_tr,
            model_variant,
            treatment_params.get("k_cyto"),
            effect_params["beta_dose"],
            treatment_params[BASELINE_DEATH_PARAMETER_NAME],
            effect_params["use_confluence_death"],
            effect_params["mu_confluence_death"],
            effect_params["confluence_death_exponent"],
            effect_params["use_hill_dose_gate"],
            effect_params["dose_gate_ec50_uM"],
            effect_params["dose_gate_hill"],
        ),
    )
    return sol[:, 0], sol[:, -1]


def simulate_joint_dfdctp_safe(
    t: np.ndarray,
    params: Sequence[float],
    N0: float,
    D0: float,
    r: float,
    K: float,
    dose_muM: float,
    dfdctp_signal_curve: DfdctpSignalSurface,
    n_tr: int,
    fit_config: JointFitConfig,
    model_variant: Optional[str] = None,
    beta_dose: Optional[float] = None,
    use_confluence_death: Optional[bool] = None,
    mu_confluence_death: Optional[float] = None,
    confluence_death_exponent: Optional[float] = None,
    use_hill_dose_gate: Optional[bool] = None,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
) -> SimulationResult:
    model_variant = model_variant or fit_config.model_variant
    model_variant = validate_model_variant(model_variant)
    fit_config = resolve_joint_fit_config(fit_config, model_variant=model_variant)
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    t_arr = np.asarray(t, dtype=float)
    if t_arr.ndim != 1 or len(t_arr) == 0:
        return SimulationResult(np.array([]), np.array([]), False, "Empty or non-1D time grid")
    if not np.all(np.isfinite(t_arr)):
        return SimulationResult(np.array([]), np.array([]), False, "Non-finite time grid")
    if np.any(np.diff(t_arr) < 0):
        return SimulationResult(np.array([]), np.array([]), False, "Time grid must be sorted ascending")

    params_arr = np.asarray(params, dtype=float)
    expected_params = len(get_treatment_parameter_names_for_config(fit_config))
    if len(params_arr) != expected_params or not np.all(np.isfinite(params_arr)):
        return SimulationResult(np.array([]), np.array([]), False, "Invalid treatment params")
    if not np.all(np.isfinite([N0, D0, r, K, dose_muM])) or not np.isfinite(n_tr):
        return SimulationResult(np.array([]), np.array([]), False, "Non-finite simulation inputs")
    if N0 < 0 or D0 < 0 or r < 0 or K <= 0 or dose_muM < 0 or int(n_tr) < 0:
        return SimulationResult(np.array([]), np.array([]), False, "Inputs outside valid ranges")
    try:
        treatment_params = unpack_treatment_params_for_config(params_arr, fit_config)
    except ValueError as exc:
        return SimulationResult(np.array([]), np.array([]), False, str(exc))
    try:
        effect_params = resolve_simulation_effect_parameters(
            treatment_params=treatment_params,
            fit_config=fit_config,
            beta_dose=beta_dose,
            use_confluence_death=use_confluence_death,
            mu_confluence_death=mu_confluence_death,
            confluence_death_exponent=confluence_death_exponent,
            use_hill_dose_gate=use_hill_dose_gate,
            dose_gate_ec50_uM=dose_gate_ec50_uM,
            dose_gate_hill=dose_gate_hill,
        )
    except (KeyError, ValueError) as exc:
        return SimulationResult(np.array([]), np.array([]), False, str(exc))

    bounds_to_check = {
        "k_tr": treatment_params["k_tr"],
        "k_kill": treatment_params["k_kill"],
        "k_clear": treatment_params["k_clear"],
        "beta_dose": effect_params["beta_dose"],
        "r": float(r),
        "K": float(K),
    }
    if "k_cyto" in treatment_params:
        bounds_to_check["k_cyto"] = treatment_params["k_cyto"]
    if treatment_params[BASELINE_DEATH_PARAMETER_NAME] > 0:
        bounds_to_check[BASELINE_DEATH_PARAMETER_NAME] = treatment_params[BASELINE_DEATH_PARAMETER_NAME]
    if effect_params["use_confluence_death"] and effect_params["mu_confluence_death"] > 0:
        bounds_to_check[CONFLUENCE_DEATH_PARAMETER_NAME] = effect_params["mu_confluence_death"]
    if effect_params["use_hill_dose_gate"]:
        bounds_to_check["dose_gate_ec50_uM"] = effect_params["dose_gate_ec50_uM"]
        bounds_to_check["dose_gate_hill"] = effect_params["dose_gate_hill"]
    for name, value in bounds_to_check.items():
        if name == CONFLUENCE_DEATH_PARAMETER_NAME and not fit_config.use_confluence_death:
            continue
        lower = fit_config.lower_bounds.get(name, -np.inf)
        upper = fit_config.upper_bounds.get(name, np.inf)
        if value < (lower - 1e-12) or value > (upper + 1e-12):
            return SimulationResult(np.array([]), np.array([]), False, f"{name}={value} outside bounds")
    if treatment_params[BASELINE_DEATH_PARAMETER_NAME] < 0:
        return SimulationResult(np.array([]), np.array([]), False, "mu_base_death must be >= 0")
    if effect_params["use_confluence_death"]:
        if effect_params["mu_confluence_death"] < 0:
            return SimulationResult(np.array([]), np.array([]), False, "mu_confluence_death must be >= 0")
        if not np.isfinite(effect_params["confluence_death_exponent"]) or effect_params["confluence_death_exponent"] <= 0:
            return SimulationResult(np.array([]), np.array([]), False, "confluence_death_exponent must be > 0")

    y0 = np.array([N0] + [0.0] * int(n_tr) + [D0], dtype=float)

    def rhs(t_eval, state):
        return np.asarray(
            single_clone_dfdctp_ode_joint(
                state,
                t_eval,
                r,
                K,
                dose_muM,
                dfdctp_signal_curve,
                treatment_params["k_tr"],
                treatment_params["k_kill"],
                treatment_params["k_clear"],
                int(n_tr),
                model_variant,
                treatment_params.get("k_cyto"),
                effect_params["beta_dose"],
                treatment_params[BASELINE_DEATH_PARAMETER_NAME],
                effect_params["use_confluence_death"],
                effect_params["mu_confluence_death"],
                effect_params["confluence_death_exponent"],
                effect_params["use_hill_dose_gate"],
                effect_params["dose_gate_ec50_uM"] if effect_params["use_hill_dose_gate"] else None,
                effect_params["dose_gate_hill"] if effect_params["use_hill_dose_gate"] else None,
            ),
            dtype=float,
        )

    try:
        result = solve_ivp(
            rhs,
            t_span=(float(t_arr[0]), float(t_arr[-1])),
            y0=y0,
            t_eval=t_arr,
            method=fit_config.solver_method,
            rtol=fit_config.solver_rtol,
            atol=fit_config.solver_atol,
        )
    except Exception as exc:
        return SimulationResult(np.array([]), np.array([]), False, f"solve_ivp exception: {exc}")

    if not result.success:
        return SimulationResult(np.array([]), np.array([]), False, result.message)
    if result.y.shape[1] != len(t_arr):
        return SimulationResult(np.array([]), np.array([]), False, "Solver returned unexpected shape")

    alive = np.asarray(result.y[0], dtype=float)
    dead = np.asarray(result.y[-1], dtype=float)
    if not np.all(np.isfinite(alive)) or not np.all(np.isfinite(dead)):
        return SimulationResult(alive, dead, False, "Non-finite predictions")
    if np.any(alive < -1e-8) or np.any(dead < -1e-8):
        return SimulationResult(alive, dead, False, "Materially negative predictions")
    alive = np.clip(alive, 0.0, np.inf)
    dead = np.clip(dead, 0.0, np.inf)
    return SimulationResult(alive, dead, True, result.message)


def safe_objective_value(
    objective_fn: Callable[[np.ndarray], float],
    x: np.ndarray,
    penalty: float,
) -> float:
    try:
        value = objective_fn(x)
    except Exception:
        return penalty
    if not np.isfinite(value):
        return penalty
    return float(value)


def get_partial_pooling_coordinate_names(
    parameter_names: Sequence[str],
    ploidies: Sequence[str],
    objective: str,
    fit_config: Optional[JointFitConfig] = None,
) -> List[str]:
    fit_config = resolve_joint_fit_config(fit_config)
    coordinate_names: List[str] = [f"mu_log_{name}" for name in parameter_names]
    for ploidy in ploidies:
        coordinate_names.extend(f"delta_log_{ploidy}_{name}" for name in parameter_names)
    if fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate:
        coordinate_names.extend(["log_dose_gate_ec50_uM", "log_dose_gate_hill"])
    if objective == "negative_binomial":
        coordinate_names.extend(["log_theta_alive", "log_theta_dead"])
    return coordinate_names


def coordinate_bound_flags(
    x: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    coordinate_names: Sequence[str],
    atol: float = 1e-10,
) -> Dict[str, Any]:
    x_arr = np.asarray(x, dtype=float)
    lower_arr = np.asarray(lower_bounds, dtype=float)
    upper_arr = np.asarray(upper_bounds, dtype=float)
    at_lower = np.isclose(x_arr, lower_arr, atol=atol, rtol=0.0)
    at_upper = np.isclose(x_arr, upper_arr, atol=atol, rtol=0.0)
    return {
        "any_at_lower_bound": bool(np.any(at_lower)),
        "any_at_upper_bound": bool(np.any(at_upper)),
        "lower_bound_coordinates": ",".join(
            coordinate_names[idx] for idx in np.flatnonzero(at_lower)
        ),
        "upper_bound_coordinates": ",".join(
            coordinate_names[idx] for idx in np.flatnonzero(at_upper)
        ),
    }


def flatten_partial_pooling_parameters(unpacked: Dict[str, Any]) -> Dict[str, float]:
    flat: Dict[str, float] = {}
    for key, value in unpacked.get("mu_logs", {}).items():
        flat[f"mu_{key}"] = float(np.exp(value))
    for ploidy, params in unpacked.get("ploidy_params", {}).items():
        for key, value in params.items():
            flat[f"{ploidy}_{key}"] = float(value)
    dose_gate_params = unpacked.get("dose_gate_params", {})
    if dose_gate_params:
        flat["dose_gate_ec50_uM"] = float(dose_gate_params.get("dose_gate_ec50_uM", np.nan))
        flat["dose_gate_hill"] = float(dose_gate_params.get("dose_gate_hill", np.nan))
    if unpacked.get("theta_alive") is not None:
        flat["theta_alive"] = float(unpacked["theta_alive"])
    if unpacked.get("theta_dead") is not None:
        flat["theta_dead"] = float(unpacked["theta_dead"])
    return flat


def diagnose_objective_near_start(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    unpack_fn: Callable[[np.ndarray], Dict[str, Any]],
    coordinate_names: Sequence[str],
    label: str,
    output_path: Path,
    penalty: float,
    debug_eval_fn: Optional[Callable[[np.ndarray], Dict[str, Any]]] = None,
    eps_values: Sequence[float] = (1e-6, 1e-4, 1e-2),
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    lower_arr = np.asarray(lower_bounds, dtype=float)
    upper_arr = np.asarray(upper_bounds, dtype=float)
    base_objective = safe_objective_value(objective_fn, x0, penalty)
    base_flat = flatten_partial_pooling_parameters(unpack_fn(x0))
    rows: List[Dict[str, Any]] = []
    for coord_idx, coord_name in enumerate(coordinate_names):
        for eps in eps_values:
            plus = x0.copy()
            minus = x0.copy()
            plus[coord_idx] = min(upper_arr[coord_idx], plus[coord_idx] + eps)
            minus[coord_idx] = max(lower_arr[coord_idx], minus[coord_idx] - eps)
            plus_obj = safe_objective_value(objective_fn, plus, penalty)
            minus_obj = safe_objective_value(objective_fn, minus, penalty)
            denom = plus[coord_idx] - minus[coord_idx]
            slope = np.nan if denom == 0 else float((plus_obj - minus_obj) / denom)

            row = {
                "label": label,
                "coordinate_index": coord_idx,
                "coordinate_name": coord_name,
                "x0_value": float(x0[coord_idx]),
                "eps": float(eps),
                "x_minus_value": float(minus[coord_idx]),
                "x_plus_value": float(plus[coord_idx]),
                "objective_x0": float(base_objective),
                "objective_x_minus": float(minus_obj),
                "objective_x_plus": float(plus_obj),
                "absolute_delta_minus": float(minus_obj - base_objective),
                "absolute_delta_plus": float(plus_obj - base_objective),
                "relative_delta_minus": float((minus_obj - base_objective) / max(abs(base_objective), 1.0)),
                "relative_delta_plus": float((plus_obj - base_objective) / max(abs(base_objective), 1.0)),
                "finite_difference_slope": slope,
                "minus_hit_penalty": bool(minus_obj >= penalty),
                "plus_hit_penalty": bool(plus_obj >= penalty),
            }
            try:
                plus_flat = flatten_partial_pooling_parameters(unpack_fn(plus))
                minus_flat = flatten_partial_pooling_parameters(unpack_fn(minus))
                changed_keys = sorted(
                    key for key in set(base_flat) | set(plus_flat) | set(minus_flat)
                    if (
                        not np.isclose(float(base_flat.get(key, np.nan)), float(plus_flat.get(key, np.nan)), equal_nan=True)
                        or not np.isclose(float(base_flat.get(key, np.nan)), float(minus_flat.get(key, np.nan)), equal_nan=True)
                    )
                )
                row["changed_biological_parameters"] = ",".join(changed_keys)
            except Exception as exc:
                row["changed_biological_parameters"] = f"UNPACK_ERROR:{exc}"
            if debug_eval_fn is not None:
                plus_debug = debug_eval_fn(plus)
                minus_debug = debug_eval_fn(minus)
                row["plus_simulations_attempted"] = plus_debug.get("simulations_attempted")
                row["plus_simulations_failed"] = plus_debug.get("simulations_failed")
                row["plus_failure_messages"] = plus_debug.get("failure_messages_joined")
                row["plus_penalty_hit"] = plus_debug.get("penalty_hit")
                row["minus_simulations_attempted"] = minus_debug.get("simulations_attempted")
                row["minus_simulations_failed"] = minus_debug.get("simulations_failed")
                row["minus_failure_messages"] = minus_debug.get("failure_messages_joined")
                row["minus_penalty_hit"] = minus_debug.get("penalty_hit")
            rows.append(row)
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def build_parameter_unpacking_check(
    x0: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    unpack_fn: Callable[[np.ndarray], Dict[str, Any]],
    coordinate_names: Sequence[str],
    output_path: Path,
    eps: float = 1e-4,
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    lower_arr = np.asarray(lower_bounds, dtype=float)
    upper_arr = np.asarray(upper_bounds, dtype=float)
    base_flat = flatten_partial_pooling_parameters(unpack_fn(x0))
    rows: List[Dict[str, Any]] = []
    for coord_idx, coord_name in enumerate(coordinate_names):
        plus = x0.copy()
        minus = x0.copy()
        plus[coord_idx] = min(upper_arr[coord_idx], plus[coord_idx] + eps)
        minus[coord_idx] = max(lower_arr[coord_idx], minus[coord_idx] - eps)
        plus_flat = flatten_partial_pooling_parameters(unpack_fn(plus))
        minus_flat = flatten_partial_pooling_parameters(unpack_fn(minus))
        for bio_name in sorted(set(base_flat) | set(plus_flat) | set(minus_flat)):
            rows.append({
                "coordinate_index": coord_idx,
                "coordinate_name": coord_name,
                "biological_parameter": bio_name,
                "value_x0": base_flat.get(bio_name, np.nan),
                "value_plus": plus_flat.get(bio_name, np.nan),
                "value_minus": minus_flat.get(bio_name, np.nan),
            })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def project_to_bounds(
    x: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
) -> Tuple[np.ndarray, bool, np.ndarray, np.ndarray]:
    x_arr = np.asarray(x, dtype=float).copy()
    lower_arr = np.asarray(lower_bounds, dtype=float)
    upper_arr = np.asarray(upper_bounds, dtype=float)
    clipped = np.clip(x_arr, lower_arr, upper_arr)
    any_clipped = bool(np.any(clipped != x_arr))
    at_lower = np.isclose(clipped, lower_arr, atol=1e-12, rtol=0.0)
    at_upper = np.isclose(clipped, upper_arr, atol=1e-12, rtol=0.0)
    return clipped, any_clipped, at_lower, at_upper


def manual_coordinate_descent_probe(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    coordinate_names: Sequence[str],
    output_path: Path,
    penalty: float,
    diagnostic_label: str,
    steps: Sequence[float] = (-1e-2, -3e-3, -1e-3, -3e-4, -1e-4, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2),
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    base_objective = safe_objective_value(objective_fn, x0, penalty)
    rows: List[Dict[str, Any]] = []
    for coord_idx, coord_name in enumerate(coordinate_names):
        for step in steps:
            x_trial = x0.copy()
            x_trial[coord_idx] += step
            x_trial, _, at_lower, at_upper = project_to_bounds(x_trial, lower_bounds, upper_bounds)
            trial_objective = safe_objective_value(objective_fn, x_trial, penalty)
            rows.append({
                "diagnostic_label": diagnostic_label,
                "coordinate_index": coord_idx,
                "coordinate_name": coord_name,
                "step": float(step),
                "objective_at_x0": float(base_objective),
                "objective_trial": float(trial_objective),
                "delta_objective": float(trial_objective - base_objective),
                "improved": bool(trial_objective < base_objective),
                "hit_penalty": bool(trial_objective >= penalty),
                "at_lower_bound_after_step": bool(at_lower[coord_idx]),
                "at_upper_bound_after_step": bool(at_upper[coord_idx]),
            })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def manual_gradient_probe(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    coordinate_names: Sequence[str],
    output_path: Path,
    penalty: float,
    h_values: Sequence[float] = (1e-4, 1e-3),
) -> Tuple[pd.DataFrame, Dict[float, np.ndarray]]:
    x0 = np.asarray(x0, dtype=float)
    rows: List[Dict[str, Any]] = []
    gradients: Dict[float, np.ndarray] = {}
    for h in h_values:
        gradient = np.full(len(x0), np.nan, dtype=float)
        for coord_idx, coord_name in enumerate(coordinate_names):
            x_minus = x0.copy()
            x_plus = x0.copy()
            x_minus[coord_idx] -= h
            x_plus[coord_idx] += h
            x_minus, _, _, _ = project_to_bounds(x_minus, lower_bounds, upper_bounds)
            x_plus, _, _, _ = project_to_bounds(x_plus, lower_bounds, upper_bounds)
            f_minus = safe_objective_value(objective_fn, x_minus, penalty)
            f_plus = safe_objective_value(objective_fn, x_plus, penalty)
            denom = x_plus[coord_idx] - x_minus[coord_idx]
            g_i = np.nan if denom == 0 else float((f_plus - f_minus) / denom)
            gradient[coord_idx] = g_i
            rows.append({
                "coordinate_index": coord_idx,
                "coordinate_name": coord_name,
                "h": float(h),
                "f_minus": float(f_minus),
                "f_plus": float(f_plus),
                "gradient": g_i,
                "hit_penalty_minus": bool(f_minus >= penalty),
                "hit_penalty_plus": bool(f_plus >= penalty),
            })
        gradients[float(h)] = gradient
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df, gradients


def gradient_comparison_table(
    manual_gradients: Dict[float, np.ndarray],
    optimizer_result: Any,
    coordinate_names: Sequence[str],
    output_path: Path,
) -> pd.DataFrame:
    opt_jac = getattr(optimizer_result, "jac", None)
    rows: List[Dict[str, Any]] = []
    if opt_jac is None:
        df = pd.DataFrame(rows)
        df.to_csv(output_path, sep="\t", index=False)
        return df
    opt_grad = np.asarray(opt_jac, dtype=float)
    opt_norm = float(np.linalg.norm(opt_grad))
    for h, grad in manual_gradients.items():
        grad_arr = np.asarray(grad, dtype=float)
        manual_norm = float(np.linalg.norm(grad_arr))
        denom = manual_norm * opt_norm
        cosine = np.nan
        if denom > 0 and np.all(np.isfinite(grad_arr)) and np.all(np.isfinite(opt_grad)):
            cosine = float(np.dot(grad_arr, opt_grad) / denom)
        max_manual_idx = int(np.nanargmax(np.abs(grad_arr))) if grad_arr.size > 0 else -1
        max_opt_idx = int(np.nanargmax(np.abs(opt_grad))) if opt_grad.size > 0 else -1
        rows.append({
            "h": float(h),
            "manual_gradient_norm": manual_norm,
            "optimizer_gradient_norm": opt_norm,
            "norm_ratio_manual_over_optimizer": manual_norm / opt_norm if opt_norm > 0 else np.nan,
            "cosine_similarity": cosine,
            "manual_max_abs_gradient": float(np.nanmax(np.abs(grad_arr))) if grad_arr.size > 0 else np.nan,
            "manual_max_abs_gradient_name": coordinate_names[max_manual_idx] if max_manual_idx >= 0 else "",
            "optimizer_max_abs_gradient": float(np.nanmax(np.abs(opt_grad))) if opt_grad.size > 0 else np.nan,
            "optimizer_max_abs_gradient_name": coordinate_names[max_opt_idx] if max_opt_idx >= 0 else "",
        })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def manual_line_search_probe(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    direction: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    output_path: Path,
    penalty: float,
    alphas: Sequence[float],
    debug_eval_fn: Optional[Callable[[np.ndarray], Dict[str, Any]]] = None,
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    direction = np.asarray(direction, dtype=float)
    base_objective = safe_objective_value(objective_fn, x0, penalty)
    rows: List[Dict[str, Any]] = []
    for alpha in alphas:
        x_trial = x0 + float(alpha) * direction
        x_trial, any_clipped, _, _ = project_to_bounds(x_trial, lower_bounds, upper_bounds)
        trial_objective = safe_objective_value(objective_fn, x_trial, penalty)
        row = {
            "alpha": float(alpha),
            "objective_at_x0": float(base_objective),
            "objective_trial": float(trial_objective),
            "delta_objective": float(trial_objective - base_objective),
            "improved": bool(trial_objective < base_objective),
            "hit_penalty": bool(trial_objective >= penalty),
            "any_bound_clipped": bool(any_clipped),
        }
        if debug_eval_fn is not None:
            debug_info = debug_eval_fn(x_trial)
            row["simulations_failed"] = debug_info.get("simulations_failed")
            row["failure_messages"] = debug_info.get("failure_messages_joined")
        rows.append(row)
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def manual_best_coordinate_line_search(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    coordinate_index: int,
    coordinate_name: str,
    direction_sign: float,
    output_path: Path,
    penalty: float,
    alpha_grid: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    base_objective = safe_objective_value(objective_fn, x0, penalty)
    alpha_grid = alpha_grid if alpha_grid is not None else np.logspace(-6, -0.5, 20)
    rows: List[Dict[str, Any]] = []
    for alpha in alpha_grid:
        signed_step = float(direction_sign * alpha)
        x_trial = x0.copy()
        x_trial[coordinate_index] += signed_step
        x_trial, _, _, _ = project_to_bounds(x_trial, lower_bounds, upper_bounds)
        trial_objective = safe_objective_value(objective_fn, x_trial, penalty)
        rows.append({
            "coordinate_index": int(coordinate_index),
            "coordinate_name": coordinate_name,
            "alpha": float(alpha),
            "signed_step": signed_step,
            "objective_trial": float(trial_objective),
            "delta_objective": float(trial_objective - base_objective),
            "improved": bool(trial_objective < base_objective),
            "hit_penalty": bool(trial_objective >= penalty),
        })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def diagnostic_alternative_optimizers(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    bounds: Sequence[Tuple[float, float]],
    output_path: Path,
    maxiter: int,
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    initial_objective = float(objective_fn(x0))
    diagnostic_maxiter = min(int(maxiter), 20)
    diagnostic_maxfev = min(max(200, diagnostic_maxiter * max(1, len(x0))), 1000)
    rows: List[Dict[str, Any]] = []
    for method, options in [
        ("Powell", {"maxiter": diagnostic_maxiter, "maxfev": diagnostic_maxfev}),
        ("Nelder-Mead", {"maxiter": diagnostic_maxiter, "maxfev": diagnostic_maxfev}),
    ]:
        try:
            result = minimize(
                objective_fn,
                x0=x0,
                method=method,
                bounds=bounds,
                options=options,
            )
            movement_norm = float(np.linalg.norm(np.asarray(result.x, dtype=float) - x0))
            rows.append({
                "method": method,
                "initial_objective": initial_objective,
                "final_objective": float(objective_fn(result.x)),
                "delta_objective": float(objective_fn(result.x) - initial_objective),
                "n_iter": getattr(result, "nit", np.nan),
                "n_function_eval": getattr(result, "nfev", np.nan),
                "success": bool(result.success),
                "message": str(getattr(result, "message", "")),
                "result_equals_start": bool(np.allclose(np.asarray(result.x, dtype=float), x0, rtol=0.0, atol=1e-12)),
                "best_parameter_movement_norm": movement_norm,
            })
        except Exception as exc:
            rows.append({
                "method": method,
                "initial_objective": initial_objective,
                "final_objective": np.nan,
                "delta_objective": np.nan,
                "n_iter": np.nan,
                "n_function_eval": np.nan,
                "success": False,
                "message": str(exc),
                "result_equals_start": np.nan,
                "best_parameter_movement_norm": np.nan,
            })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def diagnostic_normalized_objective_lbfgsb(
    objective_fn: Callable[[np.ndarray], float],
    x0: np.ndarray,
    bounds: Sequence[Tuple[float, float]],
    output_path: Path,
    maxiter: int,
    n_observations: Optional[int] = None,
    penalty: float = 1e30,
) -> pd.DataFrame:
    x0 = np.asarray(x0, dtype=float)
    safe_objective = lambda x: safe_objective_value(objective_fn, np.asarray(x, dtype=float), penalty)
    f0 = float(safe_objective(x0))
    normalization_cases = [
        ("abs_f0", max(abs(f0), 1.0)),
    ]
    if n_observations is not None and n_observations > 0:
        normalization_cases.append(("per_observation", float(n_observations)))
    rows: List[Dict[str, Any]] = []
    for normalization_type, scale in normalization_cases:
        def normalized_objective(x: np.ndarray, denom: float = float(scale)) -> float:
            return float(safe_objective(x) / denom)

        result = minimize(
            normalized_objective,
            x0=x0,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": min(int(maxiter), 20), "maxfun": 50},
        )
        final_raw = float(safe_objective(result.x))
        rows.append({
            "normalization_type": normalization_type,
            "initial_objective_raw": f0,
            "final_objective_raw": final_raw,
            "initial_objective_normalized": float(normalized_objective(x0)),
            "final_objective_normalized": float(normalized_objective(result.x)),
            "delta_objective_raw": float(final_raw - f0),
            "n_iter": getattr(result, "nit", np.nan),
            "n_function_eval": getattr(result, "nfev", np.nan),
            "success": bool(result.success),
            "message": str(getattr(result, "message", "")),
            "result_equals_start": bool(np.allclose(np.asarray(result.x, dtype=float), x0, rtol=0.0, atol=1e-12)),
            "movement_norm": float(np.linalg.norm(np.asarray(result.x, dtype=float) - x0)),
        })
    df = pd.DataFrame(rows)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    return df


def dfdctp_ng_per_ml_to_uM(
    dfdctp_ng_per_ml: float,
    molecular_weight_ng_per_nmol: float = DFDCTP_MOLECULAR_WEIGHT_NG_PER_NMOL,
) -> float:
    """Converts dFdCTP concentration from ng/mL to uM using its molecular weight."""
    if molecular_weight_ng_per_nmol <= 0:
        raise ValueError(
            "molecular_weight_ng_per_nmol must be > 0, "
            f"got {molecular_weight_ng_per_nmol}"
        )
    return dfdctp_ng_per_ml / molecular_weight_ng_per_nmol

#################### r, K fitting utilities ####################

def parse_gem_label(val):
    """Extracts standardized Gemcitabine concentration from raw text."""
    if pd.isna(val):
        return np.nan
    val_str = str(val)
    match = re.search(r'(?i)\b([0-9]+\.?[0-9]*)\s*(nM|uM|µM|pM)\b', val_str)
    if not match:
        return np.nan
    
    amount, unit = match.groups()
    unit = unit.lower().replace('µm', 'uM').replace('um', 'uM').replace('nm', 'nM')
    return f"{amount} {unit}"


def _standardize_platemap_frame(plate_long: pd.DataFrame) -> pd.DataFrame:
    required = ["plate_row", "plate_col", "ploidy", "gem", "condition_raw"]
    missing = [col for col in required if col not in plate_long.columns]
    if missing:
        raise ValueError(f"Platemap parser did not produce required columns: {missing}")

    standardized = plate_long[required].copy()
    standardized["plate_row"] = standardized["plate_row"].astype(str).str.upper().str.strip()
    standardized["plate_col"] = pd.to_numeric(standardized["plate_col"], errors="coerce").astype("Int64")
    standardized["ploidy"] = standardized["ploidy"].apply(parse_ploidy_label)
    standardized["gem"] = standardized["gem"].apply(parse_gem_label)

    invalid_coords = standardized[
        ~standardized["plate_row"].isin(list("ABCDEFGH")) |
        standardized["plate_col"].isna() |
        (standardized["plate_col"] < 1) |
        (standardized["plate_col"] > 12)
    ]
    if not invalid_coords.empty:
        raise ValueError(
            "Invalid platemap coordinates encountered: "
            f"{invalid_coords[['plate_row', 'plate_col']].drop_duplicates().to_dict(orient='records')}"
        )

    invalid_ploidy = standardized["ploidy"].notna() & ~standardized["ploidy"].isin(["2N", "4N"])
    if invalid_ploidy.any():
        raise ValueError(
            "Invalid platemap ploidy values encountered: "
            f"{standardized.loc[invalid_ploidy, 'ploidy'].dropna().unique().tolist()}"
        )

    duplicate_mask = standardized.duplicated(subset=["plate_row", "plate_col"], keep=False)
    if duplicate_mask.any():
        conflicting = []
        for key, group in standardized.loc[duplicate_mask].groupby(["plate_row", "plate_col"], dropna=False):
            if len(group.drop_duplicates()) > 1:
                conflicting.append((key, group.to_dict(orient="records")))
        if conflicting:
            raise ValueError(f"Conflicting duplicate platemap entries: {conflicting}")
        standardized = standardized.drop_duplicates()

    return standardized.reset_index(drop=True)

def load_and_clean_counts(
    max_days: float = 5.0,
    paths: Optional[ExperimentPaths] = None,
) -> pd.DataFrame:
    """Loads the aggregated parquet counts data, cleans columns, and filters time."""
    paths = paths or default_experiment_paths()
    validate_paths(paths, required=("counts_agg",))
    df = pd.read_parquet(paths.counts_agg)
    
    df = df.rename(columns={"Classifier.Phenotype": "phenotype", "N": "count"})
    
    if 'plate_row' not in df.columns or 'plate_col' not in df.columns:
        extracted = df['well'].str.extract(r'^([A-H])([0-9]{1,2})(?:_([0-9]+))?$')
        df['plate_row'] = extracted[0]
        df['plate_col'] = pd.to_numeric(extracted[1], errors='coerce').astype('Int64')
        
    if 'time_hours' not in df.columns:
        raise ValueError("Expected 'time_hours' in counts table.")
    
    df['time_days'] = df['time_hours'] / 24.0
    df = df[df['time_days'] <= max_days].copy()
    
    df = df[df['phenotype'].isin(['Alive', 'Dead', 'Transitional'])]
    
    df = df.dropna(subset=['plate_row', 'plate_col'])
    df = df[(df['plate_col'] >= 1) & (df['plate_col'] <= 12)]
    
    return df

def load_platemap(paths: Optional[ExperimentPaths] = None) -> pd.DataFrame:
    """Loads and parses supported platemap formats into a standardized layout."""
    paths = paths or default_experiment_paths()
    validate_paths(paths, required=("platemap",))
    pm = pd.read_excel(paths.platemap, sheet_name=0)
    pm = pm.rename(columns={col: normalize_column_name(col) for col in pm.columns})
    normalized_cols = list(pm.columns)

    well_col = next((c for c in normalized_cols if c in {"well", "wellid"}), None)
    row_col = next((c for c in normalized_cols if c == "row"), None)
    col_col = next((c for c in normalized_cols if c in {"col", "column"}), None)
    condition_col = next((c for c in normalized_cols if c in {"condition", "condition_raw", "treatment", "label"}), None)
    ploidy_col = next((c for c in normalized_cols if c == "ploidy"), None)
    gem_col = next((c for c in normalized_cols if c in {"gem", "dose"}), None)

    if well_col is not None or (row_col is not None and col_col is not None):
        plate_long = pm.copy()
        if well_col is not None:
            parsed = plate_long[well_col].apply(parse_well_id)
            plate_long["plate_row"] = parsed.apply(lambda item: item[0])
            plate_long["plate_col"] = parsed.apply(lambda item: item[1])
        else:
            plate_long["plate_row"] = plate_long[row_col].astype(str).str.upper().str.strip()
            plate_long["plate_col"] = pd.to_numeric(plate_long[col_col], errors="coerce")

        if condition_col is not None:
            plate_long["condition_raw"] = plate_long[condition_col].astype(str)
            plate_long["ploidy"] = plate_long["condition_raw"].apply(parse_ploidy_label)
            plate_long["gem"] = plate_long["condition_raw"].apply(parse_gem_label)
        elif ploidy_col is not None and gem_col is not None:
            plate_long["ploidy"] = plate_long[ploidy_col]
            plate_long["gem"] = plate_long[gem_col]
            plate_long["condition_raw"] = (
                plate_long[ploidy_col].astype(str).fillna("") + " | " +
                plate_long[gem_col].astype(str).fillna("")
            )
        else:
            raise ValueError(
                "Unsupported long-format platemap. Detected columns "
                f"{normalized_cols} but missing either a combined condition column or separate ploidy/gem columns."
            )
        return _standardize_platemap_frame(plate_long)

    first_col = pm.columns[0]
    row_labels = pm[first_col].astype(str).str.upper().str.strip()
    value_cols = [c for c in pm.columns[1:] if str(c).isdigit() and 1 <= int(str(c)) <= 12]
    if row_labels.isin(list("ABCDEFGH")).any() and value_cols:
        wide = pm.copy()
        wide = wide.rename(columns={first_col: "row"})
        wide["row"] = wide["row"].astype(str).str.upper().str.strip()
        wide = wide[wide["row"].isin(list("ABCDEFGH"))].copy()
        plate_long = pd.melt(
            wide,
            id_vars="row",
            value_vars=value_cols,
            var_name="plate_col",
            value_name="condition_raw",
        )
        plate_long["plate_row"] = plate_long["row"]
        plate_long["ploidy"] = plate_long["condition_raw"].apply(parse_ploidy_label)
        plate_long["gem"] = plate_long["condition_raw"].apply(parse_gem_label)
        return _standardize_platemap_frame(plate_long)

    raise ValueError(
        "Unsupported platemap format. Detected normalized columns "
        f"{normalized_cols}. Expected wide row/1-12 layout or long format with well or row+column fields."
    )

def assemble_modeling_dataset(
    paths: Optional[ExperimentPaths] = None,
    fit_config: Optional[JointFitConfig] = None,
) -> pd.DataFrame:
    """Combines counts and platemap into a single dataframe."""
    paths = paths or default_experiment_paths()
    fit_config = resolve_joint_fit_config(fit_config)
    counts_df = load_and_clean_counts(max_days=fit_config.max_days, paths=paths)
    platemap_df = load_platemap(paths=paths)
    
    merged = pd.merge(counts_df, platemap_df, on=['plate_row', 'plate_col'], how='left')
    return merged


def get_alive_like_phenotypes(count_transitional_as_alive: bool) -> Tuple[str, ...]:
    return ("Alive", "Transitional") if count_transitional_as_alive else ("Alive",)


def get_dead_like_phenotypes(count_transitional_as_alive: bool) -> Tuple[str, ...]:
    return ("Dead",) if count_transitional_as_alive else ("Dead", "Transitional")

def save_and_maybe_close(fig, output_path: Optional[Path] = None, close: bool = True, **savefig_kwargs) -> None:
    """Saves a Matplotlib figure and closes it by default for batch-script safety."""
    if output_path is not None:
        fig.savefig(output_path, **savefig_kwargs)
    if close and (output_path is not None or close):
        plt.close(fig)

def _build_phenotype_pivot(
    df: pd.DataFrame,
    gem_dose: str,
    ploidy: str,
    phenotype: Union[str, Sequence[str]],
    plate_row: Optional[str] = None,
) -> pd.DataFrame:
    """
    Aggregates counts by time/replicate and returns a pivoted matrix.

    Duplicate entries for the same time/replicate/phenotype are summed because
    the count column is an observed cell count.
    """
    phenotype_values = (phenotype,) if isinstance(phenotype, str) else tuple(phenotype)
    subset = df[
        (df['gem'] == gem_dose) &
        (df['ploidy'] == ploidy) &
        (df['phenotype'].isin(phenotype_values))
    ].copy()
    if plate_row is not None:
        subset = subset[subset['plate_row'] == plate_row].copy()
    if subset.empty:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No data found for {gem_dose}, {ploidy}, {phenotype_values}{scope}")

    grouped = (
        subset
        .groupby(['time_days', 'plate_row', 'plate_col'], as_index=False)['count']
        .sum()
    )
    return grouped.pivot_table(
        index='time_days',
        columns=['plate_row', 'plate_col'],
        values='count',
        aggfunc='sum',
    ).sort_index()

def get_aligned_live_dead_data(
    df: pd.DataFrame,
    gem_dose: str,
    ploidy: str,
    plate_row: Optional[str] = None,
    t_max: Optional[float] = None,
    count_transitional_as_alive: bool = False,
) -> Dict[str, Any]:
    """
    Inner-joins alive-like and dead-like trajectories by shared time and replicate identity.

    Alignment is performed on `time_days` and the replicate keys
    `plate_row`/`plate_col`. Time points or replicate columns that do not appear
    in both observed compartments are dropped explicitly and reported.
    """
    alive_pivot = _build_phenotype_pivot(
        df,
        gem_dose,
        ploidy,
        get_alive_like_phenotypes(count_transitional_as_alive),
        plate_row=plate_row,
    )
    dead_pivot = _build_phenotype_pivot(
        df,
        gem_dose,
        ploidy,
        get_dead_like_phenotypes(count_transitional_as_alive),
        plate_row=plate_row,
    )

    common_times = alive_pivot.index.intersection(dead_pivot.index).sort_values()
    common_cols = alive_pivot.columns.intersection(dead_pivot.columns)
    if len(common_times) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead time points for {gem_dose}, {ploidy}{scope}")
    if len(common_cols) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead replicate columns for {gem_dose}, {ploidy}{scope}")

    alive_time_set = set(alive_pivot.index.tolist())
    dead_time_set = set(dead_pivot.index.tolist())
    dropped_timepoints = len(alive_time_set.union(dead_time_set)) - len(common_times)
    dropped_replicates = len(set(alive_pivot.columns.tolist()).union(set(dead_pivot.columns.tolist()))) - len(common_cols)

    if t_max is not None:
        common_times = common_times[common_times <= t_max]
    if len(common_times) == 0:
        scope = f", row {plate_row}" if plate_row is not None else ""
        raise ValueError(f"No overlapping Alive/Dead time points remain after t_max for {gem_dose}, {ploidy}{scope}")

    alive_aligned = alive_pivot.loc[common_times, common_cols]
    dead_aligned = dead_pivot.loc[common_times, common_cols]

    return {
        "t": common_times.to_numpy(dtype=float),
        "y_alive": alive_aligned.to_numpy(dtype=float),
        "y_dead": dead_aligned.to_numpy(dtype=float),
        "replicate_columns": list(common_cols),
        "dropped_timepoints": int(dropped_timepoints),
        "dropped_replicates": int(dropped_replicates),
        "plate_row": plate_row,
        "count_transitional_as_alive": bool(count_transitional_as_alive),
    }

def get_fitting_data(
    df: pd.DataFrame,
    gem_dose: str = "0 nM",
    ploidy: str = "2N",
    phenotype: str = "Alive",
    count_transitional_as_alive: bool = False,
):
    """
    Backward-compatible wrapper returning one phenotype matrix after explicit
    Alive/Dead alignment.
    """
    aligned = get_aligned_live_dead_data(
        df,
        gem_dose=gem_dose,
        ploidy=ploidy,
        count_transitional_as_alive=count_transitional_as_alive,
    )
    if phenotype == "Alive":
        return aligned["t"], aligned["y_alive"]
    if phenotype == "Dead":
        return aligned["t"], aligned["y_dead"]
    raise ValueError(f"Unsupported phenotype {phenotype}; expected Alive or Dead")

def logistic_growth(t, N0, r, K):
    """Analytical solution to the logistic growth equation."""
    return (K * N0) / (N0 + (K - N0) * np.exp(-r * t))

# Canonical internal unit convention for the joint in vitro ODE:
# - time `t` is in days
# - `dose_muM` is nominal extracellular gemcitabine dose in uM
# - `C_dfdctp_uM` is intracellular dFdCTP concentration/effective concentration in uM
# - `k_tr` and `k_clear` are in day^-1
# - `k_kill` is in day^-1 per uM intracellular dFdCTP, so
#   `kappa = k_kill * delayed_C_dfdctp_uM` is a death hazard in day^-1
# - `k_cyto` is in uM^-1 and acts through a bounded proliferation multiplier

def decay_rate_from_half_life_days(half_life_days: float) -> float:
    """Converts an extracellular half-life in days to a decay rate in days^-1."""
    if half_life_days <= 0:
        raise ValueError(f"half_life_days must be > 0, got {half_life_days}")
    return np.log(2.0) / half_life_days


def cytostasis_multiplier(c_dfdctp_uM: float, k_cyto: float) -> float:
    if k_cyto < 0:
        raise ValueError("k_cyto must be >= 0")
    c = max(float(c_dfdctp_uM), 0.0)
    return 1.0 / (1.0 + float(k_cyto) * c)


def residuals_fixed_N0(params, t, y_true, N0_fixed):
    """Objective function that only fits r and K, keeping N0 fixed."""
    r, K = params
    y_pred = logistic_growth(t, N0_fixed, r, K)
    return y_true - y_pred

def logistic_growth_replicate_n0(t, n0_by_rep, r, K):
    """Evaluates logistic growth with shared r/K and replicate-specific initial counts."""
    t_arr = np.asarray(t, dtype=float)
    n0_arr = np.asarray(n0_by_rep, dtype=float)
    return logistic_growth(t_arr[:, None], n0_arr[None, :], r, K)

def residuals_shared_rk_replicate_n0(params, t_data, y_alive, n0_by_rep, t_max=None):
    """Residuals for shared r/K baseline growth with replicate-specific fixed N0."""
    r, K = params
    y_pred = logistic_growth_replicate_n0(t_data, n0_by_rep, r, K)
    valid_mask = ~np.isnan(y_alive)
    if t_max is not None:
        valid_mask = valid_mask & (np.asarray(t_data)[:, None] <= t_max)
    return (y_alive - y_pred)[valid_mask]

def fit_baseline_shared_rk_replicate_n0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None, close_fig: bool = True):
    """Fits shared cohort r/K while fixing each replicate to its observed baseline N0."""
    y_alive_matrix = y_alive[:, None] if y_alive.ndim == 1 else y_alive
    n0_by_replicate = y_alive_matrix[0, :].astype(float)
    replicate_mask = ~np.isnan(n0_by_replicate)
    if not np.any(replicate_mask):
        raise ValueError(f"No valid baseline alive counts found for {ploidy_label}")

    y_alive_matrix = y_alive_matrix[:, replicate_mask]
    n0_by_replicate = n0_by_replicate[replicate_mask]

    residual_mask = ~np.isnan(y_alive_matrix)
    if t_max is not None:
        residual_mask = residual_mask & (np.asarray(t_data)[:, None] <= t_max)
    n_observations = int(np.count_nonzero(residual_mask))
    if n_observations < 2:
        raise ValueError(f"Not enough alive observations to fit baseline growth for {ploidy_label}")

    K_guess = np.nanmax(y_alive_matrix[residual_mask])
    r_guess = 0.5

    res = least_squares(
        residuals_shared_rk_replicate_n0,
        [r_guess, K_guess],
        args=(t_data, y_alive_matrix, n0_by_replicate, t_max),
        bounds=([0, 0], [np.inf, np.inf]),
    )
    r_opt, K_opt = res.x
    n0_mean = float(np.nanmean(n0_by_replicate))
    n0_min = float(np.nanmin(n0_by_replicate))
    n0_max = float(np.nanmax(n0_by_replicate))

    print(f"\n{'='*60}")
    print(f"{ploidy_label} baseline fit")
    print("replicate-specific N0 enabled")
    if t_max is not None:
        print(f"Fitted only on t <= {t_max} Days")
    print(f"Replicates/wells used: {len(n0_by_replicate)}")
    print(f"N0 mean: {n0_mean:.1f}")
    print(f"N0 range: {n0_min:.1f} to {n0_max:.1f}")
    print(f"Observations used: {n_observations}")
    print(f"r = {r_opt:.4f} days^-1")
    print(f"K = {K_opt:.1f}")
    print("=" * 60 + "\n")

    if output_dir is not None or not close_fig:
        fig, ax1 = plt.subplots(figsize=(8, 5))
        color1 = 'tab:purple'
        ax1.set_xlabel('Time (Days)')
        ax1.set_ylabel('Alive Cell Count', color=color1)

        for rep_idx in range(y_alive_matrix.shape[1]):
            ax1.plot(
                t_data,
                y_alive_matrix[:, rep_idx],
                'o',
                color=color1,
                alpha=0.18,
                label='Alive (Replicates)' if rep_idx == 0 else "",
            )

        ax1.plot(t_data, np.nanmean(y_alive_matrix, axis=1), 's', color='indigo', label='Alive (Mean)')

        t_smooth = np.linspace(min(t_data), max(t_data), 100)
        fitted_curves = logistic_growth_replicate_n0(t_smooth, n0_by_replicate, r_opt, K_opt)
        for rep_idx in range(fitted_curves.shape[1]):
            ax1.plot(
                t_smooth,
                fitted_curves[:, rep_idx],
                '-',
                color='fuchsia',
                alpha=0.25,
                linewidth=1.2,
                label='Replicate-specific Fits' if rep_idx == 0 else "",
            )
        ax1.plot(
            t_smooth,
            logistic_growth(t_smooth, n0_mean, r_opt, K_opt),
            '-',
            color='darkmagenta',
            linewidth=2.2,
            label=f'Mean-N0 Summary Fit (r={r_opt:.3f})',
        )

        if t_max is not None:
            ax1.axvline(x=t_max, color='gray', linestyle='--', alpha=0.7, label=f'Fit Cutoff (t={t_max})')

        ax1.tick_params(axis='y', labelcolor=color1)

        plt.title(f"Phase 1 Control Fit: {ploidy_label} (Replicate-specific N0)")
        plt.legend(loc='upper left')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        if output_dir is not None:
            output_path = output_dir / f"baseline_replicate_n0_{slugify_label(ploidy_label)}.png"
            save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')

    return {
        "N0_by_replicate": n0_by_replicate,
        "N0_mean": n0_mean,
        "r": float(r_opt),
        "K": float(K_opt),
        "ploidy": ploidy_label,
        "n_replicates": int(len(n0_by_replicate)),
        "n_observations": n_observations,
    }

def fit_baseline_locked_N0(t_data, y_alive, y_dead, ploidy_label, t_max=None, output_dir: Optional[Path] = None, close_fig: bool = True):
    """Legacy two-stage fit helper: baseline-only r/K fit with replicate-specific N0."""
    return fit_baseline_shared_rk_replicate_n0(
        t_data,
        y_alive,
        y_dead,
        ploidy_label,
        t_max=t_max,
        output_dir=output_dir,
        close_fig=close_fig,
    )



#################### Cohort fitting utilities ####################

def single_clone_signal_ode_joint(
    y,
    t,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve: Callable[[Union[float, np.ndarray], float], Union[float, np.ndarray]],
    k_tr,
    k_kill,
    k_clear,
    n_tr,
    model_variant=JointFitConfig().model_variant,
    k_cyto=None,
    beta_dose=1.0,
    mu_base_death=0.0,
    use_confluence_death: bool = False,
    mu_confluence_death: float = 0.0,
    confluence_death_exponent: float = 4.0,
):
    """
    Backward-compatible name for the dFdCTP-driven live/dead ODE.

    The default model variant is now immediate cytostasis plus delayed death,
    not the older delayed-death-only behavior.
    """
    return single_clone_dfdctp_ode_joint(
        y=y,
        t=t,
        r=r,
        K=K,
        dose_muM=dose_muM,
        dfdctp_signal_curve=dfdctp_signal_curve,
        k_tr=k_tr,
        k_kill=k_kill,
        k_clear=k_clear,
        n_tr=n_tr,
        model_variant=model_variant,
        k_cyto=k_cyto,
        beta_dose=beta_dose,
        mu_base_death=mu_base_death,
        use_confluence_death=use_confluence_death,
        mu_confluence_death=mu_confluence_death,
        confluence_death_exponent=confluence_death_exponent,
    )


def simulate_joint_ext(
    t,
    params,
    N0,
    D0,
    r,
    K,
    dose_muM,
    dfdctp_signal_curve,
    n_tr,
    model_variant=JointFitConfig().model_variant,
    fit_config: Optional[JointFitConfig] = None,
    beta_dose: Optional[float] = None,
    use_confluence_death: Optional[bool] = None,
    mu_confluence_death: Optional[float] = None,
    confluence_death_exponent: Optional[float] = None,
    use_hill_dose_gate: Optional[bool] = None,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
):
    """
    Backward-compatible wrapper for the dFdCTP-driven live/dead simulation.

    The default model variant is now immediate cytostasis plus delayed death.
    """
    return simulate_joint_dfdctp(
        t=t,
        params=params,
        N0=N0,
        D0=D0,
        r=r,
        K=K,
        dose_muM=dose_muM,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr=n_tr,
        model_variant=model_variant,
        fit_config=fit_config,
        beta_dose=beta_dose,
        use_confluence_death=use_confluence_death,
        mu_confluence_death=mu_confluence_death,
        confluence_death_exponent=confluence_death_exponent,
        use_hill_dose_gate=use_hill_dose_gate,
        dose_gate_ec50_uM=dose_gate_ec50_uM,
        dose_gate_hill=dose_gate_hill,
    )


def poisson_nll(y, mu):
    """Poisson negative log-likelihood for nonnegative count observations."""
    y_arr = np.asarray(y, dtype=float)
    mu_arr = np.clip(np.asarray(mu, dtype=float), 1e-9, np.inf)
    if np.any(y_arr < 0):
        raise ValueError("Poisson observations must be nonnegative")
    return float(np.sum(mu_arr - y_arr * np.log(mu_arr) + gammaln(y_arr + 1.0)))


def negative_binomial_nll(y, mu, theta):
    """
    Negative binomial negative log-likelihood with mean/dispersion parameterization.

    Var(y) = mu + mu^2 / theta
    """
    if theta <= 0:
        raise ValueError(f"theta must be > 0, got {theta}")
    y_arr = np.asarray(y, dtype=float)
    mu_arr = np.clip(np.asarray(mu, dtype=float), 1e-9, np.inf)
    if np.any(y_arr < 0):
        raise ValueError("Negative binomial observations must be nonnegative")

    theta_val = float(theta)
    ll = (
        gammaln(y_arr + theta_val)
        - gammaln(theta_val)
        - gammaln(y_arr + 1.0)
        + theta_val * (np.log(theta_val) - np.log(theta_val + mu_arr))
        + y_arr * (np.log(mu_arr) - np.log(theta_val + mu_arr))
    )
    return float(-np.sum(ll))


def collect_alive_dead_observations(
    dose_data_list,
    treatment_params,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
):
    """
    Collects observed/predicted alive and dead counts for objective evaluation.

    `high_dose_weight` is an optional caller-specified dose weight. It is not a
    normalization term.
    """
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )
    fit_config = resolve_joint_fit_config(fit_config)
    model_variant = fit_config.model_variant
    dfdctp_signal_curve = assert_preferred_live_dead_driver(dfdctp_signal_curve)
    records = []

    for data in dose_data_list:
        t_data, y_alive_data, y_dead_data = data['t'], data['y_alive'], data['y_dead']
        N0, D0, dose_muM = data['N0'], data['D0'], data['dose_muM']

        sim_result = simulate_joint_dfdctp_safe(
            t=np.asarray(t_data, dtype=float),
            params=treatment_params,
            N0=N0,
            D0=D0,
            r=r_fixed,
            K=K_fixed,
            dose_muM=dose_muM,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=n_tr_test,
            model_variant=model_variant,
            fit_config=fit_config,
        )
        if not sim_result.success:
            raise RuntimeError(f"Simulation failed for dose {dose_muM}: {sim_result.message}")
        y_alive_pred, y_dead_pred = sim_result.alive, sim_result.dead

        if fit_means_only:
            y_alive_obs = np.nanmean(y_alive_data, axis=1) if np.ndim(y_alive_data) > 1 else np.asarray(y_alive_data, dtype=float)
            y_dead_obs = np.nanmean(y_dead_data, axis=1) if np.ndim(y_dead_data) > 1 else np.asarray(y_dead_data, dtype=float)
            y_alive_pred_expanded = np.asarray(y_alive_pred, dtype=float)
            y_dead_pred_expanded = np.asarray(y_dead_pred, dtype=float)
        else:
            y_alive_obs = y_alive_data.flatten() if np.ndim(y_alive_data) > 1 else np.asarray(y_alive_data, dtype=float)
            y_dead_obs = y_dead_data.flatten() if np.ndim(y_dead_data) > 1 else np.asarray(y_dead_data, dtype=float)
            if np.ndim(y_alive_data) > 1:
                y_alive_pred_expanded = np.repeat(np.asarray(y_alive_pred, dtype=float), y_alive_data.shape[1])
                y_dead_pred_expanded = np.repeat(np.asarray(y_dead_pred, dtype=float), y_dead_data.shape[1])
            else:
                y_alive_pred_expanded = np.asarray(y_alive_pred, dtype=float)
                y_dead_pred_expanded = np.asarray(y_dead_pred, dtype=float)

        weight = high_dose_weight if dose_muM >= 0.2 else 1.0

        alive_mask = ~np.isnan(y_alive_obs)
        dead_mask = ~np.isnan(y_dead_obs)
        records.append(
            {
                "dose_muM": dose_muM,
                "weight": float(weight),
                "alive_obs": np.asarray(y_alive_obs[alive_mask], dtype=float),
                "alive_pred": np.asarray(y_alive_pred_expanded[alive_mask], dtype=float),
                "dead_obs": np.asarray(y_dead_obs[dead_mask], dtype=float),
                "dead_pred": np.asarray(y_dead_pred_expanded[dead_mask], dtype=float),
            }
        )
    return records


def residuals_global_joint(treatment_params, dose_data_list, r_fixed, K_fixed, dfdctp_signal_curve, n_tr_test,
                           fit_means_only=True, high_dose_weight=1.0,
                           observation_channels="alive_only",
                           model_variant=JointFitConfig().model_variant,
                           fit_config: Optional[JointFitConfig] = None):
    """
    Residuals for the joint live/dead fit in raw observed cell-count units.

    `high_dose_weight` remains an optional caller-specified dose weight. It is
    not a normalization term and does not rescale residuals by within-dose
    maxima.
    """
    all_residuals = []
    fit_config = resolve_joint_fit_config(fit_config, model_variant=model_variant)

    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=n_tr_test,
            fit_config=fit_config,
            observation_channels=observation_channels,
            fit_means_only=fit_means_only,
            high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        all_residuals.extend((record["alive_obs"] - record["alive_pred"]) * weight)
        if observation_channels == "alive_dead":
            all_residuals.extend((record["dead_obs"] - record["dead_pred"]) * weight)

    return np.array(all_residuals)


def live_dead_objective_nll(
    treatment_params,
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr_test,
    objective="negative_binomial",
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
    theta_alive=None,
    theta_dead=None,
):
    """
    Likelihood objective for alive/dead counts using predicted ODE means.
    """
    if objective not in {"negative_binomial", "poisson"}:
        raise ValueError(f"Unsupported likelihood objective {objective}")
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )

    if objective == "negative_binomial":
        if theta_alive is None:
            raise ValueError("theta_alive is required for negative_binomial objective")
        if theta_alive <= 0:
            raise ValueError("theta_alive must be > 0")
        if observation_channels == "alive_dead":
            if theta_dead is None:
                raise ValueError("theta_dead is required for alive_dead negative_binomial objective")
            if theta_dead <= 0:
                raise ValueError("theta_dead must be > 0")

    total_nll = 0.0
    fit_config = resolve_joint_fit_config(fit_config)
    for record in collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr_test,
        fit_config=fit_config,
        observation_channels=observation_channels,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    ):
        weight = record["weight"]
        if objective == "poisson":
            total_nll += weight * poisson_nll(record["alive_obs"], record["alive_pred"])
            if observation_channels == "alive_dead":
                total_nll += weight * poisson_nll(record["dead_obs"], record["dead_pred"])
        else:
            total_nll += weight * negative_binomial_nll(record["alive_obs"], record["alive_pred"], theta_alive)
            if observation_channels == "alive_dead":
                total_nll += weight * negative_binomial_nll(record["dead_obs"], record["dead_pred"], theta_dead)
    return float(total_nll)


def summarize_live_dead_fit(
    dose_data_list,
    r_fixed,
    K_fixed,
    dfdctp_signal_curve,
    n_tr,
    treatment_params,
    objective,
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
    theta_alive=None,
    theta_dead=None,
):
    fit_config = resolve_joint_fit_config(fit_config)
    records = collect_alive_dead_observations(
        dose_data_list=dose_data_list,
        treatment_params=treatment_params,
        r_fixed=r_fixed,
        K_fixed=K_fixed,
        dfdctp_signal_curve=dfdctp_signal_curve,
        n_tr_test=n_tr,
        fit_config=fit_config,
        observation_channels=observation_channels,
        fit_means_only=fit_means_only,
        high_dose_weight=high_dose_weight,
    )
    residual_vector = np.concatenate(
        [
            (record["alive_obs"] - record["alive_pred"]) * record["weight"]
            for record in records
        ] + (
            [
                (record["dead_obs"] - record["dead_pred"]) * record["weight"]
                for record in records
            ] if observation_channels == "alive_dead" else []
        )
    ) if records else np.array([], dtype=float)

    n_obs = int(sum(
        len(record["alive_obs"]) + (len(record["dead_obs"]) if observation_channels == "alive_dead" else 0)
        for record in records
    ))
    rmse = float(np.sqrt(np.mean(residual_vector ** 2))) if residual_vector.size > 0 else np.nan
    summary = {
        "objective": objective,
        "observation_channels": observation_channels,
        "model_variant": fit_config.model_variant,
        "model_preset": fit_config.model_preset,
        "n_observations": n_obs,
        "rmse": rmse,
        "theta_alive": theta_alive,
        "theta_dead": theta_dead if observation_channels == "alive_dead" else None,
    }
    summary["fit_beta_dose"] = bool(fit_config.fit_beta_dose)
    summary["fixed_beta_dose"] = float(fit_config.fixed_beta_dose)
    unpacked_params = unpack_treatment_params_for_config(treatment_params, fit_config)
    summary.update(unpacked_params)
    summary[DOSE_SCALING_PARAMETER_NAME] = resolve_beta_dose_for_params(unpacked_params, fit_config)
    summary[CONFLUENCE_DEATH_PARAMETER_NAME] = resolve_mu_confluence_death_for_params(unpacked_params, fit_config)
    summary["use_hill_dose_gate"] = bool(fit_config.use_hill_dose_gate)
    summary["fit_hill_dose_gate"] = bool(fit_config.fit_hill_dose_gate)
    summary["use_confluence_death"] = bool(fit_config.use_confluence_death)
    summary["fit_mu_confluence_death"] = bool(fit_config.fit_mu_confluence_death)
    summary["fixed_mu_confluence_death"] = float(fit_config.fixed_mu_confluence_death)
    summary["confluence_death_exponent"] = float(fit_config.confluence_death_exponent)
    summary["dose_gate_ec50_uM"] = (
        float(fit_config.fixed_dose_gate_ec50_uM) if fit_config.use_hill_dose_gate else np.nan
    )
    summary["dose_gate_hill"] = (
        float(fit_config.fixed_dose_gate_hill) if fit_config.use_hill_dose_gate else np.nan
    )
    if objective == "least_squares":
        summary["objective_value"] = float(0.5 * np.sum(residual_vector ** 2))
        summary["nll"] = None
        summary["aic"] = None
        summary["bic"] = None
        summary["n_parameters"] = len(get_treatment_parameter_names_for_config(fit_config))
    else:
        nll = live_dead_objective_nll(
            treatment_params=treatment_params,
            dose_data_list=dose_data_list,
            r_fixed=r_fixed,
            K_fixed=K_fixed,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr_test=n_tr,
            objective=objective,
            fit_config=fit_config,
            observation_channels=observation_channels,
            fit_means_only=fit_means_only,
            high_dose_weight=high_dose_weight,
            theta_alive=theta_alive,
            theta_dead=theta_dead,
        )
        if objective == "negative_binomial":
            n_parameters = len(get_treatment_parameter_names_for_config(fit_config)) + (
                2 if observation_channels == "alive_dead" else 1
            )
        else:
            n_parameters = len(get_treatment_parameter_names_for_config(fit_config))
        summary["objective_value"] = float(nll)
        summary["nll"] = float(nll)
        summary["n_parameters"] = n_parameters
        summary["aic"] = float(2 * n_parameters + 2 * nll)
        summary["bic"] = float(n_parameters * np.log(max(n_obs, 1)) + 2 * nll)
    return summary


def get_legacy_hill_parameter_names_for_config(fit_config: JointFitConfig) -> Tuple[str, ...]:
    if fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate:
        return ("dose_gate_ec50_uM", "dose_gate_hill")
    return ()


def build_legacy_hill_eval_config(
    fit_config: JointFitConfig,
    dose_gate_ec50_uM: Optional[float],
    dose_gate_hill: Optional[float],
) -> JointFitConfig:
    if not fit_config.use_hill_dose_gate:
        return fit_config
    return replace(
        fit_config,
        fit_hill_dose_gate=False,
        fixed_dose_gate_ec50_uM=(
            fit_config.fixed_dose_gate_ec50_uM
            if dose_gate_ec50_uM is None else float(dose_gate_ec50_uM)
        ),
        fixed_dose_gate_hill=(
            fit_config.fixed_dose_gate_hill
            if dose_gate_hill is None else float(dose_gate_hill)
        ),
    )


def build_legacy_start_dicts_for_config(fit_config: JointFitConfig) -> List[Dict[str, float]]:
    model_variant = validate_model_variant(fit_config.model_variant)
    if model_variant == "delayed_death_only":
        starts = [
            {"k_tr": 0.2, "k_kill": 10.0, "k_clear": 0.2},
            {"k_tr": 0.5, "k_kill": 25.0, "k_clear": 0.5},
            {"k_tr": 1.0, "k_kill": 50.0, "k_clear": 1.0},
            {"k_tr": 2.0, "k_kill": 100.0, "k_clear": 1.0},
            {"k_tr": 5.0, "k_kill": 50.0, "k_clear": 0.5},
            {"k_tr": 10.0, "k_kill": 20.0, "k_clear": 0.2},
        ]
    else:
        starts = [
            {"k_tr": 0.2, "k_kill": 10.0, "k_clear": 0.2, "k_cyto": 1.0},
            {"k_tr": 0.5, "k_kill": 25.0, "k_clear": 0.5, "k_cyto": 10.0},
            {"k_tr": 1.0, "k_kill": 50.0, "k_clear": 1.0, "k_cyto": 100.0},
            {"k_tr": 2.0, "k_kill": 100.0, "k_clear": 1.0, "k_cyto": 10.0},
            {"k_tr": 5.0, "k_kill": 50.0, "k_clear": 0.5, "k_cyto": 1.0},
            {"k_tr": 10.0, "k_kill": 20.0, "k_clear": 0.2, "k_cyto": 100.0},
        ]
    for start in starts:
        start[DOSE_SCALING_PARAMETER_NAME] = 1.0
        start[BASELINE_DEATH_PARAMETER_NAME] = 0.02
        start[CONFLUENCE_DEATH_PARAMETER_NAME] = 0.05
        start["dose_gate_ec50_uM"] = 0.0125
        start["dose_gate_hill"] = 2.0
    return starts


def build_legacy_start_vector(
    start_values: Dict[str, float],
    treatment_names: Sequence[str],
    hill_names: Sequence[str],
) -> np.ndarray:
    ordered = [float(start_values[name]) for name in list(treatment_names) + list(hill_names)]
    return np.asarray(ordered, dtype=float)


def unpack_legacy_fit_vector(
    natural_params: Sequence[float],
    fit_config: JointFitConfig,
) -> Tuple[np.ndarray, JointFitConfig, Dict[str, float]]:
    params_arr = np.asarray(natural_params, dtype=float)
    treatment_names = get_treatment_parameter_names_for_config(fit_config)
    hill_names = get_legacy_hill_parameter_names_for_config(fit_config)
    n_treatment = len(treatment_names)
    treatment_params = np.asarray(params_arr[:n_treatment], dtype=float)
    dose_gate_params: Dict[str, float] = {
        "use_hill_dose_gate": bool(fit_config.use_hill_dose_gate),
        "fit_hill_dose_gate": bool(fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate),
        "dose_gate_ec50_uM": np.nan,
        "dose_gate_hill": np.nan,
    }
    if fit_config.use_hill_dose_gate:
        if hill_names:
            gate_ec50 = float(params_arr[n_treatment])
            gate_hill = float(params_arr[n_treatment + 1])
        else:
            gate_ec50 = float(fit_config.fixed_dose_gate_ec50_uM)
            gate_hill = float(fit_config.fixed_dose_gate_hill)
        dose_gate_params["dose_gate_ec50_uM"] = gate_ec50
        dose_gate_params["dose_gate_hill"] = gate_hill
    else:
        gate_ec50 = None
        gate_hill = None
    eval_config = build_legacy_hill_eval_config(fit_config, gate_ec50, gate_hill)
    return treatment_params, eval_config, dose_gate_params


def apply_legacy_summary_effective_dose_fields(
    summary: Dict[str, Any],
    fit_config: JointFitConfig,
    dose_gate_params: Dict[str, float],
) -> Dict[str, Any]:
    summary["use_hill_dose_gate"] = bool(fit_config.use_hill_dose_gate)
    summary["fit_hill_dose_gate"] = bool(fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate)
    summary["dose_gate_ec50_uM"] = dose_gate_params.get("dose_gate_ec50_uM", np.nan)
    summary["dose_gate_hill"] = dose_gate_params.get("dose_gate_hill", np.nan)
    hill_param_count = len(get_legacy_hill_parameter_names_for_config(fit_config))
    summary["n_parameters"] = int(summary.get("n_parameters", 0) + hill_param_count)
    if summary.get("nll") is not None and np.isfinite(summary.get("nll", np.nan)):
        nll = float(summary["nll"])
        n_obs = max(int(summary.get("n_observations", 0)), 1)
        n_parameters = int(summary["n_parameters"])
        summary["aic"] = float(2 * n_parameters + 2 * nll)
        summary["bic"] = float(n_parameters * np.log(n_obs) + 2 * nll)
    return summary


def fit_live_dead_model(
    dose_data_list,
    r_opt,
    K_opt,
    dfdctp_signal_curve,
    n_tr,
    objective="negative_binomial",
    fit_config: Optional[JointFitConfig] = None,
    observation_channels="alive_only",
    fit_means_only=True,
    high_dose_weight=1.0,
    max_nfev=3000,
):
    """
    Legacy two-stage fit: fits treatment parameters with baseline r/K held fixed.

    Available objectives:
    - `negative_binomial`
    - `poisson`
    - `least_squares` (legacy raw-residual objective)
    """
    if observation_channels not in {"alive_only", "alive_dead"}:
        raise ValueError(
            "observation_channels must be 'alive_only' or 'alive_dead', "
            f"got {observation_channels}"
        )
    fit_config = resolve_joint_fit_config(fit_config)
    model_variant = fit_config.model_variant
    treatment_names = get_treatment_parameter_names_for_config(fit_config)
    hill_names = get_legacy_hill_parameter_names_for_config(fit_config)
    start_dicts = build_legacy_start_dicts_for_config(fit_config)
    if objective == "least_squares":
        bounds = (
            [fit_config.lower_bounds[name] for name in list(treatment_names) + list(hill_names)],
            [fit_config.upper_bounds[name] for name in list(treatment_names) + list(hill_names)],
        )
        best_result = None
        best_summary = None
        best_cost = np.inf
        best_treatment_params = None
        best_dose_gate_params = {
            "use_hill_dose_gate": bool(fit_config.use_hill_dose_gate),
            "fit_hill_dose_gate": bool(fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate),
            "dose_gate_ec50_uM": (
                float(fit_config.fixed_dose_gate_ec50_uM) if fit_config.use_hill_dose_gate else np.nan
            ),
            "dose_gate_hill": (
                float(fit_config.fixed_dose_gate_hill) if fit_config.use_hill_dose_gate else np.nan
            ),
        }
        def residual_fn(full_params):
            treatment_params, eval_config, _ = unpack_legacy_fit_vector(full_params, fit_config)
            return residuals_global_joint(
                treatment_params,
                dose_data_list,
                r_opt,
                K_opt,
                dfdctp_signal_curve,
                n_tr,
                fit_means_only,
                high_dose_weight,
                observation_channels,
                model_variant,
                eval_config,
            )
        for start_values in start_dicts:
            guess_arr = build_legacy_start_vector(start_values, treatment_names, hill_names)
            res = least_squares(
                residual_fn,
                guess_arr,
                bounds=bounds,
                loss="linear",
                max_nfev=max_nfev,
            )
            if res.cost < best_cost:
                best_cost = float(res.cost)
                best_result = res
                best_treatment_params, eval_config, best_dose_gate_params = unpack_legacy_fit_vector(res.x, fit_config)
                best_summary = apply_legacy_summary_effective_dose_fields(summarize_live_dead_fit(
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_tr,
                    treatment_params=best_treatment_params,
                    objective=objective,
                    fit_config=eval_config,
                    observation_channels=observation_channels,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                ), fit_config, best_dose_gate_params)
        if best_result is None or best_summary is None:
            return None
        return {
            "success": bool(best_result.success),
            "treatment_params": np.asarray(best_treatment_params, dtype=float),
            "objective": objective,
            "summary": best_summary,
            "optimizer_result": best_result,
            "dose_gate_params": best_dose_gate_params,
        }

    if objective not in {"negative_binomial", "poisson"}:
        raise ValueError(f"Unsupported objective {objective}")

    theta_guess_grid = [(20.0, 20.0), (50.0, 50.0), (100.0, 30.0)]
    param_lower = np.array([fit_config.lower_bounds[name] for name in treatment_names], dtype=float)
    param_upper = np.array([fit_config.upper_bounds[name] for name in treatment_names], dtype=float)
    hill_lower = np.array([fit_config.lower_bounds[name] for name in hill_names], dtype=float)
    hill_upper = np.array([fit_config.upper_bounds[name] for name in hill_names], dtype=float)
    theta_lower = np.array([1e-3, 1e-3], dtype=float)
    theta_upper = np.array([1e6, 1e6], dtype=float)

    best_result = None
    best_summary = None
    best_nll = np.inf
    best_treatment_params = None
    best_dose_gate_params = {
        "use_hill_dose_gate": bool(fit_config.use_hill_dose_gate),
        "fit_hill_dose_gate": bool(fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate),
        "dose_gate_ec50_uM": (
            float(fit_config.fixed_dose_gate_ec50_uM) if fit_config.use_hill_dose_gate else np.nan
        ),
        "dose_gate_hill": (
            float(fit_config.fixed_dose_gate_hill) if fit_config.use_hill_dose_gate else np.nan
        ),
    }

    for start_values in start_dicts:
        theta_guesses = theta_guess_grid if objective == "negative_binomial" else [(None, None)]
        for theta_guess_alive, theta_guess_dead in theta_guesses:
            treatment_start = build_legacy_start_vector(start_values, treatment_names, hill_names)
            if objective == "negative_binomial":
                if observation_channels == "alive_dead":
                    x0 = np.log(np.concatenate([treatment_start, [theta_guess_alive, theta_guess_dead]]))
                    lower = np.log(np.concatenate([param_lower, hill_lower, theta_lower]))
                    upper = np.log(np.concatenate([param_upper, hill_upper, theta_upper]))
                else:
                    x0 = np.log(np.concatenate([treatment_start, [theta_guess_alive]]))
                    lower = np.log(np.concatenate([param_lower, hill_lower, theta_lower[:1]]))
                    upper = np.log(np.concatenate([param_upper, hill_upper, theta_upper[:1]]))
            else:
                x0 = np.log(treatment_start)
                lower = np.log(np.concatenate([param_lower, hill_lower]))
                upper = np.log(np.concatenate([param_upper, hill_upper]))

            def objective_fn(log_params):
                natural_params = np.exp(np.asarray(log_params, dtype=float))
                treatment_params, eval_config, _ = unpack_legacy_fit_vector(
                    natural_params[: len(treatment_names) + len(hill_names)],
                    fit_config,
                )
                if objective == "negative_binomial":
                    theta_alive_idx = len(treatment_names) + len(hill_names)
                    theta_alive = float(natural_params[theta_alive_idx])
                    theta_dead = float(natural_params[theta_alive_idx + 1]) if observation_channels == "alive_dead" else None
                else:
                    theta_alive = None
                    theta_dead = None
                return live_dead_objective_nll(
                    treatment_params=treatment_params,
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr_test=n_tr,
                    objective=objective,
                    fit_config=eval_config,
                    observation_channels=observation_channels,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                    theta_alive=theta_alive,
                    theta_dead=theta_dead,
                )

            res = minimize(
                objective_fn,
                x0=x0,
                method="L-BFGS-B",
                bounds=list(zip(lower, upper)),
                options={"maxiter": max_nfev},
            )
            if res.fun < best_nll:
                natural_params = np.exp(np.asarray(res.x, dtype=float))
                treatment_params, eval_config, best_dose_gate_params = unpack_legacy_fit_vector(
                    natural_params[: len(treatment_names) + len(hill_names)],
                    fit_config,
                )
                theta_alive = (
                    float(natural_params[len(treatment_names) + len(hill_names)])
                    if objective == "negative_binomial"
                    else None
                )
                theta_dead = (
                    float(natural_params[len(treatment_names) + len(hill_names) + 1])
                    if objective == "negative_binomial" and observation_channels == "alive_dead"
                    else None
                )
                best_nll = float(res.fun)
                best_result = res
                best_treatment_params = treatment_params
                best_summary = apply_legacy_summary_effective_dose_fields(summarize_live_dead_fit(
                    dose_data_list=dose_data_list,
                    r_fixed=r_opt,
                    K_fixed=K_opt,
                    dfdctp_signal_curve=dfdctp_signal_curve,
                    n_tr=n_tr,
                    treatment_params=treatment_params,
                    objective=objective,
                    fit_config=eval_config,
                    observation_channels=observation_channels,
                    fit_means_only=fit_means_only,
                    high_dose_weight=high_dose_weight,
                    theta_alive=theta_alive,
                    theta_dead=theta_dead,
                ), fit_config, best_dose_gate_params)

    if best_result is None or best_summary is None:
        return None

    return {
        "success": bool(best_result.success),
        "treatment_params": np.asarray(best_treatment_params, dtype=float),
        "objective": objective,
        "observation_channels": observation_channels,
        "summary": best_summary,
        "optimizer_result": best_result,
        "dose_gate_params": best_dose_gate_params,
    }


def plot_global_fit_subplots_joint(
    dose_data_list,
    ploidy_label,
    r,
    K,
    n_tr,
    dfdctp_signal_curve,
    k_tr,
    k_kill,
    k_clear,
    k_cyto=None,
    beta_dose=1.0,
    mu_base_death=0.0,
    mu_confluence_death=0.0,
    use_confluence_death: bool = False,
    confluence_death_exponent: float = 4.0,
    use_hill_dose_gate: bool = False,
    dose_gate_ec50_uM: Optional[float] = None,
    dose_gate_hill: Optional[float] = None,
    model_variant=JointFitConfig().model_variant,
    fit_config: Optional[JointFitConfig] = None,
    fit_summary: Optional[Dict[str, Any]] = None,
    output_dir: Optional[Path] = None,
    close_fig: bool = True,
):
    n_doses = len(dose_data_list)
    cols = 3
    rows = math.ceil(n_doses / cols)
    
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 5, rows * 4), sharex=True)
    if n_doses == 1: axes = [axes]
    elif n_doses > 1: axes = axes.flatten()
        
    model_variant = validate_model_variant(model_variant)
    plot_fit_config = resolve_joint_fit_config(fit_config) if fit_config is not None else JointFitConfig(
        model_variant=model_variant,
        use_hill_dose_gate=use_hill_dose_gate,
        fit_beta_dose=True,
        fixed_beta_dose=beta_dose,
        fit_hill_dose_gate=False,
        fixed_dose_gate_ec50_uM=(
            0.0125 if dose_gate_ec50_uM is None or not np.isfinite(dose_gate_ec50_uM)
            else float(dose_gate_ec50_uM)
        ),
        fixed_dose_gate_hill=(
            2.0 if dose_gate_hill is None or not np.isfinite(dose_gate_hill)
            else float(dose_gate_hill)
        ),
        use_confluence_death=use_confluence_death,
        fit_mu_confluence_death=False,
        fixed_mu_confluence_death=mu_confluence_death,
        confluence_death_exponent=confluence_death_exponent,
    )
    treatment_value_by_name = {
        "k_tr": k_tr,
        "k_kill": k_kill,
        "k_clear": k_clear,
        "k_cyto": 1e-8 if k_cyto is None else k_cyto,
        DOSE_SCALING_PARAMETER_NAME: beta_dose,
        BASELINE_DEATH_PARAMETER_NAME: mu_base_death,
        CONFLUENCE_DEATH_PARAMETER_NAME: mu_confluence_death,
    }
    treatment_params = [
        treatment_value_by_name[name]
        for name in get_treatment_parameter_names_for_config(plot_fit_config)
    ]
    
    for i, data in enumerate(dose_data_list):
        ax = axes[i]
        t_data = data['t']
        y_alive_raw = data['y_alive']
        y_dead_raw = data['y_dead']
        N0 = data['N0']
        D0 = data['D0']
        dose_muM = data['dose_muM']
        dose_label = data['dose_label']
        
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            A_sim, D_sim = simulate_joint_ext(
                t_data, treatment_params, N0, D0, r, K, dose_muM, dfdctp_signal_curve, n_tr,
                model_variant=plot_fit_config.model_variant,
                fit_config=plot_fit_config,
                beta_dose=beta_dose,
                use_confluence_death=use_confluence_death,
                mu_confluence_death=mu_confluence_death,
                confluence_death_exponent=confluence_death_exponent,
                use_hill_dose_gate=use_hill_dose_gate,
                dose_gate_ec50_uM=dose_gate_ec50_uM,
                dose_gate_hill=dose_gate_hill,
            )
            
        if y_alive_raw.ndim > 1:
            for j in range(y_alive_raw.shape[1]):
                ax.scatter(t_data, y_alive_raw[:, j], color='limegreen', alpha=0.3, s=20)
        else:
            ax.scatter(t_data, y_alive_raw, color='limegreen', alpha=0.5, s=20)
        ax.plot(t_data, A_sim, color='darkgreen', linewidth=2, label='Alive Model' if i==0 else "")
        
        if y_dead_raw.ndim > 1:
            for j in range(y_dead_raw.shape[1]):
                ax.scatter(t_data, y_dead_raw[:, j], color='lightcoral', alpha=0.3, s=20)
        else:
            ax.scatter(t_data, y_dead_raw, color='lightcoral', alpha=0.5, s=20)
        ax.plot(t_data, D_sim, color='darkred', linewidth=2, linestyle='--', label='Dead Model' if i==0 else "")
        
        ax.set_title(f"Dose: {dose_label}", fontsize=12, fontweight='bold')
        ax.grid(True, linestyle='--', alpha=0.6)
        
        if i % cols == 0: ax.set_ylabel("Cell Count", fontsize=11)
        if i >= (rows - 1) * cols: ax.set_xlabel("Time (Days)", fontsize=11)
        if i == 0: ax.legend(loc='upper left', fontsize=10)
            
    for j in range(n_doses, len(axes)): fig.delaxes(axes[j])
        
    param_lines = [
        f"Optima: n_tr = {n_tr}",
        f"Model: {fit_summary.get('model_variant', model_variant) if fit_summary is not None else model_variant}",
        f"$k_{{tr}}$ = {k_tr:.3f} | $k_{{kill}}$ = {k_kill:.3f} | $k_{{clear}}$ = {k_clear:.3f}",
    ]
    if validate_model_variant(model_variant) != "delayed_death_only" and k_cyto is not None:
        param_lines.append(f"$k_{{cyto}}$ = {k_cyto:.3f}")
    beta_dose = fit_summary.get("beta_dose") if fit_summary is not None else None
    if beta_dose is not None and np.isfinite(beta_dose):
        param_lines.append(f"$\\beta_{{dose}}$ = {beta_dose:.3f}")
    hill_ec50 = fit_summary.get("dose_gate_ec50_uM") if fit_summary is not None else None
    hill_h = fit_summary.get("dose_gate_hill") if fit_summary is not None else None
    if fit_summary is not None and fit_summary.get("use_hill_dose_gate"):
        if hill_ec50 is not None and np.isfinite(hill_ec50):
            param_lines.append(f"Hill EC50 = {hill_ec50 * 1000.0:.3f} nM")
        if hill_h is not None and np.isfinite(hill_h):
            param_lines.append(f"Hill h = {hill_h:.3f}")
    mu_base_death = fit_summary.get("mu_base_death") if fit_summary is not None else None
    if mu_base_death is not None and np.isfinite(mu_base_death):
        param_lines.append(f"$\\mu_{{base death}}$ = {mu_base_death:.4f} day$^{{-1}}$")
    mu_confluence_death = fit_summary.get("mu_confluence_death") if fit_summary is not None else None
    if mu_confluence_death is not None and np.isfinite(mu_confluence_death):
        param_lines.append(f"$\\mu_{{conf death}}$ = {mu_confluence_death:.4f} day$^{{-1}}$")
    if fit_summary is not None and fit_summary.get("use_confluence_death"):
        exponent = fit_summary.get("confluence_death_exponent")
        if exponent is not None and np.isfinite(exponent):
            param_lines.append(f"Confluence death exponent = {exponent:.2f}")
    if fit_summary is not None:
        param_lines.append(f"Fit objective: {fit_summary['objective']}")
        if fit_summary.get("observation_channels") is not None:
            param_lines.append(f"Observation channels: {fit_summary['observation_channels']}")
        if fit_summary.get("theta_alive") is not None and fit_summary.get("theta_dead") is not None:
            param_lines.append(
                f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f} | "
                f"$\\theta_{{dead}}$ = {fit_summary['theta_dead']:.3f}"
            )
        elif fit_summary.get("theta_alive") is not None:
            param_lines.append(f"$\\theta_{{alive}}$ = {fit_summary['theta_alive']:.3f}")
        if fit_summary.get("nll") is not None and fit_summary.get("aic") is not None and fit_summary.get("bic") is not None:
            param_lines.append(
                f"NLL = {fit_summary['nll']:.3f} | AIC = {fit_summary['aic']:.3f} | "
                f"BIC = {fit_summary['bic']:.3f}"
            )
        elif fit_summary.get("nll") is not None:
            param_lines.append(f"Data NLL = {fit_summary['nll']:.3f}")
        else:
            param_lines.append(f"Legacy raw-count least-squares cost = {fit_summary['objective_value']:.3f}")
        if np.isfinite(fit_summary.get("rmse", np.nan)):
            param_lines.append(f"RMSE = {fit_summary['rmse']:.3f}")
    param_text = "\n".join(param_lines)
    fig.suptitle(f"{ploidy_label} Cohort - Joint Live/Dead Kinetic Fit\n{param_text}", fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.85)

    if output_dir is not None:
        output_path = output_dir / f"cohort_joint_fit_{slugify_label(ploidy_label)}.png"
        save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')
    elif close_fig:
        plt.close(fig)

def get_fitting_data_one_row(
    df,
    gem_dose: str,
    ploidy: str,
    phenotype: str,
    plate_row: str,
    count_transitional_as_alive: bool = False,
):
    """
    Backward-compatible wrapper returning one phenotype vector after explicit
    Alive/Dead alignment for a single row replicate.

    Example:
        plate_row='A' for 2N replicate A
        plate_row='E' for 4N replicate E
    """
    aligned = get_aligned_live_dead_data(
        df,
        gem_dose=gem_dose,
        ploidy=ploidy,
        plate_row=plate_row,
        count_transitional_as_alive=count_transitional_as_alive,
    )
    matrix = aligned["y_alive"] if phenotype == "Alive" else aligned["y_dead"] if phenotype == "Dead" else None
    if matrix is None:
        raise ValueError(f"Unsupported phenotype {phenotype}; expected Alive or Dead")
    if matrix.shape[1] != 1:
        raise ValueError(
            f"Expected exactly one aligned replicate column for {gem_dose}, {ploidy}, row {plate_row}; "
            f"got {matrix.shape[1]}"
        )
    return aligned["t"], matrix[:, 0]

def fit_joint_one_replicate(
    dose_data_list,
    r_opt,
    K_opt,
    dfdctp_signal_curve,
    ploidy,
    rep_idx,
    objective="negative_binomial",
    observation_channels="alive_only",
    fit_config: Optional[JointFitConfig] = None,
):
    """
    Fits one global joint live/dead model across all doses for one replicate.
    """
    fit_config = resolve_joint_fit_config(fit_config)
    best_objective_value, best_fit, best_n_tr = np.inf, None, 1
    n_range = range(2, 10)

    for n_test in n_range:
        fit_result = fit_live_dead_model(
            dose_data_list=dose_data_list,
            r_opt=r_opt,
            K_opt=K_opt,
            dfdctp_signal_curve=dfdctp_signal_curve,
            n_tr=n_test,
            objective=objective,
            fit_config=fit_config,
            observation_channels=observation_channels,
            fit_means_only=False,
            high_dose_weight=1.0,
            max_nfev=3000,
        )
        if fit_result is None:
            continue
        objective_value = float(fit_result["summary"]["objective_value"])
        if objective_value < best_objective_value:
            best_objective_value, best_fit, best_n_tr = objective_value, fit_result, n_test

    if best_fit is None:
        return None

    fitted_treatment = unpack_treatment_params_for_config(best_fit["treatment_params"], fit_config)
    k_tr_opt = fitted_treatment["k_tr"]
    k_kill_opt = fitted_treatment["k_kill"]
    k_clear_opt = fitted_treatment["k_clear"]
    fit_summary = best_fit["summary"]
    
    print(f"\n--- {ploidy} Replicate {rep_idx} ---")
    print(f"Optimal n_tr: {best_n_tr}")
    print(f"Optimal k_tr: {k_tr_opt:.4f}")
    print(f"Optimal k_kill (Potency): {k_kill_opt:.4f}")
    print(f"Optimal k_clear: {k_clear_opt:.4f}")
    if "k_cyto" in fitted_treatment:
        print(f"Optimal k_cyto: {fitted_treatment['k_cyto']:.4f}")
    print(f"Fit objective: {objective}")
    print(f"Observation channels: {observation_channels}")
    if fit_summary.get("use_hill_dose_gate"):
        print(
            "Hill gate: "
            f"EC50={fit_summary.get('dose_gate_ec50_uM', np.nan):.6f} uM, "
            f"h={fit_summary.get('dose_gate_hill', np.nan):.4f} "
            f"(fit={fit_summary.get('fit_hill_dose_gate', False)})"
        )
    print(
        f"beta_dose: {fit_summary.get('beta_dose', fit_config.fixed_beta_dose):.4f} "
        f"(fit={fit_config.fit_beta_dose})"
    )
    if fit_summary.get("use_confluence_death"):
        print(
            "Confluence death: "
            f"mu={fit_summary.get('mu_confluence_death', fit_config.fixed_mu_confluence_death):.4f}, "
            f"exponent={fit_summary.get('confluence_death_exponent', fit_config.confluence_death_exponent):.2f}"
        )
    if fit_summary["nll"] is not None:
        print(f"Best negative log-likelihood: {fit_summary['nll']:.4f}")
        print(f"AIC: {fit_summary['aic']:.4f}")
        print(f"BIC: {fit_summary['bic']:.4f}")
        print(f"Alive dispersion theta: {fit_summary['theta_alive']:.4f}")
        if fit_summary["theta_dead"] is not None:
            print(f"Dead dispersion theta: {fit_summary['theta_dead']:.4f}")
    else:
        print(f"Best raw-count least-squares cost: {fit_summary['objective_value']:.4f}")
    print(f"RMSE: {fit_summary['rmse']:.4f}")

    return {
        'ploidy': ploidy,
        'rep_idx': rep_idx,
        'n_tr': best_n_tr,
        'k_tr': k_tr_opt,
        'k_kill': k_kill_opt,
        'k_clear': k_clear_opt,
        'k_cyto': fitted_treatment.get('k_cyto', np.nan),
        'cost': best_objective_value,
        'objective': objective,
        'observation_channels': observation_channels,
        'theta_alive': fit_summary['theta_alive'],
        'theta_dead': fit_summary['theta_dead'],
        'nll': fit_summary['nll'],
        'aic': fit_summary['aic'],
        'bic': fit_summary['bic'],
        'rmse': fit_summary['rmse'],
    }

def build_row_replicate_dose_data_list(
    df,
    gem_doses,
    ploidy,
    plate_row,
    t_max,
    count_transitional_as_alive: bool = False,
):
    """
    Builds a dose_data_list for one row replicate across all gemcitabine doses.
    """
    dose_data_list = []

    for gem_dose in gem_doses:
        try:
            aligned = get_aligned_live_dead_data(
                df,
                gem_dose=gem_dose,
                ploidy=ploidy,
                plate_row=plate_row,
                t_max=t_max,
                count_transitional_as_alive=count_transitional_as_alive,
            )
            if aligned["dropped_timepoints"] > 0:
                print(
                    f"  {gem_dose}, row {plate_row}: dropped {aligned['dropped_timepoints']} "
                    "non-overlapping Alive/Dead time points"
                )
            if aligned["y_alive"].shape[1] != 1:
                raise ValueError(
                    f"Expected one aligned replicate column for {gem_dose}, {ploidy}, row {plate_row}; "
                    f"got {aligned['y_alive'].shape[1]}"
                )
            t_fit = aligned["t"]
            y_alive_fit = aligned["y_alive"][:, 0]
            y_dead_fit = aligned["y_dead"][:, 0]
            
            dose_muM = float(gem_dose.split()[0]) / 1000.0
            N0 = y_alive_fit[0]
            D0 = y_dead_fit[0]
            
            dose_data_list.append({
                'dose_label': gem_dose,
                'dose_muM': dose_muM,
                't': t_fit,
                'y_alive': y_alive_fit,
                'y_dead': y_dead_fit,
                'N0': N0,
                'D0': D0,
                'plate_row': plate_row
            })
            
        except ValueError as e:
            print(f"  Skipping {gem_dose}, row {plate_row}: {e}")

    return dose_data_list


def build_joint_fit_trajectories(
    df: pd.DataFrame,
    ploidies: Sequence[str],
    gem_doses: Sequence[str],
    fit_t_max: Optional[float],
    count_transitional_as_alive: bool = False,
    max_replicates_per_ploidy_dose: Optional[int] = None,
) -> List[ReplicateTrajectory]:
    trajectories: List[ReplicateTrajectory] = []
    total_dropped_nonfinite = 0
    for ploidy in ploidies:
        for gem_dose in gem_doses:
            aligned = get_aligned_live_dead_data(
                df,
                gem_dose=gem_dose,
                ploidy=ploidy,
                t_max=fit_t_max,
                count_transitional_as_alive=count_transitional_as_alive,
            )
            replicate_items = list(enumerate(aligned["replicate_columns"]))
            if max_replicates_per_ploidy_dose is not None:
                replicate_items = replicate_items[:int(max_replicates_per_ploidy_dose)]
            for rep_idx, replicate_id in replicate_items:
                t_rep, alive, dead, dropped_nonfinite = trim_finite_live_dead_observations(
                    aligned["t"],
                    aligned["y_alive"][:, rep_idx],
                    aligned["y_dead"][:, rep_idx],
                )
                total_dropped_nonfinite += dropped_nonfinite
                trajectories.append(
                    ReplicateTrajectory(
                        ploidy=ploidy,
                        dose_label=gem_dose,
                        dose_uM=float(gem_dose.split()[0]) / 1000.0,
                        replicate_id=(str(replicate_id[0]), int(replicate_id[1])),
                        t=t_rep,
                        alive=alive,
                        dead=dead,
                        N0=float(alive[0]),
                        D0=float(dead[0]),
                    )
                )
    if total_dropped_nonfinite > 0:
        print(
            f"Joint fit trajectories: dropped {total_dropped_nonfinite} non-finite paired "
            "alive/dead observations before likelihood evaluation."
        )
    return trajectories


def fit_joint_partial_pooling_model(
    trajectories: Sequence[ReplicateTrajectory],
    dfdctp_signal_curve_by_ploidy: Dict[str, DfdctpSignalSurface],
    fit_config: JointFitConfig,
    n_tr: int,
    diagnostic_output_dir: Optional[Path] = None,
    diagnostic_label: Optional[str] = None,
    start_indices: Optional[Sequence[int]] = None,
    diagnostic_tasks: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """
    Primary MAP fit with log-scale partial pooling across ploidies.
    """
    if len(trajectories) == 0:
        return {
            "success": False,
            "posterior_objective": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "optimizer_message": "No trajectories provided",
            "optimizer_attempts": pd.DataFrame(),
            "model_variant": fit_config.model_variant,
            "dose_gate_params": {},
        }
    if fit_config.objective in {"negative_binomial", "poisson"} and fit_config.fit_means_only:
        raise ValueError("fit_means_only=True is not allowed with count-likelihood objectives")
    if fit_config.observation_channels == "alive_only":
        raise ValueError("Partial-pooling primary fit requires observation_channels='alive_dead' to identify k_clear")
    parameter_names = get_parameter_names_for_config(fit_config)

    ploidies = sorted({traj.ploidy for traj in trajectories})
    if any(ploidy not in dfdctp_signal_curve_by_ploidy for ploidy in ploidies):
        missing = [ploidy for ploidy in ploidies if ploidy not in dfdctp_signal_curve_by_ploidy]
        raise ValueError(f"Missing dFdCTP surfaces for ploidies: {missing}")

    alive_max = max(float(np.nanmax(traj.alive)) for traj in trajectories if len(traj.alive) > 0)
    base_start_grid = [
        {
            "r": 0.7,
            "K": max(alive_max, 1.0),
            "k_tr": 0.5,
            "k_kill": 25.0,
            "k_clear": 0.5,
            "beta_dose": 1.0,
            "mu_base_death": 0.02,
            "mu_confluence_death": 0.05,
            "dose_gate_ec50_uM": 0.0125,
            "dose_gate_hill": 2.0,
            "theta_alive": 20.0,
            "theta_dead": 20.0,
        },
        {
            "r": 1.2,
            "K": max(alive_max * 1.5, 1.0),
            "k_tr": 1.0,
            "k_kill": 50.0,
            "k_clear": 1.0,
            "beta_dose": 1.0,
            "mu_base_death": 0.02,
            "mu_confluence_death": 0.05,
            "dose_gate_ec50_uM": 0.0125,
            "dose_gate_hill": 2.0,
            "theta_alive": 50.0,
            "theta_dead": 50.0,
        },
        {
            "r": 0.3,
            "K": max(alive_max * 2.0, 1.0),
            "k_tr": 2.0,
            "k_kill": 100.0,
            "k_clear": 0.2,
            "beta_dose": 1.0,
            "mu_base_death": 0.02,
            "mu_confluence_death": 0.05,
            "dose_gate_ec50_uM": 0.0125,
            "dose_gate_hill": 2.0,
            "theta_alive": 100.0,
            "theta_dead": 30.0,
        },
    ]
    if fit_config.model_variant == "delayed_death_only":
        mu_start_grid = base_start_grid
    else:
        # Conservative 4-start scheme: retain low/mid/high baseline/treatment
        # regimes while reducing the previous 3x3 k_cyto expansion.
        start_templates = [
            (base_start_grid[0], 1.0),
            (base_start_grid[0], 100.0),
            (base_start_grid[1], 10.0),
            (base_start_grid[2], 10.0),
        ]
        mu_start_grid = []
        for start_values, k_cyto in start_templates:
            expanded = dict(start_values)
            expanded["k_cyto"] = k_cyto
            mu_start_grid.append(expanded)

    prior_sd_by_param = {
        "r": fit_config.prior_sd_log_r,
        "K": fit_config.prior_sd_log_K,
        "k_tr": fit_config.prior_sd_log_k_tr,
        "k_kill": fit_config.prior_sd_log_k_kill,
        "k_clear": fit_config.prior_sd_log_k_clear,
        "mu_base_death": fit_config.prior_sd_log_mu_base_death,
    }
    if fit_config.use_confluence_death and fit_config.fit_mu_confluence_death:
        prior_sd_by_param[CONFLUENCE_DEATH_PARAMETER_NAME] = fit_config.prior_sd_log_mu_confluence_death
    if fit_config.fit_beta_dose:
        prior_sd_by_param["beta_dose"] = fit_config.prior_sd_log_beta_dose
    if "k_cyto" in parameter_names:
        prior_sd_by_param["k_cyto"] = fit_config.prior_sd_log_k_cyto

    def unpack_parameters(log_x: np.ndarray) -> Dict[str, Any]:
        idx = 0
        mu_logs = {}
        for name in parameter_names:
            mu_logs[name] = float(log_x[idx])
            idx += 1
        delta_logs = {ploidy: {} for ploidy in ploidies}
        for ploidy in ploidies:
            for name in parameter_names:
                delta_logs[ploidy][name] = float(log_x[idx])
                idx += 1
        dose_gate_params = {
            "use_hill_dose_gate": bool(fit_config.use_hill_dose_gate),
            "fit_hill_dose_gate": bool(fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate),
            "dose_gate_ec50_uM": np.nan,
            "dose_gate_hill": np.nan,
        }
        if fit_config.use_hill_dose_gate:
            if fit_config.fit_hill_dose_gate:
                dose_gate_params["dose_gate_ec50_uM"] = float(np.exp(log_x[idx]))
                idx += 1
                dose_gate_params["dose_gate_hill"] = float(np.exp(log_x[idx]))
                idx += 1
            else:
                dose_gate_params["dose_gate_ec50_uM"] = float(fit_config.fixed_dose_gate_ec50_uM)
                dose_gate_params["dose_gate_hill"] = float(fit_config.fixed_dose_gate_hill)
        theta_alive = None
        theta_dead = None
        if fit_config.objective == "negative_binomial":
            theta_alive = float(np.exp(log_x[idx]))
            idx += 1
            theta_dead = float(np.exp(log_x[idx]))
            idx += 1

        ploidy_params = {}
        for ploidy in ploidies:
            params = {
                name: float(np.exp(mu_logs[name] + delta_logs[ploidy][name]))
                for name in parameter_names
            }
            if DOSE_SCALING_PARAMETER_NAME not in params:
                params[DOSE_SCALING_PARAMETER_NAME] = float(fit_config.fixed_beta_dose)
            if CONFLUENCE_DEATH_PARAMETER_NAME not in params:
                params[CONFLUENCE_DEATH_PARAMETER_NAME] = (
                    float(fit_config.fixed_mu_confluence_death)
                    if fit_config.use_confluence_death else 0.0
                )
            ploidy_params[ploidy] = params
        return {
            "mu_logs": mu_logs,
            "delta_logs": delta_logs,
            "ploidy_params": ploidy_params,
            "dose_gate_params": dose_gate_params,
            "theta_alive": theta_alive,
            "theta_dead": theta_dead,
        }

    def data_nll_and_penalty(
        log_x: np.ndarray,
        debug: bool = False,
    ) -> Union[Tuple[float, float], Tuple[float, float, Dict[str, Any]]]:
        unpacked = unpack_parameters(log_x)
        prior_penalty = 0.0
        for ploidy in ploidies:
            for name in parameter_names:
                prior_penalty += 0.5 * (unpacked["delta_logs"][ploidy][name] / prior_sd_by_param[name]) ** 2
        if fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate:
            log_ec50 = np.log(float(unpacked["dose_gate_params"]["dose_gate_ec50_uM"]))
            log_hill = np.log(float(unpacked["dose_gate_params"]["dose_gate_hill"]))
            prior_penalty += 0.5 * (
                (log_ec50 - np.log(fit_config.fixed_dose_gate_ec50_uM))
                / fit_config.prior_sd_log_dose_gate_ec50_uM
            ) ** 2
            prior_penalty += 0.5 * (
                (log_hill - np.log(fit_config.fixed_dose_gate_hill))
                / fit_config.prior_sd_log_dose_gate_hill
            ) ** 2

        data_nll = 0.0
        debug_info = {
            "simulations_attempted": 0,
            "simulations_succeeded": 0,
            "simulations_failed": 0,
            "failed_trajectories": [],
            "failure_messages": {},
            "penalty_hit": False,
        }
        for traj in trajectories:
            debug_info["simulations_attempted"] += 1
            params = unpacked["ploidy_params"][traj.ploidy]
            beta = resolve_beta_dose_for_params(params, fit_config)
            mu_conf = resolve_mu_confluence_death_for_params(params, fit_config)
            gate_params = resolve_hill_dose_gate_params_for_unpacked(unpacked, fit_config)
            sim = simulate_joint_dfdctp_safe(
                t=traj.t,
                params=[params[name] for name in get_treatment_parameter_names_for_config(fit_config)],
                N0=traj.N0,
                D0=traj.D0,
                r=params["r"],
                K=params["K"],
                dose_muM=traj.dose_uM,
                dfdctp_signal_curve=dfdctp_signal_curve_by_ploidy[traj.ploidy],
                n_tr=n_tr,
                model_variant=fit_config.model_variant,
                fit_config=fit_config,
                beta_dose=beta,
                use_confluence_death=fit_config.use_confluence_death,
                mu_confluence_death=mu_conf,
                confluence_death_exponent=fit_config.confluence_death_exponent,
                use_hill_dose_gate=gate_params["use_hill_dose_gate"],
                dose_gate_ec50_uM=gate_params["dose_gate_ec50_uM"] if gate_params["use_hill_dose_gate"] else None,
                dose_gate_hill=gate_params["dose_gate_hill"] if gate_params["use_hill_dose_gate"] else None,
            )
            if not sim.success:
                debug_info["simulations_failed"] += 1
                debug_info["failed_trajectories"].append(f"{traj.ploidy}:{traj.dose_label}:{traj.replicate_id}")
                debug_info["failure_messages"][sim.message] = debug_info["failure_messages"].get(sim.message, 0) + 1
                raise RuntimeError(f"Simulation failed for {traj.ploidy} {traj.dose_label} {traj.replicate_id}: {sim.message}")
            debug_info["simulations_succeeded"] += 1
            alive_mask = np.isfinite(traj.alive) & np.isfinite(sim.alive)
            dead_mask = np.isfinite(traj.dead) & np.isfinite(sim.dead)
            if not np.any(alive_mask):
                debug_info["failure_messages"]["No finite alive observations remain"] = (
                    debug_info["failure_messages"].get("No finite alive observations remain", 0) + 1
                )
                raise RuntimeError(
                    f"No finite alive observations remain for {traj.ploidy} {traj.dose_label} {traj.replicate_id}"
                )
            if not np.any(dead_mask):
                debug_info["failure_messages"]["No finite dead observations remain"] = (
                    debug_info["failure_messages"].get("No finite dead observations remain", 0) + 1
                )
                raise RuntimeError(
                    f"No finite dead observations remain for {traj.ploidy} {traj.dose_label} {traj.replicate_id}"
                )
            if fit_config.objective == "negative_binomial":
                data_nll += negative_binomial_nll(
                    traj.alive[alive_mask],
                    sim.alive[alive_mask],
                    unpacked["theta_alive"],
                )
                data_nll += negative_binomial_nll(
                    traj.dead[dead_mask],
                    sim.dead[dead_mask],
                    unpacked["theta_dead"],
                )
            elif fit_config.objective == "poisson":
                data_nll += poisson_nll(traj.alive[alive_mask], sim.alive[alive_mask])
                data_nll += poisson_nll(traj.dead[dead_mask], sim.dead[dead_mask])
            else:
                data_nll += 0.5 * float(
                    np.sum((traj.alive[alive_mask] - sim.alive[alive_mask]) ** 2)
                    + np.sum((traj.dead[dead_mask] - sim.dead[dead_mask]) ** 2)
                )
        if debug:
            debug_info["failure_messages_joined"] = "; ".join(
                f"{key} x{count}" for key, count in sorted(debug_info["failure_messages"].items())
            )
            return float(data_nll), float(prior_penalty), debug_info
        return float(data_nll), float(prior_penalty)

    def posterior_objective(log_x: np.ndarray) -> float:
        data_nll, prior_penalty = data_nll_and_penalty(log_x)
        return data_nll + prior_penalty

    def debug_objective_evaluation(log_x: np.ndarray) -> Dict[str, Any]:
        out = {
            "penalty_hit": False,
            "objective_value": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "simulations_attempted": 0,
            "simulations_succeeded": 0,
            "simulations_failed": 0,
            "failure_messages_joined": "",
        }
        try:
            data_nll, prior_penalty, debug_info = data_nll_and_penalty(log_x, debug=True)
            out.update(debug_info)
            out["data_nll"] = float(data_nll)
            out["prior_penalty"] = float(prior_penalty)
            out["objective_value"] = float(data_nll + prior_penalty)
            out["penalty_hit"] = not np.isfinite(out["objective_value"]) or out["objective_value"] >= fit_config.large_objective_penalty
        except Exception as exc:
            out["penalty_hit"] = True
            out["failure_messages_joined"] = str(exc)
        return out

    def build_start_vector(start_values: Dict[str, float]) -> np.ndarray:
        values: List[float] = []
        for name in parameter_names:
            values.append(np.log(start_values[name]))
        for _ploidy in ploidies:
            for _name in parameter_names:
                values.append(0.0)
        if fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate:
            values.extend(
                [
                    np.log(start_values["dose_gate_ec50_uM"]),
                    np.log(start_values["dose_gate_hill"]),
                ]
            )
        if fit_config.objective == "negative_binomial":
            values.extend([np.log(start_values["theta_alive"]), np.log(start_values["theta_dead"])])
        return np.asarray(values, dtype=float)

    def add_unpacked_attempt_parameters(row: Dict[str, Any], prefix: str, log_x: np.ndarray) -> None:
        try:
            unpacked = unpack_parameters(np.asarray(log_x, dtype=float))
        except Exception:
            for name in ALL_DYNAMIC_PARAMETER_NAMES:
                row[f"{prefix}_mu_{name}"] = np.nan
            for ploidy in ["2N", "4N"]:
                for name in ALL_DYNAMIC_PARAMETER_NAMES:
                    row[f"{prefix}_{ploidy}_{name}"] = np.nan
            row[f"{prefix}_dose_gate_ec50_uM"] = np.nan
            row[f"{prefix}_dose_gate_hill"] = np.nan
            row[f"{prefix}_fit_beta_dose"] = fit_config.fit_beta_dose
            row[f"{prefix}_use_hill_dose_gate"] = fit_config.use_hill_dose_gate
            row[f"{prefix}_fit_hill_dose_gate"] = fit_config.fit_hill_dose_gate
            row[f"{prefix}_theta_alive"] = np.nan
            row[f"{prefix}_theta_dead"] = np.nan
            return

        mu_logs = unpacked["mu_logs"]
        for name in ALL_DYNAMIC_PARAMETER_NAMES:
            if name == DOSE_SCALING_PARAMETER_NAME and not fit_config.fit_beta_dose:
                row[f"{prefix}_mu_{name}"] = float(fit_config.fixed_beta_dose)
            elif name == CONFLUENCE_DEATH_PARAMETER_NAME and not (
                fit_config.use_confluence_death and fit_config.fit_mu_confluence_death
            ):
                row[f"{prefix}_mu_{name}"] = float(
                    fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0
                )
            else:
                row[f"{prefix}_mu_{name}"] = float(np.exp(mu_logs[name])) if name in mu_logs else np.nan
        ploidy_parameters = unpacked["ploidy_params"]
        for ploidy in ["2N", "4N"]:
            params = ploidy_parameters.get(ploidy, {})
            for name in ALL_DYNAMIC_PARAMETER_NAMES:
                if name == DOSE_SCALING_PARAMETER_NAME and not fit_config.fit_beta_dose:
                    row[f"{prefix}_{ploidy}_{name}"] = float(fit_config.fixed_beta_dose)
                elif name == CONFLUENCE_DEATH_PARAMETER_NAME and not (
                    fit_config.use_confluence_death and fit_config.fit_mu_confluence_death
                ):
                    row[f"{prefix}_{ploidy}_{name}"] = float(
                        fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0
                    )
                else:
                    row[f"{prefix}_{ploidy}_{name}"] = params.get(name, np.nan)
        dose_gate_params = unpacked.get("dose_gate_params", {})
        row[f"{prefix}_dose_gate_ec50_uM"] = dose_gate_params.get("dose_gate_ec50_uM", np.nan)
        row[f"{prefix}_dose_gate_hill"] = dose_gate_params.get("dose_gate_hill", np.nan)
        row[f"{prefix}_model_preset"] = fit_config.model_preset
        row[f"{prefix}_fit_beta_dose"] = fit_config.fit_beta_dose
        row[f"{prefix}_use_hill_dose_gate"] = fit_config.use_hill_dose_gate
        row[f"{prefix}_fit_hill_dose_gate"] = fit_config.fit_hill_dose_gate
        row[f"{prefix}_theta_alive"] = unpacked["theta_alive"] if unpacked["theta_alive"] is not None else np.nan
        row[f"{prefix}_theta_dead"] = unpacked["theta_dead"] if unpacked["theta_dead"] is not None else np.nan

    def gradient_norm_from_result(result: Any) -> float:
        try:
            jac = getattr(result, "jac", None)
            if jac is None:
                return np.nan
            return float(np.linalg.norm(np.asarray(jac, dtype=float)))
        except Exception:
            return np.nan

    lower_bounds = []
    upper_bounds = []
    for name in parameter_names:
        lower_bounds.append(np.log(fit_config.lower_bounds[name]))
        upper_bounds.append(np.log(fit_config.upper_bounds[name]))
    for _ploidy in ploidies:
        for _name in parameter_names:
            lower_bounds.append(-5.0)
            upper_bounds.append(5.0)
    if fit_config.use_hill_dose_gate and fit_config.fit_hill_dose_gate:
        lower_bounds.extend(
            [
                np.log(fit_config.lower_bounds["dose_gate_ec50_uM"]),
                np.log(fit_config.lower_bounds["dose_gate_hill"]),
            ]
        )
        upper_bounds.extend(
            [
                np.log(fit_config.upper_bounds["dose_gate_ec50_uM"]),
                np.log(fit_config.upper_bounds["dose_gate_hill"]),
            ]
        )
    if fit_config.objective == "negative_binomial":
        lower_bounds.extend([np.log(fit_config.lower_bounds["theta_alive"]), np.log(fit_config.lower_bounds["theta_dead"])])
        upper_bounds.extend([np.log(fit_config.upper_bounds["theta_alive"]), np.log(fit_config.upper_bounds["theta_dead"])])
    bounds = list(zip(lower_bounds, upper_bounds))
    coordinate_names = get_partial_pooling_coordinate_names(parameter_names, ploidies, fit_config.objective, fit_config=fit_config)

    attempt_rows: List[Dict[str, Any]] = []
    best_result = None
    best_unpacked = None
    best_data_nll = np.nan
    best_prior_penalty = np.nan
    best_value = fit_config.large_objective_penalty
    n_observations = int(
        sum(np.count_nonzero(np.isfinite(traj.alive)) + np.count_nonzero(np.isfinite(traj.dead)) for traj in trajectories)
    )
    diagnostic_task_set = None if diagnostic_tasks is None else set(diagnostic_tasks)

    allowed_start_indices = None if start_indices is None else {int(idx) for idx in start_indices}

    for start_idx, start_values in enumerate(mu_start_grid):
        if allowed_start_indices is not None and start_idx not in allowed_start_indices:
            continue
        x0 = build_start_vector(start_values)

        def safe_posterior(x):
            return safe_objective_value(posterior_objective, np.asarray(x, dtype=float), fit_config.large_objective_penalty)

        initial_raw_objective = safe_posterior(x0)
        if np.isfinite(initial_raw_objective) and initial_raw_objective < fit_config.large_objective_penalty:
            objective_scale = max(abs(float(initial_raw_objective)), 1.0)
        else:
            objective_scale = 1.0

        def optimizer_objective(x):
            return safe_posterior(x) / objective_scale

        initial_objective = float(initial_raw_objective)
        initial_normalized_objective = float(optimizer_objective(x0))
        initial_debug = debug_objective_evaluation(x0)
        if diagnostic_output_dir is not None and diagnostic_label is not None and start_idx == 0:
            if diagnostic_task_set is None or "local_sensitivity" in diagnostic_task_set:
                diagnose_objective_near_start(
                    posterior_objective,
                    x0,
                    lower_bounds,
                    upper_bounds,
                    unpack_parameters,
                    coordinate_names,
                    f"{diagnostic_label}__ntr{n_tr}__start{start_idx}",
                    diagnostic_output_dir / "objective_local_sensitivity.tsv",
                    fit_config.large_objective_penalty,
                    debug_eval_fn=debug_objective_evaluation,
                )
            if diagnostic_task_set is None or "parameter_unpacking" in diagnostic_task_set:
                build_parameter_unpacking_check(
                    x0,
                    lower_bounds,
                    upper_bounds,
                    unpack_parameters,
                    coordinate_names,
                    diagnostic_output_dir / "parameter_unpacking_check.tsv",
                )
            coord_probe = None
            if diagnostic_task_set is None or "manual_coordinate_descent" in diagnostic_task_set:
                coord_probe = manual_coordinate_descent_probe(
                    posterior_objective,
                    x0,
                    lower_bounds,
                    upper_bounds,
                    coordinate_names,
                    diagnostic_output_dir / "manual_coordinate_descent.tsv",
                    fit_config.large_objective_penalty,
                    diagnostic_label,
                )
            if diagnostic_task_set == {"normalized_lbfgsb"}:
                diagnostic_normalized_objective_lbfgsb(
                    posterior_objective,
                    x0,
                    bounds,
                    diagnostic_output_dir / "diagnostic_normalized_objective_lbfgsb.tsv",
                    maxiter=20,
                    n_observations=n_observations,
                )
                return {
                    "success": False,
                    "posterior_objective": fit_config.large_objective_penalty,
                    "data_nll": np.nan,
                    "prior_penalty": np.nan,
                    "n_observations": n_observations,
                    "objective": fit_config.objective,
                    "observation_channels": fit_config.observation_channels,
                    "model_variant": fit_config.model_variant,
                    "theta_alive": np.nan,
                    "theta_dead": np.nan,
                    "optimizer_message": "Diagnostic-only normalized-objective run",
                    "optimizer_attempts": pd.DataFrame(),
                    "population_parameters": {},
                    "ploidy_parameters": {},
                    "dose_gate_params": {},
                }
        try:
            result = minimize(
                optimizer_objective,
                x0=x0,
                method=fit_config.optimizer_method,
                bounds=bounds,
                options={"maxiter": fit_config.max_nfev},
            )
            final_value = safe_posterior(result.x)
            final_normalized_objective = float(optimizer_objective(result.x))
            if np.isfinite(final_value) and final_value < fit_config.large_objective_penalty:
                data_nll, prior_penalty = data_nll_and_penalty(result.x)
            else:
                data_nll, prior_penalty = np.nan, np.nan
            jac = np.asarray(getattr(result, "jac", np.array([])), dtype=float)
            moved = np.allclose(np.asarray(result.x, dtype=float), x0, rtol=0.0, atol=1e-12)
            max_abs_grad = float(np.max(np.abs(jac))) if jac.size > 0 else np.nan
            max_abs_grad_idx = int(np.argmax(np.abs(jac))) if jac.size > 0 else -1
            bound_flags = coordinate_bound_flags(result.x, lower_bounds, upper_bounds, coordinate_names)
            final_debug = debug_objective_evaluation(result.x)
            fit_is_usable = bool(
                (result.success or fit_config.accept_incomplete_optimizer)
                and np.isfinite(final_value)
                and final_value < fit_config.large_objective_penalty
            )
            if diagnostic_output_dir is not None and diagnostic_label is not None and start_idx == 0:
                manual_gradients = None
                if diagnostic_task_set is None or "manual_gradient" in diagnostic_task_set or "gradient_comparison" in diagnostic_task_set or "line_search_gradient" in diagnostic_task_set:
                    _, manual_gradients = manual_gradient_probe(
                        posterior_objective,
                        x0,
                        lower_bounds,
                        upper_bounds,
                        coordinate_names,
                        diagnostic_output_dir / "manual_gradient.tsv",
                        fit_config.large_objective_penalty,
                    )
                if manual_gradients is not None and (diagnostic_task_set is None or "gradient_comparison" in diagnostic_task_set):
                    gradient_comparison_table(
                        manual_gradients,
                        result,
                        coordinate_names,
                        diagnostic_output_dir / "gradient_comparison.tsv",
                    )
                if manual_gradients is not None and (diagnostic_task_set is None or "line_search_gradient" in diagnostic_task_set):
                    manual_grad = manual_gradients[1e-4]
                    grad_norm = float(np.linalg.norm(manual_grad))
                    if np.isfinite(grad_norm) and grad_norm > 0:
                        neg_direction = -manual_grad / grad_norm
                        pos_direction = manual_grad / grad_norm
                        manual_line_search_probe(
                            posterior_objective,
                            x0,
                            neg_direction,
                            lower_bounds,
                            upper_bounds,
                            diagnostic_output_dir / "manual_line_search_negative_gradient.tsv",
                            fit_config.large_objective_penalty,
                            [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1],
                            debug_eval_fn=debug_objective_evaluation,
                        )
                        manual_line_search_probe(
                            posterior_objective,
                            x0,
                            pos_direction,
                            lower_bounds,
                            upper_bounds,
                            diagnostic_output_dir / "manual_line_search_positive_gradient.tsv",
                            fit_config.large_objective_penalty,
                            [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1],
                            debug_eval_fn=debug_objective_evaluation,
                        )
                if coord_probe is not None and (diagnostic_task_set is None or "line_search_best_coordinate" in diagnostic_task_set):
                    improving_moves = coord_probe[coord_probe["improved"]].copy()
                    if len(improving_moves) > 0:
                        best_move = improving_moves.sort_values("delta_objective", ascending=True).iloc[0]
                    else:
                        best_move = coord_probe.sort_values("delta_objective", ascending=True).iloc[0]
                    manual_best_coordinate_line_search(
                        posterior_objective,
                        x0,
                        lower_bounds,
                        upper_bounds,
                        int(best_move["coordinate_index"]),
                        str(best_move["coordinate_name"]),
                        float(np.sign(best_move["step"]) if best_move["step"] != 0 else 1.0),
                        diagnostic_output_dir / "manual_line_search_best_coordinate.tsv",
                        fit_config.large_objective_penalty,
                    )
                if diagnostic_task_set is None or "alternative_optimizers" in diagnostic_task_set:
                    diagnostic_alternative_optimizers(
                        posterior_objective,
                        x0,
                        bounds,
                        diagnostic_output_dir / "diagnostic_alternative_optimizers.tsv",
                        fit_config.max_nfev,
                    )
                if diagnostic_task_set is None or "normalized_lbfgsb" in diagnostic_task_set:
                    diagnostic_normalized_objective_lbfgsb(
                        posterior_objective,
                        x0,
                        bounds,
                        diagnostic_output_dir / "diagnostic_normalized_objective_lbfgsb.tsv",
                        fit_config.max_nfev,
                        n_observations=n_observations,
                    )
            attempt_row = {
                "start_idx": start_idx,
                "start_values": str(start_values),
                "objective": fit_config.objective,
                "observation_channels": fit_config.observation_channels,
                "model_variant": fit_config.model_variant,
                "objective_scale": float(objective_scale),
                "initial_raw_objective": float(initial_raw_objective),
                "final_raw_objective": float(final_value),
                "delta_raw_objective": float(initial_raw_objective - final_value),
                "initial_normalized_objective": float(initial_normalized_objective),
                "final_normalized_objective": float(final_normalized_objective),
                "delta_normalized_objective": float(initial_normalized_objective - final_normalized_objective),
                "initial_objective": initial_objective,
                "final_objective": final_value,
                "delta_objective": initial_objective - final_value,
                "n_iter": getattr(result, "nit", np.nan),
                "n_function_eval": getattr(result, "nfev", np.nan),
                "n_gradient_eval": getattr(result, "njev", np.nan),
                "gradient_norm": gradient_norm_from_result(result),
                "gradient_max_abs": max_abs_grad,
                "gradient_max_abs_index": max_abs_grad_idx,
                "gradient_max_abs_name": coordinate_names[max_abs_grad_idx] if max_abs_grad_idx >= 0 else "",
                "optimizer_success": bool(result.success),
                "optimizer_status": getattr(result, "status", np.nan),
                "optimizer_message": str(getattr(result, "message", "")),
                "result_equals_start": bool(moved),
                "recomputed_final_objective": float(final_value),
                "initial_simulations_attempted": initial_debug["simulations_attempted"],
                "initial_simulations_failed": initial_debug["simulations_failed"],
                "initial_failure_messages": initial_debug["failure_messages_joined"],
                "final_simulations_attempted": final_debug["simulations_attempted"],
                "final_simulations_failed": final_debug["simulations_failed"],
                "final_failure_messages": final_debug["failure_messages_joined"],
                "success": fit_is_usable,
                "message": result.message,
            }
            attempt_row.update(bound_flags)
            add_unpacked_attempt_parameters(attempt_row, "initial", x0)
            add_unpacked_attempt_parameters(attempt_row, "final", result.x)
            attempt_rows.append(attempt_row)
            if fit_is_usable and final_value < best_value:
                best_value = float(final_value)
                best_result = result
                best_unpacked = unpack_parameters(result.x)
                best_data_nll = data_nll
                best_prior_penalty = prior_penalty
        except Exception as exc:
            final_value = fit_config.large_objective_penalty
            attempt_row = {
                "start_idx": start_idx,
                "start_values": str(start_values),
                "objective": fit_config.objective,
                "observation_channels": fit_config.observation_channels,
                "model_variant": fit_config.model_variant,
                "objective_scale": float(objective_scale),
                "initial_raw_objective": float(initial_raw_objective),
                "final_raw_objective": float(final_value),
                "delta_raw_objective": float(initial_raw_objective - final_value),
                "initial_normalized_objective": float(initial_normalized_objective),
                "final_normalized_objective": np.nan,
                "delta_normalized_objective": np.nan,
                "initial_objective": initial_objective,
                "final_objective": final_value,
                "delta_objective": initial_objective - final_value,
                "n_iter": np.nan,
                "n_function_eval": np.nan,
                "n_gradient_eval": np.nan,
                "gradient_norm": np.nan,
                "gradient_max_abs": np.nan,
                "gradient_max_abs_index": np.nan,
                "gradient_max_abs_name": "",
                "optimizer_success": False,
                "optimizer_status": np.nan,
                "optimizer_message": str(exc),
                "result_equals_start": True,
                "recomputed_final_objective": float(final_value),
                "initial_simulations_attempted": initial_debug["simulations_attempted"],
                "initial_simulations_failed": initial_debug["simulations_failed"],
                "initial_failure_messages": initial_debug["failure_messages_joined"],
                "final_simulations_attempted": np.nan,
                "final_simulations_failed": np.nan,
                "final_failure_messages": str(exc),
                "success": False,
                "message": str(exc),
            }
            attempt_row.update({
                "any_at_lower_bound": np.nan,
                "any_at_upper_bound": np.nan,
                "lower_bound_coordinates": "",
                "upper_bound_coordinates": "",
            })
            add_unpacked_attempt_parameters(attempt_row, "initial", x0)
            add_unpacked_attempt_parameters(attempt_row, "final", x0)
            attempt_rows.append(attempt_row)

    attempts_df = pd.DataFrame(attempt_rows)
    if best_result is None or best_unpacked is None:
        return {
            "success": False,
            "posterior_objective": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "n_observations": n_observations,
            "objective": fit_config.objective,
            "observation_channels": fit_config.observation_channels,
            "model_variant": fit_config.model_variant,
            "theta_alive": np.nan,
            "theta_dead": np.nan,
            "optimizer_message": "All optimizer starts failed",
            "optimizer_attempts": attempts_df,
            "population_parameters": {},
            "ploidy_parameters": {},
            "dose_gate_params": {},
        }

    return {
        "success": True,
        "posterior_objective": float(best_value),
        "data_nll": float(best_data_nll),
        "prior_penalty": float(best_prior_penalty),
        "n_observations": n_observations,
        "objective": fit_config.objective,
        "observation_channels": fit_config.observation_channels,
        "model_variant": fit_config.model_variant,
        "theta_alive": best_unpacked["theta_alive"],
        "theta_dead": best_unpacked["theta_dead"],
        "optimizer_message": best_result.message,
        "optimizer_attempts": attempts_df,
        "population_parameters": best_unpacked["mu_logs"],
        "ploidy_parameters": best_unpacked["ploidy_params"],
        "dose_gate_params": best_unpacked.get("dose_gate_params", {}),
    }


def fit_one_n_tr_worker(
    args: Tuple[int, Sequence[ReplicateTrajectory], Dict[str, DfdctpSignalSurface], JointFitConfig]
) -> Tuple[int, Dict[str, Any]]:
    n_tr, trajectories, dfdctp_signal_curve_by_ploidy, fit_config = args
    start_time = time.monotonic()
    print(f"[n_tr={n_tr}] starting fit in pid {os.getpid()}", flush=True)
    try:
        result = fit_joint_partial_pooling_model(
            trajectories=trajectories,
            dfdctp_signal_curve_by_ploidy=dfdctp_signal_curve_by_ploidy,
            fit_config=fit_config,
            n_tr=n_tr,
            start_indices=fit_config.optimizer_start_indices,
        )
    except Exception as exc:
        result = {
            "success": False,
            "posterior_objective": fit_config.large_objective_penalty,
            "data_nll": np.nan,
            "prior_penalty": np.nan,
            "n_observations": np.nan,
            "objective": fit_config.objective,
            "observation_channels": fit_config.observation_channels,
            "model_variant": fit_config.model_variant,
            "theta_alive": np.nan,
            "theta_dead": np.nan,
            "optimizer_message": f"Worker exception for n_tr={n_tr}: {exc}",
            "optimizer_attempts": pd.DataFrame(),
            "population_parameters": {},
            "ploidy_parameters": {},
            "dose_gate_params": {},
        }
    elapsed_seconds = time.monotonic() - start_time
    posterior_objective = result.get("posterior_objective", np.nan)
    print(
        f"[n_tr={n_tr}] finished in {elapsed_seconds:.1f}s; "
        f"success={result.get('success')}; posterior_objective={posterior_objective}",
        flush=True,
    )
    return int(n_tr), result


def run_n_tr_model_selection(
    trajectories: Sequence[ReplicateTrajectory],
    dfdctp_signal_curve_by_ploidy: Dict[str, DfdctpSignalSurface],
    fit_config: JointFitConfig,
) -> Tuple[List[Dict[str, Any]], List[pd.DataFrame], Optional[Dict[str, Any]], Optional[int]]:
    tasks = [
        (int(n_tr), trajectories, dfdctp_signal_curve_by_ploidy, fit_config)
        for n_tr in fit_config.n_tr_values
    ]
    if fit_config.n_jobs == 1:
        print(f"Running n_tr model selection sequentially over {len(tasks)} values.")
        fit_results = [fit_one_n_tr_worker(task) for task in tasks]
    else:
        max_workers = fit_config.n_jobs
        if max_workers <= 0:
            max_workers = os.cpu_count() or 1
        print(
            f"Running n_tr model selection in parallel over {len(tasks)} values "
            f"using {max_workers} worker processes."
        )
        fit_results = []
        selection_start_time = time.monotonic()
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for task in tasks:
                n_tr = int(task[0])
                print(f"[n_tr={n_tr}] submitted to process pool", flush=True)
                futures[executor.submit(fit_one_n_tr_worker, task)] = n_tr
            for future in as_completed(futures):
                n_tr = futures[future]
                fit_results.append(future.result())
                elapsed_seconds = time.monotonic() - selection_start_time
                print(f"[n_tr={n_tr}] result collected after {elapsed_seconds:.1f}s", flush=True)

    summary_rows: List[Dict[str, Any]] = []
    attempt_frames: List[pd.DataFrame] = []
    best_fit = None
    best_n_tr = None
    best_value = np.inf

    for n_tr, result in sorted(fit_results, key=lambda item: item[0]):
        population = result.get("population_parameters", {})
        ploidy_params = result.get("ploidy_parameters", {})
        summary_row = {
            "model_preset": fit_config.model_preset,
            "n_tr": n_tr,
            "success": result["success"],
            "posterior_objective": result["posterior_objective"],
            "data_nll": result["data_nll"],
            "prior_penalty": result["prior_penalty"],
            "n_observations": result.get("n_observations", np.nan),
            "objective": result.get("objective", fit_config.objective),
            "observation_channels": result.get("observation_channels", fit_config.observation_channels),
            "model_variant": result.get("model_variant", fit_config.model_variant),
            "theta_alive": result.get("theta_alive", np.nan),
            "theta_dead": result.get("theta_dead", np.nan),
            "fit_beta_dose": fit_config.fit_beta_dose,
            "fixed_beta_dose": fit_config.fixed_beta_dose,
            "use_hill_dose_gate": fit_config.use_hill_dose_gate,
            "fit_hill_dose_gate": fit_config.fit_hill_dose_gate,
            "use_confluence_death": fit_config.use_confluence_death,
            "fit_mu_confluence_death": fit_config.fit_mu_confluence_death,
            "fixed_mu_confluence_death": fit_config.fixed_mu_confluence_death,
            "confluence_death_exponent": fit_config.confluence_death_exponent,
            "fit_dose_labels": ",".join(fit_config.fit_dose_labels) if fit_config.fit_dose_labels else "",
            "max_replicates_per_ploidy_dose": fit_config.max_replicates_per_ploidy_dose,
            "accept_incomplete_optimizer": fit_config.accept_incomplete_optimizer,
            "mu_log_r": population.get("r"),
            "mu_log_K": population.get("K"),
            "mu_log_k_tr": population.get("k_tr"),
            "mu_log_k_kill": population.get("k_kill"),
            "mu_log_k_clear": population.get("k_clear"),
            "mu_log_k_cyto": population.get("k_cyto"),
            "mu_log_beta_dose": (
                population.get("beta_dose")
                if fit_config.fit_beta_dose
                else float(fit_config.fixed_beta_dose)
            ),
            "mu_log_mu_base_death": population.get("mu_base_death"),
            "mu_mu_confluence_death": (
                population.get(CONFLUENCE_DEATH_PARAMETER_NAME)
                if fit_config.use_confluence_death and fit_config.fit_mu_confluence_death
                else float(fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0)
            ),
            "dose_gate_ec50_uM": result.get("dose_gate_params", {}).get("dose_gate_ec50_uM", np.nan),
            "dose_gate_hill": result.get("dose_gate_params", {}).get("dose_gate_hill", np.nan),
            "optimizer_message": result["optimizer_message"],
        }
        for ploidy in ["2N", "4N"]:
            params = ploidy_params.get(ploidy, {})
            for name in ALL_DYNAMIC_PARAMETER_NAMES:
                if name == DOSE_SCALING_PARAMETER_NAME and not fit_config.fit_beta_dose:
                    summary_row[f"{ploidy}_{name}"] = float(fit_config.fixed_beta_dose)
                elif name == CONFLUENCE_DEATH_PARAMETER_NAME and not (
                    fit_config.use_confluence_death and fit_config.fit_mu_confluence_death
                ):
                    summary_row[f"{ploidy}_{name}"] = float(
                        fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0
                    )
                else:
                    summary_row[f"{ploidy}_{name}"] = params.get(name, np.nan)
        summary_rows.append(summary_row)
        attempts = result.get("optimizer_attempts", pd.DataFrame()).copy()
        if not attempts.empty:
            attempts["n_tr"] = n_tr
            attempt_frames.append(attempts)
        if result["success"] and result["posterior_objective"] < best_value:
            best_value = float(result["posterior_objective"])
            best_fit = result
            best_n_tr = n_tr

    return summary_rows, attempt_frames, best_fit, best_n_tr


def plot_replicate_parameter_summary(replicate_fit_df, output_dir: Optional[Path] = None, close_fig: bool = True):
    """
    Plots fitted kinetic parameters by ploidy, with one point per row replicate.
    """
    params_to_plot = ['k_tr', 'k_kill', 'k_clear']
    param_labels = {
        'k_tr': r'$k_{tr}$',
        'k_kill': r'$k_{kill}$',
        'k_clear': r'$k_{clear}$'
    }
    
    x_map = {'2N': 0, '4N': 1}
    
    for param in params_to_plot:
        fig, ax = plt.subplots(figsize=(6, 5))
        
        for ploidy in ['2N', '4N']:
            subset = replicate_fit_df[replicate_fit_df['ploidy'] == ploidy].copy()
            
            if subset.empty:
                continue
            
            x_center = x_map[ploidy]
            jitter = np.linspace(-0.08, 0.08, len(subset))
            x_vals = x_center + jitter
            y_vals = subset[param].values
            
            ax.scatter(x_vals, y_vals, s=80, alpha=0.8, label=ploidy)
            
            # Optional: label each point by row replicate
            if 'plate_row' in subset.columns:
                for x, y, row_label in zip(x_vals, y_vals, subset['plate_row'].values):
                    ax.text(x + 0.02, y, str(row_label), fontsize=9, va='center')
            
            # Plot mean line for each ploidy
            y_mean = np.nanmean(y_vals)
            ax.hlines(y_mean, x_center - 0.18, x_center + 0.18, colors='black', linewidth=2)
        
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['2N', '4N'], fontsize=12)
        ax.set_ylabel(param_labels[param], fontsize=13)
        ax.set_title(f"Per-Replicate Estimates of {param_labels[param]}", fontsize=14, fontweight='bold')
        ax.grid(True, axis='y', linestyle='--', alpha=0.4)
        plt.tight_layout()

        if output_dir is not None:
            output_path = output_dir / f"replicate_parameter_summary_{param}.png"
            save_and_maybe_close(fig, output_path, close=close_fig, dpi=200, bbox_inches='tight')
        elif close_fig:
            plt.close(fig)


def save_effective_dose_scaling_summary(
    dfdctp_signal_curve_by_ploidy: Dict[str, DfdctpSignalSurface],
    ploidy_parameters: Dict[str, Dict[str, float]],
    fit_config: JointFitConfig,
    dose_gate_params: Optional[Dict[str, Any]],
    dose_uM_values: Sequence[float],
    output_path: Path,
) -> None:
    rows: List[Dict[str, Any]] = []
    unique_doses = sorted({float(d) for d in dose_uM_values if np.isfinite(d)})
    dose_gate_params = dose_gate_params or {}
    for ploidy, params in ploidy_parameters.items():
        curve = dfdctp_signal_curve_by_ploidy[ploidy]
        reference_dose = get_dose_scaling_reference_uM(curve)
        beta_dose = float(params.get("beta_dose", fit_config.fixed_beta_dose))
        use_hill = bool(fit_config.use_hill_dose_gate)
        gate_ec50 = (
            float(dose_gate_params.get("dose_gate_ec50_uM", fit_config.fixed_dose_gate_ec50_uM))
            if use_hill else np.nan
        )
        gate_hill = (
            float(dose_gate_params.get("dose_gate_hill", fit_config.fixed_dose_gate_hill))
            if use_hill else np.nan
        )
        ref_gate = (
            hill_dose_gate(reference_dose, gate_ec50, gate_hill)
            if use_hill else np.nan
        )
        for dose_uM in unique_doses:
            if dose_uM <= 0:
                power_correction_factor = 0.0
                hill_gate_raw = 0.0 if use_hill else np.nan
                hill_gate_normalized = 0.0 if use_hill else 1.0
            else:
                power_correction_factor = (dose_uM / reference_dose) ** (beta_dose - 1.0)
                if use_hill:
                    hill_gate_raw = hill_dose_gate(dose_uM, gate_ec50, gate_hill)
                    hill_gate_normalized = normalized_hill_dose_gate(
                        dose_uM=dose_uM,
                        reference_dose_uM=reference_dose,
                        ec50_uM=gate_ec50,
                        hill_coefficient=gate_hill,
                    )
                else:
                    hill_gate_raw = np.nan
                    hill_gate_normalized = 1.0
            rows.append({
                "ploidy": ploidy,
                "beta_dose": beta_dose,
                "reference_dose_uM": reference_dose,
                "dose_uM": dose_uM,
                "power_correction_factor": power_correction_factor,
                "use_hill_dose_gate": use_hill,
                "dose_gate_ec50_uM": gate_ec50,
                "dose_gate_hill": gate_hill,
                "hill_gate_raw": hill_gate_raw,
                "hill_gate_at_reference": ref_gate,
                "hill_gate_normalized": hill_gate_normalized,
                "combined_correction_factor": power_correction_factor * hill_gate_normalized,
            })
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def save_baseline_death_summary(
    ploidy_parameters: Dict[str, Dict[str, float]],
    trajectories: Sequence[ReplicateTrajectory],
    output_path: Path,
) -> None:
    rows: List[Dict[str, Any]] = []
    for ploidy, params in ploidy_parameters.items():
        mu_base_death = float(params.get("mu_base_death", np.nan))
        baseline_trajectories = [
            traj for traj in trajectories
            if traj.ploidy == ploidy and np.isclose(float(traj.dose_uM), 0.0)
        ]
        a0_values = [float(traj.N0) for traj in baseline_trajectories if np.isfinite(traj.N0)]
        dead_increase_per_day_values = []
        for traj in baseline_trajectories:
            valid = np.isfinite(traj.t) & np.isfinite(traj.dead)
            if np.count_nonzero(valid) < 2:
                continue
            t_valid = traj.t[valid]
            dead_valid = traj.dead[valid]
            baseline_t = float(t_valid[0])
            candidate_idx = np.flatnonzero(t_valid > baseline_t)
            if len(candidate_idx) == 0:
                continue
            idx = int(candidate_idx[0])
            dt = float(t_valid[idx] - baseline_t)
            if dt > 0:
                dead_increase_per_day_values.append(float((dead_valid[idx] - dead_valid[0]) / dt))
        rows.append(
            {
                "ploidy": ploidy,
                "mu_base_death": mu_base_death,
                "baseline_death_half_life_days": np.log(2.0) / mu_base_death if mu_base_death > 0 else np.inf,
                "A0_mean": float(np.nanmean(a0_values)) if a0_values else np.nan,
                "predicted_0nM_dead_generation_rate_at_t0": (
                    mu_base_death * float(np.nanmean(a0_values))
                    if mu_base_death >= 0 and a0_values else np.nan
                ),
                "mean_observed_0nM_dead_increase_per_day_first_interval": (
                    float(np.nanmean(dead_increase_per_day_values))
                    if dead_increase_per_day_values else np.nan
                ),
            }
        )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def save_confluence_death_summary(
    ploidy_parameters: Dict[str, Dict[str, float]],
    output_path: Path,
    fit_config: JointFitConfig,
    density_fractions: Sequence[float] = (0.25, 0.5, 0.75, 1.0, 1.25),
) -> None:
    rows: List[Dict[str, Any]] = []
    for ploidy, params in ploidy_parameters.items():
        K = float(params.get("K", np.nan))
        mu_base_death = float(params.get(BASELINE_DEATH_PARAMETER_NAME, np.nan))
        mu_confluence_death = float(
            params.get(
                CONFLUENCE_DEATH_PARAMETER_NAME,
                fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
            )
        )
        for density_fraction in density_fractions:
            alive_count = float(K * density_fraction) if np.isfinite(K) else np.nan
            confluence_hazard = (
                confluence_death_hazard(
                    alive_count=alive_count,
                    K=K,
                    mu_confluence_death=mu_confluence_death,
                    confluence_death_exponent=fit_config.confluence_death_exponent,
                )
                if fit_config.use_confluence_death and np.isfinite(K) and np.isfinite(alive_count)
                else 0.0
            )
            total_baseline_hazard = (
                total_drug_independent_death_hazard(
                    alive_count=alive_count,
                    K=K,
                    mu_base_death=mu_base_death,
                    use_confluence_death=fit_config.use_confluence_death,
                    mu_confluence_death=mu_confluence_death,
                    confluence_death_exponent=fit_config.confluence_death_exponent,
                )
                if np.isfinite(K) and np.isfinite(alive_count) and np.isfinite(mu_base_death)
                else np.nan
            )
            rows.append(
                {
                    "ploidy": ploidy,
                    "K": K,
                    "density_fraction": float(density_fraction),
                    "alive_count_at_fraction": alive_count,
                    "mu_base_death": mu_base_death,
                    "use_confluence_death": bool(fit_config.use_confluence_death),
                    "mu_confluence_death": mu_confluence_death,
                    "confluence_death_exponent": float(fit_config.confluence_death_exponent),
                    "confluence_death_hazard": confluence_hazard,
                    "total_drug_independent_death_hazard": total_baseline_hazard,
                }
            )
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False)


def main(
    paths: Optional[ExperimentPaths] = None,
    pk_config: Optional[PKConfig] = None,
    fit_config: Optional[JointFitConfig] = None,
) -> Dict[str, Any]:
    paths = paths or default_experiment_paths()
    pk_config = pk_config or PKConfig()
    fit_config = resolve_joint_fit_config(fit_config)
    validate_paths(paths, required=("counts_agg", "platemap", "pkpd_constants"))
    paths.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Saving figures to: {paths.output_dir}")
    hill_status = (
        "disabled"
        if not fit_config.use_hill_dose_gate
        else ("fitted" if fit_config.fit_hill_dose_gate else "fixed")
    )
    beta_status = "fitted" if fit_config.fit_beta_dose else f"fixed ({fit_config.fixed_beta_dose:.4f})"
    confluence_status = (
        "disabled"
        if not fit_config.use_confluence_death
        else (
            "fitted"
            if fit_config.fit_mu_confluence_death
            else f"fixed ({fit_config.fixed_mu_confluence_death:.4f})"
        )
    )
    print("Model configuration:")
    print(f"  preset: {fit_config.model_preset}")
    print(f"  model_variant: {fit_config.model_variant}")
    print(f"  beta: {beta_status}")
    print(f"  Hill gate: {hill_status}")
    print(f"  confluence death: {confluence_status}")

    pk_sheets = import_and_clean_pkpd(paths=paths, pk_config=pk_config)
    print_pk_workbook_summary(pk_sheets)
    dfdctp_signal_curve_by_ploidy = build_preferred_dfdctp_signal_surfaces(
        pk_sheets,
        ploidy_keys=("2N", "4N"),
        fallback_half_life_days=pk_config.fallback_half_life_days,
    )
    save_pk_tail_diagnostics(dfdctp_signal_curve_by_ploidy, paths.output_dir)

    df = assemble_modeling_dataset(paths=paths, fit_config=fit_config)
    ploidy_options = [p for p in df['ploidy'].unique() if pd.notna(p)]
    gem_doses = sorted([d for d in df['gem'].unique() if pd.notna(d)], key=lambda x: float(x.split()[0]))
    if fit_config.fit_dose_labels is not None:
        requested_doses = set(fit_config.fit_dose_labels)
        missing_doses = sorted(requested_doses.difference(gem_doses), key=lambda x: float(x.split()[0]))
        if missing_doses:
            raise ValueError(f"Requested fit dose labels are absent from the modeling dataset: {missing_doses}")
        gem_doses = [dose for dose in gem_doses if dose in requested_doses]
    live_dead_dose_uM_values = [float(dose.split()[0]) / 1000.0 for dose in gem_doses if dose != "0 nM"]
    dfdctp_y_limits = get_dfdctp_signal_curve_y_limits(dfdctp_signal_curve_by_ploidy)

    for ploidy_key in ["2N", "4N"]:
        print_dfdctp_signal_curve_summary(
            ploidy_key,
            dfdctp_signal_curve_by_ploidy[ploidy_key],
            live_dead_dose_uM_values=live_dead_dose_uM_values,
        )
        plot_dfdctp_signal_curve(
            ploidy_key,
            dfdctp_signal_curve_by_ploidy[ploidy_key],
            output_dir=paths.output_dir,
            y_limits=dfdctp_y_limits,
        )
        plot_dfdctp_amplitude_scaling(ploidy_key, dfdctp_signal_curve_by_ploidy[ploidy_key], output_dir=paths.output_dir)

    trajectories = build_joint_fit_trajectories(
        df=df,
        ploidies=ploidy_options,
        gem_doses=gem_doses,
        fit_t_max=fit_config.fit_t_max,
        count_transitional_as_alive=fit_config.count_transitional_as_alive,
        max_replicates_per_ploidy_dose=fit_config.max_replicates_per_ploidy_dose,
    )

    summary_rows, attempt_frames, best_fit, best_n_tr = run_n_tr_model_selection(
        trajectories=trajectories,
        dfdctp_signal_curve_by_ploidy=dfdctp_signal_curve_by_ploidy,
        fit_config=fit_config,
    )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(paths.output_dir / "joint_fit_summary.tsv", sep="\t", index=False)
    attempts_df = pd.concat(attempt_frames, ignore_index=True) if attempt_frames else pd.DataFrame()
    attempts_df.to_csv(paths.output_dir / "optimizer_attempts.tsv", sep="\t", index=False)

    if best_fit is None or best_n_tr is None:
        print("No successful joint partial-pooling fit found.")
        return {"success": False, "joint_fit_summary": summary_df, "optimizer_attempts": attempts_df}

    print(f"Best joint partial-pooling fit used n_tr={best_n_tr}")
    print(f"Posterior objective = data_nll + prior_penalty = {best_fit['posterior_objective']:.4f}")
    print(f"Data NLL: {best_fit['data_nll']:.4f}")
    print(f"Prior penalty: {best_fit['prior_penalty']:.4f}")
    print(f"Objective: {best_fit['objective']}")
    print(f"Observation channels: {best_fit['observation_channels']}")
    print(f"Model variant: {best_fit['model_variant']}")
    print(
        "Observed compartments: "
        + (
            "alive = Alive + Transitional; dead = Dead"
            if fit_config.count_transitional_as_alive
            else "alive = Alive; dead = Dead + Transitional"
        )
    )
    if fit_config.use_hill_dose_gate:
        gate_params = best_fit.get("dose_gate_params", {})
        print(
            "Hill dose gate: enabled "
            f"(fit={fit_config.fit_hill_dose_gate}, "
            f"EC50={gate_params.get('dose_gate_ec50_uM', np.nan):.6f} uM, "
            f"h={gate_params.get('dose_gate_hill', np.nan):.4f})"
        )
    else:
        print("Hill dose gate: disabled")
    if fit_config.use_confluence_death:
        print(
            "Confluence death: enabled "
            f"(fit={fit_config.fit_mu_confluence_death}, "
            f"fixed_mu={fit_config.fixed_mu_confluence_death:.4f}, "
            f"exponent={fit_config.confluence_death_exponent:.2f})"
        )
    else:
        print("Confluence death: disabled")
    save_effective_dose_scaling_summary(
        dfdctp_signal_curve_by_ploidy=dfdctp_signal_curve_by_ploidy,
        ploidy_parameters=best_fit["ploidy_parameters"],
        fit_config=fit_config,
        dose_gate_params=best_fit.get("dose_gate_params", {}),
        dose_uM_values=live_dead_dose_uM_values,
        output_path=paths.output_dir / "effective_dose_scaling_summary.tsv",
    )
    save_baseline_death_summary(
        ploidy_parameters=best_fit["ploidy_parameters"],
        trajectories=trajectories,
        output_path=paths.output_dir / "baseline_death_summary.tsv",
    )
    save_confluence_death_summary(
        ploidy_parameters=best_fit["ploidy_parameters"],
        output_path=paths.output_dir / "confluence_death_summary.tsv",
        fit_config=fit_config,
    )

    for ploidy in ["2N", "4N"]:
        params = best_fit["ploidy_parameters"][ploidy]
        param_line = (
            f"{ploidy}: r={params['r']:.4f}, K={params['K']:.4f}, "
            f"k_tr={params['k_tr']:.4f}, k_kill={params['k_kill']:.4f}, "
            f"k_clear={params['k_clear']:.4f}, "
            f"beta_dose={params.get('beta_dose', fit_config.fixed_beta_dose):.4f}, "
            f"mu_base_death={params['mu_base_death']:.4f}, "
            f"mu_confluence_death={params.get(CONFLUENCE_DEATH_PARAMETER_NAME, fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0):.4f}"
        )
        if "k_cyto" in params:
            param_line += f", k_cyto={params['k_cyto']:.4f}"
        if fit_config.use_hill_dose_gate:
            gate_params = best_fit.get("dose_gate_params", {})
            param_line += (
                f", dose_gate_ec50_uM={gate_params.get('dose_gate_ec50_uM', np.nan):.6f}, "
                f"dose_gate_hill={gate_params.get('dose_gate_hill', np.nan):.4f}"
            )
        print(param_line)

        dose_data_list = []
        for gem_dose in gem_doses:
            try:
                aligned = get_aligned_live_dead_data(
                    df,
                    gem_dose=gem_dose,
                    ploidy=ploidy,
                    t_max=fit_config.fit_t_max,
                    count_transitional_as_alive=fit_config.count_transitional_as_alive,
                )
                dose_data_list.append(
                    {
                        "dose_label": gem_dose,
                        "dose_muM": float(gem_dose.split()[0]) / 1000.0,
                        "t": aligned["t"],
                        "y_alive": aligned["y_alive"],
                        "y_dead": aligned["y_dead"],
                        "N0": float(np.nanmean(aligned["y_alive"][0, :])),
                        "D0": float(np.nanmean(aligned["y_dead"][0, :])),
                    }
                )
            except ValueError:
                continue
        if dose_data_list:
            plot_global_fit_subplots_joint(
                dose_data_list=dose_data_list,
                ploidy_label=ploidy,
                r=params["r"],
                K=params["K"],
                n_tr=best_n_tr,
                dfdctp_signal_curve=dfdctp_signal_curve_by_ploidy[ploidy],
                k_tr=params["k_tr"],
                k_kill=params["k_kill"],
                k_clear=params["k_clear"],
                k_cyto=params.get("k_cyto"),
                beta_dose=params.get("beta_dose", fit_config.fixed_beta_dose),
                mu_base_death=params["mu_base_death"],
                mu_confluence_death=params.get(
                    CONFLUENCE_DEATH_PARAMETER_NAME,
                    fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
                ),
                use_confluence_death=fit_config.use_confluence_death,
                confluence_death_exponent=fit_config.confluence_death_exponent,
                use_hill_dose_gate=fit_config.use_hill_dose_gate,
                dose_gate_ec50_uM=best_fit.get("dose_gate_params", {}).get("dose_gate_ec50_uM", np.nan),
                dose_gate_hill=best_fit.get("dose_gate_params", {}).get("dose_gate_hill", np.nan),
                model_variant=fit_config.model_variant,
                fit_config=fit_config,
                fit_summary={
                    "objective": best_fit["objective"],
                    "observation_channels": best_fit["observation_channels"],
                    "model_variant": best_fit["model_variant"],
                    "theta_alive": best_fit["theta_alive"],
                    "theta_dead": best_fit["theta_dead"],
                    "beta_dose": params.get("beta_dose", fit_config.fixed_beta_dose),
                    "mu_base_death": params["mu_base_death"],
                    "mu_confluence_death": params.get(
                        CONFLUENCE_DEATH_PARAMETER_NAME,
                        fit_config.fixed_mu_confluence_death if fit_config.use_confluence_death else 0.0,
                    ),
                    "use_confluence_death": fit_config.use_confluence_death,
                    "confluence_death_exponent": fit_config.confluence_death_exponent,
                    "use_hill_dose_gate": fit_config.use_hill_dose_gate,
                    "dose_gate_ec50_uM": best_fit.get("dose_gate_params", {}).get("dose_gate_ec50_uM", np.nan),
                    "dose_gate_hill": best_fit.get("dose_gate_params", {}).get("dose_gate_hill", np.nan),
                    "nll": best_fit["data_nll"],
                    "aic": None,
                    "bic": None,
                    "rmse": np.nan,
                },
                output_dir=paths.output_dir,
            )

    return {
        "success": True,
        "best_fit": best_fit,
        "best_n_tr": best_n_tr,
        "joint_fit_summary": summary_df,
        "optimizer_attempts": attempts_df,
    }


if __name__ == "__main__":
    main()
        
