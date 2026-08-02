from __future__ import annotations

from pathlib import Path

from figure_output_contract import write_tsv


def write_density_localization_fixture(run_root: Path) -> list[Path]:
    """Write a compact synthetic fixture satisfying the reviewed table contract."""
    critical = 2.6737373166169203
    global_p = 232 / 4900
    grid_rows = []
    for index in range(501):
        pseudotime = index / 500
        pointwise = 0.296 <= pseudotime <= 0.486
        simultaneous = 0.414 <= pseudotime <= 0.426
        if pseudotime == 0.420:
            difference = 0.271136208738626
            studentized = 2.68138128387289
        elif pseudotime == 0.452:
            difference = 0.283948126904316
            studentized = 2.49554727755694
        elif simultaneous:
            difference = 0.25
            studentized = 2.68
        elif pointwise:
            difference = 0.20
            studentized = 2.0
        else:
            difference = -0.10
            studentized = -1.0
        permutation_sd = difference / studentized
        pointwise_p = 0.01 if pointwise else 0.50
        adjusted_p = (
            global_p
            if pseudotime == 0.420
            else (0.049 if simultaneous else (0.20 if pointwise else 0.80))
        )
        grid_rows.append(
            {
                "pseudotime": pseudotime,
                "vehicle_mean_density": 1.0,
                "treated_mean_density": 1.0 + difference,
                "treated_minus_vehicle_density": difference,
                "permutation_sd": permutation_sd,
                "observed_studentized": studentized,
                "pointwise_p_two_sided": pointwise_p,
                "max_t_adjusted_p_two_sided": adjusted_p,
                "simultaneous_critical": critical,
                "simultaneous_lower_envelope": -critical * permutation_sd,
                "simultaneous_upper_envelope": critical * permutation_sd,
                "pointwise_positive_supported": str(pointwise).upper(),
                "simultaneous_positive_supported": str(simultaneous).upper(),
                "frozen_state_interval": str(
                    0.30 <= pseudotime <= 0.49
                ).upper(),
            }
        )
    grid_path = run_root / "tables/panel_7B_density_localization_grid.tsv"
    write_tsv(grid_path, grid_rows, list(grid_rows[0]))

    interval_rows = [
        {
            "support_type": "positive_pointwise_two_sided",
            "start": 0.296,
            "end": 0.486,
            "width": 0.190,
            "alpha": 0.05,
        },
        {
            "support_type": "positive_simultaneous_max_abs_t",
            "start": 0.414,
            "end": 0.426,
            "width": 0.012,
            "alpha": 0.05,
        },
    ]
    interval_path = (
        run_root / "tables/panel_7B_density_localization_intervals.tsv"
    )
    write_tsv(interval_path, interval_rows, list(interval_rows[0]))

    test_row = {
        "analysis_id": "equal_mouse_kde_exact_origin_stratified_max_t_v1",
        "cell_universe": "reviewed_qc_retained_cellcycle_2881",
        "n_cells": 2881,
        "n_samples": 16,
        "n_vehicle_samples": 8,
        "n_treated_samples": 8,
        "sample_weighting": "equal_mouse",
        "density_estimator": "stats::density Gaussian kernel",
        "bandwidth_method": "pooled label-invariant stats::bw.nrd0",
        "bandwidth": 0.0506470455660707,
        "grid_start": 0,
        "grid_end": 1,
        "grid_points": 501,
        "contrast": "treated_minus_vehicle",
        "permutation_strata": "initial_ploidy",
        "n_permutations": 4900,
        "pointwise_test": "two-sided exact permutation at each grid point",
        "pointwise_alpha": 0.05,
        "pointwise_start": 0.296,
        "pointwise_end": 0.486,
        "simultaneous_test": (
            "studentized max-absolute-T exact permutation null envelope"
        ),
        "simultaneous_alpha": 0.05,
        "simultaneous_critical": critical,
        "simultaneous_start": 0.414,
        "simultaneous_end": 0.426,
        "raw_excess_peak_pseudotime": 0.452,
        "raw_excess_peak_density_difference": 0.283948126904316,
        "max_abs_t_pseudotime": 0.420,
        "max_abs_t_observed_statistic": 2.68138128387289,
        "global_max_abs_t": 2.68138128387289,
        "global_max_t_p_two_sided": global_p,
    }
    test_path = run_root / "tables/panel_7B_density_localization_test.tsv"
    write_tsv(test_path, [test_row], list(test_row))
    return [grid_path, interval_path, test_path]
