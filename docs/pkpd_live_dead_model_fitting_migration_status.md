# PKPD Live/Dead Fitting Migration Status

Updated: 2026-06-29

## Completed Scope

The migration plan in `docs/pkpd_live_dead_model_fitting_migration_plan.md` has been executed through the compatibility milestone:

- Saved `joint_fit_summary.tsv` folders reproduce downstream plots through `Code/plot_invitro_fit_outputs.py` and `Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py`.
- `Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py` can validate inputs with `--check-inputs`.
- `run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1` runs the real fitting path serially, writes `joint_fit_summary.tsv` and `optimizer_attempts.tsv`, and does not require multiprocessing.
- A new smoke-fit summary can be passed back into the plotter for both 50 nM and 25 nM single-dose comparison plots.
- Regression tests cover saved-summary reproduction and smoke-fit-to-plotter compatibility.

## Step 6 Decision

No broad function-splitting refactor was performed after the end-to-end tests passed. The imported `src/invitro_fitting.py` snapshot is still intentionally retained as the maintained compatibility module for now.

Reasoning:

- The current completion criteria are satisfied.
- The wrapper and tests now provide a stable safety net for future extraction.
- Splitting the 5,000-line fitting snapshot into path, data-loading, PKPD, model, optimizer, and plotting modules would be a separate numerical-behavior-preservation refactor.
- Performing that split in the same migration would add risk without changing the validated user-facing workflow.

Future refactors should follow the extraction order in the plan and rerun both:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_saved_summary.py
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_fit.py
```

after each extraction.

## Validation Run

The following validations passed in the current environment:

```bash
python3 -m py_compile Code/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_saved_summary.py Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_fit.py
python3 Code/plot_invitro_fit_outputs.py --skip-cohort --skip-dose-comparison
python3 Code/plot_invitro_fit_outputs.py Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159 --skip-cohort --skip-dose-comparison
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --check-inputs
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/20260629T070102
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/20260629T070102 --comparison-dose '25 nM'
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_saved_summary.py
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_fit.py
```

Notes:

- The first unreduced smoke fit was interrupted because it used the full trajectory set and was too slow for routine validation.
- Smoke mode now uses both ploidies, control plus 25/50 nM doses, one replicate per ploidy-dose, one `n_tr` value, one optimizer start, serial execution, and a short fitting window.
- PyArrow emitted sandbox-related `sysctlbyname` warnings during smoke testing, and Matplotlib emitted tight-layout warnings. These did not fail validation.
- Timestamped validation output folders under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/` were left untracked.
