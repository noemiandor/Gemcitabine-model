# PKPD Live/Dead Fitting Migration Plan

## Objective

Migrate the in vitro gemcitabine PKPD/live-dead fitting workflow from `/Users/4470246/Repositories/miningcloneid` into this repository so that this repo can both:

1. Reproduce manuscript figures from saved `joint_fit_summary.tsv` files.
2. Regenerate new timestamped fit-output folders containing `joint_fit_summary.tsv`, `optimizer_attempts.tsv`, and all downstream plots.

The current repository already satisfies the first goal through `Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py`, the compatibility entrypoint `Code/plot_invitro_fit_outputs.py`, copied raw inputs, copied saved summaries, and regenerated figure outputs. The next migration phase should make the fitting machinery itself maintainable and testable in this repo.

## Scope

Preserve these capabilities first:

- Loading and cleaning Incucyte live/dead count data.
- Loading plate-map and count data.
- Loading and cleaning the gemcitabine PK/PD workbook.
- Building dFdCTP signal surfaces.
- Assembling aligned live/dead trajectories.
- Fitting joint 2N/4N live/dead models.
- Writing `joint_fit_summary.tsv`.
- Writing `optimizer_attempts.tsv`.
- Regenerating all downstream plots from the saved summary.

Do not stub missing functions. If a function is missing, import the actual implementation from `/Users/4470246/Repositories/miningcloneid/code/invitro_fitting.py` or copy a narrow, intact dependency block that preserves behavior.

## Current State

Committed migration:

- `Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py`
- `Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py`
- `Code/plot_invitro_fit_outputs.py`
- `Data/in-vitro/pkpd_live_dead_model/raw/`
- `Data/GemDelayKillTerm/processed/`
- `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/{alsoGoodFit_20260514T093906,bestFitSoFar_20260513T164159}/`

Current reproduction commands:

```bash
python3 -m py_compile Code/plot_invitro_fit_outputs.py
python3 Code/plot_invitro_fit_outputs.py
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159 \
  --comparison-dose '25 nM'
```

The compatibility entrypoint also accepts legacy-style positional paths containing `code/invitro_fitting_outputs/...`, but canonical target-repo documentation and tests should use `Code/` and `Data/` paths so the workflow remains portable to case-sensitive filesystems.

The copied `src/invitro_fitting.py` already contains the core fitting functions, but it is still a compatibility snapshot. The migration is not complete until there is a documented fitting entrypoint and validation that a new fit folder can be generated and fed back into the plotter. Do not broadly refactor or split this snapshot until the wrapper can run `--check-inputs`, run a smoke fit, produce `joint_fit_summary.tsv`, and feed that summary back into `plot_invitro_fit_outputs.py`.

## Visualization Contract To Preserve

The saved-summary plotter currently generates these outputs from a valid fit-output folder:

- `ploidy_parameter_comparison.tsv`
- `ploidy_parameter_log2_fold_change.png`
- `ploidy_parameter_paired_values.png`
- `pk_tail_fit_diagnostics.tsv`
- `dfdctp_signal_curve_2n.png`
- `dfdctp_signal_curve_4n.png`
- `dfdctp_signal_curve_combined_ploidy.png`
- `dfdctp_amplitude_scaling_2n.png`
- `dfdctp_amplitude_scaling_4n.png`
- `dose_response_ploidy_comparison.tsv`
- `dose_response_ploidy_comparison.png`
- `cohort_joint_fit_2n.png`
- `cohort_joint_fit_4n.png`
- `dose_50_nm_2n_vs_4n.png`
- `dose_25_nm_2n_vs_4n.png` when run with `--comparison-dose '25 nM'`

Current plot styling must remain stable unless deliberately changed:

- Alive observations are gray/black.
- Alive model lines are black solid lines.
- Dead observations are light red.
- Dead model lines are red dashed lines.
- Cohort fit plots use a 5 row x 2 column layout.
- Cohort and single-dose comparison plots use integer day x-axis ticks.
- Single-ploidy and combined dFdCTP signal-curve plots use a log-scaled x-axis with day 0 shown at a small plotting floor.

Canonical validation commands for this contract:

```bash
python3 -m py_compile Code/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py
python3 Code/plot_invitro_fit_outputs.py
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159 \
  --comparison-dose '25 nM'
```

## Target Layout

Keep the module boundaries explicit:

```text
Code/in-vitro/pkpd_live_dead_model/
  README.md
  run_invitro_fit.py
  plot_invitro_fit_outputs.py
  src/
    __init__.py
    paths.py
    data_loading.py
    pkpd.py
    model.py
    fitting.py
    plotting.py
    invitro_fitting.py        # Temporary compatibility module while refactoring
  tests/
    smoke_test_saved_summary.py
    smoke_test_fit.py

Data/in-vitro/pkpd_live_dead_model/
  raw/
  processed/
  invitro_fitting_outputs/
```

The first implementation should not aggressively split `invitro_fitting.py`. Prefer preserving the working implementation intact, adding a clear fitting entrypoint, and only then extracting reusable modules after tests are passing.

## Implementation Steps

### Step 1: Freeze Saved-Summary Reproduction Baseline

Record the current saved-summary reproduction behavior before touching the fitting path.

Actions:

- Run the existing plotter on both imported summary folders.
- Record expected output names, table schemas, nonzero file sizes, and image dimensions where stable in a small manifest under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/`.
- Do not require exact PNG checksums for Matplotlib outputs; font and renderer differences make them brittle. Use perceptual/pixel-tolerance checks only if a later review shows they are necessary.
- Add a smoke test that verifies:
  - `joint_fit_summary.tsv` exists.
  - The best row can be selected by minimum `posterior_objective`.
  - Required downstream plots are produced.

Validation:

```bash
python3 -m py_compile Code/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py
python3 Code/plot_invitro_fit_outputs.py
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159 \
  --comparison-dose '25 nM'
```

Acceptance:

- Existing saved summaries still regenerate plots.
- No fitting-code change is allowed to break this baseline.
- Every required plot listed in the visualization contract exists and has nonzero size.
- TSV outputs preserve required schemas.

### Step 2: Add a Fitting Entrypoint Without Refactoring Internals

Create `Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py` as a thin CLI wrapper around the current intact `src/invitro_fitting.py`.

Actions:

- Import the real fitting functions from `src/invitro_fitting.py`.
- Preserve the existing fitting defaults from the source workflow.
- Support a timestamped output folder under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/<timestamp>/`.
- Add options for:
  - `--output-dir`
  - `--max-days`
  - `--fit-t-max`
  - `--n-starts`
  - `--max-parallel`
  - `--n-jobs`
  - `--n-tr-values`
  - `--max-nfev`
  - `--model-preset`
  - `--smoke-test`
  - `--overwrite`
- Map every CLI option to a real existing internal field or optimizer setting, such as `JointFitConfig`, `n_jobs`, `n_tr_values`, `max_nfev`, and optimizer attempt-grid parameters. Do not add flags that are parsed but silently unused.
- Prefer CLI names that match actual internals where practical. If compatibility names such as `--max-parallel` are retained, document that they map to the internal serial/parallel execution setting.
- Define `--smoke-test` explicitly. It must:
  - force serial execution (`--max-parallel 1` or internal `n_jobs = 1`);
  - reduce `n_tr_values` to one or two values;
  - reduce optimizer starts/attempts to the minimum needed to exercise summary writing;
  - reduce `max_nfev`;
  - optionally subset doses or ploidies only if the resulting `joint_fit_summary.tsv` still exercises the same plotter contract.
- Ensure the wrapper writes `joint_fit_summary.tsv` and `optimizer_attempts.tsv`.

Validation:

```bash
python3 -m py_compile Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1
```

Acceptance:

- A new output directory is created.
- `joint_fit_summary.tsv` and `optimizer_attempts.tsv` are present.
- The output folder can be passed to `plot_invitro_fit_outputs.py`.
- Validation does not require multiprocessing; serial execution is the default for smoke tests.

### Step 3: Validate Data-Loading and PKPD Subsystems Independently

Before broad fitting runs, add targeted checks for the raw inputs.

Actions:

- Add a data-check mode or smoke test that calls:
  - `default_experiment_paths`
  - `import_and_clean_pkpd`
  - `build_preferred_dfdctp_signal_surfaces`
  - `assemble_modeling_dataset`
  - `get_aligned_live_dead_data`
- Verify that expected ploidy labels, dose labels, replicate dimensions, and time ranges are present.
- Write a short `input_manifest.tsv` with paths and row counts.

Validation:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --check-inputs
```

Acceptance:

- The command fails loudly on missing raw files.
- Counts, plate map, and PKPD workbook paths resolve only within this repo.
- No hidden dependency on `/Users/4470246/Repositories/miningcloneid` remains.

### Step 4: Feed New Fit Outputs Back Into the Plotter

Confirm the expected end-to-end workflow.

Actions:

- Run a minimal fit producing a timestamped output folder.
- Pass that folder to `plot_invitro_fit_outputs.py`.
- Confirm the plotter chooses the best row by minimum `posterior_objective`.
- Confirm the same downstream plot filenames are generated.

Validation:

```bash
OUTDIR=$(ls -td Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/* | head -1)
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py "$OUTDIR"
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py "$OUTDIR" --comparison-dose '25 nM'
```

Acceptance:

- New fitting output is plot-compatible.
- Required generated files include:
  - `ploidy_parameter_comparison.tsv`
  - `ploidy_parameter_log2_fold_change.png`
  - `ploidy_parameter_paired_values.png`
  - `pk_tail_fit_diagnostics.tsv`
  - `dfdctp_signal_curve_2n.png`
  - `dfdctp_signal_curve_4n.png`
  - `dfdctp_signal_curve_combined_ploidy.png`
  - `dfdctp_amplitude_scaling_2n.png`
  - `dfdctp_amplitude_scaling_4n.png`
  - `dose_response_ploidy_comparison.tsv`
  - `dose_response_ploidy_comparison.png`
  - `cohort_joint_fit_2n.png`
  - `cohort_joint_fit_4n.png`
  - `dose_50_nm_2n_vs_4n.png`
  - `dose_25_nm_2n_vs_4n.png`

### Step 5: Add Regression/Parity Checks Against Imported Summaries

The full optimizer may be stochastic, so do not require exact recovery of the archived best-fit summaries unless seeds and optimizer settings are fully controlled. Instead, test structural compatibility first.

Actions:

- Compare columns in new `joint_fit_summary.tsv` to imported summaries.
- Check that required columns for plotting exist:
  - `posterior_objective`
  - `objective`
  - `observation_channels`
  - `model_variant`
  - `n_tr`
  - `theta_alive`
  - `theta_dead`
  - 2N/4N fitted parameter columns used by the plotter.
- Compare output filenames and table schemas.
- Confirm required PNG outputs exist, have nonzero size, and have expected dimensions if stable.
- Confirm style/layout expectations from the visualization contract, especially black/gray alive, red/light-red dead, 5 row x 2 column cohort plots, integer day ticks, and log-scaled dFdCTP signal-curve x-axes.
- If deterministic smoke runs are introduced, record expected checksums for tables only, not PNGs.

Validation:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_saved_summary.py
python3 Code/in-vitro/pkpd_live_dead_model/tests/smoke_test_fit.py
```

Acceptance:

- New fit outputs satisfy the plotter contract.
- Saved-summary reproduction remains unchanged except where intentional plot styling changes were made.

### Step 6: Refactor Conservatively After End-to-End Success

Only after the wrapper and tests pass, split reusable functions out of the snapshot module.

Preferred extraction order:

1. Path constants and `ExperimentPaths`.
2. Data loading and cleaning.
3. PKPD cleaning and dFdCTP signal-surface construction.
4. Live/dead trajectory assembly.
5. Model simulation.
6. Optimizer and summary writers.
7. Plotting helpers.

Rules:

- Move functions intact first.
- Keep wrapper imports compatible during each move.
- After each extraction, run saved-summary and smoke-fit tests.
- Avoid changing numerical behavior while moving code.

## Validation Matrix

Required after each milestone:

```bash
python3 -m py_compile Code/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py Code/in-vitro/pkpd_live_dead_model/src/invitro_fitting.py
python3 Code/plot_invitro_fit_outputs.py --skip-cohort --skip-dose-comparison
python3 Code/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159 \
  --skip-cohort --skip-dose-comparison
```

Required once fitting entrypoint exists:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --check-inputs
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py <new-output-dir>
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py <new-output-dir> --comparison-dose '25 nM'
```

## Risks and Guardrails

- **Hidden path dependencies**: fail if any default path points outside this repo.
- **Numerical drift during refactor**: move functions intact before changing interfaces.
- **Optimizer runtime**: keep a cheap `--smoke-test` mode and do not make full fitting mandatory for every validation.
- **Multiprocessing fragility**: smoke tests must default to serial execution and should not require multiprocessing to pass in sandboxed or CI contexts.
- **Stochastic output**: compare schemas and required columns first; reserve exact checksums for deterministic fixtures.
- **Large generated outputs**: commit only manuscript-relevant saved summaries, compact smoke outputs, or explicitly requested regenerated figures.
- **Plotter contract breakage**: any new fitting output must remain readable by `plot_invitro_fit_outputs.py`.

## Completion Criteria

The fitting migration is complete when:

- `python3 Code/plot_invitro_fit_outputs.py` still reproduces saved-fit plots.
- `python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1` creates a new output folder without requiring multiprocessing.
- The new folder contains `joint_fit_summary.tsv` and `optimizer_attempts.tsv`.
- The new `joint_fit_summary.tsv` can be passed back into `plot_invitro_fit_outputs.py`.
- The README documents both saved-summary reproduction and fitting regeneration.
- No runtime path dependency on `/Users/4470246/Repositories/miningcloneid` remains.
