# PKPD Live/Dead Model Figure Reproduction

This module reproduces manuscript-style in vitro PKPD/live-dead model panels from saved fitting summaries. It was migrated from the working-tree version of `/Users/4470246/Repositories/miningcloneid/code/plot_invitro_fit_outputs.py`.

## Inputs

The plotter reads saved fit summaries and raw in vitro inputs from:

- `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/<fit-name>/joint_fit_summary.tsv`
- `Data/in-vitro/pkpd_live_dead_model/processed/counts_by_well_time_wellAggregated.parquet`
- `Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx`
- `Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx`

The local `src/invitro_fitting.py` is a dependency snapshot used to rebuild dFdCTP signal surfaces, live/dead observations, model simulations, and new fit summaries.

## Usage

```bash
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
  Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/bestFitSoFar_20260513T164159
python3 Code/in-vitro/pkpd_live_dead_model/plot_invitro_fit_outputs.py \
  code/invitro_fitting_outputs/bestFitSoFar_20260513T164159 \
  --comparison-dose '25 nM'
```

The default input folder is `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/alsoGoodFit_20260514T093906`. Old miningcloneid-style positional paths containing `code/invitro_fitting_outputs/...` are resolved to the corresponding folder under `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/`.

The default `--fit-t-max` is `4.0`, and the default `--comparison-dose` is `50 nM`.

## Fitting

The fitting entrypoint is a thin wrapper around `src/invitro_fitting.py`:

```bash
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --check-inputs
python3 Code/in-vitro/pkpd_live_dead_model/run_invitro_fit.py --smoke-test --n-starts 1 --max-parallel 1
```

`--check-inputs` validates that all raw inputs resolve inside this repository, loads the count table, platemap, PKPD workbook, dFdCTP surfaces, modeling dataset, and aligned live/dead matrices, then writes `Data/in-vitro/pkpd_live_dead_model/input_manifest.tsv`.

`--smoke-test` runs the real fitting path serially with one `n_tr` value, one optimizer start, control plus 25/50 nM doses, one replicate per ploidy-dose group, and a short fitting window. It writes a timestamped output folder containing `joint_fit_summary.tsv` and `optimizer_attempts.tsv`; that folder can be passed back into `plot_invitro_fit_outputs.py`.
