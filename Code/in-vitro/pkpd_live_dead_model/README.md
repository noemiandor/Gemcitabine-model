# PKPD Live/Dead Model Figure Reproduction

This module reproduces manuscript-style in vitro PKPD/live-dead model panels from saved fitting summaries. It was migrated from the working-tree version of `/Users/4470246/Repositories/miningcloneid/code/plot_invitro_fit_outputs.py`.

## Inputs

The plotter reads saved fit summaries and raw in vitro inputs from:

- `Data/in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/<fit-name>/joint_fit_summary.tsv`
- `Data/in-vitro/pkpd_live_dead_model/processed/counts_by_well_time_wellAggregated.parquet`
- `Data/in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx`
- `Data/in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx`

The local `src/invitro_fitting.py` is a dependency snapshot used to rebuild dFdCTP signal surfaces, live/dead observations, and model simulations from the saved summary rows. The broader fitting workflow can be formalized in this folder later without changing the figure-reproduction entrypoint.

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
