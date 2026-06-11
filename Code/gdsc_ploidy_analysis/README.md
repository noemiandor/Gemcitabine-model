# GDSC Ploidy Analysis

This folder packages the GDSC-versus-ploidy analysis used in the gemcitabine manuscript into a self-contained module.

## Contents

- `run_gdsc_ploidy_analysis.R`: main entrypoint
- `refresh_pubchem_annotations.R`: explicit network-dependent annotation refresh entrypoint
- `annotate_from_pubchem.R`: local drug-category annotation helper
- `custom_set_candidate.tsv`: reconstructed manual category overrides
- `GDSC_analysis_original.R`: original script with its historical external-path assumptions
- `GDSC_analysis_repro_reference.R`: repaired reference version developed during reproduction
- `data/`: local input files required by the analysis
- `output/`: generated results

## Inputs

The analysis uses:

- `GDSC2_fitted_dose_response_24Jul22.txt`
- `ploidyAcrossCellLines_V1.txt`
- `small_molecule_20200407234909.csv`

## Outputs

Running the analysis writes:

- `output/drugsVsPloidyCorr.pdf`
- `output/drugsVsPloidyCorr.RData`
- `output/coxIn.RData`
- `output/drugsVsPloidyCorr.xlsx`
- `output/ploidyVsDrugSensitivity.pdf`

## R package dependencies

- `EnrichIntersect`
- `xlsx`
- `plyr`
- `ChemmineR`
- `textreadr`

## Notes

- `custom_set_candidate.tsv` is a reconstructed manual override table inferred from surviving outputs and rerun behavior.
- The original script relied on an undefined in-memory object called `custom.set`; this package makes that dependency explicit.
- The main analysis uses `data/derived/pubchem_drug_annotations.tsv` by default and does not contact PubChem during normal runs.
- To refresh PubChem annotations, run `Rscript refresh_pubchem_annotations.R --output-cache=data/derived/pubchem_drug_annotations.refresh.tsv` and review the TSV diff before replacing the canonical cache.
- The PubChem helper is a near-verbatim local port of the original script, with only path localization and a local `isempty()` shim.
- `matlab` is no longer required by this module because `isempty()` is implemented locally.
