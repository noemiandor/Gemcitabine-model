# Regression Report

- Baseline: `Code/gdsc_ploidy_analysis/output_milestones/phase2_refresh_entrypoint`
- Candidate: `Code/gdsc_ploidy_analysis/output_milestones/phase3_structured_pubchem`
- Final result changed: `no scientific-result change; regenerated binary artifacts differ`
- Difference count: `1`

## Difference Summary

- `checksums.tsv`: changed (md5 38fe01ed833fba60b4cc34712014d38f -> 2be450f745dd5598261ce88b615b2fcb)

## Attribution

The first and only tracked difference is `checksums.tsv`. `drugsVsPloidyCorr.RData`, `coxIn.RData`, exported correlations, saved annotations, group counts, and workbook sheet TSV exports match Phase 2.

The changed hashes are limited to regenerated `drugsVsPloidyCorr.xlsx`, `drugsVsPloidyCorr.pdf`, and `ploidyVsDrugSensitivity.pdf`. The main analysis still consumed the same cached annotations, so the structured PubChem refresh implementation did not change the scientific result.

Local PubChem parser fixture tests passed with `Rscript Code/gdsc_ploidy_analysis/tests/test_pubchem_client.R`, covering `Drug Classes`, `Mechanism of Action` fallback, unclassified responses, and malformed JSON.
