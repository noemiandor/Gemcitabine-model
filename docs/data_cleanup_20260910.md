# Dataset organization and cleanup, 2026-09-10

Starting revision: `a5efb28749d88ca2e516f57c1505547250fda8d6` on `main`.
Pre-existing Figure 7 code/asset modifications were left in place.

## Relocation

GDSC and CCLE inputs moved from their `Code/*/data/` directories into
`Data/public_data/gdsc/` and `Data/public_data/ccle/`. The Figure 6 workbook moved
to `Data/in-vitro/metabolomics/raw/`. Their content is unchanged. Manager,
standalone analysis defaults, the GDSC workbook test, and current documentation
now use the new locations. Historical run manifests retain original paths; the
dataset guide supplies the path mapping.

Five existing local Figure 3 input/export files under
`Data/in-vitro/drug_response/` were previously untracked. They are included in
the cleanup change so fresh clones receive the dose-response measurements,
frozen ploidy export, lineage map, and matching/detail audits.

## Removed

| Path | Tracked files | Reason |
|---|---:|---|
| `Data/K01_WGDClassification_050124/` | 92 | No reference found in repository code, configuration, figure metadata, or documentation; outside the paths enumerated by the retained imaging importer. |
| `Code/gdsc_ploidy_analysis/data/raw/small_molecule_20200407234909.csv` (temporarily relocated during cleanup) | 1 | Retired drug annotation input. Current GDSC code uses the reviewed XLSX exclusively. Only an old run manifest and archived planning documents referenced this CSV. |
| `Results/public_data/gdsc_ploidy_analysis/runs/primary_secondary_candidate_gdsc/` | 76 | Superseded result snapshot; no consumer outside the snapshot and its own manager record. Current Figure 1 source runs are retained. |
| `Results/public_data/gdsc_ploidy_analysis/runs/primary_secondary_final_gdsc/` | 75 | Same; the word `final` in this historical run name does not indicate the current publication source. |
| `Results/manager/runs/primary_secondary_candidate/` | 1 | Manager record for removed candidate snapshot. |
| `Results/manager/runs/primary_secondary_final/` | 1 | Manager record for removed historical final snapshot. |

These are tracked-file removals, recoverable from the starting Git revision.
No history rewrite is performed. Other retained historical snapshots may refer
to retired paths as provenance; that is not an active analysis dependency.

## Retained deliberately / remaining deposit work

- `Data/GemDelayKillTerm/processed/counts_by_well_time*.parquet` are model inputs
  and their immediate precursors, not obsolete outputs.
- `Data/GemDelayKillTerm/raw/` (~19 GiB) and `Data/in-vitro/Archive/` (~13 GiB)
  participate in, or preserve sources for, the legacy imaging preparation path.
  That path uses recursive directory discovery and an external HALO directory.
  Their absence from Manager's immediate input list does not establish that
  they can be discarded. No configured imaging deposit was found.
- `Data/wgd_calls/`, the `RealTimePseudotimeCorrelation*.txt` files, and other
  untracked experimental exports are retained as unresolved archival data.
  The archived tracking code refers to WGD result basenames; deleting them
  solely because their directory is not named in code would be unsafe.
- `Data/M00_GemcitabinePKPD_101823/` is referenced by archived PK plotting code;
  it is retained rather than reported as unused by every script.
- The endpoint FCS/workspace dataset is used by the SI9 replay workflow. Its
  ~203 MiB should be deposited with a tested retrieval contract before removal.
- Reviewed Figure 7 references, older references used by tests/audits, the 11
  SI tables, and active `figures/*/manifest.tsv` source runs are retained.
- Optional CCLE `grbrowser` inputs still have a supported code path.
- Untracked/ignored run directories are not blanket-deleted: several are
  publication sources or contain unique exploratory analyses. New large-data
  publication or wholesale local-archive deletion is outside this cleanup.

The public guide is `Data/README.md`; the curated inventory is
`Data/figure_datasets.tsv`. Its checker reports missing and untracked registered
inputs; it intentionally does not equate no textual reference with permission
to delete a file. The cleanup is conservative where dependencies are unresolved.

## Validation

- All 35 relocated tracked files are byte-identical to their original Git blobs.
- The curated dataset checker passes for 26 registered datasets, reporting the
  missing optional SI10 table/provenance explicitly.
- Manager input checks pass for all seven default modules; shell syntax passes.
- GDSC workbook/smoke tests and CCLE smoke tests pass.
- CCLE reproduces all five baseline scientific-table/run-configuration hashes.
  The old strict baseline command reports differences only in
  `metadata/input_manifest.tsv` (relocated paths) and `metadata/session_info.txt`
  (environment); these historical baselines were not overwritten.
- Manager regenerated all three metabolomics submodules and Figure 6 assets
  into `/private/tmp/gemcitabine-dataset-cleanup-analysis` and
  `/private/tmp/gemcitabine-dataset-cleanup-figures`. All 15 published source
  CSV/TSV result tables match the previous run byte-for-byte. Their input,
  output, and figure manifests validate.
- The first metabolomics attempt used the default Homebrew Python, which lacks
  pandas. The successful run used the existing Miniconda Python environment;
  no packages or environments were added.
- The original Figure 3 measurement file has trailing empty TSV fields. Those
  represent missing measurements and were preserved, despite Git's whitespace
  warning; scientific inputs were not reformatted.
