# Endpoint tumor flow-cytometry extraction

`extract_endpoint_flow.py` turns the reviewed FlowJo workspace into portable,
deterministic TSVs. It uses only the Python standard library. It does **not**
reimplement FlowJo transformations or gates: population membership and counts
come from the reviewed workspace, while each FCS file supplies acquisition
metadata and an independent `$TOT` consistency check.

## Explicit input contract

The crosswalk is a tab-delimited file with exactly 16 data rows and these
required columns:

| column | meaning |
|---|---|
| `mouse_id` | unique manuscript mouse identifier |
| `injected_origin` | `2N` or `4N` injected lineage |
| `dose_mg_kg` | nonnegative numeric dose; vehicle is `0` |
| `fcs_file` | portable path relative to `--fcs-dir` |
| `wsp_sample_name` | exact, unique FlowJo `SampleNode` name |
| `peak_gate_name` | exact direct child of `HumanCells`, such as `1.88N` or `2N 50` |
| `reviewed_peak_annotation_n` | positive numeric transcription of the reviewed annotation; not a calibrated ploidy estimate |
| `acquisition_date` | reviewed ISO date that must agree with FCS `$DATE` |

The extractor fails rather than guessing if sample names or files are
duplicated, the FCS filename disagrees with `$FIL` or the workspace URI, FCS
`$TOT` disagrees with the workspace sample count, or any required gate is
missing. `HumanCells`, `2N`, `4N`, `mCh+\2N`, `mCh+\4N`, and the specified
peak gate must each resolve unambiguously. The extractor never guesses ploidy
from a gate label; it reports the crosswalk's explicit reviewed annotation.
The workspace hierarchy is also enforced exactly as `Human Cell Enrichment`
then `HumanCells`; both parent counts and their percentages are reported.
The named populations may
overlap; their counts are therefore reported as reviewed and are not summed to
reconstruct `HumanCells`. Only explicitly crosswalked SampleNodes are inspected
for populations; exploratory or control nodes elsewhere in the same workspace
cannot silently enter the cohort or block it because of incomplete old gates.
Each FCS file must identify `FACSCantoII` in `$CYT` and must contain exactly one
parameter named `450/50 Violet B-A` (observed as `$P7N`). Both expectations can
be changed explicitly with `--expected-cytometer` and
`--expected-dna-channel`; the observed values and parameter index are written
to the per-sample table. The pinned input manifest must match the crosswalk,
workspace, and all selected FCS files by dataset-relative path, byte size, and
SHA-256 before gate extraction proceeds. Samples below 1,000 `HumanCells` are
flagged descriptively by default but are not silently excluded.

Run from the repository root:

```bash
python3 Code/in-vivo/flow_cytometry/extract_endpoint_flow.py \
  --crosswalk Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/crosswalk.tsv \
  --paired-sensitivity Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/paired_acquisition_sensitivity.tsv \
  --workspace Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/workspace/20250129_TumorSamples.wsp \
  --expected-input-manifest Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/content_manifest.tsv \
  --fcs-dir Data/in-vivo/flow_cytometry/endpoint_tumors_20250128/fcs \
  --output-dir Results/in-vivo/flow_cytometry/endpoint_tumors
```

Outputs are:

- `endpoint_flow_per_sample.tsv`: counts, parent-relative percentages,
  mCherry-positive 2N/4N gate-count composition, acquisition metadata, QC
  flags, and the numeric peak annotation explicitly recorded in the crosswalk.
- `endpoint_flow_group_summary.tsv`: mouse-balanced summaries overall, by
  injected origin, and by origin and dose. It does not average nominal peak
  annotations as though they were calibrated measurements. The mCherry fields
  explicitly use the sum of the named 2N and 4N gate counts as denominator;
  they do not claim that those gates exhaust all mCherry-positive events.
- `input_hashes.tsv`: SHA-256 and byte size for the crosswalk, workspace, and
  all 16 FCS files, plus the paired-sensitivity selection table. Paths are
  logical/relative rather than machine-specific.
- `paired_workspace_sensitivity.tsv`: selected-versus-paired manual peak gates
  and counts for the four unresolved `_M`/`-M` workspace nodes. It explicitly
  records that the paired raw FCS files were not imported or header-validated.

Percentages are calculated from integer workspace counts and emitted to six
decimal places. No run timestamp or absolute path is embedded, so identical
inputs yield byte-identical outputs.

Run the dependency-free tests without producing bytecode artifacts:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s Code/in-vivo/flow_cytometry/tests -v
```
