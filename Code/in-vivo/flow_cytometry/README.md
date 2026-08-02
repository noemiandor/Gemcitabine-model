# Endpoint tumor flow-cytometry extraction

The endpoint-flow module has two deliberately separate layers:

1. `extract_endpoint_flow.py` turns the reviewed FlowJo workspace into
   portable, deterministic frozen-count TSVs. It uses only the Python standard
   library; population membership and counts come from the workspace, while
   each FCS file supplies acquisition metadata and an independent `$TOT`
   consistency check.
2. `reconstruct_endpoint_flow.R` reads the raw FCS event matrices and replays
   each selected sample's exact workspace geometry in parent-to-child order.
   It reports replayed counts alongside the frozen FlowJo values and renders a
   review composite; it never replaces a FlowJo count to manufacture parity.

`run_endpoint_flow.sh` runs both layers in that order.

The event renderer requires R plus `flowCore`, `xml2`, `ggplot2`, `patchwork`,
and `ragg`. It requires `ragg` rather than silently changing PNG devices, uses
`cairo_pdf` for the vector figure, and records the R/package/device/Cairo
versions in `metadata/run_config.tsv`.

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

Run the complete module from the repository root:

```bash
bash Code/in-vivo/flow_cytometry/run_endpoint_flow.sh \
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
- `endpoint_flow_count_agreement.tsv`: frozen FlowJo and raw-event-replay
  counts for the full hierarchy and the named peak, 2N, and 4N gates. The
  fixed, documented numerical tolerance accommodates small implementation-level
  boundary/statistic differences without assigning them a specific cause;
  every delta and status remains visible. The module fails rather than updating
  its successful-run pointer if any replay falls outside this tolerance.
- `endpoint_flow_gate_geometry.tsv`: the parsed numeric per-sample polygon
  vertices and rectangle bounds read from the workspace.
- `endpoint_flow_dna_histograms.tsv`: fixed 1,000-a.u. raw DNA-channel bins
  from 25,000 through 140,000. Probability mass sums to one independently
  within each mouse; fluorescence values are not centered, peak-aligned, or
  rescaled.
- `endpoint_flow_named_2n_4n_summary.tsv`: per-mouse frozen and replayed counts
  and HumanCells-parent percentages for the workspace-named 2N and 4N sibling
  gates. These narrow gates are nonexhaustive and are not called cell-state
  proportions.
- `endpoint_flow_reconstruction_per_sample.tsv` and
  `endpoint_flow_representative_selection.tsv`: event-universe, raw-range,
  low-count, and deterministic representative-selection audit tables.
- `metadata/run_config.tsv`: the linear-transform, identity-spillover,
  normalization, raw-axis, tolerance, and package-version contract.
- `figures/panel_SuppFig9_endpoint_flow_cytometry.{pdf,png}`: a 7.1-by-9-inch
  review composite containing the representative raw-event gate hierarchy,
  all 16 within-mouse-normalized DNA-content distributions, and the per-mouse
  named 2N/4N gate summaries. The 174-HumanCells sample is retained and marked
  with a dagger.

The representative scatter panels display at most 10,000 deterministic,
evenly spaced event indices for legibility; their gate counts use every event.
The histogram and named-gate panels use the full replayed HumanCells sets.

For the selected files, FSC-A, SSC-A, and `450/50 Violet B-A` all use the
workspace's linear 0--262,144 transform with gain 1, and the FCS spillover
matrix is the identity. The renderer verifies those statements for every
sample. It does not apply the FACSDiva `P7BS` display-scaling keyword as an
event offset. The common 25,000--140,000 display window contains every replayed
HumanCells DNA event; the complete acquisition range remains recorded in the
FCS data.

Percentages in the frozen-count extractor are calculated from integer
workspace counts and emitted to six decimal places. Analysis TSVs and the PNG
contain no run timestamp or absolute path and are byte-identical for identical
inputs in the pinned software environment. Cairo embeds a PDF `CreationDate`,
so independently regenerated PDFs can differ bytewise despite identical visual
and scientific content; manifests therefore hash each concrete PDF instance.

Run the dependency-free tests without producing bytecode artifacts:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest discover \
  -s Code/in-vivo/flow_cytometry/tests -v
```

Run the event-level integration test, which replays all 16 pinned files and
checks the figure/table contracts, with:

```bash
scripts/agentRrunner.sh \
  Code/in-vivo/flow_cytometry/tests/test_reconstruct_endpoint_flow.R
```
