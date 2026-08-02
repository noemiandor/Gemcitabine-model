# Endpoint-tumor flow cytometry (2025-01-28 acquisition)

This directory freezes the selected flow-cytometry evidence used to assess endpoint DNA-content state in the same 16 tumors represented in the manuscript's single-cell analysis. The portable identity contract is [`crosswalk.tsv`](crosswalk.tsv). It maps each canonical manuscript mouse to exactly one selected FCS acquisition, the exact sample name in the FlowJo workspace, the exact peak-population name to extract, and a separately reviewed numeric peak annotation.

## Layout and source

- `fcs/` contains the 16 selected FCS files. `fcs_file` in the crosswalk is deliberately a bare basename; extraction code resolves it against an explicit `--fcs-dir` rather than a machine-specific path.
- `paired_acquisition_sensitivity.tsv` explicitly maps the four unresolved `_M`/`-M` workspace pairs. Their raw FCS files are not imported, so this sensitivity is limited to counts and annotations already frozen in the reviewed workspace.
- `workspace/20250129_TumorSamples.wsp` is the reviewed FlowJo 10.10.0 workspace. Its internal URIs retain the original analyst's absolute paths and must not be used as portable file locations. Extraction joins workspace samples by the exact `wsp_sample_name` values in the crosswalk.
- `content_manifest.tsv` pins the dataset-relative path, byte size, and SHA-256 digest of the crosswalk, paired-sensitivity mapping, workspace, and all 16 selected FCS files. Extraction fails before analysis if any selected input differs.
- The read-only source was `/Volumes/Flow_Cytometry/20250128_SUM159_2Nm_4Nm_TumorSamples`. This volume path records provenance only and is not a runtime dependency.

The source directory contains 42 FCS files. Sixteen were selected for the manuscript mice. The workspace is retained intact, so it also refers to excluded files that are intentionally absent from `fcs/`.

## Selection and identifier policy

The selected universe is fixed to the 16 mice in the canonical `Data/in-vivo/sample_info.xlsx` metadata: four vehicle, two 30 mg/kg, and two 120 mg/kg tumors within each injected-origin group.

- For the four A1 tumors, the `GENOMICS` aliquot was selected because it corresponds most closely to the sequenced-tumor context. The separate `TISSUE` and `VT` acquisitions are recorded as excluded alternates.
- For A2-0 and A2-L, the initial acquisitions were selected. The later duplicate acquisitions are recorded as excluded alternates.
- Where both an unsuffixed acquisition and an `_M`/`-M` acquisition exist for a selected 4N-origin tumor, the unsuffixed acquisition was selected by a fixed filename rule. Nothing in the FCS files or workspace establishes whether `M` denotes a technical repeat, a distinct aliquot, or another specimen state, so the paired files are not called technical replicates. The scripted paired-workspace sensitivity retains the same qualitative near-2N conclusion; the selected and substituted descriptive annotation means are 2.05N and approximately 2.09N, respectively, but neither is treated as a calibrated ploidy estimate or emitted as a group measurement.
- A4-R, A4-RL, A5-L, A5-R, A6-0, and A6-RR each have one applicable acquisition in the workspace.

Twelve additional source FCS files represent tumors outside the canonical 16-mouse manuscript universe: A5-RL, A6-R, A6-RL, A7-L, A8-L, and A8-R, each with its paired acquisition where present. They were not imported and must not be substituted for a canonical mouse.

Two filename differences are resolved only through explicit crosswalk rows; extraction code must not guess them:

- Canonical `2N-A1-LR` corresponds to cytometry specimen `A1-RL` (LR/RL normalization).
- Canonical `A6-4N-O` corresponds to cytometry specimen `A6-0` (letter O/digit 0 normalization).

The differing `2N-A...`, `4N-A...`, and `A...-4N...` token orders are also preserved exactly. `mouse_id`, injected origin, and dose come from the canonical manuscript metadata, not from inference over FCS filenames.

## Dates and stored-material caveat

All selected FCS headers report `$DATE=28-JAN-2025`, represented as `2025-01-28` in the crosswalk. The workspace was last modified on 2025-02-04.

The FCS files and workspace do not record tumor harvest dates. The canonical `sample_info.xlsx` table records harvest identifiers but not calendar dates, so `harvest_date` is explicitly `NA`; a conflicting legacy absolute-date export was not used to fill this field. Harvest dates should be added only from a reviewed primary record.

The tumors were harvested before the January 2025 acquisition, but the available folder does not document storage duration, preservation, thawing, staining, or acquisition-preparation details. These flow results therefore corroborate endpoint DNA-content state subject to a stored-material caveat; they do not establish how storage affected recoverable events or fluorescence profiles.

## Gate interpretation

`peak_gate_name` is the exact name of an existing FlowJo population. The names are analyst-supplied labels, not ploidy values recalculated from FCS intensities by the repository extractor. `reviewed_peak_annotation_n` carries a numeric transcription of the reviewed annotation for per-sample display; extraction code must not infer it by parsing the gate-name string, and group-level ploidy means are intentionally not produced.

- The A1 2N-origin samples have a generic `2N` peak population. The selected A2 and A4 samples use the exact workspace population name `2N 50`. Both map to the separately reviewed value `2.00`; neither is a sample-specific numerical peak estimate computed by the extractor.
- The 4N-origin samples have manually named peak populations from `1.88N` through `2.2N`. The exact workspace spelling `2.2N` is retained even when a report formats the same value as `2.20N`.
- These labels and their gate counts are not directly interchangeable with NUMBAT cell-level ploidy estimates or chromosome-count-derived injected-state proxies.

The workspace also contains `mCh+\2N`, `mCh+\4N`, `mCh-\2N`, and `mCh-\4N` gates. Those counts may be extracted as supporting summaries, but they are distinct from the manually labelled peak populations specified by `peak_gate_name`.

The sum of the `mCh+\2N` and `mCh+\4N` gate counts is not assumed to equal all mCherry-positive cells. Any composition calculated from those two counts is labelled explicitly as a percentage of their summed gate counts, not as a percentage of all mCherry-positive events.

The selected A1-0 GENOMICS file contains only 174 `HumanCells` events (115 in its `2N` peak gate). The extractor flags samples below its predeclared 1,000-`HumanCells` descriptive-QC threshold but does not silently exclude them. Group output reports the all-sample and QC-pass-only mean peak-gate percentages separately.

Exported histogram images were not imported because they are cropped presentation derivatives and cannot replace scripted extraction from the FCS/workspace pair.

The repository now performs that scripted event-level reconstruction in
`Code/in-vivo/flow_cytometry/reconstruct_endpoint_flow.R`. The relevant
workspace transforms are linear with unit gain and the embedded FCS spillover
matrices are identity, so the renderer applies the per-sample workspace-defined
polygons and rectangles directly on recorded fluorescence values. Small replay
differences are retained in an explicit agreement table rather than being
forced to equal the workspace's frozen FlowJo counts; no specific software
mechanism is asserted as their cause.
