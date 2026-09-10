# Datasets by manuscript figure

This is the entry point for finding data. The [dataset index](figure_datasets.tsv)
provides paths, roles, availability, and Manager module names. Small experimental
inputs and frozen tables belong here; generated analysis runs belong in
`Results/`; publication images belong in `figures/`.

## Figure-to-data map

| Figure / panels | Inputs | Manager module / notes |
|---|---|---|
| **1B-C**, GDSC drug-class enrichment; public-data support for Figure 2 | [GDSC response table](public_data/gdsc/raw/GDSC2_fitted_dose_response_24Jul22.txt), [cell-line ploidy](public_data/gdsc/raw/ploidyAcrossCellLines_V1.txt), and [reviewed drug-class workbook](public_data/gdsc/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx) | `gdsc`. The workbook's `Drug counts` sheet defines the displayed collapsed classes. These are frozen local inputs; no PubChem/DrugBank query is needed. |
| **2B**, CCLE breast-cell-line drug sensitivity versus ploidy | [CCLE inputs](public_data/ccle/): per-cell-line response TSVs, cell metadata, drug aliases, curated ploidy, and expression-availability column names | `ccle`. The default response is Z score. `raw/grbrowser/` supports an optional alternative analysis and is retained because code reads it. |
| **3D-F,H**, dose-response and association with ploidy | [Gemcitabine measurements](in-vitro/drug_response/Gemcitabine.txt), [frozen CLONEID ploidy](in-vitro/drug_response/fig3h_cloneid_ploidy.tsv), plus lineage and matching audit tables in the [same directory](in-vitro/drug_response/) | `drug_response`. The frozen export allows ordinary reproduction without access to CLONEID; refreshing it requires the external CLONEID infrastructure. |
| **4-5**, intracellular active-drug measurements, live/dead curves, fitted potency and latency; related supplementary fit panels | [Well-aggregated imaging counts](GemDelayKillTerm/processed/counts_by_well_time_wellAggregated.parquet), [field-level counts](GemDelayKillTerm/processed/counts_by_well_time.parquet), [plate map](in-vitro/pkpd_live_dead_model/raw/Gemcitabine_PlateMap_20240111.xlsx), and [PK assay workbook](in-vitro/pkpd_live_dead_model/raw/drugKinetics/GemcitabineExposure_PKPD.xlsx) | `pkpd`. Both count parquets are small and versioned. [Saved fitting summaries](in-vitro/pkpd_live_dead_model/invitro_fitting_outputs/) support fast plotting; `--mode full-refit --modules pkpd` fits from the count/PK inputs. |
| **6**, metabolomics PCA, response counts, pathway enrichment, heatmap, and volcano plots | [Metabolomics workbook](in-vitro/metabolomics/raw/Metabolomics_2N_4N_Full.xlsm) | `metabolomics`. This contains the feature-abundance input supplied to the analysis, not the instrument's raw mass-spectrometry files. [Processing description](../Code/Gemcitabine_Metabolomics_Heatmap/Readme.txt). |
| **7**, tumor growth and single-cell analyses | [Tumor-growth workbook](in-vivo/dt_Gem_VT_20241223_v4.xlsx), [sample metadata](in-vivo/sample_info.xlsx), [cell-level ploidy](in-vivo/all_ploidy.tsv), [processed cell tables and frozen analysis references](in-vivo/figure7/), [NUMBAT/CBS inputs](in-vivo/scRNAseq_Numbat/), and the SI cache below | `in_vivo_figure7`. Frozen tables/references are required for routine rendering; raw recomputation uses the pinned Zenodo inputs and generates candidates for review. |
| **SI4-7**, cellular landscape, composition, endpoint ploidy, cluster pathways | [Eleven frozen plotting tables and manifest](in-vivo/SIfigures/), plus [NUMBAT/CBS matrices](in-vivo/scRNAseq_Numbat/) | `si_figures`. Some tables also support panels promoted into main Figure 7, so uncited or unplotted tables are not automatically unused. |
| **SI8**, TGI sensitivity | [Frozen sensitivity references](in-vivo/figure7/saved_tgi_sensitivity/) and Figure 7 inputs | Materialized by `in_vivo_figure7`; retain its sensitivity tables as publication dependencies. |
| **SI9**, endpoint flow cytometry | [Selected FCS files, FlowJo workspace, crosswalk, and checksum manifest](in-vivo/flow_cytometry/endpoint_tumors_20250128/) | Optional `in_vivo_endpoint_flow` module. These are currently in Git, with a combined size of about 203 MiB; transfer to a documented deposit is still needed before removing them locally. |
| **SI10**, velocity/pseudotime view | Configured frozen-table location: `Data/in-vivo/figure7/processed/cellcycle_velocity_pseudotime_umap.tsv`, plus its provenance TSV; **both are absent from this checkout** | Optional `in_vivo_velocity_pseudotime`. Routine plotting requires this reviewed bundle; raw reconstruction uses the scRNA-seq inputs below. |

The table describes data for code-generated panels. It does not imply that the
repository contains every original microscopy image, immunoblot, external
NCI-60 dataset, or final manually assembled manuscript composite. The older
[Figure-Code Map](../docs/FigureCodeMap.md) records those gaps; use this guide for
current input locations and `figures/*/manifest.tsv` for current source assets.

## Large data and Zenodo

The configured scRNA-seq deposit is
[Zenodo record 21463392](https://zenodo.org/records/21463392), DOI
`10.5281/zenodo.21463392`. The
[pinned download manifest](../Code/in-vivo/figure7/zenodo_required_files.tsv)
lists 18 loom files and one integrated Seurat RDS (about 10.61 GiB total), with
sizes and checksums. The downloader reuses verified files and retrieves missing
ones into `Results/in-vivo/figure7/raw/zenodo_21463392/`.

```sh
Rscript Code/in-vivo/figure7/download_figure7_raw_data.R --help
bash Manager.sh --mode check-only --modules in_vivo_figure7,si_figures \
  --figure7-no-download-missing-raw
```

The Seurat RDS is already a processed object; this deposit boundary is not
equivalent to FASTQs. An earlier reconstruction path accepts external Cell Ranger
H5 inputs via `--figure7-cellranger-root`. See the
[Figure 7 workflow](../Code/in-vivo/figure7/README.md) for the exact boundaries.

No imaging-data Zenodo retrieval contract was found in this repository. Local
`GemDelayKillTerm/raw/` (~19 GiB), `in-vitro/Archive/` (~13 GiB), and the
object-level imaging parquet are **not covered by the configured scRNA-seq
download manifest**. The legacy preprocessing script
[`portScript.R`](GemDelayKillTerm/portScript.R) also refers to an external HALO
directory. These inputs were retained pending provenance/deposit review;
ordinary model fitting uses the small, versioned count parquet above.

Large data should be removed from the distributed checkout only after a deposit
has a persistent record, file list, checksums, and tested retrieval instructions.
Do not remove the sole known copy merely because it is not needed for plotting.
Moving/deleting tracked files does not erase their size from Git history.

## Checks and historical paths

```sh
python3 Code/tools/check_figure_datasets.py
```

This checks registered inputs for existence and Git coverage, and reports
optional local caches separately. It does not infer that an unregistered file
is safe to delete.

Historical run manifests record the paths used when those runs were created.
For their relocated inputs, use:

| Historical prefix/path | Current prefix/path |
|---|---|
| `Code/gdsc_ploidy_analysis/data/` | `Data/public_data/gdsc/` |
| `Code/ccle_ploidy_analysis/data/` | `Data/public_data/ccle/` |
| `Code/Gemcitabine_Metabolomics_Heatmap/Metabolomics_2N_4N_Full.xlsm` | `Data/in-vitro/metabolomics/raw/Metabolomics_2N_4N_Full.xlsm` |

The [cleanup record](../docs/data_cleanup_20260910.md) lists removals and retained
archives. Old run metadata was not rewritten to pretend the relocation occurred
at the time of the original analysis.
