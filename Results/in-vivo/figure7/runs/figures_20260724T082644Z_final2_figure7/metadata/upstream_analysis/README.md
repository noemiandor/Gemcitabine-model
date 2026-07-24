# Figure 7 upstream analysis documentation

The Seurat RDS and loom files consumed by this workflow are processed artifacts with upstream data-generation and analysis steps. The reviewed descriptions, parameters, software/runtime records, and checksums for those preceding steps are archived at Zenodo DOI 10.5281/zenodo.21463392.

This Manager run begins from the processed Seurat RDS and loom artifacts. It does not rerun Cell Ranger, Seurat integration/refinement, or velocyto loom generation. The indexed Zenodo documents describe those preceding steps and remain the authoritative copies.

- Zenodo DOI: [10.5281/zenodo.21463392](https://doi.org/10.5281/zenodo.21463392)
- Record: [Resource limitation rewires chromosome instability and ploidy evolution across in vitro and in vivo cancer models](https://zenodo.org/records/21463392)
- Indexed documents: 11
- Consuming module: `in_vivo_figure7`

## Published index

| Role | File | Description |
|---|---|---|
| overview | [readme.md](https://zenodo.org/api/records/21463392/files/readme.md/content) | Reviewed description of the Seurat and loom upstream analysis steps and recorded reproduction commands. |
| checksum_manifest | [checksum.md](https://zenodo.org/api/records/21463392/files/checksum.md/content) | SHA-256 checksum inventory for the deposited data and documentation. |
| seurat_provenance | [integrated_sct_cca_seurat_final_reclustered_provenance.csv](https://zenodo.org/api/records/21463392/files/integrated_sct_cca_seurat_final_reclustered_provenance.csv/content) | Structured provenance for Seurat generation, QC, integration, annotation, filtering, clustering, and reductions. |
| loom_parameters | [loom_generation_parameters.csv](https://zenodo.org/api/records/21463392/files/loom_generation_parameters.csv/content) | Per-sample velocyto inputs, commands, parameters, logs, and generated loom metadata. |
| loom_runtime | [loom_generation_runtime_environment.csv](https://zenodo.org/api/records/21463392/files/loom_generation_runtime_environment.csv/content) | Runtime environment recorded for loom generation. |
| runtime_overview | [runtime_environment.md](https://zenodo.org/api/records/21463392/files/runtime_environment.md/content) | Human-readable software and runtime provenance for the upstream analysis. |
| runtime_software | [runtime_environment_software.csv](https://zenodo.org/api/records/21463392/files/runtime_environment_software.csv/content) | Structured software and executable version inventory. |
| runtime_packages | [runtime_environment_packages.csv](https://zenodo.org/api/records/21463392/files/runtime_environment_packages.csv/content) | Structured package inventory for the final-object environment. |
| runtime_session | [runtime_environment_sessionInfo.txt](https://zenodo.org/api/records/21463392/files/runtime_environment_sessionInfo.txt/content) | R sessionInfo for the final-object inspection environment. |
| upstream_packages | [runtime_environment_upstream_01_data_packages.csv](https://zenodo.org/api/records/21463392/files/runtime_environment_upstream_01_data_packages.csv/content) | Package inventory for the upstream 01_data Seurat integration stage. |
| upstream_session | [runtime_environment_upstream_01_data_sessionInfo.txt](https://zenodo.org/api/records/21463392/files/runtime_environment_upstream_01_data_sessionInfo.txt/content) | R sessionInfo for the upstream 01_data Seurat integration stage. |

The accompanying `zenodo_document_manifest.tsv` records the size and Zenodo MD5 for every indexed document. `provenance.tsv` records the pinned source-manifest SHA-256 and publication contract.
