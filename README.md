# Gemcitabine model and manuscript figures

Start with the [figure-by-figure dataset guide](Data/README.md) to find the
experimental inputs, public datasets, frozen analysis tables, and large-data
retrieval instructions. A machine-readable [dataset index](Data/figure_datasets.tsv)
lists the same input locations and their availability.

`Manager.sh` runs the manuscript modules and writes analysis results to
`Results/` and publication assets to `figures/`. Each figure's `manifest.tsv`
records its source outputs. See the
[module registry](docs/manuscript_figure_module_registry.tsv) for entrypoints.

Check the inputs before running an analysis:

```sh
python3 Code/tools/check_figure_datasets.py
bash Manager.sh --mode check-only --modules gdsc,ccle,drug_response,pkpd,metabolomics
```

For example, regenerate Figure 6's analysis and source panels with:

```sh
bash Manager.sh --modules metabolomics --run-id my_figure6
```

Some manuscript panels are externally assembled images, immunoblots, or
schematics; the dataset guide identifies the limits of repository coverage.
