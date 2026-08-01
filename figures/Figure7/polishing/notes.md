# Figure 7 publication polishing notes

The strict A--K identity and scientific analyses are unchanged. The display is
rebuilt natively at 7.1 x 9.7 inches instead of scaling a 24 x 28 inch canvas.
The six-row plan preserves row-major order while allocating full width to panel
I, 4.45 inches to G, and 3.46 inches each to J and K after the outer gutters.
Titles, subtitles,
formulas, and methodological prose are removed from final panel objects and
transferred to the manuscript legend.

The bundled optimizer is run as a required independent layout diagnostic. Its
candidate is retained, while `layout/layout_plan.csv` records the adopted
scientifically constrained arrangement because A--K reading order, the paired
UMAPs, the F/G relationship, and shared-width J/K comparison are part of the
panel identity contract.

Project-map decision: `docs/FigureCodeMap.md` remains the repository-wide code
map. The polishing package adds a local A--K identity map and does not renumber
or replace the project map.

No raster subpanel is loaded into the final composite. Audit PNGs are generated
for inspection only; the final PDF and PNG are assembled from reconstructed
ggplot and pheatmap grob objects in memory.
