# Figure 7 publication polishing notes

The manuscript composite is rebuilt natively at 7.1 x 10.645 inches from live
plot and grob objects. Its explicit A--L map is A=7A, B=7C, C=SI4A, D=SI4B,
E=SI4C, F=SI4E, G=SI7B, H=the ECDF component of 7B, I=7F, J=the QC-filtered
NUMBAT copy-number heatmap, K=7D, and L=7E.

The detailed treated-minus-vehicle localization formerly compressed beneath H
is regenerated as full-width Supplementary Figure 4I from the exact reviewed
2,881-cell CellCycle universe. Main H retains only the equal-mouse ECDF
comparisons. C--E are enlarged to 1.35 inches high. I and J share the
penultimate row equally; J retains all 22 chromosome columns while displaying
only chromosomes 1, 5, 9, 13, 17, and 22 as tick labels to avoid print-size
collisions.

Main J is generated from all 16 checksum-pinned CBS matrices and restricted by
exact cell identity to the 9,832-cell final-QC tumor universe, including 5,335
treated cells. Clusters 3, 4, 9, and 9c remain excluded. Its row annotations
encode injected origin, gemcitabine dose, and mouse. The former supplementary
copy-number heatmap is removed; Supplementary Figure 6 now contains panels A--E,
with the injected-reference comparison as E.

Panel-internal titles, formulas, and methodological prose remain in the
manuscript legend rather than the composite. The bundled optimizer is retained
as an independent diagnostic, while `layout/layout_plan.csv` records the
scientifically constrained A--L arrangement.

Project-map decision: `docs/FigureCodeMap.md` remains the repository-wide code
map. This package adds a local A--L identity map and does not replace the
project map. Raster subpanels are audit exports only; no raster subpanel is
loaded into the final composite.
