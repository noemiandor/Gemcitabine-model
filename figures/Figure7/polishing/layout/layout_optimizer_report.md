# Panel layout optimization report

- Input dimensions: `figures/Figure7/polishing/layout/layout_optimizer_groups.csv`
- Layout plan: `figures/Figure7/polishing/layout/layout_optimizer_candidate.csv`
- Target width: 7.10 in
- Maximum height: 10.64 in
- Gap: 0.08 in
- Coordinate convention: `x_in`/`y_in` and `x_npc`/`y_npc` are lower-left panel origins. Use `y_npc` directly; do not invert it during assembly.

## Figures

- `Figure 7`: 6.23 x 10.65 in, wasted 3.8%, tree `[[ab / [[[cde / fg] / h] / ij]] / kl]`

## Scale recommendations

Panels with `sx` or `sy` far from 1 should be resized in the subpanel-generation script, then dimensions should be regenerated and this optimizer rerun before final assembly.
