# Manuscript Figure Style

## Status and scope

This file is the repository-wide visual contract for manuscript-facing figures.
It is a journal-neutral baseline derived from the current hypoxia LTEE figures.
Record any target-journal override in the figure package README or rebuild
manifest rather than silently changing these rules.

These rules apply at the final displayed size. A panel that is readable only
when viewed substantially larger than its manuscript size does not pass visual
QC.

## Output dimensions and formats

- Build figures at their intended final width: 3.35 in for one column, 5.2 in
  for an intermediate width, or 7.1 in for two columns.
- Keep a full figure at or below 9.0 in high unless the target journal permits
  a larger page.
- Prefer vector PDF for plots, text, and schematics. Also export a white-
  background PNG at 300 dpi for review and manuscript assembly checks.
- Use 600 dpi for rasterized line art when a vector output is impossible.
- Do not resize a completed composite to compensate for undersized text.
  Regenerate it at the intended dimensions.
- Do not place a figure number, manuscript caption, or redundant figure-level
  title inside the image.

## Typography

- Use a sans-serif family consistently: Arial when available, then Helvetica,
  Liberation Sans, or the R `sans` family.
- Use regular weight for axis and legend text. Reserve bold for panel labels
  and short structural headings.
- At final size, use approximately:
  - panel labels: 11-12 pt, bold;
  - short panel headings: 9-10 pt;
  - axis titles: 8-9 pt;
  - tick and legend text: 7-8 pt;
  - annotations: at least 7 pt.
- Do not use text smaller than 7 pt in a manuscript-facing panel.
- Keep capitalization sentence-style. Avoid all-caps headings.

## Panel labels and layout

- Label panels `A`, `B`, `C`, and so on, without punctuation.
- Place labels at the upper-left of each panel, aligned across rows and outside
  the plotting data region whenever possible.
- Use one consistent label size, weight, and offset across the figure set.
- Follow left-to-right, top-to-bottom reading order.
- Use stable panel dimensions and shared axis limits when panels are intended
  for direct comparison.
- Align plotting regions, not merely image boundaries.
- Use one shared legend for repeated encodings. Place it in reserved whitespace
  rather than over data.
- Avoid nested decorative frames, shadows, and background cards.

## Plot theme

- Use a white background with dark gray or black text.
- Prefer a restrained `theme_classic()` or lightly gridded `theme_bw()` style.
- Use major grid lines only when they aid quantitative reading; draw them in a
  light neutral gray such as `#E5E5E5`.
- Remove minor grid lines unless they carry necessary scale information.
- Use axis lines and ticks in `#333333`. Do not use a full heavy border around
  every panel unless it is needed to separate a heatmap or image.
- Avoid chart subtitles that repeat information provided by the caption.

## Canonical semantic encodings

Use the same encoding across figures whenever the same biological distinction
appears. If several distinctions occur in one panel, ploidy identity takes
precedence over context or model/observation color.

### Ploidy lineage

- 2N: blue `#0072B2`.
- 4N: vermilion `#D55E00`.
- When color is unavailable, distinguish 2N and 4N with line type or point
  shape in addition to labels.

### Experimental context

Use these only when color is not already encoding ploidy:

- in vitro: sky blue `#56B4E9`;
- in vivo: bluish green `#009E73`.

### Oxygen exposure

Lower oxygen should be darker. The normoxic endpoint must remain visibly blue,
not white:

- 20.5% O2: `#8FC5E3`;
- 5% O2: `#4292C6`;
- 1% O2: `#2171B5`;
- 0% O2: `#084594`.

For continuous oxygen, interpolate through this sequence. Give pale filled
regions a visible outline when they touch a white background.

### Observations and model output

- Prefer geometry over a new color scale: observations as points or intervals;
  fitted or predicted values as lines or distributions.
- Use filled points for observed values and solid lines for the primary model
  prediction. Use dashed lines for a comparison, counterfactual, or alternate
  initial condition.
- Do not present interpolated model values with the same point geometry as
  independent measurements.

### Additional categorical colors

Use colorblind-aware colors in this order when additional categories are
unavoidable: `#009E73`, `#CC79A7`, `#E69F00`, `#56B4E9`, `#F0E442`, and
`#000000`. Do not use red and green as the sole distinction.

## Marks and uncertainty

- Primary curves: approximately 0.6-0.9 mm `ggplot2` linewidth.
- Secondary curves: approximately 0.4-0.6 mm linewidth.
- Reference lines and grid lines: approximately 0.25-0.4 mm linewidth.
- Points: normally 1.8-2.6 mm, with a visible stroke when filled.
- Uncertainty bands: the same hue as the estimate with alpha 0.15-0.25 and no
  heavy outline.
- Show individual observations when sample size and density permit. Summary
  marks must not conceal the underlying distribution.
- State in the caption whether a band shows an interval, standard error,
  standard deviation, or an empirical quantile range.

## Axes, labels, and notation

- Label every quantitative axis with the measured quantity and unit.
- Use days for experimental time unless another unit is biologically clearer.
- Use consistent limits for direct comparisons and document intentional axis
  truncation.
- Avoid dual y axes.
- Render oxygen as O2 with the 2 subscripted when the plotting system supports
  it, and always include percent units when plotting oxygen concentration.
- Define specialized model quantities in the caption on first use. Do not rely
  on internal variable names such as `p_mis` as the only reader-facing label.
- Prefer biologically interpretable labels such as "Missegregation probability
  per chromosome" and place the code variable in Methods or supplementary
  provenance.

## Accessibility and visual QC

Before a panel is accepted:

1. Inspect the final PNG directly at intended manuscript size.
2. Confirm that no labels, legends, or panel marks overlap or are clipped.
3. Confirm that all text is readable and that the plot remains interpretable in
   grayscale.
4. Confirm that critical comparisons use shape, line type, or direct labels in
   addition to color where practical.
5. Confirm that the oxygen 20.5% condition remains visible on white.
6. Confirm that panel labels follow the figure caption and manuscript order.
7. Confirm that the displayed data scope matches the source manifest.

## Regeneration and raster exceptions

- Generate manuscript-facing plots and composites from repository-local scripts
  and direct data or report-export tables.
- Assemble final composites from plot or grob objects rather than screenshots
  or previously exported panel PNGs.
- Follow `.agents/references/immutable_raster_policy.md` for any asset that
  genuinely cannot be regenerated as a plot or vector object.
- A legacy plot PNG is not an immutable scientific source merely because its
  original plotting script is currently missing.
