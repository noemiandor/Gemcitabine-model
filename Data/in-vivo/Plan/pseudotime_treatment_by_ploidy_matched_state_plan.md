# Matched-Pseudotime Treatment-by-Ploidy Program Analysis: Implementation Plan

## V2 amendment: post-hoc methodological investigation frozen on 2026-07-15

Version: `04j_v2_2026_07_15`.

This amendment freezes a second workflow, `04j_pseudotime_treatment_ploidy_programs_v2`, as a post-hoc methodological investigation. The v1 result tree remains the original primary analysis and must be retained without overwrite, movement, deletion, or reinterpretation as a v2 rerun. V2 tests why program-level conclusions may differ across score aggregation, gene-set null hypotheses, pseudotime coordinates, and technical support assumptions.

The 04i and 04j analyses estimate different quantities. The 04i accumulated-state workflow compares expression in the frozen accumulated pseudotime state against adjacent states with a common smooth and no treatment-by-pseudotime interaction estimand. The 04j workflow estimates treatment-by-baseline-initial-ploidy differences within matched pseudotime support. Evidence from 04i may be used only as an exploratory bridge projection in 04j v2; it is not an independent validation of the 04j treatment interaction.

Tier remains an interpretation layer, not a method filter. `analysis_role` determines mathematical eligibility: `directional_signature` programs may receive directional scalar-score interpretation; `broad_pathway` programs receive self-contained, competitive, and ranked GSEA evidence without assuming all genes move in one direction; `legacy_union` programs reproduce v1-style score summaries for comparison; `contextual_control` and `matched_null_control` programs calibrate specificity and false-positive behavior. All estimable programs are carried through the frozen method set, and FDR is reported separately for all-universe, focused Tier 1+2, Tier 3, and controls within each statistical method family.

Scalar scores, self-contained tests, competitive tests, and GSEA have different null hypotheses. Scalar scores test a pre-specified weighted sample-level score. Fry and mroast test coordinated set-level expression shifts under self-contained gene-set nulls. Camera-style tests are competitive and depend on inter-gene correlation assumptions. Ranked GSEA is exploratory pathway ranking evidence and is not a coefficient for the 04j interaction.

V2 freezes the program universe, score definitions, simulation scenarios, random seeds, FDR families, method eligibility, coordinate-analysis contract, and final interpretation rules before inspecting any v2 interaction result. Current 04j cameraPR results, GSEA leading edges, real interaction P values, real interaction directions, or v2 score changes must not be used to select genes, set weights, choose programs, choose coordinates, tune simulation scenarios, or decide method eligibility.

The retained v1 forty-program universe is labeled `legacy_v1_comparison`. V2 additionally decomposes broad MSigDB unions into atomic exact MSigDB set records that preserve source set ID, collection, version, response family, tier, and source program. Broad pathways are not forced to be directional signatures. Directional signatures require external pre-specified signed weights.

The frozen score set is:

1. `signed_weighted_mean_zscore`
2. `signed_pca_eigengene`
3. `signed_rank_mean`
4. `trimmed_signed_mean_zscore`
5. `within_subbin_residual_zscore`

The first four scores are calibrated before introducing the fifth score. The fifth score uses TMM-normalized mouse-subbin logCPM, fits a subbin-only model without treatment, initial ploidy, endpoint ploidy, TGI, or interaction terms, centers by the subbin-only fitted mean, scales genes with robust empirical-Bayes posterior residual SD, applies frozen signed weights normalized by `sum(abs(weights))`, and records program-observation precision weights from robust summaries of program-gene voom weights. The formula is frozen before real v2 interaction results are interpreted.

The frozen gene-set method set is fry directional mean, mroast directional mean with set-specific signed weights, mroast `set.statistic = "msq"` with Mixed P/FDR, cameraPR with `use.ranks = FALSE` and `TRUE`, cameraPR at `rho = 0.01, 0.03, 0.05, 0.10`, residual empirical inter-gene correlation, and ranked GSEA as exploratory evidence. Observed-data mroast uses a fixed seed and at least 9999 rotations when computationally available.

The frozen baseline simulation keeps the observed mouse-subbin observation count, mice, group allocation, missing-subbin pattern, design matrix, mouse block, gene residual scale, program covariance, program overlap, and inter-gene correlation summaries. Scenarios include null interaction, small/medium/large programs, 10/25/50/100 percent active genes, coherent/mixed/externally signed effects, and standardized interaction effects below, near, and above the observed minimum detectable range. Null simulations use at least 1000 replicates; each power scenario uses at least 300 replicates. Simulation outputs must include type-I error, FDR, power, effect bias, CI coverage, direction accuracy, method failure rate, and Monte Carlo SE.

The independent-coordinate contract is frozen as two sensitivity coordinates: `untreated_learned_coordinate`, learned only from control cells and applied to treated cells by a pre-specified out-of-sample projection with fixed root/orientation, and `target_genes_excluded_coordinate`, relearned after excluding nucleotide, replication, stress, repair, and tested program genes. Coordinates are compared against original pseudotime by correlation, density, ECDF, phase/coverage diagnostics, and then all mathematically estimable frozen 04j models are rerun without changing programs, weights, FDR families, or interpretation rules.

The independent-signature contract is frozen as follows. If a unique auditable independent in-vitro gemcitabine gene-level response is found in the repository or input roots, it is versioned, capped, L1-normalized, and used only as `external_invitro_projection` without thresholds chosen from 04j results. Whether or not such a source exists, `cross_fitted_04i_state_projection` is run as an exploratory bridge: each held-out mouse receives weights refit from the remaining mice, without 04j interaction labels, within response families; out-of-fold scores are then modeled with the 04j interaction design.

Final v2 interpretation is governed by pre-frozen rules. Score aggregation loss requires stable fry/mroast or independent/cross-fitted evidence plus simulation evidence that existing scores lose power under the relevant sparsity or direction pattern. Global scaling loss requires calibrated within-subbin residual-score power and stable real-data direction, leave-one-out behavior, and support. Pseudotime over-conditioning requires stronger and robust evidence on untreated/excluded coordinates across methods, not a single score. Low power means plausible-effect simulation power is low, CI/MDE are wide, and methods are jointly unstable or null. Camera-only weak evidence is flagged when rho sensitivity and controls undermine support. Robust program evidence requires at least one direct mouse-level or self-contained method passing FDR with adequate support, directional consistency, leave-one-out stability, and coordinate robustness. Failure to reject must not be interpreted as ploidy independence.

## V3 amendment: post-hoc methodological repair frozen on 2026-07-15

Version: `04j_v3_2026_07_15`.

This amendment freezes a third workflow, `04j_pseudotime_treatment_ploidy_programs_v3`, as a post-hoc methodological analysis. V3 does not overwrite V1 or V2 code or results. The V3 result root is:

```text
/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs_v3
```

The primary estimand remains the treatment by baseline initial-ploidy interaction within matched pseudotime support. Continuous endpoint tumor ploidy, thresholded end-time ploidy, and dose-specific analyses remain secondary/sensitivity analyses and are written to separate V3 directories: `secondary_etp/`, `secondary_end_time_ploidy/`, and `secondary_dose/`. V3 may use 04i only as an exploratory cross-fitted projection with a different estimand; 04i significance is not a positive control for 04j.

Before V3 real-result interpretation, the workflow must export frozen program universe, response family, tier, analysis role, eligibility, directional eligibility, control status, score definitions, permutation space, simulation scenarios, FDR families, random seeds, QC criteria, and source checksums under `frozen_definitions/` and `manifest/`. V1/V2 code, plan, and result-root checksums are recorded before and after V3 and compared in `manifest/`.

V3 method roles are:

- `legacy_non_directional`: retained legacy program scores for historical score-model comparison; these do not create directional scalar biological claims.
- `broad_pathway`: atomic MSigDB pathway records analyzed with self-contained, competitive, exact permutation, multidimensional, and trajectory-wide methods without assuming a single coherent direction.
- `mechanistic_directional_signature`: pre-frozen signed mechanistic definitions, currently limited to the custom gemcitabine handling component definition unless an independent external source is found.
- `matched_null_control`: negative/context controls analyzed in the same FDR families as biological programs.

Every eligible program must carry `program_id`, `response_family`, `tier`, `analysis_role`, `eligibility`, `directional_eligibility`, and `control_status` through score models, fry, mroast, camera, exact permutation, multidimensional, trajectory, and report tables. Any eligible program missing tier or response family is an acceptance-test failure. For each method family, V3 reports all-universe FDR, focused Tier 1 + Tier 2 FDR, Tier 3 FDR, and negative/control FDR from the same result objects used in the HTML report.

V3 replaces V2 simulation with expression-scale injection. It does not use artificial `method_scale`. Simulated effects are injected into the gene-expression matrix on the common expression scale, after preserving observed design, mouse/subbin structure, missingness pattern, residual scale, program size/overlap, and observed gene covariance as far as the local pseudobulk representation permits. The frozen scenarios include 1000 null replicates and 300 replicates per power scenario across small/medium/large program sizes, coherent/mixed/sparse direction patterns, active-gene fractions, and effect sizes. Simulation summaries include type-I error, empirical rejection/FDR proxies, power, bias, interval coverage, direction accuracy, failure rate, and Monte Carlo SE. Vectorized expansion is not used unless validated against end-to-end expression-scale simulation.

V3 self-contained gene-set tests must use the same design matrix as the primary model and pass `block = sample_id` plus the second-pass duplicateCorrelation coefficient to `fry` and `mroast`. Mroast uses a frozen rotation seed and 9999 rotations when feasible. Mixed-direction/broad pathways are interpreted primarily through Mixed or msq-style evidence rather than directional mean claims.

V3 implements true leave-one-mouse-out 04i projection. For each held-out mouse, the 04i gene-state model is refit after removing that mouse from the original 04i input cells. Fold-specific weights are generated only from training mice and then projected to held-out 04j mouse-subbin observations. Fold audit tables must record held-out mouse, training mouse IDs, training input checksum, frozen training gene universe checksum, fitted weight checksum, nonzero genes, projection checksum, and explicit failure/fallback status. Full-data 04i t-statistics cannot be reused as cross-fitted weights.

V3 directional-signature rules are conservative. Broad pathways and legacy unions are `legacy_non_directional` unless they have an independently auditable external in-vitro gemcitabine response signature or a frozen mechanistic signed definition. If no unique external in-vitro source is found, the independent compact directional score is `not_estimable` for real data; method code, search audit, and status are still exported. No 04j GSEA leading edge, gene-level interaction direction, P value, or V3 intermediate result may define a directional weight.

V3 adds six method families on the unified frozen universe:

1. Exact mouse-level ploidy-stratified permutation adaptive gene-set test. It enumerates all `choose(8,4)^2 = 4900` assignments, keeps whole mouse blocks and ploidy strata fixed, computes directional signed mean only for legal directional signatures, and computes maxmean, msq, maxabs/sparse, calibrated adaptive omnibus, BH FDR, and Westfall-Young maxT sensitivity.
2. Independent compact directional score. This runs only when a unique external in-vitro gemcitabine response source is found and audited; otherwise real-data status is `not_estimable`.
3. Control-reference residual score. Control samples define subbin/ploidy baselines and shrinkage residual SDs, and treatment observations are projected relative to those references with explicit reference-training audit.
4. Shrinkage covariance-whitened score. Legal signed weights are combined with ridge/shrinkage covariance estimated from allowed training/control observations; condition number, effective rank, inversion status, dominant genes, and visible fallback status are exported.
5. Multidimensional pathway interaction. Broad pathways are projected onto a PCA basis learned from control/training observations; PCs explaining 80% variance up to 3 PCs are tested jointly with the same 4900 mouse-level permutation space.
6. Trajectory-wide global interaction. A 20-bin pseudobulk spline model tests the global null that all treatment by initial-ploidy by spline coefficients are zero; per-bin estimates are localization/descriptive only.

V3 output directories are organized as:

```text
manifest/
frozen_definitions/
qc/
primary_initial_ploidy/score_models/
primary_initial_ploidy/exact_permutation/
primary_initial_ploidy/fry_mroast/
primary_initial_ploidy/camera/
primary_initial_ploidy/multidimensional/
primary_initial_ploidy/trajectory_global/
secondary_etp/
secondary_end_time_ploidy/
secondary_dose/
crossfit_04i/
directional_signatures/
simulation/
controls/
figures/
tables/
report/
logs/
checkpoints/
```

The V3 HTML report is self-contained with embedded images and includes scientific question/estimands, frozen definitions, cohort/design QC, program universe, simulation before real-result interpretation, tiered results, controls, score comparison, exact permutation, fry/mroast/camera/GSEA, cross-fitted 04i/external signatures, multidimensional results, trajectory-wide global interactions, ETP/end-time ploidy/dose sensitivities, cross-method concordance, integrated conclusion, limitations, not-estimable methods, and audit information. Nominal P values, exploratory GSEA, unstable single-method hits, or sensitivity-only results must not be written as confirmatory findings.

## V4 amendment: post-hoc methodological repair frozen on 2026-07-16

Version: `04j_v4_2026_07_16`.

This amendment freezes a fourth, independent workflow,
`04j_pseudotime_treatment_ploidy_programs_v4`, before inspecting any V4
P value, effect direction, adaptive component, GSEA leading edge, best window,
or coordinate-sensitivity result. V4 is a post-hoc methodological repair
motivated by review of V3 implementation and calibration. It does not overwrite
or reinterpret V1, V2, or V3 as a rerun. The V4 result root is:

```text
/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs_v4
```

The primary estimand remains:

```text
(gemcitabine - control in baseline-4N-origin tumors)
-
(gemcitabine - control in baseline-2N-origin tumors)
```

within matched common pseudotime support, using `initial_ploidy` as the primary
effect modifier. Marginal total transcriptional interaction, multiscale
pseudotime interaction, dose-specific interaction, continuous endpoint tumor
ploidy, threshold endpoint tumor ploidy, alternative coordinates, and
pseudotime-region occupancy interaction are separate secondary or sensitivity
families and are never mixed with primary FDR. The 04i state analysis is an
internal cross-fitted exploratory bridge with a different estimand; it is not an
independent validation or positive control for 04j.

No external expression cohort or external signature is used in V4. The frozen
external branch status is:

```text
disabled_by_user_no_external_cohort
```

V4 does not search GEO, use a network connection, download data, or use any
external expression queue. The supplied inputs are already restricted to the
intended cells and genes. V4 adds no species filter, human/mouse split, species
ratio QC, or species-based analysis.

### V4 preservation, checksum, and checkpoint contract

Before plan or code changes, V4 records SHA256 checksums for V1, V2, and V3 code
and result roots and the pre-amendment plan. At completion it recomputes them.
V1/V2/V3 code and result files must be byte-identical. The plan is expected to
change and is compared separately. Current CellCycle metadata, NonCellCycle
metadata, sample information, and Seurat RDS checksums are recorded and compared
with the V3 manifest wherever V3 recorded a checksum; any V3 manifest omission
is reported as not verifiable rather than silently accepted.

All intensive steps use versioned checkpoints containing the V4 config SHA256,
method version, input checksum, schema version, and completion status. A
checkpoint is reusable only when all fields match. `--overwrite=TRUE` replaces
V4-derived artifacts but preserves the original V4 preflight checksum files.
BLAS/OpenMP threads are fixed at one and the workflow uses at most eight
process-level workers.

### Mouse-level standardized-expression operator

The frozen primary common-support subbins remain those selected by the original
support-only grid. Each eligible subbin has base standardization weight
`1 / number_of_eligible_subbins`, matching the primary estimand. Before subbin
aggregation, each gene's TMM logCPM is residualized for a subbin-only nuisance
mean that uses neither treatment labels nor expression-effect selection.

Each mouse is then represented once. A mouse with all eligible subbins receives
the base weights. If an eligible subbin is absent because it fails the original
minimum-cell rule, the available base weights are renormalized to sum to one.
This missing-subbin rule depends only on the frozen support indicator and not on
expression values, treatment labels, P values, or effects. V4 exports every
mouse-by-subbin base weight, availability indicator, final weight, and the
complete 16-mouse standardization matrix.

For each retained gene, V4 computes the studentized difference-in-differences
and estimated standard error from the 16 mouse vectors. It enumerates all
`choose(8,4)^2 = 4900` ploidy-stratified whole-mouse treatment assignments.
Every assignment recomputes group means, within-group variances, the interaction
estimate, and studentization. Full gene-level permutation statistics are stored
in a chunked HDF5 artifact with a row/assignment index and checksum audit.

The reduced-null Freedman-Lane sensitivity fits a nuisance model containing
initial ploidy and a common treatment effect but no interaction. Whole-mouse
residual blocks are reassigned within the frozen ploidy-stratified assignment
space, after which the interaction and studentization are recomputed. Its
purpose is to prevent a large common treatment effect from masquerading as a
ploidy-specific interaction.

### Exact-membership adaptive gene-set tests

The exact statistical unit is the SHA256 of the sorted unique observed genes.
The 106 historical labels remain for interpretation, but duplicate labels map
to one exact membership through a frozen alias table. V4 exports exact duplicate
groups, pairwise Jaccard overlap, gene multiplicity, and family membership.

For every exact membership, the component statistics are frozen as:

1. a two-sided signed weighted mean only for a legally directional,
   independently frozen mechanistic signature;
2. true maxmean,
   `max(mean(max(z,0)), mean(max(-z,0)))`, rather than `abs(mean(z))`;
3. mean-square statistic;
4. maximum absolute statistic;
5. aSPU with gamma `1`, `2`, `4`, `8`, and `infinity`;
6. higher-criticism sparse statistic;
7. a ridge-shrinkage covariance-adjusted quadratic statistic using the
   treatment-label-independent 16-mouse covariance representation.

All components use the same 4900 assignments. For each assignment, every
component is converted to an empirical P value against the same complete null
distribution. The within-assignment minimum component P is then calibrated
against its complete permutation distribution. The reported adaptive P value is
therefore nested and selection-adjusted. Westfall-Young maxT sensitivity is
computed over unique memberships. An uncalibrated selected-component P value is
never reported as the omnibus result.

Influence outputs include dominant-gene contributions, leave-one-mouse-out
interactions and sign stability, full leave-one-gene observed-statistic changes
for small and medium sets, all gene influence values plus ranked top genes for
large sets, and permutation-calibrated leave-one-gene robustness for any set
that reaches its applicable FDR threshold. The previously observed PPP4R2/HDR
pattern is a predeclared influence audit only and cannot alter membership.

### Response-family hierarchical inference

For each response family, V4 computes a deduplicated family gene-union adaptive
test, a member-membership minP test, and a nested permutation-calibrated family
omnibus combining those two evidence routes. Family FDR is separated into the
11 focused Tier 1 + Tier 2 families, the three Tier 3 families, and controls.
Transparent all-family and all-membership FDR is also reported.

A focused program-level claim requires both family FDR below 0.10 and
within-family Westfall-Young adjusted P below 0.10. Exact aliases cannot be
counted as independent corroborating results.

### Repaired score and gene-set methods

The five frozen score methods remain:

1. `signed_weighted_mean_zscore`;
2. `signed_pca_eigengene`;
3. `signed_rank_mean`;
4. `trimmed_signed_mean_zscore`;
5. `within_subbin_residual_zscore`.

Fry and broad or mixed mroast use `PValue.Mixed` for multiplicity correction.
Mroast mean may use a directional P only for a legally directional mechanistic
signature; mroast msq always uses `PValue.Mixed`. Camera repeats the original
rank and rho grid and adds each set's residual empirical correlation. Camera is
fixed as `diagnostic_only` if the negative-control calibration rule fails.

The control-reference score uses controls from the same subbin and initial
ploidy. Treated mice use all eligible controls; a control observation uses a
leave-one-control-mouse-out reference. The reference mean and shrunk SD are
recomputed for every permuted assignment, with no self-inclusion.

The covariance-whitened score runs only for legally directional programs and
exports the real program ID, metadata, estimate, SE, CI, condition number,
effective rank, dominant genes, ridge value, and explicit fallback. Every output
row must pass metadata acceptance.

The multidimensional method uses a treatment-label-independent PCA basis learned
from all 16 standardized mouse vectors after treatment-blind centering. This
basis is frozen for all assignments, which avoids the V3 observed-control-basis
leakage. V4 exports PC variance, loadings, basis checksum, training audit, and
permutation-basis audit.

The trajectory global test and all adaptive tests use exact whole-mouse
permutation P values. Limma asymptotic statistics, when exported, are
supplementary diagnostics only.

### Corrected cross-fitted 04i bridge

Each of 16 folds deletes the held-out mouse from the original 04i cell input
before pseudobulk construction, normalization, gene filtering, centering,
scaling, state-model fitting, and weight estimation. The original 04i gene
universe is available to training and is not restricted to 04j program genes.
Held-out observations use only fold-training filtering and standardization
parameters. Full-data 04i statistics and full-data 04j gene z-scores are
forbidden. Fold audits record training cells and mice, training gene universe,
normalization/standardization checksums, weight checksum, projection checksum,
and explicit held-out exclusion.

### Corrected simulation calibration

The simulation reduced model contains subbin, initial ploidy, a common treatment
effect, and frozen nuisance terms, but no treatment-by-initial-ploidy
interaction. Its fitted mean is therefore a true interaction-null mean.
Residuals are resampled by independent whole-mouse wild-bootstrap sign flips,
preserving within-mouse subbin dependence, missing-subbin patterns, gene
covariance, and residual scale.

Each power replicate injects signal into one target program only. Target
programs and active genes are sampled with frozen seeds rather than membership
order. Program size, active-gene fraction, and coherent, mixed, or sparse
patterns alter the injected gene signal. The five score-specific true effects
are recomputed on the noiseless injected expression using each score's actual
definition; a gene-scale beta is never substituted for a score-scale truth.

The null contains at least 1000 replicates and each power scenario at least 300.
Replicate-level FDP and FDR are computed after applying the frozen q=0.10 family
within that replicate; correlated program tests are not counted as independent
Monte Carlo replicates. Outputs include type-I error, empirical FDR, power,
bias, 95% CI coverage, direction accuracy, failure rate, Monte Carlo SE, and the
score-specific true effect.

Calibration is classified without tuning:

- the binomial Monte Carlo interval for type-I error must cover 0.05, with
  0.035-0.065 also shown as a reference interval;
- 95% CI coverage is flagged in range only at 0.925-0.975;
- global-null empirical FDR at q=0.10 is flagged in range only when at most
  0.12;
- a failing score or method is `uncalibrated` and cannot support a positive
  claim.

### Multiscale, trajectory, marginal, and dose analyses

The multiscale scan uses the frozen 20-bin grid and all contiguous windows of
width 2, 3, 4, 5, and 6 bins. Window eligibility depends only on the original
minimum-cell rule and at least two contributing mice per ploidy-by-treatment
group. Every assignment repeats the complete scan. Inference uses the
across-window maximum statistic; the best window is localization only.

The trajectory-wide analysis starts at 20 bins and applies a deterministic
support-only adjacent-bin merging algorithm until all four groups contain at
least two mice per analysis bin. Spline degrees of freedom 2 and 3 are frozen.
If the stronger df is selected, selection occurs inside every permutation. V4
uses all 4900 assignments, exports simultaneous maxT bands, leaves unsupported
regions blank, and reports the `[0.30,0.49]` integrated interaction using the
same standardization reference as the primary estimand.

The marginal analysis sums raw counts over all eligible target cells within each
mouse, performs TMM normalization, and fits
`expression ~ initial_ploidy * treatment`. It runs the five scores, corrected
fry/mroast, exact adaptive membership tests, and family hierarchy in
`secondary_marginal_total_effect/`. Frozen-region mouse occupancy interactions
are reported beside marginal and matched-state evidence without causal
mediation language.

The dose analysis preserves 4 control, 2 dose-30, and 2 dose-120 mice per
initial-ploidy stratum and enumerates the complete
`[8!/(4!2!2!)]^2 = 176400` assignment space. The frozen candidate statistics
are pooled any-dose, linear trend, 30-versus-control,
120-versus-control, 120-versus-30 nonlinearity, and a joint dose-by-ploidy
statistic. Candidate selection is nested and exactly permutation-calibrated.
Dose FDR remains secondary and separate.

### Alternative-coordinate sensitivity

`original_coordinate` is the supplied frozen pseudotime. The alternative
coordinates reuse the repository's Monocle3 trajectory framework with root
cluster 6 and fixed orientation.

`untreated_learned_coordinate` learns the trajectory only from control target
cells. Treated cells are projected out of sample through a training-only PCA
representation and nearest-neighbor interpolation of the control Monocle3
pseudotime; treatment response is not used to choose root or orientation.

`target_genes_excluded_coordinate` excludes the complete frozen 04j program-gene
universe before trajectory learning. The exclusion list is exported before any
V4 result is examined. Root cluster 6, orientation, preprocessing, projection
rules, and seeds are frozen.

Each coordinate exports training cells/mice, excluded genes, root/orientation,
projection parameters and checksum, treatment-label-use audit, correlation,
density, ECDF, Phase/S.Score/G2M.Score diagnostics when available, and support
changes. Each estimable coordinate reruns the five score models, exact adaptive
tests, family hierarchy, multiscale scan, trajectory global test, and corrected
fry/mroast. A failed coordinate is reported as `not_estimable` with the exact
attempt and error rather than omitted.

### V4 outputs, report, and interpretation

The V4 tree is:

```text
00_manifest/
frozen_definitions/
qc/
simulation/
primary_initial_ploidy/
  score_models/
  fry_mroast/
  camera_diagnostic/
  mouse_level_gene_statistics/
  adaptive_gene_set/
  response_family_hierarchical/
  multiscale_pseudotime/
  trajectory_global/
  robustness/
secondary_marginal_total_effect/
secondary_dose/
secondary_etp/
secondary_end_time_ploidy/
coordinate_sensitivity/
crossfit_04i/
figures/
tables/
report/
logs/
checkpoints/
```

The self-contained HTML report embeds every PNG as base64 and contains 21
ordered sections: scientific question; frozen/post-hoc definitions; QC;
simulation; Tier 1; Tier 2; Tier 3; controls; score comparison; adaptive tests;
family hierarchy; multiscale scan; trajectory-wide interaction; marginal versus
matched; dose/ETP/end-time; coordinate sensitivity; 04i cross-fit; influence;
cross-method concordance; integrated conclusion; and limitations/final audit.
Each figure has a continuous number, formal title, stand-alone legend, and
interpretation in its result section. A forest plot must show the two simple
treatment effects and interaction with confidence intervals; estimate-versus-P
scatter plots and files named `forest_proxy` are forbidden.

A focused confirmatory statement requires family FDR below 0.10,
within-family Westfall-Young P below 0.10, a calibrated method, acceptable
negative controls, at least one mouse-level exact/adaptive or self-contained
supporting method, no complete one-mouse or one-gene domination, adequate common
support, and no clearly opposite alternative-coordinate result. Nominal,
camera-only, best-window-only, single-score, Tier-3-only, sensitivity-only, or
uncalibrated evidence is exploratory.

If no result qualifies, V4 must state that failure to reject is not evidence of
ploidy independence, identify low-power, local, sparse, or unstable patterns,
and report observed minimum detectable effects plus a mouse-number planning
range. Scientific calibration failure may coexist with successful workflow
completion; code errors, leakage, missing schemas, or corrupt metadata may not
be relabeled as scientific negatives.

## V5 amendment: five non-score gene-vector analyses frozen on 2026-07-16

Version: `04j_v5_non_score_2026_07_16`.

This amendment freezes a fifth standalone workflow,
`04j_pseudotime_treatment_ploidy_programs_v5_non_score`, before inspecting any
V5 P value, enrichment direction, leading edge, posterior estimate, state
decomposition, spline coefficient, local trajectory feature, or closed-testing
result. V5 is a post-hoc methodological extension motivated by the absence of a
qualified focused-family result in V1-V4 and by the need to answer the same
scientific question without reducing a gene set to a scalar program score.

V1-V4 code and result trees remain read-only. Their source and result-tree
checksums are captured before the V5 plan amendment and compared again at
completion. The V5 result root is:

```text
/Volumes/Protable Disk/Project/BreastCancerOrthotopicModels/Results/04j_pseudotime_treatment_ploidy_programs_v5_non_score
```

The primary estimand remains:

```text
(gemcitabine - control in baseline-4N-origin tumors)
-
(gemcitabine - control in baseline-2N-origin tumors)
```

within the frozen `[0.30,0.49]` matched-pseudotime support. The biological unit
is the mouse. The exact primary randomization space remains the complete
ploidy-stratified whole-mouse space,
`choose(8,4)^2 = 4900`, including the observed assignment.

V5 uses the V4 standardized 16-mouse gene-expression matrix, complete
gene-by-assignment HDF5 statistics, assignment registry, exact memberships,
response-family definitions, trajectory support, marginal pseudobulk
checkpoint, cell metadata, and Seurat annotations as frozen inputs. These
objects are not recomputed unless a checksum or schema validation fails. V5
adds no species filter or species split and uses no external cohort, external
signature, network query, treatment-derived state clustering, or gene selected
from V1-V4 P values, effects, leading edges, dominant-gene results, or
scientific conclusions.

Program scores are forbidden as a V5 inferential unit. Frozen gene-set and
response-family membership may define a gene vector, enrichment set,
mechanistic hierarchy, or multiplicity family, but V5 does not calculate,
model, test, or use any scalar average-expression program score to support a
scientific conclusion.

### V5 method 1: exact whole-mouse ranked enrichment

For each of the 4900 assignments, genes are ranked by the V4 studentized
gene-level treatment-by-initial-ploidy statistic. The inferential statistic is
a weighted Kolmogorov-Smirnov enrichment statistic computed separately for
positive and negative tails and combined as
`max(abs(ES_positive), abs(ES_negative))`. Tail selection is repeated inside
every assignment. The weight exponent is fixed at one; an unweighted exponent
zero result is a labeled sensitivity and is not selected after observing data.

All 86 exact memberships and all 15 deduplicated response-family unions are
tested. The exact empirical P value uses:

```text
(1 + number of null statistics at least as extreme as observed)
/
(1 + number of assignments)
```

BH FDR is reported for all memberships, focused Tier 1+2 memberships, Tier 3,
controls, all families, focused families, Tier 3 families, and controls.
Westfall-Young maxT is reported over exact memberships. Normalized enrichment
scores are display quantities only; exact inference uses the raw enrichment
statistic. Leading-edge genes are exported only after inference and cannot
define a new set, weight, mechanism, or secondary test.

### V5 method 2: matched-state versus marginal-state decomposition

This branch separates within-state expression differences from state
composition differences without claiming causal mediation. State definitions
are frozen from pre-existing treatment-independent annotations:

1. the two V4 primary matched-pseudotime subbins;
2. `cluster_final` levels `6`, `10`, and `4c`;
3. Seurat Phase levels `S` and `G2M`.

G1 is exported descriptively because the target-cell subset contains too few G1
cells for a stable four-group interaction. `cluster_final × Phase` is attempted
only as a predeclared sensitivity. A state is inferentially eligible only when
every ploidy-by-treatment group contains at least two contributing mice and
each contributing mouse has at least five cells. Unsupported states are
reported as `not_estimable`; unrelated clusters or phases are not merged to
obtain significance.

For differential abundance, each mouse contributes state counts and its total
captured target-cell count. A fixed 0.5 count replacement is used before the
centered-log-ratio transformation. The global statistic is the squared norm of
the state-level treatment-by-initial-ploidy contrast, and state-level
statistics are its studentized components. The same 4900 whole-mouse
assignments provide exact global and state-level P values and Westfall-Young
adjustment.

For expression decomposition, raw counts are summed per mouse-by-state,
TMM-normalized, and modeled gene by gene. For every gene, the marginal
interaction is decomposed using a fixed pooled reference state distribution
into:

1. within-state expression component;
2. state-composition component;
3. residual/non-additive component.

Response-family evidence is calculated from the complete gene-component vector
using exact maxmean, mean-square, and sparse maximum statistics. No component
is averaged into a scalar program expression score. All state definitions,
reference weights, support rules, and vector statistics are recomputed or
reapplied inside each exact assignment as required by their null.

### V5 method 3: Bayesian hierarchical gene-vector model

The Bayesian branch uses the frozen V4 mouse-standardized gene-expression
matrix. Each response family is fit separately using gene-level standardized
outcomes and the mouse design:

```text
expression_gene_mouse =
  gene intercept
  + gene ploidy effect
  + gene common treatment effect
  + gene treatment-by-ploidy interaction
  + treatment-blind residual factors
  + residual error
```

The gene-specific interaction is represented as a family-level coherent mean
plus regularized sparse gene deviations. A robust Student-t prior protects the
family mean from a small number of extreme genes, and a regularized horseshoe
prior permits sparse deviations without selecting genes from observed V5
effects. Up to three treatment-blind residual factors are learned from the
reduced-null 16-mouse expression matrix and are frozen before sampling.

Gene outcomes are scaled by reduced-null residual SD. The
practical-equivalence region is frozen as `[-0.20,0.20]` residual SD units. The
report includes the family posterior mean and 95% interval, posterior direction
probabilities, `P(abs(family mean) > 0.20)`, posterior coherent-gene fraction,
shrunken gene-level interactions, prior and posterior predictive checks,
leave-one-mouse-out sensitivity, and prior sensitivity.

The local RStan implementation uses four chains for the full analysis. A fit is
accepted only with `Rhat < 1.01`, bulk and tail effective sample sizes above
400 for the family-level quantities, zero divergent transitions, and passing
posterior predictive checks. Bayesian evidence is model-based supportive
evidence and cannot create a confirmatory claim without aligned evidence from
an exact whole-mouse method.

### V5 method 4: frozen mechanistic exact global/kernel test with closure

Five mechanistic families are frozen before V5 results:

1. gemcitabine handling;
2. nucleotide supply and buffering;
3. replication execution;
4. replication-stress checkpoint;
5. fork protection and repair.

Their genes and child subpanels come only from the frozen V4 role, membership,
and response-family tables. The signed burden statistic is eligible only for
the independently directional gemcitabine-handling definition. Every family
also receives an unsigned ridge-shrinkage kernel statistic.
Burden-versus-kernel adaptation is selected and calibrated inside every one of
the 4900 assignments.

Inference follows a closed hierarchy:

```text
mechanistic family
  -> frozen exact membership/subpanel
    -> individual genes
```

A child inferential statement is made only when every intersection hypothesis
on its path is rejected. Family and subpanel levels use exact permutation
P values, nested minP calibration, and Westfall-Young adjustment. Gene-level
tests use the same exact assignment statistics and closed/Holm adjustment
within an opened subpanel. Deleting the largest-contribution gene is a frozen
fragility audit; a family result that loses direction or qualification is
labeled `single_gene_sensitive`.

### V5 method 5: mouse-level gene-wise pseudotime spline interaction

This branch uses the V4 treatment-blind 20-bin grid, deterministic support-only
adjacent-bin merging, and frozen mouse support. For each gene, mouse-bin
pseudobulk expression is modeled as:

```text
expression ~ initial_ploidy * treatment * natural_spline(pseudotime)
```

with mouse blocking and observation precision weights. Candidate spline degrees
of freedom are two and three, restricted to full-rank candidates. Candidate
selection occurs inside every assignment. The gene-level global statistic is
the joint quadratic statistic for all treatment-by-initial-ploidy-by-spline
coefficients. Exact P values use the same 4900 whole-mouse assignments.

Gene-level global statistics are aggregated to exact memberships and response
families with predeclared maxmean, mean-square, and higher-criticism gene-vector
statistics, with the adaptive choice calibrated inside every assignment.
Simultaneous trajectory bands are based on permutation maxT. Pointwise curves
and local maxima are localization only.

A spline requires at least three distinct four-group-supported pseudotime bins
and a full-rank design. If `[0.30,0.49]` has fewer than three supported bins,
the primary-window spline is reported as `not_estimable`; the support threshold
is not weakened and unsupported regions are not extrapolated. The
whole-trajectory `[0,1]` support-only spline is a separate secondary estimand.
The `[0.30,0.49]` integrated curve is reported only when the fitted support
contains that interval without unsupported extrapolation.

### V5 simulation, multiplicity, and interpretation

Simulation is conducted before real-result interpretation. The observed
16-mouse design, assignment space, missing-state pattern, trajectory support,
gene residual scales, and treatment-blind covariance structure are preserved.
Frozen scenarios include global null, coherent family shift, sparse shift,
mixed-direction shift, composition-only shift, within-state-only shift, local
trajectory peak, and early/late reversal. Exact-method null simulations use at
least 1000 replicates for representative frozen families; power scenarios use
at least 300 replicates. Bayesian simulation uses prior predictive checks,
simulation-based calibration, interval coverage, direction accuracy, and
false-positive practical-equivalence summaries.

The exact frequentist reference interval for type-I error is `0.035-0.065`;
95% interval coverage is flagged in range only at `0.925-0.975`. A method that
fails its relevant calibration is labeled `uncalibrated` and cannot support a
positive biological conclusion.

Multiplicity is separated by estimand:

1. primary matched-state gene-vector evidence;
2. marginal/composition decomposition evidence;
3. whole-trajectory evidence.

P values are not pooled across these different estimands. Within the primary
matched-state estimand, method-specific FDR and exact familywise adjustments
are reported transparently. Bayesian posterior probabilities are never treated
as frequentist P values. A focused matched-state claim requires at least one
calibrated exact whole-mouse method passing the focused family threshold,
within-family adjustment below 0.10, aligned non-score evidence from a second
method, acceptable negative controls, adequate common support, and no complete
one-mouse or one-gene domination.

A V5 branch is complete when it has an inferential result or an explicit,
audited `not_estimable` status. Missing support, rank failure, or sampling
failure cannot be hidden, repaired by lowering thresholds, or relabeled as a
biological negative. Failure to reject remains inconclusive and is not evidence
of ploidy independence.

### V5 output and report contract

The V5 tree is:

```text
00_manifest/
frozen_definitions/
qc/
01_exact_ranked_enrichment/
02_state_decomposition/
03_bayesian_gene_vector/
04_exact_global_kernel_closed_testing/
05_gene_wise_trajectory/
06_simulation_calibration/
07_cross_method_synthesis/
figures/
tables/
report/
logs/
checkpoints/
```

The technical report is a self-contained HTML artifact with a technical
summary, key findings with visual evidence, data and estimand definitions,
separate sections for all five methods, calibration and robustness, integrated
cross-method interpretation, limitations, recommended next steps, and final
acceptance audit. It distinguishes the matched-state, marginal/composition,
whole-trajectory, and Bayesian model-based estimands. It does not use
vote-counting, nominal P values, leading-edge selection, or Bayesian posterior
probability alone to create a positive conclusion.

## Objective

Test whether gemcitabine is associated with different nucleotide-metabolism, DNA-replication, and replication-stress transcriptional responses in captured cells from baseline-2N-origin versus baseline-4N-origin tumors after standardizing the groups to the same inferred cell-cycle pseudotime distribution.

This is a treatment-by-initial-ploidy analysis. It is separate from the existing [accumulation-state annotation plan](pseudotime_accumulation_state_pathway_analysis_plan.md), whose estimand describes the expression state at pseudotime `0.30-0.49` relative to neighboring states and is not a treated-versus-control expression comparison.

The primary question here is:

> After standardizing pseudotime, is the gemcitabine-associated expression change in a nucleotide/replication-stress program different between cells from baseline-4N-origin and baseline-2N-origin tumors?

The sequenced cohort is outcome-selected: the manuscript states that it contains all vehicle tumors and representative high-dose tumors exhibiting treatment failure. Therefore, this analysis estimates a treatment-by-initial-ploidy expression association within the selected sequenced tumor subset and captured endpoint cells. It is hypothesis-generating and does not estimate the population-average treatment effect from the randomized in-vivo experiment.

Both pseudotime conditioning and endpoint-cell capture add further selection. If treatment or tumor origin changes entry into a state, death within that state, or capture probability, state conditioning can attenuate, create, or reverse an apparent expression interaction. The analysis therefore cannot establish a cell-autonomous effect of current ploidy, distinguish reduced physical damage from reduced signaling at equal damage, prove nucleotide buffering, measure transition speed, or describe cells that died before collection.

## Primary estimand

For every gene and pre-specified program, estimate within the frozen accumulated-state interval `[0.30, 0.49]` after all four groups have been standardized to one fixed pseudotime-bin distribution:

```text
delta_2N_origin = mean(2N-origin, gemcitabine) - mean(2N-origin, control)
delta_4N_origin = mean(4N-origin, gemcitabine) - mean(4N-origin, control)
interaction = delta_4N_origin - delta_2N_origin
```

The interaction is the primary estimand. Both simple treatment effects and their confidence intervals must accompany it. A ploidy main effect or a significant treatment effect in only one ploidy group is not evidence of interaction.

Use baseline `initial_ploidy` as the primary effect modifier. It is a tumor-of-origin assignment, not necessarily the ploidy of an endpoint cell. Endpoint tumor ploidy and single-cell ploidy are post-treatment variables and may be examined only in explicitly labeled descriptive sensitivity analyses. Results must be described as associations in cells from baseline-2N-origin or baseline-4N-origin tumors, not as cell-autonomous responses of current 2N or 4N cells.

## Treatment and dose definitions

The primary treatment variable is binary:

```text
control = 0 mg/kg
gemcitabine = 30 or 120 mg/kg
```

This gives four control and four treated biological samples within each initial-ploidy group. The resulting treatment contrast is the equally represented average of the observed 30- and 120-mg/kg arms in this selected cohort, not a generic gemcitabine effect. Pooling the two active doses is the primary analysis because it provides the best-supported treatment-by-ploidy comparison.

Secondary analyses will treat dose as a factor with levels `0`, `30`, and `120 mg/kg` and estimate separate 30-versus-control and 120-versus-control interactions. These are exploratory because each active dose contains only two mice per initial-ploidy group. A continuous-dose trend may be reported only as a sensitivity analysis because it imposes a linear dose-response assumption.

## Frozen inputs and preflight audit

Use the same raw RNA counts, cells, pseudotime values, root orientation, and cell inclusion rules used by the 04i analysis. The cell-level metadata source is:

```text
Data/in-vivo/CellCycleCells_pseudotime_distribution_per_sample_cell_level_with_ploidy_dose_tgi.csv
```

The current metadata contain 2,881 cells from 16 mice in a balanced mouse-level ploidy-by-dose design:

| Initial ploidy | 0 mg/kg | 30 mg/kg | 120 mg/kg |
|---|---:|---:|---:|
| 2N | 4 | 2 | 2 |
| 4N | 4 | 2 | 2 |

Before expression modeling:

1. Confirm exact one-to-one cell-ID matching between metadata and the raw-count matrix.
2. Require one dose and one initial-ploidy assignment per mouse.
3. Freeze pseudotime without rerooting, rescaling, or choosing new intervals after examining expression results.
4. Export mouse-by-region cell counts, library sizes, mean/median pseudotime, pseudotime quantiles, and cluster composition.
5. Audit animal identity, harvest date, tumor size/response at harvest, dissociation batch, 10x channel/library, sequencing batch, sequencing saturation, and demultiplexing scheme by ploidy and dose. Model-matrix rank cannot detect hidden or perfectly confounded processing batches. If a technical batch uniquely identifies a ploidy-by-dose group, that interaction is not biologically separable and modeling must stop with a non-estimability report.
6. Export viability/capture-quality metrics, mitochondrial fraction, UMI depth, detected genes, and any apoptotic-cell exclusion by group when available.
7. Confirm model-matrix rank and treatment/pseudotime overlap before fitting gene-level models.
8. Record input paths, checksums, code revision, package versions, and all exclusion reasons.

## Pseudotime regions

The primary expression comparison uses mouse-by-sub-bin pseudobulks from the frozen accumulated-state interval. At the whole-window level, all 16 mice currently contribute at least 10 cells; sub-bin support and library quality require the additional rules below.

| Role | Pseudotime interval | Use |
|---|---:|---|
| Primary matched state | `[0.30, 0.49]` | Primary hypothesis-generating treatment-by-ploidy comparison |
| Left state | `[0.11, 0.30)` | Localization and robustness |
| Right state | `[0.49, 0.68]` | Localization and robustness |
| Broad accumulated state | `[0.30, 0.61]` | Interval sensitivity |
| Late state | `[0.80, 0.98]` | Exploratory; coverage is weaker |

Do not use the narrow `[0.414, 0.426]` density-support interval for the primary expression analysis. It is too narrow for reliable mouse-level pseudobulk. Its role remains a density-localization result.

## Primary pseudotime-standardized pseudobulk model

A single pseudobulk from `[0.30, 0.49]` plus adjustment for its mean pseudotime is only window restriction with approximate matching: samples can have the same mean but different distributions across a nonlinear cell-cycle interval. The primary analysis will instead standardize all groups to the same fixed distribution of pseudotime sub-bins.

Start with the expression-blind, frozen sub-bins:

```text
[0.30, 0.35)
[0.35, 0.40)
[0.40, 0.45)
[0.45, 0.49]
```

For every mouse and sub-bin, sum raw counts and record the number of cells, total library size, detected genes, and pseudotime distribution. Pre-specify minimum cell, library-size, detected-gene, and voom-weight criteria after the input audit but before gene- or program-level treatment results are examined. The cell threshold alone is not a sufficient quality rule.

A sub-bin enters the standardized primary estimand only if every ploidy-by-treatment group contains at least two QC-passing mouse pseudobulks. If the four-bin grid fails this support rule, merge adjacent bins once according to the frozen pairing `[0.30, 0.395)` and `[0.395, 0.49]`; do not choose boundaries using expression effects. If fewer than two common-support sub-bins remain, the exact standardized analysis is not estimable. In that event, report the coverage failure and retain the one-pseudobulk window-restricted model below as exploratory rather than calling it matched pseudotime.

For the supported grid, apply a documented gene-expression filter, TMM normalization, and `voomWithQualityWeights`. Fit a cell-means model with a coefficient for every pseudotime-sub-bin-by-ploidy-by-treatment combination and use mouse as the repeated-observation block. Estimate within-mouse correlation with a documented two-pass `duplicateCorrelation` procedure and use robust empirical-Bayes moderation.

```text
expression ~ 0 + pseudotime_subbin:initial_ploidy:treatment
block = sample_id
```

Standardize every group to the same pre-specified reference distribution by giving the eligible sub-bins equal weight. For example, the 2N-origin treatment effect is the equally weighted average of the within-sub-bin `Gem_2N-origin - Control_2N-origin` contrasts. Apply the identical weights to all four groups and to both simple treatment effects before calculating their difference.

Pre-specified gene-level contrasts are:

```text
2N-origin treatment effect:  Gem_2N-origin - Control_2N-origin
4N-origin treatment effect:  Gem_4N-origin - Control_4N-origin
ploidy interaction:         (Gem_4N-origin - Control_4N-origin)
                          - (Gem_2N-origin - Control_2N-origin)
```

Also test one global interaction contrast spanning the complete common-support grid before localizing effects to individual sub-bins. Sub-bin-specific intervals are descriptive unless simultaneous inference is implemented.

As an approximate sensitivity analysis, create one pseudobulk per mouse over `[0.30, 0.49]` and fit:

```text
expression ~ centered_mean_pseudotime + initial_ploidy * treatment
```

This model requires density, quantile, and empirical-CDF overlap plots, not only overlap of mouse means, and must be labeled pseudotime-window-restricted. Fit analogous window-restricted models in the left, right, broad, and late regions to assess localization.

Do not include TGI, endpoint ploidy, apoptosis, cluster composition, or other plausible treatment consequences as primary-model covariates. They can block part of the treatment association or create post-treatment adjustment bias.

## Dose-specific model

On the same common-support sub-bin grid and with the same standardization weights, replace binary treatment with categorical dose:

```text
expression ~ 0 + pseudotime_subbin:initial_ploidy:dose
block = sample_id
```

Extract:

```text
(30_vs_0 in 4N-origin) - (30_vs_0 in 2N-origin)
(120_vs_0 in 4N-origin) - (120_vs_0 in 2N-origin)
```

Report estimates and confidence intervals even when not significant. Regard dose-specific estimates as descriptive because they contain only two treated tumors per dose and tumor-origin group. Directional agreement may be noted, but it is not a robustness requirement because responses may be nonlinear and replication is limited.

## Focused and exploratory response families

Freeze exact gene-set identifiers, versions, memberships, expected directions, response-family assignment, tier assignment, and multiplicity family before rerunning the updated analysis. Separate programs whose expected responses differ instead of creating one mixed replication-stress score. Because PPP, E2F, p53, homologous-recombination, Fanconi, checkpoint, and nucleotide-metabolism programs are not specific to gemcitabine, describe them as contextual response programs rather than drug-specific readouts.

### Tier 1: Original core gemcitabine-mechanism response families

These families are the primary focused families carried forward from the original plan:

- Gemcitabine uptake, activation, and inactivation: nucleoside transport; prodrug phosphorylation and active-metabolite formation; cytidine deamination; nucleotide dephosphorylation and inactivation. This descriptive axis may include independently justified genes such as `SLC29A1`, `SLC28A1`, `SLC28A3`, `DCK`, `CMPK1`, `NME1`, `NME2`, `CDA`, `DCTD`, `NT5C1A`, and other pre-specified `NT5C` family members. RNA abundance in this axis does not measure intracellular dFdCTP exposure.
- Nucleotide supply and putative buffering: de novo pyrimidine synthesis; purine synthesis; deoxyribonucleotide biosynthesis; ribonucleotide reductase activity; nucleoside salvage and transport; pentose-phosphate pathway.
- Replication execution: DNA replication initiation and elongation; replication-fork machinery; E2F targets; S-phase/DNA-synthesis programs.
- Replication-stress and checkpoint response: ATR activation in response to replication stress; CHK1/WEE1-associated checkpoint response; ATM/CHK2-associated checkpoint response; DNA-damage checkpoint.
- Fork repair and replication-associated repair: Fanconi anemia and fork repair; homologous recombination.
- Downstream fate response: p53 response; apoptosis; senescence- or cytostasis-associated expression.

### Tier 2: Upgraded focused response families from exploratory GSEA

The previous exploratory GSEA showed strong enrichment in several axes that are mechanistically close to the expected gemcitabine transcriptional response. These axes are upgraded to focused response families for the rerun, but they must be labeled separately from the original Tier 1 families because the upgrade was motivated by the first exploratory result.

- G2/M and broad cell-cycle checkpoint response: Hallmark G2M checkpoint, Reactome G2/M checkpoints, and broader cell-cycle checkpoint programs.
- Broad DNA repair, nucleotide-excision repair, double-strand-break repair, and DNA recombination: pathways expected to accompany replication-blocking nucleoside analog stress.
- TP53 transcriptional regulation and programmed cell death: TP53-regulated transcription and programmed-cell-death programs that extend the original p53/apoptosis downstream-fate readouts.
- One-carbon, purine, and nucleoside-phosphate metabolism: metabolic support programs connected to nucleotide availability and replication stress buffering.
- Mitotic progression and chromosome segregation: mitotic execution, chromosome segregation, and mitotic-checkpoint programs that capture cell-cycle consequences of replication stress.

### Tier 3: Exploratory-only contextual programs

The following axes remain exploratory-only. They are part of the frozen program universe and should be analyzed with the same program-score and competitive-enrichment machinery as Tier 1 and Tier 2. However, their evidence level remains exploratory: they should not be promoted to focused response-family evidence unless an independent justification is added in a future frozen plan.

- Translation, ribosome, RNA processing, and MYC-associated biosynthetic programs.
- Immune, interferon, antigen-presentation, T-cell, and cytotoxicity-associated programs.
- Cell-cell adhesion, extracellular organization, sensory/synapse-like, and other cell-state-organization programs.

Use independently curated Hallmark, Reactome, and GO gene sets. Keep gemcitabine handling separate from endogenous nucleotide synthesis/salvage and putative buffering. A small custom display panel may aid interpretation, but it should remain descriptive unless its membership and directionality are pre-specified independently of these data. The report must summarize results by response-family tier first, then by family, and then end with an integrated conclusion across tiers.

The analysis should use one unified frozen program universe. Tier assignment is a classification and interpretation label, not a method filter. Program-score modeling, `cameraPR` or equivalent competitive gene-set testing, dose-specific sensitivity, endpoint-ploidy sensitivity, and robustness analyses should be run on all estimable programs in Tier 1, Tier 2, Tier 3, and negative controls. The report should then interpret Tier 1 and Tier 2 as the focused response-family evidence, Tier 3 as exploratory-only program evidence, and negative controls as diagnostics.

## Program-level inference

For each short, pre-specified directional program, define a signed pseudobulk score before treatment results are examined. For example, calculate a fixed weighted mean of normalized, gene-standardized expression, splitting positive and negative regulatory components when a gene set is not directionally coherent. Fit the same sub-bin-standardized repeated-measures model to these scores. This supplies directly interpretable `delta_2N-origin`, `delta_4N-origin`, interaction estimates, and confidence intervals for the primary forest plot.

The implemented score contract is:

```text
gene_z = z-score(TMM-normalized pseudobulk logCPM for each gene across retained mouse-subbin observations)
program_score = sum(pre-specified signed_weight_g * gene_z_g) / sum(abs(pre-specified signed_weight_g))
```

The sign and weight source must be exported with the program membership audit. For replication-execution programs, the signed score is oriented so that higher values represent lower replication-execution expression, matching the expected gemcitabine response direction. For stress-response and fate-response programs, higher values represent higher response-program expression. Contextual programs must be labeled as contextual even when they use positive weights. If a custom panel contains mechanistically opposing components, such as gemcitabine uptake/activation versus inactivation/dephosphorylation, those components must have explicit opposite signs in the frozen configuration.

Retain the signed weighted mean z-score above as the primary program-score definition. To diagnose whether conclusions depend on score aggregation rather than on the matched-pseudotime model, also run the same primary model on three frozen score-definition sensitivities for every estimable program in the unified program universe:

- Signed PCA/eigengene score: multiply each observed gene z-score by the pre-specified sign of its weight, compute PC1 across retained mouse-subbin observations, orient PC1 so that it is positively correlated with the primary signed weighted mean score, and standardize the PC1 scores across retained observations. Do not orient using treatment, ploidy, interaction, `cameraPR`, or GSEA results.
- Signed rank-mean score: within each retained mouse-subbin observation, rank all retained genes by TMM-normalized logCPM, transform ranks to centered percentiles, and calculate the signed weighted mean of the program genes using the same pre-specified signs and weights. The ranking background is the retained gene universe, not genes selected by enrichment results.
- Robust trimmed signed mean z-score: calculate the signed weighted gene-z contributions and, when a program has at least ten observed genes, remove the lowest and highest 10% of signed contributions within each mouse-subbin before taking the weighted mean. For programs with fewer than ten observed genes, fall back to the primary signed weighted mean and record the fallback.

These three additional scores are score-definition sensitivity analyses, not replacements for the primary score. They must use the same mouse-subbin matched model, the same contrasts, the same response-family universe, and the same tier-stratified interpretation rules. Each score method must export program-level and mouse-bin-by-program QC, including genes used, detected program genes, detected fraction, raw program UMI, finite-score status, fallback status, score variance, interaction standard error, confidence interval width, and approximate minimum detectable interaction.

Use `cameraPR` or an equivalent competitive test on the moderated gene-level interaction statistics as complementary enrichment evidence for the same frozen program universe; `camera`-family enrichment results are not themselves program-effect coefficients or confidence intervals. Use ranked GSEA on the interaction moderated-t statistics as an exploratory genome-wide pathway-level analysis and export leading-edge genes. Genome-wide GSEA should then be mapped back to the declared response-family tiers and contextual axes for interpretation.

Define the Tier 1 original core families and Tier 2 upgraded focused families as the focused response-family interpretation set for the rerun. Tier 3 programs are analyzed in parallel with the same methods but interpreted as exploratory-only evidence. Apply FDR correction across all tested programs for a transparent all-universe screen, and also report tier-stratified or focused-set FDR summaries for the declared interpretation families. Treat simple effects as decompositions of those interactions. Correct secondary region, dose, endpoint-ploidy threshold, and continuous endpoint-ploidy families separately. Label genome-wide Hallmark, Reactome, and GO screening, and all Tier 3 contextual programs, as exploratory. A pathway significant in only one tumor-origin group is not evidence of differential response without the interaction.

All score plots must display mouse pseudobulks as the independent observations; cell-level scores must not be plotted or tested as if cells were biological replicates.

Interpret interaction signs together with simple effects:

- If gemcitabine suppresses replication execution in the 2N-origin group but the effect is smaller in the 4N-origin group, then `delta_4N-origin` is less negative and the interaction is positive.
- If gemcitabine induces checkpoint/stress-response-associated expression in the 2N-origin group but less strongly in the 4N-origin group, then `delta_4N-origin` is less positive and the interaction is negative.
- A nonzero interaction with effects in opposite directions is biologically different from a smaller response in the 4N-origin group and must be described explicitly.

Because many curated pathways contain both positive and negative regulators, inspect leading-edge direction and gene-level coherence before assigning a simple biological label to a pathway NES.

## Secondary continuous-pseudotime analysis

Use the existing mouse-by-20-bin pseudobulk data to determine whether an interaction varies along pseudotime. Fit low-complexity group-specific curves with two-pass `duplicateCorrelation` or an explicit mouse random effect, for example:

```text
expression ~ group + group:spline(pseudotime, df = 3)
group = interaction(initial_ploidy, treatment)
block = sample_id
```

First perform a global curve-interaction test. Only then estimate the treatment effect within each tumor-origin group and their difference on a common pseudotime grid. Report a grid segment only when at least two QC-passing mice from every ploidy-by-treatment group support that segment; require three per group where coverage permits. Visually shade or omit unsupported regions. Treat pointwise intervals as descriptive unless simultaneous confidence bands are implemented.

Integrate the fitted interaction over `[0.30, 0.49]` using the same fixed reference distribution used for the primary standardized analysis. Do not compare a uniformly integrated curve with an empirically cell-weighted pseudobulk estimand.

This curve analysis is secondary because many original 0.05-wide sample-bins have sparse group-level mouse coverage. It must not turn smooth extrapolation into apparent matched-pseudotime evidence.

## Robustness and diagnostic analyses

1. Leave one mouse out and refit all primary program interactions.
2. Pre-specify numerical influence flags, such as a leave-one-out sign reversal, a large change relative to the full-model standard error, or loss of common support; do not decide robustness thresholds after viewing results.
3. Compare the standardized primary model with the window-restricted model with and without mean-pseudotime adjustment.
4. Repeat in the broad, left, and right intervals.
5. Examine 30 and 120 mg/kg separately as descriptive estimates.
6. Repeat within cluster 10 as a sensitivity analysis to assess residual cluster-composition dependence.
7. Report pseudobulk cell count, library size, detected genes, voom/quality weight, residual, and leverage diagnostics. Prefer quality weights and influence analysis to equal-cell downsampling, which discards information; any repeated downsampling should be a technical sensitivity analysis with Monte Carlo variability reported.
8. Project treated cells onto a trajectory learned from untreated cells and compare with an independently assigned phase/state label where available. Seurat Phase/S-score/G2M-score provides a weaker same-transcriptome check when no orthogonal phase measurement exists. Also construct, where computationally identifiable, a sensitivity coordinate without the nucleotide/replication genes being tested. These are required diagnostics because using the same transcriptional programs to construct and test pseudotime can mechanically attenuate group differences. If no independent or treatment-invariant coordinate is available, qualify the analysis as conditional on an unverified expression-derived coordinate rather than fully matched pseudotime.
9. Model mouse-level occupancy of each frozen pseudotime region by treatment, initial ploidy, and their interaction as a descriptive companion analysis. This exposes selection/composition changes but does not correct survivor bias.
10. Include pre-specified negative-control program families to detect broad normalization or batch-driven interactions.
11. Treat endpoint-ploidy-stratified results as descriptive and never as baseline predictive evidence.
12. Report expected confidence-interval width or minimum detectable standardized interaction for the program scores using the observed mouse-level variance. With four tumors per pooled-treatment/tumor-origin cell, interaction estimates will be imprecise; with two treated tumors per active dose and tumor-origin group, dose-specific estimates are descriptive.

An interaction should be described as stable only when it is not dominated by one mouse, has adequate common pseudotime support, shows coherent gene-level direction, and survives the pre-specified technical diagnostics. Dose and interval results provide context rather than hard pass/fail criteria. Failure to reject the interaction null is inconclusive with this sample size and must not be presented as ploidy independence.

## Interpretation relative to the in-vitro mechanism

The in-vitro results motivate three distinct questions:

1. Is nucleotide metabolism less transcriptionally perturbed by gemcitabine in captured cells from baseline-4N-origin tumors than in those from baseline-2N-origin tumors?
2. Are transcriptional programs associated with replication stress and checkpoint response less induced in the baseline-4N-origin group?
3. Are transcriptional programs associated with downstream cytostasis or death different in the baseline-4N-origin group?

Matched-pseudotime RNA interactions can provide in-vivo alignment with these response patterns, but they do not measure active dFdCTP, dNTP pools, metabolic flux, physical DNA lesions, protein phosphorylation, or cells lost before capture. In particular:

- a smaller ATR/CHK1-, ATM/CHK2-, or WEE1-associated RNA-program response in the baseline-4N-origin group is compatible with either less upstream stress or a different transcriptional response, but checkpoint activation and damage sensing are primarily phosphorylation- and localization-dependent and cannot be inferred from RNA alone;
- comparatively preserved nucleotide/replication expression in the baseline-4N-origin group is compatible with buffering but does not prove larger functional nucleotide pools;
- expression of gemcitabine transport, activation, or inactivation genes does not measure intracellular dFdCTP exposure;
- fewer captured late-pseudotime cells cannot distinguish accelerated passage, arrest followed by death, or differential survival/capture.

One cross-sectional endpoint also cannot distinguish a weaker response from a delayed response or establish durable senescence/cytostasis. Those claims require longitudinal sampling, lineage tracking, or a time course.

Distinguishing reduced effective damage generation from reduced signaling requires an orthogonal matched-phase hierarchy: gemcitabine uptake and DNA-incorporated gemcitabine; dFdCTP and endogenous dCTP competition; dNTP abundance and isotope-resolved flux; DNA-fiber dynamics; appropriate lesion or strand-break assays; gamma-H2AX, pRPA, pCHK1, and pCHK2 as stress/signaling markers; and phase-resolved apoptosis or TUNEL. Pre-specify metabolite reporting per cell, DNA content, cell volume, and replication demand rather than selecting a normalization after viewing the result.

## Proposed implementation and outputs

Implement this as a standalone workflow that consumes the frozen 04i inputs rather than altering the 04i state-annotation estimand:

```text
Code/in-vivo/04j_pseudotime_treatment_ploidy_programs.R
Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_util.R
Code/in-vivo/04j_pseudotime_treatment_ploidy_programs_config.yaml
```

Minimum output set:

```text
audit/input_manifest.csv
audit/mouse_region_coverage.csv
audit/pseudotime_subbin_common_support.csv
audit/processing_batch_crosswalk.csv
audit/pseudotime_distribution_balance.csv
audit/design_matrix_rank.csv
pseudobulk/sample_subbin_metadata.csv
stats/gene_level_simple_effects_and_interactions.csv
stats/program_score_simple_effects_and_interactions.csv
stats/competitive_gene_set_tests.csv
stats/dose_specific_interactions.csv
stats/pseudotime_region_occupancy_interactions.csv
stats/program_score_precision_and_power.csv
robustness/leave_one_mouse_out.csv
robustness/interval_and_model_sensitivities.csv
robustness/pseudotime_coordinate_sensitivities.csv
figures/program_interaction_forest.pdf
figures/program_effect_heatmap.pdf
figures/continuous_pseudotime_interactions.pdf
04j_pseudotime_treatment_ploidy_programs_report.html
```

The principal figure should show, for every pre-specified program score, the 2N-origin treatment association, 4N-origin treatment association, and treatment-by-tumor-origin interaction with confidence intervals. A second heatmap may show interaction estimates across regions and doses. Continuous-pseudotime plots must include observed-support markings and mouse-level uncertainty.

## Completion criteria

The workflow is complete when:

1. inputs and frozen estimands are recorded reproducibly;
2. all primary results use mice, not cells, as biological replicates;
3. simple effects and the difference-in-differences interaction are exported together;
4. exact pathway definitions, score directions, negative controls, and multiplicity families are frozen;
5. processing-batch separability and common pseudotime support are documented;
6. unsupported pseudotime regions are not extrapolated or interpreted;
7. dose, interval, cluster, coordinate, and leave-one-mouse-out sensitivities are reported;
8. the report identifies the outcome-selected sequenced cohort and state-conditioned survivor estimand;
9. the report distinguishes transcriptional alignment from proof of active-metabolite exposure, damage burden, checkpoint activation, or nucleotide buffering.
