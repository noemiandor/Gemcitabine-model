Full 2-fold analysis was executed using the following computational steps:

1. Load metabolomics data file (e.g. Metabolomics_2N_4N_Full.xlsm).
2. Detect replicate columns matching ^(2N|4N)_[CG][1-4]$.
3. Impute missing/zero values by half of the row-wise minimum positive value.
4. Log2-transform intensities.
5. Run Welch t-tests for:
   - 2N_G vs 2N_C
   - 4N_G vs 4N_C
   - 4N_C vs 2N_C
   - 4N_G vs 2N_G
6. Run two-way ANOVA:
   log2_abundance ~ ploidy + treatment + ploidy:treatment
7. Define differential response:
   delta_log2FC = log2FC(2N_G/2N_C) - log2FC(4N_G/4N_C)
8. Define 2-fold hits using p<0.05 and |log2FC|>=1.
9. Order heatmap rows:
   2N up, 2N down, 4N up, 4N down, Common up, Common down, Opposite, Interaction only.
10. Generate:
   - PCA for samples
   - volcano plots with black insignificant and red significant points
   - row z-score heatmap
   - replicate-level fold-change heatmap
   - top interaction heatmap
   - pathway-class enrichment heatmap
