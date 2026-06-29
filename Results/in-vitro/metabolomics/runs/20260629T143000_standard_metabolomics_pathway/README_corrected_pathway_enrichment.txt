
Corrected pathway enrichment heatmap
====================================

Main correction:
The prior script used permissive raw substring matching:
    if any(k.lower() in metabolite_name for k in keywords)

This could falsely assign pathways when short tokens such as sam, sah, amp, cmp, udp, or nad
appeared inside unrelated metabolite names. This script replaces that with explicit
component-level annotations and token-aware matching for complete biochemical abbreviations.

Input features after preprocessing: 472
Features with changed annotation compared with old raw-substring approach: 57

Main outputs:
- figures/corrected_curated_pathway_enrichment_heatmap_2fold.png
- figures/corrected_curated_pathway_enrichment_heatmap_2fold.pdf
- tables/corrected_curated_pathway_enrichment_2fold.csv
- tables/corrected_curated_pathway_enrichment_heatmap_matrix.csv
- tables/annotated_metabolites_with_response_stats_corrected.csv
- tables/annotation_change_report_raw_substring_vs_corrected.csv
- tables/annotation_change_summary.csv

Interpretation:
This is a corrected curated, name-based pathway-class enrichment. It is still not a substitute
for complete KEGG/HMDB ID-mapped pathway enrichment, but it removes the most important false
positive mechanism from raw substring matching.
