# GDSC Primary/Secondary Drug-Class Simplification Plan

Status: implemented and archived.

The executed implementation plan is retained at:

```text
docs/Archive/gdsc_primary_secondary_drug_class_simplification_plan.md
```

The active Figure 1 GDSC workflow now uses this reviewed workbook as its drug-class source of truth:

```text
Code/gdsc_ploidy_analysis/data/manual/drug_class_final_used_with_primary_secondary_corrected.xlsx
```

The post-cleanup validation run regenerated the GDSC outputs here:

```text
Results/public_data/gdsc_ploidy_analysis/runs/post_cleanup_check_gdsc
```

Current entrypoints:

```sh
Rscript Code/gdsc_ploidy_analysis/run_gdsc_ploidy_analysis.R --analysis-mode=manuscript
bash Manager.sh --mode standard --modules gdsc --run-id <run_id>
```
