# Manual GDSC Drug-Class Workbook

`drug_class_final_used_with_primary_secondary_corrected.xlsx` is the reviewed
drug-class assignment workbook used by the GDSC ploidy enrichment workflow.

The workbook was compiled and reviewed from PubChem and DrugBank
drug-class/mechanism annotations. The analysis treats the workbook as a frozen
manual input: runtime code does not query PubChem or DrugBank, and the legacy
PubChem cache and fallback class-generation paths are not consumed by the main
workflow.

The enrichment grouping variable is `primary_anticancer_class`. The
manuscript-facing collapsed heatmap uses the workbook's `Drug counts` sheet:
`Combine.if.needed` is used where populated, otherwise
`primary_anticancer_class` is retained.
