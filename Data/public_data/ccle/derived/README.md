# Derived Inputs

`ccle_expression_columns.tsv` contains the CCLE sample column names from the original `CCLE.rpkm.2016-06-17f.gct` expression matrix after applying the historical script's column-name cleanup:

```r
gsub("-Tumor$", "", gsub("^fh_", "", colnames(ex)))
```

The drug-ploidy barplot only uses these names to reproduce the original expression-availability cell-line filter; it does not use expression values.
