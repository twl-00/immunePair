# immunePairMarker

`immunePairMarker` wraps the original one-step analysis script into reusable R
package functions.

## Install

```r
pkg_dir <- file.path("D:", "研究生", "免疫应答标志物", "数据分析", "immunePairMarker")
devtools::install(pkg_dir)
```

## Run the full workflow

```r
library(immunePairMarker)

result <- run_pair_marker_analysis(
  expr_file = file.path("D:", "研究生", "免疫应答标志物", "数据分析", "data", "GSE145996", "GSE145996_exp.txt"),
  clinical_file = file.path("D:", "研究生", "免疫应答标志物", "数据分析", "data", "GSE145996", "GSE145996_clinical.txt"),
  out_dir = file.path("D:", "研究生", "免疫应答标志物", "数据分析", "result"),
  response_col = "response",
  main_delta = 0.25,
  delta_list = c(0, 0.25, 0.5)
)
```
