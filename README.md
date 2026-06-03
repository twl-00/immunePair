# immunePairMarker

`immunePairMarker` is an R package for pairwise gene marker analysis in immune
response studies. It screens gene pairs using pairwise expression comparisons,
chi-square testing, Fisher's exact test, and multiple-testing correction.

## Installation

You can install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("twl-00/immunePairMarker", upgrade = "never")
```

## Quick Start

```r
library(immunePairMarker)

expr_file <- system.file(
  "extdata",
  "example_exp.txt",
  package = "immunePairMarker"
)

clinical_file <- system.file(
  "extdata",
  "example_clinical.txt",
  package = "immunePairMarker"
)

result <- run_pair_marker_analysis(
  expr_file = expr_file,
  clinical_file = clinical_file,
  out_dir = tempdir(),
  response_col = "response",
  response_label = "response",
  main_delta = 0.25,
  delta_list = c(0, 0.25, 0.5),
  dataset_name = "example"
)

result$sensitivity
head(result$pairs)
head(result$sig_pairs)
```

## Input Files

The expression file should be a tab-delimited text file. The first column should
contain gene names, and the remaining columns should contain samples.

```text
gene    S1    S2    S3
GENE_A  12    11    10
GENE_B  2     1     2
```

The clinical file should be a tab-delimited text file. The first column should
contain sample IDs, and one column should contain the response label.

```text
sample  response
S1      response
S2      non_response
```

Sample IDs in the clinical file should match sample names in the expression
matrix.

## Main Output

`run_pair_marker_analysis()` returns a list containing:

- `sensitivity`: number of gene pairs passing chi-square screening under each
  delta cutoff
- `pairs`: all screened gene pairs with chi-square p-values, Fisher p-values,
  odds ratios, and adjusted p-values
- `sig_pairs`: significant gene pairs after adjusted p-value filtering
- `output_paths`: paths of written result files

## Citation

If you use `immunePairMarker` in your research, please cite the GitHub
repository:

```text
https://github.com/twl-00/immunePairMarker
```

## License

This package is licensed under the MIT License.
