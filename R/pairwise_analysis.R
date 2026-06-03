#' Run C++ pairwise chi-square screening
#'
#' @param expr Numeric expression matrix with genes in rows and samples in columns.
#' @param resp Integer response vector with 1 for response and 0 for non-response.
#' @param min_prop Minimum proportion of pairwise comparison value 1.
#' @param max_prop Maximum proportion of pairwise comparison value 1.
#' @param chisq_cutoff Chi-square p-value cutoff.
#' @param diff_cutoff Minimum expression difference for assigning pairwise order.
#' @param min_valid_prop Minimum proportion of valid samples per gene pair.
#' @return A data frame with screened gene pairs and chi-square p-values.
#' @export
run_pairwise_chisq <- function(
    expr,
    resp,
    min_prop = 0.05,
    max_prop = 0.95,
    chisq_cutoff = 0.01,
    diff_cutoff = 0.25,
    min_valid_prop = 0.5) {

  mat_num <- log2(as.matrix(expr) + 1)
  res_cpp <- pairwise_chisq_cpp(
    mat = mat_num,
    resp = as.integer(resp),
    min_prop = min_prop,
    max_prop = max_prop,
    chisq_cutoff = chisq_cutoff,
    diff_cutoff = diff_cutoff,
    min_valid_prop = min_valid_prop
  )

  if (length(res_cpp$g1) == 0) {
    return(data.frame(
      gene1 = character(),
      gene2 = character(),
      nr0 = integer(),
      nr1 = integer(),
      r0 = integer(),
      r1 = integer(),
      chisq_p = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  data.frame(
    gene1 = rownames(expr)[res_cpp$g1],
    gene2 = rownames(expr)[res_cpp$g2],
    nr0 = res_cpp$nr0,
    nr1 = res_cpp$nr1,
    r0 = res_cpp$r0,
    r1 = res_cpp$r1,
    chisq_p = res_cpp$chisq_p,
    stringsAsFactors = FALSE
  )
}

#' Run delta sensitivity analysis
#'
#' @param expr Numeric expression matrix.
#' @param resp Integer response vector.
#' @param delta_list Delta cutoffs to test.
#' @param ... Additional arguments passed to `run_pairwise_chisq()`.
#' @return A data frame with pair counts for each delta.
#' @export
run_delta_sensitivity <- function(expr, resp, delta_list = c(0, 0.25, 0.5), ...) {
  rows <- lapply(delta_list, function(delta) {
    res <- run_pairwise_chisq(expr, resp, diff_cutoff = delta, ...)
    data.frame(delta = delta, n_pairs_chisq = nrow(res))
  })

  do.call(rbind, rows)
}
