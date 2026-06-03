#' Add Fisher exact test and adjusted p-values
#'
#' @param pair_df Pairwise chi-square result data frame.
#' @param adjust_method Method passed to [p.adjust()].
#' @return Pair result data frame with `OR`, `fisher_p`, and `adjusted_p`.
#' @export
add_fisher_test <- function(pair_df, adjust_method = "BH") {
  required_cols <- c("nr0", "nr1", "r0", "r1")
  missing_cols <- setdiff(required_cols, names(pair_df))

  if (length(missing_cols) > 0) {
    stop("pair_df is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  if (nrow(pair_df) == 0) {
    pair_df$OR <- numeric()
    pair_df$fisher_p <- numeric()
    pair_df$adjusted_p <- numeric()
    return(pair_df)
  }

  do_fisher <- function(nr0, nr1, r0, r1) {
    tab <- matrix(c(nr0, nr1, r0, r1), nrow = 2, byrow = TRUE)
    ft <- stats::fisher.test(tab)
    c(OR = unname(ft$estimate), p = ft$p.value)
  }

  fisher_mat <- t(mapply(
    do_fisher,
    pair_df$nr0,
    pair_df$nr1,
    pair_df$r0,
    pair_df$r1
  ))

  fisher_mat <- as.matrix(fisher_mat)
  colnames(fisher_mat) <- c("OR", "p")

  pair_df$OR <- fisher_mat[, "OR"]
  pair_df$fisher_p <- fisher_mat[, "p"]
  pair_df$adjusted_p <- stats::p.adjust(pair_df$fisher_p, method = adjust_method)
  pair_df[order(pair_df$fisher_p), , drop = FALSE]
}
