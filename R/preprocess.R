#' Filter expression matrix by non-zero expression
#'
#' @param expr Expression data frame or matrix with genes in rows and samples in columns.
#' @param min_nonzero_prop Minimum proportion of samples with expression greater than zero.
#' @return A numeric matrix.
#' @export
filter_expression_matrix <- function(expr, min_nonzero_prop = 0.5) {
  mat <- as.matrix(expr)
  suppressWarnings(storage.mode(mat) <- "numeric")

  if (is.null(rownames(mat))) {
    stop("Expression matrix must have gene row names.", call. = FALSE)
  }

  keep <- rowSums(mat > 0, na.rm = TRUE) > ncol(mat) * min_nonzero_prop
  mat <- mat[keep, , drop = FALSE]
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  mat
}

#' Align expression and clinical samples
#'
#' @param expr Expression matrix with samples in columns.
#' @param clinical Clinical data frame with samples in row names.
#' @return A list with aligned `expr` and `clinical`.
#' @export
align_samples <- function(expr, clinical) {
  overlap_sample <- intersect(colnames(expr), rownames(clinical))

  if (length(overlap_sample) == 0) {
    stop("No overlapping samples between expression columns and clinical row names.", call. = FALSE)
  }

  list(
    expr = expr[, overlap_sample, drop = FALSE],
    clinical = clinical[overlap_sample, , drop = FALSE]
  )
}

#' Build binary response vector
#'
#' @param clinical Clinical data frame.
#' @param response_col Column containing response labels.
#' @param response_label Label treated as response.
#' @return Integer vector: 1 for response, 0 for non-response.
#' @export
make_response_vector <- function(clinical, response_col = "response", response_label = "response") {
  if (!response_col %in% names(clinical)) {
    stop("response_col not found in clinical data: ", response_col, call. = FALSE)
  }

  response_status <- factor(clinical[[response_col]])
  if (!response_label %in% levels(response_status)) {
    stop("response_label not found in response_col: ", response_label, call. = FALSE)
  }

  ifelse(response_status == response_label, 1L, 0L)
}
