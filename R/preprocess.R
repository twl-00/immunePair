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

#' Build binary phenotype vector
#'
#' @param clinical Clinical data frame.
#' @param phenotype_col Column containing phenotype labels.
#' @param positive_label Label coded as 1.
#' @param negative_label Optional label coded as 0. If `NULL`, all samples not
#'   matching `positive_label` are coded as 0.
#' @return Integer vector: 1 for the positive class, 0 for the negative class,
#'   and `NA` for samples outside `positive_label` and `negative_label` when
#'   `negative_label` is supplied.
#' @export
make_binary_label_vector <- function(
    clinical,
    phenotype_col = "response",
    positive_label = "response",
    negative_label = NULL) {
  
  if (!phenotype_col %in% names(clinical)) {
    stop("phenotype_col not found in clinical data: ", phenotype_col, call. = FALSE)
  }
  
  phenotype_status <- as.character(clinical[[phenotype_col]])
  labels <- unique(stats::na.omit(phenotype_status))
  
  if (!positive_label %in% labels) {
    stop("positive_label not found in phenotype_col: ", positive_label, call. = FALSE)
  }
  
  if (!is.null(negative_label) && !negative_label %in% labels) {
    stop("negative_label not found in phenotype_col: ", negative_label, call. = FALSE)
  }
  
  if (is.null(negative_label)) {
    out <- ifelse(phenotype_status == positive_label, 1L, 0L)
  } else {
    out <- rep(NA_integer_, length(phenotype_status))
    out[phenotype_status == positive_label] <- 1L
    out[phenotype_status == negative_label] <- 0L
  }
  
  names(out) <- rownames(clinical)
  out
}

#' Build binary response vector
#'
#' @param clinical Clinical data frame.
#' @param response_col Column containing response labels.
#' @param response_label Label treated as response.
#' @return Integer vector: 1 for response, 0 for non-response.
#' @export
make_response_vector <- function(clinical, response_col = "response", response_label = "response") {
  unname(make_binary_label_vector(
    clinical = clinical,
    phenotype_col = response_col,
    positive_label = response_label,
    negative_label = NULL
  ))
}
