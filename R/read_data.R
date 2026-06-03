#' Read expression data
#'
#' @param file Path to a tab-delimited expression file.
#' @param gene_col Name of the gene column. If `NULL`, the first column is used
#'   when row names are not already meaningful.
#' @param sep Field separator.
#' @return A data frame with genes as row names and samples as columns.
#' @export
read_expression_data <- function(file, gene_col = NULL, sep = "\t") {
  dat <- read.table(
    file,
    sep = sep,
    row.names = NULL,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (is.null(gene_col)) {
    candidate_cols <- intersect(c("gene", "Gene", "GENE", "symbol", "Symbol"), names(dat))
    if (length(candidate_cols) > 0) {
      gene_col <- candidate_cols[1]
    } else {
      gene_col <- names(dat)[1]
    }
  }

  if (!gene_col %in% names(dat)) {
    stop("gene_col not found in expression data: ", gene_col, call. = FALSE)
  }

  genes <- dat[[gene_col]]
  keep <- !is.na(genes) & genes != "" & !duplicated(genes)
  dat <- dat[keep, , drop = FALSE]
  genes <- genes[keep]
  dat[[gene_col]] <- NULL

  rownames(dat) <- genes
  dat
}

#' Read clinical data
#'
#' @param file Path to a tab-delimited clinical file.
#' @param sep Field separator.
#' @return A data frame with samples as row names.
#' @export
read_clinical_data <- function(file, sep = "\t") {
  read.table(
    file,
    sep = sep,
    row.names = 1,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
