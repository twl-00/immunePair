#' Run the complete pair marker analysis workflow
#'
#' @param expr_file Path to expression table.
#' @param clinical_file Path to clinical table.
#' @param out_dir Output directory. If `NULL`, results are not written.
#' @param gene_col Gene column name in expression table.
#' @param response_col Clinical response column.
#' @param response_label Label treated as response.
#' @param main_delta Delta cutoff for the main analysis.
#' @param delta_list Delta cutoffs for sensitivity analysis.
#' @param dataset_name Prefix for output file names.
#' @param sig_cutoff Adjusted p-value cutoff for significant pairs.
#' @param min_nonzero_prop Minimum non-zero expression proportion.
#' @param min_prop Minimum pairwise comparison proportion.
#' @param max_prop Maximum pairwise comparison proportion.
#' @param chisq_cutoff Chi-square p-value cutoff.
#' @param min_valid_prop Minimum valid sample proportion.
#' @return A list with filtered expression matrix, aligned clinical data,
#'   response vector, sensitivity summary, all pairs, significant pairs, and
#'   output paths.
#' @export
run_pair_marker_analysis <- function(
    expr_file,
    clinical_file,
    out_dir = NULL,
    gene_col = NULL,
    response_col = "response",
    response_label = "response",
    main_delta = 0.25,
    delta_list = c(0, 0.25, 0.5),
    dataset_name = tools::file_path_sans_ext(basename(expr_file)),
    sig_cutoff = 0.05,
    min_nonzero_prop = 0.5,
    min_prop = 0.05,
    max_prop = 0.95,
    chisq_cutoff = 0.01,
    min_valid_prop = 0.5) {

  expr <- read_expression_data(expr_file, gene_col = gene_col)
  expr <- filter_expression_matrix(expr, min_nonzero_prop = min_nonzero_prop)
  clinical <- read_clinical_data(clinical_file)

  aligned <- align_samples(expr, clinical)
  expr <- aligned$expr
  clinical <- aligned$clinical
  resp <- make_response_vector(
    clinical,
    response_col = response_col,
    response_label = response_label
  )

  message("Filtered genes: ", nrow(expr), "; aligned samples: ", ncol(expr))

  sensitivity <- run_delta_sensitivity(
    expr = expr,
    resp = resp,
    delta_list = delta_list,
    min_prop = min_prop,
    max_prop = max_prop,
    chisq_cutoff = chisq_cutoff,
    min_valid_prop = min_valid_prop
  )

  pairs <- run_pairwise_chisq(
    expr = expr,
    resp = resp,
    min_prop = min_prop,
    max_prop = max_prop,
    chisq_cutoff = chisq_cutoff,
    diff_cutoff = main_delta,
    min_valid_prop = min_valid_prop
  )

  if (nrow(pairs) == 0) {
    warning("No gene pairs passed chi-square screening.", call. = FALSE)
  }
  pairs <- add_fisher_test(pairs)

  sig_pairs <- pairs[pairs$adjusted_p < sig_cutoff, , drop = FALSE]
  output_paths <- list(sensitivity = NULL, all_pairs = NULL, significant_pairs = NULL)

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE)
    }

    output_paths$sensitivity <- file.path(
      out_dir,
      paste0(dataset_name, "_delta_sensitivity_summary.txt")
    )
    output_paths$all_pairs <- file.path(
      out_dir,
      paste0(dataset_name, "_chisq_fisher_bh.txt")
    )
    output_paths$significant_pairs <- file.path(
      out_dir,
      paste0(dataset_name, "_sig_pairs_adjP", sig_cutoff, ".txt")
    )

    write.table(
      sensitivity,
      file = output_paths$sensitivity,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    write.table(
      pairs,
      file = output_paths$all_pairs,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )

    if (nrow(sig_pairs) > 0) {
      write.table(
        sig_pairs,
        file = output_paths$significant_pairs,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  list(
    expr = expr,
    clinical = clinical,
    response = resp,
    sensitivity = sensitivity,
    pairs = pairs,
    sig_pairs = sig_pairs,
    output_paths = output_paths
  )
}
