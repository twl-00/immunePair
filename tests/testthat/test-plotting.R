test_that("pair state helper encodes pairwise expression states", {
  expr <- matrix(
    c(10, 10, 1, 1, 1, 1, 10, 10),
    nrow = 2,
    byrow = TRUE
  )
  rownames(expr) <- c("A", "B")
  colnames(expr) <- paste0("S", 1:4)

  pair_state_and_score <- getFromNamespace(".pair_state_and_score", "PairMarker")
  pair <- pair_state_and_score(expr, "A", "B", delta = 0.25)

  expect_equal(unname(pair$state), c(1L, 1L, 0L, 0L))
  expect_equal(names(pair$state), colnames(expr))
})

test_that("ROC helper returns high AUC for ordered scores", {
  response <- c(1L, 1L, 0L, 0L)
  score <- c(0.9, 0.8, 0.2, 0.1)

  compute_roc <- getFromNamespace(".compute_roc", "PairMarker")
  roc <- compute_roc(response, score)

  expect_equal(roc$auc, 1)
  expect_false(roc$reversed)
})

test_that("HR labels do not print rounded zero", {
  format_hr <- getFromNamespace(".format_hr", "PairMarker")

  expect_equal(format_hr(0.0001), "< 0.01")
  expect_equal(format_hr(0.8), "= 0.80")
})

test_that("plot pair selector prefers significant pairs", {
  result <- list(
    pairs = data.frame(
      gene1 = c("A", "C"),
      gene2 = c("B", "D"),
      adjusted_p = c(0.02, 0.03)
    ),
    sig_pairs = data.frame(
      gene1 = "E",
      gene2 = "F",
      adjusted_p = 0.001
    )
  )

  select_pair_table <- getFromNamespace(".select_pair_table", "PairMarker")
  selected <- select_pair_table(result, top_n = 1)

  expect_equal(selected$gene1, "E")
  expect_equal(selected$gene2, "F")
})

test_that("integrated plot selector supports top rank and pair_id", {
  integrated <- list(
    summary = data.frame(
      pair_id = c("A|B", "C|D"),
      gene1 = c("A", "C"),
      gene2 = c("B", "D"),
      n_dataset = c(3, 2)
    )
  )

  select_integrated <- getFromNamespace(".select_integrated_plot_pairs", "PairMarker")

  top_pair <- select_integrated(integrated, top_n = 1)
  expect_equal(top_pair$gene1, "A")
  expect_equal(top_pair$gene2, "B")

  pair_by_id <- select_integrated(integrated, pair_id = "C|D")
  expect_equal(pair_by_id$gene1, "C")
  expect_equal(pair_by_id$gene2, "D")
})
