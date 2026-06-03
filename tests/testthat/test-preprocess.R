test_that("align_samples keeps common samples in expression order", {
  expr <- matrix(1:6, nrow = 2)
  rownames(expr) <- c("A", "B")
  colnames(expr) <- c("S1", "S2", "S3")
  clinical <- data.frame(response = c("response", "non_response"), row.names = c("S3", "S1"))

  aligned <- align_samples(expr, clinical)

  expect_equal(colnames(aligned$expr), c("S1", "S3"))
  expect_equal(rownames(aligned$clinical), c("S1", "S3"))
})

test_that("make_response_vector builds a binary response vector", {
  clinical <- data.frame(response = c("response", "non_response", "response"))

  expect_equal(make_response_vector(clinical), c(1L, 0L, 1L))
})
