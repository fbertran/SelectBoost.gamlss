test_that("fast_cor matches cor() behaviour for single observation", {
  df <- data.frame(a = 1, b = 2, c = 3)
  res <- SelectBoost.gamlss:::fast_cor(df, names(df))
  expect_true(all(is.na(res)))
})

test_that("fast_cor returns NA when a variable has zero variance", {
  df <- data.frame(a = c(1, 1, 1), b = c(1, 2, 3))
  res <- SelectBoost.gamlss:::fast_cor(df, names(df))
  rownames(res) <- names(df)
  colnames(res) <- names(df)
  expect_true(is.na(res["a", "a"]))
  expect_true(is.na(res["a", "b"]))
  expect_true(is.na(res["b", "a"]))
  expect_equal(res["b", "b"], 1)
})
