library(testthat)

skip_if_future_missing <- function() {
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("future.apply", quietly = TRUE)) {
    skip("future packages not available")
  }
}

test_that(".sb_parallel_lapply restores previous plan after multisession", {
  skip()
  skip_if_future_missing()
  old_plan <- try(future::plan(), silent = TRUE)
  on.exit({
    if (!inherits(old_plan, "try-error")) future::plan(old_plan) else future::plan(future::sequential)
  }, add = TRUE)
  future::plan(future::sequential)
  baseline <- future::plan()
  res <- SelectBoost.gamlss:::.sb_parallel_lapply(3, function(i) i, parallel = "multisession", workers = 1)
  expect_identical(res, list(3))
  plan_after <- future::plan()
  expect_identical(plan_after, baseline)
})
