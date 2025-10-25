library(testthat)

test_that(".glmnet_select_terms keeps factor terms under binomial family", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  n <- 120
  x1 <- rnorm(n)
  fac <- factor(sample(c("a", "b"), n, replace = TRUE))
  eta <- -1 + 2 * as.numeric(fac == "b")
  prob <- 1 / (1 + exp(-eta))
  y <- rbinom(n, size = 1, prob = prob)
  dat <- data.frame(y = factor(y), x1 = x1, fac = fac)
  terms <- c("x1", "fac")
  sel <- .glmnet_select_terms(dat, "y", terms, alpha = 1, family = "binomial")
  expect_true("fac" %in% sel)
})