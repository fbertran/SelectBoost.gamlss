# simple smoke tests (skip on CRAN)
test_that("sb_gamlss runs on a tiny example", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  set.seed(1)
  n <- 80
  x1 <- rnorm(n); x2 <- rnorm(n)
  y  <- gamlss.dist::rNO(n, mu = 1 + 1.2*x1, sigma = 1)
  dat <- data.frame(y, x1, x2)
  fit <- sb_gamlss(
    y ~ 1, data = dat, family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2, B = 10, sample_fraction = 0.8,
    pi_thr = 0.4, pre_standardize = TRUE, trace = FALSE
  )
  expect_s3_class(fit, "sb_gamlss")
  expect_true(is.list(fit$final_formula))
  expect_true(NROW(fit$selection) >= 1)
  # AICc helper
  expect_true(is.numeric(AICc_gamlss(fit$final_fit)))
})
