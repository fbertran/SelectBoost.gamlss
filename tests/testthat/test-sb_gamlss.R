# simple smoke tests (skip on CRAN)
test_that("sb_gamlss runs on a tiny example", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
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

test_that("sb_gamlss correlated resampling path produces selection metadata", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  set.seed(2)
  n <- 40
  base <- rnorm(n)
  x1 <- base + rnorm(n, sd = 0.1)
  x2 <- base + rnorm(n, sd = 0.1)
  y  <- 0.5 + 1.5 * x1 + rnorm(n)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- sb_gamlss(
    y ~ 1,
    data = dat,
    family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2,
    B = 6,
    sample_fraction = 0.75,
    use_groups = TRUE,
    c0 = 0.3,
    pi_thr = 0.4,
    trace = FALSE
  )
  tab <- selection_table(fit)
  expect_true(all(c("parameter", "term", "count", "prop") %in% names(tab)))
  expect_true(any(tab$parameter == "mu"))
  expect_true(any(tab$term %in% c("x1", "x2")))
})

test_that("sb_gamlss handles cbind binomial responses", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  set.seed(3)
  n <- 60
  x1 <- rnorm(n)
  eta <- 0.5 + 0.8 * x1
  trials <- 6L
  success <- rbinom(n, size = trials, prob = plogis(eta))
  failure <- trials - success
  dat <- data.frame(success = success, failure = failure, x1 = x1)
  fit <- sb_gamlss(
    cbind(success, failure) ~ 1,
    data = dat,
    family = gamlss.dist::BI(),
    mu_scope = ~ x1,
    B = 6,
    sample_fraction = 0.7,
    pi_thr = 0.4,
    trace = FALSE
  )
  expect_s3_class(fit, "sb_gamlss")
  expect_true(nrow(selection_table(fit)) >= 1)
})

test_that("sb_gamlss subsets weights during resampling", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  set.seed(4)
  n <- 50
  x1 <- rnorm(n)
  y <- gamlss.dist::rNO(n, mu = 1 + 0.5 * x1, sigma = 1)
  w <- seq_len(n)
  dat <- data.frame(y = y, x1 = x1)
  expect_no_error({
    sb_gamlss(
      y ~ 1,
      data = dat,
      family = gamlss.dist::NO(),
      mu_scope = ~ x1,
      B = 5,
      sample_fraction = 0.6,
      pi_thr = 0.4,
      weights = w,
      trace = FALSE
    )
  })
})

test_that("sb_gamlss works with glmnet engine", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  skip_if_not_installed("glmnet")
  set.seed(5)
  n <- 60
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- gamlss.dist::rNO(n, mu = 0.3 + 1.1 * x1 - 0.6 * x2, sigma = 1)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- sb_gamlss(
    y ~ 1,
    data = dat,
    family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2,
    B = 4,
    sample_fraction = 0.75,
    engine = "glmnet",
    pi_thr = 0.3,
    trace = FALSE
  )
  expect_s3_class(fit, "sb_gamlss")
})

test_that("sb_gamlss works with grpreg engine", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  skip_if_not_installed("grpreg")
  set.seed(6)
  n <- 60
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- gamlss.dist::rNO(n, mu = -0.2 + 0.9 * x1 + 0.7 * x2, sigma = 1)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- sb_gamlss(
    y ~ 1,
    data = dat,
    family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2,
    B = 4,
    sample_fraction = 0.75,
    engine = "grpreg",
    pi_thr = 0.3,
    trace = FALSE
  )
  expect_s3_class(fit, "sb_gamlss")
})

test_that("sb_gamlss works with sgl engine", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  skip_if_not_installed("SelectBoost")
  skip_if_not_installed("SGL")
  set.seed(7)
  n <- 60
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- gamlss.dist::rNO(n, mu = 0.1 + 1.2 * x1 - 0.5 * x2, sigma = 1)
  dat <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- sb_gamlss(
    y ~ 1,
    data = dat,
    family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2,
    B = 4,
    sample_fraction = 0.75,
    engine = "sgl",
    pi_thr = 0.3,
    trace = FALSE
  )
  expect_s3_class(fit, "sb_gamlss")
})
