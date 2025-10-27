library(testthat)

test_that(".expand_terms_to_mm covers pbm and lo smoothers", {
  terms <- c("pbm(x1)", "lo(x2, span = 0.7)", "x1:x3")
  expanded <- .expand_terms_to_mm(terms, df_smooth = 4L)
  expect_true(any(grepl("splines::bs", expanded)))
  expect_true(any(grepl("splines::ns", expanded)))
  expect_true("x1:x3" %in% expanded)
})

test_that("custom converters can override expansion", {
  conv <- list(function(term, df_smooth) {
    if (grepl("special", term)) return("I(x4^2)")
    NULL
  })
  on.exit(options(SelectBoost.gamlss.term_converters = NULL), add = TRUE)
  options(SelectBoost.gamlss.term_converters = conv)
  expanded <- .expand_terms_to_mm(c("special(x4)"), df_smooth = 5L)
  expect_identical(expanded, "I(x4^2)")
})

test_that("group designs keep column-to-term mapping", {
  skip_if_not_installed("splines")
  set.seed(123)
  n <- 15
  dat <- data.frame(
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    f = factor(sample(letters[1:3], n, replace = TRUE))
  )
  pb <- function(x, ...) x
  cs <- function(x, ...) x
  mu_terms <- attr(stats::terms(~ f + x1 + pb(x2)), "term.labels")
  gd_mu <- SelectBoost.gamlss:::.mu_group_design(dat, mu_terms, df_smooth = 4L)
  expect_length(gd_mu$groups, ncol(gd_mu$X))
  expect_identical(sort(unique(gd_mu$groups)), seq_along(mu_terms))
  
  sigma_terms <- attr(stats::terms(~ f + cs(x1)), "term.labels")
  gd_sigma <- SelectBoost.gamlss:::.param_group_design(dat, sigma_terms, df_smooth = 4L)
  expect_length(gd_sigma$groups, ncol(gd_sigma$X))
  expect_identical(sort(unique(gd_sigma$groups)), seq_along(sigma_terms))
})
