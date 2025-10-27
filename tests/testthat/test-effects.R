skip_if_not_installed <- function(pkg) {
  testthat::skip_if_not_installed(pkg)
}

library(testthat)


test_that("effect_plot handles factor predictors with ggplot2", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("ggplot2")
  set.seed(123)
  n <- 80
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f = factor(sample(letters[1:3], n, replace = TRUE))
  )
  fit <- gamlss::gamlss(y ~ x + f, 
                        data = dat, 
                        family = gamlss.dist::NO())
  plt <- effect_plot(fit, "f", data = dat, what = "mu")
  expect_s3_class(plt, "ggplot")
  df_layer <- ggplot2::layer_data(plt)
  expect_true(all(c("x", "y") %in% names(df_layer)))
})

test_that("effect_plot handles factor predictors with ggplot2", {
  skip_if_not_installed("gamlss")
  skip_if_not_installed("ggplot2")
  set.seed(123)
  n <- 80
  dat <- data.frame(
    y = rnorm(n),
    x = rnorm(n),
    f = factor(sample(letters[1:3], n, replace = TRUE))
  )
  fit <- gamlss::gamlss(y ~ x + f, 
                        data = dat, 
                        family = gamlss.dist::NO())
  plt <- effect_plot(fit, "f", data = dat, what = "not_expected")
  expect_s3_class(plt, "ggplot")
  df_layer <- ggplot2::layer_data(plt)
  expect_true(all(c("x", "y") %in% names(df_layer)))
})
