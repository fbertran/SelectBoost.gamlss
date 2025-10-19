test_that("confidence functionals exist and run", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  set.seed(1)
  n <- 120
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  y  <- gamlss.dist::rNO(n, mu = 1 + 1.1*x1 - 0.7*x3, sigma = 1)
  dat <- data.frame(y, x1, x2, x3)
  g <- sb_gamlss_c0_grid(
    y ~ 1, data = dat, family = gamlss.dist::NO(),
    mu_scope = ~ x1 + x2 + x3, sigma_scope = ~ x1 + x2,
    c0_grid = seq(0.3, 0.7, by = 0.2), B = 10, pi_thr = 0.5, pre_standardize = TRUE, trace = FALSE
  )
  cf <- confidence_functionals(g, pi_thr = 0.5, q = c(0.5, 0.8), conservative = FALSE)
  expect_s3_class(cf, "sb_confidence")
  expect_true(all(c("area","area_pos","cover") %in% names(cf)))
  # Plot functions should run without error
  png(tempfile(fileext = ".png")); plot(cf, top = 5, label_top = 3); dev.off()
  png(tempfile(fileext = ".png")); plot_stability_curves(g, terms = c("x1","x3"), parameter = "mu"); dev.off()
})
