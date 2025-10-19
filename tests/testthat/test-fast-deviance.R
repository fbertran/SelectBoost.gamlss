test_that("fast deviance basic path works for NO", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  set.seed(1)
  y <- gamlss.dist::rNO(50, mu = 0, sigma = 1)
  dat <- data.frame(y = y)
  fit <- gamlss::gamlss(y ~ 1, data = dat, family = gamlss.dist::NO())
  # newdata = same data for simplicity
  ll <- loglik_gamlss_newdata_fast(fit, dat)
  expect_true(is.numeric(ll))
})

test_that("native densities exist when available", {
  skip_on_cran()
  fams <- c("LOGLOG","DEL","ZAGA","ZIP2","ZAIG","ZALG","ZIBI","ZIBB","PARETO","SEP1","SEP2",
            "ZIPF","ZIPFmu","BCT","BCPE","SICHEL","GLG","BETA4","RS","WEI","GIG")
  for (f in fams) {
    dens_name <- paste0("d", f)
    fn <- tryCatch(getFromNamespace(dens_name, "gamlss.dist"), error = function(e) NULL)
    # If it exists, it should be a function
    if (!is.null(fn)) expect_true(is.function(fn))
  }
})
