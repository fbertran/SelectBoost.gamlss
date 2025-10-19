test_that("fast vs generic deviance agrees within tolerance on common families", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
  fams <- list(
    NO = gamlss.dist::NO(),
    PO = gamlss.dist::PO(),
    LOGNO = gamlss.dist::LOGNO(),
    GA = gamlss.dist::GA(),
    IG = gamlss.dist::IG(),
    LO = gamlss.dist::LO(),
    LOGITNO = gamlss.dist::LOGITNO(),
    GEOM = gamlss.dist::GEOM(),
    BE = gamlss.dist::BE()
  )
  gens <- list(
    NO = function(n) gamlss.dist::rNO(n, mu = 0, sigma = 1),
    PO = function(n) gamlss.dist::rPO(n, mu = 2.5),
    LOGNO = function(n) gamlss.dist::rLOGNO(n, mu = 0, sigma = 0.6),
    GA = function(n) gamlss.dist::rGA(n, mu = 2, sigma = 0.5),
    IG = function(n) gamlss.dist::rIG(n, mu = 2, sigma = 0.5),
    LO = function(n) gamlss.dist::rLO(n, mu = 0, sigma = 1),
    LOGITNO = function(n) gamlss.dist::rLOGITNO(n, mu = 0, sigma = 1),
    GEOM = function(n) gamlss.dist::rGEOM(n, mu = 3),
    BE = function(n) gamlss.dist::rBE(n, mu = 0.4, sigma = 0.2)
  )
  set.seed(123)
  n <- 1000
  tols <- .family_tolerance()
  tol_default <- tols[['.default']]
  for (nm in names(fams)) {
    fam <- fams[[nm]]
    gen <- gens[[nm]]
    y <- gen(n)
    dat <- data.frame(y = y)
    fit <- try(gamlss::gamlss(y ~ 1, data = dat, family = fam), silent = TRUE)
    if (inherits(fit, "try-error")) next
    tol_fam <- tols[[nm]] %||% tol_default
    chk <- check_fast_vs_generic(fit, dat, tol = tol_fam)
    expect_true(chk$pass, info = paste("Family", nm, "abs diff", chk$abs_diff))
  }
})
