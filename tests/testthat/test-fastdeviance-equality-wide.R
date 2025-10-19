test_that("fast vs generic deviance equality holds across a broad family sweep", {
  skip_on_cran()
  skip_if_not_installed("gamlss")
skip_long_tests()
  set.seed(999)
  fams <- c("NO","PO","LOGNO","GA","IG","LO","LOGITNO","GEOM","BE",
            "NBI","NBII","LOGLOG","DEL","ZAGA","ZIP","ZINBI","DPO","GPO",
            "ZAIG","ZALG","BCT","BCPE","ZIPF","ZIPFmu","SICHEL","GLG",
            "BETA4","RS","WEI","GIG")
  n <- 2000
  tols <- .family_tolerance()
  tol_default <- tols[['.default']]
  tol <- tol_default
  # iterate
  for (fam in fams) {
    y <- .gen_family(fam, n)
    if (is.null(y)) next
    dat <- data.frame(y = y)
    # Some families require special forms or additional args; we use default ~1 fits
    fam_obj <- try(getFromNamespace(fam, "gamlss.dist")(), silent = TRUE)
    if (inherits(fam_obj, "try-error")) next
    fit <- try(gamlss::gamlss(y ~ 1, data = dat, family = fam_obj), silent = TRUE)
    if (inherits(fit, "try-error")) next
    tol_fam <- tols[[fam]] %||% tol_default
chk <- check_fast_vs_generic(fit, dat, tol = tol_fam)
    expect_true(chk$pass, info = paste("Family", fam, "abs diff", chk$abs_diff))
  }
})
