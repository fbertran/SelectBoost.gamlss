
#' Numerical check: fast vs generic deviance log-likelihood
#'
#' Computes both \code{loglik_gamlss_newdata_fast()} and \code{loglik_gamlss_newdata()}
#' and reports the absolute difference. Useful for sanity checks.
#'
#' @param fit A \code{gamlss} fit.
#' @param newdata Data frame to evaluate on.
#' @param tol Tolerance for pass/fail (default 1e-8).
#' @return A list with fields: \code{ll_fast}, \code{ll_generic}, \code{abs_diff}, \code{pass}.
#' @export
check_fast_vs_generic <- function(fit, newdata, tol = 1e-8) {
  ll_fast <- loglik_gamlss_newdata_fast(fit, newdata)
  ll_gen  <- loglik_gamlss_newdata(fit, newdata)
  abs_diff <- abs(ll_fast - ll_gen)
  list(ll_fast = ll_fast, ll_generic = ll_gen, abs_diff = abs_diff, pass = isTRUE(abs_diff <= tol))
}
