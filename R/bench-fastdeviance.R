
#' Compare fast vs generic deviance log-likelihood evaluation
#'
#' @param fit A \code{gamlss} fit.
#' @param newdata Data frame to evaluate on.
#' @param reps Number of repetitions (default 100).
#' @param unit microbenchmark unit (default "us").
#' @return A data.frame with method, median, and relative speed.
#' @export
fast_vs_generic_ll <- function(fit, newdata, reps = 100L, unit = "us") {
  f_fast <- function() loglik_gamlss_newdata_fast(fit, newdata)
  f_gen  <- function() loglik_gamlss_newdata(fit, newdata)

  if (requireNamespace("microbenchmark", quietly = TRUE)) {
    mb <- microbenchmark::microbenchmark(
      fast    = f_fast(),
      generic = f_gen(),
      times = reps, unit = unit
    )
    sm <- summary(mb)
    out <- data.frame(method = sm$expr, median = sm$median, unit = attr(mb,"unit"), stringsAsFactors = FALSE)
  } else {
    # fallback: crude timing
    tb <- system.time(for (i in seq_len(reps)) f_fast())
    tg <- system.time(for (i in seq_len(reps)) f_gen())
    out <- data.frame(method = c("fast","generic"), median = c(tb[["elapsed"]], tg[["elapsed"]]) / reps, unit = "s/rep", stringsAsFactors = FALSE)
  }
  # relative speed: generic / fast
  if (nrow(out) == 2) {
    rel <- out$median[out$method=="generic"] / out$median[out$method=="fast"]
    out$rel_speed <- c(rel, 1/rel)[match(out$method, c("fast","generic"))]
  }
  out
}
