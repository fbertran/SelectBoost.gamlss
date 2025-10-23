
#' Stability curves over a c0 grid for sb_gamlss
#' @inheritParams sb_gamlss
#' @param c0_grid Numeric vector of `c0` thresholds in (0,1).
#' @param progress Logical; show a progress bar across `c0_grid`.
#' @return An object of class "SelectBoost_gamlss_grid" containing
#' \itemize{
#'   \item \code{results}: named list of sb_gamlss fits, names are c0 values
#'   \item \code{table}: data.frame with parameter, term, count, prop, c0
#'   \item \code{pi_thr}: the threshold used
#' }
#' @export
#' 
sb_gamlss_c0_grid <- function(
  formula, data, family,
  mu_scope, sigma_scope = NULL, nu_scope = NULL, tau_scope = NULL,
  base_sigma = ~ 1, base_nu = ~ 1, base_tau = ~ 1,
  c0_grid = seq(0.1, 0.9, by = 0.1),
  B = 60, sample_fraction = 0.7, pi_thr = 0.6, k = 2,
  direction = c("both","forward","backward"),
  pre_standardize = FALSE, trace = TRUE, progress = TRUE,
  use_groups = TRUE, corr_func = "cor", group_fun = SelectBoost::group_func_2, ...
) {
  direction <- match.arg(direction)
  stopifnot(is.numeric(c0_grid), all(c0_grid > 0), all(c0_grid < 1))
  res_list <- list()
  pb <- NULL
  if (isTRUE(progress)) pb <- utils::txtProgressBar(min = 0, max = length(c0_grid), style = 3)
  i_prog <- 0
  tabs <- list()
  for (c0 in c0_grid) {
    fit <- sb_gamlss(
      formula = formula, data = data, family = family,
      mu_scope = mu_scope, sigma_scope = sigma_scope, nu_scope = nu_scope, tau_scope = tau_scope,
      base_sigma = base_sigma, base_nu = base_nu, base_tau = base_tau,
      B = B, sample_fraction = sample_fraction, pi_thr = pi_thr, k = k,
      direction = direction, pre_standardize = pre_standardize,
      use_groups = use_groups, c0 = c0, trace = trace,
      corr_func = corr_func, group_fun = group_fun, ...
    )
    res_list[[as.character(c0)]] <- fit
    if (!is.null(pb)) { i_prog <- i_prog + 1; utils::setTxtProgressBar(pb, i_prog) }
    tab <- selection_table(fit)
    if (!is.null(tab) && NROW(tab)) {
      tab$c0 <- c0
      tabs[[length(tabs)+1L]] <- tab
    }
  }
  if (!is.null(pb)) close(pb)
  table <- if (length(tabs)) do.call(rbind, tabs) else
    data.frame(parameter=character(), term=character(), count=integer(), prop=double(), c0=double())
  out <- list(results = res_list, table = table, pi_thr = pi_thr, call = match.call())
  class(out) <- "SelectBoost_gamlss_grid"
  out
}

#' Compute SelectBoost-like confidence table across c0
#'
#' @param grid an object returned by \code{sb_gamlss_c0_grid}
#' @param pi_thr optional override of the threshold (defaults to grid$pi_thr)
#' @return data.frame with term, parameter, conf_index (mean positive excess), cover (fraction of c0 with prop>=thr)
#' @export
confidence_table <- function(grid, pi_thr = NULL) {
  stopifnot(inherits(grid, "SelectBoost_gamlss_grid"))
  tab <- grid$table
  if (is.null(tab) || !NROW(tab)) return(tab)
  thr <- if (is.null(pi_thr)) grid$pi_thr else pi_thr
  agg <- by(tab, list(tab$parameter, tab$term), function(df) {
    # mean positive part of (prop - thr)
    conf <- mean(pmax(0, df$prop - thr))
    cover <- mean(as.numeric(df$prop >= thr))
    data.frame(parameter=df$parameter[1], term=df$term[1], conf_index=conf, cover=cover)
  })
  out <- do.call(rbind, as.list(agg))
  rownames(out) <- NULL
  out[order(out$conf_index, decreasing = TRUE), , drop = FALSE]
}


#' Plot summary for sb_gamlss_c0_grid
#' @param x A `SelectBoost_gamlss_grid` object.
#' @param top Integer; how many top terms to show in the confidence barplot.
#' @param ... Ignored (reserved for future).
#' @return An invisible copy of `x`.
#' @export
#' @method plot SelectBoost_gamlss_grid 
plot.SelectBoost_gamlss_grid <- function(x, top = 15, ...) {
  tab <- x$table
  if (is.null(tab) || !NROW(tab)) { graphics::plot.new(); graphics::title("No results"); return(invisible(x)) }
  
  ## Build a clean data.frame with explicit columns (no NSE / bare symbols)
  ok  <- as.integer(tab[["prop"]] >= x$pi_thr)
  tmp <- data.frame(c0 = tab[["c0"]], parameter = tab[["parameter"]], ok = ok)
  
  ## Panel 1: number of stable terms vs c0 (summed across parameters)
  stable_by_c0 <- stats::aggregate(ok ~ c0 + parameter, data = tmp, FUN = sum)
  
  ## Panel 2: top terms by confidence
  conf <- confidence_table(x)
  if (is.finite(top) && top < nrow(conf)) conf <- utils::head(conf, top)
  
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(4, 4, 2, 1))
  
  ## Plot 1
  c0s <- sort(unique(stable_by_c0$c0))
  total_stable <- tapply(stable_by_c0$ok, stable_by_c0$c0, sum)
  graphics::plot(c0s,
                 as.numeric(total_stable[as.character(c0s)]),
                 type = "b", xlab = "c0", ylab = "# stable terms",
                 main = "Stable terms vs c0")
  graphics::abline(h = 0, lty = 2)
  
  ## Plot 2
  graphics::barplot(height = conf$conf_index,
                    names.arg = paste(conf$parameter, conf$term, sep = "::"),
                    las = 2,
                    main = "Top terms by confidence (across c0)",
                    ylab = "conf_index")
  
  invisible(x)
}


#' AutoBoost for GAMLSS (SelectBoost-style)
#'
#' Runs a c0 grid, picks the c0 that maximizes total confidence, and returns the corresponding sb_gamlss fit.
#' 
#' @inheritParams sb_gamlss
#' @param c0_grid Numeric vector of `c0` values.
#' @return A `SelectBoost_gamlss_grid` with summary plots/tables.
#' @export
autoboost_gamlss <- function(
  formula, data, family,
  mu_scope, sigma_scope = NULL, nu_scope = NULL, tau_scope = NULL,
  base_sigma = ~ 1, base_nu = ~ 1, base_tau = ~ 1,
  c0_grid = seq(0.1, 0.9, by = 0.1),
  B = 60, sample_fraction = 0.7, pi_thr = 0.6, k = 2,
  direction = c("both","forward","backward"),
  pre_standardize = FALSE, trace = TRUE, progress = TRUE,
  use_groups = TRUE, corr_func = "cor", group_fun = SelectBoost::group_func_2, ...
) {
  direction <- match.arg(direction)
  grid <- sb_gamlss_c0_grid(
    formula = formula, data = data, family = family,
    mu_scope = mu_scope, sigma_scope = sigma_scope, nu_scope = nu_scope, tau_scope = tau_scope,
    base_sigma = base_sigma, base_nu = base_nu, base_tau = base_tau,
    c0_grid = c0_grid, B = B, sample_fraction = sample_fraction, pi_thr = pi_thr, k = k,
    direction = direction, pre_standardize = pre_standardize, trace = trace,
    use_groups = use_groups, corr_func = corr_func, group_fun = group_fun, ...
  )
  conf <- confidence_table(grid)
  if (is.null(grid$table) || !NROW(grid$table) || !NROW(conf)) {
    # Fallback: pick median c0
    chosen <- as.character(stats::median(c0_grid))
  } else {
    # For each c0 compute total confidence
    tab <- grid$table
    thr <- grid$pi_thr
    score_by_c0 <- tapply(pmax(0, tab$prop - thr), tab$c0, sum)
    # choose the c0 with the highest score; tie -> middle
    best <- names(score_by_c0)[which.max(score_by_c0)]
    chosen <- if (length(best)) best else as.character(stats::median(c0_grid))
  }
  fit <- grid$results[[chosen]]
  attr(fit, "chosen_c0") <- as.numeric(chosen)
  attr(fit, "c0_grid") <- c0_grid
  attr(fit, "confidence_table") <- conf
  class(fit) <- unique(c("autoboost_gamlss", class(fit)))
  fit
}

#' FastBoost for GAMLSS (lightweight stability selection)
#'
#' A faster variant with fewer bootstraps and smaller subsamples.
#' 
#' Fast SelectBoost (single c0)
#' @inheritParams sb_gamlss
#' @return An `sb_gamlss` fit at the given `c0`.
#' @export
fastboost_gamlss <- function(
  formula, data, family,
  mu_scope, sigma_scope = NULL, nu_scope = NULL, tau_scope = NULL,
  base_sigma = ~ 1, base_nu = ~ 1, base_tau = ~ 1,
  B = 30, sample_fraction = 0.6, pi_thr = 0.6, k = 2,
  direction = c("both","forward","backward"),
  pre_standardize = FALSE, use_groups = TRUE, c0 = 0.5, trace = TRUE,
  corr_func = "cor", group_fun = SelectBoost::group_func_2, ...
) {
  sb_gamlss(
    formula = formula, data = data, family = family,
    mu_scope = mu_scope, sigma_scope = sigma_scope, nu_scope = nu_scope, tau_scope = tau_scope,
    base_sigma = base_sigma, base_nu = base_nu, base_tau = base_tau,
    B = B, sample_fraction = sample_fraction, pi_thr = pi_thr, k = k,
    direction = direction, pre_standardize = pre_standardize,
    use_groups = use_groups, c0 = c0, trace = trace,
    corr_func = corr_func, group_fun = group_fun, ...
  )
}
