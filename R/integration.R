
#' SelectBoost-style wrapper for GAMLSS
#' @description A thin wrapper around \code{sb_gamlss()} with SelectBoost-flavored arguments.
#' @param c0 Correlation threshold for grouping (as in SelectBoost::group_func_2).
#' @param use_groups Logical; enable SelectBoost grouping.
#' @inheritParams sb_gamlss
#' @return An object of class c("SelectBoost_gamlss"), with slots similar to sb_gamlss.
#' @export
SelectBoost_gamlss <- function(
  formula, data, family,
  mu_scope, sigma_scope = NULL, nu_scope = NULL, tau_scope = NULL,
  base_sigma = ~ 1, base_nu = ~ 1, base_tau = ~ 1,
  B = 100, sample_fraction = 0.7, pi_thr = 0.6, k = 2,
  direction = c("both","forward","backward"),
  pre_standardize = FALSE, use_groups = TRUE, c0 = 0.5, trace = TRUE, ...
) {
  res <- sb_gamlss(
    formula = formula, data = data, family = family,
    mu_scope = mu_scope, sigma_scope = sigma_scope, nu_scope = nu_scope, tau_scope = tau_scope,
    base_sigma = base_sigma, base_nu = base_nu, base_tau = base_tau,
    B = B, sample_fraction = sample_fraction, pi_thr = pi_thr, k = k,
    direction = direction, pre_standardize = pre_standardize,
    use_groups = use_groups, c0 = c0, trace = trace, ...
  )
  class(res) <- c("SelectBoost_gamlss", class(res))
  res
}

#' @export
summary.SelectBoost_gamlss <- function(object, prop.level = 0.6, ...) {
  sel <- selection_table(object)
  if (is.null(sel) || NROW(sel) == 0L) return(structure(list(selection = sel), class = "summary.SelectBoost_gamlss"))
  # Compute a simple "confidence index" per-term akin to SelectBoost: proportion minus threshold
  conf <- aggregate(prop ~ term + parameter, data = sel, FUN = function(x) max(0, mean(x) - prop.level))
  names(conf)[names(conf) == "prop"] <- "conf_index"
  out <- list(selection = sel, threshold = prop.level, confidence = conf)
  class(out) <- "summary.SelectBoost_gamlss"
  out
}


#' @export
plot.summary.SelectBoost_gamlss <- function(x, ...) {
  # Simple 2-panel: selection proportions and confidence index
  sel <- x$selection
  conf <- x$confidence
  if (is.null(sel) || NROW(sel) == 0L) { plot.new(); title("No results"); return(invisible()) }
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(2,1), mar = c(5, 8, 2, 1))
  # Panel 1: proportions by parameter
  params <- unique(sel$parameter)
  for (p in params) {
    sub <- sel[sel$parameter == p, , drop = FALSE]
    ord <- order(sub$prop, decreasing = TRUE); sub <- sub[ord, , drop = FALSE]
    barplot(sub$prop, names.arg = sub$term, horiz = TRUE, las = 1,
            xlab = "Selection proportion", main = paste0("Parameter: ", p))
    abline(v = x$threshold %||% 0.6, lty = 2)
  }
  # Panel 2: confidence index
  if (!is.null(conf) && NROW(conf)) {
    conf$lab <- paste(conf$parameter, conf$term, sep = "::")
    ord <- order(conf$conf_index, decreasing = TRUE); conf <- conf[ord, , drop = FALSE]
    barplot(conf$conf_index, names.arg = conf$lab, horiz = TRUE, las = 1,
            xlab = "Confidence index", main = "Confidence index (proxy)")
  }
  invisible(x)
}

#' Plot selection proportions for a single sb_gamlss
#' @param x An `sb_gamlss` object.
#' @param ... Graphical parameters.
#' @export
plot.SelectBoost_gamlss <- function(x, ...) {
  plot_sb_gamlss(x, ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
