
#' Confidence functionals from a c0 grid
#'
#' Summarize selection proportions across c0 (SelectBoost threshold) into
#' single-number confidence scores per term/parameter.
#'
#' @param x An object from \code{sb_gamlss_c0_grid()}.
#' @param pi_thr Stability threshold; defaults to \code{x$pi_thr}.
#' @param q Numeric vector of quantiles to compute (in 0..1).
#' @param weight_fun Optional function w(c0) for weighted AUSC; default uniform.
#' @param conservative If TRUE, use Wilson lower confidence bounds for proportions.
#' @param B Number of bootstraps (if not inferable when conservative = TRUE).
#' @param method Integration method: "trapezoid" (default) or "step".
#'
#' @return A data.frame with per-term summaries, classed as \code{"sb_confidence"}.
#' @export
confidence_functionals <- function(x, pi_thr = NULL,
                                   q = c(0.5, 0.8, 0.9),
                                   weight_fun = NULL,
                                   conservative = FALSE,
                                   B = NULL,
                                   method = c("trapezoid","step")) {
  stopifnot(inherits(x, "SelectBoost_gamlss_grid"))
  tab <- x$table
  if (is.null(tab) || !nrow(tab)) return(tab)
  method <- match.arg(method)
  if (is.null(pi_thr)) pi_thr <- x$pi_thr

  # Wilson lower bound for binomial proportion
  wilson_lower <- function(p_hat, B, z = 1.96) {
    if (is.null(B) || is.na(B) || B <= 0) return(p_hat)
    n <- B
    denom <- 1 + z^2/n
    center <- (p_hat + z^2/(2*n)) / denom
    half <- z * sqrt((p_hat*(1-p_hat)/n) + z^2/(4*n^2)) / denom
    pmax(0, center - half)
  }

  prop_col <- "prop"
  if (isTRUE(conservative)) {
    if (is.null(B)) {
      steps <- diff(sort(unique(tab$prop)))
      eps <- suppressWarnings(min(steps[steps > 0]))
      B <- if (is.finite(eps) && eps > 0) round(1/eps) else NA_integer_
    }
    tab$prop_lb <- wilson_lower(tab$prop, B = B)
    prop_col <- "prop_lb"
  }

  if (is.null(weight_fun)) {
    tab$w <- 1
  } else {
    tab$w <- weight_fun(tab$c0)
  }

  split_key <- interaction(tab$parameter, tab$term, drop = TRUE, lex.order = TRUE)
  spl <- split(tab, split_key)

  integrate_trap <- function(x, y) {
    if (length(x) < 2) return(mean(y))
    num <- sum( (head(y,-1) + tail(y,-1)) / 2 * diff(x) )
    den <- max(x) - min(x)
    if (isTRUE(all.equal(den, 0))) mean(y) else num / den
  }
  integrate_step <- function(x, y) {
    # left Riemann sum
    if (length(x) < 2) return(mean(y))
    num <- sum( head(y,-1) * diff(x) )
    den <- max(x) - min(x)
    if (isTRUE(all.equal(den, 0))) mean(y) else num / den
  }
  integ <- switch(method, trapezoid = integrate_trap, step = integrate_step)

  one <- function(df) {
    o   <- order(df$c0)
    c0  <- as.numeric(df$c0[o])
    p   <- as.numeric(df[[prop_col]][o])
    w   <- as.numeric(df$w[o])

    # AUSC (normalized)
    area <- integ(c0, p)

    # Thresholded positive area
    pos <- p - pi_thr; pos[pos < 0] <- 0
    area_pos <- integ(c0, pos)

    # Weighted AUSC
    # normalize by integral of weights
    if (length(c0) >= 2) {
      num <- sum( (head(p*w,-1) + tail(p*w,-1)) / 2 * diff(c0) )
      den <- sum( (head(w,-1)   + tail(w,-1))   / 2 * diff(c0) )
      w_area <- ifelse(den > 0, num / den, area)
    } else {
      w_area <- area
    }

    cover <- mean(p >= pi_thr)
    qs    <- stats::quantile(p, probs = q, names = FALSE, type = 7)
    names(qs) <- paste0("q", sprintf("%02d", as.integer(100*q)))

    sup <- max(p); inf <- min(p)

    c(area = area, area_pos = area_pos, w_area = w_area,
      cover = cover, sup = sup, inf = inf, qs)
  }

  mat <- lapply(spl, one)
  out <- as.data.frame(do.call(rbind, mat), stringsAsFactors = FALSE)
  key <- names(spl)
  par <- sub("\\..*$", "", key)
  term <- sub("^[^\\.]*\\.", "", key)

  out$parameter <- par
  out$term      <- term
  rownames(out) <- NULL

  # Composite rank score (tunable weights)
  qcols <- grep("^q\\d+$", names(out), value = TRUE)
  if (length(qcols)) {
    # pick highest quantile present
    hiq <- qcols[which.max(as.integer(sub("^q", "", qcols)))]
    out$rank_score <- 0.5*out$area_pos + 0.3*out$cover + 0.2*out[[hiq]]
  } else {
    out$rank_score <- 0.5*out$area_pos + 0.5*out$cover
  }

  attr(out, "pi_thr") <- pi_thr
  attr(out, "q")      <- q
  attr(out, "method") <- method
  class(out) <- c("sb_confidence", class(out))
  out[order(out$rank_score, decreasing = TRUE), ]
}

#' Plot confidence functionals
#'
#' Two-panel plot: (1) scatter of area_pos vs cover (size by rank), (2) barplot of top-N rank_score.
#' @param x An object from \code{confidence_functionals()}.
#' @param top Show top-N terms in the barplot (default 15).
#' @param label_top Integer; number of points to label in the scatter (default 10).
#' @param ... Graphical parameters passed to plotting backend.
#' @return An invisible copy of `x`.
#' @export
#' @method plot sb_confidence
plot.sb_confidence <- function(x, top = 15, label_top = 10, ...) {
  if (is.null(x) || !nrow(x)) { plot.new(); title("No confidence summaries"); return(invisible()) }
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(mfrow = c(2,1), mar = c(4,4,2,1))

  # Panel 1: scatter of area_pos vs cover
  sz <- 1 + 4 * (x$rank_score - min(x$rank_score)) / max(1e-9, max(x$rank_score) - min(x$rank_score))
  plot(x$area_pos, x$cover, pch = 19, cex = sz,
       xlab = "Thresholded area (area_pos)", ylab = "Coverage (prop >= pi_thr)",
       main = "Confidence landscape: area_pos vs coverage")
  lab_idx <- seq_len(min(label_top, nrow(x)))
  ord <- order(x$rank_score, decreasing = TRUE)
  lab_idx <- ord[lab_idx]
  labs <- paste(x$parameter[lab_idx], x$term[lab_idx], sep = "::")
  text(x$area_pos[lab_idx], x$cover[lab_idx], labels = labs, pos = 4, cex = 0.8)

  # Panel 2: top-N barplot of rank score
  xt <- x[ord[seq_len(min(top, nrow(x)))], ]
  barlab <- paste(xt$parameter, xt$term, sep = "::")
  barplot(xt$rank_score, names.arg = barlab, las = 2,
          ylab = "rank_score", main = "Top terms by rank_score")
  invisible(x)
}

#' Plot stability curves p(c0) for selected terms
#' @param grid An object from \code{sb_gamlss_c0_grid()}.
#' @param terms Character vector of term names to plot.
#' @param parameter Optional parameter name ('mu','sigma','nu','tau'); if NULL, all.
#' @param ncol Columns in the multi-panel layout.
#' @export
plot_stability_curves <- function(grid, terms, parameter = NULL, ncol = 2L) {
  stopifnot(inherits(grid, "SelectBoost_gamlss_grid"))
  tab <- grid$table
  
  # Filter by parameter (if provided) and by terms — avoid NSE / subset()
  if (!is.null(parameter)) {
    tab <- tab[tab[["parameter"]] %in% parameter, , drop = FALSE]
  }
  tab <- tab[tab[["term"]] %in% terms, , drop = FALSE]
  
  if (!NROW(tab)) {
    graphics::plot.new(); graphics::title("No matching terms"); return(invisible(grid))
  }
  
  terms_u <- unique(tab[["term"]])
  n       <- length(terms_u)
  nrowp   <- ceiling(n / ncol)
  
  op <- graphics::par(no.readonly = TRUE); on.exit(graphics::par(op), add = TRUE)
  graphics::par(mfrow = c(nrowp, ncol), mar = c(4, 4, 2, 1))
  
  for (t in terms_u) {
    df <- tab[tab[["term"]] == t, , drop = FALSE]
    o  <- order(df[["c0"]])
    graphics::plot(df[["c0"]][o], df[["prop"]][o],
                   type = "b",
                   xlab = "c0 (SelectBoost threshold)",
                   ylab = "selection proportion",
                   main = paste0(df[["parameter"]][o][1], "::", t))
    graphics::abline(h = grid$pi_thr, lty = 2)
  }
  
  invisible(grid)
}
