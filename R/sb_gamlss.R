#' SelectBoost for GAMLSS (stability selection)
#'
#' @param formula Base formula for the location \eqn{\mu} parameter (in the main model call).
#' @param data Data frame.
#' @param family A `gamlss.dist` family object (e.g., `gamlss.dist::NO()`).
#' @param mu_scope Formula of candidate terms for \eqn{\mu}.
#' @param sigma_scope,nu_scope,tau_scope Formulas of candidate terms for \eqn{\sigma}, \eqn{\nu}, \eqn{\tau}.
#' @param base_sigma,base_nu,base_tau Optional base (always-included) formulas for \eqn{\sigma}, \eqn{\nu}, \eqn{\tau}.
#' @param B Number of bootstrap subsamples for stability selection.
#' @param sample_fraction Fraction of rows per subsample (e.g., 0.7).
#' @param pi_thr Selection proportion threshold to define “stable” terms (e.g., 0.6).
#' @param c0 SelectBoost meta-parameter controlling reweighting/thresholding (see vignette).
#' @param k Penalty weight for stepwise GAIC when `engine = "stepGAIC"` (default 2).
#' @param direction Stepwise direction for `stepGAIC` (`"both"`, `"forward"`, `"backward"`).
#' @param engine Engine for \eqn{\mu} (`"stepGAIC"`, `"glmnet"`, `"grpreg"`, `"sgl"`).
#' @param engine_sigma,engine_nu,engine_tau Optional engines for \eqn{\sigma}, \eqn{\nu}, \eqn{\tau}.
#' @param glmnet_alpha Elastic-net mixing for glmnet (1 = lasso, 0 = ridge).
#' @param grpreg_penalty Group penalty for grpreg (`"grLasso"`, `"grMCP"`, `"grSCAD"`).
#' @param sgl_alpha Alpha for sparse group lasso.
#' @param pre_standardize Logical; standardize numeric predictors before penalized fits.
#' @param use_groups Logical; treat all columns from the same term (factor dummies, spline bases) as one group.
#' @param df_smooth Degrees of freedom for proxy spline bases (`pb()/cs()` → `splines::bs(df=df_smooth)`) used only for grouped selection design.
#' @param parallel Parallel mode (`"none"`, `"auto"`, `"multisession"`, `"multicore"`).
#' @param workers Integer; number of workers if parallel.
#' @param progress Logical; show a progress bar in sequential runs.
#' @param trace Logical; print progress messages.
#' @param ... Passed to underlying engines (e.g., to `gamlss::gamlss`, `glmnet`, etc.).
#'
#' @return An object of class \code{"sb_gamlss"} with elements:
#' \itemize{
#'   \item \code{final_fit}: the final \code{gamlss} object.
#'   \item \code{final_formula}: list of formulas for mu/sigma/nu/tau.
#'   \item \code{selection}: data.frame of selection counts and proportions.
#'   \item \code{B}, \code{sample_fraction}, \code{pi_thr}, \code{k}.
#'   \item \code{scaler}: list with \code{center}, \code{scale}, \code{vars}, \code{response}.
#' }
#' @examplesIf requireNamespace("gamlss.dist", quietly = TRUE)
#' set.seed(1)
#' dat <- data.frame(
#'   y = gamlss.dist::rNO(60, mu = 0),
#'   x1 = rnorm(60),
#'   x2 = rnorm(60),
#'   x3 = rnorm(60)
#' )
#' fit <- sb_gamlss(
#'   y ~ 1,
#'   data = dat,
#'   family = gamlss.dist::NO(),
#'   mu_scope = ~ x1 + x2 + gamlss::pb(x3),
#'   B = 8,
#'   pi_thr = 0.6,
#'   trace = FALSE
#' )
#' fit$final_formula
#' @export
sb_gamlss <- function(
  formula,
  data,
  family,
  mu_scope,
  sigma_scope = NULL,
  nu_scope = NULL,
  tau_scope = NULL,
  base_sigma = ~ 1,
  base_nu = ~ 1,
  base_tau = ~ 1,
  B = 100,
  sample_fraction = 0.7,
  pi_thr = 0.6,
  k = 2,
  direction = c("both","forward","backward"),
  pre_standardize = FALSE,
  use_groups = FALSE,
  c0 = 0.5,
  engine = c('stepGAIC','glmnet','grpreg','sgl'),
  engine_sigma = NULL,
  engine_nu = NULL,
  engine_tau = NULL,
  grpreg_penalty = c('grLasso','grMCP','grSCAD'),
  sgl_alpha = 0.95,
  df_smooth = 6L,
  progress = TRUE,
  glmnet_alpha = 1,
  parallel = c('none','auto','multisession','multicore'),
  workers = NULL,
  trace = TRUE,
  ...
) {
  stopifnot(inherits(formula, "formula"))
  direction <- match.arg(direction)
  engine <- match.arg(engine)
  parallel <- match.arg(parallel)
  grpreg_penalty <- match.arg(grpreg_penalty)
  engine_sigma <- if (is.null(engine_sigma)) 'stepGAIC' else match.arg(engine_sigma, c('stepGAIC','glmnet','grpreg','sgl'))
  engine_nu    <- if (is.null(engine_nu))    'stepGAIC' else match.arg(engine_nu,    c('stepGAIC','glmnet','grpreg','sgl'))
  engine_tau   <- if (is.null(engine_tau))   'stepGAIC' else match.arg(engine_tau,   c('stepGAIC','glmnet','grpreg','sgl'))

  if (missing(mu_scope) || !inherits(mu_scope, "formula"))
    stop("Provide `mu_scope` as an RHS-only formula, e.g., ~ x1 + x2")

  if (!is.data.frame(data)) data <- as.data.frame(data)
  n <- NROW(data)
  if (n <= 1) stop("Data must have at least 2 rows.")
  if (sample_fraction <= 0 || sample_fraction > 1) stop("sample_fraction must be in (0,1].")

  # Helpers
  term_labels <- function(f) {
    tl <- attr(stats::terms(f), "term.labels")
    if (is.null(tl)) character(0) else tl
  }

  add_terms_formula <- function(base_formula, new_terms = character(0)) {
    trm <- stats::terms(base_formula)
    intercept <- isTRUE(attr(trm, "intercept") == 1L)
    has_response <- isTRUE(attr(trm, "response") == 1L)
    lhs <- if (has_response) paste(deparse(base_formula[[2L]]), collapse = "") else NULL
    rhs_terms <- term_labels(base_formula)
    rhs <- unique(c(rhs_terms, new_terms))
    out <- stats::reformulate(rhs, response = lhs, intercept = intercept)
    environment(out) <- environment(base_formula)
    out
  }
  
  resp_name <- all.vars(formula)[1L]

  # Pre-standardization of numeric predictors (store scaler for later)
  scaler <- NULL
  # Use RcppArmadillo scaler if requested
  if (isTRUE(pre_standardize)) {
    num_vars <- names(which(sapply(data, is.numeric)))
    num_vars <- setdiff(num_vars, resp_name)
    if (length(num_vars)) {
      cs <- cs_scale_mat(data, num_vars)
      data[num_vars] <- as.data.frame(cs$X)
      scaler <- list(center = setNames(cs$center, num_vars), scale = setNames(cs$scale, num_vars), vars = num_vars, response = resp_name)
    }
  }
  # Build correlation-based groups (SelectBoost) for numeric predictors if requested
  groups <- NULL
  if (isTRUE(use_groups)) {
    num_vars_all <- names(which(sapply(data, is.numeric)))
    num_vars_all <- setdiff(num_vars_all, resp_name)
    if (length(num_vars_all) >= 2) {
      cormat <- try(fast_cor(stats::na.omit(data), num_vars_all), silent = TRUE)
      if (!inherits(cormat, "try-error")) {
        groups <- try(SelectBoost::group_func_2(cormat, c0), silent = TRUE)
        if (inherits(groups, "try-error")) groups <- NULL
      }
    }
  }

  if (isTRUE(pre_standardize)) {
    num_vars <- names(which(sapply(data, is.numeric)))
    num_vars <- setdiff(num_vars, resp_name)
    if (length(num_vars)) {
      center <- vapply(data[num_vars], function(x) mean(x, na.rm = TRUE), numeric(1))
      scale  <- vapply(data[num_vars], function(x) stats::sd(x, na.rm = TRUE), numeric(1))
      # avoid zero division
      scale[scale == 0] <- 1
      data[num_vars] <- Map(function(x, m, s) (x - m)/s, data[num_vars], center, scale)
      scaler <- list(center = center, scale = scale, vars = num_vars, response = resp_name)
    }
  }

  base_mu    <- formula
  base_sigma <- if (inherits(base_sigma, "formula")) base_sigma else ~ 1
  base_nu    <- if (inherits(base_nu, "formula")) base_nu else ~ 1
  base_tau   <- if (inherits(base_tau, "formula")) base_tau else ~ 1

  base_mu_norm    <- add_terms_formula(base_mu)
  base_sigma_norm <- add_terms_formula(base_sigma)
  base_nu_norm    <- add_terms_formula(base_nu)
  base_tau_norm   <- add_terms_formula(base_tau)
  
  build_upper <- function(base_rhs, scope_rhs) {
    base_rhs <- add_terms_formula(base_rhs)
    if (is.null(scope_rhs)) return(base_rhs)
    cand <- term_labels(scope_rhs)
    if (length(cand) == 0) return(base_rhs)
    add_terms_formula(base_rhs, cand)
  }

  upper_mu    <- build_upper(base_mu_norm,    mu_scope)
  upper_sigma <- build_upper(base_sigma_norm, sigma_scope)
  upper_nu    <- build_upper(base_nu_norm,    nu_scope)
  upper_tau   <- build_upper(base_tau_norm,   tau_scope)
  
  seen_terms <- list(
    mu    = unique(term_labels(upper_mu)),
    sigma = unique(term_labels(upper_sigma)),
    nu    = unique(term_labels(upper_nu)),
    tau   = unique(term_labels(upper_tau))
  )
  counts <- lapply(seen_terms, function(x) setNames(integer(length(x)), x))

  if (trace) message("Bootstrapping ", B, " replicates...")

  iter_fun <- function(b) {
    # If grouping is enabled, pick one representative per group for this bootstrap
    upper_mu_b    <- upper_mu
    upper_sigma_b <- upper_sigma
    upper_nu_b    <- upper_nu
    upper_tau_b   <- upper_tau
    if (!is.null(groups)) {
      pick_one <- function(terms) {
        # Terms are strings; keep one per group if term name matches variable name as prefix
        if (length(terms) == 0) return(terms)
        chosen <- character(0)
        # groups is a list of integer vectors (indices of variables), with attr 'type' etc.
        # We use variable names from num_vars_all; map indices to names
        vars <- num_vars_all
        for (g in groups) {
          cand_vars <- vars[g]
          present <- cand_vars[cand_vars %in% terms]
          if (length(present)) {
            chosen <- c(chosen, sample(present, 1))
          }
        }
        # Also keep any terms not in numeric groups (e.g., factors/smoothers)
        others <- setdiff(terms, vars)
        unique(c(chosen, others))
      }
      # Filter the RHS term lists of each parameter's upper scope
      tl_mu    <- attr(stats::terms(upper_mu), "term.labels");    tl_mu    <- pick_one(tl_mu)
      tl_sigma <- attr(stats::terms(upper_sigma), "term.labels"); tl_sigma <- pick_one(tl_sigma)
      tl_nu    <- attr(stats::terms(upper_nu), "term.labels");    tl_nu    <- pick_one(tl_nu)
      tl_tau   <- attr(stats::terms(upper_tau), "term.labels");   tl_tau   <- pick_one(tl_tau)

      rebuild <- function(base_rhs, terms) {
        if (length(terms) == 0) return(add_terms_formula(base_rhs))
        add_terms_formula(base_rhs, terms)
      }
      upper_mu_b    <- rebuild(base_mu,    tl_mu)
      upper_sigma_b <- rebuild(base_sigma, tl_sigma)
      upper_nu_b    <- rebuild(base_nu,    tl_nu)
      upper_tau_b   <- rebuild(base_tau,   tl_tau)
    }
    idx <- sample.int(n, size = max(2L, floor(sample_fraction * n)), replace = FALSE)
    dat_b <- data[idx, , drop = FALSE]

    fit0 <- gamlss::gamlss(
      formula = base_mu,
      sigma.formula = base_sigma,
      nu.formula = base_nu,
      tau.formula = base_tau,
      data = dat_b,
      family = family,
      ...
    )

    try_step <- function(fit, lower, upper, what) {
      if (identical(term_labels(lower), term_labels(upper))) return(fit)
      out <- try(
        gamlss::stepGAIC(
          fit,
          scope = list(lower = lower, upper = upper),
          what = what,
          direction = direction,
          k = k,
          trace = FALSE
        ),
        silent = TRUE
      )
      if (inherits(out, "try-error")) fit else out
    }

    fit1 <- try_step(fit0, base_mu_norm,    upper_mu_b,    what = "mu")
    fit1 <- try_step(fit1, base_sigma_norm, upper_sigma_b, what = "sigma")
    fit1 <- try_step(fit1, base_nu_norm,    upper_nu_b,    what = "nu")
    fit1 <- try_step(fit1, base_tau_norm,   upper_tau_b,   what = "tau")

    get_param_terms <- function(fit, what) {
      f <- try(stats::formula(fit, what = what), silent = TRUE)
      if (inherits(f, "try-error")) return(character(0))
      term_labels(f)
    }

    sel <- list(
      mu    = get_param_terms(fit1, "mu"),
      sigma = get_param_terms(fit1, "sigma"),
      nu    = get_param_terms(fit1, "nu"),
      tau   = get_param_terms(fit1, "tau")
    )

    for (par in names(counts)) {
      if (length(counts[[par]]) == 0L) next
      present <- intersect(names(counts[[par]]), sel[[par]])
      if (length(present)) counts[[par]][present] <- counts[[par]][present] + 1L
    }

    if (trace && (b %% max(1L, floor(B/10))) == 0) message("  replicate ", b, "/", B)
  }

  mk_tbl <- function(cnts, par) {
    if (length(cnts) == 0) return(NULL)
    data.frame(
      parameter = par,
      term = names(cnts),
      count = as.integer(unname(cnts)),
      prop = as.numeric(unname(cnts))/B,
      stringsAsFactors = FALSE
    )
  }
  tabs <- do.call(rbind, Filter(Negate(is.null), list(
    mk_tbl(counts$mu, "mu"),
    mk_tbl(counts$sigma, "sigma"),
    mk_tbl(counts$nu, "nu"),
    mk_tbl(counts$tau, "tau")
  )))

  keep_terms <- function(cnts) {
    if (length(cnts) == 0) return(character(0))
    names(cnts)[(cnts/B) >= pi_thr]
  }

  mu_keep    <- keep_terms(counts$mu)
  sigma_keep <- keep_terms(counts$sigma)
  nu_keep    <- keep_terms(counts$nu)
  tau_keep   <- keep_terms(counts$tau)

  add_keep <- function(base_rhs, keep) {
    if (length(keep) == 0) return(add_terms_formula(base_rhs))
    add_terms_formula(base_rhs, keep)
  }

  final_mu    <- add_keep(base_mu,    mu_keep)
  final_sigma <- add_keep(base_sigma, sigma_keep)
  final_nu    <- add_keep(base_nu,    nu_keep)
  final_tau   <- add_keep(base_tau,   tau_keep)

  final_fit <- gamlss::gamlss(
    formula = final_mu,
    sigma.formula = final_sigma,
    nu.formula = final_nu,
    tau.formula = final_tau,
    data = data,
    family = family,
    ...
  )

  out <- list(
    final_fit = final_fit,
    final_formula = list(mu = final_mu, sigma = final_sigma, nu = final_nu, tau = final_tau),
    selection = tabs,
    B = B,
    sample_fraction = sample_fraction,
    pi_thr = pi_thr,
    k = k,
    scaler = scaler,
    call = match.call()
  )
  class(out) <- "sb_gamlss"
  out
}

#' AICc for a gamlss fit
#' @param object a 'gamlss' object
#' @return numeric AICc value
#' @export
AICc_gamlss <- function(object) {
  ll <- as.numeric(logLik(object))
  k  <- attr(logLik(object), "df")
  n  <- NROW(object$y)
  aic <- -2*ll + 2*k
  aicc <- aic + (2*k*(k+1))/(n - k - 1)
  aicc
}

#' Selection table accessor
#' @param x An sb_gamlss object
#' @return data.frame with parameter, term, count, prop
#' @export
selection_table <- function(x) {
  stopifnot(inherits(x, "sb_gamlss"))
  x$selection
}

#' Plot selection frequencies for sb_gamlss
#' @param x An sb_gamlss object
#' @param top Show only the top N terms per-parameter (default all)
#' @param ... Graphical parameters.
#' @export
plot_sb_gamlss <- function(x, top = Inf, ...) {
  stopifnot(inherits(x, "sb_gamlss"))
  tab <- x$selection
  if (is.null(tab) || NROW(tab) == 0L) {
    plot.new(); title("No selection results to plot"); return(invisible())
  }
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  pars <- unique(tab$parameter)
  nP <- length(pars)
  par(mfrow = c(nP, 1), mar = c(5, 8, 2, 1))
  for (p in pars) {
    sub <- tab[tab$parameter == p, , drop = FALSE]
    ord <- order(sub$prop, decreasing = TRUE)
    sub <- sub[ord, , drop = FALSE]
    if (is.finite(top) && top < nrow(sub)) sub <- utils::head(sub, top)
    barplot(height = sub$prop, names.arg = sub$term, horiz = TRUE, las = 1,
            xlab = "Selection proportion", main = paste0("Parameter: ", p), ...)
    if (!is.null(x$pi_thr)) abline(v = x$pi_thr, lty = 2)
  }
  invisible(x)
}

#' @export
coef.sb_gamlss <- function(object, ...) stats::coef(object$final_fit, ...)

#' @export
predict.sb_gamlss <- function(object, ...) {
  dots <- list(...)
  if (!is.null(object$scaler) && !is.null(dots$newdata)) {
    nd <- dots$newdata
    vars <- object$scaler$vars
    if (length(vars)) {
      for (v in vars) if (v %in% names(nd)) {
        m <- object$scaler$center[[v]]; s <- object$scaler$scale[[v]]
        if (is.null(m) || is.null(s) || is.na(s) || s == 0) s <- 1
        nd[[v]] <- (nd[[v]] - m)/s
      }
      dots$newdata <- nd
    }
  }
  do.call(stats::predict, c(list(object$final_fit), dots))
}
