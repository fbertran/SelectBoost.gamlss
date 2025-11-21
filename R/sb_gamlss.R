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
#' @param glmnet_family Family passed to glmnet-based selectors ("gaussian", "binomial", "poisson").
#' @param grpreg_penalty Group penalty for grpreg (`"grLasso"`, `"grMCP"`, `"grSCAD"`).
#' @param sgl_alpha Alpha for sparse group lasso.
#' @param pre_standardize Logical; standardize numeric predictors before penalized fits.
#' @param use_groups Logical; treat SelectBoost correlation groups during resampling.
#' @param df_smooth Degrees of freedom for proxy spline bases (`pb()/cs()` → `splines::bs(df=df_smooth)`) used only for grouped selection design.
#' @param parallel Parallel mode (`"none"`, `"auto"`, `"multisession"`, `"multicore").
#' @param workers Integer; number of workers if parallel.
#' @param progress Logical; show a progress bar in sequential runs.
#' @param trace Logical; print progress messages.
#' @param corr_func Correlation function passed to `SelectBoost::boost.compcorrs`.
#' @param group_fun Grouping function passed to `SelectBoost::boost.findgroups`.
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
  glmnet_family = c("gaussian","binomial","poisson"),
  parallel = c('none','auto','multisession','multicore'),
  workers = NULL,
  trace = TRUE,
  corr_func = "cor",
  group_fun = SelectBoost::group_func_2,
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
  if (B < 1) stop("B must be >= 1")

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

  make_formula <- function(response, cols, intercept) {
    if (!length(cols)) {
      if (is.null(response)) {
        if (intercept) stats::as.formula("~ 1") else stats::as.formula("~ -1")
      } else {
        if (intercept) stats::as.formula(paste(response, "~ 1")) else stats::as.formula(paste(response, "~ -1"))
      }
    } else {
      rhs <- paste(cols, collapse = " + ")
      if (!intercept) rhs <- paste(rhs, "-1")
      if (is.null(response)) {
        stats::as.formula(paste("~", rhs))
      } else {
        stats::as.formula(paste(response, "~", rhs))
      }
    }
  }

  resp_vars <- all.vars(formula[[2L]])
  if (!length(resp_vars)) stop("Formula must have a response.")
  resp_name <- resp_vars[1L]
  response_str <- paste(deparse(formula[[2L]]), collapse = "")
  response_obj <- tryCatch({
    mf <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
    stats::model.response(mf)
  }, error = function(e) {
    env <- list2env(data, parent = environment(formula))
    eval(formula[[2L]], envir = env)
  })
  if (is.null(response_obj)) {
    if (length(resp_vars) > 1L && all(resp_vars %in% names(data))) {
      response_obj <- as.matrix(data[, resp_vars, drop = FALSE])
    } else if (resp_name %in% names(data)) {
      response_obj <- data[[resp_name]]
    }
  }

  scaler <- NULL
  if (isTRUE(pre_standardize)) {
    num_vars <- names(which(vapply(data, is.numeric, logical(1))))
    num_vars <- setdiff(num_vars, resp_vars)
    if (length(num_vars)) {
      center <- vapply(data[num_vars], function(x) mean(x, na.rm = TRUE), numeric(1))
      scale  <- vapply(data[num_vars], function(x) stats::sd(x, na.rm = TRUE), numeric(1))
      scale[scale == 0] <- 1
      data[num_vars] <- Map(function(x, m, s) (x - m)/s, data[num_vars], center, scale)
      scaler <- list(center = center, scale = scale, vars = num_vars, response = resp_vars)
    }
  }

  base_mu    <- formula
  base_sigma <- if (inherits(base_sigma, "formula")) base_sigma else ~ 1
  base_nu    <- if (inherits(base_nu, "formula")) base_nu else ~ 1
  base_tau   <- if (inherits(base_tau, "formula")) base_tau else ~ 1

  base_mu_terms    <- term_labels(base_mu)
  base_sigma_terms <- term_labels(base_sigma)
  base_nu_terms    <- term_labels(base_nu)
  base_tau_terms   <- term_labels(base_tau)

  scope_filter <- function(scope_formula, drop_terms) {
    if (is.null(scope_formula)) return(NULL)
    tl <- setdiff(term_labels(scope_formula), drop_terms)
    if (!length(tl)) return(NULL)
    f <- stats::as.formula(paste("~", paste(tl, collapse = " + ")))
    environment(f) <- environment(scope_formula)
    f
  }

  mu_scope_use    <- scope_filter(mu_scope, base_mu_terms)
  sigma_scope_use <- scope_filter(sigma_scope, base_sigma_terms)
  nu_scope_use    <- scope_filter(nu_scope, base_nu_terms)
  tau_scope_use   <- scope_filter(tau_scope, base_tau_terms)

  # build base design matrices with sanitised names
  build_base_design <- function(base_formula, prefix, response = NULL) {
    trm <- stats::terms(base_formula)
    intercept <- isTRUE(attr(trm, "intercept") == 1L)
    mm <- stats::model.matrix(base_formula, data)
    assign <- attr(mm, "assign")
    if (ncol(mm) && any(assign == 0)) {
      keep <- assign != 0
      mm <- mm[, keep, drop = FALSE]
      assign <- assign[keep]
    }
    cols <- if (ncol(mm)) {
      make.unique(paste0(prefix, "__base__", make.names(colnames(mm), unique = FALSE)), sep = "__")
    } else character(0)
    if (length(cols)) colnames(mm) <- cols
    term_map <- if (length(assign)) {
      tl <- term_labels(base_formula)
      out <- rep("", length(assign))
      nz <- assign > 0 & assign <= length(tl)
      out[nz] <- tl[assign[nz]]
      out
    } else character(0)
    list(
      matrix = mm,
      cols = cols,
      term_map = term_map,
      terms = term_labels(base_formula),
      intercept = intercept,
      response = response,
      base_formula = base_formula
    )
  }

  mu_base_design    <- build_base_design(base_mu,    "mu",    response = resp_name)
  sigma_base_design <- build_base_design(base_sigma, "sigma")
  nu_base_design    <- build_base_design(base_nu,    "nu")
  tau_base_design   <- build_base_design(base_tau,   "tau")

  base_data <- data
  if (length(mu_base_design$cols))    base_data[mu_base_design$cols]    <- mu_base_design$matrix
  if (length(sigma_base_design$cols)) base_data[sigma_base_design$cols] <- sigma_base_design$matrix
  if (length(nu_base_design$cols))    base_data[nu_base_design$cols]    <- nu_base_design$matrix
  if (length(tau_base_design$cols))   base_data[tau_base_design$cols]   <- tau_base_design$matrix

  base_formulas_internal <- list(
    mu    = make_formula(response_str, mu_base_design$cols,    mu_base_design$intercept),
    sigma = make_formula(NULL,      sigma_base_design$cols, sigma_base_design$intercept),
    nu    = make_formula(NULL,      nu_base_design$cols,    nu_base_design$intercept),
    tau   = make_formula(NULL,      tau_base_design$cols,   tau_base_design$intercept)
  )
  environment(base_formulas_internal$mu)    <- environment(base_mu)
  environment(base_formulas_internal$sigma) <- environment(base_sigma)
  environment(base_formulas_internal$nu)    <- environment(base_nu)
  environment(base_formulas_internal$tau)   <- environment(base_tau)

  make_upper_formula <- function(base_cols, cand_cols, intercept, response = NULL) {
    cols <- unique(c(base_cols, cand_cols))
    make_formula(response, cols, intercept)
  }

  mu_scope_info    <- .sb_prepare_selectboost(data, mu_scope_use,    "mu_scope",    B, corr_func = corr_func, group_fun = group_fun, corr_threshold = c0, use_groups = use_groups)
  sigma_scope_info <- .sb_prepare_selectboost(data, sigma_scope_use, "sigma_scope", B, corr_func = corr_func, group_fun = group_fun, corr_threshold = c0, use_groups = use_groups)
  nu_scope_info    <- .sb_prepare_selectboost(data, nu_scope_use,    "nu_scope",    B, corr_func = corr_func, group_fun = group_fun, corr_threshold = c0, use_groups = use_groups)
  tau_scope_info   <- .sb_prepare_selectboost(data, tau_scope_use,   "tau_scope",   B, corr_func = corr_func, group_fun = group_fun, corr_threshold = c0, use_groups = use_groups)

  upper_formulas_internal <- list(
    mu    = make_upper_formula(mu_base_design$cols,    mu_scope_info$scope$sanitized_colnames,    mu_base_design$intercept,    response_str),
    sigma = make_upper_formula(sigma_base_design$cols, sigma_scope_info$scope$sanitized_colnames, sigma_base_design$intercept),
    nu    = make_upper_formula(nu_base_design$cols,    nu_scope_info$scope$sanitized_colnames,    nu_base_design$intercept),
    tau   = make_upper_formula(tau_base_design$cols,   tau_scope_info$scope$sanitized_colnames,   tau_base_design$intercept)
  )
  environment(upper_formulas_internal$mu)    <- environment(base_mu)
  environment(upper_formulas_internal$sigma) <- environment(base_sigma)
  environment(upper_formulas_internal$nu)    <- environment(base_nu)
  environment(upper_formulas_internal$tau)   <- environment(base_tau)

  dots_original <- list(...)
  weights_full <- NULL
  if ("weights" %in% names(dots_original)) {
    weights_full <- dots_original$weights
    dots_original$weights <- NULL
    if (!is.null(weights_full)) {
      if (is.matrix(weights_full) || (!is.null(dim(weights_full)) && length(dim(weights_full)) >= 2L)) {
        if (nrow(weights_full) != n) {
          stop("`weights` must have the same number of rows as `data`.", call. = FALSE)
        }
      } else if (is.data.frame(weights_full)) {
        if (nrow(weights_full) != n) {
          stop("`weights` data frame must have the same number of rows as `data`.", call. = FALSE)
        }
      } else if (length(weights_full) != n) {
        stop("`weights` must have length equal to `nrow(data)`.", call. = FALSE)
      }
    }
  }
  if ("trace" %in% names(dots_original)) {
    dots_original$trace <- NULL
  }
  drop_final <- c("formula","sigma.formula","nu.formula","tau.formula","data","family")
  dots_original <- dots_original[setdiff(names(dots_original), drop_final)]
  
  glmnet_family <- match.arg(glmnet_family)
  
  selector_globals <- list(
    base_data = base_data,
    base_formulas = base_formulas_internal,
    upper_formulas = upper_formulas_internal,
    family = family,
    direction = direction,
    k = k,
    glmnet_alpha = glmnet_alpha,
    glmnet_family = glmnet_family,
    grpreg_penalty = grpreg_penalty,
    sgl_alpha = sgl_alpha,
    sample_fraction = sample_fraction,
    resp_name = resp_name,
    response = response_obj,
    dots = dots_original,
    weights = weights_full,
    n = n,
    mu_base_cols = mu_base_design$cols,
    sigma_base_cols = sigma_base_design$cols,
    nu_base_cols = nu_base_design$cols,
    tau_base_cols = tau_base_design$cols
  )

  clean_dots <- function(dots) {
    drop <- c("formula","sigma.formula","nu.formula","tau.formula","data","family","trace")
    dots[setdiff(names(dots), drop)]
  }
  selector_globals$dots <- clean_dots(selector_globals$dots)

  subset_dots <- function(idx) {
    dots <- selector_globals$dots
    if (!length(dots)) return(dots)
    if (!length(idx)) return(dots)
    res <- dots
    for (nm in names(res)) {
      val <- res[[nm]]
      if (is.null(val)) next
      if ((is.vector(val) || is.factor(val)) && length(val) == selector_globals$n) {
        res[[nm]] <- val[idx]
      } else if (is.matrix(val) && nrow(val) == selector_globals$n) {
        res[[nm]] <- val[idx, , drop = FALSE]
      } else if (is.data.frame(val) && nrow(val) == selector_globals$n) {
        res[[nm]] <- val[idx, , drop = FALSE]
      } else if (inherits(val, "Matrix") && nrow(val) == selector_globals$n) {
        res[[nm]] <- val[idx, , drop = FALSE]
      } else if (!is.null(dim(val)) && dim(val)[1] == selector_globals$n) {
        res[[nm]] <- val[idx, , drop = FALSE]
      }
    }
    res
  }
  
  subset_special <- function(obj, idx) {
    if (is.null(obj) || !length(idx)) return(obj)
    if (is.data.frame(obj)) {
      if (nrow(obj) == selector_globals$n) {
        return(obj[idx, , drop = FALSE])
      }
      return(obj)
    }
    if (is.matrix(obj) || (!is.null(dim(obj)) && length(dim(obj)) >= 2L)) {
      if (nrow(obj) == selector_globals$n) {
        return(obj[idx, , drop = FALSE])
      }
      return(obj)
    }
    if (length(obj) == selector_globals$n) {
      return(obj[idx])
    }
    tryCatch(obj[idx], error = function(e) obj)
  }
  
  subset_response <- function(resp, idx) {
    if (is.null(resp) || !length(idx)) return(resp)
    if (is.data.frame(resp)) {
      if (nrow(resp) == selector_globals$n) {
        return(resp[idx, , drop = FALSE])
      }
      return(resp)
    }
    if (is.matrix(resp) || (!is.null(dim(resp)) && length(dim(resp)) >= 2L)) {
      if (NROW(resp) == selector_globals$n) {
        return(resp[idx, , drop = FALSE])
      }
      return(resp)
    }
    if (length(resp) == selector_globals$n) {
      return(resp[idx])
    }
    tryCatch(resp[idx], error = function(e) resp)
  }
  
  ensure_matrix <- function(X, nrow_expected = NULL) {
    if (is.null(X)) {
      nr <- if (is.null(nrow_expected)) 0 else nrow_expected
      return(matrix(0, nrow = nr, ncol = 0))
    }
    if (is.matrix(X)) return(X)
    if (is.vector(X)) {
      if (is.null(nrow_expected)) return(matrix(X, ncol = 1))
      return(matrix(X, nrow = nrow_expected))
    }
    as.matrix(X)
  }

  make_selector <- function(param, engine_name, scope_info, base_design) {
    candidate_cols <- scope_info$scope$sanitized_colnames
    term_map <- scope_info$scope$term_map
    base_cols <- switch(param,
      mu = selector_globals$mu_base_cols,
      sigma = selector_globals$sigma_base_cols,
      nu = selector_globals$nu_base_cols,
      tau = selector_globals$tau_base_cols
    )
    function(X, Y) {
      if (!length(candidate_cols)) return(setNames(numeric(0), candidate_cols))
      X <- ensure_matrix(X, selector_globals$n)
      if (!ncol(X)) return(setNames(numeric(0), candidate_cols))
      colnames(X) <- candidate_cols
      size <- max(2L, floor(selector_globals$sample_fraction * selector_globals$n))
      idx <- if (size >= selector_globals$n) seq_len(selector_globals$n) else sample.int(selector_globals$n, size = size, replace = FALSE)
      df <- selector_globals$base_data[idx, , drop = FALSE]
      X_sub <- X[idx, , drop = FALSE]
      for (nm in candidate_cols) df[[nm]] <- X_sub[, nm]
      dots_sub <- subset_dots(idx)
      response_vec <- subset_response(selector_globals$response, idx)
      weights_sub <- subset_special(selector_globals$weights, idx)
      if (is.null(response_vec) && selector_globals$resp_name %in% names(df)) {
        response_vec <- df[[selector_globals$resp_name]]
      }
      response_as_vector <- function(y, engine_label, family = NULL) {
        if (is.null(y)) {
          stop(sprintf("engine='%s' requires a univariate response.", engine_label), call. = FALSE)
        }
        if (is.data.frame(y)) {
          if (nrow(y) != length(idx)) {
            stop(sprintf("engine='%s' received a response with incompatible dimensions.", engine_label), call. = FALSE)
          }
          y <- as.matrix(y)
        }
        if (is.matrix(y)) {
          if (nrow(y) != length(idx)) {
            stop(sprintf("engine='%s' received a response with incompatible dimensions.", engine_label), call. = FALSE)
          }
          if (ncol(y) != 1L) {
            stop(sprintf("engine='%s' currently requires a univariate response.", engine_label), call. = FALSE)
          }
          y <- y[, 1]
        } else if (!is.null(dim(y)) && length(dim(y)) > 1L) {
          stop(sprintf("engine='%s' currently requires a univariate response.", engine_label), call. = FALSE)
        }
        if (length(y) != length(idx)) {
          stop(sprintf("engine='%s' received a response of length %d but expected %d.", engine_label, length(y), length(idx)), call. = FALSE)
        }
        if (!is.null(family) && identical(engine_label, "glmnet")) {
          return(.glmnet_prepare_response(y, family))
        }
        as.numeric(y)
      }
      
      run_step <- function() {
        args <- c(list(
          formula = selector_globals$base_formulas$mu,
          sigma.formula = selector_globals$base_formulas$sigma,
          nu.formula = selector_globals$base_formulas$nu,
          tau.formula = selector_globals$base_formulas$tau,
          family = selector_globals$family,
          data = df,
          trace = FALSE
        ), dots_sub)
        if (!is.null(weights_sub)) {
          args$weights <- weights_sub
        }
        fit0 <- try(do.call(gamlss::gamlss, args), silent = TRUE)
        if (inherits(fit0, "try-error")) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        scope <- list(lower = selector_globals$base_formulas[[param]], upper = selector_globals$upper_formulas[[param]])
        fit1 <- try(gamlss::stepGAIC(fit0, scope = scope, what = param, direction = selector_globals$direction, k = selector_globals$k, trace = FALSE), silent = TRUE)
        if (inherits(fit1, "try-error")) fit1 <- fit0
        cf <- try(stats::coef(fit1, what = param), silent = TRUE)
        out <- setNames(rep(0, length(candidate_cols)), candidate_cols)
        if (!inherits(cf, "try-error") && length(cf)) {
          keep <- intersect(names(cf), candidate_cols)
          if (length(keep)) out[keep] <- cf[keep]
        }
        out
      }

      run_glmnet <- function() {
        if (!requireNamespace("glmnet", quietly = TRUE)) {
          stop("Package 'glmnet' is required for engine='glmnet'.", call. = FALSE)
        }
        X_glm <- as.matrix(df[, c(base_cols, candidate_cols), drop = FALSE])
        penalty <- c(rep(0, length(base_cols)), rep(1, length(candidate_cols)))
        y <- response_as_vector(response_vec, "glmnet", selector_globals$glmnet_family)
        cv <- try(glmnet::cv.glmnet(
          x = X_glm,
          y = y,
          weights = weights_sub,
          alpha = selector_globals$glmnet_alpha,
          family = selector_globals$glmnet_family,
          penalty.factor = penalty,
          standardize = FALSE,
          intercept = TRUE
        ), silent = TRUE)
        if (inherits(cv, "try-error")) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        lambda_min <- cv$lambda.min
        lambda_min <- cv$lambda.min
        if (!is.numeric(lambda_min) || !length(lambda_min) || !is.finite(lambda_min)) {
          return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        }
        fit <- cv$glmnet.fit
        if (is.null(fit) || is.null(fit$beta) || !length(fit$lambda)) {
          return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        }
        idx <- which.min(abs(fit$lambda - lambda_min))
        if (!length(idx) || is.na(idx)) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        beta <- as.matrix(fit$beta[, idx, drop = FALSE])
        if (is.null(rownames(beta))) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        out <- setNames(rep(0, length(candidate_cols)), candidate_cols)
        keep <- intersect(candidate_cols, rownames(beta))
        if (length(keep)) out[keep] <- beta[keep, 1]
        out
      }
      
      run_grpreg <- function() {
        if (!requireNamespace("grpreg", quietly = TRUE)) {
          stop("Package 'grpreg' is required for engine='grpreg'.", call. = FALSE)
        }
        X_mat <- as.matrix(df[, c(base_cols, candidate_cols), drop = FALSE])
        y <- response_as_vector(response_vec, "grpreg")
        term_levels <- unique(term_map[candidate_cols])
        term_levels <- term_levels[term_levels != ""]
        group_idx <- match(term_map[candidate_cols], term_levels)
        if (anyNA(group_idx)) {
          missing <- which(is.na(group_idx))
          offset <- length(term_levels)
          group_idx[missing] <- offset + seq_along(missing)
        }
        group_vec <- c(rep(0, length(base_cols)), group_idx)
        if (!is.null(weights_sub)) {
          cv <- try(grpreg::cv.grpreg(x = X_mat, y = y, group = group_vec, family = "gaussian", penalty = selector_globals$grpreg_penalty, seed = NULL, weights = weights_sub), silent = TRUE)
        } else {
          cv <- try(grpreg::cv.grpreg(x = X_mat, y = y, group = group_vec, family = "gaussian", penalty = selector_globals$grpreg_penalty, seed = NULL), silent = TRUE)
        }
        cv <- try(grpreg::cv.grpreg(x = X_mat, y = y, group = group_vec, family = "gaussian", penalty = selector_globals$grpreg_penalty, seed = NULL), silent = TRUE)
        if (inherits(cv, "try-error")) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        beta <- stats::coef(cv$fit, lambda = cv$lambda.min)
        out <- setNames(rep(0, length(candidate_cols)), candidate_cols)
        keep <- intersect(candidate_cols, rownames(beta))
        if (length(keep)) out[keep] <- beta[keep, 1]
        out
      }

      run_sgl <- function() {
        if (!requireNamespace("SGL", quietly = TRUE)) {
          stop("Package 'SGL' is required for engine='sgl'.", call. = FALSE)
        }
        X_mat <- as.matrix(df[, candidate_cols, drop = FALSE])
        if (!ncol(X_mat)) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        y <- response_as_vector(response_vec, "sgl")
        if (length(base_cols)) {
          base_mat <- as.matrix(df[, base_cols, drop = FALSE])
          fit_base <- try(stats::lm.fit(base_mat, y), silent = TRUE)
          y_work <- if (inherits(fit_base, "try-error")) y else fit_base$residuals
        } else {
          y_work <- y
        }
        term_levels <- unique(term_map[candidate_cols])
        term_levels <- term_levels[term_levels != ""]
        index <- match(term_map[candidate_cols], term_levels)
        if (anyNA(index)) {
          missing <- which(is.na(index))
          offset <- length(term_levels)
          index[missing] <- offset + seq_along(missing)
        }
        data_sgl <- list(x = X_mat, y = as.numeric(y_work))
        cv <- try(SGL::cvSGL(data = data_sgl, index = index, type = "linear", alpha = selector_globals$sgl_alpha, standardize = FALSE), silent = TRUE)
        if (inherits(cv, "try-error")) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        lam <- cv$lam[which.min(cv$lldiff)]
        fit <- try(SGL::SGL(data = data_sgl, index = index, type = "linear", alpha = selector_globals$sgl_alpha, standardize = FALSE, lam = lam), silent = TRUE)
        if (inherits(fit, "try-error")) return(setNames(rep(0, length(candidate_cols)), candidate_cols))
        beta <- as.numeric(fit$beta)
        names(beta) <- candidate_cols
        beta
      }

      switch(engine_name,
        stepGAIC = run_step(),
        glmnet   = run_glmnet(),
        grpreg   = run_grpreg(),
        sgl      = run_sgl()
      )
    }
  }

  selectors <- list(
    mu    = make_selector("mu",    engine,       mu_scope_info,    mu_base_design),
    sigma = make_selector("sigma", engine_sigma, sigma_scope_info, sigma_base_design),
    nu    = make_selector("nu",    engine_nu,    nu_scope_info,    nu_base_design),
    tau   = make_selector("tau",   engine_tau,   tau_scope_info,   tau_base_design)
  )

  run_param <- function(param, scope_info) {
    cols <- scope_info$scope$sanitized_colnames
    if (!length(cols)) {
      return(list(coefs = matrix(numeric(0), nrow = 0, ncol = B), term_counts = numeric(0), term_props = numeric(0)))
    }
    X <- ensure_matrix(scope_info$matrix, selector_globals$n)
    sims <- scope_info$simulations
    if (!is.null(dim(sims))) {
      dimnames(sims)[[1]] <- seq_len(nrow(sims))
    }
    response_full <- selector_globals$response
    if (is.null(response_full)) {
      response_full <- data[[resp_name]]
    }
    res <- SelectBoost::boost.apply(X, sims, response_full, selectors[[param]], verbose = FALSE)
    res <- ensure_matrix(res, length(cols))
    if (!is.null(cols)) {
      rownames(res) <- cols
    }
    freq <- SelectBoost::boost.select(res)
    freq <- setNames(as.numeric(freq), rownames(res))
    cnt <- freq * B
    term_map <- scope_info$scope$term_map
    names(cnt) <- names(freq)
    valid_terms <- term_map[names(cnt)]
    split_vals <- if (length(cnt)) split(cnt, valid_terms) else list()
    split_vals <- split_vals[names(split_vals) != ""]
    term_counts <- if (length(split_vals)) vapply(split_vals, max, numeric(1)) else numeric(0)
    all_terms <- scope_info$scope$term_labels
    term_counts_full <- if (length(all_terms)) setNames(rep(0, length(all_terms)), all_terms) else numeric(0)
    if (length(term_counts)) term_counts_full[names(term_counts)] <- term_counts
    term_props <- if (length(term_counts_full)) term_counts_full / B else numeric(0)
    list(coefs = res, term_counts = term_counts_full, term_props = term_props)
  }

  results <- list(
    mu    = run_param("mu",    mu_scope_info),
    sigma = run_param("sigma", sigma_scope_info),
    nu    = run_param("nu",    nu_scope_info),
    tau   = run_param("tau",   tau_scope_info)
  )

  build_selection_df <- function(param, res, base_terms) {
    term_counts <- res$term_counts
    term_props  <- res$term_props
    if (length(term_counts)) {
      df <- data.frame(parameter = param, term = names(term_counts), count = as.integer(round(term_counts)), prop = as.numeric(term_props), stringsAsFactors = FALSE)
    } else {
      df <- data.frame(parameter = character(0), term = character(0), count = integer(0), prop = numeric(0))
    }
    if (length(base_terms)) {
      base_df <- data.frame(parameter = param, term = base_terms, count = rep(B, length(base_terms)), prop = rep(1, length(base_terms)), stringsAsFactors = FALSE)
      df <- rbind(df, base_df)
    }
    df
  }

  selection_tab <- do.call(rbind, list(
    build_selection_df("mu",    results$mu,    mu_base_design$terms),
    build_selection_df("sigma", results$sigma, sigma_base_design$terms),
    build_selection_df("nu",    results$nu,    nu_base_design$terms),
    build_selection_df("tau",   results$tau,   tau_base_design$terms)
  ))

  select_terms <- function(base_formula, res, scope_formula_use) {
    if (is.null(scope_formula_use)) return(base_formula)
    term_props <- res$term_props
    if (!length(term_props)) return(base_formula)
    keep <- names(term_props)[term_props >= pi_thr]
    if (!length(keep)) return(base_formula)
    add_terms_formula(base_formula, keep)
  }

  final_mu    <- select_terms(base_mu,    results$mu,    mu_scope_use)
  final_sigma <- select_terms(base_sigma, results$sigma, sigma_scope_use)
  final_nu    <- select_terms(base_nu,    results$nu,    nu_scope_use)
  final_tau   <- select_terms(base_tau,   results$tau,   tau_scope_use)

  final_args <- c(list(
    formula = final_mu,
    sigma.formula = final_sigma,
    nu.formula = final_nu,
    tau.formula = final_tau,
    data = data,
    family = family,
    trace = trace
  ), dots_original)
  if (!is.null(weights_full)) {
    final_args$weights <- weights_full
  }
  final_fit <- do.call(gamlss::gamlss, final_args)

  out <- list(
    final_fit = final_fit,
    final_formula = list(mu = final_mu, sigma = final_sigma, nu = final_nu, tau = final_tau),
    selection = selection_tab,
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
#' 
#' @return Invisibly returns x the plotted `sb_gamlss` object.
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
