
# Utilities for deviance / log-likelihood on validation data

#' Predict distribution parameters on newdata
#' @param fit a gamlss fit
#' @param newdata data.frame
#' @return list with available components: mu, sigma, nu, tau
predict_params <- function(fit, newdata) {
  out <- list()
  for (w in c("mu","sigma","nu","tau")) {
    val <- try(stats::predict(fit, newdata = newdata, what = w, type = "response"), silent = TRUE)
    if (!inherits(val, "try-error")) out[[w]] <- as.numeric(val)
  }
  out
}

#' Get a density function for a gamlss family
#' @param fit a gamlss fit (or family name)
#' @return function(x, mu, sigma, nu, tau, log=FALSE)
get_density_fun <- function(fit) {
  fam <- try(fit$family[1], silent = TRUE)
  if (inherits(fam, "try-error") || is.null(fam)) fam <- as.character(fit$family[1])
  fname <- paste0("d", fam)
  fn <- try(getFromNamespace(fname, "gamlss.dist"), silent = TRUE)
  if (inherits(fn, "try-error")) stop("Density function ", fname, " not found in gamlss.dist")
  fn
}

#' Log-likelihood (sum) on newdata given a gamlss fit
#' @param fit gamlss object
#' @param newdata data.frame
#' @return numeric scalar: sum of log-likelihoods
loglik_gamlss_newdata <- function(fit, newdata) {
  dens <- get_density_fun(fit)
  pars <- predict_params(fit, newdata)
  y <- model.response(model.frame(stats::formula(fit, what = "mu"), data = newdata))
  # match formal args of density to pass only what is needed
  fmls <- names(formals(dens))
  args <- list(x = y, log = TRUE)
  for (nm in c("mu","sigma","nu","tau")) {
    if (nm %in% names(pars) && nm %in% fmls) args[[nm]] <- pars[[nm]]
  }
  ll <- do.call(dens, args)
  sum(ll[is.finite(ll)])
}

#' K-fold deviance for an sb_gamlss configuration
#' @param K folds
#' @param build_fit function(...) that returns an sb_gamlss object
#' @param data data.frame used inside build_fit
#' @return numeric: mean deviance across folds (-2 * mean loglik)
cv_deviance_sb <- function(K, build_fit, data) {
  n <- NROW(data)
  idx <- sample(rep_len(1:K, n))
  dev <- numeric(K)
  for (k in 1:K) {
    train <- data[idx != k, , drop = FALSE]
    valid <- data[idx == k, , drop = FALSE]
    fit <- try(build_fit(train), silent = TRUE)
    if (inherits(fit, "try-error")) { dev[k] <- Inf; next }
    # refit is inside sb_gamlss; compute loglik on validation
    # need the final gamlss fit stored at fit$final_fit
    g <- fit$final_fit
    # predict parameters on validation via g with updated data
    ll <- try(loglik_gamlss_newdata_fast(g, valid), silent = TRUE)
    dev[k] <- if (inherits(ll, "try-error")) Inf else (-2 * ll)
  }
  mean(dev)
}


# Fast-path deviance for common families (avoids gamlss.dist lookup)
# Currently supports:
# - "NO": Normal(mu, sigma) via stats::dnorm
# - "PO": Poisson(mu) via stats::dpois

loglik_gamlss_newdata_fast <- function(fit, newdata) {
  fam <- try(fit$family[1], silent = TRUE)
  fam <- if (inherits(fam, "try-error")) NULL else as.character(fam)
  fam <- gsub("[() ]", "", fam, perl = TRUE)
  pars <- predict_params(fit, newdata)
  mf <- model.frame(stats::formula(fit, what = "mu"), data = newdata)
  y <- model.response(mf)

  # helpers
  safe_sum <- function(x) sum(x[is.finite(x)], na.rm = TRUE)

  # ---- Normal (NO): mu, sigma
  if (identical(fam, "NO") && all(c("mu","sigma") %in% names(pars))) {
    ll <- stats::dnorm(y, mean = pars$mu, sd = pars$sigma, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Poisson (PO): mu
  if (identical(fam, "PO") && "mu" %in% names(pars)) {
    ll <- stats::dpois(y, lambda = pars$mu, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Log-Normal (LOGNO): mu = meanlog, sigma = sdlog
  if (identical(fam, "LOGNO") && all(c("mu","sigma") %in% names(pars))) {
    ll <- stats::dlnorm(y, meanlog = pars$mu, sdlog = pars$sigma, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Gamma (GA): Var = (sigma^2)*(mu^2) => shape = 1/sigma^2, scale = mu*sigma^2
  if (identical(fam, "GA") && all(c("mu","sigma") %in% names(pars))) {
    shape <- 1/(pars$sigma^2)
    scale <- pars$mu * (pars$sigma^2)
    ll <- stats::dgamma(y, shape = shape, scale = scale, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Inverse Gaussian (IG): Var = sigma^2 * mu^3; use closed-form density
  if (identical(fam, "IG") && all(c("mu","sigma") %in% names(pars))) {
    mu <- pars$mu; sg <- pars$sigma
    yv <- as.numeric(y)
    # domain fixes
    yv[yv <= 0] <- NA_real_
    mu[mu <= 0] <- NA_real_
    sg[sg <= 0] <- NA_real_
    ll <- -0.5*log(2*pi) - log(sg) - 1.5*log(yv) - ((yv - mu)^2) / (2 * mu^2 * sg^2 * yv)
    return(safe_sum(ll))
  }

  # ---- Negative Binomial type I (NBI): Var = mu + sigma * mu^2
  if (identical(fam, "NBI") && all(c("mu","sigma") %in% names(pars))) {
    size <- 1/pars$sigma
    ll <- stats::dnbinom(y, size = size, mu = pars$mu, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Negative Binomial type II (NBII): Var = mu * (1 + sigma) => size = mu/sigma, prob = 1/(1+sigma)
  if (identical(fam, "NBII") && all(c("mu","sigma") %in% names(pars))) {
    size <- pars$mu / pmax(pars$sigma, .Machine$double.eps)
    ll <- stats::dnbinom(y, size = size, mu = pars$mu, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Binomial (BI): y is two-column (success, failure); mu is success prob; bd = trials = rowSums(y)
  if (identical(fam, "BI") && "mu" %in% names(pars) && is.matrix(y) && ncol(y) == 2) {
    successes <- y[,1]
    trials <- rowSums(y)
    ll <- stats::dbinom(successes, size = trials, prob = pars$mu, log = TRUE)
    return(safe_sum(ll))
  }


  # ---- Logit-Normal (LOGITNO): logit(y) ~ Normal(mu, sigma); y in (0,1)
  if (identical(fam, "LOGITNO") && all(c("mu","sigma") %in% names(pars))) {
    ye <- pmin(pmax(as.numeric(y), .Machine$double.eps), 1 - .Machine$double.eps)
    z  <- log(ye/(1 - ye))  # logit
    ll <- stats::dnorm(z, mean = pars$mu, sd = pars$sigma, log = TRUE) - (log(ye) + log1p(-ye))
    return(safe_sum(ll))
  }

  # ---- Geometric (GEOM): mean mu => prob p = 1/(1+mu); support {0,1,2,...}
  if (identical(fam, "GEOM") && "mu" %in% names(pars)) {
    p <- 1/(1 + pars$mu)
    ll <- stats::dgeom(y, prob = p, log = TRUE)
    return(safe_sum(ll))
  }
  # ---- Logistic (LO): mu = location, sigma = scale
  if (identical(fam, "LO") && all(c("mu","sigma") %in% names(pars))) {
    ll <- stats::dlogis(y, location = pars$mu, scale = pars$sigma, log = TRUE)
    return(safe_sum(ll))
  }

  # ---- Beta (BE): Var = sigma^2 * mu * (1 - mu)
  # Mapping to standard Beta(alpha, beta):
  #   phi = 1/sigma^2 - 1
  #   alpha = mu * phi
  #   beta  = (1 - mu) * phi
  if (identical(fam, "BE") && all(c("mu","sigma") %in% names(pars))) {
    mu <- pmin(pmax(pars$mu, .Machine$double.eps), 1 - .Machine$double.eps)
    sg <- pmin(pmax(pars$sigma, 1e-8), 1 - 1e-8)
    phi <- 1/(sg*sg) - 1
    alpha <- mu * phi
    beta  <- (1 - mu) * phi
    ll <- stats::dbeta(y, shape1 = alpha, shape2 = beta, log = TRUE)
    return(safe_sum(ll))
  }


  # ---- Native gamlss.dist fast routes for additional families via direct density lookup
  dens_name <- paste0("d", fam)
  dn <- try(getFromNamespace(dens_name, "gamlss.dist"), silent = TRUE)
  if (!inherits(dn, "try-error")) {
    fmls <- names(formals(dn))
    args <- list(log = TRUE)
    if (is.matrix(y) && ncol(y) == 2) {
      successes <- y[, 1]
      trials <- rowSums(y)
      args$x <- successes
      if ("bd" %in% fmls) args$bd <- trials
    } else {
      args$x <- as.numeric(y)
    }
    for (nm in c("mu","sigma","nu","tau")) {
      if (nm %in% names(pars) && nm %in% fmls) args[[nm]] <- pars[[nm]]
    }
    ll <- try(do.call(dn, args), silent = TRUE)
    if (!inherits(ll, "try-error")) {
      return(safe_sum(ll))
    }
  }

  # fallback to generic density in gamlss.dist
  loglik_gamlss_newdata(fit, newdata)
}

