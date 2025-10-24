
# Internal helpers for generating synthetic data across many gamlss families

#' Reasonable defaults
#' 
#' Adjust as needed per family docs
#' 
#' @export
.family_defaults <- function() {
  # Reasonable defaults; adjust as needed per family docs
  list(
    NO     = list(mu = 0,   sigma = 1),
    PO     = list(mu = 2.5),
    LOGNO  = list(mu = 0,   sigma = 0.6),
    GA     = list(mu = 2,   sigma = 0.5),
    IG     = list(mu = 2,   sigma = 0.5),
    LO     = list(mu = 0,   sigma = 1),
    LOGITNO= list(mu = 0,   sigma = 1),
    GEOM   = list(mu = 3),
    BE     = list(mu = 0.4, sigma = 0.2),
    NBI    = list(mu = 2,   sigma = 0.5),
    NBII   = list(mu = 2,   sigma = 0.5),
    LOGLOG = list(mu = 0,   sigma = 1),
    DEL    = list(mu = 2,   sigma = 0.3, nu = 0.2),
    ZAGA   = list(mu = 2,   sigma = 0.5, nu = 0.1),
    ZIP    = list(mu = 2.5, nu = 0.2),
    ZINBI  = list(mu = 2,   sigma = 0.5, nu = 0.2),
    DPO    = list(mu = 2.5, sigma = 1.2),
    GPO    = list(mu = 2.5, sigma = 0.2),
    ZAIG   = list(mu = 2,   sigma = 0.5, nu = 0.1),
    ZALG   = list(mu = 0,   sigma = 0.6, nu = 0.1),
    BCT    = list(mu = 0,   sigma = 1,   nu = 0.5, tau = 2),
    BCPE   = list(mu = 0,   sigma = 1,   nu = 0.5, tau = 2),
    ZIPF   = list(mu = 2.5, sigma = 1),
    ZIPFmu = list(mu = 2.5),
    SICHEL = list(mu = 2.5, sigma = 0.5, nu = 1.2),
    GLG    = list(mu = 0,   sigma = 1,   nu = 1.5),
    BETA4  = list(mu = 0.5, sigma = 1,   nu = 2,   tau = 3),
    RS     = list(mu = 0,   sigma = 1,   nu = 2),
    WEI    = list(mu = 2,   sigma = 1),
    GIG    = list(mu = 2,   sigma = 1,   nu = 1)
  )
}

#' Try to generate y for a family
#' 
#' Returns NULL if generator not available or fails
#' 
#' @export
.gen_family <- function(fam, n) {
  rname <- paste0("r", fam)
  gen <- try(getFromNamespace(rname, "gamlss.dist"), silent = TRUE)
  if (inherits(gen, "try-error")) return(NULL)
  par <- .family_defaults()[[fam]]
  if (is.null(par)) return(NULL)
  args <- c(list(n = n), par)
  y <- try(do.call(gen, args), silent = TRUE)
  if (inherits(y, "try-error")) return(NULL)
  y
}

#' Per-family numeric tolerance for equality checks
#' 
#' @export
.family_tolerance <- function() {
  # Defaults are conservative; heavier-tailed / complex families get looser tol
  tol <- list(
    .default = 1e-6,
    IG = 1e-6,
    LOGNO = 1e-6,
    GA = 1e-6,
    NO = 1e-8,
    PO = 1e-8,
    LO = 1e-8,
    LOGITNO = 1e-6,
    GEOM = 1e-8,
    BE = 1e-6,
    NBI = 1e-6,
    NBII = 1e-6,
    LOGLOG = 1e-6,
    DEL = 1e-6,
    ZAGA = 1e-6,
    ZIP = 1e-6,
    ZINBI = 1e-6,
    DPO = 1e-6,
    GPO = 1e-6,
    ZAIG = 1e-6,
    ZALG = 1e-6,
    BCT = 1e-5,
    BCPE = 1e-5,
    ZIPF = 1e-6,
    ZIPFmu = 1e-6,
    SICHEL = 1e-5,
    GLG = 1e-6,
    BETA4 = 1e-6,
    RS = 1e-6,
    WEI = 1e-6,
    GIG = 1e-6
  )
  tol
}

# #' `%||%`
# #' 
# #' @export
# `%||%` <- function(a, b) if (is.null(a)) b else a
