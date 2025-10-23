## ----setup, include = FALSE---------------------------------------------------
#file.edit(normalizePath("~/.Renviron"))
LOCAL <- identical(Sys.getenv("LOCAL"), "TRUE")
#LOCAL=TRUE
knitr::opts_chunk$set(purl = LOCAL)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
set.seed(1)
library(gamlss)
library(SelectBoost.gamlss)

n <- 400
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
mu <- 1 + 1.5*x1 - 1.2*x3
y  <- gamlss.dist::rNO(n, mu = mu, sigma = 1)

dat <- data.frame(y, x1, x2, x3)

fit <- sb_gamlss(
  y ~ 1,
  mu_scope = ~ x1 + x2 + x3,
  sigma_scope = ~ x1 + x2,
  family = gamlss.dist::NO(),
  data = dat,
  B = 50,
  sample_fraction = 0.7,
  pi_thr = 0.6,
  k = 2,
  pre_standardize = TRUE,
  trace = FALSE
)

fit$final_formula
plot_sb_gamlss(fit, top = 10)

