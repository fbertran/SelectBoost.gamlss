## -----------------------------------------------------------------------------
library(gamlss)
library(SelectBoost.gamlss)

set.seed(1)
n <- 400
x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
y  <- gamlss.dist::rNO(n, mu = 1 + 1.3*x1 - 1.0*x3, sigma = 1)
dat <- data.frame(y, x1, x2, x3)

g <- sb_gamlss_c0_grid(
  y ~ 1, data = dat, family = gamlss.dist::NO(),
  mu_scope = ~ x1 + x2 + x3, sigma_scope = ~ x1 + x2,
  c0_grid = seq(0.2, 0.8, by = 0.2), B = 40, pi_thr = 0.6, pre_standardize = TRUE, trace = FALSE
)

cf <- confidence_functionals(
  g, pi_thr = 0.6, q = c(0.5, 0.8, 0.9),
  weight_fun = NULL, conservative = FALSE, method = "trapezoid"
)

plot(cf, top = 10, label_top = 6)
plot_stability_curves(g, terms = c("x1","x3"), parameter = "mu")

