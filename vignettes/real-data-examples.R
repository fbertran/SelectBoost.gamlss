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
library(SelectBoost.gamlss)
library(gamlss)
library(MASS)           # Insurance, quine
set.seed(1)

## -----------------------------------------------------------------------------
if (requireNamespace("gamlss.data", quietly = TRUE)) {
 boys <- NULL
  load_ok <- tryCatch({
    utils::data("boys", package = "gamlss.data", envir = environment())
    TRUE
  }, error = function(e) {
    cat("Unable to load 'boys' dataset:", conditionMessage(e), "— skipping Dutch boys growth example.\n")
    FALSE
  })
  if (!load_ok || is.null(boys)) {
    cat("Dataset 'boys' unavailable; skipping Dutch boys growth example.\n")
  } else {
    boys <- subset(boys, !is.na(height) & !is.na(age) & !is.na(tanner))
    boys$tanner <- factor(boys$tanner)

    fit_growth <- sb_gamlss(
      height ~ 1, data = boys, family = gamlss.dist::BCT(),
      mu_scope    = ~ pb(age) + pb(age):tanner,
      sigma_scope = ~ pb(age),
      nu_scope    = ~ pb(age),
      engine = "grpreg", engine_sigma = "grpreg", engine_nu = "grpreg",
      grpreg_penalty = "grLasso",
      B = 80, pi_thr = 0.6, pre_standardize = TRUE,
      parallel = "auto", trace = FALSE
    )

    selection_table(fit_growth)[1:12, ]

    # Effect plot for age on mu
    if (requireNamespace('ggplot2', quietly = TRUE)) {
      print(effect_plot(fit_growth, 'age', boys, what = 'mu'))
    }
    # Stability curves for the most important terms (optional)
    g <- sb_gamlss_c0_grid(height ~ 1, data = boys, family = gamlss.dist::BCT(),
    mu_scope = ~ pb(age) + pb(age):tanner, sigma_scope = ~ pb(age),
    nu_scope = ~ pb(age), c0_grid = seq(0.2, 0.8, by = 0.1),
    B = 40, pi_thr = 0.6, pre_standardize = TRUE, parallel = "auto")
    plot_stability_curves(g, terms = c("pb(age)", "pb(age):tanner"), parameter = "mu")
  }
} else {
  cat("Package 'gamlss.data' not installed; skipping Dutch boys growth example.\n")
}

## -----------------------------------------------------------------------------
data(Insurance, package = "MASS")
ins <- transform(Insurance, logH = log(Holders))

fit_po <- sb_gamlss(
  Claims ~ offset(logH), data = ins, family = gamlss.dist::PO(),
  mu_scope = ~ Age + District + Group,
  engine = "glmnet", glmnet_alpha = 1,    # lasso
  B = 100, pi_thr = 0.6, pre_standardize = FALSE,
  parallel = "auto", trace = FALSE
)

selection_table(fit_po)

## -----------------------------------------------------------------------------
data(faithful)
faith <- transform(faithful, eru2 = eruptions^2)

fit_ga <- sb_gamlss(
  waiting ~ 1, data = faith, family = gamlss.dist::GA(),
  mu_scope    = ~ pb(eruptions) + eru2,
  sigma_scope = ~ pb(eruptions),
  engine = "glmnet", glmnet_alpha = 0.5,   # elastic-net
  B = 60, pi_thr = 0.6, pre_standardize = TRUE,
  parallel = "auto", trace = FALSE
)

selection_table(fit_ga)

# Effect plot for eruptions on mu
if (requireNamespace('ggplot2', quietly = TRUE)) {
  print(effect_plot(fit_ga, 'eruptions', faith, what = 'mu'))
}

## -----------------------------------------------------------------------------
mt <- transform(mtcars,
  am = as.integer(am),
  cyl = factor(cyl), gear = factor(gear), carb = factor(carb), vs = factor(vs))

fit_bin <- sb_gamlss(
  am ~ 1, data = mt, family = gamlss.dist::BI(),
  mu_scope = ~ wt + hp + qsec + cyl + gear + carb + vs,
  engine = "grpreg", grpreg_penalty = "grLasso",
  B = 80, pi_thr = 0.6, pre_standardize = TRUE,
  parallel = "auto", trace = FALSE
)

selection_table(fit_bin)

## -----------------------------------------------------------------------------
data(quine, package = "MASS")
fit_nb2 <- sb_gamlss(
  Days ~ 1, data = quine, family = gamlss.dist::NBII(),
  mu_scope = ~ Eth + Sex + Lrn + Age,
  engine = "grpreg", grpreg_penalty = "grLasso",
  B = 80, pi_thr = 0.6, pre_standardize = FALSE,
  parallel = "auto", trace = FALSE
)

selection_table(fit_nb2)

