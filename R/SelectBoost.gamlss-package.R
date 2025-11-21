#' @keywords internal
#' @aliases SelectBoost.gamlss-package SelectBoost.gamlss NULL
#'
#' @references
#' Bertrand, F. and Maumy, M. (2024). *An Improvement for Variable Selection for
#' Generalized Additive Models for Location, Shape and Scale and Quantile
#' Regression*. Joint Statistics Meetings (JSM) 2024, Portland, OR. Presented the
#' SelectBoost approach for GAMLSS and quantile regression, emphasizing
#' correlation-aware resampling to improve recall and precision when predictors
#' are numerous and highly correlated. 
#'
#' @author This package was written by Frederic Bertrand.
#' Maintainer: Frederic Bertrand <frederic.bertrand@@lecnam.net>
#'
#' @examplesIf requireNamespace("gamlss.dist", quietly = TRUE)
#' set.seed(1)
#' dat <- data.frame(
#'   y = gamlss.dist::rNO(80, mu = 0),
#'   x1 = rnorm(80),
#'   x2 = rnorm(80)
#' )
#' fit <- SelectBoost_gamlss(
#'   y ~ 1,
#'   data = dat,
#'   family = gamlss.dist::NO(),
#'   mu_scope = ~ x1 + x2,
#'   B = 10,
#'   pi_thr = 0.6,
#'   trace = FALSE
#' )
#' fit$final_formula
#'
"_PACKAGE"


#' @importFrom graphics abline
#' @importFrom graphics barplot
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot.new
#' @importFrom graphics text
#' @importFrom graphics title
# #' @importFrom rlang sym
#' @importFrom rlang .data
#' @importFrom stats aggregate
#' @importFrom stats coef
#' @importFrom stats logLik
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @importFrom stats qlogis
#' @importFrom stats setNames
#' @importFrom stats update
#' @importFrom utils getFromNamespace
#' @importFrom utils head
#' @importFrom utils tail
#' @useDynLib SelectBoost.gamlss, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL


# useDynLib(SelectBoost.gamlss, .registration = TRUE)
# export(sb_gamlss)
# export(SelectBoost_gamlss)
# export(selection_table)
# export(plot_sb_gamlss)
# export(coef.sb_gamlss)
# export(predict.sb_gamlss)
# export(AICc_gamlss)
# importFrom(stats, terms, model.frame, update, predict, model.matrix, na.omit)
# importFrom(utils, head)
# importFrom(SelectBoost, group_func_2)
# S3method(summary,SelectBoost_gamlss)
# S3method(plot,summary.SelectBoost_gamlss)
# S3method(plot,SelectBoost_gamlss)
# export(sb_gamlss_c0_grid)
# export(confidence_table)
# export(autoboost_gamlss)
# export(fastboost_gamlss)
# S3method(plot,SelectBoost_gamlss_grid)
# export(confidence_functionals)
# export(plot_stability_curves)
# S3method(plot,sb_confidence)
# importFrom(Rcpp, sourceCpp)
# export(tune_sb_gamlss)
# export(knockoff_filter_mu)
# export(knockoff_filter_param)
# importFrom(splines, bs)
# export(fast_vs_generic_ll)
# export(check_fast_vs_generic)
# export(effect_plot)
