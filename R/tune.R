
#' Tune select engines/penalties via a small stability run
#'
#' Evaluates a grid of configurations and picks the one maximizing a stability-based score,
#' optionally penalized by complexity. Designed to be lightweight and robust.
#'
#' @param config_grid a list of named lists, each containing a subset of sb_gamlss args
#'                    (e.g., list(engine="grpreg", engine_sigma="sgl", grpreg_penalty="grLasso", sgl_alpha=0.9))
#' @param base_args a named list of arguments passed to `sb_gamlss()` common to all configs
#' @param score_lambda non-negative; penalty weight on total number of stable terms (complexity)
#' @param B_small number of bootstraps to use during tuning (defaults to 30)
#' @param metric Character; `"stability"` or `"deviance"` (K-fold CV).
#' @param K Integer; folds for deviance CV.
#' @param progress Logical; show progress bar across configs.
#' @param score_lambda Numeric; complexity penalty weight for stability metric.
#' @return a list: best_config, scores (data.frame), and the fitted sb_gamlss object for the best config.
#' @export
tune_sb_gamlss <- function(config_grid, base_args, score_lambda = 0, B_small = 30, metric = c('stability','deviance'), K = 3, progress = TRUE) {
  stopifnot(is.list(config_grid), length(config_grid) >= 1L)
  metric <- match.arg(metric)
  pb <- NULL
  if (isTRUE(progress)) pb <- utils::txtProgressBar(min = 0, max = length(config_grid), style = 3)
  i_prog <- 0
  eval_one <- function(cfg) {
    args <- base_args
    args$B <- B_small
    # merge/override with cfg
    for (nm in names(cfg)) args[[nm]] <- cfg[[nm]]
    if (metric == 'stability') {
      fit <- try(do.call(sb_gamlss, args), silent = TRUE)
      if (inherits(fit, 'try-error')) return(list(score = -Inf, fit = NULL, cfg = cfg))
      tab <- selection_table(fit)
      if (is.null(tab) || !NROW(tab)) return(list(score = -Inf, fit = fit, cfg = cfg))
      pi_thr <- fit$pi_thr %||% 0.6
      pos <- pmax(0, tab$prop - pi_thr)
      mass <- sum(pos)
      stable_terms <- sum(tab$prop >= pi_thr)
      score <- mass - score_lambda * stable_terms
      res <- list(score = score, fit = fit, cfg = cfg)
    } else {
      # deviance-based CV metric
      dat_all <- args$data
      build_fit <- function(dtrain) {
        args$data <- dtrain
        do.call(sb_gamlss, args)
      }
      dev <- cv_deviance_sb(K, build_fit, dat_all)
      # lower deviance is better; convert to score = -deviance
      fit <- try(do.call(sb_gamlss, args), silent = TRUE)
      res <- list(score = -dev, fit = if (inherits(fit, 'try-error')) NULL else fit, cfg = cfg)
    }
    if (!is.null(pb)) { i_prog <<- i_prog + 1; utils::setTxtProgressBar(pb, i_prog) }
    res
  }
  res <- lapply(config_grid, eval_one)
  scores <- data.frame(
    idx = seq_along(res),
    score = sapply(res, `[[`, "score"),
    config = I(lapply(res, `[[`, "cfg")),
    stringsAsFactors = FALSE
  )
  best_idx <- which.max(scores$score)
  if (!is.null(pb)) close(pb)
  list(best_config = scores$config[[best_idx]], scores = scores, best_fit = res[[best_idx]]$fit)
}
