## -----------------------------------------------------------------------------
library(SelectBoost.gamlss)
library(gamlss)

## -----------------------------------------------------------------------------
bc <- NULL
if (requireNamespace("pscl", quietly = TRUE)) {
  data("bioChemists", package = "pscl")
  bc <- pscl::bioChemists
  needed <- c("art", "kid5", "phd", "fem", "mar", "ment")
  missing_cols <- setdiff(needed, names(bc))
  if (length(missing_cols)) {
    cat(
      "bioChemists dataset is missing columns:",
      paste(missing_cols, collapse = ", "),
      "— skipping ZIP example.\n"
    )
    bc <- NULL
  } else {
    keep <- stats::complete.cases(bc[, needed])
    bc <- bc[keep, , drop = FALSE]
    bc$fem <- factor(bc$fem)
    bc$mar <- factor(bc$mar)
    bc$kid5 <- factor(bc$kid5)

    # ZIP engine via stepwise (mu) + grouped/elastic features
    set.seed(10)
    fit_zip <- sb_gamlss(
      art ~ 1, data = bc, family = gamlss.dist::ZIP(),
      mu_scope = ~ kid5 + fem + mar + phd + ment,
      engine = "grpreg", grpreg_penalty = "grLasso",
      B = 80, pi_thr = 0.6, pre_standardize = FALSE, parallel = "auto", trace = FALSE
    )
    head(selection_table(fit_zip))
  }
} else {
  cat("Package 'pscl' not installed; skipping bioChemists example.\n")
}

## -----------------------------------------------------------------------------
if (requireNamespace("pscl", quietly = TRUE) && !is.null(bc)) {
  set.seed(11)
  fit_zinb <- sb_gamlss(
    art ~ 1, data = bc, family = gamlss.dist::ZINBI(),
    mu_scope = ~ kid5 + fem + mar + phd + ment,
    engine = "grpreg", grpreg_penalty = "grLasso",
    B = 80, pi_thr = 0.6, pre_standardize = FALSE, parallel = "auto", trace = FALSE
  )
  head(selection_table(fit_zinb))
} else if (!requireNamespace("pscl", quietly = TRUE)) {
  cat("Package 'pscl' not installed; skipping bioChemists example.\n")
} else {
  cat("bioChemists dataset unavailable; skipping ZINBI example.\n")
}

## -----------------------------------------------------------------------------
if (requireNamespace("pscl", quietly = TRUE)) {
  if (exists("fit_zip", inherits = FALSE) && exists("fit_zinb", inherits = FALSE)) {
    zip_g <- fit_zip$final_fit
    zinb_g <- fit_zinb$final_fit
    if (inherits(zip_g, "gamlss") && inherits(zinb_g, "gamlss")) {
      safe_scalar <- function(expr) {
        out <- tryCatch(expr, error = function(e) NA_real_)
        if (length(out) == 0L) NA_real_ else as.numeric(out[1L])
      }
      met <- data.frame(
        model = c("ZIP", "ZINBI"),
        AIC = c(safe_scalar(AIC(zip_g)), safe_scalar(AIC(zinb_g))),
        GAIC_k2 = c(safe_scalar(gamlss::GAIC(zip_g, k = 2)), safe_scalar(gamlss::GAIC(zinb_g, k = 2))),
        deviance = c(safe_scalar(zip_g$deviance), safe_scalar(zinb_g$deviance))
      )
      if (any(!is.na(met$AIC)) || any(!is.na(met$GAIC_k2)) || any(!is.na(met$deviance))) {
        print(met)
      } else {
        cat("Model fit statistics unavailable; skipping metric table.\n")
      }
    } else {
      cat("Final fits not found; metrics skipped.\n")
    }
  } else {
    cat("Models not available; run the fitting chunks first.\n")
  }
} else {
  cat("Package 'pscl' not installed; skipping bioChemists example.\n")
}

## -----------------------------------------------------------------------------
aq <- subset(airquality, !is.na(Ozone) & !is.na(Wind) & !is.na(Temp) & !is.na(Solar.R))
aq$Month <- factor(aq$Month)

set.seed(20)
fit_zaga <- sb_gamlss(
  Ozone ~ 1, data = aq, family = gamlss.dist::ZAGA(),
  mu_scope    = ~ pb(Temp) + pb(Wind) + pb(Solar.R) + Month,
  sigma_scope = ~ pb(Temp),
  engine = "grpreg", grpreg_penalty = "grLasso",
  B = 80, pi_thr = 0.6, pre_standardize = TRUE, parallel = "auto", trace = FALSE
)
head(selection_table(fit_zaga))

## -----------------------------------------------------------------------------
if (requireNamespace('ggplot2', quietly = TRUE)) {
  print(effect_plot(fit_zaga, 'Temp', aq, what = 'mu'))
}

## -----------------------------------------------------------------------------
if (requireNamespace("nlme", quietly = TRUE)) {
  data_ok <- tryCatch({
    utils::data("Orthodont", package = "nlme", envir = environment())
    TRUE
  }, error = function(e) FALSE) 
  if (!data_ok || !exists("Orthodont", inherits = FALSE)) {
    cat("Dataset 'Orthodont' unavailable; skipping longitudinal example.\n")
  } else {
    od <- get("Orthodont", inherits = FALSE)
    rm(Orthodont)
    od$Subject <- factor(od$Subject)
    od$Sex <- factor(od$Sex)

    set.seed(30)
    fit_long <- sb_gamlss(
      distance ~ gamlss::random(Subject),     # random intercept
      data = od, family = gamlss.dist::NO(),
      mu_scope = ~ pb(age) + Sex + pb(age):Sex,      # select fixed effects
      engine = "grpreg", grpreg_penalty = "grLasso",
      B = 80, pi_thr = 0.6, pre_standardize = TRUE, parallel = "auto", trace = FALSE
    )
    head(selection_table(fit_long))

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      gg <- effect_plot(fit_long, "age", od, what = "mu")
      print(gg)
    }
  }
} else {
  cat("Packages 'nlme' not installed; skipping longitudinal example.\n")
}

