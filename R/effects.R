
#' One-variable effect plot from an sb_gamlss (or gamlss) fit
#'
#' Varies one variable and holds others at typical values (median/mode) to plot the
#' predicted parameter curve (default: mu). Uses ggplot2 if available, otherwise base.
#'
#' @param fit sb_gamlss object (or gamlss)
#' @param var character, name of the variable to vary
#' @param data original data.frame used to fit
#' @param what which parameter to predict ("mu","sigma","nu","tau")
#' @param grid number of grid points for numeric variable
#' @return a ggplot object if ggplot2 present; otherwise draws base plot and returns NULL
#' @export
effect_plot <- function(fit, var, data, what = "mu", grid = 100) {
  gfit <- if (!is.null(fit$final_fit)) fit$final_fit else fit
  stopifnot(inherits(gfit, "gamlss"))

  df_ref <- data[stats::complete.cases(data), , drop = FALSE]
  # typical values
  typical <- lapply(df_ref, function(col) {
    if (is.numeric(col)) stats::median(col, na.rm = TRUE)
    else if (is.factor(col)) levels(col)[which.max(tabulate(as.integer(col)))]
    else if (is.logical(col)) FALSE else col[1]
  })
  typical <- as.data.frame(typical, stringsAsFactors = FALSE)

  if (!var %in% names(df_ref)) stop("Variable '", var, "' not found in data.")

  # construct grid
  if (is.numeric(df_ref[[var]])) {
    rng <- stats::quantile(df_ref[[var]], c(0.02, 0.98), na.rm = TRUE)
    xs <- seq(rng[1], rng[2], length.out = grid)
  } else if (is.factor(df_ref[[var]])) {
    xs <- levels(df_ref[[var]])
  } else {
    stop("Only numeric or factor variables supported for effect_plot().")
  }

  newd <- typical[rep(1, length(xs)), , drop = FALSE]
  if (is.factor(df_ref[[var]])) newd[[var]] <- factor(xs, levels = levels(df_ref[[var]])) else newd[[var]] <- xs

  yhat <- try(stats::predict(gfit, newdata = newd, what = what, type = "response"), silent = TRUE)
  if (inherits(yhat, "try-error")) stop("Prediction failed for parameter '", what, "'.")
  newd$.fitted <- as.numeric(yhat)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- if (is.numeric(df_ref[[var]])) {
      ggplot2::ggplot(newd, ggplot2::aes(x = rlang::.data[[var]], y = rlang::.data[[".fitted"]])) +
        ggplot2::geom_line() +
        ggplot2::labs(x = var, y = paste0("fitted ", what), title = paste("Effect of", var))
    } else {
      ggplot2::ggplot(newd, ggplot2::aes(x = rlang::.data[[var]], y = rlang::.data[[".fitted"]])) +
        ggplot2::geom_point() + ggplot2::geom_line(group = 1) +
        ggplot2::labs(x = var, y = paste0("fitted ", what), title = paste("Effect of", var))
    }
    return(gg)
  } else {
    if (is.numeric(df_ref[[var]])) {
      plot(newd[[var]], newd$.fitted, type = "l", xlab = var, ylab = paste0("fitted ", what), main = paste("Effect of", var))
    } else {
      plot(newd[[var]], newd$.fitted, xlab = var, ylab = paste0("fitted ", what), main = paste("Effect of", var))
      lines(seq_along(xs), newd$.fitted)
    }
    return(invisible(NULL))
  }
}
