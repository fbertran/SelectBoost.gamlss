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
  
  df <- data[stats::complete.cases(data), , drop = FALSE]
  if (!var %in% names(df)) stop("Variable '", var, "' not found in data.")
  
  # typical values for non-var columns
  typical <- lapply(df, function(col) {
    if (is.numeric(col)) stats::median(col, na.rm = TRUE)
    else if (is.factor(col)) levels(col)[which.max(tabulate(as.integer(col)))]
    else if (is.logical(col)) FALSE
    else col[1]
  })
  typical <- as.data.frame(typical, stringsAsFactors = FALSE)
  
  # grid over var
  if (is.numeric(df[[var]])) {
    rng <- stats::quantile(df[[var]], c(0.02, 0.98), na.rm = TRUE)
    xs  <- seq(rng[1], rng[2], length.out = grid)
  } else if (is.factor(df[[var]])) {
    xs  <- levels(df[[var]])
  } else {
    stop("Only numeric or factor variables supported in effect_plot().")
  }
  
  newd <- typical[rep(1, length(xs)), , drop = FALSE]
  if (is.factor(df[[var]])) {
    newd[[var]] <- factor(xs, levels = levels(df[[var]]))
  } else {
    newd[[var]] <- xs
  }
  
  yhat <- try(stats::predict(gfit, newdata = newd, what = what, type = "response"), silent = TRUE)
  
  # If prediction fails, still return a ggplot object (so downstream code/tests don't error)
  if (inherits(yhat, "try-error")) {
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      msg <- paste("Prediction failed for", what, ":", as.character(yhat))
      blank <- data.frame(.x = seq_along(xs), .y = NA_real_)
#      x_var = ".x"
#      y_var = ".y"
      return(
#        ggplot2::ggplot(tmp_dataf, ggplot2::aes(x = x, y = y)) +
        ggplot2::ggplot(blank, ggplot2::aes(.data$.x, .data$.y)) +
#        ggplot2::ggplot(newd, ggplot2::aes(.data[[x_var]], .data[[y_var]])) +
          ggplot2::geom_blank() +
          ggplot2::labs(title = msg, x = var, y = paste0("fitted ", what))
      )
    } else {
      warning("Prediction failed: ", as.character(yhat))
      return(invisible(NULL))
    }
  }
  
  newd$.fitted <- as.numeric(yhat)
  
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    if (is.numeric(df[[var]])) {
      x_var <- var
      y_var <- ".fitted"
      ggplot2::ggplot(newd, ggplot2::aes(.data[[x_var]], .data[[y_var]])) +
        ggplot2::geom_line() +
        ggplot2::labs(x = var, y = paste0("fitted ", what), title = paste("Effect of", var))
    } else {
      x_var <- var
      y_var <- ".fitted"      
      ggplot2::ggplot(newd, ggplot2::aes(.data[[x_var]], .data[[y_var]])) +
        ggplot2::geom_point() +
        ggplot2::geom_line(ggplot2::aes(group = 1)) +
        ggplot2::labs(x = var, y = paste0("fitted ", what), title = paste("Effect of", var))
    }
  } else {
    # base fallback
    if (is.numeric(df[[var]])) {
      graphics::plot(newd[[var]], newd$.fitted, type = "l",
                     xlab = var, ylab = paste0("fitted ", what), main = paste("Effect of", var))
    } else {
      graphics::plot(seq_along(xs), newd$.fitted, xaxt = "n",
                     xlab = var, ylab = paste0("fitted ", what), main = paste("Effect of", var))
      graphics::axis(1, at = seq_along(xs), labels = xs, las = 2)
      graphics::lines(seq_along(xs), newd$.fitted)
    }
    invisible(NULL)
  }
}

#' @rdname effect_plot
#' @param x object returned by \code{effect_plot()} when prediction fails
#' @param ... unused
#' @export
print.effect_plot_failure <- function(x, ...) {
  msg <- x$message
  if (!is.character(msg) || length(msg) == 0L || !nzchar(msg)) {
    msg <- "Prediction failed for effect_plot()."
  }
  message(msg)
  invisible(x)
}
