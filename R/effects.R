
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
#' @return A ggplot object if ggplot2 present; otherwise draws a base plot and returns
#'   NULL. When prediction fails, returns an object of class
#'   \code{effect_plot_failure} that carries diagnostic information and prints a
#'   message describing the failure.
#' @export
effect_plot <- function(fit, var, data, what = "mu", grid = 100) {
  call <- match.call()
  gfit <- if (!is.null(fit$final_fit)) fit$final_fit else fit
  stopifnot(inherits(gfit, "gamlss"))

  df_ref <- data[stats::complete.cases(data), , drop = FALSE]
  if (nrow(df_ref) == 0L) {
    stop("No complete cases available to build reference profile for effect_plot().")
  }
  
  # typical values (preserve column classes for predict.gamlss)
  typical <- df_ref[1, , drop = FALSE]
  for (nm in names(df_ref)) {
    col <- df_ref[[nm]]
    if (is.numeric(col)) {
      typical[[nm]] <- stats::median(col, na.rm = TRUE)
    } else if (is.factor(col)) {
      lev <- levels(col)
      if (!length(lev)) next
      idx <- which.max(tabulate(as.integer(col[!is.na(col)]), nbins = length(lev)))
      if (!length(idx) || is.na(idx) || idx < 1L) idx <- 1L
      typical[[nm]] <- factor(lev[idx], levels = lev)
    } else if (is.logical(col)) {
      vals <- col[!is.na(col)]
      typical[[nm]] <- if (!length(vals)) FALSE else sum(vals) >= ceiling(length(vals) / 2)
    } else {
      nz <- col[!is.na(col)]
      if (length(nz)) typical[[nm]] <- nz[1] else typical[[nm]] <- col[1]
    }
  }
  
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
  if (inherits(yhat, "try-error")) {
    cond <- attr(yhat, "condition")
    cond_msg <- if (!is.null(cond)) conditionMessage(cond) else ""
    
    row_preds <- numeric(nrow(newd))
    row_error <- NULL
    row_fail <- NA_integer_
    for (i in seq_len(nrow(newd))) {
      one <- try(
        stats::predict(gfit, newdata = newd[i, , drop = FALSE], what = what, type = "response"),
        silent = TRUE
      )
      if (inherits(one, "try-error")) {
        row_error <- attr(one, "condition")
        row_fail <- i
        break
      }
      row_preds[i] <- as.numeric(one)
    }
    
    if (!is.null(row_error)) {
      base_extra <- if (!identical(cond_msg, "")) paste0(": ", cond_msg) else ""
      fallback_msg <- conditionMessage(row_error)
      fallback_extra <- if (!identical(fallback_msg, "")) paste0(        " (fallback failed",
        if (!is.na(row_fail)) paste0(" at row ", row_fail) else "",
        ": ",
        fallback_msg,
        ")"
      ) else if (!is.na(row_fail)) paste0(" (fallback failed at row ", row_fail, ")") else " (fallback failed)"
      
      failed_newd <- newd
      failed_newd$.fitted <- NA_real_
      
      failure <- structure(
        list(
          message = paste0("Prediction failed for parameter '", what, "'", base_extra, fallback_extra, "."),
          parameter = what,
          newdata = failed_newd,
          call = call,
          original_error = cond,
          fallback_error = row_error
        ),
        class = "effect_plot_failure"
      )
      return(failure)
    }
    yhat <- row_preds
    }
  newd$.fitted <- as.numeric(yhat)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    gg <- if (is.numeric(df_ref[[var]])) {
      ggplot2::ggplot(newd, ggplot2::aes_string(x = var, y = ".fitted")) +
        ggplot2::geom_line() +
        ggplot2::labs(x = var, y = paste0("fitted ", what), title = paste("Effect of", var))
    } else {
      ggplot2::ggplot(newd, ggplot2::aes_string(x = var, y = ".fitted")) +
        ggplot2::geom_point() + ggplot2::geom_line(ggplot2::aes(group = 1)) +
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
