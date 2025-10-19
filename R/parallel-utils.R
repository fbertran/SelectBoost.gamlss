
# Internal parallel helpers

.sb_parallel_lapply <- function(B, FUN, parallel = c("none","auto","multisession","multicore"), workers = NULL) {
  parallel <- match.arg(parallel)
  if (parallel == "none") return(lapply(seq_len(B), FUN))

  if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
    warning("future/future.apply not available; falling back to sequential lapply.")
    return(lapply(seq_len(B), FUN))
  }

  old_plan <- NULL
  undo <- function() {}
  if (parallel == "auto") {
    # set a temporary plan if none is set
    old_plan <- try(future::plan(), silent = TRUE)
    # if querying plan fails or returns something, still set a temp plan
    future::plan(future::multisession, workers = if (is.null(workers)) max(1, future::availableCores() - 1) else workers)
    undo <- function() {
      if (!inherits(old_plan, "try-error")) {
        future::plan(old_plan)
      } else {
        future::plan(future::sequential)
      }
    }
  } else if (parallel == "multisession") {
    future::plan(future::multisession, workers = workers %||% max(1, future::availableCores() - 1))
  } else if (parallel == "multicore") {
    future::plan(future::multicore, workers = workers %||% max(1, future::availableCores() - 1))
  }

  on.exit(undo(), add = TRUE)
  future.apply::future_lapply(seq_len(B), FUN)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
