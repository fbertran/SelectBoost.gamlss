# Internal parallel helpers

# internal: lapply with optional future backend, restoring plan exactly
.sb_parallel_lapply <- function(X, FUN, ..., parallel = c("none", "auto", "multisession", "multicore"), workers = NULL) {
  parallel <- match.arg(parallel)
  
  if (parallel == "none") {
    return(lapply(X, FUN, ...))
  }
  
  if (!requireNamespace("future.apply", quietly = TRUE) || !requireNamespace("future", quietly = TRUE)) {
    stop("Parallel mode requires 'future' and 'future.apply'.")
  }
  
  # Save & restore exact plan
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  # Choose strategy
  strategy <- switch(parallel,
                     auto        = if (.Platform$OS.type == "windows") future::multisession else future::multicore,
                     multisession= future::multisession,
                     multicore   = future::multicore
  )
  
  # Determine workers
  w <- tryCatch(future::availableCores(), error = function(e) 1L)
  if (!is.null(workers)) w <- as.integer(workers)
  w <- max(1L, as.integer(w))
  
  future::plan(strategy, workers = w)
  future.apply::future_lapply(X, FUN, ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
