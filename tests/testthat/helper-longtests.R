skip_long_tests <- function() {
  opt <- getOption("SelectBoost.gamlss.run_long_tests", default = FALSE)
  env <- Sys.getenv("RUN_LONG_TESTS", unset = "false")
  if (!(isTRUE(opt) || tolower(env) %in% c("1","true","yes"))) {
    testthat::skip("Long tests disabled (set options(SelectBoost.gamlss.run_long_tests=TRUE) or RUN_LONG_TESTS=true).")
  }
}
