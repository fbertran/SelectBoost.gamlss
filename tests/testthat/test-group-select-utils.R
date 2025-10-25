library(testthat)

test_that(".expand_terms_to_mm covers pbm and lo smoothers", {
  terms <- c("pbm(x1)", "lo(x2, span = 0.7)", "x1:x3")
  expanded <- .expand_terms_to_mm(terms, df_smooth = 4L)
  expect_true(any(grepl("splines::bs", expanded)))
  expect_true(any(grepl("splines::ns", expanded)))
  expect_true("x1:x3" %in% expanded)
})

test_that("custom converters can override expansion", {
  conv <- list(function(term, df_smooth) {
    if (grepl("special", term)) return("I(x4^2)")
    NULL
  })
  on.exit(options(SelectBoost.gamlss.term_converters = NULL), add = TRUE)
  options(SelectBoost.gamlss.term_converters = conv)
  expanded <- .expand_terms_to_mm(c("special(x4)"), df_smooth = 5L)
  expect_identical(expanded, "I(x4^2)")
})