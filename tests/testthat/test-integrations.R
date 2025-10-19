test_that("integration wrappers exist", {
  skip_on_cran()
  expect_true(exists("SelectBoost_gamlss"))
  expect_true(exists("sb_gamlss_c0_grid"))
  expect_true(exists("autoboost_gamlss"))
  expect_true(exists("fastboost_gamlss"))
})
