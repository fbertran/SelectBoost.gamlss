
# Internal helpers backed by RcppArmadillo

cs_scale_mat <- function(df, vars) {
  X <- as.matrix(df[vars])
  out <- cs_scale(X)
  list(center = as.numeric(out$center), scale = as.numeric(out$scale),
       X = out$X, vars = vars)
}

fast_cor <- function(df, vars) {
  X <- as.matrix(df[vars])
  cs_cor(X)
}
