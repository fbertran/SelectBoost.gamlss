#include <cmath>   // for std::isfinite
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// center/scale numeric matrix (column-wise), return list(center, scale, X)
// [[Rcpp::export]]
List cs_scale(arma::mat X) {
  arma::rowvec mu = arma::mean(X, 0);
  arma::rowvec sd = arma::stddev(X, 0, 0); // normalize by n-1
  sd.transform( [](double v){ return (v == 0.0 || !std::isfinite(v)) ? 1.0 : v; } );
  X.each_row() -= mu;
  X.each_row() /= sd;
  return List::create(_["center"]=mu, _["scale"]=sd, _["X"]=X);
}

// fast Pearson correlation (like cor(X))
// [[Rcpp::export]]
arma::mat cs_cor(const arma::mat& X) {
  if (X.n_rows < 2) {
    arma::mat S(X.n_cols, X.n_cols);
    S.fill(NA_REAL);
    return S;
  }
  arma::rowvec mu = arma::mean(X,0);
  arma::mat C = X.each_row() - mu;
  arma::mat S = (C.t() * C) / (X.n_rows - 1.0);
  arma::vec d = arma::sqrt(S.diag());
  for (arma::uword i=0;i<S.n_rows;++i)
    for (arma::uword j=0;j<S.n_cols;++j)
      S(i,j) = (d(i)==0 || d(j)==0) ? NA_REAL : S(i,j)/(d(i)*d(j));
  return S;
}

// thresholded positive trapezoidal area (normalized)
// [[Rcpp::export]]
double pos_trapz_norm(const arma::vec& c0, const arma::vec& prop, double thr) {
  if (c0.n_elem < 2) {
    arma::vec p = arma::clamp(prop - thr, 0.0, 1.0);
    return arma::mean(p);
  }
  arma::vec p = arma::clamp(prop - thr, 0.0, 1.0);
  double num = 0.0;
  for (arma::uword i=0; i < c0.n_elem-1; ++i) {
    num += 0.5 * (p[i] + p[i+1]) * (c0[i+1] - c0[i]);
  }
  double den = c0.max() - c0.min();
  return (den > 0.0) ? (num / den) : arma::mean(p);
}
