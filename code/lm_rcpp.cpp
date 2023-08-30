#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec fast_lm(const arma::vec& Zi, const arma::mat& X, const arma::mat& Y, const arma::mat& C) {
  int n = Zi.n_elem;
  arma::mat XTY = X.t() * Y;
  arma::mat XTC = X.t() * C;
  arma::mat XTX = X.t() * X;
  
  arma::mat inv_term = inv(XTX);
  
  arma::vec beta = inv_term * (XTY + XTC);
  
  return beta;
}