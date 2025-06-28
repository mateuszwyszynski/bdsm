#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

//' Residual Maker Matrix
//'
//' Create residual maker matrix from a given matrix \code{m}. See article about
//' \href{https://en.wikipedia.org/wiki/Projection_matrix}{projection matrix} on
//' the Wikipedia.
//'
//' @param m Matrix
//'
//' @return
//' M x M matrix where M is the number of rows in the \code{m} matrix.
//' @export
//'
//' @examples
//' residual_maker_matrix(matrix(c(1,2,3,4), nrow = 2))
// [[Rcpp::export]]
arma::mat residual_maker_matrix(const arma::mat &m)
{
  static Function crossprod("crossprod");
  static Function solve("solve");

  // Equivalent to R: m %*% solve(crossprod(m)) %*% t(m)
  SEXP crossprod_result = crossprod(m);
  SEXP solve_result = solve(crossprod_result);
  arma::mat solve_mat = Rcpp::as<arma::mat>(solve_result);

  arma::mat proj_matrix = m * solve_mat * m.t();
  return arma::eye(m.n_rows, m.n_rows) - proj_matrix;
}
