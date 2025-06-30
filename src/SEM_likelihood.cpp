#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Function declarations from SEM_matrices.cpp
arma::mat residual_maker_matrix(const arma::mat &m);
Rcpp::List sem_B_matrix(double alpha, int periods_n, Rcpp::Nullable<arma::vec> beta = R_NilValue);
arma::mat sem_C_matrix(double alpha, double phi_0, int periods_n,
                       Rcpp::Nullable<arma::vec> beta = R_NilValue,
                       Rcpp::Nullable<arma::vec> phi_1 = R_NilValue);
Rcpp::List sem_sigma_matrix(double err_var, const arma::vec &dep_vars,
                            Rcpp::Nullable<arma::vec> phis = R_NilValue,
                            Rcpp::Nullable<arma::vec> psis = R_NilValue);

// arma::mat operator*(const arma::mat &A, const arma::mat &B)
// {
//   static Rcpp::Function rmul("%*%");
//   return Rcpp::as<arma::mat>(rmul(A, B));
// }

// arma::mat operator-(const arma::mat &A, const arma::mat &B)
// {
//   static Rcpp::Function rsub("-");
//   return Rcpp::as<arma::mat>(rsub(A, B));
// }

inline arma::mat crossprod(const arma::mat &M)
{
  // static Rcpp::Function rcrossprod("crossprod");
  // return Rcpp::as<arma::mat>(rcrossprod(M));
  return M.t() * M;
}

// double rsum(const arma::mat &M)
// {
//   static Rcpp::Function fn("sum");
//   return Rcpp::as<arma::vec>(fn(M)).at(0);
// }

// [[Rcpp::export]]
SEXP sem_likelihood_calculate(double alpha, double phi_0, double err_var,
                              const arma::vec &dep_vars,
                              const arma::mat &Y1,
                              const arma::mat &Y2,
                              const arma::mat &cur_Z,
                              Rcpp::Nullable<arma::mat> cur_Y2 = R_NilValue,
                              Rcpp::Nullable<arma::vec> beta = R_NilValue,
                              Rcpp::Nullable<arma::vec> phi_1 = R_NilValue,
                              Rcpp::Nullable<arma::vec> phis = R_NilValue,
                              Rcpp::Nullable<arma::vec> psis = R_NilValue,
                              bool per_entity = false, bool exact_value = true)
{
  Function det("det");
  Function log("log");
  Function solve("solve");

  arma::mat res_maker_matrix = residual_maker_matrix(cur_Z);

  int n_entities = Y1.n_rows;
  int periods_n = dep_vars.n_elem;
  int tot_regressors_n = Y2.n_cols / (periods_n - 1);
  int lin_related_regressors_n = beta.isNotNull() ? as<arma::vec>(beta).n_elem : 0;

  Rcpp::List B = sem_B_matrix(alpha, periods_n, beta);
  arma::mat C = sem_C_matrix(alpha, phi_0, periods_n, beta, phi_1);
  Rcpp::List S = sem_sigma_matrix(err_var, dep_vars, phis, psis);

  arma::mat B1 = as<arma::mat>(B[0]);
  arma::mat S1 = as<arma::mat>(S[0]);
  arma::mat S2 = as<arma::mat>(S[1]);

  arma::mat U1;
  if (lin_related_regressors_n == 0)
  {
    U1 = trans(B1 * trans(Y1) - C * trans(cur_Z));
  }
  else
  {
    arma::mat B2 = as<arma::mat>(B[1]);
    U1 = trans(B1 * trans(Y1) + B2 * trans(as<arma::mat>(cur_Y2)) - C * trans(cur_Z));
  }
  // arma::mat S11_inverse = inv(S1);
  arma::mat S11_inverse = as<arma::mat>(solve(S1));
  arma::mat M = Y2 - U1 * S11_inverse * S2;
  arma::mat H = trans(M) * res_maker_matrix * M;

  double gaussian_normalization_const = as<double>(log(2 * M_PI)) *
                                        n_entities * (periods_n + (periods_n - 1) * tot_regressors_n) / 2.0;
  double trace_simplification_term =
      0.5 * n_entities * (periods_n - 1) * tot_regressors_n;

  double likelihood = -n_entities / 2.0 * as<double>(log(as<arma::vec>(det(S1)) * as<arma::vec>(det(H / n_entities))));

  if (isnan(likelihood))
  {
    return wrap(likelihood);
  }

  if (exact_value)
  {
    likelihood -= gaussian_normalization_const + trace_simplification_term;
  }

  arma::vec result;
  if (!per_entity)
  {
    likelihood -= 0.5 * sum(diagvec(S11_inverse * crossprod(U1)));
    return wrap(likelihood);
  }
  else
  {
    arma::vec per_entity_likelihood = arma::vec(n_entities);
    per_entity_likelihood.fill(likelihood / n_entities);
    arma::vec diagonal_terms = -0.5 * diagvec(U1 * S11_inverse * trans(U1));
    per_entity_likelihood += diagonal_terms;
    return wrap(per_entity_likelihood);
  }
}
