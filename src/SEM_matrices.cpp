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

//' Coefficients matrix for SEM representation
//'
//' Create coefficients matrix for Simultaneous Equations Model (SEM)
//' representation.
//'
//' @param alpha numeric
//' @param periods_n integer
//' @param beta numeric vector. Default is c() for no regressors case.
//'
//' @return List with two matrices B11 and B12
//' @export
//'
//' @examples
//' sem_B_matrix(3, 4, 4:6)
// [[Rcpp::export]]
Rcpp::List sem_B_matrix(double alpha, int periods_n, Rcpp::Nullable<arma::vec> beta = R_NilValue)
{
  // Create B11 matrix
  arma::mat B11 = arma::eye(periods_n, periods_n);

  // Create alpha_matrix (diagonal matrix with -alpha values)
  arma::mat alpha_matrix = arma::diagmat(arma::vec(periods_n - 1).fill(-alpha));

  // Add alpha_matrix to B11[2:periods_n, 1:(periods_n-1)]
  B11.submat(1, 0, periods_n - 1, periods_n - 2) =
      B11.submat(1, 0, periods_n - 1, periods_n - 2) + alpha_matrix;

  // Create B12 matrix
  arma::mat B12;
  bool has_beta = false;

  if (beta.isNotNull())
  {
    arma::vec beta_vec = Rcpp::as<arma::vec>(beta);
    if (beta_vec.n_elem > 0)
    {
      has_beta = true;
      int regressors_n = beta_vec.n_elem;
      int n_cols = regressors_n * (periods_n - 1);

      // Initialize B12 with zeros
      B12 = arma::zeros(periods_n, n_cols);

      // Fill beta_matrix for rows 2 to periods_n
      for (int row_ind = 1; row_ind < periods_n; row_ind++)
      {
        int n_zeros_front = (row_ind - 1) * regressors_n;

        // Place -beta values in the appropriate position
        for (int j = 0; j < regressors_n; j++)
        {
          B12(row_ind, n_zeros_front + j) = -beta_vec(j);
        }
      }
    }
  }

  // Return list
  if (has_beta)
  {
    return Rcpp::List::create(B11, B12);
  }
  else
  {
    return Rcpp::List::create(B11, R_NilValue);
  }
}

//' Coefficients matrix for initial conditions
//'
//' Create matrix for Simultaneous Equations Model (SEM)
//' representation with coefficients placed next to initial values
//' of regressors, dependent variable and country-specific time-invariant
//' variables.
//'
//' @param alpha numeric
//' @param phi_0 numeric
//' @param periods_n numeric
//' @param beta numeric vector. Default is c() for no regressors case.
//' @param phi_1 numeric vector. Default is c() for no regressors case.
//'
//' @return matrix
//' @export
//'
//' @examples
//' alpha <- 9
//' phi_0 <- 19
//' beta <- 11:15
//' phi_1 <- 21:25
//' periods_n <- 4
//' sem_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
// [[Rcpp::export]]
arma::mat sem_C_matrix(double alpha, double phi_0, int periods_n, 
                       Rcpp::Nullable<arma::vec> beta = R_NilValue, 
                       Rcpp::Nullable<arma::vec> phi_1 = R_NilValue)
{
  // Create C1 matrix - column vector with phi_0 repeated periods_n times
  arma::mat C1(periods_n, 1);
  C1.fill(phi_0);
  
  // Add alpha to the first element
  C1(0, 0) = C1(0, 0) + alpha;
  
  // Handle beta and phi_1 if provided
  if (beta.isNotNull()) {
    arma::vec beta_vec = Rcpp::as<arma::vec>(beta);
    if (beta_vec.n_elem > 0) {
      if (phi_1.isNotNull()) {
        arma::vec phi_1_vec = Rcpp::as<arma::vec>(phi_1);
        
        // Create col2 matrix - phi_1 repeated for each row
        arma::mat col2(periods_n, phi_1_vec.n_elem);
        for (int i = 0; i < periods_n; i++) {
          col2.row(i) = phi_1_vec.t();
        }
        
        // Add beta to the first row
        for (int j = 0; j < beta_vec.n_elem && j < phi_1_vec.n_elem; j++) {
          col2(0, j) = col2(0, j) + beta_vec(j);
        }
        
        // Combine C1 and col2
        C1 = arma::join_rows(C1, col2);
      }
    }
  }
  
  return C1;
}
