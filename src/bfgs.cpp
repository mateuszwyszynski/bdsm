#include <RcppArmadillo.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

// Custom BFGS implementation from scratch
// This doesn't rely on external optimization libraries
// [[Rcpp::export]]
List bfgs_optim(NumericVector params_no_na, Function SEM_likelihood, SEXP data, bool exact_value, List control) {
  // Extract control parameters
  double fnscale = -1.0; // Default for maximization in R
  int maxit = 1000;
  double reltol = 1e-8;
  int trace = 2;
  int report = 10;
  bool use_r = false; // Fall back to R's optim for testing?
  
  if (control.containsElementNamed("fnscale")) {
    fnscale = as<double>(control["fnscale"]);
  }
  if (control.containsElementNamed("maxit")) {
    maxit = as<int>(control["maxit"]);
  }
  if (control.containsElementNamed("reltol")) {
    reltol = as<double>(control["reltol"]);
  }
  if (control.containsElementNamed("trace")) {
    trace = as<int>(control["trace"]);
  }
  if (control.containsElementNamed("REPORT")) {
    report = as<int>(control["REPORT"]);
  }
  
  // Fall back to R's implementation for testing if needed
  if (use_r) {
    Environment stats("package:stats");
    Function optim = stats["optim"];
    
    List result = optim(
      Named("par") = params_no_na,
      Named("fn") = SEM_likelihood,
      Named("data") = data,
      Named("exact_value") = exact_value,
      Named("method") = "BFGS",
      Named("control") = control
    );
    
    return result;
  }
  
  // Initialize parameter vector
  int n = params_no_na.size();
  arma::vec x(params_no_na.begin(), n);
  
  // Counters
  int f_count = 0;
  int g_count = 0;
  
  // Objective function
  auto fn = [&](const arma::vec& pars) {
    NumericVector pars_r(pars.begin(), pars.end());
    f_count++;
    double val = as<double>(SEM_likelihood(pars_r, data, R_NilValue, R_NilValue, 
                                         R_NilValue, R_NilValue, false, exact_value));
    // For maximization problems, negate
    return fnscale < 0 ? -val : val;
  };
  
  // Numerical gradient function
  auto grad = [&](const arma::vec& pars) {
    g_count++;
    arma::vec g(n);
    double eps = std::sqrt(2.2e-16); // ~sqrt(.Machine$double.eps)
    
    double f0 = fn(pars);
    
    for (int i = 0; i < n; i++) {
      arma::vec x_plus = pars;
      double h = eps * std::max(1.0, std::abs(pars(i)));
      x_plus(i) += h;
      
      double f_plus = fn(x_plus);
      g(i) = (f_plus - f0) / h;
    }
    
    return g;
  };
  
  // Initial function value and gradient
  double f = fn(x);
  arma::vec g = grad(x);
  
  // Initial Hessian approximation (identity matrix)
  arma::mat H = arma::eye(n, n);
  
  // Report initial value
  if (trace > 0) {
    Rcout << "initial  value " << (fnscale < 0 ? -f : f) << std::endl;
  }
  
  bool converged = false;
  int iter = 0;
  
  // Main BFGS loop
  while (iter < maxit) {
    // Compute search direction
    arma::vec direction = -H * g;
    
    // Line search to find step size
    double alpha = 1.0;
    double beta = 0.5;  // Backtracking factor
    double c1 = 1e-4;   // Sufficient decrease condition constant
    
    // Initial values
    double f0 = f;
    arma::vec x0 = x;
    arma::vec g0 = g;
    
    // Backtracking line search
    double slope = arma::dot(g, direction);
    if (slope >= 0) {
      // If direction is not descent, reset Hessian to identity
      H = arma::eye(n, n);
      direction = -g;
      slope = arma::dot(g, direction);
    }
    
    // Keep reducing step size until we find a sufficient decrease
    while (true) {
      x = x0 + alpha * direction;
      f = fn(x);
      
      if (f <= f0 + c1 * alpha * slope || alpha < 1e-10) {
        break;
      }
      
      alpha *= beta;
    }
    
    // Compute new gradient
    g = grad(x);
    
    // BFGS update of Hessian approximation
    arma::vec s = alpha * direction;
    arma::vec y = g - g0;
    
    double rho = 1.0 / std::max(1e-12, arma::dot(y, s));
    
    if (std::abs(rho) < 1e10) { // Protection against numerical issues
      arma::mat identity = arma::eye(n, n);
      H = (identity - rho * s * y.t()) * H * (identity - rho * y * s.t()) + rho * s * s.t();
    }
    
    // Check convergence
    double gnorm = arma::norm(g, 2);
    
    if (gnorm < reltol || std::abs(f - f0) < reltol * (1.0 + std::abs(f0))) {
      converged = true;
      break;
    }
    
    // Reporting
    if (trace > 0 && iter % report == 0) {
      Rcout << "iter " << iter << " value " << (fnscale < 0 ? -f : f) << std::endl;
    }
    
    iter++;
  }
  
  // Final value after adjustment for maximization/minimization
  double value = fnscale < 0 ? -f : f;
  
  // Final parameters
  NumericVector par(x.begin(), x.end());
  
  // Reporting
  if (trace > 0) {
    Rcout << "final  value " << value << std::endl;
    if (converged) {
      Rcout << "converged" << std::endl;
    } else {
      Rcout << "max iterations reached" << std::endl;
    }
  }
  
  // Return in the format expected by optim()
  return List::create(
    Named("par") = par,
    Named("value") = value,
    Named("counts") = List::create(
      Named("function") = f_count,
      Named("gradient") = g_count
    ),
    Named("convergence") = converged ? 0 : 1,
    Named("message") = converged ? "converged" : "max iterations reached"
  );
}
