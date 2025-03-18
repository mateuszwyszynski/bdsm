#include <RcppArmadillo.h>
#include <roptim.h>

using namespace Rcpp;
using namespace roptim;

class SEMObjective : public Functor {
public:
  SEMObjective(SEXP data_, bool exact_value_, Function SEM_likelihood_fn_, double fnscale_, 
               arma::vec parscale_ = arma::vec()) :
    data(data_),
    exact_value(exact_value_),
    SEM_likelihood_fn(SEM_likelihood_fn_),
    fnscale(fnscale_),
    parscale(parscale_) {}

  double operator()(const arma::vec &x) {
    NumericVector x_rcpp(x.begin(), x.end());
    double result = as<double>(SEM_likelihood_fn(x_rcpp, data, R_NilValue, R_NilValue, R_NilValue,
                                       R_NilValue, false, exact_value));
    // We're maximizing but roptim minimizes, so negate if fnscale is negative
    if (fnscale < 0) {
      return -result;
    } else {
      return result;
    }
  }
  
  // Override gradient method to match R's numerical gradient behavior
  arma::vec gradient(const arma::vec &x) {
    const double eps = std::sqrt(DOUBLE_EPS); // Similar to R's .Machine$double.eps^0.5
    arma::vec grad(x.n_elem);
    
    for(unsigned int i = 0; i < x.n_elem; i++) {
      double h = eps * std::max(std::abs(x(i)), 1.0);
      if (parscale.n_elem > 0) {
        h *= parscale(i); // Adjust step size by parscale if provided
      }
      
      arma::vec x_plus = x;
      x_plus(i) += h;
      double f_plus = this->operator()(x_plus);
      
      arma::vec x_minus = x;
      x_minus(i) -= h;
      double f_minus = this->operator()(x_minus);
      
      // Central difference approximation
      grad(i) = (f_plus - f_minus) / (2 * h);
    }
    
    return grad;
  }

private:
  SEXP data;
  bool exact_value;
  Function SEM_likelihood_fn;
  double fnscale;
  arma::vec parscale;
  const double DOUBLE_EPS = 2.2204460492503131e-16; // R's .Machine$double.eps
};

// [[Rcpp::export]]
List bfgs_optim(NumericVector params_no_na, Function SEM_likelihood, SEXP data, bool exact_value, List control) {
  // Extract control parameters
  double fnscale = -1.0;
  if (control.containsElementNamed("fnscale")) {
    fnscale = as<double>(control["fnscale"]);
  }
  
  // Extract parameter scaling
  arma::vec parscale;
  if (control.containsElementNamed("parscale")) {
    parscale = as<arma::vec>(control["parscale"]);
  }
  
  SEMObjective objective(data, exact_value, SEM_likelihood, fnscale, parscale);
  Roptim<SEMObjective> optimizer;
  arma::vec params_arma(params_no_na.begin(), params_no_na.size());

  // Set appropriate defaults to better match R behavior
  optimizer.control.trace = 2;      // Match R default in optimal_model_space
  optimizer.control.maxit = 1000;   // Match R default in optimal_model_space
  optimizer.control.reltol = 1e-8;  // Default R value
  optimizer.control.abstol = -INFINITY;
  optimizer.control.alpha = 1.0;    // Line search parameters to match R
  optimizer.control.beta = 0.5;
  optimizer.control.gamma = 2.0;
  
  // Override defaults with control parameters
  if (control.containsElementNamed("maxit")) {
    optimizer.control.maxit = as<int>(control["maxit"]);
  }
  if (control.containsElementNamed("abstol")) {
    optimizer.control.abstol = as<double>(control["abstol"]);
  }
  if (control.containsElementNamed("reltol")) {
    optimizer.control.reltol = as<double>(control["reltol"]);
  }
  if (control.containsElementNamed("REPORT")) {
    optimizer.control.REPORT = as<int>(control["REPORT"]);
  }
  if (control.containsElementNamed("trace")) {
    optimizer.control.trace = as<int>(control["trace"]);
  }
  if (control.containsElementNamed("parscale")) {
    optimizer.control.parscale = parscale;
  }

  optimizer.set_method("BFGS");
  
  // Set more aggressive initial line search step (to match R behavior)
  optimizer.control.lmm = 20;  // Increase memory for L-BFGS-B approximation
  
  optimizer.minimize(objective, params_arma);

  // Get results and negate value if needed for consistency
  arma::vec result_par = optimizer.par();
  NumericVector par(result_par.begin(), result_par.end());
  double value = optimizer.value();
  
  // Need to negate back to match R's behavior if we're maximizing
  if (fnscale < 0) {
    value = -value;
  }
  
  int convergence = optimizer.convergence();
  int count = optimizer.fncount();
  int grad_count = optimizer.grcount();

  // Create and return results in the same format as stats::optim
  return List::create(
    Named("par") = par,
    Named("value") = value,
    Named("counts") = List::create(
      Named("function") = count,
      Named("gradient") = grad_count
    ),
    Named("convergence") = convergence,
    Named("message") = optimizer.message()
  );
}
