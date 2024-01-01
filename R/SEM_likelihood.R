SEM_params_to_list <- function(params, periods_n, tot_regressors_n,
                               in_regressors_n) {
  phis_n <- tot_regressors_n*(periods_n - 1)
  psis_n <- tot_regressors_n*periods_n*(periods_n - 1)/2

  alpha <- params[1]
  if (tot_regressors_n == 0) {
    beta <- c()
    phi_1 <- c()
    phis <- c()
    psis <- c()
  } else {
    if (in_regressors_n == 0) {
      beta <- c()
      phi_1 <- c()
    } else {
      beta <- params[2:(1 + in_regressors_n)]
      phi_1 <- params[(3 + in_regressors_n):(2 + 2*in_regressors_n)]
    }
    phis <-
      params[(4 + 2*in_regressors_n + periods_n):(3 + 2*in_regressors_n + periods_n + phis_n)]
    psis <-
      params[(4 + 2*in_regressors_n + periods_n + phis_n):(3 + 2*in_regressors_n + periods_n + phis_n + psis_n)]
  }
  phi_0 <- params[2 + in_regressors_n]
  err_var <- params[3 + 2*in_regressors_n]
  dep_vars <-
    params[(4 + 2*in_regressors_n):(3 + 2*in_regressors_n + periods_n)]

  list(alpha = alpha, phi_0 = phi_0, err_var = err_var, dep_vars = dep_vars,
       beta = beta, phi_1 = phi_1, phis = phis, psis = psis)
}

#' Likelihood for the SEM model
#'
#' @param params Parameters describing the model. Can be either a vector or a
#' list with named parameters. See 'Details'
#' @param data Data for the likelihood computations. Can be either a list of
#' matrices or a dataframe. If the dataframe, additional parameters are
#' required to build the matrices within the function.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestampsg
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#' @param regressors Which subset of columns should be used as regressors.
#' @param in_regressors Which subset of columns should be used as regressors
#' for the current model. In other words \code{regressors} are the total set of
#' regressors and \code{in_regressors} are the ones for which linear relation
#' is not set to zero for a given model.
#' @param per_entity Whether to compute overall likelihood or a vector of
#' likelihoods with per entity value
#' @param projection_matrix_const Wheter the residual maker matrix (and so
#' the projection matrix) should be computed for each model separately.
#' \code{TRUE} means that the matrix will be the same for all models
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}).
#' Currently \code{TRUE} adds: 1. a normalization constant coming from Gaussian
#' distribution, 2. a term disappearing during likelihood simplification in
#' Likelihood-based Estimation of Dynamic Panels with Predetermined Regressors
#' by Moral-Benito (see Appendix A.1). The latter happens when transitioning
#' from equation (47) to equation (48), in step 2: the term trace{HG_22} is
#' dropped, because it can be assumed to be constant from Moral-Benito
#' perspective. To get the exact value of the likelihood we have to take this
#' term into account.
#'
#' @details
#' The \code{params} argument is a list that should contain the following
#' components:
#'
#' \code{alpha} scalar value which determines linear dependence on lagged
#' dependent variable
#'
#' \code{phi_0} scalar value which determines linear dependence on the value
#' of dependent variable at the lowest time stamp
#'
#' \code{err_var} scalar value which determines classical error component
#' (Sigma11 matrix, sigma_epsilon^2)
#'
#' \code{dep_vars} double vector of length equal to the number of time stamps
#' (i.e. time stamps greater than or equal to the second lowest time stamp)
#'
#' \code{beta} double vector which determines the linear dependence on
#' regressors different than the lagged dependent variable; The vector should
#' have length equal to the number of regressors.
#'
#' \code{phi_1} double vector which determines the linear dependence on
#' initial values of regressors different than the lagged dependent variable;
#' The vector should have length equal to the number of regressors.
#'
#' \code{phis} double vector which together with \code{psis} determines upper
#' right and bottom left part of the covariance matrix; The vector should have
#' length equal to the number of regressors times number of time stamps minus 1,
#' i.e. \code{regressors_n * (periods_n - 1)}
#'
#' \code{psis} double vector which together with \code{psis} determines upper
#' right and bottom left part of the covariance matrix; The vector should have
#' length equal to the number of regressors times number of time stamps minus 1
#' times number of time stamps divided by 2, i.e.
#' \code{regressors_n * (periods_n - 1) * periods_n / 2}
#'
#'
#' @return
#' The value of the likelihood for SEM model (or a part of interest of the
#' likelihood)
#'
#' @export
#'
#' @examples
SEM_likelihood <- function(params, data, timestamp_col = NULL,
                           entity_col = NULL, dep_var_col = NULL,
                           regressors = NULL, in_regressors = NULL,
                           per_entity = FALSE, projection_matrix_const = TRUE,
                           exact_value = TRUE) {
  if (is.list(params) && is.list(data)) {
    alpha <- params$alpha
    phi_0 <- params$phi_0
    err_var <- params$err_var
    dep_vars <- params$dep_vars
    beta <- if(is.null(params$beta)) c() else params$beta
    phi_1 <- if(is.null(params$phi_1)) c() else params$phi_1
    phis <- if(is.null(params$phis)) c() else params$phis
    psis <- if(is.null(params$psis)) c() else params$psis

    Y1 <- data$Y1
    Y2 <- data$Y2
    cur_Y2 <- data$cur_Y2
    Z <- data$Z
    res_maker_matrix <- if (projection_matrix_const) {
      data$res_maker_matrix
    } else {
      residual_maker_matrix(Z)
    }

    n_entities <- nrow(Z)
    periods_n <- length(dep_vars)
    tot_regressors_n <- ncol(data$Y2) / (periods_n - 1)
    in_regressors_n <- length(beta)

    B <- SEM_B_matrix(alpha, periods_n, beta)
    C <- SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
    S <- SEM_sigma_matrix(err_var, dep_vars, phis, psis)

    U1 <- if (in_regressors_n == 0) {
      t(tcrossprod(B[[1]], Y1) - tcrossprod(C, Z))
    } else {
      t(tcrossprod(B[[1]], Y1) + tcrossprod(B[[2]], cur_Y2) - tcrossprod(C, Z))
    }
    S11_inverse <- solve(S[[1]])
    M <- Y2 - U1 %*% S11_inverse %*% S[[2]]
    H <- crossprod(M, res_maker_matrix) %*% M

    gaussian_normalization_const <- log(2 * pi) *
      n_entities * (periods_n + (periods_n-1) * tot_regressors_n) / 2
    trace_simplification_term <-
      1/2 * n_entities * (periods_n - 1) * tot_regressors_n

    likelihood <- -n_entities/2 * log(det(S[[1]]) * det(H/n_entities))

    if(exact_value) {
      likelihood <- likelihood -
        gaussian_normalization_const - trace_simplification_term
    }

    likelihood <- if(!per_entity) {
      likelihood - 1/2 * sum(diag(S11_inverse %*% crossprod(U1)))
    } else {
      likelihood / n_entities  - 1/2 * diag(U1 %*% S11_inverse %*% t(U1))
    }
  } else {
    if (!is.list(data)) {
      Y1 <- SEM_dep_var_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        dep_var_col = dep_var_col
      )
      Y2 <- SEM_regressors_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        regressors = regressors
      )
      cur_Y2 <- SEM_regressors_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        regressors = in_regressors
      )
      cur_Z <- exogenous_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        dep_var_col = dep_var_col, regressors_subset = in_regressors
      )
      Z <- exogenous_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        dep_var_col = dep_var_col, regressors_subset = regressors
      )
      res_maker_matrix <- residual_maker_matrix(Z)

      data = list(Y1 = Y1, Y2 = Y2, cur_Y2 = cur_Y2, Z = cur_Z,
                  res_maker_matrix = res_maker_matrix)
    }
    if (!is.list(params)) {
      periods_n <- ncol(data$Y1)
      tot_regressors_n <- ncol(data$Y2) / (periods_n - 1)
      in_regressors_n <- if(is.null(data$cur_Y2)) {
        0
      } else {
        ncol(data$cur_Y2) / (periods_n - 1)
      }
      params <- SEM_params_to_list(params, periods_n = periods_n,
                                   tot_regressors_n = tot_regressors_n,
                                   in_regressors_n = in_regressors_n)
    }
    likelihood <-
      SEM_likelihood(params = params, data = data, per_entity = per_entity,
                     projection_matrix_const = projection_matrix_const,
                     exact_value = exact_value)
  }
  likelihood
}
