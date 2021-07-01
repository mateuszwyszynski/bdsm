#' Matrix with data for regressors for SEM representation
#'
#' Create matrix which contains data for regressors used in the model in
#' Simultaneous Equations Model (SEM) representation. This matrix is then used
#' to compute likelihood for SEM analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods.
#' @param entity_col Columns which determines entities (e.g. countires, people)
#' @param start_time First time period. Only time periods greater than
#' \code{start_time} will be considered in the resulting matrix
#' @param regressors_subset Which subset of columns should be used as
#' regressors. If \code{NULL} (default) then all remaining columns regressors
#' will be used as regressors. For now columns have to be passed as list of
#' strings containg column names.
#'
#' @return
#' Matrix of size N x T*k where N is the number of entities considered, T is the
#' number of periods greater than \code{start_time} and k is the number of
#' chosen regressors
#' @importFrom rlang .data
#' @export
#'
#' @examples
SEM_regressors_matrix <- function(df, timestamp_col, entity_col, start_time,
                                  regressors_subset = NULL) {
  df %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, regressors_subset) %>%
    dplyr::filter({{ timestamp_col }} > start_time) %>%
    tidyr::pivot_wider(
      names_from = {{ timestamp_col }},
      values_from = !{{ entity_col }} & !{{ timestamp_col }}
      ) %>%
    dplyr::select(!{{ entity_col }}) %>%
    dplyr::select(order(as.numeric(gsub("[^0-9]+", "", colnames(.data))))) %>%
    as.matrix()
}

#' Coefficients matrix for SEM representation
#'
#' Create coefficients matrix for Simultaneous Equations Model (SEM)
#' representation.
#'
#' @param alpha numeric
#' @param periods_n integer
#' @param beta numeric vector. Default is c() for no regressors case.
#'
#' @return List with two matrices B11 and B12
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' SEM_B_matrix(3, 4, 4:6)
SEM_B_matrix <- function(alpha, periods_n, beta = c()) {
  alpha_matrix <- diag(rep(-alpha, periods_n-1))
  B11 <- diag(periods_n)
  B11[2:periods_n, 1:(periods_n - 1)] <-
    B11[2:periods_n, 1:(periods_n - 1)] + alpha_matrix

  B12 <- if (length(beta) != 0) {
    regressors_n <- length(beta)
    n_cols <- regressors_n*(periods_n - 1)
    beta_row <- function(row_ind, beta, n_cols, regressors_n) {
      n_zeros_front <- (row_ind-1)*regressors_n
      c(rep(0, n_zeros_front), -beta,
        rep(0, n_cols - n_zeros_front - regressors_n))
    }
    beta_matrix <- 1:(periods_n-1) %>%
      sapply(beta_row, beta = beta,
             n_cols = n_cols, regressors_n = regressors_n) %>% t()
    rbind(optimbase::zeros(1, n_cols),
          beta_matrix)
  } else {
    NULL
  }
  list(B11, B12)
}

#' Coefficients matrix for initial conditions
#'
#' Create matrix for Simultaneous Equations Model (SEM)
#' representation with coefficients placed next to initial values
#' of regressors, dependent variable and country-specific time-invariant
#' variables.
#'
#' @param alpha numeric
#' @param phi_0 numeric
#' @param periods_n numeric
#' @param beta numeric vector. Default is c() for no regressors case.
#' @param phi_1 numeric vector. Default is c() for no regressors case.
#'
#' @return matrix
#' @export
#'
#' @examples
#' alpha <- 9
#' phi_0 <- 19
#' beta <- 11:15
#' phi_1 <- 21:25
#' periods_n <- 4
#' SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
SEM_C_matrix <- function(alpha, phi_0,  periods_n, beta = c(), phi_1 = c()) {
  C1 <- matrix(rep(phi_0, periods_n))
  C1[1, 1] <- C1[1, 1] + alpha
  if (length(beta) != 0) {
    col2 <- matrix(rep(phi_1, periods_n), periods_n, byrow = TRUE)
    col2[1, 1:length(beta)] <-
      col2[1, 1:length(beta)] + beta
    C1 <- cbind(C1, col2)
  }
  C1
}

#' Covariance matrix for SEM representation
#'
#' Create covariance matrix for Simultaneous Equations Model (SEM)
#' representation. Only the part necessary to compute concentrated likelihood
#' function is computed (cf. Appendix in the Moral-Benito paper)
#'
#' @param err_var numeric
#' @param dep_vars numeric vector
#' @param phis numeric vector
#' @param psis numeric vector
#' @param psis_byrow boolean. Whether psis should be passed row-wise (i.e. they
#' will be filled into Sigma12 across rows first) or column-wise. Default is
#' TRUE i.e. row-wise.
#'
#' @return List with two matrices Sigma11 and Sigma12
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' err_var <- 1
#' dep_vars <- c(2, 2, 2, 2)
#' phis <- c(10, 10, 20, 20, 30, 30)
#' psis <- c(101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
#' SEM_sigma_matrix(err_var, dep_vars, phis, psis)
SEM_sigma_matrix <- function(err_var, dep_vars, phis = c(),
                             psis = c(), psis_byrow = TRUE) {
  periods_n <- length(dep_vars)

  O11 <- err_var^2*optimbase::ones(periods_n, periods_n) +
    diag(dep_vars^2)

  O12 <- if (length(phis) != 0) {
    regressors_n <- length(phis)/(periods_n - 1)

    phi_matrix <- matrix(rep(phis, periods_n), nrow = periods_n, byrow = TRUE)

    if (psis_byrow) {
      fill_zeros <- function(v, desired_len) {
        zeros_n <- desired_len - length(v)
        c(rep(0, zeros_n), v)
      }

      psi_matrix <- psis %>%
        split(rep(1:(periods_n-1), regressors_n*((periods_n-1):1))) %>%
        lapply(fill_zeros, desired_len = regressors_n*(periods_n-1)) %>%
        unlist() %>% matrix(nrow = periods_n - 1, byrow = TRUE) %>%
        rbind(rep(0, (periods_n - 1)*regressors_n))
    } else {
      time_fixed_psi_matrix <- function(psi, regressors_n) {
        nrows <- length(psi)/regressors_n
        t(matrix(psi, nrow = nrows, ncol = regressors_n))
      }
      psi_matrix <- psis %>%
        split(rep(1:(periods_n-1), regressors_n*(1:(periods_n-1)))) %>%
        sapply(time_fixed_psi_matrix, regressors_n = regressors_n) %>%
        plyr::rbind.fill.matrix() %>% t() %>% tidyr::replace_na(0) %>%
        rbind(rep(0, (periods_n - 1)*regressors_n))
    }

    phi_matrix + psi_matrix
  } else {
    NULL
  }
  list(O11, O12)
}

SEM_params_to_list <- function(params, periods_n, regressors_n,
                               phis_n, psis_n) {
  alpha <- params[1]
  if (regressors_n == 0) {
    beta <- c()
    phi_1 <- c()
    phis <- c()
    psis <- c()
  } else {
    beta <- params[2:(1 + regressors_n)]
    phi_1 <- params[(3 + regressors_n):(2 + 2*regressors_n)]
    phis <-
      params[(4 + 2*regressors_n + periods_n):(3 + 2*regressors_n + periods_n + phis_n)]
    psis <-
      params[(4 + 2*regressors_n + periods_n + phis_n):(3 + 2*regressors_n + periods_n + phis_n + psis_n)]
  }
  phi_0 <- params[2 + regressors_n]
  err_var <- params[3 + 2*regressors_n]
  dep_vars <-
    params[(4 + 2*regressors_n):(3 + 2*regressors_n + periods_n)]

  list(alpha = alpha, phi_0 = phi_0, err_var = err_var, dep_vars = dep_vars,
       beta = beta, phi_1 = phi_1, phis = phis, psis = psis)
}

SEM_likelihood <- function(params, n_entities,
                           cur_Y2, Y1, Z, res_maker_matrix,
                           periods_n = NULL, regressors_n = NULL,
                           phis_n = NULL, psis_n = NULL) {
  if (is.list(params)) {
    alpha <- params$alpha
    phi_0 <- params$phi_0
    err_var <- params$err_var
    dep_vars <- params$dep_vars
    beta <- if(is.null(params$beta)) c() else params$beta
    phi_1 <- if(is.null(params$phi_1)) c() else params$phi_1
    phis <- if(is.null(params$phis)) c() else params$phis
    psis <- if(is.null(params$psis)) c() else params$psis

    periods_n <- length(dep_vars)
    cur_regressors_n <- length(beta)
    B <- SEM_B_matrix(alpha, periods_n, beta)
    C <- SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
    S <- SEM_sigma_matrix(err_var, dep_vars, phis, psis)

    U1 <- if (cur_regressors_n == 0) {
      t(tcrossprod(B[[1]], Y1) - tcrossprod(C, Z))
    } else {
      t(tcrossprod(B[[1]], Y1) + tcrossprod(B[[2]], cur_Y2) - tcrossprod(C, Z))
    }
    S11_inverse <- solve(S[[1]])
    V <- if (is.null(S[[2]])) {
      cur_Y2
    } else {
      cur_Y2 - U1 %*% S11_inverse %*% S[[2]]
    }
    H <- crossprod(V, res_maker_matrix) %*% V
    likelihood <-
      -n_entities/2 * log(det(S[[1]]) * det(H/n_entities)) -
      1/2 * sum(diag(S11_inverse %*% crossprod(U1)))
  } else {
    params <- SEM_params_to_list(params, periods_n = periods_n,
                                 regressors_n = regressors_n,
                                 phis_n = phis_n, psis_n = psis_n)
    likelihood <- SEM_likelihood(params = params, n_entities = n_entities,
                                 cur_Y2 = cur_Y2, Y1 = Y1, Z = Z,
                                 res_maker_matrix = res_maker_matrix)
  }
  likelihood
}

SEM_lik_grad <- function(params, n_entities,
                         cur_Y2, Y1, Z, res_maker_matrix,
                         periods_n = NULL, regressors_n = NULL,
                         phis_n = NULL, psis_n = NULL) {
  if (is.list(params)) {
    alpha <- params$alpha
    phi_0 <- params$phi_0
    err_var <- params$err_var
    dep_vars <- params$dep_vars
    beta <- if(is.null(params$beta)) c() else params$beta
    phi_1 <- if(is.null(params$phi_1)) c() else params$phi_1
    phis <- if(is.null(params$phis)) c() else params$phis
    psis <- if(is.null(params$psis)) c() else params$psis

    periods_n <- length(dep_vars)
    cur_regressors_n <- length(beta)
    B <- SEM_B_matrix(alpha, periods_n, beta)
    C <- SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
    S <- SEM_sigma_matrix(err_var, dep_vars, phis, psis)

    U1 <- if (cur_regressors_n == 0) {
      t(tcrossprod(B[[1]], Y1) - tcrossprod(C, Z))
    } else {
      t(tcrossprod(B[[1]], Y1) + tcrossprod(B[[2]], cur_Y2) - tcrossprod(C, Z))
    }
    S11_inverse <- solve(S[[1]])
    V <- if (is.null(S[[2]])) {
      cur_Y2
    } else {
      cur_Y2 - U1 %*% S11_inverse %*% S[[2]]
    }
    H <- crossprod(V, res_maker_matrix) %*% V

    lik_vec <- optimbase::zeros(n_entities, 1)
    for (iter in 1:n_entities) {
      u10i=as.matrix(U1[iter,])
      lik_vec[iter]=-(1/2)*log(det(S[[1]]))-(1/2)*log(det(H/n_entities))-(1/2)*(t(u10i)%*%solve(S[[1]])%*%u10i)
    }
  } else {
    params <- SEM_params_to_list(params, periods_n = periods_n,
                                 regressors_n = regressors_n,
                                 phis_n = phis_n, psis_n = psis_n)
    lik_vec <- SEM_lik_grad(params = params, n_entities = n_entities,
                               cur_Y2 = cur_Y2, Y1 = Y1, Z = Z,
                               res_maker_matrix = res_maker_matrix)
  }
  lik_vec
}
