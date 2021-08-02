#' Matrix with dependent variable data for SEM representation
#'
#' Create matrix which contains dependent variable data used in the Simultaneous
#' Equations Model (SEM) representation on the left hand side of the equations.
#' The matrix contains the data for time periods greater than chosen
#' \code{start_time}. The matrix is then used to compute likelihood for SEM
#' analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#' @param start_time First time period. Only time periods greater than or equal
#' to \code{start_time} will be considered in the resulting matrix. For
#' \code{start_time=NULL} second lowest value from \code{timestamp_col} is set
#' as \code{start_time}. Default is \code{NULL}.
#'
#' @return
#' Matrix of size N x T where N is the number of entities considered and T is
#' the number of periods greater than or equal to \code{start_time}.
#'
#' @export
#'
#' @examples
SEM_dep_var_matrix <- function(df, timestamp_col, entity_col, dep_var_col,
                               start_time = NULL) {
  if (is.null(start_time)) {
    timestamps <- dplyr::select(df, {{ timestamp_col }})
    time_zero <- min(timestamps)
    start_time <- min(timestamps[timestamps != time_zero])
  }
  df %>% dplyr::filter({{ timestamp_col }} >= start_time) %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }}) %>%
    tidyr::pivot_wider(names_from = {{ timestamp_col }},
                       values_from = {{ dep_var_col }}) %>%
    dplyr::select(!{{ entity_col }}) %>% as.matrix()
}

#' Matrix with regressors data for SEM representation
#'
#' Create matrix which contains regressors data used in the Simultaneous
#' Equations Model (SEM) representation on the left hand side of the equations.
#' The matrix contains regressors data for time periods greater than chosen
#' \code{start_time}. The matrix is then used to compute likelihood for SEM
#' analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param regressors Which subset of columns should be used as regressors.
#' @param start_time First time period. Only time periods greater than
#' \code{start_time} will be considered in the resulting matrix. For
#' \code{start_time=NULL} second lowest value from \code{timestamp_col} is set
#' as \code{start_time}. Default is \code{NULL}.
#'
#' @return
#' Matrix of size N x (T-1)*k where N is the number of entities considered, T is
#' the number of periods greater than or equal to \code{start_time} and k is the
#' number of chosen regressors
#' @export
#'
#' @examples
SEM_regressors_matrix <- function(df, timestamp_col, entity_col, regressors,
                                  start_time = NULL) {
  if (is.null(start_time)) {
    timestamps <- dplyr::select(df, {{ timestamp_col }})
    time_zero <- min(timestamps)
    start_time <- min(timestamps[timestamps != time_zero])
  }
  . <- NULL
  df %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ regressors }}) %>%
    dplyr::filter({{ timestamp_col }} > start_time) %>%
    tidyr::pivot_wider(
      names_from = {{ timestamp_col }},
      values_from = !{{ entity_col }} & !{{ timestamp_col }}
      ) %>%
    dplyr::select(!{{ entity_col }}) %>%
    dplyr::select(order(as.numeric(gsub("[^0-9]+", "", colnames(.))))) %>%
    as.matrix()
}

#' Matrix with exogenous variables for SEM representation
#'
#' Create matrix which contains exogenous variables used in the Simultaneous
#' Equations Model (SEM) representation. Currently these are: dependent variable
#' from the \code{start_time} period and regressors from the next period after
#' the \code{strat_time}. The matrix is then used to compute likelihood for SEM
#' analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestampsg
#' @param start_time First time period. Only time periods greater than
#' \code{start_time} will be considered in the resulting matrix
#' @param lagged_col Column which contains lagged version of dependent variable
#' @param regressors_subset Which subset of columns should be used as
#' regressors. If \code{NULL} (default) then all remaining columns will be used
#' as regressors. For now columns have to be passed as list of column names
#' represented as strings.
#'
#' @return
#' Matrix of size N x k+1 where N is the number of entities considered and k is
#' the number of chosen regressors
#' @export
#'
#' @examples
SEM_exogenous_matrix <- function(df, timestamp_col, start_time, lagged_col,
                                 regressors_subset = NULL) {
  df %>% dplyr::filter({{ timestamp_col }} == start_time) %>%
    dplyr::select({{ lagged_col}}, regressors_subset) %>% as.matrix()
}

#' Residual Maker Matrix
#'
#' Create residual maker matrix from a given matrix \code{m}. See article about
#' \href{https://en.wikipedia.org/wiki/Projection_matrix}{projection matrix} on
#' the Wikipedia.
#'
#' @param m Matrix
#'
#' @return
#' M x M matrix where M is the number of rows in the \code{m} matrix.
#' @export
#'
#' @examples
residual_maker_matrix <- function(m) {
  proj_matrix <- m%*%solve(crossprod(m))%*%t(m)
  res_maker_matrix <- diag(nrow(m)) - proj_matrix
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

SEM_params_to_list <- function(params, periods_n, tot_regressors_n,
                               in_regressors_n, phis_n, psis_n) {
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

SEM_likelihood <- function(params, data, timestamp_col = NULL,
                           entity_col = NULL, start_time = NULL,
                           lagged_col = NULL, dep_var_col = NULL,
                           regressors = NULL, in_regressors = NULL,
                           periods_n = NULL, tot_regressors_n = NULL,
                           in_regressors_n = NULL,
                           phis_n = NULL, psis_n = NULL, grad = FALSE) {
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
    res_maker_matrix <- data$res_maker_matrix

    n_entities <- nrow(Z)
    periods_n <- length(dep_vars)
    in_regressors_n <- length(beta)
    B <- SEM_B_matrix(alpha, periods_n, beta)
    C <- SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
    S <- SEM_sigma_matrix(err_var, dep_vars, phis, psis)

    U1 <- if (in_regressors_n == 0) {
      t(tcrossprod(B[[1]], Y1) - tcrossprod(C, cur_Z))
    } else {
      t(tcrossprod(B[[1]], Y1) + tcrossprod(B[[2]], cur_Y2) - tcrossprod(C, cur_Z))
    }
    S11_inverse <- solve(S[[1]])
    V <- Y2 - U1 %*% S11_inverse %*% S[[2]]
    H <- crossprod(V, res_maker_matrix) %*% V
    likelihood <- if(!grad) {
      -n_entities/2 * log(det(S[[1]]) * det(H/n_entities)) -
        1/2 * sum(diag(S11_inverse %*% crossprod(U1)))
    } else {
      lik_vec <- optimbase::zeros(n_entities, 1)
      for (iter in 1:n_entities) {
        u10i=as.matrix(U1[iter,])
        lik_vec[iter]=-(1/2)*log(det(S[[1]]))-(1/2)*log(det(H/n_entities))-(1/2)*(t(u10i)%*%solve(S[[1]])%*%u10i)
      }
      lik_vec
    }
  } else {
    if (!is.list(params)) {
      params <- SEM_params_to_list(params, periods_n = periods_n,
                                   tot_regressors_n = tot_regressors_n,
                                   in_regressors_n = in_regressors_n,
                                   phis_n = phis_n, psis_n = psis_n)
    }
    if (!is.list(data)) {
      Y1 <- SEM_dep_var_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        dep_var_col = dep_var_col, start_time = start_time
        )
      Y2 <- SEM_regressors_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        regressors = regressors, start_time = start_time
      )
      cur_Y2 <- SEM_regressors_matrix(
        df = data, timestamp_col = timestamp_col, entity_col = entity_col,
        regressors = in_regressors, start_time = start_time
      )
      cur_Z <- SEM_exogenous_matrix(
        df = data, timestamp_col = timestamp_col, start_time = start_time,
        lagged_col = lagged_col, regressors_subset = in_regressors
      )
      Z <- SEM_exogenous_matrix(
        df = data, timestamp_col = timestamp_col, start_time = start_time,
        lagged_col = lagged_col
      )
      res_maker_matrix <- residual_maker_matrix(Z)

      data = list(Y1 = Y1, Y2 = Y2, cur_Y2 = cur_Y2, Z = cur_Z,
                  res_maker_matrix = res_maker_matrix)
    }
    likelihood <- SEM_likelihood(params = params, data = data, grad = grad)
  }
  likelihood
}
