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
#' number of chosen regressors. If there are no regressors returns \code{NULL}.
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

  df <- df %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ regressors }})

  if (length(colnames(df)) == 2) NULL else {
    . <- NULL
    df %>% dplyr::filter({{ timestamp_col }} > start_time) %>%
      tidyr::pivot_wider(
        names_from = {{ timestamp_col }},
        values_from = !{{ entity_col }} & !{{ timestamp_col }}
      ) %>%
      dplyr::select(!{{ entity_col }}) %>%
      dplyr::select(order(as.numeric(gsub("[^0-9]+", "", colnames(.))))) %>%
      as.matrix()
  }
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
#' natural numbers can be used as timestamps
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
      . <- NULL
      psi_matrix <- psis %>%
        split(rep(1:(periods_n-1), regressors_n*(1:(periods_n-1)))) %>%
        sapply(time_fixed_psi_matrix, regressors_n = regressors_n) %>%
        plyr::rbind.fill.matrix() %>% t() %>% replace(is.na(.), 0) %>%
        rbind(rep(0, (periods_n - 1)*regressors_n))
    }

    phi_matrix + psi_matrix
  } else {
    NULL
  }

  list(O11, O12)
}
