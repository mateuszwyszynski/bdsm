regressor_names <- function(df, timestamp_col, entity_col, dep_var_col) {
  df %>%
    dplyr::select(
      ! c({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }})
    ) %>% colnames()
}

determine_min_timestamps <- function(df, timestamp_col) {
  timestamps <- dplyr::select(df, {{ timestamp_col }})
  timestamp_0 <- min(timestamps)
  timestamp_1 <- min(timestamps[timestamps != timestamp_0])
  list(timestamp_0 = timestamp_0, timestamp_1 = timestamp_1)
}

#' Matrix with dependent variable data for SEM representation
#'
#' Create matrix which contains dependent variable data used in the Simultaneous
#' Equations Model (SEM) representation on the left hand side of the equations.
#' The matrix contains the data for time periods greater than or equal to the
#' second lowest time stamp. The matrix is then used to compute likelihood for
#' SEM analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#'
#' @return
#' Matrix of size N x T where N is the number of entities considered and T is
#' the number of periods greater than or equal to the second lowest time stamp.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   entities = rep(1:4, 5),
#'   times = rep(seq(1960, 2000, 10), each = 4),
#'   dep_var = stats::rnorm(20), a = stats::rnorm(20), b = stats::rnorm(20)
#' )
#' sem_dep_var_matrix(df, times, entities, dep_var)
sem_dep_var_matrix <- function(df, timestamp_col, entity_col, dep_var_col) {
  min_timestamps <-
    determine_min_timestamps(df = df, timestamp_col = {{ timestamp_col }})
  timestamp_1 <- min_timestamps$timestamp_1

  df %>% dplyr::filter({{ timestamp_col }} >= timestamp_1) %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }}) %>%
    tidyr::pivot_wider(names_from = {{ timestamp_col }},
                       values_from = {{ dep_var_col }}) %>%
    dplyr::select(!{{ entity_col }}) %>% as.matrix()
}

#' Matrix with regressors data for SEM representation
#'
#' Create matrix which contains regressors data used in the Simultaneous
#' Equations Model (SEM) representation on the left hand side of the equations.
#' The matrix contains regressors data for time periods greater than or equal to
#' the second lowest time stamp. The matrix is then used to compute likelihood
#' for SEM analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#'
#' @return
#' Matrix of size N x (T-1)*k where N is the number of entities considered, T is
#' the number of periods greater than or equal to the second lowest time stamp
#' and k is the number of chosen regressors. If there are no regressors returns
#' \code{NULL}.
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   entities = rep(1:4, 5),
#'   times = rep(seq(1960, 2000, 10), each = 4),
#'   dep_var = stats::rnorm(20), a = stats::rnorm(20), b = stats::rnorm(20)
#' )
#' sem_regressors_matrix(df, times, entities, dep_var)
sem_regressors_matrix <- function(df, timestamp_col, entity_col, dep_var_col) {
  regressors <- df %>%
    regressor_names(timestamp_col = {{ timestamp_col }},
                    entity_col = {{ entity_col }},
                    dep_var_col = {{ dep_var_col }})

  min_timestamps <-
    determine_min_timestamps(df = df, timestamp_col = {{ timestamp_col }})
  timestamp_1 <- min_timestamps$timestamp_1

  df <- df %>%
    dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ regressors }})

  if (length(colnames(df)) == 2) NULL else {
    . <- NULL
    df %>% dplyr::filter({{ timestamp_col }} > timestamp_1) %>%
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
#' from the lowest time stamp and regressors from the second lowest time stamp.
#' The matrix is then used to compute likelihood for SEM analysis.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#'
#' @return
#' Matrix of size N x k+1 where N is the number of entities considered and k is
#' the number of chosen regressors
#' @export
#'
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   entities = rep(1:4, 5),
#'   times = rep(seq(1960, 2000, 10), each = 4),
#'   dep_var = stats::rnorm(20), a = stats::rnorm(20), b = stats::rnorm(20)
#' )
#' exogenous_matrix(df, times, entities, dep_var)
exogenous_matrix <- function(df, timestamp_col, entity_col, dep_var_col) {
  regressors <- df %>%
    regressor_names(timestamp_col = {{ timestamp_col }},
                    entity_col = {{ entity_col }},
                    dep_var_col = {{ dep_var_col }})

  min_timestamps <-
    determine_min_timestamps(df = df, timestamp_col = {{ timestamp_col }})
  timestamp_1 <- min_timestamps$timestamp_1
  timestep <- timestamp_1 - min_timestamps$timestamp_0

  df_with_lagged_col <- df %>%
    dplyr::select({{ entity_col }}, {{ timestamp_col }}, {{ dep_var_col }}) %>%
    dplyr::filter({{ timestamp_col }} == (timestamp_1 - timestep)) %>%
    dplyr::mutate("{{timestamp_col}}" := {{ timestamp_col }} + timestep)

  df %>%
    dplyr::filter({{ timestamp_col }} == timestamp_1) %>%
    dplyr::select(!{{ dep_var_col }}) %>%
    dplyr::left_join(df_with_lagged_col,
              by = dplyr::join_by(
                {{ timestamp_col }} == {{ timestamp_col }},
                {{ entity_col }} == {{ entity_col }}
              )) %>%
    dplyr::select({{ dep_var_col }}, {{ regressors }}) %>% as.matrix()
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
#' sem_sigma_matrix(err_var, dep_vars, phis, psis)
sem_sigma_matrix <- function(err_var, dep_vars, phis = c(), psis = c()) {
  periods_n <- length(dep_vars)

  O11 <- err_var^2*optimbase::ones(periods_n, periods_n) +
    diag(dep_vars^2)

  O12 <- if (length(phis) != 0) {
    regressors_n <- length(phis)/(periods_n - 1)

    phi_matrix <- matrix(rep(phis, periods_n), nrow = periods_n, byrow = TRUE)
    psi_matrix <- sem_psi_matrix(psis = psis, timestamps_n = periods_n,
                                 features_n = regressors_n)

    phi_matrix + psi_matrix
  } else {
    NULL
  }

  list(O11, O12)
}
