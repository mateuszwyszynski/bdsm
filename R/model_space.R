#' Initialize model space matrix
#'
#' This function builds a representation of the model space, by creating a
#' dataframe where each column represents values of the parameters for a given
#' model. Real value means that the parameter is included in the model. A
#' parameter not present in the model is marked as \code{NA}.
#'
#' Currently the set of features is assumed to be all columns which remain after
#' excluding \code{timestamp_col}, \code{entity_col} and \code{dep_var_col}.
#'
#' A power set of all possible exclusions of linear dependence on the given
#' feature is created, i.e. if there are 4 features we end up with 2^4 possible
#' models (for each model we independently decide whether to include or not a
#' feature).
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col Column which determines time periods. For now only
#' natural numbers can be used as timestamps
#' @param entity_col Column which determines entities (e.g. countries, people)
#' @param dep_var_col Column with dependent variable
#' @param init_value Initial value for parameters present in the model. Default
#' is \code{1}.
#'
#' @return matrix of model parameters
#'
#' @examples
#' library(magrittr)
#'
#' data_prepared <- bdsm::economic_growth[, 1:5] %>%
#'   bdsm::feature_standardization(
#'     excluded_cols = c(country, year, gdp)
#'   ) %>%
#'   bdsm::feature_standardization(
#'     group_by_col  = year,
#'     excluded_cols = country,
#'     scale         = FALSE
#'   )
#'
#' init_model_space_params(data_prepared, year, country, gdp)
#' @export
init_model_space_params <- function(df, timestamp_col, entity_col,
                                    dep_var_col, init_value = 1) {
  regressors <- df %>%
    regressor_names(timestamp_col = {{ timestamp_col }},
                    entity_col = {{ entity_col }},
                    dep_var_col = {{ dep_var_col }})
  regressors_n <- length(regressors)

  counts <- df %>% dplyr::select({{ timestamp_col }}, {{ entity_col }}) %>%
    sapply(function(x) dplyr::n_distinct(x))

  timestamps_n <- counts[[1]] - 1
  entities_n <- counts[[2]]

  regressors_subsets_matrix <-
    (rje::powerSetMat(regressors_n) * init_value)

  linear_params_matrix <-
    t(cbind(regressors_subsets_matrix, regressors_subsets_matrix))

  rownames(linear_params_matrix) <-
    c(paste('beta', regressors, sep="_"), paste('phi_1', regressors, sep="_"))

  dep_var_matrix <-
    t(matrix(init_value, nrow = 2^regressors_n, ncol = 3 + timestamps_n))
  rownames(dep_var_matrix) <- c(c('alpha', 'phi_0', 'err_var'),
                                paste('dep_var', 1:timestamps_n, sep = '_'))

  phis_n <- regressors_n * (timestamps_n - 1)
  psis_n <- regressors_n * (timestamps_n - 1) * timestamps_n / 2

  psis_phis_matrix <-
    matrix(init_value, nrow = phis_n + psis_n, ncol = 2^regressors_n)

  . <- NULL
  rbind(dep_var_matrix, linear_params_matrix, psis_phis_matrix) %>%
    replace(. == 0, NA)
}


#' Helper function to extract names from a vector defining a model
#'
#' For now it is assumed that we can only exclude linear relationships between
#' regressors and the dependent variable.
#'
#' The vector needs to have named rows, i.e. it is assumed it comes from a
#' model space (see \link[bdsm]{init_model_space_params} for details).
#'
#' @param params a vector with parameters describing the model
#'
#' @return
#' Names of regressors which are assumed to be linearly connected with dependent
#' variable within the model described by the \code{params} vector.
#'
#' @examples
#' params <- c(alpha = 1, beta_gdp = 1, beta_gdp_lagged = 1, phi_0 = 1, err_var = 1)
#' regressor_names_from_params_vector(params)
#'
#' @export
regressor_names_from_params_vector <- function(params) {
  regressors_subset <-
    t(params %>% stats::na.omit()) %>% as.data.frame() %>%
    dplyr::select(tidyselect::matches('beta'))

  names(regressors_subset) <- gsub("beta_", "", names(regressors_subset))
  names(regressors_subset)
}


#' Finds MLE parameters for each model in the given model space
#'
#' Given a dataset and an initial value for parameters, initializes a model
#' space with parameters equal to the initial value for each model. Then for each
#' model performs a numerical optimization and finds parameters which maximize
#' the likelihood.
#'
#' @param df Data frame with data for the analysis.
#' @param timestamp_col The name of the column with time stamps.
#' @param entity_col Column with entities (e.g. countries).
#' @param dep_var_col Column with the dependent variable.
#' @param init_value The value with which the model space will be initialized.
#' This will be the starting point for the numerical optimization.
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[bdsm]{SEM_likelihood} for details.
#' @param cl An optional cluster object. If supplied, the function will use this
#' cluster for parallel processing. If \code{NULL} (the default),
#' \code{pbapply::pblapply} will run sequentially.
#' @param control a list of control parameters for the optimization which are
#' passed to \link[stats]{optim}. Default is
#' \code{list(trace = 2, maxit = 10000, fnscale = -1, REPORT = 100, scale = 0.05)}.
#'
#' @return
#' List (or matrix) of parameters describing analyzed models.
#'
#' @importFrom pbapply pbapply
#'
#' @export
optim_model_space_params <- function(df, timestamp_col, entity_col, dep_var_col, init_value,
                                     exact_value = TRUE, cl = NULL,
                                     control = list(trace = 2, maxit = 10000,
                                                    fnscale = -1, REPORT = 100,
                                                    scale = 0.05)) {
  matrices_shared_across_models <- df %>%
    matrices_from_df(timestamp_col = {{ timestamp_col }},
                     entity_col = {{ entity_col }},
                     dep_var_col = {{ dep_var_col }},
                     which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

  model_space <- df %>%
    init_model_space_params(timestamp_col = {{ timestamp_col }},
                            entity_col = {{ entity_col }},
                            dep_var_col = {{ dep_var_col }},
                            init_value = init_value)

  optimization_wrapper <- function(params, data) {
    regressors_subset <- regressor_names_from_params_vector(params)

    model_specific_matrices <- df %>%
      matrices_from_df(timestamp_col = {{ timestamp_col }},
                       entity_col = {{ entity_col }},
                       dep_var_col = {{ dep_var_col }},
                       lin_related_regressors = regressors_subset,
                       which_matrices = c("cur_Y2", "cur_Z"))

    data$cur_Z <- model_specific_matrices$cur_Z
    data$cur_Y2 <- model_specific_matrices$cur_Y2

    params_no_na <- stats::na.omit(params)

    # Adjust the optimization control parameters.
    control$parscale <- control$scale * params_no_na
    control$scale <- NULL

    optimized <- stats::optim(params_no_na, SEM_likelihood, data = data,
                              exact_value = exact_value,
                              method = "BFGS",
                              control = control)

    params[!is.na(params)] <- optimized[[1]]
    params
  }

  pbapply::pbapply(model_space, MARGIN = 2,  function(x) {
    optimization_wrapper(x, matrices_shared_across_models)
  }, cl = cl)
}


#' Calculation of the model_space object
#'
#' This function calculates model space, values of the maximized likelihood function, BICs, and
#' standard deviations of the parameters that will be used in Bayesian model averaging.
#'
#' @param df Data frame with data for the analysis.
#' @param timestamp_col The name of the column with time stamps
#' @param entity_col Column with entities (e.g. countries)
#' @param dep_var_col Column with the dependent variable
#' @param init_value The value with which the model space will be initialized.
#' This will be the starting point for the numerical optimization.
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[bdsm]{SEM_likelihood} for details.
#' @param cl An optional cluster object. If supplied, the function will use this
#' cluster for parallel processing. If \code{NULL} (the default),
#' \code{pbapply::pblapply} will run sequentially.
#' @param control a list of control parameters for the optimization which are
#' passed to \link[stats]{optim}. Default is
#' \code{list(trace = 2, maxit = 10000, fnscale = -1, REPORT = 100, scale = 0.05)}, but note
#' that \code{scale} is used only for adjusting the \code{parscale} element added later in the function code.
#'
#' @importFrom parallel parApply
#'
#' @return
#' List with two objects: \cr
#' 1) model_space - table with parameters of all estimated models \cr
#' 2) like_table - table with the value of maximized likelihood function, BIC, and
#' standard errors for all estimated models
#'
#' @examples
#' library(magrittr)
#'
#' data_prepared <- bdsm::economic_growth[, 1:5] %>%
#'   bdsm::feature_standardization(
#'     excluded_cols = c(country, year, gdp)
#'   ) %>%
#'   bdsm::feature_standardization(
#'     group_by_col  = year,
#'     excluded_cols = country,
#'     scale         = FALSE
#'   )
#'
#' for_bma <- model_space(
#'   df            = data_prepared,
#'   dep_var_col   = gdp,
#'   timestamp_col = year,
#'   entity_col    = country,
#'   init_value    = 0.5
#'  )
#'
#' @export
#
model_space <-
  function(df, timestamp_col, entity_col, dep_var_col, init_value,
           exact_value = TRUE, cl = NULL,
           control = list(trace = 2, maxit = 10000, fnscale = -1,
                          REPORT = 100, scale = 0.05)){
    model_space <- optim_model_space_params(
      df            = df,
      timestamp_col = {{timestamp_col}},
      entity_col    = {{entity_col}},
      dep_var_col   = {{dep_var_col}},
      init_value    = init_value,
      exact_value   = exact_value,
      cl            = cl,
      control       = control
    )

    like_table <- model_space_stats(
      df            = df,
      dep_var_col   = {{dep_var_col}},
      timestamp_col = {{timestamp_col}},
      entity_col    = {{entity_col}},
      model_space   = model_space,
      cl            = cl
    )

    for_bma <- list(model_space, like_table)
  }
