#' Calculation of the bma_prep object
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
#' @param run_parallel If \code{TRUE} the optimization is run in parallel using
#' the \link[parallel]{parApply} function. If \code{FALSE} (default value) the
#' base apply function is used. Note that using the parallel computing requires
#' setting the default cluster. See README.
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
#' \donttest{
#' library(magrittr)
#'
#' data_prepared <- economic_growth[,1:7] %>%
#'    feature_standardization(timestamp_col = year, entity_col = country) %>%
#'    feature_standardization(timestamp_col = year, entity_col = country,
#'                            time_effects = TRUE, scale = FALSE)
#'
#' for_bma <- bma_prep(df = data_prepared, dep_var_col = gdp,
#' timestamp_col = year, entity_col = country, init_value = 0.5)
#' }
#'
#' @export
#'

bma_prep <-
  function(df, timestamp_col, entity_col, dep_var_col, init_value,
           exact_value = TRUE, run_parallel = FALSE,
           control = list(trace = 2, maxit = 10000, fnscale = -1,
                          REPORT = 100, scale = 0.05)){
    model_space <- optimal_model_space(df=df, timestamp_col = {{timestamp_col}}, entity_col = {{entity_col}},
                                       dep_var_col = {{dep_var_col}}, init_value = init_value,
                                       exact_value = exact_value, run_parallel = run_parallel,
                                       control = control)

    like_table <- likelihoods_summary(df = df, dep_var_col = {{dep_var_col}}, timestamp_col = {{timestamp_col}},
                                      entity_col = {{entity_col}}, model_space = model_space,
                                      run_parallel = run_parallel)

    for_bma <- list(model_space, like_table)
  }
