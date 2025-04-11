#' Approximate standard deviations for the models
#'
#' Approximate standard deviations are computed for the models in the given
#' model space. Two versions are computed.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param dep_var_col Column with the dependent variable
#' @param timestamp_col The name of the column with timestamps
#' @param entity_col Column with entities (e.g. countries)
#' @param model_space A matrix (with named rows) with each column corresponding
#' to a model. Each column specifies model parameters. Compare with
#' \link[bdsm]{optimal_model_space}
#' @param model_prior Which model prior to use. For now there are two options:
#' \code{'uniform'} and \code{'binomial-beta'}. Default is \code{'uniform'}.
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[bdsm]{SEM_likelihood} for details.
#' @param cl An optional cluster object. If supplied, the function will use this
#' cluster for parallel processing. If \code{NULL} (the default),
#' \code{pbapply::pblapply} will run sequentially.
#'
#' @return
#' Matrix with columns describing likelihood and standard deviations for each
#' model. The first row is the likelihood for the model (computed using the
#' parameters in the provided model space). The second row is almost 1/2 * BIC_k
#' as in Raftery's Bayesian Model Selection in Social Research eq. 19 (see TODO
#' in the code below). The third row is model posterior probability. Then there
#' are rows with standard deviations for each parameter. After that we have rows
#' with robust standard deviation (not sure yet what exactly "robust" means).
#'
#' @importFrom pbapply pbapply
#' @export
#'
#' @examples
#' \donttest{
#'  library(magrittr)
#'  data_prepared <- bdsm::economic_growth[, 1:7] %>%
#'    bdsm::feature_standardization(
#'      excluded_cols = c(country, year, gdp)
#'    ) %>%
#'    bdsm::feature_standardization(
#'      group_by_col  = year,
#'      excluded_cols = country,
#'      scale         = FALSE
#'    )
#'
#'  likelihoods_summary(
#'    df            = data_prepared,
#'    dep_var_col   = gdp,
#'    timestamp_col = year,
#'    entity_col    = country,
#'    model_space   = economic_growth_ms
#'  )
#' }
#'
likelihoods_summary <- function(df, dep_var_col, timestamp_col, entity_col,
                                model_space, exact_value = TRUE,
                                model_prior = 'uniform', cl = NULL) {
  regressors <- df %>%
    regressor_names(timestamp_col = {{ timestamp_col }},
                    entity_col = {{ entity_col }},
                    dep_var_col = {{ dep_var_col }})
  regressors_n <- length(regressors)
  variables_n <- regressors_n + 1

  matrices_shared_across_models <- df %>%
    matrices_from_df(timestamp_col = {{ timestamp_col }},
                     entity_col = {{ entity_col }},
                     dep_var_col = {{ dep_var_col }},
                     which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

  n_entities <- nrow(matrices_shared_across_models$Z)
  periods_n <- nrow(df) / n_entities - 1

  prior_exp_model_size <- regressors_n / 2
  prior_inc_prob <- prior_exp_model_size / regressors_n

  # THIS WILL BE DELETED
  #print(paste("Prior Mean Model Size:", prior_exp_model_size))
  #print(paste("Prior Inclusion Probability:", prior_inc_prob))

  # parameter for beta (random) distribution of the prior inclusion probability
  b <- (regressors_n - prior_exp_model_size) / prior_exp_model_size

  std_dev_from_params <- function(params, data) {
    regressors_subset <-
      regressor_names_from_params_vector(params)

    lin_features_n <- length(regressors_subset) + 1
    features_n <- ncol(data$Z)

    model_specific_matrices <- df %>%
      matrices_from_df(timestamp_col = {{ timestamp_col }},
                       entity_col = {{ entity_col }},
                       dep_var_col = {{ dep_var_col }},
                       lin_related_regressors = regressors_subset,
                       which_matrices = c("cur_Y2", "cur_Z"))

    data$cur_Z <- model_specific_matrices$cur_Z
    data$cur_Y2 <- model_specific_matrices$cur_Y2

    params_no_na <- params %>% stats::na.omit()

    likelihood <-
      SEM_likelihood(params = params_no_na, data = data,
                     exact_value = exact_value)

    hess <- hessian(SEM_likelihood, theta = params_no_na, data = data)

    likelihood_per_entity <-
      SEM_likelihood(params_no_na, data = data, per_entity = TRUE)

    # TODO: how to interpret the Gmat and Imat
    Gmat <- rootSolve::gradient(SEM_likelihood, params_no_na, data = data,
                                per_entity = TRUE)
    Imat <- crossprod(Gmat)

    # Section 2.3.3 in Moral-Benito
    # GROWTH EMPIRICS IN PANEL DATA UNDER MODEL UNCERTAINTY AND WEAK EXOGENEITY:
    # "Finally, each model-specific posterior is given by a normal distribution
    # with mean at the MLE and dispersion matrix equal to the inverse of the
    # Fisher information."
    # This is most likely why hessian is used to compute standard errors.
    # TODO: Learn the Bernstein–von Mises theorem which explain in detail how
    # all this works
    stdr <- rep(0, features_n)
    stdh <- rep(0, features_n)

    . <- NULL
    linear_params <- t(params) %>% as.data.frame() %>%
      dplyr::select(tidyselect::matches('alpha'),
                    tidyselect::matches('beta')) %>%
      as.matrix() %>% t()

    betas_first_ind <- 4 + periods_n
    betas_last_ind <- betas_first_ind + lin_features_n - 2
    inds <- if (betas_first_ind > betas_last_ind) {
      c(1)
    } else {
      c(1, betas_first_ind:betas_last_ind)
    }

    stdr[!is.na(linear_params)] <- sqrt(diag(solve(hess) %*% Imat %*% solve(hess)))[inds]
    stdh[!is.na(linear_params)] <- sqrt(diag(solve(hess)))[inds]

    # Below we have almost 1/2 * BIC_k as in Raftery's Bayesian Model Selection
    # in Social Research eq. 19. The part with reference model M_1 is skipped,
    # because we use this formula to compute exp(logl) which is in turn used to
    # compute posterior probabilities using eqs. 34/35. Since the part connected
    # with M_1 model would be present in all posteriors it cancels out. Hence
    # the important part is the one computed below.
    #
    # TODO: Why everything is divided by n_entities?

    # Eq. 19
    loglikelihood <-
      (likelihood - (lin_features_n/2)*(log(n_entities*periods_n)))/n_entities

    # Eq. 35
    bic <- exp(loglikelihood)

    c(likelihood, bic, stdh, stdr)
  }

  pbapply::pbapply(model_space, MARGIN = 2,  function(x) {
    std_dev_from_params(x, matrices_shared_across_models)
  }, cl = cl)
}
