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
#' \link[panels]{optimal_model_space}
#' @param projection_matrix_const Whether the residual maker matrix (and so
#' the projection matrix) should be computed for each model separately.
#' \code{TRUE} means that the matrix will be the same for all models
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[panels]{SEM_likelihood} for details.
#'
#' @return
#' Matrix with columns describing likelihood and standard deviations for each
#' model. First row is the likelihood for the model (computed using the
#' parameters in the provided model space). Then there are rows with standard
#' deviations for each parameter. After that we have rows with robust standard
#' deviation (not sure yet what exactly "robust" means).
#'
#' @export
#'
#' @examples
#' data_centered_scaled <-
#'   feature_standardization(df = panels::economic_growth[,1:7],
#'                           timestamp_col = year, entity_col = country)
#' data_cross_sectional_standarized <-
#'   feature_standardization(df = data_centered_scaled, timestamp_col = year,
#'                           entity_col = country, cross_sectional = TRUE,
#'                           scale = FALSE)
#'
#' bma_result <- likelihoods_summary(df = data_cross_sectional_standarized,
#'                                   dep_var_col = gdp, timestamp_col = year,
#'                                   entity_col = country,
#'                                   model_space = economic_growth_ms,
#'                                   projection_matrix_const = TRUE)
#'
likelihoods_summary <- function(df, dep_var_col, timestamp_col, entity_col,
                                model_space, projection_matrix_const,
                                exact_value = TRUE) {
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
                     exact_value = exact_value,
                     projection_matrix_const = projection_matrix_const)

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

    stdr[!is.na(linear_params)] <- sqrt(diag(solve(hess) %*% Imat %*% solve(hess)))[1:lin_features_n]
    stdh[!is.na(linear_params)] <- sqrt(diag(solve(hess)))[1:lin_features_n]

    c(likelihood, stdh, stdr)
  }

  apply(model_space, 2, std_dev_from_params, matrices_shared_across_models)
}

#' Summary of a model space
#'
#' A summary of a given model space is prepared. This include things such as
#' posterior inclusion probability (PIP), posterior mean and so on. This is the
#' core function of the package, because it allows to make assessments and
#' decisions about the parameters and models.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param dep_var_col Column with the dependent variable
#' @param timestamp_col The name of the column with timestamps
#' @param entity_col Column with entities (e.g. countries)
#' @param model_space A matrix (with named rows) with each column corresponding
#' to a model. Each column specifies model parameters. Compare with
#' \link[panels]{optimal_model_space}
#' @param model_prior Which model prior to use. For now there are two options:
#' \code{'uniform'} and \code{'binomial-beta'}. Default is \code{'uniform'}.
#' @param projection_matrix_const Whether the residual maker matrix (and so
#' the projection matrix) should be computed for each model separately.
#' \code{TRUE} means that the matrix will be the same for all models
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[panels]{SEM_likelihood} for details.
#'
#' @return
#' List of parameters describing analysed models
#'
#' @export
bma_summary <- function(df, dep_var_col, timestamp_col, entity_col,
                        model_space, projection_matrix_const,
                        exact_value = TRUE, model_prior = 'uniform') {
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

  # parameter for beta (random) distribution of the prior inclusion probability
  b <- (regressors_n - prior_exp_model_size) / prior_exp_model_size

  mod <- optimbase::zeros(variables_n,1)
  bet <- optimbase::zeros(variables_n,1)
  pvarh <- optimbase::zeros(variables_n,1)
  pvarr <- optimbase::zeros(variables_n,1)
  fy <- optimbase::zeros(variables_n,1)
  fyt <- 0
  ppmsize <- 0
  cout <- 0

  likelihoods_info <- likelihoods_summary(
    df, dep_var_col = {{ dep_var_col }}, timestamp_col = {{ timestamp_col }},
    entity_col = {{ entity_col }}, model_space = model_space,
    projection_matrix_const = projection_matrix_const, exact_value = exact_value
    )

  regressors_subsets <- rje::powerSet(regressors)
  regressors_subsets_matrix <-
    rje::powerSetMat(regressors_n) %>% as.data.frame()

  row_ind <- 0
  for (regressors_subset in regressors_subsets) {
    row_ind <- row_ind + 1
    print(paste('Progress:', row_ind, 'out of', length(regressors_subsets)))
    mt <- as.matrix(t(regressors_subsets_matrix[row_ind, ]))
    out = (mt == 0)       # regressors out of the current model
    cur_regressors_n <- sum(mt)
    cur_variables_n <- cur_regressors_n+1

    likelihood_max <- likelihoods_info[1, row_ind]

    # Below we have almost 1/2 * BIC_k as in Raftery's Bayesian Model Selection
    # in Social Research eq. 19. The part with reference model M_1 is skipped,
    # because we use this formula to compute exp(logl) which is in turn used to
    # compute posterior probabilities using eqs. 34/35. Since the part connected
    # with M_1 model would be present in all posteriors it cancels out. Hence
    # the important part is the one computed below.
    #
    # TODO: Why everything is divided by n_entities?

    # Eq. 19
    logl <- (likelihood_max-(cur_variables_n/2)*(log(n_entities*periods_n)))/n_entities

    # Eq. 35
    bict <- exp(logl)

    if (model_prior == 'binomial-beta') {
      prior_model_prob <-
        gamma(1 + cur_regressors_n) * gamma(b + regressors_n - cur_regressors_n)
    } else if (model_prior == 'uniform') {
      prior_model_prob <-
        prior_inc_prob^cur_regressors_n *
        (1-prior_inc_prob)^(regressors_n - cur_regressors_n)
    } else {
      stop("Please specify a correct model prior!")
    }

    # posterior model probability  #
    postprob=prior_model_prob*bict

    # constructing the full vector of estimates #
    mty=rbind(1,mt)

    print(likelihoods_info[, row_ind])
    stdrt1 <- likelihoods_info[-(1:(regressors_n+2)), row_ind]
    stdht1 <- likelihoods_info[2:(regressors_n+2), row_ind]
    print(stdrt1)
    varrt1 <- stdrt1^2
    varht1 <- stdht1^2


    . <- NULL
    linear_params <- t(model_space[, row_ind]) %>% as.data.frame() %>%
      dplyr::select(tidyselect::matches('alpha'),
                    tidyselect::matches('beta')) %>% replace(is.na(.), 0) %>%
      as.matrix() %>% t()


    # calculating the percentage of significant regressions #
    ptr=linear_params/stdht1
    ntr=linear_params/stdht1
    if (row_ind==1) {
      pts=ptr; nts=ntr
    }
    else {
      pts=cbind(pts,ptr); nts=cbind(nts,ntr)
    }

    # accumulating estimates for posterior model probabilities #
    mod=mod+mty
    fy=fy+postprob*mty
    fyt=fyt+postprob
    ppmsize=ppmsize+postprob*(sum(mty))

    # storing estimates conditional on inclusion #
    bet=bet+postprob*linear_params
    pvarr=pvarr+(postprob*varrt1+postprob*(linear_params*linear_params))         # as in Leamer (1978) #
    pvarh=pvarh+(postprob*varht1+postprob*(linear_params*linear_params))         # as in Leamer (1978) #

    # here we store model-specific diagnostics and estimates (BICs, likelihoods...) #
    if (row_ind==1) {
      models_posterior_prob <- postprob
      models_prior_prob <- prior_model_prob
      liks <- exp(likelihood_max/n_entities)
      bics <- bict
      foutt <- likelihood_max
    }
    else {
      models_posterior_prob <- rbind(models_posterior_prob, postprob)
      models_prior_prob <- rbind(models_prior_prob, prior_model_prob)
      liks <- rbind(liks, exp(likelihood_max/n_entities))
      bics <- rbind(bics ,bict)
      foutt <- rbind(foutt, likelihood_max)
    }
  }

  list(prior_exp_model_size = prior_exp_model_size,
       prior_inc_prob = prior_inc_prob, variables_n = variables_n,
       models_posterior_prob = models_posterior_prob,
       models_prior_prob = models_prior_prob, liks = liks,
       bics = bics, foutt = foutt,
       bet = bet, mod = mod, pvarh = pvarh, pvarr = pvarr, fy = fy, fyt = fyt,
       ppmsize = ppmsize, cout = 0, nts = nts, pts = pts)
}
