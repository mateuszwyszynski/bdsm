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
#' @return
#' @export
#'
#' @examples
initialize_model_space <- function(df, timestamp_col, entity_col,
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
#' model space (see \link[panels]{initialize_model_space} for details).
#'
#' @param params a vector with parameters describing the model
#'
#' @return
#' Names of regressors which are assumed to be linearly connected with dependent
#' variable within the model described by the \code{params} vector.
#' @export
#'
#' @examples
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
#' space with parameters equal to initial value for each model. Then for each
#' model performs a numerical optimization and finds parameters which maximize
#' the likelihood.
#'
#' @param df Data frame with data for the SEM analysis.
#' @param timestamp_col The name of the column with time stamps
#' @param entity_col Column with entities (e.g. countries)
#' @param dep_var_col Column with the dependent variable
#' @param init_value The value with which the model space will be initialized.
#' This will be the starting point for the numerical optimization.
#' @param projection_matrix_const Whether the residual maker matrix (and so
#' the projection matrix) should be computed for each model separately.
#' \code{TRUE} means that the matrix will be the same for all models
#' @param exact_value Whether the exact value of the likelihood should be
#' computed (\code{TRUE}) or just the proportional part (\code{FALSE}). Check
#' \link[panels]{SEM_likelihood} for details.
#' @param control a list of control parameters for the optimization which are
#' passed to \link[stats]{optim}. Default is
#' \code{list(trace = 2, maxit = 10000, fnscale = -1, REPORT = 100)}, but note
#' that a \code{parscale} element is also added later in the function code.
#' For now it is hard coded with no control on the user side.
#'
#' @return
#' List of parameters describing analysed models
#'
#' @export
optimal_model_space <-
  function(df, timestamp_col, entity_col, dep_var_col, init_value,
           projection_matrix_const, exact_value = TRUE,
           control = list(trace = 2, maxit = 10000, fnscale = -1,
                          REPORT = 100)) {
  Y1 <- SEM_dep_var_matrix(
    df = df, timestamp_col = {{ timestamp_col }},
    entity_col = {{ entity_col }}, dep_var_col = {{ dep_var_col }}
  )

  Y2 <- df %>%
    SEM_regressors_matrix(timestamp_col = {{ timestamp_col }},
                          entity_col = {{ entity_col }},
                          dep_var_col = {{ dep_var_col }})

  Z <- df %>%
    exogenous_matrix(timestamp_col = {{ timestamp_col }},
                     entity_col = {{ entity_col }},
                     dep_var_col = {{ dep_var_col }})

  res_maker_matrix <- residual_maker_matrix(Z)

  model_space <- df %>%
    initialize_model_space(timestamp_col = {{ timestamp_col }},
                           entity_col = {{ entity_col}},
                           dep_var_col = {{ dep_var_col }},
                           init_value = init_value)

  constant_data <- list(Y1 = Y1, Y2 = Y2, res_maker_matrix = res_maker_matrix,
                        Z = NULL, cur_Y2 = NULL)

  optimization_wrapper <- function(params, data) {
    regressors_subset <-
      regressor_names_from_params_vector(params)

    data$Z <- df %>%
      dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }},
                    regressors_subset) %>%
      exogenous_matrix({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }})

    data$cur_Y2 <- df %>%
      dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }},
                    regressors_subset) %>%
      SEM_regressors_matrix(timestamp_col = {{ timestamp_col }},
                            entity_col = {{ entity_col }},
                            dep_var_col = {{ dep_var_col }})

    params_no_na <- params %>% stats::na.omit()

    # parscale argument somehow (don't know yet how) changes step size during
    # optimisation. Most likely optimisation methods used in Gauss are
    # scale-free and these used in R are not
    # TODO: search for methods (or implement methods) in R which are scale-free
    control$parscale = 0.05*params_no_na

    optimized <- stats::optim(params_no_na, SEM_likelihood, data = data,
                              exact_value = exact_value,
                              projection_matrix_const = projection_matrix_const,
                              method="BFGS",
                              control = control)

    params[!is.na(params)] <- optimized[[1]]
    params
  }

  model_space <- apply(model_space, 2, optimization_wrapper, constant_data)
  model_space
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

  Y1 <- SEM_dep_var_matrix(
    df = df, timestamp_col = {{ timestamp_col }},
    entity_col = {{ entity_col }}, dep_var_col = {{ dep_var_col }}
  )

  Y2 <- df %>%
    SEM_regressors_matrix(timestamp_col = {{ timestamp_col }},
                          entity_col = {{ entity_col }},
                          dep_var_col = {{ dep_var_col }})

  Z <- df %>%
    exogenous_matrix(timestamp_col = {{ timestamp_col }},
                     entity_col = {{ entity_col }},
                     dep_var_col = {{ dep_var_col }})

  n_entities <- nrow(Z)
  periods_n <- nrow(df) / n_entities - 1

  res_maker_matrix <- residual_maker_matrix(Z)

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

    cur_Z <- df %>%
      dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }},
                    regressors_subset) %>%
      exogenous_matrix({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }})

    cur_Y2 <- df %>%
      dplyr::select({{ timestamp_col }}, {{ entity_col }}, {{ dep_var_col }},
                    regressors_subset) %>%
      SEM_regressors_matrix(timestamp_col = {{ timestamp_col }},
                            entity_col = {{ entity_col }},
                            dep_var_col = {{ dep_var_col }})

    data <- list(Y1 = Y1, Y2 = Y2, cur_Y2 = cur_Y2, Z = cur_Z,
                 res_maker_matrix = res_maker_matrix)

    optimised_params <- model_space[, row_ind] %>% stats::na.omit()
    likelihood_max <-
      SEM_likelihood(params = optimised_params, data = data,
                     exact_value = exact_value,
                     projection_matrix_const = projection_matrix_const)


    hess <- hessian(SEM_likelihood, theta = optimised_params, data = data)

    likelihood_per_entity <-
      SEM_likelihood(optimised_params, data = data, per_entity = TRUE)

    Gmat <- rootSolve::gradient(SEM_likelihood, optimised_params, data = data,
                                per_entity = TRUE)
    Imat=crossprod(Gmat)
    stdr=sqrt(diag(solve(hess)%*%(Imat)%*%solve(hess)))

    # Section 2.3.3 in Moral-Benito
    # GROWTH EMPIRICS IN PANEL DATA UNDER MODEL UNCERTAINTY AND WEAK EXOGENEITY:
    # "Finally, each model-specific posterior is given by a normal distribution
    # with mean at the MLE and dispersion matrix equal to the inverse of the
    # Fisher information."
    # This is most likely why hessian is used to compute standard errors.
    # TODO: Learn the Bernstein–von Mises theorem which explain in detail how
    # all this works
    stdh=sqrt(diag((solve(hess)))) #sqrt of negative values(
    varr=stdr^2; varh=stdh^2

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

    # selecting estimates of interest (i.e. alpha and betas) #
    stdrt=stdr[1:cur_variables_n]; stdht=stdh[1:cur_variables_n]
    varht=varh[1:cur_variables_n]; varrt=varr[1:cur_variables_n]

    # constructing the full vector of estimates #
    mty=rbind(1,mt)
    stdrt1=optimbase::zeros(variables_n,1); stdht1=optimbase::zeros(variables_n,1)
    varht1=optimbase::zeros(variables_n,1); varrt1=optimbase::zeros(variables_n,1)
    it1=0
    it=1
    for (it in 1:variables_n) {
      if (mty[it]==1) {
        it1=1+it1
        stdrt1[it]=stdrt[it1]
        stdht1[it]=stdht[it1]
        varht1[it]=varht[it1]
        varrt1[it]=varrt[it1]
      }
      else {
        stdrt1[it]=0
        stdht1[it]=0
        varht1[it]=0
        varrt1[it]=0
      }
    }

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
      modprob=postprob; modelid=row_ind; modpri=prior_model_prob; liks=exp(likelihood_max/n_entities); bics=bict
      stds=stdht1; stdsr=stdrt1; foutt=likelihood_max
    }
    else {
      modprob=rbind(modprob,postprob); modelid=rbind(modelid,row_ind); modpri=rbind(modpri,prior_model_prob)
      liks=rbind(liks,exp(likelihood_max/n_entities)); bics=rbind(bics,bict);
      stds=cbind(stds,stdht1); stdsr=cbind(stdsr,stdrt1); foutt=rbind(foutt, likelihood_max)
    }
  }

  list(prior_exp_model_size = prior_exp_model_size,
       prior_inc_prob = prior_inc_prob, variables_n = variables_n,
       modprob = modprob, modelid = modelid, modpri = modpri, liks = liks,
       bics = bics, stds = stds, stdsr = stdsr, foutt = foutt,
       bet = bet, mod = mod, pvarh = pvarh, pvarr = pvarr, fy = fy, fyt = fyt,
       ppmsize = ppmsize, cout = 0, nts = nts, pts = pts)
}
