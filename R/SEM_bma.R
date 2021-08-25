SEM_bma <- function(regressors_subsets, R_df, variables_n, regressors_n,
                    periods_n, timestamp_col, year0, lagged_col, entity_col, Y1,
                    Y2, res_maker_matrix, prandom, n_entities, b, pinc) {
  mod <- optimbase::zeros(variables_n,1)
  bet <- optimbase::zeros(variables_n,1)
  pvarh <- optimbase::zeros(variables_n,1)
  pvarr <- optimbase::zeros(variables_n,1)
  fy <- optimbase::zeros(variables_n,1)
  fyt <- 0
  ppmsize <- 0
  cout <- 0

  which_regs_bin_vectors <- replicate(regressors_n, 0:1, simplify = FALSE) %>%
    expand.grid()
  which_regs_bin_vectors <-
    which_regs_bin_vectors[,order(ncol(which_regs_bin_vectors):1)]

  row_ind <- 0
  for (regressors_subset in regressors_subsets) {
    regressors_subset <- rev(regressors_subset)
    row_ind <- row_ind + 1
    mt <- as.matrix(t(which_regs_bin_vectors[row_ind, ]))
    out = (mt == 0)       # regressors out of the current model
    cur_regressors_n <- sum(mt)
    cur_variables_n <- cur_regressors_n+1

    cur_Z <- R_df %>%
      SEM_exogenous_matrix({{ timestamp_col }}, year0, {{ lagged_col }},
                           regressors_subset)

    # Initial parameter values for optimisation
    alpha <- 0.5
    phi_0 <- 0.5
    err_var <- 0.5
    dep_vars <- rep(0.5, periods_n)
    beta <- rep(0.5, cur_regressors_n)
    phi_1 <- rep(0.5, cur_regressors_n)
    phis_n <- regressors_n*(periods_n - 1)
    phis <- rep(0.5, phis_n)
    psis_n <- regressors_n*periods_n*(periods_n - 1)/2
    psis <- rep(0.5, psis_n)

    t0in <- matrix(c(alpha, beta, phi_0, phi_1, err_var, dep_vars, phis, psis))

    cur_Y2 <- R_df %>%
      SEM_regressors_matrix(timestamp_col = {{ timestamp_col }},
                            entity_col = {{ entity_col }},
                            regressors = regressors_subset, start_time = year0)

    data <- list(Y1 = Y1, Y2 = Y2, cur_Y2 = cur_Y2, Z = cur_Z,
                 res_maker_matrix = res_maker_matrix)

    # parscale argument somehow (don't know yet how) changes step size during optimisation.
    # Most likely optimisation methods used in Gauss are scale-free and these used in R are not
    # TODO: search for methods (or implement methods) in R which are scale-free
    optimized <- stats::optim(t0in, SEM_likelihood, data = data,
                              periods_n = periods_n,
                              tot_regressors_n = regressors_n,
                              in_regressors_n = cur_regressors_n,
                              phis_n = phis_n, psis_n = psis_n, method="BFGS",
                              control = list(trace=2, maxit = 10000,
                                             fnscale = -1,
                                             parscale = 0.05*t0in))
    optimised_params <- optimized[[1]]
    likelihood_max <- optimized[[2]]

    hess <- hessian(SEM_likelihood, theta = optimised_params, data = data,
                    periods_n = periods_n, tot_regressors_n = regressors_n,
                    in_regressors_n = cur_regressors_n,
                    phis_n = phis_n, psis_n = psis_n)

    likelihood_per_entity <-
      SEM_likelihood(optimised_params, data = data, per_entity = TRUE,
                     periods_n = periods_n, tot_regressors_n = regressors_n,
                     in_regressors_n = cur_regressors_n, phis_n = phis_n,
                     psis_n = psis_n)

    Gmat <- rootSolve::gradient(SEM_likelihood, optimised_params, data = data,
                                per_entity = TRUE, periods_n = periods_n,
                                tot_regressors_n = regressors_n,
                                in_regressors_n = cur_regressors_n,
                                phis_n = phis_n, psis_n = psis_n)
    Imat=crossprod(Gmat)
    stdr=sqrt(diag(solve(hess)%*%(Imat)%*%solve(hess)))
    stdh=sqrt(diag((solve(hess)))) #sqrt of negative values(
    varr=stdr^2; varh=stdh^2

    # storing results for the CURRENT model #
    logl=(likelihood_max-(cur_variables_n/2)*(log(n_entities*periods_n)))/n_entities
    bict=exp(logl)                              # integrated likelihood approximation           #
    # prior model probability (either random -Ley&Steel09- or fixed) #
    if (prandom == 1) {
      priorprobt=(gamma(1+cur_regressors_n))*(gamma(b+regressors_n-cur_regressors_n))    #optimised_params random#
    }

    if (prandom == 0) {
      priorprobt=(pinc)^(cur_regressors_n)*(1-pinc)^(regressors_n-cur_regressors_n)      #optimised_params fixed#
    }

    # posterior model probability  #
    postprob=priorprobt*bict

    # selecting estimates of interest (i.e. alpha and betas) #
    bt=optimised_params[1:cur_variables_n]; stdrt=stdr[1:cur_variables_n]; stdht=stdh[1:cur_variables_n]
    varht=varh[1:cur_variables_n]; varrt=varr[1:cur_variables_n]

    # constructing the full vector of estimates #
    mty=rbind(1,mt)
    bt1=optimbase::zeros(variables_n,1)
    stdrt1=optimbase::zeros(variables_n,1); stdht1=optimbase::zeros(variables_n,1)
    varht1=optimbase::zeros(variables_n,1); varrt1=optimbase::zeros(variables_n,1)
    it1=0
    it=1
    for (it in 1:variables_n) {
      if (mty[it]==1) {
        it1=1+it1
        bt1[it]=bt[it1]
        stdrt1[it]=stdrt[it1]
        stdht1[it]=stdht[it1]
        varht1[it]=varht[it1]
        varrt1[it]=varrt[it1]
      }
      else {
        bt1[it]=0         # if the regressor is not in the model, 0 #
        stdrt1[it]=0
        stdht1[it]=0
        varht1[it]=0
        varrt1[it]=0
      }
    }


    # calculating the percentage of significant regressions #
    ptr=bt1/stdht1
    ntr=bt1/stdht1
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
    bet=bet+postprob*bt1
    pvarr=pvarr+(postprob*varrt1+postprob*(bt1*bt1))         # as in Leamer (1978) #
    pvarh=pvarh+(postprob*varht1+postprob*(bt1*bt1))         # as in Leamer (1978) #

    # here we store model-specific diagnostics and estimates (BICs, likelihoods, betas...) #
    if (row_ind==1) {
      modprob=postprob; modelid=row_ind; modpri=priorprobt; liks=exp(likelihood_max/n_entities); bics=bict
      betas=bt1; stds=stdht1; stdsr=stdrt1; foutt=likelihood_max
    }
    else {
      modprob=rbind(modprob,postprob); modelid=rbind(modelid,row_ind); modpri=rbind(modpri,priorprobt)
      liks=rbind(liks,exp(likelihood_max/n_entities)); bics=rbind(bics,bict); betas=cbind(betas,bt1)
      stds=cbind(stds,stdht1); stdsr=cbind(stdsr,stdrt1); foutt=rbind(foutt, likelihood_max)
    }
  }

  list(modprob = modprob, modelid = modelid, modpri = modpri, liks = liks,
       bics = bics, betas = betas, stds = stds, stdsr = stdsr, foutt = foutt,
       bet = bet, mod = mod, pvarh = pvarh, pvarr = pvarr, fy = fy, fyt = fyt,
       ppmsize = ppmsize, cout = 0)
}
