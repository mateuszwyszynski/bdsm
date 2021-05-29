# this program runs BALIMLE approach
# It estimates all possible models.
library(tidyverse)
library(optimbase) #used for zeros and ones
library(rootSolve) #used for gradient function
library(numDeriv)  #used for hessian?
library(mcmc)
library(matlib)   #used for inverse
library(nlme)
library("parallel")
library(rje)

source('R/SEM_likelihood.R')
source('R/hessian.R')

no_of_cores = detectCores()
set.seed(23)
begin<-Sys.time()
prandom=0 #prandom=1 for optimised_params random Ley&Steel09. prandom = 0 for optimised_params fixed

#---------------------------------------------------------------------------------
#		   	                     LOADING THE DATASET
#---------------------------------------------------------------------------------

row=292; column=11; t=4; n=73
rawdata<-readxl::read_excel("balimle-dataset.xlsx")

#    VARIABLES IN RAWDATA
#    1.FDI  2.FDIlag  3.EI  4.LLF  5.EX  6.SW  7.RES  8.LOPW  9.INT  10.RI
#-----------------------------------------------------------------------------------------
rawdata=rawdata[,1:8]   #I select the regressors of interest @
year0 <- min(rawdata$year)
varlist<- c("FDI","EI","LLF","EX", "SW")
regressors <- rev(colnames(rawdata)[-1:-4])
regressors_n <- ncol(rawdata) - 4
variables_n <- regressors_n + 1
variables <- rev(colnames(rawdata)[-1:-3])
n_entities <- 73

#' Prepare data for LIML estimation
#'
#' @description
#' This is a function which prepares data for Limited Information Maximum
#' Likelihood (LIML) Estimation. Following operations are performed:
#'
#' 1. Data standarisation
#' 2. Cross-sectional demeaning of variables
#' 3. Organisation of data for the LIML estimation
#'
#' @param df Dataframe with data that should be prepared for LIML estimation
liml_data_prep <- function(df){
  df <- df %>% mutate(across(!(year:country), scale))

  csddata_df <- df %>% group_by(year) %>%
    summarise(country = country,
              across(!country, function(x) scale(x, scale = FALSE))) %>%
    arrange(country) %>% ungroup()
}

R_df <- liml_data_prep(rawdata)

# Dependent variables for the t periods - matrix of size N x t
Y1 <- R_df %>% select(year, country, gdp) %>%
  pivot_wider(names_from = year, values_from = gdp) %>%
  select(!country) %>% as.matrix()

# Regressors for the t periods:
# X0 - for the first year
# Y2 - for the remaining years
Y2 <- R_df %>% select(!gdp & !lag_gdp) %>% filter(year != year0) %>%
  pivot_wider(names_from = year, values_from = !country & !year) %>%
  select(!country) %>%
  select(order(as.numeric(gsub("[^0-9]+", "", colnames(.))))) %>% as.matrix()

X0 <- R_df %>% filter(year == year0) %>%
  select(!(year:lag_gdp)) %>% as.matrix()

# ---------------------------------------------------------------------------------
# 		               SOME PRELIMINAR OBJECTS BALIMLE APPROACH
# ---------------------------------------------------------------------------------
#
#***                 PRIOR STRUCTURES for THE PRIOR MODEL SIZE                  ***
# ---------------------------------------------------------------------------------

pmsize=regressors_n/2           # prior expected model size, options:
    # 1. for "SDM" priors pmsize=3
    # 2. for "FLS" priors pmsize=regressors_n/2
pinc=pmsize/regressors_n
b=(regressors_n-pmsize)/pmsize   # parameter for beta (random) distribution of the prior inclusion probability

# ***                     STORAGE OBJECTS                                       ***
# ---------------------------------------------------------------------------------

mod=zeros(variables_n,1); bet=zeros(variables_n,1);
pvarh=zeros(variables_n,1); pvarr=zeros(variables_n,1);
fy=zeros(variables_n,1); fyt=0; ppmsize=0; cout=0

#---------------------------------------------------------------------------------
#  		               LOOP COVERING FULL MODEL SPACE
#---------------------------------------------------------------------------------

regressors_subsets <- powerSet(regressors)
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

  #Z includes y0 and x0 as strictly exogenous variables
  Z <- R_df %>% filter(year == year0) %>%
    select(lag_gdp, regressors_subset) %>% as.matrix()

  proj_matrix <- Z%*%solve(crossprod(Z))%*%t(Z)
  res_maker_matrix <- diag(n) - proj_matrix

  n_params_to_estimate <- 2*cur_variables_n+t+1+(t^2+t-2)*regressors_n/2

  periods_n <- t

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

  cur_Y2 <- R_df %>% select(year, country, regressors_subset) %>%
    filter(year != year0) %>%
    pivot_wider(names_from = year, values_from = !country & !year) %>%
    select(!country) %>%
    select(order(as.numeric(gsub("[^0-9]+", "", colnames(.))))) %>% as.matrix()

  o <- orig_sigma_matrix(t0 = t0in, t = t,
                         cur_variables_n = cur_variables_n,
                         regressors_n = regressors_n)

  print(o)
}

#---------------------------------------------------------------------------------
#  		               LOOP COVERING FULL MODEL SPACE
#---------------------------------------------------------------------------------

regressors_subsets <- powerSet(regressors)
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

  #Z includes y0 and x0 as strictly exogenous variables
  Z <- R_df %>% filter(year == year0) %>%
    select(lag_gdp, regressors_subset) %>% as.matrix()

  proj_matrix <- Z%*%solve(crossprod(Z))%*%t(Z)
  res_maker_matrix <- diag(n) - proj_matrix

  n_params_to_estimate <- 2*cur_variables_n+t+1+(t^2+t-2)*regressors_n/2

  periods_n <- t

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

  cur_Y2 <- R_df %>% select(year, country, regressors_subset) %>%
    filter(year != year0) %>%
    pivot_wider(names_from = year, values_from = !country & !year) %>%
    select(!country) %>%
    select(order(as.numeric(gsub("[^0-9]+", "", colnames(.))))) %>% as.matrix()

  # parscale argument somehow (don't know yet how) changes step size during optimisation.
  # Most likely optimisation methods used in Gauss are scale-free and these used in R are not
  # TODO: search for methods (or implement methods) in R which are scale-free
  optimized <- optim(t0in, SEM_likelihood, n_entities = n_entities,
                     cur_Y2 = cur_Y2, Y1 = Y1, Y2 = Y2, Z = Z,
                     res_maker_matrix = res_maker_matrix,
                     periods_n = periods_n, regressors_n = cur_regressors_n,
                     phis_n = phis_n, psis_n = psis_n,
                     cur_variables_n = cur_variables_n,
                     tot_regressors_n = regressors_n,
                     method="BFGS",
                     control = list(trace=2, maxit = 10000, fnscale = -1,
                                    parscale = 0.05*t0in))
  optimised_params <- optimized[[1]]
  likelihood_max <- optimized[[2]]

  hess <- hessian(SEM_likelihood, theta = optimised_params,
                n_entities = n_entities,
                cur_Y2 = cur_Y2, Y1 = Y1, Y2 = Y2, Z = Z,
                res_maker_matrix = res_maker_matrix,
                periods_n = periods_n, regressors_n = cur_regressors_n,
                phis_n = phis_n, psis_n = psis_n,
                cur_variables_n = cur_variables_n,
                tot_regressors_n = regressors_n)

  likgra_val<-likgra(optimised_params)

  Gmat=gradient(likgra,optimised_params)
  Imat=crossprod(Gmat)
    stdr=sqrt(diag(solve(hess)%*%(Imat)%*%solve(hess)))
    stdh=sqrt(diag((solve(hess)))) #sqrt of negative values(
    varr=stdr^2; varh=stdh^2

    # storing results for the CURRENT model #
    logl=(likelihood_max-(cur_variables_n/2)*(log(n*t)))/n
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
    bt1=zeros(variables_n,1)
    stdrt1=zeros(variables_n,1); stdht1=zeros(variables_n,1)
    varht1=zeros(variables_n,1); varrt1=zeros(variables_n,1)
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
      modprob=postprob; modelid=row_ind; modpri=priorprobt; liks=exp(likelihood_max/n); bics=bict
      betas=bt1; stds=stdht1; stdsr=stdrt1; foutt=likelihood_max
    }
     else {
       modprob=rbind(modprob,postprob); modelid=rbind(modelid,row_ind); modpri=rbind(modpri,priorprobt)
       liks=rbind(liks,exp(likelihood_max/n)); bics=rbind(bics,bict); betas=cbind(betas,bt1)
       stds=cbind(stds,stdht1); stdsr=cbind(stdsr,stdrt1); foutt=rbind(foutt, likelihood_max)
     }
}

#                         end of loop covering full model space                        #
#--------------------------------------------------------------------------------------#
popmsize=ppmsize/fyt
modprob1=modprob/sum(modprob)

idprob=as.data.frame(cbind(modelid,modprob1,modpri,bics))
names(idprob)<-c("model","postprob","riorprob","bics")
row.names(idprob)<-NULL

# computing posterior moments CONDITIONAL on inclusion
postprobinc=fy/fyt
postmean=bet/fy
varrleamer=(pvarr/fy)-postmean^2
varhleamer=(pvarh/fy)-postmean^2
poststdr=sqrt(varrleamer)
poststdh=sqrt(varhleamer)
tr=postmean/poststdr
th=postmean/poststdh

# computing UNCONDITIONAL posterior moments
upostmean = postmean * postprobinc
uvarrleamer = (varrleamer + (postmean^2))*postprobinc - (upostmean^2)
uvarhleamer = (varhleamer + (postmean^2))*postprobinc - (upostmean^2)
upoststdr=sqrt(uvarrleamer)
upoststdh=sqrt(uvarhleamer)

# computing percentage of significant coeff estimates
nts=t(nts)
pts=t(pts)
for (jt in 1:variables_n) {
  ntss=na.omit(nts[,jt]); ptss=na.omit(pts[,jt]); nsig=ntss<(-1.96); psig=ptss>1.96
  if (jt==1) {
    negper=mean(nsig); posper=mean(psig)
  }
  else {
    negper=rbind(negper,mean(nsig)); posper=rbind(posper,mean(psig))
  }
}



result=as.data.frame(cbind(varlist,postprobinc,postmean,poststdh,poststdr,upostmean,upoststdh,upoststdr))
names(result)<-c("varname","postprob","pmean","std","stdR","unc_pmean","unc_std","unc_stdR")
the_end=Sys.time()

final<-list(result, rbind(paste("Prior Mean Model Size=",pmsize),paste("Prior Inclusion Probability=",pinc),
     paste("Posterior Mean Model Size=", popmsize)),(the_end-begin),
     t(betas),t(stds),
     t(stdsr),idprob)
names(final)<-c(" 1.- RESULTS "," 2.- FURTHER INFORMATION "," 3.- COMPUTATION TIME"," 4.- ALL BETAS (each row is a different model)",
                " 5.- ALL STD. ERRORS (each row is a different model)"," 6.- ALL ROBUST STD. ERRORS (each row is a different model)",
                " 7.- MODELS INFO ")
final
