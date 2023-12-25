# this program runs BALIMLE approach
# It estimates all possible models.
library(tidyverse)

# To install the development version of the package run:
# devtools::install()
# Then you should be able to load the package as any other package.
# More details on how to develop packages can be found here: https://r-pkgs.org
library(panels)

# Run:
# renv::restore()
# if there are any problems.
# This should recreate the most recent working environment.

set.seed(23)
begin<-Sys.time()
prandom=0 #prandom=1 for optimised_params random Ley&Steel09. prandom = 0 for optimised_params fixed
dilution <- 1
dil_power <- 1/2

#---------------------------------------------------------------------------------
#		   	                     LOADING THE DATASET
#---------------------------------------------------------------------------------

row=292; column=11;
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
periods_n <- 4

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

Y1 <- SEM_dep_var_matrix(
  df = R_df, timestamp_col = year, entity_col = country,
  dep_var_col = gdp, start_time = year0
)

Y2 <- R_df %>%
  SEM_regressors_matrix(timestamp_col = year, entity_col = country,
                        regressors = c(ish, sed, pgrw, pop),
                        start_time = year0)

Z <- R_df %>%
  SEM_exogenous_matrix(year, year0, lag_gdp,
                       regressors_subset = c('ish', 'sed', 'pgrw', 'pop'))

res_maker_matrix <- residual_maker_matrix(Z)

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

regressors_subsets <- rje::powerSet(regressors)

bma_result <- SEM_bma(regressors_subsets = regressors_subsets, R_df = R_df,
                      variables_n = variables_n, regressors_n = regressors_n,
                      periods_n = periods_n, timestamp_col = year,
                      year0 = year0, lagged_col = lag_gdp, entity_col = country,
                      Y1 = Y1, Y2 = Y2, res_maker_matrix = res_maker_matrix,
                      prandom = prandom, n_entities = n_entities, b = b,
                      pinc = pinc,
                      projection_matrix_const = TRUE)

modprob <- bma_result$modprob
modelid <- bma_result$modelid
modpri <- bma_result$modpri
liks <- bma_result$liks
bics <- bma_result$bics
betas <- bma_result$betas
stds <- bma_result$stds
stdsr <- bma_result$stdsr
foutt <- bma_result$foutt
bet <- bma_result$bet
mod <- bma_result$mod
pvarh <- bma_result$pvarh
pvarr <- bma_result$pvarr
fy <- bma_result$fy
fyt <- bma_result$fyt
ppmsize <- bma_result$ppmsize
cout <- bma_result$cout
nts <- bma_result$nts
pts <- bma_result$pts

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
