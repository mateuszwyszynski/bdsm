# The following packages are not installed automatically by renv,
# because the panels package does not depend on them.
# You should install them before running this script.
library(tidyverse)
library(readxl)
devtools::load_all()

set.seed(23)
begin<-Sys.time()
dilution <- 1
dil_power <- 1/2

# Load dataset
rawdata <- readxl::read_excel("balimle-dataset.xlsx")

#    VARIABLES IN RAWDATA
#    1.FDI  2.FDIlag  3.EI  4.LLF  5.EX  6.SW  7.RES  8.LOPW  9.INT  10.RI
#-----------------------------------------------------------------------------------------
rawdata=rawdata[,1:8]   #I select the regressors of interest @
year0 <- min(rawdata$year)
varlist<- c("FDI","EI","LLF","EX", "SW")

data_with_no_lagged_col <- rawdata %>%
  join_lagged_col(gdp, lag_gdp, year, country, 10)

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
  df <- df %>% mutate(across(!(year:country), ~ c(scale(.))))
  df
  csddata_df <- df %>% group_by(year) %>%
    summarise(country = country,
              across(!country, function(x) c(scale(x, scale = FALSE)))) %>%
    arrange(country) %>% ungroup()
}

# R_df <- liml_data_prep(rawdata)
data_prepared <- liml_data_prep(data_with_no_lagged_col)

Y1 <- SEM_dep_var_matrix(
  df = data_prepared, timestamp_col = year,
  entity_col = country, dep_var_col = gdp
)

Y2 <- data_prepared %>%
  SEM_regressors_matrix(timestamp_col = year,
                        entity_col = country,
                        regressors = c(ish, pop, pgrw, sed))
Z <- data_prepared %>%
  exogenous_matrix(timestamp_col = year, entity_col = country,
                   dep_var_col = gdp, 10,
                   regressors_subset = c(ish, pop, pgrw, sed))

bma_result <- SEM_bma(R_df = data_prepared, dep_var_col = gdp,
                      timestamp_col = year, timestep = 10,
                      entity_col = country, projection_matrix_const = TRUE)

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
for (jt in 1:bma_result$variables_n) {
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

final<-list(
  result, rbind(
    paste("Prior Mean Model Size=", bma_result$prior_exp_model_size),
    paste("Prior Inclusion Probability=", bma_result$prior_inc_prob),
    paste("Posterior Mean Model Size=", popmsize)
    ), (the_end-begin), t(betas), t(stds), t(stdsr), idprob
  )
names(final)<-c(" 1.- RESULTS "," 2.- FURTHER INFORMATION "," 3.- COMPUTATION TIME"," 4.- ALL BETAS (each row is a different model)",
                " 5.- ALL STD. ERRORS (each row is a different model)"," 6.- ALL ROBUST STD. ERRORS (each row is a different model)",
                " 7.- MODELS INFO ")
final
