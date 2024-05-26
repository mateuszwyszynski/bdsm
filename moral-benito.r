# The following packages are not installed automatically by renv,
# because the panels package does not depend on them.
# You should install them before running this script.
library(tidyverse)
library(readxl)
devtools::load_all()

set.seed(20)
begin<-Sys.time()

data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

regressors <- regressor_names(data_prepared, year, country, gdp)

bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = economic_growth_ms,
                          projection_matrix_const = TRUE)

bet <- bma_result$bet
pvarh <- bma_result$pvarh
pvarr <- bma_result$pvarr
fy <- bma_result$fy
fyt <- bma_result$fyt
ppmsize <- bma_result$ppmsize
cout <- bma_result$cout
nts <- bma_result$nts
pts <- bma_result$pts

popmsize=ppmsize/fyt
models_prob_normalized <-
  bma_result$models_posterior_prob / sum(bma_result$models_posterior_prob)

idprob <- as.data.frame(cbind(models_prob_normalized,
                              bma_result$models_prior_prob, bma_result$bics))
names(idprob)<-c("postprob", "priorprob")
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

result=as.data.frame(cbind(regressors,postprobinc,postmean,poststdh,poststdr,upostmean,upoststdh,upoststdr))
names(result)<-c("varname","postprob","pmean","std","stdR","unc_pmean","unc_std","unc_stdR")
the_end=Sys.time()

final<-list(
  result, rbind(
    paste("Prior Mean Model Size=", bma_result$prior_exp_model_size),
    paste("Prior Inclusion Probability=", bma_result$prior_inc_prob),
    paste("Posterior Mean Model Size=", popmsize)
    ), (the_end-begin), idprob
  )
names(final)<-c(" 1.- RESULTS "," 2.- FURTHER INFORMATION "," 3.- COMPUTATION TIME",
                " 4.- MODELS INFO ")
final
