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
no_of_cores = detectCores()
set.seed(23)
begin<-Sys.time()
prandom=0 #prandom=1 for theta random Ley&Steel09. prandom = 0 for theta fixed

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
regressors_n <- ncol(rawdata) - 4
variables_n <- regressors_n + 1

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

which_regs_bin_vectors <- replicate(regressors_n, 0:1, simplify = FALSE) %>%
  expand.grid()
which_regs_bin_vectors <-
  which_regs_bin_vectors[,order(ncol(which_regs_bin_vectors):1)]

for (row_ind in 1:nrow(which_regs_bin_vectors)) {
  mt <- as.matrix(t(which_regs_bin_vectors[row_ind, ]))
  out = (mt == 0)       # regressors out of the current model         
  kx=sum(mt); ky=kx+1  # number of regressors in the current model   
  
  #Z includes y0 and x0 as strictly exogenous variables
  Z <- R_df %>% filter(year == year0) %>% select(lag_gdp) %>% as.matrix()
  if (kx!=0) {
    X0j <- X0[,(mt==1)]
    Z <- cbind(Z ,X0j)
  }

  proj_matrix <- Z%*%solve(crossprod(Z))%*%t(Z)
  res_maker_matrix <- diag(n)-Z%*%solve(crossprod(Z))%*%t(Z)

  n_params_to_estimate <- 2*ky+t+1+(t^2+t-2)*regressors_n/2

  # t0in is the vector of initial values for the likelihood optimization   @
  t0in=0.5*ones(n_params_to_estimate,1)

   # ----------------------------------------------------------------------
   #                                                                        
   #                 LIKELIHOOD FUNCTION 1. BALIMLE.                       
   #                                                                        
   #  ----------------------------------------------------------------------
  
  lik<-function(t0in) {
    t0=t0in
    B0=diag(t+(t-1)*regressors_n)
    C0=zeros(t,ky)
    for (ii in 2:t) {
      B0[ii,(ii-1)]=-t0[1]      # SEM method B11 in paper
    }
    i2=1
    for (i1 in 1:regressors_n) {
      if (mt[i1]==0) {    # if x variable is not included 
        B0=B0 }
      else {              # if x variable is included     
        for (i11 in 2:t) {
          B0[i11,(t+1+(i11-2)*regressors_n+(i1-1))]=-t0[1+i2]      
        }      
    i2=i2+1
      }  
    }
  
  # C1 matrix 
    if (kx==0) {   
    C0[1,1]=t0[1]+t0[2]
      for (i3 in 2:t) {
        C0[i3,1]=t0[2]
      }
    }
    else {
    C0[1,1]=t0[1]+t0[1+ky]
    C0[1,(2:ncol(C0))]=t(t0[2:ky])+t(t0[(ky+2):(ky+1+kx)])
      for (i4 in 2:t) {
        C0[i4,]=t(t0[(ky+1):(ky+kx+1)])
      }
    }
  # C2 matrix has closed-form solutions 
  
    B110=B0[(1:t),(1:t)]
    B120=B0[(1:t),((t+1):ncol(B0))]
  
    o110=zeros(t,t)
    for (i5 in 1:t) {
        o110[i5,i5]=t0[2*ky+(i5+1)]^2
    }
    o110=o110+(t0[2*ky+1]^2)*(ones(t,t))  # Sigma 11 
    
     # Here, I split sigma 12 in the sum of two parts, phi's matrix and psi's upper triangular matrix 
    
     # phi's matrix 
    o120=zeros(t,(t-1)*regressors_n)
      for (i6 in 1:t) {
        for (i7 in 1:(t-1)) {    
          o120[i6,(1+(i7-1)*regressors_n):(i7*regressors_n)]=t(t0[(2*ky+t+2+(i7-1)*regressors_n):(2*ky+t+2+i7*regressors_n-1)])
        }
      }
    
    # psi's upper triangular matrix 
    o121=zeros(t,(t-1)*regressors_n)
    # as o121 is an upper triangular matrix, each subsequent row has 1 element less 
  
    seq=zeros(t,1)
    dseq=0
    for (iseq in 1:t) {
      seq[iseq]=dseq
      dseq=dseq+t-iseq
    }
      
    for (i8 in 1:t) {
      if (i8==t) {
        o121=o121 }
      else {
        o121[i8,((i8-1)*regressors_n+1):(ncol(o121))]=t(t0[(2*ky+t+2+(t-1+seq[i8])*regressors_n):(2*ky+t+2+(t-1+seq[i8+1])*regressors_n-1)])
      }
    }
    
    # Sigma 12  
    o120=o120+o121
    o210=t(o120)
    if (det(o110)<=0) {
      if (cout1==0) {
        cout=cout+1   
      }
      if (cout1<=250) {
        t0i=.5*ones(nrow(t0in),1)}
      if (cout1<=500) {
        t0i=t0in }
      else {
        t0i=as.vector(runif(rows(t0in))) }
      
      fact=1/(10^(round(log10(abs(t0i)),0)))
      t0in=fact*t0i
      likf=0
      cout1=cout1+1 
   }
    else {
      U10=t(B110%*%t(Y1)+B120%*%t(Y2)-C0%*%t(Z))  # Ui1 from the paper 
      H=crossprod(Y2-U10%*%solve(o110)%*%o120,res_maker_matrix)%*%(Y2-U10%*%solve(o110)%*%o120)
      likf=-(n/2)*log(det(o110))-(1/2)*sum(diag(solve(o110)%*%t(U10)%*%U10))-(n/2)*log(det(H/n))  # concentrated log-likelihood in the appendix 
    }
    return(-likf)
  }
  # parscale argument somehow (don't know yet how) changes step size during optimisation.
  # Most likely optimisation methods used in Gauss are scale-free and these used in R are not
  # TODO: search for methods (or implement methods) in R which are scale-free
  optimized<-optim(t0in,lik,method="BFGS",control = list(trace=2,maxit=10000, parscale = 0.05*t0in))
  theta<-optimized[[1]]; fout<-optimized[[2]]
  # theta returns optimized parameters, fout is the value of the function lik at the maximum
  # we now compute model-specific standard errors @
  
    #-----------------------------------------------------------------#
    #                        HESSIAN FUNCTION                         #
    #                        -----------------                        #
    #               this procedure takes as input a number            #
    #               computes the kxk (2-sided) hessian matrix         #
    #-----------------------------------------------------------------#
    
   myhess<-function(lik,theta) {
    x0<-theta
    k=nrow(x0) # x0 is theta, our model-specific optimized parameters #
    hessi=zeros(k,k)
    h=1e-3
    
    for (jc in 1:k) {
      for (jr in 1:jc) {
        x1=x0
        x2=x0
        x3=x0
        x4=x0
        x1[jr]=x0[jr]+h
        x1[jc]=x1[jc]+h
        x2[jr]=x0[jr]+h    
        x2[jc]=x2[jc]-h
        x3[jr]=x0[jr]-h
        x3[jc]=x3[jc]+h
        x4[jr]=x0[jr]-h
        x4[jc]=x4[jc]-h
        hessi[jr,jc]=(lik(x1)-lik(x2)-lik(x3)+lik(x4))/(4*h^2) # the second symmetric derivative has different formula #
        hessi[jc,jr]=hessi[jr,jc]
      }
    }
    return(hessi)
  } #for some reason, negative values are produced
  
  he=myhess(lik,theta)
  #he=hessian(lik,theta) #alternative methods
  #hess=(fdHess(theta,lik))
  #he=as.matrix(hess[[3]])
    #----------------------------------------------------------------------#
    #                                                                      #
    #                   LIKELIHOOD FUNCTION 1  (GRADIENT)                  #
    #                                                                      #
    #           this function construct the concentrated likelihood        #
    #           individual by individual and it gives the NX1 vector       #
    #           of individual log-likelihoods for the computation of       #
    #           the gradient and then the sandwich formula.                #
    #----------------------------------------------------------------------#
    
    
  likgra<-function(t0in) {
    t0=t0in
    likvec=zeros(n,1)
    
    B0=diag(t+(t-1)*regressors_n)
    C0=zeros(t,ky)
    for (ii in 2:t) {
      B0[ii,(ii-1)]=-t0[1]      # SEM method B11 in paper
    }
    
    i2=1
    for (i1 in 1:regressors_n) {
      if (mt[i1]==0) {    # if x variable is not included 
        B0=B0 }
      else {              # if x variable is included     
        for (i11 in 2:t) {
          B0[i11,(t+1+(i11-2)*regressors_n+(i1-1))]=-t0[1+i2]      
        }      
        i2=i2+1
      }  
    }
    
    
    # C1 matrix 
    if (kx==0) {   
      C0[1,1]=t0[1]+t0[2]
      for (i3 in 2:t) {
        C0[i3,1]=t0[2]
      }
    }
    else {
      C0[1,1]=t0[1]+t0[1+ky]
      C0[1,(2:ncol(C0))]=t(t0[2:ky])+t(t0[(ky+2):(ky+1+kx)])
      for (i4 in 2:t) {
        C0[i4,]=t(t0[(ky+1):(ky+kx+1)])
      }
    }
    # C2 matrix has closed-form solutions 

    B110=B0[(1:t),(1:t)]
    B120=B0[(1:t),((t+1):ncol(B0))]
    
    o110=zeros(t,t)
    for (i5 in 1:t) {
      o110[i5,i5]=t0[2*ky+(i5+1)]^2
    }
    o110=o110+(t0[2*ky+1]^2)*(ones(t,t))  # Sigma 11 
    
    # Here, I split sigma 12 in the sum of two parts, phi's matrix and psi's upper triangular matrix 
    
    # phi's matrix 
    o120=zeros(t,(t-1)*regressors_n)
    for (i6 in 1:t) {
      for (i7 in 1:(t-1)) {    
        o120[i6,(1+(i7-1)*regressors_n):(i7*regressors_n)]=t(t0[(2*ky+t+2+(i7-1)*regressors_n):(2*ky+t+2+i7*regressors_n-1)])
      }
    }
    
    # psi's upper triangular matrix 
    o121=zeros(t,(t-1)*regressors_n)
    # as o121 is an upper triangular matrix, each subsequent row has 1 element less 
    
    seq=zeros(t,1)
    dseq=0
    for (iseq in 1:t) {
      seq[iseq]=dseq
      dseq=dseq+t-iseq
    }
    
    for (i8 in 1:t) {
      if (i8==t) {
        o121=o121 }
      else {
        o121[i8,((i8-1)*regressors_n+1):(ncol(o121))]=t(t0[(2*ky+t+2+(t-1+seq[i8])*regressors_n):(2*ky+t+2+(t-1+seq[i8+1])*regressors_n-1)])
      }
    }
    
    # Sigma 12  
    o120=o120+o121
    o210=t(o120)
    
    U10=t(B110%*%t(Y1)+B120%*%t(Y2)-C0%*%t(Z))  # Ui1 from the paper 
    H=crossprod(Y2-U10%*%solve(o110)%*%o120,res_maker_matrix)%*%(Y2-U10%*%solve(o110)%*%o120)
    for (iter in 1:n) {
      u10i=as.matrix(U10[iter,])
      likvec[iter]=-(1/2)*log(det(o110))-(1/2)*log(det(H/n))-(1/2)*(t(u10i)%*%solve(o110)%*%u10i)
    }

  return(-likvec)

  }

  likgra_val<-likgra(theta)


  Gmat=gradient(likgra,theta)
  Imat=crossprod(Gmat)
    stdr=sqrt(diag(solve(he)%*%(Imat)%*%solve(he)))
    stdh=sqrt(diag((solve(he)))) #sqrt of negative values(
    varr=stdr^2; varh=stdh^2
    
    # storing results for the CURRENT model #    
    logl=(-fout-(ky/2)*(log(n*t)))/n
    bict=exp(logl)                              # integrated likelihood approximation           #
    # prior model probability (either random -Ley&Steel09- or fixed) #
    if (prandom == 1) {
      priorprobt=(gamma(1+kx))*(gamma(b+regressors_n-kx))    #theta random#
      }
    
    if (prandom == 0) {
      priorprobt=(pinc)^(kx)*(1-pinc)^(regressors_n-kx)      #theta fixed#
    }
    
    # posterior model probability  #
      postprob=priorprobt*bict
    
    # selecting estimates of interest (i.e. alpha and betas) #
      bt=theta[1:ky]; stdrt=stdr[1:ky]; stdht=stdh[1:ky]
    varht=varh[1:ky]; varrt=varr[1:ky]
    
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
      modprob=postprob; modelid=row_ind; modpri=priorprobt; liks=exp(-fout/n); bics=bict
      betas=bt1; stds=stdht1; stdsr=stdrt1; foutt=-fout
    }
     else {
       modprob=rbind(modprob,postprob); modelid=rbind(modelid,row_ind); modpri=rbind(modpri,priorprobt)
       liks=rbind(liks,exp(-fout/n)); bics=rbind(bics,bict); betas=cbind(betas,bt1)
       stds=cbind(stds,stdht1); stdsr=cbind(stdsr,stdrt1); foutt=rbind(foutt,(-fout))
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
