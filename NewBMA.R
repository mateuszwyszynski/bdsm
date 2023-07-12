# DEPENDENCIES
library(openxlsx)
library(matlib)
library(optimbase)
library(gridExtra)
library(grid)

# FUNCTIONS TO BE USED BY THE PACKAGE

## 1. comb - function that returns the number of combinations of x variables out of n variables 

comb=function(n,x){factorial(n)/(factorial(n-x)*factorial(x))} #function that calculates cimbinations 

###########################################################################################################################

## 2. ols - function that calculates all relevant ols statistics (coefficients, standard errors, R^2, Likelihood)

ols=function(y,x,const){ # function that provides ols estimates and additional statistics
  m<-nrow(y) # number of rows in the data
  if (identical(x,0)){ # when you add "0" as argument model with no variables and a constant will be calculated
    x<-ones(m,1) # creation of the vector of ones
    const=0} # we make sure another vector of ones is not gonna be added below
  # we need to calculate Dilution before we add vector of ones to the regressors data matrix
  n<-ncol(as.matrix(x))
  if (n==1){ # CONDITION for the model with one regressor
  Diluntion<-1 # we assign 1 as dilution to be used in dilution prior  
  } # the end of the CONDITION for the model with one regressor
  else{ # CONDITION for the model with more than one regressor  
    Diluntion<-(det(cor(x)))  # we assign dilution to be used in dilution prior (George 2010)
  }# end of the CONDITION for the model with more than one regressor 
  if (const==1){x<-cbind(ones(m,1),x)} # CONDITION for adding vector of ones to the model
  n<-ncol(as.matrix(x)) # number of regressors including the constant if parameter "const"=1
  if (det(t(x)%*%x)==0){stop("Determinant of X'X=0")} # # test if X'X can be inverted: Warning that X'X cannot be inverted
  if (n==1){ #test for a model with just a constant and no rgressors
    betas=(t(x)%*%x)^(-1)%*%t(x)%*%y # we estimate the model parameters - version with (no regressors and a constant) or one regressor
  }else{
    betas=inv(t(x)%*%x)%*%t(x)%*%y # we estimate the model parameters - general version
  } 
  y_hat=x%*%betas # we obtain the theoretical values
  res=y-y_hat # we calculate the residuals
  SSR=t(res)%*%res # sum of squares of the residuales
  df=m-n # we calculate the number of the degrees of freedom
  sigma2=(t(res)%*%res)/df # we calculate error variance
  sigma=sigma2^(0.5) # we obtain standard error of the regression
  if (n==1){ #test for a model with just a constant and no rgressors
    var_B=as.numeric(sigma2)*(t(x)%*%x)^-1 # we calculate variances of the coefficients  - version with (no regressors and a constant) or one regressor
  }else{
    var_B=as.numeric(sigma2)*inv(t(x)%*%x) # we calculate variances of the coefficients - for general version
  }
  se_B=diag(var_B)^(0.5) # we calculate standard erors of the coefficients
  y_m=mean(y) #we calculate the mean value of the dependent variable
  SST=t(y-y_m*ones(m,1))%*%(y-y_m*ones(m,1)) # total sum of squares of the regression
  R2<-1-(SSR/SST) # calculation of R^2
  like<-m^(-n/2)*SSR^(-m/2) # value of the likelihood function (Leamer, 1978)
  out <- list(betas,se_B,as.numeric(like),as.numeric(R2),as.numeric(df),as.numeric(Diluntion)) # creates a list of objects needed for modelSpace function
}# THE END OF THE ols FUNCTION

###########################################################################################################################

## 3. modelSpace - function that builds model space using ols statistics

#PARAMETERS TO BE DEFINED BE THE USER
M<-8 # Parameter: maximal size of the model to be considered in BMA
const=1 # Parameter: 0 - no constant in all the models
        #            1 - a constant in all the models

# loading data file
data<-read.xlsx("~/Desktop/VD determinants/Regional regression.xlsx",sheet="S_1981-2019",colNames = T, rowNames = F)

#THE modelSpace FUNCTION starts here
modelSpace=function(data,M,const,dil.Par=0.5){
dil.Par<-dil.Par  
  # collecting data characteristics
  m<-nrow(data) # number of rows in the data 
  n<-ncol(data) # number of columns in the data
  K<-n-1 # number of regressors
  
  if (M>K){# CONDITION about what to do if the user set M that is higher than K (M>K)
    # we tell the user that we are setting M=K 
    message("M>K - maximum number of regressors cannot be bigger than total number of regressors. We set M=K and continiue :)")
    M=K # we set M=K
  }# end of the CONDITION about what to do if the user set M that is higher than K (M>K)
  
  Var_names<-names(data) # names of the variables
  x_names<-Var_names[2:n] # names of the regressors
  y<-as.matrix(data[,1]) # data on the regressant 
  x<-as.matrix(data[,2:n]) # data on the regressors
  
  # s indicates the smalest model size under the consideration (0 variables or 1 variable)
  if (const==1){s=0}else {s=1}
  
  MS=0 # MS a variable representing the size of the model space (total number of models)
  for (k in s:M){# at this LOOP we add all combinations of regressors up models with M variables 
    c=comb(K,k) # number of models of the size k out of K regressors
    MS=MS+c} # this sum adds up all the models for each possible model size 
  
  # we bild a table for all ols statistics
  ols_results=zeros(MS,M+2*K+4+2*const) #THIS NEEDS TO BE FIXED
  
  # ms - model index
  ms=0 #starting value for counting the numbe of the model
  
  # k indicates the number of variables in the model
  for (k in s:M){ #at this LOOP we create all possible model sizes
    c=comb(K,k) # number of models of the size k out of K regressors
    if (k==0){ # CONDITION for the special case of a model with no variables and a constant
      ols1_model<-ols(y=y,x=0,con=1) # estimation of the model with a constant and no regressors
      ms=ms+1 #we change the number of the model
      ols_results[ms,M+1]=as.numeric(ols1_model[1]) # extraction of the coefficients 
      ols_results[ms,M+K+2]=as.numeric(ols1_model[2]) # extraction of the standard errors
      ols_results[ms,M+2*K+3]=as.numeric(ols1_model[3]) # here we extract value of the Likelihood function
      ols_results[ms,M+2*K+4]=as.numeric(ols1_model[4]) # here we extract R2
      ols_results[ms,M+2*K+5]=as.numeric(ols1_model[5]) # here we extract the number of degrees of freedom 
      ols_results[ms,M+2*K+6]=as.numeric(ols1_model[6]) # here we extract information for dilution prior
    } # end of the CONDITION for the special case of a model with no variables
    else{# CONDITION for the general case - model with regressors
      models<-as.matrix(combn(1:K,k)) #generates all possible models of the size k out of K regressors
      #i indicates the number of the model from the gour of the models with k our of K regressors
      for (i in 1:c){ # at this LOOP we are estimating the individual models with with one or more regressors (k>=0)
        mod<-models[,i] # we extraxt the indices of the regressors to be used in a given model ms
        x_ms=zeros(m,1) # we crate an artificial vector so we can perform the binding in the next LOOP
        for (t in 1:k){ # at this LOOP we colect the regressors for the model t
          h<-mod[t]  # here we extract an index of the regressor used for estimation
          x_mi=x[,h] # we extract regressor
          x_ms=cbind(x_ms,x_mi) #we bind all k regressors together 
        } # end of the LOOP that colects the regressors for the model t
        x_ms<-x_ms[,-1] # we delete an artificial vector from the regressor matrix
        ms=ms+1 # we update the index of the model
        #here I neet to introduce a condition of whether we use constant or not
        model_ms<-ols(y,x_ms,const) #estimation of the model ms
        ols_results[ms,M+2*K+2*const+1]=as.numeric(model_ms[3]) #here we extract value of the Likelihood function
        ols_results[ms,M+2*K+2*const+2]=as.numeric(model_ms[4]) #here we extract R2
        ols_results[ms,M+2*K+2*const+3]=as.numeric(model_ms[5]) #here we extract the number of degrees of freedom
        ols_results[ms,M+2*K+2*const+4]=as.numeric(model_ms[6]) # here we extract information for dilution prior
        for (p in 1:(k+const)){# LOOP performs extraction of the coefficients and the standard errors from the model
          ols_results[ms,M+p]=as.numeric(model_ms[[1]][[p]]) # extraction of the coefficients
          ols_results[ms,M+K+const+p]=as.numeric(model_ms[[2]][[p]]) # extraction of the standard errors
          if (const==1){# CONDITION - the case of a model WITH a constant
            if (p<k+const){ols_results[ms,p]=mod[p]} # here we extract indices of the used regressors
          }else if (const==0){# CONDITION - the case of a model WITHOUT a constant
            ols_results[ms,p]=mod[p]} # here we extract indices of the used regressors
        }# end of the LOOP performs extraction of the coefficients and the standard errors from the model
      } # end of the LOOP that estimates the individual models with with one or more regressors (k>=0)
    } # end of the CONDITION for the general case - model with regressors
  } # end of the LOOP that creates all model sizes

  out<-list(x_names,ols_results,ms,M,K,const) # we create a modelSpace object (mS object) - a list with:
  # 1. x_names - vector with names of the regressors
  # 2. ms - size of the mode space (the number of the last estimated model)
  # 3. M - maximum number of regressors in a model
  # 4. K- total number of regressors
   
} # end of the function modelSpace

###########################################################################################################################

modelS<-modelSpace(data,M=3,const=1) #EXAMSPLE EXECUTION OF THE modelSpace FUNCTION
x_names<-modelS[[1]][] # EXAMPLE OF THE RESULTS FROM THE modelSpace OBJECT: REGRESSORS NAMES
ols_results<-modelS[[2]][] # EXAMPLE OF THE RESULTS FROM THE modelSpace OBJECT: MAIN RESULTS 


## 4. Posterior - function that creates BMA statistics and other posterior objects using modelSpace object (mS object)

#PARAMETERS TO BE DEFINED BE THE USER
constMS=1 # Parameter: 0 - model with no variables is EXCLUDED from the model space
          #            1 - model with no variables is INCLUDED in the model space
dilution=1 # Parameter: 0 - NO application of a dilution prior
           #            1 - application of a dilution prior 
dil.Par=0.5 # Parameter associated with dilution prior (George 2010)

#THE Posterior FUNCTION starts here
Posterior=function(modelSpace,constMS,dilution=0,dil.Par=0.5){

# Extraction of the elements of the mS object
x_names<-modelS[[1]][] # extraction of the regressors names from the mS object
ols_results<-modelS[[2]][] # extraction of the ols resuls (the model space) from the mS object
MS<-modelS[[3]][1] # extraction of the total bumber of models
M<-modelS[[4]][1] # extraction of the maximum number of regressors in the model
K<-modelS[[5]][1] # extraction of the total number of regressors
const<-modelS[[6]][1] # a binary variabe indicating whether mS object incued a model with a constant:
                      # 1 - a model with just a constant is INCLUDED in model space
                      # 0 - a model with just a constant is EXCLUDED from model space

# Conditions associated with inclusion/exclusion of a model with no regressors
if(const==0){ # CONDITION preventing from inlcusion of a model with no constant if it was not included in the previoues step (in prepaering mS object)
  constMS==0 # we set constMS=0 - model with no regressors is not considered
  message("Model with regressors was not included in modelSpace object. We set constMS=0 and continiue :)")
  } # end of the CONDITION preventing from inlcusion of a model with no constant if it was not included in the previoues step (in prepaering mS object)

if (constMS==0&const==1){# CONDITION for inclusion (1) or exclusion (0) of the model with no regressors
  ols_results<-ols_results[-1,]
  MS=MS-1 # We need to lowe the total number of models as we eliminate the model with no regressors
} # end of the CONDITION for inclusion (1) or exclusion (0) of the model with no regressors

# Dividing ols results into relevant part
Reg_ID<-ols_results[,1:M] # we extract vector indices
betas<-ols_results[,(M+1):(M+K+const)] # we extract ceofficients
VAR<-ols_results[,(M+K+1+const):(M+2*K+2*const)]^2 # we extract standard errors and change to variance (VAR)
R2<-ols_results[,M+2*K+1+2*const] # we extract R^2
Like<-ols_results[,M+2*K+2+2*const] # we extract value of the likelihood function (expression proportional to likelihood)
DF<-ols_results[,M+2*K+3+2*const] # we extract the number of degrees of freedom
dilut<-ols_results[,M+2*K+4+2*const] # we extract expression associated with dilution prior
# Model Priors
uniform_models<-(1/MS)*ones(MS,1) # we create a vector of unifrom model priors (ON MODELS)
uniform_sizes<-zeros(M+constMS,1) # we create a vector to store probabilities (ON MODEL SIZES) for model sizes under uniform model prior
random_models<-zeros(MS,1) # we create a vector to store probabilities (ON MODELS) for the case of equal proabilitie in all model sizes
random_sizes<-(1/(M+constMS))*ones(M+constMS,1) # we create a vector of model size priors (ON MODEL SIZES) for the case of equal proabilitie in all model sizes

sizes<-zeros(M+constMS,1) #we create vector to store number of models in a given model size

if(constMS==0){s=1}else{s=0} # CONDITION about the inclusion/exclusion of the model with no regressors

for (k in s:M){# at this LOOP we add all combinations of regressors up models with M variables 
  sizes[k+1,1]<-comb(K,k) # number of models of the size k out of K regressors
  } # this sum adds up all the models for each possible model size 

ind<-cumsum(sizes) # we create a vector with the number of models in each model size category

h=constMS # here we create atrificial variable to perfom vector combination through rbind function
for (i in 1:(M+constMS)){ # this LOOP creates random prior for individual models and uniform prior for model sizes
    q=ones(sizes[i,1],1)*(1/(M+constMS))*(1/sizes[i,1]) # we create vectors of probabilities for models of a given size
    h=rbind(h,q) # we bind vectors with model random prior probabilities
    if (i==1){uniform_sizes[i,1]=(1/MS)} # we collect probabilities for different model sizes (UNIFORM prior): the case of the model with no regressors
    else{uniform_sizes[i,1]=sum(uniform_models[(ind[i-1]+1):ind[i],1])} # we collect probabilities for different model sizes (UNIFORM prior): the case of models with regressors
} # end of the LOOP creates random prior for individual models and uniform prior for model sizes

random_models<-as.vector(h[-1,]) # we create a vector to store probabilities (ON MODELS) for the case of equal proabilitie in all model sizes

# Posterior model probabilities (PMP) and R^2 weights of individual models
PMP_uniform<-as.matrix((uniform_models*Like)/sum(uniform_models*Like)) # calculation of PMPs based on uniform prior
PMP_R2_uniform<-as.matrix((uniform_models*R2)/sum(uniform_models*R2)) # calculation of R^2 weights based on uniform prior
PMP_random<-as.matrix((random_models*Like)/sum(random_models*Like))  # calculation of PMPs based on random prior
PMP_R2_random<-as.matrix((random_models*R2)/sum(random_models*R2)) # calculation of R^2 weights based on random prior

# CONDITION for dillution prior
if (dilution==1){# begining of the CONDITION for dilution
  # PLEASE CHECK IF THIS WORKS PROPERLY
  if (is.null(dil.Par)==1){dil.Par=0.5} # CONDITION for setting the default value of dil.Par
  PMP_uniform<-PMP_uniform*(dilut^dil.Par) # we add diltuion component to PMPs measures based on uniform prior
  PMP_R2_uniform<-PMP_R2_uniform*(dilut^dil.Par) # we add diltuion component tof R^2 weights based on uniform prior
  PMP_random<-PMP_random*(dilut^dil.Par)  # we add diltuion component to PMPs measures based on random prior
  PMP_R2_random<- PMP_R2_random*(dilut^dil.Par) # we add diltuion component tof R^2 weights based on random prior
}# THE END of the CONDITION for dilution

# Creating products of PMPs and betas and stds^2=VAR for calculation of posterior statistics
betas_PMP_uniform<-zeros(MS,K+constMS) # matrix to store products of betas and PMP_uniform
betas_R2_uniform<-zeros(MS,K+constMS) # matrix to store products of betas and PMP_R2_uniform
betas_PMP_random<-zeros(MS,K+constMS) # matrix to store products of betas and PMP_random
betas_R2_random<-zeros(MS,K+constMS) # matrix to store products of betas and PMP_R2_random
VAR_PMP_uniform<-zeros(MS,K+constMS) # matrix to store products of VAR and PMP_uniform
VAR_R2_uniform<-zeros(MS,K+constMS) # matrix to store products of VAR and PMP_R2_uniform
VAR_PMP_random<-zeros(MS,K+constMS) # matrix to store products of VAR and PMP_random
VAR_R2_random<-zeros(MS,K+constMS) # matrix to store products of VAR and PMP_R2_random

for (j in 1:(M+constMS)){ # at this LOOP we mutiply columns of betas and VAR by posterior model measures
  betas_PMP_uniform[,j]=PMP_uniform*betas[,j] # matrix to store products of betas and PMP_uniform
  betas_R2_uniform[,j]=PMP_R2_uniform*betas[,j] # matrix to store products of betas and R2_uniform
  betas_PMP_random[,j]=PMP_random*betas[,j] # matrix to store products of betas and PMP_random
  betas_R2_random[,j]=PMP_R2_random*betas[,j] # matrix to store products of betas and R2_random
  VAR_PMP_uniform[,j]=PMP_uniform*VAR[,j] # matrix to store products of VAR and PMP_uniform
  VAR_R2_uniform[,j]=PMP_R2_uniform*VAR[,j] # matrix to store products of VAR and R2_uniform
  VAR_PMP_random[,j]=PMP_random*VAR[,j] # matrix to store products of VAR and PMP_random
  VAR_R2_random[,j]=PMP_R2_random*VAR[,j] # matrix to store products of VAR and R2_random
} # end of the LOOP at which we mutiply columns of betas and VAR by posterior model measures

# POSTERIOR INCLUSION PROBABILITIES (PIP) AND POSTERIOR MEANS (PM)
PM_PMP_uniform<-zeros(K,1) # matrix to store PMs under uniform model prior
PM_R2_uniform<-zeros(K,1) # matrix to store PMs obtained with R^2 under uniform model prior
PM_PMP_random<-zeros(K,1) # matrix to store PMs under random model prior
PM_R2_random<-zeros(K,1) # matrix to store PMs obtained with R^2 under random model prior
PIP_PMP_uniform<-zeros(K,1) # matrix to store PIPs under uniform model prior
PIP_R2_uniform<-zeros(K,1) # matrix to store PIPs obtained with R^2 under uniform model prior
PIP_PMP_random<-zeros(K,1) # matrix to store PIPs under random model prior
PIP_R2_random<-zeros(K,1) # matrix to store PIPs obtained with R^2 under random model prior
Plus_PMP_uniform<-zeros(K,1) # matrix to store P(+) under uniform model prior
Plus_R2_uniform<-zeros(K,1) # matrix to store P(+) obtained with R^2 under uniform model prior
Plus_PMP_random<-zeros(K,1) # matrix to store P(+) under random model prior
Plus_R2_random<-zeros(K,1) # matrix to store P(+) obtained with R^2 under random model prior

Plus_prep<-zeros(MS,K) # matrix to store ones (1s) representing positive value of the coefficient
Denominator<-zeros(K,1) #matrix that counts total number of models with a given variable
beta_k<-zeros(MS,K)
se<-zeros(MS,K)

# Calculation of posterior inclusion probabilities (PIP), posterior means (PM), and preparations for EBA
for (k in 1:K){ # at this LOOP we move along regressors to collect elements to appropriate sums
  for (i in 1:MS){ # at this LOOP we move along all models to find the models with a given regressor k
    for (t in 1:M){ # at this LOOP we go through Reg_ID to find regressor k index
      if (Reg_ID[i,t]==k){# CONDITION for finding regressor k
        # PREPARATION OF THE OBJECTS FOR EXTREME BOUNDS ANALYSIS (EBA)
        beta_k[i,k]=betas[i,k+const] # we colect all the estimated coefficients for regressor k
        se[i,k]=(VAR[i,k+const])^(0.5) # we colect all the estimated standard errors for regressor k
        Denominator[k,1]=Denominator[k,1]+1 # we get the number of models with a given regressor - the same for every regressors
        # POSTERIOR INCLUSION PROBABILITIES (PIP)
        PIP_PMP_uniform[k,1]=PMP_uniform[i,1]+PIP_PMP_uniform[k,1] # we sum PMPs for models including a given regressor
        PIP_R2_uniform[k,1]=PMP_R2_uniform[i,1]+PIP_R2_uniform[k,1] # we sum PMPs for models including a given regressor
        PIP_PMP_random[k,1]=PMP_random[i,1]+PIP_PMP_random[k,1] # we sum PMPs for models including a given regressor
        PIP_R2_random[k,1]=PMP_R2_random[i,1]+PIP_R2_random[k,1] # we sum PMPs for models including a given regressor
        # POSTERIOR MEANS
        PM_PMP_uniform[k,1]=betas_PMP_uniform[i,k+const]+PM_PMP_uniform[k,1] # we sum products of PMPs with a given coefficient
        PM_R2_uniform[k,1]=betas_R2_uniform[i,k+const]+PM_R2_uniform[k,1] # we sum products of PMPs with a given coefficient
        PM_PMP_random[k,1]=betas_PMP_random[i,k+const]+PM_PMP_random[k,1] # we sum products of PMPs with a given coefficient
        PM_R2_random[k,1]=betas_R2_random[i,k+const]+PM_R2_random[k,1] # we sum products of PMPs with a given coefficient
        # % OF POSITVE BETAS
        if (betas[i,k+const]>0){# CONDITION for counting models with POSITIVE coefficients
          Plus_prep[i,k]=1 # we set 1 to Plus_prep for models in which regressors have positive coefficients
                           # it is latter used to calculate %(+) - percentage of positive coefficients
          # here we calculate posterior probability of a positive sign P(+) for the relevant expression
          # see Doppelhofer and Weeks (2009) p. 216 for the details of the expression
          # below is the case for the coefficient on the regressor being positive
          Plus_PMP_uniform[k,1]=Plus_PMP_uniform[k,1]+PMP_uniform[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_R2_uniform[k,1]=Plus_R2_uniform[k,1]+PMP_R2_uniform[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_PMP_random[k,1]=Plus_PMP_random[k,1]+PMP_random[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_R2_random[k,1]=Plus_R2_random[k,1]+PMP_R2_random[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
        } # end of the CONDITION for counting models with POSITIVE coefficients
        else if (betas[i,k+const]<0){# CONDITION for counting models with NEGATIVE coefficients
          # here we calculate posterior probability of a positive sign P(+) for the relevant expression
          # see Doppelhofer and Weeks (2009) p. 216 for the details of the expression
          # below is the case for the coefficient on the regressor being POSITIVE
          Plus_PMP_uniform[k,1]=Plus_PMP_uniform[k,1]+1-PMP_uniform[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_R2_uniform[k,1]=Plus_R2_uniform[k,1]+1-PMP_R2_uniform[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_PMP_random[k,1]=Plus_PMP_random[k,1]+1-PMP_random[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
          Plus_R2_random[k,1]=Plus_R2_random[k,1]+1-PMP_R2_random[i,1]*pt(betas[i,k+const]/(VAR[i,k+const])^0.5,df =DF[i])
        } # end of the CONDITION for counting models with NEGATIVE coefficients 
      } # the end of the CONDITION for finding regressor k
    } # the end of the LOOP we go through Reg_ID to find regressor k index
  } # end of the LOOP that moves along all models to find the models with a given regressor k
} # end of the LOOP that moves along regressors to collect elements to appropriate sums

# Plus denotes percentage of positive betas 
Plus<-colSums(Plus_prep)/Denominator

# POSTERIOR STANDARD DEVIATION (PSD) PREP
VAR_PMP_prep_uniform<-zeros(K,1)
VAR_R2_prep_uniform<-zeros(K,1)
VAR_PMP_prep_random<-zeros(K,1)
VAR_R2_prep_random<-zeros(K,1)

# Calculations of posterior standard deviations (PSD)
for (k in 1:K){ # at this LOOP we move along regressors to collect elements to appropriate sums
  for (i in 1:MS){ # at this LOOP we move along all models to find the models with a given regressor k
    for (t in 1:M){ # at this LOOP we go through Reg_ID to find regressor k index
      if (Reg_ID[i,t]==k){# CONDITION for finding regressor k
        # Below we add variance of regressor k multiplied by posterior model probability to the square of the difference
        # between coefficient estimated in model i and Posterior mean of regressor k also multiplied by posterior model
        # probability. As a result we get posterior variances we are doing to use to calculate posterior standard 
        # deviations (PSD) under the assumption of uniform model priorand random model prior (using PMP and R2 in both cases)
        VAR_PMP_prep_uniform[k,1]=VAR_PMP_prep_uniform[k,1]+((betas[i,k]-PM_PMP_uniform[k,1])^2)*PMP_uniform[i,1]
        VAR_R2_prep_uniform[k,1]=VAR_R2_prep_uniform[k,1]+((betas[i,k]-PM_R2_uniform[k,1])^2)*PMP_R2_uniform[i,1]
        VAR_PMP_prep_random[k,1]=VAR_PMP_prep_random[k,1]+((betas[i,k]-PM_PMP_random[k,1])^2)*PMP_random[i,1]
        VAR_R2_prep_random[k,1]=VAR_R2_prep_random[k,1]+((betas[i,k]-PM_R2_random[k,1])^2)*PMP_R2_random[i,1]
      } # the end of the CONDITION for finding regressor k
    } # the end of the LOOP we go through Reg_ID to find regressor k index
  } # end of the LOOP that moves along all models to find the models with a given regressor k
} # end of the LOOP that moves along regressors to collect elements to appropriate sums

# Calculation of posterior standard deviations (PSD)
PSD_PMP_uniform<-VAR_PMP_prep_uniform^0.5
PSD_R2_uniform<-VAR_R2_prep_uniform^0.5
PSD_PMP_random<-VAR_PMP_prep_random^0.5
PSD_R2_random<-VAR_R2_prep_random^0.5

# EXTREME BOUDS ANALYSIS (EBA)
beta_max<-as.matrix(apply(betas[,(1+const):(K+const)],2,max)) # we find the highest values of the coefficients
beta_min<-as.matrix(apply(betas[,(1+const):(K+const)],2,min)) # we find the lowest values of the coefficients
beta_mean<-as.matrix(apply(betas[,(1+const):(K+const)],2,mean)) # we find mean values of the coefficients
upper<-zeros(K,1) # matrix to store upper bounds
lower<-zeros(K,1) # matrix to store lower bounds

for (k in 1:K){ # at this LOOP we calculate upper and lower bounds
u<-which(beta_k[,k]==max(beta_k[,k])) # finding the highest value of the coefficient
if (length(u)>1){u<-1} # CONDITION TO avoid the problem of repeating values - zeros
upper[k,1]<-beta_k[u,k]+2*se[u,k] # calculation of the upper bound
l<-which(beta_k[,k]==min(beta_k[,k])) # finding the lowest value of the coefficient
if (length(l)>1){l<-1}  # CONDITION TO avoid the problem of repeating values - zeros
lower[k,1]<-beta_k[l,k]-2*se[u,k] # calculation of the lower bound
} # end of the LOOP where we calculate upper and lower bounds

# Calculation of the posterior objects (PIP, PM, and PSD) for the constant
if (const==1){ # CONDITION checking if we should include the constant in the calculations
  
  const_betas_PMP_uniform<-PMP_uniform*betas[,1] # we create a vector of products of constants and PMP_uniform
  const_betas_R2_uniform<-PMP_R2_uniform*betas[,1] # we create a vector of products of constants and R2_uniform
  const_betas_PMP_random<-PMP_random*betas[,1] # we create a vector of products of constants and PMP_random
  const_betas_R2_random<-PMP_R2_random*betas[,1] # we create a vector of products of constants and R2_random
  const_VAR_PMP_uniform<-PMP_uniform*VAR[,1] # we create a vector of products of VAR of constants and PMP_uniform
  const_VAR_R2_uniform<-PMP_R2_uniform*VAR[,1] # we create a vector of products of VAR of constants and R2_uniform
  const_VAR_PMP_random<-PMP_random*VAR[,1] # we create a vector of products of VAR of constants and PMP_random
  const_VAR_R2_random<-PMP_R2_random*VAR[,1] # we create a vector of products of VAR of constants and R2_random 
  
  # Calculation of posterior means (PM)
  const_PM_PMP_uniform<-0 # object to store PMs under uniform model prior
  const_PM_R2_uniform<-0 # object to store PMs obtained with R^2 under uniform model prior
  const_PM_PMP_random<-0 # object to store PMs under random model prior
  const_PM_R2_random<-0 # object to store PMs obtained with R^2 under random model 
  Plus_const<-0
  Plus_const_PMP_uniform<-0 # object to store P(+) under uniform model prior
  Plus_const_R2_uniform<-0 # object to store P(+) obtained with R^2 under uniform model prior
  Plus_const_PMP_random<-0 # object to store P(+) under random model prior
  Plus_const_R2_random<-0 # object to store P(+) obtained with R^2 under random model 
  
  for (i in 1:MS){ # at this LOOP we move along all models
    # POSTERIOR MEANS
    const_PM_PMP_uniform<-const_betas_PMP_uniform[i,1]+const_PM_PMP_uniform # we sum products of PMPs with a given coefficient
    const_PM_R2_uniform<-const_betas_R2_uniform[i,1]+const_PM_R2_uniform # we sum products of PMPs with a given coefficient
    const_PM_PMP_random<-const_betas_PMP_random[i,1]+const_PM_PMP_random # we sum products of PMPs with a given coefficient
    const_PM_R2_random<-const_betas_R2_random[i,1]+const_PM_R2_random # we sum products of PMPs with a given coefficient
    # % OF POSITVE BETAS
    if (betas[i,1]>0){ # CONDITION for counting models with POSITIVE coefficients 
      Plus_const<-(1/MS)+Plus_const # ath this step we add percentage associated 
      Plus_const_PMP_uniform<-Plus_const_PMP_uniform+PMP_uniform[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_R2_uniform<-Plus_const_R2_uniform+PMP_R2_uniform[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_PMP_random<-Plus_const_PMP_random+PMP_random[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_R2_random<-Plus_const_R2_random+PMP_R2_random[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    } # the end of the CONDITION for counting models with POSITIVE coefficients
    else if (betas[i,1]<0){# CONDITION for counting models with NEGATIVE coefficients
      Plus_const_PMP_uniform<-Plus_const_PMP_uniform+1-PMP_uniform[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_R2_uniform<-Plus_const_R2_uniform+1-PMP_R2_uniform[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_PMP_random<-Plus_const_PMP_random+1-PMP_random[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
      Plus_const_R2_random<-Plus_const_R2_random+1-PMP_R2_random[i,1]*pt(betas[i,1]/(VAR[i,1])^0.5,df =DF[i])
    } # the end of the CONDITION for counting models with NEGATIVE coefficients
  } # end of the LOOP that moves along all models
  # Calculations of posterior standard deviations (PSD)
  const_VAR_PMP_prep_uniform<-0 # object to store PSDs under uniform model prior
  const_VAR_R2_prep_uniform<-0 # object to store PSDs obtained with R^2 under uniform model prior
  const_VAR_PMP_prep_random<-0 # object to store PSDs under random model prior
  const_VAR_R2_prep_random<-0 # object to store PSDs obtained with R^2 under random model prior
  for (i in 1:MS){ # at this LOOP we move along all models
    # Below we add variance of the constant multiplied by posterior model probability to the square of the difference
    # between the constant estimated in model i and Posterior mean of the constantk also multiplied by posterior model
    # probability. As a result we get posterior variances we are doing to use to calculate posterior standard 
    # deviation (PSD) under the assumption of uniform model priorand random model prior (using PMP and R2 in both cases)
    const_VAR_PMP_prep_uniform<-const_VAR_PMP_prep_uniform+((betas[i,1]-const_PM_PMP_uniform)^2)*PMP_uniform[i,1]
    const_VAR_R2_prep_uniform<-const_VAR_R2_prep_uniform+((betas[i,1]-const_PM_R2_uniform)^2)*PMP_R2_uniform[i,1]
    const_VAR_PMP_prep_random<-const_VAR_PMP_prep_random+((betas[i,1]-const_PM_PMP_random)^2)*PMP_random[i,1]
    const_VAR_R2_prep_random<-const_VAR_R2_prep_random+((betas[i,1]-const_PM_R2_random)^2)*PMP_R2_random[i,1]
  } # end of the LOOP that moves along all models
  
  # Calculation of posterior standard deviations (PSD)
  const_PSD_PMP_uniform<-const_VAR_PMP_prep_uniform^0.5
  const_PSD_R2_uniform<-const_VAR_R2_prep_uniform^0.5
  const_PSD_PMP_random<-const_VAR_PMP_prep_random^0.5
  const_PSD_R2_random<-const_VAR_R2_prep_random^0.5
 
  # Extreme Bound Analysis (EBA)
  const_max<-max(betas[,1]) # finding the highest constant
  const_min<-min(betas[,1]) # finding the lowest constant
  const_mean<-mean(betas[,1]) # finding mean constant
  
  u<-which(betas[,1]==max(betas[,1])) # finding the highest value of the coefficient
  if (length(u)>1){u<-1} # CONDITION TO avoid the problem of repeating values - zeros
  const_upper<-betas[u,1]+2*se[u,1] # calculation of the upper bound
  l<-which(betas[,1]==min(betas[,1])) # finding the lowest value of the coefficient
  if (length(l)>1){l<-1}  # CONDITION TO avoid the problem of repeating values - zeros
  const_lower<-betas[l,1]-2*se[l,1] # calculation of the lower bound
  
  # At this step we add results for the constant to posterior statistics 
  PM_PMP_uniform<-as.matrix(as.numeric(rbind(const_PM_PMP_uniform,PM_PMP_uniform)))
  PM_R2_uniform<-as.matrix(as.numeric(rbind(const_PM_R2_uniform,PM_R2_uniform)))
  PM_PMP_random<-as.matrix(as.numeric(rbind(const_PM_PMP_random,PM_PMP_random)))
  PM_R2_random<-as.matrix(as.numeric(rbind(const_PM_R2_random,PM_R2_random)))
  PIP_PMP_uniform<-as.matrix(as.numeric(rbind(1,PIP_PMP_uniform)))
  PIP_R2_uniform<-as.matrix(as.numeric(rbind(1,PIP_R2_uniform)))
  PIP_PMP_random<-as.matrix(as.numeric(rbind(1,PIP_PMP_random)))
  PIP_R2_random<-as.matrix(as.numeric(rbind(1,PIP_R2_random)))
  PSD_PMP_uniform<-as.matrix(as.numeric(rbind(const_PSD_PMP_uniform,PSD_PMP_uniform)))
  PSD_R2_uniform<-as.matrix(as.numeric(rbind(const_PSD_R2_uniform,PSD_R2_uniform)))
  PSD_PMP_random<-as.matrix(as.numeric(rbind(const_PSD_PMP_random,PSD_PMP_random)))
  PSD_R2_random<-as.matrix(as.numeric(rbind(const_PSD_R2_random,PSD_R2_random)))
  lower<-as.matrix(as.numeric(rbind(const_lower,lower)))
  beta_min<-as.matrix(as.numeric(rbind(const_min,beta_min)))
  beta_mean<-as.matrix(as.numeric(rbind(const_mean,beta_mean)))
  beta_max<-as.matrix(as.numeric(rbind(const_max,beta_max)))
  upper<-as.matrix(as.numeric(rbind(const_upper,upper)))
  Plus<-as.matrix(as.numeric(rbind(Plus_const,Plus)))
  Plus_PMP_uniform<-as.matrix(as.numeric(rbind(Plus_const_PMP_uniform,Plus_PMP_uniform)))
  Plus_R2_uniform<-as.matrix(as.numeric(rbind(Plus_const_R2_uniform,Plus_R2_uniform)))
  Plus_PMP_random<-as.matrix(as.numeric(rbind(Plus_const_PMP_random,Plus_PMP_random)))
  Plus_R2_random<-as.matrix(as.numeric(rbind(Plus_const_R2_random,Plus_R2_random)))
  c_name="Constant"
  x_names<-c(c_name,x_names)
  
} # end of the CONDITION checking if we should include the constant in the calculations 

# Calculation of posterior mean to posterior standard deviation (PM/PSD) ratios 
Ratio_PMP_uniform<-PM_PMP_uniform/PSD_PMP_uniform
Ratio_R2_uniform<-PM_R2_uniform/PSD_R2_uniform
Ratio_PMP_random<-PM_PMP_random/PSD_PMP_random
Ratio_R2_random<-PM_R2_random/PSD_R2_random

# Objects to be available to the user
PMP_uniform_table<-cbind(PIP_PMP_uniform,PM_PMP_uniform,PSD_PMP_uniform,Ratio_PMP_uniform,Plus_PMP_uniform)
rownames(PMP_uniform_table)<-x_names
colnames(PMP_uniform_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
R2_uniform_table<-cbind(PIP_R2_uniform,PM_R2_uniform,PSD_R2_uniform,Ratio_R2_uniform,Plus_R2_uniform)
rownames(R2_uniform_table)<-x_names
colnames(R2_uniform_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
PMP_random_table<-cbind(PIP_PMP_random,PM_PMP_random,PSD_PMP_random,Ratio_PMP_random,Plus_PMP_random)
rownames(PMP_random_table)<-x_names
colnames(PMP_random_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
R2_random_table<-cbind(PIP_R2_random,PM_R2_random,PSD_R2_random,Ratio_R2_random,Plus_R2_random)
rownames(R2_random_table)<-x_names
colnames(R2_random_table)<-c("PIP","PM","PSD","PM/PSD","P(+)")
EBA<-cbind(lower,beta_min,beta_mean,beta_max,upper,Plus)
rownames(EBA)<-x_names
colnames(EBA)<-c("Lower bound","Minimum","Mean","Maximum","Upper bound","%(+)")

# objects to be available for other functions
PIPs<-cbind(PIP_PMP_uniform,PIP_PMP_random,PIP_R2_uniform,PIP_R2_random)
forJointness<-cbind(Reg_ID,PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random)
forBestModels<-cbind(Reg_ID,PMP_uniform,PMP_random,PMP_R2_uniform,PMP_R2_random,betas,VAR,R2)

# creation of the list of Posterior objects (Post objects)
out<-list(PMP_uniform_table,PMP_random_table,EBA,R2_uniform_table,R2_random_table,
          x_names,M,K,MS,const,PIPs,forJointness,forBestModels) # we create a Posterior object (Post object) - a list with:
# 1. OBJECT 1: Table with results with PMP under uniform model prior
# 2. OBJECT 2: Table with results with PMP under random model prior
# 3. OBJECT 3: Table with results of Extreme Bounds Analysis
# 4. OBJECT 4: Table with results with R^2 under uniform model prior
# 5. OBJECT 5: Table with results with R^2 under random model prior
# 6. OBJECT 6: x_names - vector with names of the regressors - to be used by othe functions
# 7. OBJECT 7: M - maximum number of regressors in a model
# 8. OBJECT 8: K - total number of regressors
# 9. OBJECT 9: MS - size of the mode space
# 10. OBJECT 10: const - paramter informing whether the model space includes the nostant or not
# 11. OBJECT 11: Table with PIP under different model priors for Jointness function: PIPs
# 12. OBJECT 12: Table with model IDs and PMPs for Jointness function: forJointnes
# 13. OBJECT 13: Table with model IDs, PMPs, coefficients, variances, and R^2 for bestModels function: forBestModels
# 14. OBJECT 14: HERE OBJECT FOR modelSizes FUNCTION NEEDS TO BE ADDED

}# THE END of the Posterior FUNCTION

###############################################################################################################################

Post<-Posterior(modelSpace,constMS=1)

Post[6]# x_names
Post[7]# M
Post[8]# K
Post[9]# MS
Post[10]# const
Post[11]# PIPs
Post[12]# forJointnes
Post[13]# forBestModels
Post[14]# FOR modelSizes FUNCTION

## 5. Jointness - function that calculates jointness measures for the regressors using output for the Posteriot object (Post object)

#PARAMETERS TO BE DEFINED BE THE USER

# Parameters that tell which jointness measures should be placed:
# 1. above the diagaonal
above=1
# 2. below the diagaonal
below=2
# in the output matrix
# The parameter takes the followin values:
    # 1 - for PMP_uniform
    # 2 - for PMP_random
    # 3 - for R2_uniform
    # 4 - for R2_random

# Parameters for choosing the measure of jointness
measure="HCGHM"
# The parameter takes the followin values:
    # HCGHM - for Hofmarcher et al. (2018) measure 
    # LS - for Ley & Steel (2007) measure
    # DW - for Doppelhofer & Weeks (2009) measure 
    # PPI - for posterior probability of including both variables

# If HCGHM measure is chosen then rho parameter needs to be chosen. Default is 0.5
rho=0.5 # Parameter associated with HCGHM measure of jointness (Hofmarcher et al. 2018)

#THE Jointness FUNCTION starts here
Jointness=function(Post,above=1,below=2,measure="HCGHM",rho=0.5){
  
# Extraction of the elements of the mS object

x_names<-Post[[6]]# we extract vector with names of the regressors from Posterior object (Post object)
x_names<-x_names[-1]
M<-Post[[7]]# we extract M - maximum number of regressors in a model from Posterior object (Post object)
K<-Post[[8]]# we extract K - total number of regressors from Posterior object (Post object)
MS<-Post[[9]]# we extract MS - size of the mode space from Posterior object (Post object)
const<-Post[[10]]# we extract const parameter from Posterior object (Post object)
PIPs<-Post[[11]] # we extract matrix with PIPs under different model prior assumptions
                 # from Posterior object (Post object)
if(const==1){ # CONDITION for models with a constant 
PIPs<-PIPs[-1,] # we delete a row for constants
} # the end of the CONDITION for models with a constant
forJointness<-Post[[12]] # we extract forJointness matrix from Posterior object (Post object)
Reg_ID<-forJointness[,1:M] # we extract regressor ID
PMP_uniform<-as.matrix(forJointness[,M+1]) # we extract PMP under unifrom prior
PMP_random<-as.matrix(forJointness[,M+2]) # we extract PMP under unifrom prior
PMP_R2_uniform<-as.matrix(forJointness[,M+3]) # we extract PMP under unifrom prior
PMP_R2_random<-as.matrix(forJointness[,M+4]) # we extract PMP under unifrom prior

# Information about pairs of regressors: number of pairs, list of all possible regressors pairs
c=comb(K,2) # number of pairs e.i. jointness measures
Pairs<-as.matrix(combn(1:K,2)) # a list of all possible pairs of regressors 

# introducing notation a matrices to store posterior objects
PMP_uniform_a_b<-zeros(c,1) # a_b - P(a and b) for PMP under uniform prior
PMP_uniform_Na_b<-zeros(c,1) # Na_b - P(NOT a and b) for PMP under uniform prior
PMP_uniform_a_Nb<-zeros(c,1) # a_Nb - P(a and NOT b) for PMP under uniform prior
PMP_uniform_Na_Nb<-zeros(c,1) # Na_Nb - P(NOT a and NOT b) for PMP under uniform prior
PMP_random_a_b<-zeros(c,1) # a_b - P(a and b) for PMP under random prior
PMP_random_Na_b<-zeros(c,1) # Na_b - P(NOT a and b) for PMP under random prior
PMP_random_a_Nb<-zeros(c,1) # a_Nb - P(a and NOT b) for PMP under random prior
PMP_random_Na_Nb<-zeros(c,1) # Na_Nb - P(NOT a and NOT b) for PMP under random prior
R2_uniform_a_b<-zeros(c,1) # a_b - P(a and b) for R2 under uniform prior
R2_uniform_Na_b<-zeros(c,1) # Na_b - P(NOT a and b) for R2 under uniform prior
R2_uniform_a_Nb<-zeros(c,1) # a_Nb - P(a and NOT b) for R2 under uniform prior
R2_uniform_Na_Nb<-zeros(c,1) # Na_Nb - P(NOT a and NOT b) for R2 under uniform prior
R2_random_a_b<-zeros(c,1) # a_b - P(a and b) for R2 under random prior
R2_random_Na_b<-zeros(c,1) # Na_b - P(NOT a and b) for R2 under random prior
R2_random_a_Nb<-zeros(c,1) # a_Nb - P(a and NOT b) for R2 under random prior
R2_random_Na_Nb<-zeros(c,1) # Na_Nb - P(NOT a and NOT b) for R2 under random prior

for (j in 1:c){# this LOOP goes trough all the regressors pairs
a<-Pairs[1,j] # ID of the first regressor
b<-Pairs[2,j] # ID of the second regressor
  for (i in 1:MS){# this LOOP goes through all the models
    aIN<-0 # setting the initial value before the LOOP
    bIN<-0 # setting the initial value before the LOOP
    for (t in 1:M){# this LOOP goes trough model IDs
      if (Reg_ID[i,t]==a){aIN<-1} # CONDITION for presence of the regressors a in the model
      if (Reg_ID[i,t]==b){bIN<-1} # CONDITION for presence of the regressors b in the model
    }# the end of the LOOP that goes trough model IDs
    if (aIN==1&bIN==1){# CONDNION for both variables being included in the model
      PMP_uniform_a_b[j,1]=PMP_uniform_a_b[j,1]+PMP_uniform[i,1]
      PMP_random_a_b[j,1]=PMP_random_a_b[j,1]+PMP_random[i,1]
      R2_uniform_a_b[j,1]=R2_uniform_a_b[j,1]+PMP_R2_uniform[i,1]
      R2_random_a_b[j,1]=R2_random_a_b[j,1]+PMP_R2_random[i,1]
      }# the end of the CONDNION for both variables being included in the model
    else if(aIN==0&bIN==1) {# CONDNION for only regressor b being included in the model
      PMP_uniform_Na_b[j,1]=PMP_uniform_Na_b[j,1]+PMP_uniform[i,1]
      PMP_random_Na_b[j,1]=PMP_random_Na_b[j,1]+PMP_random[i,1]
      R2_uniform_Na_b[j,1]=R2_uniform_Na_b[j,1]+PMP_R2_uniform[i,1]
      R2_random_Na_b[j,1]=R2_random_Na_b[j,1]+PMP_R2_random[i,1]
    }# the end of the CONDNION for only regressor b being included in the model
    else if(aIN==1&bIN==0) {# CONDNION for only regressor a being included in the model
      PMP_uniform_a_Nb[j,1]=PMP_uniform_a_Nb[j,1]+PMP_uniform[i,1]
      PMP_random_a_Nb[j,1]=PMP_random_a_Nb[j,1]+PMP_random[i,1]
      R2_uniform_a_Nb[j,1]=R2_uniform_a_Nb[j,1]+PMP_R2_uniform[i,1]
      R2_random_a_Nb[j,1]=R2_random_a_Nb[j,1]+PMP_R2_random[i,1]
    }# the end of the CONDNION for only regressor a being included in the model
    else if(aIN==0&bIN==0) {# CONDNION for both variables being excluded from the model
      PMP_uniform_Na_Nb[j,1]=PMP_uniform_Na_Nb[j,1]+PMP_uniform[i,1]
      PMP_random_Na_Nb[j,1]=PMP_random_Na_Nb[j,1]+PMP_random[i,1]
      R2_uniform_Na_Nb[j,1]=R2_uniform_Na_Nb[j,1]+PMP_R2_uniform[i,1]
      R2_random_Na_Nb[j,1]=R2_random_Na_Nb[j,1]+PMP_R2_random[i,1]
    }# the end of the CONDNION for both variables being excluded from the model    
  }# the end of the LOOP that goes trough all the models
}# the end of the LOOP that goes trough all the regressors pairs

# MEASURES OF JOINTNESS
PMP_uniform_HCGHM_m<-zeros(c,1) # Hofmarcher et al. (2018)
PMP_uniform_LS_m<-zeros(c,1) # Ley & Steel (2007)
PMP_uniform_DW_m<-zeros(c,1) # Doppelhofer & Weeks (2009)
PMP_random_HCGHM_m<-zeros(c,1) # Hofmarcher et al. (2018)
PMP_random_LS_m<-zeros(c,1) # Ley & Steel (2007)
PMP_random_DW_m<-zeros(c,1) # Doppelhofer & Weeks (2009)
R2_uniform_HCGHM_m<-zeros(c,1) # Hofmarcher et al. (2018)
R2_uniform_LS_m<-zeros(c,1) # Ley & Steel (2007)
R2_uniform_DW_m<-zeros(c,1) # Doppelhofer & Weeks (2009)
R2_random_HCGHM_m<-zeros(c,1) # Hofmarcher et al. (2018)
R2_random_LS_m<-zeros(c,1) # Ley & Steel (2007)
R2_random_DW_m<-zeros(c,1) # Doppelhofer & Weeks (2009)

# this LOOP calculates 3 measures of jointness
for (j in 1:c){# this LOOP goes trough all the regressors pairs
  if (identical(measure,"HCGHM")){ # CONDITION for the HCGHM measure
    PMP_uniform_HCGHM_m[j,1]=((PMP_uniform_a_b[j,1]+rho)*(PMP_uniform_Na_Nb[j,1]+rho)-(PMP_uniform_Na_b[j,1]+rho)*(PMP_uniform_a_Nb[j,1]+rho))/((PMP_uniform_a_b[j,1]+rho)*(PMP_uniform_Na_Nb[j,1]+rho)+(PMP_uniform_Na_b[j,1]+rho)*(PMP_uniform_a_Nb[j,1]+rho)-rho)
    PMP_random_HCGHM_m[j,1]=((PMP_random_a_b[j,1]+rho)*(PMP_random_Na_Nb[j,1]+rho)-(PMP_random_Na_b[j,1]+rho)*(PMP_random_a_Nb[j,1]+rho))/((PMP_random_a_b[j,1]+rho)*(PMP_random_Na_Nb[j,1]+rho)+(PMP_random_Na_b[j,1]+rho)*(PMP_random_a_Nb[j,1]+rho)-rho)
    R2_uniform_HCGHM_m[j,1]=((R2_uniform_a_b[j,1]+rho)*(R2_uniform_Na_Nb[j,1]+rho)-(R2_uniform_Na_b[j,1]+rho)*(R2_uniform_a_Nb[j,1]+rho))/((R2_uniform_a_b[j,1]+rho)*(R2_uniform_Na_Nb[j,1]+rho)+(R2_uniform_Na_b[j,1]+rho)*(R2_uniform_a_Nb[j,1]+rho)-rho)
    R2_random_HCGHM_m[j,1]=((R2_random_a_b[j,1]+rho)*(R2_random_Na_Nb[j,1]+rho)-(R2_random_Na_b[j,1]+rho)*(R2_random_a_Nb[j,1]+rho))/((R2_random_a_b[j,1]+rho)*(R2_random_Na_Nb[j,1]+rho)+(R2_random_Na_b[j,1]+rho)*(R2_random_a_Nb[j,1]+rho)-rho)
  } # the end of the CONDITION for the HCGHM measure
  else if (identical(measure,"LS")){ # CONDITION for the LS measure
    PMP_uniform_LS_m[j,1]=PMP_uniform_a_b[j,1]/(PMP_uniform_Na_b[j,1]+PMP_uniform_a_Nb[j,1])
    PMP_random_LS_m[j,1]=PMP_random_a_b[j,1]/(PMP_random_Na_b[j,1]+PMP_random_a_Nb[j,1])
    R2_uniform_LS_m[j,1]=R2_uniform_a_b[j,1]/(R2_uniform_Na_b[j,1]+R2_uniform_a_Nb[j,1])
    R2_random_LS_m[j,1]=R2_random_a_b[j,1]/(R2_random_Na_b[j,1]+R2_random_a_Nb[j,1])
  } # the end of the CONDITION for the LS measure
  else if (identical(measure,"DW")){ # CONDITION for the DW measure
    PMP_uniform_DW_m[j,1]=log((PMP_uniform_a_b[j,1]/PMP_uniform_Na_b[j,1])*(PMP_uniform_Na_Nb[j,1]/PMP_uniform_a_Nb[j,1]))
    PMP_random_DW_m[j,1]=log((PMP_random_a_b[j,1]/PMP_random_Na_b[j,1])*(PMP_random_Na_Nb[j,1]/PMP_random_a_Nb[j,1]))
    R2_uniform_DW_m[j,1]=log((R2_uniform_a_b[j,1]/R2_uniform_Na_b[j,1])*(R2_uniform_Na_Nb[j,1]/PMP_uniform_a_Nb[j,1]))
    R2_random_DW_m[j,1]=log((R2_random_a_b[j,1]/R2_random_Na_b[j,1])*(R2_random_Na_Nb[j,1]/PMP_random_a_Nb[j,1]))
  } # the end of the CONDITION for the DW measure
}# the end of the LOOP that goes trough all the regressors pairs

# here we check which combination of model priors and jointness measures has the user chosen
if (measure=="HCGHM"&above==1){first<-PMP_uniform_HCGHM_m} #above the diagonal
if (measure=="HCGHM"&above==2){first<-PMP_random_HCGHM_m} #above the diagonal
if (measure=="HCGHM"&above==3){first<-R2_uniform_HCGHM_m} #above the diagonal
if (measure=="HCGHM"&above==4){first<-R2_random_HCGHM_m} #above the diagonal
if (measure=="LS"&above==1){first<-PMP_uniform_LS_m} #above the diagonal
if (measure=="LS"&above==2){first<-PMP_random_LS_m} #above the diagonal
if (measure=="LS"&above==3){first<-R2_uniform_LS_m} #above the diagonal
if (measure=="LS"&above==4){first<-R2_random_LS_m} #above the diagonal
if (measure=="DW"&above==1){first<-PMP_uniform_DW_m} #above the diagonal
if (measure=="DW"&above==2){first<-PMP_random_DW_m} #above the diagonal
if (measure=="DW"&above==3){first<-R2_uniform_DW_m} #above the diagonal
if (measure=="DW"&above==4){first<-R2_random_DW_m} #above the diagonal
if (measure=="PPI"&above==1){first<-PMP_uniform_a_b} #above the diagonal
if (measure=="PPI"&above==2){first<-PMP_random_a_b} #above the diagonal
if (measure=="PPI"&above==3){first<-R2_uniform_a_b} #above the diagonal
if (measure=="PPI"&above==4){first<-R2_rando_a_b} #above the diagonal
if (measure=="HCGHM"&below==1){second<-PMP_uniform_HCGHM_m} #below the diagonal
if (measure=="HCGHM"&below==2){second<-PMP_random_HCGHM_m} #below the diagonal
if (measure=="HCGHM"&below==3){second<-R2_uniform_HCGHM_m} #below the diagonal
if (measure=="HCGHM"&below==4){second<-R2_random_HCGHM_m} #below the diagonal
if (measure=="LS"&below==1){second<-PMP_uniform_LS_m} #below the diagonal
if (measure=="LS"&below==2){second<-PMP_random_LS_m} #below the diagonal
if (measure=="LS"&below==3){second<-R2_uniform_LS_m} #below the diagonal
if (measure=="LS"&below==4){second<-R2_random_LS_m} #below the diagonal
if (measure=="DW"&below==1){second<-PMP_uniform_DW_m} #below the diagonal
if (measure=="DW"&below==2){second<-PMP_random_DW_m} #below the diagonal
if (measure=="DW"&below==3){second<-R2_uniform_DW_m} #below the diagonal
if (measure=="DW"&below==4){second<-R2_random_DW_m} #below the diagonal
if (measure=="PPI"&below==1){second<-PMP_uniform_a_b} #below the diagonal
if (measure=="PPI"&below==2){second<-PMP_random_a_b} #below the diagonal
if (measure=="PPI"&below==3){second<-R2_uniform_a_b} #below the diagonal
if (measure=="PPI"&below==4){second<-R2_rando_a_b} #below the diagonal

# We prepare a table to store the jointness results
JointnessTable<-zeros(K,K)

for (j in 1:c){# this LOOP goes trough all the regressors pairs
  # ABOVE THE DIAGONAL
  JointnessTable[Pairs[1,j],Pairs[2,j]]=first[j,1]
  # BELOW THE DIAGONAL
  JointnessTable[Pairs[2,j],Pairs[1,j]]=second[j,1]
}# the end of the LOOP that goes trough all the regressors pairs

# Placing X in the diagonal of the jointness table
for (i in 1:K){ # LOOP that goes through the regressors no. 1
  for (j in 1:K){ # LOOP that goes through the regressors no. 2
    if (i==j){JointnessTable[i,j]="X"} # we place "X" on the main diagonal
  }# the end of the LOOP that goes through the regressors no. 2
} # the end of the LOOP that goes through the regressors no. 1

rownames(JointnessTable)<-x_names # we add regressor names to the rows
colnames(JointnessTable)<-x_names # we add regressor names to the columns

# creation of the Jointness object (joint object)
out<-list(JointnessTable)
}# THE END of the Jointness FUNCTION 

###############################################################################################################

Joint<-Jointness(Post,above=1,below=2,measure="HCGHM",rho=0.5)
JointTable<-Joint[[1]]

# 6. BestModels - function that finds and estiomates best models according to user chosen criteria:
#                PMP_unifrom/PMP_random/R2_unifrom/R2_random - end variations within them

#PARAMETERS TO BE DEFINED BE THE USER
criterion=1 # what should be a basis of ordering:
              # 1 - PMP, uniform
              # 2 - PMP, random
              # 3 - R2, uniform
              # 4 - R2, random
best=5 # number of best models to be considered
estimate=TRUE # parameter asking if we should include model estimates
           # TRUE - we add estimation to the results
           # FALSE - we skip estimation results

#THE BestModels FUNCTION STARTS HERE
BestModels=function(Post,criterion=1,best=10,estimate=T){

#Extraction of the information form the Posterior object (Post object)
x_names<-Post[[6]]# we extract vector with names of the regressors from Posterior object (Post object)
x_names<-x_names[-1]
M<-Post[[7]]# we extract M - maximum number of regressors in a model from Posterior object (Post object)
K<-Post[[8]]# we extract K - total number of regressors from Posterior object (Post object)
MS<-Post[[9]]# we extract MS - size of the mode space from Posterior object (Post object)
const<-Post[[10]]# we extract const parameter from Posterior object (Post object)
forBestModels<-Post[[13]] # Extraction of all the information

# Division of the informtion into the relevant categories
Reg_ID<-as.matrix(forBestModels[,1:M]) # extraction of the model IDs
PMPs<-as.matrix(forBestModels[,(M+1):(M+4)]) # extraction of the model weights: PMP_uniform,PMP_random,R2_uniform,R2_random
betas<-as.matrix(forBestModels[,(M+5):(M+K+5)]) # extraction of the parameter coefficients
stds<-as.matrix((forBestModels[,(M+K+6):(M+2*K+6)])^0.5) # extraction of the parameter variances and turning them into standard errors
R2<-as.matrix(forBestModels[,M+2*K+7]) # extraction of R^2

if (best>MS){# CONDITION about what to do if the user set best that is higher than MS
  # we tell the user that we are setting best=MS 
  message("best>MS - number of best models cannot be bigger than the total number of models. We set best=MS and continiue :)")
  best=MS # we set best=MS
}# end of the CONDITION about what to do if the user set best that is higher than MS

# check for the criterion chosen by the user
if (criterion==1){ranking<-PMPs[,1]} # PMP,uniform
if (criterion==2){ranking<-PMPs[,2]} # PMP,random
if (criterion==3){ranking<-PMPs[,3]} # R2,uniform
if (criterion==4){ranking<-PMPs[,4]} # R2,random

Ranking<-cbind(ranking,Reg_ID) # we add ranking criterion based on the users choice

# here we order the models according to PMP criterion
Ranking<-Ranking[order(Ranking[,1],decreasing=T),] # ordering of the models
Best_models<-Ranking[1:best,2:(M+1)] # extraction of the matrix with ordered model IDs
Ranks<-round(Ranking[1:best,1], digits = 3) # PMPs (R^2s) of the fiest 'best' models

forRegressors<-zeros(best,K) # matrix for regressors indices

for (i in 1:best){# at this LOOP we go through all the models up to "best" chosen by the user 
  for (k in 1:M){# at this LOOP we go through all regressors ID
      if (Best_models[i,k]!=0){# CONDITION for the presence of a regressor in a model
        forRegressors[i,Best_models[i,k]]=1 #
      }# the end of the CONDITION for the presence of a regressor in a model
  }# the end of the LOOP at which we go through all regressors ID 
}# the end of the LOOP at which we go through all the models up to "best" chosen by the user 

final_names<-c(x_names,"PMP")
Rank<-zeros(best,1)
for (i in 1:best){# at this LOOP we go through all the models up to "best" chosen by the user
  Rank[i,1]=paste0("'No. ",i,"'")
}# the end of the LOOP at which we go through all the models up to "best" chosen by the user 
table<-cbind(forRegressors,Ranks) # we add a column with PMPs to the models
colnames(table)<-final_names # we add regressor names to the rows
rownames(table)<-Rank # we add regressor names to the columns
bestTable<-t(table) # we make transposition to make the table look better
BestTable<-grid.table(bestTable) # here is a more fancy version of the table (works with not many models)

# creation of the BestModels object (BS object)
out<-list(bestTable,BestTable)

}# THE END OF THE BestModels FUNCTION 

bestModels<-BestModels(Post,criterion=1,best=10,estimate=T)

modelTable<-bestModels[[1]]
modelTable

ModelTable<-bestModels[[2]] # this one does not work!!!!

# HERE WE ARE INTRODUCING A CONDITION FOR ESTIMATION OF THE FIRST best MODELS
if (estimate==T){

}

###############################################################################################################

# 8. ModelSizes - function that finds prior and posterior model probability from the perspective
#                 of the model size and makes appropriate graphs for differen types of criteria:
#                PMP_unifrom/PMP_random/R2_unifrom/R2_random










###############################################################################################################

# SOME SOLUTIONS FOR THE FUTURE

#sorted_matrix <- my_matrix[order(my_matrix[, 1], decreasing=TRUE), ] #sorting elements of a matrix

#eval(parse(text = paste0("x_", t))) # this turns string into a variable

#paste0("x_", t)

#assign(paste0("m_", i),models) # creates variable m_i with with the above
#assign(paste0("m_", i,"_",j),z)

#for (k in 1:K){assign(paste0("x_", k),x[,k])} #we assgin regressors into separate column vectors an name them x_k
