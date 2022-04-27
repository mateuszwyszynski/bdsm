/*
this program runs BALIMLE approach
It estimates all possible models.
*/

new; cls; rndseed 23;
library optmum;
#include optmum.ext;
begin=date;
prandom = 0;          @* prandom=1 for theta random Ley&Steel09. prandom = 0 for theta fixed *@
dilution = 0;
dil_power = 1/2;

@---------------------------------------------------------------------------------@
@		   	                     LOADING THE DATASET			        		  @
@---------------------------------------------------------------------------------@

row=292; column=11; t=4; n=73;
load rawdata[row,column]=data.csv;

@    VARIABLES IN RAWDATA                                                                 @
@    1.yt  2.yt1  3.ish  4.sed  5.pgrw  6.pop  7.ipr  8.opem  9.gsh  10.lnlex  11.polity  @
@-----------------------------------------------------------------------------------------@

rawdata=rawdata[.,1:6];   @ I select the regressors of interest @

let varlist[5,1] = "GDP" "ISH" "SED" "PGRW" "POP";

ktotx=cols(rawdata)-2; ktoty=ktotx+1;
{data,R}=transf1(rawdata);    @ this procedure transforms the data set:              @
                              @ first: standarization of variables                   @
                              @ second: cross-sectional de-mean (time dummies)       @
                              @ third: organization of data for the LIML estimation  @


@ dependent variable for the five periods NX4 (matrix) @
Y1=R[.,2 2+ktoty 2+2*ktoty 2+3*ktoty];

@ predetermined variables for the three periods @
Y2=R[.,4+ktotx:3+2*ktotx 5+2*ktotx:4+3*ktotx 6+3*ktotx:5+4*ktotx];
X0=R[.,3:2+ktotx];

corr_matrix=Corrx(data[.,3:2+ktotx]);

Z=R[.,1]~X0;

@---------------------------------------------------------------------------------@
@		               SOME PRELIMINAR OBJECTS BALIMLE APPROACH   			      @
@---------------------------------------------------------------------------------@

@***       PRIOR STRUCTURES for THE PRIOR MODEL SIZE                           ***@
@---------------------------------------------------------------------------------@

pmsize=ktotx/2;         @  prior expected model size, options:                    @
                     	@  1. for "SDM" priors pmsize=3	                          @
                     	@  2. for "FLS" priors pmsize=ktotx/2	                  @
pinc=pmsize/ktotx;
b=(ktotx-pmsize)/pmsize;      @ parameter for beta (random) distribution of the prior inclusion probability  @

@***                     STORAGE OBJECTS                                       ***@
@---------------------------------------------------------------------------------@

mod=zeros(ktoty,1); bet=zeros(ktoty,1); 
pvarh=zeros(ktoty,1); pvarr=zeros(ktoty,1);
fy=zeros(ktoty,1); fyt=0; ppmsize=0; cout=0;
err_var = 0;
dep_vars = zeros(t, 1);
phi_1 = zeros(ktoty, 1); @ includes phi_0 @
phis = zeros(ktotx * (t - 1), 1);
psis = zeros(ktotx * t * (t - 1) / 2, 1);

@---------------------------------------------------------------------------------@
@		               LOOP COVERING FULL MODEL SPACE           			      @
@---------------------------------------------------------------------------------@

tot=2^ktotx;    @ total number of models to be estimated @
turu=1;
do while turu<=tot;
    
    mt=msel(turu);
    out = mt.== 0;          @ regressors out of the current model         @
    kx=sumc(mt); ky=kx+1;   @ number of regressors in the current model   @
    
    @ Z includes y0 and x0 as strictly exogenous variables @
    if kx==0; cur_Z=R[.,1];
    else; X0j = (delif (X0',out))'; cur_Z=R[.,1]~X0j;
    endif;
    
    @ Q is the model specific annihilator matrix           @
    Q=eye(n)-Z*inv(Z'Z)*Z'; 
    
    @ nptbe is the Number of Parameters To Be Estimated in the current model  @
    psis=0; inp=2; do while inp <=t; psis = psis + (inp-1)*ktotx; inp=inp+1; endo;
    nptbe = t + 1 + 1 + kx + 1 + kx + ktotx*(t-1) + psis;  
    
    @ t0in is the vector of initial values for the likelihood optimization   @
    t0in=0.5*ones(nptbe,1);                                  
    
    @ we now call the optimization procedure for maximizing the likelihood function @
    optset; cout1=0; opditer=100;
    _opstmth="steep one"; _opmdmth="bfgs brent";    
    {theta,fout,grad,ret}=optmum(&lik,t0in);
    // If you want to manually check if the found MLE constructs positive definite covariance matrix uncomment two lines below:
    //print sigmaConstraint(theta);
    //wait;
    
    @ we now compute model-specific standard errors @
    he=myhess(&lik,theta);
    Gmat=gradp(&likgra,theta);
    Imat=Gmat'Gmat;
    stdr=sqrt(diag(inv(he)*(Imat)*inv(he)));
    stdh=sqrt(diag((inv(he))));
    varr=stdr.^2; varh=stdh.^2;
    
    @ storing results for the CURRENT model @    
    logl=(-fout-(ky/2)*(ln(n*t)))/n;
    bict=exp(logl);                              @ integrated likelihood approximation           @
    
    @ prior model probability (either random -Ley&Steel09- or fixed) @
    if prandom == 1;
        priorprobt=(gamma(1+kx))*(gamma(b+ktotx-kx));    @theta random@
    elseif prandom == 0;
        priorprobt=(pinc)^(kx)*(1-pinc)^(ktotx-kx);      @theta fixed@ 
    endif;
    
    if dilution == 1;
        x={};
        for i(1,ktotx,1);
            if mt[i]==1;
                x=x~i;
            endif;
        endfor;
        if cols(x)<=1;
            corr_det=1;
        else;
            cur_corr_matrix=corr_matrix[x,x];
            corr_det=det(cur_corr_matrix);
        endif;
        priorprobt=corr_det^dil_power*priorprobt;
    endif;

    @ posterior model probability  @
    postprob=priorprobt*bict;
    
    @ selecting estimates of interest (i.e. alpha and betas) @
    bt=theta[1:ky]; stdrt=stdr[1:ky]; stdht=stdh[1:ky];
    varht=varh[1:ky]; varrt=varr[1:ky];
    cur_phi_1 = theta[(ky+1):2*ky];
    cur_err_var = theta[2*ky+1:2*ky+1];
    cur_dep_vars = theta[2*ky+2:2*ky+5];
    cur_phis = theta[2*ky+6:2*ky+6+3*ktotx-1];
    cur_psis = theta[2*ky+6+3*ktotx:2*ky+6+3*ktotx+ktotx*t*(t - 1)/2-1];
    
    @ constructing the full vector of estimates @
    mty=1|mt; 
    bt1=zeros(ktoty,1);
    cur_phi_1_full = zeros(ktoty, 1); @ includes phi_0 @
    stdrt1=zeros(ktoty,1); stdht1=zeros(ktoty,1);
    varht1=zeros(ktoty,1); varrt1=zeros(ktoty,1);
    it1=0;
    it=1;
    do while it<=ktoty;
        if mty[it]==1;
            it1=1+it1;
            bt1[it]=bt[it1];
            cur_phi_1_full[it] = cur_phi_1[it1];
            stdrt1[it]=stdrt[it1];
            stdht1[it]=stdht[it1];
            varht1[it]=varht[it1];
            varrt1[it]=varrt[it1];
        else;
            bt1[it]=0;         @ if the regressor is not in the model, 0 @
            stdrt1[it]=0;
            stdht1[it]=0;
            varht1[it]=0;
            varrt1[it]=0;
        endif;
    it=it+1;
    endo;
    
    @ calculating the percentage of significant regressions @
    ptr=bt1./stdht1;
    ntr=bt1./stdht1;
    if turu==1; pts=ptr; nts=ntr;
    else; pts=pts~ptr; nts=nts~ntr;
    endif;
        
    @ here I clear some global variables to save space in the memory @
    clear _eigerr,_rtl_ver,__weight,__vtype,__vpad,__tol,__sort,__rowfac,
    __row,__range,__prec,__output,__miss,__macheps,__INFp,__INFn,__INDEFp,__INDEFn,
    __fmtnv,__fmtcv,__ff,__con,__altnam;
    
    @ accumulating estimates for posterior model probabilities @
    mod=mod+mty;
    fy=fy+postprob*mty;
    fyt=fyt+postprob;
    ppmsize=ppmsize+postprob*(sumc(mty));
    
    @ storing estimates conditional on inclusion @
    bet=bet+postprob*bt1;
    pvarr=pvarr+(postprob*varrt1+postprob*(bt1.*bt1));         @ as in Leamer (1978) @
    pvarh=pvarh+(postprob*varht1+postprob*(bt1.*bt1));         @ as in Leamer (1978) @
    phi_1 = phi_1 + postprob * cur_phi_1_full;
    err_var = err_var + postprob * cur_err_var;
    dep_vars = dep_vars + postprob * cur_dep_vars;
    phis = phis + postprob * cur_phis;
    psis = psis + postprob * cur_psis;
    
    @ here we store model-specific diagnostics and estimates (BICs, likelihoods, betas...) @   
    if turu==1;
        modprob=postprob; modelid=turu; modpri=priorprobt; liks=exp(-fout/n); bics=bict;
        betas=bt1; stds=stdht1; stdsr=stdrt1; foutt=-fout;
        phi_1s = phi_1;
        err_vars = err_var;
        dep_varss = dep_vars;
        phiss = phis;
        psiss = psis;
    else;
        modprob=modprob|postprob; modelid=modelid|turu; modpri=modpri|priorprobt; liks=liks|exp(-fout/n); bics=bics|bict;
        betas=betas~bt1; stds=stds~stdht1; stdsr=stdsr~stdrt1; foutt=foutt|(-fout);
        phi_1s = phi_1s ~ phi_1;
        err_vars = err_vars ~ err_var;
        dep_varss = dep_varss ~ dep_vars;
        phiss = phiss ~ phis;
        psiss = psiss ~ psis;
    endif; 
    cls;
    
turu=turu+1;
endo;

@                         end of loop covering full model space                        @
@--------------------------------------------------------------------------------------@


popmsize=ppmsize/fyt;
modprob1=modprob/sumc(modprob);
idprob=modelid~modprob1~modpri~bics;

@ computing posterior moments CONDITIONAL on inclusion @
postprobinc=fy/fyt;
postmean=bet./fy;
varrleamer=(pvarr./fy)-postmean.^2;
varhleamer=(pvarh./fy)-postmean.^2;
poststdr=sqrt(varrleamer);
poststdh=sqrt(varhleamer);
tr=postmean./poststdr;
th=postmean./poststdh;

@ computing UNCONDITIONAL posterior moments  @
upostmean = postmean .* postprobinc;
uvarrleamer = (varrleamer + (postmean.^2)).*postprobinc - (upostmean.^2);
uvarhleamer = (varhleamer + (postmean.^2)).*postprobinc - (upostmean.^2);
upoststdr=sqrt(uvarrleamer);
upoststdh=sqrt(uvarhleamer);

@ computing percentage of significant coeff estimates @
nts=nts'; 
pts=pts';
jt=1;
do while jt<=ktoty;
    ntss=packr(nts[.,jt]); ptss=packr(pts[.,jt]); nsig=ntss.<-1.96; psig=ptss.>1.96;
        if jt==1; negper=meanc(nsig); posper=meanc(psig);
        else; negper=negper|meanc(nsig); posper=posper|meanc(psig);
        endif;
    jt=jt+1;
endo;

@ we print final results in balimle_results.out file @
result=varlist~postprobinc~postmean~poststdh~poststdr~upostmean~upoststdh~upoststdr;
the_end=date; days=the_end[3]-begin[3]; et=(the_end[4]+8640000*days-begin[4])/100;

output file = balimle_results.out reset;
print;
print "varname -- postprob -- pmean -- std -- stdR -- unc_pmean -- unc_std -- unc_stdR";
let mask[1,8] = 0 1 1 1 1 1 1 1;
let fmt[8,3] = 
"-*.*s" 8 8 
"*.*lf" 10 4 
"*.*lf" 10 4 
"*.*lf" 10 4
"*.*lf" 10 4
"*.*lf" 10 4
"*.*lf" 10 4
"*.*lf" 10 4; 
d = printfm(result,mask,fmt);
print;
print " 1.- FURTHER INFORMATION ";
print "Prior Mean Model Size="; pmsize;
print "Prior Inclusion Probability="; pinc;
print "Posterior Mean Model Size="; popmsize;
print "minutes="; et/60;
print "hours="; (et/60)/60;
print "number of loops="; turu-1;
print "minutes per loop="; (et/60)/(turu-1);
print "----------------------------------------------------------------------";
print;
print; " 2.- ALL BETAS (each row is a different model)";
print betas';
print;
print "----------------------------------------------------------------------";
print; " 3.- ALL STD. ERRORS (each row is a different model)";
print stds';
print;
print "----------------------------------------------------------------------";
print; " 4.- ALL ROBUST STD. ERRORS (each row is a different model)";
print stdsr';
print;
print "----------------------------------------------------------------------";
print;
print " 5.- MODELS INFO ";
print " model -- postprob -- priorprob -- bics";
print idprob;
print;
print "----------------------------------------------------------------------";
print;
print " 6.- MEAN PHI_1 ";
print phi_1';
print;
print "----------------------------------------------------------------------";
print;
print " 7.- MEAN ERR VARS ";
print err_var';
print;
print "----------------------------------------------------------------------";
print;
print " 8.- MEAN DEP VARS ";
print dep_vars';
print;
print "----------------------------------------------------------------------";
print;
print " 9.- MEAN PHIS ";
print phis';
print;
print "----------------------------------------------------------------------";
print;
print " 10.- MEAN PSIS ";
print psis';
print;
print "----------------------------------------------------------------------";
print;
print " 11.- ALL PHI_1S ";
print phi_1s';
print;
print "----------------------------------------------------------------------";
print;
print " 12.- ALL ERR VARS ";
print err_vars';
print;
print "----------------------------------------------------------------------";
print;
print " 13.- ALL DEP VARS ";
print dep_varss';
print;
print "----------------------------------------------------------------------";
print;
print " 14.- ALL PHIS ";
print phiss';
print;
print "----------------------------------------------------------------------";
print;
print " 15.- ALL PSIS ";
print psiss';
output off;

#include liks1.prg;
#include oprocs.prg;
#include optims.prg;
