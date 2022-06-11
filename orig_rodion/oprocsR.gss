

@----------------------------------------------------------------------@
@                   DATA TRANSFORMATION PROCEDURE                      @
@                   -----------------------------                      @
@                this procedure transforms the data set:               @
@                first: data standarization                            @
@                second: cross-sectional de-mean                       @
@                third: organization of data for the LIML estimation   @
@----------------------------------------------------------------------@

proc(2)=transf1(rawdata);
local xdata, xmean, xstd, xdatas,
      dmtx, csddata, row_inds, limldata;

@ we first standarize the Xs regressors to facilitate convergence @
xdata = rawdata[.,3:cols(rawdata)];
xmean = ((meanc(xdata))').*.(ones(rows(xdata),1));
xstd = ((stdc(xdata))').*.(ones(rows(xdata),1));
xdatas = ( xdata - xmean ) ./ xstd;
rawdata[.,3:cols(rawdata)] = xdatas;

@ dmtx will be the matrix for cross-sectional de-meaning          @
dmtx=eye(n*t)-ones(n,n).*.((1/n)*eye(t));

@ csddata will be cross-sectional de-meaned data (time dummies)   @
csddata=dmtx*rawdata;


@ now I organize the data for the LIML parametrization @
row_inds = seqa(1, t, n);
limldata = csddata[row_inds, 2];

for period(1, t, 1);
    row_inds = seqa(period, t, n);
    limldata = limldata~csddata[row_inds, 1 3:cols(csddata)];
endfor;

retp (csddata,limldata);
endp;



@-----------------------------------------------------------------@
@                  MODEL SELECTION PROCEDURE                      @
@                  ---------------------------                    @
@               this procedure takes as input a number            @
@               and it converts the number to its binary          @
@               representation in base ktotx                      @
@               in order to be a "model".                         @
@-----------------------------------------------------------------@

proc(1) = msel(turu);
local x,z,v,i;
v=zeros(ktotx,1);
x=2^(ktotx-1);
z=turu;
i=1;
do while i<=ktotx;
    if z>x;
        v[i]=1;
        z=z-x;
    else;
        v[i]=0;
    endif;
x=x/2;
i=i+1;
endo;
retp(v);
endp;


@-----------------------------------------------------------------@
@                        HESSIAN PROCEDURE                        @
@                        -----------------                        @
@               this procedure takes as input a number            @
@               computes the kxk (2-sided) hessian matrix         @
@-----------------------------------------------------------------@

proc(1) = myhess(&f,x0);
local k, hessi, h, jc, jr, x1, x2, x3, x4;
local f:proc;

k=rows(x0);
hessi=zeros(k,k);
h=1e-3;

jc=1;
do while jc<=k;
    jr=1;
    do while jr<=jc;
        x1=x0;
        x2=x0;
        x3=x0;
        x4=x0;
        x1[jr]=x0[jr]+h;
        x1[jc]=x1[jc]+h;
        x2[jr]=x0[jr]+h;
        x2[jc]=x2[jc]-h;
        x3[jr]=x0[jr]-h;
        x3[jc]=x3[jc]+h;
        x4[jr]=x0[jr]-h;
        x4[jc]=x4[jc]-h;
        hessi[jr,jc]=(f(x1)-f(x2)-f(x3)+f(x4))/(4*h^2);
        hessi[jc,jr]=hessi[jr,jc];
    jr=jr+1;
    endo;
jc=jc+1;
endo;

retp(hessi);
endp;
