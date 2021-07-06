
@----------------------------------------------------------------------@
@                                                                      @
@                 LIKELIHOOD PROCEDURE 2. BAMLE.                       @
@                                                                      @
@----------------------------------------------------------------------@

proc(1) = likx(b0);
local iter,lik_vec,xoi,zoi,yoi,y_1oi,ymi,xmi,zmi,xzmi,b2,y_1mi,sv,a1,b1,l1,l2,su,f1,y0i,val;
lik_vec=zeros(nobs/t,1);

iter=1;
do until iter>(nobs/t);
    yoi=yo[(iter-1)*(t-1)+1:iter*(t-1),.];
    if kx > 0; xoi=xloopo[(iter-1)*(t-1)+1:iter*(t-1),.]; endif;
    y_1oi=y_1o[(iter-1)*(t-1)+1:iter*(t-1),.];  
    ymi=ym[(iter-1)*t+1,.];
    if kx > 0; xmi=xloopm[(iter-1)*t+1,.]; xzmi=xmi; endif;
    y_1mi=y_1m[(iter-1)*t+1,.];
    y0i=y0[iter,.];
    sv=b0[1,1]^2;
    a1=b0[2,1];
    f1=b0[3,1];
    su=b0[4,1]^2;
    if kx > 0; 
        b1=b0[5:klx+4,1];
        b2=b0[klx+5:klx+4+klx,1];
        lik_vec[iter,1]=-0.5*(t-1)*ln(sv)-(0.5/sv)*(yoi-a1*y_1oi-xoi*b1)'(yoi-a1*y_1oi-xoi*b1)
        -0.5*ln(su)-(0.5/su)*(ymi-a1*y_1mi-xzmi*b2-f1*y0i)^2;
    elseif kx == 0;
        lik_vec[iter,1]=-0.5*(t-1)*ln(sv)-(0.5/sv)*(yoi-a1*y_1oi)'(yoi-a1*y_1oi)
        -0.5*ln(su)-(0.5/su)*(ymi-a1*y_1mi-f1*y0i)^2;
    endif;
iter=iter+1;
endo;
val= sumc(lik_vec);
retp(-val);
endp;

@----------------------------------------------------------------------@
@                                                                      @
@                   LIKELIHOOD PROCEDURE 2  (GRADIENT)                 @
@                                                                      @
@           this procedure construct the concentrated likelihood       @
@           individual by individual and it gives the NX1 vector       @
@           of individual log-likelihoods for the computation of       @
@           the gradient and then the sandwich formula.                @
@----------------------------------------------------------------------@

proc(1) = grax(b0);
local iter,lik_vec,xoi,zoi,yoi,y_1oi,ymi,xmi,zmi,xzmi,b2,y_1mi,sv,a1,b1,l1,l2,su,f1,y0i;
lik_vec=zeros(nobs/t,1);
iter=1;
do until iter>(nobs/t);
yoi=yo[(iter-1)*(t-1)+1:iter*(t-1),.];
if kx > 0; xoi=xloopo[(iter-1)*(t-1)+1:iter*(t-1),.]; endif;
y_1oi=y_1o[(iter-1)*(t-1)+1:iter*(t-1),.];

ymi=ym[(iter-1)*t+1,.];
if kx > 0; xmi=xloopm[(iter-1)*t+1,.]; xzmi=xmi; endif;
y_1mi=y_1m[(iter-1)*t+1,.];
y0i=y0[iter,.];

sv=b0[1,1]^2;
a1=b0[2,1];
f1=b0[3,1];
su=b0[4,1]^2;
if kx > 0 ;
b1=b0[5:klx+4,1];
b2=b0[klx+5:klx+4+klx,1];
lik_vec[iter,1]=-0.5*(t-1)*ln(sv)-(0.5/sv)*(yoi-a1*y_1oi-xoi*b1)'(yoi-a1*y_1oi-xoi*b1)
-0.5*ln(su)-(0.5/su)*(ymi-a1*y_1mi-xzmi*b2-f1*y0i)^2;
elseif kx == 0;
lik_vec[iter,1]=-0.5*(t-1)*ln(sv)-(0.5/sv)*(yoi-a1*y_1oi)'(yoi-a1*y_1oi)
-0.5*ln(su)-(0.5/su)*(ymi-a1*y_1mi-f1*y0i)^2;
endif;

iter=iter+1;
endo;
retp(-lik_vec);
endp;