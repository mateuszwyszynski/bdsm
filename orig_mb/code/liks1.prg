@----------------------------------------------------------------------@
@                                                                      @
@               LIKELIHOOD PROCEDURE 1. BALIMLE.                       @
@                                                                      @
@----------------------------------------------------------------------@

proc(1)=lik(t0in);
local likf,ii,i1,i2,B0,B110,B120,C0,U10,H,o110,o120,o210,t0,t0i,fact,S11_inv,F12,M,G22_inverse,G22, likelihood_constant;
local err_var_ind, dep_vars_last_ind, o121, seq, dseq, phis_n, phis_last_ind, start_col_ind;

t0=t0in;

B0=eye(t+(t-1)*ktotx);
C0=zeros(t,ky);
ii=2;
do while ii<=t;
    B0[ii,ii-1]=-t0[1];
    ii=ii+1;
endo;

i2=1;               @ georg correction for the proper ordering of variables @
for i1(1,ktotx,1);
    if mt[i1]==0;   @ if x variable is not included @
        B0=B0;
    else;           @ if x variable is included     @
        for i11(2,t,1);
            B0[i11,t+1+(i11-2)*ktotx+(i1-1)]=-t0[1+i2];
        endfor;
        i2=i2+1;
    endif;
endfor;

if kx==0;
    C0[1,1]=t0[1]+t0[2];
    for i3(2,t,1);
        C0[i3,1]=t0[2];
    endfor;
else;
    C0[1,1]=t0[1]+t0[1+ky];
    C0[1,2:cols(C0)]=t0[2:ky]'+t0[ky+2:ky+1+kx]';
    for i4(2,t,1);
        C0[i4,.]=t0[ky+1:ky+kx+1]';
    endfor;
endif;

B110=B0[1:t,1:t];
B120=B0[1:t,(t+1):cols(B0)];

err_var_ind=2*ky+1;

o110=zeros(t,t);
for i5(1,t,1);
    o110[i5,i5]=t0[err_var_ind+i5]^2;
endfor;
o110=o110+(t0[err_var_ind]^2)*(ones(t,t));


@ phi's matrix @
dep_vars_last_ind=2*ky+t+1;

o120=zeros(t,(t-1)*ktotx);
for i6(1,t,1);
    for i7(1,t-1,1);
        for reg_ind(1,ktotx,1);
            o120[i6,(reg_ind+(i7-1)*ktotx)]=t0[dep_vars_last_ind+(i7-1)*ktotx+reg_ind];
        endfor;
    endfor;
endfor;

@ psi's upper triangular matrix @
o121=zeros(t,(t-1)*ktotx);

@ as o121 is an upper triangular matrix, each subsequent row has @
@ 1 element less @

seq=zeros(t,1);
dseq=0;
for iseq(1,t,1);
    seq[iseq]=dseq;
    dseq=dseq+t-iseq;
endfor;

phis_n=ktotx*(t-1);
phis_last_ind=dep_vars_last_ind+phis_n;

for row_ind(1,t-1,1);
    start_col_ind=row_ind;
    for col_ind(start_col_ind,t,1);
        if col_ind==t;
            o121=o121;
        else;
            for reg_ind(1,ktotx,1);
                o121[row_ind,(start_col_ind-1)*ktotx+reg_ind+(col_ind-start_col_ind)*ktotx]=t0[phis_last_ind+(col_ind-start_col_ind)*ktotx+ktotx*((row_ind-1)*(t-1)-(row_ind-2)*(row_ind-1)/2)+reg_ind];
            endfor;
        endif;
    endfor;
endfor;
@ Sigma 12 @
o120=o120+o121;

o210=o120';

if det(o110)<=0;
     if cout1==0;
        cout=cout+1;
     endif;
        if cout1<=250;
            t0i=.5*ones(rows(t0in),1);
        elseif cout1<=500;
            t0i=t0in;
        elseif cout1<=750;
            t0i=rndu(rows(t0in),1);
        elseif cout1<=1000;
            t0i=rndu(rows(t0in),1);
        else;
            cls;
            t0i=rndn(rows(t0in),1);
        endif;
        fact=1./(10.^(ceil(log(abs(t0i)))));
        t0in=fact.*t0i;
        likf=0;
        cout1=cout1+1;
    else;
        U10=(B110*Y1'+B120*Y2'-C0*cur_Z')';
        S11_inv = inv(o110);
        F12 = - S11_inv * o120;
        M = Y2 + U10 * F12;
        H = M'Q*M;
        G22_inverse = 1/n * H;
        G22 = inv(G22_inverse);
        likelihood_constant = -(n/2) * (n + (t + (t - 1)*ktotx)*ln(2*pi));
        likf = -(1/2)*(n*ln(det(o110)/det(G22)) + sumc(diag(inv(o110)*U10'U10))) + likelihood_constant;
    endif;

retp(-likf);
endp;

proc(1) = sigmaConstraint(t0in);
local ii,i1,i2,B0,B110,B120,C0,U10,H,o110,o120,o210,t0,S11_inv,F12,L,G22_inverse,G22,S22, sigma;
local err_var_ind, dep_vars_last_ind, o121, seq, dseq, phis_n, phis_last_ind, start_col_ind;
    
t0=t0in;

B0=eye(t+(t-1)*ktotx);
C0=zeros(t,ky);
ii=2;
do while ii<=t;
    B0[ii,ii-1]=-t0[1];
    ii=ii+1;
endo;

i2=1;               @ georg correction for the proper ordering of variables @
for i1(1,ktotx,1);
    if mt[i1]==0;   @ if x variable is not included @
        B0=B0;
    else;           @ if x variable is included     @
        for i11(2,t,1);
            B0[i11,t+1+(i11-2)*ktotx+(i1-1)]=-t0[1+i2];
        endfor;
        i2=i2+1;
    endif;
endfor;

if kx==0;
    C0[1,1]=t0[1]+t0[2];
    for i3(2,t,1);
        C0[i3,1]=t0[2];
    endfor;
else;
    C0[1,1]=t0[1]+t0[1+ky];
    C0[1,2:cols(C0)]=t0[2:ky]'+t0[ky+2:ky+1+kx]';
    for i4(2,t,1);
        C0[i4,.]=t0[ky+1:ky+kx+1]';
    endfor;
endif;

B110=B0[1:t,1:t];
B120=B0[1:t,(t+1):cols(B0)];

err_var_ind=2*ky+1;

o110=zeros(t,t);
for i5(1,t,1);
    o110[i5,i5]=t0[err_var_ind+i5]^2;
endfor;
o110=o110+(t0[err_var_ind]^2)*(ones(t,t));


@ phi's matrix @
dep_vars_last_ind=2*ky+t+1;

o120=zeros(t,(t-1)*ktotx);
for i6(1,t,1);
    for i7(1,t-1,1);
        for reg_ind(1,ktotx,1);
            o120[i6,(reg_ind+(i7-1)*ktotx)]=t0[dep_vars_last_ind+(i7-1)*ktotx+reg_ind];
        endfor;
    endfor;
endfor;

@ psi's upper triangular matrix @
o121=zeros(t,(t-1)*ktotx);

@ as o121 is an upper triangular matrix, each subsequent row has @
@ 1 element less @

seq=zeros(t,1);
dseq=0;
for iseq(1,t,1);
    seq[iseq]=dseq;
    dseq=dseq+t-iseq;
endfor;

phis_n=ktotx*(t-1);
phis_last_ind=dep_vars_last_ind+phis_n;

for row_ind(1,t-1,1);
    start_col_ind=row_ind;
    for col_ind(start_col_ind,t,1);
        if col_ind==t;
            o121=o121;
        else;
            for reg_ind(1,ktotx,1);
                o121[row_ind,(start_col_ind-1)*ktotx+reg_ind+(col_ind-start_col_ind)*ktotx]=t0[phis_last_ind+(col_ind-start_col_ind)*ktotx+ktotx*((row_ind-1)*(t-1)-(row_ind-2)*(row_ind-1)/2)+reg_ind];
            endfor;
        endif;
    endfor;
endfor;
@ Sigma 12 @
o120=o120+o121;

o210=o120';

U10=(B110*Y1'+B120*Y2'-C0*cur_Z')';
S11_inv = inv(o110);
F12 = - S11_inv * o120;
L = Y2 + U10 * F12;
H = L' * Q * L;
G22_inverse = 1/n * H;
S22 = G22_inverse + F12' * o110 * F12;
G22 = inv(G22_inverse);

sigma = blockDiag(o110, S22);
sigma[1:t, (t + 1):(t + (T - 1) * ktotx)] = o120;
sigma[(t + 1):(t + (T - 1) * ktotx), 1:t] = o120';

if det(G22)/det(o110)<0;
    print "hello there";
endif;

retp(minc(eigh(sigma)));

endp;




@----------------------------------------------------------------------@
@                                                                      @
@                   LIKELIHOOD PROCEDURE 1  (GRADIENT)                 @
@                                                                      @
@           this procedure construct the concentrated likelihood       @
@           individual by individual and it gives the NX1 vector       @
@           of individual log-likelihoods for the computation of       @
@           the gradient and then the sandwich formula.                @
@----------------------------------------------------------------------@


proc(1)=likgra(t0in);
local likvec,ii,i1,i2,B0,B110,B120,C0,U10,H,o110,o120,o210,iter,u10i,t0;
local err_var_ind, dep_vars_last_ind, o121, seq, dseq, phis_n, phis_last_ind, start_col_ind;


t0=t0in;


likvec=zeros(n,1);

B0=eye(t+(t-1)*ktotx);
C0=zeros(t,ky);
ii=2;
do while ii<=t;
    B0[ii,ii-1]=-t0[1];
    ii=ii+1;
endo;


i2=1;               @ georg correction for the proper ordering of variables @
for i1(1,ktotx,1);
    if mt[i1]==0;   @ if x variable is not included @
        B0=B0;
    else;           @ if x variable is included     @
        for i11(2,t,1);
            B0[i11,t+1+(i11-2)*ktotx+(i1-1)]=-t0[1+i2];
        endfor;
        i2=i2+1;
    endif;
endfor;

if kx==0;
    C0[1,1]=t0[1]+t0[2];
    for i3(2,t,1);
        C0[i3,1]=t0[2];
    endfor;
else;
    C0[1,1]=t0[1]+t0[1+ky];
    C0[1,2:cols(C0)]=t0[2:ky]'+t0[ky+2:ky+1+kx]';
    for i4(2,t,1);
        C0[i4,.]=t0[ky+1:ky+kx+1]';
    endfor;
endif;

B110=B0[1:t,1:t];
B120=B0[1:t,(t+1):cols(B0)];

err_var_ind=2*ky+1;

o110=zeros(t,t);
for i5(1,t,1);
    o110[i5,i5]=t0[err_var_ind+i5]^2;
endfor;
o110=o110+(t0[err_var_ind]^2)*(ones(t,t));


@ phi's matrix @
dep_vars_last_ind=2*ky+t+1;

o120=zeros(t,(t-1)*ktotx);
for i6(1,t,1);
    for i7(1,t-1,1);
        for reg_ind(1,ktotx,1);
            o120[i6,(reg_ind+(i7-1)*ktotx)]=t0[dep_vars_last_ind+(i7-1)*ktotx+reg_ind];
        endfor;
    endfor;
endfor;

@ psi's upper triangular matrix @
o121=zeros(t,(t-1)*ktotx);

@ as o121 is an upper triangular matrix, each subsequent row has @
@ 1 element less @

seq=zeros(t,1);
dseq=0;
for iseq(1,t,1);
    seq[iseq]=dseq;
    dseq=dseq+t-iseq;
endfor;

phis_n=ktotx*(t-1);
phis_last_ind=dep_vars_last_ind+phis_n;

for row_ind(1,t-1,1);
    start_col_ind=row_ind;
    for col_ind(start_col_ind,t,1);
        if col_ind==t;
            o121=o121;
        else;
            for reg_ind(1,ktotx,1);
                o121[row_ind,(start_col_ind-1)*ktotx+reg_ind+(col_ind-start_col_ind)*ktotx]=t0[phis_last_ind+(col_ind-start_col_ind)*ktotx+ktotx*((row_ind-1)*(t-1)-(row_ind-2)*(row_ind-1)/2)+reg_ind];
            endfor;
        endif;
    endfor;
endfor;
@ Sigma 12 @
o120=o120+o121;

o210=o120';


U10=(B110*Y1'+B120*Y2'-C0*cur_Z')';
H=(Y2-U10*inv(o110)*o120)'Q*(Y2-U10*inv(o110)*o120);

iter=1;
do while iter <= n;
    u10i=U10[iter,.]';
    likvec[iter]=-(1/2)*ln(det(o110))-(1/2)*ln(det(H/n))
                 -(1/2)*(u10i'inv(o110)*u10i)- (1/2) * (n + (t + (t - 1)*ktotx)*ln(2*pi));
    iter=iter+1;
endo;

retp(-likvec);
endp;


