@----------------------------------------------------------------------@
@                                                                      @
@               LIKELIHOOD PROCEDURE 1. BALIMLE.                       @
@                                                                      @
@----------------------------------------------------------------------@

proc(1)=lik(t0in);
local likf,ii,i1,i2,B0,B110,B120,C0,U10,H,o110,o120,o210,t0,t0i,fact,S11_inv,L,M,G22_inverse,G22;

t0=t0in;

B0=eye(t+(t-1)*ktotx);
C0=zeros(t,ky);
ii=2;
do while ii<=t;
    B0[ii,ii-1]=-t0[1];
    ii=ii+1;
endo;

i1=1;
i2=1;               @ georg correction for the proper ordering of variables @
do while i1<=ktotx;
    if mt[i1]==0;
        B0[2,5+(i1-1):4+i1]=0;
        B0[3,5+ktotx+(i1-1):4+ktotx+i1]=0;
        B0[4,5+2*ktotx+(i1-1):4+2*ktotx+i1]=0;
    elseif mt[i1]==1;
        B0[2,5+(i1-1):4+i1]=-t0[1+i2];                  @ georg @
        B0[3,5+ktotx+(i1-1):4+ktotx+i1]=-t0[1+i2];      @ georg @
        B0[4,5+2*ktotx+(i1-1):4+2*ktotx+i1]=-t0[1+i2];  @ georg @
        i2=i2+1;                                        @ georg @
    endif;
i1=i1+1;
endo;



if kx==0;
    C0[1,1]=t0[1]+t0[2];
    C0[2,1]=t0[2];
    C0[3,1]=t0[2];
    C0[4,1]=t0[2];
else;
    C0[1,1]=t0[1]+t0[1+ky];
    C0[1,2:cols(C0)]=t0[2:ky]'+t0[ky+2:ky+1+kx]';
    C0[2,.]=t0[ky+1:ky+kx+1]';
    C0[3,.]=t0[ky+1:ky+kx+1]';
    C0[4,.]=t0[ky+1:ky+kx+1]';
endif;

B110=B0[1:4,1:4];
B120=B0[1:4,5:cols(B0)];

o110=zeros(4,4);
o110[1,1]=t0[2*ky+2]^2; o110[2,2]=t0[2*ky+3]^2; o110[3,3]=t0[2*ky+4]^2; o110[4,4]=t0[2*ky+5]^2;
o110=o110+(t0[2*ky+1]^2)*(ones(4,1)*ones(4,1)');

o120=zeros(t,3*ktotx);

o120[1,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]'+t0[2*ky+6+3*ktotx:2*ky+6+3*ktotx+(ktotx-1)]';   o120[1,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]'+t0[2*ky+6+4*ktotx:2*ky+6+4*ktotx+(ktotx-1)]';   o120[1,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+5*ktotx:2*ky+6+5*ktotx+(ktotx-1)]';
o120[2,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[2,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]'+t0[2*ky+6+6*ktotx:2*ky+6+6*ktotx+(ktotx-1)]';   o120[2,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+7*ktotx:2*ky+6+7*ktotx+(ktotx-1)]';
o120[3,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[3,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]';                                                o120[3,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+8*ktotx:2*ky+6+8*ktotx+(ktotx-1)]';
o120[4,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[4,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]';                                                o120[4,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]';                          
 
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
        L = S11_inv * o120;
        M = Y2 - U10 * L;
        H = M'Q*M;
        G22_inverse = L' * o110 * L + 1/n * (M' * Q * (Y2 + U10 * L) + L' * U10' * U10 * L);
        G22 = inv(G22_inverse);
        likf=-(n/2)*ln(det(o110))-(1/2)*sumc(diag(inv(o110)*U10'U10))+(n/2)*ln(det(G22)) - 1/2 * sumc(diag(H*G22));
    endif;

retp(-likf);
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


t0=t0in;


likvec=zeros(n,1);

B0=eye(t+(t-1)*ktotx);
C0=zeros(t,ky);
ii=2;
do while ii<=t;
    B0[ii,ii-1]=-t0[1];
    ii=ii+1;
endo;


i1=1;
i2=1;               @ georg correction for the proper ordering of variables @
do while i1<=ktotx;
    if mt[i1]==0;
        B0[2,5+(i1-1):4+i1]=0;
        B0[3,5+ktotx+(i1-1):4+ktotx+i1]=0;
        B0[4,5+2*ktotx+(i1-1):4+2*ktotx+i1]=0;
    elseif mt[i1]==1;
        B0[2,5+(i1-1):4+i1]=-t0[1+i2];                  @ georg @
        B0[3,5+ktotx+(i1-1):4+ktotx+i1]=-t0[1+i2];      @ georg @
        B0[4,5+2*ktotx+(i1-1):4+2*ktotx+i1]=-t0[1+i2];  @ georg @
        i2=i2+1;                                        @ georg @
    endif;
i1=i1+1;
endo;


if kx==0;
    C0[1,1]=t0[1]+t0[2];
    C0[2,1]=t0[2];
    C0[3,1]=t0[2];
    C0[4,1]=t0[2];
else;
    C0[1,1]=t0[1]+t0[1+ky];
    C0[1,2:cols(C0)]=t0[2:ky]'+t0[ky+2:ky+1+kx]';
    C0[2,.]=t0[ky+1:ky+kx+1]';
    C0[3,.]=t0[ky+1:ky+kx+1]';
    C0[4,.]=t0[ky+1:ky+kx+1]';
endif;

B110=B0[1:4,1:4];
B120=B0[1:4,5:cols(B0)];

o110=zeros(4,4);
o110[1,1]=t0[2*ky+2]^2; o110[2,2]=t0[2*ky+3]^2; o110[3,3]=t0[2*ky+4]^2; o110[4,4]=t0[2*ky+5]^2;
o110=o110+(t0[2*ky+1]^2)*(ones(4,1)*ones(4,1)');

o120=zeros(t,3*ktotx);

o120[1,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]'+t0[2*ky+6+3*ktotx:2*ky+6+3*ktotx+(ktotx-1)]';   o120[1,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]'+t0[2*ky+6+4*ktotx:2*ky+6+4*ktotx+(ktotx-1)]';   o120[1,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+5*ktotx:2*ky+6+5*ktotx+(ktotx-1)]';
o120[2,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[2,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]'+t0[2*ky+6+6*ktotx:2*ky+6+6*ktotx+(ktotx-1)]';   o120[2,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+7*ktotx:2*ky+6+7*ktotx+(ktotx-1)]';
o120[3,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[3,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]';                                                o120[3,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]'+t0[2*ky+6+8*ktotx:2*ky+6+8*ktotx+(ktotx-1)]';
o120[4,1:ktotx]=t0[2*ky+6:2*ky+6+(ktotx-1)]';                                                o120[4,ktotx+1:2*ktotx]=t0[2*ky+6+ktotx:2*ky+6+ktotx+(ktotx-1)]';                                                o120[4,2*ktotx+1:3*ktotx]=t0[2*ky+6+2*ktotx:2*ky+6+2*ktotx+(ktotx-1)]';                          
 
o210=o120';


U10=(B110*Y1'+B120*Y2'-C0*cur_Z')';
H=(Y2-U10*inv(o110)*o120)'Q*(Y2-U10*inv(o110)*o120);

iter=1;
do while iter <= n;
    u10i=U10[iter,.]';
    likvec[iter]=-(1/2)*ln(det(o110))-(1/2)*ln(det(H/n))
                 -(1/2)*(u10i'inv(o110)*u10i);
    iter=iter+1;
endo;

retp(-likvec);
endp;


