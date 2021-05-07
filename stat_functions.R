lik <- function(t0in) {
  t0=t0in
  B0=diag(t+(t-1)*regressors_n)
  C0=zeros(t,cur_variables_n)
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
  if (cur_regressors_n==0) {
    C0[1,1]=t0[1]+t0[2]
    for (i3 in 2:t) {
      C0[i3,1]=t0[2]
    }
  }
  else {
    C0[1,1]=t0[1]+t0[1+cur_variables_n]
    C0[1,(2:ncol(C0))]=t(t0[2:cur_variables_n])+t(t0[(cur_variables_n+2):(cur_variables_n+1+cur_regressors_n)])
    for (i4 in 2:t) {
      C0[i4,]=t(t0[(cur_variables_n+1):(cur_variables_n+cur_regressors_n+1)])
    }
  }
  # C2 matrix has closed-form solutions 
  
  B110=B0[(1:t),(1:t)]
  B120=B0[(1:t),((t+1):ncol(B0))]
  
  o110=zeros(t,t)
  for (i5 in 1:t) {
    o110[i5,i5]=t0[2*cur_variables_n+(i5+1)]^2
  }
  o110=o110+(t0[2*cur_variables_n+1]^2)*(ones(t,t))  # Sigma 11 
  
  # Here, I split sigma 12 in the sum of two parts, phi's matrix and psi's upper triangular matrix 
  
  # phi's matrix 
  o120=zeros(t,(t-1)*regressors_n)
  for (i6 in 1:t) {
    for (i7 in 1:(t-1)) {    
      o120[i6,(1+(i7-1)*regressors_n):(i7*regressors_n)]=t(t0[(2*cur_variables_n+t+2+(i7-1)*regressors_n):(2*cur_variables_n+t+2+i7*regressors_n-1)])
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
      o121[i8,((i8-1)*regressors_n+1):(ncol(o121))]=t(t0[(2*cur_variables_n+t+2+(t-1+seq[i8])*regressors_n):(2*cur_variables_n+t+2+(t-1+seq[i8+1])*regressors_n-1)])
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