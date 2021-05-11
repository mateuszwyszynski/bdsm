#' Coefficients matrix for SEM representation
#'
#' Create coefficients matrix for Simultaneous Equations Model (SEM)
#' representation.
#'
#' @param params numeric vector
#' @param regressors_n integer
#' @param periods_n integer
#'
#' @return matrix
#' @export
#'
#' @examples
#' SEM_B_matrix(3:6, 3, 4)
SEM_B_matrix <- function(params, regressors_n, periods_n) {
  alpha <- params[1]
  alpha_matrix <- diag(rep(-alpha, periods_n-1))
  B11 <- diag(periods_n)
  B11[2:periods_n, 1:(periods_n - 1)] <-
    B11[2:periods_n, 1:(periods_n - 1)] + alpha_matrix

  betas <- params[-1] %>% matrix(1)
  betas_matrix <- bdiag(rep(list(-betas), periods_n - 1))
  B12 <- rbind(zeros(1, regressors_n*(periods_n - 1)), betas_matrix)

  B <- bdiag(B11, diag(regressors_n*(periods_n - 1)))
  B[2:periods_n, -1:-periods_n] <-
    B[2:periods_n, -1:-periods_n] + betas_matrix
  B
}

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
}

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

  U10=t(B110%*%t(Y1)+B120%*%t(Y2)-C0%*%t(Z))  # Ui1 from the paper
  H=crossprod(Y2-U10%*%solve(o110)%*%o120,res_maker_matrix)%*%(Y2-U10%*%solve(o110)%*%o120)
  for (iter in 1:n) {
    u10i=as.matrix(U10[iter,])
    likvec[iter]=-(1/2)*log(det(o110))-(1/2)*log(det(H/n))-(1/2)*(t(u10i)%*%solve(o110)%*%u10i)
  }

  return(-likvec)
}
