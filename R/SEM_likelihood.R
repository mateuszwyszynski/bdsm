#' Coefficients matrix for SEM representation
#'
#' Create coefficients matrix for Simultaneous Equations Model (SEM)
#' representation.
#'
#' @param alpha numeric
#' @param periods_n integer
#' @param beta numeric vector. Default is NULL for no regressors case.
#'
#' @return matrix
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' SEM_B_matrix(3, 4, 4:6)
SEM_B_matrix <- function(alpha, periods_n, beta = NULL) {
  alpha_matrix <- diag(rep(-alpha, periods_n-1))
  B <- diag(periods_n)
  B[2:periods_n, 1:(periods_n - 1)] <-
    B[2:periods_n, 1:(periods_n - 1)] + alpha_matrix

  if (!is.null(beta)) {
    regressors_n <- length(beta)
    beta <- beta %>% matrix(1)
    beta_matrix <- Matrix::bdiag(rep(list(-beta), periods_n - 1))

    B <- Matrix::bdiag(B, diag(regressors_n*(periods_n - 1)))
    B[2:periods_n, -1:-periods_n] <-
      B[2:periods_n, -1:-periods_n] + beta_matrix
  }
  B
}

#' Coefficients matrix for initial conditions
#'
#' Create matrix for Simultaneous Equations Model (SEM)
#' representation with coefficients placed next to initial values
#' of regressors, dependent variable and country-specific time-invariant
#' variables.
#'
#' @param alpha numeric
#' @param phi_0 numeric
#' @param periods_n numeric
#' @param beta numeric vector. Default is NULL for no regressors case.
#' @param phi_1 numeric vector. Default is NULL for no regressors case.
#'
#' @return matrix
#' @export
#'
#' @examples
#' alpha <- 9
#' phi_0 <- 19
#' beta <- 11:15
#' phi_1 <- 21:25
#' periods_n <- 4
#' SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
SEM_C_matrix <- function(alpha, phi_0,  periods_n, beta = NULL, phi_1 = NULL) {
  C1 <- matrix(rep(phi_0, periods_n))
  C1[1, 1] <- C1[1, 1] + alpha
  if (!is.null(beta)) {
    col2 <- matrix(rep(phi_1, periods_n), periods_n, byrow = TRUE)
    col2[1, 1:length(beta)] <-
      col2[1, 1:length(beta)] + beta
    C1 <- cbind(C1, col2)
  }
  C1
}

#' Covariance matrix for SEM representation
#'
#' Create covariance matrix for Simultaneous Equations Model (SEM)
#' representation. Only the part necessary to compute concentrated likelihood
#' function is computed (cf. Appendix in the Moral-Benito paper)
#'
#' @param err_var numeric
#' @param dep_vars numeric vector
#' @param phis numeric vector
#' @param psis list of numeric vectors
#'
#' @return matrix
#' @export
#'
#' @examples
#' err_var <- 1
#' dep_vars <- c(2, 2, 2, 2)
#' phis <- c(10, 10, 20, 20, 30, 30)
#' psis <- list(c(101, 102), c(103, 104, 105, 106), c(107, 108, 109, 110, 111, 112))
#' SEM_omega_matrix(err_var, dep_vars, phis, psis)
SEM_omega_matrix <- function(err_var, dep_vars, phis = NULL, psis = NULL) {
  periods_n <- length(dep_vars)

  O <- err_var*optimbase::ones(periods_n, periods_n) +
    diag(dep_vars)

  if (!is.null(phis)) {
    regressors_n <- length(phis)/(periods_n - 1)

    phi_matrix <- matrix(rep(phis, periods_n), nrow = periods_n, byrow = TRUE)

    time_fixed_psi_matrix <- function(psi, regressors_n) {
      nrows <- length(psi)/regressors_n
      t(matrix(psi, nrow = nrows, ncol = regressors_n))
    }
    psi_matrix <- psis %>%
      sapply(time_fixed_psi_matrix, regressors_n = regressors_n) %>%
      plyr::rbind.fill.matrix() %>% t() %>% tidyr::replace_na(0) %>%
      rbind(rep(0, (periods_n - 1)*regressors_n))

    O12 <- phi_matrix + psi_matrix
    O <- cbind(O, O12)
  }
  O
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
