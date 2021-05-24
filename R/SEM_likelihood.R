#' Coefficients matrix for SEM representation
#'
#' Create coefficients matrix for Simultaneous Equations Model (SEM)
#' representation.
#'
#' @param alpha numeric
#' @param periods_n integer
#' @param beta numeric vector. Default is c() for no regressors case.
#'
#' @return List with two matrices B11 and B12
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' SEM_B_matrix(3, 4, 4:6)
SEM_B_matrix <- function(alpha, periods_n, beta = c()) {
  alpha_matrix <- diag(rep(-alpha, periods_n-1))
  B11 <- diag(periods_n)
  B11[2:periods_n, 1:(periods_n - 1)] <-
    B11[2:periods_n, 1:(periods_n - 1)] + alpha_matrix

  B12 <- if (length(beta) != 0) {
    regressors_n <- length(beta)
    beta <- beta %>% matrix(1)
    rbind(optimbase::zeros(1, regressors_n*(periods_n - 1)),
          Matrix::bdiag(rep(list(-beta), periods_n - 1)))
  } else {
    NULL
  }
  list(B11, B12)
}

orig_B_matrix <- function(alpha, periods_n, beta = c()) {
  regressors_n <- length(beta)
  B0 <- diag(periods_n+(periods_n-1)*regressors_n)
  B11 <- B0[(1:periods_n),(1:periods_n)]
  for (ii in 2:periods_n) {
    B11[ii,(ii-1)] <- -alpha
  }
  if (length(beta) != 0) {
    B12 <- B0[(1:periods_n),((periods_n+1):ncol(B0))]
    for (row_ind in 2:periods_n) {
      B12[row_ind,(1 + (row_ind-2)*regressors_n):((row_ind-1)*regressors_n)] <-
        -beta
    }
  } else {
    B12 <- NULL
  }

  list(B11, B12)
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
#' @param beta numeric vector. Default is c() for no regressors case.
#' @param phi_1 numeric vector. Default is c() for no regressors case.
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
SEM_C_matrix <- function(alpha, phi_0,  periods_n, beta = c(), phi_1 = c()) {
  C1 <- matrix(rep(phi_0, periods_n))
  C1[1, 1] <- C1[1, 1] + alpha
  if (length(beta) != 0) {
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
#' @param psis numeric vector. Psis should be passed column-wise, i.e. they will
#' be filled into Sigma12 across columns first.
#'
#' @return List with two matrices Sigma11 and Sigma12
#' @importFrom magrittr %>%
#' @export
#'
#' @examples
#' err_var <- 1
#' dep_vars <- c(2, 2, 2, 2)
#' phis <- c(10, 10, 20, 20, 30, 30)
#' psis <- c(101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
#' SEM_sigma_matrix(err_var, dep_vars, phis, psis)
SEM_sigma_matrix <- function(err_var, dep_vars, phis = c(),
                             psis = c(), psis_byrow = TRUE) {
  periods_n <- length(dep_vars)

  O11 <- err_var^2*optimbase::ones(periods_n, periods_n) +
    diag(dep_vars^2)

  O12 <- if (length(phis) != 0) {
    regressors_n <- length(phis)/(periods_n - 1)

    phi_matrix <- matrix(rep(phis, periods_n), nrow = periods_n, byrow = TRUE)

    if (psis_byrow) {
      fill_zeros <- function(v, desired_len) {
        zeros_n <- desired_len - length(v)
        c(rep(0, zeros_n), v)
      }

      psi_matrix <- psis %>%
        split(rep(1:(periods_n-1), regressors_n*((periods_n-1):1))) %>%
        lapply(fill_zeros, desired_len = regressors_n*(periods_n-1)) %>%
        unlist() %>% matrix(nrow = periods_n - 1, byrow = TRUE) %>%
        rbind(rep(0, (periods_n - 1)*regressors_n))
    } else {
      time_fixed_psi_matrix <- function(psi, regressors_n) {
        nrows <- length(psi)/regressors_n
        t(matrix(psi, nrow = nrows, ncol = regressors_n))
      }
      psi_matrix <- psis %>%
        split(rep(1:(periods_n-1), regressors_n*(1:(periods_n-1)))) %>%
        sapply(time_fixed_psi_matrix, regressors_n = regressors_n) %>%
        plyr::rbind.fill.matrix() %>% t() %>% tidyr::replace_na(0) %>%
        rbind(rep(0, (periods_n - 1)*regressors_n))
    }

    phi_matrix + psi_matrix
  } else {
    NULL
  }
  list(O11, O12)
}

SEM_likelihood <- function(cur_Y2, Y1, Y2,
                           alpha, phi_0, err_var, dep_vars,
                           beta = c(), phi_1 = c(),
                           phis = c(), psis = c()) {
  periods_n <- length(dep_vars)
  cur_regressors_n <- length(beta)
  B <- SEM_B_matrix(alpha, periods_n, beta)
  C <- SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1)
  S <- SEM_sigma_matrix(err_var, dep_vars, phis, psis)

  if (det(S[[1]])<=0) {
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
  } else {
    Ui1 <- if (cur_regressors_n == 0) {
      t(B[[1]]%*%t(Y1)-C%*%t(Z))
    } else {
      t(B[[1]]%*%t(Y1)+B[[2]]%*%t(cur_Y2)-C%*%t(Z))
    }
    H=crossprod(Y2-Ui1%*%solve(S[[1]])%*%S[[2]],res_maker_matrix)%*%(Y2-Ui1%*%solve(S[[1]])%*%S[[2]])
    likf=-(n/2)*log(det(S[[1]]))-(1/2)*sum(diag(solve(S[[1]])%*%t(Ui1)%*%Ui1))-(n/2)*log(det(H/n))  # concentrated log-likelihood in the appendix
  }
  return(-likf)
}

lik <- function(t0in) {
  t0=t0in

  ##### Just to see if it affect execution time
  params <- t0
  alpha <- params[1]
  if (cur_regressors_n == 0) {
    beta <- c()
    phi_1 <- c()
  } else {
    beta <- params[2:(1 + cur_regressors_n)]
    phi_1 <- params[(3 + cur_regressors_n):(2 + 2*cur_regressors_n)]
  }
  phis <-
    params[(4 + 2*cur_regressors_n + periods_n):(3 + 2*cur_regressors_n + periods_n + phis_n)]
  psis <-
    params[(4 + 2*cur_regressors_n + periods_n + phis_n):(3 + 2*cur_regressors_n + periods_n + phis_n + psis_n)]
  phi_0 <- params[2 + cur_regressors_n]
  err_var <- params[3 + 2*cur_regressors_n]
  dep_vars <-
    params[(4 + 2*cur_regressors_n):(3 + 2*cur_regressors_n + periods_n)]
  #####

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
  } else {
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
