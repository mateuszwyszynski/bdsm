#-----------------------------------------------------------------#
#                        HESSIAN FUNCTION                         #
#                        -----------------                        #
#               this procedure takes as input a number            #
#               computes the kxk (2-sided) hessian matrix         #
#-----------------------------------------------------------------#

myhess <- function(lik, theta) {
  x0 <- theta
  k <- nrow(x0) # x0 is theta, our model-specific optimized parameters #
  hessi <- optimbase::zeros(k, k)
  h <- 1e-3

  for (jc in 1:k) {
    for (jr in 1:jc) {
      x1 <- x0
      x2 <- x0
      x3 <- x0
      x4 <- x0
      x1[jr] <- x0[jr] + h
      x1[jc] <- x1[jc] + h
      x2[jr] <- x0[jr] + h
      x2[jc] <- x2[jc] - h
      x3[jr] <- x0[jr] - h
      x3[jc] <- x3[jc] + h
      x4[jr] <- x0[jr] - h
      x4[jc] <- x4[jc] - h
      hessi[jr, jc] <- -(lik(x1) - lik(x2) - lik(x3) + lik(x4)) / (4 * h^2) # the second symmetric derivative has different formula #
      hessi[jc, jr] <- hessi[jr, jc]
    }
  }
  return(hessi)
}
