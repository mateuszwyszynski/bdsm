generate_test_data <- function(n_entities, n_periods) {
  structure(list(
    entities = rep(1:n_entities, n_periods),
    times = rep(1:n_periods, each=n_entities),
    dep_var = stats::rnorm(n_entities * n_periods),
    a = round(stats::rnorm(n_entities * n_periods), 5),
    b = round(stats::rnorm(n_entities * n_periods), 5),
    c = round(stats::rnorm(n_entities * n_periods), 5),
    d = round(stats::rnorm(n_entities * n_periods), 5),
    e = round(stats::rnorm(n_entities * n_periods), 5)
  ),
  class = "data.frame", row.names = c(NA, -n_entities * n_periods))
}

# One can generate a random data with a different number of regressors using the code below.
# Then we can test how sem_likelihood behaves.
# If you add some debugging messages, you can print out the eigenvalues of the matrix H and check their number.
#
# I've done multiple such tests. After going through them, I arrived at the following conclusions:
#
# 1. The symmetric matrix H awlays has dimension (n_periods - 2)*(number of regressors),
# so it will have this number of eigenvalues
# 2. The rank of H is bounded by rank(H) = rank(M^T P M) = rank( M^T P^{1/2} P^{1/2} M) = rank (P^{1/2} M) <= min(rank(P), rank(M))
# where:
#  - the P^{1/2} exists because P is a projection matrix, i.e., symmetric semi positive definite (all eigenvalues >= 0)
#  - we have rank(X^T X) = rank(X)
#  - we have rank(AB) <= min(rank(A), rank(B))
#
# Combining 1. and 2. together leads to the following necessary condition for the H matrix to be symmetric positive definite
# (all eigenvalues > 0):
#
# n_entities - n_regressors - 1 >= (n_periods - 2)*n_regressors
#
# where the RHS is simply the number of eigenvalues the H has, and the LHS is the
#
# rank(P) = rank(I) - rank(Z) = N - (n_regressors + 1)
#
# If one wants to express the condition as a bound on the number of entities needed for the approach to work:
#
# n_entities >= (n_periods - 1)*n_regressors + 1 (**)
#
# IMPORTANT: For some reason in my experiments with randomly generated data I always need to satisfy a strict inequality in (**).
# I'm not sure yet what is the reason for that.
# EDIT: I believe that the reason for the strict inequality is that the matrix Z
# contains the lagged dependent variable as one of its columns. The residual
# maker matrix maps onto the space orthogonal to the space spanned by the
# columns of Z. The matrix H is computed as M^T P M. I believe that since Z
# contains the lagged dependent variable, and M contains both the dependent
# variable and its lagged copy, one dimension is always or almost always dropped
# when H is calculated.
#
# Example:
# n_regressors = 5
# n_periods = 12
#
# so 11*5 + 1 = 56 entities are necessary to even have a chance of having non-zero eigenvalues in H.

k = 5
t = 12

set.seed(1)
df_test <- generate_test_data(n_entities = 55, n_periods = t)

matrices_shared_across_models <- df_test %>%
  matrices_from_df(timestamp_col = times,
                   entity_col = entities,
                   dep_var_col = dep_var,
                   which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

df_test_prepared <- df_test %>%
  bdsm::feature_standardization(
    excluded_cols = c(entities, times, dep_var)
  ) %>%
  bdsm::feature_standardization(
    group_by_col  = times,
    excluded_cols = entities,
    scale         = FALSE
  )

lin_rel_regs <- c("a","b","c","d","e")
params <- generate_params_vector(0.5, t, k, length(lin_rel_regs))

# for (i in 1:length(params)) {
#   params[i] = params[i]*(-1)^i
# }
#
# params <- stats::rnorm(length(params))

sem_value <- sem_likelihood(
  params,
  df_test_prepared,
  times, entities, dep_var,
  lin_related_regressors = lin_rel_regs
)

X <- cbind(
  matrices_shared_across_models$Y1,
  matrices_shared_across_models$Y2
)

check_condition <- function(X) {
  s <- svd(X, nu=0, nv=0)$d
  cond <- max(s)/min(s)
  message("Design matrix has condition number ", format(cond))
  if (cond > 1e8) {
    warning("Very high condition number -> severe collinearity")
  }
  invisible(cond)
}

check_condition(X)

diagnose_flat_direction <- function(X, tol = 1e-8) {
  svdX <- svd(X)
  small <- which(svdX$d < tol * max(svdX$d))
  if (length(small)>0) {
    cat("Found", length(small),
        "near-zero singular values.  Example flat directions:\n")
    print( svdX$v[, small, drop=FALSE] )
  } else {
    message("No extremely small singular values.")
  }
}

diagnose_flat_direction(X)
