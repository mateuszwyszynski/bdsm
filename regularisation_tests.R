generate_test_data <- function(n_entities, n_periods) {
  structure(list(
    entities = rep(1:n_entities, n_periods),
    times = rep(1:n_periods, each=n_entities),
    dep_var = stats::rnorm(n_entities * n_periods),
    a = round(stats::rnorm(n_entities * n_periods), 5),
    b = round(stats::rnorm(n_entities * n_periods), 5),
    c = round(stats::rnorm(n_entities * n_periods), 5)
  ),
  class = "data.frame", row.names = c(NA, -n_entities * n_periods))
}

set.seed(1)
# After testing different n_entities and n_periods here,
# and then evaluating the sem likelihood further in the code below,
# I believe that:
#
# 1. There are awlays (n_periods - 2)*(number of regressors) eigenvalues in the matrix H_scaled
# 2. After the first n_entities - 2 eigenvalues, the rest will be zero.
#
# So to have a working model we need:
#
# n_entities - 2 eigenvalues >= (n_periods - 2)*(number of regressors)
#
# Not sure yet why, but this seems to be the case

df_test <- generate_test_data(n_entities = 10, n_periods = 9)

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

sem_value <- sem_likelihood(
  0.5,
  df_test_prepared,
  times, entities, dep_var
)

X <- cbind(
  matrices_shared_across_models$Y1,
  matrices_shared_across_models$Y2,
  matrices_shared_across_models$Z
)

# X <- as.matrix(na.omit(stand_data1[, 3:7]))

check_condition <- function(X) {
  s <- svd(X, nu=0, nv=0)$d
  cond <- max(s)/min(s)
  message("Design matrix has condition number ", format(cond))
  if (cond > 1e8) {
    warning("Very high condition number → severe collinearity")
  }
  invisible(cond)
}

check_condition(X)

diagnose_flat_direction <- function(X, tol = 1e-8) {
  svdX <- svd(X)
  small <- which(svdX$d < tol * max(svdX$d))
  if (length(small)>0) {
    cat("Found", length(small),
        "near‑zero singular values.  Example flat directions:\n")
    print( svdX$v[, small, drop=FALSE] )
  } else {
    message("No extremely small singular values.")
  }
}

diagnose_flat_direction(X)
