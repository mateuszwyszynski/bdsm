library(parallel)
cores <- detectCores()
cl <- makeCluster(cores)
setDefaultCluster(cl)

energy <- readxl::read_excel("~/Desktop/Energy_data.xlsx")[, 1:8]

new_data1 <- join_lagged_col(energy, col = Eint, col_lagged = Eint_lag, entity_col = Country, timestamp_col = year, timestep = 1)

# mig1
stand_data1 <- feature_standardization(df = new_data1,  excluded_cols = c(year,Country))
stand_data1 <- feature_standardization(df = stand_data1, group_by_col = year, excluded_cols = Country, scale = FALSE)

matrices_shared_across_models <- stand_data1 %>%
  matrices_from_df(timestamp_col = year,
                   entity_col = Country,
                   dep_var_col = Eint,
                   which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

matrices_shared_across_models <- data_prepared %>%
  matrices_from_df(timestamp_col = year,
                   entity_col = country,
                   dep_var_col = gdp,
                   which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

generate_test_data <- function(n_entities, n_periods) {
  structure(list(
    entities = rep(1:n_entities, n_periods),
    times = rep(1:n_periods, each=n_entities),
    dep_var = stats::rnorm(n_entities * n_periods),
    a = round(stats::rnorm(n_entities * n_periods), 3),
    b = round(stats::rnorm(n_entities * n_periods), 3)
  ),
  class = "data.frame", row.names = c(NA, -n_entities * n_periods))
}

set.seed(1)
df_test <- generate_test_data(n_entities = 10, n_periods = 5)

matrices_shared_across_models <- df_test %>%
  matrices_from_df(timestamp_col = times,
                   entity_col = entities,
                   dep_var_col = dep_var,
                   which_matrices = c("Y1", "Y2", "Z", "res_maker_matrix"))

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

df_test_prepared <- df_test %>%
  bdsm::feature_standardization(
    excluded_cols = c(entities, times, dep_var)
  ) %>%
  bdsm::feature_standardization(
    group_by_col  = times,
    excluded_cols = entities,
    scale         = FALSE
  )

sem_value6 <- sem_likelihood(
  0.5,
  df_test_prepared,
  times, entities, dep_var
)

model1 <- optim_model_space(df=stand_data1,
                            timestamp_col = year,
                            entity_col = Country,
                            dep_var_col= Eint,
                            init_value = 0.5)


bma1 <- bma(model1,df = stand_data1)

bma1[[1]]
bma1[[2]]
