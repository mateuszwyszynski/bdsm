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

model1 <- optim_model_space(df=stand_data1,
                            timestamp_col = year,
                            entity_col = Country,
                            dep_var_col= Eint,
                            init_value = 0.5)


bma1 <- bma(model1,df = stand_data1)

bma1[[1]]
bma1[[2]]
