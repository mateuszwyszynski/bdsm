# code to prepare `economic_growth_bma_params` dataset
library(magrittr)

set.seed(20)

data_prepared <- bdsm::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

regressors <- regressor_names(data_prepared, year, country, gdp)

bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = economic_growth_ms)

economic_growth_bma_params <-
  parameters_summary(regressors = regressors,
                     bet = bma_result$bet, pvarh = bma_result$pvarh,
                     pvarr = bma_result$pvarr, fy = bma_result$fy,
                     fyt = bma_result$fyt, ppmsize = bma_result$ppmsize,
                     cout = bma_result$cout, nts = bma_result$nts,
                     pts = bma_result$pts, variables_n = bma_result$variables_n)

usethis::use_data(economic_growth_bma_params, overwrite = TRUE)
