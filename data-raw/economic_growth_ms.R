# code to prepare `economic_growth_ms` dataset
library(magrittr)

set.seed(23)

data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

economic_growth_ms <- optimal_model_space(df = data_prepared, dep_var_col = gdp,
                                   timestamp_col = year, entity_col = country,
                                   init_value = 0.5,
                                   projection_matrix_const = TRUE)

usethis::use_data(economic_growth_ms, overwrite = TRUE)
