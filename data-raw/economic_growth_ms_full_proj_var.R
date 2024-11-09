# code to prepare `economic_growth_ms` dataset
library(magrittr)

set.seed(23)

data_prepared <- panels::economic_growth %>%
  feature_standardization(timestamp_col = year, entity_col = gdp) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

economic_growth_ms_full_proj_var <- optimal_model_space(df = data_prepared, dep_var_col = gdp,
                                                        timestamp_col = year, entity_col = country,
                                                        init_value = 0.5)

usethis::use_data(economic_growth_ms_full_proj_var, overwrite = TRUE)
