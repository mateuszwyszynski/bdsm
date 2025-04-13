# code to prepare `economic_growth_liks` dataset
library(magrittr)

set.seed(23)

data_prepared <- bdsm::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

economic_growth_liks <-
  likelihoods_summary(df = data_prepared, dep_var_col = gdp,
                      timestamp_col = year, entity_col = country,
                      economic_growth_ms)

usethis::use_data(economic_growth_liks, overwrite = TRUE)
