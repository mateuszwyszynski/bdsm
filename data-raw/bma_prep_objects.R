"bma_prep_objects"

library(magrittr)

data_prepared <- bdsm::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          time_effects = TRUE, scale = FALSE)

bma_prep_objects <- bma_prep(df = data_prepared, dep_var_col = gdp,
                    timestamp_col = year, entity_col = country,
                    init_value = 0.5)

usethis::use_data(bma_prep_objects, overwrite = TRUE)
