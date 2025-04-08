"bma_prep_objects_full"

library(magrittr)

data_prepared <- bdsm::economic_growth %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          time_effects = TRUE, scale = FALSE)

library(parallel)
library(bdsm)
cl <- safeMakeCluster()
setDefaultCluster(cl)

bma_prep_objects_full <- bma_prep(df = data_prepared, dep_var_col = gdp,
                             timestamp_col = year, entity_col = country,
                             init_value = 0.5, run_parallel = TRUE)

usethis::use_data(bma_prep_objects_full, overwrite = TRUE)
