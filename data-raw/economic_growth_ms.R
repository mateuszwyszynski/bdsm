# code to prepare `economic_growth_ms` dataset
library(magrittr)

set.seed(23)

data_prepared <- bdsm::economic_growth[, 1:6] %>%
  bdsm::feature_standardization(
    excluded_cols = c(country, year, gdp)
  ) %>%
  bdsm::feature_standardization(
    group_by_col  = year,
    excluded_cols = country,
    scale         = FALSE
  )

economic_growth_ms <- optimal_model_space(
  df            = data_prepared,
  dep_var_col   = gdp,
  timestamp_col = year,
  entity_col    = country,
  init_value    = 0.5
)

usethis::use_data(economic_growth_ms, overwrite = TRUE)
