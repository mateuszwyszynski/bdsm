# code to prepare `small_economic_growth_liks` dataset
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

small_economic_growth_liks <- likelihoods_summary(
  df            = data_prepared,
  dep_var_col   = gdp,
  timestamp_col = year,
  entity_col    = country,
  model_space   = small_economic_growth_ms
)

usethis::use_data(small_economic_growth_liks, overwrite = TRUE)
