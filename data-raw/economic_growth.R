## code to prepare `economic_growth` dataset goes here

fpath <- system.file("extdata", "economic_growth_raw.xlsx", package = "bdsm")

economic_growth_raw <- readxl::read_excel(fpath)

economic_growth <- economic_growth_raw %>%
  join_lagged_col(gdp, lag_gdp, year, country, 10)

usethis::use_data(economic_growth, overwrite = TRUE)
