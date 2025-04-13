test_that("optim_model_space_params correctly computes small_economic_growth_ms", {
  skip_on_os(c("windows", "linux"))
  skip_on_cran()
  set.seed(23)

  data_prepared <- bdsm::economic_growth[,1:6] %>%
    bdsm::feature_standardization(
      excluded_cols = c(country, year, gdp)
    ) %>%
    bdsm::feature_standardization(
      group_by_col  = year,
      excluded_cols = country,
      scale         = FALSE
    )

  model_space <- optim_model_space_params(
    df            = data_prepared,
    dep_var_col   = gdp,
    timestamp_col = year,
    entity_col    = country,
    init_value    = 0.5
  )

  expect_equal(model_space, small_model_space[[1]])
})
