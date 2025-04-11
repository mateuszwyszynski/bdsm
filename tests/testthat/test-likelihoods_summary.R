test_that(paste("likelihoods_summary computes correct likelihoods and",
                "standard deviations based on small_economic_growth_ms"), {
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

  liks_info <- likelihoods_summary(
    df            = data_prepared,
    dep_var_col   = gdp,
    timestamp_col = year,
    entity_col    = country,
    model_space   = small_economic_growth_ms
  )

  expect_equal(liks_info, small_economic_growth_liks)
})
