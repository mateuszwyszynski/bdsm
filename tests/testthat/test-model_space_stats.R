test_that(paste("model_space_stats computes correct likelihoods and",
                "standard deviations based on small_model_space"), {
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

  liks_info <- model_space_stats(
    df            = data_prepared,
    dep_var_col   = gdp,
    timestamp_col = year,
    entity_col    = country,
    model_space   = small_model_space[[1]]
  )

  expect_equal(liks_info, small_model_space[[2]])
})
