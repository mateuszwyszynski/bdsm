test_that(paste("likelihoods_summary computes correct approximations of",
                "standard deviations based on economic_growth_ms"), {
  set.seed(23)

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    bdsm::feature_standardization(
      excluded_cols = c(country, year)
    ) %>%
    bdsm::feature_standardization(
      group_by_col  = year,
      excluded_cols = country,
      scale         = FALSE
    )

  lik_info <- likelihoods_summary(df = data_prepared, dep_var_col = gdp,
                                  timestamp_col = year, entity_col = country,
                                  model_space = economic_growth_ms)

  economic_growth_liks <- economic_growth_liks[-3,]
  expect_equal(lik_info, economic_growth_liks)
})
