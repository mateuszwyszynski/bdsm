test_that(paste("bma_prep computes correct bma_prep list"), {

  data_prepared <- bdsm::economic_growth[,1:5] %>%
    bdsm::feature_standardization(
      excluded_cols = c(country, year, gdp)
    ) %>%
    bdsm::feature_standardization(
      group_by_col  = year,
      excluded_cols = country,
      scale         = FALSE
    )

  for_bma <- bma_prep(
    df            = data_prepared,
    dep_var_col   = gdp,
    timestamp_col = year,
    entity_col    = country,
    init_value    = 0.5
  )

  expect_equal(length(for_bma), 2)
  expect_equal(class(for_bma), "list")
  expect_equal(class(for_bma[[1]]), c("matrix","array"))
  expect_equal(class(for_bma[[2]]), c("matrix","array"))
})


