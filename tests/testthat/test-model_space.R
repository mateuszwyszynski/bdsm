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


test_that(paste("model_space computes correct model_space list"), {

  data_prepared <- bdsm::economic_growth[,1:5] %>%
    bdsm::feature_standardization(
      excluded_cols = c(country, year, gdp)
    ) %>%
    bdsm::feature_standardization(
      group_by_col  = year,
      excluded_cols = country,
      scale         = FALSE
    )

  for_bma <- model_space(
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
