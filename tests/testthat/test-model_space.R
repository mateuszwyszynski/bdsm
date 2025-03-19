test_that("optimal_model_space correctly computes economic_growth_ms", {
  skip_on_os(c("windows", "linux"))
  skip_on_cran()
  set.seed(23)

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  optim_cpp <- TRUE
  model_space <- optimal_model_space(df = data_prepared, dep_var_col = gdp,
                                     timestamp_col = year, entity_col = country,
                                     init_value = 0.5, optim_cpp = optim_cpp)

  expect_equal(model_space, economic_growth_ms)
})
