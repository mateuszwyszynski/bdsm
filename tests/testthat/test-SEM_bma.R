test_that(paste("bma_stds computes correct approximations based on",
                "economic_growth_ms"), {
  set.seed(23)

  data_prepared <- panels::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            cross_sectional = TRUE, scale = FALSE)

  std_devs <- bma_stds(df = data_prepared, dep_var_col = gdp,
                       timestamp_col = year, entity_col = country,
                       model_space = economic_growth_ms,
                       projection_matrix_const = TRUE)
  expected_std_devs <- list(stds = std_devs$stds,
                            stds_robust = std_devs$stds_robust)
  expect_equal(std_devs, expected_std_devs)
})

test_that(paste("likelihoods_summary computes correct approximations of",
                "standard deviations based on economic_growth_ms"), {
  set.seed(23)

  data_prepared <- panels::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            cross_sectional = TRUE, scale = FALSE)

  lik_info <- likelihoods_summary(df = data_prepared, dep_var_col = gdp,
                                  timestamp_col = year, entity_col = country,
                                  model_space = economic_growth_ms,
                                  projection_matrix_const = TRUE)
  expected_lik_info <- rbind(economic_growth_stds, economic_growth_stds_robust)
  expect_equal(lik_info[-1,], expected_lik_info)
})
