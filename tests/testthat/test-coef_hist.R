test_that(paste("coef_hist creates correct lists with graphs"), {

  data_prepared <- bdsm::economic_growth[,1:6] %>%
    bdsm::feature_standardization(
      excluded_cols = c(country, year, gdp)
    ) %>%
    bdsm::feature_standardization(
      group_by_col  = year,
      excluded_cols = country,
      scale         = FALSE
    )

  bma_results <- bma(small_model_space, df = data_prepared, round= 3, dilution = 0)

  coef_plots <- coef_hist(bma_results, kernel = 1)

  expect_equal(class(coef_plots), "list")
  expect_true(ggplot2::is_ggplot(coef_plots[[1]]))
  expect_true(ggplot2::is_ggplot(coef_plots[[2]]))
  expect_true(ggplot2::is_ggplot(coef_plots[[3]]))
  expect_true(ggplot2::is_ggplot(coef_plots[[4]]))
})
