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
  expect_equal(class(coef_plots[[1]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[2]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[3]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[4]]), c("gg","ggplot"))
})
