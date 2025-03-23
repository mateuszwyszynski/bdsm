test_that(paste("coef_hist creates correct lists with graphs"), {

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  bma_results <- bma(bma_prep_objects, df = data_prepared, round= 3, dilution = 0)

  coef_plots <- coef_hist(bma_results, kernel = 1)

  expect_equal(class(coef_plots), "list")
  expect_equal(class(coef_plots[[1]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[2]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[3]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[4]]), c("gg","ggplot"))
  expect_equal(class(coef_plots[[5]]), c("gg","ggplot"))
})
