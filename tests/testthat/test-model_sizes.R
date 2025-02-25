test_that(paste("model_sizes creates correct lists with graphs"), {

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  bma_results <- bma(bma_prep_objects, df = data_prepared, app = 3, dilution = 0)

  size_graphs <- model_sizes(bma_results)

  expect_equal(class(size_graphs), "list")
  expect_equal(length(size_graphs), 3)
  expect_equal(class(size_graphs[[1]]), c("gg","ggplot"))
  expect_equal(class(size_graphs[[2]]), c("gg","ggplot"))
  expect_equal(class(size_graphs[[3]]), c("gg","ggplot","ggarrange"))
})
