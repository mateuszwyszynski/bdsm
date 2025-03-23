test_that(paste("model_pmp creates correct lists with graphs"), {

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  bma_results <- bma(bma_prep_objects, df = data_prepared, round= 3, dilution = 0)

  model_graphs <- model_pmp(bma_results, top = 16)

  expect_equal(class(model_graphs), "list")
  expect_equal(length(model_graphs), 3)
  expect_equal(class(model_graphs[[1]]), c("gg","ggplot"))
  expect_equal(class(model_graphs[[2]]), c("gg","ggplot"))
  expect_equal(class(model_graphs[[3]]), c("gg","ggplot","ggarrange"))
})
