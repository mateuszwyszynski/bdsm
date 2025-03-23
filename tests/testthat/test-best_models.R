test_that(paste("best_models creates correct lists with graphs"), {

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  bma_results <- bma(bma_prep_objects, df = data_prepared, round= 3, dilution = 0)

  best <- 5

  best_5_models <- best_models(bma_results, criterion = 2, best = best, estimate = TRUE, robust = TRUE)

  expect_equal(class(best_5_models), "list")
  expect_equal(length(best_5_models), 9)
  expect_equal(class(best_5_models[[1]]), c("matrix", "array"))
  expect_equal(class(best_5_models[[2]]), c("matrix", "array"))
  expect_equal(class(best_5_models[[3]]), c("matrix", "array"))
  expect_equal(ncol(best_5_models[[1]]), best)
  expect_equal(ncol(best_5_models[[2]]), best)
  expect_equal(ncol(best_5_models[[3]]), best)
  expect_equal(class(best_5_models[[4]]), "knitr_kable")
  expect_equal(class(best_5_models[[5]]), "knitr_kable")
  expect_equal(class(best_5_models[[6]]), "knitr_kable")
  expect_equal(class(best_5_models[[7]]), c("gTree", "grob",  "gDesc"))
  expect_equal(class(best_5_models[[8]]), c("gTree", "grob",  "gDesc"))
  expect_equal(class(best_5_models[[9]]), c("gTree", "grob",  "gDesc"))
})
