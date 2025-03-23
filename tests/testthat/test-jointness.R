test_that(paste("jointness computes correct jointness table"), {

  data_prepared <- bdsm::economic_growth[,1:7] %>%
    feature_standardization(timestamp_col = year, entity_col = country) %>%
    feature_standardization(timestamp_col = year, entity_col = country,
                            time_effects = TRUE, scale = FALSE)

  bma_results <- bma(bma_prep_objects, df = data_prepared, round= 3, dilution = 0)

  jointness_table <- jointness(bma_results, measure = "HCGHM", rho = 0.5, round= 3)

  expect_equal(nrow(jointness_table), bma_results[[4]])
  expect_equal(ncol(jointness_table), bma_results[[4]])
  expect_equal(class(jointness_table), c("matrix","array"))
  expect_equal(is.na(jointness_table[[1]]), TRUE)
})
