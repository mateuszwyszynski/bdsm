generate_test_data <- function() {
  n_entities <- 35
  n_periods <- 12
  structure(list(
    entities = rep(1:n_entities, n_periods),
    times = rep(1:n_periods, each=n_entities),
    dep_var = stats::rnorm(n_entities * n_periods),
    a = round(stats::rnorm(n_entities * n_periods), 5),
    b = round(stats::rnorm(n_entities * n_periods), 5),
    c = round(stats::rnorm(n_entities * n_periods), 5)
  ),
  class = "data.frame", row.names = c(NA, -n_entities * n_periods))
}

generate_test_feature_standard_data <- function(time_effects = FALSE, scale = TRUE) {
  feature_standardization(
    generate_test_data(),
    timestamp_col = times,
    entity_col = entities,
    time_effects,
    scale
  )
}

generate_test_data_prep_data <- function(standardize = TRUE, time_effects = TRUE, entity_effects = TRUE,
                                    scale = TRUE) {
  data_prep(
    generate_test_data(),
    timestamp_col = times,
    entity_col = entities,
    time_effects,
    entity_effects,
    scale
  )
}
