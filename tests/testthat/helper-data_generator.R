generate_test_data <- function() {
  structure(list(
    entities = c(
      1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L, 3L, 4L, 1L, 2L,
      3L, 4L
    ),
    times = c(
      1960, 1960, 1960, 1960, 1970, 1970, 1970, 1970, 1980, 1980, 1980, 1980,
      1990, 1990, 1990, 1990, 2000, 2000, 2000, 2000
    ),
    dep_var = c(
      -0.62645, 0.18364, -0.83562, 1.59528, 0.32950, -0.82046, 0.48742, 0.73832,
      0.57578, -0.30538, 1.51178, 0.38984, -0.62124, -2.21469, 1.12493,
      -0.04493, -0.01619, 0.94383, 0.82122, 0.59390
    ),
    a = c(
      0.91897737, 0.78213630, 0.074564983, -1.9893516, 0.6198257, -0.056128739,
      -0.15579550, -1.4707523, -0.4781500, 0.41794156, 1.3586795, -0.10278772,
      0.38767161, -0.053805040, -1.3770595, -0.4149945, -0.39428995,
      -0.059313396, 1.1000253, 0.76317574
    ),
    b = c(
      -0.16452359, -0.25336168, 0.69696337, 0.55666319, -0.6887556, -0.7074951,
      0.3645819, 0.76853292, -0.11234621, 0.88110772, 0.39810588, -0.61202639,
      0.34111969, -1.1293630, 1.4330237, 1.9803998, -0.36722147, -1.0441346,
      0.56971962, -0.13505460
    )
  ),
  class = "data.frame", row.names = c(NA, -20L))
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
