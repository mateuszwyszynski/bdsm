# microbenchmark(
#   sem_likelihood(
#     0.5,
#     generate_test_feature_standard_data(),
#     times, entities, dep_var
#   )
# )
# Unit: milliseconds
# min      lq       mean.    median   uq       max        neval
# 25.80372 26.05294 27.14582 26.22264 27.25252 36.02822   100

test_that("SEM likelihood is calculated correctly for default feature standardization parameters", {
  skip_on_cran()
  set.seed(1)
  sem_value <- sem_likelihood(
    0.5,
    feature_standardization(
      df            = generate_test_data(),
      excluded_cols = c(entities, times)
    ),
    times, entities, dep_var
  )
  expect_equal(sem_value, -2621.65382)
})

test_that("SEM likelihood is calculated correctly for time_effects TRUE", {
  skip_on_cran()
  set.seed(1)
  sem_value <- sem_likelihood(
    0.5,
    feature_standardization(
      df            = generate_test_data(),
      group_by_col  = times,
      excluded_cols = entities
    ),
    times, entities, dep_var
  )
  expect_equal(sem_value, -2566.28571)
})

test_that("SEM likelihood is calculated correctly for time_effects TRUE and scale FALSE", {
  skip_on_cran()
  set.seed(1)
  sem_value <- sem_likelihood(
    0.5,
    feature_standardization(
      df            = generate_test_data(),
      group_by_col  = times,
      excluded_cols = entities,
      scale = FALSE
    ),
    times, entities, dep_var
  )
  expect_equal(sem_value, -2520.38776)
})

test_that("SEM likelihood is calculated correctly for time_effects FALSE and scale FALSE", {
  skip_on_cran()
  sem_value <- sem_likelihood(
    0.5,
    feature_standardization(
      df            = generate_test_data(),
      excluded_cols = c(times, entities),
      scale         = FALSE
    ),
    times, entities, dep_var
  )
  expect_equal(sem_value, -2801.069937)
})
