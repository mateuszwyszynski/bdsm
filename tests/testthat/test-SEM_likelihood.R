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
  skip_on_os(c("windows", "linux"))
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
  expect_equal(sem_value, 133.223858)
})

test_that("SEM likelihood is calculated correctly for time_effects TRUE", {
  skip_on_os(c("windows", "linux"))
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
  expect_equal(sem_value, 217.693805)
})

test_that("SEM likelihood is calculated correctly for time_effects TRUE and scale FALSE", {
  skip_on_os(c("windows", "linux"))
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
  expect_equal(sem_value, 225.54665)
})

test_that("SEM likelihood is calculated correctly for time_effects FALSE and scale FALSE", {
  skip_on_os(c("windows", "linux"))
  skip_on_cran()
  set.seed(1)
  sem_value <- sem_likelihood(
    0.5,
    feature_standardization(
      df            = generate_test_data(),
      excluded_cols = c(times, entities),
      scale         = FALSE
    ),
    times, entities, dep_var
  )
  expect_equal(sem_value, 140.498138)
})

test_that("SEM likelihood is calculated incorrectly for specific data", {
  skip_on_os(c("windows", "linux"))
  skip_on_cran()
  set.seed(2)
  # TODO: That produces NaN for that particular seed.
  testthat::expect_warning(
    sem_value <- sem_likelihood(
      0.5,
      feature_standardization(
        df            = generate_test_data(),
        excluded_cols = c(times, entities),
        scale         = FALSE
      ),
      times, entities, dep_var
    )
  )
})
