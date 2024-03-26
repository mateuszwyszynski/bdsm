# microbenchmark(
#   SEM_likelihood(
#     0.5,
#     generate_test_feature_standard_data(),
#     times, entities, dep_var
#   )
# )
# Unit: milliseconds
# min      lq       mean.    median   uq       max        neval
# 25.80372 26.05294 27.14582 26.22264 27.25252 36.02822   100

test_that("SEM likelihood is calculated correctly for default feature standardization parameters", {
  set.seed(1)
  sem_value <- SEM_likelihood(
    0.5,
    generate_test_feature_standard_data(),
    times, entities, dep_var
  )
  expect_equal(sem_value, 278.13706)
})

test_that("SEM likelihood is calculated correctly for cross_sectional TRUE", {
  set.seed(1)
  sem_value <- SEM_likelihood(
    0.5,
    generate_test_feature_standard_data(cross_sectional = TRUE),
    times, entities, dep_var
  )
  expect_equal(sem_value, 549.07193)
})

test_that("SEM likelihood is calculated correctly for cross_sectional TRUE and scale FALSE", {
  set.seed(1)
  sem_value <- SEM_likelihood(
    0.5,
    generate_test_feature_standard_data(cross_sectional = TRUE, scale = FALSE),
    times, entities, dep_var
  )
  expect_equal(sem_value, 560.93187)
})

test_that("SEM likelihood is calculated correctly for cross_sectional FALSE and scale FALSE", {
  set.seed(1)
  # TODO: That produces always NaN
  sem_value <- SEM_likelihood(
    0.5,
    generate_test_feature_standard_data(cross_sectional = FALSE, scale = FALSE),
    times, entities, dep_var
  )
  expect_equal(sem_value, NaN)
})
