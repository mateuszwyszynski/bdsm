test_that("SEM_B_matrix computes proper matrix", {
  periods_n <- 4
  B <- SEM_B_matrix(3, periods_n, 4:6)
  B11_expected_data <- c(
    1, 0, 0, 0,
    -3, 1, 0, 0,
    0, -3, 1, 0,
    0, 0, -3, 1
  )
  B12_expected_data <- c(
    0, -4, 0, 0,
    0, -5, 0, 0,
    0, -6, 0, 0,
    0, 0, -4, 0,
    0, 0, -5, 0,
    0, 0, -6, 0,
    0, 0, 0, -4,
    0, 0, 0, -5,
    0, 0, 0, -6
  )
  B11_expected <- matrix(B11_expected_data, nrow = periods_n, byrow = TRUE)
  B12_expected <- matrix(B12_expected_data, nrow = periods_n)
  expect_equal(B[[1]], B11_expected, ignore_attr = TRUE)
  expect_equal(as.matrix(B[[2]]), B12_expected, ignore_attr = TRUE)
})

test_that("SEM_C_matrix computes proper matrix", {
  alpha <- 9
  phi_0 <- 19
  beta <- 11:15
  phi_1 <- 21:25
  periods_n <- 4
  C <- as.matrix(SEM_C_matrix(alpha, phi_0, periods_n, beta, phi_1))

  C_expected_data <- c(
    alpha + phi_0, rep(phi_0, periods_n-1),
    beta[1] + phi_1[1], rep(phi_1[1], periods_n-1),
    beta[2] + phi_1[2], rep(phi_1[2], periods_n-1),
    beta[3] + phi_1[3], rep(phi_1[3], periods_n-1),
    beta[4] + phi_1[4], rep(phi_1[4], periods_n-1),
    beta[5] + phi_1[5], rep(phi_1[5], periods_n-1)
  )
  C_expected <- matrix(C_expected_data, periods_n, 1 + length(beta))
  expect_equal(C, C_expected, ignore_attr = TRUE)
})

test_that("SEM_sigma_matrix computes proper matrix", {
  err_var <- 1
  dep_vars <- c(2, 2, 2, 2)
  phis <- c(10, 10, 20, 20, 30, 30)
  psis <- c(101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112)
  sigma <- as.matrix(SEM_sigma_matrix(err_var, dep_vars, phis, psis))

  sigma_11_expected_data <- c(
    err_var + dep_vars[1], err_var, err_var, err_var,
    err_var, err_var + dep_vars[2], err_var, err_var,
    err_var, err_var, err_var + dep_vars[3], err_var,
    err_var, err_var, err_var, err_var + dep_vars[4]
  )

  sigma_12_expected_data <- c(
    phis[1] + psis[1], phis[2] + psis[2], phis[3] + psis[3], phis[4] + psis[5],
    phis[5] + psis[7], phis[6] + psis[10],
    phis[1], phis[2], phis[3] + psis[4], phis[4] + psis[6],
    phis[5] + psis[8], phis[6] + psis[11],
    phis[1], phis[2], phis[3], phis[4],
    phis[5] + psis[9], phis[6] + psis[12],
    phis[1], phis[2], phis[3], phis[4], phis[5], phis[6]
  )

  sigma_11_expected <- matrix(sigma_11_expected_data, nrow = 4, byrow = TRUE)
  sigma_12_expected <- matrix(sigma_12_expected_data, nrow = 4, byrow = TRUE)
  expect_equal(sigma[[1]], sigma_11_expected, ignore_attr = TRUE)
  expect_equal(sigma[[2]], sigma_12_expected, ignore_attr = TRUE)
})
