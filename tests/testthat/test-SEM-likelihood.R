test_that("SEM_B_matrix computes proper matrix", {
  B <- as.matrix(SEM_B_matrix(3, 4, 4:6))
  B_expected_data <- c(
    1, -3, rep(0, 11),
    0, 1, -3, rep(0, 10),
    0, 0, 1, -3, rep(0, 9),
    rep(0, 3), 1, rep(0, 9),
    rep(0, 1), -4, rep(0, 2), 1, rep(0, 8),
    rep(0, 1), -5, rep(0, 3), 1, rep(0, 7),
    rep(0, 1), -6, rep(0, 4), 1, rep(0, 6),
    rep(0, 2), -4, rep(0, 4), 1, rep(0, 5),
    rep(0, 2), -5, rep(0, 5), 1, rep(0, 4),
    rep(0, 2), -6, rep(0, 6), 1, rep(0, 3),
    rep(0, 3), -4, rep(0, 6), 1, rep(0, 2),
    rep(0, 3), -5, rep(0, 7), 1, rep(0, 1),
    rep(0, 3), -6, rep(0, 8), 1, rep(0, 0)
  )
  B_expected <- matrix(B_expected_data, 13, 13)
  expect_equal(B, B_expected, ignore_attr = TRUE)
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
