test_that("SEM_B_matrix computes proper matrix", {
  B <- as.matrix(SEM_B_matrix(3:6, 3, 4))
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
