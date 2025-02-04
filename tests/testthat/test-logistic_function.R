test_that("logistic_function computes expected values", {

  # Define test parameters
  t <- 5
  R0 <- 1
  R1 <- 10
  r <- 0.5
  t0 <- 3

  # Expected calculation
  K <- R1
  L <- R0
  expected_value <- (K - L) / (1 + exp(-r * (t - t0))) + L

  # Test function output
  expect_equal(logistic_function(t, R0, R1, r, t0), expected_value)

  # Test with another set of parameters
  t <- 10
  R0 <- 2
  R1 <- 20
  r <- 1
  t0 <- 5

  expected_value <- (R1 - R0) / (1 + exp(-r * (t - t0))) + R0
  expect_equal(logistic_function(t, R0, R1, r, t0), expected_value)

  # Test edge case where t = t0
  t <- t0
  expected_value <- (R1 - R0) / 2 + R0
  expect_equal(logistic_function(t, R0, R1, r, t0), expected_value)
})
