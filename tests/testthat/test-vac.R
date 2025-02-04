test_that("vac function produces expected values", {

  # Define test parameters
  x <- c(0, 25, 50, 75, 100)
  maxv <- 0.8
  scale <- 10
  start <- 50

  # Expected calculations
  expected_values <- pmax(0, tanh((x - start) / scale) * maxv)

  # Test function output
  expect_equal(vac(x, maxv, scale, start), expected_values)

  # Test edge cases
  expect_equal(vac(start, maxv, scale, start), 0)  # At start point, should be near zero
  expect_equal(vac(start - 100, maxv, scale, start), 0)  # Far below start, should return 0
  expect_equal(vac(start + 100, maxv, scale, start), maxv)  # Far above start, should approach maxv
})
