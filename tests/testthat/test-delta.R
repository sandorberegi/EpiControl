test_that("delta function produces expected values", {

  # Define test parameters
  x <- c(0, 25, 50, 75, 100)
  scale <- 10
  start <- 50

  # Expected calculations
  expected_values <- tanh((x - start) / scale) * 0.5 + 0.5

  # Test function output
  expect_equal(delta(x, scale, start), expected_values)

  # Test edge cases
  expect_equal(delta(start, scale, start), 0.5)  # At start point, should return 0.5
  expect_equal(delta(start - 100, scale, start), 0)  # Far below start, should approach 0
  expect_equal(delta(start + 100, scale, start), 1)  # Far above start, should approach 1
})
