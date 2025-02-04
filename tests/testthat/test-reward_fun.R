test_that("reward_fun calculates correct reward values", {

  # Create a mock episimdata data frame
  episimdata <- data.frame(
    C = c(50, 100, 200),
    Re = c(1.5, 2.0, 2.5)
  )

  # Define parameters
  alpha <- 2
  ovp <- 10
  C_target <- 100
  C_target_pen <- 150
  R_target <- 2.0
  actions <- data.frame(cost_of_NPI = c(5, 15, 25))
  ii <- 2
  jj <- 1

  # Expected calculations
  C_err_pred <- abs(episimdata[ii, 'C'] - C_target) # abs(100 - 100) = 0
  R_err_pred <- abs(episimdata[ii, 'Re'] - R_target) # abs(2.0 - 2.0) = 0
  over_pen <- 0  # No over-penalty since C (100) is not above C_target_pen (150)
  expected_rew <- -alpha * C_err_pred - over_pen - actions[jj, "cost_of_NPI"] # -0 - 0 - 5 = -5

  # Test function output
  expect_equal(reward_fun(episimdata, alpha, ovp, C_target, C_target_pen, R_target, actions, ii, jj), expected_rew)

  # Test with a case exceeding the penalty threshold
  ii <- 3 # Using the third row (C = 200)
  C_err_pred <- abs(episimdata[ii, 'C'] - C_target) # abs(200 - 100) = 100
  over_pen <- ovp # Over-penalty applied since C (200) > C_target_pen (150)
  expected_rew <- -alpha * C_err_pred - over_pen - actions[jj, "cost_of_NPI"] # -200 - 10 - 5 = -215

  expect_equal(reward_fun(episimdata, alpha, ovp, C_target, C_target_pen, R_target, actions, ii, jj), expected_rew)
})
