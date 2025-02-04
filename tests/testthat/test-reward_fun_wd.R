test_that("reward_fun_wd calculates correct reward values", {

  # Create a mock episimdata data frame
  episimdata <- data.frame(
    C = c(50, 100, 200),
    Re = c(1.5, 2.0, 2.5),
    Deaths = c(5, 10, 20)
  )

  # Define parameters
  alpha <- 2
  alpha_d <- 3
  ovp <- 10
  dovp <- 15
  C_target <- 100
  C_target_pen <- 150
  D_target <- 10
  D_target_pen <- 15
  actions <- data.frame(cost_of_NPI = c(5, 15, 25))
  ii <- 2
  jj <- 1

  # Expected calculations
  C_err_pred <- abs(episimdata[ii, 'C'] - C_target) # abs(100 - 100) = 0
  D_err_pred <- abs(episimdata[ii, 'Deaths'] - D_target) # abs(10 - 10) = 0
  over_pen <- 0  # No over-penalty since C (100) is not above C_target_pen (150)
  over_pen_d <- 0  # No over-penalty since Deaths (10) is not above D_target_pen (15)
  expected_rew <- -alpha * C_err_pred - over_pen - actions[jj, "cost_of_NPI"] - alpha_d * D_err_pred - over_pen_d # -0 - 0 - 5 - 0 - 0 = -5

  # Test function output
  expect_equal(reward_fun_wd(episimdata, alpha, alpha_d, ovp, dovp, C_target, C_target_pen, D_target, D_target_pen, actions, ii, jj), expected_rew)

  # Test with a case exceeding the penalty threshold
  ii <- 3 # Using the third row (C = 200, Deaths = 20)
  C_err_pred <- abs(episimdata[ii, 'C'] - C_target) # abs(200 - 100) = 100
  D_err_pred <- abs(episimdata[ii, 'Deaths'] - D_target) # abs(20 - 10) = 10
  over_pen <- ovp # Over-penalty applied since C (200) > C_target_pen (150)
  over_pen_d <- dovp # Over-penalty applied since Deaths (20) > D_target_pen (15)
  expected_rew <- -alpha * C_err_pred - over_pen - actions[jj, "cost_of_NPI"] - alpha_d * D_err_pred - over_pen_d # -200 - 10 - 5 - 30 - 15 = -245

  expect_equal(reward_fun_wd(episimdata, alpha, alpha_d, ovp, dovp, C_target, C_target_pen, D_target, D_target_pen, actions, ii, jj), expected_rew)
})
