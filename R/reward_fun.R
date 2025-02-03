#' Compute the Reward for Epidemic Simulation
#'
#' This function calculates the reward based on case prediction errors, reproduction number deviation, and intervention costs.
#'
#' @param episimdata A data frame containing epidemic simulation data.
#' @param alpha A penalty factor for case prediction errors.
#' @param ovp A penalty value applied if the predicted cases exceed a certain threshold.
#' @param C_target The target number of cases for the given day.
#' @param C_target_pen The case threshold above which the penalty `ovp` is applied.
#' @param R_target The target effective reproduction number.
#' @param actions A matrix containing intervention actions, with a column specifying the cost of non-pharmaceutical interventions (NPI).
#' @param ii The current time step in the simulation.
#' @param jj The index of the action scenario being evaluated.
#'
#' @return A numeric value representing the reward for the given time step.
#'
#' @details The function computes the absolute error between predicted and target cases (`C_err_pred`) and the deviation of the reproduction number (`R_err_pred`). If the predicted cases exceed `C_target_pen`, an additional penalty (`ovp`) is applied. The final reward incorporates these penalties along with the intervention cost.
#'
#' @export
reward_fun <- function(episimdata,alpha,ovp,C_target,C_target_pen,R_target,actions,ii,jj) {

  C_err_pred <- abs(episimdata[ii, 'C'] - C_target)
  R_err_pred <- abs(episimdata[ii, 'Re'] - R_target)
  over_pen <- 0.0

  if (episimdata[ii, 'C'] > C_target_pen) {
    over_pen <- ovp
  }

  rew <- -alpha * C_err_pred - over_pen - actions[jj,"cost_of_NPI"]

  return(rew)
}
