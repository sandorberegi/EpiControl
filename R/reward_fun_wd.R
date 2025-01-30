#' Calculate Reward function for Model Predictive Optimal Control
#'
#' This function computes a reward value based on epidemiological simulation data,
#' penalties for exceeding target thresholds, and the cost of non-pharmaceutical interventions (NPIs).
#'
#' @param episimdata A data frame containing epidemiological simulation data. The function expects
#' specific columns: \code{"C"} for cases, \code{"Re"} for reproduction number, and \code{"Deaths"}.
#' @param alpha A numeric value representing the weight applied to case error in the reward calculation.
#' @param alpha_d A numeric value representing the weight applied to death error in the reward calculation.
#' @param ovp A numeric value for the penalty applied when cases exceed \code{C_target_pen}.
#' @param dovp A numeric value for the penalty applied when deaths exceed \code{D_target_pen}.
#' @param C_target A numeric target for the number of cases.
#' @param C_target_pen A numeric threshold for the penalty on cases.
#' @param D_target A numeric target for the number of deaths.
#' @param D_target_pen A numeric threshold for the penalty on deaths.
#' @param actions A data frame containing the costs of non-pharmaceutical interventions (NPIs).
#' It expects a column named \code{"cost_of_NPI"}.
#' @param ii An integer specifying the row index in \code{episimdata} to use for reward calculation.
#' @param jj An integer specifying the row index in \code{actions} to use for cost retrieval.
#'
#' @return A numeric value representing the calculated reward.
#'
#' @details
#' The reward is computed as:
#' \deqn{reward = -\alpha * |C - C_{target}| - penalty_{C} - NPI_{cost} - \alpha_d * |Deaths - D_{target}| - penalty_{D}}
#'
#' Where:
#' \itemize{
#'   \item \eqn{penalty_{C}} is applied if cases exceed \code{C_target_pen}.
#'   \item \eqn{penalty_{D}} is applied if deaths exceed \code{D_target_pen}.
#'   \item \code{NPI_{cost}} is retrieved from the \code{actions} data frame based on index \code{jj}.
#' }
#'
#' @examples
#' reward_fun_wd(episimdata, alpha = 0.5, alpha_d = 0.7, ovp = 50, dovp = 30,
#'               C_target = 120, C_target_pen = 140, D_target = 7, D_target_pen = 12,
#'               actions = actions, ii = 1, jj = 2)
#'
#' @export

reward_fun_wd <- function(episimdata,alpha,alpha_d,ovp,dovp,C_target,C_target_pen,D_target,D_target_pen,actions,ii,jj) {

  C_err_pred <- abs(episimdata[ii, 'C'] - C_target)
  R_err_pred <- abs(episimdata[ii, 'Re'] - R_target)
  D_err_pred <- abs(episimdata[ii, 'Deaths'] - D_target)
  over_pen <- 0.0
  over_pend_d <- 0.0

  if (episimdata[ii, 'C'] > C_target_pen) {
    over_pen <- ovp
  }

  if (episimdata[ii, 'Deaths'] > D_target_pen) {
    over_pen_d <- dovp
  }

  rew <- -alpha * C_err_pred - over_pen - actions[jj,"cost_of_NPI"] -alpha_d * D_err_pred - over_pend_d

  return(rew)
}
