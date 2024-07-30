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
