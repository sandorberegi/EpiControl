# Simulate the epidemic without control ('open-loop') pre-defined parameters.

Epi_pred <- function(episimdata, epi_par, noise_par, actions, pathogen, pred_days, kk, jj, N, ndays = nrow(episimdata), pred_susceptibles = 0, gamma = 0.95) {

  R0 <- epi_par[pathogen,"R0"]
  gen_time <- epi_par[pathogen,"gen_time"]
  gen_time_var <- epi_par[pathogen,"gen_time_var"]

  Ygen <- dgamma(1:ndays, gen_time/gen_time_var, 1/gen_time_var)
  Ygen <- Ygen/sum(Ygen)

  pred_window_end <- kk+pred_days

  if (pred_window_end > ndays){
    pred_window_end <- ndays
  }

  rew <- replicate(ndays + pred_days, 0)
  discounts <- replicate(ndays + pred_days, 1)

  for (ii in kk:pred_window_end) {

    R0act <- episimdata[kk,'R0est']*actions[jj,2]

    if (pred_susceptibles == 1) {
      Ract <- R0act * episimdata[ii-1,'S']/N
    } else {
      Ract <- R0act
    }

    episimdata[ii, 'Re'] <- Ract
    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])

    pois_input <- sum(episimdata[(ii-1):1,'C']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    #print(pois_input)
    episimdata[ii,'C'] <- rpois(1, pois_input)

    if (pred_susceptibles == 1) {
      episimdata[ii, 'S'] <- episimdata[(ii-1), 'S'] - episimdata[ii,'C']
      if (episimdata[ii, 'S'] < 0) {
        episimdata[ii, 'S'] = 0
      }
    }

    discounts[ii] <- discounts[ii-1] * gamma
    rew[ii] <- reward_fun(episimdata,alpha,ovp,C_target,C_target_pen,R_target,actions,ii,jj)

  }

  Exp_rew <- sum(rew * discounts)

  return (Exp_rew)
}
