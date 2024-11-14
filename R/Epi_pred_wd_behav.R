# Simulate the epidemic without control ('open-loop') pre-defined parameters.

Epi_pred_wd_behav <- function(episimdata, epi_par, noise_par, actions, behav, pathogen, pred_days, r_dir, kk, jj, N, ndays = nrow(episimdata), pred_susceptibles = 0, gamma = 0.95) {

  R0 <- epi_par[pathogen,"R0"]
  gen_time <- epi_par[pathogen,"gen_time"]
  gen_time_var <- epi_par[pathogen,"gen_time_var"]

  CFR <- epi_par[pathogen,"CFR"]
  mortality_mean <- epi_par[pathogen,"mortality_mean"]
  mortality_var <- epi_par[pathogen,"mortality_var"]

  BR_kappa <- behav["BR_kappa"]
  bepsilon <- behav["epsilon"]
  cost_q <- behav["cost_q"]
  cost_diff <- behav["cost_diff"]
  ind_impact <-behav["ind_impact"]
  relu_coeff <-behav["relu_coeff"]

  Ygen <- dgamma(1:ndays, gen_time/gen_time_var, 1/gen_time_var)
  Ygen <- Ygen/sum(Ygen)

  Ydeaths <- dgamma(1:ndays, mortality_mean/mortality_var, 1/mortality_var)
  Ydeaths <- Ydeaths/sum(Ydeaths)

  pred_window_end <- kk+pred_days

  if (pred_window_end > ndays){
    pred_window_end <- ndays
  }

  rew <- replicate(ndays + pred_days, 0)
  discounts <- replicate(ndays + pred_days, 1)

  for (ii in kk:pred_window_end) {

    Rcoeff <- actions[jj,2]
    Rcoeff_cost <- actions[jj,4]

    R0act <- episimdata[kk,'R0est']*actions[jj,2]

    p_max <- 1/(1.1-Rcoeff) * 0.9

    if (pred_susceptibles == 1) {
      Ract <- R0act * episimdata[ii-1,'S']/N
    } else {
      Ract <- R0act
    }

    u0S <- -(Ract/Rcoeff)*episimdata[ii-1, "Lambda_C"] - cost_diff*(episimdata[ii-1, "p"])^2
    u1S <- -Ract*episimdata[ii-1, "Lambda_C"] - cost_q*Rcoeff_cost - cost_diff*(episimdata[ii-1, "p"]-1)^2

    u0I <- -episimdata[ii-1, "p"]*(Ract/(sqrt(Rcoeff)))*episimdata[ii-1, "Lambda_C"]-(1-episimdata[ii-1, "p"])*(Ract/Rcoeff)*episimdata[ii-1, "Lambda_C"] - cost_diff*(episimdata[ii-1, "p"])^2
    u1I <- -episimdata[ii-1, "p"]*Ract*episimdata[ii-1, "Lambda_C"] - (1-episimdata[ii-1, "p"])*(Ract/(sqrt(Rcoeff)))*episimdata[ii-1, "Lambda_C"] - cost_q*Rcoeff_cost - cost_diff*(episimdata[ii-1, "p"]-1)^2

    u0 <- ind_impact*u0S + (1-ind_impact)*u0I
    u1 <- ind_impact*u1S + (1-ind_impact)*u1I

    BR <- 1/(1+exp(BR_kappa*(u0-u1))) + neg_relu(relu_coeff*(u0-u1))

    BR <- pmin(BR, p_max)

    episimdata[ii, "p"] <- bepsilon*(BR - episimdata[ii-1, "p"]) + episimdata[ii-1, "p"]

    eta <- (1 - Rcoeff)

    episimdata[ii, 'Re'] <- Ract/Rcoeff*(1-eta*episimdata[ii, "p"])

    #episimdata[ii, 'Re'] <- Ract
    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda_C'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])

    if (r_dir == 1) {
      pois_input <- episimdata[ii,'Re']*sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
    } else if (r_dir == 2) {
      Rdir <- logistic_function((ii %% rf), episimdata[(ii-(ii %% rf)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'C']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }
    #print(pois_input)
    episimdata[ii,'C'] <- rpois(1, pois_input)

    if (pred_susceptibles == 1) {
      episimdata[ii, 'S'] <- episimdata[(ii-1), 'S'] - episimdata[ii,'C']
      if (episimdata[ii, 'S'] < 0) {
        episimdata[ii, 'S'] = 0
      }
    }

    pois_input_d <- CFR*sum(episimdata[(ii-1):1,'C']*Ydeaths[1:(ii-1)])
    episimdata[ii,'Deaths'] <- rpois(1, pois_input_d)

    discounts[ii] <- discounts[ii-1] * gamma
    rew[ii] <- reward_fun_wd(episimdata,alpha,alpha_d,ovp,dovp,C_target,C_target_pen,D_target,D_target_pen,R_target,actions,ii,jj)

  }

  Exp_rew <- sum(rew * discounts)

  return (Exp_rew)
}
