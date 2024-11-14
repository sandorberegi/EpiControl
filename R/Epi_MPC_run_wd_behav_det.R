# Simulate the epidemic without control ('open-loop') pre-defined parameters.

Epi_MPC_run_wd_behav_det <- function(episimdata, epi_par, noise_par, actions, behav, pred_days, n_ens = 100, start_day = 1, ndays = nrow(episimdata), R_est_wind = 5, pathogen = 1, susceptibles = 1, delay = 0, ur = 0, r_dir = 1, N = 1e6) {

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

  Ygen <- dgamma(1:ndays, gen_time/gen_time_var, 1/gen_time_var)
  Ygen <- Ygen/sum(Ygen)

  Ydeaths <- dgamma(1:ndays, mortality_mean/mortality_var, 1/mortality_var)
  Ydeaths <- Ydeaths/sum(Ydeaths)

  number_of_actions <- nrow(actions)

  if (delay == 1) {
    repd_mean <- noise_par[1,'repd_mean'] #Reporting delay mean
    del_disp <- noise_par[1,'del_disp'] #Reporting delay variance

    Ydel <- dgamma(1:ndays, del_disp, del_disp/repd_mean)
    Ydel <- Ydel/sum(Ydel)
  }

  if (ur == 1) {
    ur_mean <- noise_par[1,'ur_mean'] #Under reporting, mean/variance
    ur_beta_a <- noise_par[1,'ur_beta_a']
    ur_beta_b <- (1-ur_mean)/ur_mean * ur_beta_a
  }

  for (ii in (start_day+1):ndays) {

    #estimate the reproduction number from data
    R_coeff_tmp <- 0.0

    if (ii-1 < R_est_wind) {
      episimdata[ii, 'Rest'] <- mean(episimdata[1:(ii-1), 'C'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
      R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff'])/sum(Ygen[1:(ii-1)])
    } else {

      episimdata[ii, 'Rest'] <- mean(episimdata[(ii-R_est_wind):(ii-1), 'C'])/mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])
      #R_coeff_tmp <- sum(Ygen[1:R_est_wind] * episimdata[(ii-1):(ii-R_est_wind), 'R_coeff'])/sum(Ygen[1:R_est_wind])

      if (r_dir == 1){
        R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
      } else {
        R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
      }
    }

    eta_tmp <- (1-R_coeff_tmp)

    c_tmp <- 1-(eta_tmp*episimdata[ii-1, 'p'])

    #print(c_tmp)

    #episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / R_coeff_tmp
    episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / c_tmp

    if (ii %% rf == 0L) {
      Rewards <- replicate(number_of_actions, 0)
      for (jj in 1:number_of_actions){
        Reward_ens <- replicate(n_ens ,0)
        for (kk in 1:n_ens){
          Reward_ens[kk] <- Epi_pred_wd_det(episimdata, epi_par, noise_par, actions, pathogen, pred_days, r_dir, ii, jj, N)
        }
        exp_reward <- mean(Reward_ens)
        Rewards[jj] <- exp_reward
      }
      #print(Rewards)
      episimdata[ii, 'policy'] <- which.max(Rewards)
    } else {
      episimdata[ii, 'policy'] <- episimdata[ii-1, 'policy']
    }

    Rcoeff <- actions[episimdata[ii, 'policy'], 'R_coeff']
    Rcoeff_cost <- actions[episimdata[ii, 'policy'], 'cost_of_NPI']
    episimdata[ii, 'R_coeff'] <- Rcoeff

    if (susceptibles == 1) {
      Ract <- Rcoeff * R0 * episimdata[ii-1,'S']/N
    } else {
      Ract <- Rcoeff * R0
    }

    u0S <- -(Ract/Rcoeff)*episimdata[ii-1, "Lambda_C"] - cost_diff*(episimdata[ii-1, "p"])^2
    u1S <- -Ract*episimdata[ii-1, "Lambda_C"] - cost_q*Rcoeff_cost - cost_diff*(episimdata[ii-1, "p"]-1)^2

    u0I <- -episimdata[ii-1, "p"]*(Ract/(sqrt(Rcoeff)))*episimdata[ii-1, "Lambda_C"]-(1-episimdata[ii-1, "p"])*(Ract/Rcoeff)*episimdata[ii-1, "Lambda_C"] - cost_diff*(episimdata[ii-1, "p"])^2
    u1I <- -episimdata[ii-1, "p"]*Ract*episimdata[ii-1, "Lambda_C"] - (1-episimdata[ii-1, "p"])*(Ract/(sqrt(Rcoeff)))*episimdata[ii-1, "Lambda_C"] - cost_q*Rcoeff_cost - cost_diff*(episimdata[ii-1, "p"]-1)^2

    u0 <- ind_impact*u0S + (1-ind_impact)*u0I
    u1 <- ind_impact*u1S + (1-ind_impact)*u1I

    BR <- 1/(1+exp(BR_kappa*(u0-u1)))

    #print(ii)
    #print(-(Ract/Rcoeff)*episimdata[ii-1, "Lambda_C"])
    #print(-Ract*episimdata[ii-1, "Lambda_C"])
    #print(cost_q*Rcoeff_cost)
    #print(episimdata[ii-1, "p"])
    #print(paste("u0:", u0))
    #print(paste("u1:", u1))

    episimdata[ii, "p"] <- bepsilon*(BR - episimdata[ii-1, "p"]) + episimdata[ii-1, "p"]

    eta <- (1 - Rcoeff)

    episimdata[ii, 'Re'] <- Ract/Rcoeff*(1-eta*episimdata[ii, "p"])
    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda_C'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])

    if (r_dir == 1) {
      pois_input <- episimdata[ii,'Re']*sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
    } else if ((r_dir == 2) && (ii > rf)) {
      Rdir <- logistic_function((ii %% rf), episimdata[(ii-(ii %% rf)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
      episimdata[ii, 'Rew'] <- Rdir
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'I']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }

    #print(pois_input)
    episimdata[ii,'I'] <- pois_input
    if (susceptibles == 1) {
      episimdata[ii, 'S'] <- episimdata[(ii-1), 'S'] - episimdata[ii,'I']
      if (episimdata[ii, 'S'] < 0) {
        episimdata[ii, 'S'] = 0
      }
    }

    pois_input_d <- CFR*sum(episimdata[(ii-1):1,'I']*Ydeaths[1:(ii-1)])

    episimdata[ii,'Deaths'] <- pois_input_d

    if (delay == 1) {
      pois_input_c <- sum(episimdata[ii:1,'I']* Ydel[1:ii])

      episimdata[ii,'C'] <- pois_input_c
    }

    if (ur == 1) {
      if (delay == 1){
        episimdata[ii,'C'] <- rbetabinom.ab(1, episimdata[ii,'C'], ur_beta_a, ur_beta_b)
      } else {
        episimdata[ii,'C'] <- rbetabinom.ab(1, episimdata[ii,'I'], ur_beta_a, ur_beta_b)
      }
    }

    if (delay+ur == 0){
      episimdata[ii,'C'] <- episimdata[ii,'I']
    }

  }
  return (episimdata)
}
