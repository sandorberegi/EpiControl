# Simulate the epidemic without control ('open-loop') pre-defined parameters.

Epi_MPC_run_wd_thr <- function(episimdata, epi_par, noise_par, actions, pred_days, n_ens = 100, start_day = 1, ndays = nrow(episimdata), R_est_wind = 5, pathogen = 1, susceptibles = 1, delay = 0, ur = 0, r_dir = 0, N = 1e6) {

  last_changed <- 2L

  R0 <- epi_par[pathogen,"R0"]
  gen_time <- epi_par[pathogen,"gen_time"]
  gen_time_var <- epi_par[pathogen,"gen_time_var"]

  CFR <- epi_par[pathogen,"CFR"]
  mortality_mean <- epi_par[pathogen,"mortality_mean"]
  mortality_var <- epi_par[pathogen,"mortality_var"]

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

    episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / R_coeff_tmp

    if (ii %% rf == 0L) {
      if (episimdata[ii-1, 'policy'] == 1L) {
        if (episimdata[ii-1, 'Deaths'] > LD_on){
          episimdata[ii, 'policy'] <- 2L
          last_changed <- ii
        } else {
          episimdata[ii, 'policy'] <- episimdata[ii-1, 'policy']
        }
      } else if (episimdata[ii-1, 'policy'] == 2L) {
        if (episimdata[ii-1, 'Deaths'] < LD_off){
          episimdata[ii, 'policy'] <- 1L
          last_changed <- ii
        } else {
          episimdata[ii, 'policy'] <- episimdata[ii-1, 'policy']
        }
      } else {
        episimdata[ii, 'policy'] <- episimdata[ii-1, 'policy']
      }
    }

    Rcoeff <- actions[episimdata[ii, 'policy'], 'R_coeff']
    episimdata[ii, 'R_coeff'] <- Rcoeff

    if (susceptibles == 1) {
      Ract <- Rcoeff * R0 * episimdata[ii-1,'S']/N
    } else {
      Ract <- Rcoeff * R0
    }

    episimdata[ii, 'Re'] <- Ract
    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda_C'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])

    if (r_dir == 1) {
      pois_input <- episimdata[ii,'Re']*sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
    } else if ((r_dir == 2) && (ii > rf)) {
      Rdir <- logistic_function((ii-last_changed), episimdata[((last_changed)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
      episimdata[ii, 'Rew'] <- Rdir
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'I']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }

    #print(pois_input)
    episimdata[ii,'I'] <- rpois(1, pois_input)
    if (susceptibles == 1) {
      episimdata[ii, 'S'] <- episimdata[(ii-1), 'S'] - episimdata[ii,'I']
      if (episimdata[ii, 'S'] < 0) {
        episimdata[ii, 'S'] = 0
      }
    }

    pois_input_d <- CFR*sum(episimdata[(ii-1):1,'I']*Ydeaths[1:(ii-1)])

    episimdata[ii,'Deaths'] <- rpois(1, pois_input_d)

    if (delay == 1) {
      pois_input_c <- sum(episimdata[ii:1,'I']* Ydel[1:ii])

      episimdata[ii,'C'] <- rpois(1, pois_input_c)
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
