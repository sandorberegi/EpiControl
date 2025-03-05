#The main constructor function for the simulations
#Epi_MPC_run_wd <- function(episimdata, epi_par, noise_par, actions, pred_days, n_ens = 100, start_day = 1, ndays = nrow(episimdata), R_est_wind = 5, pathogen = 1, susceptibles = 1, delay = 0, ur = 0, r_dir = 1, N = 1e6) {
epicontrol <- function(episim_data_ens, episettings) {
  sim_function <- episettings$sim_function
  reward_function <- episettings$reward_function
  R_estimator <- episettings$R_estimator
  noise_par <- episettings$noise_par
  epi_par <- episettings$epi_par
  actions <- episettings$actions
  sim_settings <- episettings$sim_settings
  parallel <- episettings$parallel

  pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")
  for (ii in 1:sim_ens) {
    episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
  }

  if (parallel) {


    sim_ens <- sim_settings$sim_ens

    results <- pblapply(1:sim_ens, function(idx) {
      sim_function(episim_data_ens[[idx]], episettings, epi_par, noise_par, actions, pred_days = sim_settings$pred_days, n_ens = sim_settings$n_ens, start_day = sim_settings$start_day, ndays = sim_settings$ndays, R_est_wind = sim_settings$R_est_wind, pathogen = sim_settings$pathogen, susceptibles = sim_settings$susceptibles, delay = sim_settings$delay, ur = sim_settings$ur, r_dir = sim_settings$r_dir, N = sim_settings$N)
    }, cl = cl)

  }

  else {
    for (jj in 1:sim_ens) {
      episim_data_ens[[jj]] <- sim_function(episim_data_ens[[jj]], episettings, epi_par, noise_par, actions, pred_days = sim_settings$pred_days, n_ens = sim_settings$n_ens, start_day = sim_settings$start_day, ndays = sim_settings$ndays, R_est_wind = sim_settings$R_est_wind, pathogen = sim_settings$pathogen, susceptibles = sim_settings$susceptibles, delay = sim_settings$delay, ur = sim_settings$ur, r_dir = sim_settings$r_dir, N = sim_settings$N)
      setTxtProgressBar(pb,jj)
    }
    close(pb)
    results <- episim_data_ens

  }


  return(results)
}
