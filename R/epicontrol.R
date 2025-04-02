#' Run Epidemiological Simulations with Control
#'
#' This function runs multiple ensemble simulations for epidemiological modelling,
#' incorporating noise parameters, estimated reproduction numbers, and intervention actions.
#' It supports parallel computation for efficiency.
#'
#' @param episim_data_ens A list of data frames containing initial epidemiological simulation data for each ensemble member.
#' @param episettings A list of settings controlling the simulation, including:
#'   - `sim_function`: Function to run a single epidemiological simulation.
#'   - `reward_function`: Function to calculate rewards for the control process.
#'   - `R_estimator`: Function to estimate the reproduction number.
#'   - `noise_par`: Parameters related to noise in the simulation.
#'   - `epi_par`: Epidemiological parameters.
#'   - `actions`: Possible intervention actions.
#'   - `sim_settings`: A list of simulation-specific settings. Please refer to the sim_function's documentation for the list of required parameters
#'   - `parallel`: Logical, whether to run simulations in parallel.
#'   - `cl`: Cluster object for parallel computation if `parallel = TRUE`.
#'
#' @return A list of simulation results for each ensemble member.
#'
#' @import pbapply
#' @export

epicontrol <- function(episim_data_ens, episettings) {
  sim_function <- episettings$sim_function
  reward_function <- episettings$reward_function
  R_estimator <- episettings$R_estimator
  noise_par <- episettings$noise_par
  epi_par <- episettings$epi_par
  actions <- episettings$actions
  sim_settings <- episettings$sim_settings
  parallel <- episettings$parallel

  sim_ens <- sim_settings$sim_ens
  ndays <- sim_settings$ndays

  print(sim_settings$susceptibles)

  pb <- txtProgressBar(min = 0, max = sim_ens, initial = 0, style = 3, width = 50, char = "=")
  for (ii in 1:sim_ens) {
    episim_data_ens[[ii]]$sim_id <- rep(ii, ndays)
  }


  if (parallel) {
    results <- pblapply(1:sim_ens, function(idx) {
      sim_function(episim_data_ens[[idx]], episettings, epi_par, noise_par, actions, pred_days = sim_settings$pred_days, n_ens = sim_settings$n_ens, start_day = sim_settings$start_day, ndays = sim_settings$ndays, R_est_wind = sim_settings$R_est_wind, pathogen = sim_settings$pathogen, susceptibles = sim_settings$susceptibles, delay = sim_settings$delay, ur = sim_settings$ur, r_dir = sim_settings$r_dir, N = sim_settings$N)
    }, cl = episettings$cl)

  }

  else {
    #episim_data_ens <- sim_function(episim_data_ens, episettings, epi_par, noise_par, actions, pred_days = sim_settings$pred_days, n_ens = sim_settings$n_ens, start_day = sim_settings$start_day, ndays = sim_settings$ndays, R_est_wind = sim_settings$R_est_wind, pathogen = sim_settings$pathogen, susceptibles = sim_settings$susceptibles, delay = sim_settings$delay, ur = sim_settings$ur, r_dir = sim_settings$r_dir, N = sim_settings$N)
    for (jj in 1:sim_ens) {
      episim_data_ens[[jj]] <- sim_function(episim_data_ens[[jj]], episettings, epi_par, noise_par, actions, pred_days = sim_settings$pred_days, n_ens = sim_settings$n_ens, start_day = sim_settings$start_day, ndays = sim_settings$ndays, R_est_wind = sim_settings$R_est_wind, pathogen = sim_settings$pathogen, susceptibles = sim_settings$susceptibles, delay = sim_settings$delay, ur = sim_settings$ur, r_dir = sim_settings$r_dir, N = sim_settings$N)
      setTxtProgressBar(pb,jj)
    }
    close(pb)
    results <- episim_data_ens

  }


  return(results)
}
