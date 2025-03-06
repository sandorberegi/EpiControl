#' Simulate Epidemic Dynamics with Model Predictive Control
#'
#' This function simulates an epidemic under model predictive control (MPC) using pre-defined parameters.
#' It dynamically evaluates actions to maximize expected rewards and updates epidemic states accordingly.
#'
#' @param episimdata A data frame containing simulation data. It should include columns such as
#' \code{"C"} (cases), \code{"I"} (infected individuals), \code{"Re"} (effective reproduction number),
#' \code{"S"} (susceptible individuals), \code{"Deaths"}, and \code{"Lambda"} total infectiousness.
#'@param epi_par A data frame containing epidemiological parameters for various pathogens. It should
#' have the following columns: \code{"R0"} (basic reproduction number), \code{"gen_time"}
#' (generation time), \code{"gen_time_var"} (variance of generation time), \code{"CFR"}
#' (case fatality rate), \code{"mortality_mean"} (mean mortality delay), and \code{"mortality_var"} (mortality delay variance).
#' @param noise_par A data frame containing noise parameters for under-reporting and reporting delays:
#' \itemize{
#'   \item \code{"repd_mean"}: Mean of the reporting delay.
#'   \item \code{"del_disp"}: Dispersion parameter for the reporting delay distribution.
#'   \item \code{"ur_mean"}: Mean for the under-reporting rate.
#'   \item \code{"ur_beta_a"}: Alpha parameter of the Beta distribution for under-reporting.
#' }
#' @param actions A data frame containing control actions. It should include columns like \code{"R_coeff"}.
#' @param pred_days An integer specifying the number of days to predict ahead in reward calculation.
#' @param n_ens An integer specifying the number of ensemble runs for expected reward calculation. Defaults to \code{100}.
#' @param start_day An integer specifying the start day of the simulation. Defaults to \code{1}.
#' @param ndays An integer specifying the total number of simulation days. Defaults to the number of rows in \code{episimdata}.
#' @param R_est_wind An integer specifying the window size for reproduction number estimation. Defaults to \code{5}.
#' @param pathogen An integer or string identifying the pathogen to extract corresponding parameters. Defaults to \code{1}.
#' @param susceptibles A binary value (\code{0} or \code{1}) indicating whether to update the number of susceptibles. Defaults to \code{1}.
#' @param delay A binary value (\code{0} or \code{1}) indicating whether to simulate reporting delays. Defaults to \code{0}.
#' @param ur A binary value (\code{0} or \code{1}) indicating whether to simulate under-reporting. Defaults to \code{0}.
#' @param r_dir An integer specifying the reproduction number adjustments:
#' \itemize{
#'   \item \code{1} for direct \code{Re}.
#'   \item \code{2} for logistic adjustments.
#'   \item \code{0} for using the generation time distribution.
#' }
#' @param N A numeric value representing the total population size. Defaults to \code{1e6}.
#'
#' @return A data frame containing updated simulation data with computed reproduction numbers,
#' estimated policies, daily infection incidents, cases, deaths, and other epidemic metrics.
#'
#' @details
#' The function employs a model predictive control strategy where actions are evaluated periodically
#' (based on a refresh interval, \code{rf}). Expected rewards are computed using the \code{\link{Epi_pred_wd}} function
#' over a specified prediction horizon, and the optimal action is selected. The simulation incorporates
#' effects of under-reporting, reporting delays, and susceptible depletion based on the provided parameters.
#'
#' @examples
#' # Example epidemiological data
#' episimdata <- data.frame(C = c(0, 10), I = c(5, 8), Re = c(1.5, NA), S = c(1000, 990), Deaths = c(0, 1))
#' epi_par <- data.frame(
#'   R0 = 2.5, gen_time = 5, gen_time_var = 1,
#'   CFR = 0.02, mortality_mean = 14, mortality_var = 2
#' )
#' noise_par <- data.frame(
#'   repd_mean = 2, del_disp = 1.5, ur_mean = 0.8, ur_beta_a = 2
#' )
#' actions <- data.frame(R_coeff = c(0.9, 1.1))
#' updated_data <- Epi_MPC_run_wd(
#'   episimdata = episimdata, epi_par = epi_par, noise_par = noise_par,
#'   actions = actions, pred_days = 10, n_ens = 50, start_day = 1,
#'   ndays = 20, R_est_wind = 5, pathogen = 1, susceptibles = 1,
#'   delay = 0, ur = 0, r_dir = 1, N = 1e6
#' )
#'
#' @export

Epi_MPC_run_wd <- function(episimdata, episettings, epi_par, noise_par, actions, pred_days, n_ens = 100, start_day = 1, ndays = nrow(episimdata), R_est_wind = 5, pathogen = 1, susceptibles = 1, delay = 0, ur = 0, r_dir = 1, N = 1e6) {

  R_estimator <- episettings$R_estimator
  sim_settings <- episettings$sim_settings

  rf <- sim_settings$rf
  r_trans_steep <- sim_settings$r_trans_steep
  t0 <- sim_settings$t0

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

    R_est_res <- R_estimator(episimdata, Ygen, ii, R_est_wind = R_est_wind, r_dir = r_dir)

    episimdata[ii, 'Rest'] <- R_est_res$R_est
    R_coeff_tmp <- R_est_res$R_coeff_tmp

#    if (ii-1 < R_est_wind) {
#      episimdata[ii, 'Rest'] <- mean(episimdata[1:(ii-1), 'C'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
#      R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff'])/sum(Ygen[1:(ii-1)])
#    } else {
#
#      episimdata[ii, 'Rest'] <- mean(episimdata[(ii-R_est_wind):(ii-1), 'C'])/mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])
#
#      if (r_dir == 1){
#        R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
#      } else {
#        R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
#      }
#    }

    episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / R_coeff_tmp

    if (ii %% rf == 0L) {
      Rewards <- replicate(number_of_actions, 0)
      for (jj in 1:number_of_actions){
        Reward_ens <- replicate(n_ens ,0)
        for (kk in 1:n_ens){
          Reward_ens[kk] <- Epi_pred_wd(episimdata, episettings, epi_par, noise_par, actions, pathogen, pred_days, r_dir, ii, jj, N)
        }
        exp_reward <- mean(Reward_ens)
        Rewards[jj] <- exp_reward
      }

      episimdata[ii, 'policy'] <- which.max(Rewards)
    } else {
      episimdata[ii, 'policy'] <- episimdata[ii-1, 'policy']
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
      Rdir <- logistic_function((ii %% rf), episimdata[(ii-(ii %% rf)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
      episimdata[ii, 'Rew'] <- Rdir
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'I']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }

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
        episimdata[ii,'C'] <- VGAM::rbetabinom.ab(1, episimdata[ii,'C'], ur_beta_a, ur_beta_b)
      } else {
        episimdata[ii,'C'] <- VGAM::rbetabinom.ab(1, episimdata[ii,'I'], ur_beta_a, ur_beta_b)
      }
    }

    if (delay+ur == 0){
      episimdata[ii,'C'] <- episimdata[ii,'I']
    }

  }
  return (episimdata)
}
