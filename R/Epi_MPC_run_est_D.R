#' Simulate Epidemic Dynamics with Open-loop Estimation Using Deaths
#'
#' This function simulates epidemic dynamics using a predefined set of parameters and
#' estimates reproduction numbers based on the rolling average of reported deaths.
#' It employs an open-loop strategy for policy evaluation.
#'
#' @param episimdata A data frame containing simulation data. It should include columns for:
#'   \itemize{
#'     \item \code{"I"}: Number of infected individuals.
#'     \item \code{"Deaths"}: Number of deaths.
#'     \item \code{"Lambda_C"}: Cumulative cases used for estimation.
#'     \item \code{"S"}: Number of susceptible individuals.
#'     \item \code{"R_coeff"}: Policy reproduction coefficients.
#'   }
#' @param epi_par A data frame of epidemiological parameters, including:
#'   \itemize{
#'     \item \code{"R0"}: Basic reproduction number.
#'     \item \code{"gen_time"}: Disease generation time.
#'     \item \code{"gen_time_var"}: Variance of the generation time.
#'     \item \code{"CFR"}: Case fatality rate.
#'     \item \code{"mortality_mean"}: Mean delay for mortality.
#'     \item \code{"mortality_var"}: Variance of mortality delay.
#'   }
#' @param noise_par A data frame containing noise parameters:
#'   \itemize{
#'     \item \code{"repd_mean"}: Reporting delay mean.
#'     \item \code{"del_disp"}: Dispersion parameter for reporting delays.
#'     \item \code{"ur_mean"}: Mean under-reporting rate.
#'     \item \code{"ur_beta_a"}: Alpha parameter of Beta distribution for under-reporting.
#'   }
#' @param actions A data frame containing policy actions with reproduction coefficients (\code{"R_coeff"}).
#' @param pred_days An integer specifying the number of days to predict ahead during policy evaluation.
#' @param n_ens An integer specifying the number of ensemble runs for Monte Carlo simulations. Defaults to \code{100}.
#' @param start_day An integer indicating the start day of the simulation. Defaults to \code{1}.
#' @param ndays An integer specifying the total number of simulation days. Defaults to the number of rows in \code{episimdata}.
#' @param R_est_wind An integer specifying the rolling window size for reproduction number estimation. Defaults to \code{5}.
#' @param pathogen An integer or string identifying the pathogen for parameter selection. Defaults to \code{1}.
#' @param susceptibles A binary value (\code{0} or \code{1}) indicating whether to simulate changes in susceptibles. Defaults to \code{1}.
#' @param delay A binary value (\code{0} or \code{1}) indicating whether to simulate reporting delays. Defaults to \code{0}.
#' @param ur A binary value (\code{0} or \code{1}) indicating whether to simulate under-reporting. Defaults to \code{0}.
#' @param r_dir An integer indicating reproduction number adjustments:
#' \itemize{
#' \item \code{1} for direct \code{Re}.
#' \item \code{2} for logistic adjustments.
#' \item \code{0} for using the generation time distribution.
#' @param N A numeric value representing the total population size. Defaults to \code{1e6}.
#'
#' @return A data frame containing updated simulation data with computed reproduction numbers,
#' estimated policies, rolling averages of deaths, and other epidemic metrics.
#'
#' @details
#' The function estimates reproduction numbers (\code{R0est}) based on deaths (\code{"Deaths"}) using
#' a rolling window approach. Policies are evaluated periodically using a Monte Carlo strategy,
#' and optimal policies are selected based on the highest expected rewards.
#'
#' @examples
#' # Example data and parameters
#' episimdata <- data.frame(I = c(10, 20), Deaths = c(1, 2), Lambda_C = c(10, 15), S = c(1000, 990), R_coeff = c(1.0, 0.9))
#' epi_par <- data.frame(
#'   R0 = 2.5, gen_time = 5, gen_time_var = 1, CFR = 0.02,
#'   mortality_mean = 14, mortality_var = 2
#' )
#' noise_par <- data.frame(
#'   repd_mean = 2, del_disp = 1.5, ur_mean = 0.8, ur_beta_a = 2
#' )
#' actions <- data.frame(R_coeff = c(1.0, 0.3))
#' results <- Epi_MPC_run_est_D(
#'   episimdata = episimdata, epi_par = epi_par, noise_par = noise_par,
#'   actions = actions, pred_days = 10, n_ens = 50, start_day = 1,
#'   ndays = 20, R_est_wind = 5, pathogen = 1, susceptibles = 1,
#'   delay = 0, ur = 0, r_dir = 1, N = 1e6
#' )
#'
#' @export

# Simulate the epidemic without control ('open-loop') pre-defined parameters.

Epi_MPC_run_est_D <- function(episimdata, epi_par, noise_par, actions, pred_days, n_ens = 100, start_day = 1, ndays = nrow(episimdata), R_est_wind = 5, pathogen = 1, susceptibles = 1, delay = 0, ur = 0, r_dir = 1, N = 1e6) {

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
      episimdata[ii, 'Rest'] <- mean(episimdata[1:(ii-1), 'Deaths'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
      R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'Deaths'])/sum(Ygen[1:(ii-1)])
    } else {

      episimdata[ii, 'Rest'] <- mean(episimdata[(ii-R_est_wind):(ii-1), 'Deaths'])/mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])
      #R_coeff_tmp <- sum(Ygen[1:R_est_wind] * episimdata[(ii-1):(ii-R_est_wind), 'R_coeff'])/sum(Ygen[1:R_est_wind])

      if (r_dir == 1){
        R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
      } else {
        R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
      }
    }

    episimdata[ii, 'R0est'] <- episimdata[ii, 'Rest'] / R_coeff_tmp

    if (is.infinite(episimdata[ii, 'R0est'])) {
      episimdata[ii, 'R0est'] <- 0
    }

    if (is.na(episimdata[ii, 'R0est'])) {
      episimdata[ii, 'R0est'] <- 0
    }

    #print(episimdata[ii, 'R0est'])

    if (ii %% rf == 0L) {
      Rewards <- replicate(number_of_actions, 0)
      for (jj in 1:number_of_actions){
        Reward_ens <- replicate(n_ens ,0)
        for (kk in 1:n_ens){
          Reward_ens[kk] <- Epi_pred_est_D(episimdata, epi_par, noise_par, actions, pathogen, pred_days, r_dir, ii, jj, N)
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
    episimdata[ii, 'R_coeff'] <- Rcoeff

    if (susceptibles == 1) {
      Ract <- Rcoeff * R0 * episimdata[ii-1,'S']/N
    } else {
      Ract <- Rcoeff * R0
    }

    episimdata[ii, 'Re'] <- Ract
    episimdata[ii, 'Rew'] <- sum(episimdata[ii:2,'Re']*Ygen[1:ii-1])/sum(Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'I']*Ygen[1:(ii-1)])
    episimdata[ii, 'Lambda_C'] <- sum(episimdata[(ii-1):1,'Deaths']*Ygen[1:(ii-1)])

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
