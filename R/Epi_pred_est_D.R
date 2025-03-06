#' Projection of epidemic outcomes and calculates expected reward using Deaths only
#'
#' This function projects epidemic progression based on given parameters,
#' simulating infections and deaths over a prediction window. It incorporates
#' intervention effects and computes expected reward values using Deaths only.
#' Assumes infections are not seen/considered.
#'
#' @param episimdata Data frame containing epidemic simulation data.
#' @param epi_par Data frame of epidemiological parameters indexed by pathogen.
#' @param noise_par Data frame of noise parameters (not used in function but included for consistency).
#' @param actions Data frame of intervention actions with their corresponding effects.
#' @param pathogen An integer specifying the pathogen to extract corresponding epidemiological parameters.
#' @param pred_days Integer indicating the number of days to predict ahead.
#' @param r_dir Integer controlling the reproductive number adjustment method (1, 2, or other).
#' @param kk Integer indicating the current time step in the simulation.
#' @param jj Integer indexing the intervention scenario.
#' @param N Integer representing the total population.
#' @param ndays Integer specifying the total number of days in the simulation (default: number of rows in episimdata).
#' @param pred_susceptibles Integer flag (0 or 1) determining whether to adjust Re by the susceptible population.
#' @param gamma Discount factor for calculating expected reward.
#'
#' @return Numeric value representing the expected reward over the prediction window.
#'
#' @details
#' The function simulates epidemic spread using the gamma distribution to model
#' the generation time and mortality time. It applies interventions, updates
#' the effective reproduction number (Re), and predicts deaths via a Poisson process.
#' The expected reward is calculated based on \code{\link{reward_fun_wd}}.
#'
#' @examples
#' # Example usage (assuming required inputs are available):
#' result <- Epi_pred_est_D(episimdata, epi_par, noise_par, actions, "influenza",
#'                          pred_days = 14, r_dir = 1, kk = 50, jj = 2, N = 1e6)
#'
#' @export

Epi_pred_est_D <- function(episimdata, episettings, epi_par, noise_par, actions, pathogen, pred_days, r_dir, kk, jj, N, ndays = nrow(episimdata), pred_susceptibles = 0, gamma = 0.95) {

  sim_settings <- episettings$sim_settings

  rf <- sim_settings$rf
  t0 <- sim_settings$t0
  r_trans_steep <- sim_settings$r_trans_steep

  reward_function <- episettings$reward_function

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
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'Deaths']*Ygen[1:(ii-1)])

    if (r_dir == 1) {
      pois_input <- episimdata[ii,'Re']*sum(episimdata[(ii-1):1,'Deaths']*Ygen[1:(ii-1)])
    } else if (r_dir == 2) {
      Rdir <- logistic_function((ii %% rf), episimdata[(ii-(ii %% rf)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'Deaths']*Ygen[1:(ii-1)])
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'Deaths']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }
    #print(pois_input)
    if (is.na(pois_input)) {
      pois_input <- 0L
    }
    episimdata[ii,'Deaths'] <- rpois(1, pois_input)

    discounts[ii] <- discounts[ii-1] * gamma
    rew[ii] <- reward_function(episimdata,episettings,actions,ii,jj)
  }

  Exp_rew <- sum(rew * discounts)

  return (Exp_rew)
}
