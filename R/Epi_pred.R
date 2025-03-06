#' Projection of epidemic outcomes and calculating expected reward -- use when death data is not available
#'
#' #' This function projects epidemic progression based on given parameters,
#' simulating infections and deaths over a prediction window. It incorporates
#' intervention effects and computes expected reward values.
#'
#' @param episimdata A data frame containing epidemic simulation data.
#' @param epi_par A data frame of epidemiological parameters, with rows corresponding to different pathogens.
#' @param noise_par A parameter related to stochastic noise in the epidemic simulation (not used in the function).
#' @param actions A data frame of control actions, where the second column modifies the estimated reproduction number.
#' @param pathogen An integer specifying the pathogen to extract corresponding epidemiological parameters.
#' @param pred_days The number of days to predict forward.
#' @param kk The starting day of prediction.
#' @param jj The index of the action scenario being simulated.
#' @param N The total population size.
#' @param ndays The total number of days in the epidemic simulation (default: `nrow(episimdata)`).
#' @param pred_susceptibles Logical (0 or 1), indicating whether susceptible population dynamics should be accounted for.
#' @param gamma The discount factor for future rewards (default: 0.95).
#'
#' @return The expected reward (`Exp_rew`) over the simulation period.
#'
#' @details This function updates the effective reproduction number (`Re`), computes the expected number of new cases (`C`), and applies a discount factor to rewards computed using `reward_fun`. The function assumes that cases follow a Poisson process and uses a gamma-distributed generation time to estimate the infection dynamics.
#'
#' @importFrom stats rpois dgamma
#' @export
Epi_pred <- function(episimdata, episettings, epi_par, noise_par, actions, pathogen, pred_days, kk, jj, N, ndays = nrow(episimdata), pred_susceptibles = 0, gamma = 0.95) {

  sim_settings <- episettings$sim_settings

  rf <- sim_settings$rf
  t0 <- sim_settings$t0
  r_trans_steep <- sim_settings$r_trans_steep

  reward_function <- episettings$reward_function

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
    rew[ii] <- reward_function(episimdata,episettings,actions,ii,jj)

  }

  Exp_rew <- sum(rew * discounts)

  return (Exp_rew)
}
