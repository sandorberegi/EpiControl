#' Projection of epidemic outcomes and calculating expected reward
#'
#' This function projects epidemic progression based on given parameters,
#' simulating infections and deaths over a prediction window. It incorporates
#' intervention effects and computes expected reward values.
#'
#' @param episimdata A data frame containing simulation data. It should include columns such as
#' \code{"C"} (cases), \code{"I"} (infected individuals), \code{"Re"} (effective reproduction number),
#' \code{"S"} (susceptible individuals), \code{"Deaths"}, and \code{"Lambda"}.
#' @param epi_par A data frame containing epidemiological parameters for various pathogens. It should
#' have the following columns: \code{"R0"} (basic reproduction number), \code{"gen_time"}
#' (generation time), \code{"gen_time_var"} (variance of generation time), \code{"CFR"}
#' (case fatality rate), \code{"mortality_mean"}, and \code{"mortality_var"}.
#' @param noise_par A placeholder for surveillance noise parameters. Not used in projections.
#' @param actions A data frame containing control actions. Column 2 is expected to modify the effective
#' reproduction number (\code{"Re"}).
#' @param pathogen An integer specifying the pathogen to extract corresponding epidemiological parameters.
#' @param pred_days An integer specifying the number of days to predict ahead.
#' @param r_dir An integer specifying the reproduction number adjustments:
#' \itemize{
#'   \item \code{1} for direct \code{Re}.
#'   \item \code{2} for logistic adjustments.
#'   \item \code{0} for using the generation time distribution.
#' }
#' @param kk An integer indicating the starting day for prediction within the simulation.
#' @param jj An integer specifying the row index in \code{actions} to use for control effects.
#' @param N A numeric value representing the total population size.
#' @param ndays An integer specifying the total number of days in the simulation. Defaults to the number
#' of rows in \code{episimdata}.
#' @param pred_susceptibles A binary (0 or 1) indicating whether to update the number of susceptibles
#' during the simulation. Defaults to \code{0}.
#' @param gamma A numeric value between 0 and 1 representing the discount factor for future rewards. Defaults to \code{0.95}. Smaller values will prioritise
#' immediate rewards over longer term rewards.
#'
#' @return A numeric value representing the expected discounted reward over the prediction window.
#'
#' @details
#' The function simulates the epidemic using specified parameters and computes rewards for each
#' day within the prediction window.
#' Rewards are calculated using the \code{\link{reward_fun_wd}} function and are discounted
#' exponentially using the discount factor \code{gamma}.
#'
#' @examples
#' # Example epidemiological data
#' episimdata <- data.frame(R0est = c(1.5, 1.6), C = c(0, 10), Re = c(NA, NA), S = c(1000, 990), Deaths = c(0, 1))
#' epi_par <- data.frame(
#'   R0 = c(2.5), gen_time = c(5), gen_time_var = c(1),
#'   CFR = c(0.02), mortality_mean = c(14), mortality_var = c(2)
#' )
#' actions <- data.frame(action_effect = c(0.9, 0.8))
#' Epi_pred_wd(
#'   episimdata = episimdata, epi_par = epi_par, noise_par = NULL,
#'   actions = actions, pathogen = "pathogen1", pred_days = 10,
#'   r_dir = 1, kk = 2, jj = 1, N = 1000, ndays = 20
#' )
#'
#' @export

Epi_pred_wd <- function(episimdata, episettings, epi_par, noise_par, actions, pathogen, pred_days, r_dir, kk, jj, N, ndays = nrow(episimdata), pred_susceptibles = 0, gamma = 0.95) {

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
    episimdata[ii, 'Lambda'] <- sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])

    if (r_dir == 1) {
      pois_input <- episimdata[ii,'Re']*sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
    } else if (r_dir == 2) {
      Rdir <- logistic_function((ii %% rf), episimdata[(ii-(ii %% rf)-1),'Re'], episimdata[(ii),'Re'], r_trans_steep, t0)
      pois_input <- Rdir*sum(episimdata[(ii-1):1,'C']*Ygen[1:(ii-1)])
    } else {
      pois_input <- sum(episimdata[(ii-1):1,'C']*episimdata[ii:2,'Re']*Ygen[1:(ii-1)])
    }

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
    rew[ii] <- reward_function(episimdata,episettings,actions,ii,jj)

  }

  Exp_rew <- sum(rew * discounts)

  return (Exp_rew)
}
