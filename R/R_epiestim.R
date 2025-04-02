#' Estimate Reproduction Number (R) and Reduction Coefficient Using EpiEstim
#'
#' Estimates the effective reproduction number (`R_est`) and a temporary reduction coefficient (`R_coeff_tmp`)
#' based on epidemic simulation data, using the EpiEstim package.
#'
#' @param episimdata A data frame containing epidemic simulation data. It must include the following columns:
#'   - `C`: Number of cases.
#'   - `Lambda_C`: Total infectiousness (typically computed from past cases).
#'   - `R_coeff`: Reduction coefficient representing the effect of interventions on transmission.
#' @param Ygen A numeric vector representing the generation time distribution.
#' @param ii An integer indicating the current time step in the simulation.
#' @param R_est_wind An integer specifying the window size used for estimating `R` (default is 5).
#' @param r_dir An integer indicating how the reduction coefficient is averaged:
#'   - If `1`, uses a simple mean over the estimation window.
#'   - Otherwise, uses a weighted average based on the generation distribution.
#'
#' @return A list with two elements:
#'   - `R_est`: The estimated reproduction number.
#'   - `R_coeff_tmp`: The temporary estimate of the reduction coefficient.
#'
#' @export

R_epiestim <- function(episimdata, Ygen, ii, R_est_wind = 5, r_dir = 0){
  if ((ii-3) < R_est_wind) {
    R_est <- mean(episimdata[1:(ii-1), 'C'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
    R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff'])/sum(Ygen[1:(ii-1)])
  } else {

    Ygen0 <- c(0, Ygen)

    cases_v <- episimdata[1:(ii-1), 'C']

    t_start <- seq(2, (ii-1)-R_est_wind)
    t_end <- t_start + R_est_wind

    res <- EpiEstim::estimate_R(cases_v, method = "non_parametric_si",
                      config = EpiEstim::make_config(list(si_distr = Ygen0,
                                                          t_start = t_start,
                                                          t_end = t_end)))

    epiestimR <- res$R['Mean(R)']

    yepiestim_values <- epiestimR$'Mean(R)'

    R_est <- tail(yepiestim_values, 1)

    if (r_dir == 1){
      R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
    } else {
      R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
    }
  }

  return(list(R_est = R_est, R_coeff_tmp = R_coeff_tmp))
}
