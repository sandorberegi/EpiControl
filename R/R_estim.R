#' Estimate Reproduction Number (R) and the Reduction Coefficient
#'
#' This function estimates the reproduction number (`R_est`) and an associated coefficient (`R_coeff_tmp`)
#' based on epidemic ata.
#'
#' @param episimdata A data frame containing epidemic simulation data with columns `C`, `Lambda_C`, and `R_coeff` containing, cases, total infectiousness (computed form cases), and the reduction in R by policy, respectively.
#' @param Ygen A numeric vector representing the generation time distribution.
#' @param ii An integer specifying the current time step in the simulation.
#' @param R_est_wind An integer (default = 5) defining the window size for R estimation.
#' @param r_dir An integer specifying the reproduction number adjustments.
#'
#' @return A list containing:
#'   - `R_est`: Estimated reproduction number.
#'   - `R_coeff_tmp`: Estimated reduction coefficient for `R0`.
#'
#' @export

R_estim <- function(episimdata, Ygen, ii, R_est_wind = 5, r_dir = 0){
  if ((ii-1) < R_est_wind) {
    R_est <- mean(episimdata[1:(ii-1), 'C'])/mean(episimdata[1:(ii-1), 'Lambda_C'])
    R_coeff_tmp <- sum(Ygen[1:(ii-1)] * episimdata[(ii-1):1, 'R_coeff'])/sum(Ygen[1:(ii-1)])
  } else {

    R_est <- mean(episimdata[(ii-R_est_wind):(ii-1), 'C'])/mean(episimdata[(ii-R_est_wind):(ii-1), 'Lambda_C'])

    if (r_dir == 1){
      R_coeff_tmp <-  mean(episimdata[(ii-R_est_wind):(ii-1), 'R_coeff'])
    } else {
      R_coeff_tmp <- sum(Ygen[1:(ii-R_est_wind)] * episimdata[(ii-R_est_wind):1, 'R_coeff'])/sum(Ygen[1:(ii-R_est_wind)])
    }
  }

  return(list(R_est = R_est, R_coeff_tmp = R_coeff_tmp))
}
