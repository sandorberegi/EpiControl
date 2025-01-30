#' Logistic Function for Transition Between Reproduction Numbers
#'
#' This function models a smooth transition between two reproduction numbers (\eqn{R_0} and \eqn{R_1})
#' over time using a logistic curve.
#'
#' @param t A numeric value or vector representing the time point(s) at which to calculate the reproduction number.
#' @param R0 A numeric value representing the initial reproduction number (\eqn{R_0}).
#' @param R1 A numeric value representing the final reproduction number (\eqn{R_1}).
#' @param r A numeric value representing the growth rate or steepness of the logistic transition.
#' @param t0 A numeric value representing the midpoint of the transition where the reproduction number is halfway between \eqn{R_0} and \eqn{R_1}.
#'
#' @return A numeric value or vector representing the reproduction number at each time point in \code{t}.
#'
#' @details
#' The function computes the reproduction number using the logistic equation:
#' \deqn{R(t) = \frac{R_1 - R_0}{1 + \exp(-r \cdot (t - t_0))} + R_0}
#'
#' This models a smooth, sigmoidal transition from \eqn{R_0} to \eqn{R_1} over time.
#'
#' @examples
#' # Transition from R0 = 2.5 to R1 = 1.0 over time
#' times <- seq(0, 10, by = 0.1)
#' reproduction_numbers <- logistic_function(times, R0 = 2.5, R1 = 1.0, r = 0.5, t0 = 5)
#' plot(times, reproduction_numbers, type = "l", main = "Logistic Transition of Reproduction Number")
#'
#' @export

logistic_function <- function(t, R0, R1, r, t0) {
  K <- R1
  L <- R0
  (K - L) / (1 + exp(-r * (t - t0))) + L
}
