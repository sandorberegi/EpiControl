#' Vaccination Effect Over Time
#'
#' This function models the effect of vaccination over time using a scaled hyperbolic tangent function.
#'
#' @param x Numeric vector representing time or another continuous variable.
#' @param maxv Maximum vaccination effect or uptake level.
#' @param scale Scaling factor controlling the rate of change.
#' @param start The time point at which vaccination starts to take effect.
#'
#' @return A numeric vector representing the vaccination effect, ensuring non-negative values.
#'
#' @details The function applies a hyperbolic tangent transformation to model a smooth transition of vaccination uptake. The result is capped at a minimum of zero using `pmax()` to prevent negative values.
#'
#' @examples
#' x <- seq(0, 100, by = 1)
#' maxv <- 0.8
#' scale <- 10
#' start <- 50
#' vac_effect <- vac(x, maxv, scale, start)
#' plot(x, vac_effect, type = "l", main = "Vaccination Effect Over Time")
#'
#' @export
vac <- function(x, maxv, scale, start) {
  return(pmax(0, tanh((x-start) / scale) * maxv))
}
