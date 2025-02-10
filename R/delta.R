#' Logistic Transition Function for new variant prevalence.
#'
#' Used in the COVID-19 Delta/Alpha variants example but can be used to model any generic new variant that becomes prevalent over the old one over time.
#'
#' @param x Numeric vector representing time.
#' @param scale Scaling factor controlling the rate of transition.
#' @param start The time point at which the new variant appears.
#'
#' @return A numeric vector representing the transition effect for new variant prevalence, ranging between 0 and 1.
#'
#' @details The function applies a hyperbolic tangent transformation to model a new more transmissible pathogen variant (e.g. for COVID-19 delta becoming more prevalent than alpha). The output is shifted and scaled to range between 0 and 1.
#'
#' @examples
#' x <- seq(0, 100, by = 1)
#' scale <- 10
#' start <- 50
#' delta_effect <- delta(x, scale, start)
#' plot(x, delta_effect, type = "l", main = "Logistic Transition Over Time")
#'
#' @export

delta <- function(x, scale, start) {
  return(tanh((x-start) / scale) * 0.5 + 0.5)
}
