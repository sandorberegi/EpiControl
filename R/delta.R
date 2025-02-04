#' Logistic Transition Function
#'
#' This function models a smooth logistic transition using a scaled hyperbolic tangent function.
#'
#' @param x Numeric vector representing time or another continuous variable.
#' @param scale Scaling factor controlling the rate of transition.
#' @param start The time point at which the transition begins.
#'
#' @return A numeric vector representing the transition effect, ranging between 0 and 1.
#'
#' @details The function applies a hyperbolic tangent transformation to model a smooth transition. The output is shifted and scaled to range between 0 and 1.
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
