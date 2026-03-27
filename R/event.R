#' Draw a binomial event count
#'
#' Simulate the number of individuals experiencing an event.
#'
#' @param population Integer or numeric scalar. Number of individuals at risk.
#' @param probability Numeric scalar in `[0, 1]`. Event probability.
#'
#' @return Integer scalar giving the number of events.
#'
#' @examples
#' set.seed(1)
#' event(100, 0.2)
#'
#' @export
event <- function(population, probability) {
  rbinom(1, population, probability)
}
