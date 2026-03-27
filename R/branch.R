#' Split a population into two branches
#'
#' Simulate a binomial branching process in which each individual is assigned
#' to branch A with probability `probA_over_probB`, and otherwise to branch B.
#'
#' @param population Integer or numeric scalar. Number of individuals to split.
#' @param probA_over_probB Numeric scalar in `[0, 1]`. Probability of assigning
#'   an individual to branch A.
#'
#' @return Integer vector of length 2 giving the counts in branches A and B.
#'
#' @examples
#' set.seed(1)
#' branch(100, 0.3)
#'
#' @export
branch <- function(population, probA_over_probB) {
  outcomeA <- rbinom(1, population, probA_over_probB)
  outcomeB <- population - outcomeA
  c(outcomeA, outcomeB)
}
