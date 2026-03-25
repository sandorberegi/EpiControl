#' Create a compartment
#'
#' Construct a compartment object with a name, prevalence, and a placeholder
#' for change.
#'
#' @param name Character scalar. Name of the compartment.
#' @param prevalence Numeric scalar. Initial prevalence of the compartment.
#'
#' @return An object of class `"compartment"`.
#'
#' @examples
#' s <- new_compartment("S", 999)
#' i <- new_compartment("I", 1)
#' s
#'
#' @export
new_compartment <- function(name, prevalence) {
  stopifnot(
    is.character(name), length(name) == 1,
    is.numeric(prevalence), length(prevalence) == 1
  )

  structure(
    list(
      Name = name,
      Prevalence = prevalence,
      Change = NA_real_
    ),
    class = "compartment"
  )
}
