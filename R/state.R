#' Extract the current state vector
#'
#' Return the current compartment prevalences as a named numeric vector.
#'
#' @param model A `"cm_model"` object.
#'
#' @return A named numeric vector of compartment prevalences.
#'
#' @examples
#' s <- new_compartment("S", 999)
#' i <- new_compartment("I", 1)
#' m <- build_model(s, i)
#' state(m)
#'
#' @export
state <- function(model) {
  vapply(model$compartments, `[[`, numeric(1), "Prevalence")
}
