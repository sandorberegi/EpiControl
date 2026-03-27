#' Set the model flow function
#'
#' Attach a user-defined flow function to a compartment model.
#'
#' @param model A `"cm_model"` object.
#' @param flow_fun A function describing model flows. It should accept
#'   arguments `t`, `state`, `pars`, and `dt`, and return a named vector of
#'   compartment changes. It may also include a named `"events"` attribute.
#'
#' @return The updated `"cm_model"` object.
#'
#' @examples
#' s <- new_compartment("S", 999)
#' i <- new_compartment("I", 1)
#' m <- build_model(s, i)
#'
#' flow_fun <- function(t, state, pars, dt) {
#'   c(S = -1, I = 1)
#' }
#'
#' m <- set_flows(m, flow_fun)
#'
#' @export
set_flows <- function(model, flow_fun) {
  stopifnot(is.function(flow_fun))
  model$flow_fun <- flow_fun
  model
}
