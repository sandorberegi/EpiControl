#' Advance a model by one time step
#'
#' Update all compartment prevalences using the model's flow function.
#'
#' The flow function must return a named vector with entries corresponding to
#' all compartments in the model. Optionally, it may attach an `"events"`
#' attribute containing event counts. Event counts are coerced to integers in
#' the same way as compartment deltas.
#'
#' Non-finite values (`NA`, `NaN`, `Inf`, `-Inf`) in deltas or events are
#' replaced with zero before rounding.
#'
#' @param model A `"cm_model"` object.
#' @param t Current time.
#' @param dt Time-step size.
#' @param pars A list of model parameters.
#' @param event_names Optional character vector of event names, used if the
#'   `"events"` attribute returned by the flow function is unnamed.
#'
#' @return A list with components:
#' \describe{
#'   \item{model}{Updated `"cm_model"` object.}
#'   \item{states}{Named numeric vector of updated compartment prevalences.}
#'   \item{events}{Named integer vector of event counts.}
#' }
#'
#' @examples
#' s <- new_compartment("S", 100)
#' i <- new_compartment("I", 1)
#' m <- build_model(s, i)
#'
#' flow_fun <- function(t, state, pars, dt) {
#'   n_inf <- 2
#'   deltas <- c(S = -n_inf, I = n_inf)
#'   attr(deltas, "events") <- c(infections = n_inf)
#'   deltas
#' }
#'
#' m <- set_flows(m, flow_fun)
#' out <- step_model(m, t = 1, dt = 1)
#' out$states
#' out$events
#'
#' @export
step_model <- function(model, t, dt, pars = list(), event_names = NULL) {
  if (is.null(model$flow_fun)) stop("No flow_fun set. Use set_flows().")

  st <- state(model)
  deltas <- model$flow_fun(t = t, state = st, pars = pars, dt = dt)

  needed <- names(st)
  if (is.null(names(deltas)) || !all(needed %in% names(deltas))) {
    stop(
      "flow_fun must return a named vector with entries for: ",
      paste(needed, collapse = ", ")
    )
  }

  deltas2 <- as.numeric(deltas[needed])
  deltas2[!is.finite(deltas2)] <- 0
  deltas_int <- as.integer(round(deltas2))
  names(deltas_int) <- needed

  ev <- attr(deltas, "events")
  if (is.null(ev)) ev <- numeric(0)

  ev2 <- as.numeric(ev)
  if (length(ev2) > 0) {
    ev2[!is.finite(ev2)] <- 0
    ev_int <- as.integer(round(ev2))
  } else {
    ev_int <- integer(0)
  }

  if (length(ev_int) > 0) {
    if (is.null(names(ev)) || any(names(ev) == "")) {
      if (!is.null(event_names) && length(event_names) == length(ev_int)) {
        names(ev_int) <- event_names
      } else {
        stop(
          "Events vector has no names. Supply event_names=... to step_model(), ",
          "or ensure attr(deltas, 'events') is a named vector in flow_fun."
        )
      }
    } else {
      names(ev_int) <- names(ev)
    }
  }

  for (nm in needed) {
    d <- deltas_int[[nm]]
    new_val <- model$compartments[[nm]]$Prevalence + d
    model$compartments[[nm]]$Change <- d
    model$compartments[[nm]]$Prevalence <- max(0L, new_val)
  }

  list(
    model = model,
    states = state(model),
    events = ev_int
  )
}
