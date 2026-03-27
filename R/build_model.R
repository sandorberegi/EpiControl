#' Build a compartment model
#'
#' Construct a model object from one or more compartment objects.
#'
#' @param ... Compartments supplied individually, or a single list of
#'   compartments.
#'
#' @return An object of class `"cm_model"` containing named compartments and
#'   an initially `NULL` flow function.
#'
#' @examples
#' s <- new_compartment("S", 999)
#' i <- new_compartment("I", 1)
#' m <- build_model(s, i)
#' m
#'
#' @export
build_model <- function(...) {
  comps <- list(...)

  if (length(comps) == 1 && is.list(comps[[1]]) && !inherits(comps[[1]], "compartment")) {
    comps <- comps[[1]]
  }

  nms <- vapply(comps, `[[`, "", "Name")
  if (anyDuplicated(nms)) stop("Duplicate compartment names.")
  comps <- setNames(comps, nms)

  structure(
    list(
      compartments = comps,
      flow_fun = NULL
    ),
    class = "cm_model"
  )
}
