#' Specify no time dependence of model rates
#'
#' @return A list with a single named element, `transmissibility`, which is a
#' function that returns its second argument. This matches the specification for
#' functions to be passed as the `time_dependence` argument of epidemic
#' functions.
#' @keywords internal
no_time_dependence <- function() {
  # TODO: make internal
  list(
    transmissibility = function(time, x) x
  )
}
