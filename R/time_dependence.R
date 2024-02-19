#' Specify no time dependence of model rates
#'
#' @param func_target A string for the parameter to be targeted. Defaults to the
#' transmissibility parameter.
#' @return A list with a single named element, `transmissibility`, which is a
#' function that returns its second argument. This matches the specification for
#' functions to be passed as the `time_dependence` argument of epidemic
#' functions.
#' @keywords internal
.no_time_dependence <- function(func_target = "transmissibility") {
  l <- list(
    function(time, x) x
  )
  names(l) <- func_target
  l
}
