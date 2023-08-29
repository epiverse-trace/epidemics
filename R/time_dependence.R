#' Specify no time dependence of model rates
#'
#' @return A list with a single named element, `beta`, which is a function
#' that returns its second argument. This matches the specification for
#' functions to be passed as the `time_dependence` argument of epidemic
#' functions.
#' @export
no_time_dependence <- function() {
  list(
    beta = function(time, x) x
  )
}
