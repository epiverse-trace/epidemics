#' Specify no time dependence of model rates
#'
#' @param func_target A string for the parameter to be targeted. Defaults to the
#' transmission rate parameter.
#' @return A list with a single named element, `transmission_rate`, which is a
#' function that returns its second argument. This matches the specification for
#' functions to be passed as the `time_dependence` argument of epidemic
#' functions.
#' @keywords internal
.no_time_dependence <- function(func_target = "transmission_rate") {
  l <- list(
    function(time, x) x
  )
  names(l) <- func_target
  l
}

#' Specify no population change during an epidemic
#'
#' @param population A `<population>` for which the population change list
#' should be suitable.
#' @keywords internal
.no_population_change <- function(population) {
  n_demo_groups <- length(population[["demography_vector"]])

  # return named list with 0 population change
  list(
    time = 0,
    values = list(rep(0, n_demo_groups))
  )
}
