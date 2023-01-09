
#' SIR epidemic model using 'deSolve'
#'
#' @param times A sequence of timepoints at which to solve for the proportion of
#' the population in each state (susceptible, infected, or recovered).
#' @param init The initial conditions of the population: a named vector with the
#' proportions of the population in each state ("S": susceptible, "I": infected,
#' or "R": recovered).
#' Can be easily constructed using the helper function [[initial_conditions()]].
#' @param parms The infection parameters in the form of a named vector with an
#' element "beta" for the rate of transmission \eqn{\beta}, and "gamma" for the
#' rate of recovery of infected individuals \eqn{\gamma}. Formally, both
#' \eqn{beta} and \eqn{gamma} are the rate parameters of separate exponential
#' distributions.
#'
#' @return A data table with the proportion of the population in the
#' susceptible, infected, and recovered states at each time point in `times`.
#' The data are returned in tidy format, and have \eqn{3T} rows, where \eqn{T}
#' is the number of timepoints.
#' @export
#'
#' @examples
#' library(epidemics)
#' sir_desolve(
#'   times = seq(0, 100, length.out = 100),
#'   init = initial_conditions(p_initially_infected = 0.01),
#'   parms = c(beta = 0.5, gamma = 0.05)
#' )
sir_desolve <- function(times = seq(0, 200, length.out = 2001),
                        init = initial_conditions(p_initially_infected = 0.01),
                        parms = c(beta = 0.1, gamma = 0.05)) {
  # define a function to pass to lsoda
  fn <- function(times, init, parms) {
    # define ODEs
    dS <- -parms[["beta"]] * init[["S"]] * init[["I"]]
    dI <- parms[["beta"]] * init[["S"]] * init[["I"]] - parms[["gamma"]] *
      init[["I"]]
    dR <- parms[["gamma"]] * init[["I"]]

    # return a list
    list(c(dS, dI, dR))
  }
  sir_out <- deSolve::lsoda(init, times = times, func = fn, parms = parms)
  sir_out <- data.table::as.data.table(sir_out)

  # get data in tidy format
  data.table::melt(
    sir_out,
    id.vars = "time", value.name = "proportion", variable.name = "state"
  )
}

#' Initial conditions for SIR models
#'
#' @param p_initially_infected The proportion initial infected.
#' @param p_already_recovered The proportion already recovered from infection.
#'
#' @return A named vector, ("S": susceptible, "I": infected, or "R": recovered).
#' @export
#'
#' @examples
#' initial_conditions(p_initially_infected = 0.001)
initial_conditions <- function(p_initially_infected = 0.01,
                               p_already_recovered = 0.0) {
  checkmate::assert_number(
    p_initially_infected,
    lower = 0.0, upper = 1.0
  )
  checkmate::assert_number(
    p_already_recovered,
    lower = 0.0, upper = 1.0
  )
  checkmate::assert_number(
    (p_initially_infected + p_already_recovered),
    lower = 0.0, upper = 1.0
  )
  c(
    S = 1.0 - (p_initially_infected + p_already_recovered),
    I = p_initially_infected,
    R = p_already_recovered
  )
}
