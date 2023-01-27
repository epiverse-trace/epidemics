
#' SIR epidemic model using 'deSolve'
#'
#' @param times A sequence of timepoints at which to solve for the proportion of
#' @param beta Transmission rate \eqn{beta}.
#' @param gamma Recovery rate \eqn{gamma}.
#' @param S The proportion of individuals initially susceptible to infection.
#' @param I The proportion of individuals initially infected by the pathogen.
#' @param R The proportion of individuals already recovered from infection
#' which cannot be infected again.
#' @return A list with the proportion of the population in the
#' susceptible, infected, and recovered states at each time point in `times`.
#' @export
#'
#' @examples
#' data <- sir_desolve(
#'   beta = 0.1, gamma = 0.05, S = 0.99, I = 0.01, R = 0.0,
#'   times = seq(0, 100, 0.1)
#' )
#'
#' str(data)
sir_desolve <- function(beta, gamma, S, I, R, times) {

  # make a parameter list
  params <- c(beta = beta, gamma = gamma)
  # make a list of initial conditions
  init <- c(S = S, I = I, R = R)

  # define a function to pass to lsoda
  fn <- function(times, init, parms) {
    # define ODEs
    dS <- -(params[["beta"]]) * init[["S"]] * init[["I"]]
    dI <- params[["beta"]] * init[["S"]] * init[["I"]] - params[["gamma"]] *
      init[["I"]]
    dR <- params[["gamma"]] * init[["I"]]

    # return a list
    list(c(dS, dI, dR))
  }
  sir_out <- deSolve::lsoda(init, times = times, func = fn, parms = params)

  as.list(as.data.frame(sir_out))
}
