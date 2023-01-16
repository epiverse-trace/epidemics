
#' Make convenient parameters
#'
#' @description Function to make convenient parameters.
#' @export
make_parameters_sir_stochastic <- function() {
  list(
    beta = 1e-1,
    gamma = 0.05,
    N = 1000,
    S0 = 999,
    I0 = 1,
    R0 = 0,
    tf = 200
  )
}
