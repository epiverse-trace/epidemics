
#' Make convenient parameters
#'
#' @param option Either "stochastic" or "deterministic".
#'
#' @description Function to make convenient parameters.
#' @examples 
#' make_parameters_sir(option = "stochastic")
#' make_parameters_sir(option = "deterministic")
#' @export
make_parameters_sir <- function(option = c("stochastic", "deterministic")) {
  option <- match.arg(arg = option, several.ok = FALSE)
  parameters <- switch(option,
    stochastic = list(
      beta = 1e-1,
      gamma = 0.05,
      N = 1000,
      S0 = 999,
      I0 = 1,
      R0 = 0,
      tf = 200
    ),
    deterministic = list(
      beta = 0.1,
      gamma = 0.05,
      S = 0.99,
      I = 0.01,
      R = 0.0,
      times = seq(0, 200, length.out = 2001)
    )
  )
  parameters
}
