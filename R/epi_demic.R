
#' Run an SIR epidemic model.
#'
#' @param option A string, either "deterministic" or "stochastic", for which
#' kind of model to run.
#' @param parameters A parameter list created by the convenience function
#' [make_parameters_sir()].
#'
#' @return A data frame in tidy (long) format with the columns `time`, `S`, `I`
#' and `R` for the proportions of individuals infected at each time point.
#' @examples
#' epi_demic(
#'   option = "stochastic",
#'   parameters = make_parameters_sir(option = "stochastic")
#' )
#' epi_demic(
#'   option = "deterministic",
#'   parameters = make_parameters_sir(option = "deterministic")
#' )
#' @export
epi_demic <- function(option = c("stochastic", "deterministic"),
                      parameters = list()) {
  option <- match.arg(arg = option, several.ok = FALSE)
  data <- switch(option,
    deterministic = do.call(sir_desolve, parameters),
    stochastic = run_sir_stochastic(parameters)
  )
  epi_list_to_df(data)
}

#' Tidy SIR model output
#'
#' @param x The output of an SIR model. Must have the list elements, as vectors,
#' `time`, `S`, `I`, and `R`.
#' @keywords internal
#' @return A data frame in tidy (long) format with the columns `time`, `S`, `I`
#' and `R` for the proportions of individuals infected at each time point.
epi_list_to_df <- function(x) {
  classes <- rep(c("S", "I", "R"), each = length(x$time))
  data <- data.frame(
    time = rep(x$time, 3),
    variable = classes,
    value = c(x$S, x$I, x$R)
  )

  data
}
