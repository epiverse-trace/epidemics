#' @title Model a diphtheria outbreak using a compartmental ODE model
#'
#' @name model_diphtheria
#' @rdname model_diphtheria
#'
#' @description Simulate a diphtheria outbreak using a deterministic,
#' compartmental ordinary differential equation model with the compartments
#' "susceptible", "exposed", "infectious", "hospitalised", and"recovered".
#' The model is based on Finger et al. (2019) and is intended to be used in the
#' context of internally displaced people (IDP) or refugee camps.
#' This model accommodates age or demographic structure and allows for a
#' proportion of each demographic group to be vaccinated and to not contribute
#' to the outbreak.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param transmissibility A single number for the rate at which individuals
#' move from the susceptible to the exposed compartment upon contact with an
#' infectious individual. Often denoted as \eqn{\beta}, with
#' \eqn{\beta = R_0 / \text{infectious period}}.
#' @param infectiousness_rate A single number for the rate at which individuals
#' move from the exposed to the infectious compartment. Often denoted as
#' \eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
#' This value does not depend upon the number of infectious individuals in the
#' population.
#' @param reporting_rate A single number for the proportion of infectious cases
#' that is reported; this is a precursor to hospitalisation as only reported
#' cases are hospitalised.
#' @param prop_hosp A single number for the proportion of reported cases that is
#' hospitalised.
#' @param hosp_entry_rate A single number for the rate at which reported cases
#' of infectious individuals are hospitalised.
#' This is calculated as 1 / time to hospitalisation, denoted \eqn{\tau_1}.
#' @param hosp_exit_rate A single number for the rate at which individuals are
#' discharged from hospital to enter the 'recovered' compartment.
#' This is calculated as 1 / time to discharge, denoted \eqn{\tau_2}.
#' @param recovery_rate A single number for the recovery rate, denoted
#' \eqn{\gamma}.
#' @param prop_vaccinated A numeric vector of the same length as the number of
#' demographic groups indicated the proportion of each group that is vaccinated.
#' These individuals are not included in the model dynamics.
#' @param intervention A named list of `<rate_intervention>` objects
#' representing optional pharmaceutical or non-pharmaceutical interventions
#' applied to the model parameters listed above.
#' @param time_dependence A named list where each name
#' is a model parameter, and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 100 days.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details
#'
#' ## R and Rcpp implementations
#'
#' `model_diphtheria_cpp()` is a wrapper function for [.model_diphtheria_cpp()],
#' an internal C++ function that uses Boost _odeint_ solvers for an SEIHR model.
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates transitions between
#' compartments. The default values are taken from Finger et al. (2019) where
#' possible.
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.8888889, assuming an
#' \eqn{R_0} of 4.0 and a total infectious period of 4.5 days.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.333, assuming
#' a pre-infectious period of 3 days.
#'
#' - Reporting rate (\eqn{r}, `reporting_rate`): 0.03, assuming that 3% of
#' infectious cases are detected or reported.
#'
#' - Proportion hospitalised (\eqn{\eta}, `prop_hosp`): 0.01, assuming that 1%
#' of reported cases need hospital treatment.
#'
#' - Hospital entry rate (\eqn{\tau_1}, `hosp_entry_rate`): 0.2, assuming that
#' it takes 5 days for infectious individuals to seek hospital treatment.
#'
#' - Hospital exit rate (\eqn{\tau_2}, `hosp_exit_rate`): 0.2, assuming that
#' individuals are discharged from hospital after 5 days.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.333, assuming an
#' infectious period following symptoms, of 3 days.
#'
#' @references
#' Finger, F., Funk, S., White, K., Siddiqui, M. R., Edmunds, W. J., &
#' Kucharski, A. J. (2019). Real-time analysis of the diphtheria outbreak in
#' forcibly displaced Myanmar nationals in Bangladesh. BMC Medicine, 17, 58.
#' \doi{10.1186/s12916-019-1288-7}.
#' @return A `<data.frame>` with the columns "time", "compartment", "age_group",
#' and "value".
#' @export
model_diphtheria_cpp <- function(population,
                                 transmissibility = 4.0 / 4.5,
                                 infectiousness_rate = 1.0 / 3.0,
                                 recovery_rate = 1.0 / 3.0,
                                 reporting_rate = 0.03,
                                 prop_hosp = 0.01,
                                 hosp_entry_rate = 0.2,
                                 hosp_exit_rate = 0.2,
                                 prop_vaccinated = 0.0 * get_parameter(
                                   population, "demography_vector"
                                 ),
                                 intervention = NULL,
                                 time_dependence = NULL,
                                 time_end = 100,
                                 increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(recovery_rate, lower = 0, finite = TRUE)
  # reporting rate and prop_hosp are expected to be proportions bounded 0 - 1
  checkmate::assert_number(reporting_rate, lower = 0, upper = 1)
  checkmate::assert_number(prop_hosp, lower = 0, upper = 1)
  # hospital entry and exit rate are very likely to be bounded 0 - 1 but
  # are allowed to be > 1 for now
  checkmate::assert_number(hosp_entry_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(hosp_exit_rate, lower = 0, finite = TRUE)

  # check that prop_vaccinated is the same length as demography_vector
  checkmate::assert_numeric(
    prop_vaccinated,
    lower = 0, upper = 1,
    len = length(get_parameter(population, "demography_vector"))
  )

  # allow only rate interventions in the intervention list
  checkmate::assert_list(
    intervention,
    min.len = 1, names = "unique", any.missing = FALSE,
    types = "rate_intervention", null.ok = TRUE
  )

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    names = "unique", any.missing = FALSE
  )
  invisible(
    lapply(time_dependence, checkmate::check_function,
      args = c("time", "x"),
      ordered = TRUE, null.ok = TRUE
    )
  )

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population,
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    reporting_rate = reporting_rate,
    prop_hosp = prop_hosp,
    hosp_entry_rate = hosp_entry_rate, hosp_exit_rate = hosp_exit_rate,
    prop_vaccinated = prop_vaccinated,
    intervention = intervention,
    time_dependence = time_dependence,
    time_end = time_end, increment = increment
  )

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_model_diphtheria(
    .check_args_model_diphtheria(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "recovered"
  )

  # run model over arguments
  output <- do.call(.model_diphtheria_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}
