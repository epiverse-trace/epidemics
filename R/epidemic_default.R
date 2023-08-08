#' @title Model an SEIR-V epidemic with interventions
#'
#' @description Simulate an epidemic using a deterministic, compartmental
#' epidemic model with the compartments
#' "susceptible", "exposed", "infectious", "recovered", and "vaccinated".
#' This model can accommodate heterogeneity in social contacts among demographic
#' groups, as well as differences in the sizes of demographic groups.
#'
#' The `population` and an `infection` arguments are mandatory, while passing an
#' `intervention` and `vaccination` are optional and can be used to simulate
#' scenarios with different epidemic responses or different levels of the same
#' type of response.
#' See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param infection An `infection` object created using [infection()]. Must
#' have the basic reproductive number \eqn{R_0} of the infection, the
#' infectious period, and the pre-infectious period.
#' These are used to calculate the transmission rate \eqn{\beta}, the rate
#' at which individuals move from the 'exposed' to the 'infectious' compartment,
#' \eqn{\alpha}, and the recovery rate \eqn{\gamma}.
#' @param intervention An `<intervention>` object representing an optional
#' non-pharmaceutical intervention applied to the population during the
#' epidemic. See [intervention()] for details on constructing interventions with
#' age-specific effects on social contacts, as well as for guidance on how to
#' concatenate multiple overlapping interventions into a single `<intervention>`
#' object.
#' @param vaccination A `<vaccination>` object representing an optional
#' vaccination regime with a single dose, followed during the course of the
#' epidemic, with a start and end time, and age-specific vaccination rates.
#' See [vaccination()].
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 200 days.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details This is a wrapper function for [.epidemic_default_cpp()], a C++
#' function that uses Boost _odeint_ solvers for an SEIR-V model.
#'
#' This model only allows for single, population-wide rates of
#' transition between the 'susceptible' and 'exposed' compartments, between the
#' 'exposed' and 'infectious' compartments, and in the recovery rate.
#'
#' The internal C++ function [.epidemic_default_cpp()] accepts arguments that
#' are created by processing the `population`, `infection`, `intervention` and
#' `vaccination` arguments to the wrapper function into simpler forms.
#'
#' @return A `data.table` with the columns "time", "compartment", "age_group",
#' "value". The compartments correspond to the compartments of the model
#' chosen with `model`.
#' The current default model has the compartments "susceptible", "exposed",
#' "infectious", "recovered", and "vaccinated".
#' @export
epidemic_default_cpp <- function(population,
                                 infection,
                                 intervention = NULL,
                                 vaccination = NULL,
                                 time_end = 100,
                                 increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_class(infection, "infection")

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population, infection = infection,
    time_end = time_end, increment = increment
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_class(intervention, "intervention")
    model_arguments[["intervention"]] <- intervention
  }
  if (!is.null(vaccination)) {
    checkmate::assert_class(vaccination, "vaccination")
    model_arguments[["vaccination"]] <- vaccination
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_epidemic_default(
    .check_args_epidemic_default(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )

  # run model over arguments
  output <- do.call(.epidemic_default_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}
