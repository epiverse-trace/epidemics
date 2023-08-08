#' @title Model leaky, two-dose vaccination in an epidemic using Vacamole
#'
#' @description Simulate an epidemic using the _Vacamole_ model for Covid-19
#' developed at RIVM, the National Institute for Public Health and the
#' Environment in the Netherlands.
#' This model is aimed at estimating the impact of 'leaky' vaccination on an
#' epidemic. See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param infection An `infection` object created using [infection()]. Must
#' have the basic reproductive number \eqn{R_0} of the infection, the
#' infectious period, and the pre-infectious period.
#'
#' These are used to calculate the transmission rate \eqn{\beta}, the rate
#' at which individuals move from the 'exposed' to the 'infectious' compartment,
#' \eqn{\alpha}, and the recovery rate \eqn{\gamma}.
#'
#' This model also requires the `<infection>` object to have:
#'
#'  1. A single parameter `eta` for the hospitalisation rate of infectious
#' individuals,
#'
#'  2. A single parameter `omega` for the mortality rate of infectious or
#' hospitalised individuals,
#'
#'  3. A single parameter `susc_reduction_vax`, with a value between 0.0 and 1.0
#' giving the reduction in susceptibility to infection of individuals who have
#' received two doses of the vaccine,
#'
#'  4. A single parameter `hosp_reduction_vax`, with a value between 0.0 and 1.0
#' giving the reduction in hospitalisation rates of infectious individuals who
#' have received two doses of the vaccine,
#'
#'  5. A single parameter `mort_reduction_vax`, with a value between 0.0 and 1.0
#' giving the reduction in mortality of infectious and hospitalised individuals
#' who have received two doses of the vaccine.
#'
#' @param intervention An `<intervention>` object representing an optional
#' non-pharmaceutical intervention applied to the population during the
#' epidemic. See [intervention()] for details on constructing interventions with
#' age-specific effects on social contacts, as well as for guidance on how to
#' concatenate multiple overlapping interventions into a single `<intervention>`
#' object.
#' @param vaccination A `<vaccination>` object representing an optional
#' vaccination regime with two doses followed during the course of the
#' epidemic, with a start and end time, and age-specific vaccination rates for
#' each dose.
#' See [vaccination()].
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 200 days.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details This is a wrapper function for [.epidemic_vacamole_cpp()], a C++
#' function that uses Boost _odeint_ solvers for an SEIHR-V model.
#'
#' This model allows for:
#'
#'  1. A 'hospitalised' compartment along with hospitalisation rates;
#'
#'  2. Two doses of vaccination, with 'leaky' protection, i.e., vaccination does
#' not prevent infection completely but allows for a reduction in the infection
#' rate, as well as reduced rates of moving into states considered more serious,
#' such as 'hospitalised' or 'dead'.
#'
#' The internal C++ function [.epidemic_vacamole_cpp()] accepts arguments that
#' are created by processing the `population`, `infection`, `intervention` and
#' `vaccination` arguments to the wrapper function into simpler forms.
#'
#' @return A `data.table` with the columns "time", "compartment", "age_group",
#' "value". The compartments correspond to the compartments of the model
#' chosen with `model`.
#' The current default model has the compartments "susceptible",
#' "vaccinated_one_dose", "vaccinated_two_dose", "exposed",
#' "infectious", "infectious_vaccinated", "hospitalised",
#' "hospitalised_vaccinated", "recovered",  and "dead".
#' @export
epidemic_vacamole_cpp <- function(population,
                                  infection,
                                  intervention = NULL,
                                  vaccination,
                                  time_end = 100,
                                  increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_class(infection, "infection")
  checkmate::assert_class(vaccination, "vaccination")

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population, infection = infection, vaccination = vaccination,
    time_end = time_end, increment = increment
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_class(intervention, "intervention")
    model_arguments[["intervention"]] <- intervention
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_epidemic_vacamole(
    .check_args_epidemic_vacamole(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "vaccinated_one_dose",
    "vaccinated_two_dose", "exposed",
    "exposed_vaccinated", "infectious",
    "infectious_vaccinated", "hospitalised",
    "hospitalised_vaccinated", "dead",
    "recovered"
  )

  # run model over arguments
  output <- do.call(.epidemic_vacamole_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}
