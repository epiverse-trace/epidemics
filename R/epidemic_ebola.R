#' @title Discrete probabilities for an Erlang distribution
#'
#' @name prob_discrete_erlang
#' @rdname prob_discrete_erlang
#'
#' @description A helper function that gives the probability of discrete values
#' from an Erlang distribution with a given shape and rate. The number of
#' values returned correspond to the number of discrete values over which the
#' cumulative probability reaches 0.99.
#' @param shape A single integer-like number for the shape of the Erlang
#' distribution.
#' @param rate A single number for the rate of the Erlang distribution.
#' @return A vector of variable length giving the probability of each integer
#' value for a cumulative probability of 0.99.
prob_discrete_erlang <- function(shape, rate) {
  n_bin <- 0
  factorials <- factorial(seq(0, shape))

  one_sub_cumulative_probs <- NULL
  cumulative_prob <- 0
  while (cumulative_prob <= 0.99) {
    n_bin <- n_bin + 1

    one_sub_cumulative_probs[n_bin] <- 0
    for (j in seq(0, (shape - 1))) {
      one_sub_cumulative_probs[n_bin] <-
        one_sub_cumulative_probs[n_bin] +
        (
          exp(-n_bin * rate) * ((n_bin * rate)^j) / factorials[j + 1]
        )
    }
    cumulative_prob <- 1 - one_sub_cumulative_probs[n_bin]
  }
  one_sub_cumulative_probs <- c(1, one_sub_cumulative_probs)

  density_prob <-
    utils::head(one_sub_cumulative_probs, -1) -
    utils::tail(one_sub_cumulative_probs, -1)
  density_prob <- density_prob / cumulative_prob

  return(density_prob)
}

#' @title Model a stochastic epidemic with Erlang passage times
#' @name epidemic_ebola
#' @rdname epidemic_ebola
#'
#' @description Simulate an epidemic using a discrete-time, stochastic SEIR
#' compartmental model with Erlang passage times based on Getz and Dougherty
#' (2017) in J. Biological Dynamics, and developed to model the West African
#' Ebola virus disease outbreak of 2014. See **Details** for more information.
#'
#' @param population An object of the `<population>` class, see [population()].
#'
#' This model only accepts a `<population>` without demographic structure, that
#' is, the `demography_vector` must be a single number representing the total
#' size of the affected population.
#'
#' The model also does not account for demographic differences in social
#' contacts, which means that the `contact_matrix` is ignored. For consistency,
#' the matrix must be square and have as many rows as demography groups, which
#' is one.
#'
#' @param infection An `<infection>` object created using [infection()]. Must
#' have the basic reproductive number \eqn{R_0} of the infection, the
#' infectious period, and the pre-infectious period.
#' These are used to calculate the baseline transmission rate \eqn{\beta},
#' as well as the rates \eqn{\gamma^E} and \eqn{\gamma^I} at which individuals
#' move from the 'exposed' to the 'infectious' compartment, and from the
#' 'infectious' to the 'recovered' compartment, respectively.
#' See **Details** for more information.
#'
#' @param intervention An optional `<rate_intervention>` object representing
#' pharmaceutical or non-pharmaceutical interventions applied to the infection's
#' parameters, such as the transmission rate, over the epidemic.
#' See [intervention()] for details on constructing rate interventions.
#' Defaults to `NULL`, representing no interventions on model parameters.
#' @param time_end The maximum number of timesteps over which to run the model,
#' in days. Defaults to 100 days.
#' @return
#' A `<data.table>` in
#' long format with the columns "time", "compartment", "age_group", and "value",
#' that gives the number of individuals in each model compartment over time (
#' from 1 to `time_end`).
#' @details
#'
#' ## Discrete-time ebola virus disease model following Getz & Dougherty (2017)
#'
#' The R code for this model is taken from code by Hạ Minh Lâm and initially
#' made available on _Epirecipes_ (https://github.com/epirecipes/epicookbook)
#' under the MIT licence. The model is based on Getz and Dougherty (2018); see
#' **References**.
#'
#' This model differs from Getz and Dougherty (2018) in allowing users to set
#' the basic reproductive number \eqn{R_0}, and the mean infectious
#' (\eqn{\rho^I} in Getz and Dougherty) and pre-infectious periods in days (
#' \eqn{\rho^E} in Getz and Dougherty). Getz and Dougherty instead calculate
#' these periods from other model parameters, and our change aims to make
#' the specification of the `<infection>` required for this model easier for
#' non-specialists.
#'
#' The shape of the Erlang distributions of passage times through the exposed
#' and infectious compartments (\eqn{k^E} and \eqn{k^I}) are fixed to 2 (this
#' was allowed to vary in Getz and Dougherty).
#'
#' The transition rates between the exposed and infectious, and infectious and
#' recovered compartments, \eqn{\gamma^E} and \eqn{\gamma^I} in Getz and
#' Dougherty's notation, are calculated following their equation (6).
#' \deqn{\gamma^E = \dfrac{k^E}{\rho^E} = \dfrac{2}{\rho^E} ~\text{and}~
#' \gamma^I = \dfrac{k^I}{\rho^I} = \dfrac{2}{\rho^I}}
#'
#' In this discrete time model, \eqn{\gamma^E} and \eqn{\gamma^I} are used to
#' determine the number of Erlang sub-compartments in each epidemiological
#' compartment, and the probability of newly exposed or infectious individuals
#' beginning in one of the compartments (thus allowing for variation in passage
#' times).
#'
#' @references
#'
#' Getz, W. M., & Dougherty, E. R. (2018). Discrete stochastic analogs of Erlang
#' epidemic models. Journal of Biological Dynamics, 12(1), 16–38.
#' \doi{10.1080/17513758.2017.1401677}
#'
#' @export
epidemic_ebola_r <- function(population, infection,
                             intervention = NULL, time_end = 100) {

  # input checking for the ebola R model
  assert_population(
    population,
    demography_groups = 1L,
    compartments = c("susceptible", "exposed", "infectious", "recovered")
  )
  assert_infection(
    infection,
    default_params = c(
      "r0", "infectious_period", "preinfectious_period"
    )
  )
  if (!is.null(intervention)) {
    assert_intervention(
      intervention,
      type = "rate", population = population
    )
  }

  # set Erlang shape parameters, k^E and k^I; this is a modelling decision
  shape_E <- 2L
  shape_I <- 2L

  # get Erlang rate parameters
  rate_E <- shape_E / get_parameter(infection, "preinfectious_period")
  rate_I <- shape_I / get_parameter(infection, "infectious_period")

  # prepare base transmission rate beta
  beta <- get_transmission_rate(infection = infection)

  # get initial conditions
  initial_state <- as.numeric(
    get_parameter(population, "initial_conditions")
  ) * get_parameter(population, "demography_vector")

  # round to nearest integer
  initial_state <- round(initial_state)
  names(initial_state) <- c("susceptible", "exposed", "infectious", "recovered")

  # prepare output data.frame
  population_size <- sum(initial_state)
  sim_data <- matrix(NA_integer_, nrow = time_end, ncol = 4L)
  colnames(sim_data) <- c("susceptible", "exposed", "infectious", "recovered")

  # assign initial conditions
  sim_data[1, ] <- initial_state

  # prepare probability vectors for which Erlang sub-compartment (boxcar)
  # will receive any newly exposed or infectious individuals
  exposed_boxcar_rates <- prob_discrete_erlang(
    shape = shape_E,
    rate = rate_E
  )
  infectious_boxcar_rates <- prob_discrete_erlang(
    shape = shape_I,
    rate = rate_I
  )

  # count number of exposed and infectious blocks or boxcars
  n_exposed_boxcars <- length(exposed_boxcar_rates)
  n_infectious_boxcars <- length(infectious_boxcar_rates)

  # prepare the current and past compartments
  exposed_current <- numeric(n_exposed_boxcars)
  exposed_past <- exposed_current

  infectious_current <- numeric(n_infectious_boxcars)
  infectious_past <- infectious_current

  # initialise current conditions for exposed and infectious compartments
  exposed_current[n_exposed_boxcars] <- sim_data[1, "exposed"]
  infectious_current[n_infectious_boxcars] <- sim_data[1, "infectious"]

  ## Run the simulation from time t = 2 to t = time_end
  for (time in seq(2, time_end)) {
    # get current transmission rate as base rate * p(infectious)
    transmission_rate <- beta * sim_data[time - 1, "infectious"] /
      population_size
    exposure_prob <- 1.0 - exp(-transmission_rate)

    # calculate new exposures, infectious, and recovered
    new_exposed <- stats::rbinom(
      1, sim_data[time - 1, "susceptible"], exposure_prob
    )
    new_infectious <- exposed_past[1]
    new_recovered <- infectious_past[1]

    # handle non-zero new exposures
    if (new_exposed > 0) {
      exposed_current <- as.vector(
        stats::rmultinom(1, size = new_exposed, prob = exposed_boxcar_rates)
      )
    }
    # add new exposures to seq(2, last) past boxcar compartments
    exposed_current <- exposed_current + c(exposed_past[-1], 0)

    # handle non-zero new infectious
    if (new_infectious > 0) {
      infectious_current <- as.vector(
        stats::rmultinom(
          1,
          size = new_infectious, prob = infectious_boxcar_rates
        )
      )
    }
    # add new infectious to seq(2, last) past boxcar compartments
    infectious_current <- infectious_current + c(infectious_past[-1], 0)

    # set past vectors to current vectors
    exposed_past <- exposed_current
    infectious_past <- infectious_current

    # prepare the data for output
    sim_data[time, "susceptible"] <- sim_data[time - 1, "susceptible"] -
      new_exposed
    sim_data[time, "exposed"] <- sum(exposed_current)
    sim_data[time, "infectious"] <- sum(infectious_current)
    sim_data[time, "recovered"] <- sim_data[time - 1, "recovered"] +
      new_recovered
  }

  # convert to long format
  output_to_df(
    output = list(x = sim_data, time = seq_len(time_end)),
    population = population,
    compartments = c("susceptible", "exposed", "infectious", "recovered")
  )
}

#' @title Model a stochastic epidemic with Erlang passage times using Rcpp
#' @rdname epidemic_ebola
#'
#' @param population An object of the `<population>` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param infection An `<infection>` object created using [infection()]. Must
#' have the basic reproductive number \eqn{R_0} of the infection, and the
#' infectious period.
#' These are used to calculate the transmission rate \eqn{\beta}, the rate
#' at which individuals move from the 'exposed' to the 'infectious' compartment.
#' This differs from how \eqn{\beta} is passed directly as a parameter in
#' `epidemic_ebola_r()`.
#'
#' The `<infection>` object must hold the parameters of the Erlang distributions
#' of passage times through the "exposed" and "infectious" compartments:
#' `shape_E`, `rate_E`, `shape_I`, and `rate_I`.
#'
#' @details `epidemic_ebola_cpp()` is a wrapper function for
#' [.epidemic_ebola_cpp()], a C++ function that is exposed to R via _Rcpp_.
#' The internal C++ function [.epidemic_ebola_cpp()] accepts arguments that
#' are created by processing the `population`, `infection`, `intervention` and
#' `vaccination` arguments to the wrapper function into simpler forms. This
#' processing is performed internally.
#'
#' @export
epidemic_ebola_cpp <- function(population, infection,
                               time_end = 100) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_class(infection, "infection")

  # check the time end
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population, infection = infection,
    time_end = time_end
  )

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_epidemic_ebola(
    .check_args_epidemic_ebola(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered"
  )

  # run model over arguments
  output <- do.call(.epidemic_ebola_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}
