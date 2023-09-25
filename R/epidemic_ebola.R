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
#' @description Simulate an epidemic using a stochastic an SEIR compartmental
#' model with Erlang passage times based on Getz and Dougherty (2017) in
#' J. Biological Dynamics, and developed to model the West African
#' Ebola virus disease outbreak of 2014. See **Details** for more information.
#'
#' @param initial_state A vector that contains 4 numbers corresponding to the
#' initial values of the 4 classes: S, E, I, and R.
#' @param parameters A vector that contains 5 numbers corresponding to the
#' following parameters:
#'
#' 1. `shape_E`, a single integer value for the shape parameter of the Erlang
#' distribution from which passage times are drawn for the 'exposed'
#' compartment.
#'
#'  2. `rate_E`, a single integer value for the rate parameter of the Erlang
#' distribution from which passage times are drawn for the 'exposed'
#' compartment.
#'
#'  3. `shape_I`, a single integer value for the shape parameter of the Erlang
#' distribution from which passage times are drawn for the 'infectious'
#' compartment.
#'
#'  4. `rate_I`, a single integer value for the rate parameter of the Erlang
#' distribution from which passage times are drawn for the 'infectious'
#' compartment.
#'
#'  5. `beta`, a single number for the transmission rate of the infection.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 100 days.
#' @return
#' For `epidemic_ebola_r()`, a `<data.frame>` containing the numbers of
#' individuals in each compartment, "susceptible", "exposed", "infectious", and
#' "recovered", over time (from 1 to `time_end`).
#'
#' For `epidemic_ebola_cpp()`, a `<data.table>` in long format with the columns
#' "time", "compartment", "age_group", and "value", that gives the number of
#' individuals in each model compartment over time (from 0 to `time_end`).
#' @details
#' The R code for this model is taken from code by Hạ Minh Lâm and initially
#' made available on _Epirecipes_ (https://github.com/epirecipes/epicookbook)
#' under the MIT licence. The model is based on Getz and Dougherty (2017); see
#' **References**.
#' @references
#'
#' Getz, W. M., & Dougherty, E. R. (2018). Discrete stochastic analogs of Erlang
#' epidemic models. Journal of Biological Dynamics, 12(1), 16–38.
#' \doi{10.1080/17513758.2017.1401677}
#'
#' @export
epidemic_ebola_r <- function(initial_state, parameters, time_end = 100) {

  # input checking for the ebola R model
  checkmate::assert_integer(initial_state, lower = 0, any.missing = FALSE)
  checkmate::assert_integer(initial_state, lower = 0, any.missing = FALSE)
  names(initial_state) <- c("S", "E", "I", "R")
  names(parameters) <- c(
    "erlang_shape_for_E", "erlang_rate_for_E",
    "erlang_shape_for_I", "erlang_rate_for_I",
    "base_transmission_rate"
  )

  population_size <- sum(initial_state)
  sim_data <- data.frame(
    time = seq_len(time_end),
    S = NA, E = NA, I = NA, R = NA
  )
  sim_data[1, names(initial_state)] <- initial_state

  ## Initialise a matrix to store the states of the exposed sub-blocks
  # over time.
  exposed_block_adm_rates <- prob_discrete_erlang(
    shape = parameters["erlang_shape_for_E"],
    rate = parameters["erlang_rate_for_E"]
  )
  n_exposed_blocks <- length(exposed_block_adm_rates)
  exposed_blocks <- matrix(
    data = 0, nrow = time_end,
    ncol = n_exposed_blocks
  )
  exposed_blocks[1, n_exposed_blocks] <- sim_data$E[1]

  ## Initialise a matrix to store the states of the infectious
  # sub-blocks over time.
  infectious_block_adm_rates <- prob_discrete_erlang(
    shape = parameters["erlang_shape_for_I"],
    rate = parameters["erlang_rate_for_I"]
  )
  n_infectious_blocks <- length(infectious_block_adm_rates)
  infectious_blocks <- matrix(
    data = 0, nrow = time_end,
    ncol = n_infectious_blocks
  )
  infectious_blocks[1, n_infectious_blocks] <- sim_data$I[1]

  ## Run the simulation from time t = 2 to t = time_end
  for (time in seq(2, time_end)) {
    transmission_rate <-
      parameters["base_transmission_rate"] * sim_data$I[time - 1] /
        population_size
    exposure_prob <- 1 - exp(-transmission_rate)

    new_exposed <- stats::rbinom(1, sim_data$S[time - 1], exposure_prob)
    new_infectious <- exposed_blocks[time - 1, 1]
    new_recovered <- infectious_blocks[time - 1, 1]

    if (new_exposed > 0) {
      exposed_blocks[time, ] <- t(
        stats::rmultinom(1,
          size = new_exposed,
          prob = exposed_block_adm_rates
        )
      )
    }
    exposed_blocks[time, ] <-
      exposed_blocks[time, ] +
      c(exposed_blocks[time - 1, seq(2, n_exposed_blocks)], 0)

    if (new_infectious > 0) {
      infectious_blocks[time, ] <- t(
        stats::rmultinom(1,
          size = new_infectious,
          prob = infectious_block_adm_rates
        )
      )
    }
    infectious_blocks[time, ] <-
      infectious_blocks[time, ] +
      c(infectious_blocks[time - 1, seq(2, n_infectious_blocks)], 0)

    sim_data$S[time] <- sim_data$S[time - 1] - new_exposed
    sim_data$E[time] <- sum(exposed_blocks[time, ])
    sim_data$I[time] <- sum(infectious_blocks[time, ])
    sim_data$R[time] <- sim_data$R[time - 1] + new_recovered
  }

  return(sim_data)
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
