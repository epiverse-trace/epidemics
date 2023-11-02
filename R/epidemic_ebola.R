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

#' @title Model an Ebola virus disease epidemic
#' @name epidemic_ebola
#' @rdname epidemic_ebola
#'
#' @description Simulate an epidemic using a discrete-time, stochastic SEIR
#' compartmental model with compartments based on Li et al. (2019), and with
#' Erlang passage times based on a model developed by Getz and Dougherty (2017),
#' developed to model the West African Ebola virus disease outbreak of 2014.
#' See **Details** for more information.
#'
#' `epidemic_ebola_cpp()` is an Rcpp implementation of this model that currently
#' lags behind the R implementation, and is likely to be removed.
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
#' have:
#'
#' - `r0`, a single number for the basic reproductive number \eqn{R_0} of the
#' infection,
#'
#' - `infectious_period`, a single number for the infectious period, taken to be
#' in days,
#'
#' - `preinfectious_period`, a single number for the pre-infectious period
#' before the onset of symptoms, taken to be in days,
#'
#' - `prop_community`, a single number in the range 0.0 -- 1.0 for the
#' proportion of infectious (assumed symptomatic) individuals who are not
#' hospitalised, and are instead infectious in the community,
#'
#' - `etu_risk`, a single number in the range 0.0 -- 1.0 for the relative
#' risk of EVD transmission in hospital settings (in Ebola Treatment Units; ETU)
#' between hospitalised individuals and susceptible ones.
#' This is understood to be a proportion of the baseline transmission rate
#' \eqn{\beta}, which is calculated from the `<infection>`.
#' 0.0 would indicate that hospitalisation completely eliminates transmission,
#' while 1.0 would indicate that transmission between hospitalised individuals
#' and susceptibles is the same as the baseline, which is transmission in the
#' community.
#'
#' - `funeral_risk`, a single number in the range 0.0 -- 1.0 for the relative
#' risk of funeral practices that may lead to transmission of EVD in a funeral
#' setting.
#' This is understood to be a proportion of the baseline transmission rate
#' \eqn{\beta}, which is calculated from the `<infection>`. It can alternatively
#' be interpreted as the proportion of funerals at which the risk of
#' transmission is the same as of infectious individuals in the community.
#'
#' `r0`, `infectious_period`, and `preinfectious_period` are used to calculate
#' the baseline transmission rate \eqn{\beta},
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
#'
#' @param time_dependence A named list where each name
#' is a model parameter (see `infection`), and each element is a function with
#' the first two arguments fixed as `time`, and `x`, followed by other arguments to be used by the supplied function. `time` is used internally to refer to the model time at which the named parameter will be changed.
#' See **Details** for more information.
#' @param time_end The maximum number of timesteps over which to run the model,
#' in days. Defaults to 100 days.
#' @return
#' A `<data.table>` in
#' long format with the columns "time", "compartment", "age_group", and "value",
#' that gives the number of individuals in each model compartment over time (
#' from 1 to `time_end`).
#' @details
#'
#' # Details: Discrete-time ebola virus disease model
#'
#' This model has compartments adopted from the consensus model for Ebola virus
#' disease presented in Li et al. (2019), and with transitions between
#' epidemiological compartments modelled using Erlang sub-compartments adapted
#' from Getz and Dougherty (2018); see **References**.
#'
#' The R code for this model is adapted from code by Hạ Minh Lâm and initially
#' made available on _Epirecipes_ (https://github.com/epirecipes/epicookbook)
#' under the MIT licence.
#'
#' The model implementation differs from Getz and Dougherty's (2018) in
#' allowing users to set the basic reproductive number \eqn{R_0}, and the mean
#' infectious (\eqn{\rho^I} in Getz and Dougherty) and pre-infectious periods in
#' days (\eqn{\rho^E} in Getz and Dougherty).
#' Getz and Dougherty instead calculate these from other model parameters.
#'
#' The shape of the Erlang distributions of passage times through the exposed
#' and infectious compartments (\eqn{k^E} and \eqn{k^I}) are fixed to 2 (this
#' was allowed to vary in Getz and Dougherty).
#'
#' The transition rates between the exposed and infectious, and infectious and
#' funeral compartments, \eqn{\gamma^E} and \eqn{\gamma^I} in Getz and
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
#' ## Hospitalisation, funerals, and removal
#'
#' Infectious individuals have a probability of 1.0 - `prop_community` of being
#' transferred to the hospitalised compartment, representing Ebola Treatment
#' Units (ETUs), and are considered to be infectious but no longer in the
#' community.
#' This compartment has the same number of sub-compartments as the infectious
#' compartment (i.e., infectious in the community), which means that
#' an infectious individual with \eqn{N} timesteps before exiting the
#' infectious compartment will exit the hospitalised compartment in the same
#' time.
#'
#' Hospitalised individuals can contribute to transmission of Ebola to
#' susceptibles depending on the value of `etu_risk` passed as part of the
#' `infection` argument, which scales the
#' baseline transmission rate \eqn{\beta} for hospitalised individuals.
#'
#' We assume that deaths in hospital lead to Ebola-safe funerals, and
#' individuals exiting the hospitalised compartment move to the 'removed'
#' compartment, which holds both recoveries and deaths.
#'
#' We assume that deaths outside of hospital lead to funerals that are
#' potentially unsafe burials, and the `funeral_risk` passed as part of the
#' `infection` argument scales the baseline transmission rate \eqn{\beta} for
#' funeral transmission of Ebola to susceptibles.
#'
#' Individuals are assumed to spend only a single timestep in the funeral
#' transmission compartment, before they move into the 'removed' compartment.
#'
#' Individuals in the 'removed' compartment do no affect model dynamics.
#'
#' ## Implementing vaccination
#'
#' Vaccination cannot currently be implemented in this model as it does not have
#' a "vaccinated" epidemiological compartment. This prevents the use of a
#' `<vaccination>` object.
#'
#' Instead, users can use the `time_dependence` argument to pass a function that
#' modifies model parameters --- specifically, the transmission rate --- in a
#' way that is consistent with the effect of vaccination.
#' An example is shown in the vignette about this model; run this code to open
#' the vignette: \code{vignette("ebola_model", package = "epidemics")}
#'
#' @references
#' Li, S.-L., Ferrari, M. J., Bjørnstad, O. N., Runge, M. C., Fonnesbeck, C. J.,
#' Tildesley, M. J., Pannell, D., & Shea, K. (2019). Concurrent assessment of
#' epidemiological and operational uncertainties for optimal outbreak control:
#' Ebola as a case study. Proceedings of the Royal Society B: Biological
#' Sciences, 286(1905), 20190774. \doi{10.1098/rspb.2019.0774}
#'
#' Getz, W. M., & Dougherty, E. R. (2018). Discrete stochastic analogs of Erlang
#' epidemic models. Journal of Biological Dynamics, 12(1), 16–38.
#' \doi{10.1080/17513758.2017.1401677}
#'
#' @export
epidemic_ebola_r <- function(population, infection,
                             intervention = NULL,
                             time_dependence = NULL,
                             time_end = 100) {
  # input checking for the ebola R model
  assert_population(
    population,
    compartments = c(
      "susceptible", "exposed", "infectious",
      "hospitalised", "funeral", "removed"
    )
  )
  assert_infection(
    infection,
    default_params = c(
      "name", "r0", "infectious_period", "preinfectious_period"
    ),
    extra_parameters = c(
      "prop_community", "etu_risk", "funeral_risk"
    ),
    extra_parameters_limits = list(
      prop_community = c(lower = 0, upper = 1),
      etu_risk = c(lower = 0, upper = 1),
      funeral_risk = c(lower = 0, upper = 1)
    )
  )
  if (!is.null(intervention)) {
    checkmate::assert_list(
      intervention,
      min.len = 1,
      names = "unique", any.missing = FALSE,
      types = "rate_intervention"
    )
  }
  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  if (!is.null(time_dependence)) {
    checkmate::assert_list(time_dependence, "function")
    invisible(
      lapply(time_dependence, checkmate::check_function,
        args = c("time", "x"),
        ordered = TRUE
      )
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
  names(initial_state) <- c(
    "susceptible", "exposed", "infectious",
    "hospitalised", "funeral", "removed"
  )

  # prepare output data.frame
  population_size <- sum(initial_state)
  sim_data <- matrix(NA_integer_, nrow = time_end, ncol = 6L)
  colnames(sim_data) <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "funeral", "removed"
  )

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

  # count number of infectious blocks or boxcars
  n_infectious_boxcars <- length(infectious_boxcar_rates)

  # initialise current conditions for exposed and infectious compartments
  exposed_current <- as.vector(
    stats::rmultinom(
      1,
      size = sim_data[1, "exposed"],
      prob = exposed_boxcar_rates
    )
  )
  exposed_past <- exposed_current

  infectious_current <- as.vector(
    stats::rmultinom(1,
      size = sim_data[1, "infectious"],
      prob = infectious_boxcar_rates
    )
  )
  infectious_past <- infectious_current

  hospitalised_current <- as.vector(
    stats::rmultinom(1,
      size = sim_data[1, "hospitalised"],
      prob = infectious_boxcar_rates
    )
  )
  hospitalised_past <- hospitalised_current

  funeral_trans_current <- sim_data[1, "funeral"]
  funeral_trans_past <- funeral_trans_current

  # define a fixed rounding factor for all timesteps to save function calls
  rounding_factor <- stats::rnorm(n_infectious_boxcars - 1, 0, 1e-2)

  # place parameter beta in list for rate interventions function
  parameters <- list(
    beta = beta,
    prop_community = get_parameter(infection, "prop_community"),
    etu_risk = get_parameter(infection, "etu_risk"),
    funeral_risk = get_parameter(infection, "funeral_risk")
  )

  ## Run the simulation from time t = 2 to t = time_end
  for (time in seq(2, time_end)) {
    # make a copy to assign time-dependent and intervention-affected values
    params <- parameters

    # apply time dependence before interventions
    time_dependent_params <- Map(
      parameters[names(time_dependence)],
      time_dependence,
      f = function(x, func) {
        func(time = time, x = x) # NOTE: time taken from loop index!
      }
    )
    # assign time-modified param values
    params[names(time_dependent_params)] <- time_dependent_params

    # check if an intervention is active and apply it to rates
    params <- intervention_on_rates(
      t = time,
      interventions = intervention[
        setdiff(names(intervention), "contacts")
      ],
      parameters = params
    )

    # transmission modifiers - 1.0 for baseline, user-provided for ETU risk and
    # funeral risk
    beta_modifiers <- c(1.0, params[["etu_risk"]], params[["funeral_risk"]])

    # get current transmission rate as base rate * intervention * p(infectious)
    # TODO: check if transmission rates should be summed or averaged
    transmission_rate <- sum(params[["beta"]] * beta_modifiers *
      sim_data[time - 1, c("infectious", "hospitalised", "funeral")]) /
      population_size
    exposure_prob <- 1.0 - exp(-transmission_rate)

    # calculate new exposures
    new_exposed <- stats::rbinom(
      1, sim_data[time - 1, "susceptible"], exposure_prob
    )
    # handle non-zero new exposures
    if (new_exposed > 0) {
      # distribute new exposures and add past exposures moved forward by one
      # timestep
      exposed_current <- as.vector(
        stats::rmultinom(1, size = new_exposed, prob = exposed_boxcar_rates)
      ) +
        c(exposed_past[-1], 0)
    } else {
      exposed_current <- c(exposed_past[-1], 0)
    }

    # handle hospitalisations first, as required for infectious compartment
    # new hospitalisations are a proportion of individuals from infectious
    # sub-compartments. Add a small normally distributed error to proportion
    # hospitalised to facilitate rounding to avoid fractional individuals
    # NOTE: proportion hospitalised = 1 - proportion community
    hospitalised_current <- round(
      # the SD of the normal distribution is small enough that values
      # added to zero lead to rounding to zero
      # first infectious_past compartment cannot be hospitalised and is
      # transferred to funeral compartment
      (infectious_past[-1] * (1.0 - params[["prop_community"]])) +
        rounding_factor
    )

    # calculate new infectious individuals
    new_infectious <- exposed_past[1]
    # handle non-zero new infectious
    if (new_infectious > 0) {
      infectious_current <- as.vector(
        stats::rmultinom(
          1,
          size = new_infectious, prob = infectious_boxcar_rates
        ) +
          c(infectious_past[-1] - hospitalised_current, 0)
      )
    } else {
      infectious_current <- c(infectious_past[-1] - hospitalised_current, 0)
    }

    # continue handling hospitalisations
    # concat zero to hospitalised_current at start as no infectious can go to
    # this compartment
    hospitalised_current <- c(0, hospitalised_current) +
      c(hospitalised_past[-1], 0)

    # calculate new individuals in the funeral transmission class
    funeral_trans_current <- infectious_past[1]

    # calculate new safely removed as the final hospitalised sub-compartments
    # and new burials of potentially transmitting funerals
    new_removed <- hospitalised_past[1] + funeral_trans_past

    # set past vectors to current vectors
    exposed_past <- exposed_current
    infectious_past <- infectious_current
    hospitalised_past <- hospitalised_current
    funeral_trans_past <- funeral_trans_current

    # prepare the data for output
    sim_data[time, "susceptible"] <- sim_data[time - 1, "susceptible"] -
      new_exposed
    sim_data[time, "exposed"] <- sum(exposed_current)
    sim_data[time, "infectious"] <- sum(infectious_current)
    sim_data[time, "hospitalised"] <- sum(hospitalised_current)
    sim_data[time, "funeral"] <- funeral_trans_current
    sim_data[time, "removed"] <- sim_data[time - 1, "removed"] +
      new_removed
  }

  # convert to long format
  output_to_df(
    output = list(x = sim_data, time = seq_len(time_end)),
    population = population,
    compartments = c(
      "susceptible", "exposed", "infectious",
      "hospitalised", "funeral", "removed"
    )
  )
}

#' @title Model a stochastic epidemic with Erlang passage times using Rcpp
#'
#' @name epidemic_ebola
#' @rdname epidemic_ebola
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
    "susceptible", "exposed", "infectious", "removed"
  )

  # run model over arguments
  output <- do.call(.epidemic_ebola_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}
