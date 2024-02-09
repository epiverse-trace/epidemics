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
#' @keywords internal
#' @noRd
#' @return A vector of variable length giving the probability of each integer
#' value for a cumulative probability of 0.99.
prob_discrete_erlang <- function(shape, rate) {
  n_bin <- 0
  factorials <- factorial(seq(0, shape))

  one_sub_cumulative_probs <- double(length = 1000L) # sufficient?
  cumulative_prob <- 0
  while (cumulative_prob <= 0.99) {
    n_bin <- n_bin + 1
    j <- seq_len(shape) - 1L
    total <- exp(-n_bin * rate) * ((n_bin * rate)^j) / factorials[j + 1]
    one_sub_cumulative_probs[n_bin] <- sum(total)
    cumulative_prob <- 1 - one_sub_cumulative_probs[n_bin]
  }
  one_sub_cumulative_probs <- c(1, one_sub_cumulative_probs[seq_len(n_bin)])

  density_prob <-
    utils::head(one_sub_cumulative_probs, -1) -
    utils::tail(one_sub_cumulative_probs, -1)
  density_prob <- density_prob / cumulative_prob

  return(density_prob)
}

#' @title Model an Ebola virus disease epidemic
#' @name model_ebola
#' @rdname model_ebola
#'
#' @description Simulate an epidemic using a discrete-time, stochastic SEIR
#' compartmental model with compartments based on Li et al. (2019), and with
#' Erlang passage times based on a model developed by Getz and Dougherty (2017),
#' developed to model the West African Ebola virus disease (EVD) outbreak of
#' 2013 -- 2016.
#' See **Details** for more information.
#'
#' `model_ebola_cpp()` is an Rcpp implementation of this model that currently
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
#' @param erlang_subcompartments The number of Erlang subcompartments assumed
#' for the exposed, infectious, and hospitalised compartments. Defaults to 2.
#' @inheritParams model_default
#' @param removal_rate The rate at which infectious individuals transition from
#' the infectious or hospitalised compartments to the funeral or removed
#' compartments. This model does not distinguish between recoveries and
#' deaths. Denoted in Getz and Dougherty as \eqn{\gamma^I} (see **Details**).
#' @param prop_community The proportion of infectious individuals who remain in
#' the community and are not hospitalised for treatment. Defaults to 0.5
#' @param etu_risk The relative risk of onward transmission of EVD from
#' hospitalised individuals. Must be a single value between 0.0 and 1.0, where
#' 0.0 indicates that hospitalisation completely prevents onward transmission,
#' and 1.0 indicates that hospitalisation does not prevent onward transmission
#' at all. `etu_risk` is used to scale the value of transmissibility for the
#' transmissibility \eqn{\beta}. Defaults to 0.2.
#' @param funeral_risk The relative risk of onward transmission of EVD from
#' funerals of individuals who died with EVD.
#' Must be a single value between 0.0 and 1.0, where
#' 0.0 indicates that there is no onward transmission, and 1.0 indicates that
#' funeral transmission is equivalent to transmission in the community.
#' `funeral_risk` is used to scale the value of transmissibility for the
#' transmissibility \eqn{\beta}. Defaults to 0.5.
#' @param intervention A named list of `<rate_intervention>` objects
#' representing optional pharmaceutical or non-pharmaceutical interventions
#' applied to the model parameters listed above.
#' @param time_dependence A named list where each name
#' is a model parameter, and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
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
#' The R code for this model is adapted from code by Ha Minh Lam and
#' initially made available on _Epirecipes_
#' (https://github.com/epirecipes/epicookbook) under the MIT license.
#'
#' The shape of the Erlang distributions of passage times through the exposed
#' and infectious compartments (\eqn{k^E} and \eqn{k^I}) are recommended to be
#' set to 2 as a sensible choice, which is the default value for the
#' `erlang_sbubcompartments` argument, but can be allowed to vary
#' (but not independently).
#'
#' The transition rates between the exposed and infectious, and infectious and
#' funeral compartments (and also hospitalised to removed),
#' \eqn{\gamma^E} and \eqn{\gamma^I} in Getz and Dougherty's notation, are
#' passed by the user as the `infectiousness_rate` and `removal_rate`
#' respectively.
#'
#' Getz and Dougherty's equation (6) gives the relationship between these
#' parameters and the mean pre-infectious \eqn{\rho^E} and infectious
#' \eqn{\rho^I} periods.
#' \deqn{\gamma^E = \dfrac{k^E}{\rho^E} = \dfrac{2}{\rho^E} ~\text{and}~
#' \gamma^I = \dfrac{k^I}{\rho^I} = \dfrac{2}{\rho^I}}
#'
#' In this discrete time model, \eqn{\gamma^E} and \eqn{\gamma^I} are used to
#' determine the passage times of newly exposed or infectious individuals
#' through their respective compartments (thus allowing for variation in passage
#' times).
#'
#' ## Hospitalisation, funerals, and removal
#'
#' A proportion, `1.0 - prop_community`, of infectious individuals are
#' transferred to the hospitalised compartment in each timestep,
#' This compartment represents Ebola Treatment Units (ETUs), and individuals
#' in the hospitalised compartment are considered to be infectious but no longer
#' in the community.
#'
#' The passage time of individuals in the hospitalised compartment is similar to
#' that of individuals in the infectious compartment (i.e., infectious in the
#' community), which means that an infectious individual with \eqn{N} timesteps
#' before exiting the infectious compartment will exit the hospitalised
#' compartment in the same time.
#'
#' Hospitalised individuals can contribute to transmission of Ebola to
#' susceptibles depending on the value of `etu_risk` which scales the
#' baseline transmission rate \eqn{\beta} for hospitalised individuals.
#'
#' We assume that deaths in hospital lead to Ebola-safe funerals, and
#' individuals exiting the hospitalised compartment move to the 'removed'
#' compartment, which holds both recoveries and deaths.
#'
#' We assume that deaths outside of hospital lead to funerals that are
#' potentially unsafe burials, and the `funeral_risk` argument scales the
#' baseline transmission rate \eqn{\beta} for funeral transmission of Ebola to
#' susceptibles.
#'
#' Individuals are assumed to spend only a single timestep in the funeral
#' transmission compartment, before they move into the 'removed' compartment.
#'
#' Individuals in the 'removed' compartment do no affect model dynamics.
#'
#' ## Model parameters
#'
#' The default values are:
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.125, resulting from
#' an \eqn{R_0} = 1.5 and an infectious period of 12 days.
#'
#' - Infectiousness rate (\eqn{\gamma^E}, `infectiousness_rate`): 0.4, assuming
#' a pre-infectious period of 5 days and two Erlang subcompartments.
#'
#' - Removal rate (\eqn{\gamma^I}, `recovery_rate`): 0.1667, assuming an
#' infectious period of 12 days and two Erlang subcompartments.
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
#' @examples
#' # create a population with 6 compartments
#' population <- population(
#'   contact_matrix = matrix(1),
#'   demography_vector = 14e6,
#'   initial_conditions = matrix(
#'     c(0.999998, 0.000001, 0.000001, 0, 0, 0),
#'     nrow = 1, ncol = 6L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- model_ebola_r(
#'   population = population
#' )
#'
#' # view some data
#' head(data)
#' @export
model_ebola_r <- function(population,
                          erlang_subcompartments = 2,
                          transmissibility = 1.5 / 12,
                          infectiousness_rate = erlang_subcompartments / 5,
                          removal_rate = erlang_subcompartments / 12,
                          prop_community = 0.9,
                          etu_risk = 0.7,
                          funeral_risk = 0.5,
                          intervention = NULL,
                          time_dependence = NULL, time_end = 100) {
  # input checking for the ebola R model - there is no dedicated checker fn
  # and input checking is performed here, making it different from other models
  # define compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious",
    "hospitalised", "funeral", "removed"
  )

  assert_population(
    population,
    compartments = compartments
  )

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_count(erlang_subcompartments, positive = TRUE)
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(removal_rate, lower = 0, finite = TRUE)
  # ratios are bounded 0 - 1
  checkmate::assert_number(prop_community, lower = 0, upper = 1)
  checkmate::assert_number(etu_risk, lower = 0, upper = 1)
  checkmate::assert_number(funeral_risk, lower = 0, upper = 1)

  # check time is a count
  checkmate::assert_count(time_end, positive = TRUE)

  # check all interventions are rate_interventions and target correct params
  if (!is.null(intervention)) {
    checkmate::assert_list(
      intervention,
      min.len = 1,
      names = "unique", any.missing = FALSE,
      types = "rate_intervention"
    )
    # check for model parameters targeted
    checkmate::assert_names(
      names(intervention),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "removal_rate",
        "prop_community", "etu_risk", "funeral_risk"
      )
    )
  }
  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  if (!is.null(time_dependence)) {
    checkmate::assert_list(
      time_dependence,
      types = "function",
      names = "unique", any.missing = FALSE
    )
    invisible(
      lapply(time_dependence, checkmate::check_function,
        args = c("time", "x"),
        ordered = TRUE
      )
    )
    # check for model parameters targeted
    checkmate::assert_names(
      names(time_dependence),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "removal_rate",
        "prop_community", "etu_risk", "funeral_risk"
      )
    )
  }

  # get initial conditions
  initial_state <- as.numeric(
    get_parameter(population, "initial_conditions")
  ) * get_parameter(population, "demography_vector")

  # round to nearest integer
  initial_state <- round(initial_state)
  names(initial_state) <- compartments

  # prepare output data.frame
  population_size <- sum(initial_state)
  sim_data <- matrix(NA_integer_, nrow = time_end, ncol = length(compartments))
  colnames(sim_data) <- compartments

  # assign initial conditions
  sim_data[1, ] <- initial_state

  # prepare probability vectors for which Erlang sub-compartment (boxcar)
  # will receive any newly exposed or infectious individuals
  exposed_boxcar_rates <- prob_discrete_erlang(
    shape = erlang_subcompartments,
    rate = infectiousness_rate
  )
  infectious_boxcar_rates <- prob_discrete_erlang(
    shape = erlang_subcompartments,
    rate = removal_rate
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

  # place parameter transmissibility in list for rate interventions function
  parameters <- list(
    transmissibility = transmissibility,
    prop_community = prop_community,
    etu_risk = etu_risk,
    funeral_risk = funeral_risk
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
    transmissibility_modifiers <- c(
      1.0, params[["etu_risk"]], params[["funeral_risk"]]
    )

    # get current transmissibility as base rate * intervention * p(infectious)
    # TODO: check if transmissibilities should be summed or averaged
    current_transmissibility <- sum(params[["transmissibility"]] *
      transmissibility_modifiers *
      sim_data[time - 1, c("infectious", "hospitalised", "funeral")]) /
      population_size
    exposure_prob <- 1.0 - exp(-current_transmissibility)

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
  data <- .output_to_df(
    output = list(x = sim_data, time = seq_len(time_end)),
    population = population,
    compartments = c(
      "susceptible", "exposed", "infectious",
      "hospitalised", "funeral", "removed"
    )
  )
  data.table::setDF(data)[]
}
