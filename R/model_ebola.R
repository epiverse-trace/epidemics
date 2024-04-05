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
#' @param transmission_rate A numeric vector for the rate at which individuals
#' move from the susceptible to the exposed compartment upon contact with an
#' infectious individual. Often denoted as \eqn{\beta}, with
#' \eqn{\beta = R_0 / \text{infectious period}}.
#' See **Details** for default values.
#' @param infectiousness_rate A numeric vector for the rate at which individuals
#' move from the exposed to the infectious compartment. Often denoted as
#' \eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
#' This value does not depend upon the number of infectious individuals in the
#' population.
#' See **Details** for default values.
#' @param erlang_subcompartments A numeric, integer-like vector for the number
#' of Erlang sub-compartments assumed for the exposed, infectious, and
#' hospitalised compartments. Defaults to 2.
#' @param removal_rate A numeric vector for the rate at which infectious
#' individuals transition from the infectious or hospitalised compartments to
#' the funeral or removed compartments.
#' This model does not distinguish between recoveries and deaths.
#' Denoted in Getz and Dougherty as \eqn{\gamma^I} (see **Details**).
#' @param prop_community A numeric vector for the proportion of infectious
#' individuals who remain in the community and are not hospitalised for
#' treatment. Defaults to 0.9.
#' @param etu_risk A numeric vector for the relative risk of onward transmission
#' of EVD from hospitalised individuals, with values between 0.0 and 1.0, where
#' 0.0 indicates that hospitalisation completely prevents onward transmission,
#' and 1.0 indicates that hospitalisation does not prevent onward transmission
#' at all; values are relative to the baseline transmission rate \eqn{\beta}.
#' Defaults to 0.7.
#' @param funeral_risk A numeric vector for the relative risk of onward
#' transmission of EVD from funerals of individuals who died with EVD.
#' Must be between 0.0 and 1.0, where 0.0 indicates that there is no onward
#' transmission, and 1.0 indicates that funeral transmission is equivalent to
#' the baseline transmission rate in the community \eqn{\beta}.
#' Defaults to 0.5.
#' @param intervention An optional named list of `<rate_intervention>` objects
#' representing optional pharmaceutical or non-pharmaceutical interventions
#' applied to the model parameters listed above. May also be a list of such
#' lists, in which case each set of interventions is treated as a separate
#' scenario. See **Details** below.
#' @param time_dependence An optional named list where each element is a
#' function with the first two arguments being the current simulation `time`,
#' and `x`, a value that is dependent on `time`
#' (`x` represents a model parameter).
#' List names must correspond to model parameters modified by the function.
#' Alternatively, may be a list of such lists, in which case each set of
#' functions is treated as a distinct scenario.
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end A numeric, integer-like vector for the maximum number of
#' timesteps over which to run the model, in days. Defaults to 100 days.
#' @return A `<data.table>`.
#' If the model parameters and composable elements are all scalars, a single
#' `<data.table>` with the columns "time", "compartment", "age_group", and
#' "value", giving the number of individuals per demographic group
#' in each compartment at each timestep in long (or "tidy") format is returned.
#'
#' If the model parameters or composable elements are lists or list-like,
#' a nested `<data.table>` is returned with a list column "data", which holds
#' the compartmental values described above.
#' Other columns hold parameters and composable elements relating to the model
#' run. Columns "scenario" and "param_set" identify combinations of composable
#' elements (population, interventions, vaccination regimes), and infection
#' parameters, respectively.
#' @details
#'
#' # Details: Discrete-time Ebola virus disease model
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
#' - Transmission rate (\eqn{\beta}, `transmission_rate`): 0.125, resulting from
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
#' ## Vector inputs
#'
#' ### Vector parameter inputs
#'
#' The model infection parameters and the model duration may be passed as
#' numeric or integer-like vectors (as appropriate to the parameter), to
#' simulate the effect of parameter uncertainty.
#' All parameter vectors must be of the same length, or any one parameter vector
#' may have a length > 1 while all other have a length of 1. In the first case,
#' each i-th combination of parameters is treated as a parameter set. In the
#' second case, all single value parameters (scalars) are recycled to the same
#' length as the non-scalar parameter.
#'
#' The model is run for $N$ stochastic realisations of each parameter set.
#' Random number seeds are preserved across parameter sets, so that differences
#' in outcomes in each j-th run are due to differences in parameters alone.
#'
#' ### Vector inputs for composable elements
#'
#' The `intervention` and `time_dependence` arguments also accept vectorised
#' inputs in the form of lists of intervention and time dependence sets.
#' Each combination of intervention and time-dependence sets is treated as a
#' distinct 'scenario', and realisations of each parameter set are run for each
#' scenario.
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
#' data <- model_ebola(
#'   population = population
#' )
#'
#' # view some data
#' head(data)
#' @export
model_ebola <- function(population,
                        transmission_rate = 1.5 / 12,
                        erlang_subcompartments = 2,
                        infectiousness_rate = erlang_subcompartments / 5,
                        removal_rate = erlang_subcompartments / 12,
                        prop_community = 0.9,
                        etu_risk = 0.7,
                        funeral_risk = 0.5,
                        intervention = NULL,
                        time_dependence = NULL, time_end = 100,
                        replicates = 100) {
  # input checking for the ebola R model - there is no dedicated checker fn
  # and input checking is performed here, making it different from other models
  # define compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "funeral", "removed"
  )
  assert_population(population, compartments)

  # NOTE: this relates to model structure but is allowed to be vectorised
  checkmate::assert_integerish(erlang_subcompartments, lower = 1)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(removal_rate, lower = 0, finite = TRUE)
  # ratios are bounded 0 - 1
  checkmate::assert_numeric(prop_community, lower = 0, upper = 1)
  checkmate::assert_numeric(etu_risk, lower = 0, upper = 1)
  checkmate::assert_numeric(funeral_risk, lower = 0, upper = 1)

  # check time is an integerish vector
  checkmate::assert_integerish(time_end, lower = 0)
  # replicates cannot vary
  checkmate::assert_count(replicates, positive = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    erlang_subcompartments = erlang_subcompartments,
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    removal_rate = removal_rate,
    prop_community = prop_community,
    etu_risk = etu_risk,
    funeral_risk = funeral_risk,
    time_end = time_end,
    replicates = replicates
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if `intervention` is a single intervention set or a list of such sets
  # NULL is allowed;
  # NOTE: only <rate_intervention> is allowed in this model
  is_lofints <- checkmate::test_list(
    intervention, "rate_intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) &&
    # Check that all elements of intervention sets are either
    # `<rate_intervention>` or NULL
    all(
      vapply(
        unlist(intervention, recursive = FALSE),
        FUN = function(x) {
          is_rate_intervention(x) || is.null(x)
        }, TRUE
      )
    )

  # Further input checking
  stopifnot(
    # Check if parameters can be recycled
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <rate_intervention> or a list of such
    lists" =
      is_lofints || is_lofls
  )

  # make lists if not lists
  population <- list(population) # NOTE: currently not list, but see issue #181
  if (is_lofints) {
    intervention <- list(intervention)
  }

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`, in order as the first two args
  # NOTE: this functionality is not vectorised;
  # convert to list for data.table list column
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    any.missing = FALSE, names = "unique"
  )
  # lapply on null returns an empty list
  invisible(
    lapply(time_dependence, checkmate::assert_function,
      args = c("time", "x"), ordered = TRUE
    )
  )
  # NOTE: `infectiousness_rate` and `removal_rate` are not allowed as they
  # control the number of epidemiological sub-compartments, which cannot be
  # safely changed while the model is running
  time_dependence <- list(
    .cross_check_timedep(
      time_dependence,
      c("transmission_rate", "prop_community", "etu_risk", "funeral_risk")
    )
  )

  # collect parameters and add a parameter set identifier
  params <- data.table::as.data.table(params)
  params[, "param_set" := .I]

  # this nested data.table will be returned
  model_output <- data.table::CJ(
    population = population,
    intervention = intervention,
    time_dependence = time_dependence,
    sorted = FALSE
  )

  # Send warning when > 10,000 runs (including replicates) are requested
  runs <- nrow(model_output) * nrow(params) * replicates
  warning_threshold <- 1e4
  if (runs > warning_threshold) {
    cli::cli_warn(
      sprintf(
        "Running %i scenarios and %i parameter sets with %i
        replicates each, for a total of %i model runs.",
        nrow(model_output), nrow(params), replicates, runs
      ),
      " This may take some time."
    )
  }

  # process the population, interventions, and vaccinations, after
  # cross-checking them agains the relevant population
  model_output[, args := apply(model_output, 1, function(x) {
    .check_prepare_args_ebola(c(x))
  })]
  model_output[, "scenario" := .I]

  # combine infection parameters and scenarios
  # NOTE: join X[Y] must have params as X as list cols not supported for X
  model_output <- params[, as.list(model_output), by = names(params)]

  # collect model arguments in column data, then overwrite
  model_output[, args := apply(model_output, 1, function(x) {
    c(x[["args"]], x[param_names]) # avoid including col "param_set"
  })]

  # call the internal function on all elements of the list args
  # NOTE: using internal seed preservation to ensure parameter sets use
  # identical random number streams
  # NOTE: call to `.model_ebola_internal()` wrapped in `.output_to_df_ebola()`
  # as separating them leads to unexpected list column structure when using `:=`
  model_output[, "data" := Map(args, population, f = function(l, p) {
    .output_to_df_ebola(
      output = do.call(.model_ebola_internal, l),
      population = p, # from local scope within data.table
      compartments = compartments
    )
  })]

  # remove temporary arguments
  model_output$args <- NULL

  # check for single row output, i.e., scalar arguments, and return data.table
  # do not return the parameters in this case
  if (nrow(model_output) == 1L) {
    model_output <- model_output[["data"]][[1L]] # hardcoded for special case
  }

  # return data.table
  model_output[]
}

#' @title Internal code for the Ebola model
#' @inheritParams model_ebola
#' @details
#' **NOTE** that all arguments to this internal function must be scalars for
#' infection parameters or single intervention or time-dependence sets for
#' composable elements.
#' @keywords internal
.model_ebola_internal <- function(
    initial_state, erlang_subcompartments, transmission_rate,
    infectiousness_rate, removal_rate, prop_community, etu_risk, funeral_risk,
    intervention, time_dependence, time_end, replicates) {
  withr::local_preserve_seed()
  # calculate required quantities
  population_size <- sum(initial_state)

  # prepare probability vectors for which epidemiological sub-compartment
  # will receive any newly exposed or infectious individuals
  # TODO: rename to be more informative; scheduled for follow-up PR
  exposed_boxcar_rates <- prob_discrete_erlang(
    shape = erlang_subcompartments,
    rate = infectiousness_rate
  )
  infectious_boxcar_rates <- prob_discrete_erlang(
    shape = erlang_subcompartments,
    rate = removal_rate
  )

  # compartment names required again
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "funeral", "removed"
  )

  # count number of infectious sub-compartments
  n_infectious_boxcars <- length(infectious_boxcar_rates)

  # place `transmission_rate` in list for rate_interventions function
  parameters <- list(
    transmission_rate = transmission_rate,
    prop_community = prop_community,
    etu_risk = etu_risk,
    funeral_risk = funeral_risk
  )

  # Use seed preservation to ensure that random number streams are preserved
  # for runs of intervention*parameter sets; i.e., run `i` of each set
  # should use the same random numbers
  # NOTE: Seed preservation could also be implemented one level up, but might be
  # more difficult to understand in context
  output_runs <- lapply(seq_len(replicates), function(xi_) {
    # NOTE: the original ebola model code continues here
    # Prepare data matrix and assign initial state
    sim_data <- matrix(
      NA_integer_,
      nrow = time_end, ncol = length(compartments)
    )
    colnames(sim_data) <- compartments
    sim_data[1, ] <- initial_state # at time = 1

    # Get a small random number for use with the hospitalisation functionality
    # save function calls by reusing the vector of values
    rounding_factor <- stats::rnorm(n_infectious_boxcars - 1, 0, 1e-2)

    # Distribute individuals into sub-compartments
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

    # NOTE: "funeral" has no sub-compartments
    funeral_trans_current <- sim_data[1, "funeral"]
    funeral_trans_past <- funeral_trans_current

    # Run the simulation from time t = 2 to t = time_end
    # A loop is required as conditions at each time t + 1 depend on time t.
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
      params <- .intervention_on_rates(
        t = time, interventions = intervention, parameters = params
      )

      # transmission modifiers - 1.0 for baseline, user-provided for ETU risk
      # and funeral risk
      transmission_rate_modifiers <- c(
        1.0, params[["etu_risk"]], params[["funeral_risk"]]
      )

      # get current transmission_rate as
      # base rate * intervention * p(infectious)
      # TODO: check if transmissibilities should be summed or averaged
      current_transmission_rate <- sum(params[["transmission_rate"]] *
        transmission_rate_modifiers *
        sim_data[time - 1, c("infectious", "hospitalised", "funeral")]) /
        population_size
      exposure_prob <- 1.0 - exp(-current_transmission_rate)

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
      # concat zero to hospitalised_current at start as no infectious can go
      # to this compartment
      hospitalised_current <- c(0, hospitalised_current) +
        c(hospitalised_past[-1], 0)

      # calculate new individuals in the funeral transmission class
      funeral_trans_current <- infectious_past[1]

      # calculate new safely removed as the final hospitalised
      # sub-compartments and new burials of potentially transmitting funerals
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

    # return simulated data matrix and time as a two element list
    # replicate id is handled in `.output_to_df_ebola()`
    list(x = sim_data, time = seq_len(time_end))
  })

  # return output runs
  output_runs
}
