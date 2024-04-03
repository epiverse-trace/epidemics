#' @title Internal code for the Ebola model
#' @keywords internal
.model_ebola_internal <- function(
    initial_state, erlang_subcompartments, transmission_rate,
    infectiousness_rate, removal_rate, prop_community, etu_risk, funeral_risk,
    intervention, time_dependence, time_end, replicates) {
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

  # Use seed preservation to ensure that random number streams are preserved
  # for runs of intervention*parameter sets; i.e., run `i` of each set
  # should use the same random numbers
  # NOTE: Seed preservation could also be implemented one level up, but might be
  # more difficult to understand in context
  output_runs <- withr::with_preserve_seed({
    lapply(seq_len(replicates), function(x) {
      # NOTE: the original ebola model code continues here, but is simply
      # wrapped in seed preservation code.

      # Get a small random number for use with the hospitalisation functionality
      # save function calls by reusing the vector of values
      n_infectious_boxcars <- length(infectious_boxcar_rates)
      rounding_factor <- stats::rnorm(n_infectious_boxcars - 1, 0, 1e-2)

      # Prepare data matrix and assign initial state
      sim_data <- matrix(
        NA_integer_,
        nrow = time_end, ncol = length(compartments)
      )
      colnames(sim_data) <- compartments
      sim_data[1, ] <- initial_state # at time = 1

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

      # place `transmission_rate` in list for rate_interventions function
      parameters <- list(
        transmission_rate = transmission_rate,
        prop_community = prop_community,
        etu_risk = etu_risk,
        funeral_risk = funeral_risk
      )

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
        params <- intervention_on_rates(
          t = time,
          interventions = intervention,
          parameters = params
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
  })

  # return replicates data
  output_runs
}

#' @title Replacement ebola model
#' @export
model_ebola2 <- function(population,
                         erlang_subcompartments = 2,
                         transmission_rate = 1.5 / 12,
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

  # NOTE: this relates to model structure and is not currently vectorised
  # TODO: check if this should be hard-coded and not an argument
  checkmate::assert_count(erlang_subcompartments, positive = TRUE)
  # replicates cannot vary
  checkmate::assert_count(replicates, positive = TRUE)

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
    "`intervention` must be a list of <intervention>s or a list of such lists" =
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
  model_output[, "data" := lapply(args, function(l) {
    do.call(.model_ebola_internal, l)
  })]

  # convert the raw data to output
  model_output[, "data" := Map(
    data, population,
    f = function(df, pop) {
      .output_to_df_ebola(df, pop, compartments)
    }
  )]

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
