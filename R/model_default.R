#' @title Model an SEIR-V epidemic with interventions
#'
#' @name epidemic_default
#' @rdname epidemic_default
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
#' @param time_dependence A named list where each name
#' is a model parameter (see `infection`), and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 200 days.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details
#'
#' `model_default_cpp()` is a wrapper function for [.model_default_cpp()],
#' an internal C++ function that uses Boost _odeint_ solvers for an SEIR-V model
#' .
#' [.model_default_cpp()] accepts arguments that
#' are created by processing the `population`, `infection`, `intervention` and
#' `vaccination` arguments to the wrapper function into simpler forms.
#'
#' `model_default_r()` is a wrapper around the internal function
#' `.ode_model_default()`, which is passed to [deSolve::lsoda()].
#'
#' Both models return equivalent results, but the C++ implementation is faster.
#'
#' This model only allows for single, population-wide rates of
#' transition between the 'susceptible' and 'exposed' compartments, between the
#' 'exposed' and 'infectious' compartments, and in the recovery rate.
#'
#' Model rates or parameters can be made time-dependent by passing a function
#' which modifies the parameter based on the current ODE simulation time.
#' For example, a function that modifies the transmission rate `beta` could be
#' passed as `time_dependence = list(beta = function(time, x) x + sinpi(time))`.
#' This functionality may be used to model events that are expected to have some
#' effect on model parameters, such as seasonality or annual schedules such as
#' holidays.
#'
#' @return A `data.table` with the columns "time", "compartment", "age_group",
#' "value". The compartments correspond to the compartments of the model
#' chosen with `model`.
#' The current default model has the compartments "susceptible", "exposed",
#' "infectious", "recovered", and "vaccinated".
#' @export
model_default_cpp <- function(population,
                              transmissibility = 1.3 / 7.0,
                              infectiousness_rate = 1.0 / 2.0,
                              recovery_rate = 1.0 / 7.0,
                              intervention = NULL,
                              vaccination = NULL,
                              time_dependence = NULL,
                              time_end = 100,
                              increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(recovery_rate, lower = 0, finite = TRUE)

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
    time_end = time_end, increment = increment
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_list(intervention, types = "intervention")
    model_arguments[["intervention"]] <- intervention
  }
  if (!is.null(vaccination)) {
    checkmate::assert_class(vaccination, "vaccination")
    model_arguments[["vaccination"]] <- vaccination
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
    model_arguments[["time_dependence"]] <- time_dependence
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_model_default(
    .check_args_model_default(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )

  # run model over arguments
  output <- do.call(.model_default_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}

#' Ordinary Differential Equations for the Default Model
#'
#' @description Provides the ODEs for the default SEIR-V model in a format that
#' is suitable for passing to [deSolve::lsoda()].
#' See [model_default_r()] for a list of required parameters.
#'
#' @param t A single number of the timestep at which to integrate.
#' @param y The conditions of the epidemiological compartments.
#' @param params The parameters, passed as a named list.
#'
#' @return A list with a vector with as many elements as the number of
#' demographic groups times the number of epidemiological compartments. Each
#' value gives the change in the number of individuals in that compartment.
#' @keywords internal
.ode_model_default <- function(t, y, params) {
  # no input checking, fn is expected to be called only in model_default_r()
  n_age <- nrow(params[["contact_matrix"]])

  # create a matrix
  y <- matrix(y, nrow = n_age, ncol = 5L, byrow = FALSE)

  # scale the contact matrix if within the intervention period
  contact_matrix_ <- intervention_on_cm(
    t = t,
    cm = params[["contact_matrix"]],
    time_begin = params[["npi_time_begin"]],
    time_end = params[["npi_time_end"]],
    cr = params[["npi_cr"]]
  )

  # get paramters to modify them
  infection_params <- params[c(
    "transmissibility", "infectiousness_rate", "recovery_rate"
  )]

  # apply time dependence before interventions
  time_dependent_params <- Map(
    infection_params[names(params$time_dependence)],
    params$time_dependence,
    f = function(x, func) {
      func(time = t, x = x)
    }
  )

  # assign time-modified param values
  infection_params[names(time_dependent_params)] <- time_dependent_params

  infection_params <- intervention_on_rates(
    t = t,
    interventions = params[["rate_interventions"]],
    parameters = infection_params
  )

  # modify the vaccination rate depending on the regime
  # the number of doses is already checked before passing
  current_nu <- params[["vax_nu"]] *
    ((params[["vax_time_begin"]] < t) &
      (params[["vax_time_end"]] > t))

  # calculate transitions
  sToE <- (infection_params[["transmissibility"]] * y[, 1] *
    contact_matrix_ %*% y[, 3])
  eToI <- infection_params[["infectiousness_rate"]] * y[, 2]
  iToR <- infection_params[["recovery_rate"]] * y[, 3]
  sToV <- current_nu * y[, 1]

  # define compartmental changes
  dS <- -sToE - sToV
  dE <- sToE - eToI
  dI <- eToI - iToR
  dR <- iToR
  dV <- sToV

  # return a list
  list(c(dS, dE, dI, dR, dV))
}

#' @title Model an SEIR-V epidemic with interventions
#'
#' @name epidemic_default
#' @rdname epidemic_default
#'
#' @export
model_default_r <- function(population,
                            transmissibility = 1.3 / 7.0,
                            infectiousness_rate = 1.0 / 2.0,
                            recovery_rate = 1.0 / 7.0,
                            intervention = NULL,
                            vaccination = NULL,
                            time_dependence = NULL,
                            time_end = 100,
                            increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")

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
    time_end = time_end, increment = increment,
    time_dependence = time_dependence
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_list(
      intervention,
      types = c("intervention", "list")
    )
    model_arguments[["intervention"]] <- intervention
  }
  if (!is.null(vaccination)) {
    checkmate::assert_class(vaccination, "vaccination")
    model_arguments[["vaccination"]] <- vaccination
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
    model_arguments[["time_dependence"]] <- time_dependence
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_model_default(
    .check_args_model_default(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )

  # get compartment states over timesteps
  data <- deSolve::lsoda(
    y = model_arguments[["initial_state"]],
    times = seq(0, time_end, increment),
    func = .ode_model_default,
    parms = model_arguments
  )

  # convert to long format using output_to_df() and return
  output_to_df(
    output = list(
      x = data[, setdiff(colnames(data), "time")],
      time = seq(0, time_end, increment)
    ),
    population = population,
    compartments = compartments
  )
}
