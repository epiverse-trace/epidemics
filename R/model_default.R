#' @title Model an SEIR-V epidemic with interventions
#'
#' @name model_default
#' @rdname model_default
#'
#' @description Simulate an epidemic using a deterministic, compartmental
#' epidemic model with the compartments
#' "susceptible", "exposed", "infectious", "recovered", and "vaccinated".
#' This model can accommodate heterogeneity in social contacts among demographic
#' groups, as well as differences in the sizes of demographic groups.
#'
#' The `population`, `transmissibility`, `infectiousness_rate`, and
#' `recovery_rate`
#' arguments are mandatory, while passing an `intervention` and `vaccination`
#' are optional and can be used to simulate scenarios with different epidemic
#' responses or different levels of the same type of response.
#' See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param transmissibility A numeric for the rate at which individuals
#' move from the susceptible to the exposed compartment upon contact with an
#' infectious individual. Often denoted as \eqn{\beta}, with
#' \eqn{\beta = R_0 / \text{infectious period}}.
#' @param infectiousness_rate A numeric for the rate at which individuals
#' move from the exposed to the infectious compartment. Often denoted as
#' \eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
#' This value does not depend upon the number of infectious individuals in the
#' population.
#' @param recovery_rate A numeric for the rate at which individuals move
#' from the infectious to the recovered compartment. Often denoted as
#' \eqn{\gamma}, with \eqn{\gamma = 1.0 / \text{infectious period}}.
#' @param intervention A named list of `<intervention>`s representing optional
#' non-pharmaceutical or pharmaceutical interventions applied during the
#' epidemic. Only a single intervention on social contacts of the class
#' `<contacts_intervention>` is allowed as the named element "contacts".
#' Multiple `<rate_interventions>` on the model parameters are allowed; see
#' **Details** for the model parameters for which interventions are supported.
#' @param vaccination A `<vaccination>` object representing an optional
#' vaccination regime with a single dose, followed during the course of the
#' epidemic, with a start and end time, and age-specific vaccination rates.
#' @param time_dependence A named list where each name
#' is a model parameter, and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 100 days. May be a numeric vector.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details
#'
#' ## R and Rcpp implementations
#'
#' `model_default_cpp()` is a wrapper function for the internal C++ function
#' [.model_default_cpp()] that uses Boost _odeint_ solvers, while
#' `model_default_r()` is a wrapper around [deSolve::lsoda()] which an R-only
#' implementation of the ODE system in `.ode_model_default()`.
#' Both models return equivalent results, but the C++ implementation is faster.
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates of
#' transitions between compartments per model run.
#'
#' However, model parameters may be passed as numeric vectors. These vectors
#' must follow Tidyverse recycling rules: all vectors must have the same length,
#' or, vectors of length 1 will be recycled to the length of any other vector.
#'
#' The default values are:
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.186, assuming an
#' \eqn{R_0} = 1.3 and an infectious period of 7 days.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.5, assuming
#' a pre-infectious period of 2 days.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.143, assuming an
#' infectious period of 7 days.
#'
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
#' @examples
#' # create a population
#' uk_population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0, 0),
#'     nrow = 1, ncol = 5L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' # and three discrete values of transmissibility
#' data <- model_default_cpp(
#'   population = uk_population,
#'   transmissibility = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
#' )
#'
#' # view some data
#' data
#'
#' # run epidemic simulations with differences in the end time
#' # may be useful when considering different start dates with a fixed end point
#' data <- model_default_cpp(
#'   population = uk_population,
#'   time_end = c(50, 100, 150)
#' )
#'
#' data
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
  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )
  assert_population(population, compartments)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(time_end, lower = 0, finite = TRUE)

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)
  checkmate::assert_number(increment, lower = 1e-3, finite = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if parameters can be recycled;
  # Check if `population` is a single population or a list of such
  # and convert to list for a data.table list column;
  # also check if `intervention` is a list of interventions or a list-of-lists
  # and convert to a list for a data.table list column. NULL is allowed;
  # Check if `vaccination` is a single vaccination or a list
  # and convert to a list for a data.table list column
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) && all(
    vapply(
      unlist(intervention, recursive = FALSE),
      FUN = function(x) {
        is_intervention(x) || is.null(x)
      }, TRUE
    )
  )

  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      test_recyclable(params),
    "`population` must be a <population> or a list of <population>s" =
      is_population(population) || checkmate::test_list(
        population,
        types = "population"
      ),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls,
    "`vaccination` must be a <vaccination> or a list of <vaccination>s" =
      is_vaccination(vaccination) || checkmate::test_list(
        vaccination,
        type = c("vaccination", "null"), null.ok = TRUE
      )
  )

  # make lists if not lists
  if (is_population(population)) {
    population <- list(population)
  }
  if (is_lofints) {
    intervention <- list(intervention)
  }
  if (is_vaccination(vaccination) || is.null(vaccination)) {
    vaccination <- list(vaccination)
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
  time_dependence <- list(
    .cross_check_timedep(
      time_dependence,
      c("transmissibility", "infectiousness_rate", "recovery_rate")
    )
  )

  # collect parameters and add a parameter set identifier
  params <- data.table::as.data.table(params)
  params[, "param_set" := .I]

  # this nested data.table will be returned
  model_output <- data.table::CJ(
    population = population,
    intervention = intervention,
    vaccination = vaccination,
    time_dependence = time_dependence,
    increment = increment,
    sorted = FALSE
  )

  # process the population, interventions, and vaccinations, after
  # cross-checking them agains the relevant population
  model_output[, args := apply(model_output, 1, function(x) {
    .check_prepare_args_default(c(x))
  })]
  model_output[, "scenario" := .I]

  # combine infection parameters and scenarios
  # NOTE: join X[Y] must have params as X as list cols not supported for X
  model_output <- params[, as.list(model_output), by = names(params)]

  # collect model arguments in column data, then overwrite
  model_output[, args := apply(model_output, 1, function(x) {
    c(x[["args"]], x[param_names]) # avoid including col "param_set"
  })]
  model_output[, data := Map(population, args, f = function(p, l) {
    .output_to_df(
      do.call(.model_default_cpp, l),
      population = p, # taken from local scope/env
      compartments = compartments
    )
  })]

  # check for single row output, i.e., scalar arguments, and return data.frame
  # do not return the parameters in this case
  if (nrow(model_output) == 1L) {
    model_output <- model_output[["data"]][[1L]] # hardcoded for special case
  }

  # return data.table
  model_output[]
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
  # NOTE: `model_params` refers to epdemiological parameters, while
  # `params` refers to the list passed as a function argument.
  # e.g. contact matrix, or interventions, are not included in `model_params`
  model_params <- params[c(
    "transmissibility", "infectiousness_rate", "recovery_rate"
  )]

  # apply time dependence before interventions
  time_dependent_params <- Map(
    model_params[names(params$time_dependence)],
    params$time_dependence,
    f = function(x, func) {
      func(time = t, x = x)
    }
  )

  # assign time-modified param values
  model_params[names(time_dependent_params)] <- time_dependent_params

  model_params <- intervention_on_rates(
    t = t,
    interventions = params[["rate_interventions"]],
    parameters = model_params
  )

  # modify the vaccination rate depending on the regime
  # the number of doses is already checked before passing
  current_nu <- params[["vax_nu"]] *
    ((params[["vax_time_begin"]] < t) &
      (params[["vax_time_end"]] > t))

  # calculate transitions
  sToE <- (model_params[["transmissibility"]] * y[, 1] *
    contact_matrix_ %*% y[, 3])
  eToI <- model_params[["infectiousness_rate"]] * y[, 2]
  iToR <- model_params[["recovery_rate"]] * y[, 3]
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
#' @name model_default
#' @rdname model_default
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
  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(recovery_rate, lower = 0, finite = TRUE)

  # all intervention sub-classes pass check for intervention superclass
  checkmate::assert_list(
    intervention,
    types = "intervention", null.ok = TRUE,
    names = "unique", any.missing = FALSE
  )
  # specifics of vaccination doses are checked in dedicated function
  checkmate::assert_class(vaccination, "vaccination", null.ok = TRUE)

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    names = "unique", any.missing = FALSE
  )
  # lapply on null returns an empty list
  invisible(
    lapply(time_dependence, checkmate::assert_function,
      args = c("time", "x"),
      ordered = TRUE
    )
  )

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect all model arguments
  model_arguments <- list(
    population = population,
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    intervention = intervention,
    vaccination = vaccination,
    time_dependence = time_dependence,
    time_end = time_end, increment = increment
  )

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
  data <- .output_to_df(
    output = list(
      x = data[, setdiff(colnames(data), "time")],
      time = seq(0, time_end, increment)
    ),
    population = population,
    compartments = compartments
  )
  data.table::setDF(data)[]
}
