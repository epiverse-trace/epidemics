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
#' The `population`, `transmission_rate`, `infectiousness_rate`, and
#' `recovery_rate`
#' arguments are mandatory, while passing an `intervention` and `vaccination`
#' are optional and can be used to simulate scenarios with different epidemic
#' responses or different levels of the same type of response.
#' See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param transmission_rate A numeric for the rate at which individuals
#' move from the susceptible to the exposed compartment upon contact with an
#' infectious individual. Often denoted as \eqn{\beta}, with
#' \eqn{\beta = R_0 / \text{infectious period}}. See **Details** for default
#' values.
#' @param infectiousness_rate A numeric for the rate at which individuals
#' move from the exposed to the infectious compartment. Often denoted as
#' \eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
#' This value does not depend upon the number of infectious individuals in the
#' population. See **Details** for default values.
#' @param recovery_rate A numeric for the rate at which individuals move
#' from the infectious to the recovered compartment. Often denoted as
#' \eqn{\gamma}, with \eqn{\gamma = 1.0 / \text{infectious period}}.
#' See **Details** for default values.
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
#' # Details: SEIRV model suitable for directly transmitted infections
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
#' - Transmission rate (\eqn{\beta}, `transmission_rate`): 0.186, assuming an
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
#' # and three discrete values of transmission rate
#' data <- model_default(
#'   population = uk_population,
#'   transmission_rate = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
#' )
#'
#' # view some data
#' data
#'
#' # run epidemic simulations with differences in the end time
#' # may be useful when considering different start dates with a fixed end point
#' data <- model_default(
#'   population = uk_population,
#'   time_end = c(50, 100, 150)
#' )
#'
#' data
#' @export
model_default <- function(population,
                          transmission_rate = 1.3 / 7.0,
                          infectiousness_rate = 1.0 / 2.0,
                          recovery_rate = 1.0 / 7.0,
                          intervention = NULL,
                          vaccination = NULL,
                          time_dependence = NULL,
                          time_end = 100,
                          increment = 1) {
  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )
  assert_population(population, compartments)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  checkmate::assert_integerish(time_end, lower = 0)

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)
  checkmate::assert_number(increment, lower = 1e-3, finite = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if `intervention` is a single intervention set or a list of such sets
  # NULL is allowed;
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) &&
    # Check that all elements of intervention sets are either `<intervention>`
    # or NULL
    all(
      vapply(
        unlist(intervention, recursive = FALSE),
        FUN = function(x) {
          is_intervention(x) || is.null(x)
        }, TRUE
      )
    )

  # Check if parameters can be recycled;
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls,
    # Check if `vaccination` is a single vaccination, NULL, or a list
    "`vaccination` must be a <vaccination> or a list of <vaccination>s" =
      is_vaccination(vaccination) || checkmate::test_list(
        vaccination,
        types = c("vaccination", "null"), null.ok = TRUE
      )
  )

  # make lists if not lists
  population <- list(population) # NOTE: currently not list, but see issue #181
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
      c("transmission_rate", "infectiousness_rate", "recovery_rate")
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
  model_output[, "data" := Map(population, args, f = function(p, l) {
    .output_to_df(
      do.call(.model_default_cpp, l),
      population = p, # taken from local scope/env
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

#' Ordinary Differential Equations for the Default Model
#'
#' @description Provides the ODEs for the default SEIR-V model in a format that
#' is suitable for passing to [deSolve::lsoda()].
#' See [model_default()] for a list of required parameters.
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
  # no input checking, fn is unsafe and not expected to be used
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
    "transmission_rate", "infectiousness_rate", "recovery_rate"
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

  model_params <- .intervention_on_rates(
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
  sToE <- (model_params[["transmission_rate"]] * y[, 1] *
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

#' @importFrom odin odin
#' @export
seirv_model <- odin::odin({
  # Time-dependent parameters
  # beta <- interpolate(beta_t, beta_y, "linear")
  # sigma <- interpolate(sigma_t, sigma_y, "linear")
  # gamma <- interpolate(gamma_t, gamma_y, "linear")

  # Contact matrix with time-dependent interventions
  # Then multiply across interventions - need to check this
  contact_reduction[, ] <- if (t > intervention_start[i] && t < intervention_end[i]) (1.0 - intervention_effect[i, j]) else 1 # nolint: line_length_linter.
  contact_reduction_log[, ] <- log(contact_reduction[i, j] + 1e-10) # Avoid zero
  contact_reduction_total_log[] <- sum(contact_reduction_log[, i])
  contact_reduction_total[] <- exp(contact_reduction_total_log[i])

  # Specify how transmission varies over time
  # FOI is contacts * infectious * transmission rate
  # returns a matrix, must be converted to a vector/1D array
  lambda_prod[, ] <- C[i, j] * I[j] * beta * contact_reduction_total[j] * contact_reduction_total[i] # nolint: line_length_linter.
  lambda[] <- sum(lambda_prod[i, ])

  # Vaccination - indexing over age groups
  vax_rate[] <- if (t >= vax_start[i] && t < vax_end[i]) vax_nu[i] else 0

  # ODEs
  deriv(S[]) <- -(lambda[i] * S[i]) - (vax_rate[i] * S[i])
  deriv(E[]) <- lambda[i] * S[i] - sigma * E[i]
  deriv(I[]) <- sigma * E[i] - gamma * I[i]
  deriv(R[]) <- gamma * I[i]
  deriv(V[]) <- vax_rate[i] * S[i]

  # Initial conditions
  initial(S[]) <- init_S[i]
  initial(E[]) <- init_E[i]
  initial(I[]) <- init_I[i]
  initial(R[]) <- init_R[i]
  initial(V[]) <- init_V[i]

  # User-defined parameters
  C[, ] <- user()
  n_age <- user()
  n_intervention <- user()
  beta <- user()
  sigma <- user()
  gamma <- user()

  # Not currently used (can re-add for time-dependent parameters)
  # beta_t[] <- user()
  # beta_y[] <- user()
  # sigma_t[] <- user()
  # sigma_y[] <- user()
  # gamma_t[] <- user()
  # gamma_y[] <- user()
  intervention_start[] <- user()
  intervention_end[] <- user()
  intervention_effect[, ] <- user()
  vax_start[] <- user()
  vax_end[] <- user()
  vax_nu[] <- user()
  init_S[] <- user()
  init_E[] <- user()
  init_I[] <- user()
  init_R[] <- user()
  init_V[] <- user()

  # Dimensions
  dim(C) <- c(n_age, n_age)
  dim(lambda_prod) <- c(n_age, n_age)
  dim(lambda) <- n_age
  dim(S) <- n_age
  dim(E) <- n_age
  dim(I) <- n_age
  dim(R) <- n_age
  dim(V) <- n_age
  dim(contact_reduction) <- c(n_intervention, n_age)
  dim(contact_reduction_log) <- c(n_intervention, n_age)
  dim(contact_reduction_total_log) <- c(n_age)
  dim(contact_reduction_total) <- c(n_age)
  dim(intervention_start) <- n_intervention
  dim(intervention_end) <- n_intervention
  dim(intervention_effect) <- c(n_intervention, n_age)
  dim(vax_rate) <- n_age
  dim(vax_start) <- n_age
  dim(vax_end) <- n_age
  dim(vax_nu) <- n_age
  dim(init_S) <- n_age
  dim(init_E) <- n_age
  dim(init_I) <- n_age
  dim(init_R) <- n_age
  dim(init_V) <- n_age
})

#' @export
model_default_odin <- function(population,
                               transmission_rate = 1.5 / 7.0,
                               infectiousness_rate = 1.0 / 2.0,
                               recovery_rate = 1.0 / 7.0,
                               intervention = NULL,
                               vaccination = NULL,
                               time_dependence = NULL,
                               time_end = 100,
                               increment = 1) {
  # # get compartment names
  # compartments <- c(
  #   "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  # )
  # assert_population(population, compartments)

  # # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  # checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
  # checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  # checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  # checkmate::assert_integerish(time_end, lower = 0)

  # # check the time end and increment
  # # restrict increment to lower limit of 1e-6
  # checkmate::assert_integerish(time_end, lower = 0)
  # checkmate::assert_number(increment, lower = 1e-3, finite = TRUE)

  # # Prepare model parameters
  # # Scale initial conditions for odin model
  # initial_conditions <- population$initial_conditions * population$demography_vector # nolint: line_length_linter.
  # n_age <- nrow(population$contact_matrix)
  # time_points <- seq(0, time_end, by = increment)

  # # prepare contact matrix, divide by leading eigenvalue and rowwise by popsize
  # contact_matrix_norm <- population$contact_matrix
  # contact_matrix_norm <- (contact_matrix_norm / max(Re(eigen(contact_matrix_norm)$values))) / # nolint: line_length_linter.
  #   population$demography_vector


  # # Prepare intervention parameters
  # # Add null intervention if needed as odin model requires matrix input
  # intervention_start <- intervention_end <- 0
  # intervention_effect <- rep(0, n_age)
  # if (!is.null(intervention) && length(intervention$time_begin) >= 2) {
  #   intervention_start <- as.numeric(intervention$time_begin)
  #   intervention_end <- as.numeric(intervention$time_end)
  #   intervention_effect <- t(intervention$reduction) # row is the intervention
  # }
  # # If zero or one intervention
  # if (!is.null(intervention) || length(intervention$time_begin) < 2) {
  #   null_intervention <- intervention(
  #     name = "Null closure",
  #     type = "contacts",
  #     time_begin = 1e5,
  #     time_end = 1e5 + 1,
  #     reduction = matrix(rep(0, n_age))
  #   )
  #   if (is(intervention, "list")) {
  #     intervention <- intervention[[1]]
  #   }
  #   intervention <- c(null_intervention, null_intervention, intervention)
  #   intervention_start <- as.numeric(intervention$time_begin)
  #   intervention_end <- as.numeric(intervention$time_end)
  #   intervention_effect <- t(intervention$reduction) # row is the intervention
  # }

  # n_intervention <- length(intervention_start)

  # # Prepare vaccination parameters
  # vax_start <- vax_end <- rep(0, n_age)
  # vax_nu <- rep(0, n_age)
  # if (!is.null(vaccination)) {
  #   vax_start <- as.numeric(vaccination$time_begin)
  #   vax_end <- as.numeric(vaccination$time_end)
  #   vax_nu <- as.numeric(vaccination$nu)
  # }

  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )
  assert_population(population, compartments)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  checkmate::assert_integerish(time_end, lower = 0)

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)
  checkmate::assert_number(increment, lower = 1e-3, finite = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if `intervention` is a single intervention set or a list of such sets
  # NULL is allowed;
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) &&
    # Check that all elements of intervention sets are either `<intervention>`
    # or NULL
    all(
      vapply(
        unlist(intervention, recursive = FALSE),
        FUN = function(x) {
          is_intervention(x) || is.null(x)
        }, TRUE
      )
    )

  # Check if parameters can be recycled;
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls,
    # Check if `vaccination` is a single vaccination, NULL, or a list
    "`vaccination` must be a <vaccination> or a list of <vaccination>s" =
      is_vaccination(vaccination) || checkmate::test_list(
        vaccination,
        types = c("vaccination", "null"), null.ok = TRUE
      )
  )

  # make lists if not lists
  population <- list(population) # NOTE: currently not list, but see issue #181
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
      c("transmission_rate", "infectiousness_rate", "recovery_rate")
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

  population <- model_output$population[[1]]
  C <- population$contact_matrix
  n_age <- nrow(C)
  intervention <- model_output$intervention[[1]][[1]]
  intervention_type <- strsplit(class(intervention)[1], "_", fixed = TRUE)[[1]][1]
  intervention_reduction <- ifelse(intervention_type == "rate", 0, matrix(rep(0, n_age)))
  null_intervention <- intervention(
    name = "Null closure",
    type = intervention_type,
    time_begin = 1e5,
    time_end = 1e5 + 1,
    reduction = intervention_reduction
  )
  # Always add null intervention to cater for odin
  intervention <- c(null_intervention, intervention)

  n_intervention <- length(intervention$time_begin)

  intervention_start <- as.numeric(intervention$time_begin)
  intervention_end <- as.numeric(intervention$time_end)
  if (intervention_type == "rate") {
    intervention_effect <- matrix(rep(intervention$reduction, n_age), ncol = n_age)
  } else {
    intervention_effect <- t(intervention$reduction)
  }
  # row is the intervention


  print(intervention_effect)
  #       [,1] [,2]
  # npi_1  0.5  0.0
  # npi_2  0.1  0.5
  print(class(intervention_effect))
  # [1] "matrix" "array"
  print(intervention_start)
  # [1] 0 0
  print(class(intervention_start))
  # [1] "numeric"
  print(intervention_end)
  # [1] 100 100
  print(class(intervention_end))
  # [1] "numeric"
  beta <- model_output$transmission_rate
  sigma <- model_output$infectiousness_rate
  gamma <- model_output$recovery_rate
  vaccination <- model_output$vaccination[[1]]
  if (is.null(vaccination)) {
    vax_start <- vax_end <- vax_nu <- rep(0, n_age)
    # vax_nu <- rep(0, n_age)
    # vax_start <- as.numeric(vaccination$time_begin)
    # vax_end <- as.numeric(vaccination$time_end)
    # vax_nu <- as.numeric(vaccination$nu)
  }
  initial_conditions <- population$initial_conditions
  init_S <- initial_conditions[, 1]
  init_E <- initial_conditions[, 2]
  init_I <- initial_conditions[, 3]
  init_R <- initial_conditions[, 4]
  init_V <- initial_conditions[, 5]

  # Initialize and run the model
  model <- seirv_model$new(
    C = C,
    n_age = n_age,
    n_intervention = n_intervention,
    beta = beta,
    sigma = sigma,
    gamma = gamma,
    intervention_start = intervention_start,
    intervention_end = intervention_end,
    intervention_effect = intervention_effect,
    vax_start = vax_start,
    vax_end = vax_end,
    vax_nu = vax_nu,
    init_S = init_S,
    init_E = init_E,
    init_I = init_I,
    init_R = init_R,
    init_V = init_V
  )

  time_points <- seq(0, time_end, by = increment)
  result <- model$run(time_points)

  # Add scenario information
  dt <- data.table::as.data.table(result)
  # declaring variables below to avoid data.table related lintr messages
  temp <- value <- temp_compartment <- temp_demography <-
    compartment <- demography_group <- `:=` <- NULL

  age_group_mappings <- paste0(
    seq_len(n_age),
    row.names(population$contact_matrix)
  )
  names(age_group_mappings) <- seq_len(n_age)

  mapping <- c( # prepend numbers to help during sorting. Will remove later
    S = "1susceptible", E = "2exposed", I = "3infectious",
    R = "4recovered", V = "5vaccinated", age_group_mappings
  )

  # Melt the data table to long format
  data.table::melt(dt,
    id.vars = "t",
    variable.name = "temp", # e.g. S[1], ..., V[3]
    value.name = "value"
  )[ # piping the data.table way. Possible because melt outputs a data.table
    , list(
      time = t, # alternative to using data.table::setnames(dt, "t", "time")
      temp_compartment = substring(temp, 1L, 1L), # e.g. S[1] -> S
      temp_demography = substring(temp, 3L, 3L), # e.g. S[1] -> 1
      value
    )
  ][ # |> the DT way (piping the data.table way)
    ,
    list(
      time,
      demography_group = mapping[temp_demography], # e.g. 1[0,20), 2[20,65), ...
      compartment = mapping[temp_compartment], # e.g. 1susceptible, 2exposed
      value
    )
  ][ # |> the DT way
    order(time, compartment, demography_group) # prepending numbers helps here
  ][ # |> the DT way
    ,
    `:=`( # used as the prefix form to update multiple columns
      # remove prepended numbers from `mapping`
      demography_group = substring(demography_group, 2L), # e.g. [0,20), ...
      compartment = substring(compartment, 2L) # e.g. susceptible, exposed, ...
    )
  ][ # |> the DT way
    # added because the previous operation used `:=` which doesn't output
  ]
}
