#' @title Prepare arguments to default epidemic function
#' @name prepare_default_args
#' @rdname prepare_default_args
#'
#' @description Prepare arguments to for s[.model_default_cpp()].
#'
#' @param mod_args A named list of the population, and epidemic modifiers.
#'
#' @return
#' A list of model arguments suitable
#' for [.model_default_cpp()]. This is a named list consisting of:
#'
#'  - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `transmissibility`, `infectiousness_rate`, `recovery_rate`: three numbers
#' representing the transmission rate of the infection, the rate of transition
#' from exposed to infectious, and the recovery rate, respectively;
#'
#'  - `contact_matrix`, a numeric matrix for the population contact matrix
#' scaled by the largest real eigenvalue and by the size of each groups;
#'
#'  - `npi_time_begin`, `npi_time_end`: two vectors for the start and end times
#' of any interventions applied;
#'
#'  - `npi_cr`: a matrix for the age- and intervention-specific effect on social
#' contacts;
#'
#'  - `vax_time_begin`,`vax_time_end`, `vax_nu`: three numeric matrices for the
#' age- and dose-specific start times, end times, and rates of any vaccination
#' doses implemented;
#'
#'  - `time_end`, `increment`: two numbers for the time at which to end the
#' simulation, and the value by which the simulation time
#' is incremented.
#'
#' @keywords internal
#' @details
#' `.prepare_args_model_default()` prepares arguments for
#' [.model_default_cpp()], which is the C++ function that solves the default
#' ODE system using a Boost _odeint_ solver, and for [.ode_model_default()],
#' which is passed to [deSolve::lsoda()] in [model_default_r()].
#'
#' `.prepare_args_model_default()` converts the arguments collected in
#' `mod_args` into simpler structures such as lists and numeric or integer
#' vectors that can be interpreted as C++ types such as `Rcpp::List`,
#' `Rcpp::NumericVector`, or `Eigen::MatrixXd`.
.prepare_args_model_default <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  contact_matrix <- get_parameter(mod_args[["population"]], "contact_matrix")
  contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

  # scale rows of the contact matrix by the corresponding group population
  contact_matrix <- contact_matrix /
    get_parameter(mod_args[["population"]], "demography_vector")

  # prepare initial conditions by scaling with demography
  initial_state <-
    get_parameter(mod_args[["population"]], "initial_conditions") *
      get_parameter(mod_args[["population"]], "demography_vector")

  # get NPI related times and contact reductions
  contact_interventions <- mod_args[["intervention"]][["contacts"]]
  npi_time_begin <- get_parameter(contact_interventions, "time_begin")
  npi_time_end <- get_parameter(contact_interventions, "time_end")
  npi_cr <- get_parameter(contact_interventions, "reduction")

  # get other interventions if any
  rate_interventions <- mod_args[["intervention"]][
    setdiff(names(mod_args[["intervention"]]), "contacts")
  ]

  # get vaccination related times and rates
  vax_time_begin <- get_parameter(mod_args[["vaccination"]], "time_begin")
  vax_time_end <- get_parameter(mod_args[["vaccination"]], "time_end")
  vax_nu <- get_parameter(mod_args[["vaccination"]], "nu")

  # return selected arguments for internal C++ function
  # order is important
  list(
    initial_state = initial_state,
    transmissibility = mod_args$transmissibility,
    infectiousness_rate = mod_args$infectiousness_rate,
    recovery_rate = mod_args$recovery_rate,
    contact_matrix = contact_matrix,
    npi_time_begin = npi_time_begin, npi_time_end = npi_time_end,
    npi_cr = npi_cr,
    vax_time_begin = vax_time_begin, vax_time_end = vax_time_end,
    vax_nu = vax_nu,
    rate_interventions = rate_interventions,
    time_dependence = mod_args$time_dependence,
    time_end = mod_args$time_end,
    increment = mod_args$increment
  )
}

.cross_check_intervention <- function(x, population, allowed_targets) {
  # create dummy intervention set
  tmp_intervention <- list(
    contacts = no_contacts_intervention(population),
    transmissibility = no_rate_intervention()
  )
  if (is.null(x)) {
    return(tmp_intervention)
  }

  # check that contact interventions are suitable for population
  if ("contacts" %in% names(x)) {
    assert_intervention(x[["contacts"]], "contacts", population)
  }

  # replace dummy values with user values if avaialable, and return
  tmp_intervention[names(x)] <- x
  tmp_intervention
}

.cross_check_vaccination <- function(x, population, doses) {
  if (is.null(x)) {
    no_vaccination(population)
  } else {
    assert_vaccination(mod_args[["vaccination"]], doses = doses, population)
    x
  }
}

.cross_check_timedep <- function(x, allowed_targets) {
  if (is.null(x)) {
    no_time_dependence()
  } else {
    checkmate::assert_names(names(x), subset.of = allowed_targets)
    x
  }
}

.cross_check_popchange <- function(x, population) {
  x
}
