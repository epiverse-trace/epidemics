#' @title Check arguments to default epidemic function
#' @name check_prepare_default_args
#' @rdname check_prepare_default_args
#'
#' @description Check and prepare the four main arguments to
#' [epidemic_default_cpp()] for use with [.epidemic_default_cpp()].
#'
#' `.check_args_epidemic_default()` adds an empty `<intervention>` and
#' `<vaccination>` object if these are missing from the model arguments.
#'
#' `.prepare_args_epidemics_default()` prepares arguments for
#' [.epidemic_default_cpp()], which is the C++ function that solves the default
#' ODE system using a Boost _odeint_ solver.
#' `.prepare_args_epidemics_default()` converts the arguments collected in
#' `mod_args` into simpler structures such as lists and numeric or integer
#' vectors that can be interpreted as C++ types such as `Rcpp::List`,
#' `Rcpp::NumericVector`, or `Eigen::MatrixXd`.
#'
#' @return
#'
#' `.check_args_epidemic_default()` invisibly returns the model arguments passed
#' in `mod_args`. If the model arguments did not previously contain elements
#' named `intervention` and `vaccination`, these are added as dummy objects of
#' the corresponding classes.
#'
#' `.prepare_args_epidemic_default()` returns a list of model arguments suitable
#' for [.epidemic_default_cpp()]. This is a named list consisting of:
#'
#'  1. `population`, the `<population>` object with the contact matrix
#' scaled by the largest real eigenvalue and by the size of each groups; the
#' initial conditions are also modified to represent absolute rather than
#' proportional values.
#'
#'  2. `beta`, a single number for the transmission rate of the infection.
#'
#'  3. `alpha`, a single number for the rate of transition from exposed to
#' infectious.
#'
#'  4. `gamma`, a single number for the recovery rate.
#'
#'  5. `intervention`, the `<intervention>` object,
#'
#'  6. `vaccination`, the `<vaccination>` object,
#'
#'  7. `time_end`, a single number for the time point at which to end the
#' simulation, and
#'
#'  8. `increment`,  a single number for the value by which the simulation time
#' is incremented.
#'
#' @keywords internal
.check_args_epidemic_default <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique"
  )

  # input checking on the infection object
  assert_infection(
    mod_args[["infection"]],
    extra_parameters = "preinfectious_period"
  )

  # load number of compartments to check initial conditions matrix
  compartments_default <- read_from_library(
    model_name = "default", what = "compartments"
  )
  assert_population(
    mod_args[["population"]],
    compartments = compartments_default
  )

  # add null intervention and vaccination if these are missing
  # if not missing, check that they conform to expectations
  if (!"intervention" %in% names(mod_args)) {
    mod_args[["intervention"]] <- no_intervention(
      mod_args[["population"]]
    )
  } else {
    assert_intervention(mod_args[["intervention"]], mod_args[["population"]])
  }
  if (!"vaccination" %in% names(mod_args)) {
    mod_args[["vaccination"]] <- no_vaccination(
      mod_args[["population"]]
    )
  } else {
    # default model only supports a single dose vaccination
    assert_vaccination(
      mod_args[["vaccination"]],
      doses = 1L,
      mod_args[["population"]]
    )
  }

  # return arguments invisibly
  invisible(mod_args)
}

#' @title Prepare arguments to default epidemic function
#' @name check_prepare_default_args
#' @rdname check_prepare_default_args
#' @keywords internal
.prepare_args_epidemic_default <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  contact_matrix <- mod_args[["population"]][["contact_matrix"]] /
    max(Re(eigen(mod_args[["population"]][["contact_matrix"]])$values))

  # scale rows of the contact matrix by the corresponding group population
  contact_matrix <- contact_matrix /
    mod_args[["population"]][["demography_vector"]]

  # prepare initial conditions by scaling with demography
  initial_state <-
    mod_args[["population"]][["initial_conditions"]] *
      mod_args[["population"]][["demography_vector"]]

  # calculate infection parameters
  gamma <- 1.0 / mod_args[["infection"]][["infectious_period"]]
  alpha <- 1.0 / mod_args[["infection"]][["preinfectious_period"]]
  beta <- mod_args[["infection"]][["r0"]] /
    mod_args[["infection"]][["infectious_period"]]

  # get NPI related times and contact reductions
  npi_time_begin <- mod_args[["intervention"]][["time_begin"]]
  npi_time_end <- mod_args[["intervention"]][["time_end"]]
  npi_cr <- mod_args[["intervention"]][["contact_reduction"]]

  # get vaccination related times and rates
  vax_time_begin <- mod_args[["vaccination"]][["time_begin"]]
  vax_time_end <- mod_args[["vaccination"]][["time_end"]]
  vax_nu <- mod_args[["vaccination"]][["nu"]]

  # return selected arguments for internal C++ function
  list(
    initial_state,
    beta, alpha, gamma,
    contact_matrix,
    npi_time_begin, npi_time_end, npi_cr,
    vax_time_begin, vax_time_end, vax_nu,
    mod_args$time_end,
    mod_args$increment
  )
}
