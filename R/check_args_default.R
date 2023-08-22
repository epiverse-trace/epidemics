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
#'  - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `beta`, `alpha`, `gamma`: three numbers representing the transmission rate
#' of the infection, the rate of transition from exposed to infectious, and the
#' recovery rate, respectively;
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
    # add as a list element named "contacts", and one named "beta"
    mod_args[["intervention"]] <- list(
      contacts = no_contacts_intervention(
        mod_args[["population"]]
      ),
      # a dummy intervention on the rate parameter beta
      beta = no_rate_intervention()
    )
  } else {
    # if a single intervention is passed
    stopifnot(
      "`intervention` must be a list of <intervention>s" =
        checkmate::test_list(
          mod_args[["intervention"]],
          types = c("contacts_intervention", "rate_intervention"),
          names = "unique"
        )
    )

    # check for any other intervention list element names
    checkmate::assert_names(
      names(mod_args[["intervention"]]),
      subset.of = c("beta", "gamma", "alpha", "contacts")
    )

    if ("contacts" %in% names(mod_args[["intervention"]])) {
      # check the intervention on contacts
      assert_intervention(
        mod_args[["intervention"]][["contacts"]], "contacts",
        mod_args[["population"]]
      )
    }

    # if there is only an intervention on contacts, add a dummy intervention
    # on the transmission rate beta
    if (identical(names(mod_args[["intervention"]]), "contacts")) {
      mod_args[["intervention"]]$beta <- no_rate_intervention()
    }
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
  contact_matrix <- get_parameter(mod_args[["population"]], "contact_matrix")
  contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$values))

  # scale rows of the contact matrix by the corresponding group population
  contact_matrix <- contact_matrix /
    get_parameter(mod_args[["population"]], "demography_vector")

  # prepare initial conditions by scaling with demography
  initial_state <-
    get_parameter(mod_args[["population"]], "initial_conditions") *
      get_parameter(mod_args[["population"]], "demography_vector")

  # calculate infection parameters
  gamma <- 1.0 / get_parameter(mod_args[["infection"]], "infectious_period")
  alpha <- 1.0 / get_parameter(mod_args[["infection"]], "preinfectious_period")
  beta <- get_parameter(mod_args[["infection"]], "r0") /
    get_parameter(mod_args[["infection"]], "infectious_period")

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
  list(
    initial_state = initial_state,
    beta = beta, alpha = alpha, gamma = gamma,
    contact_matrix = contact_matrix,
    npi_time_begin = npi_time_begin, npi_time_end = npi_time_end,
    npi_cr = npi_cr,
    vax_time_begin = vax_time_begin, vax_time_end = vax_time_end,
    vax_nu = vax_nu,
    rate_interventions = rate_interventions,
    time_end = mod_args$time_end,
    increment = mod_args$increment
  )
}
