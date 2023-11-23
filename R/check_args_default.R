#' @title Check arguments to default epidemic function
#' @name check_prepare_default_args
#' @rdname check_prepare_default_args
#'
#' @description Check and prepare the four main arguments to
#' [model_default_cpp()] for use with [.model_default_cpp()].
#'
#' @return
#'
#' `.check_args_model_default()` invisibly returns the model arguments passed
#' in `mod_args`. Model functionalities, such as vaccination or interventions,
#' passed as `NULL` are replaced with dummy values for internal functions.
#'
#' `.prepare_args_model_default()` returns a list of model arguments suitable
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
.check_args_model_default <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique"
  )

  # load number of compartments to check initial conditions matrix
  compartments_default <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )
  assert_population(
    mod_args[["population"]],
    compartments = compartments_default
  )

  # add null intervention and vaccination if these are missing
  # if not missing, check that they conform to expectations
  # add null rate_intervention if this is missing
  # if not missing, check that it conforms to expectations
  if (is.null(mod_args[["intervention"]])) {
    # add dummy list elements named "contacts", and one named "transmissibility"
    mod_args[["intervention"]] <- list(
      contacts = no_contacts_intervention(
        mod_args[["population"]]
      ),
      transmissibility = no_rate_intervention()
    )
  } else {
    # check intervention list names
    checkmate::assert_names(
      names(mod_args[["intervention"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "recovery_rate", "contacts"
      )
    )
    # if a contacts intervention is passed, check it
    if ("contacts" %in% names(mod_args[["intervention"]])) {
      # check the intervention on contacts
      assert_intervention(
        mod_args[["intervention"]][["contacts"]], "contacts",
        mod_args[["population"]]
      )
    } else {
      # if not contacts intervention is passed, add a dummy one
      mod_args[["intervention"]]$contacts <- no_contacts_intervention(
        mod_args[["population"]]
      )
    }

    # if there is only an intervention on contacts, add a dummy intervention
    # on the transmissibility
    if (identical(names(mod_args[["intervention"]]), "contacts")) {
      mod_args[["intervention"]]$transmissibility <- no_rate_intervention()
    }
  }

  if (is.null(mod_args[["vaccination"]])) {
    mod_args[["vaccination"]] <- no_vaccination(
      mod_args[["population"]]
    )
  } else {
    # default model only supports a single dose vaccination
    assert_vaccination(
      mod_args[["vaccination"]],
      doses = 1L, mod_args[["population"]]
    )
  }

  # handle time dependence if not present, and check targets if present
  if (is.null(mod_args[["time_dependence"]])) {
    mod_args[["time_dependence"]] <- no_time_dependence()
  } else {
    checkmate::assert_names(
      names(mod_args[["time_dependence"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "recovery_rate"
      )
    )
  }

  # return arguments invisibly
  invisible(mod_args)
}

#' @title Prepare arguments to default epidemic function
#' @name check_prepare_default_args
#' @rdname check_prepare_default_args
#' @keywords internal
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
