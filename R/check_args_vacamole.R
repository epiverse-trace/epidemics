#' @title Check arguments to the Vacamole epidemic function
#' @name check_prepare_vacamole_args
#' @rdname check_prepare_vacamole_args
#'
#' @description Check and prepare the four main arguments to
#' [model_vacamole_cpp()] for use with [.model_vacamole_cpp()].
#'
#' `.check_args_model_vacamole()` adds an empty `<intervention>` object if
#' this is missing from the model arguments.
#'
#' @return
#'
#' `.check_args_model_vacamole()` invisibly returns the model arguments
#' passed in `mod_args`. If the model arguments did not previously contain
#' elements named `intervention` this is added as dummy objects of the
#' corresponding classes.
#'
#' `.prepare_args_model_vacamole()` returns a list of model arguments
#' suitable for [.model_vacamole_cpp()]. This is a named list consisting of:
#'
#'  - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `transmissibility`, `transmissibility_vax`: two numbers representing the
#' transmission rate
#' of the infection for unvaccinated or single-dose vaccinated, and two-dose
#' vaccinated individuals, respectively;
#'
#'  - `infectiousness_rate`: a single number for the transition rate from the
#' 'exposed' and 'exposed_vaccinated' to the 'infectious' and
#' 'infectious_vaccinated' compartments;
#'
#'  - `mortality_rate`, `mortality_rate_vax`: two numbers representing the
#' mortality rate of the infection for unvaccinated or single-dose vaccinated,
#' and two-dose vaccinated individuals, respectively;
#'
#' - `hospitalisation_rate`, `hospitalisation_rate_vax`: two numbers
#' representing the hospitalisation rate of the infection for unvaccinated or
#' single-dose vaccinated, and two-dose vaccinated individuals, respectively;
#'
#' - `recovery_rate`: a single number for the recovery rate from the infection;
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
#' `.prepare_args_model_vacamole()` prepares arguments for
#' [.model_vacamole_cpp()], which is the C++ function that solves the default
#' ODE system using a Boost _odeint_ solver, and for [.ode_model_vacamole()],
#' which is passed to [deSolve::lsoda()] in [model_vacamole_r()].
#'
#' `.prepare_args_model_vacamole()` converts the arguments collected in
#' `mod_args` into simpler structures such as lists and numeric or integer
#' vectors that can be interpreted as C++ types such as `Rcpp::List`,
#' `Rcpp::NumericVector`, or `Eigen::MatrixXd`.
.check_args_model_vacamole <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique",
    must.include = "vaccination" # vaccination necessary for Vacamole
  )

  # load number of compartments to check initial conditions matrix
  compartments_vacamole <- c(
    "susceptible", "vaccinated_one_dose", "vaccinated_two_dose",
    "exposed", "exposed_vaccinated", "infectious", "infectious_vaccinated",
    "hospitalised", "hospitalised_vaccinated", "dead", "recovered"
  )
  # input checking on population parameters
  assert_population(
    mod_args[["population"]],
    compartments = compartments_vacamole
  )

  # input checking on vaccination parameters
  # the Vacamole model expects two vaccination doses
  assert_vaccination(
    mod_args[["vaccination"]],
    doses = 2L,
    population = mod_args[["population"]]
  )

  # add null intervention if this is missing
  # if not missing, check that it conforms to expectations
  if ("intervention" %in% names(mod_args)) {
    # if interventions are passed, check for the types and names
    stopifnot(
      "`intervention` must be a list of <intervention>s" =
        checkmate::test_list(
          mod_args[["intervention"]],
          types = c("contacts_intervention", "rate_intervention"),
          names = "unique"
        )
    )

    # check for any other intervention list element names - allow only
    # primary model parameters and not derived ones
    checkmate::assert_names(
      names(mod_args[["intervention"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "hospitalisation_rate",
        "mortality_rate", "recovery_rate", "contacts"
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
    # on the transmission rate transmissibility
    if (identical(names(mod_args[["intervention"]]), "contacts")) {
      mod_args[["intervention"]]$transmissibility <- no_rate_intervention()
    }
  } else {
    # add as a list element named "contacts", and one named "transmissibility"
    mod_args[["intervention"]] <- list(
      contacts = no_contacts_intervention(
        mod_args[["population"]]
      ),
      # a dummy intervention on the rate parameter transmissibility
      transmissibility = no_rate_intervention()
    )
  }

  # handle time dependence if present and check for allowed primary parameters
  if ("time_dependence" %in% names(mod_args)) {
    checkmate::assert_names(
      names(mod_args[["time_dependence"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "hospitalisation_rate",
        "mortality_rate", "recovery_rate"
      )
    )
  } else {
    mod_args[["time_dependence"]] <- no_time_dependence()
  }

  # return arguments invisibly
  invisible(mod_args)
}

#' @title Prepare arguments for the Vacamole epidemic function
#' @name check_prepare_vacamole_args
#' @rdname check_prepare_vacamole_args
#' @keywords internal
.prepare_args_model_vacamole <- function(mod_args) {
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

  # calculate derived infection parameters

  # modified parameters for the two-dose vaccinated compartment
  transmissibility_vax <- mod_args$transmissibility *
    (1.0 - mod_args[["susc_reduction_vax"]])

  hospitalisation_rate_vax <- mod_args$hospitalisation_rate *
    (1.0 - mod_args[["hosp_reduction_vax"]])

  mortality_rate_vax <- mod_args$mortality_rate *
    (1.0 - mod_args[["mort_reduction_vax"]])

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
    transmissibility = mod_args$transmissibility,
    infectiousness_rate = mod_args$infectiousness_rate,
    hospitalisation_rate = mod_args$hospitalisation_rate,
    mortality_rate = mod_args$mortality_rate,
    recovery_rate = mod_args$recovery_rate,
    transmissibility_vax = transmissibility_vax,
    hospitalisation_rate_vax = hospitalisation_rate_vax,
    mortality_rate_vax = mortality_rate_vax,
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
