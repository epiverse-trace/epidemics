#' @title Prepare arguments to default model function
#' @name prepare_default_args
#' @rdname prepare_default_args
#'
#' @description Prepare arguments to [model_default()].
#'
#' @param mod_args A named list of the population, and epidemic modifiers.
#'
#' @return
#' A list of model arguments suitable for seirv_model().
#' This is a named list consisting of:
#'
#'  - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `transmission_rate`, `infectiousness_rate`, `recovery_rate`: three numbers
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
#' `.check_prepare_args_default()` prepares arguments for
#' seirv_model(), which is the C function that solves the default
#' ODE system using odin.
.check_prepare_args_default <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  cmat_init_state <- .prepare_population(mod_args[["population"]])

  # check the interventions list against the population
  mod_args[["intervention"]] <- .cross_check_intervention(
    mod_args[["intervention"]], mod_args[["population"]],
    c("contacts", "transmission_rate", "infectiousness_rate", "recovery_rate")
  )
  # check the vaccination against the population
  mod_args[["vaccination"]] <- .cross_check_vaccination(
    mod_args[["vaccination"]], mod_args[["population"]], 1L
  )

  # get NPI related times and contact reductions - even if these are dummy
  # interventions
  contact_interventions <- mod_args[["intervention"]][["contacts"]]
  npi_time_begin <- contact_interventions[["time_begin"]]
  npi_time_end <- contact_interventions[["time_end"]]
  npi_cr <- contact_interventions[["reduction"]]

  # get other interventions if any
  rate_interventions <- mod_args[["intervention"]][
    setdiff(names(mod_args[["intervention"]]), "contacts")
  ]

  # get vaccination related times and rates
  vax_time_begin <- mod_args[["vaccination"]][["time_begin"]]
  vax_time_end <- mod_args[["vaccination"]][["time_end"]]
  vax_nu <- mod_args[["vaccination"]][["nu"]]

  # calculate vax_nu as the NUMBER of vaccinations per timestep
  # rather than a fraction/rate. See issue #198 for more details.
  vax_nu <- vax_nu * mod_args[["population"]][["demography_vector"]]

  # remove processed model arguments
  mod_args[c("population", "intervention", "vaccination")] <- NULL
  mod_args[c("param_set", "scenario")] <- NULL

  # return selected arguments for internal C++ function
  c(
    cmat_init_state,
    list(
      npi_time_begin = npi_time_begin, npi_time_end = npi_time_end,
      npi_cr = npi_cr,
      vax_time_begin = vax_time_begin, vax_time_end = vax_time_end,
      vax_nu = vax_nu,
      rate_interventions = rate_interventions
    ),
    mod_args
  )
}
