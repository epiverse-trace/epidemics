#' @title Prepare arguments to Vacamole epidemic function
#' @name prepare_vacamole_args
#' @rdname prepare_vacamole_args
#'
#' @description Prepare arguments to [model_vacamole()].
#'
#' @return
#' A list of model arguments suitable for model_vacamole().
#' This is a named list consisting of:
#'
#'  - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `transmission_rate`, `transmission_rate_vax`: two numbers representing the
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
#' `.check_prepare_args_vacamole()` prepares arguments for
#' model_vacamole(), which is the odin based C function that solves the 
#' Vacamole ODE system.
.check_prepare_args_vacamole <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  cmat_init_state <- .prepare_population(mod_args[["population"]])

  # check the interventions list against the population and allowed
  # model parameters
  mod_args[["intervention"]] <- .cross_check_intervention(
    mod_args[["intervention"]], mod_args[["population"]],
    c(
      "contacts",
      "transmission_rate", "infectiousness_rate", "recovery_rate",
      "hospitalisation_rate", "mortality_rate",
      "transmission_rate_vax", "hospitalisation_rate_vax",
      "mortality_rate_vax"
    )
  )
  # check the vaccination against the population
  mod_args[["vaccination"]] <- .cross_check_vaccination(
    mod_args[["vaccination"]], mod_args[["population"]], 2L
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
