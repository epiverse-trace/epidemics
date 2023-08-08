#' Check arguments to epidemic function for the Vacamole model
#' @keywords internal
.check_args_epidemic_vacamole <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique",
    must.include = "vaccination" # vaccination necessary for Vacamole
  )

  # input checking on infection object
  assert_infection(
    mod_args[["infection"]],
    extra_parameters = c(
      "preinfectious_period", "eta", "omega",
      "susc_reduction_vax", "hosp_reduction_vax", "mort_reduction_vax"
    )
  )

  # load number of compartments to check initial conditions matrix
  compartments_vacamole <- read_from_library(
    model_name = "vacamole", what = "compartments"
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
  if (!"intervention" %in% names(mod_args)) {
    mod_args[["intervention"]] <- no_intervention(
      mod_args[["population"]]
    )
  } else {
    assert_intervention(mod_args[["intervention"]], mod_args[["population"]])
  }

  # return arguments invisibly
  invisible(mod_args)
}

#' Prepare arguments for the Vacamole epidemic function
#' @keywords internal
.prepare_args_epidemic_vacamole <- function(mod_args) {
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

  # calculate infection parametersß
  gamma <- 1.0 / mod_args[["infection"]][["infectious_period"]]
  alpha <- 1.0 / mod_args[["infection"]][["preinfectious_period"]]
  beta <- mod_args[["infection"]][["r0"]] /
    mod_args[["infection"]][["infectious_period"]]

  eta <- mod_args[["infection"]][["eta"]]
  omega <- mod_args[["infection"]][["omega"]]

  # modified parameters for the two-dose vaccinated compartment
  beta_v <- beta *
    (1.0 - mod_args[["infection"]][["susc_reduction_vax"]])

  eta_v <- eta *
    (1.0 - mod_args[["infection"]][["hosp_reduction_vax"]])

  omega_v <- omega *
    (1.0 - mod_args[["infection"]][["hosp_reduction_vax"]])

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
    beta, beta_v, alpha, omega, omega_v, eta, eta_v, gamma,
    contact_matrix,
    npi_time_begin, npi_time_end, npi_cr,
    vax_time_begin, vax_time_end, vax_nu,
    mod_args$time_end,
    mod_args$increment
  )
}
