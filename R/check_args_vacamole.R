
#' Check arguments to epidemic function for the Vacamole model
#' @keywords internal
.check_args_epidemic_vacamole <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique",
    must.include = c(
      "population", "infection",
      "vaccination", # vaccination necessary for Vacamole
      "time_end", "increment"
    )
  )

  # add null intervention and vaccination if these are missing
  if (!"intervention" %in% names(mod_args)) {
    suppressMessages(
      mod_args$intervention <- no_intervention(
        mod_args$population
      )
    )
  }

  # some basic input checking for custom classes
  checkmate::assert_class(mod_args$population, "population")
  checkmate::assert_class(mod_args$infection, "infection")
  checkmate::assert_class(mod_args$intervention, "intervention")
  checkmate::assert_class(mod_args$vaccination, "vaccination")

  # input checking on pathogen parameters
  checkmate::assert_names(
    names(mod_args$infection),
    must.include = c(
      "r0", "infectious_period", "preinfectious_period",
      "eta", "omega", # hospitalisation and death rate
      "susc_reduction_vax", "hosp_reduction_vax", "mort_reduction_vax"
    )
  )
  checkmate::assert_number(
    mod_args$infection$r0,
    lower = 0, finite = TRUE
  )
  checkmate::assert_number(
    mod_args$infection$infectious_period,
    lower = 0, finite = TRUE
  )
  checkmate::assert_number(
    mod_args$infection$preinfectious_period,
    lower = 0, finite = TRUE
  )

  # check that compartment sizes are numerics
  checkmate::assert_matrix(
    mod_args$population$initial_conditions,
    mode = "numeric",
    ncols = 11L # hardcoded for the Vacamole model with hospitalisation no ICU
  )
  # check that compartments sum to 1.0
  checkmate::assert_numeric(
    rowSums(mod_args$population$initial_condition),
    lower = 1.0 - 1e-5, upper = 1.0 + 1e-5 # specify tolerance manually
  )

  # return arguments invisibly
  invisible(mod_args)
}

#' Prepare arguments for the Vacamole epidemic function
#' @keywords internal
.prepare_args_epidemic_vacamole <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  mod_args$population$contact_matrix <-
    mod_args$population$contact_matrix /
      max(Re(eigen(mod_args$population$contact_matrix)$value))

  # scale rows of the contact matrix by the corresponding group population
  mod_args$population$contact_matrix <-
    mod_args$population$contact_matrix /
      mod_args$population$demography_vector

  # prepare initial conditions by scaling with demography
  mod_args$population$initial_conditions <-
    mod_args$population$initial_conditions *
      mod_args$population$demography_vector

  # calculate infection parameters
  mod_args$gamma <- 1.0 / mod_args$infection$infectious_period
  mod_args$alpha <- 1.0 / mod_args$infection$preinfectious_period
  mod_args$beta <- mod_args$infection$r0 / mod_args$infection$infectious_period

  mod_args$eta <- mod_args$infection$eta
  mod_args$omega <- mod_args$infection$omega

  # modified parameters for the two-dose vaccinated compartment
  mod_args$beta_v <- mod_args$beta *
    (1.0 - mod_args$infection$susc_reduction_vax)

  mod_args$eta_v <- mod_args$infection$eta *
    (1.0 - mod_args$infection$hosp_reduction_vax)

  mod_args$omega_v <- mod_args$infection$omega *
    (1.0 - mod_args$infection$hosp_reduction_vax)

  # nu is passed through vaccination class

  # remove infection object, as parameters are passed as alpha, beta, gamma
  mod_args$infection <- NULL

  return(mod_args)
}
