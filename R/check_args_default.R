
#' Check arguments to default epidemic function
#' @keywords internal
check_args_default <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique",
    must.include = c(
      "population",
      "r0", "preinfectious_period", "infectious_period",
      "time_end", "increment"
    )
  )

  # add null intervention and vaccination if these are missing
  if (!"intervention" %in% names(mod_args)) {
    mod_args$intervention <- no_intervention(
      mod_args$population
    )
  }
  if (!"vaccination" %in% names(mod_args)) {
    mod_args$vaccination <- no_vaccination(
      mod_args$population
    )
  }

  # some basic input checking for custom classes
  checkmate::assert_class(mod_args$population, "population")
  checkmate::assert_class(mod_args$intervention, "intervention")
  checkmate::assert_class(mod_args$vaccination, "vaccination")

  # input checking on pathogen parameters
  # TODO: move to combined check function for all custom classes
  checkmate::assert_numeric(
    mod_args$r0,
    lower = 0, finite = TRUE,
    # check for as many values as age groups
    len = nrow(mod_args$population$contact_matrix)
  )
  checkmate::assert_numeric(
    mod_args$infectious_period,
    lower = 0, finite = TRUE,
    len = nrow(mod_args$population$contact_matrix)
  )

  # check that compartment sizes are numerics
  checkmate::assert_matrix(
    mod_args$population$initial_conditions,
    mode = "numeric",
    ncols = 5L # hardcoded for the present
  )
  # check that compartments sum to 1.0
  checkmate::assert_numeric(
    rowSums(mod_args$population$initial_condition),
    lower = 1.0 - 1e-6, upper = 1.0 + 1e-6 # specify tolerance manually
  )
  # TODO: more input checking to be added

  # return arguments invisibly
  invisible(mod_args)
}

#' Prepare arguments to default epidemic function
#' @keywords internal
prepare_args_default <- function(mod_args) {
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

  # calculate beta and gamma
  mod_args$gamma <- 1.0 / mod_args$infectious_period
  mod_args$alpha <- 1.0 / mod_args$preinfectious_period
  mod_args$beta <- mod_args$r0 / mod_args$infectious_period
  # nu is passed through vaccination class

  # remove r0, infectious_period, and pre-infectious period
  mod_args$r0 <- NULL
  mod_args$infectious_period <- NULL
  mod_args$preinfectious_period <- NULL

  return(mod_args)
}
