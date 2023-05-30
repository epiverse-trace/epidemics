
#' Check arguments to default epidemic function
#' @keywords internal
.check_args_epidemic_default <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique"
  )

  # add null intervention and vaccination if these are missing
  if (!"intervention" %in% names(mod_args)) {
    suppressMessages(
      mod_args$intervention <- no_intervention(
        mod_args$population
      )
    )
  }
  if (!"vaccination" %in% names(mod_args)) {
    suppressMessages(
      mod_args$vaccination <- no_vaccination(
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
    ncols = 5L # hardcoded for the present
  )
  # check that compartments sum to 1.0
  checkmate::assert_numeric(
    rowSums(mod_args$population$initial_conditions),
    lower = 1.0 - 1e-6, upper = 1.0 + 1e-6 # specify tolerance manually
  )
  # TODO: more input checking to be added

  # return arguments invisibly
  invisible(mod_args)
}

#' Prepare arguments to default epidemic function
#' @keywords internal
.prepare_args_epidemic_default <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  mod_args$population$contact_matrix <-
    mod_args$population$contact_matrix /
      max(Re(eigen(mod_args$population$contact_matrix)$values))

  # scale rows of the contact matrix by the corresponding group population
  mod_args$population$contact_matrix <-
    mod_args$population$contact_matrix /
      mod_args$population$demography_vector

  # prepare initial conditions by scaling with demography
  mod_args$population$initial_conditions <-
    mod_args$population$initial_conditions *
      mod_args$population$demography_vector

  # calculate beta and gamma
  mod_args$gamma <- 1.0 / mod_args$infection$infectious_period
  mod_args$alpha <- 1.0 / mod_args$infection$preinfectious_period
  mod_args$beta <- mod_args$infection$r0 / mod_args$infection$infectious_period
  # nu is passed through vaccination class

  # remove infection object, as parameters are passed as alpha, beta, gamma
  mod_args$infection <- NULL

  return(mod_args)
}
