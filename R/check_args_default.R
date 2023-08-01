#' Check arguments to default epidemic function
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

#' Prepare arguments to default epidemic function
#' @keywords internal
.prepare_args_epidemic_default <- function(mod_args) {
  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  mod_args[["population"]][["contact_matrix"]] <-
    mod_args[["population"]][["contact_matrix"]] /
      max(Re(eigen(mod_args[["population"]][["contact_matrix"]])$values))

  # scale rows of the contact matrix by the corresponding group population
  mod_args[["population"]][["contact_matrix"]] <-
    mod_args[["population"]][["contact_matrix"]] /
      mod_args[["population"]][["demography_vector"]]

  # prepare initial conditions by scaling with demography
  mod_args[["population"]][["initial_conditions"]] <-
    mod_args[["population"]][["initial_conditions"]] *
      mod_args[["population"]][["demography_vector"]]

  # calculate beta and gamma
  mod_args[["gamma"]] <- 1.0 / mod_args[["infection"]][["infectious_period"]]
  mod_args[["alpha"]] <- 1.0 / mod_args[["infection"]][["preinfectious_period"]]
  mod_args[["beta"]] <- mod_args[["infection"]][["r0"]] /
    mod_args[["infection"]][["infectious_period"]]
  # nu is passed through vaccination class

  # remove infection object, as parameters are passed as alpha, beta, gamma
  mod_args[["infection"]] <- NULL

  return(mod_args)
}
