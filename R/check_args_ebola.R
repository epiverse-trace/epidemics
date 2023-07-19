
#' Check arguments to epidemic ebola model from Getz and Dougherty
#' @keywords internal
.check_args_epidemic_ebola <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(
    names(mod_args),
    type = "unique"
  )

  # input checking on the infection object
  assert_infection(
    mod_args[["infection"]],
    extra_parameters = c("shape_E", "shape_I", "rate_E", "rate_I")
  )

  # load number of compartments to check initial conditions matrix
  compartments_ebola <- read_from_library(
    model_name = "ebola", what = "compartments"
  )
  assert_population(
    mod_args[["population"]],
    compartments = compartments_ebola
  )

  # return arguments invisibly
  invisible(mod_args)
}

#' Prepare arguments to epidemic ebola model from Getz and Dougherty
#' @keywords internal
.prepare_args_epidemic_ebola <- function(mod_args) {
  # prepare initial conditions by scaling with demography
  mod_args[["population"]][["initial_conditions"]] <-
    as.integer(mod_args[["population"]][["initial_conditions"]] *
      mod_args[["population"]][["demography_vector"]])
  # convert to integer values
  mod_args[["population"]][["initial_conditions"]] <- as.matrix(
    mod_args[["population"]][["initial_conditions"]],
    nrow = 1L
  )

  # calculate beta, Erlang shape and rate parameters are passed
  mod_args[["beta"]] <- mod_args[["infection"]][["r0"]] /
    mod_args[["infection"]][["infectious_period"]]
  mod_args[["shape_E"]] <- mod_args[["infection"]][["shape_E"]]
  mod_args[["rate_E"]] <- mod_args[["infection"]][["rate_E"]]
  mod_args[["shape_I"]] <- mod_args[["infection"]][["shape_I"]]
  mod_args[["rate_I"]] <- mod_args[["infection"]][["rate_I"]]

  # remove infection object, as parameters are passed as alpha, beta, gamma
  mod_args[["infection"]] <- NULL
  mod_args[["increment"]] <- NULL

  return(mod_args)
}
