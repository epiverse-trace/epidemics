#' @title Check arguments to epidemic ebola model from Getz and Dougherty
#' @name check_prepare_ebola_args
#' @rdname check_prepare_ebola_args
#'
#' @description Check and prepare the four main arguments to
#' [model_ebola_cpp()] for use with [.model_ebola_cpp()].
#' @return
#'
#' `.check_args_model_ebola()` invisibly returns the model arguments passed
#' in `mod_args`.
#'
#' `.prepare_args_model_ebola()` returns a list of model arguments suitable
#' for [.model_ebola_cpp()]. This is a named list consisting of:
#'
#' - `initial_state`: the initial conditions modified to represent absolute
#' rather than proportional values;
#'
#'  - `population_size`: the population size;
#'
#'  - `beta`: the transmission rate of the infection;
#'
#'  - `shape_E`, `rate_E`: the shape and rate parameters of an Erlang
#' distribution that represents passage times through the 'exposed' compartment;
#'
#'  - `shape_I`, `rate_I`: the shape and rate parameters of an Erlang
#' distribution that represents passage times through the 'infectious'
#' compartment;
#'
#' - `time_end`: the time at which the simulation ends.
#'
#' @keywords internal
.check_args_model_ebola <- function(mod_args) {
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
  compartments_ebola <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "funeral", "removed"
  )

  assert_population(
    mod_args[["population"]],
    compartments = compartments_ebola
  )

  # return arguments invisibly
  invisible(mod_args)
}

#' @title Prepare arguments to epidemic ebola model from Getz and Dougherty
#' @name check_prepare_ebola_args
#' @rdname check_prepare_ebola_args
#' @keywords internal
.prepare_args_model_ebola <- function(mod_args) {
  # prepare initial conditions by scaling with demography
  # we take the ceiling to ensure that models with only one
  # infectious individual can be accommodated
  initial_state <-
    ceiling(get_parameter(mod_args[["population"]], "initial_conditions") *
      get_parameter(mod_args[["population"]], "demography_vector"))

  # prepare population size
  population_size <- sum(
    get_parameter(mod_args[["population"]], "demography_vector")
  )

  # calculate beta, Erlang shape and rate parameters are passed
  beta <- get_parameter(mod_args[["infection"]], "r0") /
    get_parameter(mod_args[["infection"]], "infectious_period")
  shape_E <- get_parameter(mod_args[["infection"]], "shape_E")
  rate_E <- get_parameter(mod_args[["infection"]], "rate_E")
  shape_I <- get_parameter(mod_args[["infection"]], "shape_I")
  rate_I <- get_parameter(mod_args[["infection"]], "rate_I")

  # return only required list elements
  list(
    initial_state = initial_state, population_size = population_size,
    beta = beta, shape_E = shape_E, rate_E = rate_E,
    shape_I = shape_I, rate_I = rate_I,
    time_end = mod_args$time_end
  )
}
