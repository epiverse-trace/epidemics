#' @title Check arguments to the diphtheria model function
#' @name check_prepare_diphtheria_args
#' @rdname check_prepare_diphtheria_args
#'
#' @description Check and prepare the four main arguments to
#' [model_diphtheria_cpp()] for use with [.model_diphtheria_cpp()].
#'
#' @return
#'
#' `.check_args_model_diphtheria()` invisibly returns the model arguments
#' passed in `mod_args`. Model functionalities, such as interventions,
#' passed as `NULL` are replaced with dummy values for internal functions.
#'
#' `.prepare_args_model_diphtheria()` returns a list of model arguments
#' suitable for [.model_diphtheria_cpp()]. This is a named list consisting of:
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
#'  - `reporting_rate`: a single number for the proportion of infectious cases
#' reported;
#'
#' - `prop_hosp`: a single number for the proportion of reported cases that
#' need hospitalisation;
#'
#' - `hosp_entry_rate`, `hosp_exit_rate`: two numbers representing the rate of
#' entry and exit from the 'hospitalised' compartment;
#'
#' - `recovery_rate`: a single number for the recovery rate from the infection;
#'
#'  - `time_end`, `increment`: two numbers for the time at which to end the
#' simulation, and the value by which the simulation time
#' is incremented.
#'
#' @keywords internal
#' @details
#' `.prepare_args_model_diphtheria()` prepares arguments for
#' [.model_diphtheria_cpp()], which is the C++ function that solves the default
#' ODE system using a Boost _odeint_ solver.
#'
#' `.prepare_args_model_diphtheria()` converts the arguments collected in
#' `mod_args` into simpler structures such as lists and numeric or integer
#' vectors that can be interpreted as C++ types such as `Rcpp::List`,
#' `Rcpp::NumericVector`, or `Eigen::MatrixXd`.
.check_args_model_diphtheria <- function(mod_args) {
  # check that arguments list has expected names
  checkmate::assert_names(names(mod_args), type = "unique")

  # load number of compartments to check initial conditions matrix
  compartments_diphtheria <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "recovered"
  )

  # input checking on population parameters
  assert_population(
    mod_args[["population"]],
    compartments = compartments_diphtheria
  )

  # add null rate_intervention if this is missing
  # if not missing, check that it conforms to expectations
  if (is.null(mod_args[["intervention"]])) {
    # no need to add contact intervention as this model does not support a
    # contact matrix
    mod_args[["intervention"]] <- list(
      transmissibility = no_rate_intervention()
    )
  } else {
    # length and type of intervention is checked at the top level
    # check for any intervention list element names
    # NOTE: reporting_rate is not allowed as it cannot be meaningfully modified
    # i.e., cannot be increased in an intervention
    checkmate::assert_names(
      names(mod_args[["intervention"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "prop_hosp",
        "hosp_entry_rate", "hosp_exit_rate", "recovery_rate"
      )
    )
  }

  # handle time dependence if present and check for allowed primary parameters
  if (is.null(mod_args[["time_dependence"]])) {
    mod_args[["time_dependence"]] <- no_time_dependence()
  } else {
    checkmate::assert_names(
      names(mod_args[["time_dependence"]]),
      subset.of = c(
        "transmissibility", "infectiousness_rate", "prop_hosp",
        "reporting_rate", "hosp_entry_rate", "hosp_exit_rate", "recovery_rate"
      )
    )
  }

  # return arguments invisibly
  invisible(mod_args)
}

#' @title Prepare arguments for the Vacamole epidemic function
#' @name check_prepare_diphtheria_args
#' @rdname check_prepare_diphtheria_args
#' @keywords internal
.prepare_args_model_diphtheria <- function(mod_args) {
  # prepare initial conditions by scaling with demography
  initial_state <-
    get_parameter(mod_args[["population"]], "initial_conditions") *
      get_parameter(mod_args[["population"]], "demography_vector")

  # modify initial state by proportion vaccinated
  # NOTE: this assumes that the first column is 'susceptible'
  initial_state[, 1] <- initial_state[, 1] * (1 - mod_args[["prop_vaccinated"]])

  # return mod args without population and prop_vaccinated
  mod_args <- mod_args[!names(mod_args) %in% c("population", "prop_vaccinated")]

  # rename interventions to rate_interventions
  names(mod_args)[names(mod_args) == "intervention"] <- "rate_interventions"

  # combine with initial state
  mod_args <- c(
    list(initial_state = initial_state),
    mod_args
  )

  mod_args
}
