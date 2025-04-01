#' @title Prepare arguments to diphtheria model function
#' @name prepare_diphtheria_args
#' @rdname prepare_diphtheria_args
#'
#' @description Prepare arguments to [model_diphtheria()] for
#' .model_diphtheria_cpp().
#'
#' @param mod_args A named list of the population, and epidemic modifiers.
#'
#' @return
#' A list of model arguments suitable for .model_diphtheria_cpp().
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
#' - `recovery_rate`: a single number for the recovery rate from the infection;
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
#' - `rate_interventions`: an Rcpp List giving the interventions on model
#' parameters;
#'
#' - `time_dependence`: an Rcpp List giving the time-dependent effects on model
#' parameters in the form of R functions;
#'
#' - `pop_change_times` and `pop_change_values`: the times and values of changes
#' in the population of susceptibles;
#'
#'  - `time_end`, `increment`: two numbers for the time at which to end the
#' simulation, and the value by which the simulation time
#' is incremented.
#'
#' @keywords internal
#' @details
#' `.check_prepare_args_diphtheria()` prepares arguments for
#' .model_diphtheria_cpp(), which is the C++ function that solves the
#' ODE system using a Boost _odeint_ solver, by converting the arguments
#' collected in `mod_args` into simpler structures such as lists and numeric or
#' integer vectors that can be interpreted as C++ types such as `Rcpp::List`,
#' `Rcpp::NumericVector`, or `Eigen::MatrixXd`.
.check_prepare_args_diphtheria <- function(mod_args) {
  # prepare and use the initial state only;
  # modify by reducing susceptibles by `prop_vaccinated`
  cmat_init_state <- .prepare_population(mod_args[["population"]])
  initial_state <- cmat_init_state[["initial_state"]]

  # NOTE: this assumes that the first column is 'susceptible'
  # NOTE: there is no explicit vaccinated compartment
  initial_state[, 1] <- initial_state[, 1] * (1 - mod_args[["prop_vaccinated"]])

  # assign initial state; note prop_vaccinated have been subtracted above
  mod_args[["initial_state"]] <- initial_state

  # cross check population change and handle if NULL
  pop_change_ <- .cross_check_popchange(
    mod_args[["population_change"]], mod_args[["population"]]
  )
  # simplify population_change list
  mod_args[["pop_change_times"]] <- pop_change_[["time"]]
  mod_args[["pop_change_values"]] <- pop_change_[["values"]]

  # return mod args without population and prop_vaccinated
  mod_args <- mod_args[!names(mod_args) %in% c(
    "population", "prop_vaccinated", "population_change"
  )]

  mod_args[["intervention"]] <- .cross_check_intervention(
    mod_args[["intervention"]], mod_args[["population"]],
    c(
      "transmission_rate", "infectiousness_rate", "prop_hosp",
      "hosp_entry_rate", "hosp_exit_rate", "recovery_rate"
    )
  )

  # rename interventions to rate_interventions
  names(mod_args)[names(mod_args) == "intervention"] <- "rate_interventions"

  mod_args
}
