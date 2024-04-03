#' @title Prepare arguments to ebola model function
#' @name prepare_ebola_args
#' @rdname prepare_ebola_args
#'
#' @description Prepare arguments to [model_ebola()] for
#' [.model_ebola_internal()].
#'
#' @param mod_args A named list of composable elements; must include a
#' `<population>`.
#'
#' @return
#' A list of model arguments suitable for [.model_ebola_internal()].
#' This is a named list consisting of:
#'
#' **WIP**
#'
#' @keywords internal
#' @details
#' `.check_prepare_args_ebola()` prepares arguments for
#' [.model_ebola_internal()], by converting some of the arguments
#' collected in `mod_args` into simpler structures that are appropriate for the
#' internal model function.
.check_prepare_args_ebola <- function(mod_args) {
  # prepare the initial conditions; accessed from prepared list
  # NOTE: round to integer as this is a discrete-time model
  initial_state <- round(
    .prepare_population(
      mod_args[["population"]]
    )[["initial_state"]]
  )

  # check that rate interventions are on allowed targets
  # NOTE: type `<rate_intervention>` is already checked at top level
  # NOTE: `infectiousness_rate` and `removal_rate` are not allowed as they
  # control the number of epidemiological sub-compartments, which cannot be
  # safely changed while the model is running
  intervention <- .cross_check_intervention(
    mod_args[["intervention"]], mod_args[["population"]],
    c(
      "transmission_rate", "prop_community", "etu_risk", "funeral_risk"
    )
  )

  # extract and provide time_dependence
  time_dependence <- mod_args[["time_dependence"]]

  # return selected arguments for internal function
  list(
    initial_state = initial_state,
    intervention = intervention, # keep consistent naming
    time_dependence = time_dependence
  )
}
