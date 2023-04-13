#' Model an epidemic
#'
#' @description Simulate an epidemic using a deterministic, compartmental
#' epidemic model from the package's library of epidemic models.
#' The only option currently available is a SEIR-V model, with the compartments
#' "susceptible", "exposed", "infectious", "recovered", and "vaccinated".
#'
#' @param model A string for the epidemic model. The only currently supported
#' option is "default", for the default SEIR-V model.
#' @param ... Arguments to the model specified by `model`. See **Details** for
#' more on the supported arguments.
#'
#' @return A `data.frame` with the columns "time", "compartment", "age_group",
#' "value". The compartments correspond to the compartments of the model
#' chosen with `model`.
#' The current default model has the compartments "susceptible", "exposed",
#' "infectious", "recovered", and "vaccinated".
#' @export
#'
#' @details Arguments passed in `...` that are required for the default moedl:
#' - `population` An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' - `r0` The reproductive number of the infection, \eqn{R_0}.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' - `preinfectious_period` The mean infectious period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' - `infectious_period` The mean infections period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' - `intervention` A non-pharmaceutical intervention applied to the
#' population during the epidemic. See [intervention()].
#' - `vaccination` A vaccination regime followed during the
#' course of the epidemic, with a start and end time, and age-specific effect
#' on the transition of individuals from susceptible to vaccinated.
#' See [vaccination()].
#' - `time_end` The maximum number of timesteps over which to run the model.
#' - `increment` The size of the time increment.
#'
#' @examples
#' # create a population
#' uk_population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0, 0),
#'     nrow = 1, ncol = 5L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' epidemic_cpp(
#'   model = "default",
#'   population = uk_population,
#'   r0 = 1.5,
#'   preinfectious_period = 3,
#'   infectious_period = 7,
#'   time_end = 200,
#'   increment = 1
#' )
epidemic_cpp <- function(model = "default", ...) {

  # collect model arguments passed as `...`
  model_arguments <- list(...)

  # select epidemic model from library
  # currently supports only a single SEIRV model
  # handle the arguments check and prep functions, and the model function
  model <- match.arg(arg = model, several.ok = FALSE)

  # prepare model arguments while checking them
  args_check_fn <- switch(model,
    default = check_args_default
  )
  args_prep_fn <- switch(model,
    default = prepare_args_default
  )
  model_fn <- switch(model,
    default = .epidemic_default_cpp
  )

  # prepare and check model arguments #
  model_arguments <- args_prep_fn(args_check_fn(model_arguments))

  # RUN EPIDEMIC MODEL #
  output <- do.call(
    model_fn,
    model_arguments
  )

  # return combined output
  output_to_df(output)
}
