#' Model an epidemic
#'
#' @description Simulate an epidemic using a deterministic, compartmental
#' epidemic model from the package's library of epidemic models.
#' The only option currently available is a SEIR-V model, with the compartments
#' "susceptible", "exposed", "infectious", "recovered", and "vaccinated".
#'
#' Each call to `epidemic()` must provide a `population` and an `infection`
#' object, but all other arguments are flexible, and may be passed as `...`.
#' See **Details** for more information.
#'
#' @param model_name A string for the epidemic model. The only currently
#' supported option is "default", for the default SEIR-V model.
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param infection An `infection` object created using [infection()]. Must
#' have at least the basic reproductive number \eqn{R_0} of the infection, and
#' the infectious period. Parameters required by other models can be found in
#' the documentation for model functions.
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
#' @details Arguments passed in `...` may differ depending on the model
#' specified in `model`. The default model (`model_name = "default"`) takes the
#' following arguments.
#' - `population` An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' - `infection` An `infection` class object giving parameters appropriate to
#' an SEIR-V model. See [infection()]. The `infection` object must have the
#' parameters:
#'  - `r0` The reproductive number of the infection, \eqn{R_0}.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#'  - `preinfectious_period` The mean infectious period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#'  - `infectious_period` The mean infections period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' - `intervention` An optional non-pharmaceutical intervention applied to the
#' population during the epidemic. See [intervention()]. This is an optional
#' argument in the default model.
#' - `vaccination` An optional vaccination regime followed during the
#' course of the epidemic, with a start and end time, and age-specific effect
#' on the transition of individuals from susceptible to vaccinated.
#' See [vaccination()]. This is an optional argument in the default model.
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
#' # specify the infection parameters
#' pandemic_influenza <- infection(
#'   r0 = 1.5, infectious_period = 7, preinfectious_period = 3
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' epidemic(
#'   model_name = "default",
#'   population = uk_population,
#'   infection = pandemic_influenza,
#'   time_end = 200,
#'   increment = 1
#' )
epidemic <- function(model_name = "default",
                     population,
                     infection,
                     ...) {

  # select epidemic model from library
  # currently supports only a single SEIRV model
  # handle the arguments check and prep functions, and the model function
  model_name <- match.arg(arg = model_name, several.ok = FALSE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population, infection = infection, ...
  )

  # prepare model arguments while checking them
  args_check_fn <- read_from_library(
    model_type = "epidemic",
    model_name = model_name,
    what = "model_args_checker"
  )
  args_prep_fn <- read_from_library(
    model_type = "epidemic",
    model_name = model_name,
    what = "model_args_prepper"
  )
  model_fn <- read_from_library(
    model_type = "epidemic",
    model_name = model_name,
    what = "model_function"
  )
  compartments <- read_from_library(
    model_type = "epidemic",
    model_name = model_name,
    what = "compartments"
  )

  # prepare and check model arguments #
  model_arguments <- do.call(
    args_prep_fn,
    list(
      do.call(
        args_check_fn,
        list(model_arguments)
      )
    )
  )

  # RUN EPIDEMIC MODEL #
  output <- do.call(
    model_fn,
    model_arguments
  )

  # return combined output
  output_to_df(output, model_arguments, compartments)
}
