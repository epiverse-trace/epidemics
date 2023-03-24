#' Model an SEIR epidemic
#'
#' @description Function for a deterministic susceptible-infectious-recovered
#' model with an optional vaccination component.
#' Allows heterogeneity in social contact patterns, the demography distribution,
#' and key epidemiological parameters: the reproductive number \eqn{R_0}, and
#' the infectious period \eqn{1/\gamma}.
#' Also allows for group-specific initial proportions in each model compartment,
#' as well as group-specific vaccination start dates and vaccination rates.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param r0 The reproductive number of the infection. Must be a vector of the
#' same length as the number of demographic groups.
#' @param preinfectious_period The mean infections period. Must be a vector
#' of the same length as the number of demographic groups.
#' @param infectious_period The mean infections period. Must be a vector of the
#' same length as the number of demographic groups.
#' @param time_end The maximum number of timesteps over which to run the model.
#' @param increment The size of the time increment.
#'
#' @return A `data.frame` with the columns "time", "compartment", "age_group",
#' "value". The comparments are "susceptible", "exposed", "infectious", and
#' "recovered".
#' @export
#'
#' @examples
#' # create a population
#' population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0),
#'     nrow = 1, ncol = 4
#'   )
#' )
#'
#' # run epidemic simulation
#' epidemic_cpp(
#'   population = population,
#'   r0 = rep(1.5, nrow(population$contact_matrix)),
#'   preinfectious_period = rep(3, nrow(population$contact_matrix)),
#'   infectious_period = rep(7, nrow(population$contact_matrix)),
#'   time_end = 200,
#'   increment = 1
#' )
epidemic_cpp <- function(population,
                         r0 = 1.5,
                         preinfectious_period = 3,
                         infectious_period = 7,
                         time_end = 200,
                         increment = 1) {

  # some basic input checking
  checkmate::assert_class(population, "population")
  checkmate::assert_numeric(
    r0,
    lower = 0, finite = TRUE,
    # check for as many values as age groups
    min.len = nrow(population$contact_matrix),
    max.len = nrow(population$contact_matrix)
  )
  checkmate::assert_numeric(
    infectious_period,
    lower = 0, upper = time_end, finite = TRUE,
    min.len = nrow(population$contact_matrix),
    max.len = nrow(population$contact_matrix)
  )

  # check that compartment sizes are numerics
  checkmate::assert_matrix(population$initial_conditions,
    mode = "numeric",
    ncols = 4L # hardocoded
  )
  # check that compartments sum to 1.0
  checkmate::assert_numeric(
    apply(population$initial_condition, 1, sum),
    lower = 1.0, upper = 1.0
  )
  # more input checking to be added

  # scale contact matrix and initial conditions within the population
  population$contact_matrix <- population$contact_matrix /
    max(Re(eigen(population$contact_matrix)$value))
  population$contact_matrix <- population$contact_matrix /
    population$demography_vector
  population$initial_conditions <- population$initial_conditions *
    population$demography_vector

  # calculate beta and gamma
  gamma <- 1.0 / infectious_period
  alpha <- 1.0 / preinfectious_period
  beta <- r0 / infectious_period

  # RUN EPIDEMIC MODEL #
  output <- .epidemic_default_cpp(
    population = population,
    beta = beta,
    alpha = alpha,
    gamma = gamma,
    time_end = time_end, increment = increment
  )

  # return combined output
  output_to_df(output)
}
