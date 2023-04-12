#' Model an SIRV epidemic
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
#' @param R0 The reproductive number of the infection. May be a single number or
#' a vector of the same length as the number of demographic groups.
#' @param infectious_period The mean infections period. May be a single number
#' or a vector of the same length as the number of demographic groups.
#' @param intervention An object of the `intervention` class that specifies the
#' start and end of the intervention, and the effect of the intervention on
#' each demographic group. See [intervention()].
#' @param nu The vaccination rate \eqn{\nu}. May be a single number or
#' a vector of the same length as the number of demographic groups.
#' @param t_vax_begin The time at which vaccination begins. Must be a number,
#' less than `t_max`. Defaults to 60% of the maximum time.
#' @param t_vax_end The time at which vaccination ends. Defaults to 80% of the
#' maximum time.
#' @param t_max The maximum number of timesteps over which to run the model.
#' @param t_increment The size of the time increment.
#' @param model Epidemic model function.
#'
#' @return A `data.frame` with the columns "time", "compartment", "age_group",
#' "value". The compartments are "susceptible", "infected", "recovered", and
#' "vaccinated".
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
#' # set up intervention
#' new_intervention <- intervention(
#'   name = "close schools",
#'   time_begin = 50,
#'   time_end = 80,
#'   contact_reduction = c(0.2)
#' )
#'
#' # run epidemic simulation
#' epi_demic(
#'   population = population,
#'   R0 = 1.5, infectious_period = 7,
#'   intervention = new_intervention,
#'   nu = 0.01,
#'   t_vax_begin = 100,
#'   t_vax_end = 200,
#'   t_max = 200,
#'   t_increment = 1
#' )
epi_demic <- function(population,
                      model = "default",
                      R0 = 1.5,
                      infectious_period = 7,
                      intervention,
                      nu = 0.01,
                      t_vax_begin = t_max * 0.6,
                      t_vax_end = t_max * 0.8,
                      t_max = 200,
                      t_increment = 1) {

  # some basic input checking
  checkmate::assert_class(population, "population")
  checkmate::assert_numeric(
    t_vax_begin,
    lower = 0, upper = t_max, finite = TRUE
  )
  checkmate::assert_numeric(R0, lower = 0, finite = TRUE)
  checkmate::assert_numeric(
    infectious_period,
    lower = 0, upper = t_max, finite = TRUE
  )
  checkmate::assert_numeric(nu, lower = 0, finite = TRUE)

  # check that compartment sizes are numerics
  checkmate::assert_matrix(population$initial_conditions,
    mode = "numeric",
    ncols = 4
  )
  # check that compartments sum to 1.0
  checkmate::assert_numeric(
    apply(population$initial_condition, 1, sum),
    lower = 1, upper = 1
  )

  # check for equality of argument lengths
  stopifnot(
    "`R0` length must be the same length as the demography vector" =
      length(R0) == length(population$demography_vector)
  ) # more input checking to be added

  # scale contact matrix
  contact_matrix <- population$contact_matrix
  demography_vector <- population$demography_vector

  contact_matrix <- contact_matrix / max(Re(eigen(contact_matrix)$value))
  contact_matrix <- contact_matrix / demography_vector

  # calculate beta and gamma
  gamma <- 1 / infectious_period
  beta <- R0 / (infectious_period)

  # count age groups
  n_age_groups <- nrow(contact_matrix)

  # make a parameter list for pre vaccination
  params <- list(
    beta = beta, gamma = gamma, nu = nu,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    n_compartments = 4,
    t_vax_begin = t_vax_begin,
    t_npi_begin = intervention$time_begin,
    t_npi_end = intervention$time_end,
    npi_contact_reduction = intervention$contact_reduction
  )

  # make a vector of initial conditions times demography
  init <- population$initial_conditions * demography_vector
  init <- as.vector(init)

  # define times
  times <- seq(0, t_max, by = t_increment)

  # handle the model function
  model_fn <- switch(model,
    default = epidemic_default
  )

  # solve the equations with deSolve for times 0 to t_vax_begin
  output <- deSolve::lsoda(
    init,
    times = times, func = model_fn, parms = params
  )

  # return combined output
  output <- as.data.frame(output)
  # get column data
  compartment <- rep(c("susceptible", "infected", "recovered", "vaccinated"),
    each = nrow(output) * n_age_groups
  )
  age_group <- rep(rep(seq(n_age_groups), each = nrow(output)), 4) # hardcoded
  output_v <- unlist(as.vector(output[, -1]))

  # make long form data
  data <- data.frame(
    time = rep(times, n_age_groups * 4),
    compartment = compartment,
    age_group = age_group,
    value = output_v
  )

  # return data.frame of outputs
  data
}
