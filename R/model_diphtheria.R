#' @title Model a diphtheria outbreak using a compartmental ODE model
#'
#' @name model_diphtheria
#' @rdname model_diphtheria
#'
#' @description Simulate a diphtheria outbreak using a deterministic,
#' compartmental ordinary differential equation model with the compartments
#' "susceptible", "exposed", "infectious", "hospitalised", and"recovered".
#' The model is based on Finger et al. (2019) and is intended to be used in the
#' context of internally displaced people (IDP) or refugee camps.
#' This model accommodates age or demographic structure and allows for a
#' proportion of each demographic group to be vaccinated at the start of the
#' outbreak, and thus to not contribute to the outbreak.
#' The model also allows for changes to the number of susceptibles in each age
#' group to model influxes or evacuations from camps.
#'
#' @inheritParams model_default
#' @param reporting_rate A single number for the proportion of infectious cases
#' that is reported; this is a precursor to hospitalisation as only reported
#' cases are hospitalised.
#' @param prop_hosp A single number for the proportion of reported cases that is
#' hospitalised.
#' @param hosp_entry_rate A single number for the rate at which reported cases
#' of infectious individuals are hospitalised.
#' This is calculated as 1 / time to hospitalisation, denoted \eqn{\tau_1}.
#' @param hosp_exit_rate A single number for the rate at which individuals are
#' discharged from hospital to enter the 'recovered' compartment.
#' This is calculated as 1 / time to discharge, denoted \eqn{\tau_2}.
#' @param prop_vaccinated A numeric vector of the same length as the number of
#' demographic groups indicated the proportion of each group that is vaccinated.
#' These individuals are not included in the model dynamics.
#' @param intervention A named list of `<rate_intervention>` objects
#' representing optional pharmaceutical or non-pharmaceutical interventions
#' applied to the model parameters listed above.
#' @param population_change A two-element list, with elements named `"time"` and
#' `"values"`, giving the times of population changes, and the corresponding
#' changes in the population of each demographic group at those times.
#' `"time"` must be a numeric vector, while `"values"` must be a list of the
#' length of `"time"`, with each element a numeric vector of the same length as
#' the number of demographic groups in `population`.
#' @details
#'
#' ## Rcpp implementations
#'
#' `model_diphtheria_cpp()` is a wrapper function for [.model_diphtheria_cpp()],
#' an internal C++ function that uses Boost _odeint_ solvers for an SEIHR model.
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates transitions between
#' compartments per model run.
#'
#' However, model parameters may be passed as numeric vectors. These vectors
#' must follow Tidyverse recycling rules: all vectors must have the same length,
#' or, vectors of length 1 will be recycled to the length of any other vector.
#'
#' The default values are taken from Finger et al. (2019) where possible:
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.8888889, assuming an
#' \eqn{R_0} of 4.0 and a total infectious period of 4.5 days.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.333, assuming
#' a pre-infectious period of 3 days.
#'
#' - Reporting rate (\eqn{r}, `reporting_rate`): 0.03, assuming that 3% of
#' infectious cases are detected or reported.
#'
#' - Proportion hospitalised (\eqn{\eta}, `prop_hosp`): 0.01, assuming that 1%
#' of reported cases need hospital treatment.
#'
#' - Hospital entry rate (\eqn{\tau_1}, `hosp_entry_rate`): 0.2, assuming that
#' it takes 5 days for infectious individuals to seek hospital treatment.
#'
#' - Hospital exit rate (\eqn{\tau_2}, `hosp_exit_rate`): 0.2, assuming that
#' individuals are discharged from hospital after 5 days.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.333, assuming an
#' infectious period following symptoms, of 3 days.
#'
#' ## Modelling population changes
#'
#' This model allows changes to the number of susceptibles in each demographic
#' group, to represent influxes or evacuations from the camp as would be
#' expected in humanitarian relief situations.
#' Users can specify the times and changes (to each demographic group) of
#' changes using the `population_changes` argument, to examine the effect on
#' outbreak dynamics.
#'
#' @references
#' Finger, F., Funk, S., White, K., Siddiqui, M. R., Edmunds, W. J., &
#' Kucharski, A. J. (2019). Real-time analysis of the diphtheria outbreak in
#' forcibly displaced Myanmar nationals in Bangladesh. BMC Medicine, 17, 58.
#' \doi{10.1186/s12916-019-1288-7}.
#' @return A `data.frame` with the columns "time", "compartment", "age_group",
#' "value", and "run", giving the number of individuals per demographic group
#' in each compartment at each timestep in long (or "tidy") format, with "run"
#' indicating the unique parameter combination.
#' @examples
#' # create a dummy camp population with three age groups
#' # diphtheria model is SEIHR
#' # assume that most are susceptible, some infectious
#' # values taken from supplementary material in Finger et al. for the
#' # Kutupalong camp, rounded to the nearest 100
#' n_age_groups <- 3
#' demography_vector <- c(83000, 108200, 224600)
#' initial_conditions <- matrix(0, nrow = n_age_groups, ncol = 5)
#'
#' # set susceptibles and infectious
#' initial_conditions[, 1] <- demography_vector - 1
#' initial_conditions[, 3] <- rep(1, n_age_groups)
#'
#' camp_pop <- population(
#'   contact_matrix = matrix(1, nrow = n_age_groups, ncol = n_age_groups),
#'   demography_vector = demography_vector,
#'   initial_conditions = initial_conditions / demography_vector
#' )
#'
#' # assume younger age groups are vaccinated
#' prop_vaccinated <- c(0.2, 0.10, 0.1)
#'
#' # run model for single, default parameter set
#' data <- model_diphtheria_cpp(
#'   camp_pop,
#'   prop_vaccinated = prop_vaccinated
#' )
#' head(data)
#' tail(data)
#'
#' # run model with increase in population
#' # create population change data
#' p <- list(
#'   time = 70,
#'   values = list(
#'     c(1e4, 2e5, 1e5)
#'   )
#' )
#'
#' data <- model_diphtheria_cpp(
#'   camp_pop,
#'   prop_vaccinated = prop_vaccinated,
#'   population_change = p
#' )
#' head(data)
#' tail(data)
#' @export
model_diphtheria_cpp <- function(population,
                                 transmissibility = 4.0 / 4.5,
                                 infectiousness_rate = 1.0 / 3.0,
                                 recovery_rate = 1.0 / 3.0,
                                 reporting_rate = 0.03,
                                 prop_hosp = 0.01,
                                 hosp_entry_rate = 0.2,
                                 hosp_exit_rate = 0.2,
                                 prop_vaccinated = 0.0 * get_parameter(
                                   population, "demography_vector"
                                 ),
                                 intervention = NULL,
                                 time_dependence = NULL,
                                 population_change = NULL,
                                 time_end = 100,
                                 increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")
  checkmate::assert_numeric(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  # reporting rate and prop_hosp are expected to be proportions bounded 0 - 1
  checkmate::assert_numeric(reporting_rate, lower = 0, upper = 1)
  checkmate::assert_numeric(prop_hosp, lower = 0, upper = 1)
  # hospital entry and exit rate are very likely to be bounded 0 - 1 but
  # are allowed to be > 1 for now
  checkmate::assert_numeric(hosp_entry_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(hosp_exit_rate, lower = 0, finite = TRUE)

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    reporting_rate = reporting_rate,
    prop_hosp = prop_hosp,
    hosp_entry_rate = hosp_entry_rate,
    hosp_exit_rate = hosp_exit_rate,
    time_end = time_end
  )
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      test_recyclable(params)
  )

  # check that prop_vaccinated is the same length as demography_vector
  # TODO: treat this as a composable element for vectorisation
  checkmate::assert_numeric(
    prop_vaccinated,
    lower = 0, upper = 1,
    len = length(get_parameter(population, "demography_vector"))
  )

  # allow only rate interventions in the intervention list
  checkmate::assert_list(
    intervention,
    min.len = 1, names = "unique", any.missing = FALSE,
    types = "rate_intervention", null.ok = TRUE
  )

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    names = "unique", any.missing = FALSE
  )
  invisible(
    lapply(time_dependence, checkmate::check_function,
      args = c("time", "x"),
      ordered = TRUE, null.ok = TRUE
    )
  )

  # check population change and create
  checkmate::assert_list(
    population_change,
    null.ok = TRUE, len = 2L, types = c("numeric", "list")
  )
  if (!is.null(population_change)) {
    checkmate::assert_names(
      names(population_change),
      identical.to = c("time", "values")
    )
  }

  # collect population, infection
  # combine parameters and composable elements into a list of lists of mod args
  params <- .recycle_vectors(params)
  params <- .transpose_base(params)
  model_arguments <- lapply(
    params, function(x) {
      c(
        list(
          population = population,
          prop_vaccinated = prop_vaccinated,
          intervention = intervention,
          time_dependence = time_dependence,
          pop_change_times = population_change[["time"]],
          pop_change_values = population_change[["values"]],
          increment = increment
        ),
        x
      )
    }
  )

  # cross-check model arguments
  # TODO: simplify to only check composable elements
  model_arguments <- lapply(
    model_arguments, function(l) {
      .prepare_args_model_diphtheria(
        .check_args_model_diphtheria(l)
      )
    }
  )

  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "recovered"
  )

  # run model over arguments, prepare output, and return list
  # assign run number as list index - this is the specific combination of
  # parameters
  output <- Map(
    model_arguments, seq_along(model_arguments),
    f = function(l, i) {
      output_ <- .output_to_df(
        do.call(.model_diphtheria_cpp, l),
        population = population,
        compartments = compartments
      )
      output_$run <- i
      output_
    }
  )
  output <- data.table::rbindlist(output)

  # convert to data.frame and return
  # TODO: parameters need to be pulled along
  data.table::setDF(output)[]
}
