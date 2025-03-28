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
#' This model allows for a proportion of each demographic group to be vaccinated
#' at the start of the outbreak, and thus to not contribute to the outbreak.
#' The model also allows for changes to the number of susceptibles in each age
#' group to model influxes or evacuations from camps.
#'
#' @inheritParams model_default
#' @param reporting_rate A numeric for the proportion of infectious cases
#' that is reported; this is a precursor to hospitalisation as only reported
#' cases are hospitalised.
#' @param prop_hosp A numeric for the proportion of reported cases that is
#' hospitalised.
#' @param hosp_entry_rate A numeric for the rate at which reported cases
#' of infectious individuals are hospitalised.
#' This is calculated as 1 / time to hospitalisation, denoted \eqn{\tau_1}.
#' @param hosp_exit_rate A numeric for the rate at which individuals are
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
#' # Details: Model an infection outbreak in a humanitarian camp setting
#'
#' This model has been developed for diphtheria outbreaks in settings where
#' interventions on social contacts are difficult to implement. It it suitable
#' for application to the outbreak of similar, directly transmitted infectious
#' diseases as well.
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
#' - Transmission rate (\eqn{\beta}, `transmission_rate`): 0.8888889, assuming
#' an \eqn{R_0} of 4.0 and a total infectious period of 4.5 days.
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
#' @return A `data.table` with the columns "time", "compartment", "age_group",
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
#' data <- model_diphtheria(
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
#' data <- model_diphtheria(
#'   camp_pop,
#'   prop_vaccinated = prop_vaccinated,
#'   population_change = p
#' )
#' head(data)
#' tail(data)
#' @export
model_diphtheria_cpp <- function(population,
                                 transmission_rate = 4.0 / 4.5,
                                 infectiousness_rate = 1.0 / 3.0,
                                 recovery_rate = 1.0 / 3.0,
                                 reporting_rate = 0.03,
                                 prop_hosp = 0.01,
                                 hosp_entry_rate = 0.2,
                                 hosp_exit_rate = 0.2,
                                 prop_vaccinated = 0.0 *
                                   population[["demography_vector"]],
                                 intervention = NULL,
                                 time_dependence = NULL,
                                 population_change = NULL,
                                 time_end = 100,
                                 increment = 1) {
  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "recovered"
  )
  assert_population(population, compartments)
  
  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
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
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    reporting_rate = reporting_rate,
    prop_hosp = prop_hosp,
    hosp_entry_rate = hosp_entry_rate,
    hosp_exit_rate = hosp_exit_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)
  
  # check that prop_vaccinated is the same length as demography_vector
  # TODO: treat this as a composable element for vectorisation
  checkmate::assert_numeric(
    prop_vaccinated,
    lower = 0, upper = 1,
    len = length(population[["demography_vector"]])
  )
  # convert to list for data.table
  prop_vaccinated <- list(prop_vaccinated)
  
  # Check if `intervention` is a list of interventions or a list-of-lists
  # and convert to a list for a data.table list column. NULL is allowed;
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) && all(
    vapply(
      unlist(intervention, recursive = FALSE),
      FUN = function(x) {
        is_intervention(x) || is.null(x)
      }, TRUE
    )
  )
  
  # Check if parameters can be recycled;
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls
  )
  
  # make lists if not lists
  population <- list(population)
  if (is_lofints) {
    intervention <- list(intervention)
  }
  
  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`, in order as the first two args
  # NOTE: this functionality is not vectorised;
  # convert to list for data.table list column
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    any.missing = FALSE, names = "unique"
  )
  # lapply on null returns an empty list
  invisible(
    lapply(time_dependence, checkmate::assert_function,
           args = c("time", "x"), ordered = TRUE
    )
  )
  time_dependence <- list(
    .cross_check_timedep(
      time_dependence,
      c(
        "transmission_rate", "infectiousness_rate", "prop_hosp",
        "reporting_rate", "hosp_entry_rate", "hosp_exit_rate", "recovery_rate"
      )
    )
  )
  
  # TODO: allow vectorised input for population_change
  # convert to list for data.table
  population_change <- list(population_change)
  
  # collect parameters and add a parameter set identifier
  params <- data.table::as.data.table(params)
  params[, "param_set" := .I]
  
  # this nested data.table will be returned
  model_output <- data.table::CJ(
    population = population,
    prop_vaccinated = prop_vaccinated,
    intervention = intervention,
    time_dependence = time_dependence,
    population_change = population_change,
    increment = increment,
    sorted = FALSE
  )
  
  # process the population, interventions, and vaccinations, after
  # cross-checking them agains the relevant population
  model_output[, args := apply(model_output, 1, function(x) {
    .check_prepare_args_diphtheria(c(x))
  })]
  model_output[, "scenario" := .I]
  
  # combine infection parameters and scenarios
  # NOTE: join X[Y] must have params as X as list cols not supported for X
  model_output <- params[, as.list(model_output), by = names(params)]
  
  # collect model arguments in column data, then overwrite
  model_output[, args := apply(model_output, 1, function(x) {
    c(x[["args"]], x[param_names]) # avoid including col "param_set"
  })]
  model_output[, "data" := Map(population, args, f = function(p, l) {
    .output_to_df(
      do.call(.model_diphtheria_cpp, l),
      population = p, # taken from local scope/env
      compartments = compartments
    )
  })]
  
  # remove temporary arguments
  model_output$args <- NULL
  
  # check for single row output, i.e., scalar arguments, and return data.table
  # do not return the parameters in this case
  if (nrow(model_output) == 1L) {
    model_output <- model_output[["data"]][[1L]] # hardcoded for special case
  }
  
  # return data.table
  model_output[]
}

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
#' This model allows for a proportion of each demographic group to be vaccinated
#' at the start of the outbreak, and thus to not contribute to the outbreak.
#' The model also allows for changes to the number of susceptibles in each age
#' group to model influxes or evacuations from camps.
#'
#' @inheritParams model_default
#' @param reporting_rate A numeric for the proportion of infectious cases
#' that is reported; this is a precursor to hospitalisation as only reported
#' cases are hospitalised.
#' @param prop_hosp A numeric for the proportion of reported cases that is
#' hospitalised.
#' @param hosp_entry_rate A numeric for the rate at which reported cases
#' of infectious individuals are hospitalised.
#' This is calculated as 1 / time to hospitalisation, denoted \eqn{\tau_1}.
#' @param hosp_exit_rate A numeric for the rate at which individuals are
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
#' # Details: Model an infection outbreak in a humanitarian camp setting
#'
#' This model has been developed for diphtheria outbreaks in settings where
#' interventions on social contacts are difficult to implement. It it suitable
#' for application to the outbreak of similar, directly transmitted infectious
#' diseases as well.
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
#' - Transmission rate (\eqn{\beta}, `transmission_rate`): 0.8888889, assuming
#' an \eqn{R_0} of 4.0 and a total infectious period of 4.5 days.
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
#' @return A `data.table` with the columns "time", "compartment", "age_group",
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
#' data <- model_diphtheria(
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
#' data <- model_diphtheria(
#'   camp_pop,
#'   prop_vaccinated = prop_vaccinated,
#'   population_change = p
#' )
#' head(data)
#' tail(data)
#' @export
model_diphtheria <- function(population,
                             transmission_rate = 4.0 / 4.5,
                             infectiousness_rate = 1.0 / 3.0,
                             recovery_rate = 1.0 / 3.0,
                             reporting_rate = 0.03,
                             prop_hosp = 0.01,
                             hosp_entry_rate = 0.2,
                             hosp_exit_rate = 0.2,
                             prop_vaccinated = 0.0 *
                               population[["demography_vector"]],
                             intervention = NULL,
                             time_dependence = NULL,
                             population_change = NULL,
                             time_end = 100,
                             increment = 1) {
  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "hospitalised", "recovered"
  )
  assert_population(population, compartments)
  
  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
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
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    reporting_rate = reporting_rate,
    prop_hosp = prop_hosp,
    hosp_entry_rate = hosp_entry_rate,
    hosp_exit_rate = hosp_exit_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)
  
  # check that prop_vaccinated is the same length as demography_vector
  # TODO: treat this as a composable element for vectorisation
  checkmate::assert_numeric(
    prop_vaccinated,
    lower = 0, upper = 1,
    len = length(population[["demography_vector"]])
  )
  # convert to list for data.table
  prop_vaccinated <- list(prop_vaccinated)
  
  # Check if `intervention` is a list of interventions or a list-of-lists
  # and convert to a list for a data.table list column. NULL is allowed;
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) && all(
    vapply(
      unlist(intervention, recursive = FALSE),
      FUN = function(x) {
        is_intervention(x) || is.null(x)
      }, TRUE
    )
  )
  
  # Check if parameters can be recycled;
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls
  )
  
  # make lists if not lists
  population <- list(population)
  if (is_lofints) {
    intervention <- list(intervention)
  }
  
  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`, in order as the first two args
  # NOTE: this functionality is not vectorised;
  # convert to list for data.table list column
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    any.missing = FALSE, names = "unique"
  )
  # lapply on null returns an empty list
  invisible(
    lapply(time_dependence, checkmate::assert_function,
           args = c("time", "x"), ordered = TRUE
    )
  )
  time_dependence <- list(
    .cross_check_timedep(
      time_dependence,
      c(
        "transmission_rate", "infectiousness_rate", "prop_hosp",
        "reporting_rate", "hosp_entry_rate", "hosp_exit_rate", "recovery_rate"
      )
    )
  )
  
  # TODO: allow vectorised input for population_change
  # convert to list for data.table
  population_change <- list(population_change)
  
  # collect parameters and add a parameter set identifier
  params <- data.table::as.data.table(params)
  params[, "param_set" := .I]
  
  # this nested data.table will be returned
  model_output <- data.table::CJ(
    population = population,
    prop_vaccinated = prop_vaccinated,
    intervention = intervention,
    time_dependence = time_dependence,
    population_change = population_change,
    increment = increment,
    sorted = FALSE
  )
  
  # process the population, interventions, and vaccinations, after
  # cross-checking them agains the relevant population
  model_output[, args := apply(model_output, 1, function(x) {
    .check_prepare_args_diphtheria(c(x))
  })]
  model_output[, "scenario" := .I]
  
  # combine infection parameters and scenarios
  # NOTE: join X[Y] must have params as X as list cols not supported for X
  model_output <- params[, as.list(model_output), by = names(params)]
  
  # collect model arguments in column data, then overwrite
  model_output[, args := apply(model_output, 1, function(x) {
    c(x[["args"]], x[param_names]) # avoid including col "param_set"
  })]
  model_output[, "data" := lapply(args, function(args) {
    time_points <- seq(0, args$time_end, by = args$increment)
    n_time <- length(time_points)
    n_age <- length(args$pop_change_values[[1]])
    
    pop_change <- matrix(0, nrow = n_time, ncol = n_age)
    if(any(args$pop_change_times > 0)){ 
      pop_change[args$pop_change_times + 1,] <- do.call(rbind, args$pop_change_values)
    }
    rate_intervention_start <- as.numeric(args$rate_interventions[[1]]$time_begin)
    rate_intervention_end <- as.numeric(args$rate_interventions[[1]]$time_end)
    rate_intervention_effect <- matrix(rep(args$rate_interventions[[1]]$reduction, n_age), ncol = n_age)
    
    n_rate_intervention <- length(rate_intervention_start)
    
    time_dependent_params <- Map(
      args[names(args$time_dependence)],
      args$time_dependence,
      f = function(x, func) {
        func(time = time_points, x = x)
      }
    )
    # assign time-modified param values
    args[names(time_dependent_params)] <- time_dependent_params
    beta <- args$transmission_rate
    r <- args$reporting_rate
    eta <- args$prop_hosp
    sigma <- args$infectiousness_rate
    tau1 <- args$hosp_entry_rate
    tau2 <- args$hosp_exit_rate
    gamma <- args$recovery_rate
    
    if(length(beta) == 1) beta <- rep(beta, n_time)
    if(length(r) == 1) r <- rep(r, n_time)
    if(length(eta) == 1) eta <- rep(eta, n_time)
    if(length(sigma) == 1) sigma <- rep(sigma, n_time)
    if(length(tau1) == 1) tau1 <- rep(tau1, n_time)
    if(length(tau2) == 1) tau2 <- rep(tau2, n_time)
    if(length(gamma) == 1) gamma <- rep(gamma, n_time)
    
    initial_conditions <- args$initial_state
    init_S <- initial_conditions[, 1]
    init_E <- initial_conditions[, 2]
    init_I <- initial_conditions[, 3]
    init_H <- initial_conditions[, 4]
    init_R <- initial_conditions[, 5]
    
    # Initialize and run the model
    # model <- diphtheria_local$new(
    model <- diphtheria$new(
      time = time_points,
      n_time = n_time,
      n_age = n_age,
      pop_change = pop_change,
      n_rate_intervention = n_rate_intervention,
      beta = beta,
      r = r,
      eta = eta,
      sigma = sigma,
      tau1 = tau1,
      tau2 = tau2,
      gamma = gamma,
      rate_intervention_start = rate_intervention_start,
      rate_intervention_end = rate_intervention_end,
      rate_intervention_effect = rate_intervention_effect,
      init_S = init_S,
      init_E = init_E,
      init_I = init_I,
      init_H = init_H,
      init_R = init_R
    )
    
    result <- model$run(time_points)
    
    # Add scenario information
    dt <- data.table::as.data.table(result)
    # declaring variables below to avoid data.table related lintr messages
    temp <- value <- temp_compartment <- temp_demography <-
      compartment <- demography_group <- `:=` <- time <- NULL
    
    age_group_mappings <- paste0( # properly label demography groups
      seq_len(n_age),
      c(
        rownames(population[[1]]$contact_matrix),
        names(population[[1]]$demography_vector),
        sprintf(
          "demo_group_%i",
          seq_len(nrow(population[[1]]$contact_matrix))
        )
      )[seq_len(nrow(population[[1]]$contact_matrix))]
    )
    names(age_group_mappings) <- seq_len(nrow(population[[1]]$contact_matrix))
    
    mapping <- c( # prepend numbers to help during sorting. Will remove later
      S = "1susceptible", E = "2exposed", I = "3infectious",
      H = "4hospitalised", R = "5recovered", age_group_mappings
    )
    
    # Melt the data table to long format
    data.table::melt(dt,
                     id.vars = "t",
                     variable.name = "temp", # e.g. S[1], ..., V[3]
                     value.name = "value"
    )[ # piping the data.table way. Possible because melt outputs a data.table
      , list(
        time = t, # alternative to using data.table::setnames(dt, "t", "time")
        temp_compartment = substring(temp, 1L, 1L), # e.g. S[1] -> S
        temp_demography = substring(temp, 3L, 3L), # e.g. S[1] -> 1
        value
      )
    ][ # |> the DT way (piping the data.table way)
      ,
      list(
        time,
        demography_group = mapping[temp_demography], # e.g. 1[0,20), 2[20,65),
        compartment = mapping[temp_compartment], # e.g. 1susceptible, 2exposed
        value
      )
    ][ # |> the DT way
      order(time, compartment, demography_group) # prepending numbers helps here
    ][ # |> the DT way
      ,
      `:=`( # used as the prefix form to update multiple columns
        # remove prepended numbers from `mapping`
        demography_group = substring(demography_group, 2L), # e.g. [0,20), ...
        compartment = substring(compartment, 2L) # e.g. susceptible, exposed,
      )
    ][ # |> the DT way
      # added because the previous operation used `:=` which doesn't output
    ]
  })]
  
  # remove temporary arguments
  model_output$args <- NULL
  
  # check for single row output, i.e., scalar arguments, and return data.table
  # do not return the parameters in this case
  if (nrow(model_output) == 1L) {
    model_output <- model_output[["data"]][[1L]] # hardcoded for special case
  }
  
  # return data.table
  model_output[]
}
