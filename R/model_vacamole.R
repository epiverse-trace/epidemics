#' @title Model leaky, two-dose vaccination in an epidemic using Vacamole
#'
#' @name model_vacamole
#' @rdname model_vacamole
#'
#' @description Simulate an epidemic using the _Vacamole_ model for Covid-19
#' developed at RIVM, the National Institute for Public Health and the
#' Environment in the Netherlands.
#' This model is aimed at estimating the impact of 'leaky' vaccination on an
#' epidemic. See **Details** and **References** for more information.
#'
#' @inheritParams model_default
#' @param hospitalisation_rate A numeric for the hospitalisation rate of
#' infectious individuals.
#' @param mortality_rate A numeric for the mortality rate of
#' infectious or hospitalised individuals.
#' @param transmissibility_vax A numeric of values between 0.0 and 1.0
#' giving the transmissibility of the infection to individuals who
#' have received two doses of the vaccine. The default values is 80% of the
#' transmissibility of the infection to individuals who are not doubly
#' vaccinated.
#' @param hospitalisation_rate_vax A numeric of values between 0.0 and 1.0
#' giving the hospitalisation rate of infectious individuals who
#' have received two doses of the vaccine. The default value is 80% of the
#' hospitalisation rate of individuals who are not doubly vaccinated.
#' @param mortality_rate_vax A numeric of values between 0.0 and 1.0
#' giving the mortality of infectious and hospitalised individuals
#' who have received two doses of the vaccine. The default value is 80% of the
#' mortality rate of individuals who are not doubly vaccinated.
#' @param vaccination An optional `<vaccination>` object representing a
#' vaccination regime with **two doses** followed during the course of the
#' epidemic, with a start and end time, and age-specific vaccination rates for
#' each dose. See [vaccination()].
#' @details
#' `model_vacamole_cpp()` is a wrapper function for
#' [.model_vacamole_cpp()], an internal C++ function that uses a Boost _odeint_
#' solvers.
#'
#' The Vacamole model has the compartments "susceptible",
#' "vaccinated_one_dose", "vaccinated_two_dose", "exposed", "infectious"
#' "infectious_vaccinated", "hospitalised", "hospitalised_vaccinated",
#' "recovered",  and "dead".
#'
#' This model allows for:
#'
#'  1. A 'hospitalised' compartment along with a hospitalisation rate;
#'
#'  2. Two doses of vaccination, with 'leaky' protection, i.e., vaccination does
#' not prevent infection completely but allows for a reduction in the infection
#' rate, as well as reduced rates of moving into states considered more serious,
#' such as 'hospitalised' or 'dead'.
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates of
#' transitions between compartments per model run.
#'
#' However, model parameters may be passed as numeric vectors. These vectors
#' must follow Tidyverse recycling rules: all vectors must have the same length,
#' or, vectors of length 1 will be recycled to the length of any other vector.
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.186, resulting from
#' an \eqn{R_0} = 1.3 and an infectious period of 7 days. The transmissibility
#' for doubly vaccinated individuals (\eqn{\beta_v}) is 80% of \eqn{\beta},
#' 0.1488.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.5, assuming
#' a pre-infectious period of 2 days.
#'
#' - Hospitalisation rate (\eqn{\eta}, `hospitalistion_rate`): 1.0 / 1000,
#' assuming that one in every thousand infectious individuals is hospitalised.
#' The hospitalisation rate of doubly vaccinated individuals (\eqn{\eta_v}) is
#' 80% of \eqn{\eta}, 0.8 / 1000.
#'
#' - Mortality rate (\eqn{\omega}, `mortality_rate`): 1.0 / 1000,
#' assuming that one in every thousand infectious and hospitalised
#' individuals dies. The mortality rate of the doubly vaccinated
#' (\eqn{\omega_v}) is 80% of \eqn{\omega}, 0.8 / 1000.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.143, assuming an
#' infectious period of 7 days.
#'
#' @return A `<data.table>`.
#' If the model parameters and composable elements are all scalars, a single
#' `<data.table>` with the columns "time", "compartment", "age_group", and
#' "value", giving the number of individuals per demographic group
#' in each compartment at each timestep in long (or "tidy") format is returned.
#'
#' If the model parameters or composable elements are lists or list-like,
#' a nested `<data.table>` is returned with a list column "data", which holds
#' the compartmental values described above.
#' Other columns hold parameters and composable elements relating to the model
#' run. Columns "scenario" and "param_set" identify combinations of composable
#' elements (population, interventions, vaccination regimes), and infection
#' parameters, respectively.
#'
#' @references
#' Ainslie, K. E. C., Backer, J. A., Boer, P. T. de, Hoek, A. J. van,
#' Klinkenberg, D., Altes, H. K., Leung, K. Y., Melker, H. de, Miura, F., &
#' Wallinga, J. (2022). A scenario modelling analysis to anticipate the impact
#' of COVID-19 vaccination in adolescents and children on disease outcomes in
#' the Netherlands, summer 2021. Eurosurveillance, 27(44), 2101090.
#' \doi{10.2807/1560-7917.ES.2022.27.44.2101090}
#'
#' @examples
#' # create a population, note eleven columns for compartments
#' population <- population(
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0),
#'     nrow = 1, ncol = 11L
#'   )
#' )
#'
#' # create a vaccination regime
#' double_vax <- vaccination(
#'   nu = matrix(1e-3, ncol = 2, nrow = 1),
#'   time_begin = matrix(c(10, 30), nrow = 1),
#'   time_end = matrix(c(50, 80), nrow = 1)
#' )
#'
#' # run epidemic simulation with vaccination but no intervention
#' # with a single set of parameters
#' data <- model_vacamole_cpp(
#'   population = population,
#'   vaccination = double_vax
#' )
#'
#' # view some data
#' head(data)
#'
#' # run epidemic simulation with no vaccination or intervention
#' # and three discrete values of transmissibility
#' data <- model_default_cpp(
#'   population = uk_population,
#'   transmissibility = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
#' )
#'
#' # view some data
#' head(data)
#' tail(data)
#'
#' @export
model_vacamole_cpp <- function(population,
                               transmissibility = 1.3 / 7.0,
                               transmissibility_vax = 0.8 * transmissibility,
                               infectiousness_rate = 1.0 / 2.0,
                               hospitalisation_rate = 1.0 / 1000,
                               hospitalisation_rate_vax = 0.8 *
                                 hospitalisation_rate,
                               mortality_rate = 1.0 / 1000,
                               mortality_rate_vax = 0.8 * mortality_rate,
                               recovery_rate = 1.0 / 7.0,
                               intervention = NULL,
                               vaccination = NULL,
                               time_dependence = NULL,
                               time_end = 100,
                               increment = 1) {
  # TODO: ensure population is properly vectorised
  checkmate::assert_class(population, "population")
  # get compartment names
  compartments <- c(
    "susceptible", "vaccinated_one_dose", "vaccinated_two_dose", "exposed",
    "exposed_vaccinated", "infectious", "infectious_vaccinated", "hospitalised",
    "hospitalised_vaccinated", "dead", "recovered"
  )
  assert_population(population, compartments)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(hospitalisation_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(mortality_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)

  # parameters for rates affecting only doubly vaccinated
  checkmate::assert_numeric(transmissibility_vax, lower = 0, upper = 1)
  checkmate::assert_numeric(hospitalisation_rate_vax, lower = 0, upper = 1)
  checkmate::assert_numeric(mortality_rate_vax, lower = 0, upper = 1)

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # check all vector lengths are equal or 1L
  params <- list(
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    hospitalisation_rate = hospitalisation_rate,
    mortality_rate = mortality_rate,
    recovery_rate = recovery_rate,
    transmissibility_vax = transmissibility_vax,
    hospitalisation_rate_vax = hospitalisation_rate_vax,
    mortality_rate_vax = mortality_rate_vax,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if parameters can be recycled;
  # Check if `population` is a single population or a list of such
  # and convert to list for a data.table list column;
  # also check if `intervention` is a list of interventions or a list-of-lists
  # and convert to a list for a data.table list column. NULL is allowed;
  # Check if `vaccination` is a single vaccination or a list
  # and convert to a list for a data.table list column
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

  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      test_recyclable(params),
    "`population` must be a <population> or a list of <population>s" =
      is_population(population) || checkmate::test_list(
        population,
        types = "population"
      ),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls,
    "`vaccination` must be a <vaccination> or a list of <vaccination>s" =
      is_vaccination(vaccination) || checkmate::test_list(
        vaccination,
        types = c("vaccination", "null"), null.ok = TRUE
      )
  )

  # make lists if not lists
  if (is_population(population)) {
    population <- list(population)
  }
  if (is_lofints) {
    intervention <- list(intervention)
  }
  if (is_vaccination(vaccination) || is.null(vaccination)) {
    vaccination <- list(vaccination)
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
      c("transmissibility", "infectiousness_rate", "recovery_rate")
    )
  )

  # collect parameters and add a parameter set identifier
  params <- data.table::as.data.table(params)
  params[, "param_set" := .I]

  # this nested data.table will be returned
  model_output <- data.table::CJ(
    population = population,
    intervention = intervention,
    vaccination = vaccination,
    time_dependence = time_dependence,
    increment = increment,
    sorted = FALSE
  )

  # process the population, interventions, and vaccinations, after
  # cross-checking them agains the relevant population
  model_output[, args := apply(model_output, 1, function(x) {
    .check_prepare_args_vacamole(c(x))
  })]
  model_output[, "scenario" := .I]

  # combine infection parameters and scenarios
  # NOTE: join X[Y] must have params as X as list cols not supported for X
  model_output <- params[, as.list(model_output), by = names(params)]

  # collect model arguments in column data, then overwrite
  model_output[, args := apply(model_output, 1, function(x) {
    c(x[["args"]], x[param_names]) # avoid including col "param_set"
  })]
  model_output[, data := Map(population, args, f = function(p, l) {
    .output_to_df(
      do.call(.model_vacamole_cpp, l),
      population = p, # taken from local scope/env
      compartments = compartments
    )
  })]

  # check for single row output, i.e., scalar arguments, and return data.frame
  # do not return the parameters in this case
  if (nrow(model_output) == 1L) {
    model_output <- model_output[["data"]][[1L]] # hardcoded for special case
  }

  # return data.table
  model_output[]
}

#' Ordinary Differential Equations for the Vacamole Model
#'
#' @description Provides the ODEs for the RIVM Vacamole model in a format that
#' is suitable for passing to [deSolve::lsoda()].
#' See [model_vacamole_cpp()] for a list of required parameters.
#'
#' @param t A single number of the timestep at which to integrate.
#' @param y The conditions of the epidemiological compartments.
#' @param params The parameters, passed as a named list.
#'
#' @return A list with a vector with as many elements as the number of
#' demographic groups times the number of epidemiological compartments. Each
#' value gives the change in the number of individuals in that compartment.
#' @keywords internal
.ode_model_vacamole <- function(t, y, params) {
  # no input checking, fn is unsafe and not expected to be used
  n_age <- nrow(params[["contact_matrix"]])

  # create a matrix
  # columns correspond to the following compartments
  # 1| 2| 3|4| 5|6| 7|8| 9|10|11 # nolint
  # S|V1|V2|E|EV|I|IV|H|HV|D|R # nolint
  y <- matrix(y, nrow = n_age, ncol = 11L, byrow = FALSE)

  # scale the contact matrix if within the intervention period
  contact_matrix_ <- intervention_on_cm(
    t = t,
    cm = params[["contact_matrix"]],
    time_begin = params[["npi_time_begin"]],
    time_end = params[["npi_time_end"]],
    cr = params[["npi_cr"]]
  )

  # allow modification of only some parameters
  model_params <- params[
    c(
      "transmissibility", "transmissibility_vax", "infectiousness_rate",
      "hospitalisation_rate", "hospitalisation_rate_vax",
      "mortality_rate", "mortality_rate_vax", "recovery_rate"
    )
  ]

  # modifiable parameters
  model_params_primary <- model_params[c(
    "transmissibility", "infectiousness_rate", "hospitalisation_rate",
    "mortality_rate", "recovery_rate"
  )]

  # apply time dependence before interventions
  time_dependent_params <- Map(
    model_params_primary[names(params$time_dependence)],
    params$time_dependence,
    f = function(x, func) {
      func(time = t, x = x)
    }
  )

  # assign time-modified param values
  model_params_primary[names(time_dependent_params)] <- time_dependent_params

  model_params_primary <- intervention_on_rates(
    t = t,
    interventions = params[["rate_interventions"]],
    parameters = model_params_primary
  )

  # reassign all modified values to model_params
  model_params[names(model_params_primary)] <- model_params_primary

  # modify the vaccination rate depending on the regime
  # the number of doses is already checked before passing
  # columns are: 1: rate of first dose, 2: rate of second dose
  current_nu <- params[["vax_nu"]] *
    ((params[["vax_time_begin"]] < t) &
      (params[["vax_time_end"]] > t))

  # calculate transitions
  s_to_e <- (model_params[["transmissibility"]] * y[, 1] *
    contact_matrix_ %*% (y[, 6] + y[, 7]))
  # transitions into the vaccinated compartments
  s_to_v1 <- current_nu[, 1] * y[, 1]
  v1_to_v2 <- current_nu[, 2] * y[, 2]

  # transitions into the exposed compartment
  v1_to_e <- model_params[["transmissibility"]] * y[, 2] *
    (contact_matrix_ %*% (y[, 6] + y[, 7]))
  v2_to_ev <- model_params[["transmissibility_vax"]] * y[, 3] *
    (contact_matrix_ %*% (y[, 6] + y[, 7]))

  # transitions into the infectious compartment
  e_to_i <- model_params[["infectiousness_rate"]] * y[, 4]
  ev_to_iv <- model_params[["infectiousness_rate"]] * y[, 5]

  # transitions from infectious to hospitalised
  i_to_h <- model_params[["hospitalisation_rate"]] * y[, 6]
  iv_to_hv <- model_params[["hospitalisation_rate_vax"]] * y[, 7]

  # transitions from infectious to dead
  i_to_d <- model_params[["mortality_rate"]] * y[, 6]
  iv_to_d <- model_params[["mortality_rate_vax"]] * y[, 7]

  # transitions from hospitalied to dead
  h_to_d <- model_params[["mortality_rate"]] * y[, 8]
  hv_to_d <- model_params[["mortality_rate_vax"]] * y[, 9]

  # transitions from infectious to recovered
  i_to_r <- model_params[["recovery_rate"]] * y[, 6]
  iv_to_r <- model_params[["recovery_rate"]] * y[, 7]

  # transitions from hospitalised to recovered
  h_to_r <- model_params[["recovery_rate"]] * y[, 8]
  hv_to_r <- model_params[["recovery_rate"]] * y[, 9]

  # define compartmental changes
  dS <- -s_to_e - s_to_v1

  dV1 <- s_to_v1 - v1_to_e - v1_to_v2
  dV2 <- v1_to_v2 - v2_to_ev

  dE <- s_to_e + v1_to_e - e_to_i
  dEv <- v2_to_ev - ev_to_iv

  dI <- e_to_i - i_to_h - i_to_r - i_to_d
  dIv <- ev_to_iv - iv_to_hv - iv_to_r - iv_to_d

  dH <- i_to_h - h_to_r - h_to_d
  dHv <- iv_to_hv - hv_to_r - hv_to_d

  dD <- i_to_d + iv_to_d + h_to_d + hv_to_d
  dR <- i_to_r + iv_to_r + h_to_r + hv_to_r

  # return a list
  list(c(dS, dV1, dV2, dE, dEv, dI, dIv, dH, dHv, dD, dR))
}
