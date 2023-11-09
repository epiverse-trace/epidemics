#' @title Model leaky, two-dose vaccination in an epidemic using Vacamole
#'
#' @name model_vacamole
#' @rdname model_vacamole
#'
#' @description Simulate an epidemic using the _Vacamole_ model for Covid-19
#' developed at RIVM, the National Institute for Public Health and the
#' Environment in the Netherlands.
#' This model is aimed at estimating the impact of 'leaky' vaccination on an
#' epidemic. See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @inheritParams model_default
#' @param hospitalisation_rate A single number for the hospitalisation rate of
#' infectious individuals.
#' @param mortality_rate A single parameter for the mortality rate of
#' infectious or hospitalised individuals.
#' @param susc_reduction_vax A single value between 0.0 and 1.0
#' giving the reduction in susceptibility of infectious individuals who
#' have received two doses of the vaccine.
#' @param hosp_reduction_vax A single value between 0.0 and 1.0
#' giving the reduction in hospitalisation rate of infectious individuals who
#' have received two doses of the vaccine.
#' @param mort_reduction_vax A single value between 0.0 and 1.0
#' giving the reduction in mortality of infectious and hospitalised individuals
#' who have received two doses of the vaccine.
#' @param intervention An `<intervention>` object representing an optional
#' non-pharmaceutical intervention applied to the population during the
#' epidemic. See [intervention()] for details on constructing interventions with
#' age-specific effects on social contacts, as well as for guidance on how to
#' concatenate multiple overlapping interventions into a single `<intervention>`
#' object.
#' @param vaccination A `<vaccination>` object representing an optional
#' vaccination regime with two doses followed during the course of the
#' epidemic, with a start and end time, and age-specific vaccination rates for
#' each dose. See [vaccination()].
#' @param time_dependence A named list where each name
#' is a model parameter (see `infection`), and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 200 days.
#' @param increment The size of the time increment. Taken as days, with a
#' default value of 1 day.
#' @details
#' This model allows for:
#'
#'  1. A 'hospitalised' compartment along with a hospitalisatio rates;
#'
#'  2. Two doses of vaccination, with 'leaky' protection, i.e., vaccination does
#' not prevent infection completely but allows for a reduction in the infection
#' rate, as well as reduced rates of moving into states considered more serious,
#' such as 'hospitalised' or 'dead'.
#'
#' ## R and Rcpp implementations
#'
#' `model_vacamole_cpp()` is a wrapper function for
#' [.model_vacamole_cpp()], a C++ function that uses Boost _odeint_ solvers
#' to implement the RIVM Vacamole model.
#'
#' `model_vacamole_r()` is a wrapper around the internal function
#' `.ode_model_vacamole()`, which is passed to [deSolve::lsoda()].
#'
#' Both models return equivalent results, but the C++ implementation is faster.
#'
#' `.model_vacamole_cpp()` and `.ode_model_vacamole()` both accept
#' arguments that are created by processing the `population`, `infection`,
#' `intervention` and `vaccination` arguments to the wrapper function into
#' simpler forms.
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates of
#' transition between the 'susceptible' and 'exposed' compartments, between the
#' 'exposed' and 'infectious' compartments, and in the recovery rate.
#'
#' The default values are:
#'
#' - Transmissibility (\eqn{\beta}, `transmissibility`): 0.186, resulting from
#' an \eqn{R_0} = 1.3 and an infectious period of 7 days.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.5, assuming
#' a pre-infectious period of 2 days.
#'
#' - Hospitalisation rate (\eqn{\eta}, `hospitalistion_rate`): 1.0 / 1000,
#' assuming that one in every thousand infectious individuals is hospitalised.
#'
#' - Mortality rate (\eqn{\omega}, `mortality_rate`): 1.0 / 1000,
#' assuming that one in every thousand infectious and hospitalised
#' individuals dies.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.143, assuming an
#' infectious period of 7 days.
#'
#' - Susceptibility reduction from vaccination (`susc_reduction_vax`): 0.2,
#' assuming a 20% reduction in susceptibility for individuals who are doubly
#' vaccinated.
#'
#' - Hospitalisation reduction from vaccination (`hosp_reduction_vax`): 0.2,
#' assuming a 20% reduction in hospitalisation for individuals who are doubly
#' vaccinated.
#'
#' - Mortality reduction from vaccination (`mort_reduction_vax`): 0.2,
#' assuming a 20% reduction in mortality for individuals who are doubly
#' vaccinated.
#'
#' ## Implementing time-dependent parameters
#'
#' Model rates or parameters can be made time-dependent by passing a function
#' which modifies the parameter based on the current ODE simulation time.
#' For example, a function that increases the transmission rate
#' `transmissibility` could be passed as
#' `time_dependence = list(transmissibility = function(time, x) x * time)`.
#' This functionality may be used to model events that are expected to have some
#' effect on model parameters, such as seasonality or annual schedules such as
#' holidays.
#'
#' @return A `data.table` with the columns "time", "compartment", "age_group",
#' "value". The compartments correspond to the compartments of the model
#' chosen with `model`.
#' The current default model has the compartments "susceptible",
#' "vaccinated_one_dose", "vaccinated_two_dose", "exposed",
#' "infectious", "infectious_vaccinated", "hospitalised",
#' "hospitalised_vaccinated", "recovered",  and "dead".
#'
#' @references
#' Ainslie, K. E. C., Backer, J. A., Boer, P. T. de, Hoek, A. J. van,
#' Klinkenberg, D., Altes, H. K., Leung, K. Y., Melker, H. de, Miura, F., &
#' Wallinga, J. (2022). A scenario modelling analysis to anticipate the impact
#' of COVID-19 vaccination in adolescents and children on disease outcomes in
#' the Netherlands, summer 2021. Eurosurveillance, 27(44), 2101090.
#' \doi{10.2807/1560-7917.ES.2022.27.44.2101090}
#'
#' @export
model_vacamole_cpp <- function(population,
                               transmissibility = 1.3 / 7.0,
                               infectiousness_rate = 1.0 / 2.0,
                               hospitalisation_rate = 1.0 / 1000,
                               mortality_rate = 1.0 / 1000,
                               recovery_rate = 1.0 / 7.0,
                               susc_reduction_vax = 0.2,
                               hosp_reduction_vax = 0.2,
                               mort_reduction_vax = 0.2,
                               intervention = NULL,
                               vaccination,
                               time_dependence = NULL,
                               time_end = 100,
                               increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")

  # check model parameters
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(hospitalisation_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(mortality_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(recovery_rate, lower = 0, finite = TRUE)

  checkmate::assert_number(susc_reduction_vax, lower = 0, upper = 1)
  checkmate::assert_number(hosp_reduction_vax, lower = 0, upper = 1)
  checkmate::assert_number(mort_reduction_vax, lower = 0, upper = 1)

  checkmate::assert_class(vaccination, "vaccination")

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population,
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    hospitalisation_rate = hospitalisation_rate,
    mortality_rate = mortality_rate,
    recovery_rate = recovery_rate,
    susc_reduction_vax = susc_reduction_vax,
    hosp_reduction_vax = hosp_reduction_vax,
    mort_reduction_vax = mort_reduction_vax,
    vaccination = vaccination,
    time_end = time_end, increment = increment
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_list(
      intervention,
      types = c("contacts_intervention", "rate_intervention")
    )
    model_arguments[["intervention"]] <- intervention
  }
  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  if (!is.null(time_dependence)) {
    checkmate::assert_list(time_dependence, "function")
    invisible(
      lapply(time_dependence, checkmate::check_function,
        args = c("time", "x"),
        ordered = TRUE
      )
    )
    model_arguments[["time_dependence"]] <- time_dependence
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_model_vacamole(
    .check_args_model_vacamole(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "vaccinated_one_dose",
    "vaccinated_two_dose", "exposed",
    "exposed_vaccinated", "infectious",
    "infectious_vaccinated", "hospitalised",
    "hospitalised_vaccinated", "dead",
    "recovered"
  )

  # run model over arguments
  output <- do.call(.model_vacamole_cpp, model_arguments)

  # prepare output and return
  output_to_df(output, population, compartments)
}

#' Ordinary Differential Equations for the Vacamole Model
#'
#' @description Provides the ODEs for the RIVM Vacamole model in a format that
#' is suitable for passing to [deSolve::lsoda()].
#' See [model_vacamole_r()] for a list of required parameters.
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
  # no input checking, fn is expected to be called only in model_vacamole_r()
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

  # modify parameters
  infection_params <- params[
    c(
      "transmissibility", "transmissibility_vax", "infectiousness_rate",
      "hospitalisation_rate", "hospitalisation_rate_vax",
      "mortality_rate", "mortality_rate_vax", "recovery_rate"
    )
  ]

  infection_params <- intervention_on_rates(
    t = t,
    interventions = params[["rate_interventions"]],
    parameters = infection_params
  )

  # modify the vaccination rate depending on the regime
  # the number of doses is already checked before passing
  # columns are: 1: rate of first dose, 2: rate of second dose
  current_nu <- params[["vax_nu"]] *
    ((params[["vax_time_begin"]] < t) &
      (params[["vax_time_end"]] > t))

  # calculate transitions
  s_to_e <- (params[["transmissibility"]] * y[, 1] *
    contact_matrix_ %*% (y[, 6] + y[, 7]))
  # transitions into the vaccinated compartments
  s_to_v1 <- current_nu[, 1] * y[, 1]
  v1_to_v2 <- current_nu[, 2] * y[, 2]

  # transitions into the exposed compartment
  v1_to_e <- params[["transmissibility"]] * y[, 2] *
    (contact_matrix_ %*% (y[, 6] + y[, 7]))
  v2_to_ev <- params[["transmissibility_vax"]] * y[, 3] *
    (contact_matrix_ %*% (y[, 6] + y[, 7]))

  # transitions into the infectious compartment
  e_to_i <- params[["infectiousness_rate"]] * y[, 4]
  ev_to_iv <- params[["infectiousness_rate"]] * y[, 5]

  # transitions from infectious to hospitalised
  i_to_h <- params[["hospitalisation_rate"]] * y[, 6]
  iv_to_hv <- params[["hospitalisation_rate_vax"]] * y[, 7]

  # transitions from infectious to dead
  i_to_d <- params[["mortality_rate"]] * y[, 6]
  iv_to_d <- params[["mortality_rate_vax"]] * y[, 7]

  # transitions from hospitalied to dead
  h_to_d <- params[["mortality_rate"]] * y[, 8]
  hv_to_d <- params[["mortality_rate_vax"]] * y[, 9]

  # transitions from infectious to recovered
  i_to_r <- params[["recovery_rate"]] * y[, 6]
  iv_to_r <- params[["recovery_rate"]] * y[, 7]

  # transitions from hospitalised to recovered
  h_to_r <- params[["recovery_rate"]] * y[, 8]
  hv_to_r <- params[["recovery_rate"]] * y[, 9]

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

#' @title Model leaky, two-dose vaccination in an epidemic using Vacamole
#'
#' @name model_vacamole
#' @rdname model_vacamole
#'
#' @export
model_vacamole_r <- function(population,
                             transmissibility = 1.3 / 7.0,
                             infectiousness_rate = 1.0 / 2.0,
                             hospitalisation_rate = 1.0 / 1000,
                             mortality_rate = 1.0 / 1000,
                             recovery_rate = 1.0 / 7.0,
                             susc_reduction_vax = 0.2,
                             hosp_reduction_vax = 0.2,
                             mort_reduction_vax = 0.2,
                             intervention = NULL,
                             vaccination,
                             time_dependence = NULL,
                             time_end = 100,
                             increment = 1) {
  # check class on required inputs
  checkmate::assert_class(population, "population")

  # check model parameters
  checkmate::assert_number(transmissibility, lower = 0, finite = TRUE)
  checkmate::assert_number(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(hospitalisation_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(mortality_rate, lower = 0, finite = TRUE)
  checkmate::assert_number(recovery_rate, lower = 0, finite = TRUE)

  checkmate::assert_number(susc_reduction_vax, lower = 0, upper = 1)
  checkmate::assert_number(hosp_reduction_vax, lower = 0, upper = 1)
  checkmate::assert_number(mort_reduction_vax, lower = 0, upper = 1)

  checkmate::assert_class(vaccination, "vaccination")

  # check the time end and increment
  # restrict increment to lower limit of 1e-6
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_number(increment, lower = 1e-6, finite = TRUE)

  # collect population, infection, and model arguments passed as `...`
  model_arguments <- list(
    population = population,
    transmissibility = transmissibility,
    infectiousness_rate = infectiousness_rate,
    hospitalisation_rate = hospitalisation_rate,
    mortality_rate = mortality_rate,
    recovery_rate = recovery_rate,
    susc_reduction_vax = susc_reduction_vax,
    hosp_reduction_vax = hosp_reduction_vax,
    mort_reduction_vax = mort_reduction_vax,
    vaccination = vaccination,
    time_end = time_end, increment = increment
  )

  # check class add intervention and vaccination if not NULL
  if (!is.null(intervention)) {
    checkmate::assert_list(
      intervention,
      types = c("contacts_intervention", "rate_intervention")
    )
    model_arguments[["intervention"]] <- intervention
  }

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`
  # time must be before x, and they must be first two args
  if (!is.null(time_dependence)) {
    checkmate::assert_list(time_dependence, "function")
    invisible(
      lapply(time_dependence, checkmate::check_function,
        args = c("time", "x"),
        ordered = TRUE
      )
    )
    model_arguments[["time_dependence"]] <- time_dependence
  }

  # prepare checked arguments for function
  # this necessary as check_args adds intervention and vaccination
  # if missing
  model_arguments <- .prepare_args_model_vacamole(
    .check_args_model_vacamole(model_arguments)
  )

  # get compartment names
  compartments <- c(
    "susceptible", "vaccinated_one_dose",
    "vaccinated_two_dose", "exposed",
    "exposed_vaccinated", "infectious",
    "infectious_vaccinated", "hospitalised",
    "hospitalised_vaccinated", "dead",
    "recovered"
  )

  # get compartment states over timesteps
  data <- deSolve::lsoda(
    y = model_arguments[["initial_state"]],
    times = seq(0, time_end, increment),
    func = .ode_model_vacamole,
    parms = model_arguments
  )

  # convert to long format using output_to_df() and return
  output_to_df(
    output = list(
      x = data[, setdiff(colnames(data), "time")],
      time = seq(0, time_end, increment)
    ),
    population = population,
    compartments = compartments
  )
}
