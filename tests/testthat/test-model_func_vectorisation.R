#### Tests checking the output of multiple parameter sets ####
# Prepare contact matrix and demography vector
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# Prepare some initial objects
uk_population <- population(
  name = "UK population",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(0.9999, 0, 0.0001, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# pass multiple parameter sets
beta <- rnorm(10, 1.3 / 7, sd = 0.01)
sigma <- rnorm(10, 0.5, sd = 0.01)
gamma <- rnorm(10, 1 / 7, sd = 0.01)

test_that("Vectorised inputs: default model, parameter sets", {
  ## Parameter uncertainty, single intervention scenario (none) ##
  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      transmissibility = beta,
      infectiousness_rate = sigma,
      recovery_rate = gamma
    )
  )

  output <- model_default_cpp(
    population = uk_population,
    transmissibility = beta,
    infectiousness_rate = sigma,
    recovery_rate = gamma
  )

  expect_s3_class(output, "data.frame")

  expect_identical(
    nrow(output), length(beta)
  )
  expect_identical(
    output$transmissibility, beta
  )

  # expect parameter set column has different values
  # but scenario has single value
  expect_identical(
    output$param_set, seq_along(beta)
  )
  expect_identical(
    unique(output$scenario), 1L
  )
})

test_that("Vectorised inputs: default model, intervention sets", {
  ## Multiple intervention scenarios ##
  # create scenarios
  school_closure <- intervention(
    "school_closure", "contacts",
    time_begin = 30, time_end = 60, reduction = c(0.2, 0, 0)
  )
  work_closure <- intervention(
    "work_closure", "contacts",
    time_begin = 50, time_end = 80, reduction = c(0, 0.5, 0.1)
  )

  # intervention scenario list
  intervention_scenarios <- list(
    scenario_01 = list(contacts = school_closure),
    scenario_02 = list(contacts = c(school_closure, work_closure)),
    scenario_03 = list(contacts = work_closure)
  )

  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      intervention = intervention_scenarios
    )
  )

  output <- model_default_cpp(
    population = uk_population,
    intervention = intervention_scenarios
  )

  expect_s3_class(output, "data.frame")

  expect_identical(
    nrow(output), length(intervention_scenarios)
  )
  expect_equal(
    output$intervention, intervention_scenarios,
    ignore_attr = TRUE # ignore names
  )

  # expect parameter set column has has single value
  # but scenario has multiple values
  expect_identical(
    output$scenario, seq_along(intervention_scenarios)
  )
  expect_identical(
    unique(output$param_set), 1L
  )
})

test_that("Vectorised inputs: default model, param and intervention sets", {
  ## Multiple intervention scenarios and multiple parameter sets ##
  # create scenarios
  school_closure <- intervention(
    "school_closure", "contacts",
    time_begin = 30, time_end = 60, reduction = c(0.2, 0, 0)
  )
  work_closure <- intervention(
    "work_closure", "contacts",
    time_begin = 50, time_end = 80, reduction = c(0, 0.5, 0.1)
  )

  # intervention scenario list
  intervention_scenarios <- list(
    scenario_01 = list(contacts = school_closure),
    scenario_02 = list(contacts = c(school_closure, work_closure)),
    scenario_03 = list(contacts = work_closure)
  )

  expect_no_condition(
    model_default_cpp(
      transmissibility = beta,
      infectiousness_rate = sigma,
      recovery_rate = gamma,
      population = uk_population,
      intervention = intervention_scenarios
    )
  )
  output <- model_default_cpp(
    transmissibility = beta,
    infectiousness_rate = sigma,
    recovery_rate = gamma,
    population = uk_population,
    intervention = intervention_scenarios
  )

  expect_s3_class(output, "data.frame")

  expect_identical(
    nrow(output), length(intervention_scenarios) * length(beta)
  )

  # expect parameter set and scenario column have multiple values
  expect_identical(
    output$scenario,
    rep(seq_along(intervention_scenarios), length(beta))
  )
  expect_identical(
    output$param_set,
    rep(seq_along(beta), each = length(intervention_scenarios))
  )
})

## Tests for multi parameter, multi interventions with vaccinations ##
test_that("Vectorised inputs: default model multi-param-NPI with vaccination", {
  # create alternative vaccination regimes
  vax_regime_01 <- vaccination(
    time_begin = matrix(20, nrow(contact_matrix)),
    time_end = matrix(100, nrow(contact_matrix)),
    nu = matrix(0.01, nrow(contact_matrix))
  )
  vax_regime_02 <- vaccination(
    time_begin = matrix(10, nrow(contact_matrix)),
    time_end = matrix(40, nrow(contact_matrix)),
    nu = matrix(0.03, nrow(contact_matrix))
  )

  # create scenarios
  school_closure <- intervention(
    "school_closure", "contacts",
    time_begin = 30, time_end = 60, reduction = c(0.2, 0, 0)
  )
  work_closure <- intervention(
    "work_closure", "contacts",
    time_begin = 50, time_end = 80, reduction = c(0, 0.5, 0.1)
  )

  # intervention scenario list
  intervention_scenarios <- list(
    scenario_00 = NULL, # no response scenario
    scenario_01 = list(contacts = school_closure),
    scenario_02 = list(contacts = c(school_closure, work_closure)),
    scenario_03 = list(contacts = work_closure)
  )

  # vaccination scenarios list
  vaccination_scenarios <- list(
    scenario_00 = NULL,
    scenario_01 = vax_regime_01,
    scenario_02 = vax_regime_02
  )

  expect_no_condition(
    model_default_cpp(
      transmissibility = beta,
      infectiousness_rate = sigma,
      recovery_rate = gamma,
      population = uk_population,
      intervention = intervention_scenarios,
      vaccination = vaccination_scenarios
    )
  )

  output <- model_default_cpp(
    transmissibility = beta,
    infectiousness_rate = sigma,
    recovery_rate = gamma,
    population = uk_population,
    intervention = intervention_scenarios,
    vaccination = vaccination_scenarios
  )

  expect_s3_class(output, "data.frame")

  expect_identical(
    nrow(output),
    length(intervention_scenarios) * length(beta) *
      length(vaccination_scenarios)
  )

  # do not expect names
  expect_equal(
    unique(output$vaccination), vaccination_scenarios,
    ignore_attr = TRUE
  )

  # expect parameter set and scenario column have multiple values
  n_unique_scenarios <- length(intervention_scenarios) *
    length(vaccination_scenarios)
  expect_identical(
    output$scenario,
    rep(seq_len(n_unique_scenarios), length(beta))
  )
  expect_identical(
    output$param_set,
    rep(seq_along(beta), each = n_unique_scenarios)
  )
})
