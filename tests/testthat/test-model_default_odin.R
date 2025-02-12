#### Tests for the default model ####
# Prepare contact matrix and demography vector
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 60),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# make initial conditions - order is important
initial_conditions <- c(
  S = 1 - 1e-6, E = 0,
  I = 1e-6, R = 0, V = 0
)
initial_conditions <- rbind(
  initial_conditions,
  initial_conditions
)

# create a population
uk_population <- population(
  name = "UK population",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

# prepare a two dose vaccination regime for three age groups
single_vaccination <- vaccination(
  name = "double_vaccination",
  nu = matrix(1e-3, nrow = 2),
  time_begin = matrix(0, nrow = 2),
  time_end = matrix(100, nrow = 2)
)

# model run time
time_end <- 100L
compartments <- c(
  "susceptible", "exposed", "infectious", "recovered", "vaccinated"
)

test_that("Default odin model: basic expectations, scalar arguments", {
  # expect run with no conditions for default arguments
  expect_no_condition(model_default_odin(uk_population))

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_default_odin(uk_population)
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)
  expect_named(
    data, c("time", "demography_group", "compartment", "value"),
    ignore.order = TRUE
  )
  expect_identical(
    nrow(data),
    length(demography_vector) * (time_end + 1L) * length(compartments)
  )
  expect_identical(unique(data$compartment), compartments)
  expect_true(
    checkmate::test_numeric(
      data$value,
      upper = max(demography_vector), lower = 0, any.missing = FALSE
    )
  )
  expect_identical(
    unique(data$demography_group), rownames(contact_matrix)
  )

  # expect no individuals are vaccinated as vaccination is optional
  expect_identical(
    unique(data[grepl("vaccinated", data$compartment, fixed = TRUE), ]$value), 0
  )

  # expect constant population size overall and per demography-group
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(contact_matrix)
  )
  expect_identical(
    rowSums(final_state), uk_population$demography_vector,
    tolerance = 1e-6
  )
})

# NOTE: statistical correctness is not expected to change for vectorised input
test_that("Default odin model: statistical correctness, parameters", {
  # expect final size increases with transmission_rate
  size_beta_low <- epidemic_size(
    model_default_odin(uk_population, transmission_rate = 1.3 / 7.0)
  )
  size_beta_high <- epidemic_size(
    model_default_odin(uk_population, transmission_rate = 1.5 / 7.0)
  )
  expect_true(
    all(size_beta_high > size_beta_low)
  )

  # expect final size increases with infectiousness rate (lower incubation time)
  size_sigma_low <- epidemic_size(
    model_default_odin(uk_population, infectiousness_rate = 1 / 5)
  )
  size_sigma_high <- epidemic_size(
    model_default_odin(uk_population, infectiousness_rate = 1 / 2)
  )
  expect_true(
    all(size_sigma_high > size_sigma_low)
  )

  # expect final size increases with initial infections
  initial_conditions_high <- c(
    S = 1 - 10e-6, E = 0, I = 10e-6,
    R = 0, V = 0
  )
  initial_conditions_high <- rbind(
    initial_conditions_high,
    initial_conditions_high
  )
  uk_population_high_infections <- population(
    name = "UK population",
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions_high
  )
  size_infections_low <- epidemic_size(model_default_odin(uk_population))
  size_infections_high <- epidemic_size(
    model_default_odin(uk_population_high_infections)
  )
  expect_true(
    all(size_infections_high > size_infections_low)
  )
})

# prepare baseline for comparison of against intervention scenarios
data_baseline <- model_default_odin(uk_population)

test_that("Default odin model: contacts interventions and stats. correctness", {
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
  )
  # repeat some basic checks from default case with no intervention
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_default_odin(
      uk_population,
      intervention = list(contacts = intervention)
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_default_odin(
    uk_population,
    intervention = list(contacts = intervention)
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)

  # expect final size is lower with intervention
  expect_true(
    all(epidemic_size(data_baseline) > epidemic_size(data))
  )

  # expect model runs with multiple contacts interventions
  # expect that effect of multiple interventions is greater than single
  intervention_02 <- intervention(
    "work_closure", "contacts", 0, time_end, c(0.1, 0.5)
  )
  combined_interventions <- c(intervention, intervention_02)

  expect_no_condition(
    model_default_odin(
      uk_population,
      intervention = list(contacts = combined_interventions)
    )
  )
  data_combined <- model_default_odin(
    uk_population,
    intervention = list(contacts = combined_interventions)
  )
  # expect epidemic size is lower for combined intervention
  # Epidemic size differs depending on model_default() or model_default_odin().
  # Not sure how large of a discrepancy is expected.
  expect_true(
    all(epidemic_size(data_combined) < epidemic_size(data))
  )
})

test_that("Default odin model: rate interventions", {
  intervention_01 <- intervention(
    "mask_mandate", "rate", 0, time_end, 0.5
  )
  intervention_02 <- intervention(
    "mask_mandate", "rate", time_end / 2, time_end, 0.1
  )
  intervention <- c(intervention_01, intervention_02)
  # repeat some basic checks from default case with no intervention
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_default_odin(
      uk_population,
      intervention = list(transmission_rate = intervention)
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_default_odin(
    uk_population,
    intervention = list(transmission_rate = intervention)
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)

  # expect final size is lower with intervention
  expect_true(
    all(epidemic_size(data_baseline) > epidemic_size(data))
  )
})

# test_that("Default model: vaccination and stats. correctness", {
#   # repeat some basic checks from default case with no vaccination
#   # expect run with no conditions for default arguments
#   expect_no_condition(
#     model_default_odin(uk_population, vaccination = single_vaccination)
#   )

#   # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
#   data <- model_default_odin(uk_population, vaccination = single_vaccination)
#   expect_s3_class(data, "data.frame")
#   expect_identical(length(data), 4L)

#   # expect non-zero vaccinations towards simulation end
#   checkmate::expect_numeric(
#     tail(data[grepl("dose", data$compartment, fixed = TRUE), ]$value),
#     lower = 10
#   )

#   # expect final size is lower with intervention
#   expect_true(
#     all(epidemic_size(data_baseline) > epidemic_size(data))
#   )

#   # expect that high vaccination rates (± 1% per day, e.g. Covid-19 vax)
#   # give statistically correct results (no negative susceptibles at any time)
#   high_rate_vax <- vaccination(
#     time_begin = matrix(200, nrow(contact_matrix)),
#     time_end = matrix(200 + 150, nrow(contact_matrix)),
#     nu = matrix(0.01, nrow = nrow(contact_matrix))
#   )
#   data <- model_default_odin(
#     uk_population,
#     vaccination = high_rate_vax, time_end = 600
#   )
#   checkmate::expect_numeric(
#     data[grepl("susceptible", data$compartment, fixed = TRUE), ]$value,
#     lower = 0
#   )
# })

# test_that("Default model: time dependence", {
#   # expect time dependence is correctly handled
#   time_dependence <- list(
#     transmission_rate = function(time, x, t_change = time_end / 2) {
#       ifelse(time > t_change, x / 2, x)
#     },
#     recovery_rate = function(time, x, t_change = time_end / 2) {
#       ifelse(time > t_change, x + x / 2, x)
#     }
#   )

#   # repeat some basic checks from default case with no time_dependence
#   # expect run with no conditions for default arguments
#   expect_no_condition(
#     model_default_odin(
#       uk_population,
#       time_dependence = time_dependence
#     )
#   )

#   # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
#   data <- model_default_odin(
#     uk_population,
#     time_dependence = time_dependence
#   )
#   expect_s3_class(data, "data.frame")
#   expect_identical(length(data), 4L)

#   # expect final size is lower with intervention
#   expect_true(
#     all(epidemic_size(data_baseline) > epidemic_size(data))
#   )
# })

test_that("Default odin model: errors and warnings, scalar arguments", {
  # expect errors on basic input checking
  expect_error(
    model_default_odin(population = "population"),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  expect_error(
    model_default_odin(population = population),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  pop_wrong_compartments <- uk_population
  pop_wrong_compartments$initial_conditions <- initial_conditions[, -1]
  expect_error(
    model_default_odin(pop_wrong_compartments),
    regexp = "(Assertion on)*(initial_conditions)*failed"
  )

  # expect errors for infection parameters
  expect_error(
    model_default_odin(uk_population, transmission_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_default_odin(uk_population, infectiousness_rate = list(0.2)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_default_odin(uk_population, recovery_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )

  # expect error on time parameters
  expect_error(
    model_default_odin(uk_population, time_end = "100"),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_default_odin(uk_population, time_end = 100.5),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_default_odin(uk_population, time_end = c(100, -100, 10)),
    regexp = "(Element)*(is not >= 0)"
  )
  expect_error(
    model_default_odin(uk_population, increment = "0.1"),
    regexp = "Must be of type 'number'"
  )
  expect_error(
    model_default_odin(uk_population, increment = c(0.1, 0.2)),
    regexp = "Must have length 1"
  )

  # expect error on poorly specified interventions
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, 0.5 # needs two effects
  )
  expect_error(
    model_default_odin(
      uk_population,
      intervention = list(contacts = intervention)
    )
  )
  expect_error(
    model_default_odin(
      uk_population,
      intervention = list(transmission_rate = intervention)
    ),
    regexp = "Must inherit from class 'rate_intervention'"
  )

  # expect error on poorly specified vaccination (needs 1 dose)
  vax_double_dose <- vaccination(
    nu = matrix(1e-3, nrow = 2, ncol = 2),
    time_begin = matrix(00, nrow = 2, ncol = 2),
    time_end = matrix(100, nrow = 2, ncol = 2)
  )
  expect_error(
    model_default_odin(
      uk_population,
      vaccination = vax_double_dose
    ),
    regexp = "Must have exactly 1 cols"
  )

  # expect error on poorly specified time-dependence function list
  expect_error(
    model_default_odin(
      uk_population,
      time_dependence = function(x) x
    ),
    regexp = "Must be of type 'list'"
  )
  expect_error(
    model_default_odin(
      uk_population,
      time_dependence = list(function(x) x)
    ),
    regexp = "Must have names"
  )
  expect_error(
    model_default_odin(
      uk_population,
      time_dependence = list(transmission_rate = function(x) x)
    ),
    regexp = "Must have first formal arguments \\(ordered\\): time,x."
  )
  expect_error(
    model_default_odin(
      uk_population,
      time_dependence = list(transmission_rate = NULL)
    ),
    regexp = "Contains missing values"
  )
})

# prepare vectors of parameters
beta <- rnorm(10, 1.3 / 7, sd = 0.01)
sigma <- rnorm(10, 0.5, sd = 0.01)
gamma <- rnorm(10, 1 / 7, sd = 0.01)

test_that("Default odin model: infection parameters as vectors", {
  # expect no conditions when vectors are passed
  expect_no_condition(
    model_default_odin(
      uk_population,
      transmission_rate = beta, infectiousness_rate = sigma,
      recovery_rate = gamma
    )
  )
  # expect output structure is a nested data.table
  output <- model_default_odin(
    uk_population,
    transmission_rate = beta, infectiousness_rate = sigma,
    recovery_rate = gamma
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(beta))
  expect_identical(output$transmission_rate, beta)
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(output$param_set, seq_along(beta))
  expect_identical(unique(output$scenario), 1L)

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null")
  )
  checkmate::expect_list(
    output$vaccination,
    types = c("vaccination", "null")
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null"), any.missing = FALSE
  )
})

test_that("Default model: composable elements as lists", {
  # expect no conditions when multiple interventions or vaccinations are passed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      contacts = intervention(
        "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
      )
    ),
    scenario_02 = list(
      contacts = intervention(
        "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
      ),
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      )
    )
  )

  expect_no_condition(
    model_default_odin(uk_population, intervention = npi_list)
  )

  # expect output is a nested data.frame-like object
  output <- model_default_odin(uk_population, intervention = npi_list)
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(npi_list))
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(output$scenario, seq_along(npi_list))
  expect_identical(unique(output$param_set), 1L)

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  # some interventions may be missing
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null")
  )
  checkmate::expect_list(
    output$vaccination,
    types = c("vaccination", "null")
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null")
  )
})

test_that("Default model: multi-parameter, multi-composables", {
  # expect no conditions when multiple interventions or vaccinations are passed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      contacts = intervention(
        "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
      )
    ),
    scenario_02 = list(
      contacts = intervention(
        "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
      ),
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      )
    )
  )
  # reuse parameter sets from earlier tests
  expect_no_condition(
    model_default_odin(
      uk_population,
      transmission_rate = beta, recovery_rate = gamma,
      intervention = npi_list
    )
  )

  # expect output is a nested data.frame-like object
  output <- model_default_odin(
    uk_population,
    transmission_rate = beta, recovery_rate = gamma,
    intervention = npi_list
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(npi_list) * length(beta))
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(
    output$scenario, rep(seq_along(npi_list), length(beta))
  )
  expect_identical(unique(output$param_set), seq_along(beta))
  expect_identical(
    output$param_set, rep(seq_along(beta), each = length(npi_list))
  )

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  # some interventions or vaccinations may be missing
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null"), any.missing = TRUE
  )
  checkmate::expect_list(
    output$vaccination,
    types = c("vaccination", "null"), any.missing = TRUE
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null"), any.missing = FALSE
  )
})

test_that("Default odin model: errors on vectorised input", {
  # expect errors on poorly specified vector inputs
  expect_error(
    model_default_odin(
      uk_population,
      transmission_rate = beta[-1], recovery_rate = gamma
    ),
    regexp = "All parameters must be of the same length, or must have length 1"
  )
  expect_error(
    model_default_odin(
      uk_population,
      intervention = list(
        NULL,
        list(dummy = intervention)
      )
    ),
    regexp =
      "`intervention` must be a list of <intervention>s or a list of such lists"
  )
  expect_error(
    model_default_odin(
      uk_population,
      vaccination = list(
        NULL,
        list(vaccination) # this is a list with the function vaccination()
      )
    ),
    regexp = "`vaccination` must be a <vaccination> or a list of <vaccination>s"
  )

  # expect time-dependence cannot be vectorised
  expect_error(
    model_default_odin(
      uk_population,
      time_dependence = list(
        time_dep_01 = list(
          transmission_rate = function(x) x
        ),
        time_dep_02 = list(
          transmission_rate = function(x) x
        )
      )
    ),
    regexp = "May only contain the following types: \\{function\\}"
  )
})
