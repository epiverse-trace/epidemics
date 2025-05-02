#### Tests for the stochastic SEIR model ####

base_seed <- .Random.seed
set.seed(1)

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
  I = 1e-6, R = 0
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
  "susceptible", "exposed", "infectious", "recovered"
)
n_samples <- 100L

test_that("Stochastic SEIR model: basic expectations, scalar arguments", {
  # expect run with no conditions for default arguments
  expect_no_condition(model_stochastic_seir(uk_population))

  # expect data.frame-nheriting output with 4 cols; C++ model time begins at 0
  data <- model_stochastic_seir(uk_population, n_samples = n_samples )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L)
  expect_named(
    data, c( "sample", "time", "demography_group", "compartment", "value" ),
    ignore.order = TRUE
  )
  expect_identical(
    nrow(data),
    length(demography_vector) * (time_end + 1L) * length(compartments) * n_samples
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

  # expect constant population size overall and per demography-group
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )
  ############  
  final_state <- data[ time == max( time ), .( count = sum( value ) ), 
                       by = c("demography_group", "sample" ) ]
  dt_demography <- data.table( pop = uk_population$demography_vector, 
                               demography_group = rownames(contact_matrix))
  final_state <- dt_demography[ final_state, on = "demography_group"]
  
  expect_identical(
    final_state[, count], final_state[ ,pop],
    tolerance = 1e-6
  )
})

# NOTE: statistical correctness is not expected to change for vectorised input
test_that("Stochastic SEIR model: statistical correctness, parameters", {
  # expect final size increases with transmission_rate
  size_beta_low <- epidemic_size(
    model_stochastic_seir(uk_population, transmission_rate = 1.3 / 7.0, n_samples = n_samples )
  )
  size_beta_high <- epidemic_size(
    model_stochastic_seir(uk_population, transmission_rate = 1.5 / 7.0, n_samples = n_samples )
  )
  expect_true(
    all(size_beta_high > size_beta_low)
  )

  # expect final size increases with infectiousness rate (lower incubation time)
  size_sigma_low <- epidemic_size(
    model_stochastic_seir(uk_population, infectiousness_rate = 1 / 5, n_samples = n_samples )
  )
  size_sigma_high <- epidemic_size(
    model_stochastic_seir(uk_population, infectiousness_rate = 1 / 2, , n_samples = n_samples )
  )
  expect_true(
    all(size_sigma_high > size_sigma_low)
  )

  # expect final size increases with initial infections
  initial_conditions_high <- c(
    S = 1 - 10e-6, E = 0, I = 10e-6,
    R = 0
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
  size_infections_low <- epidemic_size(
    model_stochastic_seir(uk_population, n_samples = n_samples ) 
  )
  size_infections_high <- epidemic_size(
    model_stochastic_seir(uk_population_high_infections, n_samples = n_samples )
  )
  expect_true(
    all(size_infections_high > size_infections_low)
  )
})

# prepare baseline for comparison of against intervention scenarios
data_baseline <- model_stochastic_seir(uk_population, n_samples = n_samples )

test_that("Stochastic SEIR model: contacts interventions and stats. correctness", {
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
  )
  # repeat some basic checks from default case with no intervention
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_stochastic_seir(
      uk_population,
      intervention = list(contacts = intervention),
      n_samples = n_samples
    )
  )

  # expect data.frame-inheriting output with 5 cols
  data <- model_stochastic_seir(
    uk_population,
    intervention = list(contacts = intervention),
    n_samples = n_samples
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L)

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
    model_stochastic_seir(
      uk_population,
      intervention = list(contacts = combined_interventions),
      n_samples = n_samples
    )
  )
  data_combined <- model_stochastic_seir(
    uk_population,
    intervention = list(contacts = combined_interventions),
    n_samples = n_samples
  )
  # expect epidemic size is lower for combined intervention
  expect_true(
    all(epidemic_size(data_combined) < epidemic_size(data))
  )
})

test_that("Stochastic SEIR model: rate interventions", {
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
    model_stochastic_seir(
      uk_population,
      intervention = list(transmission_rate = intervention),
      n_samples = n_samples
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_stochastic_seir(
    uk_population,
    intervention = list(transmission_rate = intervention),
    n_samples = n_samples
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L)

  # expect final size is lower with intervention
  expect_true(
    all(epidemic_size(data_baseline) > epidemic_size(data))
  )
})

test_that("Stochastic SEIR model: errors and warnings, scalar arguments", {
  # expect errors on basic input checking
  expect_error(
    model_stochastic_seir(population = "population",n_samples = n_samples),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  expect_error(
    model_stochastic_seir(population = population,n_samples = n_samples),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  pop_wrong_compartments <- uk_population
  pop_wrong_compartments$initial_conditions <- initial_conditions[, -1]
  expect_error(
    model_stochastic_seir(pop_wrong_compartments,n_samples = n_samples),
    regexp = "(Assertion on)*(initial_conditions)*failed"
  )

  # expect errors for infection parameters
  expect_error(
    model_stochastic_seir(uk_population, transmission_rate = "0.19",n_samples = n_samples),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_stochastic_seir(uk_population, infectiousness_rate = list(0.2),n_samples = n_samples),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_stochastic_seir(uk_population, recovery_rate = "0.19",n_samples = n_samples),
    regexp = "Must be of type 'numeric'"
  )

  # expect error on time parameters
  expect_error(
    model_stochastic_seir(uk_population, time_end = "100",n_samples = n_samples),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_stochastic_seir(uk_population, time_end = 100.5,n_samples = n_samples),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_stochastic_seir(uk_population, time_end = c(100, -100, 10),n_samples = n_samples),
    regexp = "(Element)*(is not >= 0)"
  )
  # expect error on poorly specified interventions
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, 0.5 # needs two effects
  )
  expect_error(
    model_stochastic_seir(
      uk_population,
      intervention = list(contacts = intervention),
      n_samples = n_samples
    )
  )
  expect_error(
    model_stochastic_seir(
      uk_population,
      intervention = list(transmission_rate = intervention),
      n_samples = n_samples
    )
  )
})

# reset seed not to disturb other tests
.Random.seed <- base_seed

