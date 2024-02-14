#### Tests for the Vacamole model ####
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

# // 0| 1| 2|3| 4|5| 6|7| 8|9|10
# // S|V1|V2|E|EV|I|IV|H|HV|D|R

# make initial conditions - order is important
initial_conditions <- c(
  S = 1 - 1e-6,
  V1 = 0, V2 = 0,
  E = 0, EV = 0,
  I = 1e-6, IV = 0,
  H = 0, HV = 0, D = 0, R = 0
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
double_vaccination <- vaccination(
  name = "double_vaccination",
  nu = matrix(
    1e-3,
    nrow = 2, ncol = 2
  ),
  time_begin = matrix(
    00,
    nrow = 2, ncol = 2
  ),
  time_end = matrix(
    100,
    nrow = 2, ncol = 2
  )
)

# model run time
time_end <- 100L
compartments <- c(
  "susceptible", "vaccinated_one_dose", "vaccinated_two_dose", "exposed",
  "exposed_vaccinated", "infectious", "infectious_vaccinated", "hospitalised",
  "hospitalised_vaccinated", "dead", "recovered"
)

test_that("Vacamole model: basic expectations, scalar arguments", {
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_vacamole_cpp(uk_population)
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  # NOTE: increased hospitalisation and mortality for testing purposes
  data <- model_vacamole_cpp(
    uk_population,
    time_end = time_end,
    hospitalisation_rate = 1 / 100, mortality_rate = 1 / 100
  )
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
    unique(data[grepl("dose", data$compartment, fixed = TRUE), ]$value), 0
  )

  # expect individuals in hospitalised compartments
  # NOTE: vaccinated individuals are hospitalised at very low rates, not checked
  checkmate::expect_numeric(
    tail(data[data$compartment == "hospitalised", ]$value),
    lower = 1
  )

  # expect individuals in dead compartment
  checkmate::expect_numeric(
    tail(data[grepl("dead", data$compartment, fixed = TRUE), ]$value),
    lower = 1
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

test_that("Vacamole model: two dose vaccination, scalar arguments", {
  # repeat some basic checks from default case with no vaccination
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_vacamole_cpp(uk_population, vaccination = double_vaccination)
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_vacamole_cpp(uk_population, vaccination = double_vaccination)
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)

  # expect non-zero vaccinations towards simulation end
  checkmate::expect_numeric(
    tail(data[grepl("dose", data$compartment, fixed = TRUE), ]$value),
    lower = 10
  )
})

test_that("Vacamole model: contacts interventions", {
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, c(0.5, 0.0)
  )
  # repeat some basic checks from default case with no intervention
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_vacamole_cpp(
      uk_population,
      intervention = list(contacts = intervention)
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_vacamole_cpp(
    uk_population,
    intervention = list(contacts = intervention)
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)
})

test_that("Vacamole model: rate interventions", {
  # TODO: ADD TESTS
})

test_that("Vacamole model: time dependence", {
  # TODO: ADD TESTS
})

# NOTE: statistical correctness is not expected to change for vectorised input
test_that("Vacamole model: statistical correctness", {
  # expect final size increases with transmissibility
  size_beta_low <- epidemic_size(
    model_vacamole_cpp(uk_population, transmissibility = 1.3 / 7.0)
  )
  size_beta_high <- epidemic_size(
    model_vacamole_cpp(uk_population, transmissibility = 1.5 / 7.0)
  )
  expect_true(
    all(size_beta_high > size_beta_low)
  )

  # expect final size increases with infectiousness rate (lower incubation time)
  size_sigma_low <- epidemic_size(
    model_vacamole_cpp(uk_population, infectiousness_rate = 1 / 5)
  )
  size_sigma_high <- epidemic_size(
    model_vacamole_cpp(uk_population, infectiousness_rate = 1 / 2)
  )
  expect_true(
    all(size_sigma_high > size_sigma_low)
  )

  # expect final size increases with initial infections
  initial_conditions_high <- c(
    S = 1 - 10e-6, V1 = 0, V2 = 0, E = 0, EV = 0, I = 10e-6, IV = 0,
    H = 0, HV = 0, D = 0, R = 0
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
  size_infections_low <- epidemic_size(model_vacamole_cpp(uk_population))
  size_infections_high <- epidemic_size(
    model_vacamole_cpp(uk_population_high_infections)
  )
  expect_true(
    all(size_infections_high > size_infections_low)
  )

  # expect no deaths when mortality is 0
  data_no_mortality <- model_vacamole_cpp(uk_population, mortality_rate = 0)
  expect_true(
    all(data_no_mortality[data_no_mortality$compartment == "dead", ]$value == 0)
  )

  # expect final size decreases with mortality rate (reduces infectious indivs.)
  # NOTE: this may be a non-linear relationship, not tested here
  size_omega_low <- epidemic_size(
    model_vacamole_cpp(uk_population, mortality_rate = 0.001),
    include_deaths =
    )
  size_omega_high <- epidemic_size(
    model_vacamole_cpp(uk_population, mortality_rate = 0.01)
  )
  expect_true(
    all(size_omega_low > size_omega_high)
  )

  # expect no hospitalisations when hospitalisation rate = 0
  data <- model_vacamole_cpp(uk_population, hospitalisation_rate = 0.0)
  expect_identical(
    unique(data[grepl("hospitalised", data$compartment, fixed = TRUE)]$value),
    0
  )
})

test_that("Vacamole model: errors and warnings, scalar arguments", {
  # expect errors on basic input checking
  expect_error(
    model_vacamole_cpp(population = "population"),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  expect_error(
    model_vacamole_cpp(population = population),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  pop_wrong_compartments <- uk_population
  pop_wrong_compartments$initial_conditions <- initial_conditions[, -1]
  expect_error(
    model_vacamole_cpp(pop_wrong_compartments),
    regexp = "(Assertion on)*(initial_conditions)*failed"
  )

  # expect errors for infection parameters
  expect_error(
    model_vacamole_cpp(uk_population, transmissibility = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, infectiousness_rate = list(0.2)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, recovery_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, hospitalisation_rate = list(0.002)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, transmissibility = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, mortality_rate = list(0.002)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, mortality_rate_vax = "0.002"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, hospitalisation_rate_vax = "0.002"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, transmissibility = "0.002"),
    regexp = "Must be of type 'numeric'"
  )

  # expect error on time parameters
  expect_error(
    model_vacamole_cpp(uk_population, time_end = "100"),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, time_end = 100.5),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, time_end = c(100, -100, 10)),
    regexp = "(Element)*(is not >= 0)"
  )
  expect_error(
    model_vacamole_cpp(uk_population, increment = "0.1"),
    regexp = "Must be of type 'number'"
  )
  expect_error(
    model_vacamole_cpp(uk_population, increment = c(0.1, 0.2)),
    regexp = "Must have length 1"
  )

  # expect error on poorly specified interventions
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, 0.5 # needs two effects
  )
  expect_error(
    model_vacamole_cpp(
      uk_population,
      intervention = list(contacts = intervention)
    )
  )
  expect_error(
    model_vacamole_cpp(
      uk_population,
      intervention = list(transmissibility = intervention)
    ),
    regexp = "Must inherit from class 'rate_intervention'"
  )

  # TODO: ADD TESTS ON TIME-DEPENDENCE
})

# prepare vectors of parameters
beta <- rnorm(10, 1.3 / 7, sd = 0.01)
sigma <- rnorm(10, 0.5, sd = 0.01)
gamma <- rnorm(10, 1 / 7, sd = 0.01)

test_that("Vacamole model: infection parameters as vectors", {
  # expect no conditions when vectors are passed
  expect_no_condition(
    model_vacamole_cpp(
      uk_population,
      transmissibility = beta, infectiousness_rate = sigma,
      recovery_rate = gamma
    )
  )
  # expect output structure is a nested data.table
  output <- model_vacamole_cpp(
    uk_population,
    transmissibility = beta, infectiousness_rate = sigma,
    recovery_rate = gamma
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(beta))
  expect_identical(output$transmissibility, beta)
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(output$param_set, seq_along(beta))
  expect_identical(unique(output$scenario), 1L)

  # TODO: ADD TESTS for list column of interventions and vaccination
})

test_that("Vacamole model: composable elements as lists")

test_that("Vacamole model: multi-parameter, multi-composables")

#### Tests for the R implementation of the Vacamole model ####
# basic expectations
skip("Vacamole: R-only models not updated")
test_that("Output of the Vacamole epidemic model R", {
  # run epidemic model, expect no conditions
  expect_no_condition(
    model_vacamole_r(
      population = uk_population,
      intervention = list(contacts = no_contacts_intervention(uk_population)),
      vaccination = double_vaccination,
      time_end = 100, increment = 1.0
    )
  )

  data <- model_vacamole_r(
    population = uk_population,
    intervention = list(contacts = no_contacts_intervention(uk_population)),
    vaccination = double_vaccination,
    time_end = 100, increment = 1.0
  )

  # check for output type and contents
  expect_s3_class(data, "data.frame")
  expect_length(data, 4L)
  expect_named(
    data, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  expect_identical(
    unique(data$compartment),
    c(
      "susceptible", "vaccinated_one_dose", "vaccinated_two_dose",
      "exposed", "exposed_vaccinated", "infectious", "infectious_vaccinated",
      "hospitalised", "hospitalised_vaccinated", "dead", "recovered"
    )
  )

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      data$value >= 0 & data$value <= sum(uk_population$demography_vector)
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )

  # check that all age groups in the simulation are the same
  # size as the demography vector
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(contact_matrix)
  )
  expect_identical(
    rowSums(final_state),
    uk_population$demography_vector,
    tolerance = 1e-6
  )
})

# equivalence expectations
skip("Vacamole: R-only implementation not updated")
test_that("Equivalence of vacamole model R and Cpp", {
  # create an intervention and vaccination
  multi_intervention <- c(
    intervention(
      type = "contacts",
      time_begin = 50, time_end = 100,
      reduction = matrix(
        0.2, nrow(contact_matrix), 1
      )
    ),
    intervention(
      type = "contacts",
      time_begin = 70, time_end = 90,
      reduction = matrix(
        0.3, nrow(contact_matrix), 1
      )
    )
  )

  # run epidemic model, expect no conditions
  data_r <- model_vacamole_r(
    population = uk_population,
    intervention = list(contacts = multi_intervention),
    vaccination = double_vaccination,
    time_end = 100, increment = 1.0
  )

  data_cpp <- model_vacamole_cpp(
    population = uk_population,
    intervention = list(contacts = multi_intervention),
    vaccination = double_vaccination,
    time_end = 100, increment = 1.0
  )

  expect_identical(
    tail(data_r),
    tail(data_cpp),
    tolerance = 1.0 # tolerance of 1, although actual difference is around 0.1
  )

  expect_identical(
    epidemic_size(data_r),
    epidemic_size(data_cpp),
    tolerance = 1.0
  )
})
