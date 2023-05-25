
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
  contact_matrix = matrix(1, 2, 2),
  demography_vector = 67e6 * c(0.4, 0.6),
  initial_conditions = initial_conditions
)

# prepare a two dose vaccination regime for three age groups
double_vaccination <- vaccination(
  name = "double_vaccination",
  nu = matrix(
    1e-3,
    nrow = 3, ncol = 2
  ),
  time_begin = matrix(
    00,
    nrow = 3, ncol = 2
  ),
  time_end = matrix(
    200,
    nrow = 3, ncol = 2
  )
)

# make infection class for Vacamole model
# note extra arguments passed as ...
infect <- infection(
  name = "covid", r0 = 5, infectious_period = 10,
  preinfectious_period = 5,
  eta = 1 / 1000, omega = 1 / 1000,
  susc_reduction_vax = 0.5,
  hosp_reduction_vax = 0.7,
  mort_reduction_vax = 0.9
)

test_that("Vacamole model works", {
  # check model runs silently
  expect_silent(
    epidemic(
      model_name = "vacamole",
      population = uk_population,
      infection = infect,
      vaccination = double_vaccination,
      time_end = 400, increment = 1
    )
  )

  data <- epidemic(
    model_name = "vacamole",
    population = uk_population,
    infection = infect,
    vaccination = double_vaccination,
    time_end = 400, increment = 1
  )

  # expect output is a data.table
  expect_s3_class(data, "data.table")
  # expect output has correct compartments
  expect_identical(
    unique(data$compartment),
    read_from_library(model_name = "vacamole", what = "compartments")
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

#### Tests for features of Vacamole ####
# prepare a null vaccination schedule
no_vaccination <- no_vaccination(uk_population, doses = 2)

# make infection class for Vacamole model
# note extra arguments passed as ...
infect <- infection(
  name = "covid", r0 = 5, infectious_period = 10,
  preinfectious_period = 5,
  eta = 1 / 1000, omega = 1 / 1000,
  susc_reduction_vax = 0.5,
  hosp_reduction_vax = 0.7,
  mort_reduction_vax = 0.9
)

test_that("Vacamole model with no vaccination", {
  # check model runs silently
  data <- epidemic(
    model_name = "vacamole",
    population = uk_population,
    infection = infect,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are vaccinated in any compartments related to
  # vaccination - e.g. hospitalised-vaccinated etc.
  pop_vaxxed <- data[time == max(time) &
    grepl("vaccinated", data$compartment, fixed = TRUE)]$value
  expect_identical(
    unique(pop_vaxxed), 0.0,
    tolerance = 1e-6
  )
})

nonlethal_infect <- infection(
  name = "covid", r0 = 5, infectious_period = 10,
  preinfectious_period = 5,
  eta = 1 / 1000,
  omega = 0, # no deaths due to infection
  susc_reduction_vax = 0.5,
  hosp_reduction_vax = 0.7,
  mort_reduction_vax = 0.9
)

test_that("Vacamole with non-fatal infection", {
  data <- epidemic(
    model_name = "vacamole",
    population = uk_population,
    infection = nonlethal_infect,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are dead
  pop_dead <- data[time == max(time) &
    grepl("dead", data$compartment, fixed = TRUE)]$value
  expect_identical(
    unique(pop_dead), 0.0,
    tolerance = 1e-6
  )
})

no_hospitalisation <- infection(
  name = "covid", r0 = 5, infectious_period = 10,
  preinfectious_period = 5,
  eta = 0, # no hospitalisation
  omega = 1 / 1000,
  susc_reduction_vax = 0.5,
  hosp_reduction_vax = 0.7,
  mort_reduction_vax = 0.9
)

test_that("Vacamole with no hospitalisation", {
  data <- epidemic(
    model_name = "vacamole",
    population = uk_population,
    infection = no_hospitalisation,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are dead
  pop_hospitalised <- data[time == max(time) &
    grepl("hospitalised", data$compartment, fixed = TRUE)]$value
  expect_identical(
    unique(pop_hospitalised), 0.0,
    tolerance = 1e-6
  )
})
