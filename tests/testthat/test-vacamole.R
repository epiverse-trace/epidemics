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
    nrow = 2, ncol = 2
  ),
  time_begin = matrix(
    00,
    nrow = 2, ncol = 2
  ),
  time_end = matrix(
    200,
    nrow = 2, ncol = 2
  )
)

test_that("Vacamole model works", {
  # check model runs silently
  expect_no_condition(
    model_vacamole_cpp(
      population = uk_population,
      vaccination = double_vaccination,
      time_end = 400, increment = 1
    )
  )

  data <- model_vacamole_cpp(
    population = uk_population,
    vaccination = double_vaccination,
    time_end = 400, increment = 1
  )

  # expect output is a data.table
  expect_s3_class(data, "data.frame")
  # expect output has correct compartments
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

#### Tests for features of Vacamole ####
# prepare a null vaccination schedule
no_vaccination <- no_vaccination(uk_population, doses = 2)

test_that("Vacamole model with no vaccination", {
  # check model runs silently
  data <- model_vacamole_cpp(
    population = uk_population,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are vaccinated in any compartments related to
  # vaccination - e.g. hospitalised-vaccinated etc.
  pop_vaxxed <- data[data$time == max(data$time) &
    grepl("vaccinated", data$compartment, fixed = TRUE), ]$value
  expect_identical(
    unique(pop_vaxxed), 0.0,
    tolerance = 1e-6
  )
})

test_that("Vacamole with non-fatal infection", {
  data <- model_vacamole_cpp(
    population = uk_population,
    mortality_rate = 0,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are dead
  pop_dead <- data[data$time == max(data$time) &
    grepl("dead", data$compartment, fixed = TRUE), ]$value
  expect_identical(
    unique(pop_dead), 0.0,
    tolerance = 1e-6
  )
})

test_that("Vacamole with no hospitalisation", {
  data <- model_vacamole_cpp(
    population = uk_population,
    hospitalisation_rate = 0,
    vaccination = no_vaccination,
    time_end = 400, increment = 1
  )
  # test that no individuals are dead
  pop_hospitalised <- data[data$time == max(data$time) &
    grepl("hospitalised", data$compartment, fixed = TRUE), ]$value
  expect_identical(
    unique(pop_hospitalised), 0.0,
    tolerance = 1e-6
  )
})

#### Test for Vacamole model run with wrong inputs ####
test_that("Vacamole model errors correctly", {
  # error because vaccination doses are wrong
  expect_error(
    model_vacamole_cpp(
      population = uk_population,
      vaccination = no_vaccination(uk_population, doses = 3L)
    )
  )

  # error because disallowed model parameters are targeted
  expect_error(
    model_vacamole_cpp(
      population = uk_population,
      vaccination = no_vaccination(uk_population, doses = 2L),
      intervention = list(
        transmissibility = no_rate_intervention(), # works okay
        transmissibility_vax = no_rate_intervention() # should error
      )
    )
  )

  expect_error(
    model_vacamole_cpp(
      population = uk_population,
      vaccination = no_vaccination(uk_population, doses = 2L),
      time_dependence = list(
        transmissibility = no_time_dependence(), # works okay
        transmissibility_vax = no_time_dependence() # should error
      )
    )
  )
})

#### Tests for the R implementation of the Vacamole model ####
# basic expectations
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
