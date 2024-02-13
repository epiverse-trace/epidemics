# Basic tests to check for functionality
#### Tests checking the output of a single parameter set ####
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

test_that("Output of default epidemic model Cpp, scalar arguments", {
  # run epidemic model, expect no condition
  time_end <- 100L
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )

  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      intervention = list(
        contacts = no_contacts_intervention(uk_population)
      ),
      time_end = time_end, increment = 1.0
    )
  )

  output <- model_default_cpp(
    population = uk_population,
    intervention = list(
      contacts = no_contacts_intervention(uk_population)
    ),
    time_end = time_end, increment = 1.0
  )

  # check for output type and column names
  expect_s3_class(output, "data.frame")
  expect_named(
    output,
    c(
      "time", "compartment", "demography_group", "value"
    ),
    ignore.order = TRUE
  )

  # expect Nrow is time_end * compartment * demography
  # account for time = 0 as well
  expect_identical(
    nrow(output),
    (time_end + 1L) * length(demography_vector) * length(compartments)
  )

  # check data columns and types
  expect_s3_class(
    output, "data.frame"
  )
  expect_identical(
    unique(output$compartment),
    c("susceptible", "exposed", "infectious", "recovered", "vaccinated")
  )
  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      output$value >= 0 & output$value <=
        sum(uk_population$demography_vector)
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(output[output$time == min(output$time), ]$value),
    sum(output[output$time == max(output$time), ]$value),
    tolerance = 1e-6
  )

  # check that all age groups in the simulation are the same
  # size as the demography vector
  final_state <- matrix(
    unlist(output[output$time == max(output$time), ]$value),
    nrow = nrow(contact_matrix)
  )
  expect_identical(
    rowSums(final_state),
    uk_population$demography_vector,
    tolerance = 1e-6
  )
})

#### Tests for statistical correctness ####
# sense checks for variation in epidemiological parameters
test_that("Higher transmissibility gives larger final size, default model", {
  # prepare epidemic model runs with different R0 estimates
  r0_low <- 1.1
  r0_high <- 1.5
  infectious_period <- 7

  # get data
  data <- lapply(
    # transmissibility = r0 / infectious period
    c(r0_low, r0_high) / infectious_period,
    function(beta) {
      # run model on data
      data <- model_default_cpp(
        population = uk_population,
        transmissibility = beta,
        time_end = 10, increment = 1.0
      )
    }
  )

  # get final size as total recoveries
  final_sizes <- lapply(data, epidemic_size)

  # test for effect of R0
  expect_true(
    all(final_sizes[[2]] > final_sizes[[1]])
  )
})

# Tests that require uniform contact matrices
# create a dummy population with a uniform contact matrix
# Prepare some initial objects
dummy_contact_matrix <- t(matrix(
  1, 2, 2,
  byrow = TRUE
))
dummy_demography_vector <- 10e6 * c(0.5, 0.5)

dummy_population <- population(
  name = "dummy population",
  contact_matrix = dummy_contact_matrix,
  demography_vector = dummy_demography_vector,
  initial_conditions = matrix(
    c(1 - 1e-6, 1e-6 * 0.9, 1e-6 * 0.1, 0, 0),
    nrow = nrow(dummy_contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

test_that("Identical population sizes lead to identical final size", {
  data <- model_default_cpp(
    population = dummy_population,
    time_end = 200, increment = 0.1
  )

  final_sizes <- epidemic_size(data)

  # both groups have same final size
  expect_equal(
    final_sizes[1], final_sizes[2],
    tolerance = 1e-6, ignore_attr = TRUE
  )
})

test_that("Higher infectiousness rate leads to larger final size", {
  # make a temporary pre-infectious period vector
  # lower values mean quicker transition from E => I
  infectiousness_rates <- 1 / c(2, 3) # 1 / pre-infectious period in days
  data <- lapply(
    infectiousness_rates,
    function(sigma) {
      model_default_cpp(
        population = dummy_population,
        infectiousness_rate = sigma,
        time_end = 200, increment = 0.1
      )
    }
  )

  final_sizes <- lapply(data, epidemic_size)

  # both groups have same final size
  expect_true(
    all(final_sizes[[1]] > final_sizes[[2]])
  )
})

test_that("Group with more contacts has larger final size and infections", {
  # make a temporary contact matrix
  # group 1 has more contacts
  contact_matrix <- matrix(
    c(12, 2, 1, 3),
    nrow = nrow(dummy_contact_matrix),
    ncol = ncol(dummy_contact_matrix),
    byrow = TRUE
  )
  # add to dummy pop
  dummy_population$contact_matrix <- contact_matrix

  data <- model_default_cpp(
    population = dummy_population,
    time_end = 200, increment = 0.1
  )

  final_sizes <- epidemic_size(data)

  # group 1 with more contacts has higher final size
  expect_gt(
    final_sizes[1], final_sizes[2]
  )

  # calculate individuals still infected and check that
  # group with more contacts has more current infections
  current_infections <- data[data$compartment == "infectious" &
    data$time == max(data$time), ]$value
  expect_gt(
    current_infections[1], current_infections[2]
  )
})

#### Tests for the R implementation of the default model ####
# basic expectations
test_that("Output of default epidemic model R", {
  # run epidemic model, expect no conditions
  expect_no_condition(
    model_default_r(
      population = uk_population,
      time_end = 100, increment = 1.0
    )
  )

  data <- model_default_r(
    population = uk_population,
    intervention = list(
      contacts = no_contacts_intervention(uk_population)
    ),
    time_end = 100, increment = 1.0
  )

  # check for output type and contents
  expect_s3_class(data, "data.frame")
  expect_length(data, 4L)
  expect_named(
    data, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  compartments_default <- c(
    "susceptible", "exposed", "infectious", "recovered", "vaccinated"
  )
  expect_identical(
    unique(data$compartment),
    compartments_default
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
test_that("Equivalence of default model R and Cpp", {
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
  vax_regime <- vaccination(
    time_begin = matrix(10, nrow(contact_matrix), 1),
    time_end = matrix(100, nrow(contact_matrix), 1),
    nu = matrix(0.01, nrow(contact_matrix), 1)
  )

  # run epidemic model, expect no conditions
  data_r <- model_default_r(
    population = uk_population,
    intervention = list(
      contacts = multi_intervention
    ),
    vaccination = vax_regime,
    time_end = 100, increment = 1.0
  )

  data_cpp <- model_default_cpp(
    population = uk_population,
    intervention = list(
      contacts = multi_intervention
    ),
    vaccination = vax_regime,
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
