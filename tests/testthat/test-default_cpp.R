# Basic tests to check for functionality
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
    c(0.9999, 0.0001, 0, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Prepare epi parameters
pandemic <- infection(
  r0 = 3,
  preinfectious_period = 3,
  infectious_period = 7
)

test_that("Output of default epidemic model", {
  # run epidemic model
  data <- epidemic(
    model_name = "default",
    population = uk_population,
    infection = pandemic,
    intervention = no_intervention(uk_population),
    time_end = 100, increment = 1.0
  )

  # check for output type and contents
  expect_s3_class(data, "data.table")
  expect_length(data, 4L)
  expect_named(
    data, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
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

test_that("Larger R0 leads to larger final size in default epidemic model", {
  # prepare epidemic model runs with different R0 estimates
  r0 <- 1.5
  infection_list <- list(
    infection_r0_low = infection(
      r0 = r0,
      preinfectious_period = 3,
      infectious_period = 7
    ),
    infection_r0_high = infection(
      r0 = r0 + 1.0,
      preinfectious_period = 3,
      infectious_period = 7
    )
  )

  # get data
  data <- lapply(
    infection_list,
    function(infection_) {
      # run model on data
      data <- epidemic(
        population = uk_population,
        infection = infection_,
        time_end = 10, increment = 1.0
      )
    }
  )

  # get final size as total recoveries
  final_sizes <- lapply(data, epidemic_size)

  # test for effect of R0
  expect_true(
    all(final_sizes[["r0_high"]] > final_sizes[["r0_low"]])
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

# prepare epidemiological parameters
pandemic <- infection(
  r0 = 1.5,
  preinfectious_period = 3,
  infectious_period = 7
)

test_that("Identical population sizes lead to identical final size", {
  data <- epidemic(
    population = dummy_population,
    infection = pandemic,
    time_end = 200, increment = 0.1
  )

  final_sizes <- epidemic_size(data)

  # both groups have same final size
  expect_identical(
    final_sizes[1], final_sizes[2],
    tolerance = 1e-6
  )
})

test_that("Lower preinfectious period leads to larger final size", {
  # make a temporary pre-infectious period vector
  # lower values mean quicker transition from E => I
  infection_list <- list(
    infection(
      r0 = 1.5,
      preinfectious_period = 1.2,
      infectious_period = 7
    ),
    infection(
      r0 = 1.5,
      preinfectious_period = 5,
      infectious_period = 7
    )
  )

  data <- lapply(
    infection_list,
    function(infection_) {
      epidemic(
        population = dummy_population,
        infection = infection_,
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

test_that("Lower infectious period leads to larger final size", {
  # make a temporary infectious period vector
  # lower values mean quicker transition from I => R
  infection_list <- list(
    infection(
      r0 = 1.5,
      preinfectious_period = 3,
      infectious_period = 5
    ),
    infection(
      r0 = 1.5,
      preinfectious_period = 3,
      infectious_period = 7
    )
  )

  data <- lapply(
    infection_list,
    function(infection_) {
      epidemic(
        population = dummy_population,
        infection = infection_,
        time_end = 200, increment = 0.1
      )
    }
  )

  final_sizes <- lapply(data, epidemic_size)

  # group 2 must have a larger final size
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

  data <- epidemic(
    population = dummy_population,
    infection = pandemic,
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
