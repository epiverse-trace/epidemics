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
r0 <- 1.5
preinfectious_period <- 3
infectious_period <- 7

test_that("Output of default epidemic model", {
  # run epidemic model
  data <- epidemic(
    model_name = "default",
    population = uk_population,
    r0 = r0,
    intervention = no_intervention(uk_population),
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    time_end = 100, increment = 1.0
  )

  # check for output type and contents
  expect_s3_class(data, "data.frame")
  expect_length(data, length(uk_population$initial_conditions) + 1L)

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      vapply(subset(data, select = -time),
        FUN = function(x) {
          all(x >= 0 & x <= sum(uk_population$demography_vector))
        },
        FUN.VALUE = TRUE
      )
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(subset(data[1, ], select = -time)),
    sum(subset(tail(data, 1L), select = -time)),
    tolerance = 1e-6
  )

  # check that all age groups in the simulation are the same
  # size as the demography vector
  final_state <- matrix(
    unlist(subset(tail(data, 1), select = -time)),
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
  r0_list <- list(
    r0_low = r0,
    r0_high = r0 + 1.0
  )

  # get data
  data <- lapply(
    r0_list,
    function(r0_) {
      # run model on data
      data <- epidemic(
        population = uk_population,
        r0 = r0_,
        preinfectious_period = preinfectious_period,
        infectious_period = infectious_period,
        time_end = 10, increment = 1.0
      )
    }
  )

  # get final size as total recoveries
  final_sizes <- lapply(data, tail, 1)
  final_sizes <- vapply(final_sizes, `[[`, FUN.VALUE = 1, "recovered_1")

  # test for effect of R0
  expect_gt(
    final_sizes["r0_high"],
    final_sizes["r0_low"]
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
r0 <- 1.5
preinfectious_period <- 3
infectious_period <- 7

test_that("Identical population sizes lead to identical final size", {
  data <- epidemic(
    population = dummy_population,
    r0 = r0,
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    time_end = 200, increment = 0.1
  )

  final_sizes <- unname(
    unlist(
      subset(
        tail(data, 1),
        select = grepl("recovered", colnames(data), fixed = TRUE)
      )
    )
  )

  # both groups have same final size
  expect_identical(
    final_sizes[1], final_sizes[2],
    tolerance = 1e-6
  )
})

test_that("Lower preinfectious period leads to larger final size", {
  # make a temporary pre-infectious period vector
  # lower values mean quicker transition from E => I
  preinfectious_period_low <- 1.2
  preinfectious_period_high <- 5.0

  data <- lapply(
    c(preinfectious_period_low, preinfectious_period_high),
    function(x) {
      epidemic(
        population = dummy_population,
        r0 = r0,
        preinfectious_period = x,
        infectious_period = infectious_period,
        time_end = 200, increment = 0.1
      )
    }
  )

  final_sizes <- lapply(
    data,
    function(df) {
      sum(
        unlist(
          subset(
            tail(df, 1),
            select = grepl("recovered", colnames(df), fixed = TRUE)
          )
        )
      )
    }
  )

  # both groups have same final size
  expect_gt(
    final_sizes[[1]], final_sizes[[2]]
  )
})

test_that("Lower infectious period leads to larger final size", {
  # make a temporary infectious period vector
  # lower values mean quicker transition from I => R
  infectious_period_low <- 5
  infectious_period_high <- 7

  data <- lapply(
    c(infectious_period_low, infectious_period_high),
    function(x) {
      epidemic(
        population = dummy_population,
        r0 = r0,
        preinfectious_period = preinfectious_period,
        infectious_period = x,
        time_end = 200, increment = 0.1
      )
    }
  )

  final_sizes <- lapply(
    data,
    function(df) {
      sum(
        unlist(
          subset(
            tail(df, 1),
            select = grepl("recovered", colnames(df), fixed = TRUE)
          )
        )
      )
    }
  )

  # group 2 must have a larger final size
  expect_gt(
    final_sizes[[1]], final_sizes[[2]]
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
    r0 = r0,
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    time_end = 200, increment = 0.1
  )

  final_sizes <- unlist(
    subset(
      tail(data, 1),
      select = grepl("recovered", colnames(data), fixed = TRUE)
    )
  )

  # group 1 with more contacts has higher final size
  expect_gt(
    final_sizes[1], final_sizes[2]
  )

  # calculate individuals still infected and check that
  # group with more contacts has more current infections
  current_infections <- unlist(
    subset(
      tail(data, 1),
      select = grepl("infect", colnames(data), fixed = TRUE)
    )
  )
  expect_gt(
    current_infections[1], current_infections[2]
  )
})
