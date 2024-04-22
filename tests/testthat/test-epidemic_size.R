# Test function for epidemic size
# create initial conditions
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 2, ncol = 5L,
  byrow = TRUE
)
colnames(initial_conditions) <- c(
  "susceptible", "exposed", "infectious", "recovered", "vaccinated"
)

# create a population where there are no contacts between groups
# this helps test the expectation that the final size proportions
# should be the same as the demographic group proportions
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(c(1, 0, 0, 1), 2, 2), # within group contacts only
  demography_vector = 67e6 * c(0.4, 0.6),
  initial_conditions = initial_conditions
)

time_end <- 200
# run epidemic simulation with no vaccination or intervention
data <- model_default(
  population = uk_population,
  time_end = time_end,
  increment = 1
)

test_that("Epidemic size functions", {
  # test for the initial size = 0.0
  epidemic_initial_size <- epidemic_size(data, stage = 0.0)
  expect_equal(
    epidemic_initial_size,
    initial_conditions[, "infectious"],
    ignore_attr = TRUE
  )
  expect_equal(
    epidemic_size(data, time = 0),
    epidemic_initial_size,
    ignore_attr = TRUE
  )

  # test the final size
  epidemic_final_size <- epidemic_size(data)
  expect_equal(
    epidemic_final_size,
    data[data$compartment == "recovered" & data$time == max(data$time), ]$value,
    ignore_attr = TRUE
  )
  expect_equal(
    epidemic_size(data, time = time_end),
    epidemic_final_size,
    ignore_attr = TRUE
  )
  expect_equal(
    epidemic_size(data, time = time_end),
    epidemic_size(data)
  )

  # expect return types and contents
  expect_s3_class(
    epidemic_size(data, simplify = FALSE),
    "data.table"
  )
  expect_s3_class(
    epidemic_size(data, by_group = FALSE, simplify = FALSE),
    "data.table"
  )
  expect_s3_class(
    epidemic_size(data, time = c(1, 2), simplify = FALSE),
    "data.table"
  )

  # expect that the final size proportion is the same as the demography prop.
  expect_equal(
    epidemic_final_size / sum(epidemic_final_size),
    uk_population$demography_vector / sum(uk_population$demography_vector),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  # test that final size is greater than size at 50% epidemic time
  epidemic_half_size <- epidemic_size(data, stage = 0.5)
  expect_true(
    all(epidemic_final_size > epidemic_half_size)
  )

  # test that by group FALSE returns a single value
  expect_length(
    epidemic_size(data, by_group = FALSE), 1
  )
})

#### Test that epidemic size with no deaths is same as include_deaths = FALSE
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

test_that("Epidemic size with no deaths is correct", {
  data <- model_vacamole(
    population = uk_population,
    mortality_rate = 0,
    time_end = 400, increment = 1
  )
  expect_identical(
    epidemic_size(data, include_deaths = FALSE),
    epidemic_size(data)
  )
})

test_that("Epidemic size for ebola model with replicates", {
  # prepare data
  demography_vector <- 67000
  replicates <- 100L

  pop <- population(
    contact_matrix = matrix(1),
    demography_vector = demography_vector,
    initial_conditions = matrix(
      c(1 - 1e-3, 1e-3 / 2, 1e-3 / 2, 0, 0, 0),
      nrow = 1
    )
  )

  # expect that function returns data.table when replicates > 1
  output <- model_ebola(pop, replicates = replicates)
  # throws a warning as expected
  suppressWarnings(
    data <- epidemic_size(output, simplify = TRUE)
  )
  expect_s3_class(data, "data.table")
  expect_identical(
    nrow(data), replicates
  )

  # expect that output can be simplified when replicates = 1
  output <- model_ebola(pop, replicates = 1L)
  data <- epidemic_size(output, simplify = TRUE)
  expect_vector(data, numeric())
})

test_that("Epidemic size warnings and messages", {
  # expect message when time is specified
  expect_message(epidemic_size(data, time = 50))

  # expect warning when multiple time points are requested and simplify is TRUE
  expect_warning(epidemic_size(data, stage = c(0.1, 1.0)))
  expect_warning(epidemic_size(data, time = c(10, 50)))
  expect_error(
    epidemic_size(data, time = NULL, stage = NULL),
    regexp = "One of `stage` or `time` must be provided; both are NULL!"
  )

  # expect message when dead are requested but not present
  expect_warning(epidemic_size(data, time = 50, include_deaths = TRUE))

  # expect input checking errors
  expect_error(epidemic_size("data"))
  expect_error(epidemic_size(data[, c("time")]))
  data_test <- data
  colnames(data_test) <- letters[1:4]
  expect_error(
    epidemic_size(data_test),
    regexp = "Names must include the elements"
  )

  expect_error(epidemic_size(data, by_group = 1))
  expect_error(epidemic_size(data, include_deaths = 1))
  expect_error(epidemic_size(data, stage = "0.1"))
  expect_error(epidemic_size(data, time = "0.1"))
  expect_error(epidemic_size(data, time = 0.1))
})
