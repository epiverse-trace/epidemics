# Basic tests to check for functionality
# Prepare contact matrix and demography vector
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 40),
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
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = 2, ncol = 5L,
    byrow = TRUE
  )
)

# Prepare epi parameters
r0 <- rep(1.5, nrow(contact_matrix))
preinfectious_period <- rep(3, nrow(contact_matrix))
infectious_period <- rep(7, nrow(contact_matrix))

# prepare a basic intervention
close_schools <- intervention(
  name = "close_schools",
  time_begin = 100, time_end = 150,
  contact_reduction = c(0.2, 0.0)
)

# test the intervention has expected structure
test_that("Intervention is correctly initialised", {
  expect_s3_class(close_schools, "intervention")
  expect_named(
    close_schools, c("name", "time_begin", "time_end", "contact_reduction")
  )
  expect_type(
    close_schools$name, "character"
  )
  expect_length(
    close_schools$name, 1L
  )
  expect_type(
    close_schools$time_begin, "double"
  )
  expect_length(
    close_schools$time_begin, 1L
  )
  expect_type(
    close_schools$time_end, "double"
  )
  expect_length(
    close_schools$time_end, 1L
  )
  expect_type(
    close_schools$contact_reduction, "double"
  )
  expect_length(
    close_schools$contact_reduction,
    nrow(contact_matrix)
  )
})

test_that("Intervention reduces final size", {
  # run model with intervention
  data_intervention <- epidemic_cpp(
    population = uk_population,
    r0 = r0,
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    intervention = close_schools,
    time_end = 200, increment = 1.0
  )

  # run model without intervention
  data <- epidemic_cpp(
    population = uk_population,
    r0 = r0,
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    intervention = no_intervention(uk_population),
    time_end = 200, increment = 1.0
  )

  # expect that intervention reduces epidemic final size
  # test only for a single group, to which intervention is applied
  final_size_intervention <- subset(tail(data_intervention, 1L),
    select = recovered_1
  )
  final_size <- subset(tail(data, 1L), select = recovered_1)

  expect_lt(
    final_size_intervention,
    final_size
  )
})

# Test for errors on wrong intervention formulation
# prepare a basic intervention
badly_formed_intervention <- intervention(
  name = "close_schools",
  time_begin = 100, time_end = 150,
  contact_reduction = c(0.2, 0.0, 0.2) # too many values, 2 needed, 3 given
)

test_that("Error on poorly specified intervention", {
  # expect failure for poorly specified intervention
  expect_error(
    epidemic_cpp(
      population = uk_population,
      r0 = r0,
      preinfectious_period = preinfectious_period,
      infectious_period = infectious_period,
      intervention = badly_formed_intervention,
      time_end = 200, increment = 1.0
    ),
    regexp = "(Intervention must have)|(number of elements)|(contact matrix)"
  )
})

test_that("Null intervention is correctly initialised", {
  null_intervention <- no_intervention(uk_population)
  # expect message as time end and time start are identical
  expect_message(
    {
      no_intervention(uk_population)
    },
    regexp = "(time_end)*(not greater than)*(time_begin)"
  )
  expect_identical(
    null_intervention$time_begin, null_intervention$time_end
  )

  expect_type(
    null_intervention$contact_reduction, "double"
  )
  expect_length(
    null_intervention$contact_reduction,
    nrow(uk_population$contact_matrix)
  )
})
