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

# prepare a basic intervention
close_schools <- intervention(
  name = "close_schools",
  type = "contacts",
  time_begin = 100, time_end = 150,
  reduction = matrix(c(0.2, 0.0))
)

# snapshot test for printing
test_that("Printing <contacts_intervention> class", {
  expect_snapshot(close_schools)
})

# test the intervention has expected structure
test_that("Contacts intervention is correctly initialised", {
  expect_s3_class(close_schools, "contacts_intervention")
  expect_named(
    close_schools, c("name", "time_begin", "time_end", "reduction")
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
    close_schools$reduction, "double"
  )
  expect_length(
    close_schools$reduction,
    nrow(contact_matrix)
  )
})

test_that("Contacts intervention reduces final size", {
  # run model with intervention
  data_intervention <- model_default(
    population = uk_population,
    intervention = list(
      contacts = close_schools
    ),
    time_end = 200, increment = 1.0
  )

  # run model without intervention
  data <- model_default(
    population = uk_population,
    time_end = 200, increment = 1.0
  )

  # expect that intervention reduces epidemic final size
  # test only for a single group, to which intervention is applied
  final_size_intervention <- epidemic_size(data_intervention, by_group = FALSE)
  final_size <- epidemic_size(data, by_group = FALSE)

  expect_lt(
    final_size_intervention,
    final_size
  )
})

# Test for errors on wrong intervention formulation
# prepare a basic intervention
badly_formed_intervention <- intervention(
  name = "close_schools",
  type = "contacts",
  time_begin = 100, time_end = 150,
  reduction = rep(0.2, 3) # too many values, 2 needed, 3 given
)

test_that("Error on poorly specified contacts intervention", {
  # expect failure for poorly specified intervention
  expect_error(
    model_default(
      population = uk_population,
      intervention = badly_formed_intervention,
      time_end = 200, increment = 1.0
    )
  )
})

test_that("Null contacts intervention is correctly initialised", {
  null_intervention <- .no_contacts_intervention(uk_population)
  # expect no message using helper function no_contacts_intervention()
  expect_no_condition(
    .no_contacts_intervention(uk_population)
  )
  expect_identical(
    null_intervention$time_begin, null_intervention$time_end
  )

  expect_type(
    null_intervention$reduction, "double"
  )
  expect_length(
    null_intervention$reduction,
    nrow(uk_population$contact_matrix)
  )
})

# Tests for multiple stacked contact interventions
npi_1 <- intervention(
  type = "contacts",
  time_begin = 30,
  time_end = 60,
  reduction = rep(0.15, 3)
)

# second dose regime
npi_2 <- intervention(
  type = "contacts",
  time_begin = 45,
  time_end = 75,
  reduction = rep(0.1, 3)
)

# Tests for basic expectations of c.intervention
test_that("Concatenating `intervention`s works", {
  expect_no_condition(
    c(npi_1, npi_2)
  )
  multi_npi <- c(npi_1, npi_2)
  # expect that rows sum to expected values
  expect_identical(
    rowSums(multi_npi$reduction),
    rep(0.25, 3L)
  )
  # expect that there are only two columns
  expect_identical(
    ncol(multi_npi$reduction),
    2L
  )
  # snapshot of the multi-NPI
  expect_snapshot(
    multi_npi
  )
})
