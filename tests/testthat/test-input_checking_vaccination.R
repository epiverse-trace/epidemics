# Tests for input checking on vaccination objects
# prepare a population and a vaccination
contact_matrix <- matrix(1, 2, 2)
demography_vector <- c(1e6, 1e6)

# Prepare a test population
test_population <- population(
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  # initial conditions for the default model
  initial_conditions = matrix(
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Prepare a vaccination
test_vaccination <- vaccination(
  time_begin = matrix(0, nrow = nrow(contact_matrix)),
  time_end = matrix(100, nrow = nrow(contact_matrix)),
  nu = matrix(1e-4, nrow = nrow(contact_matrix))
)

# prepare a bad vaccination that does not match the population
test_vaccination_bad <- vaccination(
  time_begin = matrix(0),
  time_end = matrix(100),
  nu = matrix(1e-4)
)

test_that("Vaccinations are checked correctly", {
  # check for no conditions on a well formed vaccination
  expect_no_condition(
    assert_vaccination(
      test_vaccination,
      doses = 1L, population = test_population
    )
  )
  expect_no_condition(
    assert_vaccination(
      test_vaccination,
      doses = 1L # with population missing
    )
  )
  expect_error(
    assert_vaccination(
      test_vaccination,
      doses = 1L, population = "population" # population is wrong class
    )
  )

  # expect error when the number of doses is incorrect
  expect_error(
    assert_vaccination(
      test_vaccination,
      doses = 2L, population = test_population
    )
  )
  # expect error when the number of demographic groups is incorrect
  expect_error(
    assert_vaccination(
      test_vaccination_bad,
      doses = 2L, population = test_population
    )
  )
})
