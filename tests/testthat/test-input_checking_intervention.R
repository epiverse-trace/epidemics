# Tests for input checking on intervention objects
# prepare a population and a intervention
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

# Prepare a intervention
test_intervention <- intervention(
  time_begin = 60,
  time_end = 100,
  contact_reduction = matrix(0.2, nrow(contact_matrix))
)

# prepare a bad intervention that does not match the population
test_intervention_bad <- intervention(
  time_begin = 60,
  time_end = 100,
  contact_reduction = matrix(0.2)
)

test_that("Interventions are checked correctly", {
  # check for no conditions on a well formed intervention
  expect_no_condition(
    assert_intervention(
      test_intervention,
      population = test_population
    )
  )
  expect_no_condition(
    assert_intervention(
      test_intervention # with population missing
    )
  )
  expect_error(
    assert_intervention(
      test_intervention,
      population = "population" # population is wrong class
    )
  )

  # expect error when the number of demographic groups is incorrect
  expect_error(
    assert_intervention(
      test_intervention_bad,
      population = test_population
    )
  )
})
