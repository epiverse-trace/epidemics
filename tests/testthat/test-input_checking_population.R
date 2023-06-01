# Tests for input checking on population objects
# prepare a population
contact_matrix <- matrix(1, 2, 2)
demography_vector <- c(1e6, 1e6)

# Prepare some initial objects
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

test_that("Population are checked correctly", {
  # check for no conditions on a well formed population
  expect_no_condition(
    assert_population(
      test_population,
      compartments = read_from_library(
        model_name = "default", what = "compartments"
      )
    )
  )

  # expect failure when there is a mismatch in compartments
  # check for no conditions on a well formed population
  expect_error(
    assert_population(
      test_population,
      compartments = read_from_library(
        model_name = "vacamole", what = "compartments"
      )
    )
  )
})
