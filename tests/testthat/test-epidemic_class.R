# Tests for the <epidemic> class
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

test_that("<epidemic> class basic expectations", {
  # construction through a model run
  output <- model_default_cpp(uk_population)
  expect_snapshot(
    output
  )
  # check parameters
  expect_type(
    output[["parameters"]], "list"
  )
  # NOTE: hash during testing locally is NULL

  # manual construction
  data <- get_parameter(output, "data")
  expect_no_condition(
    epidemic(
      model_fn = "dummy_fn",
      data = data,
      parameters = list(
        dummy = "dummy_param"
      )
    )
  )

  # badly specified epidemic objects
  expect_error(
    epidemic(
      model_fn = model_default_cpp,
      data = data,
      parameters = list(
        dummy = "dummy_param"
      )
    ),
    regexp = "Must be of type 'string', not 'closure'"
  )

  expect_error(
    epidemic(
      model_fn = "dummy_name",
      data = data[, c(1, 2, 3)],
      parameters = list(
        dummy = "dummy_param"
      )
    ),
    regexp = "Must have exactly 4 cols, but has 3 cols"
  )

  expect_error(
    epidemic(
      model_fn = "dummy_name",
      data = data,
      parameters = "dummy_params"
    ),
    regexp = "Must be of type 'list'"
  )
})
