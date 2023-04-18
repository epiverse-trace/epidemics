# Tests for adding functions to the model library provided by model_library()

# read the model library using the model library function
ml_manual <- model_library()

test_that("Add model function and details to the library", {

  # expect message when writing
  expect_message(
    add_to_library(
      model_type = "test-type",
      model_name = "test-name"
    ),
    regexp = paste(
      "(Adding 'test-type' type model)*('test-name')*(model library)*",
      "(Please check that this is correct before commiting your changes.)",
      sep = "\n"
    )
  )

  # check that the test model functions exist
  test_fn <- read_from_library(
    model_type = "test-type",
    model_name = "test-name",
    which_function = "model_function"
  )
  expect_identical(
    test_fn,
    "test-type_test-name_cpp"
  )
})

# IMPORTANT: write original model library back to file
jsonlite::write_json(
  ml_manual,
  system.file("extdata", "model_library.json", package = "epidemics")
)
