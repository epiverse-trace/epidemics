# Tests for reading model library JSON and for
# reading functions from the model library provided by model_library()

ml_manual <- jsonlite::read_json(
  system.file(
    "extdata", "model_library.json",
    package = "epidemics"
  ),
  simplifyVector = TRUE
)

test_that("Reading model library JSON", {
  ml <- model_library()

  expect_identical(
    ml, ml_manual
  )
  expect_snapshot(ml)
})

test_that("Reading functions from library", {
  model_function <- read_from_library(
    model_type = "epidemic",
    model_name = "default",
    what = "model_function"
  )

  args_checker <- read_from_library(
    model_type = "epidemic",
    model_name = "default",
    what = "model_args_checker"
  )

  model_fn_manual <- ml_manual[
    ml_manual$model_type == "epidemic" &
      ml_manual$model_name == "default",
  ][["model_function"]]

  args_checker_manual <- ml_manual[
    ml_manual$model_type == "epidemic" &
      ml_manual$model_name == "default",
  ][["model_args_checker"]]

  # test that the default epidemic model function is correctly returned
  expect_identical(
    model_function,
    model_fn_manual
  )

  # test that the default epidemic model argument checker is correctly returned
  expect_identical(
    args_checker,
    args_checker_manual
  )

  # test that the default model compartments are returned correctly
  compartments_default <- read_from_library(what = "compartments")
  expect_identical(
    compartments_default,
    unlist(ml_manual$compartments[ml_manual$model_name == "default"])
  )

  # Test that requesting a model not in the library errors
  expect_error(
    read_from_library(
      model_type = "epidemic",
      model_name = "not-the-default",
      what = "model_function"
    ),
    regexp = "(No model named)*(not-the-default)*(epidemic)*(check type-name)"
  )
})
