# Tests for library helpers
test_that("multiplication works", {
  ml <- model_library()
  ml_names_manual <- ml[ml$model_type == "epidemic", ]$model_name

  ml_names <- get_model_names(model_type = "epidemic")

  expect_type(ml_names, "character")

  expect_identical(ml_names, ml_names_manual)

  # expect error when bad model_type specified
  expect_error(
    get_model_names(model_type = "some_type")
  )
})
