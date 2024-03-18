#### Tests for passing time-dependence functions ####

# Basic tests fo .no_time_dependence()
test_that("Time dependence: Basic expectations", {
  expect_identical(
    .no_time_dependence(),
    list(transmission_rate = function(time, x) x),
    ignore_function_env = TRUE
  )

  expect_identical(
    .no_time_dependence("recovery_rate"),
    list(recovery_rate = function(time, x) x),
    ignore_function_env = TRUE
  )
})

# NOTE: time-dependence functionality for models is tested in model test files
