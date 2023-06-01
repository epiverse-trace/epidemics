# Tests for input checking on infection objects
test_that("Infections are checked correctly", {
  # prepare a well formed infection object for the default model
  infection_default <- infection(
    name = "influenza", r0 = 1.3, infectious_period = 10,
    preinfectious_period = 3, mortality_rate = 1e-4
  )

  # expect no errors for well formed infection-assertion matches
  expect_no_condition(
    assert_infection(
      infection_default,
      extra_parameters = c("preinfectious_period", "mortality_rate")
    )
  )

  expect_no_condition(
    assert_infection(
      infection_default,
      extra_parameters = c("preinfectious_period", "mortality_rate"),
      extra_parameters_limits = list(
        preinfectious_period = c(
          "lower" = 0, "upper" = 5
        )
      )
    )
  )

  # expect errors when the infection does not pass assertions
  # 1. expect error when the extra arguments are missing from the object
  # dummy extra argument that is missing
  missing_arg <- "some_extra_arg"
  # exact regex is not added here as it is too long
  expect_error(
    assert_infection(
      infection_default,
      extra_parameters = c(
        "preinfectious_period", "mortality_rate", missing_arg
      )
    )
  )

  # expect error also when the `infection` has more parameters than expected
  expect_error(
    assert_infection(
      infection_default,
      extra_parameters = "preinfectious_period"
    )
  )

  # 2. expect error when argument is present but has objectively bad value
  infection_malformed <- infection(
    name = "influenza", r0 = 1.3, infectious_period = 10,
    preinfectious_period = -5 # time cannot be negative
  )
  expect_error(
    assert_infection(
      infection_malformed,
      extra_parameters = "preinfectious_period"
    )
  )

  # 3. expect error when argument present but outside limit
  expect_error(
    assert_infection(
      infection_default,
      extra_parameters = "preinfectious_period",
      extra_parameters_limits = list(
        preinfectious_period = c(
          lower = 4, upper = 5
        )
      )
    )
  )

  # 4. expect error when argument limits are passed as a single vector
  expect_error(
    assert_infection(
      infection_default,
      extra_parameters = "preinfectious_period",
      extra_parameters_limits = c(
        lower = 4, upper = 5
      )
    )
  )

  # 5. expect error when argument limits are poorly named
  expect_error(
    assert_infection(
      infection_default,
      extra_parameters = "preinfectious_period",
      extra_parameters_limits = list(
        preinfectious_period = c(
          lower_limit = 4, higher_limit = 5
        )
      )
    )
  )
})
