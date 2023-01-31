#### epi_demic() with sir_stochastic() ####
test_that("epi_demic runs without errors", {
  expect_silent(
    epi_demic(
      parameters = make_parameters_sir(option = "stochastic"),
      option = "stochastic"
    )
  )
})

# check statistical correctness
test_that("epi_demic stochastic is statistically correct", {
  outcome <- epi_demic(
    parameters = make_parameters_sir(option = "stochastic"),
    option = "stochastic"
  )

  expect_lte(
    max(outcome$value), 1.0
  )
  expect_gte(
    min(outcome$value), 0.0
  )

  # check that all elements sum to 1.0 at the final timepoint
  max_time <- outcome[outcome$time == max(outcome$time), ]
  expect_identical(
    sum(max_time$value), 1.0,
    tolerance = 1e-6
  )
})

# check failure for other options for now
test_that("epi_demic() with invalid option", {
  expect_error(
    epi_demic(
      parameters = make_parameters_sir(option = "stochastic"),
      option = "invalid option"
    ),
    regexp = "('arg' should be one of)*(stochastic)*(deterministic)"
  )
})

#### epi_demic() with sir_desolve() ####
test_that("epi_demic deterministic runs without errors", {
  expect_silent(
    epi_demic(
      parameters = make_parameters_sir(option = "deterministic"),
      option = "deterministic"
    )
  )
})

# check statistical correctness
test_that("epi_demic determinsitic is statistially correct", {
  outcome <- epi_demic(
    parameters = make_parameters_sir(option = "deterministic"),
    option = "deterministic"
  )

  expect_lte(
    max(outcome$value), 1.0
  )
  expect_gte(
    min(outcome$value), 0.0
  )

  # check that all elements sum to 1.0 at the final timepoint
  max_time <- outcome[outcome$time == max(outcome$time), ]
  expect_identical(
    sum(max_time$value), 1.0,
    tolerance = 1e-6
  )
})
