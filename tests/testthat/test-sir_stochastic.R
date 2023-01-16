#### epi_demic() with sir_stochastic() ####
test_that("epi_demic runs withour errors", {
  expect_silent(
    epi_demic(
      parameters = make_parameters_sir_stochastic(),
      option = "sir_stochastic"
    )
  )
})

# check statistical correctness
test_that("epi_demic runs withour errors", {
  outcome <- epi_demic(
    parameters = make_parameters_sir_stochastic(),
    option = "sir_stochastic"
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
      parameters = make_parameters_sir_stochastic(),
      option = "invalid option"
    ),
    regexp = "Option not found, option must be 'sir_stochastic'"
  )
})
