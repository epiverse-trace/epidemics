#### Test for output and statistical correctness of SIR model using deSolve ####
test_that("SIR deSolve model returns data.table with correct form", {
  parameters <- make_parameters_sir(option = "deterministic")
  output <- do.call(sir_desolve, parameters)
  expect_vector(
    output, list()
  )
  expect_named(
    output,
    c("time", "S", "I", "R")
  )
  expect_snapshot(
    str(output)
  )
})

#### Some basic tests of statistical correctness ####
# not that these are not exhaustive
test_that("SIR deSolve model is statistically correct", {
  parameters <- make_parameters_sir(option = "deterministic")
  parameters$I <- 0.0

  # check that no initial infections lead to a final size of 0.0
  output <- do.call(sir_desolve, parameters)
  final_size <- sum(output$R)
  expect_identical(
    final_size, 0.0
  )

  # check that setting gamma = 0 leads to all individuals in the I class, i.e., no recoveries
  parameters$I <- 0.01
  parameters$gamma <- 0.0
  output <- do.call(sir_desolve, parameters)
  recovered <- tail(output$R, 1)
  infected <- tail(output$I, 1)
  expect_identical(
    recovered, 0.0,
    tolerance = 1e-5
  )
  expect_identical(
    infected, 1.0,
    tolerance = 1e-5
  )
})
