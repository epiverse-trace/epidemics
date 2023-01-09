#### Test for output and statistical correctness of SIR model using deSolve ####
test_that("SIR deSolve model returns data.table with correct form", {
  time_seq <- seq(0, 100, length.out = 101)
  output <- sir_desolve(times = time_seq)
  expect_s3_class(
    output, "data.table"
  )
  expect_identical(
    colnames(output),
    c("time", "state", "proportion")
  )
  expect_snapshot(
    head(output)
  )
  # expect correct number of rows
  expect_identical(
    nrow(output), length(time_seq) * 3L
  )
})

#### Some basic tests of statistical correctness ####
# not that these are not exhaustive
test_that("SIR deSolve model is statistically correct", {
  time_seq <- seq(0, 1000)

  # check that no initial infections lead to a final size of 0.0
  output <- sir_desolve(
    times = time_seq,
    initial_conditions(p_initially_infected = 0.0)
  )
  final_size <- sum(output[state == "R", ][["proportion"]])
  expect_identical(
    final_size, 0.0
  )

  # check that setting gamma = 0 leads to all individuals in the I class
  output <- sir_desolve(
    times = time_seq,
    initial_conditions(p_initially_infected = 0.01),
    parms = c(beta = 0.1, gamma = 0.0)
  )
  recovered <- data.table::last(output[state == "R", ][["proportion"]])
  infected <- data.table::last(output[state == "I", ][["proportion"]])
  expect_identical(
    infected, 1.0,
    tolerance = 1e-5
  )
})
