test_that("Basic tests for default epidemic model", {
  data <- epidemic_default_cpp(
    init = c(0.98, 0.01, 0.01, 0.0),
    beta = 0.5,
    alpha = 0.1,
    gamma = 0.05
  )
  expect_type(data, "list")
  expect_named(
    data, c("x", "time")
  )
  # for an SEIR model
  expect_length(
    data[["x"]][[1]],
    4L
  )
})

test_that("Basic tests for data collection helper", {
  data <- epidemic_default_cpp(
    init = c(0.98, 0.01, 0.01, 0.0),
    beta = 0.5,
    alpha = 0.1,
    gamma = 0.05
  )
  data <- output_to_df(data)
  expect_named(
    data, c("time", "S", "E", "I", "R"),
    ignore.order = TRUE
  )
})
