# Basic tests to check for functionality
# Prepare some initial objects
init <- rbind(
  c(0.98, 0.01, 0.01, 0.0),
  c(0.9, 0.09, 0.01, 0.0),
  c(0.95, 0.04, 0.01, 0.0)
)
contact_matrix <- matrix(1.0, nrow = 3, ncol = 3)
contact_matrix <- contact_matrix / rowSums(contact_matrix)
diag(contact_matrix) <- rep(1.0, nrow(contact_matrix))

p <- population(
  contact_matrix = contact_matrix,
  demography_vector = rep(10, 3),
  initial_conditions = init
)

test_that("Basic tests for default epidemic model", {
  # create initial conditions

  data <- epidemic_default_cpp(
    population = p,
    beta = 0.5,
    alpha = 0.1,
    gamma = 0.05,
    time_end = 10, increment = 1.0
  )
  expect_type(data, "list")
  expect_named(
    data, c("x", "time")
  )
  # for an SEIR model
  expect_length(
    data[["x"]][[1]],
    length(init)
  )
})

test_that("Basic tests for data collection helper", {
  data <- epidemic_default_cpp(
    population = p,
    beta = 0.5,
    alpha = 0.1,
    gamma = 0.05,
    time_end = 200, increment = 1.0
  )
  data <- output_to_df(data)
  expect_named(
    data,
    c("time", do.call(paste0, expand.grid(c("S", "E", "I", "R"), seq(3)))),
    ignore.order = TRUE
  )
})
