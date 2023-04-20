# Test function for epidemic size
# create initial conditions
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 1, ncol = 5L
)
colnames(initial_conditions) <- read_from_library(what = "compartments")

# create a population
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = initial_conditions
)

# run epidemic simulation with no vaccination or intervention
data <- epidemic(
  model_name = "default",
  population = uk_population,
  r0 = 1.5,
  preinfectious_period = 3,
  infectious_period = 7,
  time_end = 200,
  increment = 1
)

test_that("Epidemic size functions", {
  # test for the initial size = 0.0
  epidemic_initial_size <- epidemic_size(data, stage = 0.0)
  expect_equal(
    epidemic_initial_size,
    initial_conditions[, "infectious"],
    ignore_attr = TRUE
  )

  # test the final size
  epidemic_final_size <- epidemic_size(data)
  expect_equal(
    epidemic_final_size,
    data[data$compartment == "recovered" & data$time == max(data$time), ]$value,
    ignore_attr = TRUE
  )

  # test that final size is greater than size at 50% epidemic time
  epidemic_half_size <- epidemic_size(data, stage = 0.5)
  expect_gt(
    epidemic_final_size, epidemic_half_size
  )

  # test that by group FALSE returns a single value
  expect_length(
    epidemic_size(data, by_group = FALSE), 1
  )
})
