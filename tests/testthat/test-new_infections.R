# Tests for calculation of new infections
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 1, ncol = 5L,
  byrow = TRUE
)
colnames(initial_conditions) <- read_from_library(what = "compartments")

# create a population where there are no contacts between groups
# this helps test the expectation that the final size proportions
# should be the same as the demographic group proportions
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1), # within group contacts only
  demography_vector = 67e6,
  initial_conditions = initial_conditions
)

# add a vaccination regime
vax_time_begin <- 50
vaccination <- vaccination(
  time_begin = vax_time_begin, time_end = 200, nu = 0.0005
)

# run epidemic simulation with vaccination
time_end <- 200
increment <- 1
data <- epidemic(
  model_name = "default",
  population = uk_population,
  r0 = 2.0,
  preinfectious_period = 3,
  infectious_period = 7,
  vaccination = vaccination,
  time_end = time_end,
  increment = increment
)

test_that("New infections are correctly caculated", {
  data_ <- new_infections(data, compartments_from_susceptible = "vaccinated")

  # expect correct number of rows
  expect_identical(
    nrow(data_),
    nrow(data) + (length(seq(0, time_end, by = increment)))
  )

  # expect columns with specific names
  expect_named(
    data_,
    c("time", "demography_group", "compartment", "value"),
    ignore.order = TRUE
  )

  # expect new variable
  expect_true(
    "new_infections" %in% unique(data_$compartment)
  )

  # test for correctness
  # test that new infections are exactly the same as change in susceptibles
  # - change in vaccinations
  new_infections <- data_[
    data_$compartment == "new_infections",
  ][["value"]]

  # note use of original data.table
  delta_susc <- data[data$compartment == "susceptible", ][["value"]]
  delta_susc <- c(0, -diff(delta_susc))

  delta_vax <- data[data$compartment == "vaccinated", ][["value"]]
  delta_vax <- c(0, diff(delta_vax))

  expect_identical(
    new_infections,
    delta_susc - delta_vax
  )
})

test_that("New infections without accounting for vaccination", {
  data_ <- new_infections(data)
  new_infections <- data_[
    data_$compartment == "new_infections",
  ][["value"]]

  # note use of original data.table
  delta_susc <- data[data$compartment == "susceptible", ][["value"]]
  delta_susc <- c(0, -diff(delta_susc))

  delta_vax <- data[data$compartment == "vaccinated", ][["value"]]
  delta_vax <- c(0, diff(delta_vax))

  # expect that new infections are exactly identical to change in susceptibles
  expect_identical(
    new_infections,
    delta_susc
  )
})
