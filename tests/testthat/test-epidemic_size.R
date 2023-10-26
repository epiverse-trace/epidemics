# Test function for epidemic size
# create initial conditions
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 2, ncol = 5L,
  byrow = TRUE
)
colnames(initial_conditions) <- c(
  "susceptible", "exposed", "infectious", "recovered", "vaccinated"
)

# create a population where there are no contacts between groups
# this helps test the expectation that the final size proportions
# should be the same as the demographic group proportions
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(c(1, 0, 0, 1), 2, 2), # within group contacts only
  demography_vector = 67e6 * c(0.4, 0.6),
  initial_conditions = initial_conditions
)

# prepare an infection
pandemic <- infection(
  r0 = 1.5,
  preinfectious_period = 3,
  infectious_period = 7
)

# run epidemic simulation with no vaccination or intervention
data <- model_default_cpp(
  population = uk_population,
  infection = pandemic,
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

  # expect that the final size proportion is the same as the demography prop.
  expect_equal(
    epidemic_final_size / sum(epidemic_final_size),
    uk_population$demography_vector / sum(uk_population$demography_vector),
    tolerance = 1e-6,
    ignore_attr = TRUE
  )

  # test that final size is greater than size at 50% epidemic time
  epidemic_half_size <- epidemic_size(data, stage = 0.5)
  expect_true(
    all(epidemic_final_size > epidemic_half_size)
  )

  # test that by group FALSE returns a single value
  expect_length(
    epidemic_size(data, by_group = FALSE), 1
  )
})

#### Test that epidemic size with no deaths is same as deaths = FALSE
nonlethal_infect <- infection(
  name = "covid", r0 = 5, infectious_period = 10,
  preinfectious_period = 5,
  eta = 1 / 1000,
  omega = 0, # no deaths due to infection
  susc_reduction_vax = 0.5,
  hosp_reduction_vax = 0.7,
  mort_reduction_vax = 0.9
)

# make initial conditions - order is important
initial_conditions <- c(
  S = 1 - 1e-6,
  V1 = 0, V2 = 0,
  E = 0, EV = 0,
  I = 1e-6, IV = 0,
  H = 0, HV = 0, D = 0, R = 0
)
initial_conditions <- rbind(
  initial_conditions,
  initial_conditions
)

# create a population
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1, 2, 2),
  demography_vector = 67e6 * c(0.4, 0.6),
  initial_conditions = initial_conditions
)

test_that("Epidemic size with no deaths is correct", {
  data <- model_vacamole_cpp(
    population = uk_population,
    infection = nonlethal_infect,
    vaccination = no_vaccination(uk_population, doses = 2L),
    time_end = 400, increment = 1
  )
  expect_identical(
    epidemic_size(data, deaths = FALSE),
    epidemic_size(data)
  )
})
