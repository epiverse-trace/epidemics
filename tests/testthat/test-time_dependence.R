#### Tests for passing time-dependence functions ####

# Basic tests fo no_time_dependence()
test_that("Basic test for no_time_dependence", {
  expect_identical(
    no_time_dependence(),
    list(transmissibility = function(time, x) x),
    ignore_function_env = TRUE
  )
})

# set up initial conditions
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# Prepare some initial objects
uk_population <- population(
  name = "UK population",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(0.9999, 0, 0.0001, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# prepare function
mod_transmissibility <- function(time, x, tmax = 365 / 4) {
  x + (x * sinpi(time / tmax))
}

#### Tests for functioning and equivalence ####
test_that("Basic expectations for time dependence functions", {
  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      time_dependence = list(
        transmissibility = mod_transmissibility
      ),
      time_end = 365, increment = 1
    )
  )

  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      time_dependence = list(
        transmissibility = mod_transmissibility
      ),
      time_end = 365, increment = 1
    )
  )

  output_cpp <- model_default_cpp(
    population = uk_population,
    time_dependence = list(
      transmissibility = mod_transmissibility
    ),
    time_end = 365, increment = 1
  )
  output_r <- model_default_cpp(
    population = uk_population,
    time_dependence = list(
      transmissibility = mod_transmissibility
    ),
    time_end = 365, increment = 1
  )

  expect_identical(
    tail(get_parameter(output_cpp, "data"), 20),
    tail(get_parameter(output_r, "data"), 20)
  )
})
