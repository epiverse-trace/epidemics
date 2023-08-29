#### Tests for passing time-dependence functions ####

# Basic tests fo no_time_dependence()
test_that("Basic test for no_time_dependence", {
  expect_identical(
    no_time_dependence(),
    list(beta = function(time, x) x),
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

# Prepare epi parameters
pandemic <- infection(
  r0 = 1.1,
  preinfectious_period = 3,
  infectious_period = 7
)

# prepare function
mod_beta <- function(time, x, tmax = 365 / 4) x + (x * sinpi(time / tmax))

#### Tests for functioning and equivalence ####
test_that("Basic expectations for time dependence functions", {
  expect_no_condition(
    epidemic_default_cpp(
      population = uk_population,
      infection = pandemic,
      time_dependence = list(
        beta = mod_beta
      ),
      time_end = 365, increment = 1
    )
  )

  expect_no_condition(
    epidemic_default_cpp(
      population = uk_population,
      infection = pandemic,
      time_dependence = list(
        beta = mod_beta
      ),
      time_end = 365, increment = 1
    )
  )

  data_cpp <- epidemic_default_cpp(
    population = uk_population,
    infection = pandemic,
    time_dependence = list(
      beta = mod_beta
    ),
    time_end = 365, increment = 1
  )
  data_r <- epidemic_default_cpp(
    population = uk_population,
    infection = pandemic,
    time_dependence = list(
      beta = mod_beta
    ),
    time_end = 365, increment = 1
  )

  expect_identical(
    tail(data_r, 20),
    tail(data_cpp, 20),
  )
})
