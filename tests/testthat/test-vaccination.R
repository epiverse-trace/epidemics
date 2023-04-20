# Basic tests to check for functionality of vaccination class
# Prepare contact matrix and demography vector
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 40, 65),
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
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Prepare epi parameters
r0 <- 1.5
preinfectious_period <- 3
infectious_period <- 7

# prepare a basic vaccination regime
elder_vaccination <- vaccination(
  name = "elder_vaccination",
  time_begin = c(0, 0, 0),
  time_end = c(200, 200, 200),
  nu = c(0, 0, 1e-4)
)

# test the vaccination has expected structure
test_that("Vaccination is correctly initialised", {
  expect_s3_class(elder_vaccination, "vaccination")
  expect_named(
    elder_vaccination,
    c("name", "time_begin", "time_end", "nu")
  )
  expect_type(
    elder_vaccination$name, "character"
  )
  expect_length(
    elder_vaccination$name, 1L
  )
  expect_type(
    elder_vaccination$time_end, "double"
  )
  expect_type(
    elder_vaccination$time_begin, "double"
  )
  expect_type(
    elder_vaccination$nu, "double"
  )
  expect_length(
    elder_vaccination$time_end,
    length(elder_vaccination$time_begin)
  )
  expect_length(
    elder_vaccination$nu,
    length(elder_vaccination$time_begin)
  )
})

# run model with vaccination
data_vaccination <- epidemic(
  population = uk_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  vaccination = elder_vaccination,
  time_end = 200, increment = 1.0
)

# run model without vaccination
data <- epidemic(
  population = uk_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  time_end = 200, increment = 1.0
)

test_that("Epidemic model with vaccination", {

  # expect that only the last age group is vaccinated
  total_vaccinated <- data_vaccination[data_vaccination$compartment ==
    "vaccinated" & data_vaccination$time ==
    max(data_vaccination$time), ]$value

  expect_identical(
    total_vaccinated[seq(2)], rep(0.0, 2),
    ignore_attr = TRUE
  )
  expect_gt(
    total_vaccinated[3L], 0.0
  )
  expect_true(
    all(total_vaccinated[3L] > total_vaccinated[seq(2)])
  )

  # expect that vaccination reduces epidemic final size
  # test for the overall population
  final_size_vaccination <- epidemic_size(data_vaccination)
  final_size_default <- epidemic_size(data)

  expect_true(
    all(final_size_vaccination < final_size_default)
  )
})

# expect that no vaccination prints a message
test_that("Bad vaccination schedule prints a message", {
  expect_message(
    no_vaccination(uk_population),
    regexp = "(time_end)*(not greater than)*(time_begin)"
  )
})
