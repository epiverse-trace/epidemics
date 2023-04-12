# Basic tests to check for functionality of population class
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

# a population with only 50% the size
half_population <- population(
  name = "UK population",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector * 0.5,
  initial_conditions = matrix(
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Prepare epi parameters
r0 <- rep(1.5, nrow(contact_matrix))
preinfectious_period <- rep(3, nrow(contact_matrix))
infectious_period <- rep(7, nrow(contact_matrix))

# run model with different populations
data_full_pop <- epidemic_cpp(
  population = uk_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  time_end = 200, increment = 1.0
)

data_half_pop <- epidemic_cpp(
  population = half_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  time_end = 200, increment = 1.0
)

# expect that smaller population has a smaller final size
test_that("Effect of population on final size", {
  # get final sizes
  final_size_full <- unlist(
    subset(
      tail(data_full_pop, 1),
      select = grepl("recovered", colnames(data_full_pop), fixed = TRUE)
    )
  )
  final_size_half <- unlist(
    subset(
      tail(data_half_pop, 1),
      select = grepl("recovered", colnames(data_half_pop), fixed = TRUE)
    )
  )

  # expect that the smaller population has a smaller final size
  expect_true(
    all(final_size_half < final_size_full)
  )
})