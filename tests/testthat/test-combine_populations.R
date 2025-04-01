# Basic tests to check for functionality of combining population
# Prepare contact matrix and demography vector of the first population
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 40, 65),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

# Create the first population
population1 <- population(
  name = "Population 1",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Create the second population
population2 <- population(
  name = "Population 2",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(0.9999, 0.00005, 0.00005, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)
# Set connectivity matrix between population1 and population2
prop_matrix <- matrix(c(1, 0.05, 0.05, 1), nrow = 2, ncol = 2)

combined_population <- combine_populations(
  populations = list(population1, population2),
  connectivity_matrix = prop_matrix,
  method = "linear", name = "combine"
)

# snapshot test for printing
test_that("Printing population class", {
  expect_snapshot(combined_population)
  
  # population with no rownames for the contact matrix
  names(combined_population$demography_vector) <- NULL
  rownames(combined_population$contact_matrix) <- NULL
})

# test the population has expected structure
test_that("Combined population is correctly initialised", {
  expect_s3_class(combined_population, "population")
  expect_named(
    combined_population,
    c("name", "contact_matrix", "demography_vector", "initial_conditions")
  )
  expect_type(
    combined_population$name, "character"
  )
  expect_length(
    combined_population$name, 1L
  )
  expect_type(
    combined_population$contact_matrix, "double"
  )
  expect_identical(
    nrow(combined_population$contact_matrix),
    nrow(population1$contact_matrix) + nrow(population2$contact_matrix)
  )
  expect_identical(
    nrow(combined_population$contact_matrix),
    length(combined_population$demography_vector)
  )
  expect_identical(
    nrow(combined_population$initial_conditions),
    nrow(combined_population$contact_matrix)
  )
})

# Set connectivity matrix between population1 and population2
dist_matrix <- matrix(c(1, 0.05, 0.05, 1), nrow = 2, ncol = 2)

combined_population_gravity <- combine_populations(
  populations = list(population1, population2),
  connectivity_matrix = dist_matrix,
  method = "gravity", name = "combine"
)

# expect that setting the method with "gravity", or with the gravity contact function
# leads to the same result
test_that("Using the function directly leads to the same result", {
  combined_population_function <- combine_populations(
    populations = list(population1, population2),
    connectivity_matrix = dist_matrix,
    method = gravity_contact, name = "combine"
  )
  
  expect_equal(
    combined_population_gravity,
    combined_population_function
  )
})
