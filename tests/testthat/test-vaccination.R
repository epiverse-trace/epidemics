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
r0 <- rep(1.5, nrow(contact_matrix))
preinfectious_period <- rep(3, nrow(contact_matrix))
infectious_period <- rep(7, nrow(contact_matrix))

# prepare a basic vaccination regime
elder_vaccination <- vaccination(
  name = "elder_vaccination",
  time_begin = c(0, 0, 0),
  time_end = c(200, 200, 200),
  nu = c(0, 0, 1e-4)
)

# run model with vaccination
data_vaccination <- epidemic_cpp(
  population = uk_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  vaccination = elder_vaccination,
  time_end = 200, increment = 1.0
)

# run model without vaccination
data <- epidemic_cpp(
  population = uk_population,
  r0 = r0,
  preinfectious_period = preinfectious_period,
  infectious_period = infectious_period,
  time_end = 200, increment = 1.0
)

test_that("Epidemic model with vaccination", {

  # expect that only the last age group is vaccinated
  total_vaccinated <- unlist(
    subset(
      tail(data_vaccination, 1),
      select = grepl("vaccinated", colnames(data_vaccination), fixed = TRUE)
    )
  )
  expect_identical(
    total_vaccinated[seq(2)], rep(0.0, 2),
    ignore_attr = TRUE
  )
  expect_gt(
    total_vaccinated["vaccinated_3"], 0.0
  )
  expect_true(
    all(total_vaccinated["vaccinated_3"] > total_vaccinated[seq(2)])
  )

  # expect that vaccination reduces epidemic final size
  # test for the overall population
  final_size_vaccination <- unlist(
    subset(
      tail(data_vaccination, 1),
      select = grepl("recovered", colnames(data_vaccination), fixed = TRUE)
    )
  )
  final_size_default <- unlist(
    subset(
      tail(data, 1),
      select = grepl("recovered", colnames(data_vaccination), fixed = TRUE)
    )
  )

  expect_true(
    all(final_size_vaccination < final_size_default)
  )
})
