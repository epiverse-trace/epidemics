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

# Prepare epi parameters as an infection object
pandemic <- infection(
  r0 = 1.5,
  preinfectious_period = 3,
  infectious_period = 7
)

# prepare a basic vaccination regime
elder_vaccination <- vaccination(
  name = "elder_vaccination",
  time_begin = matrix(0, 3, 1),
  time_end = matrix(200, 3, 1),
  nu = matrix(c(0, 0, 1e-4), 3, 1)
)

# prepare a three dose vaccination regime for a single age group
triple_vaccination <- vaccination(
  name = "triple_vaccination",
  nu = matrix(
    1e-4,
    nrow = 1, ncol = 3
  ),
  time_begin = matrix(
    seq(0, 30, 60),
    nrow = 1, ncol = 3
  ),
  time_end = matrix(
    seq(31, 61, 101),
    nrow = 1, ncol = 3
  )
)

# snapshot test for printing
test_that("Printing vaccination class", {
  expect_snapshot(elder_vaccination)
  expect_snapshot(triple_vaccination)
})

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
  expect_identical(
    nrow(elder_vaccination$time_begin),
    nrow(elder_vaccination$nu)
  )
  expect_identical(
    nrow(elder_vaccination$time_end),
    nrow(elder_vaccination$nu)
  )
})

# test that building a multi-dose vaccination works with c()
test_that("Multi-dose vaccination using `c()`", {
  vax_1 <- vaccination(
    name = "vax_regime",
    time_begin = matrix(1),
    time_end = matrix(100),
    nu = matrix(0.001)
  )

  # second dose regime
  vax_2 <- vaccination(
    name = "vax_regime",
    time_begin = matrix(101),
    time_end = matrix(200),
    nu = matrix(0.001)
  )

  expect_s3_class(
    c(vax_1, vax_2),
    "vaccination"
  )
  expect_snapshot(
    c(vax_1, vax_2)
  )

  # test for combining a two dose regime with another dose
  double_vax <- c(vax_1, vax_2)
  expect_s3_class(
    c(double_vax, vax_1),
    "vaccination"
  )
})

# run model with vaccination
data_vaccination <- epidemic_default_cpp(
  population = uk_population,
  infection = pandemic,
  vaccination = elder_vaccination,
  time_end = 200, increment = 1.0
)

# run model without vaccination
data <- epidemic_default_cpp(
  population = uk_population,
  infection = pandemic,
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

#### Test two dose no vaccination ####
doses <- 2L
no_vax_two_dose <- no_vaccination(uk_population, doses = doses)

test_that("No vaccination with two dose regime", {
  expect_identical(
    ncol(no_vax_two_dose$nu), doses
  )
  expect_identical(
    ncol(no_vax_two_dose$time_begin), doses
  )
  expect_identical(
    ncol(no_vax_two_dose$time_end), doses
  )
  expect_identical(
    unique(as.vector(no_vax_two_dose$nu)), 0.0,
    tolerance = 1e-6
  )
})
