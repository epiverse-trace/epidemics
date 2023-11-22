# Tests for a model with multiple stacked interventions
# Prepare population
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
    c(0.9999, 0.0001, 0, 0, 0),
    nrow = nrow(contact_matrix), ncol = 5L,
    byrow = TRUE
  )
)

# Prepare a stacked intervention
npi_1 <- intervention(
  type = "contacts",
  time_begin = 30,
  time_end = 60,
  reduction = matrix(0.15, nrow = nrow(contact_matrix))
)
npi_2 <- intervention(
  type = "contacts",
  time_begin = 45,
  time_end = 75,
  reduction = matrix(0.1, nrow = nrow(contact_matrix))
)
multi_npi <- c(npi_1, npi_2)

# Tests for the value of contact reduction from overlapping NPIs
test_that("Cumulative effect of NPIs", {
  cumulative_cr <- cumulative_contacts_intervention(
    t = 50, time_begin = multi_npi$time_begin, time_end = multi_npi$time_end,
    reduction = multi_npi$reduction
  )
  expect_identical(
    cumulative_cr,
    rep(0.25, nrow(contact_matrix))
  )

  # when no intervention is active
  expect_identical(
    cumulative_contacts_intervention(
      t = 10, time_begin = multi_npi$time_begin, time_end = multi_npi$time_end,
      reduction = multi_npi$reduction
    ),
    numeric(nrow(contact_matrix))
  )

  # expect contact matrix has scaled values
  expect_identical(
    intervention_on_cm(
      t = 50,
      cm = contact_matrix,
      time_begin = multi_npi$time_begin,
      time_end = multi_npi$time_end,
      cr = multi_npi$reduction
    ),
    (contact_matrix * 0.75) %*% diag(rep(0.75, nrow(contact_matrix)))
  )
})

# Initial test that a model with two interventions works
test_that("Default model with multiple interventions", {
  # run epidemic model
  expect_no_condition(
    model_default_cpp(
      population = uk_population,
      intervention = list(contacts = multi_npi),
      time_end = 100
    )
  )

  # expect final size is smaller with stacked interventions
  data_1_npi <- model_default_cpp(
    population = uk_population,
    intervention = list(contacts = npi_2),
    time_end = 100
  )
  data_2_npi <- model_default_cpp(
    population = uk_population,
    intervention = list(contacts = multi_npi),
    time_end = 100
  )
  expect_lte(
    epidemic_size(data_2_npi, by_group = FALSE, include_deaths = FALSE),
    epidemic_size(data_1_npi, by_group = FALSE, include_deaths = FALSE)
  )
})
