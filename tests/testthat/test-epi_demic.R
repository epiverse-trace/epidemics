#### Test for output and statistical correctness of an SIRV model ####

# with a single unstructured population
population <- population(
  contact_matrix = matrix(1),
  initial_conditions = matrix(
    c(0.999, 0.001, 0, 0),
    nrow = 1, ncol = 4
  )
)

# with a single intervention
intervention <- intervention(
  time_begin = 50,
  time_end = 80,
  contact_reduction = 0.2
)

# no vaccination in this model
output <- epi_demic(
  population = population,
  intervention = intervention,
  nu = 0.0,
  t_increment = 1
)

# check for correct form of output
test_that("SIRV model returns data.table with correct form", {
  expect_s3_class(
    output, c("data.frame", "data.table")
  )
  expect_named(
    output,
    c("time", "compartment", "age_group", "value")
  )
  expect_snapshot(
    head(output)
  )
})

#### Some basic tests of statistical correctness ####
# note that these are not exhaustive
test_that("SIRV model is statistically correct", {
  # FIXME: add more tests for correctness
  # add test for equality of population sizes

  # check that a zero vaccination rate leads to no vaccinated individuals
  expect_identical(
    unique(output$value[output$compartment == "vaccinated"]),
    0.0,
    tolerance = 1e-5
  )
})

#### with an age structured population ####
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)
contact_matrix <- t(contact_data$matrix)
demography_vector <- contact_data$demography$population

population <- population(
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = rbind(
    c(0.999, 0.001, 0, 0),
    c(0.9999, 0.0001, 0, 0), # lowest infection among 20-40
    c(0.999, 0.001, 0, 0)
  )
)

# with a single intervention
intervention <- intervention(
  time_begin = 50,
  time_end = 80,
  contact_reduction = c(0.2, 0.001, 0.0)
)

# R0 varies by age
r0 <- c(1.5, 1.3, 1.2)

# no vaccination in this model
output <- epi_demic(
  R0 = r0,
  population = population,
  intervention = intervention,
  nu = 0.0,
  t_increment = 1
)

#### Basic tests on multiple age groups ####
test_that("SIRV model with multiple age groups", {
  # multiple age groups are returned correctly
  expect_identical(
    unique(output$age_group),
    seq_along(demography_vector)
  )

  # FIXME: add more tests for correctness
  # add test for equality of population sizes

  # check that a zero vaccination rate leads to no vaccinated individuals
  expect_identical(
    unique(output$value[output$compartment == "vaccinated"]),
    0.0,
    tolerance = 1e-5
  )
})
