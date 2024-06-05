#' Run the default model using UK population data from polymod and
#' demography groups (0, 20, 40)
test_data_default_model <- function() {
  polymod <- socialmixr::polymod
  contact_data <- socialmixr::contact_matrix(
    polymod,
    countries = "United Kingdom",
    age.limits = c(0, 20, 40),
    symmetric = TRUE
  )

  # prepare contact matrix
  contact_matrix <- t(contact_data[["matrix"]])

  # prepare the demography vector
  demography_vector <- contact_data[["demography"]][["population"]]
  names(demography_vector) <- rownames(contact_matrix)


  # initial conditions: one in every 1 million is infected
  initial_i <- 1e-6
  initial_conditions <- c(
    S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
  )

  # build for all age groups
  initial_conditions <- rbind(
    initial_conditions,
    initial_conditions,
    initial_conditions
  )
  rownames(initial_conditions) <- rownames(contact_matrix)

  # prepare the population to model as affected by the epidemic
  uk_population <- population(
    name = "UK",
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions
  )

  # run an epidemic model using `epidemic()`
  output <- model_default(
    population = uk_population,
    transmission_rate = 1.5 / 7.0,
    infectiousness_rate = 1.0 / 3.0,
    recovery_rate = 1.0 / 7.0,
    # intervention = list(contacts = close_schools),
    time_end = 600, increment = 1.0
  )

  output
}

#' Run the vacamole model with no demography groups and a double vaccination regime
test_data_vacamole <- function() {
  # create a population, note eleven columns for compartments
  population <- population(
    contact_matrix = matrix(1),
    demography_vector = 67e6,
    initial_conditions = matrix(
      c(0.9999, 0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0),
      nrow = 1, ncol = 11L
    )
  )

  # create a vaccination regime
  double_vax <- vaccination(
    nu = matrix(1e-3, ncol = 2, nrow = 1),
    time_begin = matrix(c(10, 30), nrow = 1),
    time_end = matrix(c(50, 80), nrow = 1)
  )

  # data <-
  model_vacamole(
    population = population,
    vaccination = double_vax
  )
}

#' Run the ebola model with no demography groups and a double vaccination regime
test_data_ebola <- function() {
  # create a population with 6 compartments
  population <- population(
    contact_matrix = matrix(1),
    demography_vector = 14e6,
    initial_conditions = matrix(
      c(0.999998, 0.000001, 0.000001, 0, 0, 0),
      nrow = 1, ncol = 6L
    )
  )

  model_ebola(
    population = population
  )
}

test_that("Get the epidemic peak and size for the default model with 3 demography groups", {
  # run epidemic model, expect no condition
  data <- test_data_default_model()
  expect_no_condition(
    epidemic_peak(data)
  )

  out <- epidemic_peak(data)


  # check for output type and contents
  expect_s3_class(out, "data.frame")


  expect_length(data, 4L)
  expect_named(
    out, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )

  expect_identical(unique(out$compartment), "infectious")

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      out[["value"]] <= data[data[["compartment"]] == "susceptible" & data[["time"]] == 0, value]
    )
  )
})

test_that("Get the epidemic peak and size for the vacamole model", {
  # run epidemic model, expect no condition
  data <- test_data_vacamole()
  expect_no_condition(
    epidemic_peak(data)
  )
  out <- epidemic_peak(data)

  # check for output type and contents
  expect_s3_class(out, "data.frame")
  expect_length(data, 4L)
  expect_named(
    out, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )

  expect_identical(unique(out$compartment), "infectious")

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      out[["value"]] <= data[data[["compartment"]] == "susceptible" & data[["time"]] == 0, value]
    )
  )
})

test_that("Get the epidemic peak and size for the ebola model", {
  # run epidemic model, expect no condition
  data <- test_data_ebola()
  expect_no_condition(
    epidemic_peak(data)
  )
  out <- epidemic_peak(data)

  # check for output type and contents
  expect_s3_class(out, "data.frame")
  expect_length(data, 5L)
  expect_named(
    out, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )

  expect_identical(unique(out$compartment), "infectious")

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      out[["value"]] <= data[data[["compartment"]] == "susceptible" & data[["time"]] == 0, value]
    )
  )
})
