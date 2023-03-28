# Basic tests to check for functionality
# Prepare contact matrix and demography vector
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
    c(0.9999, 0.0001, 0, 0),
    nrow = nrow(contact_matrix), ncol = 4,
    byrow = TRUE
  )
)

# Prepare epi parameters
r0 <- rep(1.5, nrow(uk_population$contact_matrix))
preinfectious_period <- rep(3, nrow(uk_population$contact_matrix))
infectious_period <- rep(7, nrow(uk_population$contact_matrix))

test_that("Output of default epidemic model", {
  # run model on data
  data <- epidemic_cpp(
    population = uk_population,
    r0 = r0,
    intervention = no_intervention(uk_population),
    preinfectious_period = preinfectious_period,
    infectious_period = infectious_period,
    time_end = 100, increment = 1.0
  )

  # check for output type and contents
  expect_s3_class(data, "data.table")
  expect_length(data, length(uk_population$initial_conditions) + 1L)

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      vapply(subset(data, select = -time),
        FUN = function(x) {
          all(x >= 0 & x <= sum(uk_population$demography_vector))
        },
        FUN.VALUE = TRUE
      )
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(subset(data[1, ], select = -time)),
    sum(subset(tail(data, 1L), select = -time)),
    tolerance = 1e-6
  )
})

test_that("Larger R0 leads to larger final size in default epidemic model", {
  # prepare epidemic model runs with different R0 estimates
  r0_list <- list(
    r0_low = r0,
    r0_high = r0 + 1.0
  )

  # get data
  data <- lapply(
    r0_list,
    function(r0_) {
      # run model on data
      data <- epidemic_cpp(
        population = uk_population,
        r0 = r0_,
        preinfectious_period = preinfectious_period,
        infectious_period = infectious_period,
        intervention = no_intervention(uk_population),
        time_end = 10, increment = 1.0
      )
    }
  )

  # get final size as total recoveries
  final_sizes <- lapply(data, tail, 1)
  final_sizes <- vapply(final_sizes, `[[`, FUN.VALUE = 1, "recovered_1")

  # test for effect of R0
  expect_gt(
    final_sizes["r0_high"],
    final_sizes["r0_low"]
  )
})
