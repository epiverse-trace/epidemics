#### Tests for the ebola model ####

# prepare data
demography_vector <- 67e6

uk_pop <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(1 - 1e-6, 0, 1e-6, 0),
    nrow = 1, ncol = 4
  )
)

ebola <- infection(
  name = "ebolavirus disease",
  r0 = 1.7,
  infectious_period = 5,
  shape_E = 5L, rate_E = 1,
  shape_I = 5L, rate_I = 1
)

# basic expectations
test_that("Ebola model: basic expectations", {
  # runs without issues
  expect_no_condition(
    epidemic_ebola_cpp(
      population = uk_pop,
      infection = ebola,
      time_end = 100
    )
  )

  set.seed(1)
  # returns a data.table
  data <- epidemic_ebola_cpp(
    population = uk_pop,
    infection = ebola,
    time_end = 200
  )
  expect_s3_class(
    data, "data.table"
  )
  expect_length(data, 4L)
  expect_named(
    data, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  expect_identical(
    unique(data$compartment),
    read_from_library(model_name = "ebola", what = "compartments")
  )

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      data$value >= 0 & data$value <= sum(uk_pop$demography_vector)
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )

  # check that all age groups in the simulation are the same
  # size as the demography vector --- here, only one age group
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(uk_pop$contact_matrix)
  )
  expect_identical(
    rowSums(final_state),
    uk_pop$demography_vector,
    tolerance = 1e-6
  )

  # snaphshot test
  expect_snapshot(
    head(data)
  )
})

test_that("Larger R0 leads to larger final size in ebola model", {
  # prepare epidemic model runs with different R0 estimates
  r0_low <- 1.3
  r0_high <- 1.7
  infection_list <- list(
    ebola_r0_low = infection(
      r0 = r0_low, infectious_period = 5,
      shape_E = 5L, rate_E = 1,
      shape_I = 5L, rate_I = 1
    ),
    ebola_r0_high = infection(
      r0 = r0_high + 1.0, infectious_period = 5,
      shape_E = 5L, rate_E = 1,
      shape_I = 5L, rate_I = 1
    )
  )

  # get data
  data <- lapply(
    infection_list,
    function(infection_) {
      # run model on data
      data <- epidemic_ebola_cpp(
        population = uk_pop,
        infection = infection_,
        time_end = 100
      )
    }
  )

  # get final size as total recoveries
  final_sizes <- lapply(data, epidemic_size)

  # test for effect of R0
  expect_true(
    all(final_sizes[["r0_high"]] > final_sizes[["r0_low"]])
  )
})

# Equivalence with R model
test_that("Ebola model equivalence in R-only and RCpp", {
  set.seed(1)
  max_time <- 100
  ebola_r <- epidemic_ebola_r(
    initial_state = as.integer(uk_pop$demography_vector *
      uk_pop$initial_conditions),
    parameters = c(
      shape_E = 5L, rate_E = 1,
      shape_I = 5L, rate_I = 1,
      beta = ebola$r0 / ebola$infectious_period
    ),
    max_time = max_time
  )
  # values at last timestep
  ebola_r_values <- as.integer(
    t(
      as.matrix(ebola_r[max_time, c("S", "E", "I", "R")])
    )
  )

  # set seed for Rcpp run
  set.seed(1)
  ebola_cpp <- epidemic_ebola_cpp(
    population = uk_pop,
    infection = ebola,
    time_end = max_time - 1 # one less for C++ implementation due to zero index
  )
  # get last timestep values
  ebola_cpp_values <- tail(ebola_cpp$value, 4L)

  # expect identical
  expect_identical(
    ebola_cpp_values,
    ebola_r_values
  )
})
