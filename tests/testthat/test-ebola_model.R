#### Tests for the ebola model ####

# prepare data
demography_vector <- 67000

pop <- population(
  contact_matrix = matrix(1),
  demography_vector = demography_vector,
  initial_conditions = matrix(
    c(1 - 1e-3, 1e-3 / 2, 1e-3 / 2, 0, 0, 0),
    nrow = 1
  )
)

ebola <- infection(
  name = "ebolavirus disease",
  r0 = 1.7,
  infectious_period = 7,
  preinfectious_period = 5,
  prop_hospitalised = 0.5,
  etu_safety = 0.8,
  funeral_safety = 0.5
)

# basic expectations
test_that("Ebola model: basic expectations", {
  # runs without issues
  expect_no_condition(
    epidemic_ebola_r(
      population = pop,
      infection = ebola,
      time_end = 100
    )
  )

  set.seed(1)
  # returns a data.table
  data <- epidemic_ebola_r(
    population = pop,
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
  expect_setequal(
    unique(data$compartment),
    read_from_library(model_name = "ebola", what = "compartments")
  )

  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      data$value >= 0 & data$value <= sum(pop$demography_vector)
    )
  )

  # check for identical numbers of individuals at start and end
  # Note only valid for models without births and deaths
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1
  )

  # check that all age groups in the simulation are the same
  # size as the demography vector --- here, only one age group
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(pop$contact_matrix)
  )
  expect_identical(
    rowSums(final_state),
    pop$demography_vector,
    tolerance = 1
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
      preinfectious_period = 5,
      prop_hospitalised = 0.5,
      etu_safety = 0.8,
      funeral_safety = 0.5
    ),
    ebola_r0_high = infection(
      r0 = r0_high + 1.0, infectious_period = 5,
      preinfectious_period = 5,
      prop_hospitalised = 0.5,
      etu_safety = 0.8,
      funeral_safety = 0.5
    )
  )

  # get data
  data <- lapply(
    infection_list,
    function(infection_) {
      # run model on data
      data <- epidemic_ebola_r(
        population = pop,
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

#### Basic test of ebola model C++ version ####
# this model is not exported and is likely to be removed
test_that("Basic expectations for ebola model C++ version", {
  ebola_for_cpp <- infection(
    name = "ebolavirus disease",
    r0 = 1.7, infectious_period = 5,
    shape_E = 5, rate_E = 1, shape_I = 5, rate_I = 1
  )

  data_cpp <- epidemic_ebola_cpp(
    population = pop,
    infection = ebola_for_cpp
  )

  expect_s3_class(data_cpp, "data.table")
  expect_length(data_cpp, 4L)
  expect_named(
    data_cpp, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  # remove checks for Cpp version having same compartments as the library
})
