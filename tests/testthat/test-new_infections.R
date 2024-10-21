#### Tests for calculation of new infections ####
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 1, ncol = 5L,
  byrow = TRUE
)
colnames(initial_conditions) <- c(
  "susceptible", "exposed", "infectious", "recovered", "vaccinated"
)

# create a population where there are no contacts between groups
# this helps test the expectation that the final size proportions
# should be the same as the demographic group proportions
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1), # within group contacts only
  demography_vector = 67e6,
  initial_conditions = initial_conditions
)

# add a vaccination regime
vax_time_begin <- 50
vaccination <- vaccination(
  time_begin = matrix(vax_time_begin),
  time_end = matrix(200), nu = matrix(0.0005)
)

# run epidemic simulation with vaccination
time_end <- 100
increment <- 1
data <- model_default(
  population = uk_population,
  vaccination = vaccination,
  time_end = time_end,
  increment = increment
)

test_that("New infections: basic expectations and statistical correctness", {
  # expect run with no conditions
  expect_no_condition(
    new_infections(data, compartments_from_susceptible = "vaccinated")
  )
  data_ <- new_infections(data, compartments_from_susceptible = "vaccinated")

  # expect correct number of rows for two demography groups
  expect_identical(
    data_$time, seq(0, time_end, increment)
  )
  # expect columns with specific names
  expect_named(
    data_, c("time", "demography_group", "new_infections"),
    ignore.order = TRUE
  )

  # test for correctness
  # test that new infections are exactly the same as change in susceptibles
  # - change in vaccinations
  new_infections <- data_$new_infections

  # note use of original data.table
  delta_susc <- data[data$compartment == "susceptible", ][["value"]]
  delta_susc <- c(0, -diff(delta_susc))

  delta_vax <- data[data$compartment == "vaccinated", ][["value"]]
  delta_vax <- c(0, diff(delta_vax))

  expect_identical(
    new_infections,
    delta_susc - delta_vax
  )
})

test_that("New infections: statistical correctness without vaccination", {
  data_ <- new_infections(data)
  new_infections <- data_$new_infections

  # note use of original data.table
  delta_susc <- data[data$compartment == "susceptible", ][["value"]]
  delta_susc <- c(0, -diff(delta_susc))

  delta_vax <- data[data$compartment == "vaccinated", ][["value"]]
  delta_vax <- c(0, diff(delta_vax))

  # expect that new infections are exactly identical to change in susceptibles
  expect_identical(
    new_infections,
    delta_susc
  )
})

#### Tests for multiple demographic groups ####
# Prepare data with multiple demographic groups
initial_conditions <- matrix(
  c(0.9999, 0.0001, 0, 0, 0), # note "infectious" are zero
  nrow = 1, ncol = 5L,
  byrow = TRUE
)
initial_conditions <- rbind(initial_conditions, initial_conditions)
colnames(initial_conditions) <- c(
  "susceptible", "exposed", "infectious", "recovered", "vaccinated"
)

# create a dummy population
demography_vector <- 67e6 * c(0.4, 0.6)
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1, 2, 2),
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

# run epidemic simulation with vaccination
time_end <- 100L
data <- model_default(
  population = uk_population,
  time_end = time_end,
  increment = increment
)

test_that("New infections: aggregation over demographic group", {
  # expect that structure and contents for multiple groups is correct
  data_new_infections <- new_infections(data)
  expect_identical(
    nrow(data_new_infections),
    (time_end + 1L) * length(demography_vector)
  )
  # NOTE: use `expect_equal()` to avoid multiple integer-float conversions
  expect_equal(
    data_new_infections$time,
    rep(seq(0, time_end), each = length(demography_vector))
  )
  # snapshot test
  expect_snapshot(
    data_new_infections
  )

  # expect that setting `by_group = FALSE` is identical to manual aggregation
  data_new_infections_agg <- new_infections(data, by_group = FALSE)
  data_infections_manual_agg <- aggregate(
    new_infections ~ time,
    data_new_infections, 
    FUN = sum
  )
  expect_identical(
    data_new_infections_agg$new_infections,
    data_infections_manual_agg$new_infections
  )
})
