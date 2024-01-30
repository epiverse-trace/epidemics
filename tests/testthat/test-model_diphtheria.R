#### Tests for the diphtheria model ####

# create a dummy camp population with three age groups
# diphtheria model is SEIHR
# assume that most are susceptible, some infectious
# values taken from supplementary material in Finger et al. for the
# Kutupalong camp, rounded to the nearest 100
n_age_groups <- 3
demography_vector <- c(83000, 108200, 224600)
initial_conditions <- matrix(0, nrow = n_age_groups, ncol = 5)

# set susceptibles and infectious
initial_conditions[, 1] <- demography_vector - 1
initial_conditions[, 3] <- rep(1, n_age_groups)

camp_pop <- population(
  contact_matrix = matrix(1, nrow = n_age_groups, ncol = n_age_groups),
  demography_vector = demography_vector,
  initial_conditions = initial_conditions / demography_vector
)

# Basic expectations for the diphtheria model
test_that("Diptheria model, basic expectations", {
  # model runs with default arguments
  expect_no_condition(
    model_diphtheria_cpp(
      population = camp_pop
    )
  )

  # expectations on output
  prop_vaccinated <- c(0.2, 0.10, 0.1) # vaccinated not included in model
  output <- model_diphtheria_cpp(
    population = camp_pop,
    prop_vaccinated = prop_vaccinated
  )
  expect_s3_class(output, "epidemic")

  # check for snapshot
  # NOTE: this snapshot shows fractional individuals not equal to
  # initial conditions - this is because initial susceptibles are reduced
  # by prop_vaccinated
  expect_snapshot(
    output
  )

  data <- get_parameter(output, "data")
  expect_named(
    data, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  expect_identical(
    get_parameter(output, "compartments"),
    c("susceptible", "exposed", "infectious", "hospitalised", "recovered")
  )
  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      data$value >= 0 & data$value <= sum(demography_vector)
    )
  )
  # check for identical numbers of individuals at start and end
  # NOTE: high tolerance because hospitalised compartment is not directly
  # linked to infectious compartment per Finger et al. model structure.
  # leads to different individuals at final state than initial state
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )
  # check that all age groups in the simulation are the same
  # size as the demography vector
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = n_age_groups
  )
  # NOTE: no checks for final state equal to demography vector as
  # vaccinated individuals are removed from model
})

# Expectations for the diphtheria model with changed population sizes
test_that("Diphtheria model with population size changes", {
  # model runs with population change functionality
  # create population change data
  p <- list(
    time = 70,
    values = list(
      c(1e4, 2e5, 1e5)
    )
  )

  expect_no_condition(
    model_diphtheria_cpp(
      population = camp_pop,
      population_change = p
    )
  )

  # NOTE: expected final population size is larger than the initial
  # but identical to the original + added population
  output <- model_diphtheria_cpp(
    population = camp_pop,
    prop_hosp = 0.08,
    population_change = p,
    time_end = 200
  )
  data <- get_parameter(output, "data")

  last_value <- aggregate(
    value ~ demography_group,
    data = data[data$time == max(data$time), ], FUN = "sum"
  )
  last_value

  expect_identical(
    last_value$value,
    camp_pop$demography_vector + p$values[[1]],
    tolerance = 1e-6
  )

  ## Multiple population changes including decreases
  p <- list(
    time = c(70, 80),
    values = list(
      c(1e4, 2e5, 1e5), # influx to camp
      c(-9e3, -1.5e5, -0.5e5) # evacuation from camp
    )
  )

  expect_no_condition(
    model_diphtheria_cpp(
      population = camp_pop,
      population_change = p
    )
  )

  # NOTE: expected final population size is larger than the initial
  # but identical to the original + net added population
  output <- model_diphtheria_cpp(
    population = camp_pop,
    prop_hosp = 0.08,
    population_change = p,
    time_end = 200
  )
  data <- get_parameter(output, "data")

  last_value <- aggregate(
    value ~ demography_group,
    data = data[data$time == max(data$time), ], FUN = "sum"
  )
  last_value

  expect_identical(
    last_value$value,
    camp_pop$demography_vector + Reduce(x = p$values, f = `+`),
    tolerance = 1e-6
  )
})

# Test for poor specification of the population change mechanic
test_that("Diphtheria model handles population_change errors", {
  # badly specified pop change
  p <- list(
    time = "some time",
    values = list(
      c(1, 2, 3)
    )
  )
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      population_change = p
    ),
    regexp = "May only contain the following types: \\{numeric,list\\}"
  )

  # wrong name of first element
  p <- list(
    timesteps = 10,
    values = list(
      c(1, 2, 3)
    )
  )
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      population_change = p
    ),
    regexp = "(Names must be identical)*(time)*(values)"
  )

  # wrong length of values - must be same length as demography groups
  p <- list(
    time = 10,
    values = list(
      c(1, 2)
    )
  )
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      population_change = p
    ),
    regexp = "`population_change` `values` must be same length as demography"
  )
})

# Tests for expected failures/input checks on interventions and time-dependence
test_that("Diphtheria model input checks", {
  # expect no contacts interventions allowed
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      intervention = list(
        contacts = no_contacts_intervention(camp_pop) # a dummy intervention
      )
    ),
    regexp = "May only contain the following types: \\{rate_intervention\\}"
  )

  # expect failure on rate interventions targeting disallowed parameters
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      intervention = list(
        beta = intervention(
          type = "rate",
          time_begin = 10, time_end = 20,
          reduction = 0.1
        )
      )
    ),
    regexp = "(must be a subset of)*(but has additional elements)"
  )

  # expect failure on time dependence targeting disallowed parameters
  expect_error(
    model_diphtheria_cpp(
      camp_pop,
      time_dependence = list(
        beta = function(time, x) x
      )
    ),
    regexp = "(must be a subset of)*(but has additional elements)"
  )
})
