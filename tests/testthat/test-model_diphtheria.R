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
  expect_s3_class(output, "data.frame")
  expect_named(
    output, c("compartment", "demography_group", "value", "time"),
    ignore.order = TRUE
  )
  expect_identical(
    unique(output$compartment),
    c("susceptible", "exposed", "infectious", "hospitalised", "recovered")
  )
  # check for all positive values within the range 0 and total population size
  expect_true(
    all(
      output$value >= 0 & output$value <= sum(demography_vector)
    )
  )
  # check for identical numbers of individuals at start and end
  # NOTE: high tolerance because hospitalised compartment is not directly
  # linked to infectious compartment per Finger et al. model structure.
  # leads to more individuals at final state than initial state
  expect_identical(
    sum(output[output$time == min(output$time), ]$value),
    sum(output[output$time == max(output$time), ]$value),
    tolerance = 100
  )
  # check that all age groups in the simulation are the same
  # size as the demography vector
  final_state <- matrix(
    unlist(output[output$time == max(output$time), ]$value),
    nrow = n_age_groups
  )
  # NOTE: no checks for final state equal to demography vector as
  # vaccinated individuals are removed from model
})
