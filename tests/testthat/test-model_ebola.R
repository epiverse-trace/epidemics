#### Tests for the ebola model ####

# Prepare population and parameters
demography_vector <- 67000 # small population
contact_matrix <- matrix(1)

# manual case counts divided by pop size rather than proportions as small sizes
# introduce errors when converting to counts in the model code; extra
# individuals may appear
infectious <- 1
exposed <- 10
initial_conditions <- matrix(
  c(demography_vector - infectious - exposed, exposed, infectious, 0, 0, 0) /
    demography_vector,
  nrow = 1
)
rownames(contact_matrix) <- "full_pop"
pop <- population(
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

compartments <- c(
  "susceptible", "exposed", "infectious", "hospitalised", "funeral", "removed"
)

# prepare integer values for expectations on data length
time_end <- 100L
replicates <- 10L

test_that("Ebola model: basic expectations, scalar arguments", {
  # expect run with no conditions for default arguments
  expect_no_condition(model_ebola(pop, replicates = replicates))

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- withr::with_seed(1, model_ebola(pop, replicates = replicates))
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L) # extra column for replicate identifier
  expect_named(
    data, c("time", "demography_group", "compartment", "value", "replicate"),
    ignore.order = TRUE
  )
  expect_identical(
    nrow(data),
    length(demography_vector) * (time_end) * length(compartments) * replicates
  )
  expect_identical(unique(data$compartment), compartments)
  expect_true(
    checkmate::test_numeric(
      data$value,
      upper = max(demography_vector), lower = 0, any.missing = FALSE
    )
  )
  expect_identical(
    unique(data$demography_group), rownames(pop$contact_matrix)
  )

  # expect constant population size overall and per demography-group
  # NOTE: testing expectation for a single replicate
  data <- data[data$replicate == 1, ]
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(pop$contact_matrix)
  )
  expect_identical(
    rowSums(final_state), pop$demography_vector,
    tolerance = 1e-6
  )

  # expect snapshot is preserved; note this is a single replicate
  expect_snapshot(
    head(data, 20L)
  )
})

# NOTE: statistical correctness is not expected to change for vectorised input
test_that("Ebola model: statistical correctness, parameters", {
  # expect final size increases with transmission_rate
  size_beta_low <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, transmission_rate = 1.3 / 12, replicates = 1)
    )
  )
  size_beta_high <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, transmission_rate = 1.5 / 12, replicates = 1)
    )
  )
  expect_gt(size_beta_high, size_beta_low)

  # expect final size increases with infectiousness rate (lower incubation time)
  size_sigma_low <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, infectiousness_rate = 1 / 5, replicates = 1)
    )
  )
  size_sigma_high <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, infectiousness_rate = 1 / 2, replicates = 1)
    )
  )
  expect_gt(size_sigma_high, size_sigma_low)

  # expect final size decreases with removal rate (fewer infection events)
  size_gamma_low <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, removal_rate = 1 / 5, replicates = 1)
    )
  )
  size_gamma_high <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, removal_rate = 1 / 2, replicates = 1)
    )
  )
  expect_lt(size_gamma_high, size_gamma_low) # NOTE: expect less than!

  # expect that increased ETU safety reduces final size
  size_etu_risk_high <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, etu_risk = 0.9, replicates = 1)
    )
  )
  size_etu_risk_low <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, etu_risk = 0.1, replicates = 1)
    )
  )
  expect_gt(size_etu_risk_high, size_etu_risk_low)

  # expect that increased hospitalisation safety reduces final size
  size_prop_comm_high <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, prop_community = 0.9, etu_risk = 0.1, replicates = 1)
    )
  )
  size_prop_comm_low <- epidemic_size(
    withr::with_seed(
      1, model_ebola(pop, prop_community = 0.3, etu_risk = 0.3, replicates = 1)
    )
  )
  expect_gt(size_prop_comm_high, size_prop_comm_low)

  # NOTE: not testing the effect of `erlang_subcompartments` as this is
  # equivalent to testing the effect of `transmission_rate` or `removal_rate`
  # with less reasonable variation (lower bound 1, only integer values)

  # Further tests to check model mechanics: expectations on final size when
  # specific parameters are zero
  popsize <- 10e3
  total_cases <- 20

  population <- population(
    contact_matrix = matrix(1),
    demography_vector = popsize,
    initial_conditions = matrix(
      c(popsize - total_cases, 0, total_cases, 0, 0, 0) / popsize,
      nrow = 1L, ncol = 6L, byrow = TRUE
    )
  )

  # Expect hospitalisations are zero when all infections are in the community
  data <- model_ebola(
    population = population, prop_community = 1.0, replicates = 1
  )
  expect_identical(
    unique(data[data$compartment == "hospitalised", ]$value), 0
  )

  # Expect that full ETU safety leads to a fixed final size == total_cases
  # create a dummy population where all individuals are hospitalised
  popsize <- 10e3
  total_cases <- 20

  population <- population(
    contact_matrix = matrix(1),
    demography_vector = popsize,
    initial_conditions = matrix(
      c(popsize - total_cases, 0, 0, total_cases, 0, 0) / popsize,
      nrow = 1L, ncol = 6L, byrow = TRUE
    )
  )
  data <- model_ebola(
    population = population, prop_community = 0, etu_risk = 0, replicates = 1
  )
  expect_equal(
    epidemic_size(data), total_cases,
    ignore_attr = TRUE
  )

  # Expect that full funeral safety leads to a fixed final size == total_cases
  # create a population with only funeral infections possible
  population <- population(
    contact_matrix = matrix(1),
    demography_vector = popsize,
    initial_conditions = matrix(
      c(popsize - total_cases, 0, 0, 0, total_cases, 0) / popsize,
      nrow = 1L, ncol = 6L, byrow = TRUE
    )
  )
  data <- model_ebola(
    population = population, prop_community = 1,
    funeral_risk = 0, replicates = 1
  )
  expect_equal(
    epidemic_size(data), total_cases,
    ignore_attr = TRUE
  )
})

# prepare baseline for comparison of against intervention scenarios
data_baseline <- withr::with_seed(1, model_ebola(pop, replicates = 1))

test_that("Ebola model: rate interventions", {
  intervention_01 <- intervention(
    "mask_mandate", "rate", 0, time_end, 0.5
  )
  intervention_02 <- intervention(
    "mask_mandate", "rate", time_end / 2, time_end, 0.1
  )
  intervention <- c(intervention_01, intervention_02)

  # repeat some basic checks from default case with no intervention
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_ebola(
      pop,
      intervention = list(transmission_rate = intervention)
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- withr::with_seed(
    1,
    model_ebola(
      pop,
      intervention = list(transmission_rate = intervention),
      replicates = 1
    )
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L) # replicates column expected

  # expect final size is lower with intervention
  expect_gt(epidemic_size(data_baseline), epidemic_size(data))
})

test_that("Ebola model: time dependence", {
  # expect time dependence is correctly handled
  time_dependence <- list(
    transmission_rate = function(time, x, t_change = time_end / 2) {
      ifelse(time > t_change, x / 2, x)
    },
    # hospitalisation increases, and prop_community decreases by 50%
    prop_community = function(time, x, t_change = time_end / 2) {
      ifelse(time > t_change, x * 0.5, x)
    }
  )

  # repeat some basic checks from default case with no time_dependence
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_ebola(pop, time_dependence = time_dependence)
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- withr::with_seed(
    1,
    model_ebola(pop, time_dependence = time_dependence, replicates = 1)
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 5L)

  # expect final size is lower with intervention
  expect_gt(epidemic_size(data_baseline), epidemic_size(data))
})

test_that("Ebola model: errors and warnings, scalar arguments", {
  # expect errors on basic input checking
  expect_error(
    model_ebola(population = "population"),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  expect_error(
    model_ebola(population = population),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  pop_wrong_compartments <- pop
  pop_wrong_compartments$initial_conditions <- initial_conditions[, -1]
  expect_error(
    model_ebola(pop_wrong_compartments),
    regexp = "(Assertion on)*(initial_conditions)*failed"
  )

  # expect errors for infection parameters
  expect_error(
    model_ebola(pop, transmission_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_ebola(pop, infectiousness_rate = list(0.2)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_ebola(pop, removal_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_ebola(pop, erlang_subcompartments = "2"),
    regexp = "Must be of type 'integerish"
  )
  expect_error(
    model_ebola(pop, erlang_subcompartments = 2.5),
    regexp = "Must be of type 'integerish"
  )
  expect_error(
    model_ebola(pop, prop_community = 1.01),
    regexp = "not <= 1"
  )
  expect_error(
    model_ebola(pop, funeral_risk = "1.01"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_ebola(pop, etu_risk = 1.01),
    regexp = "not <= 1"
  )
  # NOTE: further checks on parameter type are ommitted here.

  # expect error on time parameters
  expect_error(
    model_ebola(pop, time_end = "100"),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_ebola(pop, time_end = 100.5),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_ebola(pop, time_end = c(100, -100, 10)),
    regexp = "(Element)*(is not >= 0)"
  )
  expect_error(
    model_ebola(pop, replicates = "0.1"),
    regexp = "Must be of type 'count'"
  )
  expect_error(
    model_ebola(pop, replicates = c(10, 20)),
    regexp = "Must have length 1"
  )

  # expect error on poorly specified interventions
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, 0.5 # no contacts_interventions
  )
  expect_error(
    model_ebola(
      pop,
      intervention = list(contacts = intervention)
    )
  )

  # expect error on poorly specified time-dependence function list
  expect_error(
    model_ebola(
      pop,
      time_dependence = function(x) x
    ),
    regexp = "Must be of type 'list'"
  )
  expect_error(
    model_ebola(
      pop,
      time_dependence = list(function(x) x)
    ),
    regexp = "Must have names"
  )
  expect_error(
    model_ebola(
      pop,
      time_dependence = list(transmission_rate = function(x) x)
    ),
    regexp = "Must have first formal arguments \\(ordered\\): time,x."
  )
  expect_error(
    model_ebola(
      pop,
      time_dependence = list(transmission_rate = NULL)
    ),
    regexp = "Contains missing values"
  )
})

# prepare vectors of parameters
beta <- rnorm(10, 1.3 / 7, sd = 0.01)
sigma <- rnorm(10, 0.5, sd = 0.01)
gamma <- rnorm(10, 1 / 7, sd = 0.01)

test_that("Ebola model: infection parameters as vectors", {
  # expect no conditions when vectors are passed
  expect_no_condition(
    model_ebola(
      pop,
      transmission_rate = beta, infectiousness_rate = sigma,
      removal_rate = gamma, replicates = 1
    )
  )
  # expect output structure is a nested data.table
  output <- model_ebola(
    pop,
    transmission_rate = beta, infectiousness_rate = sigma,
    removal_rate = gamma, replicates = 1
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(beta))
  expect_identical(output$transmission_rate, beta)
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(output$param_set, seq_along(beta))
  expect_identical(unique(output$scenario), 1L)

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null")
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null"), any.missing = FALSE
  )
})

test_that("Ebola model: composable elements as lists", {
  # expect no conditions when multiple interventions are passed
  # NOTE: only rate interventions are allowed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      )
    ),
    scenario_02 = list(
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      ),
      etu_risk = intervention(
        "better_isolation", "rate", time_end / 2, time_end, 0.5
      )
    )
  )

  expect_no_condition(
    model_ebola(pop, intervention = npi_list, replicates = 10)
  )

  # expect output is a nested data.frame-like object
  output <- model_ebola(pop, intervention = npi_list, replicates = 10)
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(npi_list))
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(output$scenario, seq_along(npi_list))
  expect_identical(unique(output$param_set), 1L)

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  # some interventions may be missing
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null")
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null")
  )
})

test_that("Ebola model: multi-parameter, multi-composables", {
  # expect no conditions when multiple interventions or vaccinations are passed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      )
    ),
    scenario_02 = list(
      transmission_rate = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      ),
      etu_risk = intervention(
        "better_isolation", "rate", time_end / 2, time_end, 0.5
      )
    )
  )

  # reuse parameter sets from earlier tests
  expect_no_condition(
    model_ebola(
      pop,
      transmission_rate = beta, removal_rate = gamma,
      intervention = npi_list,
      replicates = 10
    )
  )

  # expect output is a nested data.frame-like object
  output <- model_ebola(
    pop,
    transmission_rate = beta, removal_rate = gamma,
    intervention = npi_list, replicates = 10
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(npi_list) * length(beta))
  checkmate::expect_list(output$data, types = "data.frame", any.missing = FALSE)

  # expect `parameter_set` and `scenario` are correctly filled
  expect_identical(
    output$scenario, rep(seq_along(npi_list), length(beta))
  )
  expect_identical(unique(output$param_set), seq_along(beta))
  expect_identical(
    output$param_set, rep(seq_along(beta), each = length(npi_list))
  )

  # expect list column of interventions and vaccination
  checkmate::expect_list(
    output$population,
    types = "population", any.missing = FALSE
  )
  # some interventions or vaccinations may be missing
  checkmate::expect_list(
    output$intervention,
    types = c("list", "null"), any.missing = TRUE
  )
  checkmate::expect_list(
    output$time_dependence,
    types = c("list", "null"), any.missing = FALSE
  )
})

test_that("Ebola model: errors on vectorised input", {
  # expect errors on poorly specified vector inputs
  expect_error(
    model_ebola(
      pop,
      transmission_rate = beta[-1], removal_rate = gamma
    ),
    regexp = "All parameters must be of the same length, or must have length 1"
  )
  expect_error(
    model_ebola(
      pop,
      intervention = list(
        NULL,
        list(dummy = intervention)
      )
    ),
    regexp =
      "`intervention` must be a list of <rate_intervention> or a list of such"
  )

  # expect time-dependence cannot be vectorised
  expect_error(
    model_ebola(
      pop,
      time_dependence = list(
        time_dep_01 = list(
          transmission_rate = function(x) x
        ),
        time_dep_02 = list(
          transmission_rate = function(x) x
        )
      )
    ),
    regexp = "May only contain the following types: \\{function\\}"
  )
})
