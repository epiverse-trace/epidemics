#### Tests for the diphtheria model ####

# create a dummy camp population with three age groups
# diphtheria model is SEIHR
# assume that most are susceptible, some infectious
# values taken from supplementary material in Finger et al. for the
# Kutupalong camp, rounded to the nearest 100
n_age_groups <- 3
demography_vector <- c(83000, 108200, 224600)
# a dummy contact matrix
contact_matrix <- matrix(1, nrow = n_age_groups, ncol = n_age_groups)
rownames(contact_matrix) <- c("0-4", "5-14", "15+") # only for testing
initial_conditions <- matrix(0, nrow = n_age_groups, ncol = 5)

# set susceptibles and infectious
initial_conditions[, 1] <- demography_vector - 1
initial_conditions[, 3] <- rep(1, n_age_groups)

camp_population <- population(
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions / demography_vector
)

# model run time
time_end <- 100L
compartments <- c(
  "susceptible", "exposed", "infectious", "hospitalised", "recovered"
)

test_that("Diptheria model: basic expectations, scalar arguments", {
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_diphtheria_cpp(camp_population)
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_diphtheria_cpp(camp_population, time_end = time_end)
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)
  expect_named(
    data, c("time", "demography_group", "compartment", "value"),
    ignore.order = TRUE
  )
  expect_identical(
    nrow(data),
    length(demography_vector) * (time_end + 1L) * length(compartments)
  )
  expect_identical(unique(data$compartment), compartments)
  expect_true(
    checkmate::test_numeric(
      data$value,
      upper = max(demography_vector), lower = 0, any.missing = FALSE
    )
  )
  # NOTE: diphtheria model contact matrix is a dummy
  expect_identical(
    unique(data$demography_group), rownames(contact_matrix)
  )

  # expect (very small, non zero) individuals in hospitalised compartments
  checkmate::expect_numeric(
    tail(data[data$compartment == "hospitalised", ]$value),
    lower = 0
  )

  # expect constant population size overall and per demography-group
  expect_identical(
    sum(data[data$time == min(data$time), ]$value),
    sum(data[data$time == max(data$time), ]$value),
    tolerance = 1e-6
  )
  final_state <- matrix(
    unlist(data[data$time == max(data$time), ]$value),
    nrow = nrow(contact_matrix)
  )
  expect_identical(
    rowSums(final_state), camp_population$demography_vector,
    tolerance = 1e-6
  )
})

# NOTE: statistical correctness is not expected to change for vectorised input
test_that("Diptheria model: statistical correctness, parameters", {
  # expect final size increases with transmissibility
  size_beta_low <- epidemic_size(
    model_diphtheria_cpp(camp_population, transmissibility = 4 / 4.5)
  )
  size_beta_high <- epidemic_size(
    model_diphtheria_cpp(camp_population, transmissibility = 4.5 / 4.5)
  )
  expect_true(
    all(size_beta_high > size_beta_low)
  )

  # expect final size increases with infectiousness rate (lower incubation time)
  size_sigma_low <- epidemic_size(
    model_diphtheria_cpp(camp_population, infectiousness_rate = 1 / 5)
  )
  size_sigma_high <- epidemic_size(
    model_diphtheria_cpp(camp_population, infectiousness_rate = 1 / 2)
  )
  expect_true(
    all(size_sigma_high > size_sigma_low)
  )

  # expect final size increases with initial infections
  initial_conditions_high <- initial_conditions
  initial_conditions_high[, 1] <- demography_vector - 10
  initial_conditions_high[, 3] <- rep(10, n_age_groups)
  camp_population_high_infections <- population(
    name = "UK population",
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions_high / demography_vector
  )
  size_infections_low <- epidemic_size(model_diphtheria_cpp(camp_population))
  size_infections_high <- epidemic_size(
    model_diphtheria_cpp(camp_population_high_infections)
  )
  # NOTE: difference is very small, as >90% individuals are infected
  # even in the low initial infections scenario
  expect_true(
    all(size_infections_high > size_infections_low)
  )

  # expect no hospitalisations when hospitalisation rate = 0
  data <- model_diphtheria_cpp(camp_population, hosp_entry_rate = 0.0)
  expect_identical(
    unique(data[grepl("hospitalised", data$compartment, fixed = TRUE)]$value),
    0
  )

  # expect no hospitalisations when infection reporting rate = 0
  data <- model_diphtheria_cpp(camp_population, reporting_rate = 0.0)
  expect_identical(
    unique(data[grepl("hospitalised", data$compartment, fixed = TRUE)]$value),
    0
  )

  # expect no hospitalisations when proportion needing/given hospitalisation = 0
  data <- model_diphtheria_cpp(camp_population, prop_hosp = 0.0)
  expect_identical(
    unique(data[grepl("hospitalised", data$compartment, fixed = TRUE)]$value),
    0
  )
})

# prepare baseline for comparison of against intervention scenarios
data_baseline <- model_diphtheria_cpp(camp_population)

test_that("Diphtheria model: pre-exisiting immunity and stats. correctness", {
  # expect no condition on pre-existing immunity user-supplied value
  expect_no_condition(
    model_diphtheria_cpp(camp_population, prop_vaccinated = c(0.5, 0.5, 0.5))
  )
  # exepct that pre-existing immunity reduces final size
  data <- model_diphtheria_cpp(
    camp_population,
    prop_vaccinated = c(0.5, 0.5, 0.5)
  )
  size_baseline <- epidemic_size(data_baseline)
  size_pre_immunity <- epidemic_size(data)
  expect_true(
    all(size_baseline > size_pre_immunity)
  )
})

test_that("Diptheria model: rate interventions and stats. correctness", {
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
    model_diphtheria_cpp(
      camp_population,
      intervention = list(transmissibility = intervention)
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_diphtheria_cpp(
    camp_population,
    intervention = list(transmissibility = intervention)
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)

  # expect final size is lower with intervention
  expect_true(
    all(epidemic_size(data_baseline) > epidemic_size(data))
  )

  # expect intervention on proportion needing hospital reduces hospitalised
  intervention <- intervention(
    "early_detection", "rate", time_end / 2, time_end,
    reduction = 0.5
  )
  data_high_hosp <- model_diphtheria_cpp(
    camp_population,
    prop_hosp = 0.7,
    reporting_rate = 1.0
  )
  data_hosp_intervention <- model_diphtheria_cpp(
    camp_population,
    prop_hosp = 0.7,
    reporting_rate = 1.0,
    intervention = list(
      prop_hosp = intervention
    )
  )
  expect_true(
    all(
      data_high_hosp[data_high_hosp$compartment ==
        "hospitalised"]$value >=
        data_hosp_intervention[data_hosp_intervention$compartment ==
          "hospitalised"]$value
    )
  )
})

test_that("Diptheria model: time dependence", {
  # expect time dependence is correctly handled
  time_dependence <- list(
    transmissibility = function(time, x, t_change = time_end / 2) {
      ifelse(time > t_change, x / 2, x)
    },
    recovery_rate = function(time, x, t_change = time_end / 2) {
      ifelse(time > t_change, x + x / 2, x)
    }
  )

  # repeat some basic checks from default case with no time_dependence
  # expect run with no conditions for default arguments
  expect_no_condition(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = time_dependence
    )
  )

  # expect data.frame-inheriting output with 4 cols; C++ model time begins at 0
  data <- model_diphtheria_cpp(
    camp_population,
    time_dependence = time_dependence
  )
  expect_s3_class(data, "data.frame")
  expect_identical(length(data), 4L)

  # expect final size is lower with intervention
  expect_true(
    all(epidemic_size(data_baseline) > epidemic_size(data))
  )
})

# Expectations for the diphtheria model with changed population sizes
test_that("Diphtheria model: population size changes", {
  # test population change object
  p <- list(
    time = 70,
    values = list(
      c(1e4, 2e5, 1e5)
    )
  )

  # expect no conditions
  expect_no_condition(
    model_diphtheria_cpp(
      population = camp_population,
      population_change = p
    )
  )

  # NOTE: expected final population size is larger than the initial
  # but identical to the original + added population
  data <- model_diphtheria_cpp(
    population = camp_population,
    prop_hosp = 0.08,
    population_change = p,
    time_end = 200
  )

  last_value <- aggregate(
    value ~ demography_group,
    data = data[data$time == max(data$time), ], FUN = "sum"
  )
  expect_identical(
    sort(last_value$value), # as aggregate() re-orders DF by age group labels
    camp_population$demography_vector + p$values[[1]],
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
  # expect no conditions when running
  expect_no_condition(
    model_diphtheria_cpp(
      population = camp_population,
      population_change = p
    )
  )

  # NOTE: expected final population size is larger than the initial
  # but identical to the original + net added population
  data <- model_diphtheria_cpp(
    population = camp_population,
    prop_hosp = 0.08,
    population_change = p,
    time_end = 200
  )

  last_value <- aggregate(
    value ~ demography_group,
    data = data[data$time == max(data$time), ], FUN = "sum"
  )

  expect_identical(
    sort(last_value$value), # as aggregate() re-orders DF by age group labels
    camp_population$demography_vector + Reduce(x = p$values, f = `+`),
    tolerance = 1e-6
  )
})

test_that("Diptheria model: errors and warnings, scalar arguments", {
  # expect errors on basic input checking
  expect_error(
    model_diphtheria_cpp(population = "population"),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  expect_error(
    model_diphtheria_cpp(population = population),
    regexp = "(Assertion on 'population' failed)*(Must inherit)*(population)"
  )
  pop_wrong_compartments <- camp_population
  pop_wrong_compartments$initial_conditions <- initial_conditions[, -1]
  expect_error(
    model_diphtheria_cpp(pop_wrong_compartments),
    regexp = "(Assertion on)*(initial_conditions)*failed"
  )

  # expect errors for infection parameters
  expect_error(
    model_diphtheria_cpp(camp_population, transmissibility = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, infectiousness_rate = list(0.2)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, recovery_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, prop_hosp = list(0.002)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, reporting_rate = "0.19"),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, hosp_entry_rate = list(0.002)),
    regexp = "Must be of type 'numeric'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, hosp_exit_rate = "0.002"),
    regexp = "Must be of type 'numeric'"
  )

  # expect error on time parameters
  expect_error(
    model_diphtheria_cpp(camp_population, time_end = "100"),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, time_end = 100.5),
    regexp = "Must be of type 'integerish'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, time_end = c(100, -100, 10)),
    regexp = "(Element)*(is not >= 0)"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, increment = "0.1"),
    regexp = "Must be of type 'number'"
  )
  expect_error(
    model_diphtheria_cpp(camp_population, increment = c(0.1, 0.2)),
    regexp = "Must have length 1"
  )

  # expect error on contacts intervention as they are not allowed
  intervention <- intervention(
    "school_closure", "contacts", 0, time_end, c(0.5, 0.2, 0.1)
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      intervention = list(contacts = intervention)
    ),
    regexp = "has additional elements \\{'contacts'\\}"
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      intervention = list(transmissibility = intervention)
    ),
    regexp = "Must inherit from class 'rate_intervention'"
  )

  # expect error if vaccination is passed
  vax_single_dose <- vaccination(
    nu = matrix(1e-3, nrow = 2),
    time_begin = matrix(00, nrow = 2),
    time_end = matrix(100, nrow = 2)
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      vaccination = vax_single_dose
    ),
    regexp = "unused argument"
  )

  # expect error on poorly specified time-dependence function list
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = function(x) x
    ),
    regexp = "Must be of type 'list'"
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = list(function(x) x)
    ),
    regexp = "Must have names"
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = list(transmissibility = function(x) x)
    ),
    regexp = "Must have first formal arguments \\(ordered\\): time,x."
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = list(transmissibility = NULL)
    ),
    regexp = "Contains missing values"
  )

  # expect error on prop_vaccinated


  # exepect error on badly specified pop change
  p <- list(
    time = "some time",
    values = list(
      c(1, 2, 3)
    )
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
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
      camp_population,
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
      camp_population,
      population_change = p
    ),
    regexp = "`population_change` `values` must be same length as demography"
  )
})

# prepare vectors of parameters
beta <- rnorm(10, 1.3 / 7, sd = 0.01)
sigma <- rnorm(10, 0.5, sd = 0.01)
gamma <- rnorm(10, 1 / 7, sd = 0.01)

test_that("Diptheria model: infection parameters as vectors", {
  # expect no conditions when vectors are passed
  expect_no_condition(
    model_diphtheria_cpp(
      camp_population,
      transmissibility = beta, infectiousness_rate = sigma,
      recovery_rate = gamma
    )
  )
  # expect output structure is a nested data.table
  output <- model_diphtheria_cpp(
    camp_population,
    transmissibility = beta, infectiousness_rate = sigma,
    recovery_rate = gamma
  )
  expect_s3_class(output, c("data.frame", "data.table"))
  expect_identical(nrow(output), length(beta))
  expect_identical(output$transmissibility, beta)
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

test_that("Diptheria model: composable elements as lists", {
  # expect no conditions when multiple interventions or vaccinations are passed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      transmissibility = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      ),
      prop_hosp = intervention(
        "better_care", "rate", time_end / 2, time_end, 0.5
      )
    )
  )

  # artificially high hospitalisation rate
  expect_no_condition(
    model_diphtheria_cpp(
      camp_population,
      intervention = npi_list, prop_hosp = 1.0
    )
  )

  # expect output is a nested data.frame-like object
  output <- model_diphtheria_cpp(camp_population, intervention = npi_list)
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

test_that("Diptheria model: multi-parameter, multi-composables", {
  # expect no conditions when multiple interventions or vaccinations are passed
  npi_list <- list(
    scenario_baseline = NULL,
    scenario_01 = list(
      transmissibility = intervention(
        "mask_mandate", "rate", 0, time_end, 0.5
      )
    )
  )
  # reuse parameter sets from earlier tests
  expect_no_condition(
    model_diphtheria_cpp(
      camp_population,
      transmissibility = beta, recovery_rate = gamma,
      intervention = npi_list
    )
  )

  # expect output is a nested data.frame-like object
  output <- model_diphtheria_cpp(
    camp_population,
    transmissibility = beta, recovery_rate = gamma,
    intervention = npi_list
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

test_that("Vacamole: errors on vectorised input", {
  # expect errors on poorly specified vector inputs
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      transmissibility = beta[-1], recovery_rate = gamma
    ),
    regexp = "All parameters must be of the same length, or must have length 1"
  )
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      intervention = list(
        NULL,
        list(dummy = intervention)
      )
    ),
    regexp =
      "`intervention` must be a list of <intervention>s or a list of such lists"
  )

  # expect time-dependence cannot be vectorised
  expect_error(
    model_diphtheria_cpp(
      camp_population,
      time_dependence = list(
        time_dep_01 = list(transmissibility = function(x) x),
        time_dep_02 = list(transmissibility = function(x) x)
      )
    ),
    regexp = "(May only contain the following types:)*(function)"
  )
})
