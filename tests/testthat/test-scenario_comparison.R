#### Tests for scenario comparison features ####

# set up models
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)

# prepare contact matrix
contact_matrix <- t(contact_data$matrix)

# prepare the demography vector
demography_vector <- contact_data$demography$population
names(demography_vector) <- rownames(contact_matrix)

# initial conditions
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

# assign rownames for clarity
rownames(initial_conditions) <- rownames(contact_matrix)

# create population object
uk_population <- population(
  name = "UK",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

# create vector of parameters
beta <- withr::with_seed(
  1,
  rnorm(100, mean = 1.3 / 7, sd = 0.005)
)

# run the baseline model
baseline <- model_default(
  population = uk_population,
  transmission_rate = beta
)

max_time <- 100
# prepare durations as starting at 25% of the way through an epidemic
# and ending halfway through
time_begin <- max_time / 4
time_end <- max_time / 2

# create three distinct contact interventions
# prepare an intervention that models school closures for 180 days
close_schools <- intervention(
  name = "School closure",
  type = "contacts",
  time_begin = time_begin,
  time_end = time_end,
  reduction = matrix(c(0.3, 0.01, 0.01))
)

# prepare an intervention which mostly affects adults 20 -- 65
close_workplaces <- intervention(
  name = "Workplace closure",
  type = "contacts",
  time_begin = time_begin,
  time_end = time_end,
  reduction = matrix(c(0.01, 0.3, 0.01))
)

# define intervention sets
intervention_sets <- list(
  list(
    contacts = close_schools
  ),
  list(
    contacts = close_workplaces
  )
)

# run comparator scenarios
scenarios <- model_default(
  population = uk_population,
  transmission_rate = beta,
  intervention = intervention_sets
)

test_that("`outcomes_averted()`: Basic expectations", {
  # expect runs without conditions
  expect_no_condition(
    outcomes_averted(
      baseline = baseline, scenarios = scenarios
    )
  )
  expect_no_condition(
    outcomes_averted(
      baseline = baseline, scenarios = scenarios, summarise = FALSE
    )
  )
  expect_no_condition(
    outcomes_averted(
      baseline = baseline, scenarios = scenarios, by_group = FALSE
    )
  )
  expect_no_condition(
    outcomes_averted(
      baseline = baseline, scenarios = scenarios, by_group = FALSE,
      summarise = FALSE
    )
  )

  # expect return type, structure, and basic statistical correctness
  # for default options
  averted_default <- outcomes_averted(
    baseline = baseline, scenarios = scenarios, by_group = TRUE
  )
  expect_s3_class(averted_default, "data.table")
  expect_named(
    averted_default,
    c(
      "scenario", "demography_group",
      sprintf("averted_%s", c("median", "lower", "upper"))
    )
  )
  expect_true(
    all(averted_default$averted_median > averted_default$averted_lower &
      averted_default$averted_median < averted_default$averted_upper)
  )

  # expectations for aggregation over full population
  averted_aggregate <- outcomes_averted(
    baseline = baseline, scenarios = scenarios, by_group = FALSE
  )
  expect_s3_class(averted_aggregate, "data.table")
  expect_named(
    averted_aggregate,
    c("scenario", sprintf("averted_%s", c("median", "lower", "upper")))
  )
  expect_true(
    all(averted_default$averted_median > averted_default$averted_lower &
      averted_default$averted_median < averted_default$averted_upper)
  )

  # expectations when aggregation is not applied
  averted_raw <- outcomes_averted(
    baseline = baseline, scenarios = scenarios, summarise = FALSE
  )
  expect_s3_class(averted_aggregate, "data.table")
  expect_named(
    averted_raw,
    c("scenario", "demography_group", "param_set", "outcomes_averted"),
    ignore.order = TRUE
  )

  # expectations when aggregation is not applied
  averted_agg_raw <- outcomes_averted(
    baseline = baseline, scenarios = scenarios,
    summarise = FALSE, by_group = FALSE
  )
  expect_s3_class(averted_agg_raw, "data.table")
  expect_named(
    averted_agg_raw,
    c("scenario", "param_set", "outcomes_averted"),
    ignore.order = TRUE
  )
})

test_that("`outcomes_averted()`: Errors and messages", {
  # expect error on bad baseline
  expect_error(
    outcomes_averted(baseline = "baseline", scenarios),
    regexp = "Must be a data.table"
  )
  expect_error(
    outcomes_averted(baseline[, -1], scenarios),
    regexp = "(Must have)*(12 cols)"
  )
  baseline_bad_names <- baseline[, -1]
  baseline_bad_names$dummy <- NA_real_
  expect_error(
    outcomes_averted(baseline_bad_names, scenarios),
    regexp = "must be a nested <data.table> with expected column names"
  )

  # expect error on bad scenarios
  expect_error(
    outcomes_averted(baseline, split(scenarios, by = "scenario")),
    regexp = "Must be a data.table"
  )
  expect_error(
    outcomes_averted(baseline, scenarios[, -1]),
    regexp = "(Must have)*(12 cols)"
  )
  scenarios$extra_col <- NA_real_
  expect_error(
    outcomes_averted(baseline, scenarios),
    regexp = "Names must be a identical to set"
  )
  scenarios$extra_col <- NULL

  # expect that scenarios are comparable with baseline on parameter uncertainty
  scenarios_bad_vals <- data.table::copy(scenarios)
  scenarios_bad_vals[, transmission_rate := runif(nrow(scenarios))]
  expect_error(
    outcomes_averted(baseline, scenarios_bad_vals),
    regexp = "`scenarios` must have common infection parameter sets"
  )
})
