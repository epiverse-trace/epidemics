#' Calculate outcomes averted by interventions
#'
#' @param baseline A nested `<data.table>` of a single model outcome with
#' parameter uncertainty. This is expected to be the output of a single call to
#' model functions such as [model_default()] or [model_ebola()].
#' @param scenarios A nested `<data.table>` of any number of model outcomes with
#' infection parameter sets identical to those in `baseline` for comparability,
#' i.e., the `scenarios` must differ from the `baseline` only in any
#' interventions applied.
#' @param by_group A single logical value that controls whether outcomes averted
#' are calculated separately for each demographic group in the model outputs.
#' This is passed to the `by_group` argument in [epidemic_size()]. Defaults to
#' `TRUE`.
#' @param summarise A single logical value that controls whether outcomes
#' averted are summarised (returning the median and 95% uncertainty intervals),
#' aggregating over parameter sets for each scenario, or whether the raw
#' differences between the baseline and comparator scenarios are returned while
#' matching by parameter set.
#'
#' @return A `<data.table>`.
#' When `summarise = TRUE`, a `<data.table>` of the same number of
#' rows as `scenarios`, with the columns `"scenario"`, `"averted_median"`,
#' `"averted_lower"`, and `"averted_upper"`.
#'
#' When `summarise = FALSE`, a `<data.table>` with one row per `"scenario"`,
#' parameter set (`"param_set"`), and demography group (`"demography_group`),
#' with the additional column `"outcomes_averted"` giving the difference between
#' the baseline and the comparator scenario for each parameter set for each
#' demography group.
#'
#' @details
#' Both deterministic and stochastic models (currently only the Ebola model) are
#' supported.
#'
#' When comparing deterministic model scenarios, users are expected to ensure
#' that outputs comparable in terms of demographic groups and parameters.
#' The output is expected to have parameter uncertainty, and differences between
#' each scenario and the baseline are calculated after matching on parameter
#' sets.
#'
#' When comparing stochastic model scenarios, each scenario is matched against
#' the baseline on the replicate number as well as the parameter set to reduce
#' the effect of initial conditions on differences in outcomes.
#'
#' @export
#'
#' @examples
#' polymod <- socialmixr::polymod
#' contact_data <- socialmixr::contact_matrix(
#'   polymod,
#'   countries = "United Kingdom",
#'   age.limits = c(0, 20, 40),
#'   symmetric = TRUE
#' )
#'
#' # prepare contact matrix
#' contact_matrix <- t(contact_data$matrix)
#'
#' # prepare the demography vector
#' demography_vector <- contact_data$demography$population
#' names(demography_vector) <- rownames(contact_matrix)
#'
#' # initial conditions
#' initial_i <- 1e-6
#' initial_conditions <- c(
#'   S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
#' )
#'
#' # build for all age groups
#' initial_conditions <- rbind(
#'   initial_conditions,
#'   initial_conditions,
#'   initial_conditions
#' )
#'
#' # create population object
#' uk_population <- population(
#'   name = "UK",
#'   contact_matrix = contact_matrix,
#'   demography_vector = demography_vector,
#'   initial_conditions = initial_conditions
#' )
#'
#' # create vector of parameters
#' beta <- withr::with_seed(
#'   1,
#'   rnorm(100, mean = 1.3 / 7, sd = 0.005)
#' )
#'
#' baseline <- model_default(
#'   population = uk_population,
#'   transmission_rate = beta
#' )
#'
#' max_time <- 100
#' # prepare durations as starting at 25% of the way through an epidemic
#' # and ending halfway through
#' time_begin <- max_time / 4
#' time_end <- max_time / 2
#'
#' # create three distinct contact interventions
#' # prepare an intervention that models school closures for 180 days
#' close_schools <- intervention(
#'   name = "School closure",
#'   type = "contacts",
#'   time_begin = time_begin,
#'   time_end = time_end,
#'   reduction = matrix(c(0.3, 0.01, 0.01))
#' )
#'
#' # prepare an intervention which mostly affects adults 20 -- 65
#' close_workplaces <- intervention(
#'   name = "Workplace closure",
#'   type = "contacts",
#'   time_begin = time_begin,
#'   time_end = time_end,
#'   reduction = matrix(c(0.01, 0.3, 0.01))
#' )
#'
#' intervention_sets <- list(
#'   list(
#'     contacts = close_schools
#'   ),
#'   list(
#'     contacts = close_workplaces
#'   )
#' )
#'
#' scenarios <- model_default(
#'   population = uk_population,
#'   transmission_rate = beta,
#'   intervention = intervention_sets
#' )
#'
#' # Defaults to summarise = TRUE
#' outcomes_averted(
#'   baseline = baseline,
#'   scenarios = scenarios
#' )
#'
#' # Set summarise = FALSE to get raw difference data
#' outcomes_averted(
#'   baseline = baseline,
#'   scenarios = scenarios,
#'   summarise = FALSE
#' )
outcomes_averted <- function(baseline,
                             scenarios,
                             by_group = TRUE,
                             summarise = TRUE) {
  # Assign `outcome` as NULL to avoid global variable messages
  outcome <- NULL
  # TODO: reconsider whether min.cols should be checked
  checkmate::assert_data_table(
    baseline,
    min.rows = 1, min.cols = 12L # known number of columns in output
  )
  # expect column names correspond to nested output with common infection params
  # NOTE: `recovery_rate` not expected as it is replaced by `removal_rate` for
  # ebola model
  stopifnot(
    "`baseline` must be a nested <data.table> with expected column names" =
      checkmate::test_names(
        colnames(baseline),
        must.include = c(
          "transmission_rate", "infectiousness_rate", "time_end",
          "param_set", "population", "intervention",
          "time_dependence", "scenario", "data"
        )
      )
  )
  checkmate::assert_data_table(scenarios, min.rows = 1, min.cols = 12L)
  checkmate::assert_names(
    colnames(scenarios),
    identical.to = colnames(baseline)
  )
  checkmate::assert_logical(summarise, any.missing = FALSE, len = 1)
  checkmate::assert_logical(by_group, any.missing = FALSE, len = 1)

  # Check that scenarios have comparable parameters
  # First collect scenario parameter names that are not checked for equality
  names_scenario_params <- c(
    "population", "intervention", "vaccination",
    "time_dependence", "increment", "scenario", "data"
  )
  # Expect that parameter sets are shared across baseline and scenarios
  stopifnot(
    "`baseline` and `scenarios` must have common infection parameter sets" =
      identical(
        baseline[, !colnames(baseline) %in% names_scenario_params,
          with = FALSE
        ],
        unique(
          scenarios[, !colnames(scenarios) %in% names_scenario_params,
            with = FALSE
          ]
        )
      )
  )

  # get epidemic size for baseline and response scenarios
  # NOTE: make copy of data.table as assignment by reference will modify
  # original data.table
  baseline_outcomes <- data.table::copy(baseline)
  # suppressed as most models do not track deaths
  suppressWarnings(
    baseline_outcomes[, "outcome" := lapply(
      baseline_outcomes$data, epidemic_size,
      simplify = FALSE, by_group = by_group, include_deaths = TRUE
    )]
  )

  # set flag for whether data has replicates; applies only to stochastic model
  has_reps <- if ("replicate" %in% colnames(baseline[["data"]][[1L]])) {
    "replicate"
  } else {
    NULL
  }

  baseline_outcomes <- baseline_outcomes[, unlist(outcome, recursive = FALSE),
    by = c("scenario", "param_set")
  ]
  data.table::setnames(baseline_outcomes, "value", "baseline_value")

  # get epidemic sizes for response scenarios; some code repetition here
  # NOTE: work with a copy
  scenario_outcomes <- data.table::copy(scenarios)
  suppressWarnings(
    scenario_outcomes[, "outcome" := lapply(
      scenario_outcomes$data, epidemic_size,
      simplify = FALSE, by_group = by_group, include_deaths = TRUE
    )]
  )

  scenario_outcomes <- scenario_outcomes[, unlist(outcome, recursive = FALSE),
    by = c("scenario", "param_set")
  ]

  # variables to merge on; # cannot use ifelse()
  merge_variables <- c(
    "param_set", if (by_group) "demography_group" else NULL, has_reps
  )

  # merge and get difference in outcomes
  scenario_outcomes <- data.table::merge.data.table(
    scenario_outcomes,
    baseline_outcomes[, c(merge_variables, "baseline_value"), with = FALSE],
    by = merge_variables
  )
  # assume baseline has more cases
  scenario_outcomes$outcomes_averted <- scenario_outcomes$baseline_value -
    scenario_outcomes$value

  # select relevant columns
  scenario_outcomes <- scenario_outcomes[, c(
    "scenario", merge_variables, "outcomes_averted"
  ), with = FALSE]

  # return summarised values (median, 95% uncertainty) or raw differences
  if (summarise) {
    grouping_vars <- c("scenario", if (by_group) "demography_group" else NULL)
    averted_summary <- scenario_outcomes[, list(
      averted_median = stats::median(outcomes_averted),
      averted_lower = stats::quantile(outcomes_averted, probs = 0.025),
      averted_upper = stats::quantile(outcomes_averted, probs = 0.975)
    ), by = grouping_vars]

    # return difference summary
    data.table::setorderv(averted_summary, grouping_vars)
    averted_summary[]
  } else {
    data.table::setorderv(
      scenario_outcomes,
      cols = c("scenario", merge_variables)
    )
    scenario_outcomes[, colnames(scenario_outcomes) != "time", with = FALSE]
  }
}
