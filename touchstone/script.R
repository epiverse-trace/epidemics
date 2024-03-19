# see `help(run_script, package = 'touchstone')` on how to run this
# interactively

# TODO OPTIONAL Add directories you want to be available in this file or during
# the benchmarks.
# touchstone::pin_assets("some/dir")

# installs branches to benchmark
touchstone::branch_install()

#### Models without interventions ####
# NOTE: all use default SCALAR parameters
# Default model
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default_ode = {
    model_default(
      population = uk_population,
      time_end = max_time
    )
  },
  n = 25
)

# Default model with parameter vectors
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default_ode_param_vec = {
    model_default(
      population = uk_population,
      transmission_rate = beta,
      time_end = max_time
    )
  },
  n = 25
)

# Default model with intervention scenarios
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default_ode_interventions = {
    model_default(
      population = uk_population,
      intervention = intervention_scenarios,
      time_end = max_time
    )
  },
  n = 25
)

# Default model with parameter vector and intervention scenarios
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default_ode_paramvec_intervs = {
    model_default(
      population = uk_population,
      transmission_rate = beta,
      intervention = intervention_scenarios,
      time_end = max_time
    )
  },
  n = 25
)

# Ebola model
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ebola.R")
  },
  default = {
    withr::with_seed(
      seed = 1,
      {
        model_ebola(
          population = guinea_population,
          time_end = 100
        )
      }
    )
  },
  n = 25
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
