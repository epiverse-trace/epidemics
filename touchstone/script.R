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
  default = {
    model_default(
      population = uk_population,
      time_end = 600, increment = 1.0
    )
  },
  n = 50
)

# Vacamole model
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default = {
    model_vacamole(
      population = uk_population,
      time_end = 600, increment = 1.0
    )
  },
  n = 50
)

# Diphtheria model
touchstone::benchmark_run(
  expr_before_benchmark = {
    source("touchstone/setup_ode.R")
  },
  default = {
    model_diphtheria(
      population = uk_population,
      time_end = 600, increment = 1.0
    )
  },
  n = 50
)

# Ebola model
model_ebola(
  population = guinea_population,
  time_end = 100
)
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
  n = 50
)

# create artifacts used downstream in the GitHub Action
touchstone::benchmark_analyze()
