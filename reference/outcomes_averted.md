# Calculate outcomes averted by interventions

Calculate outcomes averted by interventions

## Usage

``` r
outcomes_averted(baseline, scenarios, by_group = TRUE, summarise = TRUE)
```

## Arguments

- baseline:

  A nested `<data.table>` of a single model outcome with parameter
  uncertainty. This is expected to be the output of a single call to
  model functions such as
  [`model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md)
  or
  [`model_ebola()`](https://epiverse-trace.github.io/epidemics/reference/model_ebola.md).

- scenarios:

  A nested `<data.table>` of any number of model outcomes with infection
  parameter sets identical to those in `baseline` for comparability,
  i.e., the `scenarios` must differ from the `baseline` only in any
  interventions applied.

- by_group:

  A single logical value that controls whether outcomes averted are
  calculated separately for each demographic group in the model outputs.
  This is passed to the `by_group` argument in
  [`epidemic_size()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_size.md).
  Defaults to `TRUE`.

- summarise:

  A single logical value that controls whether outcomes averted are
  summarised (returning the median and 95% uncertainty intervals),
  aggregating over parameter sets for each scenario, or whether the raw
  differences between the baseline and comparator scenarios are returned
  while matching by parameter set.

## Value

A `<data.table>`. When `summarise = TRUE`, a `<data.table>` of the same
number of rows as `scenarios`, with the columns `"scenario"`,
`"averted_median"`, `"averted_lower"`, and `"averted_upper"`.

When `summarise = FALSE`, a `<data.table>` with one row per
`"scenario"`, parameter set (`"param_set"`), and demography group
(`"demography_group`), with the additional column `"outcomes_averted"`
giving the difference between the baseline and the comparator scenario
for each parameter set for each demography group.

## Details

Both deterministic and stochastic models (currently only the Ebola
model) are supported.

When comparing deterministic model scenarios, users are expected to
ensure that outputs comparable in terms of demographic groups and
parameters. The output is expected to have parameter uncertainty, and
differences between each scenario and the baseline are calculated after
matching on parameter sets.

When comparing stochastic model scenarios, each scenario is matched
against the baseline on the replicate number as well as the parameter
set to reduce the effect of initial conditions on differences in
outcomes.

## Examples

``` r
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age_limits = c(0, 20, 40),
  symmetric = TRUE,
  return_demography = TRUE
)
#> Warning: Automatic country population lookup in `contact_matrix()` was deprecated in
#> socialmixr 0.6.0.
#> When `countries` is given (or a `country` column is present) without
#> `survey_pop`, contact_matrix() currently calls the soft-deprecated `wpp_age()`
#> to look up population data. This automatic lookup will be removed in a future
#> release: callers will then have to supply `survey_pop` whenever `symmetric`,
#> `split`, `per_capita`, `weigh_age`, or `return_demography` is TRUE.
#> ℹ Pass `survey_pop` explicitly to silence this warning, e.g. `survey_pop =
#>   survey_country_population(survey, countries)` or a data frame from the
#>   wpp2024 package.

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

intervention_sets <- list(
  list(
    contacts = close_schools
  ),
  list(
    contacts = close_workplaces
  )
)

scenarios <- model_default(
  population = uk_population,
  transmission_rate = beta,
  intervention = intervention_sets
)

# Defaults to summarise = TRUE
outcomes_averted(
  baseline = baseline,
  scenarios = scenarios
)
#>    scenario demography_group averted_median averted_lower averted_upper
#>       <int>           <char>          <num>         <num>         <num>
#> 1:        1           [0,20)       821.4310      501.3149     1314.9918
#> 2:        1          [20,40)       609.9044      362.8535      996.3781
#> 3:        1         [40,Inf)       731.7277      433.9860     1198.3892
#> 4:        2           [0,20)       560.9640      332.7719      918.8973
#> 5:        2          [20,40)       533.9183      326.5801      853.6345
#> 6:        2         [40,Inf)       597.1504      357.3972      971.1120

# Set summarise = FALSE to get raw difference data
outcomes_averted(
  baseline = baseline,
  scenarios = scenarios,
  summarise = FALSE
)
#>      scenario param_set demography_group outcomes_averted
#>         <int>     <int>           <char>            <num>
#>   1:        1         1           [0,20)         668.7940
#>   2:        1         1          [20,40)         491.6037
#>   3:        1         1         [40,Inf)         589.0750
#>   4:        1         2           [0,20)         837.4592
#>   5:        1         2          [20,40)         622.3674
#>  ---                                                     
#> 596:        2        99          [20,40)         368.9478
#> 597:        2        99         [40,Inf)         406.1099
#> 598:        2       100           [0,20)         472.3195
#> 599:        2       100          [20,40)         453.8384
#> 600:        2       100         [40,Inf)         504.1885
```
