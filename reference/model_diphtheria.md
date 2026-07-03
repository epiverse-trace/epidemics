# Model a diphtheria outbreak using a compartmental ODE model

Simulate a diphtheria outbreak using a deterministic, compartmental
ordinary differential equation model with the compartments
"susceptible", "exposed", "infectious", "hospitalised", and"recovered".
The model is based on Finger et al. (2019) and is intended to be used in
the context of internally displaced people (IDP) or refugee camps. This
model allows for a proportion of each demographic group to be vaccinated
at the start of the outbreak, and thus to not contribute to the
outbreak. The model also allows for changes to the number of
susceptibles in each age group to model influxes or evacuations from
camps.

## Usage

``` r
model_diphtheria(
  population,
  transmission_rate = 4/4.5,
  infectiousness_rate = 1/3,
  recovery_rate = 1/3,
  reporting_rate = 0.03,
  prop_hosp = 0.01,
  hosp_entry_rate = 0.2,
  hosp_exit_rate = 0.2,
  prop_vaccinated = 0 * population[["demography_vector"]],
  intervention = NULL,
  time_dependence = NULL,
  population_change = NULL,
  time_end = 100,
  increment = 1
)
```

## Arguments

- population:

  An object of the `population` class, which holds a population contact
  matrix, a demography vector, and the initial conditions of each
  demographic group. See
  [`population()`](https://epiverse-trace.github.io/epidemics/reference/population.md).

- transmission_rate:

  A numeric for the rate at which individuals move from the susceptible
  to the exposed compartment upon contact with an infectious individual.
  Often denoted as \\\beta\\, with \\\beta = R_0 / \text{infectious
  period}\\. See **Details** for default values.

- infectiousness_rate:

  A numeric for the rate at which individuals move from the exposed to
  the infectious compartment. Often denoted as \\\sigma\\, with \\\sigma
  = 1.0 / \text{pre-infectious period}\\. This value does not depend
  upon the number of infectious individuals in the population. See
  **Details** for default values.

- recovery_rate:

  A numeric for the rate at which individuals move from the infectious
  to the recovered compartment. Often denoted as \\\gamma\\, with
  \\\gamma = 1.0 / \text{infectious period}\\. See **Details** for
  default values.

- reporting_rate:

  A numeric for the proportion of infectious cases that is reported;
  this is a precursor to hospitalisation as only reported cases are
  hospitalised.

- prop_hosp:

  A numeric for the proportion of reported cases that is hospitalised.

- hosp_entry_rate:

  A numeric for the rate at which reported cases of infectious
  individuals are hospitalised. This is calculated as 1 / time to
  hospitalisation, denoted \\\tau_1\\.

- hosp_exit_rate:

  A numeric for the rate at which individuals are discharged from
  hospital to enter the 'recovered' compartment. This is calculated as 1
  / time to discharge, denoted \\\tau_2\\.

- prop_vaccinated:

  A numeric vector of the same length as the number of demographic
  groups indicated the proportion of each group that is vaccinated.
  These individuals are not included in the model dynamics.

- intervention:

  A named list of `<rate_intervention>` objects representing optional
  pharmaceutical or non-pharmaceutical interventions applied to the
  model parameters listed above.

- time_dependence:

  A named list where each name is a model parameter, and each element is
  a function with the first two arguments being the current simulation
  `time`, and `x`, a value that is dependent on `time` (`x` represents a
  model parameter). See **Details** for more information, as well as the
  vignette on time- dependence
  `vignette("time_dependence", package = "epidemics")`.

- population_change:

  A two-element list, with elements named `"time"` and `"values"`,
  giving the times of population changes, and the corresponding changes
  in the population of each demographic group at those times. `"time"`
  must be a numeric vector, while `"values"` must be a list of the
  length of `"time"`, with each element a numeric vector of the same
  length as the number of demographic groups in `population`.

- time_end:

  The maximum number of timesteps over which to run the model. Taken as
  days, with a default value of 100 days. May be a numeric vector.

- increment:

  The size of the time increment. Taken as days, with a default value of
  1 day.

## Value

A `data.table` with the columns "time", "compartment", "age_group",
"value", and "run", giving the number of individuals per demographic
group in each compartment at each timestep in long (or "tidy") format,
with "run" indicating the unique parameter combination.

## Details: Model an infection outbreak in a humanitarian camp setting

This model has been developed for diphtheria outbreaks in settings where
interventions on social contacts are difficult to implement. It it
suitable for application to the outbreak of similar, directly
transmitted infectious diseases as well.

### Model parameters

This model only allows for single, population-wide rates transitions
between compartments per model run.

However, model parameters may be passed as numeric vectors. These
vectors must follow Tidyverse recycling rules: all vectors must have the
same length, or, vectors of length 1 will be recycled to the length of
any other vector.

The default values are taken from Finger et al. (2019) where possible:

- Transmission rate (\\\beta\\, `transmission_rate`): 0.8888889,
  assuming an \\R_0\\ of 4.0 and a total infectious period of 4.5 days.

- Infectiousness rate (\\\sigma\\, `infectiousness_rate`): 0.333,
  assuming a pre-infectious period of 3 days.

- Reporting rate (\\r\\, `reporting_rate`): 0.03, assuming that 3% of
  infectious cases are detected or reported.

- Proportion hospitalised (\\\eta\\, `prop_hosp`): 0.01, assuming that
  1% of reported cases need hospital treatment.

- Hospital entry rate (\\\tau_1\\, `hosp_entry_rate`): 0.2, assuming
  that it takes 5 days for infectious individuals to seek hospital
  treatment.

- Hospital exit rate (\\\tau_2\\, `hosp_exit_rate`): 0.2, assuming that
  individuals are discharged from hospital after 5 days.

- Recovery rate (\\\gamma\\, `recovery_rate`): 0.333, assuming an
  infectious period following symptoms, of 3 days.

### Modelling population changes

This model allows changes to the number of susceptibles in each
demographic group, to represent influxes or evacuations from the camp as
would be expected in humanitarian relief situations. Users can specify
the times and changes (to each demographic group) of changes using the
`population_changes` argument, to examine the effect on outbreak
dynamics.

## References

Finger, F., Funk, S., White, K., Siddiqui, M. R., Edmunds, W. J., &
Kucharski, A. J. (2019). Real-time analysis of the diphtheria outbreak
in forcibly displaced Myanmar nationals in Bangladesh. BMC Medicine, 17,
58.
[doi:10.1186/s12916-019-1288-7](https://doi.org/10.1186/s12916-019-1288-7)
.

## Examples

``` r
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

# assume younger age groups are vaccinated
prop_vaccinated <- c(0.2, 0.10, 0.1)

# run model for single, default parameter set
data <- model_diphtheria(
  camp_pop,
  prop_vaccinated = prop_vaccinated
)
head(data)
#>     time demography_group compartment    value
#>    <num>           <char>      <char>    <num>
#> 1:     0     demo_group_1 susceptible  66399.2
#> 2:     0     demo_group_2 susceptible  97379.1
#> 3:     0     demo_group_3 susceptible 202139.1
#> 4:     0     demo_group_1     exposed      0.0
#> 5:     0     demo_group_2     exposed      0.0
#> 6:     0     demo_group_3     exposed      0.0
tail(data)
#>     time demography_group  compartment        value
#>    <num>           <char>       <char>        <num>
#> 1:   100     demo_group_1 hospitalised 2.596850e-02
#> 2:   100     demo_group_2 hospitalised 3.808464e-02
#> 3:   100     demo_group_3 hospitalised 7.905592e-02
#> 4:   100     demo_group_1    recovered 6.052887e+04
#> 5:   100     demo_group_2    recovered 8.876938e+04
#> 6:   100     demo_group_3    recovered 1.842660e+05

# run model with increase in population
# create population change data
p <- list(
  time = 70,
  values = list(
    c(1e4, 2e5, 1e5)
  )
)

data <- model_diphtheria(
  camp_pop,
  prop_vaccinated = prop_vaccinated,
  population_change = p
)
head(data)
#>     time demography_group compartment    value
#>    <num>           <char>      <char>    <num>
#> 1:     0     demo_group_1 susceptible  66399.2
#> 2:     0     demo_group_2 susceptible  97379.1
#> 3:     0     demo_group_3 susceptible 202139.1
#> 4:     0     demo_group_1     exposed      0.0
#> 5:     0     demo_group_2     exposed      0.0
#> 6:     0     demo_group_3     exposed      0.0
tail(data)
#>     time demography_group  compartment        value
#>    <num>           <char>       <char>        <num>
#> 1:   100     demo_group_1 hospitalised 2.476344e-01
#> 2:   100     demo_group_2 hospitalised 2.891514e+00
#> 3:   100     demo_group_3 hospitalised 1.702774e+00
#> 4:   100     demo_group_1    recovered 6.624751e+04
#> 5:   100     demo_group_2    recovered 1.743403e+05
#> 6:   100     demo_group_3    recovered 2.306429e+05
```
