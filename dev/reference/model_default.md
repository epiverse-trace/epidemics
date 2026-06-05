# Model an SEIR-V epidemic with interventions

Simulate an epidemic using a deterministic, compartmental epidemic model
with the compartments "susceptible", "exposed", "infectious",
"recovered", and "vaccinated". This model can accommodate heterogeneity
in social contacts among demographic groups, as well as differences in
the sizes of demographic groups.

The `population`, `transmission_rate`, `infectiousness_rate`, and
`recovery_rate` arguments are mandatory, while passing an `intervention`
and `vaccination` are optional and can be used to simulate scenarios
with different epidemic responses or different levels of the same type
of response. See **Details** for more information.

## Usage

``` r
model_default(
  population,
  transmission_rate = 1.3/7,
  infectiousness_rate = 1/2,
  recovery_rate = 1/7,
  intervention = NULL,
  vaccination = NULL,
  time_dependence = NULL,
  time_end = 100,
  increment = 1
)
```

## Arguments

- population:

  An object of the `population` class, which holds a population contact
  matrix, a demography vector, and the initial conditions of each
  demographic group. See
  [`population()`](https://epiverse-trace.github.io/epidemics/dev/reference/population.md).

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

- intervention:

  A named list of `<intervention>`s representing optional
  non-pharmaceutical or pharmaceutical interventions applied during the
  epidemic. Only a single intervention on social contacts of the class
  `<contacts_intervention>` is allowed as the named element "contacts".
  Multiple `<rate_interventions>` on the model parameters are allowed;
  see **Details** for the model parameters for which interventions are
  supported.

- vaccination:

  A `<vaccination>` object representing an optional vaccination regime
  with a single dose, followed during the course of the epidemic, with a
  start and end time, and age-specific vaccination rates.

- time_dependence:

  A named list where each name is a model parameter, and each element is
  a function with the first two arguments being the current simulation
  `time`, and `x`, a value that is dependent on `time` (`x` represents a
  model parameter). See **Details** for more information, as well as the
  vignette on time- dependence
  `vignette("time_dependence", package = "epidemics")`.

- time_end:

  The maximum number of timesteps over which to run the model. Taken as
  days, with a default value of 100 days. May be a numeric vector.

- increment:

  The size of the time increment. Taken as days, with a default value of
  1 day.

## Value

A `<data.table>`. If the model parameters and composable elements are
all scalars, a single `<data.table>` with the columns "time",
"compartment", "age_group", and "value", giving the number of
individuals per demographic group in each compartment at each timestep
in long (or "tidy") format is returned.

If the model parameters or composable elements are lists or list-like, a
nested `<data.table>` is returned with a list column "data", which holds
the compartmental values described above. Other columns hold parameters
and composable elements relating to the model run. Columns "scenario"
and "param_set" identify combinations of composable elements
(population, interventions, vaccination regimes), and infection
parameters, respectively.

## Details: SEIRV model suitable for directly transmitted infections

### Model parameters

This model only allows for single, population-wide rates of transitions
between compartments per model run.

However, model parameters may be passed as numeric vectors. These
vectors must follow Tidyverse recycling rules: all vectors must have the
same length, or, vectors of length 1 will be recycled to the length of
any other vector.

The default values are:

- Transmission rate (\\\beta\\, `transmission_rate`): 0.186, assuming an
  \\R_0\\ = 1.3 and an infectious period of 7 days.

- Infectiousness rate (\\\sigma\\, `infectiousness_rate`): 0.5, assuming
  a pre-infectious period of 2 days.

- Recovery rate (\\\gamma\\, `recovery_rate`): 0.143, assuming an
  infectious period of 7 days.

## Examples

``` r
# create a population
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0, 0),
    nrow = 1, ncol = 5L
  )
)

# run epidemic simulation with no vaccination or intervention
# and three discrete values of transmission rate
data <- model_default(
  population = uk_population,
  transmission_rate = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
)

# view some data
data
#>    transmission_rate infectiousness_rate recovery_rate time_end param_set
#>                <num>               <num>         <num>    <num>     <int>
#> 1:         0.1857143                 0.5     0.1428571      100         1
#> 2:         0.2000000                 0.5     0.1428571      100         2
#> 3:         0.2142857                 0.5     0.1428571      100         3
#>         population intervention vaccination time_dependence increment scenario
#>             <list>       <list>      <list>          <list>     <num>    <int>
#> 1: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#> 2: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#> 3: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#>                   data
#>                 <list>
#> 1: <data.table[505x4]>
#> 2: <data.table[505x4]>
#> 3: <data.table[505x4]>

# run epidemic simulations with differences in the end time
# may be useful when considering different start dates with a fixed end point
data <- model_default(
  population = uk_population,
  time_end = c(50, 100, 150)
)

data
#>    transmission_rate infectiousness_rate recovery_rate time_end param_set
#>                <num>               <num>         <num>    <num>     <int>
#> 1:         0.1857143                 0.5     0.1428571       50         1
#> 2:         0.1857143                 0.5     0.1428571      100         2
#> 3:         0.1857143                 0.5     0.1428571      150         3
#>         population intervention vaccination time_dependence increment scenario
#>             <list>       <list>      <list>          <list>     <num>    <int>
#> 1: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#> 2: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#> 3: <population[4]>       [NULL]      [NULL]       <list[1]>         1        1
#>                   data
#>                 <list>
#> 1: <data.table[255x4]>
#> 2: <data.table[505x4]>
#> 3: <data.table[755x4]>
```
