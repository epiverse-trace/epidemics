# Model leaky, two-dose vaccination in an epidemic using Vacamole

Simulate an epidemic using the *Vacamole* model for Covid-19 developed
at RIVM, the National Institute for Public Health and the Environment in
the Netherlands. This model is aimed at estimating the impact of 'leaky'
vaccination on an epidemic. See **Details** and **References** for more
information.

## Usage

``` r
model_vacamole(
  population,
  transmission_rate = 1.3/7,
  transmission_rate_vax = 0.8 * transmission_rate,
  infectiousness_rate = 1/2,
  hospitalisation_rate = 1/1000,
  hospitalisation_rate_vax = 0.8 * hospitalisation_rate,
  mortality_rate = 1/1000,
  mortality_rate_vax = 0.8 * mortality_rate,
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

- transmission_rate_vax:

  A numeric of values between 0.0 and 1.0 giving the transmission_rate
  of the infection to individuals who have received two doses of the
  vaccine. The default values is 80% of the transmission_rate of the
  infection to individuals who are not doubly vaccinated.

- infectiousness_rate:

  A numeric for the rate at which individuals move from the exposed to
  the infectious compartment. Often denoted as \\\sigma\\, with \\\sigma
  = 1.0 / \text{pre-infectious period}\\. This value does not depend
  upon the number of infectious individuals in the population. See
  **Details** for default values.

- hospitalisation_rate:

  A numeric for the hospitalisation rate of infectious individuals.

- hospitalisation_rate_vax:

  A numeric of values between 0.0 and 1.0 giving the hospitalisation
  rate of infectious individuals who have received two doses of the
  vaccine. The default value is 80% of the hospitalisation rate of
  individuals who are not doubly vaccinated.

- mortality_rate:

  A numeric for the mortality rate of infectious or hospitalised
  individuals.

- mortality_rate_vax:

  A numeric of values between 0.0 and 1.0 giving the mortality of
  infectious and hospitalised individuals who have received two doses of
  the vaccine. The default value is 80% of the mortality rate of
  individuals who are not doubly vaccinated.

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

  An optional `<vaccination>` object representing a vaccination regime
  with **two doses** followed during the course of the epidemic, with a
  start and end time, and age-specific vaccination rates for each dose.
  See
  [`vaccination()`](https://epiverse-trace.github.io/epidemics/dev/reference/vaccination.md).

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

## Details: Vacamole Covid-19 model with leaky, two-dose vaccination

The Vacamole model has the compartments "susceptible",
"vaccinated_one_dose", "vaccinated_two_dose", "exposed", "infectious"
"infectious_vaccinated", "hospitalised", "hospitalised_vaccinated",
"recovered", and "dead".

This model allows for:

1.  A 'hospitalised' compartment along with a hospitalisation rate;

2.  Two doses of vaccination, with 'leaky' protection, i.e., vaccination
    does not prevent infection completely but allows for a reduction in
    the infection rate, as well as reduced rates of moving into states
    considered more serious, such as 'hospitalised' or 'dead'.

### Model parameters

This model only allows for single, population-wide rates of transitions
between compartments per model run.

However, model parameters may be passed as numeric vectors. These
vectors must follow Tidyverse recycling rules: all vectors must have the
same length, or, vectors of length 1 will be recycled to the length of
any other vector.

- Transmission rate (\\\beta\\, `transmission_rate`): 0.186, resulting
  from an \\R_0\\ = 1.3 and an infectious period of 7 days. The
  transmission rate for doubly vaccinated individuals (\\\beta_v\\) is
  80% of \\\beta\\, 0.1488.

- Infectiousness rate (\\\sigma\\, `infectiousness_rate`): 0.5, assuming
  a pre-infectious period of 2 days.

- Hospitalisation rate (\\\eta\\, `hospitalisation_rate`): 1.0 / 1000,
  assuming that one in every thousand infectious individuals is
  hospitalised. The hospitalisation rate of doubly vaccinated
  individuals (\\\eta_v\\) is 80% of \\\eta\\, 0.8 / 1000.

- Mortality rate (\\\omega\\, `mortality_rate`): 1.0 / 1000, assuming
  that one in every thousand infectious and hospitalised individuals
  dies. The mortality rate of the doubly vaccinated (\\\omega_v\\) is
  80% of \\\omega\\, 0.8 / 1000.

- Recovery rate (\\\gamma\\, `recovery_rate`): 0.143, assuming an
  infectious period of 7 days.

## References

Ainslie, K. E. C., Backer, J. A., Boer, P. T. de, Hoek, A. J. van,
Klinkenberg, D., Altes, H. K., Leung, K. Y., Melker, H. de, Miura, F., &
Wallinga, J. (2022). A scenario modelling analysis to anticipate the
impact of COVID-19 vaccination in adolescents and children on disease
outcomes in the Netherlands, summer 2021. Eurosurveillance, 27(44),
2101090.
[doi:10.2807/1560-7917.ES.2022.27.44.2101090](https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090)

## Examples

``` r
# create a population, note eleven columns for compartments
population <- population(
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0, 0, 0, 0, 0.0001, 0, 0, 0, 0, 0),
    nrow = 1, ncol = 11L
  )
)

# create a vaccination regime
double_vax <- vaccination(
  nu = matrix(1e-3, ncol = 2, nrow = 1),
  time_begin = matrix(c(10, 30), nrow = 1),
  time_end = matrix(c(50, 80), nrow = 1)
)

# run epidemic simulation with vaccination but no intervention
# with a single set of parameters
data <- model_vacamole(
  population = population,
  vaccination = double_vax
)

# view some data
head(data)
#>     time demography_group         compartment    value
#>    <num>           <char>              <char>    <num>
#> 1:     0     demo_group_1         susceptible 66993300
#> 2:     0     demo_group_1 vaccinated_one_dose        0
#> 3:     0     demo_group_1 vaccinated_two_dose        0
#> 4:     0     demo_group_1             exposed        0
#> 5:     0     demo_group_1  exposed_vaccinated        0
#> 6:     0     demo_group_1          infectious     6700

# run epidemic simulation with no vaccination or intervention
# and three discrete values of transmission_rate
data <- model_vacamole(
  population = population,
  transmission_rate = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
)

# view some data
head(data)
#>    transmission_rate infectiousness_rate hospitalisation_rate mortality_rate
#>                <num>               <num>                <num>          <num>
#> 1:         0.1857143                 0.5                0.001          0.001
#> 2:         0.2000000                 0.5                0.001          0.001
#> 3:         0.2142857                 0.5                0.001          0.001
#>    recovery_rate transmission_rate_vax hospitalisation_rate_vax
#>            <num>                 <num>                    <num>
#> 1:     0.1428571             0.1485714                    8e-04
#> 2:     0.1428571             0.1600000                    8e-04
#> 3:     0.1428571             0.1714286                    8e-04
#>    mortality_rate_vax time_end param_set      population intervention
#>                 <num>    <num>     <int>          <list>       <list>
#> 1:              8e-04      100         1 <population[4]>       [NULL]
#> 2:              8e-04      100         2 <population[4]>       [NULL]
#> 3:              8e-04      100         3 <population[4]>       [NULL]
#>    vaccination time_dependence increment scenario                 data
#>         <list>          <list>     <num>    <int>               <list>
#> 1:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
#> 2:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
#> 3:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
tail(data)
#>    transmission_rate infectiousness_rate hospitalisation_rate mortality_rate
#>                <num>               <num>                <num>          <num>
#> 1:         0.1857143                 0.5                0.001          0.001
#> 2:         0.2000000                 0.5                0.001          0.001
#> 3:         0.2142857                 0.5                0.001          0.001
#>    recovery_rate transmission_rate_vax hospitalisation_rate_vax
#>            <num>                 <num>                    <num>
#> 1:     0.1428571             0.1485714                    8e-04
#> 2:     0.1428571             0.1600000                    8e-04
#> 3:     0.1428571             0.1714286                    8e-04
#>    mortality_rate_vax time_end param_set      population intervention
#>                 <num>    <num>     <int>          <list>       <list>
#> 1:              8e-04      100         1 <population[4]>       [NULL]
#> 2:              8e-04      100         2 <population[4]>       [NULL]
#> 3:              8e-04      100         3 <population[4]>       [NULL]
#>    vaccination time_dependence increment scenario                 data
#>         <list>          <list>     <num>    <int>               <list>
#> 1:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
#> 2:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
#> 3:      [NULL]       <list[1]>         1        1 <data.table[1111x4]>
```
