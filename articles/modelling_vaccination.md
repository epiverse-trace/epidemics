# Modelling the effect of a vaccination campaign

``` r

library(epidemics)
```

## Prepare population and initial conditions

Prepare population and contact data.

``` r

# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod

# demography data from the wpp2024 package
data("popAge1dt", package = "wpp2024")
uk_pop <- popAge1dt |>
  dplyr::filter(name == "United Kingdom", year == 2006) |>
  dplyr::select(lower.age.limit = age, population = pop) |>
  dplyr::mutate(population = population * 1000)

contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  survey_pop = uk_pop,
  age_limits = c(0, 20, 65),
  symmetric = TRUE,
  return_demography = TRUE
)
#> Warning in normalise_weighted_matrix(survey_pop = survey_pop, weighted_matrix = weighted.matrix, : Large differences in the size of the sub-populations with the current age
#> breaks are likely to result in artefacts after making the matrix symmetric.
#> ! Please reconsider the age breaks to obtain more equally sized
#>   sub-populations.
#> ℹ Normalization factors: [0.5 and 2.2]

# prepare contact matrix
contact_matrix <- t(contact_data$matrix)

# prepare the demography vector
demography_vector <- contact_data$demography$population
names(demography_vector) <- rownames(contact_matrix)
```

Prepare initial conditions for each age group.

``` r

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
```

Prepare a population as a `population` class object.

``` r

uk_population <- population(
  name = "UK",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)
```

## Prepare a vaccination campaign

Prepare a vaccination campaign targeting individuals aged over 65.

``` r

# prepare a vaccination object
vaccinate_elders <- vaccination(
  name = "vaccinate elders",
  time_begin = matrix(100, nrow(contact_matrix)),
  time_end = matrix(250, nrow(contact_matrix)),
  nu = matrix(c(0, 0, 0.0001))
)

# view vaccination object
vaccinate_elders
#> <vaccination> object
#> 
#>  Vaccination name: 
#> "vaccinate elders"
#> 
#>  Begins at: 
#>      dose_1
#> [1,]    100
#> [2,]    100
#> [3,]    100
#> 
#>  Ends at: 
#>      dose_1
#> [1,]    250
#> [2,]    250
#> [3,]    250
#> 
#>  Vaccination rate: 
#>      dose_1
#> [1,]  0e+00
#> [2,]  0e+00
#> [3,]  1e-04
```

**Note that** when the vaccination rate is high, the number of
individuals who could *potentially* transition from ‘susceptible’ to
‘vaccinated’ may be greater than the number of individuals in the
‘susceptible’ compartment:
$`dS_{i_{t}} = \nu_{i_t}S_{i_{0}} > S_{i_{t}}`$ for each demographic
group $`i`$ at time $`t`$.

This could realistically occur when there are more doses available than
individuals eligible to receive them, for instance towards the end of an
epidemic.

*epidemics* automatically handles such situations by setting
$`dS_{i_{t}}`$ to be the minimum of doses or eligible individuals: :
$`dS_{i_t} = \text{Min}(\nu_{i_t}S_{i_0}, S_{i_t})`$, such that
$`S_{i_{t+1}}`$ does not take a negative value.

The `<vaccination>` object can be passed to the `vaccination` argument
of a model function call.

``` r

# run an epidemic model using `epidemic`
output <- model_default(
  population = uk_population,
  vaccination = vaccinate_elders,
  time_end = 600, increment = 1.0
)
```

**Note that** vaccination modelling is currently only supported for the
‘default’
([`model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md))
and ‘Vacamole’
([`model_vacamole()`](https://epiverse-trace.github.io/epidemics/reference/model_vacamole.md))
models.
