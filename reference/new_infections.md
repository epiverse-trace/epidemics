# Get new infections over model time

Get new infections over model time

## Usage

``` r
new_infections(data, exclude_compartments = NULL, by_group = TRUE)
```

## Arguments

- data:

  A table of model output, typically the output of
  [`model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md)
  or similar functions.

- exclude_compartments:

  An optional argument, for a character vector of the names of model
  compartments into which individuals transition from the "susceptible"
  compartment, and which are not related to infection. A common example
  is a compartment for "vaccinated" individuals who are no longer
  susceptible, but who should also not be counted as infected.

- by_group:

  A logical representing whether the epidemic size should be returned by
  demographic group, or whether a single population-wide value is
  returned.

## Value

A table with the same columns as `data`, but with the additional
variable under `compartment`, "new_infections", resulting in additional
rows.

## Examples

``` r
# create a population
uk_population <- population(
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0, 0),
    nrow = 1, ncol = 5L
  )
)


# run epidemic simulation with no vaccination or intervention
data <- model_default(
  population = uk_population,
  time_end = 200,
  increment = 1
)

new_infections(data)
#> Key: <time, demography_group>
#>       time demography_group new_infections
#>      <num>           <char>          <num>
#>   1:     0     demo_group_1         0.0000
#>   2:     1     demo_group_1       254.2911
#>   3:     2     demo_group_1       597.4234
#>   4:     3     demo_group_1       787.2866
#>   5:     4     demo_group_1       899.5670
#>  ---                                      
#> 197:   196     demo_group_1    194151.2386
#> 198:   197     demo_group_1    196154.1154
#> 199:   198     demo_group_1    198089.5621
#> 200:   199     demo_group_1    199954.4072
#> 201:   200     demo_group_1    201745.5572
```
