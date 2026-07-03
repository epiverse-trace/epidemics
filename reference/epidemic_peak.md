# Get the time and size of a compartment's highest peak

Get the time and size of a compartment's highest peak for all demography
groups.

## Usage

``` r
epidemic_peak(data, compartments = "infectious")
```

## Arguments

- data:

  A `<data.frame>` or `<data.table>` of model output, typically the
  output of a compartmental model.

- compartments:

  A character vector for the compartments of interest.

## Value

A `<data.table>` with columns "demography_group", "compartment", "time"
and "value"; these specify the name of the demography group, the
epidemiological compartment, and the peak time and value for each
compartment in `compartments`.

## Details

This is used for epidemics with a single peak. It is useful from a
public health policy point of view to determine how bad an epidemic will
be and when that happens.

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
data <- model_default(
  population = uk_population,
  time_end = 600
)

# get the timing and peak of the exposed and infectious compartment
epidemic_peak(data, c("exposed", "infectious"))
#>    demography_group compartment  time   value
#>              <char>      <char> <num>   <num>
#> 1:     demo_group_1     exposed   220  436788
#> 2:     demo_group_1  infectious   227 1511761
```
