# Get the epidemic size

Get the size of the epidemic at any stage between the start and the end.
This is calculated as the number of individuals *recovered* from
infection at that stage of the epidemic.

## Usage

``` r
epidemic_size(
  data,
  stage = 1,
  time = NULL,
  by_group = TRUE,
  include_deaths = FALSE,
  simplify = TRUE
)
```

## Arguments

- data:

  A table of model output, typically the output of
  [`model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md)
  or similar functions.

- stage:

  A numeric vector for the stage of the epidemic at which to return the
  epidemic size; here 0.0 represents the start time of the epidemic
  i.e., the initial conditions of the epidemic simulation, while 1.0
  represents the end of the epidemic simulation model (100% of model
  time). Defaults to 1.0, at which stage returned values represent the
  *final size* of the epidemic. This value is overridden by any values
  passed to the `time` argument.

- time:

  Alternative to `stage`, an integer-like vector for the timepoint of
  the epidemic at which to return the epidemic size. Overrides any
  values passed to `stage`.

- by_group:

  A logical representing whether the epidemic size should be returned by
  demographic group, or whether a single population-wide value is
  returned. Defaults to `TRUE`.

- include_deaths:

  A logical value that indicates whether to count dead individuals in
  the epidemic size calculation. Defaults to `FALSE`. Setting
  `include_deaths = TRUE` makes the function look for a `"dead"`
  compartment in the data. If there is no such column, the function
  returns only the final number of recovered or removed individuals in
  each demographic group.

- simplify:

  A logical determining whether the epidemic size data should be
  simplified to a vector with one element for each demographic group. If
  the length of `stage` or `time` is \$\>\$ 1, this argument is
  overridden and the data are returned as a `<data.table>`.

## Value

If `simplify == TRUE` and a single timepoint is requested, returns a
vector of epidemic sizes of the same length as the number of demographic
groups. If `by_group == FALSE`, sums the epidemic size to return an
overall value for the full population.

If multiple timepoints are requested, or if multiple replicates are
present under a specially named column "replicate" (only from the Ebola
model), no simplification to a vector is possible; returns a
`<data.table>` of timepoints and epidemic sizes at each timepoint.

All options return the absolute sizes and not proportions.

## Details

This function can be used to calculate the *final size* of the epidemic,
by setting `stage = 1.0` (100% of model time; the default).

The function allows for the calculation of epidemic sizes by demographic
group as well as the total epidemic size.

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
  population = uk_population
)

# get the final epidemic size if no other arguments are specified
epidemic_size(data)
#> [1] 481200.3

# get the epidemic size at the halfway point
epidemic_size(data, stage = 0.5)
#> [1] 81921.9

# alternatively, get the epidemic size at `time = 50`
epidemic_size(data, time = 50)
#> epidemic_size(): `time` provided will override any `stage` provided
#> [1] 81921.9
```
