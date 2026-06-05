# Create a population object combining several populations

Create a population object combining several populations

## Usage

``` r
combine_populations(
  populations,
  connectivity_matrix,
  method = "linear",
  name = NA_character_
)
```

## Arguments

- populations:

  A list of `<population>` objects

- connectivity_matrix:

  A numeric matrix for the contact matrix between the elements of the
  populations

- method:

  The method to combine the contact matrices of `populations`, can be a
  character chain (`linear` or `gravity`) which will call an internal
  function, or a user-defined function with two arguments (`populations`
  and `connectivity_matrix`) and returning a numeric matrix.

- name:

  Optional string for the combined population name.

## Value

An object of the `<population>` class.

## Examples

``` r
pop1 <- population(
  name = "Population 1",
  contact_matrix = matrix(c(1, .5, .5, 1), nrow = 2),
  demography_vector = c("0-20" = 2e7, "20+" = 4e7),
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0,
      0.9999, 0.0001, 0, 0),
    nrow = 2, ncol = 4, byrow = TRUE
  )
)

pop2 <- population(
  name = "Population 2",
  contact_matrix = matrix(c(1, .5, .5, 1), nrow = 2),
  demography_vector = c("0-20" = 1e7, "20+" = 2e7),
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0,
      0.9999, 0.0001, 0, 0),
    nrow = 2, ncol = 4, byrow = TRUE
  )
)

prop_matrix <- matrix(c(1, .1, .05, 1), nrow = 2)

combined_population <- combine_populations(
      populations = list(pop1, pop2),
      connectivity_matrix = prop_matrix,
      method = "linear",
      name = "combine"
)


# check for class <population>
is_population(combined_population)
#> [1] TRUE
```
