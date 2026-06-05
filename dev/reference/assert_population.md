# Assert properties of a `population` object

Assert that objects of the `population` class have the parameters
expected by an epidemic model. See
[`population()`](https://epiverse-trace.github.io/epidemics/dev/reference/population.md)
and specific epidemic functions to check the population parameters
required by each model. This function is for internal use in argument
checking functions.

## Usage

``` r
assert_population(x, compartments, demography_vector = NULL)
```

## Arguments

- x:

  A `<population>` object.

- compartments:

  A character vector giving the names of model compartments whose length
  is taken as the reference for the number of columns in the
  `initial_conditions` matrix in `x`.

- demography_vector:

  An optional numeric vector whose length is used to check the length of
  the demography vector present in `x`.

## Value

Silently returns the `<population>` object `x`. Primarily called for its
side effects of throwing errors when `x` does not meet certain
requirements.
