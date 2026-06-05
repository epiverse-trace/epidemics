# Prepare a `<population>` for an epidemic model

Prepare a `<population>` for an epidemic model

## Usage

``` r
.prepare_population(x)
```

## Arguments

- x:

  A `<population>`.

## Value

A named list of "contact_matrix", the population social contacts matrix
scaled by the largest real eigenvalue and the demography vector, and
"initial_state", the proportional initial state multiplied by the
demography vector.
