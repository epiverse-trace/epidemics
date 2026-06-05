# Assert properties of a `vaccination` object

Assert that objects of the `vaccination` class have the parameters
expected by an epidemic model. See
[`vaccination()`](https://epiverse-trace.github.io/epidemics/dev/reference/vaccination.md)
and specific epidemic functions to check the vaccination properties
required by each model. This function is for internal use in argument
checking functions.

## Usage

``` r
assert_vaccination(x, doses, population = NULL)
```

## Arguments

- x:

  A
  [vaccination](https://epiverse-trace.github.io/epidemics/dev/reference/vaccination.md)
  object.

- doses:

  The number of doses expected in the vaccination object.

- population:

  An optional argument which is a `<population>` object. When present,
  this is used to check whether the vaccination object `x` has
  corresponding values for each demographic group in `population`.

## Value

Silently returns the `<vaccination>` object `x`. Primarily called for
its side effects of throwing errors when `x` does not meet certain
requirements.
