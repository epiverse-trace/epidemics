# Cross-check model elements

Check model elements for compatibility with the population in an
epidemic model, returning compatible dummy values when model elements
are not applied, and erroring appropriately when model elements are not
compatible with the population characteristics.

## Usage

``` r
.cross_check_intervention(x, population, allowed_targets)

.cross_check_vaccination(x, population, doses)

.cross_check_timedep(x, allowed_targets)

.cross_check_popchange(x, population)
```

## Arguments

- x:

  Model input to be checked. The expected value of `x` depends on the
  function:

  - `.cross_check_intervention()`: A named list of `<intervention>`
    objects;

  - `.cross_check_vaccination()`: A `<vaccination>` object;

  - `.cross_check_timedep()`: A named list of functions with two
    arguments, `time` and `x`, typically returning `x` as a function of
    `time`;

  - `.cross_check_popchange()`: A named list with two elements, `time`
    and `values`, describing the times and values by which the number of
    susceptibles changes in an epidemic model.

- population:

  An object of the `population` class, which holds a population contact
  matrix, a demography vector, and the initial conditions of each
  demographic group. See
  [`population()`](https://epiverse-trace.github.io/epidemics/reference/population.md).

- allowed_targets:

  The model components, or infection parameters, that the model input
  `x` affects.

- doses:

  The expected number of vaccination doses.

## Value

- `.cross_check_intervention()` returns a named list with at least the
  elements "contacts" describing a `<contacts_intervention>` on
  `population` (if this is among the allowed targets), and a
  `<rate_intervention>` on the transmission rate parameter. If these are
  present in `x`, they are returned as is, or substituted if missing.
  Any other interventions are also returned. If `x` is `NULL`, dummy
  contact and rate interventions are returned in a list.

- `.cross_check_vaccination()` returns `x` after checking that it is
  suitable for `population`, or a dummy vaccination regime with `doses`
  number of doses for each age group.

- `.cross_check_timedep()` returns `x` if `x` is not `NULL`, otherwise
  returns a dummy function operating on the transmission rate parameter
  by default; see
  [`.no_time_dependence()`](https://epiverse-trace.github.io/epidemics/reference/dot-no_time_dependence.md);

- `.cross_check_popchange()` returns `x` after checks against
  `population` if `x` is not `NULL`, otherwise returns a dummy list with
  no population change; see
  [`.no_population_change()`](https://epiverse-trace.github.io/epidemics/reference/dot-no_population_change.md).
