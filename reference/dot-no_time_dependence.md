# Specify no time dependence of model rates

Specify no time dependence of model rates

## Usage

``` r
.no_time_dependence(func_target = "transmission_rate")
```

## Arguments

- func_target:

  A string for the parameter to be targeted. Defaults to the
  transmission rate parameter.

## Value

A list with a single named element, `transmission_rate`, which is a
function that returns its second argument. This matches the
specification for functions to be passed as the `time_dependence`
argument of epidemic functions.
