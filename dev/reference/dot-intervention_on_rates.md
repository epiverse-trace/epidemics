# Apply interventions to rate parameters

Apply interventions to rate parameters

## Usage

``` r
.intervention_on_rates(t, interventions, parameters)
```

## Arguments

- t:

  A single number for the simulation time.

- interventions:

  A named list of `list`-like objects that each have at least the three
  members `"time_begin"`, `"time_end"`, and `"reduction"`. These are
  used to calculate the effect on each of the named parameters in the
  simulation.

- parameters:

  A named list of numeric parameters affected by `interventions`. This
  represents the model parameters, such as the transmission rate,
  \\\beta\\, or the recovery rate, \\\gamma\\.

## Value

A named list of the same length as `parameters`, with the same names.
These parameters can then be used in a timestep of an ODE model.
