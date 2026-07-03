# Ordinary Differential Equations for the Vacamole Model

Provides the ODEs for the RIVM Vacamole model in a format that is
suitable for passing to
[`deSolve::lsoda()`](https://rdrr.io/pkg/deSolve/man/lsoda.html). See
[`model_vacamole()`](https://epiverse-trace.github.io/epidemics/reference/model_vacamole.md)
for a list of required parameters.

## Usage

``` r
.ode_model_vacamole(t, y, params)
```

## Arguments

- t:

  A single number of the timestep at which to integrate.

- y:

  The conditions of the epidemiological compartments.

- params:

  The parameters, passed as a named list.

## Value

A list with a vector with as many elements as the number of demographic
groups times the number of epidemiological compartments. Each value
gives the change in the number of individuals in that compartment.
