# Prepare arguments to default model function

Prepare arguments to
[`model_default()`](https://epiverse-trace.github.io/epidemics/dev/reference/model_default.md).

## Usage

``` r
.check_prepare_args_default(mod_args)
```

## Arguments

- mod_args:

  A named list of the population, and epidemic modifiers.

## Value

A list of model arguments suitable for seirv_model(). This is a named
list consisting of:

- `initial_state`: the initial conditions modified to represent absolute
  rather than proportional values;

- `transmission_rate`, `infectiousness_rate`, `recovery_rate`: three
  numbers representing the transmission rate of the infection, the rate
  of transition from exposed to infectious, and the recovery rate,
  respectively;

- `contact_matrix`, a numeric matrix for the population contact matrix
  scaled by the largest real eigenvalue and by the size of each groups;

- `npi_time_begin`, `npi_time_end`: two vectors for the start and end
  times of any interventions applied;

- `npi_cr`: a matrix for the age- and intervention-specific effect on
  social contacts;

- `vax_time_begin`,`vax_time_end`, `vax_nu`: three numeric matrices for
  the age- and dose-specific start times, end times, and rates of any
  vaccination doses implemented;

- `time_end`, `increment`: two numbers for the time at which to end the
  simulation, and the value by which the simulation time is incremented.

## Details

`.check_prepare_args_default()` prepares arguments for seirv_model(),
which is the C function that solves the default ODE system using odin.
