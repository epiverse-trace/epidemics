# Prepare arguments to diphtheria model function

Prepare arguments to
[`model_diphtheria()`](https://epiverse-trace.github.io/epidemics/reference/model_diphtheria.md).

## Usage

``` r
.check_prepare_args_diphtheria(mod_args)
```

## Arguments

- mod_args:

  A named list of the population, and epidemic modifiers.

## Value

A list of model arguments suitable for model_diphtheria(). This is a
named list consisting of:

- `initial_state`: the initial conditions modified to represent absolute
  rather than proportional values;

- `transmission_rate`, `transmission_rate_vax`: two numbers representing
  the transmission rate of the infection for unvaccinated or single-dose
  vaccinated, and two-dose vaccinated individuals, respectively;

- `infectiousness_rate`: a single number for the transition rate from
  the 'exposed' and 'exposed_vaccinated' to the 'infectious' and
  'infectious_vaccinated' compartments;

- `recovery_rate`: a single number for the recovery rate from the
  infection;

- `reporting_rate`: a single number for the proportion of infectious
  cases reported;

- `prop_hosp`: a single number for the proportion of reported cases that
  need hospitalisation;

- `hosp_entry_rate`, `hosp_exit_rate`: two numbers representing the rate
  of entry and exit from the 'hospitalised' compartment;

- `rate_interventions`: a List giving the interventions on model
  parameters;

- `time_dependence`: a List giving the time-dependent effects on model
  parameters in the form of R functions;

- `pop_change_times` and `pop_change_values`: the times and values of
  changes in the population of susceptibles;

- `time_end`, `increment`: two numbers for the time at which to end the
  simulation, and the value by which the simulation time is incremented.

## Details

`.check_prepare_args_diphtheria()` prepares arguments for
model_diphtheria(), which uses an *odin* C function to solve the ODE
system.
