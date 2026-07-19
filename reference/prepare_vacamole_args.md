# Prepare arguments to Vacamole epidemic function

Prepare arguments to
[`model_vacamole()`](https://epiverse-trace.github.io/epidemics/reference/model_vacamole.md).

## Usage

``` r
.check_prepare_args_vacamole(mod_args)
```

## Value

A list of model arguments suitable for model_vacamole(). This is a named
list consisting of:

- `initial_state`: the initial conditions modified to represent absolute
  rather than proportional values;

- `transmission_rate`, `transmission_rate_vax`: two numbers representing
  the transmission rate of the infection for unvaccinated or single-dose
  vaccinated, and two-dose vaccinated individuals, respectively;

- `infectiousness_rate`: a single number for the transition rate from
  the 'exposed' and 'exposed_vaccinated' to the 'infectious' and
  'infectious_vaccinated' compartments;

- `mortality_rate`, `mortality_rate_vax`: two numbers representing the
  mortality rate of the infection for unvaccinated or single-dose
  vaccinated, and two-dose vaccinated individuals, respectively;

- `hospitalisation_rate`, `hospitalisation_rate_vax`: two numbers
  representing the hospitalisation rate of the infection for
  unvaccinated or single-dose vaccinated, and two-dose vaccinated
  individuals, respectively;

- `recovery_rate`: a single number for the recovery rate from the
  infection;

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

`.check_prepare_args_vacamole()` prepares arguments for
model_vacamole(), which is the *odin* based C function that solves the
Vacamole ODE system.
