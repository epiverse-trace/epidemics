# Package index

## Model functions

- [`model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md)
  : Model an SEIR-V epidemic with interventions
- [`model_diphtheria()`](https://epiverse-trace.github.io/epidemics/reference/model_diphtheria.md)
  : Model a diphtheria outbreak using a compartmental ODE model
- [`model_ebola()`](https://epiverse-trace.github.io/epidemics/reference/model_ebola.md)
  : Model an Ebola virus disease epidemic
- [`model_vacamole()`](https://epiverse-trace.github.io/epidemics/reference/model_vacamole.md)
  : Model leaky, two-dose vaccination in an epidemic using Vacamole

## Classes and related functions

### Intervention

- [`as.intervention()`](https://epiverse-trace.github.io/epidemics/reference/as.intervention.md)
  : Convert a list to a intervention object

- [`.intervention_on_rates()`](https://epiverse-trace.github.io/epidemics/reference/dot-intervention_on_rates.md)
  : Apply interventions to rate parameters

- [`intervention()`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  [`is_intervention()`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  [`is_contacts_intervention()`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  [`is_rate_intervention()`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  [`c(`*`<contacts_intervention>`*`)`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  [`c(`*`<rate_intervention>`*`)`](https://epiverse-trace.github.io/epidemics/reference/intervention.md)
  : Create an intervention for an epidemic model

- [`print(`*`<intervention>`*`)`](https://epiverse-trace.github.io/epidemics/reference/print_intervention.md)
  [`print(`*`<contact_intervention>`*`)`](https://epiverse-trace.github.io/epidemics/reference/print_intervention.md)
  [`print(`*`<rate_intervention>`*`)`](https://epiverse-trace.github.io/epidemics/reference/print_intervention.md)
  :

  Print an object of the `<intervention>` super-class

### Population

- [`combine_populations()`](https://epiverse-trace.github.io/epidemics/reference/combine_populations.md)
  : Create a population object combining several populations

- [`population()`](https://epiverse-trace.github.io/epidemics/reference/population.md)
  [`is_population()`](https://epiverse-trace.github.io/epidemics/reference/population.md)
  : Construct a new population for an epidemic model

- [`print(`*`<population>`*`)`](https://epiverse-trace.github.io/epidemics/reference/print.population.md)
  :

  Print a `<population>` object

### Vaccination

- [`as.vaccination()`](https://epiverse-trace.github.io/epidemics/reference/as.vaccination.md)
  : Convert a list to a vaccination object

- [`print(`*`<vaccination>`*`)`](https://epiverse-trace.github.io/epidemics/reference/print.vaccination.md)
  :

  Print a `<vaccination>` object

- [`vaccination()`](https://epiverse-trace.github.io/epidemics/reference/vaccination.md)
  [`is_vaccination()`](https://epiverse-trace.github.io/epidemics/reference/vaccination.md)
  [`c(`*`<vaccination>`*`)`](https://epiverse-trace.github.io/epidemics/reference/vaccination.md)
  : Construct a new vaccination regime for an epidemic model

## Helper functions

Functions to help with handling package classes and outputs.

- [`epidemic_size()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_size.md)
  : Get the epidemic size
- [`epidemic_peak()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_peak.md)
  : Get the time and size of a compartment's highest peak
- [`new_infections()`](https://epiverse-trace.github.io/epidemics/reference/new_infections.md)
  : Get new infections over model time

## Scenario comparison

Functions to help with comparing intervention scenarios.

- [`outcomes_averted()`](https://epiverse-trace.github.io/epidemics/reference/outcomes_averted.md)
  : Calculate outcomes averted by interventions
