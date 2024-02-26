# epidemics 0.1.0

This is an initial GitHub release of _epidemics_, an R package that ships a library of compartmental epidemic model structures that can be used, along with supplied classes that help define population characteristics and epidemic response interventions including vaccinations, to compose and model epidemic scenarios.

_epidemics_ is still being actively developed, with major changes planned for the near future. This release is aimed at supporting reproducibility projects that used _epidemics_ which would be subject to breaking changes due to planned package development.
The sections below describe the contents of this release.

## Model structures

This release of _epidemics_ includes four model structures supporting a range of composable elements to modify epidemic trajectories.

1. "Default" model: A deterministic SEIR-V model allowing heterogeneity in social contacts between demographic groups, with optional, single-dose non-leaky vaccination;

2. "Vacamole" model: A deterministic SEI-HRD-V2 implementation of a model allowing heterogeneity in social contacts between demographic groups, with a two-dose leaky vaccination, supporting different infection trajectories through the infectious and hospitalised (H) compartments for doubly vaccinated individuals, which tracks deaths (D), and which was initially developed by the Dutch public health agency RIVM for vaccine impact modelling during the Covid-19 pandemic, and published as Ainslie et al. 2022 <https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090>;

3. "Diphtheria" model: A deterministic SEIHR model tracking outcomes for different demographic groups, but not including heterogeneity in social contacts, adapted from Finger et al. 2019 <https://doi.org/10.1186/s12916-019-1288-7> and intended for application to disease outbreaks in a humanitarian camp setting;

4. "Ebola" model: A discrete time stochastic SEIHFR model suitable for modelling Ebola virus disease and other haemorrhagic fevers, and which allows varying the efficacy of isolation in a hospital setting (H), and allows modelling transmission in a funeral context (F), as adapted from a consensus Ebola virus disease model in Li et al. 2019 <https://doi.org/10.1098/rspb.2019.0774> and using simulation methods from Getz and Dougherty 2018 <https://doi.org/10.1080/17513758.2017.1401677>.

## Solving ODE systems using Boost _odeint_

_epidemics_ uses Boost's _odeint_ <https://www.boost.org/doc/libs/1_84_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html> to treat the deterministic models' ordinary differential equations (ODEs) as initial value problems and solve them.

Model ODEs are defined as structs in the package headers, and exposed to R as internal Rcpp functions. The 'default', 'Vacamole', and 'diphtheria' models are ODE models defined in this way. This is intended to help reduce overheads associated with passing ODE systems written in R back and forth from a solver (such as those provided by {deSolve}), and is an easier way to define feature-rich models than writing C code for solvers provided by {deSolve} that accept compiled code.

_epidemics_ headers include tools for handling the C++ representations of R objects used in the package (see below), and can be imported by other Rcpp packages.

The 'default' and 'Vacamole' models have equivalent R-only implementations as well which use the {deSolve} package; these are intended to be made unavailable in future releases.

## Composable elements as classes

_epidemics_ provides classes that help to organise the components of an epidemic scenario model.

1. `<population>`: An S3 class to store population characteristics including the size of demographic groups, a social contacts matrix, and initial conditions for a model;

2. `<intervention>`: An S3 abstract class and super-class that allows the definition of events that modify the epidemic trajectory:

   a. `<rate_intervention>`: A sub-class of `<intervention>` that allows the reduction of transition rates between model compartments to simulate the effect of policy interventions over a specific period;

   b. `<contacts_intervention>`: A sub-class of `<intervention>` that allows the reduction of social contacts to simulate the effect of policy interventions over a specific period;

3. `<vaccination>`: An S3 class that holds the intervals and group-specific rates at which individuals transition into the 'vaccinated' compartment(s) of a model, if available;

## Other composable elements

_epidemics_ allows models to include elements that affect an epidemic trajectory, but which are not custom classes.

1. Time-depedence: All models can be passed a list of functions with two arguments, `time` and `x` which are expected to return `x` as a function of `time`, and which may be used to model the effect of seasonality in model parameters;

2. Population changes: Applicable only to the diphtheria model, a two element list of `time` and `values`, which allow the definition of changes to the number of susceptible individuals in the model, and which may be used to model influxes and evacuations of individuals from humanitarian camps.

## Output processing functions

_epidemics_ provides functions to help process the output of an epidemic model run, to calculate the size of the epidemic in each demographic group at any stage (`epidemic_size()`), and to calculate the number of new infections in each demographic group at each timepoint in the model (`new_infections()`).

## Usage vignettes

_epidemics_ includes a range of usage vignettes that demonstrate how to:

1. Get started with the package;

2. Get started with modelling interventions on social contacts to control outbreaks;

3. Model overlapping and sequential interventions on social contacts;

4. Model interventions that modify transition rates between model compartments;

5. Get started with modelling a vaccination campaign;

6. Model time-dependence and seasonality in disease transmission dynamics;

7. Generate and model uncertainty in model parameters;

8. Reduce the number of parameters required for final size estimation;

9. Use the 'Vacamole' model for scenarios of leaky vaccination and vaccine impact on hospitalisation;

10. Use the 'Ebola' model for scenarios of responses to an Ebola virus disease outbreak;

11. Use the 'diptheria' model for scenarios of outbreaks in a humanitarian camp setting.

## Miscellaneous

1. Workflows to render the vignettes and README as a website;

2. Test code coverage of 93%.
