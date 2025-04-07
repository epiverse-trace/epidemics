# epidemics (development version)

## Helper functions

1. Added the `combine_populations()` function to combine age-structured populations into a new population object.

## Documentation

1. Added a `modelling-populations` vignette on how to combine age-structured populations.

# epidemics 0.4.0

Maintainer is changing to @rozeggo.

## Model functions

1. Internal model functions for the models which allow vaccination have been corrected to prevent vaccination introducing negative values of susceptibles; tests added to check for this (#235, initially reported by @avallecam).

## Helper functions

1. Added the `epidemic_peak()` function to calculate the timing and size of the largest peak in each compartment in an scenario model (#240) by @bahadzie.

2. Added the `outcomes_averted()` function to compare epidemic scenarios (e.g. with and without interventions or vaccination) (#225, #230).

## Documentation

1. Adds a developer-focused vignette on how to modify epidemics and model structures to address potential modelling requests or tasks (#210).

2. Splits up the 'Modelling uncertainty and scenarios' vignette into separate vignettes on uncertainty and scenario comparisons (#225).

3. Removed unnecessary plots from the vignette on modelling vaccination (#235).

4. Fixed link to _socialmixr_ package in the 'Get started' and 'Modelling interventions' vignettes.

5. Updated and added documentation for all new or modified functions.

6. Updated references JSON file.

## Package

1. Updated Codecov GitHub Actions workflow to restore code coverage reporting.

2. Updated package title and citation file. 

3. Updated `_pkgdown.yaml` with new vignette and updated section titles.

4. Updated WORDLIST.

# epidemics 0.3.0

This is a minor version release of _epidemics_ is the end point of this project that allows [vector inputs to the Ebola model](https://github.com/orgs/epiverse-trace/projects/36) (#211, #212).

## Breaking changes

1. The Ebola model run with a single parameter set and a single set of composable (i.e., functionality before this vignette) now runs 100 replicates of the model by default, and the output returns the additional column name "replicate" (#211).

2. The Ebola model disallows rate interventions on, and time-dependence for, the 'infectiousness rate' and the 'removal rate' as these control the number of epidemiological sub-compartments and cannot be allowed to vary during run time.

1. The default behaviour of `epidemic_size()` is to exclude the 'dead' compartment from epidemic size calculations; this has changed from including it by default, as most models don't have a 'dead' compartment (#212);

## Model functions

1. The Ebola model code has been split into a two-level structure similar to the ODE models. The user facing function `model_ebola()` now handles input checking and some cross-checking, and makes combinations of infection parameter sets and scenarios, and handles output. The internal model function `.model_ebola_internal()` is now called over combinations of parameters and scenario elements (#211).

2. `.model_ebola_internal()` relies on functionality from the {withr} package to preserve the random number seed between parameter set-scenario combinations, and to ensure that each $i$-th replicate of each scenario starts with the same random number stream.

## Model structures

## Classes

1. The `<intervention>` helper functions `.cumulative_rate_intervention()` and `.intervention_on_rates()` have been renamed to include the prefix `.` to indicate they are internal functions (#211).

2. Small updates to `.cross_check_intervention()` to better handle intervention sets where contacts interventions are optional (the Ebola model) (#211).

## Helper functions

1. `.check_prepare_args_ebola()` is a new argument cross-checking and preparation function for the Ebola model (#212).

2. `epidemic_size()` is substantially updated (#212):

   - Added option for `time` which returns epidemic size at a specific time point, overriding the `stage` argument, defaults to `NULL` as the intended use of the function is to return the final size;

   - Added option to return epidemic sizes at multiple stages or time points (`stage` and `time` can be vectors);

   - Added option to simplify the output to a vector, which is `TRUE` by default to keep consistency with previous functionality;

   - Added functionality to handle replicates from the Ebola model;

   - Added tests for new functionality.

3. Added function `.output_to_df_ebola()` to handle output from the Ebola model (#211).

## Documentation

1. Updates to the Ebola model vignette showing new functionality.

2. Updates to the design decisions vignette documenting decisions taken for the Ebola model.

3. Updated sections on handling social contacts data in the 'Getting started' and 'Modelling interventions on social contacts' vignettes, added reference to Wallinga et al. 2006 <https://doi.org/10.1093/aje/kwj317> (#217).

4. Enables auto development mode for the website in `_pkgdown.yml`, and manually specify Bootstrap version 5 (#213).

5. Updates WORDLIST.

6. Updates installation instructions in the Readme; link under 'Contribute' correctly directs to pull requests page (#217).

## Package

1. {withr} moved from Suggests to Imports due to use in seed management.

2. Added @bahadzie as contributor and @jamesmbaazam as reviewer.

3. Updates to the GitHub Actions workflows and linter config file.

4. Updates the license year to 2024 in all files with the year.

# epidemics 0.2.0

This is a second GitHub release of _epidemics_ which makes substantial additions to the functionality in v0.1.0, and introduces significant breaking changes (#176).
This release is the end point of [this project to ship vectorised ODE models.](https://github.com/orgs/epiverse-trace/projects/35/)

This release focuses on the ODE models in _epidemics_.

## Breaking changes

1. All model functions have been renamed to `model_<NAME>()`, removing the language suffix (#176).

2. The wrappers around R-only implementations of the 'default' and 'Vacamole' models have been removed, but the ODE system functions have been retained for potential future use (#176).

3. The "Vacamole" model has been refactored with the arguments `*_reduction_vax` for the effect of double vaccination on compartmental transition rates renamed to `*_vax`, where `*` may be one of "susceptibility", "hospitalisation", and "mortality". Previously, these parameters were implemented in the internal C++ code, but presented as inverse values to users (i.e., `susceptibility_vax = susceptibility * (1 - (susceptibility_reduction_vax)))`). This change brings the user-facing representation in line with the internal implementation, and allows these parameters to be targeted by rate interventions and time-dependence, which was not possible earlier (#176, #203, and overriding #144).

4. The function `get_parameter()` has been removed (#176).

5. The infection parameter "transmissibility" has been renamed to "transmission rate" and the corresponding function argument is `transmission_rate` (#196).

## Model functions

- ODE model functions can now be passed numeric vectors of infection parameters, and lists of intervention sets (a list of `<intervention>`s), and lists of `<vaccination>` to the `intervention` and `vaccination` argument respectively. This is the **main change in this minor version.** (#176)

- ODE model functions take on {data.table} to make combinations of interventions and vaccinations (together called a 'scenario'), and infection parameter sets to run each scenario for each parameter set (#176).

- ODE model functions all return `<data.table>`s - these may be nested with the model arguments as identity columns if vectors of parameters or multiple scenarios are passed. A simple, unnested `<data.table>` is passed if the model function is called with scalar arguments (#176).

## Model structures

- There is a fix to how vaccination is implemented in the default and Vacamole models: the group-specific vaccination rates passed in a `<vaccination>` are internally converted to a count by multiplication with the group-specific population size, and this value is subtracted from any susceptibles, while vaccination is active. The previous implementation made the number of vaccinations dependent on the number of susceptibles, which was not in line with a public health understanding of vaccination (details: #198, fix: #202).

- All C++ ODE model implementation now access values from the `std::unordered_map` of model parameters using the `at` operator on keys (`model_params.at("parameter_name")`) rather than using `[` (`model_params["parameter_name"]`); the latter method will introduce a key-value pair of `parameter_name : 0` for numeric value types when the key is missing from the map. This led to issues with the Vacamole model where the `*_vax` parameters were not correctly transferred to C++ and were substituted with 0s without throwing an error (#203).

## Classes

No substantial changes to classes; small additions of input checking to `<population>` class.

## Helper functions

1. Internal helper functions `.check_args_model_*()` and `.prepare_args_model_*()` have been combined into single functions `.check_prepare_args_*()` that are called both for their output and for input checking side effects.

2. The new internal functions `.prepare_population()`, `.cross_check_intervention()`, `.cross_check_vaccination()`, `.cross_check_timedep()`, and `.cross_check_popchange()` check and prepare a model population and check other inputs for compatibility with it. These are used in `.check_prepare_args_*()`.

3. New internal functions have been added to check and recycle lists of vectors; original implementations by @TimTaylor.

4. `output_to_df()` is renamed to `.output_to_df()`.

## Documentation

1. The benchmarking vignette has been removed as the R-only model implementations are no longer provided to users (#176).

2. The vignette on parameter uncertainty has been rewritten to show how to pass vectors of infection parameters and model composable elements to model functions, and renamed to "Modelling parameter uncertainty and epidemic scenarios" (#176).

3. A design decisions vignette has been added to help developers and contributors understand the package architecture and design choices made in _epidemics_ development; this includes a conceptual design diagram as well as an architecture diagram (#188).

4. All function documentation has been updated to reflect name changes and other minor improvements (#176).

5. The README has been reorganised to shift the model list to the end, to make the installation instructions easier to find (#176).

6. Corrected the website URL in `_pkgdown.yml`; this allows search functionality in the package website (#195).

7. Updated WORDLIST.

## Package

1. All ODE model functions have received a more extensive and more standardised (as much as possible) test suite.

2. Filenames have been standardised to show which files are related, e.g. `R/model_default.R`, `src/model_default.cpp`, and `inst/include/model_default.h`; references to filenames such as in the package header have been updated.

3. Removed {deSolve} from dependencies.

4. Added {ggdist}, {withr}, {purr} and {tidyr} to Suggests; `CITATION.cff` updated to match.

5. Added basic infrastructure for continuous relative benchmarking (#206).

# epidemics 0.1.0

This is an initial GitHub release of _epidemics_, an R package that ships a library of compartmental epidemic model structures that can be used, along with supplied classes that help define population characteristics and epidemic response interventions including vaccinations, to compose and model epidemic scenarios.

_epidemics_ is still being actively developed, with major changes planned for the near future. This release is aimed at supporting the reproducibility of projects that used _epidemics_ which would be subject to breaking changes due to planned package development.
The sections below describe the contents of this release.

## Model structures

This release of _epidemics_ includes four model structures supporting a range of composable elements to modify epidemic trajectories.

1. "Default" model: A deterministic SEIR-V model allowing heterogeneity in social contacts between demographic groups, with optional, single-dose non-leaky vaccination;

2. "Vacamole" model: A deterministic SEI-HRD-V2 implementation of a model allowing heterogeneity in social contacts between demographic groups, with a two-dose leaky vaccination, supporting different infection trajectories through the infectious and hospitalised (H) compartments for doubly vaccinated individuals, which tracks deaths (D), and which was initially developed by the Dutch public health agency RIVM for vaccine impact modelling during the Covid-19 pandemic, and published as Ainslie et al. 2022 <https://doi.org/10.2807/1560-7917.ES.2022.27.44.2101090>;

3. "Diphtheria" model: A deterministic SEIHR model tracking outcomes for different demographic groups, but not including heterogeneity in social contacts, adapted from Finger et al. 2019 <https://doi.org/10.1186/s12916-019-1288-7> and intended for application to disease outbreaks in a humanitarian camp setting;

4. "Ebola" model: A discrete time stochastic SEIHFR model suitable for modelling Ebola virus disease and other haemorrhagic fevers, and which allows varying the efficacy of isolation in a hospital setting (H), and allows modelling transmission in a funeral context (F), as adapted from a consensus Ebola virus disease model in Li et al. 2019 <https://doi.org/10.1098/rspb.2019.0774> and using simulation methods from Getz and Dougherty 2018 <https://doi.org/10.1080/17513758.2017.1401677>.

## Solving ODE systems using Boost _odeint_

_epidemics_ uses Boost's _odeint_ <https://www.boost.org/doc/libs/1_84_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html> to treat the deterministic models' ordinary differential equations (ODEs) as initial value problems and solve them.

Model ODEs are defined as `structs` with operators in the package headers, and exposed to R as internal Rcpp functions. The 'default', 'Vacamole', and 'diphtheria' models are ODE models defined in this way. This is intended to help reduce overheads associated with passing ODE systems written in R back and forth from a solver (such as those provided by {deSolve}), and is an easier way to define feature-rich models than writing C code for solvers provided by {deSolve} that accept compiled code.

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

1. Time-dependence: All models can be passed a list of functions with two arguments, `time` and `x` which are expected to return `x` as a function of `time`, and which may be used to model the effect of seasonality in model parameters;

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

11. Use the 'diphtheria' model for scenarios of outbreaks in a humanitarian camp setting.

## Miscellaneous

1. Workflows to render the vignettes and README as a website;

2. Test code coverage of 93%.
