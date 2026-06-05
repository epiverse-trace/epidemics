# Model an Ebola virus disease epidemic

Simulate an epidemic using a discrete-time, stochastic SEIR
compartmental model with compartments based on Li et al. (2019), and
with Erlang passage times based on a model developed by Getz and
Dougherty (2017), developed to model the West African Ebola virus
disease (EVD) outbreak of 2013 – 2016. See **Details** for more
information.

`model_ebola_cpp()` is an Rcpp implementation of this model that
currently lags behind the R implementation, and is likely to be removed.

## Usage

``` r
model_ebola(
  population,
  transmission_rate = 1.5/12,
  erlang_subcompartments = 2,
  infectiousness_rate = erlang_subcompartments/5,
  removal_rate = erlang_subcompartments/12,
  prop_community = 0.9,
  etu_risk = 0.7,
  funeral_risk = 0.5,
  intervention = NULL,
  time_dependence = NULL,
  time_end = 100,
  replicates = 100
)
```

## Arguments

- population:

  An object of the `<population>` class, see
  [`population()`](https://epiverse-trace.github.io/epidemics/dev/reference/population.md).

  This model only accepts a `<population>` without demographic
  structure, that is, the `demography_vector` must be a single number
  representing the total size of the affected population.

  The model also does not account for demographic differences in social
  contacts, which means that the `contact_matrix` is ignored. For
  consistency, the matrix must be square and have as many rows as
  demography groups, which is one.

- transmission_rate:

  A numeric vector for the rate at which individuals move from the
  susceptible to the exposed compartment upon contact with an infectious
  individual. Often denoted as \\\beta\\, with \\\beta = R_0 /
  \text{infectious period}\\. See **Details** for default values.

- erlang_subcompartments:

  A numeric, integer-like vector for the number of Erlang
  sub-compartments assumed for the exposed, infectious, and hospitalised
  compartments. Defaults to 2.

- infectiousness_rate:

  A numeric vector for the rate at which individuals move from the
  exposed to the infectious compartment. Often denoted as \\\sigma\\,
  with \\\sigma = 1.0 / \text{pre-infectious period}\\. This value does
  not depend upon the number of infectious individuals in the
  population. See **Details** for default values.

- removal_rate:

  A numeric vector for the rate at which infectious individuals
  transition from the infectious or hospitalised compartments to the
  funeral or removed compartments. This model does not distinguish
  between recoveries and deaths. Denoted in Getz and Dougherty as
  \\\gamma^I\\ (see **Details**).

- prop_community:

  A numeric vector for the proportion of infectious individuals who
  remain in the community and are not hospitalised for treatment.
  Defaults to 0.9.

- etu_risk:

  A numeric vector for the relative risk of onward transmission of EVD
  from hospitalised individuals, with values between 0.0 and 1.0, where
  0.0 indicates that hospitalisation completely prevents onward
  transmission, and 1.0 indicates that hospitalisation does not prevent
  onward transmission at all; values are relative to the baseline
  transmission rate \\\beta\\. Defaults to 0.7.

- funeral_risk:

  A numeric vector for the relative risk of onward transmission of EVD
  from funerals of individuals who died with EVD. Must be between 0.0
  and 1.0, where 0.0 indicates that there is no onward transmission, and
  1.0 indicates that funeral transmission is equivalent to the baseline
  transmission rate in the community \\\beta\\. Defaults to 0.5.

- intervention:

  An optional named list of `<rate_intervention>` objects representing
  optional pharmaceutical or non-pharmaceutical interventions applied to
  the model parameters listed above. May also be a list of such lists,
  in which case each set of interventions is treated as a separate
  scenario. See **Details** below.

- time_dependence:

  An optional named list where each element is a function with the first
  two arguments being the current simulation `time`, and `x`, a value
  that is dependent on `time` (`x` represents a model parameter). List
  names must correspond to model parameters modified by the function.
  Alternatively, may be a list of such lists, in which case each set of
  functions is treated as a distinct scenario. See **Details** for more
  information, as well as the vignette on time- dependence
  `vignette("time_dependence", package = "epidemics")`.

- time_end:

  A numeric, integer-like vector for the maximum number of

- replicates:

  A single number for replicates to run. Defaults to 100. timesteps over
  which to run the model, in days. Defaults to 100 days.

## Value

A `<data.table>`. If the model parameters and composable elements are
all scalars, a single `<data.table>` with the columns "time",
"compartment", "age_group", and "value", giving the number of
individuals per demographic group in each compartment at each timestep
in long (or "tidy") format is returned.

If the model parameters or composable elements are lists or list-like, a
nested `<data.table>` is returned with a list column "data", which holds
the compartmental values described above. Other columns hold parameters
and composable elements relating to the model run. Columns "scenario"
and "param_set" identify combinations of composable elements
(population, interventions, vaccination regimes), and infection
parameters, respectively.

## Details: Discrete-time Ebola virus disease model

This model has compartments adopted from the consensus model for Ebola
virus disease presented in Li et al. (2019), and with transitions
between epidemiological compartments modelled using Erlang
sub-compartments adapted from Getz and Dougherty (2018); see
**References**.

The R code for this model is adapted from code by Ha Minh Lam and
initially made available on *Epirecipes*
(https://github.com/epirecipes/epicookbook) under the MIT license.

The shape of the Erlang distributions of passage times through the
exposed and infectious compartments (\\k^E\\ and \\k^I\\) are
recommended to be set to 2 as a sensible choice, which is the default
value for the `erlang_sbubcompartments` argument, but can be allowed to
vary (but not independently).

The transition rates between the exposed and infectious, and infectious
and funeral compartments (and also hospitalised to removed),
\\\gamma^E\\ and \\\gamma^I\\ in Getz and Dougherty's notation, are
passed by the user as the `infectiousness_rate` and `removal_rate`
respectively.

Getz and Dougherty's equation (6) gives the relationship between these
parameters and the mean pre-infectious \\\rho^E\\ and infectious
\\\rho^I\\ periods. \$\$\gamma^E = \dfrac{k^E}{\rho^E} =
\dfrac{2}{\rho^E} ~\text{and}~ \gamma^I = \dfrac{k^I}{\rho^I} =
\dfrac{2}{\rho^I}\$\$

In this discrete time model, \\\gamma^E\\ and \\\gamma^I\\ are used to
determine the passage times of newly exposed or infectious individuals
through their respective compartments (thus allowing for variation in
passage times).

### Hospitalisation, funerals, and removal

A proportion, `1.0 - prop_community`, of infectious individuals are
transferred to the hospitalised compartment in each timestep, This
compartment represents Ebola Treatment Units (ETUs), and individuals in
the hospitalised compartment are considered to be infectious but no
longer in the community.

The passage time of individuals in the hospitalised compartment is
similar to that of individuals in the infectious compartment (i.e.,
infectious in the community), which means that an infectious individual
with \\N\\ timesteps before exiting the infectious compartment will exit
the hospitalised compartment in the same time.

Hospitalised individuals can contribute to transmission of Ebola to
susceptibles depending on the value of `etu_risk` which scales the
baseline transmission rate \\\beta\\ for hospitalised individuals.

We assume that deaths in hospital lead to Ebola-safe funerals, and
individuals exiting the hospitalised compartment move to the 'removed'
compartment, which holds both recoveries and deaths.

We assume that deaths outside of hospital lead to funerals that are
potentially unsafe burials, and the `funeral_risk` argument scales the
baseline transmission rate \\\beta\\ for funeral transmission of Ebola
to susceptibles.

Individuals are assumed to spend only a single timestep in the funeral
transmission compartment, before they move into the 'removed'
compartment.

Individuals in the 'removed' compartment do no affect model dynamics.

### Model parameters

The default values are:

- Transmission rate (\\\beta\\, `transmission_rate`): 0.125, resulting
  from an \\R_0\\ = 1.5 and an infectious period of 12 days.

- Infectiousness rate (\\\gamma^E\\, `infectiousness_rate`): 0.4,
  assuming a pre-infectious period of 5 days and two Erlang
  subcompartments.

- Removal rate (\\\gamma^I\\, `recovery_rate`): 0.1667, assuming an
  infectious period of 12 days and two Erlang subcompartments.

### Implementing vaccination

Vaccination cannot currently be implemented in this model as it does not
have a "vaccinated" epidemiological compartment. This prevents the use
of a `<vaccination>` object.

Instead, users can use the `time_dependence` argument to pass a function
that modifies model parameters — specifically, the transmission rate —
in a way that is consistent with the effect of vaccination. An example
is shown in the vignette about this model; run this code to open the
vignette: `vignette("ebola_model", package = "epidemics")`

### Vector inputs

#### Vector parameter inputs

The model infection parameters and the model duration may be passed as
numeric or integer-like vectors (as appropriate to the parameter), to
simulate the effect of parameter uncertainty. All parameter vectors must
be of the same length, or any one parameter vector may have a length \>
1 while all other have a length of 1. In the first case, each i-th
combination of parameters is treated as a parameter set. In the second
case, all single value parameters (scalars) are recycled to the same
length as the non-scalar parameter.

The model is run for \$N\$ stochastic realisations of each parameter
set. Random number seeds are preserved across parameter sets, so that
differences in outcomes in each j-th run are due to differences in
parameters alone.

#### Vector inputs for composable elements

The `intervention` and `time_dependence` arguments also accept
vectorised inputs in the form of lists of intervention and time
dependence sets. Each combination of intervention and time-dependence
sets is treated as a distinct 'scenario', and realisations of each
parameter set are run for each scenario.

## References

Li, S.-L., Ferrari, M. J., Bjørnstad, O. N., Runge, M. C., Fonnesbeck,
C. J., Tildesley, M. J., Pannell, D., & Shea, K. (2019). Concurrent
assessment of epidemiological and operational uncertainties for optimal
outbreak control: Ebola as a case study. Proceedings of the Royal
Society B: Biological Sciences, 286(1905), 20190774.
[doi:10.1098/rspb.2019.0774](https://doi.org/10.1098/rspb.2019.0774)

Getz, W. M., & Dougherty, E. R. (2018). Discrete stochastic analogs of
Erlang epidemic models. Journal of Biological Dynamics, 12(1), 16–38.
[doi:10.1080/17513758.2017.1401677](https://doi.org/10.1080/17513758.2017.1401677)

## Examples

``` r
# create a population with 6 compartments
population <- population(
  contact_matrix = matrix(1),
  demography_vector = 14e6,
  initial_conditions = matrix(
    c(0.999998, 0.000001, 0.000001, 0, 0, 0),
    nrow = 1, ncol = 6L
  )
)

# run epidemic simulation with no vaccination or intervention
data <- model_ebola(
  population = population
)

# view some data
head(data)
#>     time demography_group  compartment    value replicate
#>    <int>           <char>       <char>    <num>     <int>
#> 1:     0     demo_group_1  susceptible 13999972         1
#> 2:     0     demo_group_1      exposed       14         1
#> 3:     0     demo_group_1   infectious       14         1
#> 4:     0     demo_group_1 hospitalised        0         1
#> 5:     0     demo_group_1      funeral        0         1
#> 6:     0     demo_group_1      removed        0         1
```
