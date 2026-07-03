# Reducing parameters required for final size estimation

**New to *epidemics*?** It may help to read the [“Get
started”](https://epiverse-trace.github.io/epidemics/articles/epidemics.md)
vignette first!

Additionally, read how to quickly [calculate the final size of an
epidemic](https://epiverse-trace.github.io/finalsize/) - with fewer
epidemiological inputs than *epidemics* required - the *finalsize*
package, which uses an analytical solution to avoid simulating the full
epidemic dynamics.

The total number of individuals expected to have been infected over an
epidemic is called the epidemic final size ([Miller
2012](#ref-miller2012)). The R package
[*finalsize*](https://cran.r-project.org/package=finalsize) can help to
calculate the epidemic final sizes using the final size equation
([Miller 2012](#ref-miller2012)), while accounting for heterogeneity in
population contacts as well as susceptibility to infection, such as due
to vaccination.

*epidemics* can also be used to calculate the final size of an outbreak
while accounting for heterogeneity in contacts and the rollout of
vaccination campaigns.

This vignette lays out which of the two is more suitable for specific
tasks, and how to convert between when seeking more insight into the
effect of policy decisions made during epidemic response.

``` r

library(epidemics)
library(finalsize)

library(dplyr)
library(tibble)
library(tidyr)
library(bench)
```

## Different use cases of *finalsize* and *epidemics*

While both *finalsize* and *epidemics* can be used to calculate the
epidemic final size, they are aimed at different use cases.

***finalsize*** assumes that both infection and population
characteristics are fixed to their initial conditions. For the
infection, this includes properties such as the transmission rate, while
for the population, it includes social contacts between demographic
groups, and the proportion of demographic groups that have specific
levels of susceptibility to infection.

An advantage is that *finalsize* only requires that we define the
initial transmission rate of the infection and susceptibility of the
population, rather than all the time-dependent processes that drive the
epidemic shape, such as the duration of infectiousness or latent period
of the infection. For questions relating to epidemic shape, rather than
shape, this results in a much simpler set of inputs.

However, it also means that *finalsize* cannot be used to model temporal
dynamics of epidemic response, or be used to answer policy questions
with a temporal component, such as when to implement interventions.

***epidemics*** includes a number of scenario models, each with its own
assumptions. Most models allow for some modification of the initial
characteristics of the outbreak due to a range of events. For the
infection, this includes interventions (such as masking or treatments)
that reduce the number of forward transmissions or deaths, but also
seasonal effects which may increase or decrease the transmission rate.
For the population, initial conditions of social contacts can be
influenced by interventions as well.

This makes it much easier to model the temporal dynamics of
public-health policy decisions which are taken during epidemic response
as *epidemics* has many more features that allow for such modelling.

However, it requires more inputs to be defined, including time-dependent
infections processes. It can also be more difficult to model scenarios
in which more complicated susceptibility structure is required, such as
when some demographic groups have underlying immunity to infection due
to past exposure or vaccination. Thus *epidemics* is likely to be
especially useful for outbreaks of novel pathogens (such as the Covid-19
pandemic) where there is little population immunity to infection.

It is easier to configure *finalsize* out-of-the-box for scenarios with
complex demographic patterns of underlying susceptibility to (or
immunity against) infection, potentially due to a history of previous
outbreaks and the policy responses (such as vaccination). Thus, while it
cannot model temporal dynamics, it can quickly provide useful initial
estimates of the final size of outbreaks, without having to write
compartmental models which implement multiple policy decisions.

## Converting scenarios between *finalsize* and *epidemics*

Here, we show an example in which we show how to model a similar
scenario using both *finalsize* and *epidemics*. For example, we might
want to study the full epidemic dynamics in future, but start with a
simpler final size estimate in the meantime. Here we show how to build
equivalent models, while allowing extensions to model epidemic temporal
dynamics (using *epidemics*) later.

As an illustration, we use *epidemics* to model the effect of
vaccination on the trajectory of an epidemic. While *finalsize* does not
allow for the implementation of a *dynamic* vaccination calendar, we can
model reduced susceptibility to infection in the population, such as due
to prior vaccination.

The two methods are comparable if we model vaccination in *epidemics* as
occurring before the main wave of the epidemic, as this sets up
underlying susceptibility.

### Prepare population and model parameters

We first prepare the population (modelled on the U.K.), initial
conditions, and model parameters, before passing them to
[`epidemics::model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md).
This example does not include any interventions or vaccination regimes.

Code to prepare model inputs is folded below for brevity. See the [“Get
started”](https://epiverse-trace.github.io/epidemics/articles/epidemics.md)
vignette for an explanation of how to prepare these inputs.

``` r

# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age_limits = c(0, 20, 40),
  symmetric = TRUE,
  return_demography = TRUE
)

# prepare contact matrix
contact_matrix <- t(contact_data$matrix)

# prepare the demography vector
demography_vector <- contact_data$demography$population
names(demography_vector) <- rownames(contact_matrix)

# initial conditions: one in every 1 million is infected
initial_i <- 1e-6
initial_conditions <- c(
  S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
)

# build for all age groups
initial_conditions <- rbind(
  initial_conditions,
  initial_conditions,
  initial_conditions
)
rownames(initial_conditions) <- rownames(contact_matrix)
```

``` r

# prepare the population to model as affected by the epidemic
# the contact_matrix, demography_vector, and initial_conditions
# have been prepared in the folded code above
uk_population <- population(
  name = "UK",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

# view the population
uk_population
#> 
#>  Population name: 
#> 
#>  Demography 
#> [0,20): 14,799,290 (20%)
#> [20,40): 16,526,302 (30%)
#> [40,Inf): 28,961,159 (50%)
#> 
#>  Contact matrix 
#>                  age.group
#> contact.age.group   [0,20)  [20,40) [40,Inf)
#>          [0,20)   7.883663 2.794154 1.565665
#>          [20,40)  3.120220 4.854839 2.624868
#>          [40,Inf) 3.063895 4.599893 5.005571
#> 
#>  Initial Conditions 
#>                 S E     I R V
#> [0,20)   0.999999 0 1e-06 0 0
#> [20,40)  0.999999 0 1e-06 0 0
#> [40,Inf) 0.999999 0 1e-06 0 0
```

We model the spread of influenza with pandemic potential, assuming \$R_0
= \$ 1.5, a pre-infectious period of 3 days, and an infectious period of
7 days.

This leads to the following model parameters for transmission rate
$`\beta`$, the rate of transition from exposed to infectious $`\sigma`$,
and the recovery rate $`\gamma`$.

``` r

# simulate pandemic parameters
transmission_rate <- 1.5 / 7
infectiousness_rate <- 1 / 3
recovery_rate <- 1 / 7
```

### Implementing vaccination in *epidemics*

For simplicity, we implement a vaccination regime in which vaccines are
delivered to individuals over the age of 40 over the course of roughly
six months (150 days) before the main epidemic peak. We assume that 1 in
every 1000 individuals in this age group is vaccinated every day; this
translates to approximately 28,900 individuals per day.

We assume that the vaccination is non-leaky, protecting vaccinated
individuals from infection, and this due to the model structure of the
‘default’ SEIR-V model in *epidemics* that we use here.

**New to modelling vaccination using *epidemics*?** It may help to read
the [“Modelling a vaccination
campaign”](https://epiverse-trace.github.io/epidemics/articles/modelling_vaccination.md)
vignette first.

We first create a `<vaccination>` object to represent this vaccination
regime.

``` r

# prepare a vaccination object
vaccinate_elders <- vaccination(
  name = "vaccinate elders",
  time_begin = matrix(1, nrow(contact_matrix)),
  time_end = matrix(150, nrow(contact_matrix)),
  nu = matrix(c(0.0, 0.0, 0.001))
)
```

We then pass this vaccination regime to the epidemic function in the
optional `vaccination` argument, and model 600 days of the epidemic. We
obtain the ‘final’ epidemic size using
[`epidemics::epidemic_size()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_size.md).

``` r

# model epidemic with vaccination prior to the main epidemic wave
output <- model_default(
  population = uk_population,
  transmission_rate = transmission_rate,
  infectiousness_rate = infectiousness_rate,
  recovery_rate = recovery_rate,
  vaccination = vaccinate_elders,
  time_end = 600, increment = 1.0
)

# Calculate the epidemic size using the helper function
finalsize_dat <- tibble(
  demography_group = names(demography_vector),
  value = epidemic_size(output) / demography_vector
)

# View the data
finalsize_dat
#> # A tibble: 3 × 2
#>   demography_group value
#>   <chr>            <dbl>
#> 1 [0,20)           0.617
#> 2 [20,40)          0.525
#> 3 [40,Inf)         0.343
```

### Calculating individuals vaccinated in epidemic model

To compare the results of the ODE model against the analytical method,
it is necessary to set up the population’s susceptibility and
demography-in-susceptibility matrices to reflect the proportion of
individuals in each age group vaccinated before the outbreak.

To do this, we need to calculate the total proportion of individuals in
each age group vaccinated at the end of the epidemic model run.

``` r

# Proportion of the individuals vaccinated
p_vacc <- filter(output, compartment == "vaccinated", time == max(time))
p_vacc <- p_vacc$value / demography_vector
```

### Implementing vaccination in *finalsize*

We use the proportion of individuals vaccinated in each age group to set
up matrices to pass to the `susceptibility` and `p_susceptibility`
arguments of
[`finalsize::final_size()`](https://epiverse-trace.github.io/finalsize/reference/final_size.html).

**New to final size estimation with heterogeneity in susceptibility to
infection?** It may help to read the [“Modelling heterogeneous
susceptibility”
vignette](https://epiverse-trace.github.io/finalsize/articles/varying_susceptibility.html),
and the [“Guide to constructing susceptibility
matrices”](https://epiverse-trace.github.io/finalsize/articles/susceptibility_matrices.html)
for *finalsize* first.

We create the susceptibility matrix to have two groups, ‘unvaccinated’,
with full susceptibility, and ‘vaccinated’, who are immune to infection.

``` r

# create the susceptibility matrix
susceptibility <- matrix(
  data = 1,
  nrow = length(demography_vector),
  ncol = 2,
  dimnames = list(
    names(demography_vector),
    c("unvaccinated", "vaccinated")
  )
)
```

We then create the demography-in-susceptibility matrix to reflect that
only some individuals in the \>40 age group are vaccinated. Since
vaccination has not been implemented for other age groups, we assume
that all individuals in those groups are fully susceptible.

``` r

# Second column holds the vaccinated (who are protected fully)
susceptibility[, "vaccinated"] <- 0

# Assume susceptibility varies within age groups
p_susceptibility <- matrix(
  data = 1.0,
  nrow = length(demography_vector),
  ncol = 2,
  dimnames = list(
    names(demography_vector),
    c("unvaccinated", "vaccinated")
  )
)
p_susceptibility[, "vaccinated"] <- p_vacc
p_susceptibility[, "unvaccinated"] <- 1 - p_vacc
```

We then calculate the expected final size using the modified
susceptibility matrices, and combine this data with the output of the
ODE model as in the previous example.

We need to scale the contact matrix appropriately.

``` r

# Define population in each age group
scalar <- max(eigen(contact_data$matrix)$values)
contact_matrix <- (contact_data$matrix / demography_vector) / scalar
```

``` r

# Calculate the proportion of individuals infected in each age group
dat1 <- final_size(
  r0 = 1.5,
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  susceptibility = susceptibility,
  p_susceptibility = p_susceptibility
)

# Final size returns the proportion infected in each susceptibility group
# (i.e. non vaccinated and vaccinated)
# here we calculate the proportion of infected for the age group as a whole
dat3 <- dat1 |>
  select(
    -susceptibility,
    final_size = p_infected
  ) |>
  filter(susc_grp == "unvaccinated") |>
  select(-susc_grp)

fs <- dat3$unvaccinated * p_susceptibility[, "unvaccinated"]

finalsize_dat <- finalsize_dat |>
  select(demo_grp = demography_group, seir_v = value) |>
  left_join(dat3)
```

We can print the data to compare values, which are similar.

``` r

finalsize_dat
#> # A tibble: 3 × 3
#>   demo_grp seir_v final_size
#>   <chr>     <dbl>      <dbl>
#> 1 [0,20)    0.617      0.660
#> 2 [20,40)   0.525      0.543
#> 3 [40,Inf)  0.343      0.272
```

However, note that while
[`epidemics::epidemic_size()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_size.md)
can be used for any time point in the epidemic,
[`finalsize::final_size()`](https://epiverse-trace.github.io/finalsize/reference/final_size.html)
only returns the size at the end of the epidemic, showing how
*epidemics* is more suitable to examine temporal dynamics.

``` r

# note that epidemic_size() returns the absolute values
# epidemic size after 10% of the epidemic model run
epidemic_size(data = output, stage = 0.1)
#> [1] 574.5195 488.4438 585.7900

# epidemic size at 50%
epidemic_size(data = output, stage = 0.5)
#> [1] 5569912 4894521 5262117
```

## Consideration of computational speed

An important reason to use *finalsize* for final size calculations may
be speed —
[`finalsize::final_size()`](https://epiverse-trace.github.io/finalsize/reference/final_size.html)
is much faster than
[`epidemics::epidemic_size()`](https://epiverse-trace.github.io/epidemics/reference/epidemic_size.md)
applied to
[`epidemics::model_default()`](https://epiverse-trace.github.io/epidemics/reference/model_default.md).
This makes it much more suitable for high-performance applications such
as fitting to data.

The benchmarking shows that the total time taken for 1,000 runs using
both methods is very low, and both methods could be used to scan across
multiple values of the input parameters, such as when dealing with
parameter uncertainty.

``` r

# run benchmarks using {bench}
benchmark <- mark(
  analytical_method = final_size(
    r0 = 1.5,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    susceptibility = susceptibility,
    p_susceptibility = p_susceptibility
  ),
  ode_model = epidemic_size(
    model_default(
      population = uk_population,
      transmission_rate = transmission_rate,
      infectiousness_rate = infectiousness_rate,
      recovery_rate = recovery_rate,
      vaccination = vaccinate_elders,
      time_end = 600, increment = 1.0
    )
  ),
  iterations = 1000,
  time_unit = "s",
  check = FALSE
)

# view the total time for 1000 runs in seconds
select(as_tibble(benchmark), expression, total_time)
#> # A tibble: 2 × 2
#>   expression        total_time
#>   <bch:expr>             <dbl>
#> 1 analytical_method      0.479
#> 2 ode_model             12.5
```

**Note** that some model runs using *epidemics* that implement more
complex compartmental structure, or multiple interventions and
vaccination regimes are likely to be slower than the example shown here.

However, users should choose *finalsize* only when the assumptions
underlying a final size calculation are met. For example, in cases where
vaccinations are concurrent with a large number of infections, or when
interventions are applied to reduce transmission, the final size
assumptions are not met. In these cases, users are advised to use a
dynamical model such as those in *epidemics* instead.

## References

Miller, Joel C. 2012. “A Note on the Derivation of Epidemic Final
Sizes.” *Bulletin of Mathematical Biology* 74 (9): 2125–41.
<https://doi.org/10.1007/s11538-012-9749-6>.
