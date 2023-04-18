
# *epidemics*: A library of compartmental epidemic scenario models

<!-- <img src="man/figures/logo.png" align="right" width="130"/> -->

*epidemics* is an R package that provides an easy interface to a library
of compartmental models that can help to model epidemic scenarios for
directly transmitted infections such as influenza, Covid-19, or
respiratory syncytial virus (RSV).

*epidemics* currently provides a single model with susceptible, exposed,
infectious, recovered, and vaccinated compartments (SEIR-V), allowing
for heterogeneity in social contacts, the implementation of a
group-specific non-pharmaceutical intervention that reduces social
contacts, and a vaccination regime with group-specific start and end
dates.

*epidemics* implements methods outlined in Bjørnstad et al.
([2020b](#ref-bjornstad2020)) and Bjørnstad et al.
([2020a](#ref-bjornstad2020a)).

The models in *epidemics* can help provide rough estimates of the course
of epidemics, and the effectiveness of pharmaceutical and
non-pharmaceutical interventions

*epidemics* relies on [Eigen](https://gitlab.com/libeigen/eigen) via
[{RcppEigen}](https://cran.r-project.org/web/packages/RcppEigen/index.html),
and on [Boost
Odeint](https://www.boost.org/doc/libs/1_82_0/libs/numeric/odeint/doc/html/index.html)
via [{BH}](https://cran.r-project.org/web/packages/BH/index.html), and
is developed at the [Centre for the Mathematical Modelling of Infectious
Diseases](https://www.lshtm.ac.uk/research/centres/centre-mathematical-modelling-infectious-diseases)
at the London School of Hygiene and Tropical Medicine as part of the
[Epiverse-TRACE initiative](https://data.org/initiatives/epiverse/).

<!-- badges: start -->

[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R-CMD-check](https://github.com/epiverse-trace/epidemics/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/epiverse-trace/epidemics/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/epiverse-trace/epidemics/branch/main/graph/badge.svg)](https://app.codecov.io/gh/epiverse-trace/epidemics?branch=main)
[![Project Status: WIP – Initial development is in progress, but there
has not yet been a stable, usable release suitable for the
public.](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![CRAN
status](https://www.r-pkg.org/badges/version/epidemics)](https://CRAN.R-project.org/package=epidemics)
<!-- badges: end -->

## Installation

The current development version of *finalsize* can be installed from
*epidemics* from [GitHub](https://github.com/) using the `pak` package:

``` r
if(!require("pak")) install.packages("pak")
pak::pak("epiverse-trace/epidemics")
```

## Quick start

*epidemics* provides the single function …

Here, an example using …

## Package vignettes

More details on how to use *epidemics* can be found in the [online
documentation as package
vignettes](https://epiverse-trace.github.io/epidemics/), under
“Articles”.

## Help

To report a bug please open an
[issue](https://github.com/epiverse-trace/epidemics/issues/new/choose).

## Contribute

Contributions to *epidemics* are welcomed. Please follow the [package
contributing
guide](https://github.com/epiverse-trace/epidemics/blob/main/.github/CONTRIBUTING.md).

## Code of conduct

Please note that the *epidemics* project is released with a [Contributor
Code of
Conduct](https://github.com/epiverse-trace/.github/blob/main/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-bjornstad2020a" class="csl-entry">

Bjørnstad, Ottar N., Katriona Shea, Martin Krzywinski, and Naomi Altman.
2020a. “Modeling Infectious Epidemics.” *Nature Methods* 17 (5): 455–56.
<https://doi.org/10.1038/s41592-020-0822-z>.

</div>

<div id="ref-bjornstad2020" class="csl-entry">

———. 2020b. “The SEIRS Model for Infectious Disease Dynamics.” *Nature
Methods* 17 (6): 557–58. <https://doi.org/10.1038/s41592-020-0856-2>.

</div>

</div>
