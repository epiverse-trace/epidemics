% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_ebola.R
\name{.model_ebola_internal}
\alias{.model_ebola_internal}
\title{Internal code for the Ebola model}
\usage{
.model_ebola_internal(
  initial_state,
  erlang_subcompartments,
  transmission_rate,
  infectiousness_rate,
  removal_rate,
  prop_community,
  etu_risk,
  funeral_risk,
  intervention,
  time_dependence,
  time_end,
  replicates,
  local_seeds
)
}
\arguments{
\item{erlang_subcompartments}{A numeric, integer-like vector for the number
of Erlang sub-compartments assumed for the exposed, infectious, and
hospitalised compartments. Defaults to 2.}

\item{transmission_rate}{A numeric vector for the rate at which individuals
move from the susceptible to the exposed compartment upon contact with an
infectious individual. Often denoted as \eqn{\beta}, with
\eqn{\beta = R_0 / \text{infectious period}}.
See \strong{Details} for default values.}

\item{infectiousness_rate}{A numeric vector for the rate at which individuals
move from the exposed to the infectious compartment. Often denoted as
\eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
This value does not depend upon the number of infectious individuals in the
population.
See \strong{Details} for default values.}

\item{removal_rate}{A numeric vector for the rate at which infectious
individuals transition from the infectious or hospitalised compartments to
the funeral or removed compartments.
This model does not distinguish between recoveries and deaths.
Denoted in Getz and Dougherty as \eqn{\gamma^I} (see \strong{Details}).}

\item{prop_community}{A numeric vector for the proportion of infectious
individuals who remain in the community and are not hospitalised for
treatment. Defaults to 0.9.}

\item{etu_risk}{A numeric vector for the relative risk of onward transmission
of EVD from hospitalised individuals, with values between 0.0 and 1.0, where
0.0 indicates that hospitalisation completely prevents onward transmission,
and 1.0 indicates that hospitalisation does not prevent onward transmission
at all; values are relative to the baseline transmission rate \eqn{\beta}.
Defaults to 0.7.}

\item{funeral_risk}{A numeric vector for the relative risk of onward
transmission of EVD from funerals of individuals who died with EVD.
Must be between 0.0 and 1.0, where 0.0 indicates that there is no onward
transmission, and 1.0 indicates that funeral transmission is equivalent to
the baseline transmission rate in the community \eqn{\beta}.
Defaults to 0.5.}

\item{intervention}{An optional named list of \verb{<rate_intervention>} objects
representing optional pharmaceutical or non-pharmaceutical interventions
applied to the model parameters listed above. May also be a list of such
lists, in which case each set of interventions is treated as a separate
scenario. See \strong{Details} below.}

\item{time_dependence}{An optional named list where each element is a
function with the first two arguments being the current simulation \code{time},
and \code{x}, a value that is dependent on \code{time}
(\code{x} represents a model parameter).
List names must correspond to model parameters modified by the function.
Alternatively, may be a list of such lists, in which case each set of
functions is treated as a distinct scenario.
See \strong{Details} for more information, as well as the vignette on time-
dependence \code{vignette("time_dependence", package = "epidemics")}.}

\item{time_end}{A numeric, integer-like vector for the maximum number of}

\item{replicates}{A single number for replicates to run. Defaults to 100.
timesteps over which to run the model, in days. Defaults to 100 days.}
}
\description{
Internal code for the Ebola model
}
\details{
\strong{NOTE} that all arguments to this internal function must be scalars for
infection parameters or single intervention or time-dependence sets for
composable elements.
}
\keyword{internal}
