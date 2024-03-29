% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_vacamole.R
\name{.ode_model_vacamole}
\alias{.ode_model_vacamole}
\title{Ordinary Differential Equations for the Vacamole Model}
\usage{
.ode_model_vacamole(t, y, params)
}
\arguments{
\item{t}{A single number of the timestep at which to integrate.}

\item{y}{The conditions of the epidemiological compartments.}

\item{params}{The parameters, passed as a named list.}
}
\value{
A list with a vector with as many elements as the number of
demographic groups times the number of epidemiological compartments. Each
value gives the change in the number of individuals in that compartment.
}
\description{
Provides the ODEs for the RIVM Vacamole model in a format that
is suitable for passing to \code{\link[deSolve:lsoda]{deSolve::lsoda()}}.
See \code{\link[=model_vacamole]{model_vacamole()}} for a list of required parameters.
}
\keyword{internal}
