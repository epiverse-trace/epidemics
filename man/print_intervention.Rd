% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intervention.R
\name{print_intervention}
\alias{print_intervention}
\alias{print.intervention}
\alias{print.contact_intervention}
\alias{print.rate_intervention}
\title{Print an object of the \verb{<intervention>} super-class}
\usage{
\method{print}{intervention}(x, ...)

\method{print}{contact_intervention}(x, ...)

\method{print}{rate_intervention}(x, ...)
}
\arguments{
\item{x}{An object that inherits from the \verb{<intervention>} class.
For \code{print.contacts_intervention()}, an object of the
\verb{<contacts_intervention>} class. For \code{print.rate_intervention()}, an object
of the \verb{<rate_intervention>} class.}

\item{...}{Other parameters passed to \code{\link[=print]{print()}}.}
}
\value{
Invisibly returns the object \code{x}.
Called for printing side-effects.
}
\description{
Print an object of the \verb{<intervention>} super-class

Print objects of the \verb{<contact_intervention>} class

Print objects of the \verb{<rate_intervention>} class
}
