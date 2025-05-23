% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_population.R
\name{gravity_contact}
\alias{gravity_contact}
\title{Create the connectivity matrix between all groups of \code{populations} using
a gravity model}
\usage{
gravity_contact(populations, connectivity_matrix)
}
\arguments{
\item{populations}{A list of \verb{<population>} objects}

\item{connectivity_matrix}{A numeric matrix for the contact matrix between
the elements of the populations}
}
\value{
A numeric matrix combining the contact matrices in \verb{<population>} and
the connectivity matrix using a gravity model
}
\description{
Create the connectivity matrix between all groups of \code{populations} using
a gravity model
}
\keyword{internal}
