% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/population.R
\name{population}
\alias{population}
\alias{is_population}
\title{Construct a new population for an epidemic model}
\usage{
population(
  name = NA_character_,
  contact_matrix,
  demography_vector,
  initial_conditions
)

is_population(x)
}
\arguments{
\item{name}{Optional string for the population name.}

\item{contact_matrix}{A matrix giving the contacts between the demographic
groups in the population. Must be a square matrix.}

\item{demography_vector}{A vector of the sizes of each demographic group.
Must have the same length as the dimensions of the contact matrix.}

\item{initial_conditions}{Matrix representing the initial proportions of each
demographic group in the four model compartments: 'susceptible',
'infected/infectious', 'recovered', and 'vaccinated'. Must have as many rows
as the number of demographic groups. Each compartment is represented in the
columns of the matrix, so that the element \eqn{M_{ij}} represents the
proportion of individuals of demographic group \eqn{i} in compartment \eqn{j}
.}

\item{x}{An object to be checked as a valid population.}
}
\value{
An object of the \verb{<population>} S3 class.

\code{is_population()} returns a logical for whether the object is a
\verb{<population>}.
}
\description{
Construct a new population for an epidemic model

Check whether an object is a \verb{<population>}
}
\examples{
uk_pop <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0),
    nrow = 1, ncol = 4
  )
)

# print to check
uk_pop

# check for class <population>
is_population(uk_pop)
}
