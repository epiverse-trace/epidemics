% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{epidemic_peak}
\alias{epidemic_peak}
\title{Get the timing and size of a compartment's peak}
\usage{
epidemic_peak(data, compartment = "infectious")
}
\arguments{
\item{data}{A \verb{<data.frame>} of model output, typically
the output of \code{\link[=model_default_cpp]{model_default_cpp()}} or similar functions.}

\item{compartment}{The compartment of interest ;}
}
\value{
A \verb{<data.frame>} of the timing and peak of the selected compartment
for each of the demographic groups in \code{data}.
}
\description{
Get the timing and size of a compartment's peak for all demographic groups.
}
\examples{
# create a population
uk_population <- population(
  name = "UK population",
  contact_matrix = matrix(1),
  demography_vector = 67e6,
  initial_conditions = matrix(
    c(0.9999, 0.0001, 0, 0, 0),
    nrow = 1, ncol = 5L
  )
)

# run epidemic simulation with no vaccination or intervention
data <- model_default_cpp(
  population = uk_population
)

# get the timing and peak of the infectious compartment
epidemic_peak(data)

# get the timing and peak of the exposed compartment
epidemic_peak(data, compartment = "exposed")
}