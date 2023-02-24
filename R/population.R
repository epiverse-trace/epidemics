
#' Construct a new population
#'
#' @param name Optional string for the population name.
#' @param contact_matrix A matrix giving the contacts between the demographic
#' groups in the population. Must be a square matrix.
#' @param demography_vector A vector of the sizes of each demographic group.
#' Must have the same length as the dimensions of the contact matrix.
#' @param initial_conditions Matrix representing the initial proportions of each
#' demographic group in the four model compartments: 'susceptible',
#' 'infected/infectious', 'recovered', and 'vaccinated'. Must have as many rows
#' as the number of demographic groups. Each compartment is represented in the
#' columns of the matrix, so that the element \eqn{M_{ij}} represents the
#' proportion of individuals of demographic group \eqn{i} in compartment \eqn{j}
#' .
#' @keywords internal
#' @return An object of the `population` class.
#' @noRd
new_population <- function(name = NA_character_,
                           contact_matrix = matrix(),
                           demography_vector = numeric(),
                           initial_conditions = matrix()) {
  # create and return scenario class
  structure(
    list(
      name = name,
      contact_matrix = contact_matrix,
      demography_vector = demography_vector,
      initial_conditions = initial_conditions
    ),
    class = "population"
  )
}

#' Construct a new population for an epidemic model
#'
#' @param name Optional string for the population name.
#' @param contact_matrix A matrix giving the contacts between the demographic
#' groups in the population. Must be a square matrix.
#' @param demography_vector A vector of the sizes of each demographic group.
#' Must have the same length as the dimensions of the contact matrix.
#' @param initial_conditions Matrix representing the initial proportions of each
#' demographic group in the four model compartments: 'susceptible',
#' 'infected/infectious', 'recovered', and 'vaccinated'. Must have as many rows
#' as the number of demographic groups. Each compartment is represented in the
#' columns of the matrix, so that the element \eqn{M_{ij}} represents the
#' proportion of individuals of demographic group \eqn{i} in compartment \eqn{j}
#' .
#'
#' @return An object of the `population` S3 class.
#' @export
#'
#' @examples
#' population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0),
#'     nrow = 1, ncol = 4
#'   )
#' )
population <- function(name = NA_character_,
                       contact_matrix = matrix(1),
                       demography_vector = 67e6,
                       initial_conditions = matrix()) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_matrix(contact_matrix)
  checkmate::assert_numeric(demography_vector)
  checkmate::assert_matrix(initial_conditions)
  stopifnot(
    "Initial conditions must have same number of rows as contact_matrix" =
      nrow(initial_conditions) == nrow(contact_matrix)
  )

  # call scenario constructor
  population_ <- new_population(
    name = name,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions
  )

  # call scenario validator
  validate_population(object = population_)

  # return scenario object
  population_
}

#' Validate a `population` object
#'
#' @param object A `population` object for validation.
#'
#' @return A validated `population` object.
#' @keywords internal
#' @noRd
validate_population <- function(object) {
  # check for class and class invariants
  stopifnot(
    "Object should be of class `population`" =
      (is_population(object)),
    "`population` does not contain the correct attributes" =
      (c(
        "name", "contact_matrix", "demography_vector"
      ) %in% attributes(object)$names),
    "`contact_matrix` must have as many rows as length of `demography_vector`" =
      (nrow(object$contact_matrix) == length(object$demography_vector)),
    "`initial_conditions` must have as many rows as `demography_vector`" =
      (nrow(object$initial_conditions) == length(object$demography_vector)),
    "`initial_conditions` rows must always sum to 1.0" =
      (all(abs(rowSums(object$initial_conditions) - 1) < 1e-6))
  )
  invisible(object)
}

#' Check whether an object is a `population`
#'
#' @param object An object to be checked as a valid population.
#'
#' @return A logical for whether the object is a `population`.
#' @export
#'
#' @examples
#' # for the UK
#' new_pop <- population(
#'   initial_conditions = matrix(
#'     c(0.999, 0.001, 0, 0),
#'     nrow = 1, ncol = 4
#'   )
#' )
#' is_population(new_pop)
is_population <- function(object) {
  inherits(object, "population")
}
