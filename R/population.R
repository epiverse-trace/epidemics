
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
  # create and return population class
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

  # call population constructor
  population_ <- new_population(
    name = name,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions
  )

  # call population validator
  validate_population(object = population_)

  # return population object
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

#' Print a `population` object
#'
#' @param x A `population` object.
#' @param ... Other parameters passed to [print()].
#' @return None. Prints output.
#' @export
print.population <- function(x, ...) {
  format(x, ...)
}

#' Format a `population` object
#'
#' @param x A `population` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisbily returns the [`population`] object `x`. Called for printing
#' side-effects.
#' @keywords internal
#' @noRd
format.population <- function(x, ...) {

  # validate the population object
  validate_population(x)

  # header
  header <- "<population>"

  # collect information on name
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  name <- glue::glue("Population name: {name}")

  # copy demography vector
  demography_print <- prettyNum(
    x$demography_vector,
    big.mark = ",", scientific = FALSE
  )
  demography_proportions <- round(
    x$demography_vector / sum(x$demography_vector), 1
  ) * 100.0

  demography_print <- glue::glue(
    "{demography_print} ({demography_proportions}%)"
  )
  names(demography_print) <- names(x$demography_vector)

  contact_matrix <- x$contact_matrix

  # demographic group names
  if (is.null(names(x$demography_print))) {
    if (is.null(rownames(x$contact_matrix))) {
      # assign names to demography vector and contacts
      demography_print <-
        glue::glue(
          "Dem. grp. {seq_along(demography_print)}: {demography_print}"
        )
      rownames(contact_matrix) <-
        glue::glue("Dem. grp. {seq_along(demography_print)}:")

      if (is.null(colnames(contact_matrix))) {
        colnames(contact_matrix) <-
          glue::glue("Dem. grp. {seq_along(demography_print)}:")
      }
    } else {
      demography_print <-
        glue::glue("{rownames(x$contact_matrix)}: {demography_print}")
      colnames(contact_matrix) <- rownames(x$contact_matrix)
    }
  }


  # print to screen
  writeLines(
    c(
      header,
      name,
      glue::glue(
        "

        Demography:
        "
      )
    )
  )
  print(
    demography_print
  )
  print(
    glue::glue(
      "

    Contact matrix:
    "
    )
  )
  print(round(contact_matrix, 1))

  invisible(x)
}
