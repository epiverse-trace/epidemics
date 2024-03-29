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
#' @return An object of the `<population>` class.
#' @noRd
new_population <- function(name = NA_character_,
                           contact_matrix = matrix(),
                           demography_vector = numeric(),
                           initial_conditions = matrix()) {
  # create and return population class
  x <- list(
    name = name,
    contact_matrix = contact_matrix,
    demography_vector = demography_vector,
    initial_conditions = initial_conditions
  )
  class(x) <- "population"

  x
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
#' @param x An object to be checked as a valid population.
#'
#' @return An object of the `<population>` S3 class.
#'
#' `is_population()` returns a logical for whether the object is a
#' `<population>`.
#' @export
#'
#' @examples
#' uk_pop <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0),
#'     nrow = 1, ncol = 4
#'   )
#' )
#'
#' # print to check
#' uk_pop
#'
#' # check for class <population>
#' is_population(uk_pop)
population <- function(name = NA_character_,
                       contact_matrix,
                       demography_vector,
                       initial_conditions) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_matrix(contact_matrix, mode = "numeric")
  checkmate::assert_numeric(
    demography_vector,
    lower = 0, finite = TRUE, any.missing = FALSE,
    len = nrow(contact_matrix)
  )
  checkmate::assert_matrix(
    initial_conditions,
    mode = "numeric", nrows = length(demography_vector)
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

#' Validate a `<population>` object
#'
#' @param object A `<population>` object for validation.
#'
#' @return A validated `<population>` object.
#' @keywords internal
#' @noRd
validate_population <- function(object) {
  # check for class and class invariants
  stopifnot(
    "Object should be of class <population>" =
      (is_population(object)),
    "<population> does not contain the correct attributes" =
      (c(
        "name", "contact_matrix", "demography_vector"
      ) %in% attributes(object)$names),
    "`name` must be a string" =
      checkmate::test_string(object$name, na.ok = TRUE),
    "`contact_matrix` must be a numeric matrix with same rows as demography" =
      checkmate::test_matrix(
        object$contact_matrix,
        mode = "numeric",
        nrows = length(object$demography_vector)
      ),
    "`initial_conditions` must be a numeric matrix" =
      checkmate::test_matrix(
        object$initial_conditions,
        mode = "numeric",
        nrows = length(object$demography_vector)
      ),
    "`initial_conditions` rows must always sum to 1.0" =
      (all(abs(rowSums(object$initial_conditions) - 1.0) < 1e-3))
  )
  invisible(object)
}

#' Check whether an object is a `<population>`
#'
#' @name population
#' @rdname population
#'
#' @export
#'
is_population <- function(x) {
  inherits(x, "population")
}

#' Print a `<population>` object
#'
#' @param x A `<population>` object.
#' @param ... Other parameters passed to [print()].
#' @return Invisibly returns the `<population>` object `x`.
#' Called for printing side-effects.
#' @export
print.population <- function(x, ...) {
  format(x, ...)
}

#' Format a `<population>` object
#'
#' @param x A `<population>` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisibly returns the `<population>` object `x`. Called for printing
#' side-effects.
#' @keywords internal
#' @noRd
format.population <- function(x, ...) {
  # validate the population object
  validate_population(x)

  # header
  header <- class(x) # nolint: object_usage_linter

  # collect information on name
  # nolint start: object_usage_linter
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  # nolint end
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
  cat(
    cli::cli_text(
      "{.cls {header}} object"
    )
  )
  cat(
    "\n",
    cli::col_blue(
      "Population name: "
    )
  )
  cli::cli_text(
    "{cli::cli_format({name}, style = list(string_quote = \"\"))}"
  )
  cat(
    "\n",
    cli::col_blue(
      "Demography"
    ),
    "\n"
  )
  print(demography_print)
  cat(
    "\n",
    cli::col_blue(
      "Contact matrix"
    ),
    "\n"
  )
  print(contact_matrix)

  invisible(x)
}
