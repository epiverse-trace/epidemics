
#' Construct a new intervention
#'
#' @param name String for the name of the intervention.
#' @param time_begin Single number for the start time of the intervention.
#' @param time_end Single number for the end time of the intervention.
#' @param contact_reduction A vector of the same length as the
#' number of demographic groups in the target population, which gives the
#' overall or group-specific proportion reduction in contacts respectively.
#'
#' @return An `intervention` class object.
#' @keywords internal
#' @noRd
new_intervention <- function(name = NA_character_,
                             time_begin,
                             time_end,
                             contact_reduction) {
  # create and return intervention class
  structure(
    list(
      name = name,
      time_begin = time_begin,
      time_end = time_end,
      contact_reduction = contact_reduction
    ),
    class = "intervention"
  )
}

#' Construct a new intervention for an epidemic model
#'
#' @param name String for the name of the intervention.
#' @param time_begin Single number for the start time of the intervention.
#' @param time_end Single number for the end time of the intervention.
#' @param contact_reduction A vector of the same length as the
#' number of demographic groups in the target population, which gives the
#' overall or group-specific proportion reduction in contacts respectively.
#'
#' @return An object of the `intervention` S3 class.
#' @export
#'
#' @examples
#' # assuming a population with two age groups, 0 -- 18, and 18+
#' # an example in which schools are closed for 30 days (or other time units)
#' close_schools <- intervention(
#'   name = "close schools",
#'   time_begin = 50,
#'   time_end = 80,
#'   contact_reduction = c(0.5, 0.01) # reduces contacts differentially
#' )
#' close_schools
intervention <- function(name = NA_character_,
                         time_begin,
                         time_end,
                         contact_reduction) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_number(time_begin, lower = 0, finite = TRUE)
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_numeric(contact_reduction)

  # print message if time end is not before time begin
  # but allow it nonetheless
  if (!(time_end > time_begin)) {
    message(
      "Intervention: `time_end` is not greater than `time_begin`"
    )
  }

  # call intervention constructor
  intervention_ <- new_intervention(
    name = name,
    time_begin = time_begin,
    time_end = time_end,
    contact_reduction = contact_reduction
  )

  # call intervention validator
  validate_intervention(object = intervention_)

  # return intervention object
  intervention_
}

#' Validate an intervention
#'
#' @param object An object to be validated as an `intervention`.
#'
#' @return No return.
#' @noRd
#' @keywords internal
validate_intervention <- function(object) {
  # check for class and class invariants
  stopifnot(
    "Object should be of class `intervention`" =
      (is_intervention(object)),
    "`intervention` does not contain the correct attributes" =
      (c(
        "name", "time_begin", "time_end", "contact_reduction"
      ) %in% attributes(object)$names)
  )

  # check intervention class members
  checkmate::assert_string(object$name, na.ok = TRUE)
  checkmate::assert_number(object$time_begin, lower = 0, finite = TRUE)
  checkmate::assert_number(object$time_end, lower = 0, finite = TRUE)
  checkmate::assert_numeric(object$contact_reduction)

  invisible(object)
}

#' Check whether an object is an `intervention`
#'
#' @param object An object to be checked as being an `intervention`.
#'
#' @return A logical for whether the object is of the `intervention` class.
#' @export
#'
#' @examples
#' close_schools <- intervention(
#'   name = "close schools",
#'   time_begin = 50,
#'   time_end = 80,
#'   contact_reduction = c(0.5, 0.01) # reduces contacts differentially
#' )
#' is_intervention(close_schools)
is_intervention <- function(object) {
  inherits(object, "intervention")
}

#' Generate a null intervention
#'
#' @param population A `population` object with a `contact_matrix` member.
#'
#' @return An intervention that has no effect on contacts, with start and end
#' times set to 0.0
#' @export
no_intervention <- function(population) {
  checkmate::assert_class(population, "population")
  intervention(
    name = "no_intervention", time_begin = 0, time_end = 0,
    contact_reduction = rep(0.0, times = nrow(population$contact_matrix))
  )
}

#' Print a `intervention` object
#'
#' @param x A `intervention` object.
#' @param ... Other parameters passed to [print()].
#' @return Invisibly returns the [`intervention`] object `x`.
#' Called for printing side-effects.
#' @export
print.intervention <- function(x, ...) {
  format(x, ...)
}

#' Format a `intervention` object
#'
#' @param x A `intervention` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisibly returns the [`intervention`] object `x`.
#' Called for printing side-effects.
#' @keywords internal
#' @noRd
format.intervention <- function(x, ...) {

  # validate the intervention object
  validate_intervention(x)

  # header
  header <- "<intervention>"

  # collect information on name
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  name <- glue::glue("Intervention name: {name}")

  # print to screen
  writeLines(
    c(
      header,
      name,
      glue::glue(
        "

        Time begin: {x$time_begin}
        Time end: {x$time_end}
        "
      )
    )
  )
  print(
    glue::glue(
      "Contact reduction:
      {glue::glue_collapse(x$contact_reduction, sep = ', ')}
      "
    )
  )

  invisible(x)
}
