
#' Construct a new intervention
#'
#' @param name String for the name of the intervention.
#' @param time_begin Single number for the start time of the intervention.
#' @param time_end Single number for the end time of the intervention.
#' @param contact_reduction Single number, or a vector of the same length as the
#' number of demographic groups, which gives the overall or group-specific
#' proportion reduction in contacts respectively.
#'
#' @return An `intervention` class object.
#' @keywords internal
#' @noRd
new_intervention <- function(name = NA_character_,
                             time_begin,
                             time_end,
                             contact_reduction) {
  # create and return scenario class
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
#' @param contact_reduction Single number, or a vector of the same length as the
#' number of demographic groups, which gives the overall or group-specific
#' proportion reduction in contacts respectively.
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
intervention <- function(name = NA_character_,
                         time_begin,
                         time_end,
                         contact_reduction) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_number(time_begin, lower = 0, finite = TRUE)
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)
  checkmate::assert_numeric(contact_reduction)

  # call scenario constructor
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
