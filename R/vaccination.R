
#' Construct a new vaccination regime
#'
#' @param name String for the name of the vaccination regime.
#' @param time_begin Vector for the start time of the vaccination for each
#' demographic group \eqn{i}.
#' @param time_end Vector for the end time of the vaccination for each
#' demographic group \eqn{i}.
#' @param nu Vector of the same length as the
#' number of demographic groups, which gives the group-specific rate of
#' vaccination, expressed as the rate parameter \eqn{nu}.
#'
#' @return An `vaccination` class object.
#' @keywords internal
#' @noRd
new_vaccination <- function(name = NA_character_,
                            time_begin,
                            time_end,
                            nu) {
  # create and return scenario class
  structure(
    list(
      name = name,
      time_begin = time_begin,
      time_end = time_end,
      nu = nu
    ),
    class = "vaccination"
  )
}

#' Construct a new vaccination regime for an epidemic model
#'
#' @param name String for the name of the vaccination regime.
#' @param time_begin Vector for the start time of the vaccination for each
#' demographic group \eqn{i}.
#' @param time_end Vector for the end time of the vaccination for each
#' demographic group \eqn{i}.
#' @param nu Vector of the same length as the
#' number of demographic groups, which gives the group-specific rate of
#' vaccination, expressed as the rate parameter \eqn{nu}.
#'
#' @return An object of the `vaccination` S3 class.
#' @export
#'
#' @examples
#' # assuming a population with two age groups, children 0 -- 5, and others 5+
#' # an example for childhood vaccination only
#' childhood_vaccination <- vaccination(
#'   name = "childhood_vaccination",
#'   time_begin = c(0, 100), # assuming a simulation over 100 days
#'   time_end = c(100, 100),
#'   nu = c(0.0001, 0.0) # over 5s never vaccinated
#' )
vaccination <- function(name = NA_character_,
                        time_begin,
                        time_end,
                        nu) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_numeric(time_begin, lower = 0, finite = TRUE)
  checkmate::assert_numeric(time_end, lower = 0, finite = TRUE)
  checkmate::assert_numeric(nu, finite = TRUE)

  # call scenario constructor
  vaccination_ <- new_vaccination(
    name = name,
    time_begin = time_begin,
    time_end = time_end,
    nu = nu
  )

  # call vaccination validator
  validate_vaccination(object = vaccination_)

  # return vaccination object
  vaccination_
}

#' Validate a `vaccination` object
#'
#' @param object An object to be validated as a `vaccination`.
#'
#' @return No return.
#' @noRd
#' @keywords internal
validate_vaccination <- function(object) {
  # check for class and class invariants
  stopifnot(
    "Object should be of class `vaccination`" =
      (is_vaccination(object)),
    "`vaccination` does not contain the correct attributes" =
      (c(
        "name", "time_begin", "time_end", "nu"
      ) %in% attributes(object)$names)
  )
  # other checks to be added here
  invisible(object)
}

#' Check whether an object is a `vaccination`
#'
#' @param object An object to be checked as being a `vaccination`.
#'
#' @return A logical for whether the object is of the `vaccination` class.
#' @export
#'
#' @examples
#' # an example for childhood vaccination only
#' childhood_vaccination <- vaccination(
#'   name = "childhood_vaccination",
#'   time_begin = c(0, 100), # assuming a simulation over 100 days
#'   time_end = c(100, 100),
#'   nu = c(0.0001, 0.0) # over 5s never vaccinated
#' )
#' is_vaccination(childhood_vaccination)
is_vaccination <- function(object) {
  inherits(object, "vaccination")
}

#' Generate a null vaccination
#'
#' @param population A `population` object with a `contact_matrix` member.
#'
#' @return An vaccination that has no effect on the population, with start and
#' end times set to 0.0, and the rate of vaccination \eqn{nu} also set to 0.0.
#' @export
no_vaccination <- function(population) {
  checkmate::assert_class(population, "population")
  vaccination(
    name = "no_vaccination", time_begin = 0, time_end = 0,
    nu = rep(0.0, times = nrow(population$contact_matrix))
  )
}