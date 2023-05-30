
#' Construct a new vaccination regime
#'
#' @param name String for the name of the vaccination regime.
#' @param time_begin Matrix for the start time of delivering vaccination dose
#' \eqn{j} to demographic group \eqn{i}.
#' demographic group \eqn{i}.
#' @param time_end Matrix for the end time of delivering vaccination dose
#' \eqn{j} to demographic group \eqn{i}.
#' @param nu Matrix for the group-specific rate of vaccination, expressed as the
#' rate parameter \eqn{nu}. Each element of the matrix \eqn{nu_{ij}} represents
#' the rate of delivering vaccine dose \eqn{j} to demographic group \eqn{i}.
#'
#' @return A `vaccination` class object.
#' @keywords internal
#' @noRd
new_vaccination <- function(name = NA_character_,
                            time_begin,
                            time_end,
                            nu) {
  # create and return vaccination class
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
#' @param time_begin Matrix for the start time of delivering vaccination dose
#' \eqn{j} to demographic group \eqn{i}.
#' demographic group \eqn{i}.
#' @param time_end Matrix for the end time of delivering vaccination dose
#' \eqn{j} to demographic group \eqn{i}.
#' @param nu Matrix for the group-specific rate of vaccination, expressed as the
#' rate parameter \eqn{nu}. Each element of the matrix \eqn{nu_{ij}} represents
#' the rate of delivering vaccine dose \eqn{j} to demographic group \eqn{i}.
#'
#' @return An object of the `vaccination` S3 class.
#' @export
#'
#' @examples
#' # assuming a population with two age groups, children 0 -- 5, and others 5+
#' # an example for childhood vaccination only
#' childhood_vaccination <- vaccination(
#'   name = "childhood_vaccination",
#'   time_begin = matrix(c(0, 100)), # assuming a simulation over 100 days
#'   time_end = matrix(c(100, 100)),
#'   nu = matrix(c(0.0001, 0.0)) # over 5s never vaccinated
#' )
#' childhood_vaccination
vaccination <- function(name = NA_character_,
                        nu,
                        time_begin,
                        time_end) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_matrix(nu)
  checkmate::assert_matrix(time_begin, nrows = nrow(nu), ncols = ncol(nu))
  checkmate::assert_matrix(time_end, nrows = nrow(nu), ncols = ncol(nu))

  # message if any vaccinations' intervals are badly formed
  if (any(time_end <= time_begin)) {
    message(
      "Vaccination: some `time_end`s are not greater than `time_begin`s"
    )
  }

  # assign dose names
  dose_names <- glue::glue("dose_{seq(ncol(nu))}")
  colnames(nu) <- dose_names
  colnames(time_begin) <- dose_names
  colnames(time_end) <- dose_names

  # call vaccination constructor
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

  # other checks for the vaccination object
  checkmate::assert_string(object$name, na.ok = TRUE)
  checkmate::assert_matrix(object$nu)
  checkmate::assert_matrix(
    object$time_begin,
    nrows = nrow(object$nu), ncols = ncol(object$nu)
  )
  checkmate::assert_matrix(
    object$time_end,
    nrows = nrow(object$nu), ncols = ncol(object$nu)
  )

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
#'   time_begin = matrix(c(0, 100)), # assuming a simulation over 100 days
#'   time_end = matrix(c(100, 100)),
#'   nu = matrix(c(0.0001, 0.0)) # over 5s never vaccinated
#' )
#' is_vaccination(childhood_vaccination)
is_vaccination <- function(object) {
  inherits(object, "vaccination")
}

#' Generate a null vaccination
#'
#' @param population A `population` object with a `contact_matrix` member.
#' @param doses A number, defaulting to 1, to indicate the number of doses in
#' the vaccination regime.
#' @return An vaccination that has no effect on the population, with start and
#' end times set to 0.0, and the rate of vaccination \eqn{nu} also set to 0.0.
#' @export
no_vaccination <- function(population, doses = 1L) {
  checkmate::assert_class(population, "population")
  vaccination(
    name = "no_vaccination",
    time_begin = matrix(
      0.0,
      nrow = nrow(population$contact_matrix), ncol = doses
    ),
    time_end = matrix(
      0.0,
      nrow = nrow(population$contact_matrix), ncol = doses
    ),
    nu = matrix(
      0.0,
      nrow = nrow(population$contact_matrix), ncol = doses
    )
  )
}

#' Print a `vaccination` object
#'
#' @param x A `vaccination` object.
#' @param ... Other parameters passed to [print()].
#' @return Invisibly returns the [`vaccination`] object `x`.
#' Called for printing side-effects.
#' @export
print.vaccination <- function(x, ...) {
  format(x, ...)
}

#' Format a `vaccination` object
#'
#' @param x A `vaccination` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisibly returns the [`vaccination`] object `x`.
#' Called for printing side-effects.
#' @keywords internal
#' @noRd
format.vaccination <- function(x, ...) {

  # validate the vaccination object
  validate_vaccination(x)

  # header
  header <- "<vaccination>"

  # collect information on name
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  name <- glue::glue("Vaccination name: {name}")

  # print to screen
  writeLines(
    c(
      header,
      name,
      glue::glue(
        "

        Time begin:
        "
      )
    )
  )
  print(x$time_begin)

  print(glue::glue(
    "

    Time end:
    "
  ))
  print(x$time_end)

  print(
    glue::glue(
      "

      Vaccination rate:
      "
    )
  )
  print(x$nu)

  invisible(x)
}

#' Convert a list to a vaccination object
#'
#' @param x A list, or an object that inherits from a list.
#' @return A [vaccination] class object.
#' @export
#' @examples
#' # prepare a list
#' vax <- list(
#'   name = "vax_regime",
#'   time_begin = matrix(1),
#'   time_end = matrix(100),
#'   nu = matrix(0.001)
#' )
#'
#' as.vaccination(vax)
as.vaccination <- function(x) {
  # check that input is a list or vaccination
  stopifnot(
    "Input must inherit from `list`" =
      is.list(x)
  )

  x <- vaccination(
    name = x[["name"]],
    time_begin = x[["time_begin"]],
    time_end = x[["time_end"]],
    nu = x[["nu"]]
  )

  # return x
  x
}
