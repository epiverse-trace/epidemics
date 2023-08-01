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
#'   contact_reduction = matrix(c(0.5, 0.01)) # reduces contacts differentially
#' )
#' close_schools
intervention <- function(name = NA_character_,
                         time_begin,
                         time_end,
                         contact_reduction) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_matrix(contact_reduction, mode = "numeric")
  checkmate::assert_number(time_begin, lower = 0, finite = TRUE)
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)

  # message if any vaccinations' intervals are badly formed
  if (any(time_end <= time_begin)) {
    message(
      "Vaccination: some `time_end`s are not greater than `time_begin`s"
    )
  }

  # call intervention constructor
  intervention_ <- new_intervention(
    name = name,
    time_begin = matrix(time_begin),
    time_end = matrix(time_end),
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
  checkmate::assert_matrix(object$contact_reduction, mode = "numeric")
  checkmate::assert_matrix(
    object$time_begin,
    ncols = ncol(object$contact_reduction), nrows = 1L
  )
  checkmate::assert_matrix(
    object$time_end,
    ncols = ncol(object$contact_reduction), nrows = 1L
  )

  # stricter initialisation of interventions so that negative values and
  # npi intervals are not allowed, and cumulative contact reductions > 1 are
  # not allowed
  stopifnot(
    "`nu` should have positive or zero values" =
      all(object$contact_reduction >= 0.0),
    "`time_begin` should have positive or zero values" =
      all(object$time_begin >= 0.0),
    "`time_end` should have values greater-than or equal-to `time_begin`" =
      all(object$time_end >= object$time_begin),
    "Rows of `contact_reduction` must sum to <= 1.0" =
      all(rowSums(object$contact_reduction) <= 1.0)
  )

  # message if any intervention intervals are badly formed
  # tackles the case of mistakenly setting all values the same
  # this is explicitly used in no_intervention(), with message suppressed
  # also accounts for eventual extension to group-specific start and end times
  if (any(object$time_end <= object$time_begin)) {
    message(
      "Intervention: some `time_end`s are not greater than `time_begin`s"
    )
  }

  # checks on length of `contact_reduction` can only be made in the context
  # of a population, see assert_intervention()

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
#'   contact_reduction = matrix(c(0.5, 0.01)) # reduces contacts differentially
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
  # message on identical value of time_begin and time_end suppressed
  # as this is a valid use case.
  suppressMessages(
    intervention(
      name = "no_intervention", time_begin = 0, time_end = 0,
      contact_reduction = matrix(0.0, nrow(population$contact_matrix))
    )
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
      name
    )
  )
  writeLines("Time begin:")
  print(x$time_begin)

  writeLines("Time end:")
  print(x$time_end)

  writeLines("Contact reduction:")
  print(x$contact_reduction)

  invisible(x)
}

#' Convert a list to a intervention object
#'
#' @param x A list, or an object that inherits from a list.
#' @return A [intervention] class object.
#' @export
#' @examples
#' # prepare a list
#' npi <- list(
#'   name = "npi",
#'   time_begin = 30,
#'   time_end = 60,
#'   contact_reduction = matrix(0.1, 3)
#' )
#'
#' as.intervention(npi)
as.intervention <- function(x) {
  # check that input is a list or intervention
  stopifnot(
    "Input must inherit from `list`" =
      is.list(x)
  )

  x <- intervention(
    name = x[["name"]],
    time_begin = x[["time_begin"]],
    time_end = x[["time_end"]],
    contact_reduction = x[["contact_reduction"]]
  )

  # return x
  x
}

#' Concatenate intervention doses into a multi-dose intervention
#'
#' @param x An `intervention` object.
#' @param ... intervention objects to combine with `x` to create a multi-dose
#' `intervention` object.
#' @return An `intervention` object with as many doses as the overall number of
#' doses specified in `x` and in the objects passed to `...`.
#' @export
#' @examples
#' # create first dose regime
#' npi_1 <- intervention(
#'   time_begin = 30,
#'   time_end = 60,
#'   contact_reduction = matrix(0.1)
#' )
#'
#' # second dose regime
#' npi_2 <- intervention(
#'   time_begin = 45,
#'   time_end = 75,
#'   contact_reduction = matrix(0.1)
#' )
#'
#' c(npi_1, npi_2)
c.intervention <- function(x, ...) {
  # collect inputs
  multi_npi <- list(x, ...)
  invisible(
    lapply(multi_npi, validate_intervention)
  )

  # check that all intervention regimes have the same dimensions
  # of intervention rates --- these are identical to dims of start and end times
  stopifnot(
    "All `intervention`s must have identical dimensions for c, start, and end" =
      all(
        vapply(multi_npi, function(npi) {
          identical(nrow(npi$contact_reduction), nrow(x$contact_reduction))
        }, FUN.VALUE = logical(1))
      )
  )

  # strip class and `name` member from `multi_npi`
  multi_npi <- lapply(multi_npi, function(z) {
    z <- unclass(z)
    z$name <- NULL
    z
  })

  # modify x to return a list object of multiple start and end times and nu-s
  multi_npi <- do.call(
    Map, c(f = cbind, multi_npi)
  )

  # add name parameter --- take "name" of `x`
  multi_npi$name <- x$name

  # generate intervention dose names as dose_1 ... dose_n
  # get total number of doses from the sum of all columns
  npi_names <- glue::glue("npi_{seq_len(ncol(multi_npi$contact_reduction))}")

  # add names to doses for comprehension when printed
  for (i in c("time_begin", "time_end", "contact_reduction")) {
    if (is.matrix(multi_npi[[i]])) {
      colnames(multi_npi[[i]]) <- npi_names
    } else {
      names(multi_npi[[i]]) <- npi_names
    }
  }

  # convert resulting object to intervention
  multi_npi <- new_intervention(
    multi_npi$name,
    multi_npi$time_begin,
    multi_npi$time_end,
    multi_npi$contact_reduction
  )

  # validate new object
  validate_intervention(multi_npi)

  # return object
  multi_npi
}
