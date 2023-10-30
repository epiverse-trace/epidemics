#' Construct a new infection
#'
#' @param name String for the name of the infection. The default value is `NA`.
#' @param r0 The basic reproductive number of the infection, \eqn{R_0}.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' @param infectious_period The mean infections period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' Currently, only the option to pass a single infectious period common to all
#' demographic groups is available.
#' @param ... Other arguments representing other infection parameters. This
#' argument is intentionally left flexible to allow for single values or vectors
#' of values.
#' @return An `<infection>` class object.
#' @keywords internal
#' @noRd
new_infection <- function(name = NA_character_,
                          r0,
                          infectious_period,
                          extra_arguments) {
  # create and return infection class
  structure(
    append(
      list(
        name = name,
        r0 = r0,
        infectious_period = infectious_period
      ),
      extra_arguments
    ),
    class = "infection"
  )
}

#' Prepare an infection for an epidemic model
#'
#' @name infection
#' @rdname infection
#'
#' @description Prepare an S3 object of the `<infection>` class, which is used
#' to hold some basic parameters of the infection that is to be modelled as
#' causing an epidemic.
#'
#' `<infection>`s can hold any number of extra parameters as required by
#' specific epidemic models.
#'
#' @param name String for the name of the infection. The default value is `NA`.
#' @param r0 The basic reproductive number of the infection, \eqn{R_0}.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' The 'default' model requires a single number.
#' @param infectious_period The mean infections period.
#' Must be a number or a vector of numbers of the same length as the number of
#' demographic groups, depending on the model being implemented.
#' Currently, only the option to pass a single infectious period common to all
#' demographic groups is available.
#' @param ... Other arguments representing other infection parameters. This
#' argument is intentionally left flexible to allow for single values or vectors
#' of values.
#'
#' @param x An object to be checked as being an `<infection>`.
#' @return
#' `infection()` returns an object of the `<infection>` S3 class.
#'
#' `is_infection()` returns a logical for whether the object is of the
#' `<infection>` class.
#' @export
#'
#' @examples
#' pandemic_influenza <- infection(
#'   name = "pandemic influenza",
#'   r0 = 1.3,
#'   infectious_period = 5,
#'   preinfectious_period = 2
#' )
#'
#' # print to examine
#' pandemic_influenza
#'
#' # check for the infection class
#' is_infection(pandemic_influenza)
infection <- function(name = NA_character_,
                      r0,
                      infectious_period,
                      ...) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_number(r0, lower = 0, finite = TRUE)
  # currently allowing only single infectious period
  checkmate::assert_number(infectious_period, lower = 0, finite = TRUE)

  # collect extra infection arguments
  extra_args <- list(...)
  # check types
  checkmate::assert_list(
    extra_args,
    types = numeric()
  )
  # expect that all parameters have the same length
  for (i in seq_along(extra_args)) {
    stopifnot(
      "Error: All infection parameters must be the same length!" =
        length(extra_args[[i]]) == length(infectious_period)
    )
  }

  # call infection constructor
  infection_ <- new_infection(
    name = name,
    r0 = r0,
    infectious_period = infectious_period,
    extra_arguments = extra_args
  )

  # call infection validator
  validate_infection(object = infection_)

  # return infection object
  infection_
}

#' Validate an infection
#'
#' @param object An object to be validated as an `<infection>`.
#'
#' @return No return.
#' @noRd
#' @keywords internal
validate_infection <- function(object) {
  # check for class and class invariants
  stopifnot(
    "Object should be of class `<infection>`" =
      (is_infection(object)),
    "`<infection>` does not contain the minimum correct attributes" =
      (c(
        "name", "r0", "infectious_period"
      ) %in% attributes(object)$names)
  )

  # check infection class members
  checkmate::assert_string(object$name, na.ok = TRUE)
  checkmate::assert_number(object$r0, lower = 0, finite = TRUE)
  checkmate::assert_number(object$infectious_period, lower = 0, finite = TRUE)

  # check for extra arguments
  extra_args <- object[!names(object) %in% c(
    "name", "r0", "infectious_period"
  )]

  # expect that all parameters have the same length
  for (i in seq_along(extra_args)) {
    stopifnot(
      "Error: All infection parameters must be the same length!" =
        length(extra_args[[i]]) == length(object$infectious_period)
    )
  }

  invisible(object)
}

#' Check whether an object is an `<infection>`
#' @name infection
#' @rdname infection
#' @export
is_infection <- function(x) {
  inherits(x, "infection")
}

#' Print an `<infection>` object
#'
#' @param x A `<infection>` object.
#' @param ... Other parameters passed to [print()].
#' @return Invisibly returns the `<infection>` object `x`.
#' Called for printing side-effects.
#' @export
print.infection <- function(x, ...) {
  format(x, ...)
}

#' Format an `<infection>` object
#'
#' @param x A `<infection>` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisibly returns the [`<infection>`] object `x`.
#' Called for printing side-effects.
#' @keywords internal
#' @noRd
format.infection <- function(x, ...) {
  # validate the infection object
  validate_infection(x)

  # header
  header <- class(x) #nolint

  # collect information on name
  #nolint start
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  #nolint end
  # other parameters
  extra_argument_names <- names(x)[!names(x)
  %in% c("name", "r0", "infectious_period")]
  extra_argument_names <- if (length(extra_argument_names) == 0) {
    "None"
    } else {
      glue::glue_collapse(
      glue::double_quote(extra_argument_names),
      sep = ", "
      )
    }

  # print to screen
  cli::cli_text(
    cli::col_red("{.cls {header}} object")
  )
  cat(
    "\n",
    cli::col_red(
      "Name: "
      )
    )
  cli::cli_text(
    "{cli::cli_format({name}, style = list(string_quote = \"\"))}"
  )
  cat(
    "\n",
    cli::col_red(
      "R0: "
      )
  )
  cli::cli_text(
    "{.val {x$r0}}"
    )
  cat(
    "\n",
    cli::col_red(
      "Infectious period: "
      )
  )
  cli::cli_text(
    "{.val {x$infectious_period}}"
    )
  cat(
    "\n",
    cli::col_red(
      "Other parameters: "
      )
  )
  cli::cli_text(
    "{cli::cli_format({extra_argument_names},
    style = list(string_quote = \"\"))}
    "
    )

  invisible(x)
}

#' Calculate the transmission rate of an `<infection>`
#'
#' @param infection An `<infection>` object, with the parameters `r0` for the
#' basic reproductive number, and `infectious_period`, for the mean infectious
#' period.
#'
#' @return A single number representing the transmission rate of the infection,
#' \eqn{\beta}.
#' @export
#'
#' @examples
#' # create an infection
#' ebola <- infection(r0 = 1.3, infectious_period = 5, preinfectious_period = 5)
#' get_transmission_rate(infection = ebola)
get_transmission_rate <- function(infection) {
  checkmate::assert_class(infection, "infection")
  get_parameter(infection, "r0") / get_parameter(infection, "infectious_period")
}
