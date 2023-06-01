
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
  checkmate::assert_matrix(nu, mode = "numeric")
  checkmate::assert_matrix(
    time_begin,
    nrows = nrow(nu), ncols = ncol(nu), mode = "numeric"
  )
  checkmate::assert_matrix(
    time_end,
    nrows = nrow(nu), ncols = ncol(nu), mode = "numeric"
  )

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
  checkmate::assert_matrix(object$nu, mode = "numeric")
  checkmate::assert_matrix(
    object$time_begin,
    nrows = nrow(object$nu), ncols = ncol(object$nu),
    mode = "numeric"
  )
  checkmate::assert_matrix(
    object$time_end,
    nrows = nrow(object$nu), ncols = ncol(object$nu),
    mode = "numeric"
  )
  stopifnot(
    "`nu` should have positive or zero values" =
      all(object$nu >= 0.0),
    "`time_begin` should have positive or zero values" =
      all(object$time_begin >= 0.0),
    "`time_end` should have positive or zero values" =
      all(object$time_end >= 0.0)
  )

  # message if any vaccinations' intervals are badly formed
  if (any(object$time_end <= object$time_begin)) {
    message(
      "Vaccination: some `time_end`s are not greater than `time_begin`s"
    )
  }

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
  # suppress messages because we are sure we want both time begin and end as 0
  suppressMessages(
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

  writeLines(
    glue::glue(
      "

    Time end:
    "
    )
  )
  print(x$time_end)

  writeLines(
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

#' Concatenate vaccination doses into a multi-dose vaccination
#'
#' @param x A `vaccination` object.
#' @param ... Vaccination objects to combine with `x` to create a multi-dose
#' `vaccination` object.
#' @return A `vaccination` object with as many doses as the overall number of
#' doses specified in `x` and in the objects passed to `...`.
#' @export
#' @examples
#' # create first dose regime
#' vax_1 <- vaccination(
#'   name = "vax_regime",
#'   time_begin = matrix(1),
#'   time_end = matrix(100),
#'   nu = matrix(0.001)
#' )
#'
#' # second dose regime
#' vax_2 <- vaccination(
#'   name = "vax_regime",
#'   time_begin = matrix(101),
#'   time_end = matrix(200),
#'   nu = matrix(0.001)
#' )
#'
#' c(vax_1, vax_2)
c.vaccination <- function(x, ...) {
  # collect inputs
  multi_vacc <- list(x, ...)
  invisible(
    lapply(multi_vacc, validate_vaccination)
  )

  # check that all vaccination regimes have the same dimensions
  # of vaccination rates --- these are identical to dims of start and end times
  stopifnot(
    "All `vaccination`s must have identical dimensions for Nu, start, and end" =
      all(
        vapply(multi_vacc, function(vx) {
          identical(nrow(vx$nu), nrow(x$nu))
        }, FUN.VALUE = logical(1))
      )
  )

  # strip class and `name` member from `multi_vacc`
  multi_vacc <- lapply(multi_vacc, function(z) {
    z <- unclass(z)
    z$name <- NULL
    z
  })

  # modify x to return a list object of multiple start and end times and nu-s
  multi_vacc <- do.call(
    Map, c(f = cbind, multi_vacc)
  )

  # add name parameter --- take "name" of `x`
  multi_vacc$name <- x$name

  # generate vaccination dose names as dose_1 ... dose_n
  # get total number of doses from the sum of all columns
  vacc_dose_names <- glue::glue("dose_{seq_len(ncol(multi_vacc$nu))}")

  # add names to doses for comprehension when printed
  for (i in c("time_begin", "time_end", "nu")) {
    colnames(multi_vacc[[i]]) <- vacc_dose_names
  }

  # convert resulting object to vaccination
  multi_vacc <- as.vaccination(multi_vacc)

  # validate new object
  validate_vaccination(multi_vacc)

  # return object
  multi_vacc
}
