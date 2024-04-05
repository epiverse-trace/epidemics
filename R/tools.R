#' Vector recyclability checks and vector recycling
#' @name vector_recycling
#' @rdname vector_recycling
#'
#' @description Internal functions to check whether vectors can be recycled and
#' to recycle vectors.
#'
#' @author Tim Taylor
#' @param x A list of vectors to be checked for compliance with Tidyverse
#' recycling rules, or to be recycled.
#' When `x` is a list to be recycled, two cases are envisaged: all
#' elements of `x` are expected to be of the same length, or, all but one
#' element of `x` are scalars.
#' @return
#' `.test_recyclable` returns a single logical value.
#' `.recycle_vectors` returns a list of vectors recycled to the same length
#' following Tidyverse recycling rules.
#' @details
#' Note that `.test_recyclable()` will return vectors of unequal lengths if `x`
#' does not comply with length rules. This compliance is not enforced as this
#' function is only expected to be used internally, after a call to
#' `.test_recyclable()` (e.g. in a `stopifnot()`).
#' @keywords internal
.test_recyclable <- function(x) {
  # basic input checking only
  stopifnot(
    "`x` must be a list with vector elements" =
      is.list(x) && all(vapply(x, is.vector, TRUE))
  )
  lens <- lengths(x)
  not_scalar <- lens != 1L
  if (any(not_scalar)) {
    lens <- lens[not_scalar]
    # return logical
    length(unique(lens)) == 1L
  } else {
    TRUE
  }
}

#' Recycle vectors with Tidyverse recycling rules
#' @name vector_recycling
#' @keywords internal
.recycle_vectors <- function(x) {
  # basic input checking only
  stopifnot(
    "`x` must be a list with vector elements" =
      is.list(x) && all(vapply(x, is.vector, TRUE))
  )
  # vector lengths not checked as this fn is expected to be used
  # after a call to `.test_recyclable()`
  vec_lengths <- lengths(x)
  longest <- max(vec_lengths)
  x[vec_lengths != longest] <- lapply(
    x[vec_lengths != longest], rep, longest
  )

  x
}

#' Prepare a `<population>` for an epidemic model
#'
#' @param x A `<population>`.
#'
#' @return A named list of "contact_matrix", the population social contacts
#' matrix scaled by the largest real eigenvalue and the demography vector, and
#' "initial_state", the proportional initial state multiplied by the demography
#' vector.
#' @keywords internal
.prepare_population <- function(x) {
  # no input checks on internal function

  # prepare the contact matrix and the initial conditions
  # scale the contact matrix by the maximum real eigenvalue
  # scale rows of the contact matrix by the corresponding group population
  contact_matrix <- x[["contact_matrix"]]
  contact_matrix <- (contact_matrix / max(Re(eigen(contact_matrix)$values))) /
    x[["demography_vector"]]

  # prepare initial conditions by scaling with demography
  initial_state <- x[["initial_conditions"]] * x[["demography_vector"]]

  # return list
  list(contact_matrix = contact_matrix, initial_state = initial_state)
}

#' Cross-check model elements
#'
#' @name cross_checking_inputs
#' @rdname cross_checking_inputs
#'
#' @description
#' Check model elements for compatibility with the population in an epidemic
#' model, returning compatible dummy values when model elements are not applied,
#' and erroring appropriately when model elements are not compatible with the
#' population characteristics.
#'
#' @inheritParams model_default
#' @param x Model input to be checked. The expected value of `x` depends on the
#' function:
#'
#' - `.cross_check_intervention()`: A named list of `<intervention>` objects;
#'
#' - `.cross_check_vaccination()`: A `<vaccination>` object;
#'
#' - `.cross_check_timedep()`: A named list of functions with two arguments,
#' `time` and `x`, typically returning `x` as a function of `time`;
#'
#' - `.cross_check_popchange()`: A named list with two elements, `time` and
#' `values`, describing the times and values by which the number of susceptibles
#' changes in an epidemic model.
#' @param allowed_targets The model components, or infection parameters, that
#' the model input `x` affects.
#' @param doses The expected number of vaccination doses.
#'
#' @return
#' - `.cross_check_intervention()` returns a named list with at least the
#' elements "contacts" describing a `<contacts_intervention>` on `population`
#' (if this is among the allowed targets), and a `<rate_intervention>` on the
#' transmission rate parameter. If these are present in `x`, they are
#' returned as is, or substituted if missing. Any other interventions
#' are also returned. If `x` is `NULL`, dummy contact and rate interventions
#' are returned in a list.
#'
#' - `.cross_check_vaccination()` returns `x` after checking that it is suitable
#' for `population`, or a dummy vaccination regime with `doses` number of doses
#' for each age group.
#'
#' - `.cross_check_timedep()` returns `x` if `x` is not `NULL`, otherwise
#' returns a dummy function operating on the transmission rate parameter by
#' default; see [.no_time_dependence()];
#'
#' - `.cross_check_popchange()` returns `x` after checks against `population` if
#' `x` is not `NULL`, otherwise returns a dummy list with no population change;
#' see [.no_population_change()].
#' @keywords internal
.cross_check_intervention <- function(x, population, allowed_targets) {
  # create dummy intervention set
  tmp_intervention <- list(
    transmission_rate = .no_rate_intervention()
  )
  # Ebola and Diphtheria models do not allow contact interventions
  if ("contacts" %in% allowed_targets) {
    tmp_intervention[["contacts"]] <- .no_contacts_intervention(population)
  }
  if (is.null(x)) {
    tmp_intervention
  } else {
    # check that interventions are on allowed targets
    # these are typically "contacts" and some infection parameters
    checkmate::assert_names(names(x), subset.of = allowed_targets)

    # check that contact interventions are suitable for population
    if ("contacts" %in% names(x)) {
      assert_intervention(x[["contacts"]], "contacts", population)
    }
    # check class of other intervention objects
    invisible(
      lapply(x[names(x) != "contacts"], assert_intervention, type = "rate")
    )

    # replace dummy values with user values if avaialable, and return
    tmp_intervention[names(x)] <- x
    tmp_intervention
  }
}

#' Cross-check a `<vaccination>`
#' @name cross_checking_inputs
#' @rdname cross_checking_inputs
#' @keywords internal
.cross_check_vaccination <- function(x, population, doses) {
  # no input checking as this is an internal function
  if (is.null(x)) {
    .no_vaccination(population, doses = doses)
  } else {
    assert_vaccination(x, doses = doses, population)
    x
  }
}

#' Cross-check time-dependence
#' @name cross_checking_inputs
#' @rdname cross_checking_inputs
#' @keywords internal
.cross_check_timedep <- function(x, allowed_targets) {
  if (is.null(x)) {
    .no_time_dependence()
  } else {
    checkmate::assert_names(names(x), subset.of = allowed_targets)
    x
  }
}

#' Cross-check population change
#' @name cross_checking_inputs
#' @rdname cross_checking_inputs
#' @keywords internal
.cross_check_popchange <- function(x, population) {
  if (is.null(x)) {
    .no_population_change(population)
  } else {
    checkmate::assert_list(
      x,
      any.missing = FALSE, names = "unique",
      len = 2L, types = c("numeric", "list")
    )
    checkmate::assert_names(
      names(x),
      identical.to = c("time", "values")
    )
    # check that time vector and values list have identical lengths
    checkmate::assert_numeric(
      x[["time"]],
      lower = 0, finite = TRUE, min.len = 1
    )
    checkmate::assert_list(
      x[["values"]],
      any.missing = FALSE, len = length(x[["time"]])
    )
    # check that values elements (vecs) are compatible with population
    invisible(
      lapply(
        x[["values"]],
        FUN = function(le) {
          stopifnot(
            "`population_change` `values` must be same length as demography" =
              checkmate::test_numeric(
                le,
                len = length(population[["demography_vector"]]),
                any.missing = FALSE, finite = TRUE
              )
          )
        }
      )
    )

    # return x
    x
  }
}
