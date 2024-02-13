#' Transpose a list
#' @description
#' This is intended to be a base R implementation of the `transpose()`
#' function from the `{purrr}` package.
#' @author Tim Taylor
#' @param x A list of vectors or list-like elements, which have the same length.
#' @return A list of the same length as each element of `x`, where each element
#' is of the length of `x`. When `x` is named, each element of the returned list
#' takes the name of the original element of `x`.
#' @keywords internal
.transpose_base <- function(x) {
  stopifnot(
    "List elements must all be vectors or lists" =
      all(vapply(x, is.vector, TRUE)),
    "List elements must be of the same length" =
      length(unique(lengths(x))) == 1
  )
  lapply(seq_along(x[[1]]), function(y) lapply(x, `[[`, y))
}

#' Check compliance with Tidyverse vector recycling rules
#' @author Tim Taylor
#' @param x A list of vectors to be checked for compliance with Tidyverse
#' recycling rules.
#' @return A single logical value.
#' @keywords internal
test_recyclable <- function(x) {
  # basic input checking only
  stopifnot(
    "List elements must all be vectors or lists" =
      all(vapply(x, is.vector, TRUE))
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
#'
#' @param x A list of vectors to be recycled. Two cases are envisaged: all
#' elements of `x` are expected to be of the same length, or, all but one
#' element of `x` are scalars.
#' @return A list of vectors with all elements recycled to the same length.
#' @details
#' Note that this function will return vectors of unequal lengths if `x` does
#' not comply with length rules. This compliance is not enforced as this
#' function is only expected to be used internally.
#' @keywords internal
.recycle_vectors <- function(x) {
  stopifnot(
    "List elements must all be vectors" =
      all(vapply(x, is.vector, TRUE))
  )
  # vector lengths not checked as this fn is expected to be used
  # after a call to `test_recyclable()`
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

#' Cross-check an intervention list
#'
#' @param x A named list of `<intervention>` objects.
#' @param population A `<population>` to which any `<contact_intervention>` in
#' `x` applies, and with which it must be compatible.
#' @param allowed_targets A character vector of targets allowed for the
#' `<intervention>`s in `x`. May contain "contacts", and the names of any
#' infection parameters.
#'
#' @return A named list with at least the elements "contacts" describing a
#' `<contacts_intervention>` on `population`, and a `<rate_intervention>` on
#' the transmissibility parameter. If these are present in `x`, they are
#' returned as is, or substituted if missing such that the function output is
#' suitable for processing for a C++ model function. Any other interventions
#' are also returned. If `x` is `NULL`, dummy
#' contact and rate interventions are returned in a list.
#' @keywords internal
.cross_check_intervention <- function(x, population, allowed_targets) {
  # create dummy intervention set
  tmp_intervention <- list(
    contacts = no_contacts_intervention(population),
    transmissibility = no_rate_intervention()
  )
  if (is.null(x)) {
    return(tmp_intervention)
  }

  # check that interventions are on allowed targets
  # these are typically "contacts" and some infection parameters
  checkmate::assert_names(names(x), subset.of = allowed_targets)

  # check that contact interventions are suitable for population
  if ("contacts" %in% names(x)) {
    assert_intervention(x[["contacts"]], "contacts", population)
  }

  # replace dummy values with user values if avaialable, and return
  tmp_intervention[names(x)] <- x
  tmp_intervention
}

#' Cross-check a `<vaccination>`
#'
#' @param x A `<vaccination>`
#' @param population A `<population>`
#' @param doses A single number for the expected number of doses.
#'
#' @return Returns `x` after checking that it is suitable for `population`, or
#' a dummy vaccination regime with `doses` number of doses for each age group.
#' @keywords internal
.cross_check_vaccination <- function(x, population, doses) {
  # no input checking as this is an internal function
  if (is.null(x)) {
    no_vaccination(population, doses = doses)
  } else {
    assert_vaccination(x, doses = doses, population)
    x
  }
}

#' Cross-check time-dependence
#'
#' @param x A named list of functions that describe how an infection parameter
#' changes with time.
#' @param allowed_targets A character vector of targets for `x`; these are
#' typically infection parameters.
#'
#' @return If `x` is not `NULL`, returns `x`; otherwise a dummy function is
#' returned.
#' @keywords internal
.cross_check_timedep <- function(x, allowed_targets) {
  if (is.null(x)) {
    no_time_dependence()
  } else {
    checkmate::assert_names(names(x), subset.of = allowed_targets)
    x
  }
}

.cross_check_popchange <- function(x, population) {
  x
}
