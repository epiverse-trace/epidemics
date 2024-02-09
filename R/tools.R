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
