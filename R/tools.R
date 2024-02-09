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
