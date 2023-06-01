#' Assert properties of an `infection` object
#'
#' @description
#' Assert that objects of the `infection` class have the parameters expected by
#' an epidemic model. See [infection()] and [epidemic()], as well as model
#' details to check the infection parameters required by each model. This
#' function is for internal use in argument checking functions.
#'
#' @param x An [infection] object.
#' @param extra_parameters A character vector giving the names of any extra
#' parameters included in `x`. `x` may already be expected to have the members
#' `r0` and `infectious_period`.
#' @param extra_parameters_limits An optional named list of the numeric limits
#' of the acceptable limits of the parameters specifies in `extra_parameters`.
#' Each list element must be a named, two element vector of 'lower' and 'upper'
#' limits. Not all extra parameters need to have their limits passed explicitly.
#' Parameters whose limits are not explicitly specified will have their values
#' checked against a lower bound of 0, with the upper limit expected to be
#' finite.
#'
#' @keywords internal
#'
#' @return Silently returns the `infection` object `x`. Primarily called for its
#' side effects of throwing errors when `x` does not meet certain requirements.
assert_infection <- function(x, extra_parameters, extra_parameters_limits) {
  # check for input class and expected names
  checkmate::assert_class(x, "infection")
  checkmate::assert_names(
    names(x),
    must.include = c(
      "r0", "infectious_period", extra_parameters
    )
  )
  # checks on the r0 and infectious period
  checkmate::assert_number(
    x$r0,
    lower = 0, finite = TRUE
  )
  checkmate::assert_number(
    x$infectious_period,
    lower = 0, finite = TRUE
  )

  # checks on the limits of the extra parameters
  # by default check for finite non-negative parameters
  invisible(
    lapply(
      x[[extra_parameters]], checkmate::assert_number,
      lower = 0.0, finite = TRUE
    )
  )

  # extra checks on limits for extra parameters
  # most development will not use this but there may be
  # specialised models that require limits on some parameters
  # in this case both the upper and lower limits should be passed
  if (!missing(extra_parameters_limits)) {
    # check that extra params limits are a list
    checkmate::assert_list(
      extra_parameters_limits,
      min.len = 1L,
      max.len = length(extra_parameters),
      types = "numeric",
      names = "unique"
    )
    checkmate::assert_names(
      names(extra_parameters_limits),
      subset.of = extra_parameters
    )

    # check that extra_parameters_limits is a two element numeric
    # with names "lower" and "upper"
    invisible(
      lapply(
        extra_parameters_limits, checkmate::assert_numeric,
        len = 2
      )
    )
    invisible(
      lapply(
        names(extra_parameters_limits), checkmate::assert_names,
        c("lower", "upper")
      )
    )

    # now check the elements of x, the infection object, using the limits
    invisible(
      Map(
        x[[names(extra_parameters_limits)]], extra_parameters_limits,
        f = function(z, lims) {
          checkmate::assert_number(
            z,
            lower = lims[["lower"]], upper = lims[["upper"]]
          )
        }
      )
    )
  }
  # invisibly return x
  invisible(x)
}
