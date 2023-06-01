#' Assert properties of an `infection` object
#'
#' @description
#' Assert that objects of the `infection` class have the parameters expected by
#' an epidemic model. See [infection()] and [epidemic()], as well as model
#' details to check the infection parameters required by each model. This
#' function is for internal use in argument checking functions.
#'
#' @param x An [infection] object.
#' @param default_params Default parameter names present in an
#' `infection` object.
#' This argument is provided to account for potential changes to the default
#' elements of an `infection` object, and is not expected to be changed.
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
#' @examples \dontrun{
#' # prepare a well formed infection object for the default model
#' infection_default <- infection(
#'   name = "influenza", r0 = 1.3, infectious_period = 10,
#'   preinfectious_period = 3, mortality_rate = 1e-4
#' )
#'
#' # expect no errors for well formed infection-assertion matches
#' assert_infection(
#'   infection_default,
#'   extra_parameters = c("preinfectious_period", "mortality_rate")
#' )
#' }
assert_infection <- function(x,
                             default_params = c(
                               "name", "r0", "infectious_period"
                             ),
                             extra_parameters, extra_parameters_limits) {
  # check for input class and expected names
  checkmate::assert_class(x, "infection")

  # check that there are no extra parameters other than the ones specified
  # collect names other than default names
  infection_extra_params <- setdiff(names(x), default_params)
  checkmate::assert_names(
    infection_extra_params,
    identical.to = extra_parameters
  )
  # name, r0, and infectious period are checked when initialising
  # an infection object, via validate_infection

  # checks on the limits of the extra parameters
  # by default check for finite non-negative parameters
  invisible(
    lapply(
      x[extra_parameters], checkmate::assert_number,
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
        extra_parameters_limits, function(le) {
          checkmate::assert_names(
            names(le),
            identical.to = c("lower", "upper")
          )
        }
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
