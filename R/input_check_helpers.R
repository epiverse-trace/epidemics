#' Assert properties of a `population` object
#'
#' @description
#' Assert that objects of the `population` class have the parameters expected by
#' an epidemic model. See [population()] and specific epidemic functions
#' to check the population parameters required by each model. This
#' function is for internal use in argument checking functions.
#'
#' @param x A `<population>` object.
#' @param compartments A character vector giving the names of model compartments
#' whose length is taken as the reference for the number of columns in the
#' `initial_conditions` matrix in `x`.
#' @param demography_vector An optional numeric vector whose length is used to
#' check the length of the demography vector present in `x`.
#'
#' @keywords internal
#'
#' @return Silently returns the `<population>` object `x`.
#' Primarily called for its side effects of throwing errors when `x` does not
#' meet certain requirements.
assert_population <- function(x, compartments, demography_vector = NULL) {
  # check for input class
  checkmate::assert_class(x, "population")

  # check that population has a set number of demography groups
  checkmate::assert_numeric(
    get_parameter(x, "demography_vector"),
    len = demography_vector # checked against NULL in most models
  )

  # check that population has as many compartments in initial conditions
  # matrix as the length of `compartments`
  checkmate::assert_matrix(
    get_parameter(x, "initial_conditions"),
    mode = "numeric", # this is also checked when initialising a population
    nrows = demography_vector, # this is NULL for most models
    ncols = length(compartments)
  )

  # invisibly return x
  invisible(x)
}

#' Assert properties of a `vaccination` object
#'
#' @description
#' Assert that objects of the `vaccination` class have the parameters expected
#' by an epidemic model. See [vaccination()] and specific epidemic functions
#' to check the vaccination properties required by each model. This
#' function is for internal use in argument checking functions.
#'
#' @param x A [vaccination] object.
#' @param doses The number of doses expected in the vaccination object.
#' @param population An optional argument which is a `<population>` object.
#' When present, this is used to check whether the vaccination object `x` has
#' corresponding values for each demographic group in `population`.
#'
#' @keywords internal
#'
#' @return Silently returns the `<vaccination>` object `x`.
#' Primarily called for its side effects of throwing errors when `x` does not
#' meet certain requirements.
assert_vaccination <- function(x, doses, population = NULL) {
  # check for input class
  checkmate::assert_class(x, "vaccination")
  checkmate::assert_number(doses, finite = TRUE, lower = 1L)
  checkmate::assert_integerish(doses, lower = 1L)

  # check that x has as many cols in `nu` as `doses`
  # all other elements are identical dims as `nu`
  checkmate::assert_matrix(
    get_parameter(x, "nu"),
    ncols = doses
  )

  # if a population is provided, check that the rows of `nu`
  # match the number of demography groups
  if (!is.null(population)) {
    checkmate::assert_class(population, "population")
    checkmate::assert_matrix(
      get_parameter(x, "nu"),
      nrows = length(get_parameter(population, "demography_vector"))
    )
  }

  # invisibly return x
  invisible(x)
}

#' Assert properties of a `intervention` object
#'
#' @description
#' Assert that objects of the `intervention` class have the properties expected
#' by an epidemic model. See [intervention()] and specific model functions
#' to check the intervention properties required by each model. This
#' function is for internal use in argument checking functions.
#'
#' @param x A [intervention] object.
#' @param type A string for the type of intervention to check for. May be one of
#' `"contacts"` or `"rate"`.
#' @param population An optional argument which is a [population] object.
#' When present, this is used to check whether the intervention object `x` has
#' corresponding values of `reduction` for each demographic group in
#' `population`.
#'
#' @keywords internal
#'
#' @return Silently returns the `<intervention>`-superclass object `x` of the
#' same class as `x`, i.e., `<rate_intervention>` or `<contacts_intervention>`.
#' Primarily called for its side effects of throwing errors when `x` does not
#' meet certain requirements.
assert_intervention <- function(x, type = c("contacts", "rate"),
                                population = NULL) {
  # check for input class
  type <- match.arg(type, several.ok = FALSE)
  switch(type,
    contacts = {
      checkmate::assert_class(x, "contacts_intervention")
      # check that the length of reduction is the same as the number
      # of demographic groups
      n_demo_groups <- NULL
      if (!is.null(population)) {
        n_demo_groups <- length(get_parameter(population, "demography_vector"))
      }
      checkmate::assert_matrix(
        get_parameter(x, "reduction"),
        mode = "numeric",
        nrows = n_demo_groups
      )
    },
    rate = {
      checkmate::assert_class(x, "rate_intervention")
      checkmate::assert_numeric(
        get_parameter(x, "reduction"),
        lower = 0.0, upper = 1.0,
        len = length(get_parameter(x, "time_begin"))
      )
    }
  )

  # invisibly return x
  invisible(x)
}
