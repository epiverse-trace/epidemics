#' Constructor for the <intervention> super-class and sub-classes
#'
#' @name intervention_constructor
#' @rdname intervention_constructor
#'
#' @param name String for the name of the intervention.
#' @param time_begin A matrix with a single element giving the start time of the
#' intervention.
#' @param time_end A matrix with a single element giving the end time of the
#' intervention.
#' @param reduction A matrix with a single column and as many rows as there are
#' demographic groups to be targeted by the intervention. See [intervention()]
#' for details on the types of intervention that can be created and the
#' requirements for the type of `reduction`.
#' @param ... Any other parameters to be passed to the constructor.
#' @param class A string giving the type of the intervention; used to generate
#' intervention sub-classes.
#'
#' @return
#' `new_intervention()` returns an object that inherits from the
#' `<intervention>` class.
#' `new_contacts_intervention()` returns an object that is a sub-class of
#' `<intervention>` called `<contacts_intervention>`.
#' `new_rate_intervention()` returns an object that is a sub-class of
#' `<intervention>` called `<rate_intervention>`.
#' @keywords internal
new_intervention <- function(name = NA_character_,
                             time_begin,
                             time_end,
                             reduction,
                             ...,
                             class) {
  # argument class is empty to force interventions to have a sub-class

  # create and return intervention class
  x <- list(
    name = name,
    time_begin = time_begin,
    time_end = time_end,
    reduction = reduction,
    ...
  )
  class(x) <- c(class, "intervention")

  x
}

#' Constructor for a new <contacts_intervention>
#'
#' @name intervention_constructor
#' @rdname intervention_constructor
new_contacts_intervention <- function(name, time_begin, time_end,
                                      reduction) {
  new_intervention(
    name, time_begin, time_end,
    reduction,
    class = "contacts_intervention"
  )
}

#' Constructor for a new <rate_intervention>
#'
#' @name intervention_constructor
#' @rdname intervention_constructor
new_rate_intervention <- function(name, time_begin, time_end,
                                  reduction) {
  new_intervention(
    name, time_begin, time_end,
    reduction,
    class = "rate_intervention"
  )
}

#' Create an intervention for an epidemic model
#'
#' @name intervention
#' @rdname intervention
#'
#' @description
#' Prepare an object of the `<intervention>` super-class that specifies a
#' modification of the model parameters.
#'
#' A `<contacts_intervention>` is used to simulate a non-pharmaceutical
#' intervention (NPI) regime that reduces the population's social contacts.
#'
#' A `<rate_intervention>` is used to simulate a reduction in the model's rate
#' parameters (such as the transmission rate \eqn{\beta}), and can be used to
#' represent pharmaceutical interventions such as improved treatment, but also
#' NPIs such as wearing masks.
#'
#' Interventions have a single start and end time that applies to all
#' demographic groups in the population, but can have groups-specific effects on
#' the reduction of contacts.
#'
#' Combine `<intervention>`-inheriting objects to create sequential or
#' overlapping intervention regimes using `c()` on two or more
#' `<intervention>`-inheriting objects.
#'
#' @param name String for the name of the intervention.
#' @param type String for the type of intervention. May be one of `"contacts"`
#' or `"rate"`, for a `<contacts_intervention>` or `<rate_intervention>`
#' respectively.
#' @param time_begin Single number for the start time of the intervention.
#' @param time_end Single number for the end time of the intervention.
#' @param reduction
#'
#' For `<contacts_intervention>`s, a matrix with as many rows as the
#' number of demographic groups in the type population, and a single column.
#' Each element gives the group-specific proportion reduction in contacts.
#'
#' For `<rate_intervention>`s, a single number giving the proportion reduction
#' in a model parameter contacts.
#'
#' See details for how `c()` can be used to combine interventions of the same
#' sub-class.
#'
#' @param x An `<intervention>` object, or an object to be checked as an
#' `<intervention>` object.
#' @param ... intervention objects to combine with `x` to create a multi-dose
#' `<intervention>` object.
#'
#'
#' @details
#' Epidemic models that can accommodate interventions on contacts are able to
#' accommodate any number of interventions with different start and end times
#' and different group-specific effects.
#'
#' Epidemic models that can accommodate interventions on rates are also able to
#' accommodate any number of interventions with different start and end times,
#' but with only a uniform effect on the relevant rate.
#'
#' When multiple contact interventions are combined using `c()`, the reduction
#' in contacts is stacked column wise to form a matrix \eqn{[i, j]}.
#'
#' When multiple rate interventions are combined using `c()`, the reduction
#' in the rate is concatenated into a vector of the same length as the number of
#' interventions.
#'
#' Models such as model_default() are set up to treat interventions
#' with overlapping periods (i.e., overlap between the time when they are active
#' ) as having an _additive effect_ on contact or rate reductions.
#'
#' For contact reductions, the group-specific effect of \eqn{J} overlapping
#' interventions is thus a vector \eqn{\sum_{j = 1}^J x_{ij}}, for each
#' demographic group \eqn{i}. This is handled internally by the epidemic model
#' code.
#' For example, a contact reduction matrix for two perfectly overlapping
#' interventions (\eqn{J = 2}) with different effects across three demographic
#' groups (\eqn{I = 3}) would be represented as:
#' \eqn{\begin{bmatrix}0.1 & 0.05\\0.1 & 0.1\\0.1 & 0.0\end{bmatrix}}
#' In epidemic models, the cumulative group-specific effect when both
#' interventions are active would be \eqn{(0.15, 0.2, 0.1)}.
#'
#' For rate reductions, the effect of overlapping interventions that reduce a
#' particular rate is also considered to be additive.
#' @return An object of the `<intervention>` S3 super-class, with possible
#' sub-classes `<contact_intervention>` and `<rate_intervention>`.
#'
#' Concatenating two or more `<intervention>`-inheriting objects using `c()`
#' also returns a `<intervention>`-inheriting object of the same sub-class.
#' This object holds the intervention-specific start and end times, and
#' reductions specified by all the constituent intervention actions (by
#' demographic group if an intervention on contacts).
#'
#' The combined effect of these actions on the population is handled internally
#' by epidemic model functions.
#'
#' A "null" intervention generated using `.no_contacts_intervention(population)`
#' or `.no_rate_intervention()` returns a `<intervention>` of the appropriate
#' sub-class that has its start and end times, and its effect all
#' set to 0.0.
#'
#' `is_intervention()`, `is_contacts_intervention()`, and
#' `is_rate_intervention()` each return a logical value for whether the object
#' is of the `<intervention>`, `<contacts_intervention>`, or
#' `<rate_intervention>` class, respectively.
#' @export
#'
#' @examples
#' # assuming a population with two age groups, 0 -- 18, and 18+
#' # an example in which schools are closed for 30 days (or other time units)
#' close_schools <- intervention(
#'   name = "close schools",
#'   type = "contacts",
#'   time_begin = 50,
#'   time_end = 80,
#'   reduction = matrix(c(0.5, 0.01)) # reduces contacts differentially
#' )
#' close_schools
#'
#' # Check for intervention class
#' is_contacts_intervention(close_schools)
#'
#' # Concatenating interventions
#' # create first intervention
#' npi_1 <- intervention(
#'   type = "contacts",
#'   time_begin = 30,
#'   time_end = 60,
#'   reduction = matrix(0.1)
#' )
#'
#' # second intervention
#' npi_2 <- intervention(
#'   type = "contacts",
#'   time_begin = 45,
#'   time_end = 75,
#'   reduction = matrix(0.1)
#' )
#'
#' c(npi_1, npi_2)
intervention <- function(name = NA_character_,
                         type,
                         time_begin,
                         time_end,
                         reduction) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_numeric(reduction, lower = 0, upper = 1.0, finite = TRUE)
  checkmate::assert_number(time_begin, lower = 0, finite = TRUE)
  checkmate::assert_number(time_end, lower = 0, finite = TRUE)

  # check type argument
  type <- match.arg(
    type,
    choices = c("contacts", "rate"),
    several.ok = FALSE
  )

  # message if any intervention intervals are badly formed
  if (time_end <= time_begin) {
    message(
      "Intervention: `time_end` is not greater than `time_begin`"
    )
  }

  # call intervention constructor depending on type
  # validate while returning
  intervention_ <-
    switch(type,
      contacts = new_contacts_intervention(
        name = name,
        time_begin = time_begin,
        time_end = time_end,
        reduction = matrix(reduction) # matrix for contacts reduction
      ),
      rate = new_rate_intervention(
        name = name,
        time_begin = matrix(time_begin),
        time_end = matrix(time_end),
        reduction = reduction
      )
    )

  # return intervention object
  intervention_
}

#' Validate objects that inherit from the <intervention> class
#'
#' @name validate_intervention
#' @rdname validate_intervention
#'
#' @param x For `validate_contacts_intervention()`, an object to be validated as
#' a `<contacts_intervention>`. For `validate_rate_intervention()`, an object to
#' be validated as a `<rate_intervention>`.
#'
#' @return Invisibly returns the input `x`. Called primarily for its
#' input checking as a validator for objects of the `<intervention>` superclass.
#' @keywords internal
validate_contacts_intervention <- function(x) {
  # check for class and class invariants
  stopifnot(
    "`x` should be of class <contacts_intervention>" =
      (is_contacts_intervention(x)),
    "<contacts_intervention> does not contain the correct attributes" =
      (c(
        "name", "time_begin", "time_end", "reduction"
      ) %in% attributes(x)$names)
  )

  # check intervention class members
  checkmate::assert_string(x$name, na.ok = TRUE)

  # checks for conformity of reduction and start and end times
  checkmate::assert_matrix(x$reduction, mode = "numeric")
  checkmate::assert_numeric(
    x$time_begin,
    len = ncol(x$reduction)
  )
  checkmate::assert_numeric(
    x$time_end,
    len = ncol(x$reduction)
  )

  # Checks intervention members so that single effect values are in the range
  # 0 - 1 and negative npi intervals are not allowed
  stopifnot(
    "`reduction` can only have values in the range 0.0 -- 1.0" =
      x$reduction >= 0.0 & x$reduction <= 1.0,
    "`time_begin` should have positive or zero values" =
      x$time_begin >= 0.0,
    "`time_end` should have values greater-than or equal-to `time_begin`" =
      x$time_end >= x$time_begin
  )

  # message if any intervention intervals are badly formed
  # tackles the case of mistakenly setting all values the same
  # this is explicitly used in .no_contacts_intervention(), with message
  # suppressed
  # also accounts for eventual extension to group-specific start and end times
  if (any(x$time_end <= x$time_begin)) {
    message(
      "Intervention: some `time_end`s are not greater than `time_begin`s"
    )
  }

  # checks on the number of rows of `reduction` can only be made in
  # the context of a population, see assert_intervention()

  # invisibly return x
  invisible(x)
}

#' Validate objects that inherit from the <intervention> class
#' @name validate_intervention
#' @rdname validate_intervention
validate_rate_intervention <- function(x) {
  # check for class and class invariants
  stopifnot(
    "`x` should be of class <rate_intervention>" =
      (is_rate_intervention(x)),
    "<rate_intervention> does not contain the correct attributes" =
      (c(
        "name", "time_begin", "time_end", "reduction"
      ) %in% attributes(x)$names)
  )

  # check intervention class members
  checkmate::assert_string(x$name, na.ok = TRUE)

  # checks for conformity of reduction and start and end times
  checkmate::assert_numeric(
    x$reduction,
    lower = 0.0, upper = 1.0, any.missing = FALSE
  )
  checkmate::assert_numeric(
    x$time_begin,
    len = length(x$reduction)
  )
  checkmate::assert_numeric(
    x$time_end,
    len = length(x$reduction)
  )

  # Checks intervention members so that single effect values are in the range
  # 0 - 1 and negative npi intervals are not allowed
  stopifnot(
    "`reduction` can only have values in the range 0.0 -- 1.0" =
      x$reduction >= 0.0 & x$reduction <= 1.0,
    "`time_begin` should have positive or zero values" =
      x$time_begin >= 0.0,
    "`time_end` should have values greater-than or equal-to `time_begin`" =
      x$time_end >= x$time_begin
  )

  # message if any intervention intervals are badly formed
  # tackles the case of mistakenly setting all values the same
  # this is explicitly used in .no_contacts_intervention(), with message
  # suppressed
  # also accounts for eventual extension to group-specific start and end times
  if (any(x$time_end <= x$time_begin)) {
    message(
      "Intervention: some `time_end`s are not greater than `time_begin`s"
    )
  }

  # invisibly return x
  invisible(x)
}

#' Check whether an object is an `<intervention>`
#'
#' @name intervention
#' @rdname intervention
#' @export
is_intervention <- function(x) {
  inherits(x, "intervention")
}

#' Check whether an object is a `<contacts_intervention>`
#' @name intervention
#' @rdname intervention
#'
#' @export
is_contacts_intervention <- function(x) {
  inherits(x, "contacts_intervention")
}

#' Check whether an object is a `<rate_intervention>`
#' @name intervention
#' @rdname intervention
#'
#' @export
is_rate_intervention <- function(x) {
  inherits(x, "rate_intervention")
}

#' Generate a null intervention on contacts
#'
#' @param population A `<population>` for which the dummy contacts intervention
#' is suitable.
#' @return A dummy `<contacts_intervention>` object.
#' @noRd
#' @keywords internal
.no_contacts_intervention <- function(population) {
  checkmate::assert_class(population, "population")
  # message on identical value of time_begin and time_end suppressed
  # as this is a valid use case.
  suppressMessages(
    intervention(
      name = "no_contacts_intervention", type = "contacts",
      time_begin = 0, time_end = 0,
      reduction = rep(0.0, nrow(population$contact_matrix))
    )
  )
}

#' Print an object of the `<intervention>` super-class
#'
#' @name print_intervention
#' @rdname print_intervention
#'
#' @param x An object that inherits from the `<intervention>` class.
#' For `print.contacts_intervention()`, an object of the
#' `<contacts_intervention>` class. For `print.rate_intervention()`, an object
#' of the `<rate_intervention>` class.
#' @param ... Other parameters passed to [print()].
#' @return Invisibly returns the object `x`.
#' Called for printing side-effects.
#' @export
print.intervention <- function(x, ...) {
  format(x, ...)
}

#' Print objects of the `<contact_intervention>` class
#'
#' @name print_intervention
#' @rdname print_intervention
print.contact_intervention <- function(x, ...) {
  validate_contacts_intervention(x)

  format(x, ...)
}

#' Print objects of the `<rate_intervention>` class
#'
#' @name print_intervention
#' @rdname print_intervention

print.rate_intervention <- function(x, ...) {
  validate_rate_intervention(x)

  format(x, ...)
}

#' Format an `<intervention>` object
#'
#' @param x A `<intervention>` object.
#' @param ... Other arguments passed to [format()].
#'
#' @return Invisibly returns the [`<intervention>`] object `x`.
#' Called for printing side-effects.
#' @keywords internal
#' @noRd
format.intervention <- function(x, ...) {
  # header
  header <- class(x)[1] # nolint: object_usage_linter

  # collect information on name
  # nolint start: object_usage_linter
  name <- ifelse(
    is.na(x$name),
    "NA",
    glue::double_quote(x$name)
  )
  # nolint end

  # Prepare reduction effect for printing
  effect <- x$reduction
  if (is.matrix(effect)) {
    colnames(effect) <- glue::glue("Interv. {seq(ncol(effect))}")
    rownames(effect) <- glue::glue("Demo. grp. {seq(nrow(effect))}")
  } else if (is.vector(effect, mode = "numeric")) {
    names(effect) <- glue::glue("Interv. {seq_along(effect)}")
  }

  # print to screen
  cat(
    cli::cli_text(
      "{.cls {header}} object"
    )
  )
  # intervention name
  cat(
    "\n",
    cli::col_magenta(
      "Intervention name: "
    )
  )
  cli::cli_text(
    "{cli::cli_format({name}, style = list(string_quote = \"\"))}"
  )
  # intervention time begin
  cat(
    "\n",
    cli::col_magenta("Begins at:"),
    "\n"
  )
  print(x$time_begin)
  # intervention time end
  cat(
    "\n",
    cli::col_magenta(
      "Ends at:"
    ),
    "\n"
  )
  print(x$time_end)
  # intervention impact
  cat(
    "\n",
    cli::col_magenta(
      "Reduction:"
    ),
    "\n"
  )
  print(effect)

  invisible(x)
}

#' Convert a list to a intervention object
#'
#' @param x A list, or an object that inherits from a list.
#' @param type A string for the type of intervention: `"contacts"` for a
#' `<contact_intervention>` or `"rate"` for a `<rate_intervention>`.
#' @return A [intervention] class object.
#' @export
#' @examples
#' # prepare a list
#' npi <- list(
#'   name = "npi",
#'   type = "contacts",
#'   time_begin = 30,
#'   time_end = 60,
#'   reduction = rep(0.1, 3)
#' )
#'
#' as.intervention(npi)
as.intervention <- function(x, type = c("contacts", "rate")) {
  # check that input is a list or intervention
  stopifnot(
    "Input must inherit from `list`" =
      is.list(x)
  )
  # check type argument
  type <- match.arg(type)

  x <- intervention(
    name = x[["name"]],
    type = type,
    time_begin = x[["time_begin"]],
    time_end = x[["time_end"]],
    reduction = x[["reduction"]]
  )

  # return x
  x
}

#' Concatenate contact interventions for use in an epidemic model
#'
#' @name intervention
#' @rdname intervention
#'
#' @export
c.contacts_intervention <- function(x, ...) {
  # collect inputs
  multi_npi <- list(x, ...)
  invisible(
    lapply(multi_npi, validate_contacts_intervention)
  )

  # check that all intervention regimes have the same dimensions
  # of intervention rates --- these are identical to dims of start and end times
  stopifnot(
    "All <contact_intervention> `reduction`s must have identical dimensions" =
        vapply(multi_npi, function(npi) {
          identical(nrow(npi$reduction), nrow(x$reduction))
        }, FUN.VALUE = logical(1))
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
  npi_names <- glue::glue("npi_{seq_len(ncol(multi_npi$reduction))}")

  # add names to doses for comprehension when printed
  for (i in c("time_begin", "time_end", "reduction")) {
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
    multi_npi$reduction,
    class = "contacts_intervention"
  )

  # validate new object
  validate_contacts_intervention(multi_npi)

  # return object
  multi_npi
}

#' Concatenate interventions for use in an epidemic model
#'
#' @name intervention
#' @rdname intervention
#'
#' @export
c.rate_intervention <- function(x, ...) {
  # collect inputs
  multi_npi <- list(x, ...)
  invisible(
    lapply(multi_npi, validate_rate_intervention)
  )

  # no checking on reduction dimensions

  # strip class and `name` member from `multi_npi`
  multi_npi <- lapply(multi_npi, function(z) {
    z <- unclass(z)
    z$name <- NULL
    z
  })

  # modify x to return a list object of multiple start and end times and nu-s
  multi_npi <- do.call(
    Map, c(f = c, multi_npi)
  )

  # add name parameter --- take "name" of `x`
  multi_npi$name <- x$name

  # generate intervention dose names as dose_1 ... dose_n
  # get total number of doses from the sum of all columns
  interv_names <- glue::glue("interv_{seq_along(multi_npi$reduction)}")

  # add names to doses for comprehension when printed
  for (i in c("time_begin", "time_end", "reduction")) {
    names(multi_npi[[i]]) <- interv_names
  }

  # convert resulting object to intervention
  multi_npi <- new_intervention(
    multi_npi$name,
    multi_npi$time_begin,
    multi_npi$time_end,
    multi_npi$reduction,
    class = "rate_intervention"
  )

  # validate new object
  validate_rate_intervention(multi_npi)

  # return object
  multi_npi
}

#' Calculate the Cumulative Effect of Interventions on Social Contacts
#'
#' @name cumulative_contacts_intervention
#' @rdname cumulative_contacts_intervention
#'
#' @param t The current time.
#' @param time_begin A numeric vector of the start times of all interventions
#' being modelled.
#' @param time_end A numeric vector of the end times of all interventions being
#' modelled. Must be the same length as `time_begin`.
#' @param reduction A numeric matrix where rows give the effect of interventions
#' on each demographic group, and the columns give the proportion reduction in
#' contacts. When two interventions overlap, the proportions are _added_, for a
#' maximum possible value of 1.0 (i.e., no contacts).
#'
#' @keywords internal
#'
#' @return
#' `cumulative_contacts_intervention()` returns a numeric vector of the
#' proportion reduction in contacts for each demographic group.
#'
#' `intervention_on_cm()` returns the contact matrix `cm` scaled by the
#' cumulative effect of any active interventions.
#'
cumulative_contacts_intervention <- function(t, # nolint: object_length_linter
                                             time_begin, time_end, reduction) {
  # determine which interventions are active, promote to numeric
  interventions_active <- as.vector(t > time_begin & t < time_end)

  # check for null interventions where all values are zero
  # first, define value when no interventions are active
  cumulative_cr <- numeric(nrow(reduction))
  if (any(interventions_active)) {
    # multiply the columns of the contact reduction matrix by the
    # active interventions logical vector --- matrix multiplication is faster
    # than two uses of t()
    # then sum the rows for a vector of the cumulative effect
    cumulative_cr <- rowSums(reduction %*% diag(
      as.numeric(interventions_active)
    ))
  }

  # fix values > 1.0
  cumulative_cr[cumulative_cr > 1.0] <- 1.0

  # return cumulative contact reduction
  cumulative_cr
}

#' Scale a Contact Matrix by all Active Interventions on Social Contacts
#'
#' @name cumulative_contacts_intervention
#' @rdname cumulative_contacts_intervention
#'
#' @param cm A numeric matrix of social contacts between demographic groups.
#'
#' @keywords internal
intervention_on_cm <- function(t, cm, time_begin, time_end, cr) {
  # cumulative cr
  contact_scaling <- 1.0 -
    cumulative_contacts_intervention(t, time_begin, time_end, cr)
  # first multiply rows by reduction
  cm_mod <- cm * contact_scaling
  # then multiply cols by reduction and return
  cm_mod %*% diag(contact_scaling)
}

#' Calculate the Cumulative Effect of Interventions on Rate Parameters
#'
#' @name cumulative_rate_intervention
#' @rdname cumulative_rate_intervention
#'
#' @param t The current time.
#' @param time_begin A numeric vector of the start times of all interventions
#' being modelled.
#' @param time_end A numeric vector of the end times of all interventions being
#' modelled. Must be the same length as `time_begin`.
#' @param reduction A numeric vector where each element gives the effect of the
#' corresponding intervention on model rate parameters. When two interventions
#' overlap, the proportions are _added_, for a maximum possible value of 1.0
#' (i.e., rate set to zero).
#'
#' @keywords internal
#'
#' @return
#' `.cumulative_rate_intervention()` returns a number of the proportion
#' reduction in a model rate parameter.
#'
#' `intervention_on_cm()` returns the contact matrix `cm` scaled by the
#' cumulative effect of any active interventions.
#'
.cumulative_rate_intervention <- function(t, time_begin, time_end, reduction) {
  # determine which interventions are active, promote to numeric
  interventions_active <- as.vector(t > time_begin & t < time_end)

  # sum active effects --- if none are active the sum is automatically 0
  cumulative_effect <- sum(reduction[interventions_active])

  # fix values > 1.0
  cumulative_effect[cumulative_effect > 1.0] <- 1.0

  # return cumulative contact reduction
  cumulative_effect
}

#' Apply interventions to rate parameters
#'
#' @param t A single number for the simulation time.
#' @param interventions A named list of `list`-like objects that each have at
#' least the three members `"time_begin"`, `"time_end"`, and `"reduction"`.
#' These are used to calculate the effect on each of the named parameters in the
#' simulation.
#' @param parameters A named list of numeric parameters affected by
#' `interventions`.
#' This represents the model parameters, such as the transmission rate,
#' \eqn{\beta}, or the recovery rate, \eqn{\gamma}.
#' @return A named list of the same length as `parameters`, with the same names.
#' These parameters can then be used in a timestep of an ODE model.
.intervention_on_rates <- function(t, interventions, parameters) {
  new_values <- Map(
    interventions, names(interventions),
    f = function(interv, name) {
      effect <- .cumulative_rate_intervention(
        t = t,
        time_begin = interv[["time_begin"]], time_end = interv[["time_end"]],
        reduction = interv[["reduction"]]
      )
      parameters[[name]] * (1 - effect)
    }
  )

  parameters[names(new_values)] <- new_values

  # return parameters
  parameters
}

#' Generate a null intervention on rates
#'
#' @return A dummy `<rate_intervention>`.
#' @noRd
#' @keywords internal
.no_rate_intervention <- function() {
  suppressMessages(
    intervention(
      type = "rate",
      time_begin = 0, time_end = 0, reduction = 0
    )
  )
}
