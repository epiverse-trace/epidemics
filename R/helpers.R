#' Return ODE model output as a data.table
#'
#' @name output_to_df
#' @rdname output_to_df
#' @param output The model output, which must be a two element list (for
#' epidemic) models, with the names "x" and "time", where "x" represents the
#' condition of each compartment at each timestep in "time".
#' @param population A `<population>` object corresponding to the population
#' used in the epidemic model.
#' The `<population>` object is used to generate the names of the demographic
#' groups, if these are named.
#' @param compartments A vector for the model compartment names.
#' @keywords internal
#' @return A `<data.table>` with the columns "compartment", "demography_group",
#' "value", and "time"; these specify the epidemiological compartment, the
#' name of the demography group, the number of individuals of that group in the
#' compartment, and the model timestep, respectively.
#' Names for the demographic groups are generated if no names are provided in
#' the `population` object; these are of the form "demo_group_X".
.output_to_df <- function(output, population, compartments) {
  # get demographic group names if any
  groupnames_from_contacts <- rownames(population$contact_matrix)
  groupnames_from_demography <- names(population$demography_vector)
  names_demo_groups <- Find(
    function(x) !is.null(x),
    list(groupnames_from_contacts, groupnames_from_demography)
  )

  if (is.null(names_demo_groups)) {
    names_demo_groups <- sprintf(
      "demo_group_%i",
      seq_len(nrow(population$contact_matrix))
    )
  }

  # count groups and timesteps to generate vectors of compartment and demo
  # group names
  n_groups <- length(names_demo_groups)
  n_timesteps <- length(output[["time"]])

  vec_compartments <- rep(compartments, each = n_groups)
  vec_compartments <- rep(vec_compartments, times = n_timesteps)

  vec_demo_groups <- rep(names_demo_groups, length(compartments) * n_timesteps)

  # handle values as potential list, otherwise matrix
  # matrix will need to be transposed for conversion to long format
  values <- output$x
  if (is.list(values)) {
    values <- unlist(values)
  } else if (is.matrix(values)) {
    values <- t(values)
  }

  # return a data.table - this will typically be one among a list,
  # and passed to data.table::rbindlist()
  data <- data.table::data.table(
    time = rep(output$time, each = n_groups * length(compartments)),
    demography_group = vec_demo_groups,
    compartment = vec_compartments,
    value = as.vector(values)
  )

  data
}

#' Return Ebola model output as a table
#' @name output_to_df
#' @rdname output_to_df
.output_to_df_ebola <- function(output, population, compartments) {
  # calculate replicates, format output, and assign replicate numbers
  replicates <- seq_along(output)
  output <- lapply(output, .output_to_df, population, compartments)

  output <- Map(output, replicates, f = function(df, i) {
    df[["replicate"]] <- i
    df
  })

  # return single data.table
  data.table::rbindlist(output)
}

#' Get the epidemic size
#'
#' Get the size of the epidemic at any stage between the start and the end.
#' This is calculated as the number of individuals _recovered_ from infection
#' at that stage of the epidemic.
#'
#' This function can be used to calculate the
#' _final size_ of the epidemic, by setting `stage = 1.0` (100% of model time;
#' the default).
#'
#' The function allows for the calculation of epidemic sizes by demographic
#' group as well as the total epidemic size.
#'
#' @param data A table of model output, typically
#' the output of [model_default()] or similar functions.
#' @param stage A numeric vector for the stage of the epidemic at which to
#' return the epidemic size; here 0.0 represents the start time of the epidemic
#' i.e., the initial conditions of the epidemic simulation, while 1.0 represents
#' the end of the epidemic simulation model (100% of model time).
#' Defaults to 1.0, at which stage returned values represent the _final size_ of
#' the epidemic.
#' This value is overridden by any values passed to the `time` argument.
#' @param time Alternative to `stage`, an integer-like vector for the timepoint
#' of the epidemic at which to return the epidemic size.
#' Overrides any values passed to `stage`.
#' @param by_group A logical representing whether the epidemic size should be
#' returned by demographic group, or whether a single population-wide value is
#' returned. Defaults to `TRUE`.
#' @param include_deaths A logical value that indicates whether to count dead
#' individuals in the epidemic size calculation.
#' Defaults to `FALSE`. Setting `include_deaths = TRUE` makes the function look
#' for a `"dead"` compartment in the data. If there is no such column, the
#' function returns only the final number of recovered or removed individuals in
#' each demographic group.
#' @param simplify A logical determining whether the epidemic size data should
#' be simplified to a vector with one element for each demographic group.
#' If the length of `stage` or `time` is $>$ 1, this argument is overridden and
#' the data are returned as a `<data.table>`.
#' @return
#' If `simplify == TRUE` and a single timepoint is requested, returns a vector
#' of epidemic sizes of the same length as the number of demographic groups.
#' If `by_group == FALSE`, sums the epidemic size to return an overall value for
#' the full population.
#'
#' If multiple timepoints are requested, or if multiple replicates are present
#' under a specially named column "replicate" (only from the Ebola model), no
#' simplification to a vector is possible; returns a `<data.table>` of
#' timepoints and epidemic sizes at each timepoint.
#'
#' All options return the absolute sizes and not proportions.
#' @export
#'
#' @examples
#' # create a population
#' uk_population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0, 0),
#'     nrow = 1, ncol = 5L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- model_default(
#'   population = uk_population
#' )
#'
#' # get the final epidemic size if no other arguments are specified
#' epidemic_size(data)
#'
#' # get the epidemic size at the halfway point
#' epidemic_size(data, stage = 0.5)
#'
#' # alternatively, get the epidemic size at `time = 50`
#' epidemic_size(data, time = 50)
epidemic_size <- function(
    data, stage = 1.0, time = NULL, by_group = TRUE,
    include_deaths = FALSE, simplify = TRUE) {
  # input checking for data - this allows data.tables as well
  checkmate::assert_data_frame(
    data,
    min.cols = 4L, min.rows = 1, any.missing = FALSE
  )
  checkmate::assert_names(
    colnames(data),
    must.include = c("time", "demography_group", "compartment", "value")
  )
  checkmate::assert_logical(by_group, len = 1L)
  checkmate::assert_logical(include_deaths, len = 1L)
  checkmate::assert_numeric(
    stage,
    lower = 0.0, upper = 1.0, finite = TRUE,
    null.ok = TRUE, any.missing = FALSE
  )
  checkmate::assert_integerish(
    time,
    lower = 0, upper = max(data[["time"]]), # not suitable for Ebola model
    null.ok = TRUE, any.missing = FALSE
  )

  stopifnot(
    "No 'recovered' or 'removed' compartment in `data`, check compartments" =
      any(c("removed", "recovered") %in% unique(data$compartment)),
    "`data` should have only one of 'recovered' or 'removed' compartments" =
      !all(c("removed", "recovered") %in% unique(data$compartment)),
    "One of `stage` or `time` must be provided; both are NULL!" =
      !all(is.null(c(stage, time)))
  )
  # if deaths are requested to be counted, but no "dead" compartment exists
  # throw a message
  if (include_deaths && (!"dead" %in% unique(data$compartment))) {
    warning(
      "epidemic_size(): No 'dead' compartment found in `data`; counting only",
      " 'recovered' or 'removed' individuals in the epidemic size."
    )
  }
  # add include_deaths to compartments to search
  size_compartments <- ifelse(
    "recovered" %in% unique(data$compartment),
    "recovered", "removed"
  )
  if (include_deaths) {
    size_compartments <- c(size_compartments, "dead")
  }

  # calculate time to get and override stage if provided
  times_to_get <- round(max(data$time) * stage, 2)
  if (!is.null(time)) {
    cli::cli_inform(
      "epidemic_size(): `time` provided will override any `stage` provided"
    )
    times_to_get <- time
  }

  # determine grouping columns to handle ebola model special case
  grouping_cols <- "time"
  if (by_group) {
    grouping_cols <- c(grouping_cols, "demography_group")
  }
  n_replicates <- 1 # set dummy value
  if ("replicate" %in% colnames(data)) {
    grouping_cols <- c(grouping_cols, "replicate")
    n_replicates <- max(data[["replicate"]])
  }

  if ((length(times_to_get) > 1L || n_replicates > 1) && simplify) {
    warning(
      "Returning epidemic size at multiple time points, or for multiple",
      " replicates; cannot simplify output to vector; returning `<data.table>`"
    )
    simplify <- FALSE
  }

  # get final numbers recovered - operate on data.table as though data.table
  epidemic_size_ <- data[data$compartment %in% size_compartments &
    data$time %in% times_to_get, ]

  # set data.table if not already, reove after #211 is merged
  data.table::setDT(epidemic_size_)

  # NOTE: requires data.table
  epidemic_size_ <- epidemic_size_[,
    list(value = sum(.SD)),
    .SDcols = "value",
    by = grouping_cols
  ]
  if (simplify) {
    epidemic_size_ <- epidemic_size_[["value"]]
  }

  # return epidemic size
  epidemic_size_
}

#' Get new infections over model time
#'
#' @param data A table of model output, typically
#' the output of [model_default()] or similar functions.
#' @param compartments_from_susceptible An optional argument, for a character
#' vector of the names of model compartments into which individuals transition
#' from the "susceptible" compartment, and which are not related to infection.
#' A common example is a compartment for "vaccinated" individuals who are no
#' longer susceptible, but who should also not be counted as infected.
#' @param by_group A logical representing whether the epidemic size should be
#' returned by demographic group, or whether a single population-wide value is
#' returned.
#' @return A table with the same columns as `data`, but with the
#' additional variable under `compartment`, "new_infections", resulting in
#' additional rows.
#' @export
#' @examples
#' # create a population
#' uk_population <- population(
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0, 0),
#'     nrow = 1, ncol = 5L
#'   )
#' )
#'
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- model_default(
#'   population = uk_population,
#'   time_end = 200,
#'   increment = 1
#' )
#'
#' new_infections(data)
#'
new_infections <- function(data,
                           compartments_from_susceptible = NULL,
                           by_group = TRUE) {
  # input checking for class and susceptible compartment
  # input checking for data - this allows data.tables as well
  checkmate::assert_data_frame(
    data,
    min.cols = 4, min.rows = 1, any.missing = FALSE
  )
  checkmate::assert_names(
    colnames(data),
    must.include = c("time", "demography_group", "compartment", "value")
  )
  checkmate::assert_character(
    compartments_from_susceptible,
    min.len = 1, null.ok = TRUE,
    any.missing = FALSE, unique = TRUE
  )
  checkmate::assert_logical(by_group, len = 1L)
  stopifnot(
    "Compartment 'susceptible' not found in data, check compartment names." =
      "susceptible" %in% unique(data$compartment),
    "Compartments from 'susceptible' not all found in data, check names." =
      all(compartments_from_susceptible %in% unique(data$compartment)) ||
        is.null(compartments_from_susceptible)
  )

  # set data to a data.table for internal operations
  data.table::setDT(data)
  # cast data wide, this makes a copy
  data <- data.table::dcast(
    data,
    time + demography_group ~ compartment,
    value.var = "value"
  )

  # check for compartments deriving from susceptible and
  # calculate new infections as the change in susceptibles -
  # the change in susceptibles due to non-infection related transitions
  # such as vaccination
  if (is.null(compartments_from_susceptible)) {
    data[, new_infections := c(0, -diff(get("susceptible"))),
      by = "demography_group"
    ]
  } else {
    data[, new_infections := c(0, -diff(get("susceptible"))) -
      Reduce(`+`, lapply(.SD, function(x) {
        c(0, diff(x))
      })),
    .SDcols = compartments_from_susceptible,
    by = "demography_group"
    ]
  }

  # return data in long format, by demographic group by default,
  # or aggregated otherwise; do not return other compartments
  # subsetting below creates a copy - input `data` is safe
  data <- data[, c("time", "demography_group", "new_infections")]
  if (!by_group) {
    data <- data[, list(new_infections = sum(new_infections)), by = "time"]
  }

  # return data
  data
}

#' Get the time and size of a compartment's highest peak
#'
#' Get the time and size of a compartment's highest peak for all
#' demography groups.
#'
#' @details
#' This is used for epidemics with a single peak. It is useful from
#' a public health policy point of view to determine how bad an epidemic will
#' be and when that happens.
#'
#' @param data A `<data.frame>` or `<data.table>` of model output, typically
#' the output of a compartmental model.
#' @param compartments A character vector for the compartments of interest.
#' @return A `<data.table>` with columns "demography_group", "compartment",
#' "time" and "value"; these specify the name of the demography group,
#' the epidemiological compartment, and the peak time and value for each
#' compartment in `compartments`.
#'
#' @export
#'
#' @examples
#' # create a population
#' uk_population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0, 0),
#'     nrow = 1, ncol = 5L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- model_default(
#'   population = uk_population,
#'   time_end = 600
#' )
#'
#' # get the timing and peak of the exposed and infectious compartment
#' epidemic_peak(data, c("exposed", "infectious"))
epidemic_peak <- function(data, compartments = "infectious") {
  # solves "no visible binding for global variable"
  compartment <- NULL
  value <- NULL
  time <- NULL

  # basic input checks
  checkmate::assert_data_frame(
    data,
    min.cols = 4L, min.rows = 1, any.missing = FALSE
  )
  checkmate::assert_names(
    colnames(data),
    must.include = c("time", "demography_group", "compartment", "value")
  )
  checkmate::assert_character(
    compartments,
    min.len = 1,
    any.missing = FALSE, unique = TRUE,
    null.ok = FALSE
  )

  data.table::setDT(data)
  out <- data[compartment %in% compartments,
    list(
      time = data.table::first(time[value == max(value)]),
      value = data.table::first(max(value))
    ),
    by = c("demography_group", "compartment")
  ]

  # return a data.table
  out
}
