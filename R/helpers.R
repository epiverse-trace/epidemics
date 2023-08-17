
#' Get parameters from \{epidemics\} classes
#'
#' @param x An object of one of the classes provided by _epidemics_:
#' `<population>`, `<infection>`, `<intervention>`, or `<vaccination>`.
#' @param parameter A string for the parameter to access from `x`.
#'
#' @return An object of the class of `parameter` from `x`.
#' @export
#'
#' @examples
#' # prepare an infection and access the r0 member
#' pandemic <- infection(
#'   r0 = 1.5,
#'   preinfectious_period = 3,
#'   infectious_period = 7
#' )
#'
#' get_parameter(pandemic, "r0")
get_parameter <- function(x, parameter) {
  stopifnot(
    "`x` must be <population>, <infection>, <intervention>, or <vaccination>" =
      is_population(x) || is_infection(x) || is_intervention(x) ||
        is_vaccination(x)
  )
  checkmate::assert_string(parameter)

  # return list element
  x[[parameter]]
}

#' Return ODE model output as a data.table
#' @param output The model output, which must be a two element list (for
#' epidemic) models, with the names "x" and "time", where "x" represents the
#' condition of each compartment at each timestep in "time".
#' @param population A `<population>` object corresponding to the population
#' used in the epidemic model; see [population()].
#' The `<population>` object is used to generate the names of the demographic
#' groups, if these are named.
#' @param compartments A vector for the model compartment names.
#' @keywords internal
#' @return A `data.table` with the columns "compartment", "demography_group",
#' "value", and "time"; these specify the epidemiological compartment, the
#' name of the demography group, the number of individuals of that group in the
#' compartment, and the model timestep, respectively.
#' Names for the demographic groups are generated if no names are provided in
#' the `population` object; these are of the form "demo_group_X".
output_to_df <- function(output, population, compartments) {
  # input checking
  checkmate::assert_list(output,
    any.missing = FALSE, all.missing = FALSE,
    len = 2L, names = "unique"
  )
  checkmate::assert_names(
    names(output),
    must.include = c("x", "time")
  )
  checkmate::assert_character(compartments,
    any.missing = FALSE,
    all.missing = FALSE, unique = TRUE
  )
  checkmate::assert_class(population, "population")

  # get demographic group names if any
  names_demo_groups <- rownames(population$contact_matrix)
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

  # return a data.table
  data.table::data.table(
    time = rep(output$time, each = n_groups * length(compartments)),
    demography_group = vec_demo_groups,
    compartment = vec_compartments,
    value = as.vector(values)
  )
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
#' @param data A `data.table` (or `data.frame`) of model output, typically
#' the output of [epidemic_default_cpp()] or similar functions.
#' @param stage The stage of the epidemic at which to return the epidemic size;
#' here, 0.0 represents the initial conditions of the epidemic (0% of model time
#' ), while 1.0 represents the end of the epidemic model (100% of model time).
#' The values returned at `stage = 1.0` represent the _final size_ of the
#' epidemic.
#' @param by_group A logical representing whether the epidemic size should be
#' returned by demographic group, or whether a single population-wide value is
#' returned.
#' @param deaths A logical value that indicates whether to count individuals in
#' the epidemic size calculation. Setting `deaths = TRUE` looks for a `"dead"`
#' compartment in the data. If there is no such column, the function returns
#' only the final number of recovered individuals in each demographic group.
#'
#' @return A single number when `by_group = FALSE`, or a vector of numbers of
#' the same length as the number of demographic groups when `by_group = TRUE`.
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
#' # create an infection object
#' pandemic_influenza <- infection(
#'   r0 = 1.5, infectious_period = 7, preinfectious_period = 3
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- epidemic_default_cpp(
#'   population = uk_population,
#'   infection = pandemic_influenza,
#'   time_end = 200,
#'   increment = 1
#' )
#'
#' # get the final epidemic size
#' epidemic_size(data)
#'
#' # get the epidemic size at the halfway point
#' epidemic_size(data, stage = 0.5)
epidemic_size <- function(data, stage = 1.0, by_group = TRUE, deaths = TRUE) {
  # input checking for data
  checkmate::assert_data_table(data)
  checkmate::assert_logical(by_group)
  checkmate::assert_logical(deaths)
  checkmate::assert_number(stage, lower = 0.0, upper = 1.0, finite = TRUE)

  stopifnot(
    "No 'recovered' compartment found in `data`, check model compartments" =
      "recovered" %in% unique(data$compartment)
  )
  # if deaths are requested to be counted, but no "dead" compartment exists
  # throw a message
  if (deaths && (!"dead" %in% unique(data$compartment))) {
    message(
      "No 'dead' compartment found in `data`; counting only 'recovered'",
      " individuals in the epidemic size."
    )
  }
  # add deaths to compartments to search
  size_compartments <- "recovered"
  if (deaths) {
    size_compartments <- c(size_compartments, "deaths")
  }

  # get final numbers recovered - operate on data.table as though data.frame
  epidemic_size_ <- data[data$compartment %in% size_compartments &
    data$time == round(max(data$time) * stage, 2), ]

  if (by_group) {
    epidemic_size_ <- epidemic_size_[["value"]]
    names(epidemic_size_) <- unique(data$demography_group)
  } else {
    epidemic_size_ <- sum(epidemic_size_[["value"]])
    names(epidemic_size_) <- "total_population"
  }

  # return epidemic size
  epidemic_size_
}

#' Get new infections
#'
#' @param data A `data.table` (or `data.frame`) of model output, typically
#' the output of [epidemic_default_cpp()] or similar functions.
#' @param compartments_from_susceptible An optional argument, for a character
#' vector of the names of model compartments into which individuals transition
#' from the "susceptible" compartment, and which are not related to infection.
#' A common example is a compartment for "vaccinated" individuals who are no
#' longer susceptible, but who should also not be counted as infected.
#' @param by_group A logical representing whether the epidemic size should be
#' returned by demographic group, or whether a single population-wide value is
#' returned.
#' @return A `data.table` with the same columns as `data`, but with the
#' additional variable under `compartment`, "new_infections", resulting in
#' additional rows.
#' @export
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
#' # create an infection
#' pandemic_influenza <- infection(
#'   r0 = 1.5, infectious_period = 7, preinfectious_period = 3
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- epidemic_default_cpp(
#'   population = uk_population,
#'   infection = pandemic_influenza,
#'   time_end = 200,
#'   increment = 1
#' )
#'
#' new_infections(data)
#'
new_infections <- function(data,
                           compartments_from_susceptible,
                           by_group = TRUE) {
  # input checking for class and susceptible compartment
  checkmate::expect_data_table(data)
  checkmate::assert_logical(by_group, len = 1L)
  stopifnot(
    "Compartment 'susceptible' not found in data, check compartment names." =
      "susceptible" %in% unique(data$compartment)
  )

  # check for compartments deriving from susceptible
  if (!missing(compartments_from_susceptible)) {
    checkmate::assert_character(
      compartments_from_susceptible,
      any.missing = FALSE
    )
    stopifnot(
      "Compartments from 'susceptible' not all found in data, check names." =
        all(compartments_from_susceptible %in% unique(data$compartment))
    )
  }

  # cast data wide, this makes a copy
  data <- data.table::dcast(
    data,
    time + demography_group ~ compartment,
    value.var = "value"
  )

  # calculate new infections as the change in susceptibles -
  # the change in susceptibles due to non-infection related transitions
  # such as vaccination
  if (missing(compartments_from_susceptible)) {
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
  # or aggregated otherwise
  # do not return other compartments
  data <- data[, c("time", "demography_group", "new_infections")]
  if (!by_group) {
    data <- data[, list(new_infections = sum(new_infections)), by = "time"]
  }

  # return data
  data
}
