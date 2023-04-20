#' Return ODE model output as a data.table
#' @param output The model output, which must be a two element list (for
#' epidemic) models, with the names "x" and "time", where "x" represents the
#' condition of each compartment at each timestep in "time".
#' @param model_arguments A list containing the model arguments passed to
#' [epidemic()]. This is scanned for information on the population passed to
#' the model, which must be a `population` object; see [population()].
#' The `population` object is used to generate the names of the demographic
#' groups, if these are named.
#' @param compartments A vector for the model compartment names.
#' @keywords internal
#' @return A `data.table` with the columns "compartment", "demography_group",
#' "value", and "time"; these specify the epidemiological compartment, the
#' name of the demography group, the number of individuals of that group in the
#' compartment, and the model timestep, respectively.
#' Names for the demographic groups are generated if no names are provided in
#' the `population` object; these are of the form "demo_group_X".
output_to_df <- function(output, model_arguments, compartments) {
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

  if ("population" %in% names(model_arguments)) {
    names_demo_groups <- rownames(model_arguments$population$contact_matrix)
    if (is.null(names_demo_groups)) {
      names_demo_groups <- sprintf(
        "demo_group_%i",
        seq_len(nrow(model_arguments$population$contact_matrix))
      )
    }
  }

  # count groups and timesteps to generate vectors of compartment and demo
  # group names
  n_groups <- length(names_demo_groups)
  n_timesteps <- length(output[["time"]])

  vec_compartments <- rep(compartments, each = n_groups)
  vec_compartments <- rep(vec_compartments, times = n_timesteps)

  vec_demo_groups <- rep(names_demo_groups, length(compartments) * n_timesteps)

  # return a data.table
  data.table::data.table(
    compartment = vec_compartments,
    demography_group = vec_demo_groups,
    value = unlist(output$x),
    time = rep(output$time, each = n_groups * length(compartments))
  )
}

#' Get the epidemic size
#'
#' Gets the size of the epidemic at any stage between the start and the end.
#' This is calculated as the number of individuals _recovered_ from infection
#' at that stage of the epidemic. This function can be used to calculate the
#' _final size_ of the epidemic, by setting `stage = 1.0` (the default).
#' The function allows for the calculation of epidemic sizes by demographic
#' group as well as the total epidemic size.
#'
#' @param data A `data.table` (or `data.frame`) of model output, typically
#' the output of [epidemic()].
#' @param stage The stage of the epidemic at which to return the epidemic size;
#' here, 0.0 represents the initial conditions of the epidemic, while 1.0
#' represents the end of the epidemic and the values returned represent the
#' _final size_ of the epidemic.
#' @param by_group A logical representing whether the epidemic size should be
#' returned by demographic group, or whether a single population-wide value is
#' returned.
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
#' # run epidemic simulation with no vaccination or intervention
#' data <- epidemic(
#'   model_name = "default",
#'   population = uk_population,
#'   r0 = 1.5,
#'   preinfectious_period = 3,
#'   infectious_period = 7,
#'   time_end = 200,
#'   increment = 1
#' )
#'
#' # get the final epidemic size
#' epidemic_size(data)
#'
#' # get the epidemic size at the halfway point
#' epidemic_size(data, stage = 0.5)
epidemic_size <- function(data, stage = 1.0, by_group = TRUE) {
  # input checking for data
  checkmate::assert_data_table(data)
  checkmate::assert_logical(by_group)
  checkmate::assert_number(stage, lower = 0.0, upper = 1.0, finite = TRUE)

  stopifnot(
    "No 'recovered' compartment found in `data`, check model compartments" =
      "recovered" %in% unique(data$compartment)
  )

  # get final numbers recovered - operate on data.table as though data.frame
  final_recovered <- data[data$compartment == "recovered" &
    data$time == round(max(data$time) * stage, 2), ]

  if (by_group) {
    final_recovered[["value"]]
  } else {
    sum(final_recovered[["value"]])
  }
}

#' Get new infections
#'
#' @param data
#' @param compartments_from_susceptible
#' @param by_group
#'
#' @importFrom data.table :=
#' @importFrom data.table .SD
#'
#' @return
#' @export
new_infections <- function(data,
                           compartments_from_susceptible = "vaccinated",
                           by_group = TRUE) {
  # input checking
  # check that the column "susceptible" is found in data
  # check that the columns <compartments from susceptible> are found

  # cast data wide
  data_ <- data.table::dcast(
    data.table::copy(data),
    time + demography_group ~ compartment,
    value.var = "value"
  )

  data_[, new_infections := c(0, -diff(susceptible)) -
    Reduce(`+`, lapply(.SD, function(x) {
      c(0, diff(x))
    })),
  .SDcols = compartments_from_susceptible,
  by = "demography_group"
  ]

  # return data
  data_[]
}
