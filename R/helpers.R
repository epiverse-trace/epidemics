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
  )
  data$time <- l[["time"]]

  # return data
  data
}
