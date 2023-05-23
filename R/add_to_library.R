#' Read and display the model library
#'
#' @export
#' @return `A data.frame` representing the model library.
model_library <- function() {
  jsonlite::read_json(
    path = system.file(
      "extdata", "model_library.json",
      package = "epidemics",
      mustWork = TRUE
    ),
    simplifyVector = TRUE
  )
}

#' Add model functions to the model library
#'
#' @param model_type The model type, as a single word, lower case string.
#' The only currently supported type is "epidemic", with future additions likely
#' to include "outbreak".
#' @param model_name The model name, as a single word, lower case string.
#' Please name your model informatively yet concisely - such as a pathogen name,
#' author or institution name, or the name of a previous model implementation.
#' Overall, the name should help users identify your model quickly and reliably.
#' @param compartments The model compartments' names, as a single word, lower
#' case string. Please use names that correspond to prevailing epidemiological
#' naming conventions, and please check the model library for the names used
#' in other models and consider re-using these names when naming the
#' compartments in your own model.
#'
#' @return No return type; updates the model library JSON file with the model
#' type and name provided by the user.
#' @export
add_to_library <- function(model_type = "epidemic",
                           model_name = "default",
                           compartments) {
  # input checking
  checkmate::assert_string(model_type)
  checkmate::assert_string(model_name)

  # read in model library
  model_library <- model_library()

  # check whether there is already a model of this name
  stopifnot(
    "Model library already has a model of this name - choose another name" =
      !(model_name %in% model_library$model_name)
  )

  # create a data.frame of model name, type, and allied functions
  temp_library <- data.frame(
    model_type = model_type,
    model_name = model_name,
    model_function = sprintf(".%s_%s_cpp", model_type, model_name),
    model_args_checker = sprintf(
      ".check_args_%s_%s", model_type, model_name
    ),
    model_args_prepper = sprintf(
      ".prepare_args_%s_%s", model_type, model_name
    ),
    stringsAsFactors = FALSE
  )
  temp_library$compartments <- list(compartments)

  # read in the model library as a JSON file
  model_library <- jsonlite::read_json(
    path = system.file(
      "extdata", "model_library.json",
      package = "epidemics",
      mustWork = TRUE
    ),
    simplifyVector = TRUE
  )

  # bind library and data.frame
  model_library <- rbind(
    model_library,
    temp_library
  )

  # sort by type and then name
  model_library <- model_library[
    order(model_library$model_type, model_library$model_name),
  ]

  # write model library to file with a message
  jsonlite::write_json(
    model_library,
    system.file("extdata", "model_library.json", package = "epidemics"),
    pretty = TRUE
  )
  message(
    sprintf(
      "Adding '%s' type model named '%s' to the model library.\n",
      model_type, model_name
    ),
    "Please check that this is correct before commiting your changes."
  )
}
