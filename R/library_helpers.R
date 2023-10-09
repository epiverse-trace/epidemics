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

#' Get names of available models
#'
#' @param model_type The model type. The only option currently available is
#' "epidemic". Future development may see the addition of more model types,
#' including endemic disease compartmental models, and branching process
#' outbreak models.
#'
#' @export
get_model_names <- function(model_type = "epidemic") {
  checkmate::assert_string(model_type)
  model_type <- match.arg(arg = model_type, several.ok = FALSE)
  # read in the model library using the predefined helper
  model_lib <- model_library()
  model_lib[model_lib$model_type == model_type, ]$model_name
}