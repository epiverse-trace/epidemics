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