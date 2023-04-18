#' Get functions from model library
#'
#' @param model_type The model type. The only option currently available is
#' "epidemic". Future development may see the addition of more model types,
#' including endemic disease compartmental models, and branching process
#' outbreak models.
#' @param model_name The model name. The only option currently available is
#' "default", for the default epidemic compartmental SEIR-V model.
#' @param which_function A string for the function requested, which may be one
#' of "model_function", "model_args_checker", and "model_args_prepper". These
#' correspond to functions that run the model, check its arguments, and prepare
#' the arguments for the model function, respectively.
#'
#' @keywords internal
read_from_library <- function(model_type = "epidemic",
                              model_name = "default",
                              which_function = "model_function") {

  # input checking
  checkmate::assert_string(model_type)
  checkmate::assert_string(model_name)
  checkmate::assert_string(which_function)

  # read model library
  model_library <- model_library()

  # check whether there are models of the specified type
  stopifnot(
    "Model library has no models of the specified type" =
      model_type %in% unique(model_library$model_type),
    "Model library has no functions of the specified type" =
      which_function %in% grep(
        "function|args", colnames(model_library),
        value = TRUE
      )
  )

  fn <- model_library[
    model_library$model_type == model_type &
      model_library$model_name == model_name,
  ][[which_function]]

  if (length(fn) == 0L) {
    stop(
      sprintf(
        "No model named '%s' of the type '%s' found; ",
        model_name, model_type
      ),
      "check type-name combination."
    )
  } else {
    fn
  }
}
