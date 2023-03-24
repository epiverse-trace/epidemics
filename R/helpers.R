#' @title Convert epidemic model output to a data.frame
#'
#' @param l A list, usually the output of a call to [.epidemic_default_cpp()].
#' @param compartments A character vector of the compartment names.
#' @return A data.frame with N+1 columns, where N is the number of compartments,
#' with the additional column being "time".
#' @export
output_to_df <- function(l, compartments = c(
                           "susceptible", "exposed",
                           "infectious", "recovered"
                         )) {
  # input checking
  stopifnot(
    "Error: `l` must be a list" =
      is.list(l),
    "Error: `l` must have names 'x' and 'time'" =
      names(l) %in% c("x", "time")
  )
  n_col <- length(l[["x"]][[1]])
  n_groups <- nrow(l[["x"]][[1]])

  # data are in the form S_1->N, E_1->N, etc.
  data <- data.table::as.data.table(
    matrix(
      unlist(l[["x"]]),
      ncol = n_col,
      byrow = TRUE
    )
  )
  colnames(data) <- do.call(
    paste0, expand.grid(compartments, seq(n_groups))
  )
  data$time <- l[["time"]]

  # return data
  data
}
