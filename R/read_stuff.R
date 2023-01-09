#' Example of a function to read stuff
#'
#' This function reads stuff into R. More details...
#'
#' @param x the path to the file to read stuff from
#'
#' @author Thibaut Jombart
#'
#' @export
#'
#' @examples
#' some_path <- "path to your file"
#' read_stuff(some_path)

read_stuff <- function(x) {
  # check inputs
  checkmate::assertCharacter(x, len = 1L)

  
  # do stuff
  x
}
