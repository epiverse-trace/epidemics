#' Create a population object combining several populations
#'
#' @param populations A list of `<population>` objects
#' @param connectivity_matrix A numeric matrix for the contact matrix between
#' the elements of the populations
#' @param method The method to combine the contact matrices of `populations`,
#' can be a character chain (`linear` or `gravity`) which will call an internal
#' function, or a user-defined function with two arguments (`populations` and
#' `connectivity_matrix`) and returning a numeric matrix.
#' @param name Optional string for the combined population name.
#'
#' @return An object of the `<population>` class.
#' @export
#' @examples
#' pop1 <- population(
#'   name = "Population 1",
#'   contact_matrix = matrix(c(1, .5, .5, 1), nrow = 2),
#'   demography_vector = c("0-20" = 2e7, "20+" = 4e7),
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0,
#'       0.9999, 0.0001, 0, 0),
#'     nrow = 2, ncol = 4, byrow = TRUE
#'   )
#' )
#'
#' pop2 <- population(
#'   name = "Population 2",
#'   contact_matrix = matrix(c(1, .5, .5, 1), nrow = 2),
#'   demography_vector = c("0-20" = 1e7, "20+" = 2e7),
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0,
#'       0.9999, 0.0001, 0, 0),
#'     nrow = 2, ncol = 4, byrow = TRUE
#'   )
#' )
#'
#' prop_matrix <- matrix(c(1, .1, .05, 1), nrow = 2)
#'
#' combined_population <- combine_populations(
#'       populations = list(pop1, pop2),
#'       connectivity_matrix = prop_matrix,
#'       method = "linear",
#'       name = "combine"
#' )
#'
#'
#' # check for class <population>
#' is_population(combined_population)
combine_populations <- function(populations, connectivity_matrix,
                                method = "linear",
                                name = NA_character_) {
  # check input
  checkmate::assert_string(name, na.ok = TRUE)
  checkmate::assert_matrix(
    connectivity_matrix, mode = "numeric",
    ncols = length(populations), nrows = length(populations))
  if (is.character(method)) {
    checkmate::assert_choice(method, choices = c("linear", "gravity"))
  } else {
    checkmate::assert_function(
      method, args = c("populations", "connectivity_matrix"))
  }
  # Check that all populations contain the same demographic groups
  names_groups <- rownames(populations[[1]]$contact_matrix)
  if (!is.null(names_groups)) {
    for (i in seq_along(populations)) {
      checkmate::assert_names(rownames(populations[[i]]$contact_matrix),
                              identical.to = names_groups)
    }
  } else {
    for (i in seq_along(populations)) {
      checkmate::assert_vector(populations[[i]]$demography_vector,
                               len = length(populations[[1]]$demography_vector))
    }
  }


  # Create demography vector for the combined population
  combined_demography_vector <- do.call(c, lapply(populations,
                                                  extract_demography))
  # Create initial condition vector for the combined population
  combined_initial_condition <- do.call(rbind, lapply(populations,
                                                      extract_initial))

  # Create contact matrix for the combined population
  if (is.function(method)) {
    combined_contact_matrix <- method(populations, connectivity_matrix)
  } else if (method == "linear") {
    combined_contact_matrix <- combine_contact(populations, connectivity_matrix)
  } else if (method == "gravity") {
    combined_contact_matrix <- gravity_contact(populations, connectivity_matrix)
  }

  # Create combined population
  combined_populations <- population(
    name = name,
    contact_matrix = combined_contact_matrix,
    demography_vector = combined_demography_vector,
    initial_conditions = combined_initial_condition)

  return(combined_populations)
}

#' Extract demography vector from a `<population>` object
#'
#' @param population A `<population>` object
#'
#' @return A named vector corresponding to the demography vector, named after
#' the `name` element of the population object
#' @keywords internal
#' @noRd
extract_demography <- function(population) {
  # Extract name from population
  name_pop <- population$name
  # Extract demography
  demography <- population$demography_vector
  # Rename demography by combining the population name and the name of each
  # element of demography
  names(demography) <- paste(name_pop, names(demography), sep = "_")
  demography
}

#' Extract initial conditions from a `<population>` object
#'
#' @param population A `<population>` object
#'
#' @return A named vector corresponding to the initial conditions, named after
#' the `name` element of the population object
#' @keywords internal
#' @noRd
extract_initial <- function(population) {
  # Extract name from population
  name_pop <- population$name
  # Extract values
  initial <- population$initial_conditions
  # Rename initial by combining the population name and the name of each element
  # of initial
  rownames(initial) <- paste(name_pop,
                             c(rownames(initial),
                               rownames(population$contact_matrix),
                               sprintf(
                                 "demo_group_%i",
                                 seq_len(nrow(initial))
                               ))[seq_len(nrow(initial))],
                             sep = "_")
  initial
}

#' Extract the contact matrix from a `<population>` object
#'
#' @param population A `<population>` object
#'
#' @return A numeric matrix corresponding to the contact vector, named after the
#' `name` element of the population object
#' @keywords internal
#' @noRd
extract_contact <- function(population) {
  # Extract name from population
  name_pop <- population$name
  # Extract the contact matrix
  contact <- population$contact_matrix
  # Rename contact by combining the population name and the name of each element
  # of the contact matrix
  rownames(contact) <- colnames(contact) <-
    paste(name_pop, c(rownames(population$contact_matrix),
                      sprintf(
                        "demo_group_%i",
                        seq_len(nrow(population$contact_matrix))
                      ))[seq_len(nrow(population$contact_matrix))],
          sep = "_")
  contact
}

#' Create the connectivity matrix between all groups of `populations`
#'
#' @param populations A list of `<population>` objects
#' @param connectivity_matrix A numeric matrix for the proportion of contacts
#' between the elements of the populations
#'
#' @return A numeric matrix combining the contact matrices in `<population>` and
#' the connectivity matrix.
#' @keywords internal
combine_contact <- function(populations, connectivity_matrix) {
  # Extract contact matrices from each element of the "populations" list
  contact_matrices <- lapply(populations, extract_contact)

  # Initial combined matrix, an squared matrix of size equal to the sum of the
  # size of each contact matrix
  combined_matrix <-
    matrix(0, nrow = do.call(sum, lapply(contact_matrices, nrow)),
           ncol = do.call(sum, lapply(contact_matrices, ncol)))
  rownames(combined_matrix) <- colnames(combined_matrix) <-
    unlist(lapply(contact_matrices, rownames))

  # compute the number of groups for each population, will be used to
  # compute the elements of the combined matrix
  n_group_per_pop <- unlist(lapply(contact_matrices, nrow))

  for (i in seq_along(populations)) {
    for (j in seq_along(populations)) {
      index_row <-
        1 + cumsum(n_group_per_pop)[i] - rev(seq_len(n_group_per_pop[i]))
      index_col <-
        1 + cumsum(n_group_per_pop)[j] - rev(seq_len(n_group_per_pop[j]))

      # Each element of combined matrix is computed as
      # contacts from x (column of combined_matrix) to y (row) =
      #  nb of contacts from x group to y group in the population x belongs to *
      #  connectivity from the pop x belongs to the pop y belongs to
      combined_matrix[index_row, index_col] <-
        contact_matrices[[j]] * connectivity_matrix[i, j]
    }
  }
  combined_matrix
}

#' Create the connectivity matrix between all groups of `populations` using
#' a gravity model
#'
#' @param populations A list of `<population>` objects
#' @param connectivity_matrix A numeric matrix for the contact matrix between
#' the elements of the populations
#'
#' @return A numeric matrix combining the contact matrices in `<population>` and
#' the connectivity matrix using a gravity model
#' @export
#' @keywords internal
gravity_contact <- function(populations, connectivity_matrix) {
  # Get the number of inhabitants per population by extracting the demography
  # vector and summing the number of inhabitants in all groups
  population_all <- unlist(lapply(lapply(populations, extract_demography),
                                  sum))

  # Extract contact matrix
  contact_matrices <- lapply(populations, extract_contact)

  # Initial combined matrix, an squared matrix of size equal to the sum of the
  # size of each contact matrix
  combined_matrix <-
    matrix(0, nrow = do.call(sum, lapply(contact_matrices, nrow)),
           ncol = do.call(sum, lapply(contact_matrices, ncol)))
  rownames(combined_matrix) <- colnames(combined_matrix) <-
    unlist(lapply(contact_matrices, rownames))

  # compute the number of groups for each population, will be used to
  # compute the elements of the combined matrix
  n_group_per_pop <- unlist(lapply(contact_matrices, nrow))

  # Set connectivity matrix in diagonal to 1 (to avoid Inf)
  diag(connectivity_matrix) <- 1
  connectivity_matrix[connectivity_matrix == 0] <- 1

  # Compute connectivity matrix:
  for (i in seq_along(populations)) {
    for (j in seq_along(populations)) {
      connectivity_matrix[i, j] <-
        connectivity_matrix[i, j]**(-1) * population_all[i] * population_all[j]
    }
  }
  # Set the maximum value of connectivity matrix to 1
  connectivity_matrix <- connectivity_matrix / max(connectivity_matrix)

  # Set diagonal values to 1
  diag(connectivity_matrix) <- 1

  for (i in seq_along(populations)) {
    for (j in seq_along(populations)) {
      index_row <-
        1 + cumsum(n_group_per_pop)[i] - rev(seq_len(n_group_per_pop[i]))
      index_col <-
        1 + cumsum(n_group_per_pop)[j] - rev(seq_len(n_group_per_pop[j]))

      # Each element of combined matrix is computed as
      # contacts from x (column of combined_matrix) to y (row) =
      #  nb of contacts from x group to y group in the population x belongs to *
      #  connectivity from the pop x belongs to the pop y belongs to
      combined_matrix[index_row, index_col] <-
        contact_matrices[[j]] * connectivity_matrix[i, j]
    }
  }
  combined_matrix
}
