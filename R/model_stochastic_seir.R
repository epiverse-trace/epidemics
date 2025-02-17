#' @title Model an stochastic SEIR epidemic with interventions
#'
#' @name model_stochastic_seir
#' @rdname model_stochastic_seir
#'
#' @description Simulate an epidemic using a stochastic, compartmental
#' epidemic model with the compartments
#' "susceptible", "exposed", "infectious", and "recovered"
#' This model can accommodate heterogeneity in social contacts among demographic
#' groups, as well as differences in the sizes of demographic groups.
#'
#' The `population`, `transmission_rate`, `infectiousness_rate`, and
#' `recovery_rate`
#' arguments are mandatory, while passing an `intervention`
#' is optional and can be used to simulate scenarios with different epidemic
#' responses or different levels of the same type of response.
#' See **Details** for more information.
#'
#' @param population An object of the `population` class, which holds a
#' population contact matrix, a demography vector, and the initial conditions
#' of each demographic group. See [population()].
#' @param transmission_rate A numeric for the rate at which individuals
#' move from the susceptible to the exposed compartment upon contact with an
#' infectious individual. Often denoted as \eqn{\beta}, with
#' \eqn{\beta = R_0 / \text{infectious period}}. See **Details** for default
#' values.
#' @param infectiousness_rate A numeric for the rate at which individuals
#' move from the exposed to the infectious compartment. Often denoted as
#' \eqn{\sigma}, with \eqn{\sigma = 1.0 / \text{pre-infectious period}}.
#' This value does not depend upon the number of infectious individuals in the
#' population. See **Details** for default values.
#' @param recovery_rate A numeric for the rate at which individuals move
#' from the infectious to the recovered compartment. Often denoted as
#' \eqn{\gamma}, with \eqn{\gamma = 1.0 / \text{infectious period}}.
#' See **Details** for default values.
#' @param intervention A named list of `<intervention>`s representing optional
#' non-pharmaceutical or pharmaceutical interventions applied during the
#' epidemic. Only a single intervention on social contacts of the class
#' `<contacts_intervention>` is allowed as the named element "contacts".
#' Multiple `<rate_interventions>` on the model parameters are allowed; see
#' **Details** for the model parameters for which interventions are supported.
#' @param time_dependence A named list where each name
#' is a model parameter, and each element is a function with
#' the first two arguments being the current simulation `time`, and `x`, a value
#' that is dependent on `time` (`x` represents a model parameter).
#' See **Details** for more information, as well as the vignette on time-
#' dependence \code{vignette("time_dependence", package = "epidemics")}.
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 100 days. May be a numeric vector.
#' @details
#'
#' # Details: Stochastic SEIR model suitable for directly transmitted infections
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates of
#' transitions between compartments per model run.
#'
#' However, model parameters may be passed as numeric vectors. These vectors
#' must follow Tidyverse recycling rules: all vectors must have the same length,
#' or, vectors of length 1 will be recycled to the length of any other vector.
#'
#' The default values are:
#'
#' - Transmission rate (\eqn{\beta}, `transmission_rate`): 0.186, assuming an
#' \eqn{R_0} = 1.3 and an infectious period of 7 days.
#'
#' - Infectiousness rate (\eqn{\sigma}, `infectiousness_rate`): 0.5, assuming
#' a pre-infectious period of 2 days.
#'
#' - Recovery rate (\eqn{\gamma}, `recovery_rate`): 0.143, assuming an
#' infectious period of 7 days.
#'
#' @return A `<data.table>`.
#' If the model parameters and composable elements are all scalars, a single
#' `<data.table>` with the columns "time", "compartment", "age_group", and
#' "value", giving the number of individuals per demographic group
#' in each compartment at each timestep in long (or "tidy") format is returned.
#'
#' If the model parameters or composable elements are lists or list-like,
#' a nested `<data.table>` is returned with a list column "data", which holds
#' the compartmental values described above.
#' Other columns hold parameters and composable elements relating to the model
#' run. Columns "scenario" and "param_set" identify combinations of composable
#' elements (population, interventions), and infection
#' parameters, respectively.
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
#' # and three discrete values of transmission rate
#' data <- model_stochastic_seir(
#'   population = uk_population,
#'   transmission_rate = c(1.3, 1.4, 1.5) / 7.0, # uncertainty in R0
#' )
#'
#' # view some data
#' data
#'
#' # run epidemic simulations with differences in the end time
#' # may be useful when considering different start dates with a fixed end point
#' data <- model_stochastic_seir(
#'   population = uk_population,
#'   time_end = c(50, 100, 150)
#' )
#'
#' data
#' @export
model_stochastic_seir <- function(population,
                          transmission_rate = 1.3 / 7.0,
                          infectiousness_rate = 1.0 / 2.0,
                          recovery_rate = 1.0 / 7.0,
                          intervention = NULL,
                          time_dependence = NULL,
                          time_end = 100 ,
                          n_samples = 1000 ) {
  # get compartment names
  compartments <- c(
    "susceptible", "exposed", "infectious", "recovered"
  )
  assert_population(population, compartments)

  # NOTE: model rates very likely bounded 0 - 1 but no upper limit set for now
  checkmate::assert_numeric(transmission_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(infectiousness_rate, lower = 0, finite = TRUE)
  checkmate::assert_numeric(recovery_rate, lower = 0, finite = TRUE)
  checkmate::assert_integerish(time_end, lower = 0)

  # check the time end 
  # restrict increment to lower limit of 1e-6
  checkmate::assert_integerish(time_end, lower = 0)

  # check all vector lengths are equal or 1L
  params <- list(
    transmission_rate = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate = recovery_rate,
    time_end = time_end
  )
  # take parameter names here as names(DT) updates by reference!
  param_names <- names(params)

  # Check if `intervention` is a single intervention set or a list of such sets
  # NULL is allowed;
  is_lofints <- checkmate::test_list(
    intervention, "intervention",
    all.missing = FALSE, null.ok = TRUE
  )
  # allow some NULLs (a valid no intervention scenario) but not all NULLs
  is_lofls <- checkmate::test_list(
    intervention,
    types = c("list", "null"), all.missing = FALSE
  ) &&
    # Check that all elements of intervention sets are either `<intervention>`
    # or NULL
    all(
      vapply(
        unlist(intervention, recursive = FALSE),
        FUN = function(x) {
          is_intervention(x) || is.null(x)
        }, TRUE
      )
    )

  # Check if parameters can be recycled;
  stopifnot(
    "All parameters must be of the same length, or must have length 1" =
      .test_recyclable(params),
    "`intervention` must be a list of <intervention>s or a list of such lists" =
      is_lofints || is_lofls
  )

  # make lists if not lists
  if (is_lofints) {
    intervention <- list(intervention)
  }

  # check that time-dependence functions are passed as a list with at least the
  # arguments `time` and `x`, in order as the first two args
  # NOTE: this functionality is not vectorised;
  # convert to list for data.table list column
  checkmate::assert_list(
    time_dependence, "function",
    null.ok = TRUE,
    any.missing = FALSE, names = "unique"
  )
  # lapply on null returns an empty list
  invisible(
    lapply(time_dependence, checkmate::assert_function,
      args = c("time", "x"), ordered = TRUE
    )
  )
  time_dependence <- list(
    .cross_check_timedep(
      time_dependence,
      c("transmission_rate", "infectiousness_rate", "recovery_rate")
    )
  )
  
  # set up matrices for the state
  n_groups <- length( population$demography_vector )
  n_cells  <- n_groups * n_samples
  pop      <- population$demography_vector

  # TO MIMIC BEHAVIOUR OF model_default(), NEED TO DEAL WITH THE INITIAL CONDITIONS
  # BEING UNNAMED OR NAMED S,E,I,R (I.E. NOT THE COMPARTMENT NAMES)
  initial_conds <- population$initial_conditions
  name_overlap  <- intersect( compartments, colnames( population$initial_conditions ) )
  if( length( name_overlap ) == 0 ) {
    initial_conds <- `colnames<-`( initial_conds, compartments )
  } else checkmate::assert( length( name_overlap ) == ncol( initial_conds), name = "initial_conditions column names")
  
  # put in initial conditions
  states  <- list()
  for( state in compartments ) {
    states[[ state ]] <- matrix( ceiling( initial_conds[ , state ] * pop ), 
                                 nrow = n_groups, ncol = n_samples )
  }
  
  # set up rates
  rates <- list(
    SE = transmission_rate,
    EI = infectiousness_rate,
    IR = recovery_rate
  )
  
  # adjust the contact to includes rates and population information 
  contact_matrix <- population$contact_matrix
  contact_matrix <- diag( 1 / mean( contact_matrix) / n_groups / pop ) %*% contact_matrix
  
  # add time dependence to rates and contat matrix
  time_dep_rate <- .process_time_dependent_rates( contact_matrix, rates, time_end, interventions )
  rates <- time_dep_rate$rate
  contact_matrix <- time_dep_rate$contact_matrix
    
  # set up flows
  flows <- list()
  flow_names <- names( rates )
  for( flow in flow_names ) {
    flows[[ flow ]] <- matrix( NA, nrow = n_groups, ncol = n_samples )
  }
  
  # set up output 
  outputs <- vector( mode = "list", length = time_end + 1 )
  output_template <- data.table(
    sample           = rep( 1:n_samples, each = n_groups, length( compartments ) ),
    demography_group = rep( names( pop ), n_samples * length( compartments ) ),
    compartment      = rep( compartments, each = n_groups * n_samples )
  )
  outputs[[ 1 ]] <- copy( output_template )
  outputs[[ 1 ]][ , time := 0 ]
  outputs[[ 1 ]][ , value := unlist( lapply( states, as.vector ) )]
  
  for( tdx in 1:time_end ) {
    # calculate the transisiotns between compartments
    infectious_contacts <- ( ( contact_matrix[[ tdx ]] * rates[[ tdx ]]$SE ) %*% states[[ "infectious" ]] ) 
    flows[[ "SE" ]] <- matrix( rbinom( n_cells, states[[ "susceptible" ]], infectious_contacts ), nrow = n_groups )
    flows[[ "EI" ]] <- matrix( rbinom( n_cells, states[[ "exposed" ]], rates[[ tdx ]]$EI ), nrow = n_groups )
    flows[[ "IR" ]] <- matrix( rbinom( n_cells, states[[ "infectious" ]], rates[[ tdx ]]$IR ), nrow = n_groups )
    
    # update the number in each state
    states[[ "susceptible" ]] <- states[[ "susceptible" ]] - flows[[ "SE" ]]
    states[[ "exposed" ]]     <- states[[ "exposed" ]] + flows[[ "SE" ]] - flows[[ "EI" ]]
    states[[ "infectious" ]]  <- states[[ "infectious" ]] + flows[[ "EI" ]] - flows[[ "IR" ]]
    states[[ "recovered" ]]   <- states[[ "recovered" ]] + flows[[ "IR" ]] 
    
    # record output in a new data table
    outputs[[ tdx+1 ]] <- copy( output_template )
    outputs[[ tdx+1 ]][ , time := tdx ]
    outputs[[ tdx+1 ]][ , value := unlist( lapply( states, as.vector ) )]
  }
  outputs <- rbindlist( outputs )

  return( outputs )
}

#' .process_time_dependent_rates
#'
#' @description Converts interventions in to time dependent changes in parameters
#' and transmission rates
#' @param contact_matrix The base contact matrix
#' @param rates The base transition rates
#' @param time_end The length of the simulation
#' @param interventions The interverntions applied
#'
#' @return A list of the contact matrix and transition rates applicable for 
#' each step of the simulation
#' @keywords internal
.process_time_dependent_rates = function( contact_matrix, rates, time_end, interventions ) {
  # start with simply copying them
  contact_matrix <- lapply( 1:time_end, function( x ) contact_matrix )
  rates          <- lapply( 1:time_end, function( x ) rates )
  
  return( list( contact_matrix = contact_matrix, rates = rates ) )
}

