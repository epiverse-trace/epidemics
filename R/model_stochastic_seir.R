#' @title Model an stochastic SEIR epidemic with interventions
#'
#' @name model_stochastic_seir
#' @rdname model_stochastic_seir
#'
#' @description Simulate an epidemic using a stochastic compartmental
#' epidemic model with the compartments
#' "susceptible", "exposed", "infectious", and "recovered".
#' The model can accommodate heterogeneity in social contacts among demographic
#' groups, as well as differences in the sizes of demographic groups.
#' Each individual within a compartment is assumed to be independent of the
#' others, with inter-compartment transition times being geometrically distributed.
#' For exposures to infections, the compartments are assumed to be well-mixed,  
#' thus the level of exposure felt be each individual within a compartment is
#' the same. Again, we assume individuals within a compartment are independent 
#' of others in the compartment, thus the number of new infections is binomially
#' distributed. 
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
#' @param time_end The maximum number of timesteps over which to run the model.
#' Taken as days, with a default value of 100 days. 
#' @param n_samples The number of stochastic replicates of the model (default = 1,000).
#' @details
#'
#' # Details: Stochastic SEIR model suitable for directly transmitted infections
#'
#' ## Model parameters
#'
#' This model only allows for single, population-wide rates of
#' transitions between compartments per model run.
#' Additionally, the transmission, infectiousness and recovery rates must
#' be scalars.
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
#' `<data.table>` with the columns "sample", time", "compartment", "age_group", 
#' and"value", giving the number of individuals per demographic group
#' in each compartment at each timestep in long (or "tidy") format is returned.
#'
#
#' @examples
#' # create a population
#' uk_population <- population(
#'   name = "UK population",
#'   contact_matrix = matrix(1),
#'   demography_vector = 67e6,
#'   initial_conditions = matrix(
#'     c(0.9999, 0.0001, 0, 0),
#'     nrow = 1, ncol = 4L
#'   )
#' )
#'
#' # run epidemic simulation with no vaccination or intervention
#' data <- model_stochastic_seir(
#'   population = uk_population,
#'   transmission_rate = 1.5 / 7.0
#' )
#'
#' # view some data
#' data
#' @export
model_stochastic_seir <- function(population,
                          transmission_rate = 1.3 / 7.0,
                          infectiousness_rate = 1.0 / 2.0,
                          recovery_rate = 1.0 / 7.0,
                          intervention = NULL,
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
  checkmate::assert_integerish(n_samples, lower = 1)
  
  # only support scalar parameters
  checkmate::assert_scalar(transmission_rate)
  checkmate::assert_scalar(infectiousness_rate)
  checkmate::assert_scalar(recovery_rate)
  checkmate::assert_scalar(n_samples)
  
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
    transmission_rate   = transmission_rate,
    infectiousness_rate = infectiousness_rate,
    recovery_rate       = recovery_rate
  )
  
  # adjust the contact to includes rates and population information 
  contact_matrix <- population$contact_matrix
  group_names    <- colnames( contact_matrix )
  contact_matrix <- diag( 1 / mean( contact_matrix) / n_groups / pop, nrow = n_groups ) %*% contact_matrix
  
  # add time dependence to rates and contat matrix
  time_dep_rate  <- .prepare_interventions( contact_matrix, rates, time_end, intervention )
  rates          <- time_dep_rate$rates
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
    demography_group = rep( group_names, n_samples * length( compartments ) ),
    compartment      = rep( compartments, each = n_groups * n_samples )
  )
  outputs[[ 1 ]] <- data.table::copy( output_template )
  time  <- NA # HACK TO PREVENT CHECK PACKAGE COMPLAINING ABOUT 
  value <- NA # data.table SYNTAX
  outputs[[ 1 ]][ , time := 0 ]
  outputs[[ 1 ]][ , value := unlist( lapply( states, as.vector ) )]
  
  for( tdx in 1:time_end ) {
    # calculate the transisiotns between compartments
    infectious_contacts <- ( ( contact_matrix[[ tdx ]] * rates$transmission_rate[ tdx ] ) %*% states[[ "infectious" ]] ) 
    flows[[ "SE" ]] <- matrix( stats::rbinom( n_cells, states[[ "susceptible" ]], infectious_contacts ), nrow = n_groups )
    flows[[ "EI" ]] <- matrix( stats::rbinom( n_cells, states[[ "exposed" ]], rates$infectiousness_rate[ tdx ] ), nrow = n_groups )
    flows[[ "IR" ]] <- matrix( stats::rbinom( n_cells, states[[ "infectious" ]], rates$recovery_rate[ tdx ] ), nrow = n_groups )
    
    # update the number in each state
    states[[ "susceptible" ]] <- states[[ "susceptible" ]] - flows[[ "SE" ]]
    states[[ "exposed" ]]     <- states[[ "exposed" ]] + flows[[ "SE" ]] - flows[[ "EI" ]]
    states[[ "infectious" ]]  <- states[[ "infectious" ]] + flows[[ "EI" ]] - flows[[ "IR" ]]
    states[[ "recovered" ]]   <- states[[ "recovered" ]] + flows[[ "IR" ]] 
    
    # record output in a new data table
    outputs[[ tdx+1 ]] <- data.table::copy( output_template )
    outputs[[ tdx+1 ]][ , time := tdx ]
    outputs[[ tdx+1 ]][ , value := unlist( lapply( states, as.vector ) )]
  }
  outputs <- data.table::rbindlist( outputs )

  return( outputs )
}

#' .prepare_interventions
#'
#' @description Converts interventions in to time dependent changes in parameters
#' and transmission rates
#' @param contact_matrix The base contact matrix
#' @param rates The base transition rates
#' @param time_end The length of the simulation
#' @param interventions The interventions to be applied
#'
#' @return A list of the contact matrix and transition rates applicable for 
#' each step of the simulation
#' @keywords internal
.prepare_interventions <- function( contact_matrix, rates, time_end, interventions ) {
  # contact reductions
  n_groups <- nrow( contact_matrix )
  contact_reduction <- lapply( 1:time_end, function( x ) rep( 1, nrow( contact_matrix ) ) )
  intervention      <- interventions[[ "contacts" ]] 
  
  if( !is.null( intervention ) ) {
    # check it is of the correct type
    type <- class( intervention )[1]
    
    checkmate::assert( type == "contacts_intervention", name = "contacts_inteverntion expected" )
    checkmate::assert( nrow( intervention$reduction ) == nrow( contact_matrix), name = "incorrect rows in reduction of contact_intervention")
    
    for( jdx in seq( 1, length( intervention$time_begin ) ) ) {
      if( intervention$time_begin[jdx] > intervention$time_end[jdx] ) 
        next;
      
      # apply them
      times    <- seq( intervention$time_begin[jdx], intervention$time_end[jdx] ) 
      red_func <- function( x ) x - as.vector( intervention$reduction[,jdx] )
      contact_reduction[ times ] <- lapply( contact_reduction[ times ], red_func )
    } 
    
    # the cumulative effect of interventions is capped at 100%
    contact_reduction <- lapply( contact_reduction, function( x ) pmax( x, 0 ) )
  }
  # apply reductions
  contact_matrix <- lapply( contact_reduction, function( x ) diag( x, nrow = n_groups ) %*% contact_matrix )
  
  # rate reductions
  for( rate_type in names( rates ) ) {
    rate_reduction <- rep( 1, time_end )
    intervention   <- interventions[[ rate_type ]] 
    if( !is.null( intervention) ) {
      # check it is of the correct type
      type <- class( intervention )[1]
      checkmate::assert( type == "rate_intervention", name = "rate_inteverntion expected" )
      
      for( jdx in seq( 1, length( intervention$time_begin ) ) ) {
        if( intervention$time_begin[jdx] > intervention$time_end[jdx] ) 
          next;
        
        # apply them
        times    <- seq( intervention$time_begin[jdx], intervention$time_end[jdx] ) 
        rate_reduction[ times ] <- rate_reduction[ times ] - intervention$reduction[jdx]
      }   
      
      # the cumulative effect of interventions is capped at 100%
      rate_reduction <- pmax( rate_reduction, 0 )
    }
    # apply reductions
    rates[[ rate_type ]] <- rates[[ rate_type ]] * rate_reduction
  }
  
  return( list( contact_matrix = contact_matrix, rates = rates ) )
}

