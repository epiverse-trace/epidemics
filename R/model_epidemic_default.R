
#' SIR-V epidemic model.
#'
#' @param t The times.
#' @param init The initial conditions.
#' @param params The parameters.
#'
#' @return A list with numeric elements, giving the status of each compartment.
#' @keywords internal
#'
epidemic_default <- function(t, init, params) {
  n_age <- nrow(params[["contact_matrix"]])

  # operate only on the recovered and vaccinated
  S_ <- init[seq(1, n_age)]
  I_ <- init[seq(n_age + 1, 2 * n_age)]

  # scale the contact matrix if within the intervention period
  if (t >= params[["t_npi_begin"]] && t <= params[["t_npi_end"]]) {
    contact_matrix_ <- params[["contact_matrix"]] *
      (1 - params[["npi_contact_reduction"]])
  } else {
    contact_matrix_ <- params[["contact_matrix"]]
  }

  # define ODEs
  # change in susceptibles
  dS <- -(params[["beta"]] * S_) *
    as.vector(contact_matrix_ %*% I_) - # infection
    (params[["nu"]] * (t > (params[["t_vax_begin"]])) * S_) # vaccination

  # change in infecteds/infectious-es
  dI <- (params[["beta"]] * S_) *
    as.vector(contact_matrix_ %*% I_) - # new infections
    (params[["gamma"]] * I_) # new recoveries

  # change in recovereds
  dR <- params[["gamma"]] * I_ # new recoveries

  # change in vaccinateds
  dV <- params[["nu"]] * (t > params[["t_vax_begin"]]) * S_

  # return a list
  list(c(dS, dI, dR, dV))
}
