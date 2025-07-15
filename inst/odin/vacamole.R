# Time-dependent parameters
beta_t <- interpolate(time, beta, "linear")
beta_vax_t <- interpolate(time, beta_vax, "linear")
eta_t <- interpolate(time, eta, "linear") # hospitalisation
eta_vax_t <- interpolate(time, eta_vax, "linear")
sigma_t <- interpolate(time, sigma, "linear")
gamma_t <- interpolate(time, gamma, "linear")
omega_t <- interpolate(time, omega, "linear")
omega_vax_t <- interpolate(time, omega_vax, "linear")

## Compute reduction in transmission from contacts_intervention
# Contact matrix with time-dependent interventions
# Then add across interventions - need to check this
contact_reduction[, ] <-
  if (t > contact_intervention_start[i] && t < contact_intervention_end[i]) {
    contact_intervention_effect[i, j]
  } else {
    0
  }
contact_reduction_sum[] <- sum(contact_reduction[, i])
contact_reduction_total[] <- 1 - min(contact_reduction_sum[i], 1)


## Compute reduction in transmission from rate_intervention
# Contact matrix with time-dependent interventions
# Then add across interventions - rate does not seem to be age-dependent,
# so could remove age dimension for rate_reduction and rate_reduction_total
rate_reduction[, ] <-
  if (t > rate_intervention_start[i] && t < rate_intervention_end[i]) {
    rate_intervention_effect[i, j]
  } else {
    0
  }
rate_reduction_sum[] <- sum(rate_reduction[, i])
rate_reduction_total[] <- 1 - min(rate_reduction_sum[i], 1)

# Specify how transmission varies over time
# FOI is contacts * infectious * transmission rate
# returns a matrix, must be converted to a vector/1D array
# multiply FOI by  reduction in contact origin * reduction in contact destination * reduction in rate # nolint: line_length_linter.
lambda_prod[, ] <- C[i, j] * (I[j] + IV[j]) * beta_t *
  contact_reduction_total[j] * contact_reduction_total[i] *
  rate_reduction_total[i]
lambda[] <- sum(lambda_prod[i, ])

lambda_prod_vax[, ] <- C[i, j] * (I[j] + IV[j]) * beta_vax_t *
  contact_reduction_total[j] * contact_reduction_total[i] *
  rate_reduction_total[i]
lambda_vax[] <- sum(lambda_prod_vax[i, ])

# Vaccination - indexing over age groups
vax_rate_S[] <- if (t >= vax_start[i, 1] && t <= vax_end[i, 1] && S[i] > 0) {
  vax_nu[i, 1] / S[i]
} else {
  0
}
vax_rate_V1[] <- if (t >= vax_start[i, 2] && t <= vax_end[i, 2] && V1[i] > 0) {
  vax_nu[i, 2] / V1[i]
} else {
  0
}

vax_rate_S[] <- if (vax_rate_S[i] + lambda[i] > 1) {
  1 - lambda[i]
} else {
  vax_rate_S[i]
}
vax_rate_V1[] <- if (vax_rate_V1[i] + lambda[i] > 1) {
  1 - lambda[i]
} else {
  vax_rate_V1[i]
}

# ODEs
deriv(S[]) <- -(lambda[i] + vax_rate_S[i]) * S[i]

deriv(E[]) <- lambda[i] * (S[i] + V1[i]) - sigma_t * E[i]
deriv(EV[]) <- lambda_vax[i] * V2[i] - sigma_t * EV[i]

deriv(I[]) <- sigma_t * E[i] - (eta_t + gamma_t + omega_t) * I[i]
deriv(IV[]) <- sigma_t * EV[i] - (eta_vax_t + gamma_t + omega_vax_t) * IV[i]

deriv(H[]) <- eta_t * I[i] - (gamma_t + omega_t) * H[i]
deriv(HV[]) <- eta_vax_t * IV[i] - (gamma_t + omega_vax_t) * HV[i]

deriv(D[]) <- omega_t * (I[i] + H[i]) + omega_vax_t * (IV[i] + HV[i])
deriv(R[]) <- gamma_t * (I[i] + H[i] + IV[i] + HV[i])

deriv(V1[]) <- vax_rate_S[i] * S[i] - (vax_rate_V1[i] + lambda[i]) * V1[i]

deriv(V2[]) <- vax_rate_V1[i] * V1[i] - lambda_vax[i] * V2[i]

# Initial conditions
initial(S[]) <- init_S[i]
initial(E[]) <- init_E[i]
initial(EV[]) <- init_EV[i]
initial(I[]) <- init_I[i]
initial(IV[]) <- init_IV[i]
initial(H[]) <- init_H[i]
initial(HV[]) <- init_HV[i]
initial(D[]) <- init_D[i]
initial(R[]) <- init_R[i]
initial(V1[]) <- init_V1[i]
initial(V2[]) <- init_V2[i]

# User-defined parameters
C[, ] <- user()
n_age <- user()
n_contact_intervention <- user()
n_rate_intervention <- user()
beta[] <- user()
beta_vax[] <- user()
eta[] <- user()
eta_vax[] <- user()
sigma[] <- user()
gamma[] <- user()
omega[] <- user()
omega_vax[] <- user()

# Not currently used (can re-add for time-dependent parameters)
dim(time) <- n_time
dim(beta) <- n_time
dim(beta_vax) <- n_time
dim(eta) <- n_time
dim(eta_vax) <- n_time
dim(sigma) <- n_time
dim(gamma) <- n_time
dim(omega) <- n_time
dim(omega_vax) <- n_time
n_time <- user()
time[] <- user()
contact_intervention_start[] <- user()
contact_intervention_end[] <- user()
contact_intervention_effect[, ] <- user()
rate_intervention_start[] <- user()
rate_intervention_end[] <- user()
rate_intervention_effect[, ] <- user()
vax_start[, ] <- user()
vax_end[, ] <- user()
vax_nu[, ] <- user()
init_S[] <- user()
init_E[] <- user()
init_EV[] <- user()
init_I[] <- user()
init_IV[] <- user()
init_H[] <- user()
init_HV[] <- user()
init_D[] <- user()
init_R[] <- user()
init_V1[] <- user()
init_V2[] <- user()

# Dimensions
dim(C) <- c(n_age, n_age)
dim(lambda_prod) <- c(n_age, n_age)
dim(lambda_prod_vax) <- c(n_age, n_age)
dim(lambda) <- n_age
dim(lambda_vax) <- n_age
dim(S) <- n_age
dim(E) <- n_age
dim(EV) <- n_age
dim(I) <- n_age
dim(IV) <- n_age
dim(H) <- n_age
dim(HV) <- n_age
dim(D) <- n_age
dim(R) <- n_age
dim(V1) <- n_age
dim(V2) <- n_age
dim(contact_reduction) <- c(n_contact_intervention, n_age)
dim(contact_reduction_sum) <- c(n_age)
dim(contact_reduction_total) <- c(n_age)
dim(contact_intervention_start) <- n_contact_intervention
dim(contact_intervention_end) <- n_contact_intervention
dim(contact_intervention_effect) <- c(n_contact_intervention, n_age)
dim(rate_reduction) <- c(n_rate_intervention, n_age)
dim(rate_reduction_sum) <- c(n_age)
dim(rate_reduction_total) <- c(n_age)
dim(rate_intervention_start) <- n_rate_intervention
dim(rate_intervention_end) <- n_rate_intervention
dim(rate_intervention_effect) <- c(n_rate_intervention, n_age)
dim(vax_rate_S) <- n_age
dim(vax_rate_V1) <- n_age
dim(vax_start) <- c(n_age, 2)
dim(vax_end) <- c(n_age, 2)
dim(vax_nu) <- c(n_age, 2)
dim(init_S) <- n_age
dim(init_E) <- n_age
dim(init_EV) <- n_age
dim(init_I) <- n_age
dim(init_IV) <- n_age
dim(init_H) <- n_age
dim(init_HV) <- n_age
dim(init_D) <- n_age
dim(init_R) <- n_age
dim(init_V1) <- n_age
dim(init_V2) <- n_age
