# Time-dependent parameters
beta_t <- interpolate(time, beta, "linear")
r_t <- interpolate(time, r, "linear") # report
eta_t <- interpolate(time, eta, "linear") # hospitalisation
sigma_t <- interpolate(time, sigma, "linear")
tau1_t <- interpolate(time, tau1, "linear") # hospital entry
tau2_t <- interpolate(time, tau2, "linear") # hospital exit
gamma_t <- interpolate(time, gamma, "linear")
pop_change_t[] <- interpolate(time, pop_change, "linear")

## Compute reduction in transmission from rate_intervention
# Contact matrix with time-dependent interventions
# Then add across interventions - rate does not seem to be age-dependent,
# so could remove age dimension for rate_reduction and rate_reduction_total
rate_reduction[, ] <-
  if (t > rate_intervention_start[i] && t < rate_intervention_end[i]) {
    rate_intervention_effect[i, j]
  } else {
    0
  } # nolint: line_length_linter.
rate_reduction_sum[] <- sum(rate_reduction[, i])
rate_reduction_total[] <- 1 - min(rate_reduction_sum[i], 1)

# Specify how transmission varies over time
# FOI is contacts * infectious * transmission rate
# returns a matrix, must be converted to a vector/1D array
# multiply FOI by  reduction in contact origin * reduction in contact destination * reduction in rate # nolint: line_length_linter.
# NOTE: division by population size - this is typically included in the
# contact matrix scaling in other models
lambda_prod[, ] <- beta_t * I[j] * rate_reduction_total[i] /
  (sum(S[]) + sum(E[]) + sum(I[]) + sum(H[]) + sum(R[]))
lambda[] <- sum(lambda_prod[i, ])

# ODEs
deriv(S[]) <- -lambda[i] * S[i] + pop_change_t[i]
deriv(E[]) <- lambda[i] * S[i] - sigma_t * E[i]
deriv(I[]) <- sigma_t * E[i] - eta_t * tau1_t * r_t * I[i] - gamma_t * I[i]
deriv(H[]) <- tau1_t * eta_t * r_t * I[i] - tau2_t * H[i]
deriv(R[]) <- gamma_t * I[i] + tau2_t * H[i]

# Initial conditions
initial(S[]) <- init_S[i]
initial(E[]) <- init_E[i]
initial(I[]) <- init_I[i]
initial(H[]) <- init_H[i]
initial(R[]) <- init_R[i]

# User-defined parameters
n_age <- user()
n_rate_intervention <- user()
beta[] <- user()
eta[] <- user()
r[] <- user()
sigma[] <- user()
gamma[] <- user()
tau1[] <- user()
tau2[] <- user()
pop_change[, ] <- user()

# Not currently used (can re-add for time-dependent parameters)
dim(time) <- n_time
dim(beta) <- n_time
dim(eta) <- n_time
dim(r) <- n_time
dim(sigma) <- n_time
dim(gamma) <- n_time
dim(tau1) <- n_time
dim(tau2) <- n_time
dim(pop_change) <- c(n_time, n_age)
dim(pop_change_t) <- n_age
n_time <- user()
time[] <- user()
rate_intervention_start[] <- user()
rate_intervention_end[] <- user()
rate_intervention_effect[, ] <- user()
init_S[] <- user()
init_E[] <- user()
init_I[] <- user()
init_H[] <- user()
init_R[] <- user()

# Dimensions
dim(lambda_prod) <- c(n_age, n_age)
dim(lambda) <- n_age
dim(S) <- n_age
dim(E) <- n_age
dim(I) <- n_age
dim(H) <- n_age
dim(R) <- n_age
dim(rate_reduction) <- c(n_rate_intervention, n_age)
dim(rate_reduction_sum) <- c(n_age)
dim(rate_reduction_total) <- c(n_age)
dim(rate_intervention_start) <- n_rate_intervention
dim(rate_intervention_end) <- n_rate_intervention
dim(rate_intervention_effect) <- c(n_rate_intervention, n_age)
dim(init_S) <- n_age
dim(init_E) <- n_age
dim(init_I) <- n_age
dim(init_H) <- n_age
dim(init_R) <- n_age
