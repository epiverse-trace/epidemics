# Time-dependent parameters
beta_t <- interpolate(time, beta, "linear")
sigma_t <- interpolate(time, sigma, "linear")
gamma_t <- interpolate(time, gamma, "linear")

## Compute reduction in transmission from contacts_intervention
# Contact matrix with time-dependent interventions
# Then add across interventions - need to check this
contact_reduction[, ] <- if (t > contact_intervention_start[i] && t < contact_intervention_end[i]) 
  contact_intervention_effect[i, j] else 0 # nolint: line_length_linter.
contact_reduction_sum[] <- sum(contact_reduction[, i])
contact_reduction_total[] <- 1 - min(contact_reduction_sum[i],1)


## Compute reduction in transmission from rate_intervention
# Contact matrix with time-dependent interventions
# Then add across interventions - rate does not seem to be age-dependent, 
# so could remove age dimension for rate_reduction and rate_reduction_total
rate_reduction[,] <- if (t > rate_intervention_start[i] && t < rate_intervention_end[i]) 
  rate_intervention_effect[i,j] else 0 # nolint: line_length_linter.
rate_reduction_sum[] <- sum(rate_reduction[,i])
rate_reduction_total[] <- 1 - min(rate_reduction_sum[i], 1)

# Specify how transmission varies over time
# FOI is contacts * infectious * transmission rate 
# returns a matrix, must be converted to a vector/1D array
# multiply FOI by  reduction in contact origin * reduction in contact destination * reduction in rate
lambda_prod[, ] <- C[i, j] * I[j] * beta_t * contact_reduction_total[j] * contact_reduction_total[i] * 
  rate_reduction_total[i]
lambda[] <- sum(lambda_prod[i, ])

# Vaccination - indexing over age groups
vax_rate[] <- if (t >= vax_start[i] && t <= vax_end[i]) vax_nu[i]/S[i] else 0

# ODEs
deriv(S[]) <- if(vax_rate[i] + lambda[i] < 1) -(lambda[i] + vax_rate[i]) * S[i] else - S[i]
deriv(E[]) <- lambda[i] * S[i] - sigma_t * E[i]
deriv(I[]) <- sigma_t * E[i] - gamma_t * I[i]
deriv(R[]) <- gamma_t * I[i]
deriv(V[]) <- if(vax_rate[i] + lambda[i] < 1) vax_rate[i] * S[i] else S[i]*(1-lambda[i])

# Initial conditions
initial(S[]) <- init_S[i]
initial(E[]) <- init_E[i]
initial(I[]) <- init_I[i]
initial(R[]) <- init_R[i]
initial(V[]) <- init_V[i]

# User-defined parameters
C[, ] <- user()
n_age <- user()
n_contact_intervention <- user()
n_rate_intervention <- user()
beta[] <- user()
sigma[] <- user()
gamma[] <- user()

# Not currently used (can re-add for time-dependent parameters)
dim(time) <- n_time
dim(beta) <- n_time
dim(sigma) <- n_time
dim(gamma) <- n_time
n_time <- user()
time[] <- user()
contact_intervention_start[] <- user()
contact_intervention_end[] <- user()
contact_intervention_effect[, ] <- user()
rate_intervention_start[] <- user()
rate_intervention_end[] <- user()
rate_intervention_effect[, ] <- user()
vax_start[] <- user()
vax_end[] <- user()
vax_nu[] <- user()
init_S[] <- user()
init_E[] <- user()
init_I[] <- user()
init_R[] <- user()
init_V[] <- user()

# Dimensions
dim(C) <- c(n_age, n_age)
dim(lambda_prod) <- c(n_age, n_age)
dim(lambda) <- n_age
dim(S) <- n_age
dim(E) <- n_age
dim(I) <- n_age
dim(R) <- n_age
dim(V) <- n_age
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
dim(vax_rate) <- n_age
dim(vax_start) <- n_age
dim(vax_end) <- n_age
dim(vax_nu) <- n_age
dim(init_S) <- n_age
dim(init_E) <- n_age
dim(init_I) <- n_age
dim(init_R) <- n_age
dim(init_V) <- n_age
