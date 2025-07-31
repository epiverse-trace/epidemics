library("epidemics")
library("socialmixr")
library("withr")

#### Set up population characteristics ####
# load contact and population data from socialmixr::polymod
polymod <- socialmixr::polymod
contact_data <- socialmixr::contact_matrix(
  polymod,
  countries = "United Kingdom",
  age.limits = c(0, 20, 40),
  symmetric = TRUE
)

# prepare contact matrix
contact_matrix <- contact_data[["matrix"]]

# prepare the demography vector
demography_vector <- contact_data$demography$population
names(demography_vector) <- colnames(contact_matrix)

# initial conditions
initial_i <- 1e-6
initial_conditions <- c(
  S = 1 - initial_i, E = 0, I = initial_i, R = 0, V = 0
)

# build for all age groups
initial_conditions <- rbind(
  initial_conditions,
  initial_conditions,
  initial_conditions
)

# assign rownames for clarity
rownames(initial_conditions) <- rownames(contact_matrix)

# create population object
uk_population <- population(
  name = "UK",
  contact_matrix = contact_matrix,
  demography_vector = demography_vector,
  initial_conditions = initial_conditions
)

# create vector of parameters
beta <- with_seed(
  1,
  rnorm(100, mean = 1.3 / 7, sd = 0.005)
)

# create list of interventions from vignette
max_time <- 600
# prepare durations as starting at 25% of the way through an epidemic
# and ending halfway through
time_begin <- max_time / 4
time_end <- max_time / 2

# create three distinct contact interventions
# prepare an intervention that models school closures for 180 days
close_schools <- intervention(
  name = "School closure",
  type = "contacts",
  time_begin = time_begin,
  time_end = time_end,
  reduction = matrix(c(0.3, 0.01, 0.01))
)

# prepare an intervention which mostly affects adults 20 -- 65
close_workplaces <- intervention(
  name = "Workplace closure",
  type = "contacts",
  time_begin = time_begin,
  time_end = time_end,
  reduction = matrix(c(0.01, 0.3, 0.01))
)

# prepare a combined intervention
combined_intervention <- c(close_schools, close_workplaces)

# create a mask-mandate rate intervention
# prepare an intervention that models mask mandates for 300 days
mask_mandate <- intervention(
  name = "mask mandate",
  type = "rate",
  time_begin = time_begin,
  time_end = time_end,
  reduction = 0.1
)

# create intervention sets, which are combinations of contacts and rate
# interventions
intervention_scenarios <- list(
  baseline = NULL,
  scenario_01 = list(
    contacts = close_schools
  ),
  scenario_02 = list(
    contacts = close_workplaces
  ),
  scenario_03 = list(
    contacts = combined_intervention
  ),
  scenario_04 = list(
    transmission_rate = mask_mandate
  ),
  scenario_05 = list(
    contacts = close_schools,
    transmission_rate = mask_mandate
  ),
  scenario_06 = list(
    contacts = close_workplaces,
    transmission_rate = mask_mandate
  ),
  scenario_07 = list(
    contacts = combined_intervention,
    transmission_rate = mask_mandate
  )
)
