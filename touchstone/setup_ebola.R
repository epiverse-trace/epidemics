library("epidemics")
library(withr)
library(future)

#### Set up population characteristics ####
# for Guinea
population_size <- 14e6

# prepare initial conditions as proportions
initial_conditions <- c(
  S = population_size - 11, E = 10, I = 1, H = 0, F = 0, R = 0
) / population_size

guinea_population <- population(
  name = "Guinea",
  contact_matrix = matrix(1), # note dummy value
  demography_vector = 14e6, # 14 million, no age groups
  initial_conditions = matrix(
    initial_conditions,
    nrow = 1
  )
)

# allow multithreading
future::plan("multicore", workers = 4)
