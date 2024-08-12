###### Variation pi dans bras expérimental

setwd("/home/ckouadio/Documents/")
source("functions1.R")

library(cuRe)
library(dplyr)
library(parallel)  # Charger le package parallel

n<-280
nb_simulations <- 10000
gamma <- exp(0.16)
end_accrual_after_period <- 40
end_study_after_period <- 84
ts <- 20000
linkfunction <- "given"

print("------------------------------------Scenario ------------------------------------\n")
# parameters
lambdas <- c(-3.36,0.06)
pi_values <- seq(0.025,0.8,0.025)

# Fonction de simulation pour un pi
simulate_for_pi <- function(pi_values){
  betas <- c(0, pi_values)
  
  results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
  null_results <- count_null_results(results)
  test_results <- count_test_results(results)
  
  return(list(
    pi = pi_values,
    results = results,
    null_results = null_results,
    test_results = test_results
  ))
}

# Nombre de cœurs à utiliser
num_cores <- 10

# Utiliser mclapply pour paralléliser les simulations
all_results <- mclapply(pi_values, simulate_for_pi, mc.cores = num_cores)

saveRDS(all_results, "scenario_variation.rds")
