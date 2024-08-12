setwd("/home/ckouadio/Documents/")
source("functions1.R")

findbeta1 <- function(OR) {
  beta1 <- log(OR)
  return(beta1)
}

#-------------------- Fixed parameters
library(cuRe)
library(dplyr)
library(parallel)  # Charger le package parallel

n <- 280
nb_simulations <- 10000
gamma <- exp(0.16)
end_accrual_after_period <- 40
end_study_after_period <- 84
ts <- 20000

#------------------- Simulating variation of HR
lambdas <- c(-3.36, 0.06)
OR_values <- seq(0.5, 1.5, 0.05)

# Fonction de simulation pour un OR donné
simulate_for_OR <- function(OR_value) {
  beta1 <- findbeta1(OR_value)
  betas <- c(1.05, beta1)
  
  results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts)
  null_results <- count_null_results(results)
  test_results <- count_test_results(results)
  
  return(list(
    OR = OR_value,
    results = results,
    null_results = null_results,
    test_results = test_results
  ))
}

# Nombre de cœurs à utiliser
num_cores <- 10

# Utiliser mclapply pour paralléliser les simulations
all_results <- mclapply(OR_values, simulate_for_OR, mc.cores = num_cores)

saveRDS(all_results, "OR_variation.rds")
