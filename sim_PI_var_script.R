# Configuration initiale
setwd("/home/ckouadio/Documents/")
source("functions1.R")

findbeta1 <- function(pi,beta0){
  beta1 <-  log(pi/(1-pi))-beta0
  return(beta1)
}

#-------------------- Fixed parameters
library(cuRe)
library(dplyr)
n<-280
nb_simulations <- 10000
gamma <- exp(0.16)
end_accrual_after_period <- 40
end_study_after_period <- 84
ts <- 20000


pi_values <- seq(0.01, 0.8, 0.01)
lambdas <- c(-3.36, 0.06)
all_results <- list()

#------------------- Simulating variation of pi
for(i in seq_along(pi_values)) {
  beta1 <- findbeta1(pi_values[i], 1.05)
  betas <- c(1.05, beta1)
  
  results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts)
  null_results <- count_null_results(results)
  test_results <- count_test_results(results)
  
  all_results[[i]] <- list(
    pi = pi_values[i],
    results = results,
    null_results = null_results,
    test_results = test_results
  )
}

# Sauvegarde des résultats
saveRDS(all_results, "cure_variation1.rds")
