setwd("/home/ckouadio/Documents/")
source("functions1.R")


findlambda1 <- function(HR){
  lambda1 <- log(HR)
  return(lambda1)
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


#------------------- Simulating variation of HR
betas <- c(1.05,0.06)
all_results <- list()
HR_values <- seq(0.5, 1.5,0.05)

for(i in seq_along(HR_values)){
  
  lambda1 <- findlambda1(HR_values[i])
  lambdas <- c(-3.36,lambda1)
  
  
  results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts)
  null_results <- count_null_results(results)
  test_results <- count_test_results(results)
  
  all_results[[i]] <- list(
    pi = HR_values[i],
    results = results,
    null_results = null_results,
    test_results = test_results
  )
}

saveRDS(all_results, "HR_variation1.rds")


