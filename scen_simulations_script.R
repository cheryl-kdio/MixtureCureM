setwd("/home/ckouadio/Documents/")
source("functions1.R")

library(cuRe)
library(dplyr)
n<-280
nb_simulations <- 10000
gamma <- exp(0.16)
end_accrual_after_period <- 40
end_study_after_period <- 84
ts <- 20000
linkfunction <- "logit"

print("------------------------------------Scenario 0------------------------------------\n")
# parameters
lambdas <- c(-3.36,0)
betas <- c(1.05,0)

results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
null_results <- count_null_results(results)
test_results <- count_test_results(results)

save_data <- list(
  results = results,
  null_results = null_results,
  test_results = test_results
)

# Save to an RDS file
saveRDS(save_data, file = "scenario0.rds")

print("------------------------------------Scenario 1-bis------------------------------------\n")
# parameters
lambdas <- c(0.28,0.76)
betas <- c(1.05,-0.77)

results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
null_results <- count_null_results(results)
test_results <- count_test_results(results)

save_data <- list(
  results = results,
  null_results = null_results,
  test_results = test_results
)

# Save to an RDS file
saveRDS(save_data, file = "scenario1bis.rds")


# print("------------------------------------Scenario 1------------------------------------\n")
# # parameters
# lambdas <- c(-3.36,0.06)
# betas <- c(1.05,-0.77)
# 
# results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
# null_results <- count_null_results(results)
# test_results <- count_test_results(results)
# 
# save_data <- list(
#   results = results,
#   null_results = null_results,
#   test_results = test_results
# )
# 
# # Save to an RDS file
# saveRDS(save_data, file = "scenario1.rds")
# 
# 
# print("------------------------------------Scenario 2------------------------------------\n")
# lambdas <- c(-3.36,-0.6)
# betas <- c(1.05,0.06)
# 
# results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
# null_results <- count_null_results(results)
# test_results <- count_test_results(results)
# 
# save_data <- list(
#   results = results,
#   null_results = null_results,
#   test_results = test_results
# )
# 
# # Save to an RDS file
# saveRDS(save_data, file = "scenario2.rds")
# 
# print("------------------------------------Scenario 3------------------------------------\n")
# # parameters
# lambdas <- c(-3.36,0.06)
# betas <- c(1.05,0.06)
# 
# results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
# null_results <- count_null_results(results)
# test_results <- count_test_results(results)
# 
# save_data <- list(
#   results = results,
#   null_results = null_results,
#   test_results = test_results
# )
# 
# # Save to an RDS file
# saveRDS(save_data, file = "scenario3.rds")
# 
# print("------------------------------------Scenario 4------------------------------------\n")
# lambdas <- c(-3.36,-0.6)
# betas <- c(1.05,-0.77)
# 
# results <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,linkfunction)
# null_results <- count_null_results(results)
# test_results <- count_test_results(results)
# 
# save_data <- list(
#   results = results,
#   null_results = null_results,
#   test_results = test_results
# )
# 
# # Save to an RDS file
# saveRDS(save_data, file = "scenario4.rds")