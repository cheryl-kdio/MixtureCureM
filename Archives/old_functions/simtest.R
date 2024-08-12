#------------------------ 
source("functions.R")

library(cuRe)
library(dplyr)
n<-280
nb_simulations <- 10000
gamma <- exp(0.16)
end_accrual_after_period <- 40
end_study_after_period <- 84
ts <- 20000
#------------------- Scenario 1
# parameters
lambdas <- c(-3.36,0.06)

print("------------------------------------Scenario 1, logit")
betas <- c(1.05,-0.77)
result_logit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "logit")
saveRDS(result_logit, "/home/ckouadio/Documents/Scenarios_res/result_logit_scenario1.rds")


print("------------------------------------Scenario 1, loglog")
betas <- c(-1.20,0.62)
result_loglog <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "loglog")
saveRDS(result_loglog, "/home/ckouadio/Documents/Scenarios_res/result_loglog_scenario1.rds")

print("------------------------------------Scenario 1, probit")
betas <- c(0.64,-0.46)
result_probit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "probit")
saveRDS(result_probit, "/home/ckouadio/Documents/Scenarios_res/result_probit_scenario1.rds")

# results_df<-cbind(result_logit,result_loglog[,2],result_probit[,2])
# colnames(results_df)[c(2:4)]<-c("logit","loglog","probit")
# 
# saveRDS(result_logit, "/home/ckouadio/Documents/Scenarios_res/result_logit_scenario1.rds")

#------------------- Scenario 2
# parameters
lambdas <- c(-3.36,-0.6)

print("------------------------------------Scenario 2, logit")
betas <- c(1.05,0.06)
result_logit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "logit")
saveRDS(result_logit, "/home/ckouadio/Documents/Scenarios_res/result_logit_scenario2.rds")


print("------------------------------------Scenario 2, loglog")
betas <- c(-1.20,-0.04)
result_loglog <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "loglog")
saveRDS(result_loglog, "/home/ckouadio/Documents/Scenarios_res/result_loglog_scenario2.rds")

print("------------------------------------Scenario 2, probit")
betas <- c(0.64,0.03)
result_probit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "probit")
saveRDS(result_probit, "/home/ckouadio/Documents/Scenarios_res/result_probit_scenario2.rds")

# results_df<-cbind(result_logit,result_loglog[,2],result_probit[,2])
# colnames(results_df)[c(2:4)]<-c("logit","loglog","probit")
# 
# saveRDS(results_df, "/home/ckouadio/Documents/scenario2.rds")

#------------------- Scenario 3
# parameters
lambdas <- c(-3.36,0.06)

print("------------------------------------Scenario 3, logit")
betas <- c(1.05,0.06)
result_logit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "logit")
saveRDS(result_logit, "/home/ckouadio/Documents/Scenarios_res/result_logit_scenario3.rds")

print("------------------------------------Scenario 3, loglog")
betas <- c(-1.20,-0.04)
result_loglog <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "loglog")
saveRDS(result_loglog, "/home/ckouadio/Documents/Scenarios_res/result_loglog_scenario3.rds")

print("------------------------------------Scenario 3, probit")
betas <- c(0.64,0.03)
result_probit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "probit")
saveRDS(result_probit, "/home/ckouadio/Documents/Scenarios_res/result_probit_scenario3.rds")

# results_df<-cbind(result_logit,result_loglog[,2],result_probit[,2])
# colnames(results_df)[c(2:4)]<-c("logit","loglog","probit")
# 
# saveRDS(results_df, "/home/ckouadio/Documents/scenario3.rds")

#------------------- Scenario 4
# parameters
lambdas <- c(-3.36,-0.6)

print("------------------------------------Scenario 4, logit")
betas <- c(1.05,-0.77)
result_logit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "logit")
saveRDS(result_logit, "/home/ckouadio/Documents/Scenarios_res/result_logit_scenario4.rds")

print("------------------------------------Scenario 4, loglog")
betas <- c(-1.20,0.62)
result_loglog <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "loglog")
saveRDS(result_loglog, "/home/ckouadio/Documents/Scenarios_res/result_loglog_scenario4.rds")

print("------------------------------------Scenario 4, probit")
betas <- c(0.64,-0.46)
result_probit <- run_simulations(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, "probit")
saveRDS(result_probit, "/home/ckouadio/Documents/Scenarios_res/result_probit_scenario4.rds")

# results_df<-cbind(result_logit,result_loglog[,2],result_probit[,2])
# colnames(results_df)[c(2:4)]<-c("logit","loglog","probit")
# 
# saveRDS(results_df, "/home/ckouadio/Documents/scenario4.rds")

