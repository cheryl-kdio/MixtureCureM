#-------------------- Functions to generate mixture cure model survival data 
data_generation <- function(n, betas, lambdas, 
                            gamma, end_accrual_after_period, end_study_after_period,ts) {
  #' @description
  #' Generate survival data for a mixture cure model
  #' @param n: (int) Nombre de participants dans l'étude
  #' @param betas: (numeric vector) Coefficients des covariables pour le modèle de cure
  #' @param lambdas: (numeric vector) Coefficients des covariables pour 
  #' la fonction de survie de Weibull
  #' @param gamma: (numeric) Paramètre de forme pour la distribution de Weibull
  #' @param end_accrual_after_period: (numeric) Durée de l'inclusion
  #' @param end_study_after_period: (numeric) Durée totale de l'étude en mois
  #' 
  #' @export
  
  # Covariate
  X <- cbind(rep(1, n), stats::rbinom(n, 1, prob = 0.5))
  
  # Proportion of cured
  Xbetas <- as.numeric(X %*% betas)
  pi <- exp(Xbetas) / (1 + exp(Xbetas))
  
  # Generation of susceptibility status
  U <- stats::rbinom(n, size = 1, prob = 1 - pi) # Cure status U = 1 --> uncured
  
  # Generation of survival times for uncured (U=1) subject from Weibull PH
  lambda <- exp(X %*% lambdas)
  
  # Generate survival times
  # Tlat <- rweibull(n, shape = gamma, scale = (1 / lambda)^(1/gamma))
  # or
  Unif <- runif(n,0,1)
  Tlat <- (-log(Unif)/lambda)^(1 / gamma)
  
  # Very large survival time for cured subjects
  Tlat[which(U == 0)] <- ts
  tobs <- Tlat
  
  # Accrual time
  accrual <- runif(n, 0, end_accrual_after_period)
  time_since_study_beginning <- tobs + accrual 
  
  # Generate censoring
  evt <- ifelse(time_since_study_beginning > end_study_after_period, 0, 1) # censoring
  tte <- ifelse(evt == 1, tobs, end_study_after_period - accrual)
  
  # Create data frame
  df <- data.frame(tte = tte,evt = evt, group = X[,2])
  
  # Return the data frame
  return(df)
}


#-------------------- Formula to compute ratio statistic for mixture cure model
compute_rn <- function (whichTau, shape, scale, pi) {
  #' @description
  #' Formula to compute ratio statistic for mixture cure model
  #' @param whichTau: (numeric) Time of the last event
  #' @param shape: (numeric) Shape parameter of the Weibull distribution
  #' @param scale: (numeric) Scale parameter of the Weibull distribution
  #' @param pi: (numeric) Proportion of cured subjects
  
  S_u<- pweibull(whichTau,shape=shape,scale=(1/scale)**(1/shape),lower.tail=FALSE)
  S <- (pi+(1-pi)*S_u)
  r_n<- S_u / S
  return(r_n)
}

#-------------------- Compute ratio statistic given real data

compute_ratio_test <- function(data, cm) {
  #' @description
  #' Compute ratio statistic for mixture cure model
  #' @param data: (data.frame) Data frame containing the survival data in form of data(tte, evt, group)
  #' @param cm: (list) objects of class cuRe (from the cuRe package)
  
  # Compute ratio test for Group 0 (control group)
  whichTau_A <- max(data$tte[data$group == 0])
  pi_A <- exp(unlist(cm$coefs)[1]) / (1 + exp(unlist(cm$coefs)[1]))
  scale_A <- exp(unlist(cm$coefs)[3])
  shape <- exp(unlist(cm$coefs)[5])
  
  ratio_test_A <- compute_rn(whichTau_A, shape, scale_A, pi_A)
  
  # Compute ratio test for Group 1 (treatment group)
  whichTau_B <- max(data$tte[data$group == 1])
  pi <- sum(c(1, 1, 0, 0, 0) * unlist(cm$coefs))
  pi_B <- exp(pi) / (1 + exp(pi))
  scale_B <- exp(unlist(cm$coefs)[3] + unlist(cm$coefs)[4])
  
  ratio_test_B <- compute_rn(whichTau_B, shape, scale_B, pi_B)
  
  
  return(list(rn_a = ratio_test_A, pi_a = pi_A ,rn_b = ratio_test_B, pi_b = pi_B))
}

#--------------------  Fit mixture cure model and test appropriateness of the model through RECeUS
run_test <- function(data) {
  cm <- fit.cure.model(Surv(tte,evt)~group,formula.surv =list(~group,~1) , data=data, dist = "weibull",link="logit")
  results <- compute_ratio_test(data, cm)
  
  # Both conditions respected
  test_conclusion1 <- with(results, as.integer((pi_a > 0.025 & rn_a < 0.05) & (pi_b > 0.05 & rn_b < 0.05)))
  
  # At least one condition respected
  test_conclusion2 <- with(results, as.integer((pi_a > 0.025 & rn_a < 0.05) | (pi_b > 0.05 & rn_b < 0.15)))
  
  # Only one condition respected
  test_conclusion3 <- with(results, as.integer((pi_a > 0.025 & rn_a < 0.05) != (pi_b > 0.05 & rn_b < 0.05)))
  
  return(list(both = test_conclusion1, at_least_one = test_conclusion2, only_one = test_conclusion3))
}

#-------------------- Simulation step
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
betas <- c(1.05,-0.77)
lambdas <- c(-3.36,0.06)


set.seed(123)
simulations <- lapply(1:nb_simulations, function(i) {
  data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period,ts)
  run_test(data)
})

sum_results <- function(simulations, test_name) {
  sapply(simulations, function(res) res[[test_name]]) %>% sum()
}

both_sum <- sum_results(simulations, "both")
at_least_one_sum <- sum_results(simulations, "at_least_one")
only_one_sum <- sum_results(simulations, "only_one")

results_df <- data.frame(
  test_conclusion = c("both", "at_least_one", "only_one"),
  sum = c(both_sum, at_least_one_sum, only_one_sum)
)

saveRDS(results_df, "/home/ckouadio/Documents/scenario1.rds")

#------------------- Scenario 2
# parameters
betas <- c(-3.36,0.06)
lambdas <- c(1.05,-0.6)


simulations <- lapply(1:nb_simulations, function(i) {
  data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period,ts)
  run_test(data)
})

sum_results <- function(simulations, test_name) {
  sapply(simulations, function(res) res[[test_name]]) %>% sum()
}

both_sum <- sum_results(simulations, "both")
at_least_one_sum <- sum_results(simulations, "at_least_one")
only_one_sum <- sum_results(simulations, "only_one")

results_df <- data.frame(
  test_conclusion = c("both", "at_least_one", "only_one"),
  sum = c(both_sum, at_least_one_sum, only_one_sum)
)

saveRDS(results_df, "/home/ckouadio/Documents/scenario2.rds")


#------------------- Scenario 3
# parameters
betas <- c(-3.36,0.06)
lambdas <- c(-3.36,0.06)


simulations <- lapply(1:nb_simulations, function(i) {
  data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period,ts)
  run_test(data)
})

sum_results <- function(simulations, test_name) {
  sapply(simulations, function(res) res[[test_name]]) %>% sum()
}

both_sum <- sum_results(simulations, "both")
at_least_one_sum <- sum_results(simulations, "at_least_one")
only_one_sum <- sum_results(simulations, "only_one")

results_df <- data.frame(
  test_conclusion = c("both", "at_least_one", "only_one"),
  sum = c(both_sum, at_least_one_sum, only_one_sum)
)

saveRDS(results_df, "/home/ckouadio/Documents/scenario3.rds")

#------------------- Scenario 4
# parameters
betas <- c(1.05,-0.77)
lambdas <- c(1.05,-0.6)

simulations <- lapply(1:nb_simulations, function(i) {
  data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period,ts)
  run_test(data)
})

sum_results <- function(simulations, test_name) {
  sapply(simulations, function(res) res[[test_name]]) %>% sum()
}

both_sum <- sum_results(simulations, "both")
at_least_one_sum <- sum_results(simulations, "at_least_one")
only_one_sum <- sum_results(simulations, "only_one")

results_df <- data.frame(
  test_conclusion = c("both", "at_least_one", "only_one"),
  sum = c(both_sum, at_least_one_sum, only_one_sum)
)

saveRDS(results_df, "/home/ckouadio/Documents/scenario4.rds")
