#-------------------- Functions to generate mixture cure model survival data 
data_generation <- function(n, betas, lambdas, 
                            gamma, end_accrual_after_period, end_study_after_period,ts,link_function) {
  #' @param n: taille d'echantillon
  #' @param betas: (vector) Coefficients des covariables pour le modèle de cure
  #' @param lambdas: (vector) Coefficients des covariables pour 
  #' la fonction de survie de Weibull
  #' @param gamma: shape
  #' @param end_accrual_after_period: Durée de l'inclusion
  #' @param end_study_after_period: Durée totale de l'étude 
  #' @export
  
  # Covariate
  X <- cbind(rep(1, n), stats::rbinom(n, 1, prob = 0.5))
  
  # Proportion of cured
  Xbetas <- as.numeric(X %*% betas)
  pi <- switch(link_function,
               "logit" = exp(Xbetas) / (1 + exp(Xbetas)),
               "loglog" = exp(-exp(Xbetas)),
               "probit" = pnorm(Xbetas),
               stop("Link function not supported"))
  
  
  # Susceptibility status
  U <- stats::rbinom(n, size = 1, prob = 1 - pi) # Cure status U = 1 --> uncured
  
  #-----Survival times for uncured (U=1) subject from Weibull PH
  lambda <- exp(X %*% lambdas)
  
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
  df <- data.frame(tte = tte,evt = evt, group = X[,2],cure=U)
  
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
  
  S_u<- pweibull(whichTau,shape=shape,scale= scale^{-1/shape},lower.tail=FALSE)
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
  
  lambda0<-unlist(cm$coefs)[3]
  lambda1<-unlist(cm$coefs)[4]
  beta0 <- unlist(cm$coefs)[1]
  beta1 <- unlist(cm$coefs)[2]
  
  # Compute ratio test for Group 0 (control group)
  whichTau_A <- max(data$tte[data$group == 0])
  pi_A <- exp(beta0) / (1 + exp(beta0))
  scale_A <- exp(lambda0)
  shape <- exp(unlist(cm$coefs)[5])
  
  ratio_test_A <- compute_rn(whichTau_A, shape, scale_A, pi_A)
  
  # Compute ratio test for Group 1 (treatment group)
  whichTau_B <- max(data$tte[data$group == 1])
  pi_B <- exp(beta0 + beta1) / (1 + exp(beta0 + beta1))
  scale_B <- exp(lambda0 + lambda1)
  
  ratio_test_B <- compute_rn(whichTau_B, shape, scale_B, pi_B)
  
  return(list(rn_a = ratio_test_A, pi_a = pi_A ,rn_b = ratio_test_B, pi_b = pi_B))
}

#--------------------  Fit mixture cure model and test appropriateness of the model through RECeUS
run_test <- function(data,link_function) {
  cm <- fit.cure.model(Surv(tte,evt)~group,formula.surv =list(~group,~1) , data=data, dist = "weibull",link=link_function)
  results <- compute_ratio_test(data, cm)
  
  # Both conditions respected
  test_conclusion1 <- with(results, as.integer( (pi_a > 0.025 & rn_a < 0.05) & (pi_b > 0.05 & rn_b < 0.05)) )
  
  # At least one condition respected
  test_conclusion2 <- with(results, as.integer( (pi_a > 0.025 & rn_a < 0.05) | (pi_b > 0.05 & rn_b < 0.15)) )
  
  # On average 
  test_conclusion3 <- with(results, as.integer( (mean(c(pi_a,pi_b)) > 0.025 & mean(c(rn_a,rn_b)) < 0.05)) )
  
  # Overall conclusion
  cm_overall <- fit.cure.model(Surv(tte,evt)~1, data=data, dist = "weibull",link=link_function)
  # Compute ratio test for overall population
  whichTau <- max(data$tte)
  beta<-unlist(cm_overall$coefs)[1]
  pi <- exp(beta) / (1 + exp(beta))
  scale <- exp(unlist(cm_overall$coefs)[2])
  shape <- exp(unlist(cm_overall$coefs)[3])
  
  r_n <- compute_rn(whichTau, shape, scale, pi)
  test_conclusion4 <- as.integer(pi > 0.025 & r_n < 0.05)
  
  return(list(both = test_conclusion1, at_least_one = test_conclusion2, on_avg = test_conclusion3, overall=test_conclusion4))
}

run_simulations <- function(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, link_function) {
  set.seed(123)
  simulations <- lapply(1:nb_simulations, function(i) {
    data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts, link_function)
    run_test(data, link_function)
  })
  
  sum_results <- function(simulations, test_name) {
    sapply(simulations, function(res) res[[test_name]]) %>% sum(na.rm = T)
  }
  
  both_sum <- sum_results(simulations, "both")
  at_least_one_sum <- sum_results(simulations, "at_least_one")
  on_avg_sum <- sum_results(simulations, "on_avg")
  overall_sum <- sum_results(simulations, "overall")
  
  results_df <- data.frame(
    test_conclusion = c("both", "at_least_one", "on_avg","overall"),
    sum = c(both_sum, at_least_one_sum, on_avg_sum,overall_sum)
  )
  
  
  return(results_df)
}

#To compute cure for different settings
#Settings
# pi0<-c(1.05, -0.77) %*% c(1, 0)
# exp(pi0)/(1+exp(pi0))
# exp(-exp(pi0))
# pnorm(pi0)
# 
# pi1<-c(1.05, -0.77) %*% c(1, 1)
# exp(pi1)/(1+exp(pi1))
# exp(-exp(pi1))
# pnorm(pi1)
# 
# pi0<-c(1.05, 0.06) %*% c(1, 0)
# exp(pi0)/(1+exp(pi0))
# exp(-exp(pi0))
# pnorm(pi0)
# 
# pi1<-c(1.05, 0.06) %*% c(1, 1)
# exp(pi1)/(1+exp(pi1))
# exp(-exp(pi1))
# pnorm(pi1)

#--------------------- Compute required sample size in graaph dataset
# library(NPHMC)
# library(dplyr)
# 
# #covidtest<-covid %>%mutate(X=as.numeric(as.factor(covid$DXM))-1) %>% select(survie,dc,X)
# graaphtest<-qaly %>%mutate(X=as.numeric(as.factor(qaly$R1))-1) %>% select(delpfs,pfs,X)
# 
# # /!\ lambda0=lambda**(1/gamma)
# 
# NPHMC(power=0.8,alpha = 0.05,accrualtime = 84, followuptime = 84, p=0.5,
#       accrualdist = "uniform", hazardratio = 1.06, oddsratio = .46, pi0 = 0.74,
#       survdist = "weib", k=0.16,lambda0 = exp(-3.36+0.06)**(1/0.16))
# 
# NPHMC(power=0.80,alpha=0.05,accrualtime=40,followuptime=84,p=0.5,accrualdist="uniform",survdist = "weib",
#        data=graaphtest)
