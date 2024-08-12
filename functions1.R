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
               "given" = as.numeric(X %*% betas),
               stop("Link function not supported"))
  
  # Susceptibility status
  U <- stats::rbinom(n, size = 1, prob = 1 - pi) # Cure status U = 1 --> uncured
  
  #-----Survival times for uncured (U=1) subject from Weibull PH
  lambda <- exp(X %*% lambdas)
  Unif <- runif(n,0,1)
  Tlat <- (-log(Unif)/lambda)^(1 / gamma)
  Tlat[which(U == 0)] <- ts
  tobs <- Tlat
  
  # Accrual time
  accrual <- runif(n, 0, end_accrual_after_period)
  time_since_study_beginning <- tobs + accrual 
  
  # Generate censoring
  evt <- ifelse(time_since_study_beginning > end_study_after_period, 0, 1) # censoring
  tte <- ifelse(evt == 1, tobs, end_study_after_period - accrual)
  
  df <- data.frame(tte = tte,evt = evt, group = X[,2],cure=U)
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
  if(!is.null(cm$covariance) && cm$optim$convergence==0 ){
    results <- compute_ratio_test(data, cm)
    
    test_conclusion1 <- with(results, as.integer( (pi_a > 0.025 & rn_a < 0.05) & (pi_b > 0.05 & rn_b < 0.05)) )
    test_conclusion2 <- with(results, as.integer( (pi_a > 0.025 & rn_a < 0.05) | (pi_b > 0.05 & rn_b < 0.15)) )
    test_conclusion3 <- with(results, as.integer( (mean(c(pi_a,pi_b)) > 0.025 & mean(c(rn_a,rn_b)) < 0.05)) )
    
    cm_overall <- fit.cure.model(Surv(tte,evt)~1, data=data, dist = "weibull",link=link_function)
    whichTau <- max(data$tte)
    beta<-unlist(cm_overall$coefs)[1]
    pi <- exp(beta) / (1 + exp(beta))
    scale <- exp(unlist(cm_overall$coefs)[2])
    shape <- exp(unlist(cm_overall$coefs)[3])
    r_n <- compute_rn(whichTau, shape, scale, pi)
    test_conclusion4 <- as.integer(pi > 0.025 & r_n < 0.05)
    
    return(list(both = test_conclusion1, at_least_one = test_conclusion2, on_avg = test_conclusion3, overall=test_conclusion4))
  }else {
    return(NULL)
  }
}

run_simulations <- function(nb_simulations, n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,link_function="logit") {
  set.seed(123)
  print("----------- Running simulations --------------- \n")
  
  simulations <- lapply(1:nb_simulations, function(i) {
    data <- data_generation(n, betas, lambdas, gamma, end_accrual_after_period, end_study_after_period, ts,link_function)
    logit_results <- run_test(data, "logit")
    loglog_results <- run_test(data, "loglog")
    probit_results <- run_test(data, "probit")
    
    list(data = data,logit = logit_results, loglog = loglog_results, probit = probit_results)
  })
  print("----------- End of simulations --------------- \n")
  
  return(simulations)
}


count_null_results <- function(simulations) {
  ###### Function to summarize the null results
  count_null_results <- function(simulations, test_name) {
    sum(sapply(simulations, function(res) is.null(res[[test_name]])))
  }
  
  print("----------- Starting counting null results --------------- \n")
  logit_null_count <- count_null_results(simulations, "logit")
  loglog_null_count <- count_null_results(simulations, "loglog")
  probit_null_count <- count_null_results(simulations, "probit")
  
  
  null_list <- list(logit = logit_null_count, loglog = loglog_null_count, probit = probit_null_count)
  return(null_list)
  
}

count_test_results <- function(simulations) {
  ###### Function to summarize the test results
  print("----------- Starting counting test results --------------- \n")
  
  sum_results <- function(simulations, link, conclusions) {
    sums <- setNames(rep(0, length(conclusions)), conclusions)
    for (conclusion in conclusions) {
      sums[conclusion] <- sum(sapply(simulations, function(res) {
        if (!is.null(res[[link]]) && !is.null(res[[link]][[conclusion]])) {
          return(res[[link]][[conclusion]])
        } else {
          return(NA)
        }
      }), na.rm = TRUE)
    }
    return(sums)
  }
  
  conclusions_to_sum <- c("both", "on_avg", "at_least_one", "overall")
  sum_logit <- sum_results(simulations, "logit", conclusions_to_sum)
  sum_loglog <- sum_results(simulations, "loglog", conclusions_to_sum)
  sum_probit <- sum_results(simulations, "probit", conclusions_to_sum)
  
  summary_list <- list(
    logit = sum_logit,
    loglog = sum_loglog,
    probit = sum_probit
  )
  return(summary_list)
}


### Alpha test for sufficent follow-up 
alphatest <- function(data){
  Tn <- max(data$tte) #tn
  Tn_a <- max(data$tte[data$evt==1]) #tn_tilde
  
  Nn <- sum(data$evt==1 & data$tte>(2*Tn_a-Tn) & data$tte <=Tn)
  n <- nrow(data)
  alpha <- (1-Nn/n)**n
  return(alpha)
}


## Deviance test for presence of immune
deviance_test <- function(data,cm,conf_level){
  
  
  ll.mix <- function(pi_=NULL,tte,evt,group,cm) {
    if (cm$all.formulas[[1]]!="Surv(tte, evt) ~ 1"){
      lambda0<-unlist(cm$coefs)[3]
      lambda1<-unlist(cm$coefs)[4]
      beta0 <- unlist(cm$coefs)[1]
      beta1 <- unlist(cm$coefs)[2]
      
      if(is.null(pi_)){
        pi_ <- exp(beta0 + beta1*group) / (1 + exp(beta0 + beta1*group))
      }
      shape <- exp(unlist(cm$coefs)[5])
      scale <- exp(lambda0 + lambda1*group)
      
    }else{
      if(is.null(pi_)){
        pi_ <- exp(unlist(cm$coefs)[1]) / (1 + exp(unlist(cm$coefs)[1]))
      }
      scale <- exp(unlist(cm$coefs)[2])
      shape <- exp(unlist(cm$coefs)[3])
    }
    
    f_t <- dweibull(tte,shape=shape,scale= scale^{-1/shape})
    S_t <- pweibull(tte,shape=shape,scale= scale^{-1/shape},lower.tail=FALSE)
    
    ## Formulae (see paper)
    s_ <- pi_ + (1 - pi_) * S_t
    ll <- sum ( evt*log((1-pi_) * f_t) + (1-evt)*log(s_),na.rm=T )
    ## Return object
    return(ll)
  }
  
  ll<-ll.mix(tte=data$tte,evt=data$evt,group=data$group,cm=cm)
  ll_0<-ll.mix(tte=data$tte,evt=data$evt,group=data$group,cm=cm,pi_ = 0)
  
  dn <- -2*(ll_0-ll)
  
  quantile <- qchisq((conf_level - 1/2 )*2 ,1)
  return( list(dn,quantile))
}

count_default <- function(link_function,scenario) {
  nb_default_pi <- 0
  nb_default_rn <- 0
  nb_default_pi_o <- 0
  nb_default_rn_o <- 0
  
  for (i in 1:10000) {
    data <- scenario$result[[i]]$data
    cm <- fit.cure.model(Surv(tte, evt) ~ group, formula.surv = list(~group, ~1), data = data, dist = "weibull", link = link_function)
    
    if (!is.null(cm$covariance) && cm$optim$convergence == 0) {
      results <- compute_ratio_test(data, cm)
      if (with(results, pi_a < 0.025 | pi_b < 0.025)) {
        nb_default_pi <- nb_default_pi + 1
      }
      
      if (with(results, rn_a > 0.05 | rn_b > 0.05)) {
        nb_default_rn <- nb_default_rn + 1
      }
      
      cm_overall <- fit.cure.model(Surv(tte,evt)~1, data=data, dist = "weibull",link=link_function)
      whichTau <- max(data$tte)
      beta<-unlist(cm_overall$coefs)[1]
      pi <- exp(beta) / (1 + exp(beta))
      scale <- exp(unlist(cm_overall$coefs)[2])
      shape <- exp(unlist(cm_overall$coefs)[3])
      r_n <- compute_rn(whichTau, shape, scale, pi)
      
      if (pi < 0.025) {
        nb_default_pi_o <- nb_default_pi_o + 1
      }
      if(r_n > 0.05){
        nb_default_rn_o <- nb_default_rn_o + 1
      }
    }
  }
  
  return(list(nb_default_pi = nb_default_pi, nb_default_rn = nb_default_rn, nb_default_pi_overall = nb_default_pi_o, nb_default_rn_overall = nb_default_rn_o))
}


run_count_default <- function(scenario) {
  logit_results <- count_default("logit",scenario)
  loglog_results <- count_default("loglog",scenario)
  probit_results <- count_default("probit",scenario)
  
  return(list(logit = logit_results, loglog = loglog_results, probit = probit_results))

}



# setwd("/home/ckouadio/Documents/")
# source("functions1.R")
# library(cuRe)
# scenario <- readRDS("scenario1.rds")
# result <- run_count_default(scenario)
# saveRDS(result, "count_default_1.rds")
# 
# 
# setwd("/home/ckouadio/Documents/")
# source("functions1.R")
# library(cuRe)
# scenario <- readRDS("scenario2.rds")
# result <- run_count_default(scenario)
# saveRDS(result, "count_default_2.rds")
# 
# setwd("/home/ckouadio/Documents/")
# source("functions1.R")
# library(cuRe)
# scenario <- readRDS("scenario3.rds")
# result <- run_count_default(scenario)
# saveRDS(result, "count_default_3.rds")
# 
# setwd("/home/ckouadio/Documents/")
# source("functions1.R")
# library(cuRe)
# scenario <- readRDS("scenario4.rds")
# result <- run_count_default(scenario)
# saveRDS(result, "count_default_4.rds")


